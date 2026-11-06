"""Final campaign, study S3 (continued): CzeV343 - a fine scan of the period of pair A across the 3:2 resonance.

The present system of configs/campaign_oq/s3_czev343_resonance.yaml (full physics) is integrated for 1000 yr for
every P_A of a grid fine enough to resolve a resonance of width ~1e-5 in the period ratio (Tremaine 2020). For every
run the three resonant angles

    phi_1 = 3 lambda_A - 2 lambda_B,   phi_2 = phi_1 - varpi_A,   phi_3 = phi_1 - varpi_B

(lambda = M + omega + Omega, varpi = omega + Omega, osculating elements in the frame of the configuration) are followed:
their net number of turns over the run, their mean rate and, if they do not circulate, the range of their libration.
In the open-questions run (P_A = 1.209341 d) the angles turned at 3.87, 3.71 and 3.21 rad/yr; they become stationary
near P_A = 1.21016, 1.21013 and 1.21002 d (d(3 n_A)/dP_A = -4707 rad/yr per day), all inside the grid below.

    python tools/campaign_final/s3_resonance_scan.py output/campaign_final/s3_resonance_scan -j 0
        [--p-range 1.20995 1.21030] [--step 2.5e-6] [--t-max 1000]

Writes rows.jsonl (one line per P_A, appended: a restart skips the finished ones), summary.csv and angles.npz (the three
angles every 1 yr of every run) into the output folder.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

CONFIG = ROOT / "configs/campaign_oq/s3_czev343_resonance.yaml"
NAMES = ("phi1", "phi2", "phi3")


def _angles(res) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    o = res.orbits
    lam = {k: o[k]["M"] + o[k]["omega"] + o[k]["Omega"] for k in ("A", "B")}
    pa = o["A"]["omega"] + o["A"]["Omega"]
    pb = o["B"]["omega"] + o["B"]["Omega"]
    p1 = 3.0 * lam["A"] - 2.0 * lam["B"]
    raw = {"phi1": p1, "phi2": p1 - pa, "phi3": p1 - pb}
    return np.asarray(res.times, float), {k: np.unwrap(np.mod(v, 2 * np.pi)) for k, v in raw.items()}


def run_one(task: dict) -> dict:
    from swingersim.backends.nbody import run_nbody
    from swingersim.config import Config

    t0 = time.perf_counter()
    row = {"P_A": task["P_A"]}
    try:
        c = Config.model_validate(task["config"])
        c._base_dir = task["base_dir"]
        res = run_nbody(c)
        t, ang = _angles(res)
        row["t_end_yr"] = float(t[-1])
        row["event"] = res.event
        row["max_energy_error"] = float(np.nanmax(np.abs(res.energy_error)))
        grid = np.arange(0.0, t[-1] + 1e-9, 1.0)
        series = {}
        for k, u in ang.items():
            turns = (u[-1] - u[0]) / (2 * np.pi)
            rate = float(np.polyfit(t, u, 1)[0])
            dev = u - u[0]
            row[f"{k}_turns"] = float(turns)
            row[f"{k}_rate_rad_yr"] = rate
            row[f"{k}_span_deg"] = float(np.degrees(dev.max() - dev.min()))
            # turning points of the angle averaged over 20 yr (short-period terms removed): a slow
            # circulation is monotonic, a libration turns back again and again
            s1 = np.interp(grid, t, u)
            dt = float(np.median(np.diff(t)))
            w = max(int(round(20.0 / dt)), 1)
            sm = np.convolve(u, np.ones(w) / w, mode="valid") if len(u) > 2 * w else u
            sm = sm[:: max(int(round(1.0 / dt)), 1)]  # one value per year
            d = np.diff(sm)
            d = d[np.abs(d) > 1e-12]
            rev = int(np.sum(np.sign(d[1:]) != np.sign(d[:-1]))) if len(d) > 1 else 0
            row[f"{k}_reversals"] = rev
            # librating: the angle never makes a full turn and turns back at least twice
            row[f"{k}_librates"] = bool(dev.max() - dev.min() < 2 * np.pi and rev >= 2)
            series[k] = s1.tolist()
        row["series"] = series
    except Exception as exc:  # noqa: BLE001
        row.update(event="error", error=f"{type(exc).__name__}: {exc}"[:300])
    row["wall_s"] = time.perf_counter() - t0
    return row


def _lost(task: dict, message: str) -> dict:
    return {"P_A": task["P_A"], "event": "error", "error": message}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out", type=Path)
    ap.add_argument("-j", "--processes", type=int, default=0)
    ap.add_argument("--p-range", type=float, nargs=2, default=(1.20995, 1.21030))
    ap.add_argument("--step", type=float, default=2.5e-6)
    ap.add_argument("--t-max", type=float, default=1000.0)
    ap.add_argument("--n-output", type=int, default=20001)
    a = ap.parse_args(argv)

    from swingersim.config import load_config
    from swingersim.ensemble import resolve_processes, row_line
    from swingersim.parallel import run_tasks

    cfg = load_config(CONFIG)
    base = json.loads(json.dumps(cfg.model_dump(), default=str))
    n = int(round((a.p_range[1] - a.p_range[0]) / a.step)) + 1
    grid = np.round(a.p_range[0] + a.step * np.arange(n), 10)
    a.out.mkdir(parents=True, exist_ok=True)
    rows_file = a.out / "rows.jsonl"
    done = set()
    if rows_file.exists():
        for line in rows_file.read_text().splitlines():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if r.get("event") != "error":
                done.add(round(float(r["P_A"]), 10))
    tasks = []
    for p in grid:
        if round(float(p), 10) in done:
            continue
        c = json.loads(json.dumps(base))
        c["orbits"]["A"]["P"] = float(p)
        c["orbits"]["A"]["a"] = None
        c["run"]["t_max_yr"] = float(a.t_max)
        c["run"]["n_output"] = int(a.n_output)
        tasks.append({"P_A": float(p), "config": c, "base_dir": cfg.base_dir})
    (a.out / "scan.yaml").write_text(
        yaml.safe_dump({"config": str(CONFIG), "p_range": list(a.p_range), "step": a.step,
                        "t_max_yr": a.t_max, "n_output": a.n_output, "points": int(n)})  # fmt: skip
    )
    nproc = resolve_processes(a.processes)
    print(f"S3 scan: {n} periods, {n - len(tasks)} done, {len(tasks)} to run on {nproc} processes", flush=True)
    t0 = time.perf_counter()
    with open(rows_file, "a", encoding="utf-8") as f:
        for k, row in enumerate(run_tasks(run_one, tasks, nproc, _lost), 1):
            f.write(row_line(row))
            f.flush()
            lib = [nm for nm in NAMES if row.get(f"{nm}_librates")]
            print(
                f"  {k}/{len(tasks)} P_A = {row['P_A']:.7f}: {row.get('event') or 'ok'}, "
                f"turns {', '.join(f'{row.get(nm + "_turns", float("nan")):+.2f}' for nm in NAMES)}"
                f"{', librating: ' + ' '.join(lib) if lib else ''}  [{time.perf_counter() - t0:.0f} s]",
                flush=True,
            )
    rows = {}
    for line in rows_file.read_text().splitlines():
        try:
            r = json.loads(line)
        except json.JSONDecodeError:
            continue
        rows[round(float(r["P_A"]), 10)] = r
    rows = [rows[k] for k in sorted(rows)]
    n_t = max((len(r.get("series", {}).get("phi1", [])) for r in rows), default=0)

    def table(nm: str) -> np.ndarray:  # one row per P_A, NaN where a run ended early
        out = np.full((len(rows), n_t), np.nan)
        for i, r in enumerate(rows):
            s = np.asarray(r.get("series", {}).get(nm, []), float)
            out[i, : len(s)] = s
        return out

    np.savez_compressed(
        a.out / "angles.npz",
        P_A=np.array([r["P_A"] for r in rows]),
        t_yr=np.arange(n_t, dtype=float),
        **{nm: table(nm) for nm in NAMES},
    )
    pd.DataFrame([{k: v for k, v in r.items() if k != "series"} for r in rows]).to_csv(
        a.out / "summary.csv", index=False
    )
    bad = sum(1 for r in rows if r.get("event") == "error")
    print(f"S3 scan: {len(rows)} rows, {bad} errors -> {a.out}", flush=True)
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
