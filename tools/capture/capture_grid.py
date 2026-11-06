"""Capture of a 2+2 quadruple into the 3:2 mean-motion resonance of its two binaries: a grid of N-body runs against the
theory of Tremaine (2020, arXiv:2003.01141, Sect. 4 and Appendix C).

Two binaries of equal stars (all four 1 Msun) on near-circular orbits, the short one with P_S = 1 d, the long one with
P_L near 1.5 d, on a circular outer orbit coplanar with them; Newtonian gravity only. The period of the LONG binary is
made to shrink at a constant rate (orbits.A.Pdot, angular-momentum loss, swingersim.pdot), so that the ratio
P_L / P_S falls through 3:2 from above. In Tremaine's notation (his binary A = the short one, B = the long one) the
detuning Delta = 3 n_B / (2 n_A) - 1 then grows through zero, the direction in which capture into the resonance with
the argument

    phi = 3 lambda_L - 2 lambda_S - varpi_L        (first order in the eccentricity e_L of the long binary)

is possible. Tremaine's conditions (equal masses; epsilon = 2^(-1/3) (P_L / P_out)^(2/3), his eq. 3):

  * adiabatic passage, his eq. 24: 0 < dDelta/dt < c n_L epsilon^(20/3), with c of order unity (about 0.17 from his
    eqs. C8-C11 and the numerical f_tau = 2.75);
  * capture certain if the initial e_L is small, his eq. 23: e_L0 < epsilon^(5/3); otherwise probabilistic.

The grid varies epsilon (through P_out), the migration rate r = (dDelta/dt) / (n_L epsilon^(20/3)) and e_L0 / epsilon^(5/3),
with random phases. Every run starts 15 resonance widths w = epsilon^(10/3) / 4.645 (his eq. C11) above the
resonance and migrates over 30 widths. The position of the resonance in the configured (osculating) P_L is measured
first for every run, with its own phases, by two short runs without migration (the configured periods are
osculating; their offsets from the mean periods depend on the phases and are much larger than the resonance).

    python tools/capture/capture_grid.py output/capture/grid -j 0            # all runs, restartable
    python tools/capture/capture_grid.py output/capture/grid --summary       # only the table of captures

Writes rows.jsonl (one line per run: classification and a compressed history) and summary.csv.
A run counts as captured when, in the last quarter of the run, the resonant angle librates (its 4 sub-windows each
span less than 300 deg) and the period ratio has stopped following the migration.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

P_S = 1.0  # period of the short binary [d]
P_L = 1.5  # nominal period of the long binary [d]
EPS = (0.02, 0.03, 0.05)
RATES = (0.01, 0.03, 0.1, 0.3, 1.0)
E_FACTORS = (0.1, 1.0, 2.0, 4.0, 8.0)
SEEDS = 6
CONTROL = {"eps": 0.03, "rate": -0.03, "e_factor": 1.0}  # divergent passage: no capture expected
DAYS = 365.25


def p_out(eps: float) -> float:
    """Outer period [d] for Tremaine's epsilon (equal masses, his eq. 3)."""
    return P_L / (eps * 2.0 ** (1.0 / 3.0)) ** 1.5


def n_long() -> float:
    """Mean motion of the long binary [rad/yr]."""
    return 2.0 * math.pi * DAYS / P_L


def width(eps: float) -> float:
    """Natural width of the resonance in Delta (his eq. C11: delta of order unity)."""
    return eps ** (10.0 / 3.0) / 4.645


def config(p_l: float, eps: float, e_l: float, pdot: float, t_max: float, n_output: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    ang = rng.uniform(0.0, 360.0, 4)
    return {
        "topology": "2+2",
        "masses": {"m1": 1.0, "m2": 1.0, "m3": 1.0, "m4": 1.0},
        "orbits": {
            "A": {"P": p_l, "e": e_l, "i": 90.0, "omega": float(ang[0]), "Omega": 0.0, "M0": float(ang[1]),
                  "Pdot": pdot},  # fmt: skip
            "B": {"P": P_S, "e": 1e-5, "i": 90.0, "omega": 0.0, "Omega": 0.0, "M0": float(ang[2])},
            "out": {"P": p_out(eps), "e": 0.0, "i": 90.0, "omega": 0.0, "Omega": 0.0, "M0": float(ang[3])},
        },
        "physics": {"gr": False, "tides": False, "rotation": False},
        "run": {"backend": "nbody", "t_max_yr": t_max, "n_output": n_output, "extremes_every_orbits": 0},
    }


def _run(d: dict):
    from swingersim.backends.nbody import run_nbody
    from swingersim.config import Config

    return run_nbody(Config.model_validate(d))


def _series(res):
    o = res.orbits
    t = np.asarray(res.times, float)
    lam = {k: np.asarray(o[k]["M"]) + np.asarray(o[k]["omega"]) + np.asarray(o[k]["Omega"]) for k in "AB"}
    varpi = np.asarray(o["A"]["omega"]) + np.asarray(o["A"]["Omega"])
    phi1 = np.unwrap(np.mod(3.0 * lam["A"] - 2.0 * lam["B"], 2 * np.pi))
    phi = np.unwrap(np.mod(3.0 * lam["A"] - 2.0 * lam["B"] - varpi, 2 * np.pi))
    per = {k: np.sqrt(np.asarray(o[k]["a"]) ** 3 / 2.0) * DAYS for k in "AB"}  # 2 Msun per binary
    return t, phi1, phi, np.unwrap(np.mod(varpi, 2 * np.pi)), per["A"] / per["B"], np.asarray(o["A"]["e"])


def calibrate(eps: float, e_l: float, seed: int) -> dict:
    """P_L of the resonance for this epsilon and these initial phases (the configured, osculating periods differ from
    the mean ones by short-period terms of order 1e-4 that depend on the phases): phi = phi1 - varpi is stationary
    where the rate of phi1 equals that of varpi; the rate of phi1 falls by 3 * 2 pi DAYS / P_L^2 per day of P_L.
    Two runs of 60 yr without migration (e_L at least 0.002, so that varpi is defined)."""
    rows = []
    e_c = max(e_l, 0.002)
    for p in (P_L, P_L * 1.0002):
        res = _run(config(p, eps, e_c, 0.0, 60.0, 12001, seed))
        t, phi1, _, varpi, _, _ = _series(res)
        rows.append((p, float(np.polyfit(t, phi1, 1)[0]), float(np.polyfit(t, varpi, 1)[0])))
    (p1, r1, w1), (p2, r2, w2) = rows
    slope = (r2 - r1) / (p2 - p1)
    return {"P_L_res": p1 + (0.5 * (w1 + w2) - r1) / slope, "slope_rad_yr_per_d": slope,
            "varpi_rate_rad_yr": 0.5 * (w1 + w2)}  # fmt: skip


def plan(eps: float, rate: float, e_factor: float, seed: int, p_res: float) -> dict:
    """The migration run: start 15 resonance widths before the resonance in the direction of travel (Delta ~
    (3/2) P_S / P_L - 1 near it, so a width w in Delta is a relative change w of P_L) and migrate over 30 widths."""
    w = width(eps)
    speed = abs(rate) * n_long() * eps ** (20.0 / 3.0)  # |dDelta/dt| [1/yr]
    sgn = 1.0 if rate > 0 else -1.0
    p0 = p_res * (1.0 + sgn * 15.0 * w)
    t_max = 30.0 * w / speed
    pdot = -sgn * speed * p0  # dP_L/dt [d/yr]
    n_output = int(min(max(4.0 * t_max, 4001), 400001))
    return config(p0, eps, e_factor * eps ** (5.0 / 3.0), pdot, t_max, n_output, seed)


def run_one(task: dict) -> dict:
    t0 = time.perf_counter()
    row = {k: task[k] for k in ("eps", "rate", "e_factor", "seed")}
    try:
        e_l = task["e_factor"] * task["eps"] ** (5.0 / 3.0)
        cal = calibrate(task["eps"], e_l, 1000 + task["seed"])
        row["calibration"] = cal
        cfg = plan(task["eps"], task["rate"], task["e_factor"], 1000 + task["seed"], cal["P_L_res"])
        row["P_L0"] = cfg["orbits"]["A"]["P"]
        row["Pdot_d_yr"] = cfg["orbits"]["A"]["Pdot"]
        row["t_max_yr"] = cfg["run"]["t_max_yr"]
        res = _run(cfg)
        t, _, phi, _, ratio, e_l = _series(res)
        row["event"] = res.event
        row["max_energy_error"] = float(np.nanmax(np.abs(res.energy_error)))
        n = len(t)
        q = n // 4
        last = slice(n - q, n)
        spans = [float(np.degrees(np.ptp(phi[s : s + q // 4]))) for s in range(n - q, n - q // 4 + 1, q // 4)][:4]
        row["phi_spans_last_quarter_deg"] = spans
        # the ratio follows the migration unless it is held: compare its drift in the last quarter with the imposed
        drift = float(np.polyfit(t[last], ratio[last], 1)[0])
        row["ratio_drift_per_yr"] = drift
        row["imposed_ratio_drift_per_yr"] = cfg["orbits"]["A"]["Pdot"] / P_S
        row["e_L_start"] = float(np.mean(e_l[: max(q // 4, 1)]))
        row["e_L_end"] = float(np.mean(e_l[-max(q // 4, 1) :]))
        row["captured"] = bool(all(s < 300.0 for s in spans) and abs(drift) < 0.5 * abs(row["imposed_ratio_drift_per_yr"]))
        k = max(n // 200, 1)  # compressed history, 200 points
        row["history"] = {"t": t[::k].tolist(), "ratio": ratio[::k].tolist(), "e_L": e_l[::k].tolist(),
                          "phi_deg": np.degrees(np.mod(phi[::k], 2 * np.pi)).tolist()}  # fmt: skip
    except Exception as exc:  # noqa: BLE001
        row.update(event="error", error=f"{type(exc).__name__}: {exc}"[:300])
    row["wall_s"] = time.perf_counter() - t0
    return row


def tasks_for(eps_list, rates, e_factors, seeds, control: bool = True) -> list[dict]:
    combos = list(itertools.product(eps_list, rates, e_factors, range(seeds)))
    if control:
        combos += [(CONTROL["eps"], CONTROL["rate"], CONTROL["e_factor"], s) for s in range(seeds)]
    out = []
    for eps, rate, ef, seed in combos:
        speed = abs(rate) * n_long() * eps ** (20.0 / 3.0)
        out.append({"eps": eps, "rate": rate, "e_factor": ef, "seed": seed, "t_max_yr": 30.0 * width(eps) / speed})
    return out


def key(r: dict) -> tuple:
    return (round(r["eps"], 6), round(r["rate"], 6), round(r["e_factor"], 6), int(r["seed"]))


def summary(out: Path) -> pd.DataFrame:
    rows = [json.loads(x) for x in (out / "rows.jsonl").read_text().splitlines() if x.strip()]
    df = pd.DataFrame([{k: v for k, v in r.items() if k != "history"} for r in rows])
    ok = df[df.event != "error"]
    tab = ok.groupby(["eps", "rate", "e_factor"]).captured.agg(["mean", "count"]).reset_index()
    tab.to_csv(out / "summary.csv", index=False)
    return tab


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out", type=Path)
    ap.add_argument("-j", "--processes", type=int, default=0)
    ap.add_argument("--eps", type=float, nargs="+", default=list(EPS))
    ap.add_argument("--rates", type=float, nargs="+", default=list(RATES))
    ap.add_argument("--e-factors", type=float, nargs="+", default=list(E_FACTORS))
    ap.add_argument("--seeds", type=int, default=SEEDS)
    ap.add_argument("--no-control", action="store_true", help="leave out the divergent control runs")
    ap.add_argument("--summary", action="store_true", help="only write and print summary.csv")
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    if a.summary:
        print(summary(a.out).to_string(index=False))
        return 0

    from swingersim.ensemble import resolve_processes, row_line
    from swingersim.parallel import run_tasks

    tasks = tasks_for(a.eps, a.rates, a.e_factors, a.seeds, not a.no_control)
    rows_file = a.out / "rows.jsonl"
    done = set()
    if rows_file.exists():
        for line in rows_file.read_text().splitlines():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if r.get("event") != "error":
                done.add(key(r))
    todo = [t for t in tasks if key(t) not in done]
    # the longest runs first, so that the pool is not left waiting for them at the end
    todo.sort(key=lambda t: -t["t_max_yr"])
    nproc = resolve_processes(a.processes)
    years = sum(t["t_max_yr"] for t in todo)
    print(f"capture grid: {len(tasks)} runs, {len(tasks) - len(todo)} done, {len(todo)} to run "
          f"({years:.3g} yr in total) on {nproc} processes", flush=True)  # fmt: skip
    t0 = time.perf_counter()

    def lost(task, msg):
        return {k: task[k] for k in ("eps", "rate", "e_factor", "seed")} | {"event": "error", "error": msg}

    with open(rows_file, "a", encoding="utf-8") as f:
        for k, row in enumerate(run_tasks(run_one, todo, nproc, lost), 1):
            f.write(row_line(row))
            f.flush()
            print(f"  {k}/{len(todo)} eps {row['eps']:g} rate {row['rate']:g} e/eps^5/3 {row['e_factor']:g} "
                  f"seed {row['seed']}: {row.get('event') or ('CAPTURED' if row.get('captured') else 'passed')}"
                  f"  [{time.perf_counter() - t0:.0f} s]", flush=True)  # fmt: skip
    print(summary(a.out).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
