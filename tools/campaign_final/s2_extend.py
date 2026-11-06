"""Final campaign, study S2 (continued): the stability criteria for a 2+2 with one widened binary.

The ensemble of the open-questions campaign (configs/campaign_oq/s2_stability_criteria_2p2.yaml, 1200 members,
800 yr = 200 outer orbits) is re-drawn member by member (the same seed gives the same members; every member is
checked against summary.csv of that ensemble) and

  ghost  every member: the N-body ghost test of Vynatheya et al. (2023, Sect. 3.1) over 100 outer orbits, the
         label the MLP classifier was trained on (stable / unbound / diverged);
  long   the members that stayed intact in 800 yr and were chaotic (MEGNO > 4) or below the Mardling-Aarseth
         limit (ratio < 1): integrated again for 10^4 outer orbits (3.98e4 yr), with MEGNO, to see which of them
         break up later.

    python tools/campaign_final/s2_extend.py ghost output/campaign_oq/s2_stability output/campaign_final/s2_ghost -j 0
    python tools/campaign_final/s2_extend.py long  output/campaign_oq/s2_stability output/campaign_final/s2_long  -j 0

Rows go to <out>/rows.jsonl (one line per member, appended as they finish: a restart skips the finished ones) and
<out>/summary.csv at the end.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

CONFIG = ROOT / "configs/campaign_oq/s2_stability_criteria_2p2.yaml"
N_OUTER_LONG = 1.0e4  # outer orbits of the long runs
N_OUTER_GHOST = 100.0  # outer orbits of the ghost test (Vynatheya et al. 2023)
MEGNO_CHAOTIC = 4.0


def _members(ens_dir: Path):
    """The members of the S2 ensemble, re-drawn and checked against its summary.csv."""
    from swingersim.config import load_config
    from swingersim.ensemble import draw_members

    cfg = load_config(CONFIG)
    members = draw_members(cfg)
    summ = pd.read_csv(ens_dir / "summary.csv").set_index("member")
    if len(members) != len(summ):
        raise SystemExit(f"{len(members)} members drawn, {len(summ)} in {ens_dir}/summary.csv")
    for m in members:
        row = summ.loc[m.index]
        for p, v in m.params.items():
            if not math.isclose(float(row[p]), float(v), rel_tol=1e-9, abs_tol=1e-9):
                raise SystemExit(f"member {m.index}: {p} = {v} drawn, {row[p]} in summary.csv")
    return cfg, members, summ


def _ghost(task: dict) -> dict:
    from swingersim.config import Config
    from swingersim.stability import ghost_test

    t0 = time.perf_counter()
    row = {"member": task["index"]}
    try:
        c = Config.model_validate(task["config"])
        c._base_dir = task["base_dir"]
        g = ghost_test(c, n_outer=N_OUTER_GHOST)
        row.update(
            ghost_stable=float(g.stable),
            ghost_reason=g.reason,
            ghost_t_end_yr=g.t_end_yr,
            ghost_outer_orbits=g.outer_orbits,
            ghost_delta_max=g.delta_max,
        )
    except Exception as exc:  # noqa: BLE001 - one failing member must not stop the rest
        row.update(ghost_reason="error", error=f"{type(exc).__name__}: {exc}"[:300])
    row["wall_s"] = time.perf_counter() - t0
    return row


def _long(task: dict) -> dict:
    from swingersim.ensemble import run_member

    return run_member(task)


def _lost(task: dict, message: str) -> dict:
    return {"member": task["index"], "fate": "error", "ghost_reason": "error", "error": message}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("mode", choices=("ghost", "long"))
    ap.add_argument("ensemble", type=Path, help="folder of the S2 ensemble (summary.csv)")
    ap.add_argument("out", type=Path)
    ap.add_argument("-j", "--processes", type=int, default=0, help="0: all cores")
    ap.add_argument("--limit", type=int, default=0, help="only the first N selected members (tests)")
    ap.add_argument("--n-outer", type=float, default=None, help="long: outer orbits (default 1e4)")
    a = ap.parse_args(argv)

    from swingersim.ensemble import resolve_processes, row_line, write_summary
    from swingersim.parallel import run_tasks
    from swingersim.units import sma_to_period

    cfg, members, summ = _members(a.ensemble)
    if a.mode == "ghost":
        sel = members
    else:
        s = summ
        pick = (s.fate == "intact") & ((s.megno > MEGNO_CHAOTIC) | (s.ma01_min < 1.0))
        sel = [m for m in members if bool(pick.loc[m.index])]
    if a.limit:
        sel = sel[: a.limit]
    a.out.mkdir(parents=True, exist_ok=True)
    rows_file = a.out / "rows.jsonl"
    done = set()
    if rows_file.exists():
        for line in rows_file.read_text().splitlines():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if r.get("fate") != "error" and r.get("ghost_reason") != "error":
                done.add(int(r["member"]))
    n_outer = a.n_outer or N_OUTER_LONG
    tasks = []
    for m in sel:
        if m.index in done:
            continue
        conf = json.loads(json.dumps(m.config, default=str))
        if a.mode == "long":
            from swingersim.config import Config
            from swingersim.hierarchy import build_tree

            c = Config.model_validate(conf)
            c._base_dir = cfg.base_dir
            root, masses = build_tree(c)
            p_out = sma_to_period(root.elements.a, float(np.sum(masses)))
            conf["run"]["t_max_yr"] = float(n_outer * p_out)
            conf["run"]["n_output"] = 2001
            conf["run"]["megno"] = True
        tasks.append(
            {
                "index": m.index,
                "direction": "forward",
                "config": conf,
                "params": m.params,
                "base_dir": cfg.base_dir,
                "mlp": False,
                "ghost": False,
                "ghost_max_wall_s": None,
                "ghost_orbits": N_OUTER_GHOST,
                "save_dir": "",
            }
        )
    print(
        f"S2 {a.mode}: {len(sel)} members selected, {len(sel) - len(tasks)} already done, "
        f"{len(tasks)} to run on {resolve_processes(a.processes)} processes",
        flush=True,
    )
    t0 = time.perf_counter()
    fn = _ghost if a.mode == "ghost" else _long
    with open(rows_file, "a", encoding="utf-8") as f:
        for k, row in enumerate(run_tasks(fn, tasks, resolve_processes(a.processes), _lost), 1):
            f.write(row_line(row))
            f.flush()
            what = row.get("ghost_reason") if a.mode == "ghost" else row.get("fate")
            print(
                f"  {k}/{len(tasks)} member {row['member']}: {what}  "
                f"[{time.perf_counter() - t0:.0f} s]",
                flush=True,
            )
    rows = []
    for line in rows_file.read_text().splitlines():
        try:
            rows.append(json.loads(line))
        except json.JSONDecodeError:
            pass
    last = {}
    for r in rows:  # the last line of a member wins (a re-run after an error)
        last[int(r["member"])] = r
    rows = [last[k] for k in sorted(last)]
    if a.mode == "long":
        write_summary(rows, a.out)
    else:
        pd.DataFrame(rows).to_csv(a.out / "summary.csv", index=False)
    bad = sum(1 for r in rows if r.get("fate") == "error" or r.get("ghost_reason") == "error")
    print(f"S2 {a.mode}: {len(rows)} rows, {bad} errors -> {a.out}", flush=True)
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
