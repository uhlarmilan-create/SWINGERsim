"""V907 Sco: integrate posterior samples of the MCMC fit (output/squadra_test3/v907_sco/fit_v2_mcmc) over 1890-2070:
inclination curve, eclipse switches and nodal period of every sample, for tools/squadra3_report.py.

    python tools/squadra3_posterior.py 100 output/squadra_test3/v907_sco/fit_v2_mcmc/posterior_inc.json
        [--fit-dir output/squadra_test3/v907_sco/fit_v2_mcmc] [--processes N] [--seed 7]

The epoch (BJD of t = 0) is read from the fitted configuration; the samples run on all cores (--processes). Every
sample is still integrated twice (swingersim.incfit.model_inclination for i(t) over 1890-2070 and the axis track
2021-2087 of squadra3_report.axis_track), as before.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import squadra3_report as S
from _common import atomic_write_text, eclipse_switches, fit_epoch, nodal_period

YEARS = np.arange(1890.0, 2071.0, 0.5)
AXIS_YEARS = np.arange(2021.0, 2087.0, 0.5)
_STATE: dict = {}  # base configuration, samples and epoch of a worker


def cfg_from(base: dict, row) -> dict:
    """The configuration of one posterior sample: the fitted configuration with the sampled parameters (orbits.*,
    masses.*, stars.<star>.<key>; ecosw / esinw give e and omega; rv.* are not part of the dynamics)."""
    d = copy.deepcopy(base)
    for k, v in row.items():
        p = k.split(".")
        if p[0] in ("rv", "log_prob"):  # not parameters of the dynamics
            continue
        if p[0] == "masses":
            d["masses"][p[1]] = float(v)
            continue
        if p[0] == "stars":
            d.setdefault("stars", {}).setdefault(p[1], {})[p[2]] = float(v)
            continue
        if p[0] != "orbits" or len(p) != 3:
            raise KeyError(f"posterior column {k!r} is not a parameter of the configuration")
        o = d["orbits"][p[1]]
        if p[2] in ("ecosw", "esinw"):
            o["_" + p[2]] = float(v)
        else:
            o[p[2]] = float(v)
    for o in d["orbits"].values():
        if "_ecosw" in o:
            ec, es = o.pop("_ecosw"), o.pop("_esinw")
            o["e"] = math.hypot(ec, es)
            o["omega"] = math.degrees(math.atan2(es, ec)) % 360
    return d


def _init(fit_dir: str) -> None:
    fit_dir = Path(fit_dir)
    with open(fit_dir / "stage2_fitted.yaml", encoding="utf-8") as f:
        _STATE["base"] = yaml.safe_load(f)
    _STATE["samples"] = pd.read_csv(fit_dir / "stage2_samples.csv")
    _STATE["t0"] = fit_epoch(fit_dir / "stage2_fitted.yaml")


def run(k: int) -> dict:
    d = cfg_from(_STATE["base"], _STATE["samples"].iloc[k])
    t0 = _STATE["t0"]
    c = S.inc_curve(d, t0, YEARS)
    tr = S.axis_track(d, t0, AXIS_YEARS)
    return {
        "k": int(k),
        "i": c["i"],
        "limit": c["limit"],
        "Pnod": nodal_period(AXIS_YEARS, tr["h"], np.array(tr["L"])[0]),
        "changes": eclipse_switches(YEARS, c["i"], c["limit"]),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("n", type=int, help="number of posterior samples")
    ap.add_argument("output", type=Path)
    ap.add_argument(
        "--fit-dir", type=Path, default=HERE.parent / "output/squadra_test3/v907_sco/fit_v2_mcmc"
    )
    ap.add_argument("--processes", type=int, default=0, help="worker processes (0: all cores)")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args(argv)
    _init(str(a.fit_dir))
    n_samples = len(_STATE["samples"])
    ks = np.random.default_rng(a.seed).choice(n_samples, min(a.n, n_samples), replace=False)
    t = time.time()
    with Pool(
        a.processes or os.cpu_count() or 1, initializer=_init, initargs=(str(a.fit_dir),)
    ) as p:
        out = p.map(run, [int(k) for k in ks])
    atomic_write_text(a.output, json.dumps(out))
    print(f"{len(out)} samples in {time.time() - t:.0f} s -> {a.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
