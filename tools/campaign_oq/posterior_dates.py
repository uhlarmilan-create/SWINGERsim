"""Campaign "open questions", study S4: eclipse dates of a pair from the posterior of an MCMC fit.

Every drawn posterior sample (stage2_samples.csv of `swingersim fit` with sampler emcee) is integrated over a span of
calendar years; its inclination i(t) is compared with the pair's own geometric eclipse limit (swingersim.incfit), and
the years at which the eclipses start or stop are recorded.

    python tools/campaign_oq/posterior_dates.py output/campaign_oq/s4_qx_cas_mcmc 300 \
        --years 1900 2150 --pair inner [--processes 0]

Writes posterior_dates.json (all samples: i(t), limit(t), switches) and posterior_dates.txt (16-50-84 % of the end of
the past eclipse window and of the next return) into the fit folder.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import squadra3_report as S
from _common import atomic_write_text, eclipse_switches, fit_epoch
from squadra3_posterior import cfg_from

_STATE: dict = {}


def _init(fit_dir: str, years: list, pair: str) -> None:
    fit_dir = Path(fit_dir)
    with open(fit_dir / "stage2_fitted.yaml", encoding="utf-8") as f:
        _STATE["base"] = yaml.safe_load(f)
    _STATE["samples"] = pd.read_csv(fit_dir / "stage2_samples.csv")
    _STATE["t0"] = fit_epoch(fit_dir / "stage2_fitted.yaml")
    _STATE["years"] = np.asarray(years, float)
    _STATE["pair"] = pair


def run(k: int) -> dict:
    row = _STATE["samples"].iloc[k].drop(labels=["log_prob"], errors="ignore")
    d = cfg_from(_STATE["base"], row)
    yrs = _STATE["years"]
    try:
        c = S.inc_curve(d, _STATE["t0"], yrs, pair=_STATE["pair"])
    except Exception as exc:  # noqa: BLE001 - a sample that cannot be integrated
        return {"k": int(k), "error": f"{type(exc).__name__}: {exc}"[:200]}
    return {
        "k": int(k),
        "i": c["i"],
        "limit": c["limit"],
        "changes": eclipse_switches(yrs, c["i"], c["limit"]),
    }


def summary(out: list, now: float) -> str:
    ok = [r for r in out if "error" not in r]
    ends, returns, never = [], [], 0
    for r in ok:
        past = [y for y in r["changes"] if y <= now]
        fut = [y for y in r["changes"] if y > now]
        if past:
            ends.append(past[-1])
        if fut:
            returns.append(fut[0])
        else:
            never += 1

    def q(x):
        return "n/a" if not x else " - ".join(f"{v:.1f}" for v in np.percentile(x, [16, 50, 84]))

    lines = [
        f"{len(ok)} posterior samples integrated ({len(out) - len(ok)} failed)",
        f"last eclipse switch before {now:.0f} (16-50-84 %): {q(ends)}",
        f"next eclipse switch after {now:.0f} (16-50-84 %): {q(returns)}",
        f"samples without a switch after {now:.0f} within the span: {never}",
    ]
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "fit_dir", type=Path, help="folder of the MCMC fit (stage2_fitted.yaml, stage2_samples.csv)"
    )
    ap.add_argument("n", type=int, help="number of posterior samples")
    ap.add_argument("--years", type=float, nargs=2, default=(1900.0, 2150.0))
    ap.add_argument("--step", type=float, default=0.5, help="output step [yr]")
    ap.add_argument("--pair", default="inner")
    ap.add_argument("--now", type=float, default=2026.0)
    ap.add_argument("--processes", type=int, default=0, help="worker processes (0: all cores)")
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args(argv)
    years = list(np.arange(a.years[0], a.years[1] + 1e-9, a.step))
    _init(str(a.fit_dir), years, a.pair)
    n_samples = len(_STATE["samples"])
    ks = np.random.default_rng(a.seed).choice(n_samples, min(a.n, n_samples), replace=False)
    t = time.time()
    with Pool(
        a.processes or os.cpu_count() or 1,
        initializer=_init,
        initargs=(str(a.fit_dir), years, a.pair),
    ) as p:
        out = p.map(run, [int(k) for k in ks])
    atomic_write_text(
        a.fit_dir / "posterior_dates.json", json.dumps({"years": years, "samples": out})
    )
    text = summary(out, a.now)
    atomic_write_text(a.fit_dir / "posterior_dates.txt", text + "\n")
    print(text)
    print(f"{len(out)} samples in {time.time() - t:.0f} s -> {a.fit_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
