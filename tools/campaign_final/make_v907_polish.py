"""Final campaign: V907 Sco - the fit of Part III (fit v2) repeated from the best point of its MCMC.

The MCMC of the fit v2 (output/campaign/mcmc1000/v907_sco_fit_v2_mcmc_1000, 48 walkers x 1000 steps) has samples with
a chi^2 2.3 below the least-squares optimum (220.46 against 222.74, with the final jitters of the fit; chi^2 of a sample
is -2 log_prob minus the Jacobian term 2 ln e of the ecosw/esinw parametrisation). This script writes
configs/campaign_final/v907_sco_fit_v2_polished.yaml: configs/squadra_test3/v907_sco_fit_v2.yaml started at the MCMC
sample with the highest probability (no anchoring, which would move the start). Result of the run of 2026-11-06: the
fit, which re-estimates the jitters from the quoted errors, returned to the optimum of the fit v2 (chi^2 222.73, every
parameter within 0.03 sigma). The optimum does not depend on the start; Delta chi^2 = 2.3 is far inside the 1-sigma
region of 14 parameters (Delta chi^2 = 15.9). Part III uses this repeated fit.

    python tools/campaign_final/make_v907_polish.py [--samples <stage2_samples.csv>]
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "configs/squadra_test3/v907_sco_fit_v2.yaml"
SAMPLES = ROOT / "output/campaign/mcmc1000/v907_sco_fit_v2_mcmc_1000/stage2_samples.csv"
OUT = ROOT / "configs/campaign_final/v907_sco_fit_v2_polished.yaml"

HEADER = """# Final campaign: V907 Sco, fit v2 of Part III repeated from the best sample of its MCMC (chi^2 220.46 against the
# 222.74 of the least-squares optimum, with the final jitters). Written by
# tools/campaign_final/make_v907_polish.py from configs/squadra_test3/v907_sco_fit_v2.yaml: the same data, free
# parameters, errors and priors; start = sample with the largest log_prob of {samples}; anchoring off.
"""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--samples", type=Path, default=SAMPLES)
    ap.add_argument("-o", "--output", type=Path, default=OUT)
    a = ap.parse_args(argv)
    s = pd.read_csv(a.samples)
    best = s.loc[s.log_prob.idxmax()]
    d = yaml.safe_load(BASE.read_text(encoding="utf-8"))
    inner, out = d["orbits"]["inner"], d["orbits"]["out"]
    ec, es = float(best["orbits.inner.ecosw"]), float(best["orbits.inner.esinw"])
    inner.update(
        P=float(best["orbits.inner.P"]),
        M0=float(best["orbits.inner.M0"]) % 360.0,
        i=float(best["orbits.inner.i"]),
        e=math.hypot(ec, es),
        omega=math.degrees(math.atan2(es, ec)) % 360.0,
    )
    out.update(
        P=float(best["orbits.out.P"]),
        M0=float(best["orbits.out.M0"]) % 360.0,
        e=float(best["orbits.out.e"]),
        omega=float(best["orbits.out.omega"]) % 360.0,
        i=float(best["orbits.out.i"]),
        Omega=float(best["orbits.out.Omega"]) % 360.0,
    )
    d["masses"]["m3"] = float(best["masses.m3"])
    d["fit"]["rv"]["gamma"] = float(best["rv.gamma"])
    d["fit"]["rv"]["offsets"]["Lacy1999"] = float(best["rv.offset.Lacy1999"])
    d["fit"]["stage2"]["anchor"] = False
    d["fit"]["output"]["dir"] = "../../output/campaign_final/v907_sco_fit_v2_polished"
    a.output.parent.mkdir(parents=True, exist_ok=True)
    try:
        src = a.samples.resolve().relative_to(ROOT)
    except ValueError:
        src = a.samples
    a.output.write_text(
        HEADER.format(samples=src) + yaml.safe_dump(d, sort_keys=False, width=120), encoding="utf-8"
    )
    print(f"wrote {a.output} (start: sample with log_prob {best.log_prob:.3f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
