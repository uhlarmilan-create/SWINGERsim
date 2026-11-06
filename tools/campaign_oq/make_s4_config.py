"""Campaign "open questions", study S4: write the configuration of the QX Cas MCMC.

The fit v3 of Part III (configs/squadra_test3/qx_cas_fit_v3.yaml: 96 minima 1933-1974 and the inclination history with
the one-sided limits of inclination_limits_v3.csv) is restarted at its converged solution (stage2_fitted.yaml of the
earlier run, if present on this computer) and sampled with emcee. The posterior then goes to posterior_dates.py for the
range of the dates at which the eclipses of QX Cas ended and when they return.

    python tools/campaign_oq/make_s4_config.py [--start output/special/qx_cas_fit_v3/stage2_fitted.yaml]
        [--walkers 24] [--steps 1500] [--tau-target 0] [--max-steps 5000] [--iterations 5]
        [--output-dir ../../output/campaign_oq/s4_qx_cas_mcmc] [-o configs/campaign_oq/s4_qx_cas_mcmc.yaml]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "configs/squadra_test3/qx_cas_fit_v3.yaml"
START = ROOT / "output/special/qx_cas_fit_v3/stage2_fitted.yaml"

HEADER = """# Campaign "open questions", study S4: QX Cas - MCMC of the dynamical fit with the inclination history as data.
# Written by tools/campaign_oq/make_s4_config.py from configs/squadra_test3/qx_cas_fit_v3.yaml{start}.
# Question (Part III, open questions): a random ensemble weighted afterwards keeps only one or two systems consistent
# with the whole inclination history, so it gives no range of dates. Here the fit itself is sampled: emcee with
# {walkers} walkers x {steps} steps{tau} (first half discarded); then tools/campaign_oq/posterior_dates.py integrates
# posterior samples over 1900-2150 for the end of the eclipses and their return.
"""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--start", type=Path, default=START)
    ap.add_argument("--walkers", type=int, default=24)
    ap.add_argument("--steps", type=int, default=1500)
    ap.add_argument(
        "--iterations", type=int, default=5, help="optimiser iterations before the MCMC"
    )
    ap.add_argument(
        "--tau-target",
        type=float,
        default=0.0,
        help="extend the chain until it is this many autocorrelation times long (0: exactly --steps)",
    )
    ap.add_argument("--max-steps", type=int, default=5000, help="--tau-target: the longest chain")
    ap.add_argument("--output-dir", default="../../output/campaign_oq/s4_qx_cas_mcmc")
    ap.add_argument(
        "-o", "--output", type=Path, default=ROOT / "configs/campaign_oq/s4_qx_cas_mcmc.yaml"
    )
    a = ap.parse_args(argv)
    raw = yaml.safe_load(BASE.read_text(encoding="utf-8"))
    note = ""
    if a.start.exists():
        fit = yaml.safe_load(a.start.read_text(encoding="utf-8"))
        for k in ("masses", "orbits", "stars"):
            raw[k] = fit[k]
        note = (
            f", started at {a.start.relative_to(ROOT) if a.start.is_relative_to(ROOT) else a.start}"
        )
    else:
        # without the earlier result the optimiser converges first (as fit v3: 40 iterations)
        a.iterations = max(a.iterations, 40)
        note = " (no earlier solution found: the optimiser converges first)"
    s2 = raw["fit"]["stage2"]
    s2.update(sampler="emcee", walkers=a.walkers, mcmc_steps=a.steps, max_iterations=a.iterations)
    if a.tau_target > 0:
        s2.update(mcmc_tau_target=a.tau_target, mcmc_max_steps=a.max_steps)
    raw["fit"]["output"] = {"dir": a.output_dir, "predict": {"years": 50, "samples": 100}}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(
        HEADER.format(
            start=note,
            walkers=a.walkers,
            steps=a.steps,
            tau=(
                f", extended until {a.tau_target:g} autocorrelation times (at most {a.max_steps})"
                if a.tau_target > 0
                else ""
            ),
        )
        + yaml.safe_dump(raw, sort_keys=False, width=120),
        encoding="utf-8",
    )
    print(f"wrote {a.output}{note}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
