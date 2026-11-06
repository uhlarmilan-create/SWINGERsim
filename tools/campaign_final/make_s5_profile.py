"""Final campaign, study S5 (continued): BU CMi - profile of chi^2 in the apsidal constant of pair A.

In the open-questions fit (configs/campaign_oq/s5_bu_cmi_k2_free.yaml, output/campaign_oq/s5_bu_cmi_k2) the apsidal
constant of m1 (pair A) went to its lower bound 0.0003 (stellar models: 0.0042), where the linear error means nothing.
Here k2 of m1 is fixed at the values of a grid and everything else of that fit is fitted again (k2 of m3 free as
there). chi^2 against k2 is the profile; its minimum and where it rises by 1 and 4 give the 1- and 2-sigma range.

chi^2 of different fits can only be compared with the same errors. The fit of S5 re-estimated the extra error of every
source (jitter) and the red-noise factor of the light curves from its own residuals; a profile with those estimates
redone in every fit would compare chi^2 with different denominators. So the errors are frozen: this script writes the
data the S5 fit kept (177 minima, 305 velocities, 4147 light-curve points; its 2 outlying minima left out) with the
final errors of that fit (quoted error with the jitter of its source, light curves times the red-noise factor) to
output/campaign_final/s5_profile/data/, and the profile fits use them as they are (error_scale none, no clipping, no
red-noise factor, one fit). At k2 = 0.0003 the profile fit must give back the chi^2 of the S5 fit (a check).

    python tools/campaign_final/make_s5_profile.py      # data + configs/campaign_final/s5_k2A_<value>.yaml
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "configs/campaign_oq/s5_bu_cmi_k2_free.yaml"
S5 = ROOT / "output/campaign_oq/s5_bu_cmi_k2"
DATA = ROOT / "output/campaign_final/s5_profile/data"
GRID = (0.0003, 0.001, 0.002, 0.003, 0.0042, 0.006, 0.009)

HEADER = """# Final campaign, study S5 (profile): BU CMi with the apsidal constant of m1 (pair A) FIXED at {k2} (stellar models:
# 0.0042); k2 of m3 free; all other parameters, priors and start values as configs/campaign_oq/s5_bu_cmi_k2_free.yaml.
# Data: the points kept by that fit with its final errors frozen (output/campaign_final/s5_profile/data, written by
# tools/campaign_final/make_s5_profile.py), fitted once without re-estimating the errors, so that chi^2 of the
# profile fits can be compared.
"""


def write_data(s5: Path, out: Path) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    m = pd.read_csv(s5 / "stage2_residuals.csv")
    m = m[m.status == "fitted"]
    kinds = {"p": "primary", "s": "secondary"}
    pd.DataFrame(
        {
            "binary": m.pair,
            "kind": m["type"].map(kinds),
            "bjd_tdb": m.bjd_tdb.map(lambda v: f"{v:.6f}"),
            "sigma_d": m.sigma_d.map(lambda v: f"{v:.8f}"),
            "ref": m.source,
        }
    ).to_csv(out / "minima.csv", index=False)
    r = pd.read_csv(s5 / "stage2_rv.csv")
    r = r[r.status == "fitted"]
    pd.DataFrame(
        {"star": r.star, "time": r.bjd_tdb, "rv": r.rv, "sigma": r.sigma, "source": r.source}
    ).to_csv(out / "rv.csv", index=False)
    lc = pd.read_csv(s5 / "stage2_lc.csv")
    lc = lc[lc.status == "fitted"]
    pd.DataFrame(
        {"time": lc.bjd_tdb.map(lambda v: f"{v:.7f}"), "flux": lc.flux, "sigma": lc.sigma, "segment": lc.segment}
    ).to_csv(out / "lc.csv", index=False)
    if m["type"].map(kinds).isna().any():
        raise SystemExit("unknown minimum type in stage2_residuals.csv")
    return {"minima": len(m), "rv": len(r), "lc": len(lc)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--grid", type=float, nargs="+", default=list(GRID))
    ap.add_argument("--s5", type=Path, default=S5, help="output folder of the S5 fit")
    ap.add_argument("-d", "--config-dir", type=Path, default=ROOT / "configs/campaign_final")
    a = ap.parse_args(argv)
    if a.config_dir.resolve().parent.parent != ROOT.resolve():
        raise SystemExit("the configurations must lie two levels below the repository (relative paths)")
    n = write_data(a.s5, DATA)
    print(f"frozen data: {n['minima']} minima, {n['rv']} velocities, {n['lc']} light-curve points -> {DATA}")
    raw = yaml.safe_load(BASE.read_text(encoding="utf-8"))
    rel = "../../output/campaign_final/s5_profile/data"
    a.config_dir.mkdir(parents=True, exist_ok=True)
    for k2 in a.grid:
        d = yaml.safe_load(yaml.safe_dump(raw))
        d["stars"]["m1"]["apsidal_k2"] = float(k2)
        f = d["fit"]
        f["minima"] = f"{rel}/minima.csv"
        f["rv"]["file"] = f"{rel}/rv.csv"
        lc = f["lc"][0]
        lc.update(file=f"{rel}/lc.csv", red_noise=False)
        s2 = f["stage2"]
        s2["free"] = [p for p in s2["free"] if p != "stars.m1.apsidal_k2"]
        s2.get("bounds", {}).pop("stars.m1.apsidal_k2", None)
        # the frozen errors include the jitter: the selection by the quoted error (max_sigma_d) was made by S5
        s2.update(error_scale="none", clip_sigma=None, refits=1, max_sigma_d=None)
        s2.pop("error_model", None)
        tag = f"{k2:.4f}"
        f["output"]["dir"] = f"../../output/campaign_final/s5_profile/k2A_{tag}"
        out = a.config_dir / f"s5_k2A_{tag}.yaml"
        out.write_text(HEADER.format(k2=tag) + yaml.safe_dump(d, sort_keys=False, width=120), encoding="utf-8")
        print(out.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
