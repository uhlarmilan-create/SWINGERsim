"""Cut a light curve to windows around the model eclipses and bin it, for fit.lc of swingersim fit.

    python tools/make_lc_windows.py CONFIG BJD_T0 OUT.csv LABEL=TESS.csv [LABEL=TESS.csv ...]
        [--half-window 0.3] [--bin-min 10] [--pairs A B]

Every input (columns bjd_tdb or time, flux, flux_err or sigma) becomes one segment named LABEL;
the eclipse times come from swingersim timings of CONFIG (t = 0 at BJD_T0; data before BJD_T0 use
the eclipses integrated backward), so the configuration should already follow the minima
(anchored). The flux is divided by its median.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from swingersim import units
from swingersim.config import load_config
from swingersim.timing import eclipse_times


def read(path: Path) -> pd.DataFrame:
    d = pd.read_csv(path, comment="#")
    d.columns = [c.strip().lower() for c in d.columns]
    t = d["bjd_tdb"] if "bjd_tdb" in d else d["time"]
    s = d["flux_err"] if "flux_err" in d else d["sigma"]
    out = pd.DataFrame({"time": t, "flux": d["flux"], "sigma": s}).dropna()
    return out.sort_values("time")


def binned(d: pd.DataFrame, minutes: float) -> pd.DataFrame:
    if minutes <= 0:
        return d
    w = minutes / 1440.0
    k = np.floor((d["time"] - d["time"].min()) / w).astype(int)
    g = d.groupby(k)
    out = pd.DataFrame(
        {
            "time": g["time"].mean(),
            "flux": g["flux"].mean(),
            "sigma": g["sigma"].apply(lambda s: np.sqrt(np.sum(s**2)) / len(s)),
            "n": g.size(),
        }
    )
    return out.drop(
        columns="n"
    )  # a bin of one point is kept (cadence >= the bin, e.g. 10-min FFIs)


def model_eclipses(cfg, bjd_t0: float, b0: float, b1: float, pairs=None) -> np.ndarray:
    """BJD of the model eclipses between b0 and b1 (sorted): integrated backward from t = 0 for
    the part before bjd_t0 (swingersim.timing.eclipse_times supports negative times)."""
    y0 = (b0 - bjd_t0) / units.DAYS_PER_YEAR - 0.01
    y1 = (b1 - bjd_t0) / units.DAYS_PER_YEAR + 0.01
    ecl = []
    if y0 < 0:
        ecl += eclipse_times(cfg, y0, pairs)[0]
    if y1 > 0:
        ecl += eclipse_times(cfg, y1, pairs)[0]
    tm = np.sort([bjd_t0 + e.t_obs_yr * units.DAYS_PER_YEAR for e in ecl])
    return tm[(tm > b0) & (tm < b1)]


def distance_to_nearest(t: np.ndarray, tm: np.ndarray) -> np.ndarray:
    """|t - nearest of the sorted times tm| for every t (inf without tm): binary search instead
    of an N x M distance matrix (a TESS sector against thousands of eclipses)."""
    t = np.asarray(t, float)
    if len(tm) == 0:
        return np.full(len(t), np.inf)
    if len(tm) == 1:
        return np.abs(t - tm[0])
    j = np.clip(np.searchsorted(tm, t), 1, len(tm) - 1)
    return np.minimum(np.abs(t - tm[j - 1]), np.abs(t - tm[j]))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    ap.add_argument("bjd_t0", type=float)
    ap.add_argument("output", type=Path)
    ap.add_argument("inputs", nargs="+", help="LABEL=file.csv")
    ap.add_argument("--half-window", type=float, default=0.3, help="[d] around each eclipse")
    ap.add_argument("--bin-min", type=float, default=10.0, help="binning [min], 0 = none")
    ap.add_argument("--pairs", nargs="+", help="pairs whose eclipses are kept (default all)")
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    rows = []
    for item in a.inputs:
        label, _, path = item.partition("=")
        d = read(Path(path))
        med = float(np.median(d["flux"]))
        d["flux"] = d["flux"] / med
        d["sigma"] = d["sigma"] / med
        t0, t1 = float(d["time"].min()), float(d["time"].max())
        tm = model_eclipses(cfg, a.bjd_t0, t0 - a.half_window, t1 + a.half_window, a.pairs)
        near = distance_to_nearest(d["time"].to_numpy(), tm)
        w = binned(d[near <= a.half_window], a.bin_min)
        w["segment"] = label
        rows.append(w)
        print(f"{label}: {len(tm)} eclipses, {len(w)} points")
    out = pd.concat(rows)
    with open(a.output, "w") as f:
        f.write(
            f"# light curve windows of +-{a.half_window} d around the model eclipses of "
            f"{a.config.name}, binned to {a.bin_min:g} min (tools/make_lc_windows.py)\n"
        )
        out[["time", "flux", "sigma", "segment"]].to_csv(f, index=False, float_format="%.7f")
    print(f"{len(out)} points in {a.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
