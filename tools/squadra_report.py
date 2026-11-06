"""Numbers and figures of the SQUADRA test document (docs/<lang>/squadra_test.md).

Three published multiple stars - BU CMi (2+2), TIC 120362137 (3+1), V994 Her (2+2+2) - are
modelled by swingersim and compared with the published solutions, the published times of minima
and the TESS light curves.

    python tools/squadra_report.py --root .            # everything (needs the job outputs)
    python tools/squadra_report.py --root . --figures  # figures only, from the cached models

Inputs (under --root): configs/squadra_test/*, data/observed/<system>/*, output/squadra_test/
<system>/* (the jobs sq1-sq3). Model light curves and eclipse times of the three
configurations of every system (published, anchored, fitted) are cached in
output/squadra_test/<system>/report/ and recomputed when the configuration, the data or the
SWINGERsim version change (a .key file beside every cached CSV). The epoch (BJD of t = 0) of every
configuration is read from the configuration itself (tools/_common.fit_epoch). Figures go to
docs/figures/squadra_test/, the numbers to docs/figures/squadra_test/numbers.json.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import DPY, file_cached, fit_epoch, read_table

SYSTEMS = {
    "bu_cmi": {
        "name": "BU CMi",
        "pairs": ["A", "B"],
        "minima": "minima_pribulla2023.csv",
        "sectors": {
            "S7": ["tess_s07_tess-spoc.csv"],
            "S34": ["tess_s34_spoc.csv"],
            "S61": ["tess_s61_spoc.csv"],
        },
    },
    "tic120362137": {
        "name": "TIC 120362137",
        "pairs": ["inner"],
        "minima": "minima_borkovits2026.csv",
        "sectors": {
            "S14": ["tess_s14_tess-spoc.csv"],
            "S26": ["tess_s26_tess-spoc.csv"],
            "S40-41": ["tess_s40_tess-spoc.csv", "tess_s41_tess-spoc.csv"],
            "S53-54": ["tess_s53_tess-spoc.csv", "tess_s54_tess-spoc.csv"],
            "S74": ["tess_s74_tess-spoc.csv"],
        },
    },
    "v994_her": {
        "name": "V994 Her",
        "pairs": ["A", "B", "C"],
        "minima": "minima_zasche2023.csv",
        "sectors": {
            "S26": ["tess_s26_spoc.csv"],
            "S40": ["tess_s40_qlp.csv"],
            "S53": ["tess_s53_qlp.csv"],
        },
    },
}
MODELS = ("published", "anchored", "fitted")


def config_path(root: Path, sys_: str, model: str) -> Path:
    if model == "fitted":
        return root / "output" / "squadra_test" / sys_ / "fit_mcmc" / "stage2_fitted.yaml"
    return root / "configs" / "squadra_test" / f"{sys_}_{model}.yaml"


def model_epoch(root: Path, sys_: str, model: str) -> float:
    """BJD of t = 0 of a model configuration (published, anchored or fitted)."""
    return fit_epoch(config_path(root, sys_, model))


def tess_files(root: Path, sys_: str, files: list[str]) -> list[Path]:
    return [root / "data" / "observed" / sys_ / f for f in files]


def read_tess(root: Path, sys_: str, files: list[str]) -> pd.DataFrame:
    """TESS light curves of a sector (files may start with '#' comment lines)."""
    d = pd.concat([read_table(f) for f in tess_files(root, sys_, files)])
    d = d[np.isfinite(d.flux) & np.isfinite(d.bjd_tdb)].sort_values("bjd_tdb")
    return d.reset_index(drop=True)


# ------------------------------------------------------------------ cached model computations
def _run_swingersim(cmd: list[str], label: str) -> None:
    print(label, flush=True)
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        print(f"FAILED: {label}\n{exc.stderr}", file=sys.stderr, flush=True)
        raise


def lightcurve_jobs(root: Path, sys_: str):
    """(output CSV, input files, request, command builder) of every model light curve of a
    system: three models x the TESS sectors."""
    s = SYSTEMS[sys_]
    out = root / "output" / "squadra_test" / sys_ / "report"
    for sec, files in s["sectors"].items():
        d = read_tess(root, sys_, files)
        cad = float(np.median(np.diff(d.bjd_tdb))) * 1440.0
        cad = 2.0 if cad < 2.5 else (10.0 if cad < 15 else 30.0)
        for m in MODELS:
            cfg = config_path(root, sys_, m)
            t0 = model_epoch(root, sys_, m)
            start = float(d.bjd_tdb.min() - t0) - 0.2
            days = float(d.bjd_tdb.max() - d.bjd_tdb.min()) + 0.4
            args = [
                "--start-day",
                f"{start:.4f}",
                "--days",
                f"{days:.4f}",
                "--cadence-min",
                f"{min(cad, 10.0):g}",
                "--exposure-min",
                f"{cad:g}",
                "--supersample",
                str(max(1, int(cad // 2))),
                "--band",
                "TESS",
            ]
            yield (
                out / f"lc_{m}_{sec}.csv",
                [cfg, *tess_files(root, sys_, files)],
                {"t0": t0, "args": args},
                lambda tmp, cfg=cfg, args=args, sec=sec, m=m: _run_swingersim(
                    [
                        sys.executable,
                        "-m",
                        "swingersim",
                        "lightcurve",
                        str(cfg),
                        *args,
                        "-o",
                        str(tmp),
                    ]
                    + ["--no-html"],
                    f"lightcurve {sys_} {sec} {m}",
                ),
            )


def model_lightcurves(root: Path, systems, force=False, processes: int = 0) -> None:
    """Model light curves of all systems, in parallel subprocesses (SWINGERsim lightcurve)."""
    jobs = [j for s in systems for j in lightcurve_jobs(root, s)]
    if force:
        for f, *_ in jobs:
            f.with_name(f.name + ".key").unlink(missing_ok=True)
    n = processes or os.cpu_count() or 1

    def one(job):
        f, files, params, write = job
        return file_cached(f, write, files, params)

    with ThreadPoolExecutor(max_workers=n) as ex:  # each job is a subprocess
        list(ex.map(one, jobs))


def model_eclipses(root: Path, sys_: str, force=False) -> None:
    """All eclipses of the three configurations over the span of the measured minima."""
    from swingersim.config import load_config
    from swingersim.timing import eclipse_times

    s = SYSTEMS[sys_]
    mins = read_minima(root, sys_)
    for model in MODELS:
        f = root / "output" / "squadra_test" / sys_ / "report" / f"timings_{model}.csv"
        if force:
            f.with_name(f.name + ".key").unlink(missing_ok=True)
        cfg_path = config_path(root, sys_, model)
        t0 = model_epoch(root, sys_, model)
        y0 = (mins.bjd_tdb.min() - t0) / DPY - 0.05
        y1 = (mins.bjd_tdb.max() - t0) / DPY + 0.05

        def write(tmp, cfg_path=cfg_path, t0=t0, y0=y0, y1=y1, model=model):
            cfg = load_config(cfg_path)
            ecl = (eclipse_times(cfg, y0)[0] if y0 < 0 else []) + (
                eclipse_times(cfg, y1)[0] if y1 > 0 else []
            )
            rows = [
                {"binary": e.binary, "kind": e.kind, "bjd_obs": t0 + e.t_obs_yr * DPY} for e in ecl
            ]
            pd.DataFrame(rows, columns=["binary", "kind", "bjd_obs"]).to_csv(tmp, index=False)
            print("eclipses", sys_, model, len(rows), flush=True)

        minima_file = root / "data" / "observed" / sys_ / s["minima"]
        file_cached(f, write, [cfg_path, minima_file], {"t0": t0, "y0": y0, "y1": y1})


def read_minima(root: Path, sys_: str) -> pd.DataFrame:
    s = SYSTEMS[sys_]
    d = read_table(root / "data" / "observed" / sys_ / s["minima"])
    return d[np.isfinite(d.bjd_tdb)].reset_index(drop=True)


def model_times_table(root: Path, sys_: str, model: str) -> pd.DataFrame:
    return pd.read_csv(
        root / "output" / "squadra_test" / sys_ / "report" / f"timings_{model}.csv", comment="#"
    )


def match_minima(mins: pd.DataFrame, tab: pd.DataFrame) -> np.ndarray:
    """Model time of the nearest model eclipse of the same pair and type, for every minimum."""
    out = np.full(len(mins), np.nan)
    for (b, k), g in mins.groupby(["binary", "kind"]):
        t = np.sort(tab[(tab.binary == b) & (tab.kind == k)].bjd_obs.values)
        if len(t) == 0:
            continue
        j = np.clip(np.searchsorted(t, g.bjd_tdb.values), 1, len(t) - 1)
        near = np.where(
            np.abs(t[j - 1] - g.bjd_tdb.values) < np.abs(t[j] - g.bjd_tdb.values), t[j - 1], t[j]
        )
        # a minimum farther than half a period from every model eclipse has no partner
        half = 0.5 * np.median(np.diff(t)) if len(t) > 1 else np.inf
        near = np.where(np.abs(near - g.bjd_tdb.values) <= half, near, np.nan)
        out[g.index.values] = near
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--figures", action="store_true", help="figures only (models cached)")
    ap.add_argument("--force", action="store_true", help="recompute the cached models")
    ap.add_argument("--processes", type=int, default=0, help="parallel light curves (0: all cores)")
    ap.add_argument("--no-png", action="store_true", help="do not fail when a PNG is not written")
    ap.add_argument("--out", type=Path, default=None, help="figure folder")
    a = ap.parse_args(argv)
    if not a.figures:
        model_lightcurves(a.root, list(SYSTEMS), a.force, a.processes)
        for s in SYSTEMS:
            model_eclipses(a.root, s, a.force)
    from squadra_figures import make_all

    make_all(
        a.root,
        a.out or Path(__file__).resolve().parents[1] / "docs" / "figures" / "squadra_test",
        png=not a.no_png,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
