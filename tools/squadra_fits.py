"""Figures and numbers of the extended mini-fit (jobs sq1-sq3 9-14) for the
SQUADRA test document (docs/<lang>/squadra_test.md, section "The extended fit").

    python tools/squadra_fits.py --root . [--output output] [--out docs/figures/squadra_test]

Inputs: output/squadra_test/<system>/<variant>/stage2_report.txt, stage2_params.csv,
stage2_lc.csv, stage2_rv.csv (written by `swingersim fit` stage 2). Numbers go to
docs/figures/squadra_test/numbers_fits.json.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import C_DATA, FigureWriter, atomic_write_text, mpl_pyplot

plt = mpl_pyplot()
C_FIT = "#1F4FB0"
C_PUB = "#C2571A"

SYSTEMS = {
    "bu_cmi": {"name": "BU CMi", "pairs": ["A", "B"]},
    "tic120362137": {"name": "TIC 120362137", "pairs": ["inner"]},
    "v994_her": {"name": "V994 Her", "pairs": ["A", "B", "C"]},
}
VARIANTS = ["fit", "fit_robust", "fit_priors", "fit_auto", "fit_rv", "fit_lc", "fit_lc_roche"]
VLABEL = {
    "fit": "first fit (0.30)",
    "fit_robust": "+ errors/outliers",
    "fit_priors": "+ priors",
    "fit_auto": "auto (from published)",
    "fit_rv": "+ RV",
    "fit_lc": "+ light curve",
    "fit_lc_roche": "+ LC, Roche",
}
# published values (+- sigma) of the orbits that the minima constrain only weakly
PUBLISHED = {
    "bu_cmi": {"orbits.out.e": (0.2728, 0.0037), "orbits.out.omega": (141.0, 3.6)},
    "tic120362137": {"orbits.out.e": (0.274, 0.013), "orbits.out.omega": (289.0, 4.8)},
    "v994_her": {"orbits.AB.e": (0.687, 0.044), "orbits.AB.omega": (59.6, 3.8)},
}


FIGURES: FigureWriter | None = None  # set by make_all


def _save(fig, out: Path, name: str) -> str:
    w = FIGURES if FIGURES is not None and FIGURES.out == Path(out) else FigureWriter(out)
    return w.save_mpl(fig, name)


# ------------------------------------------------------------------ parsing
NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"  # a number, also -1.2e-05


def _section(text: str, title: str) -> list[str]:
    """The indented lines after a section title of the report (empty when it is missing)."""
    if title not in text:
        return []
    return text.split(title, 1)[1].splitlines()[1:]


def parse_report(path: Path) -> dict:
    """The numbers of a stage2_report.txt: chi2, rms per pair, RV and LC rms, jitter (minima
    [min], velocities [km/s], light curves [ppt]), outliers. Sections that a fit does not have
    (no minima, no velocities, no light curves) are left empty."""
    t = path.read_text()
    r: dict = {"rms": {}, "rv_rms": {}, "lc": {}, "jitter": {}, "n_outliers": 0}
    m = re.search(rf"chi\^2 ({NUM}) -> ({NUM}) \(reduced ({NUM})", t)
    if m:
        r["chi2_start"], r["chi2"], r["chi2_red"] = map(float, m.groups())
    m = re.search(r"(\d+) minima", t)
    r["n_minima"] = int(m.group(1)) if m else None
    for line in _section(t, "rms of the residuals [min], start -> fitted:"):
        m = re.match(rf"\s+(\w+)\s+(\d+) minima\s+({NUM}) ->\s+({NUM})", line)
        if not m:
            break
        r["rms"][m.group(1)] = {
            "n": int(m.group(2)),
            "start": float(m.group(3)),
            "fit": float(m.group(4)),
        }
    for line in _section(t, "radial velocities, rms of the residuals [km/s], start -> fitted:"):
        m = re.match(rf"\s+(m\d+)\s+(\d+) velocities\s+({NUM}) ->\s+({NUM})", line)
        if not m:
            break
        r["rv_rms"][m.group(1)] = {"n": int(m.group(2)), "fit": float(m.group(4))}
    lc_title = "light curves, rms of measured - calibrated model [ppt], start -> fitted, scale a:"
    for line in _section(t, lc_title):
        m = re.match(rf"\s+(\S+)\s+(\d+) points\s+({NUM}) ->\s+({NUM})\s+a = ({NUM})", line)
        if not m:
            break
        r["lc"][m.group(1)] = {
            "n": int(m.group(2)),
            "start": float(m.group(3)),
            "fit": float(m.group(4)),
            "a": float(m.group(5)),
        }
    m = re.search(r"jitter added in quadrature to the quoted errors: (.*)", t)
    if m:
        for part in m.group(1).split(";"):
            mm = re.match(rf"\s*(.+?) ({NUM}) (min|km/s|ppt)\s*$", part)
            if mm:
                r["jitter"][mm.group(1).strip()] = float(mm.group(2))
    m = re.search(r"(\d+) outlying minima left out", t)
    r["n_outliers"] = int(m.group(1)) if m else 0
    return r


def read_params(path: Path) -> dict:
    if not path.exists():
        return {}
    d = pd.read_csv(path, comment="#")
    return {
        row.parameter: {
            "value": float(row.value),
            "error": float(row.error),
            "start": float(row.start),
        }
        for row in d.itertuples()
    }


# ------------------------------------------------------------------ figures
def fig_minima_rms(numbers: dict, out: Path) -> str:
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.3), gridspec_kw={"width_ratios": [2, 1.2, 3]})
    for ax, (s, info) in zip(axes, SYSTEMS.items(), strict=True):
        vs = [v for v in VARIANTS if v in numbers[s]]
        pairs = info["pairs"]
        w = 0.8 / len(vs)
        cmap = plt.get_cmap("Blues")
        for k, v in enumerate(vs):
            vals = [numbers[s][v]["report"]["rms"].get(p, {}).get("fit", np.nan) for p in pairs]
            x = np.arange(len(pairs)) + (k - (len(vs) - 1) / 2) * w
            ax.bar(
                x, vals, w * 0.92, color=cmap(0.3 + 0.7 * k / max(len(vs) - 1, 1)), label=VLABEL[v]
            )
        ax.set_xticks(
            np.arange(len(pairs)), [f"pair {p}" if p != "inner" else "Aa–Ab" for p in pairs]
        )
        ax.set_title(info["name"])
        ax.set_ylabel("rms of measured − model [min]")
        ax.grid(axis="x", visible=False)
    axes[-1].legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)
    fig.tight_layout()
    return _save(fig, out, "fits_minima_rms")


def fig_weak_elements(numbers: dict, out: Path) -> str:
    fig, axes = plt.subplots(3, 2, figsize=(8.6, 6.4))
    for row, (s, info) in enumerate(SYSTEMS.items()):
        vs = [v for v in VARIANTS if v in numbers[s]]
        for col, (par, (pv, ps)) in enumerate(PUBLISHED[s].items()):
            ax = axes[row, col]
            ax.axhspan(pv - ps, pv + ps, color=C_PUB, alpha=0.15, lw=0)
            ax.axhline(pv, color=C_PUB, lw=1, ls=(0, (5, 3)), label="published ± 1σ")
            for k, v in enumerate(vs):
                p = numbers[s][v]["params"].get(par)
                if p is None:
                    continue
                ax.errorbar(k, p["value"], yerr=p["error"], fmt="o", ms=4.5, color=C_FIT, capsize=2)
            ax.set_xticks(
                range(len(vs)), [VLABEL[v] for v in vs], rotation=35, ha="right", fontsize=7
            )
            orbit = par.split(".")[1]
            name = {"e": "e", "omega": "ω [°]"}[par.split(".")[2]]
            ax.set_title(f"{info['name']}: {name} of the {orbit} orbit", fontsize=8.5)
            if row == 0 and col == 0:
                ax.legend(loc="best", frameon=False)
    fig.tight_layout()
    return _save(fig, out, "fits_weak_elements")


def _windows(t: np.ndarray, gap: float = 0.08) -> list[np.ndarray]:
    br = np.where(np.diff(t) > gap)[0]
    return np.split(np.arange(len(t)), br + 1)


def fig_lc(root_out: Path, s: str, variant: str, out: Path, n_win: int = 4) -> str:
    d = pd.read_csv(root_out / s / variant / "stage2_lc.csv", comment="#")
    segs = list(dict.fromkeys(d.segment))
    picks = []  # deepest windows spread over the segments
    for seg in segs:
        g = d[d.segment == seg].sort_values("bjd_tdb").reset_index(drop=True)
        wins = [g.iloc[ix] for ix in _windows(g.bjd_tdb.values) if len(ix) > 15]
        wins.sort(key=lambda w: w.flux.min())
        picks += wins[: max(1, n_win // len(segs))]
    picks = picks[:n_win]
    fig, axes = plt.subplots(
        2,
        len(picks),
        figsize=(2.6 * len(picks), 3.8),
        sharey="row",
        gridspec_kw={"height_ratios": [3, 1.3]},
    )
    axes = np.atleast_2d(axes).reshape(2, -1)
    for k, w in enumerate(picks):
        tm = w.bjd_tdb.values
        x = (tm - tm[np.argmin(w.model.values)]) * 24
        a, b = axes[0, k], axes[1, k]
        a.plot(x, w.flux, ".", ms=2.5, color=C_DATA, label="TESS")
        a.plot(x, w.model_start, color=C_PUB, lw=0.9, ls=(0, (5, 3)), label="start")
        a.plot(x, w.model, color=C_FIT, lw=1.1, label="fit")
        b.plot(x, 1e3 * w.o_minus_model, ".", ms=2.5, color=C_DATA)
        b.axhline(0, color=C_FIT, lw=0.8)
        a.set_title(f"{w.segment.iloc[0]}  BJD {tm.mean():.1f}", fontsize=7.5)
        b.set_xlabel("hours from minimum")
    axes[0, 0].set_ylabel("relative flux")
    axes[1, 0].set_ylabel("O − model [ppt]")
    axes[0, 0].legend(loc="lower left", frameon=False)
    fig.tight_layout()
    return _save(fig, out, f"{s}_{variant}")


STAR_ORBIT = {
    "bu_cmi": {"m1": "A", "m2": "A", "m3": "B", "m4": "B"},
    "tic120362137": {"m1": "inner", "m2": "inner", "m3": "mid", "m4": "out"},
    "v994_her": {"m1": "A", "m2": "A", "m3": "B", "m4": "B"},
}


def fig_rv(root_out: Path, s: str, variant: str, out: Path) -> str:
    """Radial velocities against the phase of the orbit that dominates them, with the model."""
    d = pd.read_csv(root_out / s / variant / "stage2_rv.csv", comment="#")
    par = read_params(root_out / s / variant / "stage2_params.csv")
    stars = sorted(d.star.unique())
    fig, axes = plt.subplots(
        2,
        len(stars),
        figsize=(2.6 * len(stars), 4.0),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1.3]},
    )
    axes = np.atleast_2d(axes).reshape(2, -1)
    srcs = list(dict.fromkeys(d.source))
    cols = ["#1F4FB0", "#C2571A", "#2E7D4F", "#7A4FB0", "#8E1F1F"]
    for k, st in enumerate(stars):
        orb = STAR_ORBIT[s][st]
        P = par.get(f"orbits.{orb}.P", {}).get("value")
        g = d[d.star == st]
        ph = ((g.bjd_tdb - g.bjd_tdb.min()) / P) % 1.0
        a, b = axes[0, k], axes[1, k]
        fit = g.status == "fitted"
        a.plot(ph[fit], g.model[fit], "_", ms=5, color="#9AA5B5", label="model")
        for j, src in enumerate(srcs):
            h = fit & (g.source == src)
            a.plot(ph[h], g.rv[h], "o", ms=2.5, color=cols[j % 5], label=src)
            b.plot(ph[h], g.o_minus_model[h], "o", ms=2.5, color=cols[j % 5])
        b.plot(ph[~fit], g.o_minus_model[~fit], "x", ms=4, color="k")
        rms = np.sqrt(np.mean(g[fit].o_minus_model ** 2))
        a.set_title(f"{st} (orbit {orb}, P = {P:.3f} d): rms {rms:.1f} km/s", fontsize=7.5)
        b.axhline(0, color="k", lw=0.6)
        b.set_xlabel("phase of the orbit")
    axes[0, 0].set_ylabel("RV [km/s]")
    axes[1, 0].set_ylabel("O − model [km/s]")
    axes[0, 0].legend(loc="best", frameon=False, fontsize=6.5)
    fig.tight_layout()
    return _save(fig, out, f"{s}_{variant}_rv")


def make_all(root: Path, root_out: Path, out: Path, png: bool = True) -> dict:
    global FIGURES
    FIGURES = FigureWriter(out, png=png)
    numbers: dict = {}
    for s in SYSTEMS:
        numbers[s] = {}
        for v in VARIANTS:
            rp = root_out / s / v / "stage2_report.txt"
            if not rp.exists():
                continue
            numbers[s][v] = {
                "report": parse_report(rp),
                "params": read_params(root_out / s / v / "stage2_params.csv"),
            }
    figs = [fig_minima_rms(numbers, out), fig_weak_elements(numbers, out)]
    for s in SYSTEMS:
        if (root_out / s / "fit_lc" / "stage2_lc.csv").exists():
            figs.append(fig_lc(root_out, s, "fit_lc", out))
        v = "fit_lc" if (root_out / s / "fit_lc" / "stage2_rv.csv").exists() else "fit_rv"
        if (root_out / s / v / "stage2_rv.csv").exists():
            figs.append(fig_rv(root_out, s, v, out))
    atomic_write_text(out / "numbers_fits.json", json.dumps(numbers, indent=1))
    print("\n".join(figs))
    FIGURES.finish()
    return numbers


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument(
        "--output",
        type=Path,
        default=None,
        help="folder with the job outputs (default root/output)",
    )
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--no-png", action="store_true", help="do not fail when a PNG is not written")
    a = ap.parse_args(argv)
    root_out = (a.output or a.root / "output") / "squadra_test"
    make_all(
        a.root, root_out, a.out or a.root / "docs" / "figures" / "squadra_test", png=not a.no_png
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
