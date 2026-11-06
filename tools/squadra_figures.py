"""Figures and numbers of the SQUADRA test document (called by tools/squadra_report.py)."""

from __future__ import annotations

import csv
import itertools
import json
import re
from pathlib import Path

import h5py
import matplotlib.ticker
import numpy as np
import pandas as pd
import yaml
from _common import (
    C_DATA,
    FigureWriter,
    atomic_write_text,
    ensemble_observed,
    ensemble_rows,
    ensemble_weights,
    mpl_pyplot,
    nested_value,
    read_table,
    wrap_deg,
)
from _common import year as YEAR
from squadra_report import (
    MODELS,
    SYSTEMS,
    config_path,
    match_minima,
    model_epoch,
    model_times_table,
    read_minima,
    read_tess,
)

from swingersim.ensemble import effective_size

plt = mpl_pyplot()

# colours: data dark grey; published orange dashed; anchored grey-blue dotted; fitted blue solid
C_MODEL = {"published": "#C2571A", "anchored": "#7C8AA0", "fitted": "#1F4FB0"}
LS_MODEL = {"published": (0, (5, 3)), "anchored": (0, (1, 2)), "fitted": "-"}
LABEL = {
    "published": "published solution",
    "anchored": "anchored (P, M0 of the pairs)",
    "fitted": "SWINGERsim fit",
}
MI_COLS = [
    "#1F4FB0",
    "#C2571A",
    "#2E7D4F",
    "#7A4FB0",
    "#8E1F1F",
    "#3A3F45",
    "#B08A1F",
    "#4FA0B0",
    "#B04F8A",
    "#6B7480",
]
FIGURES: FigureWriter | None = None  # set by make_all


def _save(fig, out: Path, name: str) -> str:
    w = FIGURES if FIGURES is not None and FIGURES.out == Path(out) else FigureWriter(out)
    return w.save_mpl(fig, name)


# ------------------------------------------------------------------ light curves
def _pair_stars(sys_: str) -> dict[str, list[str]]:
    if sys_ == "tic120362137":
        return {
            "Aa+Ab (3.3-d pair)": ["flux_m1", "flux_m2"],
            "B (51-d star)": ["flux_m3"],
            "C (1046-d star)": ["flux_m4"],
        }
    if sys_ == "v994_her":
        return {
            "A": ["flux_m1", "flux_m2"],
            "B": ["flux_m3", "flux_m4"],
            "C": ["flux_m5", "flux_m6"],
        }
    return {"A": ["flux_m1", "flux_m2"], "B": ["flux_m3", "flux_m4"]}


def _design(t: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Columns: the model flux (scale a) and, for every continuous segment of the sector (gaps
    longer than 0.5 d split it), a constant, a slope and a curvature: the extra light and the slow
    instrumental trends of the TESS photometry. These are the only adjusted numbers; the eclipses
    themselves come from the model unchanged."""
    cols = [f]
    cuts = np.concatenate([[0], np.nonzero(np.diff(t) > 0.5)[0] + 1, [len(t)]])
    for a, b in itertools.pairwise(cuts):
        seg = np.zeros_like(t)
        seg[a:b] = 1.0
        x = np.zeros_like(t)
        x[a:b] = (t[a:b] - t[a:b].mean()) / max(t[b - 1] - t[a], 1e-3)
        cols += [seg, seg * x, seg * x**2]
    return np.vstack(cols).T


def lc_compare(root: Path, sys_: str, sec: str, model: str, coef=None):
    """Observed TESS flux and the model at the observed times, calibrated linearly (scale and a
    smooth per-segment baseline); with coef the calibration of another model is used."""
    s = SYSTEMS[sys_]
    d = read_tess(root, sys_, s["sectors"][sec])
    m = pd.read_csv(
        root / "output" / "squadra_test" / sys_ / "report" / f"lc_{model}_{sec}.csv", comment="#"
    )
    tm = model_epoch(root, sys_, model) + m.t_day.values
    f = np.interp(d.bjd_tdb.values, tm, m.flux.values)
    A = _design(d.bjd_tdb.values, f)
    if coef is None:
        w = 1.0 / np.maximum(d.flux_err.values, 1e-5)
        coef, *_ = np.linalg.lstsq(A * w[:, None], d.flux.values * w, rcond=None)
    fit = A @ coef
    res = d.flux.values - fit
    mad = 1.4826 * np.median(np.abs(res - np.median(res)))
    ok = np.abs(res) < 5 * mad
    base = A[:, 1:] @ coef[1:]
    return {
        "t": d.bjd_tdb.values,
        "obs": d.flux.values,
        "model": fit,
        "res": res,
        "coef": coef,
        "scale": float(coef[0]),
        "rms_ppt": float(np.sqrt(np.mean(res**2)) * 1e3),
        "rms_clip_ppt": float(np.sqrt(np.mean(res[ok] ** 2)) * 1e3),
        "extra_light": float(np.median(base) / (coef[0] + np.median(base))),
        "noise_ppt": float(np.median(d.flux_err.values) * 1e3),
        "mtab": m,
        "tm": tm,
    }


def fig_lightcurves(root: Path, sys_: str, out: Path, numbers: dict) -> str:
    s = SYSTEMS[sys_]
    secs = list(s["sectors"])
    fig, axes = plt.subplots(
        2 * len(secs),
        2,
        figsize=(11, 2.7 * len(secs)),
        gridspec_kw={"height_ratios": [3, 1.1] * len(secs), "width_ratios": [2.1, 1]},
    )
    nums = {}
    for k, sec in enumerate(secs):
        r = {m: lc_compare(root, sys_, sec, m) for m in MODELS}
        nums[sec] = {
            m: {
                kk: r[m][kk]
                for kk in ("rms_ppt", "rms_clip_ppt", "extra_light", "scale", "noise_ppt")
            }
            for m in MODELS
        }
        nums[sec]["n_points"] = len(r["fitted"]["t"])
        # the plot shows every model with the calibration of the fitted one
        shown = {m: lc_compare(root, sys_, sec, m, coef=r["fitted"]["coef"]) for m in MODELS}
        f = r["fitted"]
        t = f["t"] - 2450000
        ax, axr = axes[2 * k, 0], axes[2 * k + 1, 0]
        ax.plot(t, f["obs"], ".", ms=1.2, color=C_DATA, alpha=0.5, label="TESS", rasterized=True)
        gi = np.nonzero(np.diff(t) > 0.3)[0] + 1
        ax.plot(
            np.insert(t, gi, np.nan),
            np.insert(f["model"], gi, np.nan),
            "-",
            lw=0.6,
            color=C_MODEL["fitted"],
            label=LABEL["fitted"],
        )
        ax.set_title(f"{s['name']} – TESS sector {sec}: whole sector", loc="left")
        ax.set_ylabel("relative flux")
        axr.plot(t, f["res"] * 1e3, ".", ms=1.0, color=C_DATA, alpha=0.5, rasterized=True)
        axr.axhline(0, color=C_MODEL["fitted"], lw=0.8)
        axr.set_ylabel("O − C [ppt]")
        axr.set_xlabel("BJD − 2450000")
        lim = max(4 * f["rms_clip_ppt"], 2)
        axr.set_ylim(-lim, lim)
        # zoom: 2.5 days around the deepest model minimum in the middle of the sector
        mid = len(t) // 2
        j0 = np.argmin(f["model"][max(0, mid - len(t) // 6) : mid + len(t) // 6]) + max(
            0, mid - len(t) // 6
        )
        zc = t[j0]
        z = (t > zc - 1.25) & (t < zc + 1.25)
        az, azr = axes[2 * k, 1], axes[2 * k + 1, 1]
        hz = (t[z] - zc) * 24.0
        az.plot(hz, f["obs"][z], ".", ms=2.0, color=C_DATA, alpha=0.6)
        for m in ("published", "fitted"):
            zg = np.nonzero(np.diff(hz) > 7.0)[0] + 1
            az.plot(
                np.insert(hz, zg, np.nan),
                np.insert(shown[m]["model"][z], zg, np.nan),
                color=C_MODEL[m],
                ls=LS_MODEL[m],
                lw=1.1,
                label=LABEL[m],
            )
            azr.plot(hz, shown[m]["res"][z] * 1e3, ".", ms=1.6, color=C_MODEL[m], alpha=0.7)
        azr.axhline(0, color="#999999", lw=0.6)
        az.set_title(f"zoom around BJD {zc + 2450000:.2f}", loc="left")
        azr.set_xlabel("hours from the middle of the zoom")
        zl = max(np.nanmax(np.abs(shown["published"]["res"][z])) * 1e3 * 1.1, lim)
        azr.set_ylim(-zl, zl)
        if k == 0:
            ax.legend(loc="lower left", ncol=2, frameon=False)
            az.legend(loc="lower left", frameon=False)
    numbers.setdefault(sys_, {})["lightcurves"] = nums
    fig.tight_layout(h_pad=0.6)
    return _save(fig, out, f"{sys_}_lightcurves")


def fig_lc_components(root: Path, sys_: str, out: Path, sec: str, span: float = 4.0) -> str:
    """Light of the system and of each pair (or star) separately, from the fitted model."""
    s = SYSTEMS[sys_]
    m = pd.read_csv(
        root / "output" / "squadra_test" / sys_ / "report" / f"lc_fitted_{sec}.csv", comment="#"
    )
    t = model_epoch(root, sys_, "fitted") + m.t_day.values - 2450000
    t_start = t[0] + 2.0
    z = (t > t_start) & (t < t_start + span)
    d = read_tess(root, sys_, s["sectors"][sec])
    td = d.bjd_tdb.values - 2450000
    zd = (td > t_start) & (td < t_start + span)
    groups = _pair_stars(sys_)
    fig, axes = plt.subplots(
        len(groups) + 1, 1, figsize=(9, 1.6 * (len(groups) + 1) + 0.6), sharex=True
    )
    lc = lc_compare(root, sys_, sec, "fitted")
    axes[0].plot(td[zd], d.flux.values[zd], ".", ms=2.2, color=C_DATA, alpha=0.6, label="TESS")
    axes[0].plot(
        td[zd],
        lc["model"][zd],
        "-",
        color=C_MODEL["fitted"],
        lw=1.1,
        label="SWINGERsim: whole system",
    )
    axes[0].set_ylabel("system")
    axes[0].legend(loc="lower left", ncol=2, frameon=False)
    axes[0].set_title(
        f"{s['name']} – sector {sec}: the light of the whole system and of its parts (SWINGERsim)",
        loc="left",
    )
    cols = ["#1F4FB0", "#C2571A", "#2E7D4F"]
    for k, (g, cs) in enumerate(groups.items()):
        fl = m[cs].sum(axis=1).values
        axes[k + 1].plot(t[z], fl[z], "-", color=cols[k % 3], lw=1.1)
        axes[k + 1].set_ylabel(g, fontsize=8)
        share = np.nanmax(fl)
        axes[k + 1].text(
            0.99,
            0.08,
            f"share of the light {share * 100:.0f} %",
            transform=axes[k + 1].transAxes,
            ha="right",
            fontsize=7.5,
            color="#3D4651",
        )
    axes[-1].set_xlabel("BJD − 2450000")
    fig.tight_layout(h_pad=0.4)
    return _save(fig, out, f"{sys_}_components")


# ------------------------------------------------------------------ eclipse timings
def timing_residuals(root: Path, sys_: str) -> pd.DataFrame:
    mins = read_minima(root, sys_)
    for m in MODELS:
        tab = model_times_table(root, sys_, m)
        mins[f"t_{m}"] = match_minima(mins, tab)
        mins[f"r_{m}"] = (mins.bjd_tdb - mins[f"t_{m}"]) * 1440.0
    return mins


def fig_oc(root: Path, sys_: str, out: Path, numbers: dict) -> str:
    s = SYSTEMS[sys_]
    mins = timing_residuals(root, sys_)
    fitted = model_times_table(root, sys_, "fitted")
    pairs = s["pairs"]
    views = [(p, None) for p in pairs]
    if sys_ == "tic120362137":  # one season, where the 51-d variation can be seen
        views.append(("inner", (2021.30, 2021.80)))
    fig, axes = plt.subplots(
        2,
        len(views),
        figsize=(4.6 * len(views) + 0.6, 6.2),
        squeeze=False,
        gridspec_kw={"height_ratios": [2.2, 1.6]},
    )
    nums = {}
    for k, (p, xl) in enumerate(views):
        g = mins[mins.binary == p]
        ftab = fitted[(fitted.binary == p)].sort_values("bjd_obs")
        pri = ftab[ftab.kind == "primary"].bjd_obs.values
        P0 = np.median(np.diff(pri))
        n = np.concatenate([[0.0], np.cumsum(np.round(np.diff(pri) / P0))])
        P, T0 = np.polyfit(n, pri, 1)
        oc = lambda t, kind, T0=T0, P=P: (
            (t - (T0 + P * (np.round((t - T0) / P - 0.5 * kind) + 0.5 * kind))) * 1440.0
        )
        ax, axr = axes[0, k], axes[1, k]
        for m in ("anchored", "fitted"):
            tab = model_times_table(root, sys_, m)
            tab = tab[tab.binary == p].sort_values("bjd_obs")
            for kind, ls_extra in (("primary", 1.0), ("secondary", 0.6)):
                tt = tab[tab.kind == kind].bjd_obs.values
                if len(tt) == 0:
                    continue
                stp = max(1, len(tt) // 1500)
                tt = tt[::stp]
                ax.plot(
                    YEAR(tt),
                    oc(tt, 0 if kind == "primary" else 1),
                    color=C_MODEL[m],
                    ls=LS_MODEL[m],
                    lw=0.9 * ls_extra + 0.3,
                    alpha=0.95 if m == "fitted" else 0.8,
                    label=f"{LABEL[m]} ({kind[:3]}.)",
                )
        if sys_ == "tic120362137":  # the published ETV model of the paper (Source Data of Fig. 2)
            etv = read_table(root / "data" / "observed" / sys_ / "etv_model_borkovits2026.csv")
            for kind in ("primary", "secondary"):
                e = etv[etv.kind == kind]
                ax.plot(
                    YEAR(e.bjd_tdb),
                    oc(e.bjd_tdb.values, 0 if kind == "primary" else 1),
                    color="#8E1F1F",
                    lw=0.6,
                    alpha=0.6,
                    label="Borkovits et al. (2026) ETV model" if kind == "primary" else None,
                )
        for kind, mk in (("primary", "o"), ("secondary", "s")):
            gg = g[g.kind == kind]
            ax.errorbar(
                YEAR(gg.bjd_tdb),
                oc(gg.bjd_tdb.values, 0 if kind == "primary" else 1),
                yerr=gg.sigma_d * 1440,
                fmt=mk,
                ms=3.0,
                lw=0.5,
                color=C_DATA,
                mfc=C_DATA if kind == "primary" else "white",
                label=f"measured {kind}",
                zorder=5,
            )
        ax.set_title(f"pair {p}: O − C against T = {T0:.5f} + {P:.7f} E", loc="left", fontsize=8.5)
        if xl is not None:
            ax.set_xlim(xl)
            axr.set_xlim(xl)
            ax.set_title(
                f"pair {p}: one season ({xl[0]:.1f}–{xl[1]:.1f})", loc="left", fontsize=8.5
            )
        ax.set_ylabel("O − C [min]")
        for m in MODELS:
            axr.plot(
                YEAR(g.bjd_tdb),
                g[f"r_{m}"],
                ".",
                ms=3,
                color=C_MODEL[m],
                label=LABEL[m],
                alpha=0.85,
            )
        axr.set_yscale("symlog", linthresh=10)
        axr.axhline(0, color="#999999", lw=0.6)
        axr.set_ylabel("measured − model [min]")
        axr.set_xlabel("year")
        if xl is not None:
            continue
        nums[p] = {
            "n": len(g),
            **{f"rms_{m}": float(np.sqrt(np.nanmean(g[f"r_{m}"] ** 2))) for m in MODELS},
            **{f"median_abs_{m}": float(np.nanmedian(np.abs(g[f"r_{m}"]))) for m in MODELS},
            "T0": float(T0),
            "P": float(P),
        }
    h, lab = axes[0, 0].get_legend_handles_labels()
    h2, lab2 = axes[1, 0].get_legend_handles_labels()
    fig.legend(
        h + h2,
        lab + [f"residual: {x}" for x in lab2],
        loc="lower center",
        ncol=3,
        frameon=False,
        fontsize=7.5,
    )
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    # residuals of the fitted model per source
    src = {}
    for ref, g in mins.groupby("ref"):
        src[str(ref)] = {
            "n": len(g),
            "years": f"{int(YEAR(g.bjd_tdb.min()))}–{int(YEAR(g.bjd_tdb.max()))}",
            "median_sigma_s": float(np.median(g.sigma_d) * 86400),
            "rms_fitted_min": float(np.sqrt(np.nanmean(g.r_fitted**2))),
            "mean_fitted_min": float(np.nanmean(g.r_fitted)),
        }
    numbers.setdefault(sys_, {})["timings"] = nums
    numbers[sys_]["timing_sources"] = src
    return _save(fig, out, f"{sys_}_oc")


# ------------------------------------------------------------------ evolution
def fig_evolution(root: Path, sys_: str, out: Path, numbers: dict) -> str:
    s = SYSTEMS[sys_]
    base = root / "output" / "squadra_test" / sys_
    runs = []
    for d in ("backward", "forward"):
        with h5py.File(base / f"published_{d}.h5") as f:
            t = f["time_yr"][:]
            orb = {
                k: {q: f["orbits"][k][q][:] for q in ("e", "inc", "omega", "Omega", "a")}
                for k in f["orbits"]
            }
            mi = {k: f["mutual_inclination"][k][:] for k in f["mutual_inclination"]}
            runs.append((t, orb, mi))
    orbits = list(runs[0][1])
    show = [o for o in orbits if o != "out" or sys_ != "v994_her"]
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 6.4))
    cols = {
        o: c
        for o, c in zip(show, ["#1F4FB0", "#C2571A", "#2E7D4F", "#7A4FB0", "#3A3F45"], strict=False)
    }
    for t, orb, mi in runs:
        o_ = np.argsort(t)
        for o in show:
            lab = o if t[-1] > 0 else None
            axes[0, 0].plot(
                t[o_], orb[o]["e"][o_], color=cols[o], lw=0.35, label=lab, rasterized=True
            )
            w = np.degrees(np.unwrap(orb[o]["omega"][o_] + orb[o]["Omega"][o_]))
            w = w - w[np.argmin(np.abs(t[o_]))]
            axes[0, 1].plot(t[o_], w, color=cols[o], lw=0.8, label=lab)
            axes[1, 1].plot(t[o_], np.degrees(orb[o]["inc"][o_]), color=cols[o], lw=0.8, label=lab)
        for k, (key, v) in enumerate(sorted(mi.items())):
            axes[1, 0].plot(
                t[o_],
                np.degrees(v[o_]),
                lw=0.8,
                color=MI_COLS[k % len(MI_COLS)],
                label=key if t[-1] > 0 else None,
            )
    for ax in axes.flat:
        ax.axvline(0, color="#999999", lw=0.6)
        ax.set_xlabel("time from the epoch of the published solution [yr]")
    axes[0, 0].set_ylabel("eccentricity e")
    axes[0, 1].set_ylabel("longitude of periastron ϖ − ϖ(0) [deg]")
    axes[1, 0].set_ylabel("mutual inclination [deg]")
    axes[1, 1].set_ylabel("inclination to the sky [deg]")
    for ax in axes.flat:
        leg = ax.legend(fontsize=7, frameon=False, ncol=2)
        for line in leg.get_lines():
            line.set_linewidth(1.6)
    axes[0, 0].set_title(
        f"{s['name']}: published solution integrated backward and forward (N-body, GR, tides)",
        loc="left",
    )
    fig.tight_layout()
    return _save(fig, out, f"{sys_}_evolution")


# ------------------------------------------------------------------ ensemble
ENS_DIR = {"bu_cmi": "ensemble", "tic120362137": "ensemble2", "v994_her": "ensemble"}


def ensemble_config(root: Path, sys_: str) -> Path:
    """The configuration of the ensemble run: the copy in its output folder (what was run),
    else configs/squadra_test/<system>_ensemble.yaml."""
    run = root / "output" / "squadra_test" / sys_ / ENS_DIR[sys_] / "ensemble.yaml"
    return run if run.exists() else root / "configs" / "squadra_test" / f"{sys_}_ensemble.yaml"


def fig_ensemble(root: Path, sys_: str, out: Path, numbers: dict) -> str:
    s = SYSTEMS[sys_]
    summary = root / "output" / "squadra_test" / sys_ / ENS_DIR[sys_] / "summary.csv"
    d = pd.read_csv(summary, comment="#")
    fw = d[d.direction == "forward"].set_index("member")
    # the measured values ({value, sigma}) of the ensemble configuration; its hard conditions
    # (eclipsing, intact) are not applied here: the histograms show the measured quantities
    gauss = ensemble_observed(ensemble_config(root, sys_), gaussian_only=True)
    obs = {c: (o["value"], o["sigma"]) for c, o in gauss.items()}
    w = ensemble_weights(ensemble_rows(summary), gauss, fw.index.values)
    neff = effective_size(w)
    w = w / w.sum() if w.sum() > 0 else w
    fig, axes = plt.subplots(1, len(obs), figsize=(4.0 * len(obs), 3.0), squeeze=False)
    nums = {"n": len(fw), "n_eff": neff, "intact_fw": float((fw.fate == "intact").mean())}
    bw = d[d.direction == "backward"]
    nums["intact_bw"] = float((bw.fate == "intact").mean())
    for k, (c, (v, sg)) in enumerate(obs.items()):
        ax = axes[0, k]
        x = fw[c].values
        ok = np.isfinite(x)
        lo, hi = np.nanpercentile(x, [0.5, 99.5])
        lo, hi = min(lo, v - 3 * sg), max(hi, v + 3 * sg)
        bins = np.linspace(lo, hi, 40)
        ax.hist(x[ok], bins=bins, color="#B8C4D9", label="all members (prior)")
        ax.hist(
            x[ok],
            bins=bins,
            weights=w[ok] * ok.sum(),
            histtype="step",
            color="#1F4FB0",
            lw=1.4,
            label="weighted by the observations",
        )
        ax.axvspan(v - sg, v + sg, color="#C2571A", alpha=0.18)
        ax.axvline(v, color="#C2571A", lw=1.2, label="published")
        name = c.replace("_U_measured_yr", "")
        ax.set_xlabel(f"apsidal period of orbit {name} in SWINGERsim [yr]")
        ax.set_ylabel("members")
        q = np.nanpercentile(x, [16, 50, 84])
        nums[c] = {
            "median": float(q[1]),
            "p16": float(q[0]),
            "p84": float(q[2]),
            "published": v,
            "frac_finite": float(ok.mean()),
        }
        if k == 0:
            ax.legend(frameon=False, fontsize=7)
    axes[0, 0].set_title(
        f"{s['name']}: ensemble of {len(fw)} systems within the published errors", loc="left"
    )
    numbers.setdefault(sys_, {})["ensemble"] = nums
    fig.tight_layout()
    return _save(fig, out, f"{sys_}_ensemble")


# ------------------------------------------------------------------ mini-fit
_LINEAR = re.compile(r"^\s+(\S+)\s+.*\(linear ([-\d.eE+]+)\)")


def linear_errors(report: Path) -> dict[str, float]:
    """Linear (covariance) errors of the parameters, as stage2_report.txt of an MCMC run gives
    them in brackets: the same run as the emcee errors."""
    out = {}
    if report.exists():
        for line in report.read_text(encoding="utf-8").splitlines():
            m = _LINEAR.match(line)
            if m:
                out[m.group(1)] = float(m.group(2))
    return out


ANGLES = ("omega", "Omega", "M0", "i", "inc")


def fit_table(root: Path, sys_: str) -> list[dict]:
    base = root / "output" / "squadra_test" / sys_
    pub = yaml.safe_load(config_path(root, sys_, "published").read_text())
    rows = []
    lin = linear_errors(base / "fit_mcmc" / "stage2_report.txt")
    for r in csv.DictReader((base / "fit_mcmc" / "stage2_params.csv").read_text().splitlines()):
        name = r["parameter"]
        try:
            pv = nested_value(pub, name)  # orbits.*, masses.*, stars.*.*
        except (KeyError, TypeError):
            pv = float("nan")  # not in the published configuration (rv.*, offsets)
        rows.append(
            {
                "parameter": name,
                "published": pv,
                "start": float(r["start"]),
                "fitted": float(r["value"]),
                "err_mcmc": float(r["error"]),
                "err_linear": lin.get(name, float("nan")),
                "unit": r["unit"],
            }
        )
    return rows


def _offset(x: float, ref: float, name: str) -> float:
    """x - ref, wrapped to (-180, 180] for angles."""
    d = x - ref
    return float(wrap_deg(d)) if name.split(".")[-1] in ANGLES else d


def fig_fit(root: Path, sys_: str, out: Path, numbers: dict) -> str:
    s = SYSTEMS[sys_]
    rows = fit_table(root, sys_)
    base = root / "output" / "squadra_test" / sys_
    samp = pd.read_csv(base / "fit_mcmc" / "stage2_samples.csv", comment="#")
    fig, axes = plt.subplots(
        1, 2, figsize=(11, 0.32 * len(rows) + 1.8), gridspec_kw={"width_ratios": [1.4, 1]}
    )
    ax = axes[0]
    y = np.arange(len(rows))[::-1]
    for yy, r in zip(y, rows, strict=True):
        e = r["err_mcmc"] if r["err_mcmc"] > 0 else np.nan
        p = r["parameter"]
        ax.errorbar(
            _offset(r["fitted"], r["published"], p) / e,
            yy,
            xerr=1.0,
            fmt="o",
            color=C_MODEL["fitted"],
            ms=4,
        )
        ax.plot(
            _offset(r["start"], r["published"], p) / e,
            yy,
            "|",
            color=C_MODEL["anchored"],
            ms=10,
            mew=2,
        )
    ax.axvline(0, color=C_MODEL["published"], lw=1.2)
    ax.set_yticks(y)
    ax.set_yticklabels([r["parameter"].replace("orbits.", "") for r in rows], fontsize=7.5)
    ax.set_xscale("symlog", linthresh=3)
    ax.set_xlabel("(SWINGERsim − published) / σ of the fit")
    ax.set_title(f"{s['name']}: fitted elements against the published ones", loc="left")
    ax.plot([], [], "o", color=C_MODEL["fitted"], label="SWINGERsim fit ± 1σ")
    ax.plot([], [], "|", color=C_MODEL["anchored"], mew=2, ms=10, label="start (anchored)")
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    # the most strongly correlated pair of parameters in the samples
    c = samp.corr().values.copy()
    np.fill_diagonal(c, 0)
    i, j = np.unravel_index(np.argmax(np.abs(c)), c.shape)
    a2 = axes[1]
    a2.plot(
        samp.iloc[:, i],
        samp.iloc[:, j],
        ".",
        ms=1.2,
        alpha=0.25,
        color=C_MODEL["fitted"],
        rasterized=True,
    )
    a2.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    a2.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
    a2.set_xlabel(samp.columns[i].replace("orbits.", ""))
    a2.set_ylabel(samp.columns[j].replace("orbits.", ""))
    a2.set_title(f"emcee samples, strongest correlation r = {c[i, j]:+.2f}", loc="left")
    numbers.setdefault(sys_, {})["fit"] = rows
    numbers[sys_]["fit_corr"] = {"pair": [samp.columns[i], samp.columns[j]], "r": float(c[i, j])}
    fig.tight_layout()
    return _save(fig, out, f"{sys_}_fit")


def fig_summary(numbers: dict, out: Path) -> str:
    labels, vals = [], {m: [] for m in MODELS}
    for sys_, s in SYSTEMS.items():
        for p, v in numbers[sys_]["timings"].items():
            labels.append(f"{s['name']}\n{p}")
            for m in MODELS:
                vals[m].append(v[f"rms_{m}"])
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(10, 3.4))
    for k, m in enumerate(MODELS):
        ax.bar(x + (k - 1) * 0.27, vals[m], 0.25, color=C_MODEL[m], label=LABEL[m])
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("rms of measured − model [min]")
    ax.legend(frameon=False, ncol=3, fontsize=8)
    ax.set_title("Times of minima: how far the three models are from the measurements", loc="left")
    fig.tight_layout()
    return _save(fig, out, "summary_timings")


def make_all(root: Path, out: Path, png: bool = True) -> dict:
    global FIGURES
    FIGURES = FigureWriter(out, png=png)
    numbers: dict = {}
    for sys_ in SYSTEMS:
        print("figures", sys_, flush=True)
        fig_lightcurves(root, sys_, out, numbers)
        sec = {"bu_cmi": "S34", "tic120362137": "S40-41", "v994_her": "S26"}[sys_]
        fig_lc_components(root, sys_, out, sec)
        fig_oc(root, sys_, out, numbers)
        fig_evolution(root, sys_, out, numbers)
        fig_ensemble(root, sys_, out, numbers)
        fig_fit(root, sys_, out, numbers)
    fig_summary(numbers, out)
    atomic_write_text(out / "numbers.json", json.dumps(numbers, indent=1, default=float))
    FIGURES.finish()
    return numbers
