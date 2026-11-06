"""Figures and numbers of part III of the SQUADRA test (docs/<lang>/squadra_test3.md): three eclipsing binaries
whose inclination changes - V907 Sco, QX Cas and AH Cep.

    python tools/squadra3_report.py --root . [--output output] [--out docs/figures/squadra_test3]

Every figure is a plotly figure: written as PNG (for the PDF; needs kaleido and a Chrome/Chromium, e.g.
BROWSER_PATH=/opt/pw-browsers/chromium-1194/chrome-linux/chrome) and as JSON (for the interactive HTML version
made by tools/build_report_html.py; a PNG that cannot be written is fatal unless --no-png). Model curves
(inclination against time, eclipse times, the precession of the pair's orbit, light curves) are integrated here from
the configurations and cached in output/squadra_test3/<system>/report/; a cached curve is recomputed when its
configuration (the dict actually integrated), its data files, its request (epoch, grid of years) or the swingersim
version change. The epoch (BJD of t = 0) of every configuration is read from the configuration. Inputs: configs/squadra_test3/, data/observed/<system>/, the outputs of
the jobs in output/squadra_test3/<system>/ (final fits: v907_sco/fit_v2, v907_sco/fit_lc_v2, qx_cas/fit_v2,
qx_cas/ensemble, ah_cep/fit; control jobs when present: v907_sco/fit_p99, v907_sco/fit_v2_mcmc,
ah_cep/ensemble_4body, ah_cep/fit_lc). Numbers: docs/figures/squadra_test3/numbers.json.
"""

from __future__ import annotations

import argparse
import copy
import glob
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import yaml
from plotly.subplots import make_subplots

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import (
    C_DATA,
    DPY,
    FigureWriter,
    JsonCache,
    atomic_write_text,
    band_luminosity,
    digest,
    eclipse_switches,
    ensemble_observed,
    ensemble_rows,
    ensemble_weights,
    file_cached,
    fit_epoch,
    nodal_period,
    read_table,
    time_bins,
)
from _common import bjd as BJD
from _common import year as YEAR

C_ART = "#C2571A"  # the published (article) system
C_Q = "#1F4FB0"  # swingersim
C_Q2 = "#2E7D4F"  # swingersim, second model (light-curve fit)
C_GREY = "rgba(140,140,140,0.22)"
FONT = dict(family="DejaVu Sans, Arial", size=13, color="#222")


def _layout(fig, h=420, w=1000, **kw):
    fig.update_layout(
        template="plotly_white",
        font=FONT,
        height=h,
        width=w,
        margin=dict(l=70, r=20, t=50, b=60),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0, font=dict(size=11)),
        **kw,
    )
    fig.update_xaxes(
        showline=True,
        linecolor="#888",
        mirror=False,
        gridcolor="#eee",
        zeroline=False,
        tickformat=".6~g",
    )
    fig.update_yaxes(
        showline=True, linecolor="#888", mirror=False, gridcolor="#eee", zeroline=False
    )
    return fig


class Out(FigureWriter):
    """Plotly figures: JSON (interactive HTML) and PNG (PDF)."""

    def save(self, fig, name: str, scale: float = 2.0):
        return self.save_plotly(fig, name, scale)


# ------------------------------------------------------------------ model computations (cached)
def _cfg(d):
    from swingersim.config import Config

    d = copy.deepcopy(d)
    d.pop("fit", None)
    d.pop("ensemble", None)
    return Config.model_validate(d)


def _load(path) -> dict:
    return yaml.safe_load(Path(path).read_text())


def _data_dirs(root: Path, out_root: Path | None = None) -> list[Path]:
    """Folders of observed data: the repository, and the data folder next to the output folder (the TESS light
    curves are not in the repository)."""
    dirs = [root / "data/observed"]
    if out_root is not None:
        extra = out_root.parent.parent / "data" / "observed"
        if extra.resolve() != dirs[0].resolve():
            dirs.append(extra)
    return dirs


DATA_DIRS: list[Path] = []  # extra data folders (set once by main)


def _data(root: Path, system_: str, pattern: str) -> Path:
    """The folder of observed data of a system that has files matching the pattern."""
    for d in [root / "data/observed" / system_] + [x / system_ for x in DATA_DIRS]:
        if glob.glob(str(d / pattern)):
            return d
    return root / "data/observed" / system_


def cached(cache: Path, key: str, compute, *inputs):
    """A model curve cached in cache/<key>.json, recomputed when one of its inputs changes: configuration dicts and
    request parameters (hashed), data files (Path: size and modification time), and the SWINGERsim version."""
    files = [x for x in inputs if isinstance(x, Path)]
    params = [digest(x) for x in inputs if not isinstance(x, Path)]
    return JsonCache(cache).get(key, compute, files, params)


def inc_curve(d, t0, years, pair="inner"):
    from swingersim.incfit import model_inclination

    m = model_inclination(_cfg(d), pair, BJD(years), t0, with_limit=True)
    return {"year": list(map(float, years)), "i": m[0].tolist(), "limit": m[1].tolist()}


def unfold(i_folded, years, ref_year=2021.5):
    """The inclination folded to 0-90 deg (light curves measure only sin i: i and 180 - i give the same eclipses)
    unfolded to a continuous curve on 0-180 deg. Every point can be i or 180 - i; the sequence with the smallest sum
    of squared second differences (the smoothest curve) is chosen by dynamic programming, so that a cusp at 90 deg
    (the orbit passes through edge-on) is unfolded while a smooth maximum below 90 deg stays as it is. The overall
    sign is chosen so that the inclination at ref_year is below 90 deg."""
    f = np.asarray(i_folded, float)
    n_ = len(f)
    if n_ < 3:
        return f.copy()
    v = np.stack([f, 180.0 - f], axis=1)  # value of state s at point k
    # cost[k][a, b]: best cost of a path ending with states (a at k-1, b at k)
    cost = np.zeros((2, 2))
    back = np.zeros((n_, 2, 2), dtype=int)
    for k in range(2, n_):
        new = np.full((2, 2), np.inf)
        for b_ in range(2):
            for c_ in range(2):
                cand = [
                    cost[a_, b_] + (v[k, c_] - 2 * v[k - 1, b_] + v[k - 2, a_]) ** 2
                    for a_ in range(2)
                ]
                j = int(np.argmin(cand))
                new[b_, c_] = cand[j]
                back[k, b_, c_] = j
        cost = new
    b_, c_ = np.unravel_index(int(np.argmin(cost)), (2, 2))
    st = np.zeros(n_, dtype=int)
    st[-1], st[-2] = c_, b_
    for k in range(n_ - 1, 1, -1):
        st[k - 2] = back[k, st[k - 1], st[k]]
    x = v[np.arange(n_), st]
    k0 = int(np.argmin(np.abs(np.asarray(years) - ref_year)))
    return 180.0 - x if x[k0] > 90.0 else x


def align_branch(i_folded, ref):
    """A folded inclination curve put on the branches of a reference unfolded curve: between the moments when it is
    near 90 deg (where i and 180 - i coincide) every stretch is taken as i or as 180 - i, whichever is closer on
    average to the reference."""
    f = np.asarray(i_folded, float)
    ref = np.asarray(ref, float)
    x = f.copy()
    far = np.abs(f - 90.0) >= 1.5
    k = 0
    while k < len(f):
        if not far[k]:
            k += 1
            continue
        j = k
        while j < len(f) and far[j]:
            j += 1
        a, b = f[k:j], 180.0 - f[k:j]
        x[k:j] = a if np.mean(np.abs(a - ref[k:j])) <= np.mean(np.abs(b - ref[k:j])) else b
        k = j
    return x


def fold(x):
    """An inclination as eclipses measure it: i and 180 - i alike, folded to 0-90 deg."""
    x = np.asarray(x, float)
    return np.where(x > 90.0, 180.0 - x, x)


def mirror_to(model_i, i_obs):
    """A measured inclination (<= 90 deg) drawn on the side of 90 deg where the model is."""
    return np.where(np.asarray(model_i) > 90.0, 180.0 - np.asarray(i_obs), np.asarray(i_obs))


def eclipse_band(fig, years, lim, name="no eclipses (outside the eclipse window)"):
    """Grey: the inclinations without eclipses, i < limit or i > 180 - limit."""
    lim = np.asarray(lim)
    fig.add_trace(
        go.Scatter(
            x=years,
            y=lim,
            mode="lines",
            line=dict(width=0),
            fill="tozeroy",
            fillcolor=C_GREY,
            name=name,
            hoverinfo="skip",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=np.r_[years, years[::-1]],
            y=np.r_[180.0 - lim, np.full(len(lim), 180.0)],
            fill="toself",
            fillcolor=C_GREY,
            line=dict(width=0),
            showlegend=False,
            hoverinfo="skip",
        )
    )
    fig.add_hline(y=90, line=dict(color="#999", width=1, dash="dot"))


def axis_track(d, t0, years, pair="inner"):
    """Unit vectors of the pair's orbital angular momentum and of the total angular momentum at the given calendar
    years (x, y on the sky, z away from the observer)."""
    from swingersim.backends.nbody import build_simulation

    cfg = _cfg(d)
    T = (BJD(years) - t0) / DPY
    out = np.full((len(T), 3), np.nan)
    tot = np.full((len(T), 3), np.nan)
    for sign in (1.0, -1.0):
        idx = np.nonzero(T >= 0 if sign > 0 else T < 0)[0]
        if not len(idx):
            continue
        idx = idx[np.argsort(sign * T[idx], kind="stable")]
        sim, _, _, root, _ = build_simulation(cfg)
        sim.collision = "none"
        if sign < 0:
            sim.dt = -abs(sim.dt)
        node = next(n for n in root.nodes() if n.name == pair)
        i, j = node.left, node.right
        ps = sim.particles
        for k in idx:
            sim.integrate(T[k], exact_finish_time=1)
            r = np.array(ps[j].xyz) - np.array(ps[i].xyz)
            v = np.array(ps[j].vxyz) - np.array(ps[i].vxyz)
            h = np.cross(r, v)
            out[k] = h / np.linalg.norm(h)
            L = sum(p.m * np.cross(np.array(p.xyz), np.array(p.vxyz)) for p in ps)
            tot[k] = L / np.linalg.norm(L)
    return {"year": list(map(float, years)), "h": out.tolist(), "L": tot.tolist()}


def oc_model(d, t0, y0, y1):
    """BJD of every model eclipse of the inner pair between the calendar years y0 and y1."""
    from swingersim.timing import eclipse_times

    cfg = _cfg(d)
    a = (BJD(y0) - t0) / DPY
    b = (BJD(y1) - t0) / DPY
    ecl = []
    if a < 0:
        ecl += [e for e in eclipse_times(cfg, float(a))[0] if e.t_obs_yr >= a and e.t_obs_yr <= b]
    if b > 0:
        ecl += [e for e in eclipse_times(cfg, float(b))[0] if e.t_obs_yr >= max(a, 0.0)]
    out = {"primary": [], "secondary": []}
    for e in ecl:
        out[e.kind].append(float(t0 + e.t_obs_yr * DPY))
    return {k: sorted(v) for k, v in out.items()}


def moved(m, to_model, pair="inner"):
    """A copy of the system model m moved so that its pair has the centre of mass of the same pair in to_model (to
    draw a published binary where the fitted triple has it)."""
    m = copy.deepcopy(m)
    a = np.array(next(o for o in m["orbits"] if o["name"] == pair)["com"])
    b = np.array(next(o for o in to_model["orbits"] if o["name"] == pair)["com"])
    d = b - a
    for s_ in m["stars"]:
        s_["pos"] = (np.array(s_["pos"]) + d).tolist()
    for o in m["orbits"]:
        o["com"] = (np.array(o["com"]) + d).tolist()
        o["path_left"] = (np.array(o["path_left"]) + d).tolist()
        o["path_right"] = (np.array(o["path_right"]) + d).tolist()
    return m


def system(d):
    from swingersim.systemview import system_model

    return system_model(_cfg(d))


# ------------------------------------------------------------------ the three views of a system
def _ellipse_traces(m, view, color, dash, name, show_legend=True, only=None, centre=None):
    pick = (0, 1) if view == "sky" else (0, 2)
    sgn = 1.0 if view == "sky" else -1.0  # top view: -z up, i.e. the observer at the top (z points away)
    c0 = centre or [0, 0, 0]
    tr = []
    for k, o in enumerate(m["orbits"]):
        if only and o["name"] not in only:
            continue
        for side, path in ((0, o["path_left"]), (1, o["path_right"])):
            p = np.array(path)
            tr.append(
                go.Scatter(
                    x=p[:, pick[0]] - c0[pick[0]],
                    y=sgn * (p[:, pick[1]] - c0[pick[1]]),
                    mode="lines",
                    line=dict(
                        color=color, width=2 if not dash else 1.6, dash="dash" if dash else "solid"
                    ),
                    name=name,
                    legendgroup=name,
                    showlegend=show_legend and k == 0 and side == 0,
                    hovertemplate=f"{o['children'][side]} around the centre of orbit {o['name']}<extra>{name}</extra>",
                )
            )
    return tr


def _stars(m, view, color, only=None, centre=None, label=True):
    pick = (0, 1) if view == "sky" else (0, 2)
    sgn = 1.0 if view == "sky" else -1.0
    c0 = centre or [0, 0, 0]
    sel = [s for s in m["stars"] if (only is None or s["name"] in only)]
    return go.Scatter(
        x=[s["pos"][pick[0]] - c0[pick[0]] for s in sel],
        y=[sgn * (s["pos"][pick[1]] - c0[pick[1]]) for s in sel],
        mode="markers+text" if label else "markers",
        text=[s["name"] for s in sel],
        textposition="bottom right",
        textfont=dict(size=11, color=color),
        marker=dict(size=8, color=color, line=dict(color="#333", width=0.6)),
        showlegend=False,
        hoverinfo="skip",
    )


def _top_labels(m):
    """One label per body in the view from above: the pair as a whole ('pair m1+m2') and every other star, so that the
    names of the two close stars do not overlap."""
    o_in = next(o for o in m["orbits"] if o["name"] == "inner")
    pair = set(o_in["children"])
    xs, ys, txt = [], [], []
    com = o_in["com"]
    xs.append(com[0])
    ys.append(-com[2])
    txt.append("pair " + "+".join(sorted(pair)))
    for s in m["stars"]:
        if s["name"] not in pair:
            xs.append(s["pos"][0])
            ys.append(-s["pos"][2])
            txt.append(s["name"])
    return xs, ys, txt


def views_figure(
    models,
    title_left="From above (the observer is at the top)",
    zoom_pair="inner",
    cone=None,
    i_lim=None,
    cone_note="",
    label_years=(),
    sky_pad=1.35,
):
    """Three panels: the whole system from above, the eclipsing pair on the sky (stars to scale), the precession of
    the pair's orbital axis (zoomed on the track). models: [(label, system_model dict, colour, dash)]. Star names are
    written only for the last model (the swingersim one), so that the labels of nearly identical models do not overlap."""
    fig = make_subplots(
        rows=1,
        cols=3,
        subplot_titles=(
            title_left,
            "The eclipsing pair on the sky (stars to scale)",
            "Where the pair's orbit points",
        ),
        horizontal_spacing=0.09,
    )
    sky_extent = 0.0
    for n_m, (label, m, col, dash) in enumerate(models):
        last = n_m == len(models) - 1
        for tr in _ellipse_traces(m, "top", col, dash, label):
            fig.add_trace(tr, 1, 1)
        fig.add_trace(_stars(m, "top", col, label=False), 1, 1)
        if last:
            xs, ys, txt = _top_labels(m)
            fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=ys,
                    mode="text",
                    text=txt,
                    textposition="bottom right",
                    textfont=dict(size=11, color="#222"),
                    showlegend=False,
                    hoverinfo="skip",
                ),
                1,
                1,
            )
        o = next(o for o in m["orbits"] if o["name"] == zoom_pair)
        for tr in _ellipse_traces(
            m, "sky", col, dash, label, show_legend=False, only=[zoom_pair], centre=o["com"]
        ):
            fig.add_trace(tr, 1, 2)
            sky_extent = max(sky_extent, float(np.nanmax(np.abs(np.r_[tr.x, tr.y]))))
        for s in m["stars"]:
            if s["name"] in o["children"] and s["R_au"]:
                x = s["pos"][0] - o["com"][0]
                y = s["pos"][1] - o["com"][1]
                sky_extent = max(sky_extent, abs(x) + s["R_au"], abs(y) + s["R_au"])
                fig.add_shape(
                    type="circle",
                    x0=x - s["R_au"],
                    x1=x + s["R_au"],
                    y0=y - s["R_au"],
                    y1=y + s["R_au"],
                    line=dict(color=col, width=1.2, dash="dash" if dash else "solid"),
                    fillcolor="rgba(255,236,190,0.55)" if not dash else "rgba(0,0,0,0)",
                    layer="below",
                    row=1,
                    col=2,
                )
                if last:
                    fig.add_annotation(
                        x=x,
                        y=y,
                        text=s["name"],
                        showarrow=False,
                        font=dict(size=11, color="#222"),
                        row=1,
                        col=2,
                    )
    fig.update_xaxes(title_text="au", row=1, col=1)
    fig.update_yaxes(title_text="au (↓ towards Earth)", scaleanchor="x", scaleratio=1, row=1, col=1)
    ext = sky_extent * 1.08
    fig.update_xaxes(title_text="au (on the sky)", range=[-ext, ext], row=1, col=2)
    fig.update_yaxes(
        title_text="au", range=[-ext, ext], scaleanchor="x2", scaleratio=1, row=1, col=2
    )
    # the cone: the pair's orbital axis projected on the sky; inside the ring of radius sin(i_lim) it points too
    # close to the line of sight for eclipses, between the ring and the rim it eclipses. The panel is zoomed on the
    # track of the axis; the rim (r = 1, axis in the sky plane = edge-on orbit) and the inner edge of the eclipse ring
    # are drawn as arcs.
    th = np.linspace(0, 2 * np.pi, 1441)
    allx, ally = [], []
    for _, tr_, _, _ in cone or []:
        h = np.array(tr_["h"], float)
        allx += list(h[:, 0])
        ally += list(h[:, 1])
        L = np.array(tr_["L"], float)[0]
        allx.append(L[0])
        ally.append(L[1])
    if allx:
        cx, cy = 0.5 * (min(allx) + max(allx)), 0.5 * (min(ally) + max(ally))
        half = 0.5 * max(max(allx) - min(allx), max(ally) - min(ally)) * sky_pad + 0.03
        hspan = max(
            max(np.ptp(np.array(tr_["h"], float)[:, 0]), np.ptp(np.array(tr_["h"], float)[:, 1]))
            for _, tr_, _, _ in cone
        )
        if hspan < 0.05:  # the axis does not move: show the whole sphere of directions
            cx, cy, half = 0.0, 0.0, 1.05
    else:
        cx, cy, half = 0.0, 0.0, 1.05
    fig.add_trace(
        go.Scatter(
            x=np.cos(th),
            y=np.sin(th),
            mode="lines",
            line=dict(color="#777", width=1),
            showlegend=False,
            hoverinfo="skip",
        ),
        1,
        3,
    )
    if i_lim:
        r = math.sin(math.radians(i_lim))
        fig.add_trace(
            go.Scatter(
                x=np.r_[np.cos(th), r * np.cos(th[::-1])],
                y=np.r_[np.sin(th), r * np.sin(th[::-1])],
                fill="toself",
                fillcolor="rgba(255,196,0,0.22)",
                line=dict(width=0),
                name="eclipses (axis close to the sky plane)",
                hoverinfo="skip",
            ),
            1,
            3,
        )
        fig.add_trace(
            go.Scatter(
                x=r * np.cos(th),
                y=r * np.sin(th),
                mode="lines",
                line=dict(color="#C9A227", width=1, dash="dot"),
                showlegend=False,
                hoverinfo="skip",
            ),
            1,
            3,
        )
    for label, tr_, col, dash in cone or []:
        h = np.array(tr_["h"], float)
        yr = np.array(tr_["year"])
        hx, hy = (
            h[:, 0],
            h[:, 1],
        )  # the axis projected on the sky (no folding: the track stays continuous)
        fig.add_trace(
            go.Scatter(
                x=hx,
                y=hy,
                mode="lines",
                line=dict(color=col, width=2, dash="dash" if dash else "solid"),
                name=f"{label}: axis of the pair",
                text=[f"{y:.0f}" for y in yr],
                hovertemplate="%{text}<extra></extra>",
                showlegend=False,
            ),
            1,
            3,
        )
        if not dash:  # year labels only on the swingersim track
            marks = [
                int(np.argmin(np.abs(yr - y))) for y in label_years if yr.min() <= y <= yr.max()
            ]
            cyt = 0.5 * (hy.min() + hy.max())
            fig.add_trace(
                go.Scatter(
                    x=hx[marks],
                    y=hy[marks],
                    mode="markers+text",
                    text=[f"{yr[k]:.0f}" for k in marks],
                    textposition=["bottom center" if hy[k] < cyt else "top center" for k in marks],
                    textfont=dict(size=10, color=col),
                    marker=dict(size=5, color=col),
                    showlegend=False,
                    hoverinfo="skip",
                ),
                1,
                3,
            )
            k0 = int(np.argmin(np.abs(yr - 2025.0)))
            fig.add_trace(
                go.Scatter(
                    x=[hx[k0]],
                    y=[hy[k0]],
                    mode="markers",
                    marker=dict(size=10, color=col, line=dict(color="white", width=1.5)),
                    name=f"{label}: axis of the pair today",
                    hovertext=f"axis of the pair, {yr[k0]:.0f}",
                    hoverinfo="text",
                    showlegend=False,
                ),
                1,
                3,
            )
        L = np.array(tr_["L"], float)[0]
        fig.add_trace(
            go.Scatter(
                x=[L[0]],
                y=[L[1]],
                mode="markers",
                marker=dict(symbol="x-thin", size=12, color=col, line=dict(width=2, color=col)),
                name=f"{label}: total angular momentum",
                showlegend=False,
            ),
            1,
            3,
        )
    fig.update_xaxes(range=[cx - half, cx + half], title_text="axis on the sky, x", row=1, col=3)
    fig.update_yaxes(
        range=[cy - half, cy + half],
        scaleanchor="x3",
        scaleratio=1,
        title_text="axis on the sky, y",
        row=1,
        col=3,
    )
    if cone_note:
        fig.add_annotation(
            text=cone_note,
            xref="x3 domain",
            yref="y3 domain",
            x=0.5,
            y=-0.30,
            showarrow=False,
            font=dict(size=10, color="#555"),
        )
    _layout(fig, h=540, w=1350)
    fig.update_layout(
        legend=dict(orientation="h", yanchor="top", y=-0.22, x=0.0),
        margin=dict(l=70, r=30, t=50, b=130),
    )
    return fig


# ------------------------------------------------------------------ V907 Sco
def v907(root: Path, out_root: Path, O: Out, N: dict):
    s, o = "v907_sco", out_root / "v907_sco"
    cache = o / "report"
    pub_f, fit_f = (
        root / "configs/squadra_test3/v907_sco_published.yaml",
        o / "fit_v2/stage2_fitted.yaml",
    )
    pub, fit = _load(pub_f), _load(fit_f)
    t0, t0_pub = fit_epoch(fit_f), fit_epoch(pub_f)
    n = N.setdefault(s, {})
    # systems
    m_pub = cached(cache, "sys_published", lambda: system(pub), pub)
    m_fit = cached(cache, "sys_fit_v2", lambda: system(fit), fit)
    yrs = np.arange(2021.0, 2021.0 + 66.0, 0.5)
    cone_pub = cached(
        cache, "axis_published", lambda: axis_track(pub, t0_pub, yrs), pub, t0_pub, yrs
    )
    cone_fit = cached(cache, "axis_fit_v2", lambda: axis_track(fit, t0, yrs), fit, t0, yrs)
    fig = views_figure(
        [("Zasche et al. (2023)", m_pub, C_ART, True), ("SWINGERsim fit", m_fit, C_Q, False)],
        cone=[
            ("Zasche et al. (2023)", cone_pub, C_ART, True),
            ("SWINGERsim fit", cone_fit, C_Q, False),
        ],
        i_lim=m_fit["pairs"][0]["i_limit_deg"],
        cone_note="one turn of the axis, 2021-2087; yellow: eclipses",
        label_years=(2030, 2040, 2050, 2060, 2070, 2080),
    )
    O.save(fig, "v907_models")
    n["system"] = {
        "published": {
            "i_mut": m_pub["mutual"][0]["deg"],
            "triples": m_pub["triples"],
            "pairs": m_pub["pairs"],
        },
        "fit": {
            "i_mut": m_fit["mutual"][0]["deg"],
            "triples": m_fit["triples"],
            "pairs": m_fit["pairs"],
            "orbits": [
                {k: x[k] for k in ("name", "P_d", "a_au", "e", "i_deg")} for x in m_fit["orbits"]
            ],
        },
    }
    # inclination history
    years = np.arange(1890.0, 2071.0, 0.5)
    c_pub = cached(
        cache, "inc_published", lambda: inc_curve(pub, t0_pub, years), pub, t0_pub, years
    )
    c_fit = cached(cache, "inc_fit_v2", lambda: inc_curve(fit, t0, years), fit, t0, years)
    inc_f = root / "data/observed/v907_sco/inclination_zasche2023.csv"
    inc = read_table(inc_f)
    post_f = o / "fit_v2_mcmc" / "posterior_inc.json"
    post = json.loads(post_f.read_text()) if post_f.exists() else []
    lim = np.array(c_fit["limit"])
    u_fit = unfold(c_fit["i"], years)
    u_pub = unfold(c_pub["i"], years)
    # top: the inclination as eclipses measure it (folded to 0-90 deg, i and 180 - i alike: the "M" shape where the
    # orbit passes edge-on); bottom: the same curve unfolded to 0-180 deg (the side of the sky the N-body model knows)
    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.07,
        row_heights=[0.5, 0.5],
        subplot_titles=(
            "as the light curve measures it: 0–90°",
            "the same orbit unfolded: 0–180° (which side of edge-on)",
        ),
    )
    fig.add_trace(
        go.Scatter(
            x=years,
            y=lim,
            mode="lines",
            line=dict(width=0),
            fill="tozeroy",
            fillcolor=C_GREY,
            name="no eclipses (below the eclipse limit)",
            hoverinfo="skip",
        ),
        1,
        1,
    )
    fig.add_trace(
        go.Scatter(
            x=years,
            y=lim,
            mode="lines",
            line=dict(width=0),
            fill="tozeroy",
            fillcolor=C_GREY,
            showlegend=False,
            hoverinfo="skip",
        ),
        2,
        1,
    )
    fig.add_trace(
        go.Scatter(
            x=np.r_[years, years[::-1]],
            y=np.r_[180.0 - lim, np.full(len(lim), 180.0)],
            fill="toself",
            fillcolor=C_GREY,
            line=dict(width=0),
            showlegend=False,
            hoverinfo="skip",
        ),
        2,
        1,
    )
    for r_ in (1, 2):
        fig.add_hline(y=90, line=dict(color="#999", width=1, dash="dot"), row=r_, col=1)
    if post:
        Pi = np.array([align_branch(p_["i"], u_fit) for p_ in post])
        Pf = fold(Pi)
        for r_, P_ in ((1, Pf), (2, Pi)):
            lo, hi = np.percentile(P_, 2.5, axis=0), np.percentile(P_, 97.5, axis=0)
            fig.add_trace(
                go.Scatter(
                    x=np.r_[years, years[::-1]],
                    y=np.r_[hi, lo[::-1]],
                    fill="toself",
                    fillcolor="rgba(31,79,176,0.20)",
                    line=dict(width=0),
                    name=f"SWINGERsim fit: 95 % range of {len(post)} posterior samples",
                    hoverinfo="skip",
                    showlegend=r_ == 1,
                ),
                r_,
                1,
            )
    for r_, f_ in ((1, fold), (2, np.asarray)):
        fig.add_trace(
            go.Scatter(
                x=years,
                y=f_(u_pub),
                mode="lines",
                line=dict(color=C_ART, width=2, dash="dash"),
                name="system of Zasche et al. (2023), integrated by SWINGERsim",
                showlegend=r_ == 1,
            ),
            r_,
            1,
        )
        fig.add_trace(
            go.Scatter(
                x=years,
                y=f_(u_fit),
                mode="lines",
                line=dict(color=C_Q, width=2.5),
                name="SWINGERsim fit (minima + velocities + inclinations)",
                showlegend=r_ == 1,
            ),
            r_,
            1,
        )
    yo = YEAR(inc.mid_hjd)
    for r_, y_ in ((1, inc.i_deg), (2, mirror_to(np.interp(yo, years, u_fit), inc.i_deg))):
        fig.add_trace(
            go.Scatter(
                x=yo,
                y=y_,
                mode="markers",
                error_y=dict(type="data", array=inc.sigma_deg, thickness=1.2),
                marker=dict(color=C_DATA, size=8),
                name="measured inclination (Zasche et al. 2023, Table 2)",
                text=inc.source,
                hovertemplate="%{x:.1f}: %{text}<extra></extra>",
                showlegend=r_ == 1,
            ),
            r_,
            1,
        )
    for x0, x1, txt in (
        (1899, 1918, "eclipses seen<br>1899–1918"),
        (1963, 1986, "eclipses seen<br>1963–1986"),
        (2019, 2025, "eclipses<br>2019–"),
    ):
        for r_ in (1, 2):
            fig.add_vrect(
                x0=x0, x1=x1, fillcolor="rgba(46,125,79,0.10)", line_width=0, row=r_, col=1
            )
        fig.add_annotation(
            x=0.5 * (x0 + x1),
            y=53.5,
            text=txt,
            showarrow=False,
            font=dict(size=10, color="#2E5D3F"),
            row=1,
            col=1,
        )
    fig.update_xaxes(range=[1890, 2070])
    fig.update_xaxes(title_text="year", row=2, col=1)
    fig.update_yaxes(title_text="inclination [°]", range=[50, 92], row=1, col=1)
    fig.update_yaxes(title_text="inclination [°]", range=[50, 100], row=2, col=1)
    _layout(fig, h=820)
    fig.update_layout(legend=dict(y=1.03), margin=dict(t=130))
    O.save(fig, "v907_inclination")
    n["eclipse_changes"] = {
        "published": eclipse_switches(years, c_pub["i"], c_pub["limit"]),
        "fit": eclipse_switches(years, c_fit["i"], c_fit["limit"]),
    }
    if post:
        ch = [p_["changes"] for p_ in post]

        def pick(lo_, hi_):
            v = [next((x for x in c if lo_ <= x <= hi_), np.nan) for c in ch]
            return [float(np.nanpercentile(v, q)) for q in (2.5, 50, 97.5)]

        n["eclipse_changes_posterior"] = {
            "off_1920s": pick(1912, 1935),
            "on_1950s": pick(1940, 1965),
            "off_1980s": pick(1975, 2000),
            "on_2010s": pick(2005, 2025),
            "end": pick(2035, 2070),
        }
        pn = [p_["Pnod"] for p_ in post]
        n["nodal_period_posterior_yr"] = [
            float(np.percentile(pn, q)) for q in (2.5, 16, 50, 84, 97.5)
        ]
        n["posterior_n"] = len(post)
    # the 14 inclinations: chi^2 of the models at the exact epochs
    from swingersim.incfit import model_inclination

    def chi2_exact(d, tt0, key):
        mi = cached(
            cache,
            f"inc_epochs_{key}",
            lambda: model_inclination(_cfg(d), "inner", inc.mid_hjd.values, tt0).tolist(),
            d,
            tt0,
            inc_f,
        )
        return float(
            np.sum(((inc.i_deg.values - np.array(mi)) / np.maximum(inc.sigma_deg.values, 0.1)) ** 2)
        )

    n["chi2_inclination"] = {
        "published": chi2_exact(pub, t0_pub, "published"),
        "fit": chi2_exact(fit, t0, "fit_v2"),
    }
    # nodal period from the axis track
    n["nodal_period_fit_yr"] = nodal_period(yrs, cone_fit["h"], np.array(cone_fit["L"])[0])
    n["nodal_period_published_yr"] = nodal_period(yrs, cone_pub["h"], np.array(cone_pub["L"])[0])
    # light-time semi-amplitude of the pair on the outer orbit (to separate it from the dynamical eclipse-time signal)
    oo = next(x for x in m_fit["orbits"] if x["name"] == "out")
    mt = sum(fit["masses"].values())
    w2 = math.radians(fit["orbits"]["out"]["omega"])
    n["ltte_semi_amplitude_min"] = float(
        oo["a_au"]
        * fit["masses"]["m3"]
        / mt
        * math.sin(math.radians(oo["i_deg"]))
        * 499.005
        / 60.0
        * math.sqrt(1 - (oo["e"] * math.cos(w2)) ** 2)
    )
    # minima: O - C with the model
    res = pd.read_csv(o / "fit_v2/stage2_residuals.csv")
    oc = cached(
        cache, "oc_fit_v2", lambda: oc_model(fit, t0, 2021.2, 2025.7), fit, t0, (2021.2, 2025.7)
    )
    fig = oc_figure(res, oc, "inner", "TESS minima 2021–2025")
    O.save(fig, "v907_minima")

    # what the eclipse-time signal is made of: the light-travel time across the outer orbit and the dynamical part
    # (the third star changes the period and the orientation of the pair's orbit), one outer period 2021.6-2022.1
    def oc_parts():
        from swingersim.timing import eclipse_times

        ecl, _ = eclipse_times(_cfg(fit), (BJD(2022.25) - t0) / DPY)
        pr = [e for e in ecl if e.kind == "primary" and e.t_obs_yr * DPY + t0 > BJD(2021.55)]
        tg = np.array([e.t_yr for e in pr]) * DPY + t0
        to_ = np.array([e.t_obs_yr for e in pr]) * DPY + t0
        E = np.round((to_ - to_[0]) / fit["orbits"]["inner"]["P"])
        Pl, T0l = np.polyfit(E, to_, 1)
        return {
            "t": to_.tolist(),
            "total": ((to_ - (T0l + Pl * E)) * 1440).tolist(),
            "dyn": ((tg - (T0l + Pl * E)) * 1440).tolist(),
            "lt": ((to_ - tg) * 1440).tolist(),
        }

    pa = cached(cache, "oc_parts", oc_parts, fit, t0)
    tt_ = np.array(pa["t"])
    lt_ = np.array(pa["lt"])
    lt_ -= lt_.mean()
    dy = np.array(pa["dyn"])
    dy -= dy.mean()
    tot = np.array(pa["total"])
    tot -= tot.mean()
    Pout = fit["orbits"]["out"]["P"]
    tper = t0 - fit["orbits"]["out"]["M0"] / 360.0 * Pout + np.arange(-3, 4) * Pout
    fig = go.Figure()
    for tp in tper:
        if tt_.min() < tp < tt_.max():
            fig.add_vline(x=float(YEAR(tp)), line=dict(color=C_ART, width=1, dash="dot"))
            fig.add_annotation(
                x=float(YEAR(tp)),
                y=1.0,
                yref="paper",
                text="periastron of the third star",
                showarrow=False,
                yanchor="bottom",
                font=dict(size=10, color=C_ART),
            )
    fig.add_trace(
        go.Scatter(
            x=YEAR(tt_),
            y=tot,
            mode="lines+markers",
            marker=dict(size=4),
            line=dict(color=C_Q, width=2.2),
            name="model: total (as in Fig. 3)",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=YEAR(tt_),
            y=dy,
            mode="lines",
            line=dict(color=C_Q2, width=1.6, dash="dash"),
            name="dynamical part (the third star perturbs the pair)",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=YEAR(tt_),
            y=lt_,
            mode="lines",
            line=dict(color=C_ART, width=1.6, dash="dot"),
            name="light-travel time across the outer orbit",
        )
    )
    fig.update_xaxes(title_text="year")
    fig.update_yaxes(title_text="O − C of the primary eclipses [min]")
    _layout(fig, h=420)
    fig.update_layout(legend=dict(y=1.08), margin=dict(t=110))
    O.save(fig, "v907_oc_parts")
    n["oc_parts_ptp_min"] = {
        "total": float(np.ptp(tot)),
        "dynamical": float(np.ptp(dy)),
        "light_time": float(np.ptp(lt_)),
    }
    r = res[res.status == "fitted"]
    n["minima"] = {
        "n": len(r),
        "rms_min": float(np.sqrt(np.mean(r.o_minus_model_min**2))),
        "mean_p": float(r[r.type == "p"].o_minus_model_min.mean()),
        "mean_s": float(r[r.type == "s"].o_minus_model_min.mean()),
    }
    # light curve: what the reference fit (which never saw the light curves) predicts for the TESS eclipses
    pred_f = cache / "lc_pred_fit_v2.csv"
    win_f = root / "data/observed/v907_sco/lc_tess_windows.csv"

    def predict(tmp):
        from swingersim.photometry import light_curve

        dlc = read_table(win_f)
        dlc["model"] = light_curve(
            _cfg(fit), dlc.time.values - t0, "TESS", exposure_min=10.0, n_sub=3
        ).flux
        dlc.to_csv(tmp, index=False)

    file_cached(pred_f, predict, [win_f], [digest(fit), t0])
    dlc = pd.read_csv(pred_f)
    tt = dlc.time.values
    win = np.concatenate([[0], np.cumsum(np.diff(tt) > 0.1)])
    cal = np.empty(len(dlc))
    dep = []
    for w_ in np.unique(win):
        k = win == w_
        X = np.column_stack([np.ones(k.sum()), tt[k] - tt[k].mean()])
        c = np.linalg.lstsq(X, dlc.flux.values[k] - dlc.model.values[k], rcond=None)[0]
        cal[k] = dlc.model.values[k] + X @ c
        dep.append((dlc.segment.values[k][0], 1 - dlc.flux.values[k].min(), 1 - cal[k].min()))
    lcp = pd.DataFrame(
        {
            "bjd_tdb": tt,
            "flux": dlc.flux.values,
            "model": cal,
            "o_minus_model": dlc.flux.values - cal,
            "segment": dlc.segment.values,
        }
    )
    O.save(
        lc_figure(lcp, 4, model_name="SWINGERsim fit, predicted (light curves not fitted)"),
        "v907_lightcurve",
    )
    dd = pd.DataFrame(dep, columns=["seg", "obs", "mod"]).groupby("seg").median()
    n["lc_prediction"] = {
        seg: {"depth_obs": float(r_.obs), "depth_model": float(r_["mod"])}
        for seg, r_ in dd.iterrows()
    }
    n["lc_prediction_rms_ppt"] = {
        seg: float(1e3 * np.sqrt(np.mean(g.o_minus_model**2))) for seg, g in lcp.groupby("segment")
    }
    lc = pd.read_csv(o / "fit_lc_v2/stage2_lc.csv")
    n["lc_rms_ppt"] = {
        seg: float(1e3 * np.sqrt(np.mean(g.o_minus_model**2))) for seg, g in lc.groupby("segment")
    }
    st = fit["stars"]

    def lum(k):
        return band_luminosity(st[k]["R"], st[k]["Teff"], "TESS")

    n["third_light_model_tess"] = float(lum("m3") / sum(lum(k) for k in ("m1", "m2", "m3")))
    # radial velocities
    rv = pd.read_csv(o / "fit_v2/stage2_rv.csv")
    prim = cached(
        cache,
        "oc_rv_epochs",
        lambda: {
            "p": oc_model(fit, t0, 1980.0, 1990.0)["primary"]
            + oc_model(fit, t0, 2021.8, 2023.0)["primary"]
        },
        fit,
        t0,
    )["p"]
    O.save(rv_figure(rv, fit["orbits"]["inner"]["P"], prim), "v907_rv")
    n["rv_rms"] = {
        k: float(np.sqrt(np.mean(g[g.status == "fitted"].o_minus_model ** 2)))
        for k, g in rv.groupby("source")
    }
    # parameters
    for name, f in (
        ("fit", o / "fit_v2/stage2_params.csv"),
        ("fit_lc", o / "fit_lc_v2/stage2_params.csv"),
    ):
        p = pd.read_csv(f)
        n[f"params_{name}"] = {
            r.parameter: [float(r.value), float(r.error)] for r in p.itertuples()
        }
    for extra in ("fit_p99", "fit_v2_mcmc"):
        rp = o / extra / "stage2_report.txt"
        if rp.exists():
            n[extra] = rp.read_text()[:4000]
            pp = o / extra / "stage2_params.csv"
            if pp.exists():
                p = pd.read_csv(pp)
                n[f"params_{extra}"] = {
                    r.parameter: [float(r.value), float(r.error)] for r in p.itertuples()
                }


def oc_figure(res, oc, pair, title):
    """Measured minima and the model eclipses against one linear ephemeris (from the model primaries)."""
    tp = np.array(oc["primary"])
    P0 = np.median(np.diff(tp))
    E = np.round((tp - tp[0]) / P0)
    P, T0 = np.polyfit(E, tp, 1)

    def ocv(t, half):
        e = np.round((np.asarray(t) - T0) / P - half) + half
        return (np.asarray(t) - (T0 + P * e)) * 1440.0

    fig = make_subplots(
        rows=2, cols=1, shared_xaxes=True, row_heights=[0.68, 0.32], vertical_spacing=0.05
    )
    for kind, half, dash, nm in (
        ("primary", 0.0, "solid", "model: primary eclipses"),
        ("secondary", 0.5, "dot", "model: secondary eclipses"),
    ):
        t = np.array(oc[kind])
        if len(t):
            fig.add_trace(
                go.Scatter(
                    x=YEAR(t),
                    y=ocv(t, half),
                    mode="lines",
                    line=dict(color=C_Q, width=1.5, dash=dash),
                    name=nm,
                ),
                1,
                1,
            )
    r = res[(res.pair == pair) & (res.status == "fitted")]
    for kind, half, sym, nm in (
        ("p", 0.0, "circle", "measured primary"),
        ("s", 0.5, "circle-open", "measured secondary"),
    ):
        g = r[r.type == kind]
        fig.add_trace(
            go.Scatter(
                x=YEAR(g.bjd_tdb),
                y=ocv(g.bjd_tdb.values, half),
                mode="markers",
                marker=dict(symbol=sym, size=7, color=C_DATA, line=dict(width=1.2, color=C_DATA)),
                error_y=dict(type="data", array=g.sigma_d * 1440, thickness=0.8, width=0),
                name=nm,
                text=g.source,
                hovertemplate="%{x:.3f}: %{y:.1f} min %{text}<extra></extra>",
            ),
            1,
            1,
        )
        fig.add_trace(
            go.Scatter(
                x=YEAR(g.bjd_tdb),
                y=g.o_minus_model_min,
                mode="markers",
                marker=dict(symbol=sym, size=6, color=C_DATA, line=dict(width=1.1, color=C_DATA)),
                error_y=dict(type="data", array=g.sigma_d * 1440, thickness=0.8, width=0),
                showlegend=False,
            ),
            2,
            1,
        )
    fig.add_hline(y=0, line=dict(color="#888", width=1), row=2, col=1)
    fig.update_yaxes(title_text="O − C [min]", row=1, col=1)
    fig.update_yaxes(title_text="measured − model [min]", row=2, col=1)
    fig.update_xaxes(title_text="year", row=2, col=1)
    _layout(fig, h=540)
    fig.update_layout(
        title=dict(text=title, x=0.01, y=0.985, font=dict(size=13)), margin=dict(t=90)
    )
    return fig


def lc_figure(lc, nwin=4, model_name="SWINGERsim model"):
    lc = lc.sort_values("bjd_tdb").reset_index(drop=True)
    t = lc.bjd_tdb.values
    br = np.nonzero(np.diff(t) > 0.08)[0] + 1
    wins = [w for w in np.split(np.arange(len(t)), br) if len(w) > 15]
    wins.sort(key=lambda w: lc.flux.values[w].min())
    segs = list(dict.fromkeys(lc.segment.values[w[0]] for w in sorted(wins, key=lambda w: t[w[0]])))
    per = max(
        1, nwin // len(segs)
    )  # the deepest eclipses of every segment (sector), so that the change of depth shows
    picks = sorted(
        [w for sg in segs for w in [w for w in wins if lc.segment.values[w[0]] == sg][:per]],
        key=lambda w: t[w[0]],
    )
    fig = make_subplots(
        rows=2,
        cols=len(picks),
        shared_yaxes="rows",
        row_heights=[0.7, 0.3],
        vertical_spacing=0.06,
        horizontal_spacing=0.03,
        subplot_titles=[f"{lc.segment.values[w[0]]}, {YEAR(t[w[0]]):.2f}" for w in picks],
    )
    for k, w in enumerate(picks):
        g = lc.iloc[w]
        x = (g.bjd_tdb.values - g.bjd_tdb.values[np.argmin(g.model.values)]) * 24
        fig.add_trace(
            go.Scatter(
                x=x,
                y=g.flux,
                mode="markers",
                marker=dict(size=4, color=C_DATA),
                name="TESS",
                showlegend=k == 0,
            ),
            1,
            k + 1,
        )
        fig.add_trace(
            go.Scatter(
                x=x,
                y=g.model,
                mode="lines",
                line=dict(color=C_Q, width=2),
                name=model_name,
                showlegend=k == 0,
            ),
            1,
            k + 1,
        )
        fig.add_trace(
            go.Scatter(
                x=x,
                y=1e3 * g.o_minus_model,
                mode="markers",
                marker=dict(size=3.5, color=C_DATA),
                showlegend=False,
            ),
            2,
            k + 1,
        )
        fig.update_xaxes(title_text="hours from minimum", row=2, col=k + 1)
    fig.update_yaxes(title_text="relative flux", row=1, col=1)
    fig.update_yaxes(title_text="O − model [ppt]", row=2, col=1)
    _layout(fig, h=500, w=1250)
    fig.update_layout(legend=dict(y=1.09), margin=dict(t=90))
    return fig


def rv_figure(rv, P, primaries):
    fig = make_subplots(
        rows=2,
        cols=2,
        shared_xaxes=True,
        row_heights=[0.7, 0.3],
        vertical_spacing=0.06,
        subplot_titles=("star 1 (m1)", "star 2 (m2)"),
    )
    # phase from the model's own primary eclipses (light-time effect and precession shift the phase between 1983
    # and 2022, so one linear ephemeris does not fold both epochs)
    cols = {"CHIRON2022": C_Q, "Lacy1999": C_ART}
    tp = np.sort(np.asarray(primaries, float))

    def phase(t):
        out = np.empty(len(t))
        for n_, ti in enumerate(t):
            near = tp[
                np.argsort(np.abs(tp - ti))[:60]
            ]  # a local linear ephemeris from the nearest model primaries
            E = np.round((near - near.min()) / P)
            Pl, T0l = np.polyfit(E, near, 1)
            out[n_] = ((ti - T0l) / Pl) % 1.0
        return out

    for k, st in enumerate(("m1", "m2")):
        g = rv[(rv.star == st) & (rv.status == "fitted")]
        ph = pd.Series(phase(g.bjd_tdb.values), index=g.index)
        fig.add_trace(
            go.Scatter(
                x=ph,
                y=g.model,
                mode="markers",
                marker=dict(symbol="line-ew-open", size=10, color="#999"),
                name="model at the same times",
                showlegend=k == 0,
            ),
            1,
            k + 1,
        )
        for src, c in cols.items():
            h = g.source == src
            fig.add_trace(
                go.Scatter(
                    x=ph[h],
                    y=g.rv[h],
                    mode="markers",
                    marker=dict(size=6, color=c),
                    name={
                        "CHIRON2022": "CHIRON 2022 (Zasche et al. 2023)",
                        "Lacy1999": "1983–1994 (Lacy et al. 1999)",
                    }[src],
                    showlegend=k == 0,
                ),
                1,
                k + 1,
            )
            fig.add_trace(
                go.Scatter(
                    x=ph[h],
                    y=g.o_minus_model[h],
                    mode="markers",
                    marker=dict(size=5, color=c),
                    error_y=dict(type="data", array=g.sigma[h], thickness=0.7, width=0),
                    showlegend=False,
                ),
                2,
                k + 1,
            )
        fig.update_xaxes(title_text="phase of the 3.78-day orbit", row=2, col=k + 1)
    fig.update_yaxes(title_text="radial velocity [km/s]", row=1, col=1)
    fig.update_yaxes(title_text="O − model", row=2, col=1)
    _layout(fig, h=540, w=1150)
    fig.update_layout(legend=dict(y=1.09), margin=dict(t=90))
    return fig


# ------------------------------------------------------------------ QX Cas
def qx_cas(root: Path, out_root: Path, O: Out, N: dict):
    s, o = "qx_cas", out_root / "qx_cas"
    cache = o / "report"
    fit_f = o / "fit_v2/stage2_fitted.yaml"
    fit = _load(fit_f)
    t0 = fit_epoch(fit_f)
    pub = _load(
        root / "configs/squadra_test3/qx_cas_published.yaml"
    )  # the binary of Guinan et al. (2012)
    n = N.setdefault(s, {})
    m_fit = cached(cache, "sys_fit_v2", lambda: system(fit), fit)
    art = {
        "topology": "custom",
        "hierarchy": {"inner": ["m1", "m2"]},
        "masses": {k: pub["masses"][k] for k in ("m1", "m2")},
        "orbits": {"inner": dict(fit["orbits"]["inner"], i=pub["orbits"]["inner"]["i"])},
        "stars": {k: fit["stars"][k] for k in ("m1", "m2")},
        "stellar_model": fit.get("stellar_model"),
        "physics": fit.get("physics"),
        "run": fit.get("run"),
    }
    m_art = moved(cached(cache, "sys_article", lambda: system(art), art), m_fit)
    yrs = np.arange(1900.0, 1900.0 + 420.0, 2.0)
    cone_fit = cached(cache, "axis_fit_v2", lambda: axis_track(fit, t0, yrs), fit, t0, yrs)
    fig = views_figure(
        [
            ("binary of Guinan et al. (2012)", m_art, C_ART, True),
            ("SWINGERsim fit", m_fit, C_Q, False),
        ],
        cone=[("SWINGERsim fit", cone_fit, C_Q, False)],
        i_lim=m_fit["pairs"][0]["i_limit_deg"],
        cone_note="axis 1900-2320 (three turns of ~127 yr); yellow: eclipses",
        label_years=(1965, 2000, 2029, 2060),
    )
    O.save(fig, "qx_models")
    n["nodal_period_fit_yr"] = nodal_period(yrs, cone_fit["h"], np.array(cone_fit["L"])[0])
    n["system"] = {
        "triples": m_fit["triples"],
        "i_mut": m_fit["mutual"][0]["deg"],
        "pairs": m_fit["pairs"],
        "orbits": [
            {k: x[k] for k in ("name", "P_d", "a_au", "e", "i_deg")} for x in m_fit["orbits"]
        ],
    }
    years = np.arange(1900.0, 2101.0, 0.5)
    c_fit = cached(cache, "inc_fit_v2", lambda: inc_curve(fit, t0, years), fit, t0, years)
    lim = np.array(c_fit["limit"])
    n["eclipse_changes"] = eclipse_switches(years, c_fit["i"], lim)
    inc = read_table(root / "data/observed/qx_cas/inclination_limits.csv")
    fig = go.Figure()
    eclipse_band(fig, years, lim)
    fig.add_trace(
        go.Scatter(
            x=years,
            y=unfold(c_fit["i"], years, 1980.0),
            mode="lines",
            line=dict(color=C_Q, width=2.5),
            name="SWINGERsim fit",
        )
    )
    meas = inc[inc.limit.isna()]
    fig.add_trace(
        go.Scatter(
            x=YEAR(meas.mid_hjd),
            y=meas.i_deg,
            mode="markers",
            marker=dict(size=9, color=C_DATA),
            error_y=dict(type="data", array=meas.sigma_deg),
            name="light-curve solutions (Guinan et al. 2012)",
        )
    )
    up = inc[inc.limit == "upper"]
    fig.add_trace(
        go.Scatter(
            x=YEAR(up.mid_hjd),
            y=up.i_deg,
            mode="markers",
            marker=dict(size=13, symbol="triangle-down-open", color=C_DATA, line=dict(width=2)),
            name="upper limit 2003–2008 (Guinan et al. 2012)",
        )
    )
    ne = inc[inc.limit == "no_eclipse"]
    fig.add_trace(
        go.Scatter(
            x=YEAR(ne.mid_hjd),
            y=np.interp(YEAR(ne.mid_hjd), years, lim),
            mode="markers",
            marker=dict(
                size=11, symbol="x-thin", color="#B3261E", line=dict(width=2.5, color="#B3261E")
            ),
            name="no eclipses seen (FCAPT 2003–2008, TESS 2019–2024)",
        )
    )
    ec = inc[inc.limit == "eclipse"]
    fig.add_trace(
        go.Scatter(
            x=YEAR(ec.mid_hjd),
            y=np.interp(YEAR(ec.mid_hjd), years, lim),
            mode="markers",
            marker=dict(
                size=11, symbol="cross-thin", color="#2E7D4F", line=dict(width=2.5, color="#2E7D4F")
            ),
            name="photographic minima 1932–1935 (fit: 'the pair eclipsed', geometrically)",
        )
    )
    fig.update_xaxes(title_text="year")
    fig.update_yaxes(title_text="inclination of the 6.0-day pair [°]", range=[45, 95])
    O.save(_layout(fig, h=460), "qx_inclination")
    ex = o / "extra"
    if (ex / "depth_history_fit.csv").exists():
        O.save(qx_depth(ex, fit, n), "qx_depth")
    # ensemble: weights of the observed values of its configuration (swingersim.ensemble.observation_weights)
    from swingersim.ensemble import effective_size

    summary = o / "ensemble/summary.csv"
    d = pd.read_csv(summary)
    f = d[d.direction == "forward"].set_index("member")
    b = d[d.direction == "backward"].set_index("member")
    m = f.join(b, rsuffix="_b", how="inner")
    ens_cfg = o / "ensemble/ensemble.yaml"
    ens_cfg = ens_cfg if ens_cfg.exists() else root / "configs/squadra_test3/qx_cas_ensemble.yaml"
    w = ensemble_weights(ensemble_rows(summary), ensemble_observed(ens_cfg), m.index.values)

    def wq(x, q):
        o_ = np.argsort(x)
        c = np.cumsum(w[o_])
        return float(np.interp(q, c / c[-1], x[o_]))

    P = m["orbits.out.P"].values
    M3 = m["masses.m3"].values
    n["ensemble"] = {
        "n": len(m),
        "n_eff": effective_size(w),
        "P3_d": [wq(P, q) for q in (0.16, 0.5, 0.84)],
        "m3": [wq(M3, q) for q in (0.16, 0.5, 0.84)],
        "J": [wq(m.imut_out_inner_0.values, q) for q in (0.16, 0.5, 0.84)],
        "frac_P3_100_1000": float(w[(P > 100) & (P < 1000)].sum() / w.sum()),
        "collisions_all": float(((m.fate != "intact") | (m.fate_b != "intact")).mean()),
    }
    fig = make_subplots(
        rows=1, cols=2, subplot_titles=("period of the third star", "mass of the third star")
    )
    fig.add_annotation(
        text="dashed: the SWINGERsim fit",
        xref="paper",
        yref="paper",
        x=1.0,
        y=1.13,
        showarrow=False,
        font=dict(size=11, color=C_ART),
    )
    for k, (x, bins, lab) in enumerate(
        (
            (P, np.logspace(np.log10(30), np.log10(11000), 31), "days"),
            (M3, np.logspace(np.log10(0.3), np.log10(6), 21), "M☉"),
        )
    ):
        hist_all, e = np.histogram(x, bins=bins)
        hist_w, _ = np.histogram(x, bins=bins, weights=w * len(x) / w.sum())
        c = np.sqrt(e[1:] * e[:-1])
        fig.add_trace(
            go.Bar(
                x=c,
                y=hist_all,
                width=np.diff(e),
                marker_color="#d5dbe6",
                name="all 2000 tried systems",
                showlegend=k == 0,
            ),
            1,
            k + 1,
        )
        fig.add_trace(
            go.Scatter(
                x=c,
                y=hist_w,
                mode="lines",
                line=dict(color=C_Q, width=2.5, shape="hvh"),
                name="weighted by the measured i and di/dt",
                showlegend=k == 0,
            ),
            1,
            k + 1,
        )
        tv = [30, 100, 300, 1000, 3000, 10000] if k == 0 else [0.3, 0.5, 1, 2, 3, 5]
        fig.update_xaxes(
            type="log",
            title_text=lab,
            tickvals=tv,
            ticktext=[f"{v:g}" for v in tv],
            row=1,
            col=k + 1,
        )
    fig.add_vline(x=fit["orbits"]["out"]["P"], line=dict(color=C_ART, dash="dash"), row=1, col=1)
    fig.add_vline(x=fit["masses"]["m3"], line=dict(color=C_ART, dash="dash"), row=1, col=2)
    fig.update_yaxes(title_text="systems", row=1, col=1)
    _layout(fig, h=440)
    fig.update_layout(legend=dict(y=1.1), margin=dict(t=90))
    O.save(fig, "qx_ensemble")
    # TESS: the binary seen without eclipses
    tess = qx_tess(root, out_root, cache, fit, t0, n)
    O.save(tess, "qx_tess")
    res = pd.read_csv(o / "fit_v2/stage2_residuals.csv")
    n["minima_rms_by_class"] = {
        k: float(np.sqrt(np.mean(g.o_minus_model_min**2)))
        for k, g in res[res.status == "fitted"].groupby("source")
    }
    p = pd.read_csv(o / "fit_v2/stage2_params.csv")
    n["params_fit"] = {r.parameter: [float(r.value), float(r.error)] for r in p.itertuples()}


def qx_depth(ex: Path, fit: dict, n: dict):
    """QX Cas: was a model eclipse deep enough to be seen? Left: depth of the primary eclipse in B against the
    inclination (the fitted stars; band: all orientations of the eccentric orbit); right: the depth along the fitted
    inclination history, for several third lights, with the detection limits of photographic and photoelectric
    photometry (tools/qx_cas_extra.py)."""
    dv = pd.read_csv(ex / "depth_vs_inclination.csv")
    h = pd.read_csv(ex / "depth_history_fit.csv")
    sm = json.loads((ex / "summary.json").read_text())
    n["depth"] = {k: sm[k] for k in ("detectable_inclination_deg", "third_light_needed")}
    n["depth"]["fit"] = sm["fit"]
    n["depth"]["ensemble"] = sm.get("ensemble")
    mag = lambda d: -2.5 * np.log10(1.0 - np.asarray(d, float))
    fig = make_subplots(
        rows=1,
        cols=2,
        column_widths=[0.4, 0.6],
        horizontal_spacing=0.1,
        subplot_titles=(
            "primary eclipse depth in B against the inclination",
            "depth along the fitted inclination, 1900–2100",
        ),
    )
    b = dv[dv.band == "B"]
    w_fit = float(fit["orbits"]["inner"].get("omega", 0.0))
    g = b.groupby("i").drop_primary
    fig.add_trace(
        go.Scatter(
            x=np.r_[g.min().index, g.max().index[::-1]],
            y=np.r_[mag(g.min().values), mag(g.max().values[::-1])],
            fill="toself",
            fillcolor="rgba(31,79,176,0.18)",
            line=dict(width=0),
            name="all orientations ω of the eccentric orbit",
            hoverinfo="skip",
        ),
        1,
        1,
    )
    bf = b[np.isclose(b.omega, w_fit)].sort_values("i")
    fig.add_trace(
        go.Scatter(
            x=bf.i,
            y=mag(bf.drop_primary),
            mode="lines",
            line=dict(color=C_Q, width=2.5),
            name=f"ω of the fit ({w_fit % 360:.0f}°), no third light",
        ),
        1,
        1,
    )
    cols = {"0.00": C_Q, "0.10": "#5B7FD0", "0.20": "#8EA6E0", "0.37": "#B9C7EC"}
    for k, c in cols.items():
        fig.add_trace(
            go.Scatter(
                x=h.year,
                y=h[f"depth_B_l3_{k}"],
                mode="lines",
                line=dict(color=c, width=2.5 if k == "0.00" else 1.8),
                name=f"third light {float(k):.0%} of the total",
                legendgroup="l3",
            ),
            1,
            2,
        )
    for thr, txt in (
        (0.10, "≈ limit of photographic plates (0.10 mag)"),
        (0.02, "≈ limit of photoelectric photometry (0.02 mag)"),
    ):
        for c_ in (1, 2):
            fig.add_hline(y=thr, line=dict(color="#666", width=1, dash="dash"), row=1, col=c_)
        fig.add_annotation(
            x=2100,
            y=thr,
            text=txt,
            showarrow=False,
            xanchor="right",
            yanchor="bottom",
            font=dict(size=10, color="#555"),
            row=1,
            col=2,
        )
    fig.add_vrect(x0=1933, x1=1974, fillcolor="rgba(46,125,79,0.10)", line_width=0, row=1, col=2)
    fig.add_annotation(
        x=1953,
        y=0.62,
        text="photographic minima<br>1933–1974 (VarAstro)",
        showarrow=False,
        font=dict(size=10, color="#2E5D3F"),
        row=1,
        col=2,
    )
    fig.update_xaxes(title_text="inclination of the 6.0-day pair [°]", row=1, col=1)
    fig.update_xaxes(title_text="year", range=[1900, 2100], row=1, col=2)
    for c_ in (1, 2):
        fig.update_yaxes(
            title_text="eclipse depth in B [mag]" if c_ == 1 else None,
            range=[0.68, -0.01],
            row=1,
            col=c_,
        )
    _layout(fig, h=470)
    fig.update_layout(legend=dict(y=1.12), margin=dict(t=120))
    return fig


def qx_tess(root, out_root, cache, fit, t0, n):
    """Search the TESS light curves of QX Cas for the light variation of the pair that swingersim predicts without
    eclipses (tidal distortion and reflection near periastron): harmonics of the orbital period."""
    from astropy.timeseries import LombScargle

    from swingersim.photometry import light_curve

    D = _data(root, "qx_cas", "tess_s*_spoc.csv")
    files = [Path(f) for f in sorted(glob.glob(str(D / "tess_s*_spoc.csv")))]
    ts, fs = [], []
    for f in files:
        d = read_table(f).dropna()
        tb, fb = time_bins(d.bjd_tdb.values, d.flux.values, 48)  # 30-min bins
        ts.append(tb)
        fs.append(fb / np.median(fb))
    if not ts:
        return go.Figure()
    t = np.concatenate(ts)
    f = np.concatenate(fs)
    mm = np.array(
        cached(
            cache,
            "qx_tess_model",
            lambda: light_curve(
                _cfg(fit), t - t0, "TESS", exposure_min=30.0, n_sub=3
            ).flux.tolist(),
            fit,
            t0,
            *files,
        )
    )
    seg = np.concatenate([[0], np.cumsum(np.diff(t) > 1)])

    def detrend(y):
        y = y.copy()
        for s_ in np.unique(seg):
            k = seg == s_
            tc = t[k] - t[k].mean()
            y[k] -= np.polyval(np.polyfit(tc, y[k], 2), tc)
        return y

    y = detrend(f - 1)
    m = detrend(mm - 1)
    fr = np.linspace(1 / 15, 1 / 0.3, 60000)
    P = 1 / fr
    removed = []
    for _ in range(
        6
    ):  # the strongest signals that are not harmonics of the binary (other stars in the aperture)
        p = LombScargle(t, y).power(fr)
        ok = (
            (np.abs(P - 6.0) > 0.3)
            & (np.abs(P - 3.0) > 0.15)
            & (np.abs(P - 2.0) > 0.08)
            & (np.abs(P - 1.5) > 0.05)
        )
        Pk = float(P[np.argmax(np.where(ok, p, 0))])
        removed.append(Pk)
        X = np.column_stack(
            [fn(2 * np.pi * h * t / Pk) for h in (1, 2) for fn in (np.sin, np.cos)]
            + [np.ones_like(t)]
        )
        y = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
        m = m - X @ np.linalg.lstsq(X, m, rcond=None)[0]
    p1 = LombScargle(t, y).power(fr)
    H = 4
    # the period of the pair's light variation in 2019-2024 as the fitted model predicts it (the mean apparent
    # period of the model light curve near the fitted P, which is an osculating period at t = 0)
    Pb = model_period(t, m, float(fit["orbits"]["inner"]["P"]), H)
    n["binary_period_model_d"] = Pb
    X = np.column_stack(
        [fn(2 * np.pi * h * t / Pb) for h in range(1, H + 1) for fn in (np.sin, np.cos)]
        + [np.ones_like(t)]
    )
    cd = np.linalg.lstsq(X, y, rcond=None)[0]
    cm = np.linalg.lstsq(X, m, rcond=None)[0]
    resid = y - X @ cd
    rng = np.random.default_rng(1)
    boot = []
    for _ in range(300):
        c = np.linalg.lstsq(
            X, X @ cd + np.roll(resid, int(rng.integers(200, len(resid) - 200))), rcond=None
        )[0]
        boot.append([math.hypot(c[2 * k], c[2 * k + 1]) for k in range(H)])
    boot = np.array(boot)
    rnd = []
    for Pr in np.random.default_rng(2).uniform(2.2, 12, 400):
        if min(abs(Pr - Pb / k) for k in (1, 2, 3, 4)) < 0.05:
            continue
        Xr = np.column_stack(
            [np.sin(2 * np.pi * t / Pr), np.cos(2 * np.pi * t / Pr), np.ones_like(t)]
        )
        c = np.linalg.lstsq(Xr, y, rcond=None)[0]
        rnd.append(math.hypot(c[0], c[1]))
    n["tess"] = {
        "removed_periods_d": removed,
        "noise_99_ppt": float(np.percentile(rnd, 99) * 1e3),
        "noise_median_ppt": float(np.median(rnd) * 1e3),
        "harmonics": [
            {
                "k": k + 1,
                "P_d": Pb / (k + 1),
                "data_ppt": 1e3 * math.hypot(cd[2 * k], cd[2 * k + 1]),
                "err_ppt": 1e3 * float(np.std(boot[:, k])),
                "model_ppt": 1e3 * math.hypot(cm[2 * k], cm[2 * k + 1]),
            }
            for k in range(H)
        ],
        "points": len(t),
        "sectors": len(ts),
    }
    fig = make_subplots(
        rows=1,
        cols=2,
        column_widths=[0.55, 0.45],
        subplot_titles=(
            "periodogram of TESS 2019–2024 (other stars' signals removed)",
            "folded on the binary period: TESS and SWINGERsim",
        ),
    )
    fig.add_trace(
        go.Scatter(
            x=P, y=p1, mode="lines", line=dict(color=C_DATA, width=1), name="TESS periodogram"
        ),
        1,
        1,
    )
    for k in range(1, H + 1):
        fig.add_vline(x=Pb / k, line=dict(color=C_ART, width=1, dash="dash"), row=1, col=1)
    fig.add_annotation(
        text="dashed: the orbital period P and P/2, P/3, P/4",
        xref="x domain",
        yref="y domain",
        x=0.98,
        y=0.96,
        showarrow=False,
        font=dict(size=10),
        row=1,
        col=1,
    )
    fig.update_xaxes(
        type="log", range=[np.log10(1.2), np.log10(12)], title_text="period [days]", row=1, col=1
    )
    fig.update_yaxes(title_text="power", row=1, col=1)
    ph = (t / Pb) % 1.0
    bins = np.linspace(0, 1, 41)
    idx = np.digitize(ph, bins) - 1
    cen = 0.5 * (bins[1:] + bins[:-1])
    yb = np.array([np.mean(y[idx == k]) for k in range(40)])
    eb = np.array(
        [np.std(y[idx == k]) / math.sqrt(max((idx == k).sum() / 12, 1)) for k in range(40)]
    )
    mb = np.array([np.mean(m[idx == k]) for k in range(40)])
    fig.add_trace(
        go.Scatter(
            x=cen,
            y=1e3 * yb,
            mode="markers",
            marker=dict(size=7, color=C_DATA),
            error_y=dict(type="data", array=1e3 * eb, thickness=1),
            name="TESS (binned)",
        ),
        1,
        2,
    )
    fig.add_trace(
        go.Scatter(
            x=cen,
            y=1e3 * mb,
            mode="lines",
            line=dict(color=C_Q, width=2.5),
            name="SWINGERsim prediction (same processing)",
        ),
        1,
        2,
    )
    fig.update_xaxes(title_text=f"phase of the {Pb:.4f}-day period", row=1, col=2)
    fig.update_yaxes(title_text="brightness change [ppt]", row=1, col=2)
    _layout(fig, h=460, w=1250)
    fig.update_layout(legend=dict(y=1.12), margin=dict(t=95))
    return fig


def model_period(t, m, p_guess: float, harmonics: int = 4, rel: float = 2e-3) -> float:
    """Period that fits the model light curve m(t) best with `harmonics` harmonics, searched within +-rel of
    p_guess (the residual of a linear least-squares fit, minimised on a grid and refined around the best point)."""

    def resid(P):
        X = np.column_stack(
            [fn(2 * np.pi * h * t / P) for h in range(1, harmonics + 1) for fn in (np.sin, np.cos)]
            + [np.ones_like(t)]
        )
        return float(np.sum((m - X @ np.linalg.lstsq(X, m, rcond=None)[0]) ** 2))

    grid = p_guess * (1.0 + np.linspace(-rel, rel, 401))
    k = int(np.argmin([resid(P) for P in grid]))
    fine = np.linspace(grid[max(k - 1, 0)], grid[min(k + 1, len(grid) - 1)], 101)
    return float(fine[int(np.argmin([resid(P) for P in fine]))])


# ------------------------------------------------------------------ AH Cep
def ah_cep(root: Path, out_root: Path, O: Out, N: dict):
    s, o = "ah_cep", out_root / "ah_cep"
    cache = o / "report"
    fit_f = o / "fit/stage2_fitted.yaml"
    fit = _load(fit_f)
    t0 = fit_epoch(fit_f)
    fit["orbits"]["out"]["Omega"] = (
        fit["orbits"]["out"]["Omega"] - fit["orbits"]["inner"]["Omega"]
    ) % 360.0
    fit["orbits"]["inner"]["Omega"] = 0.0  # only the node difference is defined
    fit["stars"]["m3"] = {"R": 3.23, "Teff": 21050}  # MIST 7 Msun at 6 Myr
    n = N.setdefault(s, {})
    art = {
        "topology": "custom",
        "hierarchy": {"inner": ["m1", "m2"]},
        "masses": {"m1": 16.14, "m2": 13.69},
        "orbits": {
            "inner": {
                "P": 1.7747420,
                "e": 0.0,
                "i": 69.81,
                "Omega": 0.0,
                "omega": 90.0,
                "M0": 15.89,
            }
        },
        "stars": {"m1": {"R": 6.51, "Teff": 30700}, "m2": {"R": 5.64, "Teff": 28800}},
        "stellar_model": fit.get("stellar_model"),
        "physics": fit.get("physics"),
        "run": fit.get("run"),
    }
    m_fit = cached(cache, "sys_fit", lambda: system(fit), fit)
    m_art = moved(cached(cache, "sys_article", lambda: system(art), art), m_fit)
    yrs = np.arange(1900.0, 2301.0, 5.0)
    cone_fit = cached(cache, "axis_fit", lambda: axis_track(fit, t0, yrs), fit, t0, yrs)
    fig = views_figure(
        [
            ("binary of Pavlovski et al. (2018)", m_art, C_ART, True),
            ("SWINGERsim fit", m_fit, C_Q, False),
        ],
        cone=[("SWINGERsim fit", cone_fit, C_Q, False)],
        i_lim=m_fit["pairs"][0]["i_limit_deg"],
        cone_note="the pair's axis 1900-2300: it does not move visibly",
        label_years=(),
        sky_pad=1.0,
    )
    O.save(fig, "ah_models")
    n["system"] = {
        "triples": m_fit["triples"],
        "i_mut": m_fit["mutual"][0]["deg"],
        "pairs": m_fit["pairs"],
        "orbits": [
            {k: x[k] for k in ("name", "P_d", "a_au", "e", "i_deg")} for x in m_fit["orbits"]
        ],
    }
    h = np.array(cone_fit["h"])
    i_t = np.degrees(np.arccos(np.abs(h[:, 2])))
    n["i_range_1900_2300"] = [float(i_t.min()), float(i_t.max())]
    res = pd.read_csv(o / "fit/stage2_residuals.csv")
    oc = cached(
        cache, "oc_fit", lambda: oc_model(fit, t0, 1900.5, 2026.3), fit, t0, (1900.5, 2026.3)
    )
    res2 = res.copy()
    res2["source"] = res2.source.str.replace(r"TESS S\d+", "TESS", regex=True)
    O.save(oc_figure(res2, oc, "inner", "times of minima 1901–2026"), "ah_oc")
    r = res2[res2.status == "fitted"]
    n["minima_rms_by_class"] = {
        k: float(np.sqrt(np.mean(g.o_minus_model_min**2))) for k, g in r.groupby("source")
    }
    n["minima_n_by_class"] = {k: len(g) for k, g in r.groupby("source")}
    p = pd.read_csv(o / "fit/stage2_params.csv")
    n["params_fit"] = {rr.parameter: [float(rr.value), float(rr.error)] for rr in p.itertuples()}
    # the TESS light-curve fit (inclination, radii and T2 free; the depth of every segment scaled by a free factor a)
    rp = o / "fit_lc" / "stage2_report.txt"
    if rp.exists():
        import re

        txt = rp.read_text()
        a_ = [float(x) for x in re.findall(r"a = ([0-9.]+)", txt)]
        fl = _load(o / "fit_lc" / "stage2_fitted.yaml")
        pq = pd.read_csv(o / "fit_lc" / "stage2_params.csv").set_index("parameter")
        st = fl["stars"]

        def L(k, band):
            return band_luminosity(st[k]["R"], st[k]["Teff"], band)

        l3m = L("m3", "TESS") / (L("m1", "TESS") + L("m2", "TESS") + L("m3", "TESS"))
        l3 = [
            1 - a * (1 - l3m) for a in a_
        ]  # observed depth = a x model depth; the model has the third star's light l3m
        sin3 = (
            math.sin(math.radians(69.81)) / math.sin(math.radians(fl["orbits"]["inner"]["i"]))
        ) ** 3
        n["fit_lc_derived"] = {
            "i": float(fl["orbits"]["inner"]["i"]),
            "i_err": float(pq.loc["orbits.inner.i", "error"]),
            "R1": float(st["m1"]["R"]),
            "R1_err": float(pq.loc["stars.m1.R", "error"]),
            "R2": float(st["m2"]["R"]),
            "R2_err": float(pq.loc["stars.m2.R", "error"]),
            "T2": float(st["m2"]["Teff"]),
            "T2_err": float(pq.loc["stars.m2.Teff", "error"]),
            "scale_a": a_,
            "l3_model": float(l3m),
            "l3_eff": l3,
            "l3_eff_mean": float(np.mean(l3)),
            "l3_eff_spread": float(0.5 * (max(l3) - min(l3))),
            "light_ratio_B": float(L("m2", "B") / L("m1", "B")),
            "light_ratio_B_published": band_luminosity(5.64, 28800, "B")
            / band_luminosity(6.51, 30700, "B"),
            "masses_scaled_to_i": [16.14 * sin3, 13.69 * sin3],
            "mass_factor": float(sin3),
        }
        n["_fit_lc_table"] = pd.read_csv(o / "fit_lc" / "stage2_lc.csv")
    # third light in TESS
    O.save(ah_third_light(root, cache, fit, t0, n), "ah_third_light")
    n.pop("_fit_lc_table", None)
    # control jobs (when present)
    e4 = o / "ensemble_4body" / "summary.csv"
    if e4.exists():
        O.save(
            ah_fourth(
                e4,
                fit,
                t0,
                n,
                o / "fit" / "stage2_residuals.csv",
                rate_limit=round(n.get("tess_constancy", {}).get("rate_deg_yr", 0.06), 2),
            ),
            "ah_fourth_body",
        )
    rp = o / "fit_lc" / "stage2_report.txt"
    if rp.exists():
        n["fit_lc"] = rp.read_text()[:4000]
        pp = o / "fit_lc" / "stage2_params.csv"
        if pp.exists():
            pq = pd.read_csv(pp)
            n["params_fit_lc"] = {
                rr.parameter: [float(rr.value), float(rr.error)] for rr in pq.itertuples()
            }


def ah_third_light(root, cache, fit, t0, n):
    """TESS eclipse depths of AH Cep against the published pair: for every inclination the extra (third) light
    that the calibration of the depths needs, with the published radii and temperatures."""
    from swingersim.photometry import light_curve
    from swingersim.shapes import apply_policy

    D = _data(root, "ah_cep", "tess_s*_spoc.csv")
    sec, files = {}, []
    for s in ("s77", "s78", "s18"):
        f = D / f"tess_{s}_spoc.csv"
        if not f.exists():
            continue
        files.append(f)
        x = read_table(f).dropna()
        tb, fb = time_bins(x.bjd_tdb.values, x.flux.values, 144)  # 10-min bins
        sec[s] = (tb, fb / np.median(fb))
    pair = copy.deepcopy(fit)
    pair["stars"]["m3"] = {
        "R": 0.01,
        "Teff": 3500,
    }  # the pair alone; the third light is the free constant
    scan = []
    keep_curves = {}
    for inc in (69.8, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0):
        key = f"tl_{inc:.1f}"

        def compute(inc=inc):
            d = copy.deepcopy(pair)
            d["orbits"]["inner"]["i"] = inc
            d, _ = apply_policy(d, "auto", 0.5)
            cfg = _cfg(d)
            return {
                s: light_curve(cfg, t - t0, "TESS", exposure_min=10.0, n_sub=3).flux.tolist()
                for s, (t, y) in sec.items()
            }

        mods = cached(cache, key, compute, pair, inc, t0, *files)
        row = {"i": inc}
        for s, (t, y) in sec.items():
            F = np.array(mods[s])
            seg = np.concatenate([[0], np.cumsum(np.diff(t) > 0.5)])
            cols = [F]
            for k in np.unique(seg):
                mm = (seg == k).astype(float)
                cols += [mm, (t - t[seg == k].mean()) * mm]
            X = np.column_stack(cols)
            keep = np.ones(len(t), bool)
            for _ in range(3):
                c = np.linalg.lstsq(X[keep], y[keep], rcond=None)[0]
                r = y - X @ c
                sd = 1.4826 * np.median(abs(r[keep]))
                keep = abs(r) < 5 * sd
            const = np.mean([c[1 + 2 * j] for j in range(len(np.unique(seg)))])
            row[s] = {
                "l3": float(const / (const + c[0] * F.max())),
                "rms_ppt": float(np.std(r[keep]) * 1e3),
            }
            if s == "s77":
                keep_curves[inc] = (t, y, X @ c, F)
        scan.append(row)
    n["third_light_scan"] = scan
    if all(k in scan[0] for k in ("s77", "s78", "s18")):
        ii_ = np.array([r_["i"] for r_ in scan])
        k74 = int(np.argmin(np.abs(ii_ - 74.0)))
        slope = np.polyfit(ii_[2:], [r_["s77"]["l3"] for r_ in scan][2:], 1)[
            0
        ]  # extra light per degree near 72-76
        dl = scan[k74]["s18"]["l3"] - 0.5 * (scan[k74]["s77"]["l3"] + scan[k74]["s78"]["l3"])
        dt_ = float(np.mean(YEAR(sec["s77"][0])) - np.mean(YEAR(sec["s18"][0])))
        n["tess_constancy"] = {
            "slope_l3_per_deg": float(slope),
            "delta_l3_s18_minus_s77s78": float(dl),
            "delta_i_deg": float(dl / slope),
            "years": dt_,
            "rate_deg_yr": float(abs(dl / slope) / dt_),
        }
    # the light of a 7-Msun main-sequence star of the age of the pair (MIST, 6 Myr) and of two 3.4-Msun stars
    from swingersim.mist import mist_track

    def lum(mass):
        st = mist_track(mass, 0.0, 0.0).at(6.0)
        return band_luminosity(st.radius_rsun, st.teff_k, "TESS")

    pairL = lum(16.14) + lum(13.69)
    n["third_light_expected"] = {
        "one_7Msun": lum(7.0) / (pairL + lum(7.0)),
        "one_6.6Msun": lum(6.6) / (pairL + lum(6.6)),
        "two_3.5_3.3Msun": (lum(3.5) + lum(3.3)) / (pairL + lum(3.5) + lum(3.3)),
    }
    fig = make_subplots(
        rows=2,
        cols=2,
        column_widths=[0.56, 0.44],
        row_heights=[0.72, 0.28],
        vertical_spacing=0.07,
        horizontal_spacing=0.1,
        specs=[[{}, {"rowspan": 2}], [{}, None]],
        subplot_titles=(
            "TESS sector 77 (2024), one orbit: data and fits",
            "extra light the TESS eclipses need",
        ),
    )
    lcf = n.get("_fit_lc_table")
    if lcf is not None:
        g = lcf[lcf.segment == "S77a"].sort_values("bjd_tdb")
        tt = g.bjd_tdb.values
        w = (tt > tt[0] + 1.0) & (tt < tt[0] + 1.0 + 1.78)
        x = (tt[w] - tt[w][0]) * 24
        x = np.where(
            np.r_[False, np.diff(x) > 0.5], np.nan, x
        )  # break the lines at gaps in the data
        t, y, _, F = keep_curves[69.8]
        Fw = np.interp(
            tt[w], t, F
        )  # the published pair as it is: only a zero point and a slope fitted, no extra light
        Xw = np.column_stack([np.ones(w.sum()), tt[w] - tt[w].mean()])
        m_pub = Fw + Xw @ np.linalg.lstsq(Xw, g.flux.values[w] - Fw, rcond=None)[0]
        fig.add_trace(
            go.Scatter(
                x=x,
                y=g.flux.values[w],
                mode="markers",
                marker=dict(size=4, color=C_DATA),
                name="TESS",
            ),
            1,
            1,
        )
        fig.add_trace(
            go.Scatter(
                x=x,
                y=m_pub,
                mode="lines",
                line=dict(color=C_ART, width=2, dash="dash"),
                name="published pair, i = 69.8°, no extra light",
            ),
            1,
            1,
        )
        fig.add_trace(
            go.Scatter(
                x=x,
                y=g.model.values[w],
                mode="lines",
                line=dict(color=C_Q, width=2),
                name="SWINGERsim fit: i, radii, T₂ free (Roche shapes)",
            ),
            1,
            1,
        )
        fig.add_trace(
            go.Scatter(
                x=x,
                y=1e3 * (g.flux.values[w] - m_pub),
                mode="markers",
                marker=dict(size=3, color=C_ART),
                showlegend=False,
            ),
            2,
            1,
        )
        fig.add_trace(
            go.Scatter(
                x=x,
                y=1e3 * g.o_minus_model.values[w],
                mode="markers",
                marker=dict(size=3, color=C_Q),
                showlegend=False,
            ),
            2,
            1,
        )
        fig.add_hline(y=0, line=dict(color="#888", width=1), row=2, col=1)
        fig.update_yaxes(title_text="O − model [ppt]", row=2, col=1)
        fig.update_xaxes(title_text="hours", row=2, col=1)
    fig.update_yaxes(title_text="relative flux", row=1, col=1)
    ii = [r["i"] for r in scan]
    for s, dash in (("s77", "solid"), ("s78", "dash"), ("s18", "dot")):
        if s in scan[0]:
            fig.add_trace(
                go.Scatter(
                    x=ii,
                    y=[100 * r[s]["l3"] for r in scan],
                    mode="lines+markers",
                    line=dict(color=C_Q, dash=dash),
                    name=f"sector {s[1:]}",
                    text=[f"rms {r[s]['rms_ppt']:.2f} ppt" for r in scan],
                    hovertemplate="i %{x}°: %{y:.1f} %, %{text}<extra></extra>",
                ),
                1,
                2,
            )
    e = n["third_light_expected"]
    fig.add_hrect(
        y0=100 * e["one_6.6Msun"],
        y1=100 * e["one_7Msun"],
        fillcolor="rgba(194,87,26,0.18)",
        line_width=0,
        row=1,
        col=2,
    )
    fig.add_annotation(
        text="third star of 6.6–7 M☉<br>(MIST, 6 Myr)",
        x=75.6,
        y=100 * e["one_6.6Msun"] - 1.0,
        showarrow=False,
        yanchor="top",
        font=dict(size=10, color=C_ART),
        row=1,
        col=2,
    )
    fig.add_hrect(y0=-25, y1=2.0, fillcolor="rgba(0,0,0,0.05)", line_width=0, row=1, col=2)
    fig.add_annotation(
        text="spectra: third light < 2 %<br>(Pavlovski et al. 2018)",
        x=75.6,
        y=-6,
        showarrow=False,
        font=dict(size=10, color="#444"),
        row=1,
        col=2,
    )
    fig.add_hline(y=0, line=dict(color="#888"), row=1, col=2)
    fig.add_vrect(x0=69.2, x1=70.4, fillcolor="rgba(0,0,0,0.08)", line_width=0, row=1, col=2)
    fig.add_annotation(
        text="published i<br>69.8 ± 0.6°",
        x=69.8,
        y=19,
        showarrow=False,
        font=dict(size=10),
        row=1,
        col=2,
    )
    fl = n.get("fit_lc_derived")
    if fl:
        fig.add_trace(
            go.Scatter(
                x=[fl["i"]],
                y=[100 * fl["l3_eff_mean"]],
                mode="markers",
                marker=dict(size=13, symbol="star", color=C_Q2, line=dict(width=1, color="#333")),
                error_x=dict(type="data", array=[fl["i_err"]]),
                error_y=dict(type="data", array=[100 * fl["l3_eff_spread"]]),
                name="SWINGERsim fit of the TESS light curves (i, radii free)",
            ),
            1,
            2,
        )
    fig.update_xaxes(title_text="inclination of the pair [°]", range=[69, 76.6], row=1, col=2)
    fig.update_yaxes(title_text="extra light [% of the total]", range=[-21, 21], row=1, col=2)
    _layout(fig, h=560, w=1300)
    fig.update_layout(legend=dict(y=1.06), margin=dict(t=120))
    return fig


def ah_fourth(path, fit, t0, n, res_path, rate_limit=0.06):
    """Possible fourth bodies on an orbit around the pair (hierarchy ((pair, m4), m3)), stable by Mardling & Aarseth:
    which of them would have been seen in the eclipse times (their light-time effect added to the residuals of the
    three-body fit makes chi^2 worse by more than 25, after refitting a quadratic ephemeris), and how fast would the
    hidden ones tilt the pair (to be compared with the limit from the constant TESS eclipse depths 2019-2024)."""
    d = pd.read_csv(path)
    f = d[d.direction == "forward"].set_index("member")
    b = d[d.direction == "backward"].set_index("member")
    m = f.join(b, rsuffix="_b", how="inner")
    m = m[(m.fate == "intact") & (m.fate_b == "intact")]
    rate = (0.5 * (m.inner_didt0_deg_yr + m.inner_didt0_deg_yr_b)).abs().values
    res = pd.read_csv(res_path)
    res = res[res.status == "fitted"]
    t = res.bjd_tdb.values
    r = res.o_minus_model_min.values
    W = 1.0 / (res.sigma_d.values * 1440.0)
    mpair = fit["masses"]["m1"] + fit["masses"]["m2"]
    X = np.column_stack([np.ones_like(t), (t - t0) / 3e4, ((t - t0) / 3e4) ** 2])

    def chi2(y):
        c = np.linalg.lstsq(X * W[:, None], y * W, rcond=None)[0]
        return float(np.sum(((y - X @ c) * W) ** 2))

    c0 = chi2(r)
    dchi, amp = [], []
    for _, row in m.iterrows():
        P = row["orbits.mid.P"]
        e = row["orbits.mid.e"]
        m4 = row["masses.m4"]
        a = (
            ((mpair + m4) * (P / DPY) ** 2) ** (1 / 3) * m4 / (mpair + m4)
        )  # au, orbit of the pair around the pair-m4 centre
        M = np.mod(np.radians(row["orbits.mid.M0"]) + 2 * np.pi * (t - t0) / P, 2 * np.pi)
        E = M.copy()
        for _ in range(40):
            E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
        z = (
            a
            * (1 - e * np.cos(E))
            * np.sin(nu + np.radians(row["orbits.mid.omega"]))
            * np.sin(np.radians(row["orbits.mid.i"]))
            * 499.005
            / 60.0
        )
        dchi.append(chi2(r - z) - c0)
        amp.append(0.5 * float(np.ptp(z)))
    dchi = np.array(dchi)
    amp = np.array(amp)
    P = m["orbits.mid.P"].values
    m4 = m["masses.m4"].values
    hidden = dchi < 25.0
    fast = rate > rate_limit

    def frac(sel):  # fraction of hidden members in a selection (NaN for an empty one)
        return float(hidden[sel].mean()) if sel.any() else float("nan")

    def top(x):  # largest value among the hidden members (NaN if none is hidden)
        return float(x[hidden].max()) if hidden.any() else float("nan")

    if not len(m):
        n["fourth_body"] = {"n_stable": 0}
        return go.Figure()
    n["fourth_body"] = {
        "n_stable": len(m),
        "P_range_d": [float(P.min()), float(P.max())],
        "m4_range": [float(m4.min()), float(m4.max())],
        "hidden": int(hidden.sum()),
        "hidden_fraction": float(hidden.mean()),
        "hidden_fraction_m4_gt3": frac(m4 > 3),
        "hidden_fraction_m4_lt1": frac(m4 < 1),
        "max_amp_hidden_min": top(amp),
        "rate_limit_deg_yr": rate_limit,
        "hidden_fast": int((hidden & fast).sum()),
        "hidden_slow": int((hidden & ~fast).sum()),
        "hidden_fast_fraction_of_hidden": float((hidden & fast).sum() / max(hidden.sum(), 1)),
        "max_rate_hidden_deg_yr": top(rate),
    }
    fig = go.Figure()
    for sel, nm, sym, col in (
        (~hidden, "would have been seen in the eclipse times", "x-thin", "#9AA0A6"),
        (
            hidden & fast,
            f"hidden in the eclipse times, but tilts the pair faster than {rate_limit}°/yr (ruled out by TESS)",
            "circle-open",
            C_ART,
        ),
        (
            hidden & ~fast,
            f"hidden, tilts the pair slower than {rate_limit}°/yr: allowed",
            "circle",
            C_Q,
        ),
    ):
        fig.add_trace(
            go.Scatter(
                x=P[sel] / DPY,
                y=m4[sel],
                mode="markers",
                marker=dict(size=6, symbol=sym, color=col, line=dict(width=1.2, color=col)),
                name=nm,
                text=[
                    f"|di/dt| {q:.3f}°/yr, light-time ±{aa:.2f} min"
                    for q, aa in zip(rate[sel], amp[sel], strict=True)
                ],
                hovertemplate="%{text}<extra></extra>",
            )
        )
    fig.update_xaxes(
        type="log",
        title_text="period of the fourth body around the pair [yr]",
        tickvals=[0.3, 0.5, 1, 2, 3],
        ticktext=["0.3", "0.5", "1", "2", "3"],
    )
    fig.update_yaxes(
        type="log",
        title_text="mass of the fourth body [M☉]",
        tickvals=[0.3, 0.5, 1, 2, 3, 5, 8],
        ticktext=["0.3", "0.5", "1", "2", "3", "5", "8"],
    )
    _layout(fig, h=500, w=1000)
    fig.update_layout(
        legend=dict(orientation="v", y=1.02, yanchor="bottom", x=0), margin=dict(t=120)
    )
    return fig


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument(
        "--no-png",
        action="store_true",
        help="JSON figures only (no kaleido): missing PNGs are not an error",
    )
    a = ap.parse_args(argv)
    out_root = (a.output or a.root / "output") / "squadra_test3"
    DATA_DIRS[:] = _data_dirs(a.root, out_root)[1:]
    O = Out(a.out or a.root / "docs" / "figures" / "squadra_test3", png=not a.no_png)
    nf = O.out / "numbers.json"
    N = json.loads(nf.read_text()) if nf.exists() else {}
    for name, fn in (("v907_sco", v907), ("qx_cas", qx_cas), ("ah_cep", ah_cep)):
        if a.only and name not in a.only:
            continue
        print(name, flush=True)
        N[name] = {}
        fn(a.root, out_root, O, N)
    atomic_write_text(nf, json.dumps(N, indent=1, default=float))
    O.finish()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
