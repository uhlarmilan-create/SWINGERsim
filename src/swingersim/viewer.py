"""Interactive HTML viewer for a swingersim run (HDF5 -> one self-contained HTML file).

    swingersim view output/run.h5 -o output/run.html

The viewer reconstructs the star positions and the orbital ellipses at each output from the
stored osculating elements (exact, because the elements were computed from those positions),
and shows them together with the time series, with a time slider and a play button.

Requires plotly:  pip install -e ".[viz]"
"""

from __future__ import annotations

import json
from itertools import pairwise
from pathlib import Path

import h5py
import numpy as np
import yaml

from . import units
from .backends.nbody import default_radius_rsun
from .config import Config
from .hierarchy import Node, mass_of, skeleton
from .kepler import Elements, elements_to_state

COLORS = ["#e8a33d", "#d9603b", "#2a6fdb", "#3b9c5a", "#8e5bd0", "#a07a2a", "#cc4f95"]  # stars
ORBIT_COLORS = {
    "out": "#7f7f7f",
    "A": "#d9603b",
    "B": "#2a6fdb",
    "mid": "#2a6fdb",
    "inner": "#d9603b",
    "out5": "#b0b0b0",
    "Bt": "#5aa0e0",
}
N_ELLIPSE = 64
STATIC_MAX_POINTS = 2000  # points of each static time series (min/max decimation)
SIG_DIGITS = 6  # significant digits of the plotted data (HTML size)


def _compact(a) -> np.ndarray:
    """Data rounded to SIG_DIGITS significant digits for the HTML: float32 where plotly
    stores arrays in binary (plotly >= 6), short decimal numbers otherwise."""
    import plotly

    a = np.array(a, dtype=float, copy=True)
    ok = np.isfinite(a) & (a != 0.0)
    if ok.any():
        e = np.floor(np.log10(np.abs(a[ok])))
        scale = 10.0 ** (SIG_DIGITS - 1 - e)
        a[ok] = np.round(a[ok] * scale) / scale
    if int(plotly.__version__.split(".")[0]) >= 6:
        return a.astype(np.float32)
    return a


def _decimate(n: int, *ys, n_max: int = STATIC_MAX_POINTS) -> np.ndarray:
    """Indices of at most ~n_max points of series of length n that keep the envelope: the
    first and last point and the minimum and maximum of every bin of each series ys."""
    if n <= n_max:
        return np.arange(n)
    n_bins = max(n_max // (2 * max(len(ys), 1)), 1)
    edges = np.linspace(0, n, n_bins + 1).astype(int)
    keep = {0, n - 1}
    for y in ys:
        y = np.asarray(y, float)
        for lo, hi in pairwise(edges):
            seg = y[lo:hi]
            if hi <= lo or not np.isfinite(seg).any():
                continue
            keep.add(lo + int(np.nanargmin(seg)))
            keep.add(lo + int(np.nanargmax(seg)))
    return np.array(sorted(keep))


def _load(path: Path):
    with h5py.File(path) as f:
        cfg = Config.model_validate(yaml.safe_load(f.attrs["config_yaml"]))
        times = f["time_yr"][:]
        orbits = {n: {k: f["orbits"][n][k][:] for k in f["orbits"][n]} for n in f["orbits"]}
        mutual = {k: f["mutual_inclination"][k][:] for k in f["mutual_inclination"]}
        transitions = (
            json.loads(f["hierarchy"].attrs["transitions_json"]) if "hierarchy" in f else []
        )
        event = str(f.attrs.get("event", ""))
        energy_error = f["energy_error"][:]
        radii, mass0, mass_series = {}, {}, {}
        for k in f.get("stars", {}):
            g = f["stars"][k]
            if "radius_rsun" in g.attrs:
                radii[k] = float(g.attrs["radius_rsun"])
            if "mass_msun" in g.attrs:  # dynamical mass at t = 0
                mass0[k] = float(g.attrs["mass_msun"])
            if "mass_msun" in g and len(g["mass_msun"]) == len(times):
                mass_series[k] = g["mass_msun"][:]
        try:
            system = json.loads(f.attrs.get("system_json", "{}") or "{}")
        except (TypeError, ValueError):
            system = {}
    return (
        cfg, times, orbits, mutual, transitions, event, energy_error, radii,
        {"mass0": mass0, "mass_series": mass_series, "system": system},
    )  # fmt: skip


def _tree_from(mapping: dict) -> Node | None:
    """Tree of a mapping orbit -> [child, child] with stars named m<k> (None when invalid)."""
    kids = {c for v in mapping.values() for c in v}
    roots = [o for o in mapping if o not in kids]
    if len(roots) != 1:
        return None

    def make(name, depth=0):
        if name in mapping and depth < 16:
            a, b = mapping[name]
            return Node(name, make(a, depth + 1), make(b, depth + 1))
        return int(str(name)[1:]) - 1

    try:
        return make(roots[0])
    except (ValueError, TypeError):
        return None


class _Epochs:
    """The structure of the system over the run: the hierarchy (it changes when stars merge,
    physics.mergers), the stars still present, their labels and masses."""

    def __init__(self, cfg, times, extra):
        names = list(cfg.star_names)
        n = len(names)
        m0 = np.array([extra["mass0"].get(k, cfg.masses[k]) for k in names], float)
        self.times = times
        self.base = np.tile(m0, (len(times), 1))
        for j, k in enumerate(names):
            ser = extra["mass_series"].get(k)
            if ser is not None:
                ok = np.isfinite(ser)
                self.base[ok, j] = ser[ok]
        root = skeleton(cfg.tree)
        self.epochs = [(-np.inf if times[-1] >= times[0] else np.inf, root, set(), {})]
        hidden, labels = set(), {}
        mergers = [
            m
            for m in (extra["system"].get("_mergers") or [])
            if isinstance(m, dict) and "stop" not in m and len(m.get("stars", [])) == 2
        ]
        for mg in sorted(mergers, key=lambda m: m["time_yr"]):
            i, j = (int(str(x)[1:]) - 1 for x in mg["stars"][:2])
            if not (0 <= i < n and 0 <= j < n):
                continue
            hidden = hidden | {j}
            group = sorted({i, j} | self._group(labels, i, j))
            labels = {**labels, i: "+".join(str(x + 1) for x in group)}
            tree = _tree_from(mg.get("orbits_tree") or {}) or self.epochs[-1][1]
            t = float(mg["time_yr"])
            self.epochs.append((t, tree, hidden, labels))
            after = (times >= t) if times[-1] >= times[0] else (times <= t)
            if mg.get("mass_msun") is not None and names[i] not in extra["mass_series"]:
                self.base[after, i] = float(mg["mass_msun"])
            self.base[after, j] = 0.0
        self.mergers = mergers

    @staticmethod
    def _group(labels, i, j) -> set:
        out = set()
        for k in (i, j):
            if k in labels:
                out |= {int(x) - 1 for x in labels[k].split("+")}
        return out

    def at(self, k: int):
        """(tree, hidden stars, labels, masses) at output k."""
        t = self.times[k]
        fwd = self.times[-1] >= self.times[0]
        cur = self.epochs[0]
        for ep in self.epochs[1:]:
            if (t >= ep[0]) if fwd else (t <= ep[0]):
                cur = ep
        return cur[1], cur[2], cur[3], self.base[k]

    def orbit_names(self) -> list[str]:
        out = []
        for _, tree, _, _ in self.epochs:
            for node in tree.nodes():
                if node.name not in out:
                    out.append(node.name)
        return out


def _elements(orbits, name, k) -> Elements | None:
    o = orbits[name]
    vals = [o[key][k] for key in ("a", "e", "inc", "Omega", "omega", "M")]
    if not np.all(np.isfinite(vals)):
        return None
    return Elements(*map(float, vals))


def _place(root: Node, masses, orbits, k):
    """Star positions and, for every orbit, the barycentric paths of its two components."""
    pos = np.full((len(masses), 3), np.nan)
    paths: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    def rec(child, R):
        if not isinstance(child, Node):
            pos[child] = R
            return
        m1, m2 = mass_of(child.left, masses), mass_of(child.right, masses)
        M = m1 + m2
        el = _elements(orbits, child.name, k) if child.name in orbits and M > 0 else None
        if el is None:  # orbit not followed here (e.g. its stars merged): both at its centre
            rec(child.left, R)
            rec(child.right, R)
            return
        r, _ = elements_to_state(units.G * M, el)
        # sample uniformly in eccentric anomaly: smooth also near periastron of eccentric orbits
        ecc_anom = np.linspace(0, 2 * np.pi, N_ELLIPSE)
        mean_anom = ecc_anom - el.e * np.sin(ecc_anom)
        mu = units.G * M
        ring = np.array(
            [
                elements_to_state(mu, Elements(el.a, el.e, el.inc, el.Omega, el.omega, m))[0]
                for m in mean_anom
            ]
        )
        paths[child.name] = (
            (R - (m2 / M) * ring).astype(np.float32),
            (R + (m1 / M) * ring).astype(np.float32),
            R,
        )
        rec(child.left, R - (m2 / M) * r)
        rec(child.right, R + (m1 / M) * r)

    rec(root, np.zeros(3))
    return pos, paths


def _frames_index(n: int, max_frames: int) -> np.ndarray:
    return np.unique(np.linspace(0, n - 1, min(n, max_frames)).astype(int))


def _zoom_title(nodes, k: int, mergers=()) -> str:
    """Title of zoom panel k (a triple has one inner orbit, so its second panel is empty)."""
    if k >= len(nodes):
        return ""
    title = f"Orbit {nodes[k].name}: sky view around its centre of mass [AU]"
    members = {f"m{s + 1}" for s in nodes[k].stars()}
    for mg in mergers:
        if set(mg.get("stars", [])[:2]) <= members:
            title += f"<br><sub>its stars merged at t = {mg['time_yr']:.4g} yr</sub>"
            break
    return title


def build_figure(path: str | Path, max_frames: int = 150):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    cfg, times, orbits, mutual, transitions, event, energy_error, radii, extra = _load(Path(path))
    epochs = _Epochs(cfg, times, extra)
    masses = epochs.base[0]
    root = skeleton(cfg.tree)
    all_orbits = epochs.orbit_names()
    inner_nodes = [n for n in root.nodes() if n is not root]  # zoom panels
    idx = _frames_index(len(times), max_frames)

    fig = make_subplots(
        rows=4,
        cols=2,
        specs=[
            [{"type": "scene"}, {"type": "xy"}],
            [{"type": "xy"}, {"type": "xy"}],
            [{"type": "xy"}, {"type": "xy"}],
            [{"type": "xy"}, {"type": "xy"}],
        ],
        subplot_titles=(
            "System in 3D [AU]",
            "System on the sky (x, y) [AU]",
            _zoom_title(inner_nodes, 0, epochs.mergers),
            _zoom_title(inner_nodes, 1, epochs.mergers),
            "Eccentricity",
            "Mutual inclination [deg]",
            "Inclination to the sky [deg]; shaded: eclipses possible",
            "Relative energy error",
        ),
        row_heights=[0.36, 0.28, 0.18, 0.18],
        vertical_spacing=0.07,
        horizontal_spacing=0.08,
    )

    # ---- static time series (rows 3) ----
    def series(y):
        """(x, y) of a static time series: decimated and rounded (HTML size)."""
        y = np.asarray(y, float)
        k = _decimate(len(times), y)
        return _compact(times[k]), _compact(y[k])

    for name, o in orbits.items():
        xs, ys = series(o["e"])
        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                name=f"e {name}",
                mode="lines",
                line={"color": ORBIT_COLORS.get(name, "#555"), "width": 1.5},
                legendgroup="ts",
            ),
            row=3,
            col=1,
        )
    for key, s in mutual.items():
        xs, ys = series(np.degrees(s))
        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                name=f"i_mut {key}",
                mode="lines",
                line={"width": 1.5},
                legendgroup="ts",
            ),
            row=3,
            col=2,
        )
    # sky inclination of the innermost binaries with eclipse windows
    for node in root.nodes():
        if isinstance(node.left, Node) or isinstance(node.right, Node):
            continue
        r_sum = 0.0
        for s in node.stars():
            spec = cfg.stars.get(f"m{s + 1}")
            if f"m{s + 1}" in radii:  # radius used in the run (stellar model)
                r_sum += radii[f"m{s + 1}"]
            else:  # files written before v0.4
                r_sum += spec.R if spec and spec.R else default_radius_rsun(masses[s])
        a_mean = float(np.nanmean(orbits[node.name]["a"]))
        hw = 90.0 - np.degrees(np.arccos(min(1.0, r_sum * units.RSUN_AU / a_mean)))
        color = ORBIT_COLORS.get(node.name, "#555")
        fig.add_trace(
            go.Scatter(
                x=[times[0], times[-1], times[-1], times[0]],
                y=[90 - hw, 90 - hw, 90 + hw, 90 + hw],
                fill="toself",
                mode="none",
                fillcolor=color,
                opacity=0.12,
                showlegend=False,
                hoverinfo="skip",
            ),
            row=4,
            col=1,
        )
        xs, ys = series(np.degrees(orbits[node.name]["inc"]))
        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines",
                name=f"i {node.name} (eclipses if |i-90| < {hw:.1f})",
                line={"color": color, "width": 1.5},
                legendgroup="ts",
            ),
            row=4,
            col=1,
        )
    xs, ys = series(np.maximum(energy_error, 1e-17))
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=ys,
            mode="lines",
            name="|dE/E|",
            line={"color": "#555", "width": 1},
            showlegend=False,
        ),
        row=4,
        col=2,
    )
    fig.update_yaxes(type="log", exponentformat="power", row=4, col=2)
    n_static = len(fig.data)

    # ---- animated traces: one builder used for the first frame and for every frame ----
    def animated(k):
        tree, hidden, labels, m_k = epochs.at(k)
        pos, paths = _place(tree, m_k, orbits, k)
        for h in hidden:  # a star absorbed by a merger
            pos[h] = np.nan
        text = [labels.get(i, str(i + 1)) for i in range(len(masses))]
        # every orbit in every frame (the number of traces must not change between frames):
        # an orbit not followed at this output (its stars merged) is drawn empty
        nan = np.full((1, 3), np.nan, dtype=np.float32)
        paths = {n: paths.get(n, (nan, nan, np.full(3, np.nan))) for n in all_orbits}
        out = []
        # 3D: orbits + stars
        for name, (p1, p2, _) in paths.items():
            for p in (p1, p2):
                out.append(
                    go.Scatter3d(
                        x=_compact(p[:, 0]),
                        y=_compact(p[:, 1]),
                        z=_compact(p[:, 2]),
                        mode="lines",
                        line={"color": ORBIT_COLORS.get(name, "#555"), "width": 3},
                        showlegend=False,
                        hoverinfo="skip",
                    )
                )
        out.append(
            go.Scatter3d(
                x=_compact(pos[:, 0]),
                y=_compact(pos[:, 1]),
                z=_compact(pos[:, 2]),
                mode="markers+text",
                text=text,
                marker={"size": 5, "color": COLORS[: len(masses)]},
                showlegend=False,
            )
        )
        # sky view
        for name, (p1, p2, _) in paths.items():
            for p in (p1, p2):
                out.append(
                    go.Scatter(
                        x=_compact(p[:, 0]),
                        y=_compact(p[:, 1]),
                        mode="lines",
                        line={"color": ORBIT_COLORS.get(name, "#555"), "width": 1.5},
                        showlegend=False,
                        hoverinfo="skip",
                        xaxis="x",
                        yaxis="y",
                    )
                )
        out.append(
            go.Scatter(
                x=_compact(pos[:, 0]),
                y=_compact(pos[:, 1]),
                mode="markers+text",
                text=text,
                textposition="top center",
                marker={"size": 9, "color": COLORS[: len(masses)]},
                showlegend=False,
                xaxis="x",
                yaxis="y",
            )
        )
        # zoom panels (relative to the orbit's centre of mass)
        for panel, node in zip(((2, 1), (2, 2)), inner_nodes):
            p1, p2, com = paths.get(
                node.name, (np.full((1, 3), np.nan), np.full((1, 3), np.nan), np.full(3, np.nan))
            )
            members = node.stars()
            for p in (p1, p2):
                out.append(
                    go.Scatter(
                        x=_compact(p[:, 0] - com[0]),
                        y=_compact(p[:, 1] - com[1]),
                        mode="lines",
                        line={"color": ORBIT_COLORS.get(node.name, "#555"), "width": 1.5},
                        showlegend=False,
                        hoverinfo="skip",
                    )
                )
            out.append(
                go.Scatter(
                    x=_compact(pos[members, 0] - com[0]),
                    y=_compact(pos[members, 1] - com[1]),
                    mode="markers+text",
                    text=[text[i] for i in members],
                    textposition="top center",
                    marker={"size": 9, "color": [COLORS[i] for i in members]},
                    showlegend=False,
                )
            )
        return out

    first = animated(idx[0])
    # subplot placement of the animated traces (same order as in `animated`)
    n_paths = 2 * len(all_orbits)
    placement = [(1, 1)] * (n_paths + 1) + [(1, 2)] * (n_paths + 1)
    for panel, node in zip(((2, 1), (2, 2)), inner_nodes):
        placement += [panel] * 3
    for tr, (r, c) in zip(first, placement):
        fig.add_trace(tr, row=r, col=c)
    anim_ids = list(range(n_static, len(fig.data)))

    def vline(t):
        return [
            {
                "type": "line",
                "xref": f"x{ax}",
                "yref": f"y{ax} domain",
                "x0": t,
                "x1": t,
                "y0": 0,
                "y1": 1,
                "line": {"color": "black", "width": 1, "dash": "dot"},
            }
            for ax in (4, 5, 6, 7)
        ]

    def extent(traces, dims=("x", "y")):
        vals = [np.abs(np.asarray(getattr(t, d), float)) for t in traces for d in dims]
        vals = np.concatenate([v[np.isfinite(v)] for v in vals if v.size]) if vals else []
        return float(np.max(vals)) * 1.08 if len(vals) else 1.0

    n_sys = n_paths + 1
    lim_3d = lim_sky = 0.0
    lim_zoom = [0.0, 0.0]
    zoom_k: list[list[float]] = [[], []]  # extent of each zoom panel in every frame
    frames = []
    for k in idx:
        data = animated(k)
        lim_3d = max(lim_3d, extent(data[:n_sys], ("x", "y", "z")))
        lim_sky = max(lim_sky, extent(data[n_sys : 2 * n_sys]))
        for z in range(2):
            ez = extent(data[2 * n_sys + 3 * z : 2 * n_sys + 3 * z + 3])
            lim_zoom[z] = max(lim_zoom[z], ez)
            zoom_k[z].append(ez)
        # assign the same subplot axes as the first frame
        for tr, ref in zip(data, [fig.data[i] for i in anim_ids]):
            if hasattr(ref, "xaxis") and ref.xaxis:
                tr.xaxis, tr.yaxis = ref.xaxis, ref.yaxis
            if hasattr(ref, "scene") and ref.scene:
                tr.scene = ref.scene
        frames.append(
            go.Frame(
                data=data,
                traces=anim_ids,
                name=f"{times[k]:.4g}",
                layout={"shapes": vline(times[k])},
            )
        )
    # a zoom panel whose orbit changes its size by more than a factor 2 during the run (mass
    # transfer, a binary circularised by tides, a merger): its range follows the orbit, with
    # the largest extent of the frames of a similar size
    for z, ax in enumerate(("2", "3")):
        e = np.array(zoom_k[z], float)
        if not len(e) or e.min() <= 0 or e.max() / e.min() <= 2.0:
            continue
        for f, ek in zip(frames, e):
            lim = float(e[(e > 0.5 * ek) & (e < 2.0 * ek)].max())
            f.layout.update(
                {f"xaxis{ax}": {"range": [-lim, lim]}, f"yaxis{ax}": {"range": [-lim, lim]}}
            )
    fig.frames = frames

    # fixed, symmetric ranges so the view does not jump during the animation
    rng3 = [-lim_3d, lim_3d]
    fig.update_scenes(xaxis_range=rng3, yaxis_range=rng3, zaxis_range=rng3, aspectmode="cube")
    fig.update_xaxes(range=[-lim_sky, lim_sky], row=1, col=2)
    fig.update_yaxes(range=[-lim_sky, lim_sky], row=1, col=2)
    for z, col in enumerate((1, 2)):
        fig.update_xaxes(range=[-lim_zoom[z], lim_zoom[z]], row=2, col=col)
        fig.update_yaxes(range=[-lim_zoom[z], lim_zoom[z]], row=2, col=col)
    # equal aspect for spatial panels
    fig.update_yaxes(scaleanchor="x", scaleratio=1, row=1, col=2)
    fig.update_yaxes(scaleanchor="x2", scaleratio=1, row=2, col=1)
    fig.update_yaxes(scaleanchor="x3", scaleratio=1, row=2, col=2)
    fig.update_xaxes(title_text="t [yr]", row=4, col=1)
    fig.update_xaxes(title_text="t [yr]", row=4, col=2)
    fig.update_scenes(xaxis_title="x", yaxis_title="y", zaxis_title="z (line of sight)")

    notes = []
    for t in transitions:
        if t.get("kind") == "merger":
            notes.append(f"{t['time_yr']:.4g} yr: merger {t.get('label', '')}, hierarchy {t['to']}")
        elif t.get("kind") == "mass_transfer_onset":
            notes.append(f"{t['time_yr']:.4g} yr: {t.get('label', 'mass transfer started')}")
        else:
            notes.append(
                f"{t['time_yr']:.4g} yr: {t['kind'].replace('_', ' ')} {t['from']} → {t['to']}"
            )
    hidden = [n.name for n in inner_nodes[2:]]  # two zoom panels: the first two inner orbits
    if hidden:
        notes.append(
            f"zoom panels: orbits {inner_nodes[0].name}, {inner_nodes[1].name}; without a zoom "
            f"panel: {', '.join(hidden)} (shown in the 3D and sky views)"
        )
    extra_orbits = [n for n in all_orbits if n not in {x.name for x in root.nodes()}]
    if extra_orbits:
        notes.append(f"orbits after the mergers: {', '.join(extra_orbits)} (3D and sky views)")
    p = Path(path)
    name = p.parent.name if p.stem in ("run", "output", "out") and p.parent.name else p.stem
    notes.insert(
        0,
        "drag the 3D view to turn it, scroll to zoom; ▶ play or the slider move in time; "
        "click a legend entry to hide a curve; the camera icon saves a picture",
    )
    title = (
        f"<b>{name}</b> — swingersim, {cfg.topology}, {len(cfg.star_names)} stars — "
        f"{len(times)} outputs to {times[-1]:.4g} yr"
        + (f" — event: {event.replace('_', ' ')}" if event else "")
    )
    sub = "<br>".join(notes[:6]) + ("<br>…" if len(notes) > 6 else "")
    n_lines = min(len(notes), 7)
    margin_t, height = 170 + 16 * n_lines, 1700 + 16 * n_lines
    plot_h = height - margin_t - 80
    fig.update_layout(
        title={
            "text": title,
            "x": 0.01,
            "xanchor": "left",
            "y": 1.0,
            "yanchor": "top",
            "yref": "container",
            "pad": {"t": 12},
        },
        annotations=[
            *fig.layout.annotations,
            *(
                [
                    {
                        "text": sub,
                        "xref": "paper",
                        "yref": "paper",
                        "x": 0.0,
                        "y": 1.0 + (margin_t - 40) / plot_h,
                        "xanchor": "left",
                        "yanchor": "top",
                        "align": "left",
                        "showarrow": False,
                        "font": {"size": 11, "color": "#555"},
                    }
                ]
                if sub
                else []
            ),
        ],
        height=height,
        margin={"t": margin_t, "b": 80},
        template="plotly_white",
        shapes=vline(times[idx[0]]),
        legend={"orientation": "h", "y": -0.04},
        updatemenus=[
            {
                "type": "buttons",
                "x": 0.0,
                "y": 1.055,
                "yanchor": "top",
                "xanchor": "left",
                "direction": "left",
                "buttons": [
                    {
                        "label": "▶ play",
                        "method": "animate",
                        "args": [
                            None,
                            {"frame": {"duration": 60, "redraw": True}, "fromcurrent": True},
                        ],
                    },
                    {
                        "label": "❚❚ pause",
                        "method": "animate",
                        "args": [
                            [None],
                            {"frame": {"duration": 0, "redraw": False}, "mode": "immediate"},
                        ],
                    },
                ],
            }
        ],
        sliders=[
            {
                "active": 0,
                "y": 1.075,
                "yanchor": "top",
                "x": 0.16,
                "len": 0.84,
                "font": {"size": 1, "color": "rgba(0,0,0,0)"},
                "ticklen": 0,
                "minorticklen": 0,  # no row of tick marks across the panel titles below
                "tickcolor": "rgba(0,0,0,0)",
                "currentvalue": {
                    "prefix": "t = ",
                    "suffix": " yr",
                    "visible": True,
                    "font": {"size": 14, "color": "#222"},
                },
                "steps": [
                    {
                        "label": f.name,
                        "method": "animate",
                        "args": [
                            [f.name],
                            {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"},
                        ],
                    }
                    for f in frames
                ],
            }
        ],
    )
    return fig


def is_run_file(path: Path) -> bool:
    """Is this HDF5 file a saved run (swingersim run -o), with orbits over time?"""
    try:
        with h5py.File(path, "r") as f:
            return "config_yaml" in f.attrs and "time_yr" in f and "orbits" in f
    except (OSError, ValueError):
        return False


def run_files(path: Path):
    """The run files of a path: the file itself, or every run file below a folder (not the GUI
    cache), sorted."""
    path = Path(path)
    if path.is_file():
        if is_run_file(path):
            yield path
        return
    for p in sorted(path.rglob("*.h5")):
        if ".gui_cache" in p.parts or p.name.endswith(".partial.h5"):
            continue
        if is_run_file(p):
            yield p


def write_html(
    h5_path: str | Path, html_path: str | Path, max_frames: int = 150, offline: bool = True
) -> Path:
    fig = build_figure(h5_path, max_frames=max_frames)
    html_path = Path(html_path)
    html_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(html_path, include_plotlyjs=True if offline else "cdn", auto_play=False)
    # a title for the browser tab
    p = Path(h5_path)
    name = p.parent.name if p.stem in ("run", "output", "out") and p.parent.name else p.stem
    text = html_path.read_text(encoding="utf-8")
    safe = name.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    text = text.replace("<head>", f"<head>\n    <title>{safe} – swingersim viewer</title>", 1)
    html_path.write_text(text, encoding="utf-8")
    return html_path


__all__ = ["build_figure", "default_radius_rsun", "is_run_file", "run_files", "write_html"]
