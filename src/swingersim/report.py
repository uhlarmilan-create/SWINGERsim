"""A self-contained HTML report of one SWINGERsim command (one file, opens in any web browser).

Every command that produces results writes, next to its output, a report with the command
line, the configuration, the printed summary, interactive plots of the results (plotly, embedded
so that the file works offline and can be sent by e-mail) and the PNG figures the command made:

    swingersim run cfg.yaml -o out/run.h5            ->  out/run.html
    swingersim ensemble cfg.yaml -o out/ens          ->  out/ens/report.html
    swingersim lightcurve ... -o out/lc.csv          ->  out/lc.html
    swingersim validate                              ->  output/validation/report.html

``--html PATH`` chooses another file, ``--no-html`` switches the report off. Without plotly
(``pip install -e ".[viz]"``) the report has no interactive plots but everything else.
"""

from __future__ import annotations

import base64
import csv
import datetime as dt
import html
import io
import itertools
import math
import sys
from pathlib import Path

import numpy as np

# commands without a report (servers, the viewer itself, help-type commands)
NO_REPORT = {"gui", "worker", "view", "pulsation", "applegate", "mesa-grid", "mesa-reduce"}
MAX_LOG_LINES = 400
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#8e5bd0", "#c9a227", "#d4507a", "#4aa3a2", "#7a7a7a"]

CSS = """
:root { --bg: #fcfcfb; --fg: #1d1d1b; --muted: #5d5d58; --line: #dedcd6; --card: #ffffff;
        --pass: #1b8a5a; --open: #b86e00; --fail: #c0392b; --accent: #2a78d6; }
@media (prefers-color-scheme: dark) {
  :root { --bg: #17181a; --fg: #ececea; --muted: #a5a5a0; --line: #34353a; --card: #1f2023;
          --pass: #46c08a; --open: #e0a040; --fail: #ef6f5f; --accent: #6aa8f0; } }
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--fg);
       font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1180px; margin: 0 auto; padding: 24px 16px 64px; }
h1 { font-size: 24px; margin: 0 0 4px; } h2 { font-size: 18px; margin: 32px 0 8px; }
.meta { color: var(--muted); font-size: 13px; }
.card { background: var(--card); border: 1px solid var(--line); border-radius: 8px;
        padding: 12px 16px; margin: 12px 0; overflow-x: auto; }
pre { margin: 0; font: 12.5px/1.45 ui-monospace, "DejaVu Sans Mono", monospace;
      white-space: pre; }
table { border-collapse: collapse; font-size: 13.5px; }
th, td { border-bottom: 1px solid var(--line); padding: 4px 10px; text-align: left;
         vertical-align: top; }
th { color: var(--muted); font-weight: 600; }
td.num { text-align: right; font-variant-numeric: tabular-nums; }
.PASS { color: var(--pass); font-weight: 600; } .EXPLAINED { color: var(--open); font-weight: 600; }
.FAIL, .ERROR { color: var(--fail); font-weight: 600; } .SKIP { color: var(--muted); }
img { max-width: 100%; height: auto; border-radius: 4px; background: #fff; }
details summary { cursor: pointer; color: var(--accent); }
nav a { color: var(--accent); margin-right: 14px; font-size: 14px; text-decoration: none; }
"""


# --------------------------------------------------------------------------- page assembly


class Report:
    def __init__(self, title: str, base: Path | None = None):
        self.title = title
        # folder of the report: PNG files there are referenced, not embedded (see image)
        self.base = None if base is None else Path(base)
        self.blocks: list[tuple[str, str]] = []  # (section title, html)
        self.n_fig = 0

    def text(self, heading: str, text: str) -> None:
        if text and text.strip():
            body = f'<div class="card"><pre>{html.escape(text.rstrip())}</pre></div>'
            self.blocks.append((heading, body))

    def table(self, heading: str, header: list[str], rows: list[list], note: str = "") -> None:
        if not rows:
            return
        head = "".join(f"<th>{html.escape(str(h))}</th>" for h in header)
        body = []
        for r in rows:
            cells = []
            for v in r:
                if isinstance(v, tuple):  # (text, css class)
                    cells.append(f'<td class="{v[1]}">{html.escape(str(v[0]))}</td>')
                elif isinstance(v, (int, float, np.floating, np.integer)):
                    cells.append(f'<td class="num">{_fmt(v)}</td>')
                else:
                    cells.append(f"<td>{html.escape(str(v))}</td>")
            body.append("<tr>" + "".join(cells) + "</tr>")
        extra = f'<p class="meta">{html.escape(note)}</p>' if note else ""
        self.blocks.append(
            (
                heading,
                (f'<div class="card"><table><tr>{head}</tr>{"".join(body)}</table>{extra}</div>'),
            )
        )

    def figure(self, heading: str, fig) -> None:
        if fig is None:
            return
        self.n_fig += 1
        div = fig.to_html(
            full_html=False,
            include_plotlyjs=False,
            auto_play=False,
            config={"responsive": True, "displaylogo": False},
        )
        self.blocks.append((heading, f'<div class="card">{div}</div>'))

    def image(self, heading: str, path: Path) -> None:
        """A PNG figure: referenced by its file name when it lies next to the report (the
        products of a run), embedded (base64) otherwise."""
        path = Path(path)
        if self.base is not None and path.exists():
            try:
                same = path.parent.resolve() == self.base.resolve()
            except OSError:
                same = False
            if same:
                self.blocks.append(
                    (
                        heading,
                        (
                            f'<div class="card"><img alt="{html.escape(heading)}" '
                            f'src="{html.escape(path.name)}"></div>'
                        ),
                    )
                )
                return
        try:
            data = base64.b64encode(path.read_bytes()).decode("ascii")
        except OSError:
            return
        self.blocks.append(
            (
                heading,
                (
                    f'<div class="card"><img alt="{html.escape(heading)}" '
                    f'src="data:image/png;base64,{data}"></div>'
                ),
            )
        )

    def link(self, heading: str, href: str, text: str, note: str = "") -> None:
        extra = f'<p class="meta">{html.escape(note)}</p>' if note else ""
        self.blocks.append(
            (
                heading,
                (
                    f'<div class="card"><a href="{html.escape(href)}">{html.escape(text)}</a>'
                    f"{extra}</div>"
                ),
            )
        )

    def details(self, heading: str, summary: str, text: str) -> None:
        if text:
            self.blocks.append(
                (
                    heading,
                    (
                        f'<div class="card"><details><summary>{html.escape(summary)}'
                        f"</summary><pre>{html.escape(text)}</pre></details></div>"
                    ),
                )
            )

    def render(self, meta: str) -> str:
        plotly_js = ""
        if self.n_fig:
            from plotly.offline import get_plotlyjs

            plotly_js = f"<script>{get_plotlyjs()}</script>"
        sections, nav, last = [], [], None
        for k, (heading, body) in enumerate(self.blocks):
            if heading != last:
                sections.append(f'<h2 id="s{k}">{html.escape(heading)}</h2>')
                nav.append(f'<a href="#s{k}">{html.escape(heading)}</a>')
                last = heading
            sections.append(body)
        return (
            '<!doctype html>\n<html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width, initial-scale=1">'
            f"<title>{html.escape(self.title)}</title><style>{CSS}</style>{plotly_js}</head>"
            f'<body><main><h1>{html.escape(self.title)}</h1><p class="meta">{meta}</p>'
            f"<nav>{''.join(nav)}</nav>{''.join(sections)}</main></body></html>\n"
        )


def _fmt(v) -> str:
    v = float(v)
    if math.isnan(v):
        return "–"
    if v != 0 and (abs(v) >= 1e6 or abs(v) < 1e-3):
        return f"{v:.4g}"
    return f"{v:.6g}"


# --------------------------------------------------------------------------- plotly helpers


def _plotly():
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        return go, make_subplots
    except ImportError:
        return None, None


def _layout(fig, title: str = "", height: int = 420, **kw):
    fig.update_layout(
        title=title,
        height=height,
        template="plotly_white",
        margin={"l": 60, "r": 20, "t": 50 if title else 20, "b": 50},
        legend={"orientation": "h", "y": -0.18},
        colorway=PALETTE,
        **kw,
    )
    return fig


def _read_csv(path: Path) -> dict[str, np.ndarray | list]:
    with open(path, newline="", encoding="utf-8") as f:
        # leading '#' lines (provenance, units) are skipped; files without them work as before
        from .tables import data_lines

        rows = list(csv.DictReader(data_lines(f)))
    if not rows:
        return {}
    out = {}
    for k in rows[0]:
        vals = [r[k] for r in rows]
        try:
            out[k] = np.array([float(v) if v not in ("", None) else np.nan for v in vals])
        except ValueError:
            out[k] = vals
    return out


# --------------------------------------------------------------------------- per command

MAX_PLOT_POINTS = 5000  # points per plotted series of a run (min/max envelope beyond)

# the columns of the "Stars" table of a run report: attribute -> header
STAR_COLUMNS = {
    "mass_msun": "M [Msun]",
    "mass_msun_config": "M config [Msun]",
    "radius_rsun": "R [Rsun]",
    "teff_k": "Teff [K]",
    "luminosity_lsun": "L [Lsun]",
    "phase": "phase",
    "age_myr": "age [Myr]",
    "apsidal_k2": "k2 (apsidal)",
    "lag_s": "tidal lag [s]",
    "roche_fill_max": "max R/R_L",
}


def _cell(v):
    if isinstance(v, bytes):
        return v.decode()
    if isinstance(v, (np.floating, np.integer)):
        return float(v)
    return v


def decimate(x, y, n_max: int = MAX_PLOT_POINTS):
    """At most about n_max points of a series for plotting: the minimum and the maximum of y
    in each of n_max / 2 equal bins of the index (keeps the peaks a stride would miss)."""
    x, y = np.asarray(x), np.asarray(y, float)
    n = len(y)
    if n <= n_max:
        return x, y
    nb = max(n_max // 2, 1)
    edges = np.linspace(0, n, nb + 1).astype(int)
    idx = []
    for lo, hi in itertools.pairwise(edges):
        if hi <= lo:
            continue
        seg = y[lo:hi]
        ok = np.isfinite(seg)
        if not ok.any():
            idx.append(lo)  # keep a gap (NaN) visible
            continue
        i1 = lo + int(np.nanargmin(seg))
        i2 = lo + int(np.nanargmax(seg))
        idx.extend(sorted({i1, i2}))
    idx = np.unique(np.concatenate([[0], idx, [n - 1]]).astype(int))
    return x[idx], y[idx]


def _run_sections(rep: Report, h5: Path) -> None:
    import h5py

    go, make_subplots = _plotly()
    with h5py.File(h5) as f:
        t = f["time_yr"][:]
        orbits = {n: {k: f["orbits"][n][k][:] for k in f["orbits"][n]} for n in f["orbits"]}
        mutual = {k: f["mutual_inclination"][k][:] for k in f["mutual_inclination"]}
        energy = f["energy_error"][:] if "energy_error" in f else None
        event = str(f.attrs.get("event", "") or "")
        stars = {}
        if "stars" in f:
            for k in f["stars"]:
                stars[k] = {a: f["stars"][k].attrs[a] for a in f["stars"][k].attrs}
    rows = []
    for n, o in orbits.items():
        row = [n]
        for key in ("a", "e"):
            v = o.get(key)
            row += [float(np.nanmin(v)), float(np.nanmax(v))] if v is not None else ["", ""]
        rows.append(row)
    rep.table(
        "Orbits",
        ["orbit", "a min [AU]", "a max [AU]", "e min", "e max"],
        rows,
        note=f"{len(t)} outputs over {t[-1] - t[0]:.4g} yr"
        + (f"; event: {event}" if event else ""),
    )
    if stars:
        present = {k for st in stars.values() for k in st}
        keys = [k for k in STAR_COLUMNS if k in present]
        labels = [STAR_COLUMNS[k] for k in keys]
        rep.table(
            "Stars",
            ["star", *labels],
            [[s, *[_cell(stars[s].get(k, "")) for k in keys]] for s in stars],
        )
    viewer = h5.with_name(f"{h5.stem}_viewer.html")
    viewer_ok = viewer.exists()
    if viewer_ok:  # the animated orbit viewer of swingersim view: a separate page of the run
        rep.link(
            "Orbits in the sky (animation)",
            viewer.name,
            f"open the 3D viewer of the run ({viewer.name})",
            note="animated orbits, inclinations and conservation (written by swingersim run)",
        )
    if go is None:
        return
    tt, unit = (t / 1e6, "Myr") if abs(t[-1]) > 1e5 else (t, "yr")
    note = "" if len(t) <= MAX_PLOT_POINTS else f" (envelope of {len(t)} outputs)"
    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.06,
        subplot_titles=("eccentricity", "semi-major axis [AU]"),
    )
    for k, (n, o) in enumerate(orbits.items()):
        c = PALETTE[k % len(PALETTE)]
        x, y = decimate(tt, o["e"])
        fig.add_trace(
            go.Scattergl(x=x, y=y, name=f"{n}", line={"color": c, "width": 1.4}, legendgroup=n),
            1,
            1,
        )
        x, y = decimate(tt, o["a"])
        fig.add_trace(
            go.Scattergl(
                x=x,
                y=y,
                name=f"{n}",
                line={"color": c, "width": 1.4},
                legendgroup=n,
                showlegend=False,
            ),
            2,
            1,
        )
    fig.update_yaxes(type="log", row=2, col=1)
    fig.update_xaxes(title_text=f"time [{unit}]", row=2, col=1)
    rep.figure("Orbital elements" + note, _layout(fig, height=560))
    if viewer_ok:
        return
    if mutual:
        fig = go.Figure()
        for n, v in mutual.items():
            x, y = decimate(tt, np.degrees(v))
            fig.add_trace(go.Scattergl(x=x, y=y, name=n, line={"width": 1.4}))
        fig.update_xaxes(title_text=f"time [{unit}]")
        fig.update_yaxes(title_text="mutual inclination [deg]")
        rep.figure("Mutual inclinations", _layout(fig))
    if energy is not None and np.any(np.isfinite(energy)):
        x, y = decimate(tt, np.abs(energy) + 1e-17)
        fig = go.Figure(go.Scattergl(x=x, y=y, name="|dE/E|", line={"width": 1.2}))
        fig.update_yaxes(type="log", title_text="|dE/E| (or |dL/L| with dissipation)")
        fig.update_xaxes(title_text=f"time [{unit}]")
        rep.figure("Conservation check", _layout(fig, height=320))


def _ensemble_sections(rep: Report, out_dir: Path) -> None:
    from collections import Counter

    import yaml

    from .ensemble import load_rows

    rows = load_rows(out_dir)
    if not rows:
        return
    meta = (
        yaml.safe_load((out_dir / "ensemble.yaml").read_text())
        if (out_dir / "ensemble.yaml").exists()
        else {}
    )
    vary = list(((meta.get("config") or {}).get("ensemble") or {}).get("vary") or {})
    dirs = sorted({r.get("direction", "forward") for r in rows})
    fates = Counter((r.get("direction", "forward"), r.get("fate", "?")) for r in rows)
    names = sorted({f for _, f in fates})
    rep.table(
        "Fates",
        ["fate", *dirs],
        [[f, *[fates.get((d, f), 0) for d in dirs]] for f in names],
        note=f"{len(rows)} runs ({len({r['member'] for r in rows})} members)",
    )
    go, make_subplots = _plotly()
    if go is None:
        return
    fig = go.Figure()
    for k, d in enumerate(dirs):
        fig.add_trace(
            go.Bar(
                x=names, y=[fates.get((d, f), 0) for f in names], name=d, marker_color=PALETTE[k]
            )
        )
    fig.update_yaxes(title_text="runs")
    rep.figure("Fates", _layout(fig, height=340, barmode="group"))
    first = "forward" if "forward" in dirs else dirs[0]
    fwd = [r for r in rows if r.get("direction", "forward") == first]
    ok = [r.get("fate") == "intact" for r in fwd]
    params = [p for p in vary if any(isinstance(r.get(p), (int, float)) for r in fwd)][:12]
    if params:
        ncol = 3
        nrow = math.ceil(len(params) / ncol)
        fig = make_subplots(
            rows=nrow, cols=ncol, subplot_titles=params, vertical_spacing=0.12 / max(nrow / 3, 1)
        )
        for k, p in enumerate(params):
            x = np.array([float(r.get(p, np.nan)) for r in fwd])
            for good, c, lab in ((True, PALETTE[0], "intact"), (False, PALETTE[1], "not intact")):
                sel = np.array(ok) == good
                fig.add_trace(
                    go.Histogram(
                        x=x[sel],
                        name=lab,
                        marker_color=c,
                        opacity=0.75,
                        legendgroup=lab,
                        showlegend=k == 0,
                        nbinsx=30,
                    ),
                    k // ncol + 1,
                    k % ncol + 1,
                )
        rep.figure(
            f"Drawn parameters ({first}): intact and not intact",
            _layout(fig, height=260 * nrow, barmode="stack"),
        )
    emax = sorted({k for r in fwd[:50] for k in r if k.endswith("_e_max")})
    if emax:
        fig = go.Figure()
        for k in emax:
            fig.add_trace(
                go.Histogram(
                    x=[r.get(k) for r in fwd], name=k.removesuffix("_e_max"), opacity=0.6, nbinsx=40
                )
            )
        fig.update_xaxes(title_text="largest eccentricity in the run")
        fig.update_yaxes(title_text="members")
        rep.figure("Largest eccentricities", _layout(fig, height=360, barmode="overlay"))


def _lightcurve_sections(rep: Report, path: Path) -> None:
    d = _read_csv(path)
    if not d:
        return
    go, _ = _plotly()
    t, f = d["t_day"], d["flux"]
    mag = -2.5 * np.log10(np.clip(f, 1e-12, None) / np.nanmax(f))
    rep.table(
        "Light curve",
        ["points", "days", "min flux", "max flux", "largest dip [mag]"],
        [
            [
                len(t),
                float(t[-1] - t[0]),
                float(np.nanmin(f)),
                float(np.nanmax(f)),
                float(np.nanmax(mag)),
            ]
        ],
    )
    if go is None:
        return
    fig = go.Figure()
    fig.add_trace(go.Scattergl(x=t, y=f, name="total", mode="lines", line={"width": 1.2}))
    for k in d:
        if k.startswith("flux_"):
            fig.add_trace(
                go.Scattergl(
                    x=t, y=d[k], name=k[5:], mode="lines", line={"width": 1.0}, visible="legendonly"
                )
            )
    fig.update_xaxes(title_text="time [d after t = 0]", rangeslider={"visible": True})
    fig.update_yaxes(title_text="relative flux")
    rep.figure("Light curve (click the legend for the stars)", _layout(fig, height=520))


def _oc_sections(rep: Report, path: Path, kind: str) -> None:
    d = _read_csv(path)
    if not d:
        return
    go, _ = _plotly()
    x = d["t_min_d"] / 365.25 if kind == "minima" else d["t_obs_yr"]
    y = d["o_c_d"] * 86400.0
    groups = sorted({(b, k) for b, k in zip(d["binary"], d["kind"], strict=True)})
    rows = []
    for b, k in groups:
        sel = np.array([(bb, kk) == (b, k) for bb, kk in zip(d["binary"], d["kind"], strict=True)])
        rows.append(
            [
                b,
                k,
                int(sel.sum()),
                float(np.nanmin(y[sel])),
                float(np.nanmax(y[sel])),
                float(np.nanstd(y[sel])),
            ]
        )
    rep.table("Eclipses", ["binary", "kind", "n", "O−C min [s]", "O−C max [s]", "rms [s]"], rows)
    if go is None:
        return
    fig = go.Figure()
    for j, (b, k) in enumerate(groups):
        sel = np.array([(bb, kk) == (b, k) for bb, kk in zip(d["binary"], d["kind"], strict=True)])
        fig.add_trace(
            go.Scattergl(
                x=x[sel],
                y=y[sel],
                mode="markers",
                name=f"{b} {k}",
                marker={"size": 4, "color": PALETTE[j % len(PALETTE)]},
            )
        )
    fig.update_xaxes(title_text="time [yr]")
    fig.update_yaxes(title_text="O − C [s]")
    rep.figure("O − C", _layout(fig, height=480))
    if kind == "minima" and "depth_mag" in d:
        fig = go.Figure()
        for j, (b, k) in enumerate(groups):
            sel = np.array(
                [(bb, kk) == (b, k) for bb, kk in zip(d["binary"], d["kind"], strict=True)]
            )
            fig.add_trace(
                go.Scattergl(
                    x=x[sel],
                    y=d["depth_mag"][sel],
                    mode="markers",
                    name=f"{b} {k}",
                    marker={"size": 4, "color": PALETTE[j % len(PALETTE)]},
                )
            )
        fig.update_xaxes(title_text="time [yr]")
        fig.update_yaxes(title_text="depth [mag]", autorange="reversed")
        rep.figure("Depths of the minima", _layout(fig, height=380))


def _validate_sections(rep: Report, outcomes) -> None:
    from collections import Counter

    verdicts = ("PASS", "EXPLAINED", "FAIL", "ERROR", "SKIP")
    rows, count = [], Counter()
    for o in outcomes:
        b = o.benchmark
        err = getattr(o, "error", "")
        if o.skipped:
            verdict = "SKIP"
        elif err:
            verdict = "ERROR"
        elif o.passed:
            verdict = "PASS"
        else:
            verdict = "EXPLAINED" if b.explanation else "FAIL"
        count[(b.tier, verdict)] += 1
        if err:
            obtained, dev = "-", err
        else:
            obtained = f"{o.obtained:.6g} {b.unit}"
            dev = f"{100 * o.deviation:+.2f} %" if b.relative else f"{o.deviation:+.3g} {b.unit}"
        tol = f"{100 * b.tolerance:.3g} %" if b.relative else f"{b.tolerance:.3g} {b.unit}"
        rows.append(
            [
                b.id,
                b.tier,
                b.quantity,
                f"{b.expected:.6g} {b.unit}",
                obtained,
                dev,
                tol,
                (verdict, verdict),
                f"{getattr(o, 'seconds', float('nan')):.1f}",
                b.source,
            ]
        )
    tiers = sorted({t for t, _ in count})
    rep.table(
        "Summary",
        ["tier", *verdicts],
        [[t, *(count[(t, v)] for v in verdicts)] for t in tiers],
    )
    rep.table(
        "Benchmarks",
        [
            "id",
            "tier",
            "quantity",
            "expected",
            "swingersim",
            "difference",
            "tolerance",
            "verdict",
            "s",
            "source",
        ],
        rows,
    )
    opens = [
        [o.benchmark.id, getattr(o, "error", "") or o.skipped or o.benchmark.explanation]
        for o in outcomes
        if not o.passed and (getattr(o, "error", "") or o.benchmark.explanation or o.skipped)
    ]
    rep.table("Explained differences, errors and skipped benchmarks", ["id", "why"], opens)


# --------------------------------------------------------------------------- entry points


def default_path(args) -> Path | None:
    """Where the report of a command goes (None: no report)."""
    if args.cmd in NO_REPORT or getattr(args, "no_html", False):
        return None
    if getattr(args, "html", None):
        return Path(args.html)
    out = getattr(args, "output", None)
    if args.cmd in ("ensemble", "chaosmap", "scatter") and out:
        return Path(out) / "report.html"
    if args.cmd == "ensemble-report":
        return Path(args.directory) / "report.html"
    if args.cmd == "fit":
        return Path(out) / "report.html" if out else None
    if args.cmd == "validate":
        return Path("output") / "validation" / "report.html"
    if out:
        return Path(out).with_suffix(".html")
    plot = getattr(args, "plot", None)
    return Path(plot).with_suffix(".html") if plot else None


def write(args, argv: list[str], log: str, path: Path, context: dict | None = None) -> Path:
    """Build and write the report of one finished command."""
    from . import __version__

    context = context or {}
    cmd = args.cmd
    title = f"SWINGERsim {cmd}"
    cfg_path = getattr(args, "config", None)
    if isinstance(cfg_path, (str, Path)):
        title += f" – {Path(cfg_path).stem}"
    rep = Report(title, base=Path(path).parent)
    lines = log.splitlines()
    if len(lines) > MAX_LOG_LINES:
        lines = [
            *lines[:40],
            f"... ({len(lines) - MAX_LOG_LINES} lines left out) ...",
            *lines[-(MAX_LOG_LINES - 40) :],
        ]
    out = getattr(args, "output", None)
    try:
        if cmd == "run" and out and Path(out).exists():
            _run_sections(rep, Path(out))
        elif cmd in ("ensemble", "ensemble-report"):
            _ensemble_sections(rep, Path(out or args.directory))
        elif cmd == "lightcurve" and out:
            _lightcurve_sections(rep, Path(out))
        elif cmd in ("minima", "timings") and out:
            _oc_sections(rep, Path(out), cmd)
        elif cmd == "validate" and "outcomes" in context:
            _validate_sections(rep, context["outcomes"])
    except Exception as exc:  # noqa: BLE001 - a broken plot must not lose the report
        rep.text("Plots", f"the plots could not be made: {type(exc).__name__}: {exc}")
    rep.text("Printed summary", "\n".join(lines))
    pngs = [getattr(args, k, None) for k in ("plot", "seasons_plot")]
    for p in pngs:
        if p and Path(p).suffix.lower() == ".png" and Path(p).exists():
            rep.image("Figures", Path(p))
    if cmd == "run" and out:
        from .runproducts import product_paths

        for key in ("orbits", "inclinations", "conservation", "stars"):
            p = product_paths(Path(out))[key]
            if p.exists():
                rep.image("Figures", p)
    out_dir = context.get("out_dir") or out
    if out_dir and Path(out_dir).is_dir():
        for p in sorted(Path(out_dir).glob("*.png"))[:20]:
            rep.image("Figures", p)
    if isinstance(cfg_path, (str, Path)) and Path(cfg_path).exists():
        rep.details("Configuration", str(cfg_path), Path(cfg_path).read_text(encoding="utf-8"))
    when = dt.datetime.now(dt.UTC).strftime("%Y-%m-%d %H:%M UTC")
    meta = (
        f"SWINGERsim {html.escape(__version__)} · {when} · "
        f"<code>{html.escape(' '.join(['swingersim', *argv]))}</code>"
    )
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(rep.render(meta), encoding="utf-8")
    if cmd == "run" and out and Path(out).exists():
        from .runproducts import add_report_to_summary

        add_report_to_summary(Path(out), path)
    return path


class Tee(io.TextIOBase):
    """Writes to the real stdout and keeps a copy for the report."""

    def __init__(self, stream):
        self.stream = stream
        self.buffer_ = io.StringIO()

    def write(self, s):
        self.buffer_.write(s)
        try:
            return self.stream.write(s)
        except ValueError:  # the real stream was closed (e.g. by a test runner)
            return len(s)

    def flush(self):
        try:
            self.stream.flush()
        except ValueError:
            pass

    def getvalue(self) -> str:
        return self.buffer_.getvalue()

    @property
    def encoding(self):
        return getattr(self.stream, "encoding", "utf-8")

    def isatty(self):
        try:
            return self.stream.isatty()
        except ValueError:
            return False

    def fileno(self):
        return self.stream.fileno()


def capture():
    """Start copying stdout; returns the Tee (sys.stdout is restored by release())."""
    tee = Tee(sys.stdout)
    sys.stdout = tee
    return tee


def release(tee: Tee) -> None:
    if sys.stdout is tee:
        sys.stdout = tee.stream


__all__ = ["Report", "capture", "default_path", "release", "write"]
