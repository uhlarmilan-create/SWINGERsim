"""The standard products of one run, written next to its HDF5 file.

``swingersim run cfg.yaml -o out/sys.h5`` (and every run job of the worker and the web
application) leaves, besides the HDF5 file and the HTML report ``out/sys.html``:

* ``sys_summary.txt`` - the printed summary of the run (orbits, precession, stability, events)
  and the list of the files below;
* ``sys_timeseries.csv`` - one row per output time: a, e, sky and invariable-plane inclination,
  node and periastron of every orbit, the mutual inclinations, the energy error and, for
  evolving stars, mass, radius and effective temperature of every star; for spreadsheets,
  TOPCAT or own scripts. The first lines, starting with '#', describe the run, the units and
  the columns (pandas: ``read_csv(path, comment="#")``; the csv module: skip the lines that
  start with '#'; TOPCAT reads the HDF5 file or the CSV without them). Runs with more than CSV_MAX_ROWS outputs get every k-th
  output (and the last one); the HDF5 file has all of them;
* ``sys_viewer.html`` - the interactive 3D viewer of the run (``swingersim view``; needs plotly);
* ``sys_orbits.png`` - semi-major axis, eccentricity and sky inclination of every orbit over time;
* ``sys_inclinations.png`` - the mutual inclinations of the orbits over time;
* ``sys_conservation.png`` - relative errors of energy and angular momentum (quality of the run);
* ``sys_stars.png`` - mass, radius and temperature of the stars (only for evolving stars).

``--no-viewer`` and ``--no-plots`` switch off the viewer and the PNG figures. The HTML report
(``sys.html``, swingersim.report) is written after these products and adds itself to the list of
files in the summary.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np

from .report import PALETTE

CSV_MAX_ROWS = 20000  # larger runs: every k-th output in the CSV (all of them in the HDF5)


def product_paths(h5: Path) -> dict[str, Path]:
    h5 = Path(h5)
    s = h5.with_suffix("")
    return {
        "summary": Path(f"{s}_summary.txt"),
        "timeseries": Path(f"{s}_timeseries.csv"),
        "viewer": Path(f"{s}_viewer.html"),
        "orbits": Path(f"{s}_orbits.png"),
        "inclinations": Path(f"{s}_inclinations.png"),
        "conservation": Path(f"{s}_conservation.png"),
        "stars": Path(f"{s}_stars.png"),
    }


def _read(h5: Path) -> dict:
    d: dict = {"orbits": {}, "mutual": {}, "stars": {}, "meta": {}}
    with h5py.File(h5) as f:
        for k in ("swingersim_version", "rebound_version", "reboundx_version", "config_yaml"):
            if k in f.attrs:
                v = f.attrs[k]
                d["meta"][k] = v.decode() if isinstance(v, bytes) else str(v)
        if "swingersim_version" not in d["meta"] and "quadsim_version" in f.attrs:
            v = f.attrs["quadsim_version"]  # files written before the rename
            d["meta"]["swingersim_version"] = v.decode() if isinstance(v, bytes) else str(v)
        d["meta"]["event"] = str(f.attrs.get("event", "") or "")
        d["t"] = f["time_yr"][:]
        for name in f.get("orbits", {}):
            g = f[f"orbits/{name}"]
            d["orbits"][name] = {k: g[k][:] for k in g}
        for name in f.get("mutual_inclination", {}):
            d["mutual"][name] = np.degrees(f[f"mutual_inclination/{name}"][:])
        for key in ("energy_error", "angmom_error"):
            d[key] = f[key][:] if key in f else None
        for name in f.get("stars", {}):
            g = f[f"stars/{name}"]
            d["stars"][name] = {k: g[k][:] for k in g if g[k].dtype.kind in "fiu"}
    return d


COLUMN_HELP = {
    "t_yr": "time from the start of the run [yr] (the last row of a run stopped by an event: "
    "the time of the event)",
    "a": "semi-major axis of the orbit [AU] (osculating, Jacobi)",
    "e": "eccentricity",
    "i": "inclination to the sky plane [deg]",
    "i_inv": "inclination to the invariable plane at t = 0 [deg]",
    "Omega": "longitude of the ascending node on the sky [deg]",
    "Omega_inv": "node on the invariable plane at t = 0 [deg]",
    "omega": "argument of periastron [deg]",
    "imut": "mutual inclination of two orbits [deg]",
    "dE_over_E": "relative energy error (quality of the run)",
    "dL_over_L": "relative angular-momentum error",
    "M": "stellar mass [Msun]",
    "R": "stellar radius [Rsun]",
    "Teff": "effective temperature [K]",
}


def _csv_header(d: dict, cols: list[str], step: int, n: int) -> list[str]:
    """The '#' lines of the CSV: provenance, units, time system and the meaning of the columns."""
    from . import __version__

    meta = d.get("meta", {})
    lines = [
        f"SWINGERsim {meta.get('swingersim_version', __version__)} timeseries ("
        + (f"REBOUND {meta['rebound_version']}" if meta.get("rebound_version") else "")
        + (f", REBOUNDx {meta['reboundx_version']}" if meta.get("reboundx_version") else "")
        + ")",
        (
            "time: t_yr in Julian years (365.25 d) from the start of the run (t = 0 = the epoch"
            " of the configured elements); negative for backward runs"
        ),
        (
            "units: lengths AU, masses Msun, radii Rsun, angles deg; elements osculating, Jacobi;"
            " sky frame x-y = sky plane"
        ),
    ]
    if step > 1:
        lines.append(
            f"rows: every {step}-th of {n} outputs and the last one (all in the HDF5 file)"
        )
    if meta.get("event"):
        lines.append(f"event: {meta['event']} (the last row is the state at the event)")
    lines.append("columns (<quantity>_<orbit or star>_<unit>):")
    seen = set()
    for c in cols:
        key = c.split("_")[0] if c not in COLUMN_HELP else c
        if c.startswith("i_inv"):
            key = "i_inv"
        elif c.startswith("Omega_inv"):
            key = "Omega_inv"
        if key in COLUMN_HELP and key not in seen:
            seen.add(key)
            lines.append(f"  {key}: {COLUMN_HELP[key]}")
    cfg = meta.get("config_yaml", "")
    if cfg:
        lines.append("configuration:")
        lines += [f"  {line}" for line in cfg.rstrip().splitlines()]
    return [f"# {line}" for line in lines]


def write_timeseries(d: dict, path: Path, max_rows: int | None = CSV_MAX_ROWS) -> Path:
    """The CSV table of the run. With more outputs than ``max_rows`` (None or 0: no limit) only
    every k-th output and the last one are written."""
    cols: dict[str, np.ndarray] = {"t_yr": d["t"]}
    units = {
        "a": "au",
        "e": "",
        "inc": "deg",
        "inc_inv": "deg",
        "Omega": "deg",
        "Omega_inv": "deg",
        "omega": "deg",
    }
    for name, g in d["orbits"].items():
        for k, unit in units.items():
            if k not in g:
                continue
            v = g[k]
            if unit == "deg":
                v = np.degrees(v)
                if k != "inc" and k != "inc_inv":
                    v = v % 360.0
            label = {"inc": "i", "inc_inv": "i_inv"}.get(k, k)
            cols[f"{label}_{name}" + (f"_{unit}" if unit else "")] = v
    for name, v in d["mutual"].items():
        cols[f"imut_{name}_deg"] = v
    if d["energy_error"] is not None:
        cols["dE_over_E"] = d["energy_error"]
    if d["angmom_error"] is not None:
        cols["dL_over_L"] = d["angmom_error"]
    for star, g in d["stars"].items():
        for k, label in (("mass_msun", "M"), ("radius_rsun", "R"), ("teff_k", "Teff")):
            if k in g:
                cols[f"{label}_{star}_{k.split('_', 1)[1]}"] = g[k]
    n = len(d["t"])
    cols = {k: v for k, v in cols.items() if len(v) == n}
    step = int(np.ceil(n / max_rows)) if max_rows and n > max_rows else 1
    rows = np.arange(0, n, step)
    if step > 1 and rows[-1] != n - 1:
        rows = np.append(rows, n - 1)
    header = "\n".join([*_csv_header(d, list(cols), step, n), ",".join(cols)])
    np.savetxt(
        path,
        np.column_stack([np.asarray(v)[rows] for v in cols.values()]),
        delimiter=",",
        header=header,
        comments="",
        fmt="%.10g",
    )
    return path


def _time_axis(t: np.ndarray) -> tuple[np.ndarray, str]:
    span = float(np.nanmax(np.abs(t))) if len(t) else 0.0
    if span >= 2e6:
        return t / 1e6, "time [Myr]"
    if span >= 2e3:
        return t / 1e3, "time [kyr]"
    return t, "time [yr]"


def write_plots(d: dict, paths: dict[str, Path], title: str) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = []
    x, xl = _time_axis(d["t"])
    names = list(d["orbits"])
    fig, axes = plt.subplots(3, 1, figsize=(8, 7.5), sharex=True)
    for k, name in enumerate(names):
        g, c = d["orbits"][name], PALETTE[k % len(PALETTE)]
        if "a" in g:
            a0 = g["a"][0] if g["a"][0] else 1.0
            axes[0].plot(x, g["a"] / a0, lw=0.9, color=c, label=name)
        if "e" in g:
            axes[1].plot(x, g["e"], lw=0.9, color=c, label=name)
        if "inc" in g:
            axes[2].plot(x, np.degrees(g["inc"]), lw=0.9, color=c, label=name)
    axes[0].set_ylabel("a / a(0)")
    axes[1].set_ylabel("eccentricity")
    axes[2].set_ylabel("inclination to the sky [deg]")
    axes[2].set_xlabel(xl)
    axes[0].legend(loc="best", fontsize=8, ncol=min(4, len(names)) or 1)
    axes[0].set_title(f"{title}: orbits over time")
    for ax in axes:
        ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(paths["orbits"], dpi=130)
    plt.close(fig)
    out.append(paths["orbits"])

    if d["mutual"]:
        fig, ax = plt.subplots(figsize=(8, 3.6))
        for k, (name, v) in enumerate(d["mutual"].items()):
            ax.plot(x, v, lw=0.9, color=PALETTE[k % len(PALETTE)], label=name)
        ax.set_ylabel("mutual inclination [deg]")
        ax.set_xlabel(xl)
        ax.set_title(f"{title}: mutual inclinations")
        ax.legend(loc="best", fontsize=8)
        ax.grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(paths["inclinations"], dpi=130)
        plt.close(fig)
        out.append(paths["inclinations"])

    errs = [
        (lbl, v)
        for lbl, v in (
            ("energy |dE/E|", d["energy_error"]),
            ("angular momentum |dL/L|", d["angmom_error"]),
        )
        if v is not None and np.any(np.isfinite(v))
    ]
    if errs:
        fig, ax = plt.subplots(figsize=(8, 3.4))
        for k, (lbl, v) in enumerate(errs):
            ax.semilogy(x, np.maximum(np.abs(v), 1e-17), lw=0.9, color=PALETTE[k], label=lbl)
        ax.set_ylabel("relative error")
        ax.set_xlabel(xl)
        ax.set_title(f"{title}: conservation (below 1e-6 is a clean run)")
        ax.legend(loc="best", fontsize=8)
        ax.grid(alpha=0.25, which="both")
        fig.tight_layout()
        fig.savefig(paths["conservation"], dpi=130)
        plt.close(fig)
        out.append(paths["conservation"])

    evolving = {
        s: g
        for s, g in d["stars"].items()
        if "mass_msun" in g and len(g["mass_msun"]) == len(d["t"])
    }
    if evolving:
        fig, axes = plt.subplots(3, 1, figsize=(8, 7.5), sharex=True)
        for k, (s, g) in enumerate(evolving.items()):
            c = PALETTE[k % len(PALETTE)]
            axes[0].plot(x, g["mass_msun"], lw=0.9, color=c, label=s)
            if "radius_rsun" in g:
                axes[1].semilogy(x, g["radius_rsun"], lw=0.9, color=c, label=s)
            if "teff_k" in g:
                axes[2].plot(x, g["teff_k"], lw=0.9, color=c, label=s)
        axes[0].set_ylabel("mass [Msun]")
        axes[1].set_ylabel("radius [Rsun]")
        axes[2].set_ylabel("Teff [K]")
        axes[2].set_xlabel(xl)
        axes[0].legend(loc="best", fontsize=8, ncol=min(5, len(evolving)))
        axes[0].set_title(f"{title}: the stars over time")
        for ax in axes:
            ax.grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(paths["stars"], dpi=130)
        plt.close(fig)
        out.append(paths["stars"])
    return out


def write_products(
    h5: Path,
    summary: str,
    viewer: bool = True,
    plots: bool = True,
    timeseries: bool = True,
    csv_max_rows: int | None = CSV_MAX_ROWS,
) -> list[Path]:
    """Write the products of a saved run; returns the files written. A failing optional product
    (viewer without plotly, figures without matplotlib) is reported in the summary, not raised.
    The HTML report is not listed here: swingersim.report adds it once it has been written."""
    h5 = Path(h5)
    paths = product_paths(h5)
    d = _read(h5)
    written, notes = [], []
    if timeseries:
        written.append(write_timeseries(d, paths["timeseries"], csv_max_rows))
    if plots:
        try:
            written += write_plots(d, paths, h5.stem)
        except ImportError:
            notes.append("no PNG figures: matplotlib is not installed")
    if viewer:
        try:
            from .viewer import write_html

            written.append(write_html(h5, paths["viewer"]))
        except ImportError:
            notes.append("no 3D viewer: plotly is not installed (pip install -e '.[viz]')")
    lines = [
        summary.rstrip(),
        "",
        "files of this run:",
        f"  {h5.name}  (all data, HDF5)",
    ]
    lines += [f"  {p.name}" for p in written]
    lines += [f"  ({n})" for n in notes]
    paths["summary"].write_text("\n".join(lines) + "\n", encoding="utf-8")
    return [paths["summary"], *written]


def add_report_to_summary(h5: Path, report: Path) -> None:
    """List the HTML report among the files of the run (called once the report exists)."""
    summary = product_paths(Path(h5))["summary"]
    if not summary.exists():
        return
    report = Path(report)
    name = report.name if report.parent.resolve() == Path(h5).parent.resolve() else str(report)
    line = f"  {name}  (report, opens in a web browser)"
    text = summary.read_text(encoding="utf-8")
    if line in text:
        return
    anchor = f"  {Path(h5).name}  (all data, HDF5)\n"
    text = text.replace(anchor, anchor + line + "\n", 1) if anchor in text else text + line + "\n"
    summary.write_text(text, encoding="utf-8")


__all__ = ["add_report_to_summary", "product_paths", "write_products", "write_timeseries"]
