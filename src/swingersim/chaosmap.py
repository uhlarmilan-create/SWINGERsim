"""Chaos maps: MEGNO on a grid of two (or one) configuration parameters.

MEGNO (Mean Exponential Growth factor of Nearby Orbits, Cincotta & Simo 2000) tends to 2 for
quasi-periodic (regular) motion, grows linearly with time for chaotic motion and is 0 for
stable periodic orbits. A map of its value after a fixed number of outer orbits shows the
stable regions, the chaotic zones and the islands of mean-motion resonances (e.g. the 3:2
resonance of the two binaries of a 2+2 quadruple).

Each grid point is a Newtonian point-mass N-body integration of the configuration with the two
parameters replaced (paths as in ensemble.vary, e.g. orbits.out.P, orbits.B.e, masses.m3).
Relativity and the structure of the stars are left out because REBOUND's variational equations
cover only Newtonian gravity. A point stops early (and counts as unstable) when a body escapes
beyond `escape_factor` times the largest initial semi-major axis or two stars collide.

Classes of the grid points (summarize_map): regular and chaotic (MEGNO after the full time),
escape or collision (a physical instability), wall_time (stopped by the computing-time limit
before t_max: its MEGNO is not that of the full time, NaN in megno_array) and invalid / error
(the configuration could not be set up or the integration failed).
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
import rebound

from . import units
from .config import Config
from .ensemble import _set, read_jsonl, resolve_processes, rewrite_jsonl, row_line
from .hierarchy import build_tree
from .parallel import run_tasks

MEGNO_CHAOTIC = 4.0  # usual threshold: MEGNO well above 2 = chaotic


def periods(cfg: Config) -> dict[str, float]:
    """Orbital period [yr] of every orbit of the configuration."""
    root, masses = build_tree(cfg)
    from .hierarchy import mass_of

    out = {}
    for node in root.nodes():
        m = mass_of(node.left, masses) + mass_of(node.right, masses)
        out[node.name] = units.sma_to_period(node.elements.a, m)
    return out


def megno_point(cfg_dict: dict, t_max_yr: float, max_wall_s: float, escape_factor: float):
    """MEGNO of one configuration after t_max_yr (Newtonian, point masses)."""
    from . import massloss
    from .backends.nbody import MEGNO_SEED, build_simulation, lyapunov_time_yr

    t0 = time.time()
    try:
        # initial masses as in run_nbody (mass_loss {source: track, masses: current})
        cfg = massloss.resolve_masses(Config.model_validate(cfg_dict))
        sim, _, _, root, _ = build_simulation(cfg, gr=False, structure=False)
    except ValueError as exc:  # invalid drawn configuration
        return {"megno": math.nan, "lyapunov_yr": math.nan, "status": f"invalid: {exc}"}
    a_max = max(n.elements.a for n in root.nodes())
    sim.collision = "direct" if cfg.run.collisions else "none"
    sim.exit_max_distance = escape_factor * a_max
    sim.init_megno(seed=MEGNO_SEED)  # fixed: the same map on every run
    status = "ok"
    steps = 200
    try:
        for k in range(1, steps + 1):
            sim.integrate(t_max_yr * k / steps, exact_finish_time=0)
            if time.time() - t0 > max_wall_s:
                status = "wall_time"
                break
    except rebound.Escape:
        status = "escape"
    except rebound.Collision:
        status = "collision"
    megno = sim.megno() if status in ("ok", "wall_time") else math.nan
    return {
        "megno": float(megno),
        "lyapunov_yr": lyapunov_time_yr(sim),
        "t_end_yr": float(sim.t),
        "status": status,
    }


_JOB: dict = {}  # the map being computed, set once per worker process by _init


def _init(base: dict, x: str, y: str | None, t_max: float, wall: float, esc: float) -> None:
    """Pool initializer: the base configuration and the settings reach every worker process
    once; the tasks carry only (i, j, x, y)."""
    _JOB.update(base=base, x=x, y=y, t_max=t_max, wall=wall, esc=esc)


def _point_config(base: dict, x: str, xv: float, y: str | None, yv: float) -> dict:
    data = json.loads(json.dumps(base, default=str))
    _set(data, x, float(xv))
    if y is not None:
        _set(data, y, float(yv))
    return data


def _task(args):
    i, j, xv, yv = args
    job = _JOB
    try:
        cfg_dict = _point_config(job["base"], job["x"], xv, job["y"], yv)
        r = megno_point(cfg_dict, job["t_max"], job["wall"], job["esc"])
    except Exception as exc:  # noqa: BLE001  (keep the map going)
        r = {"megno": math.nan, "lyapunov_yr": math.nan, "status": f"error: {exc!r}"}
    return {"i": i, "j": j, "x": xv, "y": yv, **r}


def _lost(task, message: str) -> dict:
    i, j, xv, yv = task
    return {"i": i, "j": j, "x": xv, "y": yv, "megno": math.nan, "lyapunov_yr": math.nan,
            "status": f"error: {message}"}  # fmt: skip


def grid_values(lo: float, hi: float, n: int, log: bool) -> np.ndarray:
    if n == 1:
        return np.array([lo])
    return np.geomspace(lo, hi, n) if log else np.linspace(lo, hi, n)


def chaos_map(
    cfg: Config,
    out: Path,
    x: str,
    xvals: np.ndarray,
    y: str | None = None,
    yvals: np.ndarray | None = None,
    t_orbits: float = 1000.0,
    t_max_yr: float | None = None,
    processes: int = 1,
    max_wall_s: float = 300.0,
    escape_factor: float = 20.0,
    resume: bool = False,
    progress=None,
) -> list[dict]:
    """Compute the map; rows are written to out/grid.jsonl as they finish. A resumed map must
    have the same grid, time, escape factor and base configuration."""
    import platform

    from . import __version__

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    base = cfg.model_dump()
    base.pop("ensemble", None)
    if t_max_yr is None:
        t_max_yr = t_orbits * max(periods(cfg).values())
    yvals = np.array([math.nan]) if y is None else np.asarray(yvals)
    meta = {
        "x": x,
        "y": y,
        "xvals": [float(v) for v in xvals],
        "yvals": [float(v) for v in yvals],
        "t_max_yr": float(t_max_yr),
        "escape_factor": escape_factor,
        "config": json.loads(json.dumps(base, default=str)),
        "swingersim": __version__,
        "rebound": rebound.__version__,
        "python": platform.python_version(),
    }
    rows_file, meta_file = out / "grid.jsonl", out / "chaosmap.json"
    done = {}
    if rows_file.exists():
        if not resume:
            raise FileExistsError(f"{rows_file} exists; use --resume or another folder")
        old = json.loads(meta_file.read_text())
        keys = ("x", "y", "xvals", "yvals", "t_max_yr")
        if json.dumps([old[k] for k in keys]) != json.dumps([meta[k] for k in keys]):  # NaN-safe
            raise ValueError("the grid changed since the map was started")
        if old.get("escape_factor", escape_factor) != escape_factor:
            raise ValueError("escape_factor changed since the map was started")
        if "config" in old and json.dumps(old["config"], sort_keys=True) != json.dumps(
            meta["config"], sort_keys=True
        ):
            raise ValueError("the configuration changed since the map was started")
        raw, bad = read_jsonl(rows_file)
        for r in raw:
            done[(r["i"], r["j"])] = r
        if bad or len(done) < len(raw):  # a torn last line (or duplicates): rewrite
            rewrite_jsonl(rows_file, list(done.values()))
    meta_file.write_text(json.dumps(meta, indent=1))
    tasks = [
        (i, j, float(xv), float(yv))
        for i, xv in enumerate(xvals)
        for j, yv in enumerate(yvals)
        if (i, j) not in done
    ]
    rows = list(done.values())
    total = len(xvals) * len(yvals)
    init = (base, x, y, float(t_max_yr), float(max_wall_s), float(escape_factor))
    with rows_file.open("a") as fh:
        for r in run_tasks(_task, tasks, resolve_processes(processes), _lost, _init, init):
            fh.write(row_line(r))
            fh.flush()
            rows.append(r)
            if progress:
                progress(
                    f"  [{len(rows)}/{total}] x = {r['x']:.6g}, y = {r['y']:.6g}: "
                    f"MEGNO {r['megno']:.3g} ({r['status']})"
                )
    rows.sort(key=lambda r: (r["i"], r["j"]))
    return rows


def point_class(r: dict) -> str:
    """regular, chaotic, escape (escape or collision), wall_time or error (invalid set-up or a
    failed integration) of one grid point."""
    st = str(r.get("status", ""))
    if st in ("escape", "collision"):
        return "escape"
    if st == "wall_time":
        return "wall_time"
    m = r.get("megno")
    if st != "ok" or m is None or not math.isfinite(m):
        return "error"
    return "regular" if m < MEGNO_CHAOTIC else "chaotic"


def megno_array(rows: list[dict]) -> np.ndarray:
    """MEGNO on the grid: inf for an escape or collision, NaN for points stopped by the
    computing-time limit (their MEGNO is not that of the full time), failed set-ups and errors."""
    ni = max(r["i"] for r in rows) + 1
    nj = max(r["j"] for r in rows) + 1
    a = np.full((ni, nj), np.nan)
    for r in rows:
        c = point_class(r)
        if c in ("regular", "chaotic"):
            a[r["i"], r["j"]] = r["megno"]
        elif c == "escape":
            a[r["i"], r["j"]] = np.inf
    return a


def summarize_map(meta: dict, rows: list[dict]) -> str:
    from collections import Counter

    a = megno_array(rows)
    n = a.size
    count = Counter(point_class(r) for r in rows)
    lines = [
        f"chaos map: {meta['x']}"
        + (f" x {meta['y']}" if meta["y"] else "")
        + f", {n} points, {meta['t_max_yr']:.6g} yr each",
        f"  regular (MEGNO < {MEGNO_CHAOTIC:g}): {count['regular']}",
        f"  chaotic (MEGNO >= {MEGNO_CHAOTIC:g}): {count['chaotic']}",
        f"  escape or collision: {count['escape']}",
    ]
    if count["wall_time"]:
        lines.append(
            f"  stopped by the computing-time limit (--max-wall) before t_max: "
            f"{count['wall_time']} (not classified)"
        )
    if count["error"]:
        lines.append(f"  invalid configuration or failed integration: {count['error']}")
    xs = np.array(meta["xvals"])
    x_st = stable_from(xs, a)
    if not math.isnan(x_st):
        lines.append(
            f"  smallest {meta['x']} with regular motion from there outwards"
            + (f" for every {meta['y']}" if meta["y"] else "")
            + f": {x_st:.6g}"
        )
    return "\n".join(lines)


def stable_from(xs: np.ndarray, megno: np.ndarray) -> float:
    """Smallest x beyond which every grid point is regular; for a 2-D map, regular for every
    value of the second parameter (e.g. all initial phases)."""
    m = megno.reshape(len(xs), -1)
    reg = (np.isfinite(m) & (m < MEGNO_CHAOTIC)).all(axis=1)
    for k in range(len(xs)):
        if reg[k:].all():
            return float(xs[k])
    return math.nan


def plot_map(meta: dict, rows: list[dict], path: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    a = megno_array(rows)
    xs = np.array(meta["xvals"])
    fig, ax = plt.subplots(figsize=(7, 5))
    if meta["y"] is None:
        m = np.where(np.isinf(a[:, 0]), np.nan, a[:, 0])
        ax.plot(xs, np.clip(m, 0, 10), "o-", color="#3b6ea5", ms=3)
        bad = np.isinf(a[:, 0])
        ax.plot(xs[bad], np.full(bad.sum(), 10), "x", color="#c0392b", label="escape / collision")
        ax.axhline(2, color="0.5", lw=0.8, ls="--")
        ax.axhline(MEGNO_CHAOTIC, color="0.5", lw=0.8, ls=":")
        ax.set_ylabel("MEGNO (clipped at 10)")
        if bad.any():
            ax.legend(frameon=False)
    else:
        ys = np.array(meta["yvals"])
        m = np.where(np.isinf(a), 10.0, a)
        im = ax.pcolormesh(
            xs, ys, np.clip(m.T, 0, 10), cmap="RdYlBu_r", vmin=0, vmax=10, shading="nearest"
        )
        fig.colorbar(im, ax=ax, label="MEGNO (2 = regular, > 4 chaotic, 10 = escape)")
        ax.set_ylabel(meta["y"])
    ax.set_xlabel(meta["x"])
    ax.set_title(f"MEGNO after {meta['t_max_yr']:.4g} yr")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return Path(path)
