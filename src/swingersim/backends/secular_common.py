"""Helpers shared by the two secular backends (backends/secular.py for 2+2 quadruples and
backends/secular_nested.py for any nested hierarchy): extremes from a dense solution, first
threshold crossings in the direction of the run, the end of a run at a Roche-lobe overflow,
the invariable plane and the common post-processing of the result."""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from .. import massloss, masstransfer, units
from .. import masstransfer as _mt
from ..diagnostics.binary import binary_diagnostics
from ..diagnostics.extremes import THRESHOLDS_1ME, threshold_key
from ..diagnostics.zlk import analyse, system_diagnostics
from . import evolving as ev
from .chunked import merge_extremes, merge_roche
from .nbody import parent_of, spin_vector


def first_crossing(t: np.ndarray, f: np.ndarray, fun, level: float, forward: bool = True) -> float:
    """First time, in the direction of the run, at which f(t) exceeds level (NaN if never),
    refined on the dense solution fun(t). t ascending; forward = False for a backward run
    (t < 0), whose first crossing is the one closest to the start t = 0, i.e. the latest."""
    above = np.flatnonzero(f > level)
    if not len(above):
        return _first_peak_crossing(t, f, fun, level, forward)
    if forward:
        k, prev = int(above[0]), int(above[0]) - 1
    else:
        k, prev = int(above[-1]), int(above[-1]) + 1
    if prev < 0 or prev >= len(t):  # already above at the start of the interval
        return float(t[k])
    lo, hi = sorted((t[prev], t[k]))
    try:
        return float(brentq(lambda x: fun(x) - level, lo, hi, xtol=1e-9 * abs(t[k]) + 1e-12))
    except ValueError:
        return float(t[k])


def _first_peak_crossing(t, f, fun, level: float, forward: bool, n_peaks: int = 20) -> float:
    """No sample above level: the highest local maxima of the samples are refined on the dense
    solution, and the crossing before the earliest one (in the direction of the run) that
    exceeds level is returned (NaN if none does)."""
    if len(t) < 3:
        return float("nan")
    c, lft, rgt = f[1:-1], f[:-2], f[2:]
    inner = np.flatnonzero((c >= lft) & (c >= rgt) & ((c > lft) | (c > rgt))) + 1  # not flat
    if not len(inner):
        return float("nan")
    cands = inner[np.argsort(f[inner])[::-1][:n_peaks]]
    hits = []
    for k in cands:
        lo, hi = float(t[k - 1]), float(t[k + 1])
        r = minimize_scalar(
            lambda x: -fun(x), bounds=(lo, hi), method="bounded",
            options={"xatol": 1e-9 * (hi - lo) + 1e-12},
        )  # fmt: skip
        if -r.fun > level:
            a, b = (lo, float(r.x)) if forward else (float(r.x), hi)
            try:
                tc = brentq(lambda x: fun(x) - level, a, b, xtol=1e-9 * abs(r.x) + 1e-12)
            except ValueError:
                tc = float(r.x)
            hits.append(float(tc))
    if not hits:
        return float("nan")
    return min(hits) if forward else max(hits)


def refine(t: np.ndarray, f: np.ndarray, g, n_peaks: int = 5) -> tuple[float, float]:
    """Global maximum of a sampled function f(t): the n_peaks highest local maxima are refined
    by minimising g (= -f evaluated on the dense solution) between their neighbours.
    Returns (t, g(t)) of the best one."""
    if len(t) < 3:
        k = int(np.argmax(f))
        return float(t[k]), float(g(t[k]))
    inner = np.flatnonzero((f[1:-1] >= f[:-2]) & (f[1:-1] >= f[2:])) + 1
    cands = np.concatenate([inner, [0, len(t) - 1]])
    cands = cands[np.argsort(f[cands])[::-1][:n_peaks]]
    best_t, best_g = float(t[cands[0]]), float(g(t[cands[0]]))
    for k in cands:
        lo, hi = t[max(k - 1, 0)], t[min(k + 1, len(t) - 1)]
        if hi > lo:
            r = minimize_scalar(
                g, bounds=(lo, hi), method="bounded", options={"xatol": 1e-9 * (hi - lo) + 1e-12}
            )
            if r.fun < best_g:
                best_t, best_g = float(r.x), float(r.fun)
        gk = float(g(t[k]))
        if gk < best_g:
            best_t, best_g = float(t[k]), gk
    return best_t, best_g


def dense_extremes(
    sol,
    times,
    orbits: list[tuple[str, int]],
    a_of,
    lobes: dict[str, list[tuple[str, float, float]]],
    n_fine: int | None = None,
    t_stop: float | None = None,
    overflow: dict | None = None,
):
    """Extremes from the dense solution of a run (or of a chunk of it): every integrator step
    plus a uniform grid, then the highest local maxima are refined on the dense interpolant (a
    high-eccentricity spike can be much shorter than the output spacing).

    orbits: (name, index of the first component of the Laplace vector in the state) of every
    orbit; a_of(Y, q): semi-major axis of orbit q for the state Y (one state, or states in the
    columns of a 2-D array); lobes[name]: (star, Roche-lobe radius in units of the periastron
    distance, stellar radius [AU]) of every star with a lobe in that orbit. t_stop: only the
    part from times[0] to t_stop. overflow: when a dict, it receives for every star the first
    time (in the direction of the run) at which it fills its Roche lobe (NaN if never)."""
    t_a = float(times[0])
    t_b = float(times[-1]) if t_stop is None else float(t_stop)
    forward = float(times[-1]) >= float(times[0])
    lo, hi = min(t_a, t_b), max(t_a, t_b)
    n_fine = n_fine or max(20 * len(times), 2000)
    steps = np.asarray(sol.t, float)
    fine = np.union1d(np.linspace(t_a, t_b, n_fine), steps[(steps >= lo) & (steps <= hi)])
    Yf = sol.sol(fine)
    extremes, roche = {}, {}
    for q, (name, sl) in enumerate(orbits):

        def ecc(t, sl=sl):
            return float(np.linalg.norm(sol.sol(t)[sl : sl + 3]))

        def peri(t, sl=sl, q=q):
            y = sol.sol(t)
            return float(a_of(y, q) * (1.0 - np.linalg.norm(y[sl : sl + 3])))

        e_f = np.linalg.norm(Yf[sl : sl + 3], axis=0)
        a_f = np.broadcast_to(np.asarray(a_of(Yf, q), float), e_f.shape)
        q_f = a_f * (1.0 - e_f)
        t_e, e_max = refine(fine, e_f, lambda t, ecc=ecc: -ecc(t))
        t_q, q_min = refine(fine, -q_f, peri)
        ex = {"e_max": -e_max, "t_e_max_yr": t_e, "q_min_au": q_min, "t_q_min_yr": t_q}
        for thr in THRESHOLDS_1ME:
            t_cross = first_crossing(fine, e_f, ecc, 1.0 - thr, forward)
            if np.isnan(t_cross) and -e_max > 1.0 - thr:
                t_cross = t_e  # a peak shorter than the sampling, found by the refinement
            ex[threshold_key(thr)] = t_cross
        extremes[name] = ex
        for star, rl, radius in lobes.get(name, ()):
            fill = radius / (rl * q_f)
            kmax = int(np.argmax(fill))
            fill_q = radius / (rl * q_min)
            roche[star] = {
                "fill_max": float(max(fill[kmax], fill_q)),
                "t_fill_max_yr": float(t_q if fill_q >= fill[kmax] else fine[kmax]),
            }
            if overflow is not None:

                def fill_at(t, peri=peri, radius=radius, rl=rl):
                    return radius / (rl * peri(t))

                t_ov = first_crossing(fine, fill, fill_at, 1.0 - 1e-15, forward)
                if np.isnan(t_ov) and roche[star]["fill_max"] >= 1.0:
                    t_ov = roche[star]["t_fill_max_yr"]  # a peak between the samples
                overflow[star] = t_ov
    return extremes, roche


def _first_overflow(overflow: dict, sign: float):
    """(time, stars) of the first Roche-lobe overflow of a chunk, or None."""
    hits = {s: t for s, t in overflow.items() if not math.isnan(t)}
    if not hits:
        return None
    t_ov = min(hits.values(), key=lambda x: sign * x)
    return t_ov, sorted(s for s, t in hits.items() if sign * t <= sign * t_ov)


def integrate_until_overflow(integrate, y0, times, extremes_of, stop: bool, chunk: int = 10):
    """Chunked integration as backends/chunked.integrate_chunked (bounded memory) that, with
    stop = True (run.stop_on_roche_overflow), ends the run at the first Roche-lobe overflow
    like the N-body backend: the outputs after it are dropped (at least two are kept) and the
    extremes cover the run up to it.

    extremes_of(sol, t_chunk, t_stop, overflow) gives the (extremes, roche) of a chunk (see
    dense_extremes). Returns (times, Y, extremes, roche, overflow) with overflow = (time,
    stars) or None."""
    times = np.asarray(times, float)
    sign = 1.0 if times[-1] >= times[0] else -1.0
    Y = np.empty((len(y0), len(times)))
    Y[:, 0] = y0
    extremes: dict = {}
    roche: dict = {}
    y = np.array(y0, float)
    for i0 in range(0, len(times) - 1, chunk):
        i1 = min(i0 + chunk, len(times) - 1)
        t_chunk = times[i0 : i1 + 1]
        sol = integrate(y, t_chunk)
        if sol.status < 0:
            raise RuntimeError(f"secular integration failed: {sol.message}")
        Y[:, i0 + 1 : i1 + 1] = sol.sol(t_chunk[1:])
        Y[:, i1] = sol.y[:, -1]  # the chunk ends exactly at the output
        over: dict | None = {} if stop else None
        ex, ro = extremes_of(sol, t_chunk, None, over)
        hit = _first_overflow(over, sign) if stop else None
        if hit is not None:
            keep = max(int(np.sum(sign * times <= sign * hit[0])), 2)
            t_last = times[keep - 1] if sign * times[keep - 1] > sign * hit[0] else hit[0]
            ex, ro = extremes_of(sol, t_chunk, t_last, None)
            merge_extremes(extremes, ex)
            merge_roche(roche, ro)
            return times[:keep], Y[:, :keep], extremes, roche, hit
        merge_extremes(extremes, ex)
        merge_roche(roche, ro)
        y = sol.y[:, -1].copy()
        del sol
    return times, Y, extremes, roche, None


def extremes_until_overflow(sol, times, Y, extremes_of, stop: bool):
    """The same for a run integrated in one piece (dense solution sol of the whole run):
    (times, Y, extremes, roche, overflow)."""
    times = np.asarray(times, float)
    sign = 1.0 if times[-1] >= times[0] else -1.0
    over: dict | None = {} if stop else None
    ex, ro = extremes_of(sol, times, None, over)
    hit = _first_overflow(over, sign) if stop else None
    if hit is None:
        return times, Y, ex, ro, None
    keep = max(int(np.sum(sign * times <= sign * hit[0])), 2)
    t_last = times[keep - 1] if sign * times[keep - 1] > sign * hit[0] else hit[0]
    ex, ro = extremes_of(sol, times, t_last, None)
    return times[:keep], Y[:, :keep], ex, ro, hit


def truncate_at_stream_stop(times, Y, stream):
    """Evolving stars: the run ends at the Roche-lobe overflow found by the stream; the last row
    holds the state at that moment (the stream stores it in the rows after the overflow)."""
    sign = 1.0 if times[-1] >= 0 else -1.0
    keep = max(int(np.sum(sign * times <= sign * stream.t_stop)), 1)
    if keep < len(times) and (
        keep == 1 or sign * (stream.t_stop - times[keep - 1]) > 1e-9 * abs(times[-1])
    ):
        times = np.append(times[:keep], stream.t_stop)
        keep += 1
    stream.snapshots, stream.h_ref = stream.snapshots[:keep], stream.h_ref[:keep]
    return times[:keep], Y[:, :keep]


def invariable_basis(L0: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Unit vectors (x, y, z) of the invariable plane: z along the total angular momentum L0,
    x along the line of nodes on the reference plane (or the reference x for z = e_z)."""
    z = L0 / np.linalg.norm(L0)
    x = np.cross([0.0, 0.0, 1.0], z)
    x = x / np.linalg.norm(x) if np.linalg.norm(x) > 1e-12 else np.array([1.0, 0.0, 0.0])
    return x, np.cross(z, x), z


def invariable_angles(k: np.ndarray, basis) -> tuple[float, float]:
    """(inclination, longitude of the node) of the orbit normal k in the invariable plane."""
    x, yv, z = basis
    return (
        float(np.arccos(np.clip(k @ z, -1.0, 1.0))),
        float(np.mod(np.arctan2(k @ x, -(k @ yv)), 2.0 * np.pi)),
    )


def finish_result(
    cfg,
    res,
    *,
    root,
    masses,
    m0,
    props,
    pos0,
    vel0,
    evolving: bool,
    mass_loss: bool,
    stream=None,
    stop_event: str | None = None,
    psi_max: dict | None = None,
    overflow=None,
) -> None:
    """Post-processing common to the secular backends: series of the evolving stars and of
    the masses, the events (end of the stellar evolution, Roche-lobe overflow), the
    diagnostics of the pairs, of the system and of the binaries."""
    if evolving:
        res.stellar = ev.stellar_series(cfg, masses, res.times) if cfg.physics.evolution else {}
        if mass_loss:
            ms = massloss.mass_series(cfg, masses, res.times)
            for k in range(len(masses)):
                res.stellar.setdefault(f"m{k + 1}", {})["mass_msun"] = ms[k]
            res.mass_loss = True
        if stop_event is not None:
            res.event = stop_event
            res.event_time = float(stream.t_stop if stream.t_stop is not None else res.times[-1])
            res.event_detail = {"age_myr": cfg.stellar_model.age_myr + res.event_time / 1e6}
            if stream.t_stop is not None:
                res.event_detail["stars"] = list(getattr(stream, "stop_stars", None) or []) or [
                    s for s, r in stream.roche.items() if r["fill_max"] >= 1.0
                ]
        over = [s for s, r in res.roche.items() if r["fill_max"] >= 1.0]
        # (a stream checks every overflow itself, sparing those of a modelled mass transfer)
        if over and _mt.stop_on_overflow(cfg) and stream is None:
            res.event = "roche_overflow"
            res.event_time = min(res.roche[s]["t_fill_max_yr"] for s in over)
            res.event_detail = {"stars": over}
    elif overflow is not None:  # run.stop_on_roche_overflow: the run ended there
        res.event, res.event_time = "roche_overflow", float(overflow[0])
        res.event_detail = {"stars": list(overflow[1])}
    radii_au = np.array([p.radius_rsun * units.RSUN_AU for p in props])
    res.diagnostics = analyse(root, m0, cfg.physics.gr, res, (pos0, vel0), radii_au)
    res.system = system_diagnostics(root, m0, res.diagnostics, (pos0, vel0))
    if psi_max is not None:
        res.system.update({f"mass_loss_psi_max_{n}": v for n, v in psi_max.items()})
    res.system.update(masstransfer.system_entries(cfg))
    res.system.update(masstransfer.eccentric_transfer(cfg, res))
    spins0 = [
        spin_vector(cfg, f"m{k + 1}", parent_of(root, k), pos0, vel0, m0)
        if cfg.physics.tides
        else np.zeros(3)
        for k in range(len(masses))
    ]
    res.binaries = binary_diagnostics(
        root, m0, props, spins0, pos0, vel0, res, cfg.physics.gr, cfg.physics.tides, cfg=cfg
    )
    res.backend = "secular"
