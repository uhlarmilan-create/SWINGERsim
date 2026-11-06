"""Evolving stars in the secular backends (physics.evolution).

The stars change on nuclear time scales (Myr-Gyr), much more slowly than the secular dynamics in
most systems, so the coupling uses operator splitting as in MSE (Hamers et al. 2021): the run is
cut into segments within which no stellar radius changes by more than ``run.evolution_max_dr``
(default 1 %); in each segment the stellar structure (radius, apsidal constant, moment of
inertia, tidal lag) is held at its value in the middle of the segment, and the orbital state is
carried over from segment to segment. Spins conserve their angular momentum across a change of
the moment of inertia. The run stops at helium ignition of any star (end of the track grid,
event ``stellar_evolution_end``) and, backwards in time, at the birth of the stars (age 0,
event ``stellar_birth``).
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, replace
from itertools import pairwise

import numpy as np

from .. import massloss, masstransfer, units
from ..evolution import EvolvedBeyondGrid
from ..massloss import YR_PER_MYR
from ..stellar import StarProperties, age_offset, star_properties, stellar_track


def time_varying(cfg) -> bool:
    """Evolving stars or changing masses: the secular backends integrate in segments."""
    return bool(cfg.physics.evolution) or massloss.enabled(cfg)


def _age0(cfg) -> float:
    return float(cfg.stellar_model.age_myr)


def props_at(cfg, masses, t_yr: float) -> list[StarProperties]:
    """Stellar properties at run time t_yr (age = age_myr + t)."""
    from copy import copy

    c = copy(cfg)
    c.stellar_model = cfg.stellar_model.model_copy(
        update={"age_myr": _age0(cfg) + t_yr / YR_PER_MYR}
    )
    mt_on = masstransfer.enabled(cfg)
    out = []
    for k, m in enumerate(masses):
        tp = masstransfer.track_point(cfg, k, t_yr) if mt_on else None
        if tp is None:
            out.append(star_properties(c, f"m{k + 1}", m))
            continue
        # a star re-modelled after a mass transfer: the structure of the track it follows
        m_ev, age = tp
        ck = copy(cfg)
        ck.stellar_model = cfg.stellar_model.model_copy(
            update={"age_myr": age - age_offset(cfg, f"m{k + 1}")}
        )
        out.append(star_properties(ck, f"m{k + 1}", m_ev))
    if mt_on:
        for k, p in enumerate(out):  # radius and luminosity of the re-modelled stars
            s = masstransfer.star_state(cfg, k, t_yr)
            if s is not None:
                out[k] = p = replace(
                    p,
                    radius_rsun=s.radius_rsun,
                    luminosity_lsun=s.luminosity_lsun,
                    teff_k=s.teff_k,
                )
            # a donor in transfer is held at its Roche lobe
            r = masstransfer.donor_radius(cfg, k, t_yr)
            if r is not None and r < p.radius_rsun:
                out[k] = replace(p, radius_rsun=r)
    return out


def structure_series(cfg, masses, t_yr: np.ndarray) -> np.ndarray:
    """Radius, k2 and I/(MR^2) of all stars at the run times t_yr (shape 3 n_star x n_t)."""
    rows = []
    for t in np.asarray(t_yr, float):
        ps = props_at(cfg, masses, float(t))
        rows.append([v for p in ps for v in (p.radius_rsun, p.apsidal_k2, p.rg2)])
    return np.array(rows).T


def _t_end(tr) -> float:
    """End of the evolution a track covers: helium ignition (SSE), end of core helium burning
    or of the model ages (MIST)."""
    return float(getattr(tr, "t_end", None) or tr.t_hei)


def _wolf_rayet(tr, age_myr: float) -> bool:
    """Wolf-Rayet flag of the initial-mass track; a star whose mass transfer moved it to another
    track (a donor or gainer) can outlive its initial track, and is then not flagged."""
    try:
        return bool(tr.is_wolf_rayet(age_myr))
    except EvolvedBeyondGrid:
        return False


def time_limit(cfg, masses, t_max: float) -> tuple[float, str | None]:
    """Run time at which the evolution leaves the grid in the direction of t_max, or a mass
    transfer ends in an event swingersim does not follow (common envelope, contact, a donor
    stripped to its core)."""
    if not cfg.physics.evolution:
        return t_max, None
    if masstransfer.enabled(cfg) and t_max >= 0:
        st = masstransfer.stop(cfg)
        t_lim, why = _evolution_limit(cfg, masses, t_max)
        if st is not None and st[0] < t_lim:
            return st[0], "mass_transfer: " + st[1]
        return t_lim, why
    return _evolution_limit(cfg, masses, t_max)


def _evolution_limit(cfg, masses, t_max: float) -> tuple[float, str | None]:
    age0 = _age0(cfg)
    offs = [age_offset(cfg, f"m{k + 1}") for k in range(len(masses))]
    if t_max >= 0:
        ends = []
        for k, (m, o) in enumerate(zip(masses, offs, strict=True)):
            lim = masstransfer.track_limit(cfg, k) if masstransfer.enabled(cfg) else None
            if lim is None:
                lim = (_t_end(stellar_track(cfg.stellar_model, m)) - o - age0) * YR_PER_MYR
            ends.append(lim)
        t_end = min(ends) * (1.0 - 1e-9)
        return (t_max, None) if t_end >= t_max else (max(t_end, 0.0), "stellar_evolution_end")
    t_start = max(
        getattr(stellar_track(cfg.stellar_model, m), "t_start", 0.0) - o
        for m, o in zip(masses, offs)
    )
    t_end = -(age0 - t_start) * YR_PER_MYR * (1.0 - 1e-9)
    return (t_max, None) if t_end <= t_max else (t_end, "stellar_birth")


def segment_edges(cfg, masses, t_max: float, max_dr: float, extra_edges=None) -> np.ndarray:
    """Segment boundaries from 0 to t_max such that no radius, apsidal constant, moment of
    inertia or mass changes by more than max_dr (relative); the times extra_edges (within the
    run, e.g. encounters) are boundaries as well."""
    n = 4000
    sgn = -1.0 if t_max < 0 else 1.0

    def in_run_order(x):  # unique times ordered from 0 towards t_max (backward runs: descending)
        x = np.unique(x)
        return x[::-1] if sgn < 0 else x

    lo, hi = min(0.0, t_max), max(0.0, t_max)
    t = np.linspace(0.0, t_max, n)
    if masstransfer.enabled(cfg):  # resolve the transfer episodes
        tg = masstransfer.time_grid(cfg)
        t = in_run_order(np.concatenate([t, tg[(tg > lo) & (tg < hi)]]))
        n = len(t)
    if cfg.physics.evolution:
        # the nodes of the stellar tracks (MIST: equivalent evolutionary points), between which
        # the radius is interpolated linearly: every maximum of a thermal pulse is a node, so
        # with the nodes in the grid no pulse can fall between two grid points
        nodes = []
        for k, m in enumerate(masses):
            tr = stellar_track(cfg.stellar_model, m)
            ages = getattr(tr, "_age", None)
            if ages is None:
                continue
            tk = (np.asarray(ages) - age_offset(cfg, f"m{k + 1}") - _age0(cfg)) * YR_PER_MYR
            nodes.append(tk[(tk > lo) & (tk < hi)])
        if nodes:
            t = in_run_order(np.concatenate([t, *nodes]))
            n = len(t)

    def series(tt):
        parts = []
        if cfg.physics.evolution:
            parts.append(np.log(structure_series(cfg, masses, tt)))
        if massloss.enabled(cfg):
            parts.append(np.log(massloss.mass_series(cfg, masses, tt)))
        return np.vstack(parts)

    lr = series(t)
    # refine the grid where a quantity changes by more than max_dr between two grid points
    # (fast phases: the tip of the RGB, thermal pulses on the AGB, which the uniform grid of
    # t_max/4000 can step over - a giant of 500 Rsun was seen at 313 Rsun before 0.27)
    min_dt = max(abs(t_max) * 1e-9, 1e-3)
    for _ in range(14):
        jump = np.max(np.abs(np.diff(lr, axis=1)), axis=0)
        bad = np.nonzero((jump > max_dr) & (np.abs(np.diff(t)) > 2 * min_dt))[0]
        if bad.size == 0:
            break
        mids = 0.5 * (t[bad] + t[bad + 1])
        lm = series(mids)
        t = np.concatenate([t, mids])
        order = np.argsort(sgn * t, kind="stable")
        t = t[order]
        lr = np.hstack([lr, lm])[:, order]
    n = len(t)
    edges = [0]
    start = 0
    for j in range(1, n):
        if np.max(np.abs(lr[:, j] - lr[:, start])) > max_dr:
            edges.append(j)
            start = j
    if edges[-1] != n - 1:
        edges.append(n - 1)
    out = t[edges]
    if extra_edges is not None and len(extra_edges):
        x = np.asarray(extra_edges, float)
        out = in_run_order(np.concatenate([out, x[(x > lo) & (x < hi)]]))
    mb = getattr(cfg.physics, "magnetic_braking", None)
    if mb is not None:  # braking: segments no longer than max_step_yr
        fine = [out[0]]
        for a_, b_ in pairwise(out):
            k = max(1, int(np.ceil(abs(b_ - a_) / mb.max_step_yr)))
            fine.extend(np.linspace(a_, b_, k + 1)[1:])
        out = np.array(fine)
    return out


@dataclass
class Segment:
    t0: float
    t1: float
    sol: object  # dense solution of solve_ivp
    props: list[StarProperties]
    state: dict  # backend-specific structure snapshot


class SegmentedSolution:
    """Concatenation of dense solutions of consecutive segments (same interface as the
    solve_ivp result used by the backends: .sol(t), .t, .status, .message)."""

    def __init__(self, segments: list[Segment]):
        self.segments = segments
        self.status = min((s.sol.status for s in segments), default=0)
        self.message = "; ".join(s.sol.message for s in segments if s.sol.status < 0) or "ok"
        self.t = np.concatenate([s.sol.t for s in segments])
        self._forward = segments[-1].t1 >= segments[0].t0

    def index(self, t) -> np.ndarray:
        t = np.atleast_1d(np.asarray(t, float))
        if self._forward:
            starts = np.array([s.t0 for s in self.segments])
            idx = np.searchsorted(starts, t, side="right") - 1
        else:
            starts = np.array([-s.t0 for s in self.segments])
            idx = np.searchsorted(starts, -t, side="right") - 1
        return np.clip(idx, 0, len(self.segments) - 1)

    def sol(self, t):
        scalar = np.ndim(t) == 0
        tt = np.atleast_1d(np.asarray(t, float))
        idx = self.index(tt)
        first = self.segments[0].sol.sol(tt[:1])
        out = np.empty((first.shape[0], len(tt)))
        for k in np.unique(idx):
            sel = idx == k
            out[:, sel] = self.segments[k].sol.sol(tt[sel])
        return out[:, 0] if scalar else out


SEGMENT = {"dt": 0.0}  # length [yr] of the operator-splitting segment being installed


def brake_spin(cfg, k: int, omega: float, props: StarProperties, mass: float, dt_yr: float):
    """Spin rate [rad/yr] of star k after dt_yr of magnetic braking (physics.magnetic_braking,
    Matt et al. 2015), with the structure props in force; in a backward run the spin dt_yr
    earlier (faster)."""
    spec = getattr(cfg.physics, "magnetic_braking", None)
    if spec is None or dt_yr == 0 or props.teff_k is None or omega == 0.0:
        return omega
    if spec.stars is not None and f"m{k + 1}" not in spec.stars:
        return omega
    from .. import magnetic

    w = magnetic.spin_down(
        omega / magnetic.SEC_PER_YR * 1.0,
        mass,
        props.radius_rsun,
        props.teff_k,
        props.rg2,
        -abs(dt_yr) if cfg.run.t_max_yr < 0 else abs(dt_yr),  # a backward run: spin earlier
        spec.chi,
        spec.p,
        spec.t0_erg,
    )
    return w * magnetic.SEC_PER_YR


def scale_spin(
    omega: float,
    old: StarProperties,
    new: StarProperties,
    m_old: float | None = None,
    m_new: float | None = None,
) -> float:
    """Spin rate after a change of the moment of inertia at constant spin angular momentum.
    A star losing mass (wind, mass transfer) loses the spin angular momentum of
    massloss.wind_spin_factor (Hurley et al. 2002, eq. 11); an accreting star gains the
    specific angular momentum of a Keplerian orbit at its surface, dJ = dM sqrt(G M R), and
    spins at most at the break-up rate sqrt(G M / R^3) (Hurley et al. 2002, Sect. 2.6.4: the
    accretor spun up beyond it would shed the excess)."""
    i_old = old.rg2 * old.radius_rsun**2
    i_new = new.rg2 * new.radius_rsun**2
    if i_new <= 0:
        return omega
    if m_old is None or m_new is None or m_new == m_old:
        return omega * i_old / i_new
    if m_new < m_old:
        rg2 = 0.5 * (old.rg2 + new.rg2)
        f = massloss.wind_spin_factor(m_old, m_new, rg2) * m_old / m_new
        w = omega * f * i_old / i_new
    else:
        r_new = new.radius_rsun * units.RSUN_AU
        j = m_old * i_old * units.RSUN_AU**2 * omega
        j += (m_new - m_old) * np.sqrt(units.G * m_new * r_new) * np.sign(omega)
        w = j / (m_new * i_new * units.RSUN_AU**2)
    w_crit = np.sqrt(units.G * m_new / (new.radius_rsun * units.RSUN_AU) ** 3)
    return float(np.clip(w, -w_crit, w_crit))


class Stream:
    """Result of a streamed integration: states at the outputs and extremes collected chunk by
    chunk (the dense solution of a Gyr run does not fit in memory)."""

    def __init__(self, Y, extremes, roche, h_ref, snapshots, n_steps, t_stop=None):
        self.t_stop = t_stop  # run time of a Roche-lobe overflow that ended the run, or None
        self.Y = Y
        self.extremes = extremes
        self.roche = roche
        self.h_ref = h_ref  # energy reference of every output (start of its segment)
        self.snapshots = snapshots  # stellar structure in force at every output
        self.n_steps = n_steps
        self.stop_stars: list = []  # the stars whose overflow ended the run


@dataclass
class OrbitSpec:
    name: str
    e_index: int  # first index of the Laplace vector (3 components) in the state
    a_index: int | None  # index of a in the state, or None (a0, or a_of() when given)
    a0: float
    stars: tuple[int, int] | None  # the two stars of a two-star binary (Roche lobes)
    # (star, companions) for every star with a Roche lobe in this orbit: both stars of a
    # two-star binary, or a star orbiting a subsystem (a tertiary around an inner binary, whose
    # lobe is taken against the total mass of the subsystem; Eggleton 1983)
    lobes: tuple[tuple[int, tuple[int, ...]], ...] | None = None
    a_of: object = None  # callable returning the semi-major axis in force (masses changing)

    def __post_init__(self):
        if self.lobes is None and self.stars is not None:
            i, j = self.stars
            self.lobes = ((i, (j,)), (j, (i,)))


def stream_evolving(
    cfg,
    masses,
    props,
    times,
    y0,
    rhs,
    energy,
    apply_props,
    snapshot,
    orbits: list[OrbitSpec],
    thresholds,
    threshold_key,
    refine,
    roche_lobe,
    angmom=None,
    scale_orbits=None,
    extra_edges=None,
) -> Stream:
    """Integrate with evolving stars: operator-splitting segments (``segment_edges``) cut into
    chunks; each chunk's dense solution is reduced to running extremes and discarded.

    apply_props(pm, prev, y, mm, m_prev) installs the structure pm and the stellar masses mm
    (scaling spins and, with mass loss, the semi-major axes in y or in the system; prev and
    m_prev are those of the previous segment), snapshot() returns the structure in force,
    energy(y) the conserved energy of the segment and angmom(y) (optional) its angular
    momentum. Without physics.evolution the structure stays ``props``; without mass loss the
    masses stay ``masses``. scale_orbits(y, factors) multiplies the semi-major axes of the
    named orbits (mass transfer: the change of a beyond a M = const). extra_edges are further
    segment boundaries (encounters: applied at their own time)."""
    from scipy.integrate import solve_ivp

    t_end = float(times[-1])
    edges = segment_edges(cfg, masses, t_end, cfg.run.evolution_max_dr, extra_edges)
    sign = 1.0 if t_end >= 0 else -1.0
    s_times = sign * np.asarray(times, float)  # increasing
    rtol = cfg.run.secular_rtol
    y = np.array(y0, float)
    n_out = len(times)
    Y = np.empty((len(y), n_out))
    Y[:, 0] = y
    h_ref = np.empty(n_out)
    l_ref = np.zeros((n_out, 3))
    snaps = [None] * n_out
    evolving = bool(cfg.physics.evolution)
    ext = {}
    for o in orbits:
        ext[o.name] = {"e_max": -1.0, "t_e_max_yr": 0.0, "q_min_au": np.inf, "t_q_min_yr": 0.0}
        for thr in thresholds:
            ext[o.name][threshold_key(thr)] = float("nan")
    rl, roche = {}, {}
    e_fill = {}  # star -> eccentricity of its orbit at the largest Roche-lobe filling

    def lobes(mm):
        for o in orbits:
            for s, others in o.lobes or ():
                rl[s] = roche_lobe(mm[s] / sum(mm[k] for k in others))

    m_prev = massloss.initial_masses(cfg, masses)
    lobes(m_prev)
    for s in rl:
        roche[f"m{s + 1}"] = {"fill_max": 0.0, "t_fill_max_yr": 0.0}
    n_steps, j_out, prev, h0, L0 = 0, 1, props, 0.0, np.zeros(3)
    # chunk length: about 10^4 integrator steps per chunk (bounded memory of the dense
    # solution); the first chunk probes the step size
    h_chunk = (times[1] - times[0]) / 1000.0 if n_out > 1 else t_end / 1000.0
    mt_on = masstransfer.enabled(cfg)
    mt_stars = masstransfer.transferring_stars(cfg) if mt_on else set()
    t_prev_mid = 0.0
    for ta, tb in pairwise(edges):
        tm = 0.5 * (ta + tb)
        pm = props_at(cfg, masses, tm) if evolving else list(props)
        mm = massloss.masses_at(cfg, masses, tm)
        SEGMENT.update(dt=abs(tb - ta), t0=ta, t1=tb)
        apply_props(pm, prev, y, mm, m_prev)
        if mt_on and scale_orbits is not None:
            f = masstransfer.a_factors(cfg, t_prev_mid, tm)
            scale_orbits(y, {n: v for n, v in f.items() if v != 1.0})
        t_prev_mid = tm
        for o in orbits:  # semi-major axes in force in this segment (changing masses)
            if o.a_of is not None:
                o.a0 = float(o.a_of())
        lobes(mm)
        h0 = energy(y)
        L0 = np.asarray(angmom(y)) if angmom is not None else L0
        if ta == 0.0:
            for o in orbits:
                if not o.lobes:
                    continue
                a_now = y[o.a_index] if o.a_index is not None else o.a0
                q_now = a_now * (1.0 - np.linalg.norm(y[o.e_index : o.e_index + 3]))
                for s_, _ in o.lobes:
                    if mt_on and s_ in mt_stars and masstransfer.roche_excluded(cfg, s_, 0.0):
                        continue  # a transfer from t = 0 (e.g. physics.mass_transfer.dynamic_onset)
                    if pm[s_].radius_rsun * units.RSUN_AU >= rl[s_] * q_now:
                        raise ValueError(
                            f"star m{s_ + 1} ({pm[s_].phase}, R = {pm[s_].radius_rsun:.2f} Rsun) "
                            f"overfills its Roche lobe at t = 0 in orbit {o.name}"
                        )
        if snaps[0] is None:
            snaps[0], h_ref[0], l_ref[0] = snapshot(), h0, L0
        radii = {s: pm[s].radius_rsun * units.RSUN_AU for s in rl}
        if evolving and rl:
            try:
                pa, pb = props_at(cfg, masses, ta), props_at(cfg, masses, tb)
                radii = _SegmentRadii(
                    {
                        s: (pa[s].radius_rsun * units.RSUN_AU, radii[s], pb[s].radius_rsun * units.RSUN_AU)
                        for s in rl
                    },
                    ta, tm, tb,
                )  # fmt: skip
            except (ValueError, RuntimeError):  # an edge beyond the end of a track
                pass
        k_stop = j_out + int(np.searchsorted(s_times[j_out:], sign * tb, side="left"))
        stops = [*times[j_out:k_stop], tb]
        t_now = ta
        for stop in stops:
            while sign * (stop - t_now) > 1e-12 * max(abs(stop), 1.0):
                full = sign * (t_now + h_chunk - stop) < 0
                t_next = t_now + h_chunk if full else stop
                sol = solve_ivp(
                    rhs,
                    (t_now, t_next),
                    y,
                    method="DOP853",
                    rtol=rtol,
                    atol=rtol * 1e-2,
                    dense_output=True,
                )
                if sol.status < 0:
                    raise RuntimeError(f"secular integration failed: {sol.message}")
                n_steps += len(sol.t) - 1
                first: dict = {}  # star -> first contact with its lobe in this chunk
                saved = (deepcopy(ext), deepcopy(roche), dict(e_fill))
                _chunk_extremes(
                    sol, orbits, ext, roche, radii, rl, thresholds, threshold_key, refine, e_fill,
                    first,
                )  # fmt: skip
                y = sol.y[:, -1].copy()
                t_now = t_next
                over_by = {
                    k_: first.get(int(k_[1:]) - 1, r["t_fill_max_yr"])
                    for k_, r in roche.items()
                    if r["fill_max"] >= 1.0
                    and not (
                        int(k_[1:]) - 1 in mt_stars
                        and masstransfer.roche_excluded(
                            cfg,
                            int(k_[1:]) - 1,
                            r["t_fill_max_yr"],
                            e_fill.get(int(k_[1:]) - 1, 0.0),
                        )
                    )
                }
                over = list(over_by.values())
                if over:  # mass transfer is not modelled: the run ends at the overflow
                    t_stop = min(over, key=lambda x: sign * x)
                    stop_stars = sorted(k_ for k_, t in over_by.items() if t == t_stop)
                    # the state at the first contact with the lobe within this chunk, and the
                    # extremes up to that moment only
                    if min(sol.t[0], sol.t[-1]) <= t_stop <= max(sol.t[0], sol.t[-1]):
                        y = sol.sol(t_stop).copy()
                        ext, roche, e_fill = saved
                        _chunk_extremes(
                            _Until(sol, t_stop), orbits, ext, roche, radii, rl, thresholds,
                            threshold_key, refine, e_fill,
                        )  # fmt: skip
                    while j_out < n_out:
                        Y[:, j_out], h_ref[j_out], snaps[j_out] = y, h0, snapshot()
                        l_ref[j_out] = L0
                        j_out += 1
                    for k_ in stop_stars:  # at its lobe at t_stop (to the bisection tolerance)
                        if roche[k_]["fill_max"] < 1.0:
                            roche[k_]["fill_max"], roche[k_]["t_fill_max_yr"] = 1.0, float(t_stop)
                    out = Stream(Y, ext, roche, h_ref, snaps, n_steps, t_stop)
                    out.l_ref = l_ref
                    out.stop_stars = stop_stars
                    return out
                dt = abs(sol.t[-1] - sol.t[0])
                if dt > 0 and (full or len(sol.t) > 10000):
                    h_chunk = sign * dt * min(10000.0 / max(len(sol.t) - 1, 1), 100.0)
            while j_out < n_out and sign * times[j_out] <= sign * stop * (1 + 1e-12) + 1e-9:
                Y[:, j_out], h_ref[j_out], snaps[j_out] = y, h0, snapshot()
                l_ref[j_out] = L0
                j_out += 1
        prev, m_prev = pm, mm
    while j_out < n_out:
        Y[:, j_out], h_ref[j_out], snaps[j_out] = y, h0, snapshot()
        l_ref[j_out] = L0
        j_out += 1
    out = Stream(Y, ext, roche, h_ref, snaps, n_steps)
    out.l_ref = l_ref
    return out


class _Until:
    """The dense solution of a chunk cut at time t_end (its steps up to t_end, then t_end)."""

    def __init__(self, sol, t_end: float):
        fwd = sol.t[-1] >= sol.t[0]
        keep = sol.t[sol.t < t_end] if fwd else sol.t[sol.t > t_end]
        self.t = np.append(keep, t_end)
        self.sol = sol.sol


class _SegmentRadii:
    """Radii [AU] of the stars with a Roche lobe through one segment: log-linear between the
    start, the middle and the end of the segment (the dynamics uses the structure of the
    middle; the overflow check follows the radius, so that the moment of contact does not
    depend on run.evolution_max_dr)."""

    def __init__(self, nodes: dict, ta: float, tm: float, tb: float):
        self.nodes = nodes  # star -> (R(ta), R(tm), R(tb))
        self.t = np.array(sorted((ta, tm, tb)))
        self.rev = ta > tb

    def at(self, s, t):
        ra, rm, rb = self.nodes[s]
        r = (rb, rm, ra) if self.rev else (ra, rm, rb)
        return np.exp(np.interp(t, self.t, np.log(r)))

    def __getitem__(self, s):  # the radius of the middle (as before)
        return self.nodes[s][1]


def _radius(radii, s, t):
    return radii.at(s, t) if hasattr(radii, "at") else radii[s]


def _first_fill(sol, fine, fill, o, sl, r_over_rl, peak: float | None = None) -> float:
    """First time of a chunk (in the direction of integration) at which the star fills its lobe
    at periastron, a (1 - e) = R / r_L: bisection on the dense interpolant. peak: the time of a
    refined maximum of the filling >= 1 between the samples (no sample fills the lobe)."""
    forward = sol.t[-1] >= sol.t[0]
    if peak is None:
        idx = np.flatnonzero(fill >= 1.0)
        k = int(idx[0] if forward else idx[-1])
        j = k - 1 if forward else k + 1
        if j < 0 or j >= len(fine):
            return float(fine[k])
        t_in, t_out = float(fine[k]), float(fine[j])
    else:  # the sample before the peak in the direction of the run
        before = fine[fine < peak] if forward else fine[fine > peak]
        if not len(before):
            return float(peak)
        t_in, t_out = float(peak), float(before[-1] if forward else before[0])

    def f(x):
        y = sol.sol(x)
        a = y[o.a_index] if o.a_index is not None else o.a0
        rr = r_over_rl(x) if callable(r_over_rl) else r_over_rl
        return rr / (a * (1.0 - np.linalg.norm(y[sl : sl + 3]))) - 1.0

    lo, hi = t_out, t_in  # f(lo) < 0 <= f(hi)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if f(mid) >= 0.0:
            hi = mid
        else:
            lo = mid
        if abs(hi - lo) <= 1e-12 * max(abs(hi), 1.0):
            break
    return hi


def _chunk_extremes(
    sol, orbits, ext, roche, radii, rl, thresholds, threshold_key, refine, e_fill=None, first=None
):
    """Running extremes from one chunk: every integrator step plus 50 uniform points; the
    largest eccentricity is refined on the dense interpolant. e_fill (optional) receives the
    eccentricity at the largest Roche-lobe filling of every star; first (optional) the first time
    each star fills its lobe."""
    t = sol.t
    if len(t) < 2:
        return
    fine = np.union1d(t, np.linspace(min(t[0], t[-1]), max(t[0], t[-1]), 50))
    Yf = sol.sol(fine)
    for o in orbits:
        sl = o.e_index
        e_f = np.linalg.norm(Yf[sl : sl + 3], axis=0)
        a_f = Yf[o.a_index] if o.a_index is not None else np.full(len(fine), o.a0)
        q_f = a_f * (1.0 - e_f)
        e = ext[o.name]
        peak_e = None
        if e_f.max() > e["e_max"]:

            def ecc(x, sl=sl):
                return float(np.linalg.norm(sol.sol(x)[sl : sl + 3]))

            t_e, g = refine(fine, e_f, lambda x: -ecc(x))
            e["e_max"], e["t_e_max_yr"] = -g, t_e
            peak_e = (-g, t_e)
        k = int(np.argmin(q_f))
        t_q, q_ref = float(fine[k]), float(q_f[k])
        near_lobe = (
            bool(o.lobes)
            and max(float(np.max(_radius(radii, s, fine) / (rl[s] * q_f))) for s, _ in o.lobes)
            > 0.8
        )
        if q_ref < 1.05 * e["q_min_au"] or near_lobe:
            # a periastron passage closer than the samples (a short eccentricity spike)

            def peri(x, sl=sl, o=o):
                y = sol.sol(x)
                a = y[o.a_index] if o.a_index is not None else o.a0
                return float(a * (1.0 - np.linalg.norm(y[sl : sl + 3])))

            t_r, q_r = refine(fine, -q_f, peri)
            if q_r < q_ref:
                t_q, q_ref = float(t_r), float(q_r)
        if q_ref < e["q_min_au"]:
            e["q_min_au"], e["t_q_min_yr"] = q_ref, t_q
        for thr in thresholds:
            key = threshold_key(thr)
            if np.isnan(e[key]) and e_f.max() > 1.0 - thr:
                above = np.flatnonzero(e_f > 1.0 - thr)
                e[key] = float(fine[above[0] if t[-1] >= t[0] else above[-1]])
            elif np.isnan(e[key]) and peak_e is not None and peak_e[0] > 1.0 - thr:
                e[key] = float(peak_e[1])  # a peak shorter than the sampling
        if o.lobes:
            for s, _ in o.lobes:
                fill = _radius(radii, s, fine) / (rl[s] * q_f)
                kk = int(np.argmax(fill))
                best, t_b, e_b = float(fill[kk]), float(fine[kk]), float(e_f[kk])
                fill_r = float(_radius(radii, s, t_q)) / (rl[s] * q_ref)
                if fill_r > best:  # the refined periastron fills the lobe more than any sample
                    best, t_b = float(fill_r), t_q
                    e_b = float(np.linalg.norm(sol.sol(t_q)[sl : sl + 3]))
                r = roche[f"m{s + 1}"]
                if best > r["fill_max"]:
                    r["fill_max"], r["t_fill_max_yr"] = best, t_b
                    if e_fill is not None:
                        e_fill[s] = e_b
                if first is not None and best >= 1.0 and s not in first:
                    first[s] = _first_fill(
                        sol, fine, fill, o, sl,
                        lambda x, s=s: float(_radius(radii, s, x)) / rl[s],
                        peak=t_q if fill[kk] < 1.0 else None,
                    )  # fmt: skip


__all__ = [
    "OrbitSpec",
    "Segment",
    "SegmentedSolution",
    "Stream",
    "brake_spin",
    "props_at",
    "replace",
    "scale_spin",
    "segment_edges",
    "stream_evolving",
    "time_limit",
    "time_varying",
]


def _mist_series(tr, ages: np.ndarray) -> dict[str, np.ndarray] | None:
    """Radius, luminosity, Teff and phase of a MIST track at many ages at once (the same
    interpolation as MistTrack.at); None when an age lies outside the track (the caller then
    evaluates point by point and raises as before)."""
    from ..evolution import TEFF_SUN
    from ..mist import EEP_RGB_TIP, EEP_TACHEB, EEP_TAMS, EEP_TPAGB, EEP_ZAMS

    if not len(ages) or ages.min() < tr.t_start or ages.max() > tr.t_end:
        return None
    if "surface_h1" not in tr._cols:
        return None
    x = np.arange(1, tr.n_eep + 1)
    eep = np.clip(np.interp(ages, tr._age, x), 1.0, tr.n_eep)
    r = 10.0 ** np.interp(eep, x, tr._t[:, tr._cols["log_R"]])
    lum = 10.0 ** np.interp(eep, x, tr._t[:, tr._cols["log_L"]])
    bounds = [EEP_ZAMS, EEP_TAMS, tr.eep_gb, EEP_RGB_TIP, EEP_TACHEB, EEP_TPAGB]  # mist._phase
    teff = TEFF_SUN * lum**0.25 / np.sqrt(r)
    h1 = np.interp(eep, x, tr._t[:, tr._cols["surface_h1"]])
    return {
        "radius_rsun": r,
        "luminosity_lsun": lum,
        "teff_k": teff,
        "phase": np.searchsorted(bounds, eep, side="right").astype(np.int8),
        "age_myr": ages,
        "wolf_rayet": (h1 < 0.4) & (teff > 1.0e4),  # MistTrack.is_wolf_rayet
    }


def stellar_series(cfg, masses, times) -> dict[str, dict[str, np.ndarray]]:
    """Radius, luminosity, Teff and phase (1 MS, 2 HG, 3 GB) of every star at the outputs."""
    out = {}
    times = np.asarray(times, float)
    ages = _age0(cfg) + times / YR_PER_MYR
    mt_on = masstransfer.enabled(cfg)
    for k, m in enumerate(masses):
        tr = stellar_track(cfg.stellar_model, m)
        off = age_offset(cfg, f"m{k + 1}")
        # re-modelled states where the star follows a mass-transfer history
        remodel = [None] * len(times)
        if mt_on and masstransfer.track_point(cfg, k, float(np.max(np.abs(times)))) is not None:
            remodel = [masstransfer.star_state(cfg, k, float(t)) for t in times]
        if cfg.stellar_model.kind == "mist" and all(s is None for s in remodel):
            fast = _mist_series(tr, ages + off)
            if fast is not None:
                out[f"m{k + 1}"] = fast
                continue
        sts = [s if s is not None else tr.at(a + off) for s, a in zip(remodel, ages, strict=True)]
        out[f"m{k + 1}"] = {
            "radius_rsun": np.array([s.radius_rsun for s in sts]),
            "luminosity_lsun": np.array([s.luminosity_lsun for s in sts]),
            "teff_k": np.array([s.teff_k for s in sts]),
            "phase": np.array([s.phase for s in sts], dtype=np.int8),
            "age_myr": ages + off,
        }
        if hasattr(tr, "is_wolf_rayet"):  # MIST: surface H < 0.4 and Teff > 1e4 K
            out[f"m{k + 1}"]["wolf_rayet"] = np.array([_wolf_rayet(tr, a + off) for a in ages])
    return out
