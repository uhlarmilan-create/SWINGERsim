"""Model observables of a configuration from one N-body integration per direction (swingersim
fit, stage 2): the times of the measured minima, the radial velocities and the inclinations of
eclipsing pairs.

Every data type needs the N-body state at its own times: a minimum near the conjunction of its
pair, a velocity at its emission time, an inclination at its epoch. Instead of one integration
per data type over the whole span of the data, all requests go into one event queue ordered by
the time to which the integration has to go next, and one simulation (per direction, forward
and backward from t = 0) serves them all in turn. Each request type keeps its own logic:

* minima (dynfit.model_times): the pair's conjunction is searched with the two-body motion from
  0.75 P before its expected time, the simulation is integrated to it and the time refined,
  t_obs = t + z_cm / c - as before, and in the same order;
* radial velocities: the emission time t_e with t_e + z(t_e)/c = t_obs is approached from
  before (a margin larger than any change of the light time) and found by Newton steps with the
  acceleration of the star until the remaining mismatch is below 1 ms; the velocity is then
  corrected to first order with the acceleration. A point whose emission time was passed (it
  should not happen) is computed in a second integration at its estimated emission time;
* inclinations: the orbit of the pair at the epoch (no light time).

Stopping the integration at the times of the other data changes the IAS15 steps; the results
agree with separate integrations to the round-off of the integrator.
"""

from __future__ import annotations

import heapq
import itertools

import numpy as np

from . import units
from .timing import roemer_in_binary

DPY = units.DAYS_PER_YEAR
KMS_PER_AUYR = units.AU_KM / units.SECONDS_PER_YEAR  # 1 AU/yr = 4.74 km/s
# radial velocities: the emission time is accepted when the light-time mismatch is below
# RV_FINE_YR; the velocity is then corrected to first order with the acceleration of the star.
# REBOUND's stored acceleration belongs to the start of the last IAS15 step (~1 % off), so the
# correction must stay tiny: 1 ms leaves < 1e-8 km/s even for a 100 km/s, 2-day orbit. A point
# whose emission time lies behind the state by more than RV_BEHIND_YR goes to the second pass.
RV_FINE_YR = 1e-3 / units.SECONDS_PER_YEAR
RV_BEHIND_YR = 1e-2 / units.SECONDS_PER_YEAR
RV_TOL_YR = 5.0 / units.SECONDS_PER_YEAR
RV_MAX_STEPS = 40


def _leaves(c) -> list[int]:
    return c.stars() if hasattr(c, "stars") else [c]


def star_indices(star, n: int) -> np.ndarray:
    """Index (0-based) of every star name m1 ... mn."""
    out = np.empty(len(star), int)
    for j, s in enumerate(star):
        try:
            k = int(s[1:]) - 1
        except ValueError:
            k = -1
        if not s.startswith("m") or not 0 <= k < n:
            raise ValueError(f"radial velocities: unknown star '{s}' (m1 ... m{n})")
        out[j] = k
    return out


def pair_stars(root, pair: str) -> tuple[int, int]:
    for node in root.nodes():
        if node.name == pair:
            if not (isinstance(node.left, int) and isinstance(node.right, int)):
                raise ValueError(f"fit.inclination: orbit {pair} is not a pair of two stars")
            return node.left, node.right
    raise ValueError(f"fit.inclination: no orbit named {pair}")


def system_reach_au(root) -> float:
    """Upper bound of the distance of any star from the centre of mass [AU]: the sum of the
    apocentre distances of all orbits of the hierarchy."""
    tot = 0.0
    for n in root.nodes():
        el = n.elements
        a, e = abs(float(el.a)), float(el.e)
        tot += a * (1.0 + e) if 0.0 <= e < 1.0 else 10.0 * a
    return tot


def eclipse_limit(G, pi, pj, r, v) -> float:
    """Inclination [deg, 0-90] below which the pair does not eclipse at the nearer conjunction:
    cos i = (R1 + R2) / r_conj, r_conj = a (1 - e^2) / (1 + e |sin omega|), omega measured from
    the node on the sky plane (z is the line of sight)."""
    r = np.asarray(r, float)
    v = np.asarray(v, float)
    mu = G * (pi.m + pj.m)
    rn = float(np.linalg.norm(r))
    h = np.cross(r, v)
    hn = float(np.linalg.norm(h))
    ev = np.cross(v, h) / mu - r / rn
    e = float(np.linalg.norm(ev))
    a = 1.0 / (2.0 / rn - float(np.dot(v, v)) / mu)
    node = np.cross([0.0, 0.0, 1.0], h)
    nn = float(np.linalg.norm(node))
    sin_w = 0.0
    if e > 1e-9 and nn > 1e-12:
        sin_w = float(np.clip(np.dot(np.cross(node, ev), h) / (nn * e * hn), -1.0, 1.0))
    r_conj = a * (1.0 - e * e) / (1.0 + e * abs(sin_w))
    return float(np.degrees(np.arccos(min(1.0, (pi.r + pj.r) / r_conj))))


def observables(
    cfg,
    bjd_t0: float,
    minima=None,
    rv=None,
    incs=(),
    separation: np.ndarray | None = None,
) -> dict:
    """Model observables of the configuration cfg from one integration per direction.

    minima: object with arrays pair, kind (0 p, 1 s), time [BJD] (dynfit.Observed) or None;
    rv: (star, time [BJD]) or None; incs: list of (pair, time [BJD], with_limit).
    Returns {"t": model BJD of the minima (NaN where no conjunction was found) or None,
    "rv": velocities [km/s, without gamma] or None, "inc": [inclination [deg, folded to 0-90]
    or array (2, n) with the eclipse limit]}. With ``separation`` (length of minima) the
    projected separation [AU] of the pair at each conjunction is stored there."""
    out_t = None
    if minima is not None:
        T_min = (np.asarray(minima.time, float) - bjd_t0) / DPY
        out_t = np.full(len(T_min), np.nan)
    else:
        T_min = np.zeros(0)
    out_rv = None
    if rv is not None:
        rv_star = np.asarray(rv[0]).astype(str)
        T_rv = (np.asarray(rv[1], float) - bjd_t0) / DPY
        out_rv = np.full(len(T_rv), np.nan)
    else:
        rv_star, T_rv = np.zeros(0, str), np.zeros(0)
    inc_req = []
    out_inc = []
    for pair, times, with_limit in incs:
        T = (np.asarray(times, float) - bjd_t0) / DPY
        inc_req.append((pair, T, bool(with_limit)))
        out_inc.append((np.full(len(T), np.nan), np.full(len(T), np.nan)))
    # the direction of an observation is that of its emission time, T - z/c (z at t = 0): an
    # observation just after t = 0 can have been emitted before it (and the reverse)
    E_min, E_rv = T_min, T_rv
    if len(T_min) or len(T_rv):
        E_min, E_rv = _emission_estimates(cfg, minima, T_min, rv_star, T_rv)
    for sign in (1.0, -1.0):

        def side(T, sign=sign):
            return np.nonzero(T >= 0.0 if sign > 0 else T < 0.0)[0]

        jm, jr = side(E_min), side(E_rv)
        ji = [side(T) for _, T, _ in inc_req]
        if not (len(jm) or len(jr) or any(len(j) for j in ji)):
            continue
        _direction(
            cfg, sign, bjd_t0, minima, T_min, jm, out_t, separation, rv_star, T_rv, jr, out_rv,
            inc_req, ji, out_inc,
        )  # fmt: skip
    incs_out = [
        np.vstack([i, lim]) if req[2] else i for (i, lim), req in zip(out_inc, inc_req, strict=True)
    ]
    return {"t": out_t, "rv": out_rv, "inc": incs_out}


def _emission_estimates(cfg, minima, T_min, rv_star, T_rv):
    """Emission times [yr] of the observations to first order: T - z(0)/c of the star (radial
    velocities) or of the centre of mass of the eclipsing pair (minima)."""
    from .backends.nbody import build_simulation

    sim, _, _, root, masses = build_simulation(cfg)
    z = np.array([sim.particles[k].z for k in range(len(masses))])
    c = units.C
    e_rv = T_rv
    if len(T_rv):
        e_rv = T_rv - z[star_indices(rv_star, len(masses))] / c
    e_min = T_min
    if len(T_min):
        delay = {}
        for n in root.nodes():
            ix = _leaves(n.left) + _leaves(n.right)
            delay[n.name] = float(np.sum(masses[ix] * z[ix]) / np.sum(masses[ix])) / c
        e_min = T_min - np.array([delay.get(str(p), 0.0) for p in np.asarray(minima.pair)])
    return e_min, e_rv


def _direction(
    cfg, sign, bjd_t0, minima, T_min, jm, out_t, separation, rv_star, T_rv, jr, out_rv,
    inc_req, ji, out_inc,
):  # fmt: skip
    from .backends.nbody import build_simulation
    from .dynfit import _conjunction
    from .pdot import integrator

    sim, _, _, root, masses = build_simulation(cfg)
    sim.collision = "none"
    if sign < 0:
        sim.dt = -abs(sim.dt)
    ps = sim.particles
    advance = integrator(cfg, sim, root, masses)  # with the kicks of orbits.<pair>.Pdot
    c = units.C
    heap: list = []
    seq = itertools.count()

    def push(t, fn, *args):
        heapq.heappush(heap, (sign * t, next(seq), fn, args))

    def go(t):
        if sign * (t - sim.t) > 0.0:
            advance(t)

    # ---------------------------------------------------------------- minima
    if len(jm):
        order = jm[np.argsort(sign * T_min[jm], kind="stable")]
        nodes = {n.name: n for n in root.nodes()}
        missing = sorted(set(np.asarray(minima.pair)[order].tolist()) - set(nodes))
        if missing:
            raise ValueError(f"fit.stage2: {missing} is not an orbit of the configuration")
        sides = {k: (_leaves(n.left), _leaves(n.right)) for k, n in nodes.items()}
        mass_of_side = {
            k: (float(np.sum(masses[a])), float(np.sum(masses[b]))) for k, (a, b) in sides.items()
        }
        period = {
            k: units.sma_to_period(n.elements.a, sum(mass_of_side[k])) for k, n in nodes.items()
        }
        activity = bool(getattr(cfg.physics, "activity", None))
        pair_of = np.asarray(minima.pair)
        kind_of = np.asarray(minima.kind)

        def state(n):
            a_idx, b_idx = sides[n.name]
            ma, mb = mass_of_side[n.name]

            def com(ix, m):
                x = np.array([[ps[k].x, ps[k].y, ps[k].z] for k in ix])
                v = np.array([[ps[k].vx, ps[k].vy, ps[k].vz] for k in ix])
                w = masses[ix] / m
                return w @ x, w @ v

            xa, va = com(a_idx, ma)
            xb, vb = com(b_idx, mb)
            zc = (ma * xa[2] + mb * xb[2]) / (ma + mb)
            vzc = (ma * va[2] + mb * vb[2]) / (ma + mb)
            return xb - xa, vb - va, zc, vzc

        # light time of every orbit's centre of mass at t = 0: the first minimum of a pair is
        # searched near its emission time (with a short pair in a wide outer orbit z/c can
        # exceed half the period, and a start at T would pick the neighbouring cycle)
        delay = {k: state(nn)[2] / c for k, nn in nodes.items()}
        queue = iter(order.tolist())

        def next_minimum():
            j = next(queue, None)
            if j is None:
                return
            n = nodes[pair_of[j]]
            P = period[n.name]
            target = T_min[j] - delay[n.name]  # dynamical time of the minimum (estimate)
            start = target - sign * 0.75 * P
            push(start, min_search, j, n, P, target, start)

        def min_search(j, n, P, target, start):
            go(start)
            r, v, zc, vzc = state(n)
            mu = units.G * sum(mass_of_side[n.name])
            rel = target - sim.t
            tau = _conjunction(
                mu, r, v, min(rel, 0.0) - 0.8 * P, max(rel, 0.0) + 0.8 * P, rel, kind_of[j]
            )
            if tau is None:
                next_minimum()
                return
            tc = sim.t + tau
            if sign * (tc - sim.t) > 0.0:  # integrate to the conjunction and refine it
                push(tc, min_refine, j, n, P, mu, tc)
            else:  # the state is already past it (closely spaced minima)
                min_done(j, n, tc, r + v * (tc - sim.t), zc + vzc * (tc - sim.t))

        def min_refine(j, n, P, mu, tc):
            go(tc)
            r, v, zc, vzc = state(n)
            fine = _conjunction(mu, r, v, -0.1 * P, 0.1 * P, 0.0, kind_of[j])
            if fine is not None:
                tc = sim.t + fine
                zc = zc + vzc * fine
                r = r + v * fine
            min_done(j, n, tc, r, zc)

        def min_done(j, n, tc, r, zc):
            delay[n.name] = zc / c
            ma, mb = mass_of_side[n.name]
            # the light time of the pair's centre of mass and across the pair (Kaplan 2010)
            t_obs = tc + zc / c + roemer_in_binary(float(r[2]), ma, mb)
            if activity and len(sides[n.name][0]) == 1 and len(sides[n.name][1]) == 1:
                from .magnetic import activity_shift_yr

                t_obs += float(activity_shift_yr(cfg, (n.left, n.right), tc))
            out_t[j] = bjd_t0 + t_obs * DPY
            if separation is not None:
                separation[j] = float(np.hypot(r[0], r[1]))
            next_minimum()

        next_minimum()

    # ---------------------------------------------------------------- radial velocities
    redo = []
    if len(jr):
        kk = star_indices(rv_star, len(masses))
        d0 = np.array([ps[k].z / c for k in range(len(masses))])  # light time at t = 0
        # the light time of a star changes by less than 2 x its largest distance from the
        # centre of mass / c: approaching from that far before, t_e is never passed
        margin = 2.4 * system_reach_au(root) / c + RV_TOL_YR

        def rv_step(j, it, target):
            go(target)
            p = ps[int(kk[j])]
            f = sim.t + p.z / c - T_rv[j]
            s1 = -f / (1.0 + p.vz / c)
            behind = sign * s1 < 0.0
            if abs(s1) <= RV_FINE_YR or (behind and abs(s1) <= RV_BEHIND_YR):
                out_rv[j] = (p.vz + p.az * s1) * KMS_PER_AUYR
                return
            if behind or it >= RV_MAX_STEPS:
                redo.append((j, sim.t + s1))
                return
            # Newton step with the acceleration; land short of the estimate by a guard (twice
            # the second-order term + 1e-3 of the step) until the guard is below RV_FINE_YR
            s2 = -(f + p.az * s1 * s1 / (2.0 * c)) / (1.0 + p.vz / c)
            s = s2 if sign * s2 > 0.0 else s1
            guard = 2.0 * abs(s2 - s1) + 1e-3 * abs(s)
            step = abs(s) if guard < RV_FINE_YR else max(abs(s) - guard, 0.5 * abs(s))
            nxt = sim.t + sign * step
            push(nxt, rv_step, j, it + 1, nxt)

        for j in jr:
            a = T_rv[j] - d0[kk[j]] - sign * margin
            push(a, rv_step, int(j), 0, a)

    # ---------------------------------------------------------------- inclinations
    for q, ((pair, T, with_limit), idx) in enumerate(zip(inc_req, ji, strict=True)):
        if not len(idx):
            continue
        i, j = pair_stars(root, pair)

        def inc_at(k, T=T, i=i, j=j, q=q, with_limit=with_limit):
            go(T[k])
            r = np.array(ps[j].xyz) - np.array(ps[i].xyz)
            v = np.array(ps[j].vxyz) - np.array(ps[i].vxyz)
            h = np.cross(r, v)
            out_inc[q][0][k] = np.degrees(np.arccos(min(1.0, abs(h[2]) / np.linalg.norm(h))))
            if with_limit:
                out_inc[q][1][k] = eclipse_limit(sim.G, ps[i], ps[j], r, v)

        for k in idx[np.argsort(sign * T[idx], kind="stable")]:
            push(T[k], inc_at, int(k))

    while heap:
        _, _, fn, args = heapq.heappop(heap)
        fn(*args)

    if redo:  # emission times passed (should not happen): a second integration, as before
        sim, _, _, _, _ = build_simulation(cfg)
        sim.collision = "none"
        if sign < 0:
            sim.dt = -abs(sim.dt)
        ps = sim.particles
        advance = integrator(cfg, sim, root, masses)
        for j, te in sorted(redo, key=lambda r: sign * r[1]):
            if sign * (te - sim.t) > 0.0:
                advance(te)
            p = ps[int(kk[j])]
            f = sim.t + p.z / c - T_rv[j]
            dt = -f / (1.0 + p.vz / c)
            out_rv[j] = (p.vz + p.az * dt) * KMS_PER_AUYR


__all__ = ["eclipse_limit", "observables", "pair_stars", "star_indices"]
