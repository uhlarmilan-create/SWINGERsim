"""Secular period change of a two-star binary (orbits.<pair>.Pdot).

Semi-detached and contact binaries change their period by mass transfer and by the loss of
angular momentum (magnetised wind, the outer Lagrangian point); the eclipse times then lie on
a parabola, O - C = (1/2) (dP/dt) E^2 P (e.g. Kwee 1958; Hilditch 2001, Sect. 4.7). SWINGERsim
models this as an imposed, constant dP/dt of the binary, P(t) = P0 + Pdot t (t from t = 0 of the
configuration), inside the N-body integration, so that every observable (minima, radial
velocities, light curves, inclinations) and the third-body dynamics see the same orbit:

* the time axis is cut into intervals [(k - 1/2) dt, (k + 1/2) dt) around the grid times k dt;
  in interval k the binary has the period P(k dt) (interval 0 holds t = 0 with the configured
  P0). Where the integration crosses a boundary, in either direction, the relative orbit of the
  pair is scaled to the period of the interval it enters, a -> a (P_new / P_old)^(2/3), at fixed
  eccentricity vector and phase (the centre of mass of the pair and its velocity unchanged).
  The eclipse times then deviate from the parabola by at most (dP/dt / P) dt^2 / 8, and dt is
  chosen so that this stays below KICK_TOL_S;
* mode angular_momentum: the masses stay (loss of angular momentum, e.g. magnetic braking);
  mode conservative: conservative mass transfer between the two stars, dP/P = 3 dM1 (M1 - M2) /
  (M1 M2) (total mass and orbital angular momentum constant), so the masses change linearly with
  dM1/dt = (dP/dt / P) M1 M2 / (3 (M1 - M2)) and dM2/dt = -dM1/dt. The relative orbit depends
  on the total mass only, so the masses are changed in place and the separation scaled.

N-body backend only (the secular backends average over the orbits whose period is changed).
"""

from __future__ import annotations

import math

import numpy as np

from . import units

KICK_TOL_S = 0.05  # largest deviation of an eclipse time from the parabola [s]
KICK_MAX_YR = 0.25  # longest interval between two kicks [yr]


class PeriodDrift:
    """Kicks of the binaries with orbits.<pair>.Pdot of a simulation (one per simulation)."""

    def __init__(self, items: list[dict]):
        self.items = items
        self.changes_masses = any(it["mode"] == "conservative" for it in items)
        self.cur = {id(it): 0 for it in items}  # interval of every binary (0 holds t = 0)

    @classmethod
    def of(cls, cfg, root, masses) -> PeriodDrift | None:
        """The drifting binaries of the configuration (None without any Pdot)."""
        items = []
        for node in root.nodes():
            spec = cfg.orbits[node.name]
            pdot = float(getattr(spec, "Pdot", 0.0) or 0.0)
            if pdot == 0.0:
                continue
            i, j = node.left, node.right
            if not (isinstance(i, int) and isinstance(j, int)):
                raise ValueError(  # noqa: TRY004 - a configuration error
                    f"orbits.{node.name}.Pdot: only an orbit of two stars can drift (this orbit "
                    "contains a binary)"
                )
            m_i, m_j = float(masses[i]), float(masses[j])
            p0 = units.sma_to_period(node.elements.a, m_i + m_j)  # [yr]
            rate = pdot / units.DAYS_PER_YEAR / p0  # (dP/dt) / P [1/yr]
            mode = getattr(spec, "Pdot_mode", "angular_momentum")
            mdot = 0.0
            if mode == "conservative":
                if abs(m_i - m_j) < 1e-9 * (m_i + m_j):
                    raise ValueError(
                        f"orbits.{node.name}.Pdot_mode conservative: equal masses do not change "
                        "the period by a transfer (use angular_momentum)"
                    )
                mdot = rate * m_i * m_j / (3.0 * (m_i - m_j))  # dM_i/dt [Msun/yr]
            tol = KICK_TOL_S / units.SECONDS_PER_YEAR
            dt = min(max(math.sqrt(8.0 * tol / abs(rate)), p0), KICK_MAX_YR)
            items.append(
                {
                    "name": node.name,
                    "i": i,
                    "j": j,
                    "p0": p0,
                    "pdot": pdot / units.DAYS_PER_YEAR,  # [yr / yr]
                    "mode": mode,
                    "m0": (m_i, m_j),
                    "mdot": mdot,
                    "dt": dt,
                }
            )
        return cls(items) if items else None

    def period(self, it: dict, t: float) -> float:
        return it["p0"] + it["pdot"] * t

    def _kick(self, sim, it: dict, k_new: int) -> None:
        """Scale the binary from the period of its current interval to that of interval k_new."""
        k_old = self.cur[id(it)]
        if k_new == k_old:
            return
        p_prev, p_now = self.period(it, k_old * it["dt"]), self.period(it, k_new * it["dt"])
        if p_now <= 0.0 or p_prev <= 0.0:
            raise ValueError(
                f"orbits.{it['name']}.Pdot: the period reaches zero at t = {k_new * it['dt']:.4g} yr"
            )
        i, j = it["i"], it["j"]
        pi, pj = sim.particles[i], sim.particles[j]
        mi, mj = pi.m, pj.m
        M = mi + mj
        ri, rj = np.array(pi.xyz), np.array(pj.xyz)
        vi, vj = np.array(pi.vxyz), np.array(pj.vxyz)
        rc, vc = (mi * ri + mj * rj) / M, (mi * vi + mj * vj) / M
        g = (p_now / p_prev) ** (2.0 / 3.0)
        r, v = (rj - ri) * g, (vj - vi) / math.sqrt(g)
        if it["mode"] == "conservative":
            t = k_new * it["dt"]
            mi = it["m0"][0] + it["mdot"] * t
            mj = it["m0"][1] - it["mdot"] * t
            if mi <= 0.0 or mj <= 0.0:
                raise ValueError(f"orbits.{it['name']}.Pdot: a star loses all its mass")
            pi.m, pj.m = mi, mj
        pi.x, pi.y, pi.z = rc - mj / M * r
        pj.x, pj.y, pj.z = rc + mi / M * r
        pi.vx, pi.vy, pi.vz = vc - mj / M * v
        pj.vx, pj.vy, pj.vz = vc + mi / M * v
        self.cur[id(it)] = k_new

    def integrate(self, sim, t: float, exact_finish_time: int = 1) -> None:
        """sim.integrate(t) with the kicks of every drifting binary at the interval boundaries
        it crosses (in either direction)."""
        t0 = float(sim.t)
        for it in self.items:  # a previous call stopped past a boundary (exact_finish_time = 0)
            x = t0 / it["dt"]
            if abs(x - self.cur[id(it)]) > 0.5 + 1e-9:  # not just at the boundary itself
                self._kick(sim, it, round(x))
        if t == t0:
            return
        s = 1 if t > t0 else -1
        kicks = []  # (boundary time, item, interval entered)
        for it in self.items:
            dt = it["dt"]
            # boundaries (m + 1/2) dt strictly between t0 and t
            lo, hi = min(t0, t), max(t0, t)
            m0 = math.ceil(lo / dt - 0.5)
            m1 = math.floor(hi / dt - 0.5)
            for m in range(m0, m1 + 1):
                tb = (m + 0.5) * dt
                if lo < tb < hi or (tb == hi and s > 0) or (tb == lo and s < 0):
                    kicks.append((tb, it, m + 1 if s > 0 else m))
        kicks.sort(key=lambda x: s * x[0])
        for tb, it, k_new in kicks:
            if s * (tb - sim.t) > 0.0:
                sim.integrate(tb, exact_finish_time=1)
            self._kick(sim, it, k_new)
        if s * (t - sim.t) > 0.0:
            sim.integrate(t, exact_finish_time=exact_finish_time)


def integrator(cfg, sim, root, masses):
    """Function t, exact -> integrate sim to t (with the kicks of orbits.<pair>.Pdot, if any)."""
    drift = PeriodDrift.of(cfg, root, masses)
    if drift is None:
        return lambda t, exact=1: sim.integrate(t, exact_finish_time=exact)
    return lambda t, exact=1: drift.integrate(sim, t, exact)


__all__ = ["KICK_TOL_S", "PeriodDrift", "integrator"]
