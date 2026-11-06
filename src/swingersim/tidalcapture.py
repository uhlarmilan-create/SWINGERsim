"""Energy dissipated by tides in a close stellar encounter (tidal capture).

During a close periastron passage each star is raised a tide that is left oscillating after
the passage; the energy taken from the orbit is (Press & Teukolsky 1977)

    E_tide = sum over the two stars of  G M_other^2 / R [ (R/R_p)^6 T_2(eta) + (R/R_p)^8 T_3(eta) ],
    eta = (M / (M + M_other))^(1/2) (R_p / R)^(3/2),

with the periastron distance R_p and the dimensionless functions T_l of the stellar structure,
here the fifth-order fits of Portegies Zwart & Meinen (1993, A&A 280, 174, Table 1 and eqs. 4,
7) to the polytropes n = 1.5, 2 and 3 of Lee & Ostriker (1986) and Ray et al. (1987):
log10 T_l = A + B x + C x^2 + D x^3 + E x^4 + F x^5 with x = log10 eta. The fits cover
0 <= x <= 1; below x = 0 the stars touch, above x = 1 the passage is slow and T_l falls steeply
(continued here as a power law with the slope of the fit at x = 1).

The oscillation energy is assumed to be dissipated before the next periastron (Portegies Zwart &
Meinen 1993, Sect. 4). swingersim scatter applies E_tide at every periastron passage of two stars
(physics of the scattering experiments, ``tidal_capture``) by reducing their relative velocity
at fixed centre of mass.
"""

from __future__ import annotations

import math

import numpy as np

from . import units

# Table 1 of Portegies Zwart & Meinen (1993): coefficients A..F of log10 T_l(log10 eta)
PZM93 = {
    1.5: {
        2: (-0.397, 1.678, 1.277, -12.42, 9.446, -5.550),
        3: (-0.909, 1.574, 12.37, -57.40, 80.10, -46.43),
    },
    2.0: {
        2: (-0.517, -0.906, 23.88, -93.49, 112.3, -44.15),
        3: (-1.040, -1.354, 37.64, -139.9, 168.2, -66.53),
    },
    3.0: {
        2: (-1.124, 0.877, -13.37, 21.55, -16.48, 4.124),
        3: (-1.703, 2.653, -14.34, 12.85, -0.492, -3.600),
    },
}


def _poly(c, x):
    return sum(ck * x**k for k, ck in enumerate(c))


def _dpoly(c, x):
    return sum(k * ck * x ** (k - 1) for k, ck in enumerate(c) if k)


def t_l(eta: float, n: float, l: int) -> float:
    """T_l(eta) of a polytrope of index n (1.5, 2 or 3) for the multipole l (2 or 3)."""
    c = PZM93[float(n)][l]
    x = math.log10(max(eta, 1e-12))
    x = max(x, 0.0)  # the stars overlap there; the run stops at the collision anyway
    if x <= 1.0:
        return 10.0 ** _poly(c, x)
    slope = min(_dpoly(c, 1.0), -2.0)
    return 10.0 ** (_poly(c, 1.0) + slope * (x - 1.0))


def tidal_energy(m1, r1, n1, m2, r2, n2, rp) -> float:
    """Orbital energy [Msun AU^2 / yr^2] dissipated at one periastron passage at distance rp
    [AU] (radii r1, r2 in AU, polytropic indices n1, n2): Portegies Zwart & Meinen (1993),
    eq. 4."""
    G = units.G
    e = 0.0
    for m, r, n, mo in ((m1, r1, n1, m2), (m2, r2, n2, m1)):
        if r <= 0 or n is None:
            continue
        eta = math.sqrt(m / (m + mo)) * (rp / r) ** 1.5
        e += G * mo**2 / r * ((r / rp) ** 6 * t_l(eta, n, 2) + (r / rp) ** 8 * t_l(eta, n, 3))
    return e


def capture_radius(m1, r1, n1, m2, r2, n2, v_inf_auyr) -> float:
    """Largest periastron distance [AU] at which one passage removes the orbital energy
    (1/2) mu v_inf^2 of an unbound encounter (bisection between the touching distance and
    100 stellar radii). Reference implementation (tests/validation)."""
    mu = m1 * m2 / (m1 + m2)
    e_orb = 0.5 * mu * v_inf_auyr**2
    lo, hi = r1 + r2, 100.0 * max(r1, r2)
    if tidal_energy(m1, r1, n1, m2, r2, n2, lo) < e_orb:
        return float("nan")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if tidal_energy(m1, r1, n1, m2, r2, n2, mid) >= e_orb:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def default_polytrope(mass: float) -> float:
    """n = 1.5 for fully convective stars (below 0.4 Msun) and giants are given explicitly;
    n = 3 for main-sequence stars otherwise."""
    return 1.5 if mass < 0.4 else 3.0


class PeriastronDissipation:
    """REBOUND heartbeat that removes E_tide from every pair of stars at its periastron
    (the sign of r.v changes from - to +). The pair's centre of mass is kept."""

    def __init__(self, sim, masses, radii, polytropes):
        self.m = np.asarray(masses, float)
        self.r = np.asarray(radii, float)
        self.n = list(polytropes)
        k = len(self.m)
        self.pairs = [
            (i, j) for i in range(k) for j in range(i + 1, k) if self.r[i] > 0 and self.r[j] > 0
        ]
        self.prev = {}
        self.dissipated = 0.0
        self.events = 0
        self.sim = sim

    def __call__(self, sim_ptr):
        sim = self.sim
        ps = sim.particles
        for i, j in self.pairs:
            pi, pj = ps[i], ps[j]
            dx, dy, dz = pj.x - pi.x, pj.y - pi.y, pj.z - pi.z
            d = math.sqrt(dx * dx + dy * dy + dz * dz)
            reach = 12.0 * (self.r[i] + self.r[j])
            if d > reach:
                self.prev[(i, j)] = None
                continue
            dvx, dvy, dvz = pj.vx - pi.vx, pj.vy - pi.vy, pj.vz - pi.vz
            rv = dx * dvx + dy * dvy + dz * dvz
            before = self.prev.get((i, j))
            self.prev[(i, j)] = rv
            if before is None or not (before < 0.0 <= rv):
                continue
            # periastron of the osculating two-body orbit of the pair
            mu_g = units.G * (self.m[i] + self.m[j])
            v2 = dvx * dvx + dvy * dvy + dvz * dvz
            hx, hy, hz = dy * dvz - dz * dvy, dz * dvx - dx * dvz, dx * dvy - dy * dvx
            h2 = hx * hx + hy * hy + hz * hz
            energy = 0.5 * v2 - mu_g / d
            ecc = math.sqrt(max(1.0 + 2.0 * energy * h2 / mu_g**2, 0.0))
            q = h2 / mu_g / (1.0 + ecc)
            de = tidal_energy(self.m[i], self.r[i], self.n[i], self.m[j], self.r[j], self.n[j], q)
            mu = self.m[i] * self.m[j] / (self.m[i] + self.m[j])
            v2_new = v2 - 2.0 * de / mu
            if v2_new <= 0.0:
                v2_new = 1e-6 * v2
            f = math.sqrt(v2_new / v2) - 1.0
            M = self.m[i] + self.m[j]
            # change the relative velocity by f (dv), keeping the centre-of-mass velocity
            ddx, ddy, ddz = f * dvx, f * dvy, f * dvz
            pi.vx -= self.m[j] / M * ddx
            pi.vy -= self.m[j] / M * ddy
            pi.vz -= self.m[j] / M * ddz
            pj.vx += self.m[i] / M * ddx
            pj.vy += self.m[i] / M * ddy
            pj.vz += self.m[i] / M * ddz
            self.dissipated += 0.5 * mu * (v2 - v2_new)
            self.events += 1


__all__ = ["PZM93", "PeriastronDissipation", "capture_radius", "t_l", "tidal_energy"]
