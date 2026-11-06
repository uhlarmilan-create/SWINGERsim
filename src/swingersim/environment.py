"""The environment of a multiple star: the tide of the Galactic disc and passing stars.

Galactic tide (physics.galactic_tide). In the solar neighbourhood the vertical tide of the
disc dominates (Heisler & Tremaine 1986). Its potential is 2 pi G rho0 z^2 (z along the
Galactic pole; Grishin & Perets 2022, eq. 3) and, averaged over an orbit of semi-major axis
a, the interaction energy of the orbit (reduced mass mu) is

    H = pi G rho0 mu a^2 [ (1 - e^2) - (j . z)^2 + 5 (e . z)^2 ]          (eq. 4)

with the vectorial elements j = sqrt(1 - e^2) k and e. It drives ZLK-like cycles of the
eccentricity with e_max = sqrt(1 - (5/4) cos^2 i) above the critical inclination
cos^2 i < 4/5 (eq. 6), on the time scale tau_GT = 3 / (8 pi G rho0) sqrt(G M / a^3)
= 0.89 Gyr (M / Msun)^(1/2) (a / 10^4 AU)^(-3/2) (eq. 7). swingersim adds the corresponding
terms to the secular equations of every orbit.

Encounters (physics.encounters). Stars passing at a periapsis distance Q much larger than an
orbit change its eccentricity and orientation but hardly its energy. To first order in the
perturbation the changes of e and j after one passage are (Hamers & Samsing 2019, eq. 9; in
the frame of the perturber's orbit, which lies in the x-y plane with its periapsis on x)

    eps = [M^2 / (m (m + M)) (a / Q)^3]^(1/2) (1 + E)^(-3/2),   L = arccos(-1/E),

with the perturber mass M, the mass of the orbit m and the eccentricity of the perturber's
orbit E. swingersim draws a Poisson sequence of encounters (rate n pi b_max^2 <v>, relative
speeds from a Maxwellian with the given dispersion weighted by v, impact parameters uniform
in area, isotropic orientations) and applies the kicks to every orbit for which the passage is
distant (Q > 3 a (1 + e)); closer passages are counted, not applied.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from . import units

PC_AU = 206264.806
KMS_AUYR = units.SECONDS_PER_YEAR / units.AU_KM  # 1 km/s in AU/yr (Julian year)


def rho_au(rho_msun_pc3: float) -> float:
    return rho_msun_pc3 / PC_AU**3


def pole_vector(i_deg: float, omega_deg: float) -> np.ndarray:
    """Unit vector of the Galactic pole (same convention as an orbit normal)."""
    i, o = math.radians(i_deg), math.radians(omega_deg)
    return np.array([math.sin(i) * math.sin(o), -math.sin(i) * math.cos(o), math.cos(i)])


def tide_time_scale_yr(m_tot: float, a_au: float, rho_msun_pc3: float = 0.185) -> float:
    """tau_GT = 3 / (8 pi G rho0) sqrt(G M / a^3) [yr] (Grishin & Perets 2022, eq. 7).
    Reference implementation (tests/validation)."""
    rho = rho_au(rho_msun_pc3)
    return 3.0 / (8.0 * math.pi * units.G * rho) * math.sqrt(units.G * m_tot / a_au**3)


def tide_emax(inc_deg: float) -> float:
    """Largest eccentricity reached from a circular orbit inclined by inc to the disc.
    Reference implementation (tests/validation)."""
    c2 = math.cos(math.radians(inc_deg)) ** 2
    return math.sqrt(max(1.0 - 1.25 * c2, 0.0))


def tide_derivatives(K, p, c, z):
    """(dK/dt, dp/dt) of one orbit from H = c [ K.K - (K.z)^2 + 5 (p.z)^2 ] (K.K = 1 - e^2)
    with the equations of motion of swingersim.secular.nested, divided by L."""
    gK = 2.0 * c * (K - (K @ z) * z)
    gp = 10.0 * c * (p @ z) * z
    return -(np.cross(p, gp) + np.cross(K, gK)), -(np.cross(K, gp) + np.cross(p, gK))


def tide_energy(K, p, c, z) -> float:
    return float(c * (1.0 - p @ p - (K @ z) ** 2 + 5.0 * (p @ z) ** 2))


# ------------------------------------------------------------------ encounters
def first_order_kick(e, j, a, m, m_per, q, ecc_per):
    """(Delta e, Delta j) of an orbit (vectors e, j in the frame of the perturber's orbit,
    x along its periapsis, z along its angular momentum) for one passage: Hamers & Samsing
    (2019) eq. 9. a, q in AU, masses in Msun, ecc_per >= 1."""
    E = max(ecc_per, 1.0)
    eps = math.sqrt(m_per**2 / (m * (m + m_per)) * (a / q) ** 3) * (1.0 + E) ** -1.5
    s = math.sqrt(max(1.0 - 1.0 / E**2, 0.0))
    L = math.acos(-1.0 / E)
    ex, ey, ez = e
    jx, jy, jz = j
    de = np.array(
        [
            eps
            / (2 * E)
            * (
                s * (ez * jy * (1 - 10 * E**2) + ey * jz * (-5 + 2 * E**2))
                - 3 * E * (3 * ez * jy + ey * jz) * L
            ),
            eps
            / (2 * E)
            * (
                s * (ez * jx * (1 + 8 * E**2) + ex * jz * (-5 + 8 * E**2))
                + 3 * E * (3 * ez * jx + ex * jz) * L
            ),
            eps
            / E
            * (
                s * (ey * jx * (2 + E**2) + ex * jy * (2 - 5 * E**2))
                + 3 * E * L * (ey * jx - ex * jy)
            ),
        ]
    )
    dj = np.array(
        [
            -eps / (2 * E) * (5 * ey * ez - jy * jz) * (s * (1 + 2 * E**2) + 3 * E * L),
            eps / (2 * E) * (5 * ex * ez - jx * jz) * (s * (-1 + 4 * E**2) + 3 * E * L),
            eps / E * (5 * ex * ey - jx * jy) * s * (1 - E**2),
        ]
    )
    return de, dj


@dataclass
class Encounter:
    t_yr: float
    q_au: float
    ecc: float
    rotation: np.ndarray  # rows: x (periapsis), y, z (normal) of the perturber orbit
    m_per: float


def encounter_sequence(spec, m_sys: float, a_out: float, t_max: float) -> list[Encounter]:
    """Poisson sequence of encounters over 0..t_max [yr] (sign of t_max kept)."""
    rng = np.random.default_rng(spec.seed)
    b_max = spec.b_max_au or 30.0 * a_out
    n = spec.n_star_pc3 / PC_AU**3
    sig = spec.sigma_kms * KMS_AUYR
    mean_v = sig * math.sqrt(8.0 / math.pi)
    rate = n * math.pi * b_max**2 * mean_v
    out = []
    t = 0.0
    T = abs(t_max)
    while True:
        t += rng.exponential(1.0 / rate)
        if t >= T:
            break
        # speed weighted by v: v^3 exp(-v^2 / 2 sigma^2) -> chi distribution with 4 dof
        v = sig * math.sqrt(rng.chisquare(4))
        b = b_max * math.sqrt(rng.random())
        mu = units.G * (m_sys + spec.m_perturber)
        q = -mu / v**2 + math.sqrt((mu / v**2) ** 2 + b**2)
        ecc = 1.0 + q * v**2 / mu
        # isotropic orientation: random rotation (QR of a Gaussian matrix)
        qm, rm = np.linalg.qr(rng.normal(size=(3, 3)))
        qm = qm * np.sign(np.diag(rm))
        if np.linalg.det(qm) < 0:
            qm[:, 0] = -qm[:, 0]
        out.append(Encounter(math.copysign(t, t_max), q, ecc, qm.T, spec.m_perturber))
    return out


def apply_encounter(enc: Encounter, K, p, a: float, m: float):
    """Kicked (K, p) of one orbit, or None when the passage is not distant."""
    e_mag = math.sqrt(p @ p)
    if enc.q_au < 3.0 * a * (1.0 + e_mag):
        return None
    R = enc.rotation
    de, dj = first_order_kick(R @ p, R @ K, a, m, enc.m_per, enc.q_au, enc.ecc)
    p2, K2 = p + R.T @ de, K + R.T @ dj
    # back onto the constraint surface: K.p = 0, |K|^2 + |p|^2 = 1
    p2 = p2 - (K2 @ p2) / (K2 @ K2) * K2
    norm = math.sqrt(K2 @ K2 + p2 @ p2)
    return K2 / norm, p2 / norm


__all__ = [
    "Encounter",
    "apply_encounter",
    "encounter_sequence",
    "first_order_kick",
    "pole_vector",
    "tide_derivatives",
    "tide_emax",
    "tide_energy",
    "tide_time_scale_yr",
]


# ------------------------------------------------------------------ N-body backend
def hyperbolic_state(mu: float, q: float, ecc: float, f: float):
    """Position and velocity (orbital plane: periapsis on x, motion towards +y) and the time
    since periapsis [yr] on a hyperbolic (or parabolic) orbit at true anomaly f."""
    p = q * (1.0 + ecc)
    r = p / (1.0 + ecc * math.cos(f))
    pos = np.array([r * math.cos(f), r * math.sin(f), 0.0])
    vel = math.sqrt(mu / p) * np.array([-math.sin(f), ecc + math.cos(f), 0.0])
    if abs(ecc - 1.0) < 1e-9:  # parabola (Barker's equation)
        D = math.tan(f / 2)
        t = math.sqrt(p**3 / mu) * 0.5 * (D + D**3 / 3.0)
    else:
        a = q / (1.0 - ecc)  # negative
        F = 2.0 * math.atanh(math.sqrt((ecc - 1.0) / (ecc + 1.0)) * math.tan(f / 2))
        t = math.sqrt((-a) ** 3 / mu) * (ecc * math.sinh(F) - F)
    return pos, vel, t


class NBodyEnvironment:
    """Galactic tide (velocity kicks of the tidal acceleration relative to the barycentre of the
    system, -4 pi G rho0 (dr . z) z, every ``kick_dt``) and passing stars (each perturber is
    integrated as a particle from distance ``r_start`` before to ``r_start`` after its
    periapsis) for the N-body backend."""

    def __init__(self, cfg, n_star: int, m_sys: float, a_out: float, p_out_yr: float, t_max):
        self.n = n_star
        self.sign = 1.0 if t_max >= 0 else -1.0
        gt = cfg.physics.galactic_tide
        self.tide = None
        if gt is not None:
            self.tide = (rho_au(gt.rho0_msun_pc3), pole_vector(gt.pole_i_deg, gt.pole_Omega_deg))
            self.kick_dt = p_out_yr / 400.0
        self.events: list[tuple[float, str, Encounter]] = []
        self._next = 0  # index of the first event not yet due
        self.counts = {"n": 0}
        spec = cfg.physics.encounters
        if spec is not None:
            mu = units.G * (m_sys + spec.m_perturber)
            for enc in encounter_sequence(spec, m_sys, a_out, t_max):
                r_start = max(20.0 * enc.q_au, 20.0 * a_out)
                f0 = math.acos(
                    max(min((enc.q_au * (1 + enc.ecc) / r_start - 1.0) / enc.ecc, 1.0), -1.0)
                )
                _, _, half = hyperbolic_state(mu, enc.q_au, enc.ecc, f0)
                t0, t1 = enc.t_yr - self.sign * half, enc.t_yr + self.sign * half
                if self.sign * t0 < 0:  # already approaching at the start: begin at t = 0
                    t0 = 0.0
                if self.sign * t1 > self.sign * t_max:
                    t1 = t_max
                self.events.append((t0, "add", enc))
                self.events.append((t1, "remove", enc))
            self.events.sort(key=lambda e: self.sign * e[0])
        self.last_kick = 0.0
        self.active: dict[int, int] = {}  # id(enc) -> particle index

    @property
    def enabled(self) -> bool:
        return self.tide is not None or bool(self.events)

    def boundaries(self, t0: float, t1: float) -> list[float]:
        s = self.sign
        out = []
        if self.tide is not None:
            k = math.floor(s * t0 / self.kick_dt) + 1
            while k * self.kick_dt < s * t1:
                out.append(s * k * self.kick_dt)
                k += 1
        out += [t for t, _, _ in self.events if s * t0 < s * t <= s * t1]
        out = sorted(set(out), key=lambda x: s * x)
        if not out or out[-1] != t1:
            out.append(t1)
        return out

    def _com(self, sim):
        ps = sim.particles
        m = np.array([ps[k].m for k in range(self.n)])
        r = np.array([[ps[k].x, ps[k].y, ps[k].z] for k in range(self.n)])
        return m, r, m @ r / m.sum()

    def tide_energy(self, sim) -> float:
        if self.tide is None:
            return 0.0
        m, r, com = self._com(sim)
        z = (r - com) @ self.tide[1]
        return float(2.0 * math.pi * units.G * self.tide[0] * np.sum(m * z * z))

    def at(self, sim, t: float) -> bool:
        """Apply what is due at time t; True when particles were added or removed."""
        changed = False
        if self.tide is not None and t != self.last_kick:
            dt = t - self.last_kick
            _m, r, com = self._com(sim)
            zhat = self.tide[1]
            acc = -4.0 * math.pi * units.G * self.tide[0] * ((r - com) @ zhat)
            for k in range(self.n):
                p = sim.particles[k]
                p.vx += acc[k] * zhat[0] * dt
                p.vy += acc[k] * zhat[1] * dt
                p.vz += acc[k] * zhat[2] * dt
            self.last_kick = t
        # the events are sorted in the run direction and fall on the boundaries (t == te)
        ev, j, s = self.events, self._next, self.sign
        while j < len(ev) and s * ev[j][0] < s * t:
            j += 1  # passed without a boundary on it (not applied, as before)
        while j < len(ev) and ev[j][0] == t:
            changed = True
            if ev[j][1] == "add":
                self._add(sim, ev[j][2])
            else:
                self._remove(sim, ev[j][2])
            j += 1
        self._next = j
        return changed

    def _add(self, sim, enc: Encounter):
        m, _r, com = self._com(sim)
        ps = sim.particles
        v = np.array([[ps[k].vx, ps[k].vy, ps[k].vz] for k in range(self.n)])
        vcom = m @ v / m.sum()
        mu = units.G * (m.sum() + enc.m_per)
        # position at time sim.t relative to the periapsis time enc.t_yr
        dt = sim.t - enc.t_yr
        lo, hi = -math.pi * 0.999999, math.pi * 0.999999
        if enc.ecc > 1:
            lim = math.acos(-1.0 / enc.ecc) * 0.999999
            lo, hi = -lim, lim
        for _ in range(200):  # bisection on the true anomaly for the time since periapsis
            mid = 0.5 * (lo + hi)
            if hyperbolic_state(mu, enc.q_au, enc.ecc, mid)[2] < dt:
                lo = mid
            else:
                hi = mid
        pos, vel, _ = hyperbolic_state(mu, enc.q_au, enc.ecc, 0.5 * (lo + hi))
        # the relative orbit belongs to the perturber about the system's barycentre
        R = enc.rotation.T
        rp, vp = R @ pos, R @ vel
        M = m.sum()
        frac = M / (M + enc.m_per)
        sim.add(
            m=enc.m_per,
            x=com[0] + frac * rp[0],
            y=com[1] + frac * rp[1],
            z=com[2] + frac * rp[2],
            vx=vcom[0] + frac * vp[0],
            vy=vcom[1] + frac * vp[1],
            vz=vcom[2] + frac * vp[2],
            r=0.0,
        )
        shift = enc.m_per / (M + enc.m_per)
        for k in range(self.n):  # keep the barycentre of everything fixed
            ps[k].x -= shift * rp[0]
            ps[k].y -= shift * rp[1]
            ps[k].z -= shift * rp[2]
            ps[k].vx -= shift * vp[0]
            ps[k].vy -= shift * vp[1]
            ps[k].vz -= shift * vp[2]
        self.active[id(enc)] = len(sim.particles) - 1
        self.counts["n"] += 1

    def _remove(self, sim, enc: Encounter):
        idx = self.active.pop(id(enc), None)
        if idx is not None:
            sim.remove(idx)
            # indices of later perturbers shift down by one
            self.active = {k: (v - 1 if v > idx else v) for k, v in self.active.items()}
