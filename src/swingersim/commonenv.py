"""Common-envelope evolution with the energy formalism (physics.mass_transfer.common_envelope).

When mass transfer becomes dynamically unstable (q above q_crit, swingersim.masstransfer) the
donor's envelope engulfs both stars. swingersim does not follow the hydrodynamics of the envelope;
it computes the outcome with the energy formalism of BSE (Hurley, Tout & Pols 2002,
Sect. 2.7.1, after Webbink 1984 and Tout et al. 1997):

* binding energy of the envelopes (eq. 69): E_bind = -(G / lambda) (M1 M_env1 / R1 +
  M2 M_env2' / R2), with effective values for the companion: a main-sequence companion counts
  as all core (M_c2' = M2, R_c2' = R2, M_env2' = 0);
* orbital energy of the cores before and after (eqs. 70, 72): E_orb = -G M_c1 M_c2' / (2 a);
* the cores spiral in, giving orbital energy to the envelope with efficiency alpha_CE
  (eq. 71): E_bind = alpha_CE (E_orb,f - E_orb,i), which gives the final separation a_f;
* if a core would fill its Roche lobe (Eggleton 1983) at a_f, the cores merge at the
  separation a_L where that happens, and the merger keeps the part of the envelope the
  released energy could not unbind (eqs. 73-77; R ~ M^-x on the dynamical time scale,
  x = 0.3).

Core radii: a degenerate helium core (stars below the helium-flash mass, ~2 Msun) or a
carbon-oxygen core on the AGB takes the white-dwarf radius (Hurley, Pols & Tout 2000,
eq. 91), a non-degenerate helium core the radius of a zero-age helium star (eq. 78).

The envelope leaves the binary isotropically, in swingersim over ``duration_yr``; when that is
much shorter than an outer period the outer orbits receive an impulsive mass loss
(Hadjidemetriou 1963): they widen, and become unbound when the system loses half of its mass
or more (circular orbit) - followed exactly by the N-body backend. After the envelope is gone
the core is kept as a point mass of constant mass (swingersim does not evolve remnants); the
companion continues on the track of its present mass.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .stellar import roche_lobe_eggleton
from .units import RSUN_AU

X_GIANT = 0.3  # R ~ M^-x of giants on the dynamical time scale (Hurley et al. 2002)
M_CH = 1.44  # Chandrasekhar mass [Msun]
M_HEF = 2.0  # helium-flash mass at Z = 0.02 [Msun]: lighter stars have degenerate cores


def r_white_dwarf(m: float) -> float:
    """White-dwarf radius [Rsun], Hurley et al. (2000) eq. 91 (neutron-star floor 1.4e-5)."""
    m = min(max(m, 1e-3), M_CH * 0.999)
    return max(1.4e-5, 0.0115 * math.sqrt((M_CH / m) ** (2 / 3) - (m / M_CH) ** (2 / 3)))


def r_helium_zams(m: float) -> float:
    """Radius [Rsun] of a zero-age helium main-sequence star, Hurley et al. (2000) eq. 78."""
    return 0.2391 * m**4.6 / (m**4 + 0.162 * m**3 + 0.0065)


M_MASSIVE = 8.0  # [Msun] heavier stars do not form a degenerate core on the AGB phases


def core_radius(m_core: float, m_initial: float, phase: int) -> float:
    """Radius [Rsun] of the core of a giant (Hurley et al. 2000, Sect. 6.1): degenerate
    (white-dwarf radius) on the giant branch of stars below the helium-flash mass and on the
    AGB of stars below M_MASSIVE; a naked helium star (zero-age helium main sequence, eq. 78)
    otherwise - in core helium burning (the core is not degenerate after the helium flash) and
    in the late phases of massive stars."""
    if m_core <= 0.0:
        return 0.0
    if phase in (2, 3) and m_initial < M_HEF:
        return r_white_dwarf(m_core)
    if phase >= 5 and m_initial < M_MASSIVE:
        return r_white_dwarf(m_core)
    return r_helium_zams(m_core)


@dataclass
class CeOutcome:
    outcome: str  # "survive" or "merge"
    a_f: float  # final separation [AU] (survive) or separation of the merger a_L
    m1_f: float  # donor remnant (core) mass [Msun]; merger: merged mass in m1_f, m2_f = 0
    m2_f: float
    e_bind: float  # binding energy of the envelopes [Msun^2 / AU, G = 1]
    e_orb_i: float
    # survive: final orbital energy of the cores; merge: binding energy of the envelope left
    # around the merger (Hurley et al. 2002, eq. 73), not an orbital energy
    e_orb_f: float
    rc1: float  # core radii [Rsun]
    rc2: float


def outcome(
    m1: float,
    mc1: float,
    r1: float,
    rc1: float,
    m2: float,
    mc2: float,
    r2: float,
    rc2: float,
    a_i: float,
    alpha: float = 1.0,
    lam: float = 0.5,
    x: float = X_GIANT,
) -> CeOutcome:
    """Outcome of a common envelope (Hurley et al. 2002, eqs. 69-77). Star 1 is the donor
    (giant: mass m1, core mc1, radius r1, core radius rc1); star 2 the companion, given with
    its effective core (a main-sequence star: mc2 = m2, rc2 = r2). Masses in Msun, radii in
    Rsun, separations in AU; energies with G = 1."""
    menv1, menv2 = m1 - mc1, m2 - mc2
    r1_au, r2_au = r1 * RSUN_AU, r2 * RSUN_AU
    e_bind = -(m1 * menv1 / r1_au + (m2 * menv2 / r2_au if menv2 > 0 else 0.0)) / lam
    e_orb_i = -0.5 * mc1 * mc2 / a_i
    e_orb_f = e_orb_i + e_bind / alpha
    a_f = -0.5 * mc1 * mc2 / e_orb_f
    # Roche lobes of the cores at a_f
    rl1 = roche_lobe_eggleton(mc1 / mc2) * a_f / RSUN_AU
    rl2 = roche_lobe_eggleton(mc2 / mc1) * a_f / RSUN_AU
    if rc1 < rl1 and rc2 < rl2:
        return CeOutcome("survive", a_f, mc1, mc2, e_bind, e_orb_i, e_orb_f, rc1, rc2)
    # coalescence where the first core fills its lobe
    a_l = max(
        rc1 * RSUN_AU / roche_lobe_eggleton(mc1 / mc2),
        rc2 * RSUN_AU / roche_lobe_eggleton(mc2 / mc1),
    )
    e_bind_f = e_bind + alpha * (0.5 * mc1 * mc2 / a_l + e_orb_i)  # eq. 73
    mt = m1 + m2
    # core of the merger: both cores; a main-sequence companion is absorbed into the envelope
    mc3 = min(mc1 + (mc2 if menv2 > 0 else 0.0), mt)
    ratio = e_bind_f / e_bind if e_bind != 0 else 0.0
    if ratio <= 0.0 or mt - mc3 <= 0.0:
        mf = mc3
    else:
        ratio = min(ratio, 1.0)

        def f(mf):
            return (mf / mt) ** (1 + x) * (mf - mc3) / (mt - mc3) - ratio

        lo, hi = mc3, mt
        for _ in range(200):  # bisection: f increases monotonically from -ratio to 1 - ratio
            mid = 0.5 * (lo + hi)
            if f(mid) > 0:
                hi = mid
            else:
                lo = mid
        mf = 0.5 * (lo + hi)
    return CeOutcome("merge", a_l, mf, 0.0, e_bind, e_orb_i, e_bind_f, rc1, rc2)


def impulsive_orbit(m_before: float, m_after: float, a: float, e: float, true_anomaly: float):
    """Relative orbit after an instantaneous mass loss from a two-body orbit (positions and
    velocities unchanged, Hadjidemetriou 1963): (a', e') or (inf, e') when unbound. Masses in
    Msun (total of the orbit), a in AU, true anomaly in radians. Reference implementation
    (tests/validation)."""
    p = a * (1 - e * e)
    r = p / (1 + e * math.cos(true_anomaly))
    v2 = m_before * (2 / r - 1 / a)  # G = 1 (AU, Msun, yr / 2 pi)
    eps = 0.5 * v2 - m_after / r
    h2 = m_before * p  # specific angular momentum squared, unchanged
    e_new = math.sqrt(max(1 + 2 * eps * h2 / m_after**2, 0.0))
    if eps >= 0:
        return math.inf, e_new
    return -m_after / (2 * eps), e_new


__all__ = [
    "CeOutcome",
    "core_radius",
    "impulsive_orbit",
    "outcome",
    "r_helium_zams",
    "r_white_dwarf",
]
