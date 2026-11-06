"""Double-averaged quadrupole-order evolution of one inner orbit perturbed by a point mass.

Equations (A27), (A29) and (A36) of Naoz et al. (2013) MNRAS 431, 2155, valid for any mass
ratio (not only the test-particle limit): the inner eccentricity e1 and the argument of
periastron g1 (measured from the line of nodes in the invariable plane) evolve, the outer
orbit keeps e2 and G2 constant, and the mutual inclination follows from conservation of the
total angular momentum.

This is used for orientation: it predicts the eccentricity and inclination range and the
period of the ZLK / precession cycle of each orbit pair, which the N-body result is compared
with.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .. import units


@dataclass
class QuadrupolePrediction:
    e_min: float
    e_max: float
    i_min_deg: float
    i_max_deg: float
    period_yr: float  # period of the e / i oscillation (NaN if not found)
    g1_0_deg: float
    g1_ratio: float  # G1 / G2 at t = 0


def _angular_momenta(m_a, m_b, m_c, a1, a2):
    G = units.G
    m_in = m_a + m_b
    L1 = m_a * m_b * np.sqrt(G * a1 / m_in)
    L2 = m_in * m_c * np.sqrt(G * a2 / (m_in + m_c))
    return L1, L2


def c2_coefficient(m_a, m_b, m_c, a1, a2, e2):
    """C2 of Naoz et al. (2013) eq. (21), with k^2 = G."""
    G = units.G
    m_in = m_a + m_b
    L1, L2 = _angular_momenta(m_a, m_b, m_c, a1, a2)
    G2 = L2 * np.sqrt(1.0 - e2**2)
    return (
        G**2
        / 16.0
        * m_in**7
        * m_c**7
        / ((m_in + m_c) ** 3 * (m_a * m_b) ** 3)
        * L1**4
        / (L2**3 * G2**3)
    )


def predict(
    m_a: float,
    m_b: float,
    m_c: float,
    a1: float,
    a2: float,
    e1: float,
    e2: float,
    i_tot_deg: float,
    g1_deg: float,
    n_cycles: int = 3,
) -> QuadrupolePrediction:
    """Integrate the quadrupole equations over a few cycles and return ranges and period."""
    L1, L2 = _angular_momenta(m_a, m_b, m_c, a1, a2)
    G2 = L2 * np.sqrt(1.0 - e2**2)
    C2 = c2_coefficient(m_a, m_b, m_c, a1, a2, e2)
    e1 = max(e1, 1e-6)
    G1_0 = L1 * np.sqrt(1.0 - e1**2)
    cos_i0 = np.cos(np.radians(i_tot_deg))
    Gtot2 = G1_0**2 + G2**2 + 2.0 * G1_0 * G2 * cos_i0

    def cos_i(e):
        G1 = L1 * np.sqrt(1.0 - e**2)
        return np.clip((Gtot2 - G1**2 - G2**2) / (2.0 * G1 * G2), -1.0, 1.0)

    def rhs(_t, y):
        e, g = y
        e = min(max(e, 1e-9), 1.0 - 1e-12)
        G1 = L1 * np.sqrt(1.0 - e**2)
        c = cos_i(e)
        s2 = 1.0 - c**2
        de = 30.0 * C2 * e * (1.0 - e**2) / G1 * s2 * np.sin(2.0 * g)
        dg = (
            6.0
            * C2
            * (
                (4.0 * c**2 + (5.0 * np.cos(2.0 * g) - 1.0) * (1.0 - e**2 - c**2)) / G1
                + c / G2 * (2.0 + e**2 * (3.0 - 5.0 * np.cos(2.0 * g)))
            )
        )
        return [de, dg]

    # fastest secular time scale: g1 turns at a rate ~ C2 (1/G1 + 1/G2)
    t_scale = 1.0 / (C2 * (1.0 / G1_0 + 1.0 / G2))
    e_max_event = lambda t, y: rhs(t, y)[0]
    e_max_event.direction = -1.0
    e_max_event.terminal = n_cycles + 1  # stop after enough eccentricity maxima
    sol = solve_ivp(
        rhs,
        (0.0, 2000.0 * t_scale),
        [e1, np.radians(g1_deg)],
        method="DOP853",
        rtol=1e-10,
        atol=1e-12,
        events=e_max_event,
        dense_output=True,
        max_step=t_scale / 20.0,
    )
    t_ev = sol.t_events[0]
    if len(t_ev) >= n_cycles + 1:
        t_end = t_ev[n_cycles]
        if n_cycles > 1:  # full cycles between the maxima after the first one
            period = float((t_ev[n_cycles] - t_ev[1]) / (n_cycles - 1))
        else:  # one cycle: the interval between the two maxima found
            period = float(t_ev[1] - t_ev[0])
    else:
        t_end = sol.t[-1]
        period = float(np.diff(t_ev).mean()) if len(t_ev) >= 2 else float("nan")
    tt = np.linspace(0.0, t_end, 20000)
    e = sol.sol(tt)[0]
    if e.max() - e.min() < 1e-4:  # no eccentricity cycle (e.g. outside the ZLK window)
        period = float("nan")
    i = np.degrees(np.arccos(cos_i(e)))
    return QuadrupolePrediction(
        e_min=float(e.min()),
        e_max=float(e.max()),
        i_min_deg=float(i.min()),
        i_max_deg=float(i.max()),
        period_yr=period,
        g1_0_deg=float(g1_deg),
        g1_ratio=float(G1_0 / G2),
    )
