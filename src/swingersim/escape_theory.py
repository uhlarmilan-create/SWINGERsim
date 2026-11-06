"""Statistical (ergodic) theory of the break-up of a non-hierarchical three-body system.

Stone & Leigh (2019), Nature 576, 406: when three stars interact resonantly ("democratic"
phase), the outcome - which star escapes and the energy, angular momentum and orientation of
the binary left behind - is distributed according to the phase-space volume allowed by the
conservation of energy E0 and angular momentum L0. Their eq. (2), for escaper s and binary
(a, b) with M = m_a + m_b + m_s, m_B = m_a + m_b, m = m_B m_s / M:

    sigma_s = 2 pi^4 G^2 M^(5/2) m_B / (m_a m_b m_s)^(3/2)
              * Int L_B dE_B dL_B dC_B / [L_s (-E_B)^(3/2) (E0 - E_B)^(3/2)] * B,
    B = sqrt(2 M (E0-E_B) / (G^2 m_s^3 m_B^3)) sqrt(2 m (E0-E_B) R^2 + 2 G M m^2 R - L_s^2)
        - acosh[(1 + 2 (E0-E_B) R / (G m_s m_B)) / sqrt(1 + 2 M (E0-E_B) L_s^2/(G^2 m_s^3 m_B^3))],
    L_s^2 = L_B^2 (1 - C_B^2) + (L_B C_B - L0)^2,

with the strong-interaction radius R = alpha a_B (1 + e_B) ("apocentric escape", their fiducial
choice, alpha = 2) or R = alpha a_B. C_B is the cosine of the angle between the binary's angular
momentum and L0. swingersim evaluates the integral on a grid for each of the three possible
escapers; the totals give the escape probabilities, and the weights give the marginal
distributions of the binary's energy, eccentricity and inclination.

The theory describes resonant encounters only (about half of all low-velocity binary-single
scatterings in their experiments); prompt exchanges and fly-bys are not covered. Units: Msun,
AU, yr (G = swingersim.units.G); E0 < 0 (a bound three-body system) is required.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np

from . import units


@dataclass
class EscapeResult:
    probabilities: dict[int, float]  # escaper index -> probability
    e_grid: np.ndarray  # eccentricity bin centres
    e_density_all: np.ndarray  # dP/de over all escapers
    u_grid: np.ndarray  # E_B / E0 (>= 1) bin centres (log-spaced)
    u_density_all: np.ndarray  # dP/du over all escapers
    mean_e: float
    mean_u: float  # mean E_B / E0: how much more bound the final binary is


def _escaper_weights(masses, s, E0, L0, alpha, criterion, u, du, nl, nc):
    """Weights of the cells (u, e, C_B) of the grid for escaper s: u = E_B / E0 (bin centres
    u with widths du, a slice of the log-spaced grid), e and C_B uniform with nl and nc bins.
    Returns (w, e, C_B), arrays of shape (len(u), nl, nc)."""
    G = units.G
    a_i, b_i = [k for k in range(3) if k != s]
    ma, mb, ms = masses[a_i], masses[b_i], masses[s]
    mB = ma + mb
    M = mB + ms
    m = mB * ms / M
    muB = ma * mb / mB
    # grid: eccentricity e (L_B = L_c sqrt(1 - e^2)) and C_B
    e_edges = np.linspace(0.0, 1.0, nl + 1)
    ecc = 0.5 * (e_edges[1:] + e_edges[:-1])
    de = np.diff(e_edges)
    c_edges = np.linspace(-1.0, 1.0, nc + 1)
    C = 0.5 * (c_edges[1:] + c_edges[:-1])
    dc = np.diff(c_edges)
    U, e, CC = np.meshgrid(u, ecc, C, indexing="ij")
    EB = E0 * U  # < 0
    Es = E0 - EB  # escaper energy > 0
    a = G * ma * mb / (-2.0 * EB)
    Lc = muB * np.sqrt(G * mB * a)
    LB = Lc * np.sqrt(1.0 - e**2)
    R = alpha * a * (1.0 + e) if criterion == "AE" else alpha * a
    Ls2 = LB**2 * (1.0 - CC**2) + (LB * CC - L0) ** 2
    Ls = np.sqrt(np.maximum(Ls2, 1e-300))
    k = G**2 * ms**3 * mB**3
    arg = 2.0 * m * Es * R**2 + 2.0 * G * M * m**2 * R - Ls2
    first = np.sqrt(2.0 * M * Es / k) * np.sqrt(np.clip(arg, 0.0, None))
    ac = (1.0 + 2.0 * Es * R / (G * ms * mB)) / np.sqrt(1.0 + 2.0 * M * Es * Ls2 / k)
    second = np.arccosh(np.clip(ac, 1.0, None))
    B = np.where(arg > 0, np.clip(first - second, 0.0, None), 0.0)
    pref = 2.0 * np.pi**4 * G**2 * M**2.5 * mB / (ma * mb * ms) ** 1.5
    w = pref * LB / (Ls * (-EB) ** 1.5 * Es**1.5) * B
    # measure dE_B dL_B dC_B = |E0| du * Lc e / sqrt(1 - e^2) de * dC
    w = w * abs(E0) * du[:, None, None] * Lc * e / np.sqrt(1.0 - e**2)
    w = w * de[None, :, None] * dc[None, None, :]
    return w, e, CC


CHUNK_CELLS = 200_000  # grid cells per slice in u (bounds the memory: ~2 MB per array)


def escape_outcomes(
    masses,
    E0: float,
    L0: float,
    alpha: float = 2.0,
    criterion: str = "AE",
    nu: int = 160,
    nl: int = 200,
    nc: int = 64,
    u_max: float = 1.0e3,
    ne: int = 20,
) -> EscapeResult:
    """Outcome distribution of the break-up of a bound three-body system (E0 < 0)."""
    if E0 >= 0:
        raise ValueError("the statistical escape theory needs a bound system (E0 < 0)")
    masses = np.asarray(masses, float)
    e_edges = np.linspace(0.0, 1.0, ne + 1)
    e_c = 0.5 * (e_edges[1:] + e_edges[:-1])
    # u = E_B / E0 log-spaced in [1, u_max], integrated in slices of u (the full grid of
    # nu x nl x nc cells held ~20 temporary arrays: 377 MB for the default grid)
    lu = np.linspace(0.0, np.log(u_max), nu + 1)
    u_c = np.exp(0.5 * (lu[1:] + lu[:-1]))
    du = np.diff(np.exp(lu))
    step = max(1, CHUNK_CELLS // (nl * nc))
    totals = {}
    all_e = np.zeros(ne)
    u_hist = np.zeros(nu)
    sum_e = sum_u = 0.0
    for s in range(3):
        totals[s] = 0.0
        for k in range(0, nu, step):
            sl = slice(k, min(k + step, nu))
            w, e, _cc = _escaper_weights(
                masses, s, E0, L0, alpha, criterion, u_c[sl], du[sl], nl, nc
            )
            totals[s] += float(w.sum())
            all_e += np.histogram(e.ravel(), bins=e_edges, weights=w.ravel())[0]
            u_hist[sl] += w.sum(axis=(1, 2))
            sum_e += float((w * e).sum())
            sum_u += float((w * u_c[sl, None, None]).sum())
    tot = sum(totals.values())
    de = np.diff(e_edges)
    return EscapeResult(
        probabilities={s: totals[s] / tot for s in totals},
        e_grid=e_c,
        e_density_all=all_e / all_e.sum() / de,
        u_grid=u_c,
        u_density_all=u_hist / u_hist.sum() / du,
        mean_e=sum_e / tot,
        mean_u=sum_u / tot,
    )


def binary_single_integrals(m_bin, m_single, a_au, e_bin, v_inf_kms, b_au, cos_angle=1.0):
    """E0 and |L0| of a binary-single encounter: binary (m1, m2, a, e) plus a single star
    with velocity at infinity v_inf and impact parameter b; cos_angle = cosine of the angle
    between the binary's and the encounter's angular momentum."""
    G = units.G
    m1, m2 = m_bin
    mB = m1 + m2
    mu_inf = mB * m_single / (mB + m_single)
    v = v_inf_kms * units.SECONDS_PER_YEAR / units.AU_KM
    E0 = 0.5 * mu_inf * v**2 - G * m1 * m2 / (2.0 * a_au)
    Lb = m1 * m2 / mB * np.sqrt(G * mB * a_au * (1.0 - e_bin**2))
    Lo = mu_inf * b_au * v
    L0 = float(np.sqrt(Lb**2 + Lo**2 + 2.0 * Lb * Lo * cos_angle))
    return float(E0), L0


def angular_momentum_scale(masses, E0: float) -> float:
    """L_c(E0): the largest circular-orbit angular momentum of a pair of the stars at the
    total energy E0 (the unit of L0 in Stone & Leigh 2019)."""
    G = units.G
    best = 0.0
    for i, j in combinations(range(3), 2):
        ma, mb = masses[i], masses[j]
        mB = ma + mb
        a = G * ma * mb / (-2.0 * E0)
        best = max(best, ma * mb / mB * np.sqrt(G * mB * a))
    return best


def compare_with_scatter(rows: list[dict], masses, min_scrambles: int = 2, max_rows: int = 400):
    """Compare the theory with binary-single scattering runs (`swingersim scatter`): for the
    resonant encounters (at least `min_scrambles` phases without any bound pair) that ended as
    a binary plus an escaping star, the observed escape fractions and mean eccentricity against
    the theory averaged over the encounters' own E0 and L0."""
    import re

    sel = []
    for r in rows:
        if r.get("status") != "resolved" or r.get("scrambles", 0) < min_scrambles:
            continue
        if r.get("E0") is None or r["E0"] >= 0:
            continue
        parts = [p.strip() for p in r["hierarchy"].split("|")]
        singles = [p for p in parts if re.fullmatch(r"\d+", p)]
        pairs = [p for p in parts if re.fullmatch(r"\[\d+,\d+\]", p)]
        if len(parts) != 2 or len(singles) != 1 or len(pairs) != 1:
            continue
        orb = r["orbits"].get(pairs[0])
        if orb is None:
            continue
        sel.append((int(singles[0]) - 1, orb["e"], r["E0"], r["L0"]))
    if not sel:
        raise ValueError("no resonant binary + single outcomes in these rows")
    n = len(sel)
    obs = {s: sum(1 for x in sel if x[0] == s) / n for s in range(3)}
    obs_e = float(np.mean([x[1] for x in sel]))
    pred = {s: 0.0 for s in range(3)}
    pred_e = 0.0
    sub = sel[:: max(1, n // max_rows)]
    for _, _, E0, L0 in sub:
        res = escape_outcomes(masses, E0, L0, nu=80, nl=100, nc=32)
        for s in range(3):
            pred[s] += res.probabilities[s] / len(sub)
        pred_e += res.mean_e / len(sub)
    return {
        "n": n,
        "observed": obs,
        "predicted": pred,
        "observed_mean_e": obs_e,
        "predicted_mean_e": pred_e,
    }
