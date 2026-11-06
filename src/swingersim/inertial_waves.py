"""Frequency-averaged tidal dissipation of inertial waves in convective regions.

Inertial waves are waves restored by the Coriolis force; a tidal forcing with frequency
|omega| < 2 Omega (Omega = spin of the star) excites them in a convective envelope and their
dissipation, averaged over the tidal frequency, is the "dynamical tide" of convective zones
(Ogilvie 2013; Mathis 2015; Barker 2020, 2022). It scales with eps^2 = Omega^2 R^3 / (G M):

    1 / <Q'_IW> = eps^2 * C(structure).

This module computes C in two ways:

* ``barker_profile``: the impulsive formalism for a continuous radial profile rho(r),
  d rho / d p (Barker 2022, ApJL 927, L36, Appendix A; Barker 2020, MNRAS 498, 2270, eqs. 12-38,
  l = m = 2): the
  irrotational equilibrium tide X (A4) with the self-gravity of the perturbation (A7), the
  pressure perturbation W of the waves (A10) with vanishing radial velocity at the boundaries
  of the convective zone (A11), and the kinetic energy of the impulsively excited waves
  (A17-A20). Units G = M = R = 1.
* ``two_layer``: the closed formula of Ogilvie (2013) for a homogeneous core (radius alpha R,
  mass beta M) inside a homogeneous envelope, 1/<Q'> = (2/3) <D> (Barker 2020, eq. 39;
  <Q'> = 3/(2 <D>) as in Gallet et al. 2017). The numerical method reproduces it for a
  smoothed two-layer profile to < 0.2 % independently of alpha and beta (test and
  benchmark), which checks the solver and the normalisation of eqs. 30-38 of Barker (2020).

A fully convective star (no radiative core) is an n = 1.5 polytrope: C = 1 / 130.9.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

L_TIDE = 2
M_TIDE = 2


def _laplacian_rows(r: np.ndarray, rho: np.ndarray):
    """Tridiagonal rows of (1/r^2) d/dr (r^2 rho du/dr) - l(l+1) rho u / r^2 (conservative
    finite differences on a non-uniform grid); boundary rows are set by ``_solve``."""
    n = len(r)
    dr = np.diff(r)
    flux = (0.5 * (r[1:] + r[:-1])) ** 2 * 0.5 * (rho[1:] + rho[:-1]) / dr
    vol = 0.5 * (r[2:] - r[:-2]) * r[1:-1] ** 2
    lo, di, up = np.zeros(n), np.zeros(n), np.zeros(n)
    lo[1:-1] = flux[:-1] / vol
    up[1:-1] = flux[1:] / vol
    di[1:-1] = -(flux[:-1] + flux[1:]) / vol - L_TIDE * (L_TIDE + 1) * rho[1:-1] / r[1:-1] ** 2
    return lo, di, up


def _solve(r, lo, di, up, rhs, bc0, bc1):
    """Solve the tridiagonal system; bc = (c_u, c_du, value) means c_u u + c_du u' = value at
    that end (second-order one-sided derivative)."""
    from scipy.sparse import diags
    from scipy.sparse.linalg import spsolve

    n = len(r)
    a = diags([lo[1:], di, up[:-1]], [-1, 0, 1], format="lil")
    b = np.array(rhs, dtype=float)
    for k, (cu, cd, v), idx in ((0, bc0, [0, 1, 2]), (n - 1, bc1, [n - 1, n - 2, n - 3])):
        x0, x1, x2 = r[idx]
        w1 = (x0 - x2) / ((x1 - x0) * (x1 - x2))
        w2 = (x0 - x1) / ((x2 - x0) * (x2 - x1))
        w0 = -(w1 + w2)
        a[k, :] = 0.0
        a[k, idx[0]] = cu + cd * w0
        a[k, idx[1]] = cd * w1
        a[k, idx[2]] = cd * w2
        b[k] = v
    return spsolve(a.tocsr(), b)


def _mass(r, rho):
    seg = 0.5 * (rho[1:] * r[1:] ** 2 + rho[:-1] * r[:-1] ** 2) * np.diff(r) * 4.0 * np.pi
    return np.concatenate([[0.0], np.cumsum(seg)])


def barker_profile(
    r: np.ndarray, rho: np.ndarray, drho_dp: np.ndarray, convective: np.ndarray
) -> float:
    """C = 1 / (<Q'_IW> eps^2) of a star given on the grid r (r[0] > 0 small; surface at
    r = 1, units G = M = R = 1, total mass 1). ``convective`` marks one contiguous convective
    zone (the envelope, or the whole star when it starts at r[0])."""
    l, m = L_TIDE, M_TIDE
    g = _mass(r, rho) / r**2
    psi = r**l
    # perturbed gravitational potential (A7), whole star; regular centre (A8), vacuum (A9)
    lo, di, up = _laplacian_rows(r, np.ones_like(r))
    di = di + 4.0 * np.pi * rho * drho_dp
    phi = _solve(
        r,
        lo,
        di,
        up,
        -4.0 * np.pi * rho * drho_dp * psi,
        (l / r[0], -1.0, 0.0),
        ((l + 1) / r[-1], 1.0, 0.0),
    )
    idx = np.flatnonzero(convective)
    i0, i1 = idx[0], idx[-1] + 1
    rc, rh = r[i0:i1], rho[i0:i1]
    forcing = phi[i0:i1] + psi[i0:i1]
    full = i0 == 0
    lo, di, up = _laplacian_rows(rc, rh)
    # equilibrium tide (A4) with the boundary conditions (A5, A6)
    bc0 = (l / rc[0], -1.0, 0.0) if full else (0.0, 1.0, -forcing[0] / g[i0])
    bc1 = (0.0, 1.0, -forcing[-1] / g[i1 - 1])
    x = _solve(rc, lo, di, up, rh * drho_dp[i0:i1] * forcing, bc0, bc1)  # eq. 17: + sign
    dx = np.gradient(x, rc)
    # inertial waves, W = i w (A10, A11); Omega = 1
    rhs = 2.0 * m * np.gradient(rh, rc) * x / rc
    bc0 = (l / rc[0], -1.0, 0.0) if full else (0.0, 1.0, 2.0 * m * x[0] / rc[0])
    bc1 = (0.0, 1.0, 2.0 * m * x[-1] / rc[-1])
    w = _solve(rc, lo, di, up, rhs, bc0, bc1)
    dw = np.gradient(w, rc)

    def q(k):
        return np.sqrt((k * k - m * m) / (4.0 * k * k - 1.0)) / k

    a = 2.0 * m * x / rc - dw  # (A12), times -i
    b = 2.0 * m / (l * (l + 1) * rc**2) * (rc * dx + x) - w / rc**2  # (A13)
    c_lo = 2.0 * q(l) / rc**2 * (rc * dx + (l + 1) * x)  # (A14)
    c_hi = -2.0 * q(l + 1) / rc**2 * (rc * dx - l * x)  # (A15)
    integrand = (
        rh
        * rc**2
        * (
            a**2
            + l * (l + 1) * rc**2 * b**2
            + l * (l - 1) * rc**2 * c_lo**2
            + (l + 1) * (l + 2) * rc**2 * c_hi**2
        )
    )
    energy = 0.25 * np.trapezoid(integrand, rc)
    return float(32.0 * np.pi**2 / (3.0 * (2 * l + 1)) * energy)  # (A17), |A| = 1


def ogilvie_d(alpha: float, beta: float) -> float:
    """<D>_omega / eps^2 of the two-layer model (Ogilvie 2013; eq. 1 of Mathis 2015)."""
    if not (0.0 < alpha < 1.0 and 0.0 < beta < 1.0):
        return 0.0
    g = alpha**3 * (1.0 - beta) / (beta * (1.0 - alpha**3))
    if not 0.0 < g < 1.0:
        return 0.0
    a3, a5 = alpha**3, alpha**5
    num = (
        100.0
        * np.pi
        / 63.0
        * a5
        / (1.0 - a5)
        * (1.0 - g) ** 2
        * (1.0 - alpha) ** 4
        * (1.0 + 2.0 * alpha + 3.0 * alpha**2 + 1.5 * a3) ** 2
        * (1.0 + (1.0 - g) / g * a3)
    )
    den = 1.0 + 1.5 * g + 2.5 / g * (1.0 + 0.5 * g - 1.5 * g**2) * a3 - 2.25 * (1.0 - g) * a5
    return float(num / den**2)


def two_layer(alpha: float, beta: float) -> float:
    """C = 1/(<Q'_IW> eps^2) of the two-layer model, (2/3) <D> (Barker 2020, eq. 39)."""
    return 2.0 / 3.0 * ogilvie_d(alpha, beta)


def two_layer_profile(alpha: float, beta: float, width: float = 2.5e-4, n: int = 20000):
    """Smoothed two-layer star (tanh steps of the given width at the interface and at the
    surface) for checking ``barker_profile`` against the closed formula."""
    base = np.linspace(1e-3, 1.0 + 10.0 * width, n)
    r = np.unique(
        np.concatenate(
            [
                base,
                alpha + np.linspace(-12.0 * width, 12.0 * width, n // 4),
                1.0 + np.linspace(-12.0 * width, 10.0 * width, n // 4),
            ]
        )
    )
    rho_c = beta / (4.0 / 3.0 * np.pi * alpha**3)
    rho_e = (1.0 - beta) / (4.0 / 3.0 * np.pi * (1.0 - alpha**3))
    step = 0.5 * (1.0 - np.tanh((r - alpha) / width))
    surface = 0.5 * (1.0 - np.tanh((r - 1.0) / width))
    rho = (rho_e + (rho_c - rho_e) * step) * surface
    g = _mass(r, rho) / r**2
    with np.errstate(divide="ignore", invalid="ignore"):
        drho_dp = np.nan_to_num(np.gradient(rho, r) / (-rho * g))
    convective = (r > alpha + 6.0 * width) & (r < 1.0 - 6.0 * width)
    return r, rho, drho_dp, convective


def polytrope_profile(n_index: float, n: int = 16000):
    """Lane-Emden polytrope of index n_index normalised to G = M = R = 1. Reference
    implementation (tests/validation)."""
    from scipy.integrate import solve_ivp

    def f(x, y):
        return [y[1], -(max(y[0], 0.0) ** n_index) - 2.0 * y[1] / x]

    def surface(x, y):
        return y[0]

    surface.terminal, surface.direction = True, -1
    x0 = 1e-6
    sol = solve_ivp(
        f,
        [x0, 20.0],
        [1.0 - x0**2 / 6.0, -x0 / 3.0],
        events=surface,
        rtol=1e-12,
        atol=1e-14,
        dense_output=True,
    )
    x1 = sol.t_events[0][0]
    x = np.linspace(x0, x1, n)
    theta = np.clip(sol.sol(x)[0], 0.0, None)
    r = x / x1
    rho = theta**n_index
    rho = rho / _mass(r, rho)[-1]
    g = _mass(r, rho) / r**2
    dpdr = -rho * g
    p = np.concatenate([[0.0], np.cumsum(-0.5 * (dpdr[1:] + dpdr[:-1]) * np.diff(r))])
    p = p[-1] - p
    with np.errstate(divide="ignore", invalid="ignore"):
        drho_dp = np.where(p > 0, n_index * rho / ((n_index + 1.0) * p), 0.0)
    keep = slice(0, n - 2)  # drop the zero-density surface point
    return r[keep], rho[keep], drho_dp[keep]


@lru_cache(maxsize=4)
def polytrope_c(n_index: float = 1.5) -> float:
    """C = 1/(<Q'_IW> eps^2) of a fully convective polytrope (n = 1.5: 1/130.9). Reference
    implementation (tests/validation): FULLY_CONVECTIVE_C holds its value."""
    r, rho, drho_dp = polytrope_profile(n_index)
    return barker_profile(r, rho, drho_dp, np.ones(len(r), bool))


FULLY_CONVECTIVE_C = 1.0 / 130.9  # polytrope_c(1.5), checked in the tests


def structure_c(alpha: float, beta: float) -> float:
    """C for a star with a radiative core (two-layer) or fully convective (alpha, beta ~ 0)."""
    if alpha <= 0.02 or beta <= 0.02:
        return FULLY_CONVECTIVE_C
    return two_layer(alpha, beta)
