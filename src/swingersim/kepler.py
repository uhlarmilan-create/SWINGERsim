"""Two-body Keplerian conversions (bound orbits, e < 1).

Reference frame convention: the x-y plane is the reference plane (the sky
plane when the elements are given relative to the sky), z is its normal.
Angles in radians. Rotation: R = Rz(Omega) Rx(inc) Rz(omega).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

TWO_PI = 2.0 * np.pi


@dataclass(frozen=True)
class Elements:
    a: float  # semi-major axis [AU]
    e: float  # eccentricity
    inc: float  # inclination [rad]
    Omega: float  # longitude of ascending node [rad]
    omega: float  # argument of pericentre [rad]
    M: float  # mean anomaly [rad]


def solve_kepler(M: float, e: float, max_iter: int = 100) -> float:
    """Eccentric anomaly E from mean anomaly M (Newton iteration, to machine precision)."""
    M = np.mod(M, TWO_PI)
    E = M if e < 0.8 else np.pi
    tol = 8.0 * np.finfo(float).eps * TWO_PI
    for _ in range(max_iter):
        f = E - e * np.sin(E) - M
        dE = -f / (1.0 - e * np.cos(E))
        E += dE
        if abs(dE) < tol:
            return E
    if abs(E - e * np.sin(E) - M) < 1e-12:
        return E
    raise RuntimeError(f"Kepler equation did not converge (M={M}, e={e})")


def solve_kepler_array(M, e, max_iter: int = 100) -> np.ndarray:
    """Eccentric anomaly E from mean anomaly M, vectorised (arrays of M and e, broadcast).

    The same Newton iteration as solve_kepler (start at M for e < 0.8, at pi otherwise),
    iterated until every element converged to machine precision."""
    M = np.mod(np.asarray(M, float), TWO_PI)
    e = np.broadcast_to(np.asarray(e, float), M.shape)
    E = np.where(e < 0.8, M, np.pi)
    tol = 8.0 * np.finfo(float).eps * TWO_PI
    for _ in range(max_iter):
        dE = -(E - e * np.sin(E) - M) / (1.0 - e * np.cos(E))
        E = E + dE
        if np.all(np.abs(dE) < tol):
            return E
    if np.all(np.abs(E - e * np.sin(E) - M) < 1e-12):
        return E
    raise RuntimeError("Kepler equation did not converge for some elements")


def cross3(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Cross product of arrays of 3-vectors (..., 3): the arithmetic of np.cross (the same
    numbers) without its per-call overhead, which dominates for a few vectors."""
    a0, a1, a2 = a[..., 0], a[..., 1], a[..., 2]
    b0, b1, b2 = b[..., 0], b[..., 1], b[..., 2]
    return np.stack([a1 * b2 - a2 * b1, a2 * b0 - a0 * b2, a0 * b1 - a1 * b0], axis=-1)


def eccentricity_vector(mu, r, v) -> np.ndarray:
    """Eccentricity (Laplace-Runge-Lenz) vector e = v x h / mu - r / |r| of relative orbits;
    r, v of shape (..., 3), mu a scalar or of shape (...)."""
    r = np.asarray(r, float)
    v = np.asarray(v, float)
    h = cross3(r, v)
    rn = np.sqrt(np.einsum("...i,...i->...", r, r))
    return cross3(v, h) / np.asarray(mu, float)[..., None] - r / rn[..., None]


def states_to_elements(mu, r, v) -> dict[str, np.ndarray]:
    """Osculating elements of many relative orbits at once (r, v of shape (n, 3), mu of shape
    (n,)): a dict of arrays a, e, inc, Omega, omega, M (keys of Elements). The formulas and
    conventions are those of state_to_elements; orbits that are not bound get NaN."""
    r = np.atleast_2d(np.asarray(r, float))
    v = np.atleast_2d(np.asarray(v, float))
    mu = np.broadcast_to(np.asarray(mu, float), r.shape[:1])
    rn = np.sqrt(np.einsum("ij,ij->i", r, r))
    h = cross3(r, v)
    hn = np.sqrt(np.einsum("ij,ij->i", h, h))
    energy = 0.5 * np.einsum("ij,ij->i", v, v) - mu / rn
    bound = energy < 0.0
    with np.errstate(invalid="ignore", divide="ignore"):
        a = -mu / (2.0 * energy)
        e_vec = cross3(v, h) / mu[:, None] - r / rn[:, None]
        e = np.sqrt(np.einsum("ij,ij->i", e_vec, e_vec))
        inc = np.arccos(np.clip(h[:, 2] / hn, -1.0, 1.0))
        n = np.column_stack([-h[:, 1], h[:, 0], np.zeros(len(h))])  # z x h
        nn = np.sqrt(np.einsum("ij,ij->i", n, n))
        eps = 1e-12
        equatorial = nn < eps * hn
        n_hat = np.where(equatorial[:, None], np.array([1.0, 0.0, 0.0]), n / nn[:, None])
        Omega = np.where(equatorial, 0.0, np.arctan2(n[:, 1], n[:, 0]))
        h_hat = h / hn[:, None]
        circular = e < eps
        p_hat = np.where(circular[:, None], n_hat, e_vec / e[:, None])
        omega = np.where(
            circular,
            0.0,
            np.arctan2(
                np.einsum("ij,ij->i", cross3(n_hat, p_hat), h_hat),
                np.einsum("ij,ij->i", n_hat, p_hat),
            ),
        )
        nu = np.arctan2(
            np.einsum("ij,ij->i", cross3(p_hat, r), h_hat), np.einsum("ij,ij->i", p_hat, r)
        )
        E = 2.0 * np.arctan2(
            np.sqrt(1.0 - e) * np.sin(nu / 2.0), np.sqrt(1.0 + e) * np.cos(nu / 2.0)
        )
        M = E - e * np.sin(E)
    out = {
        "a": a,
        "e": e,
        "inc": inc,
        "Omega": np.mod(Omega, TWO_PI),
        "omega": np.mod(omega, TWO_PI),
        "M": np.mod(M, TWO_PI),
    }
    return {k: np.where(bound, x, np.nan) for k, x in out.items()}


def rotation_matrix(inc: float, Omega: float, omega: float) -> np.ndarray:
    cO, sO = np.cos(Omega), np.sin(Omega)
    ci, si = np.cos(inc), np.sin(inc)
    cw, sw = np.cos(omega), np.sin(omega)
    return np.array(
        [
            [cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
            [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
            [sw * si, cw * si, ci],
        ]
    )


def elements_to_state(mu: float, el: Elements) -> tuple[np.ndarray, np.ndarray]:
    """Relative position [AU] and velocity [AU/yr]; mu = G * (m1 + m2)."""
    if not 0.0 <= el.e < 1.0:
        raise ValueError(f"only bound orbits supported, got e={el.e}")
    E = solve_kepler(el.M, el.e)
    cE, sE = np.cos(E), np.sin(E)
    b = np.sqrt(1.0 - el.e**2)
    r_pf = np.array([el.a * (cE - el.e), el.a * b * sE, 0.0])
    r = el.a * (1.0 - el.e * cE)
    v_pf = np.sqrt(mu * el.a) / r * np.array([-sE, b * cE, 0.0])
    R = rotation_matrix(el.inc, el.Omega, el.omega)
    return R @ r_pf, R @ v_pf


def state_to_elements(mu: float, r: np.ndarray, v: np.ndarray) -> Elements:
    """Osculating elements of a bound relative orbit.

    Degenerate cases (inc = 0 or e = 0) return Omega = 0 or omega = 0 by convention.
    """
    r = np.asarray(r, float)
    v = np.asarray(v, float)
    rn = np.linalg.norm(r)
    h = np.cross(r, v)
    hn = np.linalg.norm(h)
    energy = 0.5 * v @ v - mu / rn
    if energy >= 0.0:
        raise ValueError("orbit is not bound")
    a = -mu / (2.0 * energy)
    e_vec = np.cross(v, h) / mu - r / rn
    e = np.linalg.norm(e_vec)
    inc = np.arccos(np.clip(h[2] / hn, -1.0, 1.0))

    n = np.array([-h[1], h[0], 0.0])  # z x h
    nn = np.linalg.norm(n)
    eps = 1e-12
    if nn < eps * hn:  # equatorial orbit
        n_hat = np.array([1.0, 0.0, 0.0])
        Omega = 0.0
    else:
        n_hat = n / nn
        Omega = np.arctan2(n[1], n[0])

    h_hat = h / hn
    if e < eps:  # circular: measure from the node
        omega = 0.0
        p_hat = n_hat
    else:
        p_hat = e_vec / e
        omega = np.arctan2(np.cross(n_hat, p_hat) @ h_hat, n_hat @ p_hat)

    nu = np.arctan2(np.cross(p_hat, r) @ h_hat, p_hat @ r)
    E = 2.0 * np.arctan2(np.sqrt(1.0 - e) * np.sin(nu / 2.0), np.sqrt(1.0 + e) * np.cos(nu / 2.0))
    M = E - e * np.sin(E)
    return Elements(
        a=a,
        e=e,
        inc=inc,
        Omega=np.mod(Omega, TWO_PI),
        omega=np.mod(omega, TWO_PI),
        M=np.mod(M, TWO_PI),
    )


def angular_momentum_unit(r: np.ndarray, v: np.ndarray) -> np.ndarray:
    h = np.cross(r, v)
    return h / np.linalg.norm(h)
