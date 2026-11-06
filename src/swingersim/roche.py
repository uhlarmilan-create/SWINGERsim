"""Roche geometry: distorted stars as surfaces of constant potential (tier A+ photometry,
contact binaries).

Two models:

* **Binary Roche model** (Kopal 1959; Wilson & Devinney 1971): a circular, synchronised binary
  with mass ratio q = M_2 / M_1, separation 1, star 1 at the origin and star 2 at (1, 0, 0),
  orbit normal along z. The dimensionless potential is

      Omega(x, y, z) = 1/r_1 + q (1/r_2 - x) + (1 + q)/2 (x^2 + y^2),

  larger inside the stars. The inner Lagrangian point L1 lies between the stars; the inner
  critical surface Omega_in = Omega(L1) bounds the Roche lobes, the outer critical surface
  Omega_out = max(Omega(L2), Omega(L3)) the common envelope of a contact binary. The fill-out
  factor is f = (Omega_in - Omega) / (Omega_in - Omega_out) (Lucy & Wilson 1979; Rucinski
  1973): f < 0 detached, 0 <= f <= 1 contact. A star's lobe is split from its companion's at
  the plane through L1 (x = x_L1).
* **Rotating star** (Roche model of a rigidly rotating star, the mass concentrated in the
  centre): in units of the polar radius, 1/r + (4/27) w^2 r^2 sin^2(theta) = 1 with
  w = Omega / Omega_crit, Omega_crit = sqrt(8 G M / (27 R_pol^3)) (Maeder 2009, Sect. 4.2); the
  equator reaches R_eq / R_pol = 1.5 at w = 1 (exact).

Every surface is sampled on a grid of directions from the centre of its star (equal-area rings
in cos theta, uniform in phi); the radius along each direction is found by bisection, the local
normal and effective gravity from the gradient of the potential. Volumes give the
volume-equivalent radius that the stellar models describe.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

# --------------------------------------------------------------------------- binary model


def potential(x, y, z, q: float):
    """Dimensionless binary Roche potential (star 1 at the origin, star 2 at x = 1)."""
    r1 = np.sqrt(x * x + y * y + z * z)
    r2 = np.sqrt((x - 1.0) ** 2 + y * y + z * z)
    return 1.0 / r1 + q * (1.0 / r2 - x) + 0.5 * (1.0 + q) * (x * x + y * y)


def gradient(x, y, z, q: float):
    """Gradient of the binary potential (components)."""
    r1 = np.sqrt(x * x + y * y + z * z)
    r2 = np.sqrt((x - 1.0) ** 2 + y * y + z * z)
    gx = -x / r1**3 - q * (x - 1.0) / r2**3 - q + (1.0 + q) * x
    gy = -y / r1**3 - q * y / r2**3 + (1.0 + q) * y
    gz = -z / r1**3 - q * z / r2**3
    return gx, gy, gz


def _dphi_dx_axis(x, q):
    """dOmega/dx on the x axis (y = z = 0)."""
    return -np.sign(x) / x**2 - q * np.sign(x - 1.0) / (x - 1.0) ** 2 - q + (1.0 + q) * x


def _root(f, lo, hi, n=200):
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        if np.sign(f(mid)) == np.sign(f(lo)):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


@lru_cache(maxsize=512)
def lagrange_points(q: float) -> tuple[float, float, float]:
    """x coordinates of L1 (0 < x < 1), L2 (x > 1) and L3 (x < 0)."""
    eps = 1e-9
    x1 = _root(lambda x: _dphi_dx_axis(x, q), eps, 1.0 - eps)
    x2 = _root(lambda x: _dphi_dx_axis(x, q), 1.0 + eps, 3.0)
    x3 = _root(lambda x: _dphi_dx_axis(x, q), -3.0, -eps)
    return x1, x2, x3


def critical_potentials(q: float) -> tuple[float, float]:
    """(Omega_in, Omega_out): inner critical potential (L1) and outer critical potential
    (the larger of L2 and L3)."""
    x1, x2, x3 = lagrange_points(q)
    o1 = float(potential(x1, 0.0, 0.0, q))
    o2 = float(potential(x2, 0.0, 0.0, q))
    o3 = float(potential(x3, 0.0, 0.0, q))
    return o1, max(o2, o3)


def fillout_potential(q: float, f: float) -> float:
    """Potential of the common envelope for the fill-out factor f."""
    oi, oo = critical_potentials(q)
    return oi - f * (oi - oo)


def _directions(n_theta: int, n_phi: int):
    """Equal-area grid of directions: (unit vectors (N, 3), solid angles (N,))."""
    mu = (np.arange(n_theta) + 0.5) / n_theta * 2.0 - 1.0  # cos theta
    phi = (np.arange(n_phi) + 0.5) / n_phi * 2.0 * math.pi
    mu_g, phi_g = np.meshgrid(mu, phi, indexing="ij")
    s = np.sqrt(1.0 - mu_g**2)
    d = np.stack([s * np.cos(phi_g), s * np.sin(phi_g), mu_g], axis=-1).reshape(-1, 3)
    w = np.full(len(d), 4.0 * math.pi / (n_theta * n_phi))
    return d, w


_GRID = np.linspace(0.0, 1.0, 301)[1:]  # sampling of each ray for the first crossing


RAY_CACHE_MAX_RAYS = 4096  # grids up to this size (the volume radii) are cached (8 MB each)


def _ray_potentials(q: float, star: int, n_theta: int, n_phi: int):
    """Potential on the sampling grid of every ray of a star (independent of the surface
    potential, so that the potentials of a root search reuse it): (dirs, solid angles, t_max,
    potentials (N, 300)). Cached for the small grids of volume_radius."""
    if n_theta * n_phi <= RAY_CACHE_MAX_RAYS:
        return _ray_potentials_cached(q, star, n_theta, n_phi)
    return _ray_potentials_cached.__wrapped__(q, star, n_theta, n_phi)


@lru_cache(maxsize=4)
def _ray_potentials_cached(q: float, star: int, n_theta: int, n_phi: int):
    x_l1 = lagrange_points(q)[0]
    dirs, dw = _directions(n_theta, n_phi)
    cx = 0.0 if star == 0 else 1.0
    dx = dirs[:, 0]
    # distance to the L1 plane along the direction (inf when moving away from it)
    with np.errstate(divide="ignore", invalid="ignore"):
        t_plane = (x_l1 - cx) / dx
    t_plane = np.where(t_plane > 0, t_plane, np.inf)
    t_max = np.minimum(t_plane, 1.5)
    t = _GRID[None, :] * t_max[:, None]
    pot = potential(cx + t * dirs[:, 0:1], t * dirs[:, 1:2], t * dirs[:, 2:3], q)  # (N, 300)
    for a in (dirs, dw, t_max, pot):
        a.setflags(write=False)
    return dirs, dw, t_max, pot


def _radii_binary(q: float, omega0: float, star: int, n_theta: int, n_phi: int):
    """Radius along each direction from the centre of star 1 (star = 0) or 2 (star = 1) to the
    surface Omega = omega0 (or to the plane x = x_L1, then flagged as not on the surface)."""
    dirs, dw, t_max, pot = _ray_potentials(q, star, n_theta, n_phi)
    cx = 0.0 if star == 0 else 1.0

    def phi(t):
        return potential(cx + t * dirs[:, 0], t * dirs[:, 1], t * dirs[:, 2], q) - omega0

    # first crossing below omega0 on a grid along each ray (the centrifugal term makes the
    # potential rise again far from the stars, so the bracket must be the first crossing)
    grid = _GRID
    below = pot < omega0
    surface = below.any(axis=1)  # the potential drops below omega0 before the plane
    first = np.argmax(below, axis=1)
    hi = np.where(surface, grid[first] * t_max, t_max)
    lo = np.where(surface & (first > 0), grid[np.maximum(first - 1, 0)] * t_max, 1e-6)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        inside = phi(mid) > 0.0
        lo = np.where(inside, mid, lo)
        hi = np.where(inside, hi, mid)
    r = np.where(surface, 0.5 * (lo + hi), t_max)
    return dirs, dw, r, surface


@dataclass
class Surface:
    """Sampled surface of one star in its body frame (binary separation or polar radius = 1)."""

    points: np.ndarray  # (N, 3) relative to the frame origin (star 1 or the star's centre)
    normals: np.ndarray  # (N, 3) outward unit normals
    areas: np.ndarray  # (N,) surface element areas
    gravity: np.ndarray  # (N,) effective gravity (relative units)
    volume: float
    centre: np.ndarray  # centre of the star in the frame

    @property
    def r_vol(self) -> float:
        return (3.0 * self.volume / (4.0 * math.pi)) ** (1.0 / 3.0)

    @property
    def area(self) -> float:
        return float(self.areas.sum())


def binary_surface(q: float, omega0: float, star: int, n_theta: int = 48, n_phi: int = 96):
    """Surface of star 1 (star = 0) or star 2 (star = 1) at the potential omega0, in the frame
    of star 1 (separation 1)."""
    dirs, dw, r, surf = _radii_binary(float(q), float(omega0), star, n_theta, n_phi)
    volume = float(np.sum(r**3 / 3.0 * dw))
    cx = 0.0 if star == 0 else 1.0
    p = np.column_stack([cx + r * dirs[:, 0], r * dirs[:, 1], r * dirs[:, 2]])[surf]
    gx, gy, gz = gradient(p[:, 0], p[:, 1], p[:, 2], q)
    g = np.sqrt(gx * gx + gy * gy + gz * gz)
    n = -np.column_stack([gx, gy, gz]) / g[:, None]  # the potential grows inwards
    d = dirs[surf]
    cosg = np.clip(np.sum(n * d, axis=1), 0.05, 1.0)
    areas = r[surf] ** 2 * dw[surf] / cosg
    return Surface(p, n, areas, g, volume, np.array([cx, 0.0, 0.0]))


@lru_cache(maxsize=4096)
def volume_radius(q: float, omega0: float, star: int) -> float:
    """Volume-equivalent radius (in units of the separation) of a star's part of the surface
    Omega = omega0."""
    return binary_surface(q, omega0, star, 40, 80).r_vol


def roche_lobe_volume_radius(q: float, star: int = 0) -> float:
    """Volume-equivalent radius of the Roche lobe (Omega_in), the quantity Eggleton (1983)
    fits: r_L(q_E) with q_E = M_star / M_companion."""
    return volume_radius(q, critical_potentials(q)[0], star)


def _guess_potential(q: float, r_vol: float, star: int) -> float:
    """Potential at the pole of a sphere of radius r_vol around the star (a starting value)."""
    if star == 0:
        return float(potential(0.0, 0.0, r_vol, q))
    return float(potential(1.0, 0.0, r_vol, q))


@lru_cache(maxsize=1024)
def potential_for_radius(q: float, r_vol: float, star: int, tol: float = 1e-13) -> float:
    """Potential whose surface around the given star has the volume-equivalent radius r_vol
    (detached: r_vol below the lobe; contact: above, up to the outer critical surface).

    The root of volume_radius(q, Omega) = r_vol is found by a safeguarded secant iteration
    (bisection whenever a step leaves the bracket), started from the potential at the pole of a
    sphere of radius r_vol and a Newton step with dr/dOmega ~ -r^2 (sphere), to |r - r_vol| <
    tol r_vol; volume_radius decreases monotonically with Omega. A radius beyond the outer
    critical surface is an error."""
    q, r_vol = float(q), float(r_vol)
    oo = critical_potentials(q)[1]
    lo, hi = oo, 50.0 + 2.0 / max(r_vol, 1e-3)  # V(lo) >= r_vol > V(hi)
    r_out = volume_radius(q, oo, star)
    if r_vol > r_out * (1.0 + 1e-12):
        raise ValueError(
            f"Roche star {star + 1} (q = {q:.4g}): volume radius {r_vol:.5g} (units of the "
            f"separation) exceeds the outer critical surface ({r_out:.5g})"
        )

    def f(om):
        return volume_radius(q, om, star) - r_vol

    x0 = min(max(_guess_potential(q, r_vol, star), lo), hi)
    f0 = f(x0)
    if f0 > 0.0:
        lo = x0
    else:
        hi = x0
    x1 = min(max(x0 + f0 / r_vol**2, lo), hi)
    for _ in range(100):
        if x1 in (lo, hi) or not lo < x1 < hi:
            x1 = 0.5 * (lo + hi)
        f1 = f(x1)
        if abs(f1) < tol * r_vol or hi - lo < 1e-14 * hi:
            return x1
        if f1 > 0.0:
            lo = x1
        else:
            hi = x1
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0) if f1 != f0 else 0.5 * (lo + hi)
        x0, f0 = x1, f1
        x1 = x2 if lo < x2 < hi else 0.5 * (lo + hi)
    return x1


def contact_radii(q: float, f: float) -> tuple[float, float]:
    """Volume-equivalent radii of the two stars of a contact binary with fill-out f (units of
    the separation)."""
    o = fillout_potential(q, f)
    return volume_radius(q, o, 0), volume_radius(q, o, 1)


def fillout_of(q: float, r1_vol: float) -> float:
    """Fill-out factor of a contact binary whose primary has the volume radius r1_vol."""
    oi, oo = critical_potentials(q)
    return (oi - potential_for_radius(q, r1_vol, 0)) / (oi - oo)


# --------------------------------------------------------------------------- rotating star


def rotating_radius(theta, w: float):
    """Radius (units of the polar radius) of a Roche-model rotating star at colatitude theta,
    w = Omega / Omega_crit: root of 1/r + (4/27) w^2 r^2 sin^2 theta = 1."""
    s2 = np.sin(theta) ** 2
    lo = np.ones_like(s2)
    hi = np.full_like(s2, 1.5)
    k = 4.0 / 27.0 * w * w
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        f = 1.0 / mid + k * mid * mid * s2 - 1.0
        lo = np.where(f > 0, mid, lo)
        hi = np.where(f > 0, hi, mid)
    return 0.5 * (lo + hi)


def rotating_surface(w: float, n_theta: int = 48, n_phi: int = 96) -> Surface:
    """Surface of a rotating star (polar radius 1, rotation axis z)."""
    dirs, dw = _directions(n_theta, n_phi)
    theta = np.arccos(np.clip(dirs[:, 2], -1.0, 1.0))
    r = rotating_radius(theta, w)
    p = dirs * r[:, None]
    k = 8.0 / 27.0 * w * w  # gradient of 1/r + (4/27) w^2 (x^2 + y^2)
    rr = np.linalg.norm(p, axis=1)
    gx = -p[:, 0] / rr**3 + k * p[:, 0]
    gy = -p[:, 1] / rr**3 + k * p[:, 1]
    gz = -p[:, 2] / rr**3
    g = np.sqrt(gx * gx + gy * gy + gz * gz)
    n = -np.column_stack([gx, gy, gz]) / g[:, None]
    cosg = np.clip(np.sum(n * dirs, axis=1), 0.05, 1.0)
    areas = r**2 * dw / cosg
    volume = float(np.sum(r**3 / 3.0 * dw))
    return Surface(p, n, areas, g, volume, np.zeros(3))


def omega_crit_fraction(p_rot_day: float, m_msun: float, r_vol_rsun: float) -> float:
    """w = Omega / Omega_crit of a star with rotation period p_rot_day, mass and
    volume-equivalent radius (iterated: the polar radius depends on w)."""
    from . import units

    G = units.GAUSS_K**2  # AU^3 / (Msun day^2)
    rsun_au = units.RSUN_AU
    om = 2.0 * math.pi / p_rot_day
    w = 0.0
    for _ in range(50):
        s = rotating_surface(w, 24, 48) if w > 0 else None
        r_pol = r_vol_rsun * rsun_au / (s.r_vol if s is not None else 1.0)
        w_new = om / math.sqrt(8.0 * G * m_msun / (27.0 * r_pol**3))
        if w_new >= 1.0:
            raise ValueError(
                f"rotation period {p_rot_day} d exceeds the break-up rate of the star "
                f"({m_msun} Msun, {r_vol_rsun} Rsun)"
            )
        if abs(w_new - w) < 1e-10:
            break
        w = w_new
    return w


__all__ = [
    "Surface",
    "binary_surface",
    "contact_radii",
    "critical_potentials",
    "fillout_of",
    "fillout_potential",
    "lagrange_points",
    "omega_crit_fraction",
    "potential",
    "potential_for_radius",
    "roche_lobe_volume_radius",
    "rotating_radius",
    "rotating_surface",
    "volume_radius",
]
