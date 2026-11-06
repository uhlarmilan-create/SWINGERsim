"""Dust discs and rings as occulters in the light curves (tier A photometry).

A disc is a slab: a cylinder of outer radius r_out, optional inner hole r_in and full thickness
H, centred on a star or on the centre of mass of a binary (or any subsystem). It has no mass;
it can shine and scatter starlight (disc_light) and can precess as a rigid
body torqued by the binary (dynamics: rigid); the stars behind it are dimmed by exp(-tau L / H), with L the
length of the line of sight inside the slab between the star and the observer (opaque if tau is
null). An edge-on slab is thus a dark band of height H, a tilted one an ellipse with a rim, and
the band gets darker towards the edge-on view (longer path). Examples: the eclipses of
epsilon Aur (a disc around the companion), EE Cep (a precessing disc, eclipse depths changing
from cycle to cycle) and KH 15D (a precessing circumbinary ring covering the binary orbit).

The stellar disk is sampled on a sunflower lattice of 4000 equal-area points with limb-darkening
weights; the occultation of a star by a sharp edge is resolved to ~1e-3 of its flux (mean
4e-4 for a straight edge at a random position).

Since 0.27 a disc can also

* precess as a rigid body under the torque of the binary it surrounds (``dynamics: rigid``):
  the disc normal l obeys dl/dt = -w [(1 - e_b^2)(l.l_b) l_b x l - 5 (l.e_b) e_b x l] (Zanazzi &
  Lai 2018, eq. 3), with w the test-particle rate (3/4)(M_1 M_2 / M^2)(a_b / r)^2 n(r) averaged
  over the disc's angular momentum (Sigma ~ r^-p); around an eccentric binary the disc
  circulates about l_b or librates about e_b (polar alignment; Farago & Laskar 2010), and a
  viscous damping time aligns it with the state it librates or circulates around (Martin &
  Lubow 2017). A circumstellar disc precesses about the orbit of its host's companion at
  w = -(3/4)(G M_c / a^3)(1 - e^2)^(-3/2) cos(delta) int(Sigma r^3) / int(Sigma Omega r^3)
  (Larwood et al. 1996; Bate et al. 2000);
* have a soft outer edge (opacity ~ 1 / (1 + exp((r - r_out) / w))), integrated along every
  line of sight with Gauss-Legendre quadrature;
* shine: thermal emission of a flat reprocessing disc heated by its host star(s),
  T = (2/3 pi)^(1/4) (R/r)^(3/4) T (Chiang & Goldreich 1997), and starlight scattered by its
  illuminated surface with single-scattering albedo omega and a Henyey-Greenstein phase function.

Sky frame of swingersim: x, y on the sky, z along the line of sight, away from the observer.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .kepler import rotation_matrix

N_POINTS = 4000


@dataclass
class Disc:
    name: str
    stars: list[int]  # stars whose centre of mass is the disc centre (one for a star host)
    r_out: float
    r_in: float
    half: float  # half thickness [AU]
    kappa: float  # 1/AU (inf = opaque)
    normal0: np.ndarray  # unit normal at t = 0
    axis: np.ndarray | None  # precession axis (unit) or None
    rate: float  # precession rate [rad/yr], signed
    edge: float | None = None  # soft-edge width [AU]
    rigid: dict | None = None  # circumbinary rigid-body dynamics (Zanazzi & Lai 2018)
    spec: object = None
    _sol: object = None
    _emission: dict | None = None  # thermal emission per (band, host radii and Teff)

    def normal(self, t_yr: float) -> np.ndarray:
        if self.rigid is not None:
            return self._rigid_normal(t_yr)
        if self.axis is None or self.rate == 0.0:
            return self.normal0
        a = self.axis
        th = self.rate * t_yr
        n = self.normal0
        # Rodrigues rotation of the normal about the precession axis
        return (
            n * math.cos(th)
            + np.cross(a, n) * math.sin(th)
            + a * float(a @ n) * (1.0 - math.cos(th))
        )

    def _rigid_normal(self, t_yr: float) -> np.ndarray:
        from scipy.integrate import solve_ivp

        rg = self.rigid
        sol = self._sol
        if sol is None or not (sol.t[0] <= t_yr <= sol.t[-1] or sol.t[-1] <= t_yr <= sol.t[0]):
            span = max(abs(t_yr) * 2.0, 10.0 * 2.0 * math.pi / max(rg["w"], 1e-30), 1.0)
            end = math.copysign(span, t_yr if t_yr != 0 else 1.0)
            sol = solve_ivp(
                lambda t, l: rigid_rhs(l, rg),
                (0.0, end),
                self.normal0,
                method="DOP853",
                rtol=1e-9,
                atol=1e-12,
                dense_output=True,
            )
            self._sol = sol
        n = np.asarray(sol.sol(t_yr), float)
        return n / np.linalg.norm(n)


def rigid_rhs(l, rg: dict) -> np.ndarray:
    """dl/dt of a rigid disc around a binary (unit vectors l, l_b, e-vector e_b): Zanazzi & Lai
    (2018), eq. 3, plus the optional viscous damping towards the nearest stationary state."""
    lb, eb, w, e2 = rg["lb"], rg["eb"], rg["w"], rg["e2"]
    out = -w * (
        (1.0 - e2) * float(l @ lb) * np.cross(lb, l) - 5.0 * float(l @ eb) * np.cross(eb, l)
    )
    td = rg.get("damping")
    if td:
        lam = (1.0 - e2) * float(l @ lb) ** 2 - 5.0 * float(l @ eb) ** 2
        if lam < 0.0 and e2 > 0.0:  # librating: polar alignment along +-e_b
            ehat = eb / math.sqrt(e2)
            target = ehat * math.copysign(1.0, float(l @ ehat))
        else:
            target = lb * math.copysign(1.0, float(l @ lb))
        out = out + (target - float(l @ target) * l) / td
    return out


def ring_average_rate(m1: float, m2: float, a_b: float, r_in: float, r_out: float, p: float):
    """Test-particle precession rate (3/4)(M_1 M_2 / M^2)(a_b / r)^2 n(r) [rad/yr] averaged
    over a disc from r_in to r_out with Sigma ~ r^-p, weighted by its angular momentum."""
    from . import units

    m = m1 + m2
    r = np.linspace(max(r_in, 1e-6 * r_out), r_out, 4001)
    n = np.sqrt(units.G * m / r**3)
    wt = r ** (-p) * r**3 * n  # Sigma r^3 Omega dr ~ dL
    wb = 0.75 * m1 * m2 / m**2 * (a_b / r) ** 2 * n
    return float(np.trapezoid(wt * wb, r) / np.trapezoid(wt, r))


def circumstellar_rate(
    m_host: float, m_comp: float, a: float, e: float, r_in: float, r_out: float, p: float, cos_d
) -> float:
    """Precession rate [rad/yr] of a rigid disc around a star perturbed by a companion on an
    orbit (a, e): -(3/4)(G M_c / a^3)(1 - e^2)^(-3/2) cos(delta) int(Sigma r^3 dr) /
    int(Sigma Omega r^3 dr) (Larwood et al. 1996; Bate et al. 2000), negative = retrograde."""
    from . import units

    r = np.linspace(max(r_in, 1e-6 * r_out), r_out, 4001)
    om = np.sqrt(units.G * m_host / r**3)
    i3 = np.trapezoid(r ** (3.0 - p), r)
    io = np.trapezoid(r ** (3.0 - p) * om, r)
    return float(-0.75 * units.G * m_comp / a**3 * (1.0 - e * e) ** -1.5 * cos_d * i3 / io)


def _unit_normal(i_deg: float, Omega_deg: float) -> np.ndarray:
    return rotation_matrix(math.radians(i_deg), math.radians(Omega_deg), 0.0)[:, 2]


def build_discs(cfg, root) -> list[Disc]:
    """The discs of a configuration (empty list without discs)."""
    from .hierarchy import _stars

    if not getattr(cfg, "discs", None):
        return []
    nodes = {n.name: n for n in root.nodes()}
    out = []
    for name, d in cfg.discs.items():
        if d.host in nodes:
            stars = [int(k) for k in _stars(nodes[d.host])]
        else:
            stars = [int(d.host[1:]) - 1]
        axis, rate = None, 0.0
        if d.relative_to is not None:
            o = cfg.orbits[d.relative_to]
            R = rotation_matrix(math.radians(o.i), math.radians(o.Omega), 0.0)
            # tilt the orbit normal by `tilt` about the line of nodes at `tilt_node` in the
            # orbit plane (orbit frame: x = ascending node, z = orbit normal)
            ph, tl = math.radians(d.tilt_node), math.radians(d.tilt)
            # rotating z by tl about the node direction (cos ph, sin ph, 0)
            perp = np.array([-math.sin(ph), math.cos(ph), 0.0])
            n_orb = -math.sin(tl) * perp + math.cos(tl) * np.array([0.0, 0.0, 1.0])
            normal0 = R @ n_orb
            axis = R[:, 2]
            if d.precession_period_yr is not None:
                rate = 2.0 * math.pi / d.precession_period_yr
        else:
            normal0 = _unit_normal(d.i, d.Omega)
        normal0 = normal0 / np.linalg.norm(normal0)
        rigid = None
        if d.dynamics == "rigid":
            rigid, axis, rate = _rigid_setup(cfg, d, nodes, stars, root, normal0)
        if d.tau is None:
            kappa = math.inf if d.edge_au is None else 50.0 / d.thickness_au
        else:
            kappa = d.tau / d.thickness_au
        out.append(
            Disc(
                name,
                stars,
                d.r_out_au,
                d.r_in_au,
                0.5 * d.thickness_au,
                kappa,
                normal0,
                axis,
                rate,
                d.edge_au,
                rigid,
                d,
            )
        )
    return out


def _rigid_setup(cfg, d, nodes, stars, root, normal0):
    """(rigid parameters or None, precession axis, rate) of a disc with dynamics: rigid."""
    from .hierarchy import build_tree, initialize, mass_of, node_elements

    _, masses = build_tree(cfg)
    pos, vel = initialize(root, masses)
    o = nodes[d.relative_to]
    el = node_elements(o, pos, vel, masses)
    o_stars = sorted(o.stars())
    m1, m2 = mass_of(o.left, masses), mass_of(o.right, masses)
    R = rotation_matrix(
        math.radians(cfg.orbits[o.name].i), math.radians(cfg.orbits[o.name].Omega), 0.0
    )
    lb = R[:, 2]
    if sorted(stars) == o_stars:  # circumbinary disc around the orbit's stars
        Rw = rotation_matrix(
            math.radians(cfg.orbits[o.name].i),
            math.radians(cfg.orbits[o.name].Omega),
            math.radians(cfg.orbits[o.name].omega),
        )
        e = float(el.e)
        w = ring_average_rate(m1, m2, float(el.a), d.r_in_au, d.r_out_au, d.sigma_index)
        rg = {"lb": lb, "eb": e * Rw[:, 0], "w": w, "e2": e * e, "damping": d.damping_yr}
        return rg, None, 0.0
    # circumstellar: the host and its companion in orbit relative_to
    host_child = (
        o.left
        if stars[0] in (o.left.stars() if not isinstance(o.left, int) else [o.left])
        else o.right
    )
    m_host = mass_of(host_child, masses)
    m_comp = (m1 + m2) - m_host
    cos_d = float(normal0 @ lb)
    rate = circumstellar_rate(
        m_host, m_comp, float(el.a), float(el.e), d.r_in_au, d.r_out_au, d.sigma_index, cos_d
    )
    return None, lb, rate


def _intervals_le(A, B, C, lim2):
    """t with A t^2 + B t + C <= lim2 (A >= 0), as (lo, hi) arrays; empty: lo > hi."""
    c = C - lim2
    lo = np.full_like(C, np.inf)
    hi = np.full_like(C, -np.inf)
    quad = A > 1e-14
    if np.any(quad):
        disc = B * B - 4.0 * A * c
        ok = quad & (disc >= 0)
        sq = np.sqrt(np.where(ok, disc, 0.0))
        a2 = np.where(quad, 2.0 * A, 1.0)
        lo = np.where(ok, (-B - sq) / a2, lo)
        hi = np.where(ok, (-B + sq) / a2, hi)
    # A ~ 0 (line of sight along the disc axis): the radius is constant along the line
    lin = ~quad
    inside = lin & (c <= 0)
    lo = np.where(inside, -np.inf, lo)
    hi = np.where(inside, np.inf, hi)
    return lo, hi


def path_length(x, y, z_max, normal, r_out, r_in, half):
    """Length of the lines of sight (x, y), running from z = -inf (observer) to z_max, inside
    the slab |u.n| <= half, r_in <= |u - (u.n) n| <= r_out (disc centre at the origin)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    nx, ny, nz = normal
    a0 = x * nx + y * ny  # axial coordinate at z = 0; a(t) = a0 + nz t
    # axial interval
    if abs(nz) > 1e-14:
        t1 = (-half - a0) / nz
        t2 = (half - a0) / nz
        ax_lo, ax_hi = np.minimum(t1, t2), np.maximum(t1, t2)
    else:
        ok = np.abs(a0) <= half
        ax_lo = np.where(ok, -np.inf, np.inf)
        ax_hi = np.where(ok, np.inf, -np.inf)
    ax_hi = np.minimum(ax_hi, z_max)
    # radial distance^2 along the line: A t^2 + B t + C
    A = np.full_like(x, 1.0 - nz * nz)
    B = -2.0 * a0 * nz
    C = x * x + y * y - a0 * a0
    o_lo, o_hi = _intervals_le(A, B, C, r_out * r_out)
    lo = np.maximum(ax_lo, o_lo)
    hi = np.minimum(ax_hi, o_hi)
    length = np.clip(hi - lo, 0.0, None)
    length = np.where(np.isfinite(length), length, 0.0)
    if r_in > 0:
        h_lo, h_hi = _intervals_le(A, B, C, r_in * r_in)
        cut_lo = np.maximum(lo, h_lo)
        cut_hi = np.minimum(hi, h_hi)
        cut = np.clip(cut_hi - cut_lo, 0.0, None)
        cut = np.where(np.isfinite(cut), cut, 0.0)
        length = np.where(length > 0, np.clip(length - cut, 0.0, None), 0.0)
    return length


def _chord(x, y, z_max, normal, r_out, half):
    """(lo, hi) of the lines of sight inside the slab |u.n| <= half, |u - (u.n) n| <= r_out,
    up to z_max (empty: hi <= lo)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    nx, ny, nz = normal
    a0 = x * nx + y * ny
    if abs(nz) > 1e-14:
        t1 = (-half - a0) / nz
        t2 = (half - a0) / nz
        ax_lo, ax_hi = np.minimum(t1, t2), np.maximum(t1, t2)
    else:
        ok = np.abs(a0) <= half
        ax_lo = np.where(ok, -1e30, 1e30)
        ax_hi = np.where(ok, 1e30, -1e30)
    ax_hi = np.minimum(ax_hi, z_max)
    A = np.full_like(x, 1.0 - nz * nz)
    B = -2.0 * a0 * nz
    C = x * x + y * y - a0 * a0
    o_lo, o_hi = _intervals_le(A, B, C, r_out * r_out)
    lo = np.maximum(ax_lo, np.maximum(o_lo, -1e30))
    hi = np.minimum(ax_hi, np.minimum(o_hi, 1e30))
    return lo, hi, a0, nz


_GL16 = np.polynomial.legendre.leggauss(16)


def soft_optical_depth(x, y, z_max, normal, d: Disc) -> np.ndarray:
    """Optical depth along the lines of sight through a disc with a soft outer edge (and its
    inner hole), Gauss-Legendre quadrature (16 points) along the chord."""
    w_edge = d.edge
    r_ext = d.r_out + 8.0 * w_edge
    lo, hi, a0, nz = _chord(x, y, z_max, normal, r_ext, d.half)
    ok = hi > lo
    tau = np.zeros(np.shape(lo))
    if not np.any(ok):
        return tau
    xg, wg = _GL16
    lo_, hi_ = lo[ok], hi[ok]
    t = 0.5 * (hi_ - lo_)[:, None] * (xg[None, :] + 1.0) + lo_[:, None]
    xx = np.asarray(x, float)[ok][:, None]
    yy = np.asarray(y, float)[ok][:, None]
    a = a0[ok][:, None] + nz * t  # axial coordinate along the line
    r2 = xx * xx + yy * yy + t * t - a * a
    r = np.sqrt(np.clip(r2, 0.0, None))
    k = d.kappa / (1.0 + np.exp(np.clip((r - d.r_out) / w_edge, -50.0, 50.0)))
    k = np.where(r < d.r_in, 0.0, k)
    tau[ok] = 0.5 * (hi_ - lo_) * (k * wg[None, :]).sum(axis=1)
    return tau


def _stellar_grid(ld_profile, a):
    """Points on the unit stellar disk and their limb-darkened weights (sum 1): a sunflower
    (Fibonacci) lattice of equal-area cells, which has no preferred direction, so a straight
    edge crossing the star is resolved without the bias of a ring x sector grid."""
    k = np.arange(N_POINTS) + 0.5
    s = np.sqrt(k / N_POINTS)
    phi = k * math.pi * (3.0 - math.sqrt(5.0))
    mu = np.sqrt(np.clip(1.0 - s * s, 0.0, 1.0))
    w = ld_profile(a, mu)
    return s * np.cos(phi), s * np.sin(phi), w / w.sum()


def transmission(discs, t_yr, pos, masses, lights, cache: dict | None = None) -> np.ndarray:
    """Fraction of the flux of every star that passes the discs (1 = not occulted)."""
    from .photometry import ld_profile

    n = len(lights)
    out = np.ones(n)
    if not discs:
        return out
    cache = {} if cache is None else cache
    for d in discs:
        m = masses[d.stars]
        c = (m[:, None] * pos[d.stars]).sum(axis=0) / m.sum()
        nrm = d.normal(t_yr)
        reach = d.r_out + d.half + (8.0 * d.edge if d.edge is not None else 0.0)
        for s in range(n):
            R = lights[s].radius_au
            dx, dy = pos[s, 0] - c[0], pos[s, 1] - c[1]
            if math.hypot(dx, dy) > reach + R:
                continue
            if s not in cache:
                cache[s] = _stellar_grid(ld_profile, lights[s].ld)
            px, py, w = cache[s]
            if d.edge is not None:
                tau = soft_optical_depth(dx + R * px, dy + R * py, pos[s, 2] - c[2], nrm, d)
                if np.any(tau > 0):
                    out[s] *= float(w @ np.exp(-tau))
                continue
            L = path_length(
                dx + R * px, dy + R * py, pos[s, 2] - c[2], nrm, d.r_out, d.r_in, d.half
            )
            if not np.any(L > 0):
                continue
            if math.isinf(d.kappa):
                tr = np.where(L > 0, 0.0, 1.0)
            else:
                tr = np.exp(-d.kappa * L)
            out[s] *= float(w @ tr)
    return out


def _hg(cos_t, g):
    """Henyey-Greenstein phase function (normalised over 4 pi)."""
    return (1.0 - g * g) / (4.0 * math.pi * (1.0 + g * g - 2.0 * g * cos_t) ** 1.5)


def _thermal_emission(d: Disc, r, dA, lights, teffs, band) -> float:
    """Sum of the passband intensity x area of the rings of a flat reprocessing disc, T(r) =
    (2/3 pi)^(1/4) (R/r)^(3/4) Teff summed over its host stars (Chiang & Goldreich 1997). It
    depends only on the rings and the hosts' radii and Teff: cached on the disc per (band,
    radii, Teff), since the 40 passband integrals cost ~1 ms per time point."""
    from .photometry import surface_brightness

    key = (band, *((float(lights[k].radius_au), float(teffs[k])) for k in d.stars))
    cache = d._emission
    if cache is None:
        cache = d._emission = {}
    v = cache.get(key)
    if v is None:
        t4 = np.zeros(len(r))
        for k in d.stars:
            rk = lights[k].radius_au
            t4 += (2.0 / (3.0 * math.pi)) * (rk / r) ** 3 * teffs[k] ** 4
        temp = t4**0.25
        w = np.array([surface_brightness(round(float(tt), 1), band) for tt in temp])
        v = float(np.sum(w * dA))
        if len(cache) >= 256:  # evolving hosts: a new key at every time point
            cache.clear()
        cache[key] = v
    return v


def disc_light(discs, t_yr, pos, masses, lights, teffs, band) -> np.ndarray:
    """Light of every disc (thermal emission and scattered starlight), in the units of the
    stellar fluxes of swingersim.photometry (pi R^2 w for a star outside eclipse)."""
    out = np.zeros(len(discs))
    nr, nphi = 40, 72
    for q, d in enumerate(discs):
        sp = d.spec
        if sp is None or (not sp.emission and sp.albedo <= 0.0):
            continue
        m = masses[d.stars]
        c = (m[:, None] * pos[d.stars]).sum(axis=0) / m.sum()
        nrm = d.normal(t_yr)
        cos_i = abs(float(nrm[2]))
        edges = np.geomspace(max(d.r_in, 1e-3 * d.r_out), d.r_out, nr + 1)
        r = np.sqrt(edges[1:] * edges[:-1])
        dA = math.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
        tau_perp = math.inf if math.isinf(d.kappa) else d.kappa * 2.0 * d.half
        if sp.emission:
            seen = 1.0 if math.isinf(tau_perp) else 1.0 - math.exp(-tau_perp / max(cos_i, 1e-3))
            out[q] += _thermal_emission(d, r, dA, lights, teffs, band) * cos_i * seen
        if sp.albedo > 0.0:
            # the face turned to the observer, sampled in rings and azimuths
            face = nrm if nrm[2] < 0 else -nrm  # its outward normal points to the observer
            e1 = np.cross(nrm, [0.0, 0.0, 1.0])
            e1 = e1 / np.linalg.norm(e1) if np.linalg.norm(e1) > 1e-9 else np.array([1.0, 0, 0])
            e2 = np.cross(nrm, e1)
            phi = (np.arange(nphi) + 0.5) / nphi * 2.0 * math.pi
            pts = (
                c[None, None, :]
                + r[:, None, None]
                * (np.cos(phi)[None, :, None] * e1 + np.sin(phi)[None, :, None] * e2)
                + d.half * face[None, None, :]
            )
            dA2 = (dA / nphi)[:, None]
            graze = np.clip(d.half / r, 0.02, 1.0)[:, None]  # flat slab lit at grazing incidence
            for k in range(len(lights)):
                v = pos[k][None, None, :] - pts  # element -> star
                dist = np.linalg.norm(v, axis=2)
                u = v / dist[..., None]
                cos_inc = np.einsum("ijk,k->ij", u, face)
                inc = np.where(cos_inc > 0, np.maximum(cos_inc, graze), graze)
                if k not in d.stars:
                    inc = np.where(cos_inc > 0, cos_inc, 0.0)  # external star: lit side only
                f_in = math.pi * lights[k].radius_au ** 2 * lights[k].weight / dist**2
                cos_t = u[..., 2]  # scattering angle: star -> element against element -> observer
                out[q] += float(np.sum(f_in * inc * dA2 * sp.albedo * _hg(cos_t, sp.g)))
    return out


__all__ = [
    "Disc",
    "build_discs",
    "circumstellar_rate",
    "disc_light",
    "path_length",
    "ring_average_rate",
    "soft_optical_depth",
    "transmission",
]
