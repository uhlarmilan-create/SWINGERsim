"""Light of distorted stars on surface meshes (tier A+ photometry; swingersim.roche geometry).

Used by the light curves when a star has ``stars.<name>.shape: rotating`` (an oblate rigidly
rotating star) or ``roche`` (a star of a two-star binary filling its Roche equipotential), or an
orbit has ``contact_fillout`` (a contact binary in a common envelope; Lucy 1968). Then every star
of the system is a mesh of surface elements (spheres too, so that all occultations are treated
alike):

* the surface of each body (a sphere, a rotating star, the two stars of a Roche binary) is
  computed once in its body frame (swingersim.roche) and placed on the sky at every instant: a Roche
  binary in the frame of its two stars (x along the separation, z along the orbit normal,
  scale = their instantaneous separation; the model assumes a circular, synchronised orbit;
  detached lobes keep their size in a slightly eccentric or perturbed orbit, each scaled about
  its centre by a / separation, as the constant-volume surfaces of Wilson 1979), a
  rotating star with its spin axis (stars.<name>.obliquity, spin_azimuth against its orbit),
* an element emits I = I_ref (g / <g>)^y L(mu): gravity darkening with the band exponent y of the
  star (von Zeipel 1924; Claret & Bloemen 2011), the four-coefficient limb darkening L(mu) and an
  area-weighted mean intensity equal to the star's surface brightness; the stars of a contact
  binary share one effective temperature, sigma T^4 (A_1 + A_2) = L_1 + L_2 (Lucy 1968),
* an element is seen when its normal points to the observer and no body lies on its line of
  sight: every element inside the sky circle of a body nearer to the observer is tested along
  its ray (48 points) against the inside of that body's equipotential - this includes the
  self-occultation of a contact binary's two lobes;
* reflection (Lambert spheres of the volume radii; not inside a contact envelope) and Doppler
  beaming as in swingersim.photometry (photometry.proximity_terms); the ellipsoidal variation of
  a Roche body comes from its geometry, the other pairs (spheres, rotating stars) get the
  analytic term of swingersim.photometry.

Speed: the Roche geometry (surfaces, intensities, phase tables) is cached per process and
reused while (q, Omega_1, Omega_2) and the resolution are the same; the phase tables lie on a
fixed grid of inclinations (TABLE_DI_DEG) and are interpolated linearly in between.

Accuracy: ~3e-4 of the flux with the default 64 x 128 directions (benchmark
``mesh-sphere-eclipse`` against the exact sphere eclipses; the edges of spheres in front are
anti-aliased, other bodies are tested along each ray with a golden-section refinement).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from . import roche, units

N_THETA, N_PHI = 64, 128
N_RAY = 48
# fast mode: unocculted Roche bodies from a table over the orbital phase, unocculted spheres
# analytically; the full mesh only where bodies overlap on the sky (eclipses between bodies)
FAST = True
TABLE_PHASES = 256
# phase tables on a fixed grid of inclinations [deg], interpolated linearly in between
TABLE_DI_DEG = 0.05
TABLE_MAX_N = 16  # tables kept per Roche body
# FAST mode: a Roche pair whose separation differs from the separation of its phase table by more
# than this (relative) is computed in full while its stars can eclipse each other - the table has
# the eclipse geometry of one separation (eccentric orbits: 1.2 % flux error at e = 0.1 otherwise)
TABLE_MAX_DK = 1e-4
SELF_ECLIPSE_MARGIN = 1.25  # sky distance below (R1 + R2) * this: an eclipse may be in progress
_SPHERE_FRAME = np.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]])


@dataclass
class Body:
    kind: str  # sphere, rotating, roche
    stars: list[int]
    surfaces: list  # roche.Surface per star of the body (body units)
    q: float = 1.0  # roche
    omegas: tuple[float, float] = (0.0, 0.0)  # roche: potentials of the two stars
    w: float = 0.0  # rotating: Omega / Omega_crit
    scale: float = 1.0  # sphere, rotating: AU per body unit (radius or polar radius)
    axis: np.ndarray | None = None  # rotating: spin axis on the sky
    extent: float = 1.0  # radius of a sphere around the body centre enclosing it (body units)
    contact: bool = False  # roche: a common envelope (contact_fillout)
    a_ref: float = 0.0  # roche, detached: separation [AU] at which the lobes have their radii
    k: float = 1.0  # roche: lobe size factor a_ref / separation of the last frame
    tables: dict | None = None  # roche, fast mode: inclination [deg] -> phase spline of fluxes

    def frame(self, pos, vel, masses):
        """(origin [AU], 3x3 matrix body -> sky, AU per body unit)."""
        if self.kind == "roche":
            i, j = self.stars
            r = pos[j] - pos[i]
            d = float(np.linalg.norm(r))
            ex = r / d
            h = np.cross(r, vel[j] - vel[i])
            ez = h / np.linalg.norm(h)
            ey = np.cross(ez, ex)
            # detached stars keep their size when the separation changes (an eccentric or
            # perturbed orbit): each lobe is scaled about its own centre by a_ref / d
            self.k = self.a_ref / d if self.a_ref > 0.0 and not self.contact else 1.0
            return pos[i], np.column_stack([ex, ey, ez]), d
        k = self.stars[0]
        if self.kind == "rotating":
            ez = self.axis
            tmp = np.array([1.0, 0.0, 0.0]) if abs(ez[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
            ex = np.cross(tmp, ez)
            ex /= np.linalg.norm(ex)
            ey = np.cross(ez, ex)
            return pos[k], np.column_stack([ex, ey, ez]), self.scale
        # a sphere: the grid's pole in the sky plane, so that its rings are not concentric with
        # the stellar disk (no aliasing with the edges of occulters)
        return pos[k], _SPHERE_FRAME, self.scale

    def centre_offset(self) -> np.ndarray:
        return np.array([0.5, 0.0, 0.0]) if self.kind == "roche" else np.zeros(3)

    def margin(self, p: np.ndarray) -> np.ndarray:
        """Continuous inside-measure of points (body units, (..., 3)): > 0 inside the body."""
        x, y, z = p[..., 0], p[..., 1], p[..., 2]
        if self.kind == "sphere":
            return 1.0 - (x * x + y * y + z * z)
        if self.kind == "rotating":
            r = np.sqrt(x * x + y * y + z * z)
            with np.errstate(divide="ignore"):
                f = 1.0 / np.maximum(r, 1e-12) + 4.0 / 27.0 * self.w**2 * (x * x + y * y)
            return f - 1.0
        x_l1 = roche.lagrange_points(self.q)[0]
        if self.k != 1.0:  # back to the lobes at their own size (see frame)
            cx = np.where(x < x_l1, 0.0, 1.0)
            x, y, z = cx + (x - cx) / self.k, y / self.k, z / self.k
        with np.errstate(divide="ignore"):
            pot = roche.potential(x, y, z, self.q)
        om = np.where(x < x_l1, self.omegas[0], self.omegas[1])
        return pot - om

    def inside(self, p: np.ndarray) -> np.ndarray:
        """Points (body units, (..., 3)) inside the body."""
        return self.margin(p) > 0.0


@dataclass
class MeshStar:
    star: int
    body: Body
    surface: roche.Surface
    rel_intensity: np.ndarray  # (g/<g>)^y normalised to an area-weighted mean of 1
    unit_area: float  # sum of the areas (body units)


def _rel_intensity(s: roche.Surface, y: float) -> np.ndarray:
    gbar = float(np.sum(s.gravity * s.areas) / np.sum(s.areas))
    rel = (s.gravity / gbar) ** y
    return rel * np.sum(s.areas) / np.sum(rel * s.areas)


def spin_axis(cfg, root, k: int, pos, vel, masses) -> np.ndarray:
    """Spin axis of star k on the sky: the normal of its orbit tilted by the obliquity."""
    from .backends.nbody import parent_of
    from .hierarchy import node_h_unit

    node = parent_of(root, k)
    h = node_h_unit(node, pos, vel, masses)
    spec = cfg.stars.get(f"m{k + 1}")
    nvec = np.cross([0.0, 0.0, 1.0], h)
    nvec = nvec / np.linalg.norm(nvec) if np.linalg.norm(nvec) > 1e-12 else np.array([1.0, 0, 0])
    ob = math.radians(spec.obliquity if spec else 0.0)
    az = math.radians(spec.spin_azimuth if spec else 0.0)
    return math.cos(ob) * h + math.sin(ob) * (
        math.cos(az) * nvec + math.sin(az) * np.cross(h, nvec)
    )


def needed(cfg) -> bool:
    """The configuration asks for distorted stars in the light curves."""
    if any(getattr(s, "shape", "sphere") != "sphere" for s in cfg.stars.values()):
        return True
    return any(getattr(o, "contact_fillout", None) is not None for o in cfg.orbits.values())


# per-process cache of the Roche geometry: surfaces, enclosing radius, relative intensities and
# phase tables of a (q, Omega_1, Omega_2, resolution); a fit whose Jacobian columns leave the
# Roche parameters unchanged (inclination, periods, ...) reuses them
SCENE_CACHE_SIZE = 8
_SCENE_CACHE: dict = {}


def _r12(x: float) -> float:
    return float(f"{float(x):.12g}")


def _roche_geometry(q: float, om: tuple[float, float]) -> dict:
    key = (_r12(q), _r12(om[0]), _r12(om[1]), N_THETA, N_PHI)
    geo = _SCENE_CACHE.pop(key, None)
    if geo is None:
        surfs = [roche.binary_surface(q, om[s], s, N_THETA, N_PHI) for s in (0, 1)]
        ext = max(float(np.max(np.linalg.norm(sf.points - [0.5, 0, 0], axis=1))) for sf in surfs)
        geo = {"surfs": surfs, "extent": ext, "rel": {}, "tables": {}}
        while len(_SCENE_CACHE) >= SCENE_CACHE_SIZE:
            _SCENE_CACHE.pop(next(iter(_SCENE_CACHE)))
    _SCENE_CACHE[key] = geo  # most recently used last
    return geo


def _cached_rel(geo: dict, s: int, y: float) -> np.ndarray:
    key = (s, _r12(y))
    if key not in geo["rel"]:
        geo["rel"][key] = _rel_intensity(geo["surfs"][s], y)
    return geo["rel"][key]


def _cached_tables(geo: dict, lights2) -> dict:
    """Phase-table store of a geometry for the intensities, limb and gravity darkening of its
    two stars and the resolution of the tables (they depend on nothing else)."""
    key: tuple = (TABLE_PHASES, TABLE_DI_DEG, N_RAY)
    for L in lights2:
        key += (_r12(L.weight), tuple(_r12(v) for v in np.atleast_1d(L.ld)), _r12(L.gd_y))
    tabs = geo["tables"]
    if key not in tabs:
        if len(tabs) >= 4:
            tabs.pop(next(iter(tabs)))
        tabs[key] = {}
    return tabs[key]


def _detached_potential(q: float, r_vol: float, s: int, k: int, orbit: str) -> float:
    """Potential of a detached Roche star of volume radius r_vol (units of the separation); a
    star larger than its Roche lobe is an error (a contact binary needs contact_fillout)."""
    r_vol = float(f"{r_vol:.13g}")  # float noise must not defeat the caches
    o_in = roche.critical_potentials(q)[0]
    r_lobe = roche.volume_radius(q, o_in, s)
    if r_vol >= r_lobe:
        raise ValueError(
            f"orbit {orbit}: star m{k + 1} (shape roche) has the volume radius {r_vol:.5f} a, "
            f"larger than its Roche lobe ({r_lobe:.5f} a, q = {q:.4f}); a detached Roche pair "
            "needs both stars inside their lobes (contact binaries: orbits.<name>."
            "contact_fillout)"
        )
    return roche.potential_for_radius(q, r_vol, s)


def build_scene(cfg, root, masses, props, lights, pos, vel, band: str):
    """Bodies and mesh stars of the system at t = 0 (pos, vel: barycentric state)."""
    from .backends.nbody import parent_of
    from .hierarchy import node_elements

    n = len(masses)
    used = set()
    bodies, stars = [], {}
    for node in root.nodes():
        if not (isinstance(node.left, int) and isinstance(node.right, int)):
            continue
        i, j = node.left, node.right
        spec = cfg.orbits.get(node.name)
        fill = getattr(spec, "contact_fillout", None) if spec else None
        shapes = [getattr(cfg.stars.get(f"m{k + 1}"), "shape", "sphere") for k in (i, j)]
        if fill is None and "roche" not in shapes:
            continue
        q = float(masses[j] / masses[i])
        a = float(node_elements(node, pos, vel, masses).a)
        if fill is not None:
            o = roche.fillout_potential(q, float(fill))
            om = (o, o)
        else:
            if any(shapes[s] != "roche" for s in (0, 1)):
                raise ValueError(
                    f"orbit {node.name}: shape roche must be given to both of its stars"
                )
            om = tuple(
                _detached_potential(q, props[k].radius_rsun * units.RSUN_AU / a, s, k, node.name)
                for s, k in enumerate((i, j))
            )
        geo = _roche_geometry(q, om)
        surfs = geo["surfs"]
        b = Body(
            "roche",
            [i, j],
            surfs,
            q=q,
            omegas=om,
            extent=geo["extent"],
            contact=fill is not None,
            a_ref=0.0 if fill is not None else a,
        )
        bodies.append(b)
        for s, k in enumerate((i, j)):
            stars[k] = MeshStar(k, b, surfs[s], _cached_rel(geo, s, lights[k].gd_y), surfs[s].area)
            used.add(k)
        if fill is not None:  # one envelope, one temperature (Lucy 1968)
            from .photometry import surface_brightness

            lum = sum(props[k].luminosity_lsun for k in (i, j))
            area_rsun2 = sum(surfs[s].area * (a / units.RSUN_AU) ** 2 for s in (0, 1))
            t_c = 5772.0 * (lum / (area_rsun2 / (4.0 * math.pi))) ** 0.25
            for k in (i, j):
                lights[k].weight = surface_brightness(round(t_c, 1), band)
            for k, s in ((i, 0), (j, 1)):  # volume radii of the envelope's lobes
                lights[k].radius_au = surfs[s].r_vol * a
        # phase tables of this geometry and these intensities (shared between calls)
        b.tables = _cached_tables(geo, [lights[k] for k in (i, j)])
    for k in range(n):
        if k in used:
            continue
        shape = getattr(cfg.stars.get(f"m{k + 1}"), "shape", "sphere")
        r_au = lights[k].radius_au
        if shape == "rotating":
            spec = cfg.stars.get(f"m{k + 1}")
            node = parent_of(root, k)
            if spec is None or spec.P_rot in ("sync", "pseudosync"):
                el = node_elements(node, pos, vel, masses)
                p_rot = units.sma_to_period(el.a, float(sum(masses[node.stars()])))
                p_rot *= units.DAYS_PER_YEAR
                if spec is not None and spec.P_rot == "pseudosync":
                    from .tides import pseudo_sync_factor

                    p_rot /= pseudo_sync_factor(float(el.e))
            else:
                p_rot = float(spec.P_rot)
            w = roche.omega_crit_fraction(p_rot, float(masses[k]), props[k].radius_rsun)
            sf = roche.rotating_surface(w, N_THETA, N_PHI)
            b = Body(
                "rotating",
                [k],
                [sf],
                w=w,
                scale=r_au / sf.r_vol,
                axis=spin_axis(cfg, root, k, pos, vel, masses),
                extent=float(np.max(np.linalg.norm(sf.points, axis=1))),
            )
        else:
            sf = roche.rotating_surface(0.0, N_THETA, N_PHI)
            b = Body("sphere", [k], [sf], scale=r_au, extent=1.0)
        bodies.append(b)
        stars[k] = MeshStar(k, b, sf, _rel_intensity(sf, lights[k].gd_y), sf.area)
    return bodies, stars


def _frames(bodies, pos, vel, masses) -> dict:
    """Frames of all bodies at one instant (Body.frame; sets the lobe size factors b.k)."""
    return {id(b): b.frame(pos, vel, masses) for b in bodies}


def mesh_fluxes(bodies, stars, pos, vel, masses, lights, frames=None, elements=None) -> np.ndarray:
    """Flux of every star (same units as photometry.flux_at before normalisation: pi R^2 w for
    an unocculted sphere). ``frames``: the frames of this instant when already computed;
    ``elements``: a dict that receives, per star, (sky positions [AU], flux of every seen
    element, its flux without occultations) of its visible elements (centre-of-light
    velocities, swingersim.rvlight)."""
    from .photometry import ld_mean, ld_profile

    n = len(lights)
    out = np.zeros(n)
    if frames is None:
        frames = _frames(bodies, pos, vel, masses)
    # sky circle and nearest depth of every body
    circles = {}
    for b in bodies:
        o, R, s = frames[id(b)]
        c = o + R @ b.centre_offset() * s
        circles[id(b)] = (c, b.extent * s * max(1.0, b.k))
    kf = {id(b): b.k for b in bodies}  # lobe size factors of this instant
    for k, ms in stars.items():
        b = ms.body
        o, R, s = frames[id(b)]
        pts = ms.surface.points
        kl = kf[id(b)]
        if kl != 1.0:  # a detached lobe at its own size, about its centre
            cx = np.array([float(b.stars.index(k)), 0.0, 0.0])
            pts = cx + (pts - cx) * kl
        sa2 = (s * kl) ** 2  # body-unit areas -> AU^2
        p = o + (pts @ R.T) * s  # sky positions
        nrm = ms.surface.normals @ R.T
        mu = -nrm[:, 2]  # +z points away from the observer
        vis = mu > 0.0
        cover = np.zeros(len(mu))  # covered fraction of each element (sphere occulters)
        re_max = math.sqrt(float(ms.surface.areas.max()) * sa2 / math.pi)
        for ob in bodies:
            c, rad = circles[id(ob)]
            pad = re_max if ob.kind == "sphere" else 3.0 * re_max  # partly covered edge elements
            cand = vis & ((p[:, 0] - c[0]) ** 2 + (p[:, 1] - c[1]) ** 2 < (rad + pad) ** 2)
            cand &= p[:, 2] > c[2] - rad  # the body lies (partly) in front of the point
            if ob is b and b.kind != "roche":
                continue  # a convex single star cannot hide its own visible elements
            if ob is b:
                # a lobe hides its own elements only through the other lobe (or, in a contact
                # envelope, near the neck): elements behind the other lobe's sky circle, or
                # within 0.15 separations of the inner Lagrangian point
                si = ob.stars.index(k)
                other = ob.surfaces[1 - si]
                oc = np.array([1.0 - si, 0.0, 0.0])  # body frame centre of the other star
                ext_o = float(np.max(np.linalg.norm(other.points - oc, axis=1))) * kl
                cs = o + R @ oc * s
                near_o = (p[:, 0] - cs[0]) ** 2 + (p[:, 1] - cs[1]) ** 2 < (ext_o * s + pad) ** 2
                xl1 = roche.lagrange_points(ob.q)[0]
                neck = np.abs(ms.surface.points[:, 0] - xl1) < 0.15
                cand &= near_o | neck
            if not np.any(cand):
                continue
            oo, RR, ss = frames[id(ob)]
            idx = np.flatnonzero(cand)
            start = p[idx]
            # the ray toward the observer crosses the body's enclosing sphere between t0, t1
            # rays of the body's own elements start a few element sizes away from the surface
            # (the margin is ~0 on the surface itself, which must not count as an occultation)
            re_el = np.sqrt(ms.surface.areas[idx] * sa2 / math.pi)
            eps = 3.0 * re_el if ob is b else 0.0
            t1 = start[:, 2] - (c[2] - rad)
            t0 = np.maximum(start[:, 2] - (c[2] + rad), eps)
            t1 = np.maximum(t1, t0)
            frac = np.linspace(0.0, 1.0, N_RAY)[None, :]
            t = t0[:, None] + frac * (t1 - t0)[:, None]  # (m, N_RAY)
            if ob.kind == "sphere":  # the sky circle of a sphere in front, anti-aliased:
                # each element counts as a small disk of its projected area
                rr = ob.scale
                dd = np.hypot(p[idx, 0] - oo[0], p[idx, 1] - oo[1])
                re = np.sqrt(ms.surface.areas[idx] * sa2 * mu[idx] / math.pi)
                frac = np.clip((rr - dd) / (2.0 * re) + 0.5, 0.0, 1.0)
                cover[idx] = np.maximum(cover[idx], np.where(oo[2] < p[idx, 2], frac, 0.0))
                continue
            ray = start[:, None, :].repeat(N_RAY, axis=1)
            ray[..., 2] -= t
            mg = ob.margin(((ray - oo) @ RR) / ss)  # sky -> body frame
            if ob is b:
                mg[:, 0] = -np.inf  # the first sample is the element itself
            jbest = np.argmax(mg, axis=1)
            rows_all = np.arange(len(jbest))
            best = mg[rows_all, jbest]
            tb = t[rows_all, jbest]
            # Smooth occultation: an element is a small disk of its projected radius re; the
            # covered fraction follows from the signed distance of its line of sight from the
            # occulter's limb, d = margin / |grad margin| at the deepest point of the ray
            # (anti-aliasing as for spheres). A hard inside/outside test would make the flux
            # a step function of the parameters and break the derivatives of the fit.
            re = np.sqrt(ms.surface.areas[idx] * sa2 * np.clip(mu[idx], 0.0, 1.0) / math.pi)

            def signed_distance(tq, rows, start=start, oo=oo, RR=RR, ss=ss, ob=ob):
                q0 = start[rows].copy()
                q0[:, 2] -= tq
                qb = ((q0 - oo) @ RR) / ss
                m0 = ob.margin(qb)
                h = 1e-5
                g = np.empty_like(qb)
                for ax in range(3):
                    dq = np.zeros(3)
                    dq[ax] = h
                    g[:, ax] = (ob.margin(qb + dq) - ob.margin(qb - dq)) / (2.0 * h)
                gn = np.maximum(np.linalg.norm(g, axis=1), 1e-12)
                return m0, m0 / gn * ss

            m_b, dist = signed_distance(tb, rows_all)
            # refine the deepest point where the element is near the limb (or the sampling
            # may have missed a grazing chord)
            near = (dist > -6.0 * re) & (dist < 6.0 * re) | ((best <= 0.0) & (best > -0.5))
            if np.any(near):
                rows = np.flatnonzero(near)
                jj = jbest[rows]
                step = (t1 - t0)[rows] / (N_RAY - 1)
                lo = np.maximum(t[rows, jj] - step, t0[rows])
                hi = np.minimum(t[rows, jj] + step, t1[rows])
                gr = (math.sqrt(5.0) - 1.0) / 2.0
                for _ in range(24):
                    c1 = hi - gr * (hi - lo)
                    c2 = lo + gr * (hi - lo)
                    p1 = start[rows].copy()
                    p2 = start[rows].copy()
                    p1[:, 2] -= c1
                    p2[:, 2] -= c2
                    f1 = ob.margin(((p1 - oo) @ RR) / ss)
                    f2 = ob.margin(((p2 - oo) @ RR) / ss)
                    left = f1 > f2
                    hi = np.where(left, c2, hi)
                    lo = np.where(left, lo, c1)
                tq = 0.5 * (lo + hi)
                m_r, d_r = signed_distance(tq, rows)
                better = m_r > m_b[rows]
                dist[rows] = np.where(better, d_r, dist[rows])
            frac = np.clip(dist / (2.0 * np.maximum(re, 1e-12)) + 0.5, 0.0, 1.0)
            cover[idx] = np.maximum(cover[idx], frac)
        L = lights[k]
        i1 = L.weight / ld_mean(L.ld)  # intensity at disk centre
        lmu = ld_profile(L.ld, np.clip(mu, 0.0, 1.0))
        area_au = ms.surface.areas * sa2
        w_free = (i1 * ms.rel_intensity * lmu * mu * area_au)[vis]
        w = w_free * (1.0 - cover[vis])
        out[k] = float(np.sum(w))
        if elements is not None:
            elements[k] = (p[vis], w, w_free)
    return out


def _free_fluxes(b: Body, stars: dict, lights, n: int, l_body: np.ndarray) -> np.ndarray:
    """Fluxes of the stars of a Roche body seen alone from the direction l_body (body frame,
    unit vector to the observer), at separation 1 (multiply by the separation squared)."""
    zs = -np.asarray(l_body, float)
    xs = np.cross([0.0, 0.0, 1.0], zs)
    if np.linalg.norm(xs) < 1e-9:
        xs = np.array([1.0, 0.0, 0.0])
    xs /= np.linalg.norm(xs)
    ys = np.cross(zs, xs)
    i, j = b.stars
    pos = np.zeros((n, 3))
    vel = np.zeros((n, 3))
    pos[j] = [xs[0], ys[0], zs[0]]  # body x axis on the sky
    vel[j] = [xs[1], ys[1], zs[1]]  # body y axis: the orbital velocity direction
    sub = {k: stars[k] for k in (i, j)}
    a_ref, k_now = b.a_ref, b.k  # the lobes at separation 1, as built; the lobe size factor
    b.a_ref = 0.0  # of the current instant is restored afterwards (frames are reused)
    try:
        f = mesh_fluxes([b], sub, pos, vel, np.ones(n), lights)
    finally:
        b.a_ref, b.k = a_ref, k_now
    return f[[i, j]]


def _phase_table(b: Body, stars: dict, lights, n: int, inc: float):
    """Periodic spline over the orbital phase of the fluxes of the two stars of an unocculted
    Roche body (separation 1) seen at the inclination inc [deg]."""
    from scipy.interpolate import CubicSpline

    # the Roche geometry is symmetric about the orbital x-z plane, f(phi) = f(-phi): half a
    # revolution is computed and mirrored
    half = TABLE_PHASES // 2
    ph = np.linspace(0.0, 2.0 * math.pi, 2 * half + 1)
    si, ci = math.sin(math.radians(inc)), math.cos(math.radians(inc))
    first = np.array(
        [
            _free_fluxes(b, stars, lights, n, np.array([si * math.cos(a), si * math.sin(a), ci]))
            for a in ph[: half + 1]
        ]
    )
    vals = np.vstack([first, first[-2::-1]])
    return CubicSpline(ph, vals, axis=0, bc_type="periodic")


def _table_fluxes(b: Body, stars: dict, lights, n: int, R: np.ndarray) -> np.ndarray:
    """Fluxes of the two stars of an unocculted Roche body (separation 1) from the phase tables:
    tables on a fixed grid of inclinations (TABLE_DI_DEG), interpolated linearly between the two
    that bracket the current inclination (a continuous function of the inclination, as the
    derivatives of a fit need)."""
    l_body = R.T @ np.array([0.0, 0.0, -1.0])
    inc = math.degrees(math.acos(max(-1.0, min(1.0, float(l_body[2])))))
    phi = math.atan2(float(l_body[1]), float(l_body[0])) % (2.0 * math.pi)
    if b.tables is None:
        b.tables = {}
    x = inc / TABLE_DI_DEG
    i0 = math.floor(x)
    w = x - i0
    out = np.zeros(2)
    for idx, wt in ((i0, 1.0 - w), (i0 + 1, w)):
        if wt == 0.0:
            continue
        if idx not in b.tables:
            numeric = [k for k in b.tables if not isinstance(k, str)]
            if len(numeric) >= TABLE_MAX_N:
                b.tables.pop(numeric[0])
            b.tables[idx] = _phase_table(b, stars, lights, n, min(idx * TABLE_DI_DEG, 180.0))
        out += wt * np.asarray(b.tables[idx](phi), float)
    return out


def _overlapping(bodies, frames) -> dict:
    """Bodies whose sky circles overlap another body's: id(body) -> bool."""
    circ = []
    for b in bodies:
        o, R, s = frames[id(b)]
        c = o + R @ b.centre_offset() * s
        circ.append((c, b.extent * s * max(1.0, b.k)))
    out = {}
    for a, b in enumerate(bodies):
        ov = False
        for c2, (c, r) in enumerate(circ):
            if c2 == a:
                continue
            ca, ra = circ[a]
            if (ca[0] - c[0]) ** 2 + (ca[1] - c[1]) ** 2 < (ra + r) ** 2 * 1.02:
                ov = True
                break
        out[id(b)] = ov
    return out


def scene_flux(bodies, stars, pos, vel, masses, lights, pairs, proximity=True):
    """Total and per-star fluxes normalised like photometry.flux_at."""
    base = np.array([math.pi * s.radius_au**2 * s.weight for s in lights])
    if FAST:
        f = np.zeros(len(lights))
        frames = _frames(bodies, pos, vel, masses)
        ov = _overlapping(bodies, frames)
        rest = {}
        for b in bodies:
            if ov[id(b)] or b.kind == "rotating":
                rest.update({k: stars[k] for k in b.stars})
            elif b.kind == "sphere":
                # the unocculted flux of its own mesh (identical to the full computation, so
                # that switching between the two modes does not make a step)
                k = b.stars[0]
                if b.tables is None or "sphere" not in b.tables:
                    b.tables = b.tables or {}
                    b.tables["sphere"] = float(
                        mesh_fluxes([b], {k: stars[k]}, pos, vel, masses, lights, frames)[k]
                    )
                f[k] = b.tables["sphere"]
            else:  # an unocculted Roche body: the phase table, scaled by the separation^2
                _, R, sc = frames[id(b)]
                kl = b.k  # lobe size factor (frame), before the table may be rebuilt
                i, j = b.stars
                d_sky = math.hypot(pos[j][0] - pos[i][0], pos[j][1] - pos[i][1])
                touch = SELF_ECLIPSE_MARGIN * (lights[i].radius_au + lights[j].radius_au)
                if abs(kl - 1.0) > TABLE_MAX_DK and d_sky < touch:
                    # the stars eclipse each other at a separation other than the one of the
                    # table (eccentric or perturbed orbit): the table holds the eclipse
                    # geometry of separation a_ref only, so this phase is computed in full
                    rest.update({k: stars[k] for k in b.stars})
                    continue
                f[b.stars] = _table_fluxes(b, stars, lights, len(lights), R) * (sc * kl) ** 2
        if rest:
            fr = mesh_fluxes(bodies, rest, pos, vel, masses, lights, frames)
            for k in rest:
                f[k] = fr[k]
    else:
        f = mesh_fluxes(bodies, stars, pos, vel, masses, lights)
    if proximity:
        from .photometry import proximity_terms

        # ellipsoidal variation: from the geometry for Roche bodies, the analytic term (as in
        # photometry.flux_at, applied before the eclipses) for all other pairs; reflection
        # (Lambert spheres of the volume radii) also between detached Roche lobes, left out
        # only inside the common envelope of a contact binary
        roche_pairs = {tuple(b.stars) for b in bodies if b.kind == "roche"}
        contact_pairs = {tuple(b.stars) for b in bodies if b.kind == "roche" and b.contact}
        ell, refl, beam = proximity_terms(
            pos,
            vel,
            masses,
            lights,
            pairs,
            base,
            ellipsoidal={p for p in pairs if tuple(p) not in roche_pairs},
            reflection={p for p in pairs if tuple(p) not in contact_pairs},
        )
        # the reflected light of a star is eclipsed with the star (its visible fraction of
        # the disk of its volume radius, as in photometry.flux_at)
        if np.any(refl):
            from .photometry import visible_fractions

            refl = refl * visible_fractions(np.asarray(pos, float), lights)
        f = (f * (1.0 + ell) + refl) * beam
    return float(f.sum() / base.sum()), f / base.sum()


__all__ = ["Body", "MeshStar", "build_scene", "mesh_fluxes", "needed", "scene_flux"]
