"""Radial velocities of the centre of light (fit.rv.centre: light).

A spectrograph measures the Doppler shift of the light of a star, i.e. the flux-weighted mean
line-of-sight velocity of its visible surface, not the velocity of its centre of mass. For
spherical, slowly rotating stars the two agree; in close binaries they differ (Wilson & Sofia
1976; Wilson 1979; Kopal 1959):

* a tidally distorted (Roche) star is brighter and larger on some sides, and its synchronous
  rotation adds velocity to its surface, so the centre of light moves on a different orbit
  than the centre of mass (the "proximity effects" in the radial velocities; corrections of
  several km/s in contact binaries, e.g. Rucinski 1992);
* during an eclipse the occulted part of a rotating star is missing from the mean, which
  shifts the velocity (the Rossiter-McLaughlin effect; Rossiter 1924, McLaughlin 1924).

swingersim computes this on the surface meshes of swingersim.mesh, the same surfaces, gravity and
limb darkening (passband `band`) and occultations as its light curves: every visible element
contributes its flux and the line-of-sight velocity of its rigid motion - stars of a Roche or
contact binary corotate with their orbit (v = v_cm + w x (r - r_cm), w the angular velocity of
the line of centres), other stars rotate with stars.<name>.P_rot (default: synchronous with
their orbit) about their spin axis (obliquity, spin_azimuth). Irradiation (reflection) and
Doppler beaming are not included in the weights. The offset light - centre of mass is computed
on the mesh and added to the exact centre-of-mass velocity (with the light time) of
swingersim.radvel; for spheres and rotating stars the offset of the unocculted star (zero by
symmetry, a few 0.01 km/s of mesh discretisation) is subtracted. Validated against the flux-weighted radial
velocities of PHOEBE 2 (benchmark rv-light-vs-phoebe).
"""

from __future__ import annotations

import numpy as np

from . import units

KMS_PER_AUYR = units.AU_KM / units.SECONDS_PER_YEAR


def _spin(cfg, root, k: int, pos, vel, masses) -> np.ndarray:
    """Angular velocity [rad/yr] of a non-Roche star: stars.<name>.P_rot, or synchronous with
    the orbit it belongs to (sync / pseudosync / not given), about its spin axis."""
    from .backends.nbody import parent_of
    from .hierarchy import node_elements, node_h_unit

    spec = cfg.stars.get(f"m{k + 1}")
    p_rot = spec.P_rot if spec is not None else "sync"
    node = parent_of(root, k)
    if p_rot in ("sync", "pseudosync"):
        el = node_elements(node, pos, vel, masses)
        m_tot = float(np.sum(masses[node.stars()]))
        omega = 2.0 * np.pi / units.sma_to_period(el.a, m_tot)
        if p_rot == "pseudosync":
            from .tides import pseudo_sync_factor

            omega *= pseudo_sync_factor(float(el.e))
    else:
        omega = 2.0 * np.pi / units.days_to_years(float(p_rot))
    h = node_h_unit(node, pos, vel, masses)
    nvec = np.cross([0.0, 0.0, 1.0], h)
    nn = np.linalg.norm(nvec)
    nvec = nvec / nn if nn > 1e-12 else np.array([1.0, 0.0, 0.0])
    ob = np.radians(spec.obliquity if spec is not None else 0.0)
    az = np.radians(spec.spin_azimuth if spec is not None else 0.0)
    axis = np.cos(ob) * h + np.sin(ob) * (np.cos(az) * nvec + np.sin(az) * np.cross(h, nvec))
    return omega * axis


def light_offset(cfg, star, t_day, band: str = "V", light_time: bool = True) -> np.ndarray:
    """Centre of light minus centre of mass: radial velocity [km/s] of the light of the stars
    `star` (names m1, m2, ... or indices, one per time) minus the velocity of their centres of
    mass, at the observing times t_day [d after t = 0]. Both from the same state, so that the
    offset is free of the light-time approximations of the light-curve states; add it to the
    exact centre-of-mass velocities (swingersim.radvel / obsmodel)."""
    from . import massloss, mesh
    from .backends.nbody import model_props
    from .hierarchy import build_tree, initialize
    from .obsmodel import star_indices
    from .photometry import SeenStates, star_lights

    cfg = massloss.resolve_masses(cfg)
    seen = SeenStates(cfg, light_time)
    m_lab = np.array(build_tree(cfg)[1], float)
    root, masses = seen.root, massloss.initial_masses(cfg, m_lab)
    p0, v0 = initialize(root, masses)
    props = model_props(cfg, root, m_lab, p0, v0, masses)
    lights = star_lights(cfg, masses, props, band, None)
    bodies, stars = mesh.build_scene(cfg, root, masses, props, lights, p0, v0, band)
    kk = star_indices(np.asarray(star).astype(str), len(masses))
    t_day = np.asarray(t_day, float)
    out = np.full(len(t_day), np.nan)
    order = np.argsort(t_day, kind="stable")
    for j in order:
        pos, vel = seen.at(t_day[j] / units.DAYS_PER_YEAR)
        k = int(kk[j])
        m_now = masses
        frames = mesh._frames(bodies, pos, vel, m_now)
        el: dict = {}
        mesh.mesh_fluxes(bodies, {k: stars[k]}, pos, vel, m_now, lights, frames, elements=el)
        p, w, w_free = el.get(k, (np.zeros((0, 3)), np.zeros(0), np.zeros(0)))
        b = stars[k].body
        if b.kind == "roche":
            i, jj = b.stars
            r = pos[jj] - pos[i]
            omega = np.cross(r, vel[jj] - vel[i]) / float(r @ r)
            mi, mj = m_now[i], m_now[jj]
            rc = (mi * pos[i] + mj * pos[jj]) / (mi + mj)
            vc = (mi * vel[i] + mj * vel[jj]) / (mi + mj)
            vz = vc[2] + np.cross(omega, p - rc)[:, 2]
        else:
            om = _spin(cfg, root, k, pos, vel, m_now)
            vz = vel[k][2] + np.cross(om, p - pos[k])[:, 2]
        sw = float(np.sum(w))
        light = float(np.sum(w * vz)) / sw if sw > 0.0 else float(vel[k][2])
        if b.kind != "roche" and np.sum(w_free) > 0.0:
            # an unocculted sphere or rotating star is symmetric about its projected spin axis:
            # its offset is zero; subtracting the free mean removes the mesh's discretisation
            light -= float(np.sum(w_free * vz)) / float(np.sum(w_free)) - float(vel[k][2])
        out[j] = (light - float(vel[k][2])) * KMS_PER_AUYR
    return out


def light_rv(cfg, star, t_day, bjd_t0: float = 0.0, band: str = "V") -> np.ndarray:
    """Centre-of-light radial velocity [km/s, + receding, without the systemic velocity]: the
    exact centre-of-mass velocity (radvel.model_rv, with the light time) plus light_offset."""
    from .radvel import model_rv

    t_day = np.asarray(t_day, float)
    rv = model_rv(cfg, star, bjd_t0 + t_day, bjd_t0)
    return rv + light_offset(cfg, star, t_day, band)


__all__ = ["light_offset", "light_rv"]
