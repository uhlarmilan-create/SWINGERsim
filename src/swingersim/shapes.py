"""Choice of the star shapes for light curves: spheres where they are accurate enough, Roche
surfaces (swingersim.mesh) where the stars fill their Roche lobes.

The error of spheres against Roche surfaces grows steeply with the filling of the lobe
(volume radius / Eggleton radius; measured with swingersim.mesh for a 1.5-d pair): about 0.5 ppt
rms at 0.3, 2 ppt at 0.5, 5 ppt at 0.7, 9 ppt at 0.85 and 17 ppt at 0.95 (maxima 2-3 times
larger), and 15-45 ppt for contact binaries. The Roche model of swingersim assumes a circular,
synchronised orbit; eccentric pairs keep spheres (with a note).
"""

from __future__ import annotations

import copy

from .config import Config

FILL_MAX_DETACHED = 0.99  # Eggleton fill above which a detached Roche model is not chosen


def fill_factors(cfg: Config) -> dict:
    """{pair: (fill of the first star, fill of the second, e)} for every binary of two stars:
    volume radius over the Roche-lobe radius of Eggleton (1983) at the semi-major axis. The
    radii are those of the light curves (photometry.light_curve: stellar models of the track
    labels, backends.nbody.model_props)."""
    import numpy as np

    from . import massloss, units
    from .backends.nbody import build_simulation, model_props
    from .hierarchy import build_tree, initialize
    from .stellar import roche_lobe_eggleton

    cfg = massloss.resolve_masses(cfg)
    _, _, _, root, masses = build_simulation(cfg, gr=False, structure=False)
    m_lab = np.array(build_tree(cfg)[1], float)
    props = model_props(cfg, root, m_lab, *initialize(root, masses), masses)
    out = {}
    for node in root.nodes():
        if not (isinstance(node.left, int) and isinstance(node.right, int)):
            continue
        i, j = node.left, node.right
        a_rsun = float(node.elements.a) / units.RSUN_AU
        ri, rj = props[i].radius_rsun, props[j].radius_rsun
        qi, qj = masses[i] / masses[j], masses[j] / masses[i]
        out[node.name] = (
            ri / (roche_lobe_eggleton(float(qi)) * a_rsun),
            rj / (roche_lobe_eggleton(float(qj)) * a_rsun),
            float(node.elements.e),
            f"m{i + 1}",
            f"m{j + 1}",
        )
    return out


def apply_policy(d: dict, policy: str, threshold: float = 0.5) -> tuple[dict, list[str]]:
    """Copy of the configuration dict d with the star shapes of the light-curve model, and
    notes on what was chosen."""
    notes: list[str] = []
    if policy == "as_configured":
        return d, notes
    out = copy.deepcopy(d)
    if policy == "spheres":
        for st in (out.get("stars") or {}).values():
            if isinstance(st, dict):
                st.pop("shape", None)
        for o in (out.get("orbits") or {}).values():
            if isinstance(o, dict):
                o.pop("contact_fillout", None)
        return out, notes
    cfg = Config.model_validate(d)
    ff = fill_factors(cfg)
    stars = out.setdefault("stars", {}) or {}
    out["stars"] = stars
    for pair, (f1, f2, e, s1, s2) in ff.items():
        orbit = (out.get("orbits") or {}).get(pair) or {}
        if orbit.get("contact_fillout") is not None:
            notes.append(f"{pair}: contact binary (contact_fillout), Roche envelope")
            continue
        shapes = [(stars.get(s) or {}).get("shape", "sphere") for s in (s1, s2)]
        if "roche" in shapes:
            notes.append(f"{pair}: Roche surfaces as configured (fill {f1:.2f}, {f2:.2f})")
            continue
        if max(f1, f2) < threshold:
            continue
        if max(f1, f2) >= FILL_MAX_DETACHED:
            notes.append(
                f"{pair}: fill {f1:.2f}, {f2:.2f} - a star fills its Roche lobe; a detached "
                "Roche pair cannot model it (contact binaries: contact_fillout); spheres kept"
            )
            continue
        if e > 0.02:
            notes.append(
                f"{pair}: fill {f1:.2f}, {f2:.2f} but e = {e:.3f} - the Roche model needs a "
                "circular orbit; spheres kept"
            )
            continue
        for s in (s1, s2):
            st = stars.get(s) or {}
            stars[s] = st
            st["shape"] = "roche"
        notes.append(f"{pair}: fill {f1:.2f}, {f2:.2f} >= {threshold} - Roche surfaces")
    return out, notes
