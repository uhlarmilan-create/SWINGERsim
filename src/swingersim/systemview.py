"""What a configuration looks like: the model system as data for pictures (swingersim gui, reports).

``system_model(cfg)`` returns, for t = 0 of the configuration, the stars (mass, radius, effective
temperature, position) and every orbit of the hierarchy (elements, period, size, the paths of its two
components around their common centre of mass), the eclipse geometry of every pair of two stars
(inclination against the eclipse limit) and the stability ratios of every nested triple. Coordinates:
x, y on the sky and z along the line of sight, away from the observer (right-handed; the inclination
is the angle between the orbital angular momentum and z, the node Omega is measured from x towards y).
Which directions on the celestial sphere x and y are (north, east) is not fixed by the program:
dynamics, eclipses and light time depend only on the orientation of the orbits relative to each
other and to z.

Everything is plain JSON (lists and floats), so that the GUI and the report scripts draw the same
picture: a sky view (x, y), a top view (x, -z: the observer at the top), a zoom on every pair of two stars
with the stars to scale, and the hierarchy as a tree.
"""

from __future__ import annotations

import math

import numpy as np

from . import units
from .config import Config
from .hierarchy import Node, build_tree, initialize, mass_of, node_elements, node_h_unit
from .kepler import Elements, rotation_matrix

RSUN_AU = units.RSUN_AU
DPY = units.DAYS_PER_YEAR


def _f(x) -> float | None:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _ellipse(el: Elements, n: int = 240) -> np.ndarray:
    """Relative orbit (second component around the first) as n x 3 points [AU]."""
    E = np.linspace(0.0, 2.0 * np.pi, n)
    xp = el.a * (np.cos(E) - el.e)
    yp = el.a * np.sqrt(max(1.0 - el.e**2, 0.0)) * np.sin(E)
    R = rotation_matrix(el.inc, el.Omega, el.omega)
    return (R @ np.vstack([xp, yp, np.zeros_like(xp)])).T


def _name(child, star_names: list[str]) -> str:
    return child.name if isinstance(child, Node) else star_names[child]


def _label(child, star_names) -> str:
    if isinstance(child, Node):
        return "(" + " + ".join(_label(c, star_names) for c in (child.left, child.right)) + ")"
    return star_names[child]


def system_model(cfg: Config, n_points: int = 240) -> dict:
    from .backends.nbody import model_props

    root, m_lab = build_tree(cfg)
    masses = np.asarray(m_lab, float)
    pos, vel = initialize(root, masses)
    star_names = [f"m{k + 1}" for k in range(len(masses))]
    try:
        props = model_props(cfg, root, masses, pos, vel)
    except Exception:  # noqa: BLE001 - a picture without radii is better than none
        props = [None] * len(masses)
    stars = []
    for k, m in enumerate(masses):
        spec = cfg.stars.get(star_names[k]) if cfg.stars else None
        p = props[k]
        R = _f(getattr(p, "radius_rsun", None)) or _f(getattr(spec, "R", None))
        T = _f(getattr(spec, "Teff", None)) or _f(getattr(p, "teff_k", None))
        stars.append(
            {
                "name": star_names[k],
                "mass": _f(m),
                "R_rsun": R,
                "R_au": R * RSUN_AU if R else None,
                "Teff": T,
                "pos": [float(v) for v in pos[k]],
            }
        )
    orbits = []
    pairs = []
    for node in root.nodes():
        el = node_elements(node, pos, vel, masses)
        m1, m2 = mass_of(node.left, masses), mass_of(node.right, masses)
        M = m1 + m2
        P_d = 2.0 * np.pi * math.sqrt(el.a**3 / (units.G * M)) * DPY if el.a > 0 else None
        sl, sr = node.left, node.right
        idx_l = sl.stars() if isinstance(sl, Node) else [sl]
        idx_r = sr.stars() if isinstance(sr, Node) else [sr]
        com_l = (masses[idx_l, None] * pos[idx_l]).sum(0) / masses[idx_l].sum()
        com_r = (masses[idx_r, None] * pos[idx_r]).sum(0) / masses[idx_r].sum()
        com = (m1 * com_l + m2 * com_r) / M
        rel = _ellipse(el, n_points)
        path_l = com - (m2 / M) * rel
        path_r = com + (m1 / M) * rel
        h = node_h_unit(node, pos, vel, masses)
        inc = float(np.degrees(el.inc))
        o = {
            "name": node.name,
            "children": [_name(sl, star_names), _name(sr, star_names)],
            "label": _label(node, star_names),
            "masses": [_f(m1), _f(m2)],
            "P_d": _f(P_d),
            "a_au": _f(el.a),
            "e": _f(el.e),
            "i_deg": _f(inc),
            "Omega_deg": _f(np.degrees(el.Omega) % 360.0),
            "omega_deg": _f(np.degrees(el.omega) % 360.0),
            "com": [float(v) for v in com],
            "h": [float(v) for v in h],
            "path_left": np.round(path_l, 8).tolist(),
            "path_right": np.round(path_r, 8).tolist(),
            "two_stars": not isinstance(sl, Node) and not isinstance(sr, Node),
        }
        orbits.append(o)
        if o["two_stars"]:
            R1 = stars[sl]["R_au"] or 0.0
            R2 = stars[sr]["R_au"] or 0.0
            # separation at the nearer conjunction (z is the line of sight)
            sin_w = math.sin(el.omega)
            r_conj = el.a * (1.0 - el.e**2) / (1.0 + el.e * abs(sin_w))
            i_fold = inc if inc <= 90.0 else 180.0 - inc
            lim = math.degrees(math.acos(min(1.0, (R1 + R2) / r_conj))) if r_conj > 0 else None
            impact = r_conj * math.cos(math.radians(i_fold)) / (R1 + R2) if R1 + R2 > 0 else None
            pairs.append(
                {
                    "orbit": node.name,
                    "stars": [star_names[sl], star_names[sr]],
                    "i_deg": _f(i_fold),
                    "i_limit_deg": _f(lim),
                    "impact": _f(impact),
                    "eclipsing": bool(lim is not None and i_fold > lim),
                    "fill": [
                        _f(stars[sl]["R_au"] / _roche(el.a, m1, m2)) if stars[sl]["R_au"] else None,
                        _f(stars[sr]["R_au"] / _roche(el.a, m2, m1)) if stars[sr]["R_au"] else None,
                    ],
                }
            )
    triples = []
    try:
        from .stability import nested_triples

        for t in nested_triples(root, masses, pos, vel):
            triples.append(
                {
                    "label": t.label,
                    "inner": t.inner,
                    "outer": t.outer,
                    "a_ratio": _f(t.a_in / t.a_out),
                    "i_mut_deg": _f(t.i_mut_deg),
                    "ma01": _f(t.ma01),
                    "v22": _f(t.v22),
                    "stable": bool(t.ma01 > 1.0 and t.v22 > 1.0),
                }
            )
    except Exception:  # noqa: BLE001, S110 - stability ratios are optional in the picture
        pass
    # mutual inclinations of all pairs of orbits
    mutual = []
    for a_ in range(len(orbits)):
        for b_ in range(a_ + 1, len(orbits)):
            ha, hb = np.array(orbits[a_]["h"]), np.array(orbits[b_]["h"])
            mutual.append(
                {
                    "orbits": [orbits[a_]["name"], orbits[b_]["name"]],
                    "deg": _f(np.degrees(np.arccos(np.clip(ha @ hb, -1, 1)))),
                }
            )
    return {
        "hierarchy": _label(root, star_names),
        "stars": stars,
        "orbits": orbits,
        "pairs": pairs,
        "triples": triples,
        "mutual": mutual,
        "stable": all(t["stable"] for t in triples) if triples else None,
    }


def _roche(a: float, m: float, m_other: float) -> float:
    """Eggleton (1983) Roche-lobe radius of the star of mass m [AU]."""
    from .stellar import roche_lobe_eggleton

    return a * float(roche_lobe_eggleton(m / m_other))


def system_model_from_file(path) -> dict:
    from .config import load_config

    return system_model(load_config(path))
