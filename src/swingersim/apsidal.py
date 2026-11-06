"""Apsidal motion of the eclipsing pairs, split into its contributions (swingersim apsidal).

For every binary of two stars: the classical advance of Sterne (1939) from the tidal and the
rotational distortion of both stars, separately, with the apsidal constants k2, radii and spins
that SWINGERsim uses (stellar model or per-star values), for synchronous, pseudo-synchronous (Hut
1981) and the configured rotation; the relativistic advance; optionally the total advance
measured in an N-body run (which adds the perturbations by the other stars). With a target
rate (e.g. a published apsidal period) it gives the factor by which k2 of both stars would have
to change to reach it.

Published photodynamical solutions often fit the apsidal advance rate of a pair as a free
parameter (Lightcurvefactory: Zasche et al. 2023 for V994 Her) instead of computing it from
stellar structure; a difference against swingersim then measures the tension between the measured
rate and stellar models, not a difference between two codes.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from . import units
from .config import Config
from .tides import apsidal_rate_gr, apsidal_rate_sterne, pseudo_sync_factor


def _deg_yr(rate_rad_yr: float) -> float:
    return math.degrees(rate_rad_yr)


def _period(rate_deg_yr: float) -> float:
    return 360.0 / rate_deg_yr if rate_deg_yr > 0 else math.inf


def apsidal_table(cfg: Config, targets: dict | None = None, nbody_years: float | None = None):
    """Rows (dicts) of the apsidal contributions of every binary of two stars, with the
    stellar properties, radii and apsidal constants (k2 corrected for rotation with
    physics.k2_rotation) of the N-body run: the stellar models of the track labels (initial
    masses), the dynamics with the current masses (backends.nbody.model_props)."""
    from . import massloss
    from .backends.nbody import build_simulation, model_props, parent_of, spin_vector
    from .hierarchy import Node, build_tree, initialize, node_elements
    from .tides import k2_rotating

    targets = targets or {}
    cfg = massloss.resolve_masses(cfg)
    _, _, _, root, masses = build_simulation(cfg, gr=False, structure=False)
    pos, vel = initialize(root, masses)
    m_lab = np.array(build_tree(cfg)[1], float)  # track labels, as in build_simulation
    props_all = model_props(cfg, root, m_lab, pos, vel, masses)
    rows = []
    measured = {}
    if nbody_years:
        from .backends.nbody import run_nbody
        from .diagnostics.binary import measured_apsidal_period

        run = cfg.model_copy(
            update={
                "run": cfg.run.model_copy(update={"t_max_yr": float(nbody_years), "n_output": 2001})
            }
        )
        res = run_nbody(run)
        for name, o in res.orbits.items():
            measured[name] = measured_apsidal_period(res.times, o["Omega"], o["omega"])
    for node in root.nodes():
        if isinstance(node.left, Node) or isinstance(node.right, Node):
            continue
        i, j = node.left, node.right
        el = node_elements(node, pos, vel, masses)
        m = (float(masses[i]), float(masses[j]))
        props = [props_all[k] for k in (i, j)]
        r = tuple(p.radius_rsun * units.RSUN_AU for p in props)
        n = math.sqrt(units.G * (m[0] + m[1]) / el.a**3)
        f_ps = pseudo_sync_factor(el.e)
        # spin_vector is zero without physics.rotation
        spin_cfg = tuple(
            float(
                np.linalg.norm(spin_vector(cfg, f"m{k + 1}", parent_of(root, k), pos, vel, masses))
            )
            for k in (i, j)
        )
        k2 = tuple(
            k2_rotating(cfg, k, p.apsidal_k2, w, float(masses[k]), p.radius_rsun)
            for k, p, w in zip((i, j), props, spin_cfg, strict=True)
        )
        tidal = _deg_yr(apsidal_rate_sterne(m, r, k2, el.a, el.e, (0.0, 0.0)))
        rot = {}
        for lab, w in (
            ("sync", (n, n)),
            ("pseudosync", (f_ps * n, f_ps * n)),
            ("configured", spin_cfg),
        ):
            rot[lab] = _deg_yr(apsidal_rate_sterne(m, r, k2, el.a, el.e, w)) - tidal
        gr = _deg_yr(apsidal_rate_gr(m[0] + m[1], el.a, el.e))
        total = tidal + rot["configured"] + gr
        row = {
            "pair": node.name,
            "stars": f"m{i + 1}+m{j + 1}",
            "P_d": el.a**1.5
            / math.sqrt(units.G * (m[0] + m[1]))
            * 2
            * math.pi
            * units.DAYS_PER_YEAR,
            "e": el.e,
            "k2_1": k2[0],
            "k2_2": k2[1],
            "R1_rsun": props[0].radius_rsun,
            "R2_rsun": props[1].radius_rsun,
            "P_rot_1_d": (2 * math.pi / spin_cfg[0] * units.DAYS_PER_YEAR)
            if spin_cfg[0]
            else math.inf,
            "P_rot_2_d": (2 * math.pi / spin_cfg[1] * units.DAYS_PER_YEAR)
            if spin_cfg[1]
            else math.inf,
            "pseudosync_factor": f_ps,
            "tidal_deg_yr": tidal,
            "rot_sync_deg_yr": rot["sync"],
            "rot_pseudosync_deg_yr": rot["pseudosync"],
            "rot_configured_deg_yr": rot["configured"],
            "gr_deg_yr": gr,
            "classical_gr_deg_yr": total,
            "U_classical_gr_yr": _period(total),
        }
        if node.name in measured and el.e < 0.02:
            row["U_nbody_yr"] = math.nan  # omega is not defined well for a nearly circular orbit
            row["dynamical_deg_yr"] = math.nan
        elif node.name in measured:
            u = measured[node.name]
            row["U_nbody_yr"] = u
            row["dynamical_deg_yr"] = (360.0 / u - total) if np.isfinite(u) and u != 0 else math.nan
        t = targets.get(node.name)
        if t is not None:
            # target as an apsidal period [yr]; the k2 factor scales tidal and rotational terms
            target = 360.0 / t
            other = (
                gr + row.get("dynamical_deg_yr", 0.0)
                if np.isfinite(row.get("dynamical_deg_yr", 0.0))
                else gr
            )
            row["U_target_yr"] = t
            row["k2_factor_for_target"] = (target - other) / (tidal + rot["configured"])
        rows.append(row)
    return rows


def report(rows) -> str:
    lines = [
        "swingersim apsidal: apsidal advance of the eclipsing pairs [deg/yr]",
        (
            "  pair  stars    P [d]     e      k2 (1, 2)          tidal  rot sync  pseudo  config"
            "   GR     total   U [yr]"
        ),
    ]
    for r in rows:
        lines.append(
            f"  {r['pair']:<5} {r['stars']:<8} {r['P_d']:8.4f} {r['e']:6.4f}  "
            f"{r['k2_1']:.5f} {r['k2_2']:.5f}  {r['tidal_deg_yr']:7.3f} {r['rot_sync_deg_yr']:8.3f}"
            f" {r['rot_pseudosync_deg_yr']:7.3f} {r['rot_configured_deg_yr']:7.3f} "
            f"{r['gr_deg_yr']:6.3f} {r['classical_gr_deg_yr']:8.3f} {r['U_classical_gr_yr']:8.1f}"
        )
        if "U_nbody_yr" in r and not np.isfinite(r["U_nbody_yr"]):
            lines.append("        N-body: no apsidal period (e < 0.02, omega not well defined)")
        elif "U_nbody_yr" in r:
            lines.append(
                f"        N-body: U = {r['U_nbody_yr']:.1f} yr, dynamical (other stars) "
                f"{r['dynamical_deg_yr']:+.3f} deg/yr"
            )
        if "U_target_yr" in r:
            lines.append(
                f"        target U = {r['U_target_yr']:.1f} yr needs k2 x "
                f"{r['k2_factor_for_target']:.2f} (both stars)"
            )
    lines.append(
        "  tidal + rotational after Sterne (1939); rotation 'config' = P_rot of the "
        "configuration (sync, pseudosync after Hut 1981, or a period); GR after Einstein"
    )
    return "\n".join(lines)


def write_csv(rows, path: Path) -> Path:
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return Path(path)
