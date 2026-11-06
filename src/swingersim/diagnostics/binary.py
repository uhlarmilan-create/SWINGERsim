"""Diagnostics of the binaries made of two stars: apsidal motion and tidal time scales.

For eclipsing binaries the apsidal motion is directly observable (O-C of primary and secondary
minima). swingersim reports the classical expectation of Sterne (1939) from the tidal and
rotational distortion of the stars, the relativistic contribution, and the rate measured in the
simulation, which also contains the perturbation by the other stars of the quadruple.
"""

from __future__ import annotations

import numpy as np

from .. import units
from ..hierarchy import Node, node_elements
from ..stellar import StarProperties
from ..tides import apsidal_rate_gr, apsidal_rate_sterne, leconte_rates


def _period(rate: float) -> float:
    return 2.0 * np.pi / rate if rate > 0 else np.inf


# below this eccentricity the line of apsides is not defined well enough to fit its advance
# (osculating omega of a nearly circular orbit is dominated by short-period terms); the same
# guard as swingersim.apsidal
APSIDAL_E_MIN = 0.02


def measured_apsidal_period(
    times: np.ndarray, Omega: np.ndarray, omega: np.ndarray, e: np.ndarray | None = None
) -> float:
    """Period [yr] of the longitude of periastron varpi = Omega + omega from a linear fit,
    negative when the line of apsides regresses.

    Pass ``Omega`` = 0 to fit the argument of periastron omega alone. The outputs where an
    angle is undefined (NaN) are dropped before unwrapping; with the eccentricities ``e`` the
    fit returns NaN when the mean eccentricity is below APSIDAL_E_MIN."""
    t = np.asarray(times, float)
    varpi = np.asarray(Omega, float) + np.asarray(omega, float)
    ok = np.isfinite(varpi) & np.isfinite(t)
    if e is not None:
        e = np.asarray(e, float)
        ok &= np.isfinite(e)
        if ok.sum() and float(np.mean(e[ok])) < APSIDAL_E_MIN:
            return float("nan")  # line of apsides not defined
    if ok.sum() < 3:
        return float("nan")
    t, varpi = t[ok], np.unwrap(varpi[ok])
    if abs(varpi[-1] - varpi[0]) < 0.2:
        return float("nan")  # too little advance for a meaningful fit
    rate = np.polyfit(t, varpi, 1)[0]
    return float(2.0 * np.pi / rate)  # negative: the line of apsides regresses


def eclipsing(orbit: dict, r_sum_au: float) -> np.ndarray:
    """At each output: does the binary eclipse for an observer along the line of sight (z)?

    Spherical stars on a circular orbit of the current semi-major axis eclipse when
    |cos i| < (R1 + R2) / a (the criterion of the viewer and of Vokrouhlicky 2016, Sect. 3.3).
    """
    with np.errstate(invalid="ignore"):
        return np.abs(np.cos(orbit["inc"])) < r_sum_au / orbit["a"]


def radius_sum_au(res, props: list[StarProperties], i: int, j: int):
    """R_i + R_j [AU]: constant, or at every output for evolving stars (res.stellar)."""
    st = getattr(res, "stellar", None) or {}
    ri = st.get(f"m{i + 1}", {}).get("radius_rsun", props[i].radius_rsun)
    rj = st.get(f"m{j + 1}", {}).get("radius_rsun", props[j].radius_rsun)
    return (np.asarray(ri) + np.asarray(rj)) * units.RSUN_AU


def both_eclipsing_fraction(res, root: Node, props: list[StarProperties]) -> float:
    """2+2: fraction of outputs at which both inner binaries eclipse."""
    flags = []
    for node in (root.left, root.right):
        r = radius_sum_au(res, props, node.left, node.right)
        flags.append(eclipsing(res.orbits[node.name], r))
    return float(np.mean(flags[0] & flags[1]))


def eclipse_summary(res, root: Node, props: list[StarProperties]) -> dict:
    """Simultaneous eclipses: for every 2+2 quadruple of the hierarchy the fraction of outputs
    at which both of its binaries eclipse (``eclipse_fraction_both`` for the first one, also
    stored in res.system["quadruples"]), and with three or more two-star binaries the fraction
    at which all of them eclipse (``eclipse_fraction_all``)."""
    from .zlk import quadruple_cores

    out = {}
    cores = [
        c for c in quadruple_cores(root) if isinstance(c.left, Node) and isinstance(c.right, Node)
    ]
    quads = {q["core"]: q for q in res.system.get("quadruples", [])}
    for k, core in enumerate(cores):
        f = both_eclipsing_fraction(res, core, props)
        if k == 0:
            out["eclipse_fraction_both"] = f
        if core.name in quads:
            quads[core.name]["eclipse_fraction_both"] = f
    pairs = [
        n for n in root.nodes() if not isinstance(n.left, Node) and not isinstance(n.right, Node)
    ]
    if len(pairs) >= 3:
        flags = [
            eclipsing(res.orbits[n.name], radius_sum_au(res, props, n.left, n.right)) for n in pairs
        ]
        out["eclipse_fraction_all"] = float(np.mean(np.logical_and.reduce(flags)))
    return out


def binary_diagnostics(
    root: Node,
    masses: np.ndarray,
    props: list[StarProperties],
    spins: list[np.ndarray],
    pos0: np.ndarray,
    vel0: np.ndarray,
    res=None,
    gr: bool = True,
    structure: bool = False,
    cfg=None,
) -> dict[str, dict[str, float]]:
    """cfg: with physics.k2_rotation the apsidal constants are corrected for the rotation as in
    the simulation (tides.k2_rotating)."""
    out = {}
    for node in root.nodes():
        if isinstance(node.left, Node) or isinstance(node.right, Node):
            continue
        i, j = node.left, node.right
        el = node_elements(node, pos0, vel0, masses)
        m = (masses[i], masses[j])
        r = (props[i].radius_rsun * units.RSUN_AU, props[j].radius_rsun * units.RSUN_AU)
        w = (float(np.linalg.norm(spins[i])), float(np.linalg.norm(spins[j])))
        k2 = (props[i].apsidal_k2, props[j].apsidal_k2)
        if cfg is not None:
            from ..tides import k2_rotating

            k2 = tuple(
                k2_rotating(cfg, k, k2[n], w[n], float(masses[k]), props[k].radius_rsun)
                for n, k in enumerate((i, j))
            )
        d = {
            "U_gr_yr": _period(apsidal_rate_gr(m[0] + m[1], el.a, el.e)),
            "U_tidal_yr": _period(apsidal_rate_sterne(m, r, k2, el.a, el.e, w))
            if structure
            else np.inf,
        }
        total = (2 * np.pi / d["U_gr_yr"] if gr else 0.0) + 2 * np.pi / d["U_tidal_yr"]
        d["U_expected_yr"] = _period(total)
        d["gr_included"] = float(gr)
        lag = (props[i].lag_s / units.SECONDS_PER_YEAR, props[j].lag_s / units.SECONDS_PER_YEAR)
        if structure and max(lag) > 0:
            rates = leconte_rates(m, r, k2, (props[i].rg2, props[j].rg2), lag, el.a, el.e, w)
            d.update(rates.timescales(el.a, el.e))
        if res is not None:
            o = res.orbits[node.name]
            d["U_measured_yr"] = measured_apsidal_period(
                res.times, o["Omega"], o["omega"], o.get("e")
            )
            # omega alone (from the node on the sky): differs from varpi when the plane of the
            # binary precesses
            d["U_measured_omega_yr"] = measured_apsidal_period(
                res.times, 0.0 * np.asarray(o["omega"]), o["omega"], o.get("e")
            )
            d["eclipse_fraction"] = float(np.mean(eclipsing(o, radius_sum_au(res, props, i, j))))
        out[node.name] = {k: float(v) for k, v in d.items()}
        out[node.name]["stars"] = f"m{i + 1}+m{j + 1}"
    return out
