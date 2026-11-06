"""ZLK, precession and stability diagnostics for each inner orbit and its parent orbit.

For every orbit that has a parent in the hierarchy (A and B in 2+2; inner and mid in 3+1)
swingersim reports analytic expectations computed from the initial elements and the values
measured in the simulation. Formulas: see docs/en/concepts.md.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction

import numpy as np

from .. import units
from ..hierarchy import Node, mass_of, node_h_unit, relative_state
from ..secular.quadrupole import predict
from ..stability import v22_ratio
from .orbits import gr_periastron_advance_rate

ZLK_CRITICAL_DEG = float(np.degrees(np.arccos(np.sqrt(3.0 / 5.0))))  # 39.23 deg


@dataclass
class PairDiagnostics:
    inner: str
    outer: str
    # analytic, from the initial state
    i_mut0_deg: float
    in_zlk_window: bool
    e_max_quadrupole: float  # test-particle quadrupole estimate (initially circular)
    t_zlk_yr: float  # order-of-magnitude ZLK time scale
    t_prec_yr: float  # nodal precession period (quadrupole, circular inner orbit)
    t_gr_yr: float  # 1PN apsidal precession period of the inner orbit
    gr_suppresses_zlk: bool  # t_gr < t_zlk
    epsilon_oct: float  # octupole strength (Naoz 2016)
    epsilon_sa: float  # single-averaging parameter (Luo, Katz & Dong 2016); large = mild hierarchy
    stability_ratio: float  # (a_out/a_in) / Mardling-Aarseth critical value; > 1 = stable
    # measured in the simulation
    e_min: float
    e_max: float
    i_mut_min_deg: float
    i_mut_max_deg: float
    n_zlk_cycles: int  # completed eccentricity cycles with amplitude > 0.1
    zlk_period_yr: float  # mean period of those cycles (NaN if fewer than 2)
    # quadrupole secular prediction (Naoz et al. 2013), perturber as a point mass
    quad_e_min: float = float("nan")
    quad_e_max: float = float("nan")
    quad_i_min_deg: float = float("nan")
    quad_i_max_deg: float = float("nan")
    quad_period_yr: float = float("nan")
    g1_0_deg: float = float("nan")  # initial argument of periastron from the line of nodes
    # extremes between outputs and precession measured in the simulation
    e_max_true: float = float("nan")
    t_e_max_yr: float = float("nan")
    q_min_au: float = float("nan")  # smallest periastron distance a(1-e)
    q_min_over_radii: float = float("nan")  # q_min / (R1 + R2) for binaries of two stars
    prec_period_measured_yr: float = float("nan")  # nodal precession about the invariable plane
    stability_ratio_v22: float = float("nan")  # Vynatheya et al. (2022) criterion; > 1 = stable

    def as_dict(self) -> dict:
        return asdict(self)


def _parent_pairs(root: Node) -> list[tuple[Node, Node, Node | int]]:
    """(inner node, parent node, sibling of the inner node) for every nested orbit."""
    out = []
    for parent in root.nodes():
        for child, sibling in ((parent.left, parent.right), (parent.right, parent.left)):
            if isinstance(child, Node):
                out.append((child, parent, sibling))
    return out


def mardling_aarseth_ratio(a_in, a_out, e_out, m_in, m_out, i_mut_rad) -> float:
    """R_p / R_p,crit of Mardling & Aarseth (2001): stability.ma01_ratio with q_out = m_out / m_in."""
    from ..stability import ma01_ratio

    return ma01_ratio(a_in, a_out, e_out, m_out / m_in, i_mut_rad)


def _cycles(t: np.ndarray, e: np.ndarray) -> tuple[int, float]:
    good = np.isfinite(e)
    t, e = t[good], e[good]
    if len(e) < 3 or e.max() - e.min() < 0.1:
        return 0, float("nan")
    mid = 0.5 * (e.max() + e.min())
    up = np.flatnonzero((e[:-1] < mid) & (e[1:] >= mid))
    if len(up) < 2:
        return 0, float("nan")
    return len(up) - 1, float(abs(t[up[-1]] - t[up[0]]) / (len(up) - 1))


def _initial_g1(inner, parent, pos, vel, masses) -> float:
    """Argument of periastron of the inner orbit from its node on the local invariable plane."""
    r_in, v_in = relative_state(inner, pos, vel, masses)
    r_out, v_out = relative_state(parent, pos, vel, masses)
    m_a, m_b = mass_of(inner.left, masses), mass_of(inner.right, masses)
    m_in = m_a + m_b
    m_tot = mass_of(parent.left, masses) + mass_of(parent.right, masses)
    mu_in = units.G * m_in
    h_in = np.cross(r_in, v_in) * (m_a * m_b / m_in)
    h_out = np.cross(r_out, v_out) * (m_in * (m_tot - m_in) / m_tot)
    z = h_in + h_out
    z /= np.linalg.norm(z)
    hh = h_in / np.linalg.norm(h_in)
    node = np.cross(z, hh)
    if np.linalg.norm(node) < 1e-9:
        return 0.0
    node /= np.linalg.norm(node)
    vrel = v_in
    h = np.cross(r_in, vrel)
    evec = np.cross(vrel, h) / mu_in - r_in / np.linalg.norm(r_in)
    if np.linalg.norm(evec) < 1e-9:
        return 0.0
    return float(np.degrees(np.arctan2(np.cross(node, evec) @ hh, node @ evec)) % 360.0)


def _measured_precession(t: np.ndarray, omega_inv: np.ndarray) -> float:
    good = np.isfinite(omega_inv)
    if good.sum() < 3:
        return float("nan")
    w = np.unwrap(omega_inv[good])
    if abs(w[-1] - w[0]) < 2.0 * np.pi:  # not even one full circulation
        return float("nan")
    rate = np.polyfit(t[good], w, 1)[0]
    return float(2.0 * np.pi / abs(rate))


def nearest_ratio(x: float, max_den: int = 5) -> tuple[str, float]:
    """Closest small-integer ratio p:q to x and the relative deviation from it."""
    f = Fraction(x).limit_denominator(max_den)
    if f.numerator == 0 or f.numerator > 5 * max_den:
        return "-", float("nan")
    return f"{f.numerator}:{f.denominator}", float(x / float(f) - 1.0)


def jab_max_vokrouhlicky(root: Node, masses: np.ndarray, pos: np.ndarray, vel: np.ndarray) -> float:
    """Upper limit of the mutual inclination J_AB of binaries A and B [deg].

    Vokrouhlicky (2016), Sect. 2.3, eqs. (26)-(36): for circular A and B and quadrupole
    order, the normals k_A, k_B, k precess about each other with the first integrals
    C1 = s' x^2 + s'' y^2 and C2 = (s2/s'') x + (s1/s') y + z (x = k_A.k, y = k_B.k,
    z = k_A.k_B). The smallest z on their intersection, and the bound from |x|, |y| <= 1,
    give cos J_AB^max. NaN if the bound is not informative (z_min < -1).
    """
    a_node, b_node = root.left, root.right
    MA, MB = mass_of(a_node, masses), mass_of(b_node, masses)
    MAB = MA + MB
    aA, aB, a = a_node.elements.a, b_node.elements.a, root.elements.a
    eta = np.sqrt(1.0 - root.elements.e**2)
    G = units.G
    n, nA, nB = (np.sqrt(G * M / x**3) for M, x in ((MAB, a), (MA, aA), (MB, aB)))
    xAa, xAb = masses[a_node.left] / MA, masses[a_node.right] / MA
    xBa, xBb = masses[b_node.left] / MB, masses[b_node.right] / MB
    s1 = -0.75 / eta**3 * n * (n / nA) * MB / MAB
    s2 = -0.75 / eta**3 * n * (n / nB) * MA / MAB
    sp = -0.75 / eta**4 * n * xAa * xAb * (aA / a) ** 2
    spp = -0.75 / eta**4 * n * xBa * xBb * (aB / a) ** 2
    kA = node_h_unit(a_node, pos, vel, masses)
    kB = node_h_unit(b_node, pos, vel, masses)
    k = node_h_unit(root, pos, vel, masses)
    x0, y0, z0 = kA @ k, kB @ k, kA @ kB
    al, be = s2 / spp, s1 / sp
    c1 = abs(sp) * x0**2 + abs(spp) * y0**2
    zmin_ellipse = z0 + al * x0 + be * y0 - np.sqrt(c1 * (al**2 / abs(sp) + be**2 / abs(spp)))
    zmin_box = z0 + al * (x0 - 1.0) + be * (y0 - 1.0)
    zmin = max(zmin_ellipse, zmin_box)
    if zmin < -1.0:
        return float("nan")
    return float(np.degrees(np.arccos(min(zmin, 1.0))))


def _shape(node) -> str | None:
    """'2+2' or '3+1' if the node is such a quadruple, else None."""
    if not isinstance(node, Node) or len(node.stars()) != 4:
        return None
    left, right = node.left, node.right
    if isinstance(left, Node) and isinstance(right, Node):
        return "2+2" if len(left.stars()) == 2 and len(right.stars()) == 2 else None
    if isinstance(left, Node) and len(left.stars()) == 3 and isinstance(left.left, Node):
        return "3+1"
    return None


def quadruple_cores(root: Node) -> list[Node]:
    """All sub-trees of exactly four stars that form a 2+2 or 3+1 quadruple (the whole tree
    for a quadruple; e.g. AB in a 2+2+2 sextuple; AB and CD in ((AB)(CD))E), outermost first."""
    out = []
    for node in root.nodes():
        if _shape(node):
            out.append(node)
    return out


def system_diagnostics(
    root: Node,
    masses: np.ndarray,
    diags: list,
    init_state: tuple[np.ndarray, np.ndarray] | None = None,
) -> dict:
    """Relations between the orbits of every 2+2 or 3+1 quadruple in the hierarchy.

    Both topologies: the parameter beta of Hamers & Lai (2017), eqs. (32)-(33), the ratio of the
    precession rate of the perturber's orbit to the ZLK rate of an inner orbit; for beta ~ 0.3-3
    the inner eccentricity can be excited strongly even at low inclination (secular chaos).
    2+2 also: period and precession commensurabilities, analytic J_AB limit.

    The keys of the first (outermost) quadruple are given flat (``period_ratio_A_B`` ...);
    ``quadruples`` lists the results of every quadruple, which differ from the first only in
    hierarchies with several of them (e.g. ((AB)(CD))E).
    """
    d = {x.inner: x for x in diags}
    cores = quadruple_cores(root)
    out: dict = {}
    for k, core in enumerate(cores):
        res = _core_diagnostics(core, masses, d, init_state)
        if k == 0:
            out.update(res)
        else:
            out.update({key: v for key, v in res.items() if key.startswith("beta_")})
        out.setdefault("quadruples", []).append({"core": core.name, "shape": _shape(core), **res})
    return out


def _core_diagnostics(core: Node, masses, d: dict, init_state) -> dict:
    if not (isinstance(core.left, Node) and isinstance(core.right, Node)):
        n_in, n_mid = core.left.left, core.left
        mid = d[n_mid.name]
        m01 = mass_of(n_in, masses)
        m2 = mass_of(n_mid.right, masses)
        m3 = mass_of(core.right, masses)
        a1, a2, a3 = n_in.elements.a, n_mid.elements.a, core.elements.a
        e2, e3 = n_mid.elements.e, core.elements.e
        beta = (
            0.75
            * np.cos(np.radians(mid.i_mut0_deg))
            * (a2**3 / (a1 * a3**2)) ** 1.5
            * np.sqrt(m01 / (m01 + m2))
            * (m3 / m2)
            * ((1 - e2**2) / (1 - e3**2)) ** 1.5
        )
        return {f"beta_hamers_lai_{n_in.name}": float(beta)}
    a_node, b_node = core.left, core.right
    pa = units.sma_to_period(a_node.elements.a, mass_of(a_node, masses))
    pb = units.sma_to_period(b_node.elements.a, mass_of(b_node, masses))
    ratio = pa / pb
    pq, dev = nearest_ratio(ratio)
    da, db = d[a_node.name], d[b_node.name]
    out = {
        "binaries": [a_node.name, b_node.name],
        "period_ratio_A_B": float(ratio),
        "nearest_period_commensurability": pq,
        "period_commensurability_deviation": dev,
        "precession_ratio_analytic": float(da.t_prec_yr / db.t_prec_yr),
        "precession_ratio_measured": float(da.prec_period_measured_yr / db.prec_period_measured_yr),
    }
    ma, mb = mass_of(a_node, masses), mass_of(b_node, masses)
    aa, ab = a_node.elements.a, b_node.elements.a
    out[f"beta_hamers_lai_{a_node.name}"] = float(
        0.75 * np.cos(np.radians(db.i_mut0_deg)) * (ab / aa) ** 1.5 * (ma / mb) ** 1.5
    )
    out[f"beta_hamers_lai_{b_node.name}"] = float(
        0.75 * np.cos(np.radians(da.i_mut0_deg)) * (aa / ab) ** 1.5 * (mb / ma) ** 1.5
    )
    if init_state is not None:
        out["J_AB_max_analytic_deg"] = jab_max_vokrouhlicky(core, masses, *init_state)
    return out


def analyse(
    root: Node,
    masses: np.ndarray,
    gr: bool,
    res,
    init_state: tuple[np.ndarray, np.ndarray] | None = None,
    radii_au: np.ndarray | None = None,
) -> list[PairDiagnostics]:
    """Diagnostics for every nested orbit; `res` is a RunResult."""
    out = []
    for inner, parent, sibling in _parent_pairs(root):
        ei, eo = inner.elements, parent.elements
        m_in = mass_of(inner, masses)
        m_pert = mass_of(sibling, masses)
        m_tot = m_in + m_pert
        p_in = units.sma_to_period(ei.a, m_in)
        p_out = units.sma_to_period(eo.a, m_tot)
        key = f"{parent.name}-{inner.name}"
        i_mut = res.mutual[key]
        i0 = float(i_mut[0])
        m_a, m_b = mass_of(inner.left, masses), mass_of(inner.right, masses)

        in_window = ZLK_CRITICAL_DEG < np.degrees(i0) < 180.0 - ZLK_CRITICAL_DEG
        e_max_q = float(np.sqrt(max(0.0, 1.0 - 5.0 / 3.0 * np.cos(i0) ** 2))) if in_window else 0.0
        t_zlk = p_out**2 / p_in * m_tot / m_pert * (1.0 - eo.e**2) ** 1.5
        cos_i = abs(np.cos(i0))
        t_prec = 4.0 / 3.0 * t_zlk / cos_i if cos_i > 1e-3 else float("inf")
        t_gr = 2.0 * np.pi / gr_periastron_advance_rate(m_in, ei.a, ei.e)
        eps = abs(m_a - m_b) / (m_a + m_b) * ei.a / eo.a * eo.e / (1.0 - eo.e**2)
        eps_sa = (ei.a / (eo.a * (1.0 - eo.e**2))) ** 1.5 * m_pert / np.sqrt(m_in * m_tot)

        e_series = res.orbits[inner.name]["e"]
        n_cyc, p_cyc = _cycles(res.times, e_series)

        extra = {}
        if init_state is not None:
            g1 = _initial_g1(inner, parent, *init_state, masses)
            q = predict(m_a, m_b, m_pert, ei.a, eo.a, ei.e, eo.e, float(np.degrees(i0)), g1)
            extra.update(
                quad_e_min=q.e_min,
                quad_e_max=q.e_max,
                quad_i_min_deg=q.i_min_deg,
                quad_i_max_deg=q.i_max_deg,
                quad_period_yr=q.period_yr,
                g1_0_deg=g1,
            )
        ext = getattr(res, "extremes", {}).get(inner.name)
        if ext:
            extra.update(
                e_max_true=ext["e_max"], t_e_max_yr=ext["t_e_max_yr"], q_min_au=ext["q_min_au"]
            )
            if (
                radii_au is not None
                and not isinstance(inner.left, Node)
                and not isinstance(inner.right, Node)
            ):
                extra["q_min_over_radii"] = ext["q_min_au"] / (
                    radii_au[inner.left] + radii_au[inner.right]
                )
        if "Omega_inv" in res.orbits[inner.name]:
            extra["prec_period_measured_yr"] = _measured_precession(
                res.times, res.orbits[inner.name]["Omega_inv"]
            )
        out.append(
            PairDiagnostics(
                inner=inner.name,
                outer=parent.name,
                i_mut0_deg=float(np.degrees(i0)),
                in_zlk_window=bool(in_window),
                e_max_quadrupole=e_max_q,
                t_zlk_yr=float(t_zlk),
                t_prec_yr=float(t_prec),
                t_gr_yr=float(t_gr),
                gr_suppresses_zlk=bool(gr and t_gr < t_zlk),
                epsilon_oct=float(eps),
                epsilon_sa=float(eps_sa),
                stability_ratio=float(mardling_aarseth_ratio(ei.a, eo.a, eo.e, m_in, m_pert, i0)),
                stability_ratio_v22=v22_ratio(ei.a, eo.a, ei.e, eo.e, m_pert / m_in, i0),
                e_min=float(np.nanmin(e_series)),
                e_max=float(np.nanmax(e_series)),
                i_mut_min_deg=float(np.degrees(np.nanmin(i_mut))),
                i_mut_max_deg=float(np.degrees(np.nanmax(i_mut))),
                n_zlk_cycles=n_cyc,
                zlk_period_yr=p_cyc,
                **extra,
            )
        )
    return out
