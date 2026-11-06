"""Dynamical stability of hierarchical multiple systems (quadruples and larger hierarchies).

Three tools, from cheap to definitive:

1. Algebraic criteria applied to the two "nested" triples of the quadruple (2+2: A with B as a
   point mass, B with A as a point mass; 3+1: inner with star 3, mid with star 4):
   * Mardling & Aarseth (2001), MA01;
   * the updated criterion of Vynatheya et al. (2022) MNRAS 516, 4146, eqs. (3)-(4), which adds
     a dependence on the inner eccentricity (including its ZLK growth) and a better
     inclination dependence.
2. The machine-learning classifiers (MLP) of Vynatheya, Mardling & Hamers (2023) MNRAS 525,
   2388, trained directly on 2+2 and 3+1 quadruples. They are distributed by the authors
   (github.com/pavanvyn/quadruple-stability) and used only if that repository is present.
   Larger hierarchies are classified sub-quadruple by sub-quadruple (deeper sub-systems as
   point masses), an approximation outside the training set.
3. A direct N-body test with the definition of stability of Vynatheya et al. (2023), Sect. 2:
   the system and a "ghost" copy with the inner semi-major axis larger by 1e-6 AU are
   integrated for 100 outer orbits (Newtonian point masses); the system is unstable if it
   becomes unbound or if the relative difference of the inner semi-major axes exceeds 1e-2.
"""

from __future__ import annotations

import importlib.util
import os
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import units
from .config import Config
from .hierarchy import Node, build_tree, initialize, mass_of, node_elements, node_h_unit

# --------------------------------------------------------------------------- algebraic criteria


def ma01_ratio(a_in, a_out, e_out, q_out, i_mut_rad) -> float:
    """R_p / R_p,crit of Mardling & Aarseth (2001), as written in Vynatheya et al. (2022) eq. 2.

    q_out = m_out / m_in; > 1 means stable.
    """
    crit = (
        2.8
        * ((1.0 + q_out) * (1.0 + e_out) / np.sqrt(1.0 - e_out)) ** 0.4
        * (1.0 - 0.3 * i_mut_rad / np.pi)
    )
    return float(a_out * (1.0 - e_out) / a_in / crit)


def v22_ratio(a_in, a_out, e_in, e_out, q_out, i_mut_rad) -> float:
    """Y~ / Y~crit of Vynatheya et al. (2022), eqs. (3)-(4); > 1 means stable."""
    c = np.cos(i_mut_rad)
    e_max = np.sqrt(max(0.0, 1.0 - 5.0 / 3.0 * c * c))
    e_t = max(e_in, 0.5 * e_max**2)
    y = a_out * (1.0 - e_out) / (a_in * (1.0 + e_t))
    y_crit = (
        2.4
        * ((1.0 + q_out) / ((1.0 + e_t) * np.sqrt(1.0 - e_out))) ** 0.4
        * ((1.0 - 0.2 * e_t + e_out) / 8.0 * (c - 1.0) + 1.0)
    )
    return float(y / y_crit)


@dataclass
class NestedTriple:
    label: str  # e.g. "A in out (+B)"
    inner: str
    outer: str
    a_in: float
    a_out: float
    e_in: float
    e_out: float
    q_out: float  # mass of the third body / mass of the inner binary
    i_mut_deg: float
    ma01: float
    v22: float


def nested_triples(root: Node, masses: np.ndarray, pos, vel) -> list[NestedTriple]:
    out = []
    for parent in root.nodes():
        for child, sibling in ((parent.left, parent.right), (parent.right, parent.left)):
            if not isinstance(child, Node):
                continue
            ei = node_elements(child, pos, vel, masses)
            eo = node_elements(parent, pos, vel, masses)
            q = mass_of(sibling, masses) / mass_of(child, masses)
            hi, ho = node_h_unit(child, pos, vel, masses), node_h_unit(parent, pos, vel, masses)
            i = float(np.arccos(np.clip(hi @ ho, -1.0, 1.0)))
            sib = sibling.name if isinstance(sibling, Node) else f"m{sibling + 1}"
            out.append(
                NestedTriple(
                    label=f"{child.name} in {parent.name} (+{sib})",
                    inner=child.name,
                    outer=parent.name,
                    a_in=ei.a,
                    a_out=eo.a,
                    e_in=ei.e,
                    e_out=eo.e,
                    q_out=q,
                    i_mut_deg=float(np.degrees(i)),
                    ma01=ma01_ratio(ei.a, eo.a, eo.e, q, i),
                    v22=v22_ratio(ei.a, eo.a, ei.e, eo.e, q, i),
                )
            )
    return out


# --------------------------------------------------------------------------- MLP classifiers


MODEL_VERSIONS = ("_v1.2.2", "_v1.0.2", "")  # file suffixes in the authors' repository


def mlp_directory() -> Path | None:
    """Directory of the Vynatheya et al. (2023) repository, if available.

    Looked up in $SWINGERSIM_VYNATHEYA_DIR, then in data/quadruple-stability next to the package
    source tree and in the current directory.
    """
    cands = []
    if os.environ.get("SWINGERSIM_VYNATHEYA_DIR"):
        cands.append(Path(os.environ["SWINGERSIM_VYNATHEYA_DIR"]))
    here = Path(__file__).resolve()
    cands += [
        here.parents[2] / "data" / "quadruple-stability",
        Path.cwd() / "data" / "quadruple-stability",
    ]
    for c in cands:
        if (c / "classify_quad_2p2.py").exists() and _model_file(c, "2p2") is not None:
            return c
    return None


def _model_file(d: Path, kind: str) -> Path | None:
    for v in MODEL_VERSIONS:
        f = d / f"mlp_model_{kind}_ghost{v}.pkl"
        if f.exists():
            return f
    return None


def _load(path: Path, module: str):
    spec = importlib.util.spec_from_file_location(module, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sub_quadruples(root: Node) -> list[tuple]:
    """The 2+2 and 3+1 quadruples contained in a hierarchy, with deeper sub-systems replaced
    by point masses at their centres of mass:

    * ("2+2", l, c1, c2): orbit l whose two components are both orbits (binaries c1, c2);
    * ("3+1", d, c, l): orbit l with a component c that is an orbit, which in turn has a
      component d that is an orbit (inner d, middle c, outer l).

    A 2+2 or 3+1 quadruple gives exactly one entry; a 2+2+2 sextuple gives 2+2 (A, B in AB)
    and 3+1 entries (A in AB in out, B in AB in out)."""
    out = []
    for l in root.nodes():
        kids = [c for c in (l.left, l.right) if isinstance(c, Node)]
        if len(kids) == 2:
            out.append(("2+2", l, kids[0], kids[1]))
        for c in kids:
            for d in (c.left, c.right):
                if isinstance(d, Node):
                    out.append(("3+1", d, c, l))
    return out


def _quad_label(q: tuple) -> str:
    return (
        f"2+2 ({q[2].name}, {q[3].name} in {q[1].name})"
        if q[0] == "2+2"
        else (f"3+1 ({q[1].name} in {q[2].name} in {q[3].name})")
    )


def _quad_parameters(q: tuple, pos, vel, masses) -> dict:
    """The 11 scale-free parameters of Vynatheya et al. (2023), Sect. 2, with their ordering
    conventions (q <= 1, heavier binary first for 2+2), for one sub-quadruple. Inclinations
    in radians."""

    def el(n):
        return node_elements(n, pos, vel, masses)

    def imut(a, b):
        return float(
            np.arccos(
                np.clip(node_h_unit(a, pos, vel, masses) @ node_h_unit(b, pos, vel, masses), -1, 1)
            )
        )

    def comps(n):
        return mass_of(n.left, masses), mass_of(n.right, masses)

    if q[0] == "2+2":
        _, l, c1, c2 = q
        b1, b2 = (c1, c2) if sum(comps(c1)) >= sum(comps(c2)) else (c2, c1)
        m1, m2 = comps(b1), comps(b2)
        return {
            "qi1": min(m1) / max(m1),
            "qi2": min(m2) / max(m2),
            "qo": sum(m2) / sum(m1),
            "ali1o": el(b1).a / el(l).a,
            "ali2o": el(b2).a / el(l).a,
            "ei1": el(b1).e,
            "ei2": el(b2).e,
            "eo": el(l).e,
            "ii1i2": imut(b1, b2),
            "ii1o": imut(b1, l),
            "ii2o": imut(b2, l),
        }
    _, d, c, l = q
    md = comps(d)
    m_in = sum(md)
    m3 = mass_of(c.right if c.left is d else c.left, masses)
    m4 = mass_of(l.right if l.left is c else l.left, masses)
    return {
        "qi": min(md) / max(md),
        "qm": m3 / m_in,
        "qo": m4 / (m_in + m3),
        "alim": el(d).a / el(c).a,
        "almo": el(c).a / el(l).a,
        "ei": el(d).e,
        "em": el(c).e,
        "eo": el(l).e,
        "iim": imut(d, c),
        "iio": imut(d, l),
        "imo": imut(c, l),
    }


def initial_state(cfg: Config) -> tuple[Node, np.ndarray, np.ndarray, np.ndarray]:
    """(tree, masses, positions, velocities) of the dynamics at t = 0, as the N-body backend
    sets them up: masses resolved (physics.mass_loss masses: current), the current masses of
    the tracks at t = 0 and the periods of the orbits with those masses."""
    from . import massloss

    cfg = massloss.resolve_masses(cfg)
    root, m_lab = build_tree(cfg)
    masses = massloss.initial_masses(cfg, m_lab)
    if not np.array_equal(masses, m_lab):
        root = build_tree(cfg, kepler_masses=masses)[0]
    pos, vel = initialize(root, masses)
    return root, masses, pos, vel


def mlp_parameters(cfg: Config) -> dict:
    """Parameters of the (first) sub-quadruple; for a 2+2 or 3+1 the quadruple itself."""
    root, masses, pos, vel = initial_state(cfg)
    quads = sub_quadruples(root)
    if not quads:
        raise ValueError("no 2+2 or 3+1 sub-quadruple in this hierarchy")
    return _quad_parameters(quads[0], pos, vel, masses)


def mlp_verdicts(cfg: Config) -> list[tuple[str, bool | None]] | None:
    """MLP verdict (True = stable) for every sub-quadruple of the hierarchy; None if the
    classifiers are not installed or the hierarchy has fewer than four stars. A verdict is
    None for a sub-quadruple whose kind (2+2 or 3+1) has no model file."""
    d = mlp_directory()
    if d is None:
        return None
    root, masses, pos, vel = initial_state(cfg)
    quads = sub_quadruples(root)
    if not quads:
        return None
    out = []
    with warnings.catch_warnings():
        # the networks were pickled with scikit-learn 1.2.2; according to the authors' README
        # newer versions also work but warn about a possible incompatibility
        warnings.simplefilter("ignore")
        for q in quads:
            p = _quad_parameters(q, pos, vel, masses)
            out.append((_quad_label(q), _mlp_predict(q[0], d, p)))
    return out


def mlp_stable(cfg: Config) -> bool | None:
    """Prediction of the Vynatheya et al. (2023) MLP (True = stable); None if unavailable.
    For hierarchies of more than four stars every 2+2 and 3+1 sub-quadruple (deeper
    sub-systems as point masses) must be classified stable; this extension beyond the training
    set of the classifiers is an approximation."""
    v = mlp_verdicts(cfg)
    if v is None or any(ok is None for _, ok in v):
        return None
    return all(ok for _, ok in v)


def _mlp_predict(kind: str, d: Path, p: dict) -> bool | None:
    """Verdict of the classifier of ``kind``; None if its script or model file is missing."""
    tag = "2p2" if kind == "2+2" else "3p1"
    model, script = _model_file(d, tag), d / f"classify_quad_{tag}.py"
    if model is None or not script.exists():
        return None
    mod = _load(script, f"classify_quad_{tag}")
    return bool(getattr(mod, f"mlp_classifier_{tag}")(str(model), *p.values()))


# --------------------------------------------------------------------------- N-body ghost test


@dataclass
class GhostResult:
    stable: bool
    reason: str  # "stable", "unbound", "diverged"
    t_end_yr: float
    outer_orbits: float
    delta_max: float
    detail: dict = field(default_factory=dict)


def _ghost_target(cfg: Config) -> str:
    """Inner orbit whose semi-major axis is perturbed: the lighter inner binary (2+2) or the
    innermost binary (3+1), as in Vynatheya et al. (2023)."""
    if cfg.topology == "3+1":
        return "inner"
    m = cfg.masses
    if cfg.topology == "2+2":
        return "A" if m["m1"] + m["m2"] < m["m3"] + m["m4"] else "B"
    # other hierarchies: the lightest binary made of two stars
    binaries = [n for n, (a, b) in cfg.tree.items() if a in m and b in m]
    return min(binaries, key=lambda n: sum(m[s] for s in cfg.tree[n]))


def ghost_test(
    cfg: Config,
    n_outer: float = 100.0,
    delta_crit: float = 1e-2,
    d_a_au: float = 1e-6,
    max_wall_s: float | None = None,
) -> GhostResult:
    """N-body stability test of Vynatheya et al. (2023), Newtonian point masses. With
    max_wall_s a test that takes longer is abandoned (reason 'timeout', stable False; such
    systems were left out of the training sets of Vynatheya et al., Sect. 3.1)."""
    import time

    wall0 = time.perf_counter()
    from . import massloss
    from .backends.nbody import build_simulation

    base = massloss.resolve_masses(cfg).model_copy(deep=True)
    base.physics.gr = False
    base.physics.tides = False
    base.physics.rotation = False
    base.run.collisions = False
    target = _ghost_target(base)
    # the tree and masses of the dynamics (current masses of the tracks, periods with them)
    _, _, _, root, masses = build_simulation(base, gr=False, structure=False)
    nodes = {n.name: n for n in root.nodes()}
    ghost = base.model_copy(deep=True)
    spec = ghost.orbits[target]
    a0 = nodes[target].elements.a
    ghost.orbits[target] = spec.model_copy(update={"a": a0 + d_a_au, "P": None})
    base.orbits[target] = base.orbits[target].model_copy(update={"a": a0, "P": None})

    sims = []
    for c in (base, ghost):
        sim, *_ = build_simulation(c, gr=False, structure=False)
        sims.append(sim)
    p_out = units.sma_to_period(root.elements.a, float(masses.sum()))
    t_end = n_outer * p_out
    checks = np.linspace(0.0, t_end, int(max(50, 10 * n_outer)) + 1)[1:]
    delta_max = 0.0
    tn = nodes[target]
    unbound_streak = 0
    for t in checks:
        if max_wall_s is not None and time.perf_counter() - wall0 > max_wall_s:
            return GhostResult(False, "timeout", t, t / p_out, delta_max)
        a = []
        unbound = None
        for sim in sims:
            sim.integrate(t, exact_finish_time=1)
            n = len(sim.particles)
            pos, vel = np.zeros((n, 3)), np.zeros((n, 3))
            sim.serialize_particle_data(xyz=pos, vxvyvz=vel)
            for node in root.nodes():
                try:
                    node_elements(node, pos, vel, masses)
                except ValueError:
                    unbound = node.name
                    break
            if unbound is None:
                a.append(node_elements(tn, pos, vel, masses).a)
        # an orbit counts as unbound only if it stays so (a passing star can make the
        # osculating two-body energy of a binary positive for a moment)
        unbound_streak = unbound_streak + 1 if unbound else 0
        if unbound_streak >= 3:
            return GhostResult(False, "unbound", t, t / p_out, delta_max, {"orbit": unbound})
        if len(a) == 2:
            delta = abs(a[0] - a[1]) / a[0]
            delta_max = max(delta_max, delta)
            if delta > delta_crit:
                return GhostResult(False, "diverged", t, t / p_out, delta_max)
    return GhostResult(True, "stable", t_end, n_outer, delta_max)
