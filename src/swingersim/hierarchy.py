"""Hierarchical (Jacobi-tree) initialisation: 2+2 and 3+1 quadruples, and any other
hierarchy of binaries (e.g. a quadruple with a fifth star, see config.PRESETS).

Each tree node is a binary made of two children (a star index or another node).
Its orbit describes the relative motion of child 2 around child 1, both taken
as their centres of mass, with mu = G * (M_child1 + M_child2):

    2+2:  out = [ A = (1, 2) ,  B = (3, 4) ]
    3+1:  out = [ mid = [ inner = (1, 2) , 3 ] , 4 ]
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import units
from .config import Config, OrbitSpec
from .kepler import Elements, angular_momentum_unit, elements_to_state, state_to_elements


@dataclass
class Node:
    name: str
    left: Node | int
    right: Node | int
    elements: Elements | None = None

    def stars(self) -> list[int]:
        out: list[int] = []
        for c in (self.left, self.right):
            out.extend(c.stars() if isinstance(c, Node) else [c])
        return out

    def nodes(self) -> list[Node]:
        """All binaries in the tree, outermost first."""
        out = [self]
        for c in (self.left, self.right):
            if isinstance(c, Node):
                out.extend(c.nodes())
        return out


def _stars(child: Node | int) -> list[int]:
    return child.stars() if isinstance(child, Node) else [child]


def mass_of(child: Node | int, masses: np.ndarray) -> float:
    return float(masses[_stars(child)].sum())


def skeleton(topology: str | dict) -> Node:
    """Tree of a preset topology name or of a hierarchy mapping (orbit -> two children)."""
    from .config import PRESETS, check_hierarchy

    spec = PRESETS[topology] if isinstance(topology, str) else topology
    root, _stars = check_hierarchy(spec)

    def make(name: str) -> Node | int:
        if name in spec:
            return Node(name, make(spec[name][0]), make(spec[name][1]))
        return int(name[1:]) - 1

    return make(root)


def _to_elements(spec: OrbitSpec, m_total: float) -> Elements:
    a = spec.a if spec.a is not None else units.period_to_sma(units.days_to_years(spec.P), m_total)
    return Elements(
        a=a,
        e=spec.e,
        inc=spec.i * units.DEG,
        Omega=spec.Omega * units.DEG,
        omega=spec.omega * units.DEG,
        M=spec.M0 * units.DEG,
    )


def build_tree(cfg: Config, kepler_masses=None) -> tuple[Node, np.ndarray]:
    """Tree with elements filled in, and the mass array (index k = star m<k+1>).

    Orbits given by their period P get a from Kepler's third law with ``kepler_masses`` (default:
    the configured masses). With mass loss along the tracks the configured masses label the
    tracks (ZAMS) while the dynamics starts with the current masses: pass those, so that the
    period at t = 0 is the configured one."""
    masses = np.array([cfg.masses[k] for k in cfg.star_names], float)
    mk = masses if kepler_masses is None else np.asarray(kepler_masses, float)
    root = skeleton(cfg.tree)
    for node in root.nodes():
        m_tot = mass_of(node.left, mk) + mass_of(node.right, mk)
        node.elements = _to_elements(cfg.orbits[node.name], m_tot)
    return root, masses


def initialize(root: Node, masses: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Barycentric positions [AU] and velocities [AU/yr] of the stars."""
    pos = np.zeros((len(masses), 3))
    vel = np.zeros((len(masses), 3))

    def place(child: Node | int, R: np.ndarray, V: np.ndarray) -> None:
        if not isinstance(child, Node):
            pos[child], vel[child] = R, V
            return
        m1, m2 = mass_of(child.left, masses), mass_of(child.right, masses)
        M = m1 + m2
        r, v = elements_to_state(units.G * M, child.elements)
        place(child.left, R - (m2 / M) * r, V - (m2 / M) * v)
        place(child.right, R + (m1 / M) * r, V + (m1 / M) * v)

    place(root, np.zeros(3), np.zeros(3))
    return pos, vel


def _com(child: Node | int, pos, vel, masses) -> tuple[np.ndarray, np.ndarray]:
    idx = _stars(child)
    m = masses[idx]
    return (m @ pos[idx]) / m.sum(), (m @ vel[idx]) / m.sum()


def relative_state(node: Node, pos, vel, masses) -> tuple[np.ndarray, np.ndarray]:
    r1, v1 = _com(node.left, pos, vel, masses)
    r2, v2 = _com(node.right, pos, vel, masses)
    return r2 - r1, v2 - v1


def node_elements(node: Node, pos, vel, masses) -> Elements:
    r, v = relative_state(node, pos, vel, masses)
    M = mass_of(node.left, masses) + mass_of(node.right, masses)
    return state_to_elements(units.G * M, r, v)


def node_h_unit(node: Node, pos, vel, masses) -> np.ndarray:
    r, v = relative_state(node, pos, vel, masses)
    return angular_momentum_unit(r, v)
