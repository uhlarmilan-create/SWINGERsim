"""Hierarchy detection: which stars are bound to which, at a given moment.

Algorithm (in the spirit of the binary-tree classification used in few-body scattering
codes): start with every star as a separate object. Repeatedly take the pair of objects whose
two-body relative orbit is bound and has the smallest semi-major axis, and replace the pair by
their centre of mass. Stop when no bound pair remains. The result is a list of top-level
trees; a stable quadruple gives exactly one tree containing all four stars.

Trees are nested tuples of star indices (0-based), e.g. ((0, 1), (2, 3)) for a 2+2 system.
"""

from __future__ import annotations

import math
from itertools import combinations

import numpy as np


def canonical(tree):
    """Sort children by their lowest star index so equal hierarchies compare equal."""
    if isinstance(tree, (int, np.integer)):
        return int(tree)
    a, b = (canonical(t) for t in tree)
    return (a, b) if _lowest(a) < _lowest(b) else (b, a)


def _lowest(tree) -> int:
    return tree if isinstance(tree, int) else min(_lowest(t) for t in tree)


def _size(tree) -> int:
    return 1 if isinstance(tree, int) else sum(_size(t) for t in tree)


def detect(pos: np.ndarray, vel: np.ndarray, masses: np.ndarray, G: float) -> tuple:
    """Top-level trees, sorted, as a tuple (one element for a bound quadruple)."""
    # plain Python floats: for a handful of stars much faster than numpy per pair
    objs = [
        (i, float(m), *map(float, p), *map(float, v))
        for i, (m, p, v) in enumerate(zip(masses, pos, vel))
    ]
    while len(objs) > 1:
        best, best_a = None, math.inf
        for i, j in combinations(range(len(objs)), 2):
            _, ma, xa, ya, za, ua, va, wa = objs[i]
            _, mb, xb, yb, zb, ub, vb, wb = objs[j]
            dx, dy, dz = xb - xa, yb - ya, zb - za
            du, dv, dw = ub - ua, vb - va, wb - wa
            mu = G * (ma + mb)
            energy = 0.5 * (du * du + dv * dv + dw * dw) - mu / math.sqrt(
                dx * dx + dy * dy + dz * dz
            )
            if energy < 0.0:
                a = -mu / (2.0 * energy)
                if a < best_a:
                    best, best_a = (i, j), a
        if best is None:
            break
        i, j = best
        ta, ma, *sa = objs[i]
        tb, mb, *sb = objs[j]
        m = ma + mb
        merged = ((ta, tb), m, *((ma * qa + mb * qb) / m for qa, qb in zip(sa, sb)))
        objs = [o for k, o in enumerate(objs) if k not in (i, j)] + [merged]
    return tuple(sorted((canonical(o[0]) for o in objs), key=_lowest))


def label(trees: tuple) -> str:
    """Short human-readable classification, e.g. '2+2', '3+1', 'triple + single'; for five or
    more stars the sizes of the two top-level subsystems, e.g. '3+2' or '4+1'."""
    if len(trees) == 1 and _size(trees[0]) >= 4:
        left, right = trees[0]
        a, b = sorted((_size(left), _size(right)), reverse=True)
        return f"{a}+{b}"
    parts = sorted((_size(t) for t in trees), reverse=True)
    names = {1: "single", 2: "binary", 3: "triple", 4: "quadruple", 5: "quintuple", 6: "sextuple"}
    return " + ".join(names.get(p, f"{p}-star system") for p in parts)


def to_string(trees: tuple) -> str:
    """Star numbers 1..N in brackets, e.g. '[[1,2],[3,4]]'; unbound parts joined by ' | '."""

    def s(t):
        return str(t + 1) if isinstance(t, int) else f"[{s(t[0])},{s(t[1])}]"

    return " | ".join(s(t) for t in trees)


def classify_change(old: tuple, new: tuple) -> str:
    """Name of a transition between two hierarchies."""
    if len(new) > len(old):
        # a star or a binary became unbound
        return "ejection" if any(_size(t) == 1 for t in new) and len(new) == 2 else "breakup"
    if len(new) < len(old):
        return "capture"
    if label(old) == label(new):
        return "exchange"
    return "reconfiguration"


def from_node(node) -> object:
    """Convert a swingersim.hierarchy.Node tree into the nested-tuple form."""
    if isinstance(node, int):
        return node
    return canonical((from_node(node.left), from_node(node.right)))
