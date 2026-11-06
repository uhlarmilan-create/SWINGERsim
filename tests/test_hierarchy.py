"""Validation test 2: elements -> vectors -> elements for both topologies."""

import numpy as np
import pytest

from swingersim.hierarchy import build_tree, initialize, node_elements


def _angle_diff(x, y):
    return abs((x - y + np.pi) % (2 * np.pi) - np.pi)


@pytest.mark.parametrize("name", ["example_2p2", "example_3p1"])
def test_roundtrip(name, request):
    cfg = request.getfixturevalue(name)
    root, masses = build_tree(cfg)
    pos, vel = initialize(root, masses)

    # barycentric frame
    assert np.allclose(masses @ pos, 0, atol=1e-14)
    assert np.allclose(masses @ vel, 0, atol=1e-12)

    for node in root.nodes():
        back = node_elements(node, pos, vel, masses)
        el = node.elements
        assert back.a == pytest.approx(el.a, rel=1e-12), node.name
        assert back.e == pytest.approx(el.e, abs=1e-12), node.name
        for k in ("inc", "Omega", "omega", "M"):
            assert _angle_diff(getattr(back, k), getattr(el, k)) < 1e-10, (node.name, k)


def test_tree_shapes(example_2p2, example_3p1):
    r22, _ = build_tree(example_2p2)
    assert [n.name for n in r22.nodes()] == ["out", "A", "B"]
    r31, _ = build_tree(example_3p1)
    assert [n.name for n in r31.nodes()] == ["out", "mid", "inner"]
    assert r31.stars() == [0, 1, 2, 3]
