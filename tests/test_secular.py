"""Secular 2+2 backend: internal consistency and agreement with the N-body backend."""

from __future__ import annotations

import numpy as np
import pytest

from swingersim.backends.nbody import run_nbody
from swingersim.backends.secular import build_secular, run_secular
from swingersim.config import Config, load_config
from swingersim.secular.twotwo import elements_from_vectors, vectors_from_elements

from .conftest import EXAMPLES


def test_vector_element_round_trip():
    e, inc, Om, om = 0.3, 1.1, 2.0, 4.0
    K, p = vectors_from_elements(e, inc, Om, om)
    assert K @ p == pytest.approx(0.0, abs=1e-14)
    assert K @ K + p @ p == pytest.approx(1.0)
    assert np.allclose(elements_from_vectors(K, p), (e, inc, Om, om))


@pytest.mark.parametrize("tides", [False, True])
def test_gradients_match_finite_differences(tides):
    cfg = load_config(EXAMPLES / "example_2p2_tides.yaml")
    cfg.physics.tides = tides
    cfg.physics.rotation = tides
    cfg.masses["m2"] = 0.7  # unequal masses: octupole terms on
    sys_, y0, *_ = build_secular(cfg)
    rng = np.random.default_rng(1)
    y = y0 + 1e-3 * rng.standard_normal(len(y0))  # off the constraint surface is fine for H
    H, g = sys_.energy_and_gradients(sys_.unpack(y))
    names = ("KA", "pA", "KB", "pB", "K", "p")
    h = 1e-7
    for k, name in enumerate(names):
        for c in range(3):
            yp, ym = y.copy(), y.copy()
            yp[3 * k + c] += h
            ym[3 * k + c] -= h
            num = (
                sys_.energy_and_gradients(sys_.unpack(yp))[0]
                - sys_.energy_and_gradients(sys_.unpack(ym))[0]
            ) / (2 * h)
            assert g[name][c] == pytest.approx(num, rel=1e-5, abs=1e-12 * abs(H))


def test_energy_and_constraints_conserved():
    cfg = load_config(EXAMPLES / "example_2p2_zlk.yaml")
    cfg.run.backend = "secular"
    res = run_secular(cfg)
    assert np.max(res.energy_error) < 1e-8
    assert np.max(res.angmom_error) < 1e-12
    assert res.backend == "secular"


def test_secular_agrees_with_nbody():
    cfg = load_config(EXAMPLES / "example_2p2.yaml")
    nb = run_nbody(cfg)
    cfg.run.backend = "secular"
    se = run_secular(cfg)
    for key in ("out-A", "out-B", "A-B"):
        for f in (np.min, np.max):
            assert np.degrees(f(se.mutual[key])) == pytest.approx(
                np.degrees(f(nb.mutual[key])), abs=0.5
            )


def test_dissipation_conserves_angular_momentum():
    cfg = load_config(EXAMPLES / "example_2p2_tides.yaml")
    cfg.run.backend = "secular"
    cfg.physics.rotation = False  # spins start at rest and are spun up by the tides
    cfg.stars["m1"].lag_s = 100.0
    cfg.run.t_max_yr, cfg.run.n_output = 20.0, 40
    res = run_secular(cfg)
    assert res.dissipative
    assert res.spins["m1"]["P_rot_d"][-1] < 1e3  # spun up from zero
    # the aligned-spin approximation transfers no angular momentum, so this stays small
    assert np.max(res.angmom_error) < 1e-6


def test_kernel_matches_reference_implementation():
    """The compiled kernel (swingersim.secular.kernel) reproduces energy_and_gradients."""
    from swingersim.backends.secular import build_secular
    from swingersim.secular.twotwo import _cross

    for tides in (False, True):
        cfg = Config.model_validate(
            {
                "topology": "2+2",
                "masses": {"m1": 1.3, "m2": 0.9, "m3": 1.1, "m4": 0.7},
                "orbits": {
                    "A": {"P": 3.0, "e": 0.2, "i": 40, "Omega": 10, "omega": 30},
                    "B": {"P": 1.5, "e": 0.1, "i": 70, "Omega": 100, "omega": 80},
                    "out": {"P": 200.0, "e": 0.4, "i": 5, "Omega": 0, "omega": 200},
                },
                "physics": {"gr": True, "tides": tides},
                "run": {"t_max_yr": 1.0, "backend": "secular"},
            }
        )
        sys_, y, *_ = build_secular(cfg)
        v = sys_.unpack(y)
        _h, g = sys_.energy_and_gradients(v)
        ref = np.empty(18)
        for sk, sp, mu, m, a, sl in (
            ("KA", "pA", sys_.muA, sys_.MA, v["aA"], 0),
            ("KB", "pB", sys_.muB, sys_.MB, v["aB"], 6),
            ("K", "p", sys_.mu, sys_.MAB, sys_.a_out, 12),
        ):
            L = sys_._L(mu, m, a)
            ref[sl : sl + 3] = -(_cross(v[sp], g[sp]) + _cross(v[sk], g[sk])) / L
            ref[sl + 3 : sl + 6] = -(_cross(v[sk], g[sp]) + _cross(v[sp], g[sk])) / L
        new = sys_.rhs(0.0, y)[:18]
        assert np.max(np.abs(new - ref)) <= 1e-12 * np.max(np.abs(ref))


def test_nested_backend_equals_twotwo_kernel():
    """The general nested secular equations reduce exactly to the 2+2 ones."""
    from swingersim.hierarchy import build_tree
    from swingersim.secular.nested import SecularNested

    cfg = Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 1.3, "m2": 0.9, "m3": 1.1, "m4": 0.7},
            "orbits": {
                "A": {"P": 3.0, "e": 0.2, "i": 40, "Omega": 10, "omega": 30},
                "B": {"P": 1.5, "e": 0.1, "i": 70, "Omega": 100, "omega": 80},
                "out": {"P": 200.0, "e": 0.4, "i": 5, "Omega": 0, "omega": 200},
            },
            "physics": {"gr": True},
            "run": {"t_max_yr": 1.0, "backend": "secular"},
        }
    )
    sys2, y, *_ = build_secular(cfg)
    root, masses = build_tree(cfg)
    nest = SecularNested(root, masses, gr=True)
    order = ["A", "B", "out"]
    yn = np.concatenate([y[6 * order.index(n) : 6 * order.index(n) + 6] for n in nest.names])
    rn = nest.rhs(0.0, yn)
    back = np.concatenate(
        [rn[6 * nest.names.index(n) : 6 * nest.names.index(n) + 6] for n in order]
    )
    ref = sys2.rhs(0.0, y)[:18]
    assert np.max(np.abs(back - ref)) <= 1e-12 * np.max(np.abs(ref))


def test_secular_3p1_follows_nbody():
    cfg = load_config(EXAMPLES / "example_3p1.yaml")
    cfg.run.t_max_yr = 60.0
    nb = run_nbody(cfg)
    cfg.run.backend = "secular"
    sec = run_secular(cfg)
    assert sec.backend == "secular" and set(sec.orbits) == {"inner", "mid", "out"}
    assert np.nanmax(sec.energy_error) < 1e-8 and np.nanmax(sec.angmom_error) < 1e-12
    # ZLK cycle of the inner binary: the same amplitude within the mild-hierarchy corrections
    assert abs(sec.extremes["inner"]["e_max"] - nb.extremes["inner"]["e_max"]) < 0.04
    for pair in ("mid-inner", "out-mid"):
        assert abs(np.degrees(np.max(sec.mutual[pair]) - np.max(nb.mutual[pair]))) < 2.0
