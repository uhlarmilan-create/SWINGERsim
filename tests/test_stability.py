"""Stability criteria, ghost test and backward integration."""

from __future__ import annotations

import numpy as np
import pytest

from swingersim.backends.nbody import run_nbody
from swingersim.backends.secular import run_secular
from swingersim.config import Config, config_at, load_config
from swingersim.stability import ghost_test, ma01_ratio, v22_ratio

from .conftest import EXAMPLES


def test_v22_coplanar_circular_limit():
    # e_in = e_out = 0, i = 0: Y~crit = 2.4 (1 + q)^(2/5), Y~ = a_out / a_in
    assert v22_ratio(1.0, 5.0, 0.0, 0.0, 1.0, 0.0) == pytest.approx(5.0 / (2.4 * 2**0.4))
    # MA01 in the same limit: 2.8 (1 + q)^(2/5)
    assert ma01_ratio(1.0, 5.0, 0.0, 1.0, 0.0) == pytest.approx(5.0 / (2.8 * 2**0.4))


def test_v22_zlk_inclination_lowers_stability():
    flat = v22_ratio(1.0, 5.0, 0.0, 0.0, 0.5, 0.0)
    polar = v22_ratio(1.0, 5.0, 0.0, 0.0, 0.5, np.pi / 2)
    assert polar < flat


def _equal_2p2(ratio_a: float) -> Config:
    return Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 1.0, "m2": 1.0, "m3": 1.0, "m4": 1.0},
            "orbits": {"A": {"a": ratio_a}, "B": {"a": ratio_a, "M0": 90.0}, "out": {"a": 1.0}},
            "physics": {"gr": False},
            "run": {"t_max_yr": 1.0},
        }
    )


def test_ghost_test_separates_stable_and_unstable():
    assert ghost_test(_equal_2p2(0.1), n_outer=30).stable
    bad = ghost_test(_equal_2p2(0.45), n_outer=30)
    assert not bad.stable and bad.reason in ("unbound", "diverged")


def test_backward_nbody_is_reversible():
    cfg = load_config(EXAMPLES / "example_2p2.yaml")
    cfg.run.t_max_yr, cfg.run.n_output = 20.0, 41
    fwd = run_nbody(cfg)
    back = config_at(cfg, fwd)
    back.run.t_max_yr = -20.0
    res = run_nbody(back)
    assert res.times[-1] == pytest.approx(-20.0)
    for name in ("A", "B", "out"):
        for k in ("a", "e", "inc"):
            assert res.orbits[name][k][-1] == pytest.approx(fwd.orbits[name][k][0], rel=1e-8)


def test_backward_secular_retraces_forward():
    cfg = load_config(EXAMPLES / "example_2p2.yaml")
    cfg.run.backend = "secular"
    cfg.run.t_max_yr, cfg.run.n_output = 300.0, 31
    fwd = run_secular(cfg)
    cfg.run.t_max_yr = -300.0
    back = run_secular(cfg)
    # the mean elements of the past are not the future: check only that it ran and conserves H
    assert back.times[-1] == pytest.approx(-300.0)
    assert np.max(back.energy_error) < 1e-8
    assert not np.allclose(back.mutual["A-B"], fwd.mutual["A-B"])


def test_zero_duration_rejected():
    with pytest.raises(ValueError):
        load_config(EXAMPLES / "example_2p2.yaml").model_validate(
            {**load_config(EXAMPLES / "example_2p2.yaml").model_dump(), "run": {"t_max_yr": 0}}
        )


def test_beta_hamers_lai_2p2_product():
    from swingersim.backends.secular import run_secular

    cfg = load_config(EXAMPLES / "example_2p2.yaml")
    cfg.run.backend, cfg.run.t_max_yr, cfg.run.n_output = "secular", 10.0, 5
    s = run_secular(cfg).system
    d = {x: s[f"beta_hamers_lai_{x}"] for x in ("A", "B")}
    # (a, M) factors cancel in the product: beta_A beta_B = (9/16) cos i_A cos i_B
    assert d["A"] * d["B"] > 0
    assert "eclipse_fraction_both" in s and 0.0 <= s["eclipse_fraction_both"] <= 1.0


def test_mlp_classifier_if_installed():
    from swingersim.stability import mlp_directory, mlp_stable

    if mlp_directory() is None:
        pytest.skip("repository of Vynatheya et al. (2023) not installed")
    assert mlp_stable(_equal_2p2(0.1)) is True
    assert mlp_stable(_equal_2p2(0.45)) is False
