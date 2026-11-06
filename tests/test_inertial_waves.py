"""Inertial-wave tidal dissipation (Ogilvie 2013; Barker 2022)."""

import numpy as np
import pytest

from swingersim.config import Config
from swingersim.inertial_waves import (
    FULLY_CONVECTIVE_C,
    barker_profile,
    ogilvie_d,
    polytrope_c,
    structure_c,
    two_layer,
    two_layer_profile,
)
from swingersim.tides import inertial_wave_k_over_t, inertial_wave_q


def test_numerical_method_reproduces_two_layer_formula():
    for alpha, beta in ((0.71, 0.98), (0.6, 0.95)):
        prof = two_layer_profile(alpha, beta)
        assert barker_profile(*prof) == pytest.approx(two_layer(alpha, beta), rel=1.5e-2)


def test_formula_limits_and_polytrope():
    assert ogilvie_d(0.0, 0.5) == 0.0 and ogilvie_d(0.7, 1.0) == 0.0
    # solar values (alpha 0.71, beta 0.98): <D>/eps^2 = 7.66e-4
    assert ogilvie_d(0.71, 0.98) == pytest.approx(7.66e-4, rel=2e-3)
    assert polytrope_c(1.5) == pytest.approx(FULLY_CONVECTIVE_C, rel=5e-3)
    assert structure_c(0.0, 0.0) == FULLY_CONVECTIVE_C


def test_q_scaling_and_frequency_window():
    m, r, menv, renv = 1.0, 1.0, 0.03, 0.27
    w10 = 2 * np.pi / (10 / 365.25)
    q10 = inertial_wave_q(m, r, menv, renv, w10)
    assert 5e6 < q10 < 5e7
    assert inertial_wave_q(m, r, menv, renv, 2 * w10) == pytest.approx(q10 / 4)
    assert inertial_wave_q(m, r, 0.0, 0.0, w10) == np.inf
    # a = 0.05 AU twin binary: n ~ 2 pi / 3.6 d; spin slower than n/2 -> no inertial waves
    a = 0.05
    n = np.sqrt(4 * np.pi**2 * 2.0 / a**3)
    assert inertial_wave_k_over_t(m, r, menv, renv, 1.0, a, 0.4 * n) == 0.0
    assert inertial_wave_k_over_t(m, r, menv, renv, 1.0, a, n) > 0.0


def test_config_needs_hurley():
    base = {
        "topology": "2+2",
        "masses": {"m1": 1.0, "m2": 1.0, "m3": 1.0, "m4": 1.0},
        "orbits": {"A": {"P": 3.0}, "B": {"P": 4.0}, "out": {"P": 2000.0}},
        "run": {"t_max_yr": 10.0},
    }
    with pytest.raises(ValueError):
        Config.model_validate({**base, "physics": {"tides": True, "inertial_waves": True}})
    c = Config.model_validate(
        {
            **base,
            "physics": {"tides": True, "tidal_lag": "hurley2002", "inertial_waves": True},
        }
    )
    assert c.physics.inertial_waves
