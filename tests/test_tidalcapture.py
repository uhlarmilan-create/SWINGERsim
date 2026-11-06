import math
import warnings

import pytest

from swingersim import tidalcapture as tc
from swingersim import units
from swingersim.scatter import ScatterSpec, run_encounter


def test_fit_values_and_monotonic_decline():
    # Portegies Zwart & Meinen (1993), Table 1: log T_2(eta = 1) = A for n = 3
    assert math.log10(tc.t_l(1.0, 3.0, 2)) == pytest.approx(-1.124)
    assert tc.t_l(30.0, 3.0, 2) < tc.t_l(10.0, 3.0, 2) < tc.t_l(3.0, 3.0, 2)


def _encounter(q_over_r):
    r = 0.65
    spec = ScatterSpec.model_validate(
        {
            "target": {"masses": [0.6], "radii_rsun": [r], "polytrope_n": [3.0]},
            "projectile": {"masses": [0.6], "radii_rsun": [r], "polytrope_n": [3.0]},
            "v_inf_kms": 10.0,
            "b_au": 1.0,
            "tidal_capture": True,
            "n": 1,
            "t_max_yr": 5.0,
        }
    )
    v = 10.0 * units.SECONDS_PER_YEAR / units.AU_KM
    q = q_over_r * r * units.RSUN_AU
    b = math.sqrt(q * q + 2 * units.G * 1.2 * q / v**2)
    enc = {"id": 0, "b_au": b, "q_au": q, "phase_t": 0.0, "phase_p": 0.0}
    from swingersim.scatter import draw_encounters

    enc = {**draw_encounters(spec)[0], "b_au": b, "q_au": q}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return run_encounter(spec, enc)


def test_capture_inside_and_flyby_outside_the_capture_radius():
    r = 0.65 * units.RSUN_AU
    v = 10.0 * units.SECONDS_PER_YEAR / units.AU_KM
    rc = tc.capture_radius(0.6, r, 3.0, 0.6, r, 3.0, v) / r
    inside, outside = _encounter(0.97 * rc), _encounter(1.05 * rc)
    assert inside["outcome"] == "bound binary" and inside["tidal_events"] >= 1
    assert outside["outcome"] == "flyby"
    assert outside["energy_error"] < 1e-8
