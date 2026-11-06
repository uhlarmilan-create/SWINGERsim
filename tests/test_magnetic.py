import math
import warnings

import pytest

from swingersim import magnetic as mg
from swingersim.config import Config


def test_turnover_time_and_saturation():
    assert mg.tau_cz_days(5772.0) == pytest.approx(12.9, abs=0.05)  # Matt et al. 2015 Table 1
    assert mg.tau_cz_days(7500.0) < 0.1  # no convective envelope
    assert mg.omega_sat(5772.0) == pytest.approx(10 * mg.OMEGA_SUN, rel=0.01)


def test_spin_down_regimes():
    rg2 = 0.109
    # saturated: exponential; unsaturated: the Skumanich power law
    w0 = 50 * mg.OMEGA_SUN
    t = 1e7
    inertia = rg2 * mg.MSUN_G * mg.RSUN_CM**2
    rate = mg.T0_SUN * 100 / (inertia * mg.OMEGA_SUN)
    assert mg.spin_down(w0, 1, 1, 5772, rg2, t) == pytest.approx(
        w0 * math.exp(-rate * t * mg.SEC_PER_YR), rel=1e-9
    )
    w = mg.spin_down(mg.OMEGA_SUN, 1, 1, 5772, rg2, 1e9)
    w2 = mg.spin_down(w, 1, 1, 5772, rg2, 1e9)  # composition of two steps = one step
    assert w2 == pytest.approx(mg.spin_down(mg.OMEGA_SUN, 1, 1, 5772, rg2, 2e9), rel=1e-9)
    assert mg.spin_down(-w0, 1, 1, 5772, rg2, t) < 0  # retrograde spins keep their sign


def test_applegate_inverse():
    k = mg.applegate_max_oc(0.18, 0.22, 3400, 0.69, 0.0868, 19.4)
    dpp = mg.dp_over_p_from_oc(k, 19.4)
    e = mg.applegate_energy_ratio(0.18, 0.22, 3400, 0.69, 0.0868, dpp, 19.4)
    assert e == pytest.approx(1.0, rel=1e-9)


def test_config_needs_rotation():
    base = {
        "topology": "2+2",
        "masses": {"m1": 1, "m2": 1, "m3": 1, "m4": 1},
        "orbits": {"A": {"P": 2}, "B": {"P": 3}, "out": {"P": 1e5}},
        "run": {"t_max_yr": 10},
    }
    with pytest.raises(ValueError):
        Config.model_validate({**base, "physics": {"magnetic_braking": {}}})


def test_radius_inflation():
    from swingersim.stellar import star_properties

    def cfg(x):
        return Config.model_validate(
            {
                "topology": "2+2",
                "masses": {"m1": 0.3, "m2": 0.3, "m3": 0.3, "m4": 0.3},
                "orbits": {"A": {"P": 2}, "B": {"P": 3}, "out": {"P": 1e5}},
                "stars": {"m1": {"radius_inflation": x}},
                "stellar_model": {"kind": "sse", "age_myr": 1000},
                "run": {"t_max_yr": 10},
            }
        )

    p0, p1 = star_properties(cfg(0.0), "m1", 0.3), star_properties(cfg(0.05), "m1", 0.3)
    assert p1.radius_rsun == pytest.approx(1.05 * p0.radius_rsun)
    assert p1.luminosity_lsun == pytest.approx(p0.luminosity_lsun)
    assert p1.teff_k < p0.teff_k


@pytest.mark.slow
def test_braking_shrinks_close_binary():
    from swingersim.backends.secular import run_secular

    def cfg(mb):
        d = {
            "topology": "2+2",
            "masses": {"m1": 1.0, "m2": 0.9, "m3": 1.0, "m4": 0.9},
            "orbits": {"A": {"P": 2.0}, "B": {"P": 2.5}, "out": {"P": 3.65e6, "e": 0.1, "i": 20}},
            "stellar_model": {"kind": "sse", "age_myr": 100.0, "z": 0.02},
            "physics": {
                "tides": True,
                "rotation": True,
                "tidal_lag": "hurley2002",
                "evolution": True,
            },
            "run": {"backend": "secular", "t_max_yr": 2e7, "n_output": 5},
        }
        if mb:
            d["physics"]["magnetic_braking"] = {}
        return Config.model_validate(d)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r0, r1 = run_secular(cfg(False)), run_secular(cfg(True))
    assert r1.orbits["A"]["a"][-1] < r0.orbits["A"]["a"][-1]
    # the spin stays (nearly) synchronised: the braking torque drains the orbit
    assert r1.spins["m1"]["P_rot_d"][-1] == pytest.approx(r0.spins["m1"]["P_rot_d"][-1], rel=0.01)


def test_activity_cycle_in_eclipse_timings():
    import warnings

    import numpy as np

    from swingersim.config import Config
    from swingersim.magnetic import activity_report
    from swingersim.timing import eclipse_times, o_minus_c, sinusoid_amplitude

    cfg = Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"out": ["A", "m3"], "A": ["m1", "m2"]},
            "masses": {"m1": 0.93, "m2": 0.84, "m3": 0.3},
            "orbits": {"A": {"P": 0.5212, "i": 89}, "out": {"P": 200000.0}},
            "stars": {"m1": {"R": 0.96, "Teff": 5040}},
            "physics": {"gr": False, "activity": {"m1": {"p_mod_yr": 1.0, "amplitude_s": 60.0}}},
            "run": {"t_max_yr": 3.0},
        }
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ecl, p = eclipse_times(cfg, 3.0)
    rows = [r for r in o_minus_c(ecl, p) if r["kind"] == "primary"]
    t = np.array([r["t_yr"] for r in rows])
    oc = np.array([r["o_c_d"] for r in rows]) * 86400.0
    assert sinusoid_amplitude(t, oc, 1.0) == pytest.approx(60.0, rel=0.02)
    rep = activity_report(cfg)[0]
    assert rep["star"] == "m1" and rep["dE_over_E"] > 0
