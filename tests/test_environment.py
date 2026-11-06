import math
import warnings

import numpy as np
import pytest

from swingersim import environment as env
from swingersim.config import Config


def test_tide_time_scale_and_emax():
    # Grishin & Perets (2022) eq. 7: 0.89 Gyr for 1 Msun at 1e4 AU
    assert env.tide_time_scale_yr(1.0, 1.0e4) == pytest.approx(0.89e9, rel=0.02)
    assert env.tide_emax(90.0) == pytest.approx(1.0)
    assert env.tide_emax(20.0) == 0.0  # below the critical inclination (cos^2 i > 4/5)


def test_parabolic_kick_matches_heggie_rasio_scalar():
    # eq. 10 of Hamers & Samsing (2019) for E = 1: De = +(15 pi / 4) eps e sqrt(1-e^2) sin 2w sin^2 i
    e, inc, om = 0.4, math.radians(60), math.radians(25)
    ehat = np.array([math.cos(om), math.sin(om) * math.cos(inc), math.sin(om) * math.sin(inc)])
    jhat = np.array([0.0, -math.sin(inc), math.cos(inc)])
    de, _ = env.first_order_kick(e * ehat, math.sqrt(1 - e * e) * jhat, 1.0, 2.0, 1.0, 30.0, 1.0)
    eps = math.sqrt(1.0 / (2.0 * 3.0) * (1.0 / 30.0) ** 3) * 2.0**-1.5
    expected = (
        15 * math.pi / 4 * eps * e * math.sqrt(1 - e * e) * math.sin(2 * om) * math.sin(inc) ** 2
    )
    assert ehat @ de == pytest.approx(expected, rel=1e-9)


def test_encounter_sequence_rate():
    spec = Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 1, "m2": 1, "m3": 1, "m4": 1},
            "orbits": {"A": {"P": 2}, "B": {"P": 3}, "out": {"P": 1e5}},
            "physics": {"encounters": {"n_star_pc3": 0.1, "sigma_kms": 30, "b_max_au": 1e4}},
            "run": {"t_max_yr": 1e9},
        }
    ).physics.encounters
    encs = env.encounter_sequence(spec, 4.0, 100.0, 1e10)
    n = 0.1 / env.PC_AU**3
    rate = n * math.pi * 1e8 * 30 * env.KMS_AUYR * math.sqrt(8 / math.pi)
    assert len(encs) == pytest.approx(rate * 1e10, rel=0.15)
    assert all(e.ecc >= 1.0 and e.q_au <= 1e4 for e in encs)


def test_encounters_in_secular_run():
    from swingersim.backends.secular import run_secular

    cfg = Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"out": ["A", "m3"], "A": ["m1", "m2"]},
            "masses": {"m1": 1.0, "m2": 1.0, "m3": 1.0},
            "orbits": {"A": {"a": 1.0}, "out": {"a": 3000.0, "e": 0.3, "i": 40.0}},
            "physics": {"gr": False, "encounters": {"n_star_pc3": 1.0, "seed": 3}},
            "run": {"backend": "secular", "t_max_yr": 1e8, "n_output": 101},
        }
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = run_secular(cfg)
    assert r.system["encounters_n"] > 10
    e = r.orbits["out"]["e"]
    assert np.all((e >= 0) & (e < 1))
    assert np.std(np.diff(e)) > 0  # the kicks changed the outer orbit


def test_encounters_and_tide_in_nbody_run():
    from swingersim.backends.nbody import run_nbody

    cfg = Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"out": ["A", "m3"], "A": ["m1", "m2"]},
            "masses": {"m1": 1.0, "m2": 1.0, "m3": 1.0},
            "orbits": {"A": {"a": 50.0}, "out": {"a": 3000.0, "e": 0.3, "i": 40.0}},
            "physics": {
                "gr": False,
                "galactic_tide": {},
                "encounters": {"n_star_pc3": 50.0, "seed": 2, "b_max_au": 3e4},
            },
            "run": {"backend": "nbody", "t_max_yr": 3e5, "n_output": 31},
        }
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = run_nbody(cfg)
    assert r.system["encounters_n"] >= 1
    assert np.all(np.isfinite(r.orbits["A"]["e"][: len(r.times)]))
    assert r.energy_error.max() < 1e-6


def test_ensemble_members_get_own_encounter_sequence():
    from swingersim.config import Config
    from swingersim.ensemble import draw_members

    cfg = Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"out": ["m1", "m2"]},
            "masses": {"m1": 1.0, "m2": 0.5},
            "orbits": {"out": {"a": 5000.0, "e": 0.3}},
            "physics": {"gr": False, "encounters": {"seed": 7}},
            "run": {"backend": "secular", "t_max_yr": 1.0e6},
            "ensemble": {"n": 3, "vary": {"orbits.out.e": {"uniform": [0.1, 0.5]}}},
        }
    )
    seeds = [m.config["physics"]["encounters"]["seed"] for m in draw_members(cfg)]
    assert seeds == [7, 8, 9]
