import math
import warnings

import numpy as np
import pytest

from swingersim import masstransfer as mt
from swingersim import units
from swingersim.config import Config


def _binary(m1, m2, p_day, age, t_max, kind="mist", **mtspec):
    sm = {"kind": kind, "age_myr": age}
    if kind == "mist":
        sm["vvcrit"] = 0.0
    return Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"A": ["m1", "m2"]},
            "masses": {"m1": m1, "m2": m2},
            "orbits": {"A": {"P": p_day}},
            "stellar_model": sm,
            "physics": {"gr": False, "evolution": True, "mass_transfer": mtspec},
            "run": {"t_max_yr": t_max, "backend": "secular"},
        }
    )


def test_orbit_modes_limits():
    # conservative: a (Md Ma)^2 const; alpha = 1 (Jeans mode): a M const
    md, ma, a = mt.orbit_path(3.0, 2.0, 1.0, 0.0, 0.0, 1.0, n=20)
    assert np.allclose(a * (md * ma) ** 2, 36.0, rtol=1e-8)
    md, ma, a = mt.orbit_path(3.0, 2.0, 1.0, 1.0, 0.0, 1.0, n=20)
    assert np.allclose(ma, 2.0)
    assert np.allclose(a * (md + ma), 5.0, rtol=1e-8)


def test_bse_rate_and_qcrit():
    assert mt.bse_rate(1.0, 1.0, 1.1) == 0.0
    assert mt.bse_rate(6.0, math.e, 1.0) == pytest.approx(3e-6 * 25.0)
    assert mt.q_crit(2, 3.0, 0.3) == 4.0
    assert mt.q_crit(1, 0.5, 0.0) == 0.695
    assert mt.q_crit(1, 2.0, 0.0) == 3.0  # zeta_ad = 4.72 (Hurley et al. 2002)
    assert mt.q_crit(4, 2.0, 0.3) == 3.0
    assert mt.q_crit(0, 1.0, 0.0) is None
    assert mt.q_crit(3, 1.0, 0.0) == pytest.approx((1.67 - 0.3) / 2.13)


def test_case_b_stable_transfer_reverses_mass_ratio():
    cfg = _binary(3.0, 2.0, 3.0, 355.0, 2.0e7)
    a0 = (units.G * 5.0 * (3.0 / units.DAYS_PER_YEAR) ** 2 / (4 * math.pi**2)) ** (1 / 3)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        h = mt.binary_history(cfg, (0, 1), [3.0, 2.0], a0, 2.0e7)
    assert h.donor == 0 and h.t_onset is not None
    assert h.t_reversal is not None and h.t_reversal > h.t_onset
    # the donor, re-modelled for its lower mass, ends as a giant with a thin envelope that
    # transfers on the nuclear time scale (Algol phase)
    assert h.event in (None, "donor_stripped")
    assert h.m[0, -1] < 1.0 and h.m[1, -1] > 4.0
    assert h.m_core[0, -1] < h.m[0, -1]
    assert h.m_ev[1, -1] == pytest.approx(h.m[1, -1])  # the accretor follows its new mass
    assert h.age[1, -1] < 355.0  # and is rejuvenated (same fractional age)
    assert h.a[-1] > 3 * a0  # the orbit widens after the reversal


def test_giant_donor_with_high_q_goes_common_envelope():
    cfg = _binary(3.0, 2.0, 10.0, 355.0, 2.0e7)
    a0 = (units.G * 5.0 * (10.0 / units.DAYS_PER_YEAR) ** 2 / (4 * math.pi**2)) ** (1 / 3)
    h = mt.binary_history(cfg, (0, 1), [3.0, 2.0], a0, 2.0e7)
    assert h.event == "common_envelope"
    assert h.m[0, -1] == pytest.approx(3.0, abs=0.01)


def test_wide_binary_stays_detached():
    cfg = _binary(1.0, 0.8, 1.0e5, 1000.0, 1.0e8)
    a0 = (units.G * 1.8 * (1.0e5 / units.DAYS_PER_YEAR) ** 2 / (4 * math.pi**2)) ** (1 / 3)
    h = mt.binary_history(cfg, (0, 1), [1.0, 0.8], a0, 1.0e8)
    assert h.t_onset is None and h.event is None


def test_scale_pair_keeps_centre_and_eccentricity():
    import rebound

    from swingersim.backends.nbody import _scale_pair

    sim = rebound.Simulation()
    sim.G = units.G
    sim.add(m=2.0)
    sim.add(m=1.0, a=0.1, e=0.3, omega=0.4, f=1.0)
    sim.add(m=0.5, a=3.0, e=0.1, f=2.0)
    sim.move_to_com()
    o0 = sim.particles[1].orbit(primary=sim.particles[0])
    com0 = np.array(sim.com().xyz)
    _scale_pair(sim, 0, 1, 1.7)
    o1 = sim.particles[1].orbit(primary=sim.particles[0])
    assert o1.a == pytest.approx(1.7 * o0.a, rel=1e-10)
    assert o1.e == pytest.approx(o0.e, rel=1e-10)
    assert o1.f == pytest.approx(o0.f, rel=1e-10)
    assert np.allclose(np.array(sim.com().xyz), com0, atol=1e-14)


def test_triple_secular_run_follows_history():
    from swingersim.backends.secular_nested import run_secular_nested

    cfg = Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"out": ["A", "m3"], "A": ["m1", "m2"]},
            "masses": {"m1": 3.0, "m2": 2.0, "m3": 1.5},
            "orbits": {"A": {"P": 3.0}, "out": {"P": 2000.0, "e": 0.2, "i": 30.0}},
            "stellar_model": {"kind": "mist", "age_myr": 358.5, "vvcrit": 0.0},
            "physics": {"gr": False, "evolution": True, "mass_transfer": {}},
            "run": {"backend": "secular", "t_max_yr": 1.0e7, "n_output": 101},
        }
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = run_secular_nested(cfg)
    h = mt.histories(cfg)["A"]
    if h.event in (None, "track_end"):
        assert res.event is None
    else:
        assert res.event.startswith(f"mass_transfer: {h.event}")
    a = res.orbits["A"]["a"]
    assert a[-1] == pytest.approx(h.a_at(res.times[-1]), rel=1e-9)
    # conservative transfer: the outer orbit keeps its size
    assert res.orbits["out"]["a"][-1] == pytest.approx(res.orbits["out"]["a"][0], rel=1e-6)
    assert res.system["mass_transfer_A_donor"] == "m1"
