"""Stellar mass loss in the dynamics (swingersim.massloss, physics.mass_loss)."""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from scipy.integrate import quad

from swingersim import massloss
from swingersim.backends.nbody import run_nbody
from swingersim.backends.secular import run_secular
from swingersim.config import Config

warnings.filterwarnings("ignore", message="REBOUNDx")


def _two_two(physics: dict, backend: str = "secular", t_max: float = 1e6, n: int = 51) -> Config:
    return Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 1.2, "m2": 1.0, "m3": 0.9, "m4": 0.8},
            "orbits": {
                "A": {"P": 6.74, "e": 0.05, "i": 85.0},
                "B": {"P": 2.68, "e": 0.02, "i": 60.0, "Omega": 40.0},
                "out": {"P": 500, "e": 0.3, "i": 70.0, "Omega": 10.0},
            },
            "physics": physics,
            "run": {"backend": backend, "t_max_yr": t_max, "n_output": n},
        }
    )


RATES = {"gr": True, "mass_loss": {"source": "rates", "rates": {"m1": -4e-7, "m3": -2e-7}}}


def test_config_checks():
    with pytest.raises(ValueError, match="needs rates"):
        _two_two({"mass_loss": {"source": "rates"}})
    with pytest.raises(ValueError, match="needs physics.evolution"):
        _two_two({"mass_loss": {"source": "track"}})
    with pytest.raises(ValueError, match="unknown stars"):
        _two_two({"mass_loss": {"source": "rates", "rates": {"m9": -1e-6}}})


def test_masses_and_rates():
    cfg = _two_two(RATES)
    m = massloss.masses_at(cfg, [1.2, 1.0, 0.9, 0.8], 1e6)
    assert np.allclose(m, [0.8, 1.0, 0.7, 0.8])
    assert np.allclose(massloss.mass_rates(cfg, [1.2, 1.0, 0.9, 0.8], 0.0), [-4e-7, 0, -2e-7, 0])
    with pytest.raises(ValueError, match="no mass left"):
        massloss.masses_at(cfg, [1.2, 1.0, 0.9, 0.8], 1e7)


def test_analytic_helpers():
    assert massloss.impulsive_circular(0.5)[0] == float("inf")
    a, e = massloss.impulsive_circular(0.2)
    assert a == pytest.approx(0.8 / 0.6) and e == pytest.approx(0.25)
    assert massloss.adiabatic_a(1.0, 2.0, 1.5) == pytest.approx(4.0 / 3.0)
    # a uniform shell (rg2 = 2/3) leaves the spin rate unchanged: J ~ M
    assert massloss.wind_spin_factor(2.0, 1.0, 2.0 / 3.0) == pytest.approx(0.5)


def test_secular_2p2_adiabatic_invariant():
    res = run_secular(_two_two(RATES))
    st = res.stellar
    m_a = st["m1"]["mass_msun"] + st["m2"]["mass_msun"]
    m_b = st["m3"]["mass_msun"] + st["m4"]["mass_msun"]
    for name, mm in (("A", m_a), ("B", m_b), ("out", m_a + m_b)):
        am = res.orbits[name]["a"] * mm
        assert np.max(np.abs(am / am[0] - 1.0)) < 1e-12
    assert res.mass_loss
    assert res.system["mass_loss_psi_max_out"] < 1e-6
    assert np.max(res.energy_error) < 1e-8  # within each segment
    assert np.max(res.angmom_error) < 1e-10


def test_secular_nested_adiabatic_invariant():
    cfg = Config.model_validate(
        {
            "topology": "3+1",
            "masses": {"m1": 1.0, "m2": 0.8, "m3": 0.7, "m4": 0.6},
            "orbits": {
                "inner": {"P": 3.0, "e": 0.01, "i": 80.0},
                "mid": {"P": 60.0, "e": 0.2, "i": 40.0, "Omega": 30.0},
                "out": {"P": 3000, "e": 0.4, "i": 20.0, "Omega": 70.0},
            },
            "physics": {"gr": True, "mass_loss": {"source": "rates", "rates": {"m1": -2e-5}}},
            "run": {"backend": "secular", "t_max_yr": 2e4, "n_output": 51},
        }
    )
    res = run_secular(cfg)
    m1 = res.stellar["m1"]["mass_msun"]
    for name, mm in (("inner", m1 + 0.8), ("mid", m1 + 1.5), ("out", m1 + 2.1)):
        am = res.orbits[name]["a"] * mm
        assert np.max(np.abs(am / am[0] - 1.0)) < 1e-12


def test_secular_zlk_is_reparametrised():
    """Test-particle quadrupole ZLK with a perturber of changing mass: the flow is the
    constant-mass flow in the time tau = int c2(t)/c2(0) dt (c2 ~ m3 / a_mid^3)."""
    from swingersim.validation import _ml_up_crossings, _ml_zlk_triple

    r0 = run_secular(_ml_zlk_triple(0.0, "secular"))
    r1 = run_secular(_ml_zlk_triple(0.5, "secular"))
    t0 = _ml_up_crossings(r0.times, r0.orbits["inner"]["e"])
    t1 = _ml_up_crossings(r1.times, r1.orbits["inner"]["e"])

    def m3(t):
        return 1.0 - 0.5 * t / 30000.0

    def c2(t):  # m3 / a_mid^3 with a_mid (1 + m3) = const
        return m3(t) * (1.0 + m3(t)) ** 3 / 8.0

    for k in range(3):
        tau = quad(c2, 0.0, t1[k])[0]
        assert tau == pytest.approx(t0[k], rel=2e-3)


def test_secular_refuses_large_psi():
    cfg = _two_two({"mass_loss": {"source": "rates", "rates": {"m1": -2.0}}}, t_max=0.3)
    with pytest.raises(ValueError, match="mass-loss index"):
        run_secular(cfg)


def test_nbody_forward_backward_and_output(tmp_path):
    from swingersim.io import save_run

    rates = {"gr": True, "mass_loss": {"source": "rates", "rates": {"m1": -1e-3}}}
    cfg = _two_two(rates, backend="nbody", t_max=50.0, n=21)
    res = run_nbody(cfg)
    assert res.mass_loss and res.event is None
    m = res.stellar["m1"]["mass_msun"]
    assert m[-1] == pytest.approx(1.2 - 1e-3 * 50.0)
    mm = res.stellar["m1"]["mass_msun"] + res.stellar["m2"]["mass_msun"]
    am = res.orbits["A"]["a"] * mm
    # osculating a of the perturbed binary wiggles by ~5e-4; without the invariant a M the
    # change would be the mass change, 2.3 %
    assert np.max(np.abs(am / am[0] - 1.0)) < 2e-3
    save_run(tmp_path / "ml.h5", cfg, res)
    back = run_nbody(_two_two(rates, backend="nbody", t_max=-50.0, n=21))
    assert back.stellar["m1"]["mass_msun"][-1] == pytest.approx(1.2 + 1e-3 * 50.0)


def test_megno_refused_with_mass_loss():
    with pytest.raises(ValueError, match="megno"):
        Config.model_validate(
            {**_two_two(RATES).model_dump(exclude_none=True), "run": {"t_max_yr": 1, "megno": True}}
        )


def test_nbody_evolving_stars():
    cfg = Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 7.5, "m2": 6.0, "m3": 5.0, "m4": 4.0},
            "orbits": {
                "A": {"P": 200, "e": 0.1, "i": 80},
                "B": {"P": 300, "e": 0.1, "i": 60, "Omega": 30},
                "out": {"P": 40000, "e": 0.3, "i": 70},
            },
            "stellar_model": {"kind": "mist", "age_myr": 40.0},
            "physics": {"gr": True, "tides": True, "evolution": True},
            "run": {"t_max_yr": 500.0, "n_output": 11},
        }
    )
    res = run_nbody(cfg)
    r = res.stellar["m1"]["radius_rsun"]
    assert r[-1] > r[0]  # the 7.5 Msun star is leaving the main sequence
    assert np.max(res.energy_error) < 1e-9  # per output step
