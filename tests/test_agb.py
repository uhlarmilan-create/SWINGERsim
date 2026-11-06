"""Giants to the end of the AGB (MIST EEPs to 1409) and their convective envelopes."""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from swingersim.config import Config
from swingersim.mist import EEP_POSTAGB, EEP_TACHEB, mist_track
from swingersim.stellar import mist_envelope, mist_properties


@pytest.mark.parametrize("mass", [1.0, 2.0, 4.0])
def test_tracks_reach_the_end_of_the_agb(mass):
    tr = mist_track(mass, 0.0, 0.0)
    assert tr.n_eep == EEP_POSTAGB
    t_cheb_end = tr._age[EEP_TACHEB - 1]
    assert tr.t_end > t_cheb_end
    end = tr.t_end * (1 - 1e-9)
    st = tr.at(end)
    assert st.phase == 6
    assert st.radius_rsun > 100.0  # an AGB giant
    m_fin = tr.current_mass(end)
    assert 0.5 < m_fin < 1.1  # the core of a future white dwarf (Cummings et al. 2018)
    ages = tr._age
    assert np.all(np.diff(ages) >= 0)


def test_massive_tracks_stop_after_carbon_burning():
    tr = mist_track(20.0, 0.0, 0.0)
    assert tr.n_eep <= 808


def test_giant_envelope_is_the_mass_above_the_core():
    tr = mist_track(1.5, 0.0, 0.0)
    age = tr._age[580]  # on the red giant branch
    st = tr.at(age)
    m_env, r_env = mist_envelope(tr, age, 0.0)
    assert st.phase == 3
    assert m_env == pytest.approx(tr.current_mass(age) - st.m_core, rel=1e-9)
    assert r_env == pytest.approx(st.radius_rsun)
    p = mist_properties(1.5, age, 0.0, 0.0)
    assert p.turnover_yr > 0.1  # giants: turnover of months to a year


def test_main_sequence_envelope_follows_sse():
    tr = mist_track(1.0, 0.0, 0.0)
    age = 0.5 * (tr._age[201] + tr._age[453])
    m_env, r_env = mist_envelope(tr, age, 0.0)
    assert 0.005 < m_env < 0.06  # the Sun: ~0.02 Msun
    assert 0.2 < r_env / tr.at(age).radius_rsun < 0.35  # the Sun: 0.29 R
    assert mist_envelope(mist_track(5.0, 0.0, 0.0), 30.0, 0.0)[0] < 1e-3  # radiative


def test_hurley2002_lags_with_mist_stars():
    cfg = Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"out": ["A", "m3"], "A": ["m1", "m2"]},
            "masses": {"m1": 1.0, "m2": 0.9, "m3": 2.0},
            "orbits": {"A": {"P": 5.0}, "out": {"P": 3000.0, "e": 0.2, "i": 30.0}},
            "stellar_model": {"kind": "mist", "age_myr": 1000.0, "vvcrit": 0.0},
            "physics": {"tides": True, "tidal_lag": "hurley2002"},
            "run": {"t_max_yr": 1.0},
        }
    )
    from swingersim.backends.nbody import build_simulation

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        build_simulation(cfg)


def test_roche_lobe_of_a_tertiary():
    """A star orbiting a binary has a Roche lobe against the total mass of the binary
    (Eggleton 1983 with q = m3 / (m1 + m2)); a giant tertiary that fills it stops the run."""
    from swingersim.backends.secular_nested import run_secular_nested
    from swingersim.stellar import roche_lobe_eggleton

    d = {
        "topology": "custom",
        "hierarchy": {"out": ["A", "m3"], "A": ["m1", "m2"]},
        "masses": {"m1": 1.0, "m2": 1.0, "m3": 2.0},
        "orbits": {"A": {"P": 5.0}, "out": {"a": 4.0, "e": 0.0, "i": 10.0}},
        "stars": {"m3": {"R": 100.0}},
        "run": {
            "t_max_yr": 1000.0,
            "n_output": 11,
            "backend": "secular",
            "stop_on_roche_overflow": True,
        },
    }
    res = run_secular_nested(Config.model_validate(d))
    from swingersim import units

    expected = 100.0 * units.RSUN_AU / (roche_lobe_eggleton(1.0) * 4.0)
    assert res.roche["m3"]["fill_max"] == pytest.approx(expected, rel=0.02)
    d["stars"]["m3"]["R"] = 400.0
    res = run_secular_nested(Config.model_validate(d))
    assert res.event == "roche_overflow"
