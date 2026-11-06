import math
import warnings

import numpy as np
import pytest

from swingersim import pulsation as pu
from swingersim.config import Config


def _cfg(puls, kind="zams", age=None, m1=1.7):
    sm = {"kind": kind} if age is None else {"kind": kind, "age_myr": age}
    return Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"A": ["m1", "m2"]},
            "masses": {"m1": m1, "m2": 0.5},
            "orbits": {"A": {"P": 50.0, "i": 60.0}},
            "stars": {"m1": {"pulsation": puls}},
            "stellar_model": sm,
            "physics": {"gr": False},
            "run": {"t_max_yr": 1.0},
        }
    )


def test_trabucchi_table2_breakpoints():
    # at M = Mbk every parameter equals k_b
    assert 10 ** pu._trabucchi_param("Rb", 2.6) == pytest.approx(421.0)
    assert 10 ** pu._trabucchi_param("Ps", 1.0) == pytest.approx(388.0)
    assert pu._trabucchi_param("alpha", 2.6) == pytest.approx(math.log10(49.7))


def test_lpv_period_continuous_and_saturates():
    m = 1.5
    rb = 10 ** pu._trabucchi_param("Rb", m)
    rs = 10 ** pu._trabucchi_param("Rs", m)
    for r in (rb, rs):
        assert pu.lpv_fundamental_period(m, r * (1 - 1e-9)) == pytest.approx(
            pu.lpv_fundamental_period(m, r * (1 + 1e-9)), rel=1e-6
        )
    assert pu.lpv_fundamental_period(m, 2 * rs) == pytest.approx(10 ** pu._trabucchi_param("Ps", m))
    # periods grow with radius below the saturation and stay below the death line (eq. 3)
    rr = np.linspace(100, rs, 20)
    pp = [pu.lpv_fundamental_period(m, r) for r in rr]
    assert np.all(np.diff(pp) > 0)


def test_mira_like_period():
    # a 1 Msun AGB star of 300 Rsun pulsates with a Mira-like period of about a year
    assert 250 < pu.lpv_fundamental_period(1.0, 300.0) < 420


def test_q_values_inside_poro_ranges():
    for mode, q in pu.Q_RADIAL.items():
        lo, hi = pu.POROS_2024_Q_RANGES[mode]
        assert lo <= q < hi


def test_q_period_sun():
    assert pu.q_period(0.033, 1.0, 1.0) == pytest.approx(0.033)


def test_rr_lyrae_typical_period():
    p = pu.rr_lyrae_fundamental_period(50.0, 6500.0, 0.65)
    assert 0.4 < p < 0.8


def test_flux_modulation_and_emission_time():
    from swingersim.photometry import light_curve

    cfg = _cfg({"type": "manual", "period_d": 0.2, "amplitude_mag": 0.05})
    t = np.linspace(0, 1.0, 400)
    lc = light_curve(cfg, t, band="V", proximity=False)
    y = -2.5 * np.log10(lc.components["m1"] / np.median(lc.components["m1"]))
    assert y.max() - y.min() == pytest.approx(0.1, rel=0.02)
    # the companion is not modulated
    assert np.ptp(lc.components["m2"]) < 1e-6
    assert np.allclose(lc.flux, lc.components["m1"] + lc.components["m2"])


def test_extra_mode_frequency_ratio():
    cfg = _cfg(
        {
            "type": "manual",
            "period_d": 0.1,
            "extra_modes": [{"frequency_ratio": 1.3, "amplitude_mag": 0.01}],
        }
    )
    from swingersim.stellar import star_properties

    p = star_properties(cfg, "m1", 1.7)
    pz = pu.pulsators(cfg, [1.7, 0.5], [p, p])[0]
    assert pz.modes[1].period_d == pytest.approx(0.1 / 1.3)


def test_relation_needs_period():
    with pytest.raises(ValueError, match="period_d"):
        from swingersim.stellar import star_properties

        cfg = _cfg({"type": "gamma_dor"})
        p = star_properties(cfg, "m1", 1.7)
        pu.pulsators(cfg, [1.7, 0.5], [p, p])


def test_extra_mode_needs_one_period():
    with pytest.raises(ValueError):
        _cfg({"type": "manual", "period_d": 0.1, "extra_modes": [{"amplitude_mag": 0.01}]})


def test_period_table_delta_scuti_mist():
    cfg = _cfg({"type": "delta_scuti"}, kind="mist", age=600.0, m1=1.8)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rows = pu.period_table(cfg)
    assert len(rows) == 1
    r = rows[0]
    assert 0.03 < r["period_d"] < 0.2
    assert r["dP_dt_s_per_yr"] is not None and r["dP_dt_s_per_yr"] > 0  # expanding on the MS


def test_lpv_warns_below_dominant_radius():
    from swingersim.stellar import StarProperties

    props = StarProperties(radius_rsun=50.0, apsidal_k2=0.01, rg2=0.1, teff_k=3500.0)

    class S:
        type, mode, period_d, Q = "lpv", "F", None, None

    _, notes = pu.main_period(S, 1.0, props)
    assert any("overtone" in n for n in notes)


def test_cepheid_reproduces_bono_models():
    # Bono et al. (2000), Table 2 (fundamental, Z = 0.02) and Table 1 (first overtone)
    assert pu.cepheid_period(10**3.07, 5800.0, 5.0, 0.02, "F") == pytest.approx(3.5140, rel=0.04)
    assert pu.cepheid_period(10**3.07, 6000.0, 5.0, 0.02, "1O") == pytest.approx(2.1847, rel=0.04)


def test_cepheid_strip_and_z_interpolation():
    red, blue = pu.cepheid_strip(10**3.5, 0.02)
    assert 4500 < red < blue < 6500
    # coefficients at a tabulated Z are the table values
    assert pu._bono_interp(pu.BONO_2000_F, 0.008) == pytest.approx(pu.BONO_2000_F[0.008])


def test_breger_q_and_beta_cep():
    from swingersim.stellar import StarProperties

    props = StarProperties(radius_rsun=5.0, apsidal_k2=0.01, rg2=0.1, teff_k=25000.0)

    class S:
        type, mode, period_d, Q = "beta_cep", "F", None, None

    p, _ = pu.main_period(S, 12.0, props)
    assert p == pytest.approx(0.033 / math.sqrt(12.0 / 125.0))
    assert pu.Q_RADIAL["F"] == 0.0329


def test_gmode_series():
    from swingersim.stellar import star_properties

    cfg = _cfg(
        {"type": "spb", "period_d": 1.0, "period_spacing_s": 6400.0, "n_modes": 4},
        m1=1.7,
    )
    p = star_properties(cfg, "m1", 1.7)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pz = pu.pulsators(cfg, [1.7, 0.5], [p, p])[0]
    per = [md.period_d for md in pz.modes]
    assert len(per) == 4
    assert np.allclose(np.diff(per), 6400.0 / 86400.0)
