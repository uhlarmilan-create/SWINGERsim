import warnings

import numpy as np
import pytest

from swingersim import massloss
from swingersim.config import Config


def test_sse_wind_rates():
    # Reimers on the giant branch: 4e-13 eta L R / M
    assert massloss.sse_wind_rate(1.0, 100.0, 20.0, 0.02, 3, 0.5) == pytest.approx(-4e-10)
    assert massloss.sse_wind_rate(1.0, 100.0, 20.0, 0.02, 1, 0.5) == 0.0
    # Nieuwenhuijzen & de Jager for L > 4000 Lsun
    r = massloss.sse_wind_rate(10.0, 1e4, 5.0, 0.02, 1)
    assert r == pytest.approx(-9.6e-15 * 5.0**0.81 * 1e4**1.24 * 10.0**0.16)


def test_sse_wind_track_loses_mass_on_giant_branch():
    tr = massloss.SSEWindTrack(1.0, 0.02, 0.5)
    t = tr.track
    assert tr.current_mass(0.5 * t.t_ms) == pytest.approx(1.0)
    m_tip = tr.current_mass(t.t_hei * 0.999999)
    assert 0.7 < m_tip < 0.85


def _wind_cfg(**phys):
    return Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"out": ["A", "m3"], "A": ["m1", "m2"]},
            "masses": {"m1": 25.0, "m2": 18.0, "m3": 10.0},
            "orbits": {"A": {"P": 6.0}, "out": {"P": 3000.0, "e": 0.2}},
            "stellar_model": {"kind": "mist", "age_myr": 1.0, "vvcrit": 0.4},
            "physics": {"gr": False, "evolution": True, **phys},
            "run": {"backend": "secular", "t_max_yr": 8.0e6},
        }
    )


def test_winds_with_mass_transfer():
    from swingersim import masstransfer as mt

    cfg = _wind_cfg(mass_transfer={}, mass_loss={"source": "track"})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        h = mt.histories(cfg)["A"]
        h0 = mt.histories(_wind_cfg(mass_transfer={}))["A"]
    # the winds start before the transfer and widen the binary (a M = const)
    i = int(np.searchsorted(h.t, h.t_onset))
    assert h.m[0, i - 1] < h.m[0, 0] < 25.0
    am = h.a[:i] * h.m[:, :i].sum(axis=0)
    assert np.ptp(am) / am[0] < 1e-8
    assert h.a[i - 1] > h0.a[0]
    # the tertiary loses mass by its own wind
    m = massloss.masses_at(cfg, [25.0, 18.0, 10.0], float(h.t[-1]))
    assert m[2] < 10.0


def test_sse_track_source_allowed():
    cfg = Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 1.0, "m2": 0.9, "m3": 1.0, "m4": 0.8},
            "orbits": {"A": {"P": 5}, "B": {"P": 7}, "out": {"P": 3000}},
            "stellar_model": {"kind": "sse", "age_myr": 100.0},
            "physics": {"evolution": True, "mass_loss": {"source": "track", "reimers_eta": 0.3}},
            "run": {"backend": "secular", "t_max_yr": 1e6},
        }
    )
    assert cfg.physics.mass_loss.reimers_eta == 0.3
    m = massloss.masses_at(cfg, [1.0, 0.9, 1.0, 0.8], 0.0)
    assert np.allclose(m, [1.0, 0.9, 1.0, 0.8])
