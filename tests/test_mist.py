"""MIST v1.2 tracks (stellar_model.kind: mist) and the pre-main-sequence / core-helium-burning
segments of the Claret (2023) structure tables."""

import numpy as np
import pytest

from swingersim.config import Config
from swingersim.evolution import EvolvedBeyondGrid
from swingersim.mist import available, mist_track
from swingersim.stellar import claret_structure, star_properties


def _cfg(**sm):
    return Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 1.0, "m2": 0.9, "m3": 1.2, "m4": 2.0},
            "orbits": {"A": {"P": 5.0}, "B": {"P": 6.0}, "out": {"P": 3000.0}},
            "run": {"t_max_yr": 10.0},
            "stellar_model": {"kind": "mist", **sm},
        }
    )


def test_grids_installed():
    grids = available()
    assert (0.0, 0.0) in grids and (0.0, 0.4) in grids
    assert min(g[0] for g in grids) <= -1.0 and max(g[0] for g in grids) >= 0.5


def test_track_phases_and_limits():
    t = mist_track(1.0, 0.0, 0.0)
    assert t.t_start < 0.01
    assert 10.0 < t.t_ms < t.t_hei < t.t_end
    assert t.at(1.0).phase == 0 and t.at(1000.0).phase == 1 and t.at(t.t_end - 1.0).phase == 6
    r = [t.at(a).radius_rsun for a in (1.0, 3.0, 10.0, 30.0)]
    assert np.all(np.diff(r) < 0)  # contraction on the pre-main sequence
    with pytest.raises(EvolvedBeyondGrid):
        t.at(t.t_end + 1.0)
    # metallicity interpolation lies between the neighbouring grids
    lo, hi = mist_track(1.0, -0.25, 0.0).t_ms, mist_track(1.0, 0.0, 0.0).t_ms
    assert min(lo, hi) <= mist_track(1.0, -0.1, 0.0).t_ms <= max(lo, hi)


def test_claret_segments():
    # pre-main sequence of 1 Msun near its start: almost an n = 1.5 polytrope (k2 = 0.143)
    k2, rg2 = claret_structure(1.0, 3.0, 0.0, "pms")
    assert 0.1 < k2 < 0.15 and 0.18 < rg2 < 0.21
    k2c, _ = claret_structure(1.0, 2.5, 0.0, "cheb")
    assert 0.01 < k2c < 0.05
    with pytest.raises(ValueError):
        claret_structure(1.0, 4.0, 0.0, "agb")


def test_mist_star_properties():
    p = star_properties(_cfg(age_myr=1.0, vvcrit=0.0), "m1", 1.0)
    assert p.phase == "pre-main sequence" and p.radius_rsun > 2.0 and p.apsidal_k2 > 0.1
    p = star_properties(_cfg(age_myr=1150.0, vvcrit=0.0), "m4", 2.0)
    assert p.phase == "core helium burning" and p.radius_rsun > 5.0
    with pytest.raises(ValueError):
        _cfg()  # age needed
    # since 0.18 the convective envelopes of MIST stars are known (giants: above the core;
    # main sequence: SSE fractions), so the envelope lags of Hurley et al. (2002) are allowed
    cfg = Config.model_validate(
        {
            **_cfg(age_myr=1000.0).model_dump(exclude={"ensemble"}),
            "physics": {"tides": True, "tidal_lag": "hurley2002"},
        }
    )
    assert star_properties(cfg, "m1", 1.0).turnover_yr > 0.0
