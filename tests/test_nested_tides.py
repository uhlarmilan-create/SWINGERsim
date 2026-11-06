import warnings

import numpy as np
import pytest

from swingersim.backends.secular import run_secular
from swingersim.backends.secular_nested import run_secular_nested
from swingersim.config import Config


def _cfg(**over):
    d = {
        "topology": "2+2",
        "masses": {"m1": 1.0, "m2": 0.8, "m3": 0.9, "m4": 0.7},
        "orbits": {
            "A": {"P": 3.0, "e": 0.3},
            "B": {"P": 4.0, "e": 0.2, "i": 40},
            "out": {"P": 800.0, "e": 0.2},
        },
        "stars": {k: {"lag_s": 3000.0} for k in ("m1", "m2", "m3", "m4")},
        "physics": {"gr": True, "tides": True},
        "run": {"backend": "secular", "t_max_yr": 1e5, "n_output": 11},
    }
    d.update(over)
    return Config.model_validate(d)


def test_nested_dissipation_matches_2p2_backend():
    cfg = _cfg()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = run_secular(cfg)
        b = run_secular_nested(cfg)
    assert b.dissipative
    for n in ("A", "B", "out"):
        assert np.allclose(a.orbits[n]["e"], b.orbits[n]["e"], rtol=1e-5, atol=1e-8)
        assert np.allclose(a.orbits[n]["a"], b.orbits[n]["a"], rtol=1e-8)
    assert b.spins["m1"]["P_rot_d"][-1] == pytest.approx(a.spins["m1"]["P_rot_d"][-1], rel=1e-5)


def test_3p1_with_evolving_stars_and_structure_lags():
    cfg = Config.model_validate(
        {
            "topology": "3+1",
            "masses": {"m1": 1.2, "m2": 1.0, "m3": 0.9, "m4": 0.7},
            "orbits": {
                "inner": {"P": 5.0, "e": 0.2},
                "mid": {"P": 600.0, "e": 0.1, "i": 50},
                "out": {"P": 30000.0, "e": 0.1, "i": 10},
            },
            "stellar_model": {"kind": "sse", "age_myr": 100.0},
            "physics": {
                "gr": True,
                "tides": True,
                "rotation": True,
                "evolution": True,
                "tidal_lag": "hurley2002",
            },
            "run": {"backend": "secular", "t_max_yr": 2e4, "n_output": 11},
        }
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = run_secular(cfg)
    assert r.dissipative and r.event is None
    assert np.all(np.isfinite(r.orbits["inner"]["a"]))
    assert "m1" in r.spins
