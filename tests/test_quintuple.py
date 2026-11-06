"""Hierarchies beyond quadruples: a fifth star (presets 2+2+1, 3+1+1, 2+(2+1)) and custom trees."""

import numpy as np
import pytest
from pydantic import ValidationError

from swingersim.backends.nbody import run_nbody
from swingersim.config import Config, check_hierarchy
from swingersim.hierarchy import build_tree, initialize, skeleton


def _cfg(topology="2+2+1", **extra):
    orbits = {
        "2+2+1": {
            "A": {"P": 5.0, "e": 0.1, "i": 80},
            "B": {"P": 3.0, "e": 0.05, "i": 85, "Omega": 30},
            "out": {"P": 400, "e": 0.2, "i": 70, "Omega": 10},
            "out5": {"P": 20000, "e": 0.3, "i": 20, "Omega": 100},
        },
        "2+(2+1)": {
            "A": {"P": 6.7, "e": 0.3, "i": 87},
            "B": {"P": 0.92, "e": 0.0, "i": 63},
            "Bt": {"P": 119.0, "e": 0.2, "i": 40, "Omega": 50},
            "out": {"P": 5000.0, "e": 0.3, "i": 80},
        },
    }[topology]
    data = {
        "topology": topology,
        "masses": {f"m{k}": 1.0 for k in range(1, 6)},
        "orbits": orbits,
        "physics": {"gr": False},
        "run": {"t_max_yr": 2.0, "n_output": 20},
        **extra,
    }
    return Config.model_validate(data)


def test_presets_and_tree():
    cfg = _cfg()
    assert cfg.star_names == ["m1", "m2", "m3", "m4", "m5"]
    root, masses = build_tree(cfg)
    assert root.name == "out5" and root.right == 4 and len(masses) == 5
    pos, vel = initialize(root, masses)
    assert np.allclose(masses @ pos, 0.0, atol=1e-12) and np.allclose(masses @ vel, 0.0)
    assert skeleton("2+(2+1)").right.name == "Bt"


def test_custom_hierarchy_and_errors():
    spec = {"top": ["A", "C"], "A": ["m1", "m2"], "C": ["m3", "T"], "T": ["m4", "m5"]}
    assert check_hierarchy(spec) == ("top", ["m1", "m2", "m3", "m4", "m5"])
    data = {
        "topology": "custom",
        "hierarchy": spec,
        "masses": {f"m{k}": 1.0 for k in range(1, 6)},
        "orbits": {
            "top": {"a": 50.0, "e": 0.1},
            "A": {"a": 0.5},
            "C": {"a": 5.0, "e": 0.1, "i": 30},
            "T": {"a": 0.3},
        },
        "run": {"t_max_yr": 1.0},
    }
    Config.model_validate(data)
    for bad in (
        {"top": ["A", "m5"], "A": ["m1", "m2"]},  # m3, m4 missing
        {"top": ["A", "B"], "A": ["m1", "m2"], "B": ["m2", "m3"]},  # m2 twice
        {"top": ["A", "x"], "A": ["m1", "m2"]},  # unknown child
        {"A": ["m1", "m2"], "B": ["m3", "m4"]},  # two roots
    ):
        with pytest.raises(ValueError):
            check_hierarchy(bad)
    with pytest.raises(ValidationError):
        Config.model_validate({**data, "topology": "2+2"})  # hierarchy only with custom
    with pytest.raises(ValidationError):
        Config.model_validate({**data, "masses": {f"m{k}": 1.0 for k in range(1, 5)}})


@pytest.mark.parametrize("topology", ["2+2+1", "2+(2+1)"])
def test_five_body_run(topology):
    res = run_nbody(_cfg(topology))
    assert res.event is None
    assert np.nanmax(res.energy_error) < 1e-9
    assert set(res.orbits) == set(_cfg(topology).tree)
    assert len(res.stars) == 5


def test_secular_backend_handles_five_bodies():
    from swingersim.backends.secular import run_secular

    cfg = _cfg()
    cfg.run.t_max_yr = 200.0
    sec = run_secular(cfg)
    nb = run_nbody(cfg)
    assert set(sec.orbits) == set(cfg.tree) and np.nanmax(sec.energy_error) < 1e-8
    for pair in ("out-A", "A-B"):
        assert abs(np.degrees(np.ptp(sec.mutual[pair]) - np.ptp(nb.mutual[pair]))) < 2.0
