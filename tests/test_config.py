import pytest
from pydantic import ValidationError

from swingersim.config import Config

BASE = {
    "topology": "2+2",
    "masses": {"m1": 1, "m2": 1, "m3": 1, "m4": 1},
    "orbits": {"A": {"P": 5}, "B": {"P": 3}, "out": {"P": 400}},
    "run": {"t_max_yr": 1},
}


def _with(**kw):
    d = {**BASE, **kw}
    return Config.model_validate(d)


def test_examples_load(example_2p2, example_3p1):
    assert example_2p2.topology == "2+2"
    assert example_3p1.topology == "3+1"


def test_minimal_ok():
    assert _with().physics.gr is True


def test_wrong_orbits_for_topology():
    with pytest.raises(ValidationError):
        _with(topology="3+1")


def test_both_P_and_a():
    with pytest.raises(ValidationError):
        _with(orbits={"A": {"P": 5, "a": 0.1}, "B": {"P": 3}, "out": {"P": 400}})


def test_unbound_eccentricity():
    with pytest.raises(ValidationError):
        _with(orbits={"A": {"P": 5, "e": 1.2}, "B": {"P": 3}, "out": {"P": 400}})


def test_unimplemented_physics():
    with pytest.raises(ValidationError):
        _with(physics={"evolution": True})


def test_rotation_needs_tides():
    with pytest.raises(ValidationError):
        _with(physics={"rotation": True})
    _with(physics={"tides": True, "rotation": True})


def test_star_and_model_specs():
    c = _with(stars={"m1": {"R": 1.2, "apsidal_k2": 0.01, "P_rot": 2.5, "lag_s": 1.0}})
    assert c.stars["m1"].P_rot == 2.5 and c.stars["m1"].rg2 is None
    with pytest.raises(ValidationError):
        _with(stars={"m1": {"P_rot": -1}})
    with pytest.raises(ValidationError):
        _with(stellar_model={"kind": "table"})


def test_typo_rejected():
    with pytest.raises(ValidationError):
        _with(run={"t_max_yr": 1, "integrtor": "ias15"})


def test_save_creates_missing_dir(tmp_path, example_2p2):
    from swingersim.backends.nbody import run_nbody
    from swingersim.io import save_run

    cfg = example_2p2.model_copy(deep=True)
    cfg.run.t_max_yr = 0.1
    cfg.run.n_output = 3
    out = tmp_path / "new" / "dir" / "run.h5"
    save_run(out, cfg, run_nbody(cfg))
    assert out.exists()
