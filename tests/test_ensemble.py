import json
import math

import numpy as np
import pytest
from pydantic import ValidationError

from swingersim.config import Config, Distribution
from swingersim.ensemble import (
    draw_members,
    draw_value,
    load_rows,
    report,
    run_ensemble,
    wilson,
)


def _cfg(**ens):
    data = {
        "topology": "2+2",
        "masses": {"m1": 1.0, "m2": 1.0, "m3": 1.0, "m4": 1.0},
        "orbits": {
            "A": {"a": 0.5, "e": 0.1},
            "B": {"a": 0.4, "e": 0.1},
            "out": {"a": 8.0, "e": 0.2},
        },
        "physics": {"gr": False},
        "run": {"t_max_yr": 50.0, "n_output": 20, "backend": "secular"},
        "ensemble": {"n": 3, **ens},
    }
    return Config.model_validate(data)


def test_distributions():
    rng = np.random.default_rng(1)
    assert draw_value(2.5, rng) == 2.5
    iso = np.array([draw_value("isotropic", rng) for _ in range(20000)])
    assert abs(np.mean(np.cos(np.radians(iso)))) < 0.02  # cos i uniform in (-1, 1)
    assert abs(np.std(np.cos(np.radians(iso))) - 1 / math.sqrt(3)) < 0.01
    ang = [draw_value("angle", rng) for _ in range(1000)]
    assert 0 <= min(ang) and max(ang) < 360
    d = Distribution(normal=(1.0, 1.0), min=0.5, max=1.2)
    v = [draw_value(d, rng) for _ in range(500)]
    assert 0.5 <= min(v) and max(v) <= 1.2
    d = Distribution(isotropic=(80.0, 90.0))
    v = [draw_value(d, rng) for _ in range(500)]
    assert 80.0 <= min(v) and max(v) <= 90.0
    d = Distribution(loguniform=(1.0, 100.0))
    v = np.log10([draw_value(d, rng) for _ in range(5000)])
    assert abs(np.mean(v) - 1.0) < 0.05


def test_distribution_needs_one_kind():
    with pytest.raises(ValidationError):
        Distribution(uniform=(0, 1), normal=(0, 1))
    with pytest.raises(ValidationError):
        Distribution()


def test_vary_paths_are_checked():
    with pytest.raises(ValidationError):
        _cfg(vary={"orbits.C.e": 0.1})  # no orbit C in 2+2
    with pytest.raises(ValidationError):
        _cfg(vary={"orbits.A.ecc": 0.1})
    with pytest.raises(ValidationError):
        _cfg(vary={"masses.m5": 1.0})
    with pytest.raises(ValidationError):
        _cfg(vary={"run.t_max_yr": 1.0})
    _cfg(vary={"stars.m1.R": 1.0, "masses.m2": 0.9, "orbits.out.P": 3000.0})


def test_members_deterministic_and_p_replaces_a():
    cfg = _cfg(vary={"orbits.out.P": {"uniform": [2000, 3000]}, "orbits.A.i": "isotropic"})
    m1, m2 = draw_members(cfg), draw_members(cfg)
    assert [m.params for m in m1] == [m.params for m in m2]
    for m in m1:
        assert m.config["orbits"]["out"]["a"] is None
        assert 2000 <= m.config["orbits"]["out"]["P"] <= 3000
    other = draw_members(_cfg(seed=2, vary={"orbits.A.i": "isotropic"}))
    assert other[0].params != m1[0].params


def test_invalid_draws_are_redrawn():
    # e beyond 1 is invalid and must be drawn again
    cfg = _cfg(vary={"orbits.out.e": {"uniform": [0.5, 1.5]}})
    for m in draw_members(cfg):
        assert m.params["orbits.out.e"] < 1.0
    # an outer orbit inside the binaries fails the stability requirement
    cfg = _cfg(vary={"orbits.out.a": {"uniform": [0.6, 20.0]}}, require=["stable_ma01"])
    for m in draw_members(cfg):
        assert m.params["orbits.out.a"] > 2.0


def test_run_resume_and_report(tmp_path):
    cfg = _cfg(vary={"orbits.A.i": "isotropic"}, directions=["forward", "backward"])
    rows = run_ensemble(cfg, tmp_path)
    assert len(rows) == 6
    assert {r["direction"] for r in rows} == {"forward", "backward"}
    assert all(r["fate"] == "intact" for r in rows)
    back = [r for r in rows if r["direction"] == "backward"]
    assert all(r["t_end_yr"] == -50.0 for r in back)
    assert (tmp_path / "summary.csv").exists()
    # resume: an extra member only runs what is missing
    with pytest.raises(FileExistsError):
        run_ensemble(cfg, tmp_path)
    lines = (tmp_path / "rows.jsonl").read_text().splitlines()
    (tmp_path / "rows.jsonl").write_text("\n".join(lines[:4]) + "\n")
    calls = []
    run_ensemble(cfg, tmp_path, resume=True, progress=lambda k, n, r: calls.append(k))
    assert len(calls) == 2
    assert len(load_rows(tmp_path)) == 6
    text = report(load_rows(tmp_path))
    assert "forward: 3 members" in text and "backward: 3 members" in text
    # a changed configuration cannot be resumed into the same directory
    cfg2 = _cfg(vary={"orbits.A.i": "isotropic"}, seed=5)
    with pytest.raises(ValueError):
        run_ensemble(cfg2, tmp_path, resume=True)
    json.loads(lines[0])


def test_failing_member_is_recorded(tmp_path):
    from swingersim.ensemble import run_member

    cfg = _cfg()
    task = {
        "index": 0,
        "direction": "forward",
        "config": {**cfg.model_dump(exclude={"ensemble"}), "topology": "5+5"},
        "params": {},
        "base_dir": None,
        "mlp": False,
        "ghost": False,
        "ghost_orbits": 10.0,
        "save_dir": "",
    }
    row = run_member(task)
    assert row["fate"] == "error" and "ValidationError" in row["error"]


def test_wilson_interval():
    p, lo, hi = wilson(50, 100)
    assert p == 0.5 and lo < 0.5 < hi and abs((hi - lo) / 2 - 0.05) < 0.002
    assert wilson(0, 10)[1] == 0.0


def test_mirror_inclination():
    rng = np.random.default_rng(3)
    d = Distribution(normal=(85.0, 1.0), max=90.0, mirror=True)
    v = np.array([draw_value(d, rng) for _ in range(2000)])
    assert set(np.unique(v > 90.0)) == {False, True}
    assert abs(np.mean(v > 90.0) - 0.5) < 0.04
    assert np.all((v <= 90.0) | (v >= 90.0))


def _zlk_nbody(**run):
    return Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 1.0, "m2": 1.0, "m3": 1.0, "m4": 1.0},
            "orbits": {
                "A": {"a": 0.573, "e": 0.3, "i": 72.5, "Omega": 40, "omega": 10},
                "B": {"a": 0.351, "e": 0.1, "i": 95.7, "Omega": 200, "omega": 70},
                "out": {"a": 6.0, "e": 0.3},
            },
            "physics": {"gr": False},
            "run": {"t_max_yr": 1000.0, "n_output": 200, "collisions": False, **run},
        }
    )


def test_stop_when_e_reached_and_wall_time():
    from swingersim.backends.nbody import run_nbody

    with pytest.raises(ValidationError):
        _zlk_nbody(stop_when_e_reached={"C": 0.9})
    with pytest.raises(ValidationError):
        _zlk_nbody(stop_when_e_reached={"B": 1.0})
    res = run_nbody(_zlk_nbody(stop_when_e_reached={"B": 0.9}, extremes_every_orbits=1.0))
    assert res.event == "e_reached"
    assert res.extremes["B"]["e_max"] >= 0.9
    assert res.times[-1] < 1000.0
    res = run_nbody(_zlk_nbody(max_wall_s=0.05))
    assert res.event == "wall_time" and res.times[-1] < 1000.0


def test_resume_with_more_members(tmp_path):
    cfg = _cfg(vary={"orbits.A.i": "isotropic"})
    run_ensemble(cfg, tmp_path)
    cfg.ensemble.n = 5
    calls = []
    rows = run_ensemble(cfg, tmp_path, resume=True, progress=lambda k, n, r: calls.append(k))
    assert len(calls) == 2 and len(rows) == 5
    # the first three members are the same draws as before
    assert [r["orbits.A.i"] for r in rows[:3]] == [
        m.params["orbits.A.i"] for m in draw_members(_cfg(vary={"orbits.A.i": "isotropic"}))
    ]


def test_conditioning_recovers_parameter():
    from swingersim.ensemble import effective_size, observation_weights

    rng = np.random.default_rng(5)
    x = rng.uniform(0, 1, 4000)
    rows = []
    for k, xi in enumerate(x):
        for d in ("forward", "backward"):
            rows.append(
                {
                    "member": k,
                    "direction": d,
                    "p.x": xi,
                    "fate": "intact",
                    "obs": 2 * xi if d == "forward" else -1.0,
                    "A_e_max": 0.1,
                }
            )
    w = observation_weights(rows, {"obs": {"value": 1.0, "sigma": 0.1}})
    ww = np.array([w[k] for k in range(len(x))])
    mean = np.sum(ww * x) / np.sum(ww)
    sd = np.sqrt(np.sum(ww * (x - mean) ** 2) / np.sum(ww))
    assert abs(mean - 0.5) < 0.01 and abs(sd - 0.05) < 0.01  # forward values are used
    assert 620 < effective_size(ww) < 800  # (sum w)^2 / sum w^2 = 709 for this case
    text = report(rows, weights=w)
    assert "conditioned on the observations" in text and "p.x" in text
    with pytest.raises(ValueError):
        observation_weights(rows, {"nonexistent": {"value": 1.0, "sigma": 0.1}})


def test_tilt_rate_column():
    from swingersim.backends.secular import run_secular
    from swingersim.ensemble import summarize

    cfg = _cfg()
    cfg.orbits["B"] = cfg.orbits["B"].model_copy(update={"i": 60.0})
    cfg.orbits["out"] = cfg.orbits["out"].model_copy(update={"i": 80.0})
    res = run_secular(cfg)
    row = summarize(cfg, res)
    # |i - 90| changes with the opposite sign of i below 90 deg
    assert row["B_dtiltdt0_deg_yr"] == pytest.approx(-row["B_didt0_deg_yr"], rel=1e-6)


def test_conditioning_limits_and_backward_columns():
    from swingersim.config import ObservedValue
    from swingersim.ensemble import observation_weights

    rows = [
        {"member": 0, "direction": "forward", "x": 1.0, "fill": 0.5},
        {"member": 0, "direction": "backward", "x": 1.0, "fill": 1.5},
        {"member": 1, "direction": "forward", "x": 1.0, "fill": 0.5},
        {"member": 1, "direction": "backward", "x": 1.0, "fill": 0.8},
    ]
    w = observation_weights(rows, {"backward:fill": ObservedValue(max=1.0)})
    assert w == {0: 0.0, 1: 1.0}
    w = observation_weights(rows, {"fill": {"min": 0.6}})
    assert w == {0: 0.0, 1: 0.0}
    with pytest.raises(ValidationError):
        ObservedValue(value=1.0)
    with pytest.raises(ValidationError):
        ObservedValue()


def test_resume_after_schema_extension(tmp_path):
    """An ensemble written by an older swingersim (without newer config keys) can be resumed."""
    import yaml

    cfg = _cfg(vary={"orbits.A.i": "isotropic"})
    run_ensemble(cfg, tmp_path)
    meta = yaml.safe_load((tmp_path / "ensemble.yaml").read_text())
    meta["fingerprint"] = "old-version"
    del meta["config"]["hierarchy"]  # a key that older versions did not have
    (tmp_path / "ensemble.yaml").write_text(yaml.safe_dump(meta))
    cfg.ensemble.n = 4
    assert len(run_ensemble(cfg, tmp_path, resume=True)) == 4
    meta = yaml.safe_load((tmp_path / "ensemble.yaml").read_text())
    meta["fingerprint"] = "x"
    meta["config"]["run"]["t_max_yr"] = 99.0  # a real change is still refused
    (tmp_path / "ensemble.yaml").write_text(yaml.safe_dump(meta))
    with pytest.raises(ValueError, match="run.t_max_yr"):
        run_ensemble(cfg, tmp_path, resume=True)


def test_conditioning_on_survival():
    """'fate' constraints: an observed system must have survived into the past (and the near
    future); members that break up in either direction are excluded."""
    from swingersim.config import ObservedValue
    from swingersim.ensemble import observation_weights

    rows = [
        {"member": 0, "direction": "forward", "fate": "intact", "x": 1.0},
        {"member": 0, "direction": "backward", "fate": "ejection", "x": 1.0},
        {"member": 1, "direction": "forward", "fate": "intact", "x": 1.0},
        {"member": 1, "direction": "backward", "fate": "intact", "x": 1.0},
        {"member": 2, "direction": "forward", "fate": "collision", "x": 1.0},
        {"member": 2, "direction": "backward", "fate": "intact", "x": 1.0},
    ]
    both = {"backward:fate": {"equals": "intact"}, "forward:fate": {"equals": "intact"}}
    assert observation_weights(rows, both) == {0: 0.0, 1: 1.0, 2: 0.0}
    past = {"backward:fate": ObservedValue(equals=["intact", "collision"])}
    assert observation_weights(rows, past) == {0: 0.0, 1: 1.0, 2: 1.0}
    with pytest.raises(ValidationError):
        ObservedValue(equals="intact", max=1.0)


def test_report_prints_the_evidence():
    from swingersim.ensemble import evidence, report

    rows = [
        {"member": i, "direction": "forward", "fate": "stable", "wall_s": 1.0} for i in range(4)
    ]
    w = {0: 1.0, 1: 0.5, 2: 0.0, 3: 0.5}
    z, dz, n = evidence(rows, w)
    assert n == 4 and abs(z - 0.5) < 1e-12 and dz > 0
    assert "evidence Z = 0.5" in report(rows, weights=w)
    none = {0: 1e-20, 1: 0.0, 2: 0.0, 3: 0.0}
    assert "quote Z < 1/n" in report(rows, weights=none)
