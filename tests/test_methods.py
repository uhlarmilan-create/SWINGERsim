"""Chaos maps, statistical escape theory, eclipse timings and simulation-based inference."""

import importlib.util
import json

import numpy as np
import pytest

from swingersim import units
from swingersim.chaosmap import chaos_map, megno_array, stable_from
from swingersim.config import Config
from swingersim.escape_theory import (
    angular_momentum_scale,
    binary_single_integrals,
    compare_with_scatter,
    escape_outcomes,
)
from swingersim.sbi_fit import training_set
from swingersim.timing import eclipse_times, linear_ephemeris, o_minus_c

CB = {
    "topology": "custom",
    "hierarchy": {"out": ["inner", "m3"], "inner": ["m1", "m2"]},
    "masses": {"m1": 1.0, "m2": 1.0, "m3": 1.0e-9},
    "orbits": {"inner": {"a": 1.0}, "out": {"a": 3.0}},
    "stars": {"m1": {"R": 0.001}, "m2": {"R": 0.001}, "m3": {"R": 0.001}},
    "physics": {"gr": False},
    "run": {"t_max_yr": 10, "collisions": False},
}


def test_chaos_map_finds_unstable_and_regular_orbits(tmp_path):
    cfg = Config.model_validate(CB)
    xs = np.array([1.6, 3.0])
    rows = chaos_map(cfg, tmp_path, "orbits.out.a", xs, t_max_yr=300.0)
    a = megno_array(rows)
    assert np.isinf(a[0, 0])  # the particle at 1.6 a_bin escapes
    assert a[1, 0] == pytest.approx(2.0, abs=0.2)  # regular at 3 a_bin
    assert stable_from(xs, a) == 3.0
    meta = json.loads((tmp_path / "chaosmap.json").read_text())
    assert meta["x"] == "orbits.out.a"
    with pytest.raises(FileExistsError):
        chaos_map(cfg, tmp_path, "orbits.out.a", xs, t_max_yr=300.0)
    rows = chaos_map(cfg, tmp_path, "orbits.out.a", xs, t_max_yr=300.0, resume=True)
    assert len(rows) == 2


def test_megno_point_is_reproducible_and_lyapunov_time_is_1_over_lambda():
    """MEGNO uses a fixed seed (the same value on every run), and the Lyapunov time is
    1 / lambda = 1 / (2 sim.lyapunov()): REBOUND's lyapunov() is the slope of <Y> ~ lambda t / 2."""
    import rebound

    from swingersim.backends.nbody import MEGNO_SEED, lyapunov_time_yr
    from swingersim.chaosmap import megno_point

    d = json.loads(json.dumps(CB))
    d["orbits"]["out"]["a"] = 3.0
    r1 = megno_point(d, 200.0, 60.0, 50.0)
    r2 = megno_point(d, 200.0, 60.0, 50.0)
    assert r1["megno"] == r2["megno"] and r1["lyapunov_yr"] == r2["lyapunov_yr"]
    sim = rebound.Simulation()
    sim.add(m=1.0)
    sim.add(m=1e-3, a=1.0, e=0.1)
    sim.add(m=1e-3, a=1.2, e=0.1, f=2.0)
    sim.move_to_com()
    sim.init_megno(seed=MEGNO_SEED)
    sim.integrate(300.0 * 2 * np.pi)
    slope = sim.lyapunov()
    assert slope > 0
    assert lyapunov_time_yr(sim) == pytest.approx(0.5 / slope)


def test_escape_theory_symmetry_and_mass_dependence():
    E0 = -units.G / 2.0
    eq = escape_outcomes([1.0, 1.0, 1.0], E0, 0.3 * angular_momentum_scale([1, 1, 1], E0))
    assert all(p == pytest.approx(1 / 3, abs=1e-6) for p in eq.probabilities.values())
    r = escape_outcomes([1.0, 1.0, 0.3], E0, 0.0)
    assert r.probabilities[2] > 0.8  # the lightest star escapes most often
    assert sum(r.probabilities.values()) == pytest.approx(1.0)
    with pytest.raises(ValueError):
        escape_outcomes([1, 1, 1], 1.0, 0.0)


def test_binary_single_integrals():
    E0, L0 = binary_single_integrals((1.0, 1.0), 1.0, 1.0, 0.0, 0.0, 0.0)
    assert E0 == pytest.approx(-units.G / 2.0)
    assert L0 == pytest.approx(0.5 * np.sqrt(2 * units.G))


def test_compare_with_scatter_picks_resonant_outcomes():
    row = {
        "status": "resolved",
        "scrambles": 3,
        "E0": -units.G / 2,
        "L0": 0.0,
        "hierarchy": "[1,3] | 2",
        "orbits": {"[1,3]": {"a_au": 0.8, "e": 0.5}},
    }
    c = compare_with_scatter([row, {**row, "scrambles": 0}], [1.0, 1.0, 1.0])
    assert c["n"] == 1 and c["observed"][1] == 1.0
    assert c["predicted"][1] == pytest.approx(1 / 3, abs=1e-3)


def _two_binaries(i_a, i_b):
    orb = {"e": 0.0, "Omega": 0.0}
    return Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 1.0, "m2": 1.0, "m3": 1.0, "m4": 1.0},
            "orbits": {
                "A": {"P": 5.0, "i": i_a, **orb},
                "B": {"P": 7.0, "i": i_b, **orb},
                "out": {"a": 20.0, "i": 90.0, **orb},
            },
            "stars": {f"m{k}": {"R": 1.0} for k in range(1, 5)},
            "physics": {"gr": False},
            "run": {"t_max_yr": 1.0},
        }
    )


def test_eclipse_timings():
    ecl, p_nom = eclipse_times(_two_binaries(90.0, 30.0), 0.5)
    a = [e for e in ecl if e.binary == "A"]
    assert not [e for e in ecl if e.binary == "B"]  # B at 30 deg never eclipses
    n_expected = 0.5 * units.DAYS_PER_YEAR / 5.0
    assert abs(len([e for e in a if e.kind == "primary"]) - n_expected) <= 1
    assert abs(len([e for e in a if e.kind == "secondary"]) - n_expected) <= 1
    assert all(e.impact < 0.01 for e in a)  # edge-on: central eclipses
    rows = o_minus_c(ecl, p_nom, observed=False)
    assert max(abs(r["o_c_d"]) for r in rows) < 1e-3  # strictly periodic over half a year
    assert rows[0]["P_d"] == pytest.approx(5.0, rel=1e-3)


def test_linear_ephemeris():
    t = 2.0 + 0.01 * np.arange(50) + 1e-6 * np.sin(np.arange(50))
    _t0, P, E = linear_ephemeris(t, 0.0101)
    assert P == pytest.approx(0.01, rel=1e-5)
    assert E[-1] == 49


def test_sbi_training_set():
    rows = []
    for m in range(30):
        rows.append(
            {
                "member": m,
                "direction": "forward",
                "orbits.out.P": 100.0 + m,
                "rate": 0.1 * m,
                "fill": 0.5,
            }
        )
        rows.append({"member": m, "direction": "backward", "fill": 1.5 if m < 10 else 0.5})
    observed = {"rate": {"value": 1.0, "sigma": 0.2}, "backward:fill": {"max": 1.0}}
    th, _x, cols, sig, val = training_set(rows, ["orbits.out.P"], observed, ["orbits.out.P"])
    assert len(th) == 20  # 10 members fail the backward limit
    assert th[0, 0] == pytest.approx(np.log10(110.0))
    assert cols == ["rate"] and sig[0] == 0.2 and val[0] == 1.0


@pytest.mark.slow
@pytest.mark.skipif(importlib.util.find_spec("sbi") is None, reason="sbi not installed")
def test_sbi_recovers_linear_gaussian():
    from swingersim.sbi_fit import train

    rng = np.random.default_rng(0)
    th = rng.uniform(-3, 3, size=(1500, 1))
    model = train(th, 2.0 * th, np.array([0.2]), ["a"], ["x"], max_epochs=200)
    s = model.sample(np.array([1.0]), 5000)
    assert s.mean() == pytest.approx(0.5, abs=0.05)
    assert s.std() == pytest.approx(0.1, rel=0.2)
    assert model.history["training_loss"] and model.history["validation_loss"]


@pytest.mark.slow
@pytest.mark.skipif(importlib.util.find_spec("sbi") is None, reason="sbi not installed")
def test_sbi_calibration(tmp_path):
    from swingersim.sbi_fit import calibration, coverage_curve, plot_diagnostics, train

    rng = np.random.default_rng(1)
    th = rng.uniform(-3, 3, size=(1500, 1))
    x = 2.0 * th
    model = train(th[:1200], x[:1200], np.array([0.2]), ["a"], ["x"], max_epochs=200)
    u = calibration(model, th[1200:], x[1200:], n_test=150, n_samples=300)
    _, cov = coverage_curve(u, [0.5, 0.9])
    assert cov[:, 0] == pytest.approx([0.5, 0.9], abs=0.12)
    assert plot_diagnostics(model, u, tmp_path / "d.png").exists()
