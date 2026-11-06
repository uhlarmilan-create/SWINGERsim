"""Direct tests of modules the code audit found tested only indirectly (diagnostics, secular
kernels, quadrupole theory, report) and of implementations that must stay identical."""

import math
import warnings

import numpy as np
import pytest

from swingersim import units
from swingersim.config import Config

warnings.simplefilter("ignore")


# ---------------------------------------------------------------- same equations, two codes
def test_leconte_numba_kernel_equals_python_rates():
    from swingersim.secular.kernel import leconte_pair
    from swingersim.tides import leconte_rates

    rng = np.random.default_rng(3)
    G = units.G
    for _ in range(200):
        m = rng.uniform(0.3, 5.0, 2)
        r = rng.uniform(0.5, 4.0, 2) * units.RSUN_AU
        k2 = rng.uniform(0.002, 0.05, 2)
        rg2 = rng.choice([0.0, 0.05, 0.1], 2)
        lag = rng.choice([0.0, 1e-7, 3e-6], 2)
        a = rng.uniform(0.03, 0.3)
        e = rng.uniform(0.0, 0.7)
        w = rng.uniform(0.5, 3.0, 2) * math.sqrt(G * m.sum() / a**3)
        py = leconte_rates(tuple(m), tuple(r), tuple(k2), tuple(rg2), tuple(lag), a, e, tuple(w))
        nb = leconte_pair(
            G,
            m[0],
            m[1],
            r[0],
            r[1],
            k2[0],
            k2[1],
            rg2[0],
            rg2[1],
            lag[0],
            lag[1],
            a,
            e,
            w[0],
            w[1],
        )
        np.testing.assert_allclose(nb, [py.da_dt, py.de_dt, *py.dspin_dt], rtol=1e-11, atol=1e-300)


def test_ma01_one_formula():
    from swingersim.diagnostics.zlk import mardling_aarseth_ratio
    from swingersim.stability import ma01_ratio

    assert mardling_aarseth_ratio(1.0, 20.0, 0.3, 2.0, 1.0, 0.5) == pytest.approx(
        ma01_ratio(1.0, 20.0, 0.3, 0.5, 0.5), rel=1e-14
    )
    # Mardling & Aarseth (2001): coplanar, equal masses, circular: R_p,crit = 2.8 * 1.5^0.4 a_in
    assert ma01_ratio(1.0, 2.8 * 1.5**0.4, 0.0, 0.5, 0.0) == pytest.approx(1.0)


# ---------------------------------------------------------------- diagnostics
def test_gr_periastron_advance_matches_textbook():
    from swingersim.diagnostics.orbits import gr_periastron_advance_rate, mutual_inclination

    # Mercury: 42.98 arcsec per century (Sun alone)
    rate = gr_periastron_advance_rate(1.0, 0.387098, 0.205630)
    assert math.degrees(rate) * 3600 * 100 == pytest.approx(42.98, rel=2e-3)
    h1, h2 = np.array([0, 0, 1.0]), np.array([0, math.sin(0.3), math.cos(0.3)])
    assert mutual_inclination(h1, h2) == pytest.approx(0.3)


def test_measured_apsidal_period_and_eclipses():
    from swingersim.diagnostics.binary import eclipsing, measured_apsidal_period

    t = np.linspace(0, 100, 201)
    omega = np.mod(2 * np.pi * t / 400.0, 2 * np.pi)
    assert measured_apsidal_period(t, 0.0 * t, omega) == pytest.approx(400.0, rel=1e-9)
    assert math.isnan(measured_apsidal_period(t, 0.0 * t, omega, e=np.full_like(t, 0.001)))
    orbit = {"inc": np.radians([90.0, 85.0, 60.0]), "a": np.array([0.05, 0.05, 0.05])}
    assert eclipsing(orbit, 0.01).tolist() == [True, True, False]  # cos 85 deg = 0.087 < 0.2


def test_extreme_tracker_on_a_kepler_orbit():
    from swingersim.diagnostics.extremes import ExtremeTracker
    from swingersim.hierarchy import build_tree, initialize
    from swingersim.kepler import Elements, elements_to_state

    cfg = Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"A": ["m1", "m2"]},
            "masses": {"m1": 1.0, "m2": 0.5},
            "orbits": {"A": {"a": 1.0, "e": 0.4}},
            "run": {"t_max_yr": 1.0},
        }
    )
    root, masses = build_tree(cfg)
    tr = ExtremeTracker(root, masses, units.G)
    mu = units.G * 1.5
    for M in np.linspace(0, 2 * np.pi, 400):
        r, v = elements_to_state(mu, Elements(1.0, 0.4, 0.0, 0.0, 0.0, M))
        pos = np.array([-0.5 / 1.5 * r, 1.0 / 1.5 * r])
        vel = np.array([-0.5 / 1.5 * v, 1.0 / 1.5 * v])
        tr.update(float(M), pos, vel)
    res = tr.result()["A"]
    assert res["e_max"] == pytest.approx(0.4, rel=1e-9)
    assert res["q_min_au"] == pytest.approx(0.6, rel=1e-9)
    assert initialize(root, masses)[0].shape == (2, 3)


# ---------------------------------------------------------------- quadrupole ZLK theory
def test_quadrupole_test_particle_limit_emax():
    """Test-particle quadrupole limit: e_max = sqrt(1 - 5/3 cos^2 i0) for e0 -> 0, omega free."""
    from swingersim.secular.quadrupole import predict

    p = predict(1.0, 1e-6, 1.0, 1.0, 100.0, 1e-3, 0.0, 65.0, 0.0)
    assert p.e_max == pytest.approx(
        math.sqrt(1 - 5 / 3 * math.cos(math.radians(65.0)) ** 2), abs=0.01
    )


# ---------------------------------------------------------------- report
def test_run_report_html_has_the_sections(tmp_path):
    from swingersim.cli import main

    cfg = tmp_path / "r.yaml"
    cfg.write_text(
        """
topology: '2+2'
masses: {m1: 1.0, m2: 0.8, m3: 0.9, m4: 0.7}
orbits: {A: {P: 5.0, i: 88}, B: {P: 7.0, i: 85}, out: {P: 900.0, e: 0.2, i: 70}}
run: {t_max_yr: 5.0, n_output: 51}
"""
    )
    assert main(["run", str(cfg), "-o", str(tmp_path / "r.h5"), "--no-viewer", "--no-plots"]) in (
        0,
        None,
    )
    html = (tmp_path / "r.html").read_text(encoding="utf-8")
    for part in ("swingersim", "Orbits", "r.yaml"):
        assert part in html
