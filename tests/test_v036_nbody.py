"""Regression tests of the v0.36 fixes of the N-body backend, run outputs and diagnostics."""

from __future__ import annotations

import csv
import io as _io
import json
import math
from itertools import combinations

import numpy as np
import pytest

from swingersim import units
from swingersim.config import Config


def _tight_2p2(**run):
    d = {
        "topology": "2+2",
        "masses": {"m1": 1.2, "m2": 1.0, "m3": 0.9, "m4": 0.7},
        "orbits": {
            "A": {"P": 3.0, "e": 0.2},
            "B": {"P": 5.0, "e": 0.1, "i": 50},
            "out": {"P": 3000.0, "e": 0.1, "i": 10},
        },
        "physics": {"gr": False},
        "run": {"backend": "nbody", "t_max_yr": 2.0, "n_output": 5, **run},
    }
    return Config.model_validate(d)


# ------------------------------------------------------------------ kepler helpers


def test_vectorised_elements_equal_scalar():
    from swingersim.kepler import (
        Elements,
        cross3,
        eccentricity_vector,
        elements_to_state,
        solve_kepler,
        solve_kepler_array,
        state_to_elements,
        states_to_elements,
    )

    rng = np.random.default_rng(3)
    rs, vs, mus = [], [], []
    for _ in range(200):
        mu = rng.uniform(1, 100)
        el = Elements(
            rng.uniform(0.01, 50),
            rng.choice([0.0, rng.uniform(0, 0.99)]),
            rng.choice([0.0, rng.uniform(0, np.pi)]),
            rng.uniform(0, 2 * np.pi),
            rng.uniform(0, 2 * np.pi),
            rng.uniform(0, 2 * np.pi),
        )
        r, v = elements_to_state(mu, el)
        rs.append(r), vs.append(v), mus.append(mu)
    rs.append(np.array([1.0, 0, 0])), vs.append(np.array([0, 100.0, 0])), mus.append(1.0)
    rs, vs, mus = np.array(rs), np.array(vs), np.array(mus)
    many = states_to_elements(mus, rs, vs)
    for k in range(len(mus) - 1):
        one = state_to_elements(mus[k], rs[k], vs[k])
        for key in ("a", "e", "inc"):
            assert many[key][k] == pytest.approx(getattr(one, key), rel=1e-12, abs=1e-12)
        for key in ("Omega", "omega", "M"):
            d = abs(many[key][k] - getattr(one, key)) % (2 * np.pi)
            assert min(d, 2 * np.pi - d) < 1e-10
    assert all(np.isnan(many[key][-1]) for key in many)  # unbound
    assert np.array_equal(cross3(rs, vs), np.cross(rs, vs))
    ev = eccentricity_vector(mus, rs, vs)
    assert np.allclose(np.linalg.norm(ev[:-1], axis=1), many["e"][:-1], atol=1e-12)
    M = rng.uniform(-10, 10, 500)
    e = rng.uniform(0, 0.999, 500)
    E = solve_kepler_array(M, e)
    assert np.allclose(E, [solve_kepler(m, x) for m, x in zip(M, e, strict=True)], atol=1e-13)


def test_units_constants():
    assert units.MSUN_G == pytest.approx(1.3271244e20 / 6.67430e-11 * 1e3, rel=1e-5)
    assert units.KMS_TO_AUYR == pytest.approx(0.2109, rel=1e-3)


def test_gr_rate_is_the_tides_one():
    from swingersim.diagnostics import gr_periastron_advance_rate
    from swingersim.tides import apsidal_rate_gr

    assert gr_periastron_advance_rate(2.5, 0.05, 0.3) == apsidal_rate_gr(2.5, 0.05, 0.3)


# ------------------------------------------------------------------ hierarchy detection


def _detect_reference(pos, vel, masses, G):
    """The pre-0.36 implementation (numpy per pair)."""
    objs = [
        (i, float(m), np.asarray(p, float), np.asarray(v, float))
        for i, (m, p, v) in enumerate(zip(masses, pos, vel, strict=True))
    ]
    from swingersim.diagnostics.hierarchy import _lowest, canonical

    while len(objs) > 1:
        best, best_a = None, np.inf
        for i, j in combinations(range(len(objs)), 2):
            a_, b_ = objs[i], objs[j]
            dr, dv = b_[2] - a_[2], b_[3] - a_[3]
            mu = G * (a_[1] + b_[1])
            energy = 0.5 * dv @ dv - mu / np.linalg.norm(dr)
            if energy < 0.0 and -mu / (2.0 * energy) < best_a:
                best, best_a = (i, j), -mu / (2.0 * energy)
        if best is None:
            break
        i, j = best
        a_, b_ = objs[i], objs[j]
        m = a_[1] + b_[1]
        merged = (
            (a_[0], b_[0]),
            m,
            (a_[1] * a_[2] + b_[1] * b_[2]) / m,
            (a_[1] * a_[3] + b_[1] * b_[3]) / m,
        )
        objs = [o for k, o in enumerate(objs) if k not in (i, j)] + [merged]
    return tuple(sorted((canonical(o[0]) for o in objs), key=_lowest))


def test_detect_matches_reference():
    from swingersim.diagnostics.hierarchy import detect

    rng = np.random.default_rng(7)
    for _ in range(2000):
        n = int(rng.integers(2, 7))
        pos = rng.normal(size=(n, 3)) * 10 ** rng.uniform(-2, 2, size=(n, 1))
        vel = rng.normal(size=(n, 3)) * 3
        m = rng.uniform(0.1, 3, n)
        assert detect(pos, vel, m, units.G) == _detect_reference(pos, vel, m, units.G)


# ------------------------------------------------------------------ run_nbody


def test_event_row_has_the_event_time(example_2p2):
    from swingersim.backends.nbody import run_nbody

    cfg = example_2p2.model_copy(deep=True)
    cfg.orbits["B"].P = None
    cfg.orbits["B"].a = 0.005  # stars overlap almost at once
    cfg.run.t_max_yr = 50.0
    cfg.run.n_output = 5
    res = run_nbody(cfg)
    assert res.event == "collision"
    assert res.times[-1] == pytest.approx(res.event_time, abs=0) and res.times[-1] < 12.5


def test_t_max_zero_is_refused():
    from swingersim.backends.nbody import run_nbody

    cfg = _tight_2p2()
    cfg.run.t_max_yr = 0.0
    with pytest.raises(ValueError, match="t_max_yr"):
        run_nbody(cfg)


def test_build_setup_matches_build_simulation():
    from swingersim.backends.nbody import build_setup, build_simulation

    cfg = _tight_2p2()
    a = build_simulation(cfg)
    b = build_setup(cfg)
    assert len(a) == 5 and len(b) == 7
    assert np.array_equal(a[4], b[4]) and np.array_equal(b[5], [1.2, 1.0, 0.9, 0.7])
    assert len(b[6]) == 4


def test_peters_substeps_follow_the_circular_decay():
    from swingersim.backends.nbody import peters_evolve
    from swingersim.gw import merger_time_circular_yr

    tau = merger_time_circular_yr(1.4, 1.4, 0.003)
    a, e = peters_evolve(1.4, 1.4, 0.003, 0.0, 0.9 * tau)
    assert a == pytest.approx(0.003 * 0.1**0.25, rel=1e-6) and e == 0.0
    a_back, _ = peters_evolve(1.4, 1.4, a, 0.0, -0.9 * tau)
    assert a_back == pytest.approx(0.003, rel=1e-6)
    assert peters_evolve(1.4, 1.4, 0.003, 0.0, 1.1 * tau)[0] == 0.0  # merged


def test_unresolved_gw_inspiral_stops_the_run():
    from swingersim.backends.nbody import run_nbody
    from swingersim.gw import merger_time_circular_yr

    tau = merger_time_circular_yr(1.4, 1.4, 0.001)
    cfg = Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"out": ["A", "m3"], "A": ["m1", "m2"]},
            "masses": {"m1": 1.4, "m2": 1.4, "m3": 1.0},
            "orbits": {"A": {"a": 0.001, "e": 0.0}, "out": {"a": 100.0, "e": 0.1, "i": 10}},
            "stars": {"m1": {"R": 1e-5}, "m2": {"R": 1e-5}},
            "physics": {"gr": False, "gw": True},
            "run": {"backend": "nbody", "t_max_yr": 0.5 * tau, "n_output": 2},
        }
    )
    res = run_nbody(cfg)
    assert res.event == "gw_unresolved"
    assert len(res.times) == 1 and res.event_detail["stars"] == ["m1", "m2"]


def test_extremes_roche_follow_mass_changes():
    from swingersim.diagnostics.extremes import ExtremeTracker
    from swingersim.hierarchy import build_tree

    root, m = build_tree(_tight_2p2())
    radii = np.full(4, 0.01)
    tr = ExtremeTracker(root, m, units.G, radii)
    before = {s: sc for _, s, sc in tr.roche}
    m2 = m.copy()
    m2[0] *= 0.5  # star 1 loses half its mass: its lobe shrinks, its companion's grows
    tr.set_stars(m2)
    after = {s: sc for _, s, sc in tr.roche}
    assert after[0] > before[0] and after[1] < before[1]


# ------------------------------------------------------------------ binary diagnostics


def test_measured_apsidal_period_nan_and_guard():
    from swingersim.diagnostics.binary import APSIDAL_E_MIN, measured_apsidal_period

    t = np.linspace(0, 100, 401)
    omega = np.mod(2 * np.pi * t / 40.0, 2 * np.pi)
    Omega = np.zeros_like(t)
    assert measured_apsidal_period(t, Omega, omega) == pytest.approx(40.0, rel=1e-9)
    gappy = omega.copy()
    gappy[1::2] = np.nan  # a NaN at every other output used to break np.unwrap
    assert measured_apsidal_period(t, Omega, gappy) == pytest.approx(40.0, rel=1e-9)
    e = np.full_like(t, 0.5 * APSIDAL_E_MIN)
    assert math.isnan(measured_apsidal_period(t, Omega, omega, e))
    assert measured_apsidal_period(t, Omega, omega, e + 0.1) == pytest.approx(40.0, rel=1e-9)


def test_binaries_report_omega_only_period():
    from swingersim.backends.nbody import run_nbody

    res = run_nbody(_tight_2p2())
    for b in res.binaries.values():
        assert "U_measured_omega_yr" in b


# ------------------------------------------------------------------ stability


def test_mlp_missing_model_kind_gives_none(tmp_path):
    from swingersim import stability

    (tmp_path / "classify_quad_2p2.py").write_text(
        "def mlp_classifier_2p2(path, *p):\n    return True\n"
    )
    (tmp_path / "mlp_model_2p2_ghost.pkl").write_text("")
    assert stability._mlp_predict("2+2", tmp_path, {"x": 1}) is True
    assert stability._mlp_predict("3+1", tmp_path, {"x": 1}) is None


def test_ghost_test_uses_resolved_masses(monkeypatch):
    """ghost_test resolves physics.mass_loss masses: current and measures the orbits with
    the masses of the simulation."""
    from swingersim import massloss, stability

    seen = []
    real = massloss.resolve_masses

    def spy(cfg):
        seen.append(cfg)
        return real(cfg)

    monkeypatch.setattr(massloss, "resolve_masses", spy)
    r = stability.ghost_test(_tight_2p2(), n_outer=0.05)
    assert seen and r.reason in ("stable", "diverged", "unbound")
    masses = stability.initial_state(_tight_2p2())[1]
    assert np.array_equal(masses, [1.2, 1.0, 0.9, 0.7])


# ------------------------------------------------------------------ HDF5 and products


def _saved(tmp_path, cfg=None):
    from swingersim.backends.nbody import run_nbody
    from swingersim.io import save_run

    cfg = cfg or _tight_2p2(n_output=200)
    res = run_nbody(cfg)
    h5 = tmp_path / "s.h5"
    save_run(h5, cfg, res)
    return h5, res


def test_hdf5_metadata_compression_and_masses(tmp_path):
    import h5py

    h5, res = _saved(tmp_path)
    with h5py.File(h5) as f:
        for k in (
            "swingersim_version",
            "rebound_version",
            "reboundx_version",
            "numpy_version",
            "created_utc",
            "wall_s",
            "units_json",
        ):
            assert k in f.attrs, k
        assert f["orbits/A/a"].compression == "gzip" and f["orbits/A/a"].shuffle
        assert f["orbits/A/a"].attrs["unit"] == "AU"
        assert f["orbits/A/omega"].attrs["unit"] == "rad"
        assert f["time_yr"].attrs["unit"] == "yr"
        assert f["stars/m1"].attrs["mass_msun"] == 1.2
        assert f["stars/m1"].attrs["mass_msun_config"] == 1.2
        assert "detected" not in f["hierarchy"]
        assert f["hierarchy"].attrs["final"] == res.hierarchy[-1]
        changes = json.loads(f["hierarchy"].attrs["detected_changes_json"])
        assert changes[0]["time_yr"] == 0.0 and changes[0]["to"] == res.hierarchy[0]
        assert np.array_equal(f["orbits/A/e"][()], res.orbits["A"]["e"])


def test_timeseries_header_and_decimation(tmp_path):
    from swingersim.runproducts import _read, write_timeseries

    h5, res = _saved(tmp_path)
    d = _read(h5)
    p = write_timeseries(d, tmp_path / "full.csv")
    text = p.read_text()
    assert (
        text.startswith("# SWINGERsim")
        and "timeseries" in text.splitlines()[0]
        and "# time: t_yr" in text
    )
    body = [line for line in text.splitlines() if not line.startswith("#")]
    rows = list(csv.DictReader(_io.StringIO("\n".join(body))))
    assert len(rows) == 200 and "e_A" in rows[0]
    try:
        import pandas as pd

        assert len(pd.read_csv(p, comment="#")) == 200
    except ImportError:
        pass
    p = write_timeseries(d, tmp_path / "dec.csv", max_rows=60)
    body = [line for line in p.read_text().splitlines() if not line.startswith("#")]
    rows = list(csv.DictReader(_io.StringIO("\n".join(body))))
    assert len(rows) <= 61 and float(rows[-1]["t_yr"]) == pytest.approx(res.times[-1])
    assert "every 4-th of 200 outputs" in p.read_text()


def test_summary_lists_report_only_when_written(tmp_path):
    from swingersim.cli import main
    from swingersim.runproducts import product_paths

    cfg = tmp_path / "sys.yaml"
    cfg.write_text(
        "topology: '2+2'\nmasses: {m1: 1.2, m2: 1.0, m3: 1.1, m4: 0.9}\n"
        "orbits: {A: {P: 2.0}, B: {P: 3.0}, out: {P: 200.0, e: 0.2}}\n"
        "run: {t_max_yr: 0.2, n_output: 21}\n"
    )
    h5 = tmp_path / "a.h5"
    assert main(["run", str(cfg), "-o", str(h5), "--no-viewer", "--no-plots", "--no-html"]) == 0
    assert "a.html" not in product_paths(h5)["summary"].read_text()
    h5 = tmp_path / "b.h5"
    assert main(["run", str(cfg), "-o", str(h5), "--no-viewer", "--no-plots"]) == 0
    assert (tmp_path / "b.html").exists()
    assert "b.html  (report" in product_paths(h5)["summary"].read_text()


def test_run_report_links_viewer_and_references_pngs(tmp_path):
    pytest.importorskip("plotly")
    pytest.importorskip("matplotlib")
    from swingersim.cli import main

    cfg = tmp_path / "sys.yaml"
    cfg.write_text(
        "topology: '2+2'\nmasses: {m1: 1.2, m2: 1.0, m3: 1.1, m4: 0.9}\n"
        "orbits: {A: {P: 2.0}, B: {P: 3.0}, out: {P: 200.0, e: 0.2}}\n"
        "run: {t_max_yr: 0.2, n_output: 21}\n"
    )
    h5 = tmp_path / "r.h5"
    assert main(["run", str(cfg), "-o", str(h5)]) == 0
    page = (tmp_path / "r.html").read_text()
    assert 'href="r_viewer.html"' in page
    assert 'src="r_orbits.png"' in page and '<img alt="Figures" src="data:' not in page
    assert page.count("Plotly.newPlot") >= 1 and page.count("plotly.js v") <= 1
    assert "Teff [K]" in page or "R [Rsun]" in page
    assert "M [Msun]" in page


def test_decimate_keeps_peaks():
    from swingersim.report import decimate

    t = np.arange(100_000.0)
    y = np.zeros_like(t)
    y[12_345] = 1.0
    x2, y2 = decimate(t, y, 5000)
    assert len(x2) <= 5002 and 1.0 in y2 and x2[0] == 0 and x2[-1] == t[-1]
    assert decimate(t[:10], y[:10])[0].tolist() == t[:10].tolist()
