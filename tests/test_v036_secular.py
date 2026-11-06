"""Fixes of the secular backends, discs, phenom, escape theory (review 2026-10-22)."""

from __future__ import annotations

import dataclasses
import math
import warnings

import numpy as np
import pytest

from swingersim.config import Config


def _nt(**over):
    """The configuration of tests/test_nested_tides.py (2+2, tides in both binaries)."""
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


def _reference_rhs(sys_, y):
    """Nested dissipative right-hand side assembled from the public pieces (coefficients,
    the pair kernel, swingersim.tides.leconte_rates, swingersim.gw.peters_rates)."""
    from swingersim.gw import peters_rates
    from swingersim.secular.kernel import secular_rhs_pairs
    from swingersim.tides import leconte_rates

    n6 = 6 * sys_.n_orb
    a = sys_.a_state(y)
    spins = sys_.spins_state(y)
    prec, c2, c3 = sys_.coefficients(a, spins)
    out = np.zeros_like(y)
    out[:n6] = secular_rhs_pairs(y[:n6].copy(), sys_.pk, sys_.pl, c2, c3, prec, sys_.octupole)[1]
    for q in sys_.diss:
        nd = sys_.nodes[q]
        i, j = nd.left, nd.right
        K, p = y[6 * q : 6 * q + 3], y[6 * q + 3 : 6 * q + 6]
        e, eta = float(np.linalg.norm(p)), float(np.linalg.norm(K))
        si, sj = sys_.structure[i], sys_.structure[j]
        r = leconte_rates(
            (sys_.masses[i], sys_.masses[j]),
            (si.r_au, sj.r_au),
            (si.apsidal_k2, sj.apsidal_k2),
            (si.rg2, sj.rg2),
            (si.lag_yr, sj.lag_yr),
            float(a[q]),
            e,
            (spins[i], spins[j]),
        )
        da, de = r.da_dt, r.de_dt
        if sys_.gw:
            ga, ge = peters_rates(sys_.masses[i], sys_.masses[j], float(a[q]), e)
            da, de = da + ga, de + ge
        out[sys_.a_index[q]] = da
        out[sys_.spin_index[q] : sys_.spin_index[q] + 2] = r.dspin_dt
        out[6 * q + 3 : 6 * q + 6] += de * p / e
        out[6 * q : 6 * q + 3] += -(e / eta) * de * K / eta
    return out


@pytest.mark.parametrize("gw", [False, True])
@pytest.mark.parametrize("topology", ["2+2", "3+1"])
def test_nested_dissipative_kernel_matches_reference(topology, gw):
    """The compiled nested right-hand side equals the one assembled from swingersim.tides."""
    from swingersim.backends.secular_nested import build_nested
    from swingersim.secular import kernel

    over = {"physics": {"gr": True, "tides": True, "gw": gw}}
    if topology == "3+1":
        over.update(
            topology="3+1",
            orbits={
                "inner": {"P": 3.0, "e": 0.3},
                "mid": {"P": 90.0, "e": 0.2, "i": 50},
                "out": {"P": 4000.0, "e": 0.3, "i": 20},
            },
        )
    sys_, y0 = build_nested(_nt(**over))[:2]
    assert sys_.dissipative
    rng = np.random.default_rng(4)
    for _ in range(10):
        y = y0 * (1.0 + 1e-2 * rng.standard_normal(len(y0)))
        ref = _reference_rhs(sys_, y)
        new = sys_.rhs(0.0, y)
        assert np.allclose(new, ref, rtol=1e-12, atol=1e-14 * np.max(np.abs(ref)))
    # the same code without numba (the package runs without it)
    py = getattr(kernel.nested_dissipative_rhs, "py_func", kernel.nested_dissipative_rhs)
    assert np.allclose(py(y0.copy(), *sys_._kernel_args()), sys_.rhs(0.0, y0), rtol=1e-13)


def test_nested_kernel_arguments_follow_the_structure():
    """A new structure (evolving stars) or masses rebuild the cached kernel arguments."""
    from swingersim.backends.secular_nested import build_nested

    sys_, y0 = build_nested(_nt())[:2]
    r0 = sys_.rhs(0.0, y0)
    st = {k: dataclasses.replace(v, lag_yr=2.0 * v.lag_yr) for k, v in sys_.structure.items()}
    sys_.set_structure(st)
    r1 = sys_.rhs(0.0, y0)
    q = sys_.diss[0]
    assert r1[sys_.a_index[q]] == pytest.approx(2.0 * r0[sys_.a_index[q]], rel=1e-12)
    sys_.scale_orbits(y0, {sys_.names[q]: 1.0})  # in-place a0 change invalidates as well
    assert sys_._kargs is None


def test_leconte_rates_lag0_and_rg2_0_are_consistent():
    """A star with lag 0 does not dissipate and keeps its spin; a dissipating star without a
    moment of inertia (rg2 = 0) gives finite rates (pseudo-synchronous limit), the same in
    the 2+2 and the nested kernel."""
    from swingersim.backends.secular import build_secular
    from swingersim.backends.secular_nested import build_nested
    from swingersim.secular.kernel import leconte_pair

    G = 39.47841760435743
    # star 0: rg2 = 0, lag > 0; star 1: rg2 > 0, lag = 0
    args = [G, 1.0, 0.8, 0.005, 0.004, 0.02, 0.02, 0.0, 0.08, 1e-4, 0.0, 0.04, 0.3, 300.0, 250.0]
    da, de, dw0, dw1 = leconte_pair(*args)
    assert all(math.isfinite(v) for v in (da, de, dw0, dw1))
    assert dw0 == 0.0 and dw1 == 0.0 and da != 0.0 and de != 0.0
    other = list(args)
    other[6], other[8], other[14] = 0.5, 0.3, 1.0  # star 1 (lag 0): k2, rg2, spin irrelevant
    assert leconte_pair(*other) == (da, de, dw0, dw1)

    cfg = _nt()
    s2, y2 = build_secular(cfg)[:2]
    sn, yn = build_nested(cfg)[:2]
    # m1: no moment of inertia, m2: no dissipation, in both models
    b = s2.binaries["A"]
    s2.set_structure(
        {**s2.binaries, "A": dataclasses.replace(b, rg2=(0.0, b.rg2[1]), lag_yr=(b.lag_yr[0], 0))}
    )
    st = dict(sn.structure)
    st[0] = dataclasses.replace(st[0], rg2=0.0)
    st[1] = dataclasses.replace(st[1], lag_yr=0.0)
    sn.set_structure(st)
    r2, rn = s2.rhs(0.0, y2), sn.rhs(0.0, yn)
    assert np.all(np.isfinite(r2)) and np.all(np.isfinite(rn))
    q = sn.names.index("A")
    assert rn[sn.a_index[q]] == pytest.approx(r2[18], rel=1e-12)
    assert rn[sn.spin_index[q]] == r2[20] == 0.0
    assert rn[sn.spin_index[q] + 1] == r2[21] == 0.0


def test_twotwo_rhs_has_no_dead_branch_and_integrate_is_shared():
    from swingersim.secular import nested, twotwo

    assert nested.integrate is twotwo.integrate
    assert not hasattr(twotwo, "sl_p") and not hasattr(twotwo, "leconte_rates")
    assert not hasattr(nested, "leconte_rates_fast")


def _zlk_cfg(t_max, **run):
    return _nt(
        stars={k: {"lag_s": 0.0} for k in ("m1", "m2", "m3", "m4")},
        physics={"gr": False, "tides": False},
        orbits={
            "A": {"P": 3.0, "e": 0.1},
            "B": {"P": 10.0, "e": 0.1, "i": 87},
            "out": {"P": 800.0, "e": 0.1},
        },
        run={"backend": "secular", "t_max_yr": t_max, "n_output": 201, **run},
    )


def test_backward_run_threshold_crossing_is_the_first_from_t0():
    """Backward runs: the first crossing of a 1 - e threshold is the one closest to t = 0 (was
    the one farthest from it within the first chunk, -9918 yr instead of -84 yr)."""
    from swingersim.backends.secular import run_secular

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = run_secular(_zlk_cfg(-2e4))
    t_x = r.extremes["B"]["t_1me_lt_1e-1_yr"]
    e, t = r.orbits["B"]["e"], r.times
    assert t_x < 0
    assert np.all(e[np.abs(t) < abs(t_x)] < 0.9)  # nothing above before the crossing
    first_above = t[np.argmax(e > 0.9)]
    assert abs(t_x) <= abs(first_above)


def test_first_crossing_direction():
    from swingersim.backends.secular_common import first_crossing

    t = np.linspace(-10.0, 0.0, 11)
    f = np.array([0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0], float)
    lin = lambda x: float(np.interp(x, t, f))
    assert first_crossing(t, f, lin, 1.0, forward=True) == pytest.approx(-9.5)
    assert first_crossing(t, f, lin, 1.0, forward=False) == pytest.approx(-3.5)
    assert first_crossing(t, np.full(11, 2.0), lin, 1.0, forward=False) == 0.0
    assert math.isnan(first_crossing(t, np.zeros(11), lin, 1.0))


def test_secular_roche_overflow_truncates_the_run():
    """Non-evolving secular runs with stop_on_roche_overflow end at the first overflow like
    the N-body backend (the series went on to t_max before)."""
    from swingersim.backends.secular import run_secular

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = run_secular(_zlk_cfg(2e5, stop_on_roche_overflow=True))
    assert r.event == "roche_overflow"
    assert r.event_detail["stars"]
    assert len(r.times) < 201
    assert r.times[-1] <= max(r.event_time, r.times[1])
    for o in r.orbits.values():
        assert len(o["e"]) == len(r.times)
    assert len(r.energy_error) == len(r.times)
    # the overflow is the first one: the filling first reaches 1 there
    assert 0.0 < r.event_time < r.extremes["B"]["t_e_max_yr"] + 1e-9


def test_nested_encounters_with_evolving_stars_at_their_time(monkeypatch):
    """Encounters with evolving stars are applied at their own time: the segment edges are cut
    there (they were applied at the start of the segment containing them)."""
    from swingersim.backends import evolving as ev
    from swingersim.backends import secular_nested as sn
    from swingersim.environment import Encounter

    t_enc = [1234.5678, 4321.0]
    seen = []
    monkeypatch.setattr(
        sn,
        "_encounters",
        lambda cfg, sys_: [Encounter(t, 1e9, 1.5, np.eye(3), 0.1) for t in t_enc],
    )
    monkeypatch.setattr(
        sn, "_kick", lambda sys_, y, encs, counts: seen.extend(ev.SEGMENT["t0"] for _ in encs)
    )
    orig = ev.segment_edges
    cfg = Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"out": ["A", "m3"], "A": ["m1", "m2"]},
            "masses": {"m1": 1.0, "m2": 0.9, "m3": 0.8},
            "orbits": {"A": {"a": 1.0}, "out": {"a": 300.0, "e": 0.3, "i": 40.0}},
            "stellar_model": {"kind": "sse", "age_myr": 100.0},
            "physics": {"gr": False, "evolution": True},
            "run": {"backend": "secular", "t_max_yr": 1e4, "n_output": 11},
        }
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sn.run_secular_nested(cfg)
    assert seen == t_enc
    assert ev.segment_edges is orig  # restored


def test_disc_light_emission_cache_gives_identical_values():
    from pathlib import Path

    import yaml

    from swingersim import discs as dm
    from swingersim.hierarchy import build_tree, initialize
    from swingersim.photometry import light_curve

    path = Path(__file__).resolve().parents[1] / "configs/catalog/eps_aur.yaml"
    d = yaml.safe_load(path.read_text(encoding="utf-8"))
    d["discs"]["disc"].update({"tau": 3.0, "edge_au": 0.2, "emission": True, "albedo": 0.0})
    cfg = Config.model_validate(d)
    lc = light_curve(cfg, np.array([100.0, 900.0]), "V")
    assert lc.components["disc_disc"][0] > 0
    root, masses = build_tree(cfg)
    disc = dm.build_discs(cfg, root)[0]
    pos = initialize(root, masses)[0]
    Light = type("Light", (), {})
    lights = []
    for k in range(len(masses)):
        x = Light()
        x.radius_au, x.weight = 0.3 / (k + 1), 1.0
        lights.append(x)
    teffs = [7000.0 + 4000.0 * k for k in range(len(masses))]
    first = dm.disc_light([disc], 0.0, pos, masses, lights, teffs, "V")
    assert len(disc._emission) == 1
    second = dm.disc_light([disc], 0.0, pos, masses, lights, teffs, "V")
    disc._emission = None  # recomputed from scratch
    fresh = dm.disc_light([disc], 0.0, pos, masses, lights, teffs, "V")
    assert first[0] > 0 and first[0] == second[0] == fresh[0]
    teffs[disc.stars[0]] = 8000.0  # another host temperature: another key
    dm.disc_light([disc], 0.0, pos, masses, lights, teffs, "V")
    assert len(disc._emission) == 2


def test_phenom_minima_csv_is_read_by_swingersim_fit(tmp_path):
    from swingersim import ocfit, phenom

    mins = [
        phenom.Minimum("A", "primary", 0.0, 1.0, 1.00001, 1e-5, 0.3, 1e-4),
        phenom.Minimum("A", "secondary", 0.5, 3.5, 3.50002, 2e-5, 0.2, 1e-4),
    ]
    res = phenom.MinimaResult(mins, {}, meta={"seed": 7, "swingersim_version": "x"})
    path = tmp_path / "m.csv"
    phenom.write_minima_csv(res, path, config="cfg.yaml", bjd_t0=2459000.5)
    m = ocfit.read_minima(path)
    assert m.header["time_system"] == "BJD_TDB"
    assert float(m.header["bjd_t0"]) == 2459000.5
    assert m.header["config"] == "cfg.yaml" and m.header["seed"] == "7"
    assert list(m.pair) == ["A", "A"] and list(m.kind) == [0, 1]
    assert m.time == pytest.approx([2459001.50001, 2459004.00002], abs=1e-8)
    assert m.sigma == pytest.approx([1e-5, 2e-5])
    # the default time origin: fit.stage2.bjd_t0 of the configuration, else J2000.0
    phenom.write_minima_csv(res, path)
    assert ocfit.read_minima(path).time[0] == pytest.approx(phenom.DEFAULT_BJD_T0 + 1.00001)
    from types import SimpleNamespace

    cfg = SimpleNamespace(fit=SimpleNamespace(stage2=SimpleNamespace(bjd_t0=2458000.0)))
    assert phenom._bjd_t0_of(cfg, None) == 2458000.0
    assert phenom._bjd_t0_of(_nt(), None) == phenom.DEFAULT_BJD_T0
    assert phenom._bjd_t0_of(cfg, 1.0) == 1.0


def test_escape_outcomes_in_slices_equal_the_full_grid(monkeypatch):
    from swingersim import escape_theory as et
    from swingersim import units

    E0 = -units.G / 2.0
    args = ((1.0, 1.0, 0.6), E0, 0.4 * et.angular_momentum_scale((1, 1, 0.6), E0))
    kw = {"nu": 40, "nl": 50, "nc": 16}
    sliced = et.escape_outcomes(*args, **kw)
    monkeypatch.setattr(et, "CHUNK_CELLS", 10**9)  # one slice: the former single grid
    full = et.escape_outcomes(*args, **kw)
    for k in ("e_density_all", "u_density_all"):
        assert np.allclose(getattr(sliced, k), getattr(full, k), rtol=1e-12, atol=0)
    assert sliced.mean_e == pytest.approx(full.mean_e, rel=1e-12)
    assert sliced.mean_u == pytest.approx(full.mean_u, rel=1e-12)
    for s in range(3):
        assert sliced.probabilities[s] == pytest.approx(full.probabilities[s], rel=1e-12)


def test_quadrupole_period_with_one_cycle():
    from swingersim.secular.quadrupole import predict

    args = (1.0, 1.0, 1.0, 1.0, 20.0, 0.05, 0.1, 70.0, 0.0)
    one = predict(*args, n_cycles=1)
    three = predict(*args, n_cycles=3)
    assert math.isfinite(one.period_yr)
    assert one.period_yr == pytest.approx(three.period_yr, rel=1e-3)


# ---------------------------------------------------------------- tools/ (report builders)
def _tools():
    import sys
    from pathlib import Path

    d = str(Path(__file__).resolve().parents[1] / "tools")
    if d not in sys.path:
        sys.path.insert(0, d)
    import _common

    return _common


def test_tools_json_cache_invalidates(tmp_path, monkeypatch):
    import os

    import swingersim

    c = _tools()
    data = tmp_path / "data.csv"
    data.write_text("a\n1\n")
    calls = []

    def compute():
        calls.append(1)
        return len(calls)

    cache = c.JsonCache(tmp_path / "cache")
    assert cache.get("x", compute, [data], {"t0": 1.0}) == 1
    assert cache.get("x", compute, [data], {"t0": 1.0}) == 1  # cached
    assert cache.get("x", compute, [data], {"t0": 2.0}) == 2  # another request
    st = data.stat()
    os.utime(data, ns=(st.st_atime_ns, st.st_mtime_ns + 10**9))  # the data changed
    assert cache.get("x", compute, [data], {"t0": 2.0}) == 3
    monkeypatch.setattr(swingersim, "__version__", "9.9.9")  # another swingersim
    assert cache.get("x", compute, [data], {"t0": 2.0}) == 4
    assert cache.get("x", compute, [data], {"t0": 2.0}) == 4
    assert not [p for p in (tmp_path / "cache").iterdir() if p.suffix == ".tmp"]
    # a configuration dict modified in memory is part of the key through its digest
    assert c.digest({"a": 1, "b": [1, 2]}) == c.digest({"b": [1, 2], "a": 1})
    assert c.digest({"a": 1}) != c.digest({"a": 2})


def test_tools_file_cache_and_epoch(tmp_path):
    c = _tools()
    src = tmp_path / "cfg.yaml"
    src.write_text("# t = 0 of this configuration = BJD 2459383.9.\nmasses: {m1: 1}\n")
    out = tmp_path / "lc.csv"
    n = []

    def write(tmp):
        n.append(1)
        tmp.write_text("t,flux\n0,1\n")

    assert c.file_cached(out, write, [src], {"t0": 1}) is True
    assert c.file_cached(out, write, [src], {"t0": 1}) is False
    src.write_text(src.read_text() + "# changed\n")
    assert c.file_cached(out, write, [src], {"t0": 1}) is True
    assert len(n) == 2 and out.read_text().startswith("t,flux")
    assert c.fit_epoch(src) == 2459383.9
    src.write_text("# t = 0 = BJD 2444500 (1980.6).\nfit: {stage2: {bjd_t0: 2458000.5}}\n")
    assert c.fit_epoch(src) == 2458000.5  # fit.stage2.bjd_t0 before the comment
    src.write_text("fit_epoch_bjd_tdb: 2450000.25\nfit: {stage2: {bjd_t0: 2458000.5}}\n")
    assert c.fit_epoch(src) == 2450000.25
    src.write_text("# t = 0 = BJD 2444500 (1980.6).\n")
    assert c.fit_epoch(src) == 2444500.0
    src.write_text("masses: {m1: 1}\n")
    with pytest.raises(ValueError):
        c.fit_epoch(src)


def test_tools_failed_png_is_removed_and_fatal(tmp_path):
    c = _tools()

    class Bad:
        def to_json(self):
            return "{}"

        def write_image(self, path, scale=2.0):
            raise RuntimeError("no kaleido")

    (tmp_path / "f.png").write_bytes(b"old")
    w = c.FigureWriter(tmp_path)
    w.save_plotly(Bad(), "f")
    assert not (tmp_path / "f.png").exists() and (tmp_path / "f.json").exists()
    with pytest.raises(c.FigureFailures):
        w.finish()
    (tmp_path / "f.png").write_bytes(b"old")
    w = c.FigureWriter(tmp_path, png=False)  # --no-png: JSON only, not an error
    w.save_plotly(Bad(), "f")
    assert not (tmp_path / "f.png").exists()
    w.finish()


def test_tools_parse_report_is_robust(tmp_path):
    _tools()
    import squadra_fits

    p = tmp_path / "stage2_report.txt"
    p.write_text(
        "  0 minima, 120 radial velocities, 400 light-curve points\n"
        "  chi^2 1.5e+03 -> 2.4e-05 (reduced 1.2e-07; priors 0)\n"
        "  errors (per_source), jitter added in quadrature to the quoted errors: "
        "Src A 1.5 min; rv:CTIO 18.6 km/s; lc:TESS:S34 3.05 ppt; lc:TESS:a 1e-3 ppt\n"
        "  light curves, rms of measured - calibrated model [ppt], start -> fitted, scale a:\n"
        "    TESS:S34     2080 points     3.58 ->    3.04   a = 0.995\n"
    )
    r = squadra_fits.parse_report(p)  # no section of minima: was an IndexError
    assert r["chi2_start"] == 1500.0 and r["chi2"] == 2.4e-05 and r["chi2_red"] == 1.2e-07
    assert r["rms"] == {}
    assert r["jitter"] == {"Src A": 1.5, "rv:CTIO": 18.6, "lc:TESS:S34": 3.05, "lc:TESS:a": 1e-3}
    assert r["lc"]["TESS:S34"] == {"n": 2080, "start": 3.58, "fit": 3.04, "a": 0.995}


def test_tools_light_curve_helpers(tmp_path):
    _tools()
    import make_lc_windows
    import tess_minima

    rng = np.random.default_rng(1)
    t = np.sort(rng.uniform(0.0, 100.0, 3000))
    tm = np.sort(rng.uniform(-5.0, 105.0, 200))
    ref = np.min(np.abs(t[:, None] - tm[None, :]), axis=1)
    assert np.array_equal(make_lc_windows.distance_to_nearest(t, tm), ref)
    assert np.all(np.isinf(make_lc_windows.distance_to_nearest(t, np.array([]))))
    f = tmp_path / "tess_s26_spoc_noosc.csv"  # '#' header lines, as remove_oscillation writes
    tt = np.linspace(2459010.0, 2459012.0, 500)
    rows = "\n".join(f"{x:.7f},{1.0 + 1e-4 * np.sin(x)},0.001" for x in tt)
    f.write_text("# oscillation removed\n# rms ...\nbjd_tdb,flux,flux_err\n" + rows + "\n")
    t1, f1 = tess_minima.load(f)  # np.genfromtxt took the first comment line as the names
    assert len(t1) == 500 and abs(f1.mean() - 1.0) < 1e-3


def test_tools_fit_table_offsets_wrap_angles():
    _tools()
    import squadra_figures

    assert squadra_figures._offset(359.0, 1.0, "orbits.A.omega") == pytest.approx(-2.0)
    assert squadra_figures._offset(2.5, 1.0, "masses.m3") == pytest.approx(1.5)
