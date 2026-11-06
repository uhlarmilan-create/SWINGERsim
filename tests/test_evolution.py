"""F4 (a): stellar evolution tracks (SSE, Hurley et al. 2000) and evolving stars in the
secular backends."""

import numpy as np
import pytest

from swingersim.config import Config
from swingersim.evolution import EvolvedBeyondGrid, track
from swingersim.stellar import star_properties, zams_radius_tout


def test_zams_matches_tout():
    for m in (0.3, 0.8, 1.0, 1.7, 2.5, 4.0):
        assert track(m, 0.02).at(0.0).radius_rsun == pytest.approx(zams_radius_tout(m), rel=2e-3)


def test_sun():
    s = track(1.0, 0.02).at(4570.0)
    assert s.radius_rsun == pytest.approx(1.0, rel=0.05)
    assert s.luminosity_lsun == pytest.approx(1.0, rel=0.05)
    assert s.teff_k == pytest.approx(5772.0, rel=0.03)
    assert 0.0 < s.convective_turnover_yr < 0.2


def test_phases_and_growth():
    t = track(2.0, 0.02)
    assert 0 < t.t_ms < t.t_bgb < t.t_hei
    r = [t.at(a).radius_rsun for a in np.linspace(0, 0.99 * t.t_ms, 20)]
    assert np.all(np.diff(r) > -1e-6) and r[-1] > 1.8 * r[0]
    assert t.at(0.5 * (t.t_ms + t.t_bgb)).phase == 2
    assert t.at(0.5 * (t.t_bgb + t.t_hei)).phase == 3
    with pytest.raises(EvolvedBeyondGrid):
        t.at(1.01 * t.t_hei)
    assert track(2.0, 0.008).t_ms < t.t_ms  # metal-poor stars live shorter


def test_kraft_break_envelopes():
    assert track(1.0, 0.02).at(100).m_env > 0.01
    assert track(2.0, 0.02).at(100).m_env < 1e-4


def _cfg(**over):
    d = {
        "topology": "2+2",
        "masses": {"m1": 1.4, "m2": 1.2, "m3": 1.1, "m4": 1.0},
        "orbits": {
            "A": {"P": 2.5, "e": 0.05, "i": 85},
            "B": {"P": 3.0, "e": 0.05, "i": 70, "Omega": 30},
            "out": {"P": 3000, "e": 0.3, "i": 80, "Omega": 10},
        },
        "physics": {"tides": True, "evolution": True},
        "stellar_model": {"kind": "sse", "age_myr": 1500},
        "run": {"t_max_yr": 2.0e5, "n_output": 50, "backend": "secular"},
    }
    for k, v in over.items():
        d[k] = {**d.get(k, {}), **v} if isinstance(v, dict) else v
    return Config.model_validate(d)


def test_config_rules():
    with pytest.raises(ValueError):
        Config.model_validate({**_cfg().model_dump(exclude={"ensemble"}), "stellar_model": {}})
    with pytest.raises(ValueError):
        _cfg(stellar_model={"kind": "sse", "age_myr": None})
    with pytest.raises(ValueError):
        _cfg(stars={"m1": {"R": 1.5}})


def test_sse_properties():
    c = _cfg(physics={"evolution": False})
    p = star_properties(c, "m1", 1.4)
    assert p.phase == "main sequence" and p.age_myr == 1500
    assert p.radius_rsun > zams_radius_tout(1.4)
    assert p.teff_k == pytest.approx(5772 * p.luminosity_lsun**0.25 / p.radius_rsun**0.5)


def test_evolving_run_2p2():
    from swingersim.backends.secular import run_secular

    c = _cfg()
    r = run_secular(c)
    assert r.event is None
    assert r.energy_error.max() < 1e-7
    rad = r.stellar["m1"]["radius_rsun"]
    assert rad[-1] > rad[0]
    static = run_secular(_cfg(physics={"evolution": False}))
    # 0.2 Myr: the stars hardly change, both runs agree
    assert r.orbits["A"]["e"][-1] == pytest.approx(static.orbits["A"]["e"][-1], abs=1e-4)


def test_end_of_grid_and_birth():
    from swingersim.backends.secular import run_secular

    t_hei = track(1.4, 0.02).t_hei
    c = _cfg(stellar_model={"kind": "sse", "age_myr": t_hei - 0.05}, run={"t_max_yr": 1e6})
    c.orbits["A"].P = 5000.0  # a wide pair, no Roche-lobe overflow of the giant (~150 Rsun)
    c.orbits["B"].P = 5000.0
    c.orbits["out"].P = 3.0e6
    r = run_secular(c)
    assert r.event in ("stellar_evolution_end", "roche_overflow")
    assert r.times[-1] <= 1e6
    b = run_secular(_cfg(stellar_model={"kind": "sse", "age_myr": 0.1}, run={"t_max_yr": -1e6}))
    assert b.event == "stellar_birth"
    assert b.times[-1] == pytest.approx(-1e5, rel=1e-6)


def test_nbody_evolving_stars_short_run():
    """Since 0.14 the N-body backend evolves the stars output step by output step."""
    from swingersim.backends.nbody import run_nbody

    res = run_nbody(_cfg(run={"backend": "nbody", "t_max_yr": 1.0, "n_output": 3}))
    assert res.event is None
    assert set(res.stellar["m1"]) >= {"radius_rsun", "age_myr"}


def test_hurley_tidal_lags():
    from swingersim.backends.nbody import model_props
    from swingersim.hierarchy import build_tree, initialize
    from swingersim.tides import hurley_k_over_t, lag_from_k_over_t

    c = _cfg(
        physics={"evolution": False, "tidal_lag": "hurley2002"},
        masses={"m1": 1.0, "m2": 0.9, "m3": 2.0, "m4": 1.8},
        stellar_model={"kind": "sse", "age_myr": 500},
    )
    root, m = build_tree(c)
    props = model_props(c, root, m, *initialize(root, m))
    assert all(p.lag_s > 0 for p in props)
    # convective envelopes dissipate far more than radiative ones (per unit structure)
    assert props[0].lag_s > 100 * props[2].lag_s
    # (k/T) of a synchronised solar-type star: 2/21 M_env/(M tau_conv), Hurley et al. (2002)
    kt = hurley_k_over_t(1.0, 1.0, 0.03, 0.06, 1.0, 0.05, 0.01, 2 * np.pi / 0.01)
    assert kt == pytest.approx(2 / 21 * 0.03 / 0.06, rel=1e-12)
    assert lag_from_k_over_t(kt, 1.0, 1.0, 0.0144) > 0
    with pytest.raises(ValueError):
        _cfg(physics={"tides": False, "tidal_lag": "hurley2002", "evolution": False})


def test_claret2023_structure():
    from swingersim.stellar import claret_covers, claret_structure, polytrope_constants

    k2, rg2 = claret_structure(2.0, 4.069, 0.0)  # table row, Claret (2023), 2 Msun, Z = 0.0134
    assert np.log10(k2) == pytest.approx(-2.478, abs=2e-3)
    assert np.sqrt(rg2) == pytest.approx(0.1907, abs=5e-4)
    # between grid masses and metallicities the values lie between the neighbours
    lo, hi = claret_structure(1.0, 4.4, 0.0)[0], claret_structure(1.2, 4.4, 0.0)[0]
    mid = claret_structure(1.1, 4.4, 0.0)[0]
    assert min(lo, hi) <= mid <= max(lo, hi)
    # the main sequence concentrates the star (k2 drops), the deep convective envelope of a
    # giant raises k2 and I/(MR^2) again (Claret 2023, Fig. 1)
    assert claret_structure(1.6, 4.0)[0] < claret_structure(1.6, 4.3)[0]
    assert claret_structure(1.6, 3.3)[0] > 10 * claret_structure(1.6, 3.8)[0]
    assert not claret_covers(0.1) and claret_covers(1.0)
    with pytest.raises(ValueError):
        claret_structure(10.0, 4.0)
    # sse uses Claret by default (k2 of a 2 Msun star ~3x below the n = 3 polytrope)
    c = _cfg(physics={"evolution": False}, stellar_model={"kind": "sse", "age_myr": 500})
    p = star_properties(c, "m1", 2.0)
    assert p.apsidal_k2 < 0.6 * polytrope_constants(3.0)[0]
    poly = _cfg(
        physics={"evolution": False},
        stellar_model={"kind": "sse", "age_myr": 500, "structure": "polytrope"},
    )
    assert star_properties(poly, "m1", 2.0).apsidal_k2 == pytest.approx(polytrope_constants(3.0)[0])
    with pytest.raises(ValueError):
        _cfg(stellar_model={"kind": "table", "file": "x.csv", "structure": "claret2023"})


def test_segments_resolve_thermal_pulses():
    """Regression (0.26.1): the segment grid must follow the radius through the thermal pulses
    of the TP-AGB, where the radius peaks at ~500 Rsun between two nodes of the MIST track; a
    uniform grid stepped over the peaks (313 Rsun seen) and the Roche overflow was missed."""
    from pathlib import Path

    from swingersim.backends.evolving import segment_edges, structure_series, time_limit
    from swingersim.config import load_config
    from swingersim.hierarchy import build_tree

    cfg = load_config(Path(__file__).parent.parent / "configs/examples/example_giant_triple.yaml")
    m = np.array(build_tree(cfg)[1], float)
    t_end, _ = time_limit(cfg, m, cfg.run.t_max_yr)
    e = segment_edges(cfg, m, t_end, cfg.run.evolution_max_dr)
    late = e[e > t_end - 4e6]
    mids = 0.5 * (late[1:] + late[:-1])
    r0 = structure_series(cfg, m, late[:-1])[6]
    rm = structure_series(cfg, m, mids)[6]
    assert np.max(np.abs(np.log(rm / r0))) <= 2 * cfg.run.evolution_max_dr
    assert structure_series(cfg, m, e)[6].max() > 450.0
