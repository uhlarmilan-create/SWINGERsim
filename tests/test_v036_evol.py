"""Regression tests of the v0.36 fixes in stellar evolution coupling and mass transfer."""

import math
import warnings
from types import SimpleNamespace

import numpy as np
import pytest

from swingersim import masstransfer as mt
from swingersim import units
from swingersim.config import Config

warnings.simplefilter("ignore")


def _binary(m1, m2, p_day, age, t_max, vvcrit=0.0, transfer=True, **phys):
    if transfer:
        phys["mass_transfer"] = {}
    return Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"A": ["m1", "m2"]},
            "masses": {"m1": m1, "m2": m2},
            "orbits": {"A": {"P": p_day}},
            "stellar_model": {"kind": "mist", "age_myr": age, "vvcrit": vvcrit},
            "physics": {"gr": False, "evolution": True, **phys},
            "run": {"t_max_yr": t_max, "backend": "secular"},
        }
    )


def _a0(m_tot, p_day):
    return (units.G * m_tot * (p_day / units.DAYS_PER_YEAR) ** 2 / (4 * math.pi**2)) ** (1 / 3)


@pytest.fixture(scope="module")
def case_b():
    """Case-B binary 3 + 2 Msun, P = 3 d, up to shortly after the onset (~4 s)."""
    cfg = _binary(3.0, 2.0, 3.0, 355.0, 5.2e6)
    return cfg, mt.binary_history(cfg, (0, 1), [3.0, 2.0], _a0(5.0, 3.0), 5.2e6)


# ---------------------------------------------------------------- onset, episodes, summary
def test_onset_is_the_root_of_the_overflow(case_b):
    """t_onset was the first solution point in transfer (up to max_step = t_end / 2000 late);
    it is now the root of max(R / R_L) = 1 (non-terminal ODE event)."""
    _, h = case_b
    first = h.t[np.argmax(h.fill > 1.0)]
    assert h.t_onset < first
    assert h.fill_at(h.t_onset) == pytest.approx(1.0, abs=1e-5)
    assert h.donor == 0 and h.donor_core_onset is not None


def test_episode_summary_is_written(case_b):
    _, h = case_b
    assert h.phases and h.phases[0][1] in ("nuclear", "thermal")
    assert h.peak_rate > 0 and h.t_peak_rate >= h.t_onset
    s = h.episode_summary()
    assert s["regimes"].startswith(h.phases[0][1] + "@")
    assert s["peak_regime"] == h.peak_regime


def test_time_grid_is_thinned_but_keeps_the_episode(case_b, monkeypatch):
    """time_grid returned every LSODA point (~15000 for case B); it keeps a point only when a
    mass, separation or radius correction moved by evolution_max_dr / 4, plus the episode
    times."""
    cfg, h = case_b
    monkeypatch.setattr(mt, "histories", lambda c: {"A": h})
    tg = mt.time_grid(cfg)
    assert len(tg) < len(h.t) / 5
    for t in (h.t[0], h.t[-1], h.t_onset, h.t_peak_rate, *(t for t, _ in h.phases)):
        assert np.any(tg == t)
    thr = 0.25 * cfg.run.evolution_max_dr
    rows = np.vstack([np.log(h.m), np.log(h.m_ev), h.r_corr, np.log(h.a)])
    kept = tg[np.isin(tg, h.t)]
    k = np.searchsorted(kept, h.t, side="right") - 1  # last kept point at or before each one
    j = np.searchsorted(h.t, kept[k])
    assert np.max(np.abs(rows - rows[:, j])) < thr


def test_roche_exclusion_follows_the_geometry(monkeypatch):
    """roche_excluded accepted every overflow within 0.5 % of the onset time (5 Myr at 1 Gyr);
    before the onset it now needs the circular history to fill the lobe to (1 - e)(1 - 2 %)."""
    t = np.linspace(0.0, 1.0e9, 1001)
    h = mt.History(pair=(0, 1), t=t, m=np.ones((2, 1001)), a=np.ones(1001))
    h.t_onset, h.donor = 1.0e9, 0
    h.fill = np.full(1001, 0.9)
    h.fill[-1] = 1.0  # 0.95 at 0.9995 Gyr, 0.99 at 0.99999 Gyr
    monkeypatch.setattr(mt, "histories", lambda c: {"A": h})
    assert mt.roche_excluded(None, 0, 1.0e9)
    assert not mt.roche_excluded(None, 0, 0.9995e9)  # was excluded (within 0.5 % of t_onset)
    assert mt.roche_excluded(None, 0, 0.99999e9)  # fill 0.99 >= 0.98
    assert mt.roche_excluded(None, 0, 0.9995e9, e=0.05)  # reaches the lobe at periastron
    assert not mt.roche_excluded(None, 2, 1.0e9)


# ---------------------------------------------------------------- failed integration
def test_failed_integration_ends_the_history(monkeypatch):
    """solve_ivp failures were ignored and the history padded to t_end."""
    real = mt._integrate

    def failing(*args, **kw):
        sol = real(*args, **kw)
        sol.status, sol.message = -1, "stand-in failure"
        return sol

    monkeypatch.setattr(mt, "_integrate", failing)
    cfg = _binary(1.0, 0.8, 1.0e5, 1000.0, 1.0e8)
    with pytest.warns(UserWarning, match="failed"):
        h = mt.binary_history(cfg, (0, 1), [1.0, 0.8], _a0(1.8, 1.0e5), 1.0e8)
    assert h.event == "integration_failed" and h.donor is None
    assert h.t[-1] == h.t_event  # not padded
    monkeypatch.setattr(mt, "histories", lambda c: {"A": h})
    t_stop, what = mt.stop(cfg)  # no donor: used to fail on h.donor + 1
    assert t_stop == h.t_event and what.startswith("integration_failed (orbit A)")


# ---------------------------------------------------------------- re-modelled stars
def test_grid_tracks_at_the_top_of_the_grid():
    """The upper neighbour of the top grid mass does not exist (ValueError before)."""
    from swingersim.evolution import grid_limits

    sm = SimpleNamespace(kind="sse", z=0.02, feh=None, vvcrit=None)
    top = grid_limits()["mass"][1]
    t0, t1, w = mt._grid_tracks(sm, top)
    assert t0 is t1 and w == 0.0
    assert mt.remodelled_state(sm, top, 1.0).radius_rsun > 0


def test_mist_wind_is_interpolated_between_the_grid_tracks():
    """The wind of a re-modelled MIST star came from a new MistTrack per 1e-4 Msun (a step
    function of m_ev that evicted the grid tracks from the cache); it is now interpolated
    between the M_STEP grid tracks at the same evolutionary coordinate."""
    from swingersim.stellar import stellar_track

    cfg = _binary(25.0, 18.0, 6.0, 1.0, 8.0e6, vvcrit=0.4, mass_loss={"source": "track"})
    sm = cfg.stellar_model
    age = 5.0
    tr = stellar_track(sm, 25.0)
    st = mt.remodelled_state(sm, 25.0, age)
    w_node = mt.wind_rate(cfg, 0, 20.0, st, age, 25.0)  # grid node: that track exactly
    assert w_node == pytest.approx(20.0 * tr.mdot(age) / tr.current_mass(age), rel=1e-10)
    assert mt.track_label(cfg, 20.0, 25.0, age) == pytest.approx(
        20.0 * 25.0 / tr.current_mass(age), rel=1e-10
    )
    # continuous in m_ev (no 1e-4 Msun steps)
    ms = 25.0 + np.linspace(0.0, 2 * mt.M_STEP, 41)
    rates = np.array([mt.wind_rate(cfg, 0, 20.0, st, age, m) for m in ms])
    assert np.max(np.abs(np.diff(rates))) < 0.01 * abs(rates[0])


# ---------------------------------------------------------------- segments and outputs
def test_segment_edges_take_extra_edges():
    from swingersim.backends.evolving import segment_edges

    cfg = _binary(1.2, 1.0, 5.0, 2000.0, 1.0e9, transfer=False)
    base = segment_edges(cfg, [1.2, 1.0], 1.0e9, 0.01)
    extra = [3.3e8, 7.77e8, 2.0e9, -1.0]  # the last two lie outside the run
    edges = segment_edges(cfg, [1.2, 1.0], 1.0e9, 0.01, extra_edges=extra)
    assert set(edges) == set(base) | {3.3e8, 7.77e8}
    assert np.all(np.diff(edges) > 0)


def test_stellar_series_vectorised_matches_the_track():
    from swingersim.backends.evolving import stellar_series
    from swingersim.stellar import stellar_track

    cfg = _binary(1.2, 1.0, 5.0, 2000.0, 3.0e9, transfer=False)
    times = np.linspace(0.0, 3.0e9, 501)
    out = stellar_series(cfg, [1.2, 1.0], times)
    tr = stellar_track(cfg.stellar_model, 1.2)
    for j in (0, 137, 400, 500):
        s = tr.at(2000.0 + times[j] / 1e6)
        assert out["m1"]["radius_rsun"][j] == pytest.approx(s.radius_rsun, rel=1e-14)
        assert out["m1"]["teff_k"][j] == pytest.approx(s.teff_k, rel=1e-14)
        assert out["m1"]["phase"][j] == s.phase
        assert out["m1"]["wolf_rayet"][j] == tr.is_wolf_rayet(2000.0 + times[j] / 1e6)
    assert list(out["m1"]) == [
        "radius_rsun",
        "luminosity_lsun",
        "teff_k",
        "phase",
        "age_myr",
        "wolf_rayet",
    ]


# ---------------------------------------------------------------- winds, tracks, constants
def test_sse_wind_ramp_and_reimers_from_the_gap():
    """SSE/BSE mlwind: Nieuwenhuijzen & de Jager ramped by min(1, (L - 4000) / 500), Reimers
    from kw = 2 (Hertzsprung gap) on."""
    from swingersim.massloss import sse_wind_rate

    full = -9.6e-15 * 5.0**0.81 * 4250.0**1.24 * 10.0**0.16
    assert sse_wind_rate(10.0, 4250.0, 5.0, 0.02, 1) == pytest.approx(0.5 * full)
    assert sse_wind_rate(10.0, 4000.0 + 1e-9, 5.0, 0.02, 1) == pytest.approx(0.0, abs=1e-15)
    assert sse_wind_rate(1.0, 100.0, 20.0, 0.02, 2, 0.5) == pytest.approx(-4e-10)
    assert sse_wind_rate(1.0, 100.0, 20.0, 0.02, 1, 0.5) == 0.0


def test_sse_kw_of_the_nearest_grid_track():
    """kw came from the lower-mass corner even with zero weight: a star exactly at the grid
    mass just above the kw 0/1 boundary was flagged deeply convective."""
    from swingersim.evolution import _grid, track

    g = _grid()
    iz = int(np.argmin(np.abs(g["z"] - 0.02)))
    sel = np.flatnonzero(g["phase"] == 1)
    kw0 = np.array([g["kw"][iz, i][sel][0] for i in range(len(g["mass"]))])
    i = int(np.flatnonzero(kw0 == 1)[0])  # first grid mass with kw 1
    t = track(float(g["mass"][i]), float(g["z"][iz]))
    assert t.at(0.1 * t.t_ms).sse_type == 1


def test_shared_constants():
    from swingersim import environment, magnetic, mesa, scatter

    assert environment.KMS_AUYR == pytest.approx(scatter.KMS_TO_AUYR, rel=1e-15)
    assert magnetic.MSUN_G == mesa.MSUN_G == 1.98841e33
    assert magnetic.SEC_PER_YR == units.SECONDS_PER_YEAR


# ---------------------------------------------------------------- MESA runs
def _mesa_run(tmp_path, n_prof=5, log=""):
    rdir = tmp_path / "run"
    (rdir / "LOGS").mkdir(parents=True)
    (rdir / "LOGS" / "history.data").write_text("x\n")
    for k in range(n_prof):
        (rdir / "LOGS" / f"profile{k + 1}.data").write_text("x\n")
    (rdir / "run.log").write_text(log)
    return rdir


def test_mesa_run_finished_only_after_a_normal_end(tmp_path):
    """A crashed MESA run with > 3 profiles was marked done and reduced."""
    from swingersim.mesa import finished

    assert finished(_mesa_run(tmp_path / "a"), 0)
    assert not finished(_mesa_run(tmp_path / "b"), 1)
    assert finished(_mesa_run(tmp_path / "c", log="termination code: max_age\n"), 1)
    crash = "terminated evolution: cannot find acceptable model\ntermination code: min_timestep"
    assert not finished(_mesa_run(tmp_path / "d", log=crash), 0)
    assert not finished(_mesa_run(tmp_path / "e", n_prof=3), 0)


def test_reduce_grid_skips_runs_without_done(tmp_path):
    import importlib.util
    from pathlib import Path

    spec_ = importlib.util.spec_from_file_location(
        "_t_mesa", Path(__file__).with_name("test_mesa.py")
    )
    mod = importlib.util.module_from_spec(spec_)
    spec_.loader.exec_module(mod)
    write_profile = mod.write_profile

    from swingersim.mesa import MesaGridSpec, grid_dir, reduce_grid, run_name

    spec = MesaGridSpec(name="g", masses=[0.3, 0.4], z=[0.02], threads_per_run=1)
    for m, done in ((0.3, True), (0.4, False)):
        rdir = grid_dir(spec, tmp_path) / "runs" / run_name(m, 0.02)
        (rdir / "LOGS").mkdir(parents=True)
        for k, age in enumerate((1e7, 1e8)):
            write_profile(rdir / "LOGS" / f"profile{k + 1}.data", 1.5, mass=m, radius=m, age=age)
        if done:
            (rdir / "done").write_text("exit code 0\n")
    with pytest.warns(UserWarning, match="left out"):
        out = reduce_grid(spec, tmp_path, processes=1)
    with np.load(out) as d:
        assert set(d["initial_mass"]) == {0.3}
