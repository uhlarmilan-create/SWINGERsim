"""Fixes from the independent code review: formulas checked against their derivation or the
literature, and the joining of runs across a merger or a dynamic onset."""

import math
import warnings

import numpy as np
import pytest

from swingersim import units

warnings.simplefilter("ignore")


def test_tidal_sync_time_is_the_e_folding_time():
    """t_sync = C n^2 / (2 K Omega(e)) from dOmega/dt = -2K/(C n) (Omega(e) Omega/n - N(e))
    (Leconte et al. 2010): equal to (Omega_eq - Omega) / (dOmega/dt), units of time."""
    from swingersim.tides import leconte_rates

    r = 1.0 * units.RSUN_AU
    a, e, n = 0.05, 0.0, math.sqrt(units.G * 2.0 / 0.05**3)
    lag = 100.0 / units.SECONDS_PER_YEAR
    spin = 0.5 * n
    tr = leconte_rates((1.0, 1.0), (r, r), (0.014, 0.014), (0.07, 0.07), (lag, lag), a, e,
                       (spin, spin))  # fmt: skip
    t_e = (tr.spin_pseudo_sync - spin) / tr.dspin_dt[0]
    assert tr.t_sync[0] == pytest.approx(t_e, rel=1e-9)
    assert 1e2 < tr.t_sync[0] < 1e4  # ~1.2e3 yr (was 1.5 yr, a factor n too short)


def _res(times, extremes=None, roche=None, stellar=None):
    from swingersim.backends.nbody import RunResult

    n = len(times)
    return RunResult(
        times=np.asarray(times, float),
        orbits={"A": {"e": np.zeros(n)}},
        mutual={},
        energy_error=np.zeros(n),
        hierarchy=["(1 2)"] * n,
        extremes=extremes or {},
        roche=roche or {},
        stellar=stellar or {},
    )


def test_join_thresholds_roche_and_absorbed_star():
    """The joined run: a 1 - e threshold first crossed after the merger is kept (at the time
    of the whole run), the Roche filling keeps the time of its larger maximum, and the star
    absorbed by the merger has no data after it (not those of another star)."""
    from swingersim.merger import _join

    nan = float("nan")
    a = _res(
        [0, 1, 2],
        extremes={"A": {"e_max": 0.5, "t_e_max_yr": 1.0, "q_min_au": 0.1, "t_q_min_yr": 1.0,
                        "t_1me_lt_0.01": nan}},
        roche={"m1": {"fill_max": 0.5, "t_fill_max_yr": 1.0}},
        stellar={"m1": {"radius_rsun": np.ones(3)}, "m2": {"radius_rsun": 2 * np.ones(3)},
                 "m3": {"radius_rsun": 3 * np.ones(3)}},
    )  # fmt: skip
    b = _res(
        [0, 1, 2],
        extremes={"A": {"e_max": 0.995, "t_e_max_yr": 1.5, "q_min_au": 0.05, "t_q_min_yr": 1.5,
                        "t_1me_lt_0.01": 1.5}},
        roche={"m1": {"fill_max": 0.9, "t_fill_max_yr": 1.5},
               "m2": {"fill_max": 0.3, "t_fill_max_yr": 0.5}},
        stellar={"m1": {"radius_rsun": 5 * np.ones(3)}, "m2": {"radius_rsun": 7 * np.ones(3)}},
    )  # fmt: skip
    info = {"names": {"m1": "m1", "m2": "m3"}, "kind": "merger", "hierarchy": "(1 3)",
            "stars": ["m1", "m2"], "mass_msun": 1.9, "label": "", "time_yr": 2.0}  # fmt: skip
    j = _join(a, b, info, 2.0)
    ex = j.extremes["A"]
    assert ex["t_1me_lt_0.01"] == pytest.approx(3.5)
    assert ex["e_max"] == 0.995 and ex["t_e_max_yr"] == pytest.approx(3.5)
    assert j.roche["m1"] == {"fill_max": 0.9, "t_fill_max_yr": pytest.approx(3.5)}
    assert j.roche["m3"]["t_fill_max_yr"] == pytest.approx(2.5)
    r2 = j.stellar["m2"]["radius_rsun"]  # absorbed into m1 at t = 2
    assert np.all(r2[:3] == 2.0) and np.all(np.isnan(r2[3:]))
    assert np.all(j.stellar["m3"]["radius_rsun"][3:] == 7.0)  # old m3 = new m2
    assert np.all(j.stellar["m1"]["radius_rsun"][3:] == 5.0)


def test_roemer_delay_across_the_binary():
    """Kaplan (2010): eclipses of the lighter star (farther from the centre of mass) are late by
    z (m_heavy - m_light) / (M c); equal masses: no shift."""
    from swingersim.timing import roemer_in_binary

    assert roemer_in_binary(0.1, 1.0, 1.0) == 0.0
    dt = roemer_in_binary(0.1, 2.0, 1.0) * units.SECONDS_PER_YEAR
    assert dt == pytest.approx(0.1 * units.AU_KM / 3.0 / 299792.458, rel=1e-6)


def test_magnetic_braking_backward_is_the_inverse():
    from swingersim import magnetic

    for w0 in (3e-5, 1e-5, 4e-6):  # saturated, crossing Omega_sat, unsaturated
        for dt in (1e6, 1e8, 3e9):
            w = magnetic.spin_down(w0, 1.0, 1.0, 5772.0, 0.07, dt)
            assert w < w0
            assert magnetic.spin_down(w, 1.0, 1.0, 5772.0, 0.07, -dt) == pytest.approx(w0, 1e-9)


def test_first_crossing_between_samples():
    """A spike shorter than the sampling crosses the level: found by refining the local
    maxima, at the right time (secular overflow detection)."""
    from swingersim.backends.secular_common import first_crossing

    def fun(x):
        return 0.5 + 0.6 * math.exp(-(((x - 3.03) / 0.01) ** 2))

    t = np.linspace(0, 10, 101)  # spacing 0.1, the spike's peak 0.03 from a sample
    f = np.array([fun(x) for x in t])
    assert f.max() < 1.0
    tc = first_crossing(t, f, fun, 1.0, True)
    assert fun(tc) == pytest.approx(1.0, abs=1e-6) and tc < 3.03
    assert math.isnan(first_crossing(t, f, fun, 1.2, True))


def test_regressing_apsides_have_a_negative_period():
    from swingersim.diagnostics.binary import measured_apsidal_period

    t = np.linspace(0, 10, 50)
    w = -0.5 * t
    assert measured_apsidal_period(t, 0 * t, w) == pytest.approx(-2 * math.pi / 0.5)
    assert measured_apsidal_period(t, 0 * t, -w) == pytest.approx(2 * math.pi / 0.5)


def test_circular_orbit_is_seeded_in_the_secular_backend():
    """e = 0 is a fixed point of the secular equations: a circular inner orbit at 80 deg
    must still develop Kozai-Lidov cycles (as in an N-body run)."""
    import yaml

    from swingersim.backends.secular import run_secular
    from swingersim.config import Config

    d = yaml.safe_load(
        """
topology: 2+2
masses: {m1: 1.0, m2: 1.0, m3: 1.0, m4: 1.0}
orbits:
  A: {a: 1.0, e: 0.0}
  B: {a: 1.0, e: 0.1, i: 10}
  out: {a: 15.0, e: 0.1, i: 80}
run: {t_max_yr: 60000, n_output: 601, backend: secular}
"""
    )
    res = run_secular(Config.model_validate(d))
    assert res.orbits["A"]["e"].max() > 0.5


def test_tides_without_rotation_keep_the_spins_after_an_onset():
    """physics.tides without physics.rotation: the continued run switches rotation on (else
    the spins set for it would be ignored); stars without a set spin do not rotate."""
    import yaml

    from swingersim.config import Config
    from swingersim.hierarchy import build_tree, initialize
    from swingersim.merger import NON_ROTATING_D, post_onset_config

    d = yaml.safe_load(
        """
topology: custom
hierarchy: {out: [A, m3], A: [m1, m2]}
masses: {m1: 1.0, m2: 0.8, m3: 1.0}
orbits:
  A: {a: 0.1, e: 0.9}
  out: {a: 1.5, e: 0.2, i: 87.0}
stellar_model: {kind: sse, age_myr: 1000, z: 0.02}
physics: {tides: true, evolution: true, mass_transfer: {dynamic_onset: angular_momentum}}
run: {t_max_yr: 300, n_output: 61, backend: secular}
"""
    )
    cfg = Config.model_validate(d)
    root, m = build_tree(cfg)
    pos, vel = initialize(root, m)
    new, info = post_onset_config(
        cfg, {"pos": pos, "vel": vel, "masses": m, "spins": None}, "A", 0, 1, 40.0, 300.0
    )
    assert new is not None, info
    assert new.physics.rotation
    assert new.stars["m1"].P_rot == "sync" and new.stars["m2"].P_rot == "sync"
    assert new.stars["m3"].P_rot == NON_ROTATING_D


def test_mist_giant_branch_starts_on_the_hayashi_line():
    """The Hertzsprung gap of a MIST track ends where the star reaches the Hayashi line (not at
    the middle of the TAMS - RGB tip EEPs): a 1 Msun red giant of 6 Rsun is on the giant
    branch, a 25 Msun star crossing the gap at 21 000 K is not."""
    pytest.importorskip("swingersim.mist")
    from swingersim.masstransfer import q_crit
    from swingersim.mist import EEP_RGB_TIP, MistTrack, available

    if not available():
        pytest.skip("no MIST grid installed")
    t1 = MistTrack(1.0)
    st = t1.at_eep(t1.eep_gb + 1.0)
    assert st.phase == 3 and 1.5 < st.radius_rsun < 3.0 and st.teff_k < 5500
    k = next(e for e in range(454, 605) if t1.at_eep(e).radius_rsun > 6.0)
    st = t1.at_eep(k)
    assert st.phase == 3
    assert q_crit(st.phase, 1.0, st.m_core) < 1.0  # BSE eq. 57: dynamically unstable at q 2.5
    t25 = MistTrack(25.0)
    assert t25.eep_gb == EEP_RGB_TIP
    assert t25.at_eep(560).phase == 2


def test_common_envelope_core_radii():
    """Hurley et al. (2000), Sect. 6.1: a core-helium-burning core is a naked helium star (not
    degenerate after the helium flash), as the core of a massive star."""
    from swingersim.commonenv import core_radius, r_helium_zams, r_white_dwarf

    assert core_radius(0.3, 1.0, 3) == r_white_dwarf(0.3)
    assert core_radius(0.47, 1.5, 4) == pytest.approx(r_helium_zams(0.47))
    assert core_radius(0.47, 1.5, 4) == pytest.approx(0.10, abs=0.01)
    assert core_radius(4.0, 15.0, 5) == pytest.approx(0.53, abs=0.01)
    assert core_radius(0.6, 2.0, 6) == r_white_dwarf(0.6)


def test_mcmc_is_reproducible_with_a_seed():
    """The emcee stretch moves draw from the sampler's own generator: seeded with the fit's
    seed, two chains are identical."""
    from swingersim import dynfit

    class _LinearProblem:  # y = a + b t (as tests/test_v036_fit.py)
        def __init__(self):
            self.t = np.linspace(-1, 1, 30)
            self.y = 0.3 + 2.0 * self.t + 0.1 * np.random.default_rng(4).normal(size=30)

        def n_data(self):
            return 30

        def models(self, X, pool, store=True):
            return list(np.atleast_2d(X))

        def resid(self, m, x):
            return (self.y - (m[0] + m[1] * self.t)) / 0.1

    out = []
    for _ in range(2):
        pb = _LinearProblem()
        x = np.array([0.3, 2.0])
        cov = np.diag([0.02**2, 0.03**2])
        mc = dynfit.run_mcmc(pb, x, np.sqrt(np.diag(cov)), cov, np.array([-5.0, -5.0]),
                             np.array([5.0, 5.0]), 1.0, 8, 40, 11, None, 1,
                             lambda *a: None)  # fmt: skip
        out.append(mc["log_prob"])
    assert np.array_equal(out[0], out[1])
