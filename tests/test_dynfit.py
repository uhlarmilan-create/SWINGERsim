"""Stage 2 of swingersim fit: the joint N-body model of the times of minima (dynfit.py)."""

import copy
import math
from pathlib import Path

import numpy as np
import pytest
import yaml

from swingersim import dynfit, ocfit
from swingersim.config import Config, load_config
from swingersim.timing import eclipse_times

ROOT = Path(__file__).resolve().parents[1]
BU_CMI = ROOT / "configs" / "squadra_test" / "bu_cmi_anchored.yaml"
T0 = 2458491.5


def _eclipses(cfg, years_back, years_fwd, spp=40):
    ecl = eclipse_times(cfg, -years_back, samples_per_period=spp)[0]
    ecl += eclipse_times(cfg, years_fwd, samples_per_period=spp)[0]
    return ecl


def _observed(ecl, sigma=1e-4):
    return dynfit.Observed(
        np.array([e.binary for e in ecl]),
        np.array([0 if e.kind == "primary" else 1 for e in ecl]),
        np.array([T0 + e.t_obs_yr * 365.25 for e in ecl]),
        np.full(len(ecl), sigma),
        np.array(["synthetic"] * len(ecl)),
    )


def test_model_times_agree_with_eclipse_times_forward_and_backward():
    cfg = load_config(BU_CMI)
    ecl = _eclipses(cfg, 0.3, 0.3, spp=400)[::3]
    obs = _observed(ecl)
    truth = obs.time.copy()
    obs.time = obs.time + 0.004 * np.sin(np.arange(len(obs)))  # measured != model: the nearest
    tm = dynfit.model_times(cfg, obs, T0)
    assert np.all(np.isfinite(tm))
    assert np.max(np.abs(tm - truth)) * 86400 < 0.5  # seconds
    assert np.any(obs.pair == "A") and np.any(obs.pair == "B")
    assert np.any(obs.time < T0) and np.any(obs.time > T0)


def test_ecosw_esinw_round_trip_and_unknown_names():
    raw = yaml.safe_load(BU_CMI.read_text())
    names = ["orbits.A.ecosw", "orbits.A.esinw", "orbits.out.P", "masses.m1"]
    x = dynfit.start_values(raw, names)
    e, w = raw["orbits"]["A"]["e"], math.radians(raw["orbits"]["A"]["omega"])
    assert x[0] == pytest.approx(e * math.cos(w)) and x[1] == pytest.approx(e * math.sin(w))
    d = dynfit.apply(raw, names, x + np.array([0.0, 0.0, 1.0, 0.1]))
    assert d["orbits"]["A"]["e"] == pytest.approx(e)
    assert d["orbits"]["A"]["omega"] == pytest.approx(raw["orbits"]["A"]["omega"])
    assert d["orbits"]["out"]["P"] == pytest.approx(raw["orbits"]["out"]["P"] + 1.0)
    assert d["masses"]["m1"] == pytest.approx(raw["masses"]["m1"] + 0.1)
    assert raw["orbits"]["out"]["P"] != d["orbits"]["out"]["P"]  # the input is not changed
    with pytest.raises(ValueError):
        dynfit.start_values(raw, ["orbits.A.a"])


def test_read_minima_accepts_the_gui_columns(tmp_path):
    p = tmp_path / "m.csv"
    p.write_text(
        "# minima\nbinary,kind,bjd_tdb,sigma_d,cycle,ref\n"
        "A,primary,2459000.1,0.001,0,TESS\nB,secondary,2459001.6,0.002,,ZTF\n"
    )
    m = ocfit.read_minima(p)
    assert list(m.pair) == ["A", "B"] and list(m.kind) == [0, 1]
    assert list(m.source) == ["TESS", "ZTF"] and m.sigma[1] == 0.002


def test_stage2_recovers_perturbed_elements(tmp_path):
    """Closure test: synthetic minima of the anchored BU CMi, fitted from perturbed values."""
    raw = yaml.safe_load(BU_CMI.read_text())
    cfg = Config.model_validate(raw)
    ecl = _eclipses(cfg, 0.5, 0.5)[::2]
    rng = np.random.default_rng(3)
    with open(tmp_path / "minima.csv", "w") as f:
        f.write("pair,type,time,sigma,source\n")
        for e in ecl:
            t = T0 + e.t_obs_yr * 365.25 + rng.normal(0.0, 1e-4)
            f.write(f"{e.binary},{e.kind[0]},{t:.7f},0.0001,synthetic\n")
    start = yaml.safe_load(BU_CMI.read_text())
    start["orbits"]["A"]["P"] += 2e-5
    start["orbits"]["A"]["M0"] += 0.05
    start["orbits"]["out"]["P"] += 0.3
    start["orbits"]["out"]["M0"] -= 0.5
    free = ["orbits.A.P", "orbits.A.M0", "orbits.out.P", "orbits.out.M0"]
    start["fit"] = {
        "minima": "minima.csv",
        "stage2": {"bjd_t0": T0, "free": free, "processes": 1, "max_iterations": 12},
        "output": {"dir": "out"},
    }
    p = tmp_path / "start.yaml"
    p.write_text(yaml.safe_dump(start, sort_keys=False))
    res = dynfit.run_stage2(p, log=lambda *a: None)
    truth = dynfit.start_values(raw, free)
    pull = (res.x - truth) / res.err
    assert np.all(np.abs(pull) < 4.0), pull
    assert res.chi2_red < 2.0
    assert max(res.rms_min(res.t_model).values()) < 0.5  # minutes
    for f in ("stage2_params.csv", "stage2_residuals.csv", "stage2_fitted.yaml", "stage2_oc.png"):
        assert (tmp_path / "out" / f).exists()
    fitted = load_config(tmp_path / "out" / "stage2_fitted.yaml")
    assert fitted.fit is None
    assert fitted.orbits["A"].P == pytest.approx(res.x[0])
    assert "stage 2" in dynfit.report(res)


def test_stage2_spec_is_validated():
    raw = yaml.safe_load(BU_CMI.read_text())
    raw["fit"] = {"minima": "x.csv", "stage2": {"bjd_t0": T0, "free": []}}
    with pytest.raises(ValueError):
        Config.model_validate(raw)
    raw["fit"]["stage2"] = {"bjd_t0": T0, "free": ["orbits.A.P"], "sampler": "emcee"}
    assert Config.model_validate(raw).fit.stage2.mcmc_steps == 300


def test_missing_emcee_fails_before_the_fit(monkeypatch, tmp_path):
    import sys

    raw = yaml.safe_load(BU_CMI.read_text())
    raw["fit"] = {
        "minima": str(ROOT / "data" / "observed" / "bu_cmi" / "minima_pribulla2023.csv"),
        "stage2": {"bjd_t0": T0, "free": ["orbits.A.P"], "sampler": "emcee"},
    }
    p = tmp_path / "c.yaml"
    p.write_text(yaml.safe_dump(raw))
    monkeypatch.setitem(sys.modules, "emcee", None)  # import emcee -> ImportError
    with pytest.raises(ValueError, match="pip install emcee"):
        dynfit.run_stage2(p, log=lambda *a: None)


# ------------------------------------------------------------------ errors and outliers (v0.31)
def test_source_groups_pool_small_sources():
    g = dynfit.source_groups(["TESS"] * 6 + ["RU"] * 5 + ["X", "Y"], 5)
    assert list(g[:11]) == ["TESS"] * 6 + ["RU"] * 5
    assert list(g[11:]) == [dynfit.OTHER] * 2


def test_error_factors_scale_and_jitter():
    rng = np.random.default_rng(1)
    n = 4000
    groups = np.array(["good"] * n + ["bad"] * n)
    sigma = np.full(2 * n, 1e-4)
    resid = np.concatenate([rng.normal(0, 0.8e-4, n), rng.normal(0, 5e-4, n)])
    f = dynfit.error_factors(resid, sigma, groups, "per_source", "scale", 0)
    assert f["good"] == 1.0  # never below the quoted errors
    assert f["bad"] == pytest.approx(5.0, rel=0.05)
    j = dynfit.error_factors(resid, sigma, groups, "per_source", "jitter", 0)
    assert j["good"] == 0.0
    assert j["bad"] == pytest.approx(math.sqrt(25 - 1) * 1e-4, rel=0.05)
    s = dynfit.apply_errors(sigma, groups, "per_source", "jitter", j)
    chi2 = np.sum((resid[groups == "bad"] / s[groups == "bad"]) ** 2) / n
    assert chi2 == pytest.approx(1.0, rel=1e-6)
    g = dynfit.error_factors(resid, sigma, groups, "global", "scale", 0)
    assert list(g) == ["all"] and g["all"] == pytest.approx(math.sqrt((0.64 + 25) / 2), rel=0.05)


def test_clip_outliers_robust_and_capped():
    rng = np.random.default_rng(2)
    z = rng.normal(0, 1, 200)
    z[[3, 50, 120]] = [40.0, -25.0, np.nan]
    g = np.array(["a"] * 200)
    keep, ratio = dynfit.clip_outliers(z, g, np.ones(200, bool), 5.0, 0.1)
    assert sorted(np.where(~keep)[0]) == [3, 50, 120] and ratio[3] > 30
    keep2, _ = dynfit.clip_outliers(z, g, np.ones(200, bool), 5.0, 0.01)  # room for 2
    assert int(np.sum(~keep2)) == 2 and not keep2[120] and not keep2[3]
    # a source that is ten times worse than quoted keeps its minima (robust scale of the group)
    zb = np.concatenate([z[:100], 10 * rng.normal(0, 1, 100)])
    gb = np.array(["a"] * 100 + ["b"] * 100)
    keep3, _ = dynfit.clip_outliers(zb, gb, np.ones(200, bool), 5.0, 0.5)
    assert np.all(keep3[100:])


def test_stage2_per_source_errors_and_outliers(tmp_path):
    """Synthetic minima of the anchored BU CMi from two sources: 'good' with honest errors and
    'survey' with errors quoted five times too small, plus three gross outliers."""
    raw = yaml.safe_load(BU_CMI.read_text())
    cfg = Config.model_validate(raw)
    ecl = _eclipses(cfg, 0.5, 0.5)[::3]
    rng = np.random.default_rng(5)
    bad = {5, 17, 29}
    with open(tmp_path / "minima.csv", "w") as f:
        f.write("pair,type,time,sigma,source\n")
        for k, e in enumerate(ecl):
            src, true = ("good", 1e-4) if k % 2 == 0 else ("survey", 5e-4)
            t = T0 + e.t_obs_yr * 365.25 + rng.normal(0.0, true)
            if k in bad:
                t += 0.03  # 43 min
            f.write(f"{e.binary},{e.kind[0]},{t:.7f},0.0001,{src}\n")
    start = yaml.safe_load(BU_CMI.read_text())
    start["orbits"]["A"]["P"] += 1e-5
    start["orbits"]["out"]["M0"] -= 0.3
    free = ["orbits.A.P", "orbits.A.M0", "orbits.out.P", "orbits.out.M0"]
    start["fit"] = {
        "minima": "minima.csv",
        "stage2": {
            "bjd_t0": T0,
            "free": free,
            "processes": 1,
            "max_iterations": 10,
            "error_scale": "per_source",
            "clip_sigma": 5.0,
        },
        "output": {"dir": "out"},
    }
    p = tmp_path / "start.yaml"
    p.write_text(yaml.safe_dump(start, sort_keys=False))
    res = dynfit.run_stage2(p, log=lambda *a: None)
    assert res.outliers is not None
    rows = (tmp_path / "minima.csv").read_text().splitlines()[1:]
    injected = sorted(float(rows[k].split(",")[2]) for k in bad)
    assert sorted(res.outliers.time.tolist()) == pytest.approx(injected, abs=1e-6)
    assert res.error_factors["good"] == pytest.approx(1.0, abs=0.25)
    assert 3.5 < res.error_factors["survey"] < 6.5
    assert 0.6 < res.chi2_red < 1.5
    assert len(res.passes) >= 2
    truth = dynfit.start_values(raw, free)
    assert np.all(np.abs((res.x - truth) / res.err) < 4.0)
    rep = dynfit.report(res)
    assert "outlying minima left out" in rep and "survey x" in rep
    txt = (tmp_path / "out" / "stage2_residuals.csv").read_text()
    assert txt.count(",outlier") == 3 and "sigma_quoted_d" in txt
    assert (tmp_path / "out" / "stage2_oc.png").exists()


def test_stage2_error_scale_inherits_fit_error_scale():
    raw = yaml.safe_load(BU_CMI.read_text())
    raw["fit"] = {
        "minima": "x.csv",
        "error_scale": "per_source",
        "stage2": {"bjd_t0": T0, "free": ["orbits.A.P"]},
    }
    c = Config.model_validate(raw)
    assert c.fit.stage2.error_scale is None and c.fit.error_scale == "per_source"
    raw["fit"]["stage2"]["clip_sigma"] = 0
    with pytest.raises(ValueError):
        Config.model_validate(raw)


# ------------------------------------------------------------------ priors (v0.31)
def test_priors_are_validated_and_accept_the_ensemble_form():
    raw = yaml.safe_load(BU_CMI.read_text())
    st = {"bjd_t0": T0, "free": ["orbits.A.P", "orbits.out.ecosw", "orbits.out.esinw"]}
    raw["fit"] = {"minima": "x.csv", "stage2": dict(st, priors={"orbits.A.P": [2.9, 0.1]})}
    assert Config.model_validate(raw).fit.stage2.priors["orbits.A.P"] == (2.9, 0.1)
    raw["fit"]["stage2"] = dict(st, priors={"orbits.out.e": {"normal": [0.27, 0.004]}})
    assert Config.model_validate(raw).fit.stage2.priors["orbits.out.e"] == (0.27, 0.004)
    for bad in (
        {"orbits.A.e": [0.2, 0.01]},  # fixed parameter
        {"orbits.A.P": [2.9, 0.0]},  # sigma
        {"orbits.A.P": {"uniform": [2.8, 3.0]}},
    ):
        raw["fit"]["stage2"] = dict(st, priors=bad)
        with pytest.raises(ValueError):
            Config.model_validate(raw)


def test_prior_residuals_wrap_angles_and_follow_ecosw():
    raw = yaml.safe_load(BU_CMI.read_text())
    raw["orbits"]["out"]["omega"] = 359.0
    names = ["orbits.out.ecosw", "orbits.out.esinw", "orbits.out.M0"]
    x = dynfit.start_values(raw, names)
    e = raw["orbits"]["out"]["e"]
    pri = {
        "orbits.out.omega": (1.0, 1.0),
        "orbits.out.e": (e + 0.01, 0.01),
        "orbits.out.M0": (raw["orbits"]["out"]["M0"] - 350.0, 5.0),
    }
    r = dynfit.prior_residuals(raw, names, x, pri)
    assert r[0] == pytest.approx(-2.0, abs=1e-6)  # 359 against 1 deg: -2 deg, not +358
    assert r[1] == pytest.approx(-1.0, abs=1e-6)
    assert r[2] == pytest.approx(-2.0, abs=1e-6)  # +350 deg = -10 deg


def _short_fit(tmp_path, priors, years=0.25):
    raw = yaml.safe_load(BU_CMI.read_text())
    cfg = Config.model_validate(raw)
    ecl = _eclipses(cfg, years, years)[::2]
    rng = np.random.default_rng(11)
    d = tmp_path / ("p" if priors else "n")
    d.mkdir()
    with open(d / "minima.csv", "w") as f:
        f.write("pair,type,time,sigma,source\n")
        for e in ecl:
            t = T0 + e.t_obs_yr * 365.25 + rng.normal(0.0, 3e-4)
            f.write(f"{e.binary},{e.kind[0]},{t:.7f},0.0003,synthetic\n")
    start = yaml.safe_load(BU_CMI.read_text())
    free = ["orbits.out.P", "orbits.out.M0", "orbits.out.e"]
    start["fit"] = {
        "minima": "minima.csv",
        "stage2": {
            "bjd_t0": T0,
            "free": free,
            "processes": 1,
            "max_iterations": 8,
            "priors": priors,
        },
        "output": {"dir": "out"},
    }
    p = d / "start.yaml"
    p.write_text(yaml.safe_dump(start, sort_keys=False))
    return dynfit.run_stage2(p, log=lambda *a: None), dynfit.start_values(raw, free), d


def test_stage2_prior_tightens_and_pulls(tmp_path):
    free_res, truth, _ = _short_fit(tmp_path, {})
    e_true = float(truth[2])
    pri = {"orbits.out.e": [e_true + 0.0005, 0.0001]}
    res, _, d = _short_fit(tmp_path, pri)
    # the prior is about six times tighter than the data alone; the fit follows it
    assert free_res.err[2] > 3e-4
    assert res.err[2] < 1.2e-4
    assert res.x[2] == pytest.approx(e_true + 0.0005, abs=1.5e-4)
    assert res.chi2_priors > 0 and len(res.priors) == 1
    assert "priors (mean +- sigma" in dynfit.report(res)
    rows = (d / "out" / "stage2_params.csv").read_text().splitlines()
    assert ",prior_mean,prior_sigma," in rows[0] and "0.0001" in rows[3]


# ------------------------------------------------------------------ anchoring and predictions (v0.31)
def test_data_ephemeris_bridges_gaps():
    P, T0 = 2.9384888, 2458000.1234
    E = np.concatenate([np.arange(0, 40), np.arange(900, 940), np.arange(2400, 2410)])
    kind = np.tile([0, 1], len(E))[: len(E)]
    t = T0 + P * (E + 0.5 * kind) + 0.003 * np.sin(E / 50.0)  # a curved O - C of 4 min
    _, Pf, Ef = dynfit.data_ephemeris(t, kind, P * (1 + 3e-4), T0 + 905 * P)
    assert np.all(Ef - Ef[0] == (E + 0.5 * kind) - (E[0] + 0.5 * kind[0]))
    pri = kind == 0
    assert Pf == pytest.approx(np.polyfit(E[pri], t[pri], 1)[0], rel=1e-9)


def test_anchor_recovers_period_and_phase(tmp_path):
    raw = yaml.safe_load(BU_CMI.read_text())
    cfg = Config.model_validate(raw)
    ecl = _eclipses(cfg, 1.0, 1.0)[::3]
    obs = _observed(ecl, sigma=2e-4)
    start = copy.deepcopy(raw)
    start["orbits"]["A"]["P"] *= 1 - 3.3e-4  # as the published BU CMi against swingersim
    start["orbits"]["A"]["M0"] += 1.5
    start["orbits"]["B"]["P"] *= 1 - 4.0e-4
    start["orbits"]["B"]["M0"] -= 2.0
    a = dynfit.anchor(start, None, obs, T0, log=lambda *x: None)
    for p in ("A", "B"):
        assert a.raw["orbits"][p]["P"] == pytest.approx(raw["orbits"][p]["P"], rel=3e-7)
        dm = (a.raw["orbits"][p]["M0"] - raw["orbits"][p]["M0"] + 180) % 360 - 180
        assert abs(dm) < 0.02
        assert a.rms_min[p] < 0.3
        assert abs(a.pairs[p]["drift_min"]) > 60  # the start drifted by more than an hour
    p = dynfit.write_anchored(a, tmp_path, "start.yaml", T0)
    assert load_config(p).orbits["A"].P == pytest.approx(a.raw["orbits"]["A"]["P"])
    assert "anchored" in (tmp_path / "anchor_report.txt").read_text()


def test_predictions_spread_grows_with_the_period_error(tmp_path):
    raw = yaml.safe_load(BU_CMI.read_text())
    cfg = Config.model_validate(raw)
    ecl = _eclipses(cfg, 0.0, 0.1)
    obs = _observed(ecl)
    names = ["orbits.A.P", "orbits.A.M0"]
    x = dynfit.start_values(raw, names)
    sP, sM = 2e-6, 0.01
    res = dynfit.Stage2Result(
        names,
        x,
        x,
        np.array([sP, sM]),
        np.diag([sP**2, sM**2]),
        obs,
        obs.time,
        obs.time,
        0,
        0,
        0,
        0,
        0,
        raw,
        T0,
    )
    pr = dynfit.predict_minima(
        res,
        None,
        years=0.3,
        pairs=["A"],
        start=float(obs.time.max()),
        samples=60,
        processes=1,
        log=lambda *a: None,
    )
    assert set(pr.pair.tolist()) == {"A"} and np.all(pr.t >= obs.time.max())
    n = (pr.t - T0) / x[0]
    expect = np.sqrt((sM / 360 * x[0]) ** 2 + (sP * n) ** 2)
    assert np.all(np.isfinite(pr.sigma))
    assert np.median(pr.sigma / expect) == pytest.approx(1.0, abs=0.3)
    assert pr.sigma[-1] > pr.sigma[0]
    assert np.all(pr.duration > 0) and np.all(pr.impact < 1)
    from swingersim.config import CoordinatesSpec

    path = dynfit.write_predictions(
        pr, tmp_path / "p.csv", CoordinatesSpec(ra_deg=119.524529, dec_deg=7.213475)
    )
    rows = path.read_text().splitlines()
    assert rows[0].startswith("pair,type,bjd_tdb,sigma_min") and len(rows) == len(pr.t) + 1
    utc = rows[1].split(",")[10]
    assert utc[:2] == "20" and pr.jd_utc is not None
    assert abs(pr.t[0] - pr.jd_utc[0]) * 1440 < 9.5  # barycentric correction + TT - UTC
    dynfit.plot_predictions(res, pr, tmp_path / "p.png")
    assert (tmp_path / "p.png").stat().st_size > 10000


# ------------------------------------------------------------------ radial velocities (v0.31)
def test_model_rv_light_time_and_sign():
    """RV = +v_z (away from the observer), the time derivative of the light time z / c, and the
    light time of the star itself is taken into account."""
    from swingersim import radvel

    cfg = load_config(BU_CMI)
    t = T0 + np.linspace(1.0, 3.0, 9)
    v = radvel.model_rv(cfg, np.full(9, "m1"), t, T0)
    h = 1e-4
    zp = radvel.model_rv(cfg, np.full(9, "m1"), t + h, T0)
    assert np.all(np.isfinite(v)) and np.max(np.abs(v)) > 50
    assert np.max(np.abs(zp - v)) < 0.1  # smooth
    # m1 and m2 of pair A move in opposite directions (momentum)
    v2 = radvel.model_rv(cfg, np.full(9, "m2"), t, T0)
    assert np.corrcoef(v, v2)[0, 1] < -0.9


def test_stage2_fits_radial_velocities(tmp_path):
    """Closure test: synthetic minima and velocities of the anchored BU CMi; the masses of
    pair A and gamma come back from perturbed start values; an offset of one source too."""
    from swingersim import radvel

    raw = yaml.safe_load(BU_CMI.read_text())
    cfg = Config.model_validate(raw)
    ecl = _eclipses(cfg, 0.2, 0.2)[::3]
    rng = np.random.default_rng(21)
    with open(tmp_path / "minima.csv", "w") as f:
        f.write("pair,type,time,sigma,source\n")
        for e in ecl:
            t = T0 + e.t_obs_yr * 365.25 + rng.normal(0.0, 2e-4)
            f.write(f"{e.binary},{e.kind[0]},{t:.7f},0.0002,synthetic\n")
    tt = np.sort(T0 + rng.uniform(-60.0, 60.0, 30))
    stars = np.repeat(["m1", "m2", "m3", "m4"], len(tt))
    times = np.tile(tt, 4)
    v = radvel.model_rv(cfg, stars, times, T0) + 23.0
    src = np.where(np.arange(len(times)) % 3 == 0, "B", "A")
    v = v + np.where(src == "B", 2.5, 0.0) + rng.normal(0.0, 1.0, len(v))
    with open(tmp_path / "rv.csv", "w") as f:
        f.write("# synthetic\nstar,time,rv,sigma,source\n")
        f.writelines(f"{stars[k]},{times[k]:.5f},{v[k]:.3f},1.0,{src[k]}\n" for k in range(len(v)))
    start = yaml.safe_load(BU_CMI.read_text())
    start["masses"]["m1"] *= 1.03
    start["masses"]["m2"] *= 0.97
    free = ["masses.m1", "masses.m2", "rv.gamma", "rv.offset.B"]
    start["fit"] = {
        "minima": "minima.csv",
        "rv": {"file": "rv.csv", "gamma": 20.0},
        "stage2": {"bjd_t0": T0, "free": free, "processes": 1, "max_iterations": 10},
        "output": {"dir": "out"},
    }
    p = tmp_path / "start.yaml"
    p.write_text(yaml.safe_dump(start, sort_keys=False))
    res = dynfit.run_stage2(p, log=lambda *a: None)
    truth = np.array([raw["masses"]["m1"], raw["masses"]["m2"], 23.0, 2.5])
    pull = (res.x - truth) / res.err
    assert np.all(np.abs(pull) < 4.0), (res.x, res.err)
    assert res.chi2_red < 1.6
    assert max(res.rms_rv(res.rv_model).values()) < 1.5
    assert "radial velocities" in dynfit.report(res)
    for f in ("stage2_rv.csv", "stage2_rv.png", "stage2_params.csv"):
        assert (tmp_path / "out" / f).exists()
    assert "gamma" in (tmp_path / "out" / "stage2_fitted.yaml").read_text().splitlines()[4]


# ------------------------------------------------------------------ light curves (v0.31)
def test_stage2_fits_a_light_curve(tmp_path):
    """Closure test: two segments of the synthetic TESS light curve of the anchored BU CMi
    (around eclipses of pair A, diluted, with a trend and noise) and the minima; the
    inclination of pair A and the radius of m1 come back from perturbed start values."""
    from swingersim import lcfit
    from swingersim.config import LCSpec
    from swingersim.photometry import light_curve

    raw = yaml.safe_load(BU_CMI.read_text())
    cfg = Config.model_validate(raw)
    ecl = [e for e in _eclipses(cfg, 0.0, 0.1) if e.binary == "A"]
    rng = np.random.default_rng(8)
    tt = []
    for e in ecl[2:4]:
        tm = T0 + e.t_obs_yr * 365.25
        tt.append(np.linspace(tm - 0.2, tm + 0.2, 70))
    tt = np.concatenate(tt)
    f = light_curve(cfg, tt - T0, "TESS").flux
    x = tt - tt.mean()
    flux = 0.7 * f + 0.3 + 2e-4 * x + rng.normal(0, 5e-4, len(tt))
    with open(tmp_path / "lc.csv", "w") as fh:
        fh.write("time,flux,sigma\n")
        fh.writelines(f"{tt[k]:.6f},{flux[k]:.6f},0.0005\n" for k in range(len(tt)))
    with open(tmp_path / "minima.csv", "w") as fh:
        fh.write("pair,type,time,sigma,source\n")
        for e in _eclipses(cfg, 0.0, 0.2)[::2]:
            fh.write(f"{e.binary},{e.kind[0]},{T0 + e.t_obs_yr * 365.25:.7f},0.0002,syn\n")
    start = yaml.safe_load(BU_CMI.read_text())
    start["orbits"]["A"]["i"] += 0.6
    start["stars"]["m1"]["R"] *= 1.02
    free = ["orbits.A.i", "stars.m1.R"]
    start["fit"] = {
        "minima": "minima.csv",
        "lc": [{"file": "lc.csv", "label": "tess"}],
        "stage2": {"bjd_t0": T0, "free": free, "processes": 1, "max_iterations": 8},
        "output": {"dir": "out"},
    }
    p = tmp_path / "start.yaml"
    p.write_text(yaml.safe_dump(start, sort_keys=False))
    res = dynfit.run_stage2(p, log=lambda *a: None)
    truth = dynfit.start_values(raw, free)
    pull = (res.x - truth) / res.err
    assert np.all(np.abs(pull) < 4.0), (res.x, truth, res.err)
    rms = res.rms_lc()
    assert max(rms.values()) < 0.7  # ppt, noise 0.5
    a = [c[0] for c in res.lc_coef[0].values()]
    assert np.allclose(a, 0.7, atol=0.03)  # the scale is partly degenerate with i and R
    assert "light curves" in dynfit.report(res)
    for fn in ("stage2_lc.csv", "stage2_lc.png"):
        assert (tmp_path / "out" / fn).exists()
    # the calibration alone: exact for a model that is a scaled copy of the data
    lc = lcfit.ObservedLC(tt, flux, np.full(len(tt), 5e-4), np.array(["x"] * len(tt)), "t")
    _, coef = lcfit.calibrate(lc, f)
    assert coef["x"][0] == pytest.approx(0.7, abs=0.01)
    assert LCSpec(file="a.csv").baseline_degree == 1


def test_pseudosync_and_apsidal_table():
    from swingersim import apsidal
    from swingersim.tides import pseudo_sync_factor

    assert pseudo_sync_factor(0.0) == pytest.approx(1.0)
    assert pseudo_sync_factor(0.2) == pytest.approx(1.2419, abs=1e-4)
    raw = yaml.safe_load(BU_CMI.read_text())
    rows = apsidal.apsidal_table(Config.model_validate(raw), {"A": 28.7})
    a = next(r for r in rows if r["pair"] == "A")
    assert a["rot_pseudosync_deg_yr"] > a["rot_sync_deg_yr"] > 0
    tot = a["tidal_deg_yr"] + a["rot_configured_deg_yr"] + a["gr_deg_yr"]
    assert a["classical_gr_deg_yr"] == pytest.approx(tot)
    assert a["k2_factor_for_target"] > 0
    raw["stars"]["m1"]["P_rot"] = "pseudosync"
    raw["stars"]["m2"]["P_rot"] = "pseudosync"
    b = next(r for r in apsidal.apsidal_table(Config.model_validate(raw)) if r["pair"] == "A")
    assert b["rot_configured_deg_yr"] == pytest.approx(b["rot_pseudosync_deg_yr"], rel=1e-6)
    assert b["P_rot_1_d"] == pytest.approx(b["P_d"] / b["pseudosync_factor"], rel=1e-6)


def test_free_stellar_parameters_round_trip():
    raw = yaml.safe_load(BU_CMI.read_text())
    names = ["stars.m1.R", "stars.m2.apsidal_k2"]
    x = dynfit.start_values(raw, names)
    assert x[0] == raw["stars"]["m1"]["R"] and 0 < x[1] < 0.1  # k2 from the stellar model
    d = dynfit.apply(raw, names, x * 1.1)
    assert d["stars"]["m1"]["R"] == pytest.approx(1.1 * x[0])
    assert d["stars"]["m2"]["apsidal_k2"] == pytest.approx(1.1 * x[1])
    with pytest.raises(ValueError):
        dynfit.start_values(raw, ["stars.m1.L"])


def test_conjunctions_of_a_wider_orbit_and_their_prediction():
    """Third-body events of TIC 120362137: conjunctions of the 51-d orbit (B against the
    inner pair) are found, alternate in type, and are predicted with impact < 1."""
    tic = ROOT / "configs" / "squadra_test" / "tic120362137_anchored.yaml"
    raw = yaml.safe_load(tic.read_text())
    t0 = 2458683.3
    seeds = dynfit.Observed(
        np.array(["mid"] * 4),
        np.array([0, 1, 0, 1]),
        t0 + 100.0 + 25.65 * np.arange(4),
        np.ones(4),
        np.array(["s"] * 4),
    )
    sep = np.full(4, np.nan)
    tc = dynfit.model_times(Config.model_validate(raw), seeds, t0, sep)
    assert np.all(np.isfinite(tc)) and np.all(np.diff(tc) > 15.0)
    assert np.all(sep * 215.032 < 5.0)  # Rsun: B passes across the inner pair
    obs = dynfit.Observed(
        np.array(["inner"] * 3),
        np.zeros(3, int),
        t0 + np.array([10.0, 13.3, 16.6]),
        np.full(3, 1e-4),
        np.array(["x"] * 3),
    )
    names = ["orbits.mid.M0"]
    x = dynfit.start_values(raw, names)
    res = dynfit.Stage2Result(
        names,
        x,
        x,
        np.array([0.05]),
        np.array([[0.0025]]),
        obs,
        obs.time,
        obs.time,
        0,
        0,
        0,
        0,
        0,
        raw,
        t0,
    )
    pr = dynfit.predict_minima(
        res,
        str(tic.parent),
        years=0.3,
        pairs=["mid"],
        start=t0 + 200.0,
        samples=12,
        processes=1,
        log=lambda *a: None,
    )
    assert set(pr.pair.tolist()) == {"mid"} and len(pr.t) >= 3
    assert np.all(pr.impact < 1.0) and np.all(np.isnan(pr.duration))
    assert set(pr.kind.tolist()) == {0, 1}


def test_mapped_mesh_fit_converges_to_the_direct_mesh_fit(tmp_path):
    """Output space mapping: the fit on a coarse mesh plus (accurate - coarse) reaches the
    solution of the fit on the accurate mesh (contact binary VW LMi, inclination free)."""
    from swingersim.photometry import light_curve

    t0 = 2459000.0
    raw = yaml.safe_load((ROOT / "configs" / "catalog" / "vw_lmi_contact.yaml").read_text())
    raw["run"] = {"t_max_yr": 1.0, "n_output": 11}
    cfg = Config.model_validate(raw)
    from swingersim.timing import eclipse_times

    ecl = eclipse_times(cfg, 0.05)[0]
    with open(tmp_path / "minima.csv", "w") as fh:
        fh.write("pair,type,time,sigma,source\n")
        fh.writelines(
            f"{e.binary},{e.kind[0]},{t0 + e.t_obs_yr * 365.25:.7f},0.0002,syn\n" for e in ecl[::3]
        )
    t = np.linspace(1.0, 1.48, 60)
    rng = np.random.default_rng(2)
    f = light_curve(cfg, t, "TESS").flux + rng.normal(0, 5e-4, len(t))
    with open(tmp_path / "lc.csv", "w") as fh:
        fh.write("time,flux,sigma\n")
        fh.writelines(f"{t0 + t[k]:.6f},{f[k]:.6f},0.0005\n" for k in range(len(t)))
    results = {}
    for fid in ("mapped", "direct"):
        start = yaml.safe_load(yaml.safe_dump(raw))
        start["orbits"]["A"]["i"] += 1.0
        start["fit"] = {
            "minima": "minima.csv",
            "lc": [
                {
                    "file": "lc.csv",
                    "fidelity": fid,
                    "mesh_resolution": [32, 64],
                    "coarse_resolution": [16, 32],
                }
            ],
            "stage2": {"bjd_t0": t0, "free": ["orbits.A.i"], "processes": 1, "max_iterations": 8},
            "output": {"dir": f"out_{fid}"},
        }
        p = tmp_path / f"{fid}.yaml"
        p.write_text(yaml.safe_dump(start, sort_keys=False))
        results[fid] = dynfit.run_stage2(p, log=lambda *a: None, plots=False)
    m, d = results["mapped"], results["direct"]
    assert m.fidelity and m.fidelity[-1]["change_sigma"] < 0.3
    assert abs(m.x[0] - d.x[0]) < 0.5 * d.err[0] + 1e-3, (m.x, d.x, d.err)
    assert abs(d.x[0] - raw["orbits"]["A"]["i"]) < 4 * d.err[0] + 0.05
    assert "mapped meshes" in dynfit.report(m)
