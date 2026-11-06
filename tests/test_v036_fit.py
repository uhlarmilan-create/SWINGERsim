"""Fitting chain: one integration for all data types, model cache, emcee
convergence, degrees of freedom, outputs (dynfit, obsmodel, radvel, incfit, lcfit, ocfit)."""

import csv
import math
import warnings

import numpy as np
import pytest
import yaml

from swingersim import dynfit, ocfit, units
from swingersim.config import Config, Stage2Spec
from swingersim.incfit import ObservedInclination, model_inclination
from swingersim.obsmodel import KMS_PER_AUYR, observables
from swingersim.radvel import ObservedRV, model_rv

T0 = 2459000.0
TOY = """
topology: '2+2'
masses: {m1: 1.2, m2: 1.0, m3: 1.1, m4: 0.9}
orbits: {A: {P: 2.0, i: 89.0, M0: 40, e: 0.1, omega: 30}, B: {P: 3.0, i: 88.5, M0: 200},
         out: {P: 300.0, e: 0.3, i: 85.0, M0: 90, Omega: 10}}
stars: {m1: {R: 1.2}, m2: {R: 1.0}, m3: {R: 1.1}, m4: {R: 0.9}}
run: {t_max_yr: 0.3, n_output: 61}
"""


def _toy():
    return yaml.safe_load(TOY)


def _cfg(raw=None):
    return Config.model_validate(raw or _toy())


def _obs(times, pairs, kinds, sigma=2e-4, source="s"):
    n = len(times)
    return dynfit.Observed(
        np.asarray(pairs),
        np.asarray(kinds, int),
        np.asarray(times, float),
        np.full(n, sigma),
        np.array([source] * n if isinstance(source, str) else source),
    )


# ---------------------------------------------------------------- one integration
def test_one_integration_gives_the_separate_models():
    cfg = _cfg()
    rng = np.random.default_rng(2)
    tt = T0 + np.sort(rng.uniform(-60, 60, 40))
    obs = _obs(tt, np.array(["A", "B"] * 20), np.array([0, 1, 1, 0] * 10))
    rs = np.array(["m1", "m3", "m2", "m4"] * 5)
    rt = T0 + np.sort(rng.uniform(-60, 60, 20))
    it = T0 + np.linspace(-60, 60, 7)
    o = observables(cfg, T0, minima=obs, rv=(rs, rt), incs=[("A", it, True), ("B", it, False)])
    tm = dynfit.model_times(cfg, obs, T0)
    assert np.all(np.isfinite(tm))
    assert np.max(np.abs(o["t"] - tm)) * 86400 < 1e-4  # seconds: integrator round-off
    rv = model_rv(cfg, rs, rt, T0)
    assert np.max(np.abs(o["rv"] - rv)) < 1e-7
    ia = model_inclination(cfg, "A", it, T0, with_limit=True)
    assert o["inc"][0].shape == (2, len(it))
    assert np.max(np.abs(o["inc"][0] - ia)) < 1e-9
    assert np.max(np.abs(o["inc"][1] - model_inclination(cfg, "B", it, T0))) < 1e-9


def _rv_reference(cfg, star, t_bjd):
    """Velocity at the converged emission time: fresh integration to t_e, iterated to 1e-6 s."""
    from swingersim.backends.nbody import build_simulation

    T = (t_bjd - T0) / units.DAYS_PER_YEAR
    k = int(star[1:]) - 1
    te = T
    for _ in range(8):
        sim = build_simulation(cfg)[0]
        sim.collision = "none"
        if te < 0:
            sim.dt = -abs(sim.dt)
        sim.integrate(te, exact_finish_time=1)
        p = sim.particles[k]
        new = T - p.z / units.C
        if abs(new - te) < 1e-6 / units.SECONDS_PER_YEAR:
            break
        te = new
    return p.vz * KMS_PER_AUYR


def test_rv_light_time_matches_the_converged_emission_time():
    raw = _toy()
    raw["orbits"]["out"]["P"] = 3000.0  # light time ~ 25 min, RV curve of the 2-d pair steep
    cfg = _cfg(raw)
    t = T0 + np.array([-40.3, -40.3, -12.1, 3.7, 3.7, 25.9, 40.0])
    s = np.array(["m1", "m3", "m2", "m1", "m4", "m3", "m2"])
    v = model_rv(cfg, s, t, T0)
    ref = np.array([_rv_reference(cfg, a, b) for a, b in zip(s, t, strict=True)])
    assert np.max(np.abs(v - ref)) < 1e-6  # km/s


def test_merged_model_isolates_a_failing_data_type():
    raw = _toy()
    raw["fit"] = {"minima": "x.csv"}
    obs = _obs([T0 + 1.0, T0 + 2.0], ["nosuchorbit", "A"], [0, 0])
    rv = ObservedRV(
        np.array(["m1"]), np.array([T0 + 1.0]), np.zeros(1), np.ones(1), np.array(["X"])
    )
    ctx = {"obs": obs, "rv": rv, "inc": [], "lc": [], "bjd_t0": T0, "base_dir": None}
    m = dynfit.model_at(raw, ctx)
    assert m["t"] is None  # unknown orbit: the minima fail ...
    assert m["rv"] is not None and np.isfinite(m["rv"][0])  # ... the velocities still count


# ---------------------------------------------------------------- degrees of freedom
def test_error_factors_share_the_parameters_over_all_data():
    rng = np.random.default_rng(1)
    r = 2.0 * rng.normal(size=20)
    s = np.ones(20)
    g = np.array(["rv:X"] * 20)
    alone = dynfit.error_factors(r, s, g, "per_source", "scale", 10)["rv:X"]
    shared = dynfit.error_factors(r, s, g, "per_source", "scale", 10, n_total=400)["rv:X"]
    chi2 = float(np.sum(r**2))
    assert alone == pytest.approx(max(math.sqrt(chi2 / 10.0), 1.0))  # 20 - 10 dof (old)
    assert shared == pytest.approx(max(math.sqrt(chi2 / (20 - 10 * 20 / 400)), 1.0))
    assert shared < alone


def test_chi2_red_counts_neither_priors_nor_satisfied_limits():
    q = ObservedInclination(
        np.array([T0, T0 + 1, T0 + 2]),
        np.array([80.0, 70.0, np.nan]),
        np.full(3, 0.5),
        np.array(["a"] * 3, dtype=object),
        "A",
        "t",
        np.array(["", "upper", "no_eclipse"], dtype=object),
    )
    obs = _obs(T0 + np.arange(10.0), ["A"] * 10, [0] * 10)
    pb = dynfit._Problem(_toy(), ["masses.m1"], None, T0, {"masses.m1": (1.2, 0.1)}, obs, None,
                         [], [q])  # fmt: skip
    r = np.concatenate([np.ones(10), [1.0, 0.0, 0.0], [5.0]])  # minima, inc rows, prior
    chi2, red, dof = dynfit._chi2(pb, r, 1)
    assert dof == 10 + 1 - 1  # 10 minima + 1 measured inclination - 1 parameter
    assert chi2 == pytest.approx(36.0) and red == pytest.approx(11.0 / 10)


def test_covariance_warns_about_dropped_directions():
    J = np.array([[1.0, 2.0, 1.0], [2.0, 4.0, 0.5], [3.0, 6.0, 0.1], [1.0, 2.0, 2.0]])
    msgs = []
    cov = dynfit.covariance(J, log=msgs.append, names=["a", "b", "c"])
    assert len(msgs) == 1 and "1 direction" in msgs[0]
    assert np.all(np.isfinite(cov))
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        dynfit.covariance(J)
    assert any("not constrained" in str(x.message) for x in w)
    assert dynfit.covariance(np.eye(3)) == pytest.approx(np.eye(3))


# ---------------------------------------------------------------- emcee
class _LinearProblem:
    """Stub of dynfit._Problem: y = a + b t with Gaussian errors (fast emcee test)."""

    def __init__(self):
        self.t = np.linspace(-1, 1, 30)
        self.y = 0.3 + 2.0 * self.t + 0.1 * np.random.default_rng(4).normal(size=30)
        self.calls = 0

    def n_data(self):
        return 30

    def models(self, X, pool, store=True):
        self.calls += 1
        return [x for x in np.atleast_2d(X)]

    def resid(self, m, x):
        return (self.y - (m[0] + m[1] * self.t)) / 0.1


def test_run_mcmc_rounds_walkers_and_reports_tau():
    pb = _LinearProblem()
    x = np.array([0.3, 2.0])
    cov = np.diag([0.02**2, 0.03**2])
    mc = dynfit.run_mcmc(pb, x, np.sqrt(np.diag(cov)), cov, np.array([-5.0, -5.0]),
                         np.array([5.0, 5.0]), 1.0, 7, 60, 3, None, 3, lambda *a: None)  # fmt: skip
    assert mc["walkers"] == 8  # max(7, 2 n + 2 = 6) -> even; independent of the 3 cores
    one = dynfit.run_mcmc(_LinearProblem(), x, np.sqrt(np.diag(cov)), cov, np.array([-5.0, -5.0]),
                          np.array([5.0, 5.0]), 1.0, 7, 60, 3, None, 1, lambda *a: None)  # fmt: skip
    assert one["walkers"] == 8 and np.array_equal(one["samples"], mc["samples"])  # same chain
    assert len(mc["samples"]) == 8 * 30 and len(mc["log_prob"]) == len(mc["samples"])
    assert len(mc["tau"]) == 2 and np.isfinite(mc["tau_max"])
    short = len(mc["samples"]) / 8 < 10 * mc["tau_max"]  # 30 steps after burn-in
    assert ("shorter than 10 autocorrelation times" in mc["warning"]) == short
    assert pb.calls <= 2 * 60 + 1  # one vectorised call per half-step
    assert np.all(np.isfinite(mc["log_prob"]))


def test_run_mcmc_extends_the_chain_to_the_tau_target():
    pb = _LinearProblem()
    x = np.array([0.3, 2.0])
    cov = np.diag([0.02**2, 0.03**2])
    lines = []
    mc = dynfit.run_mcmc(pb, x, np.sqrt(np.diag(cov)), cov, np.array([-5.0, -5.0]),
                         np.array([5.0, 5.0]), 1.0, 8, 40, 3, None, 2, lines.append,
                         tau_target=50.0, max_steps=3000)  # fmt: skip
    assert mc["steps"] > 40  # extended
    assert mc["steps"] >= 50 * mc["tau_max"] or mc["steps"] == 3000
    assert len(mc["samples"]) == mc["walkers"] * (mc["steps"] - mc["steps"] // 2)
    assert any("more steps" in ln for ln in lines) and not mc["warning"]
    capped = dynfit.run_mcmc(pb, x, np.sqrt(np.diag(cov)), cov, np.array([-5.0, -5.0]),
                             np.array([5.0, 5.0]), 1.0, 8, 40, 3, None, 2, lambda *a: None,
                             tau_target=1e6, max_steps=70)  # fmt: skip
    assert capped["steps"] == 70


def test_pool_honours_swingersim_cores(monkeypatch):
    monkeypatch.setenv("SWINGERSIM_CORES", "1")
    pool, n = dynfit._pool(0, {})
    assert pool is None and n == 1


# ---------------------------------------------------------------- fit and outputs
def _toy_data(n_min=24, n_rv=12, span=60.0, seed=5):
    raw = _toy()
    raw["fit"] = {"minima": "x.csv", "rv": {"file": "rv.csv", "gamma": 5.0}}
    cfg = _cfg(raw)
    rng = np.random.default_rng(seed)
    tt = T0 + np.sort(rng.uniform(-span, span, n_min))
    obs = _obs(tt, np.array(["A", "B"] * (n_min // 2)), rng.integers(0, 2, n_min),
               source=np.array(["s1", "s2", "s3"] * (n_min // 3)))  # fmt: skip
    tm = dynfit.model_times(cfg, obs, T0)
    ok = np.isfinite(tm)
    obs = dynfit.Observed(obs.pair[ok], obs.kind[ok], tm[ok] + 2e-4 * rng.normal(size=ok.sum()),
                          obs.sigma[ok], obs.source[ok])  # fmt: skip
    obs.time[3] += 0.02  # an outlier
    rs = np.array(["m1", "m3"] * (n_rv // 2))
    rt = T0 + np.sort(rng.uniform(-span, span, n_rv))
    v = model_rv(cfg, rs, rt, T0) + 5.0 + rng.normal(size=n_rv)
    rv = ObservedRV(rs, rt, v, np.ones(n_rv), np.array(["X"] * n_rv))
    return raw, obs, rv


def _spec(**kw):
    d = {
        "bjd_t0": T0,
        "free": ["orbits.A.P", "orbits.A.M0", "rv.gamma"],
        "max_iterations": 3,
        "processes": 1,
        "refits": 3,
    }
    d.update(kw)
    return Stage2Spec.model_validate(d)


def test_robust_fit_reuses_models_and_writes_the_outputs(tmp_path, monkeypatch):
    raw, obs, rv = _toy_data()
    start = yaml.safe_load(yaml.safe_dump(raw))
    start["orbits"]["A"]["P"] = 2.0 + 1e-5
    start["fit"]["rv"]["gamma"] = 4.5
    spec = _spec(error_scale="per_source", clip_sigma=5.0)
    calls = []
    real = dynfit.model_at

    def counted(d, ctx):
        calls.append((d["orbits"]["A"]["P"], d["orbits"]["A"]["M0"], d["fit"]["rv"]["gamma"]))
        return real(d, ctx)

    monkeypatch.setattr(dynfit, "model_at", counted)
    res = dynfit.fit_stage2_robust(start, None, obs, spec, "per_source", 1, lambda *a: None, rv=rv)
    assert res.outliers is not None and len(res.outliers) == 1
    assert len(res.passes) >= 2
    # every parameter set computed once: the start model and the Jacobian of the next fit in a
    # row (at the solution of the previous one) and the outputs come from the cache
    assert res.n_eval == len(calls) == len(set(calls))
    assert res.base_dir is None and res.dof == len(res.obs) + len(res.rv) - 3
    assert res.chi2_red == pytest.approx(
        float(np.sum(((res.obs.time - res.t_model) / res.obs.sigma) ** 2)
              + np.sum(((res.rv.rv - res.rv_model) / res.rv.sigma) ** 2)) / res.dof,
        rel=1e-6,
    )  # fmt: skip
    files = dynfit.write_outputs(res, tmp_path, "start.yaml", plots=False)
    names = {f.name for f in files}
    assert {"stage2_covariance.csv", "stage2_params.csv", "stage2_fitted.yaml"} <= names
    rows = list(csv.reader((tmp_path / "stage2_covariance.csv").read_text().splitlines()))
    assert rows[0] == ["parameter", *spec.free] and len(rows) == 4
    params = list(csv.DictReader((tmp_path / "stage2_params.csv").read_text().splitlines()))
    assert {"error_linear", "p16", "p50", "p84"} <= set(params[0])
    text = (tmp_path / "stage2_fitted.yaml").read_text()
    assert "SWINGERsim " in text and "BJD_TDB" in text and "sampler none" in text
    fitted = yaml.safe_load(text)
    if dynfit.EPOCH_KEY in Config.model_fields:
        assert fitted[dynfit.EPOCH_KEY] == T0
    else:  # until config.py has the key: in the header only
        assert dynfit.EPOCH_KEY not in fitted and f"BJD_TDB {T0}" in text
    rep = (tmp_path / "stage2_report.txt").read_text()
    assert "errors of all data scaled" in rep and "degrees of freedom" in rep


def test_mcmc_runs_once_after_the_fits_in_a_row(monkeypatch, tmp_path):
    raw, obs, rv = _toy_data()
    spec = _spec(error_scale="per_source", clip_sigma=5.0, sampler="emcee", walkers=8,
                 mcmc_steps=20, max_iterations=2)  # fmt: skip
    seen = []

    def fake(pb, x, err, cov, lo, hi, scale, walkers, steps, seed, pool, ncore, log, **kw):
        seen.append(x.copy())
        s = x + 1e-6 * np.random.default_rng(0).normal(size=(40, len(x)))
        return {
            "samples": s,
            "log_prob": -np.arange(40.0),
            "acceptance": 0.4,
            "tau": [3.0] * len(x),
            "tau_max": 3.0,
            "walkers": 8,
            "steps": 20,
            "discard": 10,
            "seed": seed,
            "warning": "",
        }

    monkeypatch.setattr(dynfit, "run_mcmc", fake)
    res = dynfit.fit_stage2_robust(raw, None, obs, spec, "per_source", 1, lambda *a: None, rv=rv)
    assert len(seen) == 1 and np.all(seen[0] == res.x)  # around the last refit, no extra fit
    assert res.err_linear is not None and res.mcmc["tau_max"] == 3.0
    dynfit.write_outputs(res, tmp_path, "x.yaml", plots=False)
    head = (tmp_path / "stage2_samples.csv").read_text().splitlines()[0].split(",")
    assert head[-1] == "log_prob" and head[:-1] == spec.free
    p = list(csv.DictReader((tmp_path / "stage2_params.csv").read_text().splitlines()))
    assert all(r["p50"] != "" for r in p)
    assert "tau max 3.0" in (tmp_path / "stage2_report.txt").read_text()


def test_fidelity_passes_keep_the_first_start_and_sample_once(monkeypatch):
    from swingersim import lcfit

    raw, obs, rv = _toy_data(n_min=12, n_rv=6)
    t = T0 + np.linspace(1.0, 1.4, 30)
    lc = lcfit.ObservedLC(t, np.ones(30), np.full(30, 1e-3), np.array(["a"] * 30), "lc")
    spec = _spec(error_scale="per_source", sampler="emcee", max_iterations=1)
    monkeypatch.setattr(lcfit, "needs_mesh", lambda d: True)  # spheres treated as a mapped mesh
    changes = iter([5.0, 2.0, 1.0])
    monkeypatch.setattr(dynfit, "FIDELITY_PASSES", 2)
    real_set = dynfit._set_deltas
    monkeypatch.setattr(
        dynfit, "_set_deltas", lambda pb, x, idx: real_set(pb, x, idx) * 0.0 + next(changes)
    )
    n_mcmc = []
    monkeypatch.setattr(dynfit, "_mcmc_info", lambda *a, **k: n_mcmc.append(1) or None)
    res = dynfit._fit_with_fidelity(raw, None, obs, spec, "per_source", 1, lambda *a: None, rv,
                                    [lc])  # fmt: skip
    assert len(res.fidelity) == 2 and "not converged" in res.note
    assert len(n_mcmc) == 1
    assert np.all(res.x0 == dynfit.start_values(raw, spec.free))
    first = dynfit.model_times(_cfg(raw), res.obs, T0)
    assert np.max(np.abs(res.t_start - first)) * 86400 < 1e-3  # start of the first pass
    assert "accurate meshes" in res.fidelity_note


# ---------------------------------------------------------------- inclinations, plots
def test_measured_inclinations_above_90_are_folded(tmp_path):
    from swingersim.incfit import load

    p = tmp_path / "i.csv"
    p.write_text("time,i_deg,sigma_deg,limit\n2459000,95.0,1.0,\n2459100,80.0,1.0,upper\n")

    class S:
        file = str(p)
        min_sigma_deg = 0.0
        pair = "A"
        label = None

    q = load(S, None)
    assert q.inc[0] == pytest.approx(85.0) and q.inc[1] == 80.0
    p.write_text("time,i_deg,sigma_deg,limit\n2459000,95.0,1.0,upper\n")
    with pytest.raises(ValueError, match="ambiguous"):
        load(S, None)


def test_plot_rv_with_seven_stars(tmp_path, monkeypatch):
    from swingersim import radvel

    raw = _toy()
    raw["fit"] = {"minima": "x.csv", "rv": {"file": "rv.csv"}}
    stars = np.array([f"m{k}" for k in range(1, 8)] * 2)
    rv = ObservedRV(stars, T0 + np.arange(14.0), np.zeros(14), np.ones(14), np.array(["X"] * 14))
    monkeypatch.setattr(radvel, "model_rv", lambda cfg, s, t, t0: np.zeros(len(t)))
    obs = _obs([T0], ["A"], [0])
    res = dynfit.Stage2Result(["orbits.A.P"], np.array([2.0]), np.array([2.0]), np.ones(1),
                              np.eye(1), obs, obs.time, obs.time, 0, 0, 0, 0, 0, raw, T0)  # fmt: skip
    res.rv, res.rv_model, res.rv_start = rv, np.zeros(14), np.zeros(14)
    assert dynfit.plot_rv(res, tmp_path / "rv.png").exists()  # KeyError for m7 before


# ---------------------------------------------------------------- stage 1
def _stage1_result(tmp_path):
    path = tmp_path / "m.csv"
    ocfit.synthetic_minima(
        path,
        {"A": {"terms": ["lite", "apsidal"], "T0": 2450000.0, "P": 2.5, "e": 0.1, "w0": 1.0,
               "dw": 0.002, "asini": 1.0}},
        {"P3": 900.0, "T3": 2450100.0, "e3": 0.3, "w3": 1.0},
        2450000.0,
        2453000.0,
        every=7,
    )  # fmt: skip
    mins = ocfit.read_minima(path)
    pm = ocfit.PairModel("A", ["lite", "apsidal"])
    model = ocfit.Model([pm], False)
    n = ocfit.assign_epochs(mins.time, mins.kind, 2450000.0, 2.5)
    t_ref = float(np.round(np.median(mins.time)))
    model.t_ref = t_ref
    d = ocfit.PairData(pm, n, mins.kind, mins.time - t_ref, mins.sigma, mins.source,
                       2450000.0 - t_ref, 2.5)  # fmt: skip
    return model, [d]


def test_stage1_predictions_vectorised_match_the_loop(tmp_path):
    model, data = _stage1_result(tmp_path)
    x = []
    for nm in model.names:
        k = nm.split(".", 1)[1]
        x.append({"T0": data[0].T0, "P": 2.5, "e": 0.1, "w0": 1.0, "dw": 0.002, "asini": 1.0,
                  "P3": 900.0, "T3": 2450100.0 - model.t_ref, "e3": 0.3, "w3": 1.0}[k])  # fmt: skip
    x = np.array(x)
    cov = np.diag((1e-4 * np.maximum(np.abs(x), 1e-3)) ** 2)
    res = ocfit.FitResult(list(model.names), x, np.sqrt(np.diag(cov)), cov, 0.0, 1, model, data)
    p = tmp_path / "pred.csv"
    ocfit.write_predictions(res, p, years=0.2)
    rows = list(csv.DictReader(p.read_text().splitlines()))
    # the loop of 0.35, one minimum at a time, with half the 16-84 % range
    rng = np.random.default_rng(2)
    samp = rng.multivariate_normal(res.x, res.cov, size=300)
    d = data[0]
    n0 = math.ceil((float(d.t.max()) - d.T0) / d.P)
    exp = []
    for n in range(n0, n0 + int(0.2 * 365.25 / d.P) + 1):
        for j in (0, 1):
            nn, kk = np.array([n], float), np.array([j])
            t = float(model.predict(res.x, d.model, nn, kk)[0])
            t = float(model.predict(res.x, d.model, nn, kk, np.array([t]))[0])
            ts = [float(model.predict(s, d.model, nn, kk, np.array([t]))[0]) for s in samp]
            exp.append(
                (n + j / 2, t + model.t_ref, 0.5 * np.subtract(*np.percentile(ts, [84, 16])))
            )
    assert len(rows) == len(exp)
    for r, (e, t, s) in zip(rows, exp, strict=True):
        assert float(r["epoch"]) == e
        assert float(r["time_bjd_tdb"]) == pytest.approx(t, abs=2e-6)
        assert float(r["sigma_d"]) == pytest.approx(s, abs=2e-6)


def test_stage1_covariance_from_the_svd(tmp_path, monkeypatch):
    model, data = _stage1_result(tmp_path)
    seen = []
    real = dynfit.covariance

    def spy(J, log=None, names=None):
        seen.append(real(J, log, names))
        return seen[-1]

    monkeypatch.setattr(dynfit, "covariance", spy)
    res = ocfit.fit_stage1(model, data)
    red = max(res.chi2 / res.dof, 1.0)
    assert len(seen) == 1 and np.allclose(res.cov, seen[0] * red)


def test_run_stage2_sampler_override(monkeypatch, tmp_path):
    seen = {}

    def stop(raw, base_dir, obs, spec, *a, **k):
        seen["sampler"] = spec.sampler
        raise RuntimeError("stop")

    monkeypatch.setattr(dynfit, "_fit_with_fidelity", stop)
    raw, obs, _ = _toy_data(n_min=6, n_rv=2)
    (tmp_path / "x.csv").write_text(
        "pair,type,time,sigma\n" + "".join(f"A,p,{t},0.001\n" for t in obs.time[:4])
    )
    (tmp_path / "rv.csv").write_text("star,time,rv,sigma\nm1,2459001.0,1.0,1.0\n")
    raw["fit"]["stage2"] = {"bjd_t0": T0, "free": ["orbits.A.P"], "processes": 1}
    p = tmp_path / "c.yaml"
    p.write_text(yaml.safe_dump(raw))
    with pytest.raises(RuntimeError, match="stop"):
        dynfit.run_stage2(p, log=lambda *a: None, sampler="emcee")
    assert seen["sampler"] == "emcee"
