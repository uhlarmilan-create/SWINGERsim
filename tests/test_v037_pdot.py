"""orbits.<pair>.Pdot: secular period change of a two-star binary (swingersim.pdot)."""

import warnings

import numpy as np
import pytest
import yaml

from swingersim import dynfit, units
from swingersim.config import Config, Stage2Spec

warnings.simplefilter("ignore")

BINARY = """
topology: custom
hierarchy: {A: [m1, m2]}
masses: {m1: 1.5, m2: 0.8}
orbits: {A: {P: 2.5, e: 0.0, i: 89.0, M0: 0}}
physics: {gr: false}
run: {t_max_yr: 30, n_output: 3}
"""
TRIPLE = """
topology: custom
hierarchy: {out: [A, m3], A: [m1, m2]}
masses: {m1: 1.5, m2: 0.8, m3: 0.7}
orbits:
  A: {P: 2.5, e: 0.0, i: 89.0, M0: 30}
  out: {P: 2000.0, e: 0.3, i: 80.0, omega: 40, M0: 100}
physics: {gr: true}
run: {t_max_yr: 30, n_output: 3}
"""
T0 = 2459000.0


def _primaries(cfg, years):
    from swingersim.timing import eclipse_times

    ecl, _ = eclipse_times(cfg, years)
    return np.array(sorted(e.t_obs_yr for e in ecl if e.binary == "A" and e.kind == "primary"))


@pytest.mark.parametrize("mode", ["angular_momentum", "conservative"])
@pytest.mark.parametrize("years", [30.0, -30.0])
def test_eclipse_times_follow_the_parabola(mode, years):
    """T(E) = T0 + P0 E + (1/2) P0 Pdot E^2 (Kwee 1958): the quadratic coefficient to 1e-4,
    forward and backward in time."""
    d = yaml.safe_load(BINARY)
    d["orbits"]["A"].update(Pdot=4e-6, Pdot_mode=mode)
    t = _primaries(Config.model_validate(d), years) * units.DAYS_PER_YEAR
    E = np.round((t - t[0]) / 2.5)
    c2, b, _ = np.polyfit(E, t, 2)
    assert c2 / (0.5 * 2.5 * 4e-6 / units.DAYS_PER_YEAR) == pytest.approx(1.0, abs=2e-4)
    # P at t = 0 is the configured one: the slope at the epoch
    e0 = np.round((0.0 - t[0]) / 2.5)
    assert b + 2 * c2 * e0 == pytest.approx(2.5, rel=1e-8)


def test_conservative_transfer_changes_the_masses():
    from swingersim.backends.nbody import build_simulation
    from swingersim.pdot import PeriodDrift

    d = yaml.safe_load(BINARY)
    d["orbits"]["A"].update(Pdot=4e-6, Pdot_mode="conservative")
    cfg = Config.model_validate(d)
    sim, _, _, root, masses = build_simulation(cfg)
    drift = PeriodDrift.of(cfg, root, masses)
    drift.integrate(sim, 30.0)
    m1, m2 = sim.particles[0].m, sim.particles[1].m
    rate = 4e-6 / units.DAYS_PER_YEAR / units.days_to_years(2.5)  # (dP/dt)/P [1/yr]
    mdot = rate * 1.5 * 0.8 / (3 * (1.5 - 0.8))  # dP/P = 3 dM1 (M1 - M2) / (M1 M2)
    assert m1 + m2 == pytest.approx(2.3, rel=1e-12)
    assert m1 - 1.5 == pytest.approx(mdot * 30.0, rel=0.02)  # the last interval's value


@pytest.mark.slow  # ~100 s
def test_fit_recovers_pdot_with_a_third_body():
    raw = yaml.safe_load(TRIPLE)
    truth = 3e-6
    raw["orbits"]["A"]["Pdot"] = truth
    cfg = Config.model_validate(raw)
    rng = np.random.default_rng(1)
    tt = T0 + np.sort(rng.uniform(-8000.0, 8000.0, 160))
    obs = dynfit.Observed(
        np.array(["A"] * 160),
        rng.integers(0, 2, 160),
        tt,
        np.full(160, 2e-4),
        np.array(["s"] * 160),
    )
    tm = dynfit.model_times(cfg, obs, T0)
    ok = np.isfinite(tm)
    obs = dynfit.Observed(obs.pair[ok], obs.kind[ok], tm[ok] + 2e-4 * rng.normal(size=ok.sum()),
                          obs.sigma[ok], obs.source[ok])  # fmt: skip
    start = yaml.safe_load(yaml.safe_dump(raw))
    start["orbits"]["A"]["Pdot"] = 0.0
    free = ["orbits.A.P", "orbits.A.M0", "orbits.A.Pdot"]
    spec = Stage2Spec.model_validate(
        {"bjd_t0": T0, "free": free, "max_iterations": 20, "processes": 1}
    )
    res = dynfit.fit_stage2(start, None, obs, spec, 1, lambda *a: None, mcmc=False)
    k = res.names.index("orbits.A.Pdot")
    assert res.x[k] == pytest.approx(truth, abs=3 * res.err[k])
    assert res.err[k] < 0.05 * truth


def test_pdot_needs_a_two_star_orbit_and_nbody():
    d = yaml.safe_load(TRIPLE)
    d["orbits"]["out"]["Pdot"] = 1e-6
    with pytest.raises(ValueError, match="only an orbit of two stars"):
        Config.model_validate(d)
    d = yaml.safe_load(TRIPLE)
    d["orbits"]["A"]["Pdot"] = 1e-6
    d["run"]["backend"] = "secular"
    from swingersim.backends.secular import run_secular

    with pytest.raises(ValueError, match="needs run.backend: nbody"):
        run_secular(Config.model_validate(d))


def test_light_curve_and_rv_see_the_drift():
    """The light curve and the radial velocities use the same drifting orbit: the eclipse in the
    light curve moves with the parabola."""
    from swingersim.photometry import light_curve

    d = yaml.safe_load(BINARY)
    d["stars"] = {"m1": {"R": 1.6}, "m2": {"R": 0.8}}
    d["orbits"]["A"]["Pdot"] = 7e-6  # a shift of ~0.2 d after 20 yr
    cfg = Config.model_validate(d)
    shift = 0.5 * (7e-6 / units.DAYS_PER_YEAR / 2.5) * (20.0 * units.DAYS_PER_YEAR) ** 2  # [d]
    t_pred = _primaries(cfg, 20.5) * units.DAYS_PER_YEAR
    tp = t_pred[np.argmin(np.abs(t_pred - 20.0 * units.DAYS_PER_YEAR))]
    t = tp + np.linspace(-0.2, 0.2, 801)
    lc = light_curve(cfg, t)
    assert abs(t[np.argmin(lc.flux)] - tp) < 2e-3
    d0 = yaml.safe_load(BINARY)
    d0["stars"] = d["stars"]
    t0 = _primaries(Config.model_validate(d0), 20.5) * units.DAYS_PER_YEAR
    assert abs((tp - t0[np.argmin(np.abs(t0 - tp))]) - shift) < 0.01 * shift
