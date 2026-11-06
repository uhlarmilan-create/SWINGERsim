"""Centre-of-light radial velocities (fit.rv.centre: light; swingersim.rvlight)."""

import warnings

import numpy as np
import pytest
import yaml

from swingersim import dynfit
from swingersim.config import Config, Stage2Spec
from swingersim.radvel import ObservedRV, model_rv
from swingersim.rvlight import light_offset, light_rv

warnings.simplefilter("ignore")

DETACHED = """
topology: custom
hierarchy: {A: [m1, m2]}
masses: {m1: 1.5, m2: 1.0}
orbits: {A: {P: 3.0, i: 89.5}}
stars: {m1: {R: 2.0, Teff: 7000}, m2: {R: 1.2, Teff: 5800}}
physics: {gr: false}
run: {t_max_yr: 0.1}
"""


def test_spheres_out_of_eclipse_and_rossiter_mclaughlin():
    cfg = Config.model_validate(yaml.safe_load(DETACHED))
    t = np.linspace(0.0, 3.0, 241)
    off = light_offset(cfg, ["m1"] * len(t), t)
    out = np.abs(off) < 1e-9
    assert out.mean() > 0.8  # a symmetric star out of eclipse: exactly the centre of mass
    # in eclipse: the Rossiter-McLaughlin effect, antisymmetric about mid-eclipse, of the
    # size v_rot (R_occulter / R)^2 ~ 34 km/s x 0.36
    ecl = np.flatnonzero(~out)
    assert 6.0 < off[ecl].max() < 12.0 and -12.0 < off[ecl].min() < -6.0
    assert off[ecl].max() == pytest.approx(-off[ecl].min(), rel=0.05)
    np.testing.assert_allclose(light_rv(cfg, ["m1"] * 5, t[:5], 2459000.0),
                               model_rv(cfg, ["m1"] * 5, 2459000.0 + t[:5], 2459000.0), atol=1e-9)  # fmt: skip


def test_contact_binary_masses_need_the_centre_of_light():
    """RVs of a contact binary made with the centre of light: fitting them with the centre of
    mass gives biased masses; with fit.rv.centre: light the masses come back."""
    from swingersim.validation_v027 import _contact_cfg

    cfg = _contact_cfg()
    raw = yaml.safe_load(yaml.safe_dump(cfg.model_dump(exclude_none=True)))
    raw["stars"] = {"m1": {"R": 1.0, "Teff": 6200}, "m2": {"R": 0.7, "Teff": 6200}}
    rng = np.random.default_rng(3)
    T0 = 2459000.0
    n = 100
    t = T0 + np.sort(rng.uniform(0.0, 3.0, n))
    stars = np.array(["m1", "m2"] * (n // 2))
    v = light_rv(cfg, stars, t - T0, T0) - 10.0 + rng.normal(0, 0.5, n)
    rv = ObservedRV(stars, t, v, np.full(n, 0.5), np.array(["s"] * n))
    empty = dynfit.Observed(
        np.zeros(0, str), np.zeros(0, int), np.zeros(0), np.zeros(0), np.zeros(0, str)
    )
    free = ["masses.m1", "masses.m2", "rv.gamma"]
    spec = Stage2Spec.model_validate(
        {"bjd_t0": T0, "free": free, "max_iterations": 15, "processes": 1}
    )
    out = {}
    for centre in ("mass", "light"):
        start = yaml.safe_load(yaml.safe_dump(raw))
        start["masses"] = {"m1": 1.0, "m2": 0.5}
        start["fit"] = {"minima": "x.csv", "rv": {"file": "x.csv", "gamma": 0.0, "centre": centre}}
        res = dynfit.fit_stage2(start, None, empty, spec, 1, lambda *a: None, mcmc=False, rv=rv)
        out[centre] = (res.x, res.err, res.chi2_red)
    x, e, chi_l = out["light"]
    assert x[0] == pytest.approx(1.2, abs=3 * e[0]) and x[1] == pytest.approx(0.504, abs=3 * e[1])
    assert chi_l < 1.5
    xm, _em, chi_m = out["mass"]
    # the centre of mass fits worse and gives masses too large by several per cent
    assert chi_m > 2 * chi_l
    assert xm[0] / 1.2 - 1 > 0.03 and xm[0] - 1.2 > abs(x[0] - 1.2) + 2 * e[0]
