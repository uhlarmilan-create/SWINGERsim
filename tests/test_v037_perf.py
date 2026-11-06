"""Speed-ups of the code audit without changing results: dense output of the light-curve
integration, scalar clips of the stellar tracks."""

import warnings

import numpy as np
import yaml

from swingersim.config import Config

warnings.simplefilter("ignore")

ECC = """
topology: '2+2'
masses: {m1: 2.5, m2: 2.4, m3: 2.6, m4: 2.3}
orbits: {A: {P: 2.94, e: 0.6, i: 86, omega: 101, M0: 280}, B: {P: 3.26, e: 0.22, i: 84, omega: 257}, out: {P: 121.7, e: 0.26, i: 85}}
stars: {m1: {R: 2.3, Teff: 9000}, m2: {R: 2.0, Teff: 8800}, m3: {R: 2.3, Teff: 9100}, m4: {R: 1.9, Teff: 8700}}
run: {t_max_yr: 1.0}
"""


def test_dense_light_curve_equals_exact_integration():
    from swingersim import photometry

    cfg = Config.model_validate(yaml.safe_load(ECC))
    t = np.concatenate([np.arange(-3.0, 3.0, 4 / 1440), np.arange(60.0, 62.0, 4 / 1440)])
    dense = photometry.light_curve(cfg, t, exposure_min=2.0, n_sub=3)
    orig = photometry.SeenStates.__init__
    photometry.SeenStates.__init__ = lambda self, c, lt=True: orig(self, c, lt, dense=False)
    try:
        exact = photometry.light_curve(cfg, t, exposure_min=2.0, n_sub=3)
    finally:
        photometry.SeenStates.__init__ = orig
    assert 1.0 - exact.flux.min() > 0.1  # the eclipses are in
    assert np.max(np.abs(dense.flux - exact.flux)) < 1e-8
    assert "Hermite" in dense.meta["positions"]


def test_scalar_clip_matches_numpy():
    from swingersim.stellar import _clipf

    for x, lo, hi in ((0.3, 0.0, 1.0), (-2.0, 0.0, 1.0), (5.0, 0.0, 1.0), (7, 0, 3)):
        assert _clipf(x, lo, hi) == float(np.clip(x, lo, hi))
