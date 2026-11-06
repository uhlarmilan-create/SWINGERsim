"""swingersim-PHOEBE hybrid light curves (F4 b, tier C); skipped without the phoebe package."""

from __future__ import annotations

import warnings

import numpy as np
import pytest

phoebe = pytest.importorskip("phoebe")

from swingersim.phoebe_lc import PASSBANDS, phoebe_light_curve
from swingersim.photometry import light_curve
from swingersim.validation import _phoebe_test_cfg


def test_passband_names():
    assert PASSBANDS["V"] == "Johnson:V" and PASSBANDS["TESS"] == "TESS:T"


def test_same_physics_as_tier_a_around_the_primary_eclipse():
    """Spheres, black bodies and the same limb darkening: PHOEBE and tier A agree; this
    checks incl = i, per0 = omega + 180 deg and the time of periastron."""
    cfg = _phoebe_test_cfg()
    a = light_curve(cfg, np.linspace(0, 5, 201), band="V", proximity=False)
    t = np.linspace(-0.12, 0.12, 49) + a.t_day[np.argmin(a.flux)]
    a = light_curve(cfg, t, band="V", proximity=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p = phoebe_light_curve(
            cfg,
            t,
            band="V",
            atm="blackbody",
            distortion="sphere",
            irradiation=False,
            swingersim_ld=True,
        )
    assert a.flux.min() < 0.8  # a real eclipse
    d = a.flux - p.flux
    assert np.max(np.abs(d)) < 1.5e-3
    assert np.max(np.abs(d - np.median(d))) < 1e-3
    assert set(p.components) == {"A", "B"}  # PHOEBE light of both binaries
