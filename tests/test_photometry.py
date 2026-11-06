"""F4 (b): synthetic photometry (tier A)."""

import numpy as np
import pytest

from swingersim.config import load_config
from swingersim.photometry import (
    bands,
    beaming_alpha,
    blocked_fraction,
    claret_coefficients,
    covered_fraction,
    ld_mean,
    light_curve,
    surface_brightness,
)

QUAD = np.array([0.0, 0.4 + 0.5, 0.0, -0.25])  # quadratic law u = (0.4, 0.25)


def test_passbands():
    assert {"TESS", "V", "B", "R", "I", "Gaia_G"} <= set(bands())
    assert surface_brightness(7000.0, "V") > surface_brightness(5000.0, "V")
    # hotter stars are relatively brighter in blue
    r_hot = surface_brightness(8000.0, "B") / surface_brightness(8000.0, "I")
    r_cool = surface_brightness(5000.0, "B") / surface_brightness(5000.0, "I")
    assert r_hot > r_cool
    assert 0.5 < 3.0 - beaming_alpha(6000.0, "TESS") < 6.0
    with pytest.raises(ValueError):
        surface_brightness(6000.0, "XYZ")


def test_claret_coefficients():
    a, y = claret_coefficients(5800.0, 4.4, 0.0, "V")
    mu = np.linspace(0, 1, 11)
    prof = 1 - sum(ak * (1 - mu ** ((k + 1) / 2)) for k, ak in enumerate(a))
    assert np.all(np.diff(prof) > 0) and 0.2 < prof[0] < 0.6
    assert 0.0 < y < 1.0
    aT, _ = claret_coefficients(5800.0, 4.4, 0.0, "TESS")  # proxy: I band
    aI, _ = claret_coefficients(5800.0, 4.4, 0.0, "I")
    assert np.allclose(aT, aI)


def test_blocked_fraction_limits():
    assert blocked_fraction(1.0, QUAD, [(0.0, 0.0, 1.5)]) == pytest.approx(1.0, abs=1e-12)
    assert blocked_fraction(1.0, QUAD, [(3.0, 0.0, 0.5)]) == 0.0
    one = blocked_fraction(1.0, QUAD, [(0.3, 0.2, 0.2)])
    twice = blocked_fraction(1.0, QUAD, [(0.3, 0.2, 0.2), (0.3, 0.2, 0.2)])
    assert twice == pytest.approx(one, rel=1e-10)
    a = blocked_fraction(1.0, QUAD, [(0.5, 0.0, 0.2)])
    b = blocked_fraction(1.0, QUAD, [(-0.5, 0.0, 0.2)])
    ab = blocked_fraction(1.0, QUAD, [(0.5, 0.0, 0.2), (-0.5, 0.0, 0.2)])
    assert ab == pytest.approx(a + b, rel=1e-9)
    # uniform disk: geometric overlap area
    assert blocked_fraction(1.0, np.zeros(1), [(0.0, 0.0, 0.5)]) == pytest.approx(0.25, rel=1e-9)
    assert ld_mean(QUAD) == pytest.approx(1 - 0.4 / 3 - 0.25 / 6)


def test_covered_fraction_union():
    r = np.array([0.5])
    f = covered_fraction(r, [(0.5, 0.0, 0.3), (0.0, 0.5, 0.3)])
    f1 = covered_fraction(r, [(0.5, 0.0, 0.3)])
    assert f[0] == pytest.approx(2 * f1[0], rel=1e-12)


def test_light_curve_quadruple(tmp_path):
    cfg = load_config("configs/squadra/czev4315_ensemble.yaml")
    t = np.arange(0, 1.0, 5 / 1440)
    lc = light_curve(cfg, t, band="V")
    assert lc.flux.shape == t.shape
    assert 0.5 < lc.flux.min() < lc.flux.max() < 1.1
    assert np.allclose(sum(lc.components.values()), lc.flux)
    flat = light_curve(cfg, t, band="V", proximity=False)
    assert flat.flux.max() <= 1.0 + 1e-12
    lc.write_csv(tmp_path / "lc.csv")
    text = (tmp_path / "lc.csv").read_text()
    assert text.startswith("# SWINGERsim")  # metadata header (since 0.36)
    cols = next(ln for ln in text.splitlines() if not ln.startswith("#"))
    assert cols.startswith("t_day,flux,flux_m1")
    ld = {f"m{k}": [0.0, 0.9, 0.0, -0.25] for k in range(1, 5)}
    lc2 = light_curve(cfg, t[:10], band="V", ld=ld)
    assert np.allclose(lc2.stars[0].ld, ld["m1"])


def test_hot_star_limb_darkening_from_tlusty():
    from swingersim.photometry import ld_mean, limb_darkening, tlusty_coefficients

    a, y = limb_darkening(40000.0, 4.0, 0.0, "V")  # beyond the Claret & Bloemen grid
    assert np.allclose(a, tlusty_coefficients(40000.0, 4.0, 0.0, "V"))
    assert 0.2 < y < 0.5
    # the switch at 15 000 K is smooth: TLUSTY and ATLAS agree to 1 % in <I>/I(1)
    a_t = tlusty_coefficients(15500.0, 4.0, 0.0, "V")
    a_c, _ = claret_coefficients(15500.0, 4.0, 0.0, "V")
    assert ld_mean(a_t) == pytest.approx(ld_mean(a_c), rel=0.01)
    with pytest.raises(ValueError, match="no TLUSTY"):
        tlusty_coefficients(45000.0, 2.0, 0.0, "V")  # beyond the Eddington limit
