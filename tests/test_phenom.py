"""Phenomenological minima (Mikulasek 2015) of synthetic light curves (F4, phase 1)."""

from __future__ import annotations

import csv
import math
import warnings

import numpy as np
import pytest

from swingersim import phenom
from swingersim.validation import _romer_cfg, phenom_romer_delay_expected


def test_profile_shape():
    dt = np.linspace(-0.2, 0.2, 401)
    y = phenom.profile(dt, 0.3, 0.02, 1.2)
    assert y[200] == pytest.approx(0.3)
    assert np.allclose(y, y[::-1])
    assert y[0] < 1e-6 and np.all(y >= 0) and np.all(y <= 0.3 + 1e-12)
    # the correction C deepens the wings, Gamma sharpens the minimum
    assert phenom.profile(0.02, 0.3, 0.02, 1.2, 0.5) > phenom.profile(0.02, 0.3, 0.02, 1.2)
    # a larger Gamma gives a flatter, box-like minimum
    assert phenom.profile(0.01, 0.3, 0.02, 2.0) > phenom.profile(0.01, 0.3, 0.02, 1.0)


def _segments(shifts, rng=None, noise=0.0):
    dt = np.arange(-0.12, 0.12, 1.0 / 1440.0)
    out = []
    for k, sh in enumerate(shifts):
        base = 10.0 + 0.01 * k - 0.02 * dt + 0.5 * dt**2
        mag = base + phenom.profile(dt - sh, 0.12, 0.018, 1.3, 0.1)
        if noise:
            mag = mag + rng.normal(0.0, noise, len(dt))
        out.append((dt, mag))
    return out


def test_shape_and_minimum_are_recovered():
    segs = _segments([0.0, 1e-3, -2e-3])
    s = phenom._fit_shape(segs, 0.15)
    assert s.depth == pytest.approx(0.12, rel=1e-4)
    assert s.width == pytest.approx(0.018, rel=1e-4)
    assert s.gamma == pytest.approx(1.3, rel=1e-3)
    assert s.corr == pytest.approx(0.1, abs=1e-3)
    off, _, depth, width, rms = phenom._fit_minimum(*segs[1], s)
    assert width == pytest.approx(0.018, rel=1e-5)
    assert off == pytest.approx(1e-3, abs=1e-8)
    assert depth == pytest.approx(0.12, rel=1e-5) and rms < 1e-8


def test_minimum_error_matches_the_scatter():
    """The covariance error of the time of minimum agrees with the scatter of noisy fits."""
    rng = np.random.default_rng(3)
    segs = _segments([0.0] * 60, rng, noise=2e-3)
    s = phenom._fit_shape(segs[:10], 0.15)
    fits = [phenom._fit_minimum(dt, m, s) for dt, m in segs]
    offs = np.array([f[0] for f in fits])
    sig = np.median([f[1] for f in fits])
    assert np.std(offs) == pytest.approx(sig, rel=0.25)
    assert abs(np.mean(offs)) < 3 * sig / math.sqrt(len(offs))


@pytest.fixture(scope="module")
def romer():
    return phenom.synthetic_minima(_romer_cfg(), 0.08, ["A"], cadence_min=2.0, proximity=False)


def test_minima_show_the_romer_delay(romer):
    """The minima of the two stars are displaced by -/+ half the Kaplan (2010) delay from the
    barycentric conjunctions; swingersim timings include this delay (timing.roemer_in_binary), so
    they agree with the minima of the light curve."""
    half = phenom_romer_delay_expected() / 2.0
    assert half > 10.0  # a delay that matters: ~14 s
    for kind in ("primary", "secondary"):
        d = [(x.t_min_d - x.t_geom_d) * 86400.0 for x in romer.minima if x.kind == kind]
        assert len(d) >= 5
        assert abs(np.mean(d)) < 0.03 * half
        assert np.std(d) < 0.1
    assert romer.ephemeris["A"][1] == pytest.approx(5.0, rel=1e-6)
    assert romer.n_blended == 0


def test_csv_outputs(romer, tmp_path):
    phenom.write_minima_csv(romer, tmp_path / "m.csv")
    with open(tmp_path / "m.csv", encoding="utf-8") as f:  # '#' header lines since 0.36
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    assert len(rows) == len(romer.minima)
    assert {r["kind"] for r in rows} == {"primary", "secondary"}
    assert "skipped" not in phenom.summarize(romer)


def test_seasonal_shapes_and_blending():
    """CzeV4315-like 2+2: in the total light, eclipses of A overlapping those of B are skipped
    and the rest masked; with isolated light (ideal disentangling) none is lost."""
    from swingersim.config import load_config

    cfg = load_config("configs/examples/example_timings_czev4315.yaml")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = phenom.synthetic_minima(
            cfg, 0.03, ["B"], cadence_min=4.0, season_days=5.0, isolate=False
        )
        iso = phenom.synthetic_minima(cfg, 0.03, ["B"], cadence_min=4.0, proximity=False)
    assert res.n_blended >= 1 and iso.n_blended == 0
    assert len(iso.minima) > len(res.minima)
    assert res.seasons and all(s["n"] >= 2 for s in res.seasons)
    widths = [s["width_d"] for s in res.seasons if s["kind"] == "primary"]
    assert np.std(widths) < 0.02 * np.mean(widths)  # no evolution in a few days
    for s in res.shapes.values():
        assert s.rms < 1e-3


def test_cli_minima(tmp_path):
    from swingersim.cli import main

    cfg = tmp_path / "c.yaml"
    import yaml

    cfg.write_text(yaml.safe_dump(_romer_cfg().model_dump(exclude_none=True, mode="json")))
    out = tmp_path / "m.csv"
    assert main(["minima", str(cfg), "--years", "0.03", "-o", str(out), "--binaries", "A"]) == 0
    assert out.exists()


def test_light_time_of_a_binary_far_from_the_centre():
    """QZ Car: the eclipsing binary Ac is ~35 AU from the centre of the system (z/c ~ 0.2 d for a
    6-d binary). Its minima must still follow swingersim timings - the barycentric conjunctions
    displaced by the Romer delay inside the binary (Kaplan 2010, +-16 s here); a first-order
    light-time correction of the positions shifted them by ~200 s."""
    import math

    from swingersim import units
    from swingersim.config import load_config

    cfg = load_config("configs/benchmarks/qzcar_minima.yaml")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = phenom.synthetic_minima(cfg, 0.06, ["A"], band="V", cadence_min=15.0)
    m1, m2 = cfg.masses["m1"], cfg.masses["m2"]
    a = units.period_to_sma(cfg.orbits["A"].P / units.DAYS_PER_YEAR, m1 + m2)
    half = a * math.sin(math.radians(93.0)) * (m2 - m1) / (m1 + m2) / units.C
    half *= units.SECONDS_PER_YEAR  # +16 s: the heavier Ac2 is behind at the primary eclipse
    assert abs(half) > 10.0
    for kind in ("primary", "secondary"):
        d = [(x.t_min_d - x.t_geom_d) * 86400.0 for x in res.minima if x.kind == kind]
        assert abs(np.mean(d)) < 8.0
