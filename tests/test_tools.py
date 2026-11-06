"""The helper scripts of tools/ produce published numbers: tests of their building blocks."""

import sys
from pathlib import Path

import numpy as np
import pytest

TOOLS = Path(__file__).resolve().parents[1] / "tools"
sys.path.insert(0, str(TOOLS))

import _common as C
import tess_minima


def test_calendar_years_and_epochs(tmp_path):
    assert C.year(2451544.5) == pytest.approx(2000.0)
    assert C.bjd(C.year(2459010.25)) == pytest.approx(2459010.25)
    p = tmp_path / "c.yaml"
    p.write_text(
        "# t = 0 of this configuration = BJD 2459000.5\nfit: {stage2: {bjd_t0: 2459001.0}}\n"
    )
    assert C.fit_epoch(p) == 2459001.0  # fit.stage2 before the comment
    p.write_text("fit_epoch_bjd_tdb: 2459002.0\nfit: {stage2: {bjd_t0: 2459001.0}}\n")
    assert C.fit_epoch(p) == 2459002.0  # the key written by swingersim fit wins
    p.write_text("# t = 0 = BJD 2459003.25\nmasses: {}\n")
    assert C.fit_epoch(p) == 2459003.25
    p.write_text("masses: {}\n")
    with pytest.raises(ValueError):
        C.fit_epoch(p)


def test_file_cache_recomputes_when_an_input_changes(tmp_path):
    src = tmp_path / "in.txt"
    src.write_text("a")
    out = tmp_path / "out.csv"
    calls = []

    def write(tmp):
        calls.append(1)
        tmp.write_text(src.read_text())

    assert C.file_cached(out, write, files=[src], params={"x": 1})
    assert not C.file_cached(out, write, files=[src], params={"x": 1})
    assert C.file_cached(out, write, files=[src], params={"x": 2})  # other parameters
    import os
    import time

    src.write_text("bb")
    os.utime(src, (time.time() + 5, time.time() + 5))
    assert C.file_cached(out, write, files=[src], params={"x": 2})  # changed input
    assert out.read_text() == "bb" and len(calls) == 3


def test_angles_switches_and_nodal_period():
    assert C.wrap_deg([190.0, -190.0, 10.0]).tolist() == pytest.approx([-170.0, 170.0, 10.0])
    years = np.arange(2000.0, 2010.0, 0.5)
    inc = np.linspace(80.0, 90.0, len(years))
    sw = C.eclipse_switches(years, inc, np.full(len(years), 85.0))
    assert len(sw) == 1 and 2004.0 < sw[0] < 2006.0
    # a unit vector precessing about z with a period of 40 yr, 30 deg from the axis
    t = np.linspace(0.0, 20.0, 41)
    ph = 2 * np.pi * t / 40.0
    h = np.c_[np.sin(0.5) * np.cos(ph), np.sin(0.5) * np.sin(ph), np.full_like(t, np.cos(0.5))]
    assert C.nodal_period(t, h, [0.0, 0.0, 1.0]) == pytest.approx(40.0, rel=1e-6)


def test_tess_minima_recovers_eclipse_times():
    rng = np.random.default_rng(2)
    P, tc, hw = 2.5, 0.37, 0.12
    t = np.arange(0.0, 25.0, 2.0 / 1440.0)
    true = tc + P * np.arange(-1, 12) + rng.normal(0, 20 / 86400, 13)  # ~20 s timing jitter

    def flux(tt):
        f = np.ones_like(tt)
        for tm in true:
            x = (tt - tm) / 0.06
            f -= 0.3 * np.exp(-0.5 * x * x)
        return f

    f = flux(t) + rng.normal(0, 5e-4, len(t))
    out = tess_minima.fit_eclipses(t, f, tc, P, hw, clip=4.0)
    assert len(out) >= 9
    for tm, err, _n, _rms in out:
        k = np.argmin(np.abs(true - tm))
        assert abs(tm - true[k]) < 5 * err + 5 / 86400  # within the stated error
        assert err < 60 / 86400
