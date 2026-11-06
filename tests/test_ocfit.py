import json
import math
import warnings

import numpy as np
import pytest

from swingersim import ocfit
from swingersim.config import Config, load_config
from swingersim.validation_v023 import gimenez_bastero_shift


def test_apsidal_matches_series_and_zero_for_circular():
    n = np.arange(0, 2000, 5.0)
    for j in (0, 1):
        k = np.full(len(n), j)
        a = ocfit.apsidal_shift(n, k, 3.0, 0.1, 0.4, 3e-3)
        b = gimenez_bastero_shift(n, k, 3.0, 0.1, 0.4, 3e-3)
        assert np.ptp(a - b) * 86400 < 1e-4
    assert np.all(ocfit.apsidal_shift(n, np.zeros(len(n), int), 3.0, 0.0, 0.4, 3e-3) == 0)


def test_inclination_shifts_minima_of_eccentric_pair():
    # for i < 90 deg the minimum of the projected separation moves away from conjunction
    a = ocfit.apsidal_shift(np.array([0.0]), np.array([0]), 3.0, 0.3, 0.5, 0.0, math.pi / 2)
    b = ocfit.apsidal_shift(np.array([0.0]), np.array([0]), 3.0, 0.3, 0.5, 0.0, math.radians(80))
    assert 0 < abs(a[0] - b[0]) * 86400 < 600


def test_time_conversion_hjd_to_bjd():
    t = np.array([2459000.5])
    bjd = ocfit.to_bjd_tdb(t, "HJD_UTC", 83.8, -5.4)
    # TDB - UTC = 69.184 s in 2020; HJD and BJD differ by a few seconds at most
    assert abs((bjd[0] - t[0]) * 86400 - 69.184) < 10


def test_read_minima_header_and_types(tmp_path):
    p = tmp_path / "m.csv"
    p.write_text(
        "# system: X\n# time_system: BJD_TDB\npair,type,time,sigma,source,method,depth,width\n"
        "A,p,2459000.1,0.001,TESS,individual,0.3,0.1\nA,s,2459001.6,0.001,ZTF,season,,\n"
    )
    m = ocfit.read_minima(p)
    assert list(m.kind) == [0, 1] and m.header["system"] == "X"
    assert math.isnan(m.depth[1])


def test_stage1_recovers_light_time_orbit(tmp_path):
    outer = {"P3": 900.0, "T3": 2457300.0, "e3": 0.4, "w3": 2.0}
    pairs = {
        "A": {"terms": ["lite"], "T0": 2457000.1, "P": 2.1, "asini": 0.8},
        "B": {"terms": ["lite"], "T0": 2457000.3, "P": 1.7, "asini": 0.6},
    }
    path = tmp_path / "m.csv"
    ocfit.synthetic_minima(path, pairs, outer, 2457000, 2459800, sigma=1e-4, every=25)
    cfg = Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 1, "m2": 1, "m3": 1, "m4": 1},
            "orbits": {"A": {"P": 2.1}, "B": {"P": 1.7}, "out": {"P": 1000}},
            "fit": {
                "minima": str(path),
                "stage1": {"terms": ["lite"], "shared_outer": True},
                "output": {
                    "dir": str(tmp_path / "out"),
                    "tables": ["csv", "latex"],
                    "predict": {"years": 1},
                },
            },
        }
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ocfit.run_fit(cfg, plots=False)
    p, e = res.params(), res.errors()
    assert p["out.P3"] == pytest.approx(900, abs=4 * e["out.P3"] + 1)
    assert p["out.e3"] == pytest.approx(0.4, abs=4 * e["out.e3"] + 0.005)
    assert p["A.asini"] == pytest.approx(0.8, rel=0.02)
    assert p["B.asini"] / p["A.asini"] == pytest.approx(0.75, rel=0.02)
    for f in (
        "stage1_parameters.csv",
        "stage1_parameters.tex",
        "stage1_oc.csv",
        "stage1_summary.json",
        "stage1_predictions.csv",
    ):
        assert (tmp_path / "out" / f).exists()


def test_fit_section_without_run(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text(
        "topology: '2+2'\nmasses: {m1: 1, m2: 1, m3: 1, m4: 1}\n"
        "orbits: {A: {P: 2}, B: {P: 3}, out: {P: 500}}\nfit: {minima: m.csv}\n"
    )
    cfg = load_config(p)
    assert cfg.fit.stage1.terms == ["lite"] and cfg.run.t_max_yr == 100.0
    with pytest.raises(ValueError):
        Config.model_validate(
            {
                "topology": "2+2",
                "masses": {"m1": 1, "m2": 1, "m3": 1, "m4": 1},
                "orbits": {"A": {"P": 2}, "B": {"P": 3}, "out": {"P": 500}},
                "fit": {"minima": "m.csv", "time_system": "HJD_UTC"},
            }
        )


def test_stage1_clips_outliers(tmp_path):
    outer = {"P3": 900.0, "T3": 2457300.0, "e3": 0.4, "w3": 2.0}
    pairs = {"A": {"terms": ["lite"], "T0": 2457000.1, "P": 2.1, "asini": 0.8}}
    path = tmp_path / "m.csv"
    ocfit.synthetic_minima(path, pairs, outer, 2457000, 2459800, sigma=1e-4, every=25)
    lines = path.read_text().splitlines()
    head = next(k for k, line in enumerate(lines) if not line.startswith("#"))
    cols = lines[head].split(",")
    it = cols.index("time")
    bad = [head + 5, head + 20, head + 33]
    for k in bad:
        c = lines[k].split(",")
        c[it] = f"{float(c[it]) + 0.02:.7f}"  # 29 min
        lines[k] = ",".join(c)
    path.write_text("\n".join(lines) + "\n")
    base = {
        "topology": "2+2",
        "masses": {"m1": 1, "m2": 1, "m3": 1, "m4": 1},
        "orbits": {"A": {"P": 2.1}, "B": {"P": 1.7}, "out": {"P": 1000}},
        "fit": {
            "minima": str(path),
            "stage1": {"terms": ["lite"], "clip_sigma": 5.0},
            "output": {"dir": str(tmp_path / "out")},
        },
    }
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = ocfit.run_fit(Config.model_validate(base), plots=False)
    assert len(res.outliers) == 3
    assert all(abs(o["residual_min"] - 28.8) < 2 for o in res.outliers)
    assert res.chi2 / res.dof < 1.6
    assert res.params()["A.asini"] == pytest.approx(0.8, rel=0.02)
    assert "outlying minima left out" in ocfit.report(res)
    s = json.loads((tmp_path / "out" / "stage1_summary.json").read_text())
    assert len(s["outliers"]) == 3


def test_times_before_1960_use_delta_t():
    """Before 1960 the 'UTC' times are Universal Time: TT = UT + Delta T (Espenak & Meeus 2006),
    not ERFA's TT - UTC = 32.184 s; from 1960 on astropy's UTC is kept."""
    jd1900 = 2415020.5 + 0.5 * 365.25  # 1900.5
    assert ocfit.delta_t_s(np.array([2451545.0 - 100 * 365.25]))[0] == pytest.approx(-2.79, abs=0.02)
    assert ocfit.delta_t_s(np.array([2451545.0 - 50 * 365.25]))[0] == pytest.approx(29.07, abs=0.02)
    assert 33.0 < ocfit.delta_t_s(np.array([ocfit.JD_UTC_1960 - 1]))[0] < 33.3
    b = ocfit.to_bjd_tdb(np.array([jd1900]), "JD_UTC", 0.0, 60.0)
    ref = ocfit.to_bjd_tdb(np.array([2459000.5]), "JD_UTC", 0.0, 60.0)
    back = ocfit.from_bjd_tdb(b, 0.0, 60.0)
    assert abs(back[0] - jd1900) * 86400 < 1e-3
    assert abs(ocfit.from_bjd_tdb(ref, 0.0, 60.0)[0] - 2459000.5) * 86400 < 1e-3
    # TT - UT at 1900.5 is Delta T ~ -2 s: the BJD differs from the geocentric JD by the light time
    # (|.| < 8.4 min) and Delta T only, not by 32.184 s more
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astropy.time import Time

    star = SkyCoord(0.0 * u.deg, 60.0 * u.deg)
    tt = Time(jd1900 + ocfit.delta_t_s(np.array([jd1900]))[0] / 86400, format="jd", scale="tt")
    with ocfit._quiet_erfa():
        expect = tt.tdb.jd + ocfit._light_time(tt.tdb, star, "barycentric")
    assert abs(b[0] - expect) * 86400 < 1e-3
    assert np.ndim(ocfit.to_bjd_tdb(2415100.0, "HJD_UTC", 0.0, 60.0)) == 0
