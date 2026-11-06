"""Swingersim gui 0.36: caches, decimation that keeps extremes, CSV/HDF5 readers, draft paths,
model curves, quick check, help texts and the server fixes of the 0.36 review."""

from __future__ import annotations

import json
import os
import threading
import time
import urllib.error
import urllib.request
from pathlib import Path

import h5py
import numpy as np
import pytest
import yaml

from swingersim.gui import assist
from swingersim.gui import digest as D
from swingersim.gui import modelviews as MV
from swingersim.gui import server as gs

APP = Path(gs.__file__).with_name("app.html")
CFG = (
    "topology: '2+2'\nmasses: {m1: 1.2, m2: 1.0, m3: 1.1, m4: 0.9}\n"
    "orbits: {A: {P: 2.0, i: 89.0}, B: {P: 3.0, i: 88.5}, out: {P: 200.0, e: 0.2, i: 89.0}}\n"
    "stars: {m1: {R: 1.2}, m2: {R: 1.0}, m3: {R: 1.1}, m4: {R: 0.9}}\n"
    "run: {t_max_yr: 10}\n"
)


# --------------------------------------------------------------------------- decimation
def test_envelope_keeps_a_narrow_peak_that_striding_misses():
    n = 100_000
    y = np.sin(np.linspace(0, 20, n))
    y[12_347] = 9.0  # one sample: a short eccentricity peak
    y[77_001] = -7.0
    stride = D._down(n, 1500)
    assert y[stride].max() < 2  # the old decimation stepped over both
    idx = D._decimate([y], 1500)
    assert len(idx) <= 1500 and idx[0] == 0 and idx[-1] == n - 1
    assert y[idx].max() == 9.0 and y[idx].min() == -7.0
    # several series share the indices and keep every one's extremes
    z = np.cos(np.linspace(0, 3, n))
    z[5] = 4.0
    idx = D._decimate([y, z], 3000)
    assert len(idx) <= 3000 and y[idx].max() == 9.0 and z[idx].max() == 4.0
    assert len(D._decimate([y], 10 * n)) == n  # short series stay whole


def test_envelope_keeps_gaps():
    y = np.arange(1000.0)
    y[500:520] = np.nan
    idx = D._envelope(y, 20)
    assert np.isnan(y[idx]).any()


# --------------------------------------------------------------------------- readers
def test_read_csv_skips_comment_headers_and_empty_cells(tmp_path):
    p = tmp_path / "t.csv"
    p.write_text(
        "# written by SWINGERsim\n# columns: ...\nbinary,kind,t_obs_yr,o_c_d,impact\n"
        "A,primary,0.1,1e-4,0.2\nA,secondary,0.2,,0.3\n# a remark\nB,primary,0.3,2e-4,\n"
    )
    d = D._read_csv(p)
    assert list(d["binary"]) == ["A", "A", "B"]
    assert np.isnan(d["o_c_d"][1]) and d["t_obs_yr"].dtype == float
    assert np.isnan(d["impact"][2])
    # without pandas the csv module gives the same
    import builtins

    real = builtins.__import__

    def no_pandas(name, *a, **k):
        if name == "pandas":
            raise ImportError
        return real(name, *a, **k)

    builtins.__import__ = no_pandas
    try:
        d2 = D._read_csv(p)
    finally:
        builtins.__import__ = real
    assert list(d2["binary"]) == ["A", "A", "B"]
    np.testing.assert_array_equal(d2["o_c_d"], d["o_c_d"])
    # the Outputs page recognises the file behind its comment lines
    assert gs.classify(p) == "timings" or gs._csv_header(p)[0] == "binary"


def test_lightcurve_with_header_comments_and_an_empty_one(tmp_path):
    p = tmp_path / "lc.csv"
    t = np.arange(0, 3, 1 / 720)
    f1 = np.where(np.abs((t % 1.0) - 0.5) < 0.002, 0.2, 0.5)  # a 6-minute eclipse each day
    rows = "\n".join(f"{a:.6f},{b + 0.5:.6f},{b:.6f},0.5" for a, b in zip(t, f1, strict=True))
    p.write_text("# swingersim lightcurve\nt_day,flux,flux_m1,flux_m2\n" + rows + "\n")
    assert gs.classify(p) == "lightcurve"
    old = D.LC_POINTS
    D.LC_POINTS = 300
    try:
        d = D.lightcurve_digest(p, None, {}, full=True)
    finally:
        D.LC_POINTS = old
    assert d["data"]["n"] == len(t) and len(d["data"]["t"]) <= 300
    assert min(d["data"]["flux"]) == pytest.approx(0.7)  # every eclipse bottom survived
    e = tmp_path / "e.csv"
    e.write_text("t_day,flux,flux_m1\n")
    c = D.lightcurve_digest(e, None, {}, full=True)
    assert c["data"]["n"] == 0 and c["facts"][0]["v"] == "–"


def test_timings_decimated_with_extremes_and_total_count(tmp_path, monkeypatch):
    n = 9000
    rows = ["binary,kind,epoch,t_obs_yr,o_c_d,P_d,T0_yr,impact"]
    for e in range(n):
        oc = 1e-3 * np.sin(e / 300.0) + (5e-3 if e == 4321 else 0.0)
        imp = "" if e == 7 else "0.1"
        rows.append(f"A,primary,{e},{e * 2 / 365.25:.9f},{oc:.9f},2.0,0.0,{imp}")
    p = tmp_path / "t.csv"
    p.write_text("\n".join(rows) + "\n")
    monkeypatch.setattr(D, "TIMING_POINTS", 1000)
    d = D.timings_digest(p, None, full=True)
    s = d["data"]["series"]["A"]["primary"]
    assert s["n_total"] == n and len(s["oc"]) <= 1000
    peak = (1e-3 * np.sin(4321 / 300.0) + 5e-3) * 1440
    assert max(s["oc"]) == pytest.approx(peak, rel=1e-6)  # the outlier is drawn
    b = d["data"]["binaries"][0]
    assert b["n_primary"] == n and b["oc_range"][1] == pytest.approx(max(s["oc"]), rel=1e-6)
    assert b["semi_amplitude"] == pytest.approx((b["oc_range"][1] - b["oc_range"][0]) / 2)
    assert "[min]" in d["spark"]["label"]


def test_minima_with_an_empty_cell(tmp_path):
    p = tmp_path / "m.csv"
    p.write_text(
        "binary,kind,epoch,t_min_d,sigma_d,t_geom_d,o_c_d,o_c_geom_d,depth_mag,width_d,rms_mag\n"
        "A,primary,0,1.0,1e-4,1.0,1e-5,0,0.1,0.2,1e-3\nA,primary,1,3.0,,3.0,,0,0.1,0.2,1e-3\n"
    )
    d = D.minima_digest(p, "", None, full=True)
    s = d["data"]["series"]["A"]["primary"]
    assert s["oc"][1] is None and s["sigma"][1] is None and s["oc"][0] == pytest.approx(0.864)
    assert "[s]" in d["spark"]["label"]


# --------------------------------------------------------------------------- run files
def _fake_run(path: Path, n: int = 50_000, hierarchy: str = "none", spike: int = 31_234):
    """A run file as io.save_run writes it (gzip-compressed datasets), with an eccentricity
    peak in one sample."""
    t = np.linspace(0, 100, n)
    with h5py.File(path, "w") as f:
        f.attrs["config_yaml"] = CFG
        f.attrs["binaries_json"] = json.dumps({})
        f.attrs["system_json"] = json.dumps({})
        f.attrs["backend"] = "nbody"
        kw = {"compression": "gzip"}
        f.create_dataset("time_yr", data=t, **kw)
        f.create_dataset("energy_error", data=1e-9 * np.ones(n), **kw)
        for name, a, inc in (("A", 0.035, 89.0), ("B", 0.046, 88.5), ("out", 0.95, 89.0)):
            g = f.create_group(f"orbits/{name}")
            e = 0.1 + 0.01 * np.sin(t)
            if name == "A" and spike < n:
                e[spike] = 0.93
            for k, v in (("a", a), ("inc", np.radians(inc)), ("Omega", 0.0), ("omega", 0.1)):
                g.create_dataset(k, data=v * np.ones(n), **kw)
            g.create_dataset("e", data=e, **kw)
            g.create_dataset("M", data=np.mod(t, 6.28), **kw)
        f.create_dataset("mutual_inclination/out-A", data=np.radians(5.0) * np.ones(n), **kw)
        for s, r in (("m1", 1.2), ("m2", 1.0), ("m3", 1.1), ("m4", 0.9)):
            f.create_group(f"stars/{s}").attrs["radius_rsun"] = r
        if hierarchy == "group":  # transitions only, no detected dataset
            h = f.create_group("hierarchy")
            h.attrs["configured"] = "((m1+m2)+(m3+m4))"
            h.attrs["transitions_json"] = "[]"


@pytest.mark.parametrize("hierarchy", ["none", "group"])
def test_run_digest_reads_compressed_files_without_hierarchy(tmp_path, hierarchy):
    p = tmp_path / "r.h5"
    _fake_run(p, hierarchy=hierarchy)
    card = D.run_digest(p, full=False)
    assert any(f["k"] == "Hierarchy changes" and f["v"] == "0" for f in card["facts"])
    d = D.run_digest(p, full=True)["data"]
    e = d["series"]["orbits"]["A"]["e"]
    assert e["range"][1] == pytest.approx(0.93)  # caption: the whole run
    assert max(e["y"]) == pytest.approx(0.93) and len(e["y"]) <= D.MAX_POINTS
    assert len(e["t"]) == len(e["y"]) and e["n_total"] == 50_000
    ecl = {x["binary"]: x for x in d["eclipses"]}
    assert len(ecl["A"]["t"]) == len(ecl["A"]["inc"]) == len(ecl["A"]["limit"])
    assert ecl["A"]["fraction"] == pytest.approx(1.0)  # i = 89 deg: always eclipsing


def test_card_reads_only_what_it_needs(tmp_path, monkeypatch):
    p = tmp_path / "r.h5"
    _fake_run(p)
    read = []
    orig = D._RunFile.get

    def spy(self, name):
        read.append(name)
        return orig(self, name)

    monkeypatch.setattr(D._RunFile, "get", spy)
    D.run_digest(p, full=False)
    assert "orbits/out/e" not in read and "orbits/A/M" not in read and "time_yr" not in read


def test_h5_series_keeps_the_extremes(tmp_path, monkeypatch):
    (tmp_path / "output").mkdir()
    _fake_run(tmp_path / "output" / "r.h5")
    monkeypatch.setattr(gs, "SERIES_LIMIT", 2000)
    s = gs.h5_series(tmp_path, "output/r.h5", ["orbits/A/e", "orbits/A/inc"])
    assert len(s["time_yr"]) <= 2000 and s["n_total"] == 50_000
    assert max(s["orbits/A/e"]) == pytest.approx(0.93)
    assert s["orbits/A/inc"][0] == pytest.approx(89.0)  # degrees


def test_run_viewer_reuses_the_page_of_the_run(tmp_path):
    out = tmp_path / "output"
    out.mkdir()
    h5 = out / "r.h5"
    _fake_run(h5, n=100)
    own = out / "r_viewer.html"
    own.write_text("<html>viewer of the run</html>")
    os.utime(own, ns=(h5.stat().st_mtime_ns + 10**9,) * 2)
    assert gs.run_viewer(tmp_path, "output/r.h5") == "output/r_viewer.html"


# --------------------------------------------------------------------------- fits
def _fit_folder(tmp_path: Path, lc: bool = True) -> Path:
    f = tmp_path / "fit"
    f.mkdir()
    (f / "stage2_report.txt").write_text(
        "  8 minima, 0 priors, 3 free parameters, 50 model evaluations, 12 s\n"
        "  chi^2 460.3 -> 189.3 (reduced 1.01)\n"
    )
    (f / "stage2_params.csv").write_text(
        "parameter,start,value,error,unit,prior_mean,prior_sigma,error_linear,p16,p50,p84\n"
        "orbits.inner.P,3.7747,3.77472,1e-08,d,,,1e-08,3.77471,3.77472,3.77473\n"
        "orbits.out.e,0.4,0.95,0.01,,,,0.01,,,\n"
        "masses.m3,3.0,5.0,0.5,Msun,,,0.4,4.5,5.0,5.5\n"
    )
    (f / "stage2_inclination.csv").write_text(
        "label,pair,time,i_obs,sigma,limit,i_start,i_model,i_eclipse_limit,z,source\n"
        "L,inner,2436800.0,85.4,1.0,,84.0,85.0,73.0,0.4,UBV\n"
        "L,B,2458777.5,80.0,0.3,,75.0,76.0,72.0,-13.3,TESS\n"
    )
    if lc:  # segment a: 2 points of residual 1e-3; segment b: 8 points of residual 0
        rows = ["bjd_tdb,flux,model,model_start,o_minus_model,source,segment"]
        rows += [f"{2459000 + k},1,1,1,0.001,TESS,a" for k in range(2)]
        rows += [f"{2459100 + k},1,1,1,0.0,TESS,b" for k in range(8)]
        (f / "stage2_lc.csv").write_text("\n".join(rows) + "\n")
    return f


def test_fit_bounds_rms_percentiles_and_pairs(tmp_path):
    f = _fit_folder(tmp_path)
    cfg = tmp_path / "sys_fit.yaml"
    cfg.write_text(
        yaml.safe_dump({"fit": {"minima": "m.csv", "stage2": {"bounds": {"masses.m3": [1, 9]}}}})
    )
    d = D.fit_digest(f, full=True, root=tmp_path, cfg_path=cfg)
    p = {x["name"]: x for x in d["data"]["params"]}
    assert not p["orbits.inner.P"]["tiny_error"]  # precise, not on a bound
    assert p["orbits.out.e"]["tiny_error"] and p["orbits.out.e"]["bounds"] == [0.0, 0.95]
    assert not p["masses.m3"]["tiny_error"] and p["masses.m3"]["bounds"] == [1.0, 9.0]
    assert p["masses.m3"]["p16"] == 4.5 and "p16" not in p["orbits.out.e"]
    assert any("orbits.out.e" in w and "on the bound" in w for w in d["data"]["warnings"])
    # the light curves: rms of all 10 points (sqrt(2/10) ppt), not of the two segment rms
    lc = next(x for x in d["facts"] if x["k"] == "light curves: rms")
    assert lc["v"] == f"{D.fmt(1e3 * np.sqrt(2e-6 / 10), 3)} ppt"
    assert d["data"]["inc_pairs"] == ["B", "inner"]
    assert {r["pair"] for r in d["data"]["inclination"]} == {"inner", "B"}
    # without the configuration of the fit only an error of 0 marks a bound
    d2 = D.fit_digest(f, full=True, root=tmp_path)
    assert not any(x["tiny_error"] for x in d2["data"]["params"])


# --------------------------------------------------------------------------- model curves
class _Ecl:
    def __init__(self, t, binary="inner", kind="primary"):
        self.t_obs_yr, self.binary, self.kind = t, binary, kind


def test_oc_curve_keeps_to_the_span(tmp_path, monkeypatch):
    from swingersim import timing

    calls = []

    def fake(cfg, years):
        calls.append(years)
        step = 0.01 if years > 0 else -0.01
        return [_Ecl(k * step) for k in range(int(abs(years) / 0.01) + 1)], {}

    monkeypatch.setattr(timing, "eclipse_times", fake)
    monkeypatch.setattr(MV, "_cfg", lambda path: None)
    cfg = tmp_path / "c.yaml"
    cfg.write_text(CFG)
    (tmp_path / "output").mkdir()
    t0 = 2459000.5  # 2020.0 is about BJD 2458849.5
    # a span that lies wholly before t = 0
    oc = MV.oc_curve(tmp_path, cfg, t0, 2019.0, 2019.5)
    yrs = (np.array(oc["eclipses"]["inner"]["primary"]) - t0) / 365.25
    a = (2451544.5 + 19 * 365.25 - t0) / 365.25
    b = (2451544.5 + 19.5 * 365.25 - t0) / 365.25
    assert yrs.min() >= a - 1e-9 and yrs.max() <= b + 1e-9 and calls == [a]
    # a span across t = 0: no eclipse twice
    oc = MV.oc_curve(tmp_path, cfg, t0, 2020.0, 2021.0)
    v = oc["eclipses"]["inner"]["primary"]
    assert len(v) == len(set(v))
    with pytest.raises(ValueError):
        MV.oc_curve(tmp_path, cfg, t0, 2021.0, 2020.0)


def test_model_cache_follows_version_and_referenced_files(tmp_path, monkeypatch):
    import swingersim

    tab = tmp_path / "stars.csv"
    tab.write_text("mass_msun,radius_rsun,apsidal_k2,rg2\n0.5,0.5,0.01,0.1\n2,2,0.01,0.1\n")
    cfg = tmp_path / "c.yaml"
    cfg.write_text(CFG + "stellar_model: {kind: table, file: stars.csv}\n")
    k1 = MV._cache(tmp_path, cfg, "system", {})
    monkeypatch.setattr(swingersim, "__version__", "9.9.9")
    k2 = MV._cache(tmp_path, cfg, "system", {})
    assert k1 != k2
    os.utime(tab, ns=(tab.stat().st_mtime_ns + 10**9,) * 2)
    assert MV._cache(tmp_path, cfg, "system", {}) != k2


def test_prune_cache(tmp_path):
    d = tmp_path / "output" / ".gui_cache" / "models"
    d.mkdir(parents=True)
    old, new = d / "old.json", d / "new.json"
    old.write_text("{}")
    new.write_text("{}")
    past = time.time() - 40 * 86400
    os.utime(old, (past, past))
    assert gs.prune_cache(tmp_path) == 1
    assert new.exists() and not old.exists()


# --------------------------------------------------------------------------- drafts
FIT_CFG = {
    "topology": "2+2",
    "stellar_model": {"kind": "table", "file": "tables/stars.csv"},
    "fit": {
        "minima": "data/minima.csv",
        "rv": {"file": "data/rv.csv"},
        "lc": [{"file": "data/lc1.csv"}, {"file": "/abs/lc2.csv"}],
        "inclination": [{"file": "data/inc.csv", "pair": "A"}],
        "output": {"dir": "out/fit"},
    },
}


def test_rebase_every_path_field(tmp_path):
    data = json.loads(json.dumps(FIT_CFG))
    changed = gs.rebase_paths(data, tmp_path / "configs" / "sys", tmp_path / "jobs" / "drafts")
    up = "../../configs/sys/"
    assert data["fit"]["minima"] == up + "data/minima.csv"
    assert data["fit"]["rv"]["file"] == up + "data/rv.csv"
    assert data["fit"]["lc"][0]["file"] == up + "data/lc1.csv"
    assert data["fit"]["lc"][1]["file"] == "/abs/lc2.csv"  # absolute: kept
    assert data["fit"]["inclination"][0]["file"] == up + "data/inc.csv"
    assert data["stellar_model"]["file"] == up + "tables/stars.csv"
    assert data["fit"]["output"]["dir"] == up + "out/fit"
    assert "fit.lc[0].file" in changed and len(changed) == 6
    e = {"stellar_model": {"kind": "python", "entry": "my/mod.py:props"}}
    gs.rebase_paths(e, tmp_path / "a", tmp_path / "b")
    assert e["stellar_model"]["entry"] == "../a/my/mod.py:props"
    e = {"stellar_model": {"kind": "python", "entry": "pkg.mod:props"}}
    assert gs.rebase_paths(e, tmp_path / "a", tmp_path / "b") == []


def test_save_draft_and_freeze_rewrite_the_paths(tmp_path):
    src = tmp_path / "configs" / "sys"
    src.mkdir(parents=True)
    (src / "sys_fit.yaml").write_text(yaml.safe_dump(FIT_CFG))
    text = (src / "sys_fit.yaml").read_text()
    r = gs.save_draft(tmp_path, "sysfit", text, False, source="configs/sys/sys_fit.yaml")
    assert r["rebased"] and r["data"]["fit"]["rv"]["file"] == "../../configs/sys/data/rv.csv"
    saved = yaml.safe_load((tmp_path / r["saved"]).read_text())
    assert saved["fit"]["lc"][0]["file"] == "../../configs/sys/data/lc1.csv"
    # saved again from the draft itself: nothing moves
    r2 = gs.save_draft(
        tmp_path, "sysfit", (tmp_path / r["saved"]).read_text(), True, "jobs/drafts/sysfit.yaml"
    )
    assert r2["rebased"] == []
    # frozen for a job (jobs/submitted/): every path still points at the same file
    dest = tmp_path / "jobs" / "submitted" / "001_x.yaml"
    dest.parent.mkdir(parents=True, exist_ok=True)
    gs._freeze_config(src / "sys_fit.yaml", dest)
    frozen = yaml.safe_load(dest.read_text().split("\n", 1)[1])
    for got, want in (
        (frozen["fit"]["rv"]["file"], src / "data/rv.csv"),
        (frozen["fit"]["inclination"][0]["file"], src / "data/inc.csv"),
        (frozen["fit"]["lc"][0]["file"], src / "data/lc1.csv"),
    ):
        assert (dest.parent / got).resolve() == want.resolve()


# --------------------------------------------------------------------------- quick check
def test_quickcheck_resolves_a_relative_stellar_table(tmp_path):
    (tmp_path / "stars.csv").write_text(
        "mass_msun,radius_rsun,apsidal_k2,rg2\n0.5,0.5,0.01,0.1\n2,2.0,0.01,0.1\n"
    )
    data = yaml.safe_load(
        CFG.replace("stars: {m1: {R: 1.2}, m2: {R: 1.0}, m3: {R: 1.1}, m4: {R: 0.9}}\n", "")
    )
    data["stellar_model"] = {"kind": "table", "file": "stars.csv"}
    r = assist.quickcheck(data, tmp_path)
    assert r["ok"], r
    # the radii come from the table next to the configuration (before 0.36 the table was not
    # found and the stars had no size: eclipse limit 90 deg)
    pair = next(i for i in r["items"] if i["k"].startswith("pair A"))
    assert "eclipse limit 90.0°" not in pair["v"] and pair["tone"] == "good"


def test_quickcheck_with_unknown_stability_ratios(monkeypatch):
    import swingersim.systemview as sv

    fake = {
        "triples": [{"label": "A in out", "ma01": None, "v22": 1.7}],
        "pairs": [],
        "orbits": [],
        "mutual": [],
        "hierarchy": "x",
        "stable": None,
    }
    monkeypatch.setattr(sv, "system_model", lambda cfg, n_points=8: fake)
    r = assist.quickcheck(yaml.safe_load(CFG))
    it = r["items"][0]
    assert it["v"] == "MA01 unknown · V22 1.70" and it["tone"] == "good"


def test_methods_text_names_every_fitted_kind_of_data(tmp_path):
    cfg = tmp_path / "f.yaml"
    d = yaml.safe_load(CFG)
    d["fit"] = {
        "minima": "m.csv",
        "rv": {"file": "rv.csv"},
        "lc": [{"file": "lc.csv"}],
        "inclination": [{"file": "i.csv", "pair": "A"}],
        "stage2": {"free": ["orbits.A.P"], "bjd_t0": 2459000.0},
    }
    cfg.write_text(yaml.safe_dump(d))
    m = assist.methods_text(cfg)
    t = m["text"]
    assert "radial velocities" in t and "light curves" in t and "inclinations" in t
    assert "jointly" in t and any("Bloemen" in r for r in m["references"])


def test_help_texts():
    h = assist.help_texts()
    assert "auto chooses" not in h["config"]["stars.*.shape"]["what"]
    assert "30 °/century" in h["config"]["physics.gr"]["what"]
    app = APP.read_text()
    for k in h["charts"]:  # every chart help is shown somewhere
        assert f'"{k}"' in app, k
    assert "fit" not in h["tabs"]


# --------------------------------------------------------------------------- server
@pytest.fixture
def root(tmp_path):
    out = tmp_path / "output"
    out.mkdir()
    (tmp_path / "configs").mkdir()
    (tmp_path / "configs" / "sys.yaml").write_text(CFG)
    return tmp_path


def test_table_cache_numpy_eviction_and_truncation(root, monkeypatch):
    p = root / "output" / "tab.csv"
    p.write_text("a,b,c\n1,x,True\n2,y,False\n,z,True\n")
    t = gs.table(p)
    assert t["a"].dtype == float and np.isnan(t["a"][2]) and t["b"].dtype == object
    info = {c["name"]: c for c in gs.output_info(root, "output/tab.csv")["columns"]}
    assert info["a"]["type"] == "number" and info["b"]["categories"] == {"x": 1, "y": 1, "z": 1}
    assert info["c"]["type"] == "category"
    monkeypatch.setattr(gs, "MAX_ROWS", 2)
    c = gs.columns(root, "output/tab.csv", ["a"])
    assert c["a"] == [1.0, 2.0] and c["_truncated"] and c["_n_rows"] == 3
    # a new version of the file replaces the old one in the cache
    time.sleep(0.01)
    p.write_text("a\n5\n")
    assert list(gs.table(p)["a"]) == [5.0]
    assert sum(1 for k in gs._TABLE_CACHE if k[0] == str(p)) == 1


def test_classify_parses_json_once(root, monkeypatch):
    p = root / "output" / "emu.json"
    p.write_text(json.dumps({"inputs": [], "targets": [], "metrics": {}}))
    assert gs.classify(p) == "emulator"
    monkeypatch.setattr(gs.json, "loads", lambda s: (_ for _ in ()).throw(AssertionError("parsed")))
    assert gs.classify(p) == "emulator"


def test_log_tail_from_the_end(tmp_path):
    p = tmp_path / "x.log"
    p.write_text("".join(f"line {k}\n" for k in range(200_000)))
    assert gs._tail(p, 3) == "line 199997\nline 199998\nline 199999"
    assert gs._tail(p, 3, block=7).splitlines()[-1] == "line 199999"


def test_parallel_submissions_get_distinct_numbers(root):
    gs.save_draft(root, "sys", CFG, True)
    ids, errors = [], []

    def go():
        try:
            ids.extend(gs.submit(root, ["sys"], "run")["jobs"])
        except Exception as exc:  # noqa: BLE001
            errors.append(exc)

    th = [threading.Thread(target=go) for _ in range(6)]
    for x in th:
        x.start()
    for x in th:
        x.join()
    assert not errors and len(ids) == 6 and len({i.split("_")[0] for i in ids}) == 6


def test_removed_job_time_is_utc(root):
    gs.save_draft(root, "sys", CFG, True)
    jid = gs.submit(root, ["sys"], "run")["jobs"][0]
    gs.job_action(root, [jid], "remove")
    st = json.loads((root / "jobs" / "cancel" / f"{jid}.status.json").read_text())
    assert st["finished"].endswith("+00:00")


def test_start_worker_does_not_hold_the_request(root, monkeypatch):
    class P:
        pid = 999999

        def poll(self):
            return None

    monkeypatch.setattr(gs.subprocess, "Popen", lambda *a, **k: P())
    t0 = time.monotonic()
    r = gs.start_worker(root)
    assert time.monotonic() - t0 < gs.START_WAIT_S + 1.0 and r["started"]


def test_full_digest_cached_until_the_output_changes(root, monkeypatch):
    lc = root / "output" / "lc.csv"
    lc.write_text("t_day,flux,flux_m1\n0,1,0.5\n1,0.9,0.4\n")
    d = root / "jobs" / "done"
    d.mkdir(parents=True)
    (d / "001_lc.yaml").write_text(
        yaml.safe_dump(
            {"command": "lightcurve", "config": "configs/sys.yaml", "output": "output/lc.csv"}
        )
    )
    (d / "001_lc.status.json").write_text(json.dumps({"state": "done", "finished": "x"}))
    calls = []
    orig = gs._digest

    def spy(*a, **k):
        calls.append(k.get("full"))
        return orig(*a, **k)

    monkeypatch.setattr(gs, "_digest", spy)
    a = gs.result(root, "001_lc")
    b = gs.result(root, "001_lc")
    assert calls.count(True) == 1 and a["data"] == b["data"]
    time.sleep(0.01)
    lc.write_text("t_day,flux,flux_m1\n0,1,0.5\n1,0.8,0.3\n2,1,0.5\n")
    assert gs.result(root, "001_lc")["data"]["n"] == 3 and calls.count(True) == 2
    # a new swingersim version rebuilds the cards
    import swingersim

    s1 = gs.card(root, "001_lc")["_stamp"]
    monkeypatch.setattr(swingersim, "__version__", "9.9.9")
    assert gs.card(root, "001_lc")["_stamp"] != s1


def test_timings_csv_read_once(root, monkeypatch):
    rows = ["binary,kind,epoch,t_obs_yr,o_c_d,P_d,T0_yr,impact,bjd_obs"]
    for e in range(10):
        rows.append(f"A,primary,{e},{e * 2 / 365.25:.9f},0.0,2.0,0.0,0.1,{2459000 + 2 * e:.6f}")
    (root / "output" / "t.csv").write_text("\n".join(rows) + "\n")
    (root / "m.csv").write_text("binary,kind,bjd_tdb,sigma_d,ref\nA,p,2459004.001,0.0005,x\n")
    (root / "configs" / "sys.yaml").write_text(CFG + "fit: {minima: ../m.csv}\n")
    d = root / "jobs" / "done"
    d.mkdir(parents=True)
    (d / "001_t.yaml").write_text(
        yaml.safe_dump(
            {"command": "timings", "config": "configs/sys.yaml", "output": "output/t.csv"}
        )
    )
    (d / "001_t.status.json").write_text(json.dumps({"state": "done", "finished": "x"}))
    n = []
    orig = D._read_csv
    monkeypatch.setattr(D, "_read_csv", lambda p: n.append(p) or orig(p))
    r = gs.result(root, "001_t")
    assert len(n) == 1 and r["data"]["measured"]["binaries"]["A"]["primary"]["rms"] > 1


def test_http_plotly_cache_and_tracebacks(root, monkeypatch, capsys):
    monkeypatch.setattr(gs, "scan_outputs", lambda *a: 1 / 0)
    srv = gs.make_server(root, "127.0.0.1", 0)
    th = threading.Thread(target=srv.serve_forever, daemon=True)
    th.start()
    try:
        base = f"http://127.0.0.1:{srv.server_address[1]}"
        if gs.plotly_js_path() is not None:
            r = urllib.request.urlopen(base + "/plotly.js")
            assert r.headers["Cache-Control"] == "max-age=86400" and r.read()
        assert urllib.request.urlopen(base + "/api/version").headers["Cache-Control"] == "no-store"
        with pytest.raises(urllib.error.HTTPError) as exc:
            urllib.request.urlopen(base + "/api/outputs")
        assert exc.value.code == 500
    finally:
        srv.shutdown()
        srv.server_close()
    assert "ZeroDivisionError" in capsys.readouterr().err


def test_page_has_no_mock_and_purges_plots():
    html = APP.read_text()
    assert "FitMock" not in html and 'data-page="fit"' not in html
    assert "function clearPlots" in html and "Plotly.purge" in html
    assert ".mono {" in html and "--border-strong:" in html
    assert "Math.min(...yy)" not in html  # the measured light curve: O(n), not O(n^2)


@pytest.mark.parametrize(
    ("head", "row", "kind"),
    [
        (
            "binary,kind,epoch,t_min_d,sigma_d,t_geom_d,o_c_d,o_c_geom_d,depth_mag,width_d,rms_mag",
            "A,primary,0,1.0,1e-4,1.0,1e-5,0,0.1,0.2,1e-3",
            "minima",
        ),
        ("binary,kind,t_mid_d,n,depth_mag", "A,primary,10.0,5,0.1", "seasons"),
        ("binary,kind,epoch,t_yr,o_c_d", "A,primary,0,0.1,1e-4", "timings"),
        ("t_day,flux,flux_m1", "0.5,0.9,0.4", "lightcurve"),
        ("time_yr,A_e,A_inc", "0.5,0.1,89.0", "table"),  # timeseries of a run
    ],
)
@pytest.mark.parametrize("comments", [True, False])
def test_commented_csv_headers_of_every_kind(tmp_path, head, row, kind, comments):
    """Minima, light curves and run time series start with '#' lines (0.36); old files do not."""
    (tmp_path / "output").mkdir()
    p = tmp_path / "output" / "x.csv"
    pre = "# SWINGERsim\n# units: days\n" if comments else ""
    p.write_text(pre + head + "\n" + row + "\n" + row + "\n")
    assert gs.classify(p) == kind
    info = gs.output_info(tmp_path, "output/x.csv")
    assert info["n_rows"] == 2 and [c["name"] for c in info["columns"]] == head.split(",")
    d = D._read_csv(p)
    assert list(d) == head.split(",") and len(d[head.split(",")[0]]) == 2
    if kind == "minima":
        assert D.minima_digest(p, "", None)["data"]["series"]["A"]["primary"]["t"] == [1.0, 1.0]
    if kind == "lightcurve":
        assert D.lightcurve_digest(p, None, {})["data"]["n"] == 2
