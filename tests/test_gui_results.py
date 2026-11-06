"""The Outputs page of swingersim gui: result cards and pages read from finished jobs (0.28)."""

import json
import subprocess
import sys

import pytest
import yaml

from swingersim.gui import digest as D
from swingersim.gui import server as gs

CFG = (
    "topology: '2+2'\nmasses: {m1: 1.2, m2: 1.0, m3: 1.1, m4: 0.9}\n"
    "orbits: {A: {P: 2.0, i: 89.0}, B: {P: 3.0, i: 88.5}, out: {P: 200.0, e: 0.2, i: 89.0}}\n"
    "stars: {m1: {R: 1.2}, m2: {R: 1.0}, m3: {R: 1.1}, m4: {R: 0.9}}\n"
    "run: {t_max_yr: 0.3, n_output: 121}\n"
)
ENS = (
    "topology: '2+2'\nmasses: {m1: 1.2, m2: 1.0, m3: 1.1, m4: 0.9}\n"
    "orbits: {A: {P: 2.0, i: 89.0}, B: {P: 3.0, i: 88.5}, out: {P: 200.0, e: 0.2, i: 89.0}}\n"
    "stars: {m1: {R: 1.2}, m2: {R: 1.0}, m3: {R: 1.1}, m4: {R: 0.9}}\n"
    "run: {t_max_yr: 2.0, n_output: 41, backend: secular}\n"
    "ensemble:\n  n: 8\n  seed: 3\n  directions: [forward, backward]\n"
    "  vary:\n    orbits.out.P: {loguniform: [100.0, 1000.0]}\n    orbits.out.i: isotropic\n"
    "  observed:\n    A_e_max: {max: 0.5}\n    out_e0: {value: 0.2, sigma: 0.05}\n"
)


def _q(root, *args):
    r = subprocess.run(
        [sys.executable, "-m", "swingersim", *args],
        cwd=root,
        capture_output=True,
        text=True,
        check=False,
    )
    assert r.returncode == 0, r.stderr[-2000:]
    return r.stdout


def _done(root, jid, job, log=""):
    d = root / "jobs" / "done"
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{jid}.yaml").write_text(yaml.safe_dump(job))
    (d / f"{jid}.status.json").write_text(
        json.dumps(
            {"id": jid, "state": "done", "finished": f"2026-10-13T10:{jid[:2]}:00", "elapsed_s": 3}
        )
    )
    (d / f"{jid}.log").write_text(log)


@pytest.fixture(scope="module")
def root(tmp_path_factory):
    root = tmp_path_factory.mktemp("q")
    (root / "configs").mkdir()
    (root / "configs" / "sys.yaml").write_text(CFG)
    (root / "configs" / "ens.yaml").write_text(ENS)
    (root / "output").mkdir()
    log = _q(root, "run", "configs/sys.yaml", "-o", "output/sys.h5")
    _done(
        root,
        "001_sys",
        {"command": "run", "config": "configs/sys.yaml", "output": "output/sys.h5"},
        log,
    )
    log = _q(
        root,
        "lightcurve",
        "configs/sys.yaml",
        "--days",
        "4",
        "--cadence-min",
        "20",
        "-o",
        "output/lc.csv",
    )
    _done(
        root,
        "002_sys_lightcurve",
        {
            "command": "lightcurve",
            "config": "configs/sys.yaml",
            "output": "output/lc.csv",
            "args": {"days": 4},
        },
        log,
    )
    log = _q(root, "timings", "configs/sys.yaml", "--years", "0.3", "-o", "output/t.csv")
    _done(
        root,
        "003_sys_timings",
        {
            "command": "timings",
            "config": "configs/sys.yaml",
            "output": "output/t.csv",
            "args": {"years": 0.3},
        },
        log,
    )
    log = _q(root, "ensemble", "configs/ens.yaml", "-o", "output/ens")
    _done(
        root,
        "004_ens",
        {"command": "ensemble", "config": "configs/ens.yaml", "output": "output/ens"},
        log,
    )
    (root / "output" / "m.csv").write_text(
        "binary,kind,epoch,t_min_d,sigma_d,t_geom_d,o_c_d,o_c_geom_d,depth_mag,width_d,rms_mag\n"
        "A,primary,0,1.0,1e-4,1.0,1e-5,0,0.1,0.2,1e-3\nA,secondary,0.5,2.0,2e-4,2.0,-1e-5,0,0.05,0.2,1e-3\n"
    )
    _done(
        root,
        "005_sys_minima",
        {"command": "minima", "config": "configs/sys.yaml", "output": "output/m.csv"},
        "  A primary  : 1 minima; shape A = 0.1 mag, D = 1.2 h, Gamma = 0.8, C = +0.1, rms 1e-03 mag; "
        "minimum - conjunction +2.0 s (scatter 3.0 s), median sigma 8.6 s\n"
        "  A: ephemeris of the minima T0 = 1.0 d, P = 2.0 d\n  skipped: 0 blended, 1 truncated eclipses\n",
    )
    _done(
        root,
        "006_validate",
        {"command": "validate"},
        "| a | x | PASS |\n| b | y | EXPLAINED |\n  EXPLAINED b: reason\n",
    )
    return root


def test_cards_grouped_by_configuration(root):
    res = gs.results(root)
    groups = {g["name"]: g["n"] for g in res["groups"]}
    assert groups["sys"] == 4 and groups["ens"] == 1 and groups["validate"] == 1
    cards = {c["id"]: gs.card(root, c["id"]) for c in res["cards"]}
    run = cards["001_sys"]
    assert run["kind"] == "run" and any(f["k"] == "Hierarchy changes" for f in run["facts"])
    assert "version" not in json.dumps(run).lower()  # no swingersim version on the pages
    assert cards["002_sys_lightcurve"]["kind"] == "lightcurve"
    assert cards["006_validate"]["facts"][0]["v"] == "1"
    # the second listing comes from the cache
    assert not any(c.get("pending") for c in gs.results(root)["cards"])


def test_run_page(root):
    d = gs.result(root, "001_sys")
    data = d["data"]
    assert data["tree"]["binaries"] == ["A", "B"] and data["tree"]["root"] == "out"
    assert {e["binary"] for e in data["eclipses"]} == {"A", "B"}
    assert all(0 < e["limit_now"] < 90 for e in data["eclipses"])
    assert len(data["view"]["t"]) == 121 and set(data["view"]["elements"]) == {"A", "B", "out"}
    assert [s["pair"] for s in data["secular"]] == ["A in out", "B in out"]
    assert {r["id"] for r in d["related"]} == {
        "002_sys_lightcurve",
        "003_sys_timings",
        "005_sys_minima",
    }


def test_lightcurve_and_timings_pages(root):
    lc = gs.result(root, "002_sys_lightcurve")["data"]
    assert [b["name"] for b in lc["binaries"]] == ["A", "B"]
    assert abs(sum(s["share"] for s in lc["star_rows"]) - 1) < 0.05
    assert lc["binaries"][0]["P_d"] == 2.0 and lc["binaries"][0]["span_periods"] > 1.9
    t = gs.result(root, "003_sys_timings")["data"]
    a = next(b for b in t["binaries"] if b["binary"] == "A")
    assert a["primary_star"] == "m1" and a["outer"] == "out" and a["T0_yr"] is not None
    p = t["series"]["A"]["primary"]
    assert abs(a["semi_amplitude"] - (max(p["oc"]) - min(p["oc"])) / 2) < 1e-6


def test_minima_page_reads_the_log(root):
    d = gs.result(root, "005_sys_minima")
    sh = d["data"]["shapes"]["A:primary"]
    assert sh["A_mag"] == 0.1 and sh["offset_s"] == 2.0
    assert d["data"]["ephemeris"]["A"]["P_d"] == 2.0
    assert d["data"]["series"]["A"]["primary"]["oc"][0] == pytest.approx(0.864)


def test_ensemble_page_and_reweighting(root):
    d = gs.result(root, "004_ens")["data"]
    assert d["n_members"] == 8 and set(d["directions"]) == {"forward", "backward"}
    assert {o["col"] for o in d["observations"]} == {"A_e_max", "out_e0"}
    assert d["conditioned"] and 0 < d["n_eff"] <= 8
    p = next(h for h in d["params"] if h["col"] == "orbits.out.P")
    assert p["log"] and p["drawn"]  # log-uniform draws are shown on a log axis
    assert not next(h for h in d["params"] if h["col"] == "orbits.out.i")["log"]
    none = gs.result(root, "004_ens", use=[])["data"]
    assert not none["conditioned"]
    hard = gs.result(root, "004_ens", use=["A_e_max"])["data"]
    assert hard["n_eff"] == pytest.approx(
        sum(
            1
            for r in D._rows(root / "output" / "ens")
            if r["direction"] == "forward" and r["A_e_max"] <= 0.5
        )
    )


def test_member_becomes_a_draft(root):
    r = gs.save_member_draft(root, "004_ens", 3)
    text = (root / r["saved"]).read_text()
    data = yaml.safe_load(text)
    assert "ensemble" not in data and "null" not in text  # unset values are dropped at every level
    row = next(
        x
        for x in D._rows(root / "output" / "ens")
        if x["member"] == 3 and x["direction"] == "forward"
    )
    assert data["orbits"]["out"]["P"] == pytest.approx(row["orbits.out.P"])
    assert gs.validate_config(data)["ok"]


def test_compute_jobs_carry_a_name_suffix(root):
    gs.save_draft(root, "mysys", CFG, overwrite=True)
    r = gs.submit(root, ["mysys"], "lightcurve", {"days": 2}, suffix="lightcurve")
    jid = r["jobs"][0]
    assert jid.endswith("_mysys_lightcurve")
    job = yaml.safe_load(
        (root / "jobs" / "inbox" / f"{jid}.yaml.hold").read_text().split("\n", 1)[1]
    )
    assert job["output"] == f"output/{jid}/lightcurve.csv"


def test_format_uses_superscripts():
    assert D.fmt(4.7e-10, 2) == "4.7·10⁻¹⁰"
    assert D.fmt(1234.5) == "1234"


def test_page_has_tools_and_no_commands_tab():
    html = (gs.APP).read_text()
    assert 'data-page="tools"' in html and 'data-page="commands"' not in html
    assert "What to compute" in html and "made with swingersim" not in html
