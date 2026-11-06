import json
import threading
import urllib.request

import pytest

from swingersim.gui import server as gs


@pytest.fixture
def root(tmp_path):
    out = tmp_path / "output"
    (out / "ens").mkdir(parents=True)
    (out / "ens" / "ensemble.yaml").write_text("n: 2\n")
    (out / "ens" / "summary.csv").write_text(
        "member,direction,fate,orbits.out.P\n0,forward,stable,100.5\n1,forward,collision,nan\n"
    )
    (out / "minima.csv").write_text(
        "binary,kind,epoch,t_min_d,sigma_d,t_geom_d,o_c_d,o_c_geom_d,depth_mag,width_d,rms_mag\n"
        "A,primary,0,1.0,1e-4,1.0,0.001,0.001,0.1,0.2,1e-3\n"
    )
    (out / "lc.csv").write_text("t_day,flux,flux_m1\n0,1,0.5\n1,0.9,0.4\n")
    (out / "emu.json").write_text(
        json.dumps({"inputs": ["a"], "targets": ["b"], "metrics": {"b": {"r2": 0.9}}})
    )
    cm = out / "map"
    cm.mkdir()
    (cm / "chaosmap.json").write_text(json.dumps({"x": "orbits.out.P", "config": {}}))
    (cm / "grid.jsonl").write_text('{"x": 1, "y": 2, "megno": 2.1}\n{"x": 2, "y": 2, "megno": 5}\n')
    (tmp_path / "configs").mkdir()
    (tmp_path / "configs" / "a.yaml").write_text(
        "topology: '2+2'\nmasses: {m1: 1, m2: 1, m3: 1, m4: 1}\n"
        "orbits: {A: {P: 3}, B: {P: 4}, out: {P: 1000}}\nrun: {t_max_yr: 10}\n"
    )
    jobs = tmp_path / "jobs" / "inbox"
    jobs.mkdir(parents=True)
    (jobs / "001_a.yaml.hold").write_text(
        "command: run\nconfig: configs/a.yaml\noutput: output/a.h5\n"
    )
    return tmp_path


def test_scan_and_classify(root):
    items = {i["path"]: i for i in gs.scan_outputs(root)}
    assert items["output/ens"]["kind"] == "ensemble"
    assert items["output/ens/summary.csv"]["parent"] == "output/ens"
    assert items["output/minima.csv"]["kind"] == "minima"
    assert items["output/lc.csv"]["kind"] == "lightcurve"
    assert items["output/emu.json"]["kind"] == "emulator"
    assert items["output/map"]["kind"] == "chaosmap"


def test_output_info_and_columns(root):
    info = gs.output_info(root, "output/ens")
    names = {c["name"]: c for c in info["columns"]}
    assert names["fate"]["type"] == "category" and names["fate"]["categories"]["stable"] == 1
    assert names["orbits.out.P"]["type"] == "number"
    assert info["presets"][0]["x"] == "fate"
    cols = gs.columns(root, "output/ens", ["orbits.out.P"])
    assert cols["orbits.out.P"][0] == 100.5 and cols["orbits.out.P"][1] is None  # NaN -> null
    m = gs.output_info(root, "output/map")
    assert m["n_rows"] == 2 and m["presets"][0]["type"] == "heat"
    assert gs.output_info(root, "output/emu.json")["emulator"]["metrics"]["b"]["r2"] == 0.9


def test_paths_stay_inside(root):
    with pytest.raises(gs.GuiError):
        gs.output_info(root, "../etc/passwd")
    with pytest.raises(gs.GuiError):
        gs.save_config(root, "../x.yaml", "a: 1\n", False)


def test_validate_and_save(root):
    good = gs.load_config(root, "configs/a.yaml")["data"]
    assert gs.validate_config(good)["ok"]
    bad = dict(good, orbits={"A": {"P": 3, "e": 1.5}, "B": {"P": 4}, "out": {"P": 1000}})
    r = gs.validate_config(bad)
    assert not r["ok"] and any(e["loc"] == "orbits.A.e" for e in r["errors"])
    gs.save_config(root, "configs/b.yaml", r["yaml"], False)
    with pytest.raises(gs.GuiError):
        gs.save_config(root, "configs/b.yaml", r["yaml"], False)  # no silent overwrite
    assert (root / "configs" / "b.yaml").exists()


def test_jobs(root):
    j = gs.jobs(root)
    item = j["folders"]["inbox"][0]
    assert item["held"] and item["command"] == "run"


def test_http_roundtrip(root):
    srv = gs.make_server(root, "127.0.0.1", 0)
    th = threading.Thread(target=srv.serve_forever, daemon=True)
    th.start()
    try:
        base = f"http://127.0.0.1:{srv.server_address[1]}"
        page = urllib.request.urlopen(base + "/").read().decode()
        assert "swingersim" in page and "Configuration" in page
        schema = json.loads(urllib.request.urlopen(base + "/api/schema").read())
        assert "masses" in schema["properties"]
        req = urllib.request.Request(
            base + "/api/validate",
            data=json.dumps({"yaml": "topology: '2+2'\n"}).encode(),
            headers={"Content-Type": "application/json"},
        )
        r = json.loads(urllib.request.urlopen(req).read())
        assert r["ok"] is False and r["errors"]
        with pytest.raises(urllib.error.HTTPError):
            urllib.request.urlopen(base + "/file?path=../../etc/passwd")
    finally:
        srv.shutdown()
        srv.server_close()


CFG = (
    "topology: '2+2'\nmasses: {m1: 1, m2: 1, m3: 1, m4: 1}\n"
    "orbits: {A: {P: 3}, B: {P: 4}, out: {P: 1000}}\nrun: {t_max_yr: 10}\n"
)


def test_drafts_submit_and_actions(root):
    gs.save_draft(root, "sys", CFG, False)
    with pytest.raises(gs.GuiError):
        gs.save_draft(root, "sys", CFG, False)  # no silent overwrite
    gs.save_draft(root, "sys", CFG.replace("t_max_yr: 10", "t_max_yr: 20"), True)
    with pytest.raises(gs.GuiError):
        gs.save_draft(root, "../x", CFG, False)
    assert [d["name"] for d in gs.list_drafts(root)] == ["sys"]
    r = gs.submit(root, ["sys"])
    jid = r["jobs"][0]
    assert jid.startswith("002_")  # 001 is taken by the job of the fixture
    frozen = root / "jobs" / "submitted" / f"{jid}.yaml"
    assert "t_max_yr: 20" in frozen.read_text()
    assert (root / "jobs" / "inbox" / f"{jid}.yaml.hold").exists()
    # editing the draft afterwards does not change the submitted job
    gs.save_draft(root, "sys", CFG, True)
    assert "t_max_yr: 20" in frozen.read_text()
    gs.job_action(root, [jid], "queue")
    assert (root / "jobs" / "inbox" / f"{jid}.yaml").exists()
    gs.job_action(root, [jid], "hold")
    gs.job_action(root, [jid], "remove")
    assert (root / "jobs" / "cancel" / f"{jid}.yaml").exists()
    gs.job_action(root, [jid], "again")  # back to the inbox, the old status kept aside
    assert (root / "jobs" / "inbox" / f"{jid}.yaml.hold").exists()
    assert (root / "jobs" / "cancel" / ".previous" / f"{jid}.status.json").exists()
    with pytest.raises(gs.GuiError):
        gs.job_action(root, [jid], "pause")  # not running
    gs.delete_draft(root, "sys")
    assert not gs.list_drafts(root)


def test_submit_job_for_any_command(root):
    r = gs.submit_job(
        root,
        {
            "command": "lightcurve",
            "config": "configs/a.yaml",
            "output": "output/{id}/lc.csv",
            "args": {"days": 2},
        },
        name="lc",
    )
    jid = r["jobs"][0]
    text = (root / "jobs" / "inbox" / f"{jid}.yaml.hold").read_text()
    assert f"output/{jid}/lc.csv" in text
    with pytest.raises(gs.GuiError):
        gs.submit_job(root, {"command": "run", "config": "../a.yaml"})
    with pytest.raises(gs.GuiError):
        gs.submit_job(root, {"command": "gui"})


def test_job_products_and_partial_series(root):
    d = root / "jobs" / "running"
    d.mkdir(parents=True)
    (d / "009_x.yaml").write_text("command: run\nconfig: configs/a.yaml\noutput: output/x/run.h5\n")
    (root / "output" / "x").mkdir()
    (root / "output" / "x" / "run.h5.partial.json").write_text(
        json.dumps(
            {
                "t_end_yr": 100.0,
                "time_yr": [0, 1, 2],
                "energy_error": [0, 1e-12, 2e-12],
                "orbits": {"A": {"e": [0.1, 0.2, 0.3], "inc": [0.0, 0.1, 0.2]}},
                "mutual_inclination": {"out-A": [0.5, 0.5, 0.5]},
            }
        )
    )
    p = gs.job_products(root, "009_x")
    assert p["partial"] == "output/x/run.h5.partial.json"
    idx = gs.h5_index(root, p["partial"])
    assert idx["partial"] and idx["n"] == 3 and idx["t_yr"] == 2
    s = gs.h5_series(root, p["partial"], ["orbits/A/e", "orbits/A/inc"])
    assert s["orbits/A/e"][-1] == 0.3
    assert abs(s["orbits/A/inc"][1] - 5.729577951308232) < 1e-9  # radians -> deg


def test_ensemble_presets_follow_the_columns(root):
    (root / "output" / "ens" / "summary.csv").write_text(
        "member,direction,fate,orbits.out.P,A_e_max,imut_out_A_0,imut_out_A_max,A_eclipsing0\n"
        "0,forward,intact,100,0.2,30,40,True\n1,forward,collision,200,0.9,80,120,False\n"
    )
    titles = [p["title"] for p in gs.output_info(root, "output/ens")["presets"]]
    assert "max e of A" in titles and "i_mut out_A: start vs max" in titles
    assert "A eclipsing now" in titles and "Roche-lobe filling of m1" not in titles


def test_commands_catalog(root):
    from swingersim.cli import command_catalog

    cat = command_catalog()
    assert "lightcurve" in cat and "worker" in cat
    days = [a for a in cat["lightcurve"]["arguments"] if a["dest"] == "days"]
    assert days and days[0]["kind"] in ("float", "int")


def test_run_only_the_selected_jobs(root):
    """A run from the GUI must not take other jobs that happen to be queued."""
    inbox = root / "jobs" / "inbox"
    (inbox / "002_b.yaml").write_text("command: run\nconfig: configs/a.yaml\noutput: output/b.h5\n")
    (inbox / "003_c.yaml").write_text("command: run\nconfig: configs/a.yaml\noutput: output/c.h5\n")
    gs.job_action(root, ["001_a"], "queue")
    held = gs.hold_others(root, ["001_a"])
    assert held == ["002_b", "003_c"]
    assert (inbox / "001_a.yaml").exists()
    assert (inbox / "002_b.yaml.hold").exists() and (inbox / "003_c.yaml.hold").exists()


def test_gui_finds_the_swingersim_folder(root):
    from swingersim.cli import find_root

    (root / "jobs" / "inbox").mkdir(parents=True, exist_ok=True)
    assert find_root(root / "jobs" / "inbox") == root.resolve()
    assert find_root(root) == root.resolve()


def test_gui_jobs_use_all_cores(root):
    gs.save_draft(root, "ens", CFG, False)
    r = gs.submit_job(
        root, {"command": "ensemble", "config": "configs/a.yaml", "output": "output/{id}"}
    )
    text = (root / "jobs" / "inbox" / f"{r['jobs'][0]}.yaml.hold").read_text()
    assert "processes: 0" in text
    r = gs.submit_job(
        root,
        {
            "command": "ensemble",
            "config": "configs/a.yaml",
            "output": "output/{id}",
            "args": {"processes": 4},
        },
    )
    assert "processes: 4" in (root / "jobs" / "inbox" / f"{r['jobs'][0]}.yaml.hold").read_text()
    r = gs.submit_job(
        root, {"command": "lightcurve", "config": "configs/a.yaml", "output": "output/{id}/lc.csv"}
    )
    assert "processes" not in (root / "jobs" / "inbox" / f"{r['jobs'][0]}.yaml.hold").read_text()


def test_stale_selection_is_skipped(root):
    """A job selected in the page that has meanwhile left the inbox does not block the others."""
    (root / "jobs" / "done").mkdir(parents=True, exist_ok=True)
    (root / "jobs" / "done" / "000_old.yaml").write_text("command: run\nconfig: configs/a.yaml\n")
    r = gs.job_action(root, ["000_old", "001_a"], "queue")
    assert r["jobs"] == ["001_a"] and r["skipped"][0]["id"] == "000_old"
    with pytest.raises(gs.GuiError):
        gs.job_action(root, ["000_old"], "queue")


def test_draft_that_is_a_job_file(root):
    """A complete job file put into jobs/drafts (command, config, args) is queued as that job,
    not wrapped as the configuration of a 'run' (job 060 on VYVAR, 2026-10-12)."""
    (root / "jobs" / "drafts").mkdir(parents=True, exist_ok=True)
    (root / "jobs" / "drafts" / "060_validate_x.yaml").write_text(
        "command: validate\nids: [disc-face-on-tau]\nnote: test\n"
    )
    (root / "jobs" / "drafts" / "063_lc.yaml").write_text(
        "command: lightcurve\nconfig: configs/a.yaml\noutput: output/x/lc.csv\nargs: {days: 2}\n"
    )
    kinds = {d["name"]: d["kind"] for d in gs.list_drafts(root)}
    assert kinds == {"060_validate_x": "validate", "063_lc": "lightcurve"}
    made = gs.submit(root, ["060_validate_x", "063_lc"])["jobs"]
    assert made[0].endswith("_validate_x") and not made[0].endswith("_060_validate_x")
    import yaml

    jobs = {
        j: yaml.safe_load((root / "jobs" / "inbox" / f"{j}.yaml.hold").read_text()) for j in made
    }
    assert jobs[made[0]]["command"] == "validate"
    assert jobs[made[0]]["ids"] == ["disc-face-on-tau"]
    assert jobs[made[1]]["command"] == "lightcurve"
    assert jobs[made[1]]["config"] == "configs/a.yaml"
