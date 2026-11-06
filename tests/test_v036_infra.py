"""Robustness of long jobs and of the configuration (worker, ensembles, chaos
maps, scattering, validation, command line, emulator, SBI training sets, stderr filter)."""

import json
import math
import os
import signal
import subprocess
import sys
import textwrap
import threading
import time
import warnings
from pathlib import Path

import numpy as np
import pytest
import yaml
from pydantic import ValidationError

from swingersim.config import Config, config_at
from swingersim.ensemble import run_member as _ORIGINAL_RUN_MEMBER

BASE = {
    "topology": "2+2",
    "masses": {"m1": 1.0, "m2": 1.0, "m3": 1.0, "m4": 1.0},
    "orbits": {
        "A": {"a": 0.5, "e": 0.1},
        "B": {"a": 0.4, "e": 0.1},
        "out": {"a": 8.0, "e": 0.2},
    },
    "physics": {"gr": False},
    "run": {"t_max_yr": 20.0, "n_output": 10, "backend": "secular"},
}


def _cfg(ens=None, **extra):
    data = json.loads(json.dumps(BASE))
    if ens is not None:
        data["ensemble"] = ens
    data.update(extra)
    return Config.model_validate(data)


# --------------------------------------------------------------------------- process pools


def _maybe_die(k):
    if k == 2:
        os._exit(139)  # like a segmentation fault or an out-of-memory kill
    time.sleep(0.02)
    return k


def test_run_tasks_survives_a_dying_worker():
    """multiprocessing.Pool waited forever for the task of a dead worker; now the task is
    reported lost and the others finish."""
    from swingersim.parallel import run_tasks

    t0 = time.time()
    # a Python SIGTERM handler of the parent (swingersim worker installs one) is inherited by the
    # forked pool workers: the pool's clean-up after the break hung on them before
    old = signal.signal(signal.SIGTERM, lambda *_: None)
    try:
        out = list(run_tasks(_maybe_die, range(8), 3, lambda t, msg: ("lost", t, msg)))
    finally:
        signal.signal(signal.SIGTERM, old)
    assert time.time() - t0 < 60
    lost = [o for o in out if isinstance(o, tuple)]
    assert len(lost) == 1 and lost[0][1] == 2 and "died" in lost[0][2]
    assert sorted(o for o in out if not isinstance(o, tuple)) == [0, 1, 3, 4, 5, 6, 7]


def _crashing_member(task):
    if task["index"] == 1:
        os._exit(139)
    return _ORIGINAL_RUN_MEMBER(task)


def test_ensemble_records_a_dead_member_and_finishes(tmp_path, monkeypatch):
    import swingersim.ensemble as ens

    monkeypatch.setattr(ens, "run_member", _crashing_member)
    cfg = _cfg({"n": 4, "vary": {"orbits.A.i": "isotropic"}})
    rows = ens.run_ensemble(cfg, tmp_path, processes=2)
    assert len(rows) == 4
    bad = [r for r in rows if r["fate"] == "error"]
    assert [r["member"] for r in bad] == [1] and "died" in bad[0]["error"]
    assert (tmp_path / "summary.csv").exists()


def test_resolve_processes(monkeypatch):
    from swingersim.ensemble import resolve_processes

    monkeypatch.delenv("SWINGERSIM_CORES", raising=False)
    assert resolve_processes(3) == 3
    assert resolve_processes(0) == (os.cpu_count() or 1)
    monkeypatch.setenv("SWINGERSIM_CORES", "2")
    assert resolve_processes(0) == 2 and resolve_processes(-1) == 2 and resolve_processes(5) == 5


# --------------------------------------------------------------------------- ensemble outputs


def test_rows_torn_line_duplicates_nan_and_shrinking_n(tmp_path):
    from swingersim.ensemble import load_rows, run_ensemble

    cfg = _cfg({"n": 3, "vary": {"orbits.A.i": "isotropic"}})
    run_ensemble(cfg, tmp_path)
    f = tmp_path / "rows.jsonl"
    text = f.read_text()
    assert "NaN" not in text and "Infinity" not in text  # null instead of bare NaN
    for line in text.splitlines():
        json.loads(line)  # strict JSON
    first = text.splitlines()[0]
    f.write_text(text + first + "\n" + '{"member": 7, "direction": "forw')  # dup + torn line
    rows = load_rows(tmp_path)
    assert sorted((r["member"], r["direction"]) for r in rows) == [
        (0, "forward"),
        (1, "forward"),
        (2, "forward"),
    ]
    assert any(isinstance(v, float) and math.isnan(v) for v in rows[0].values())  # null -> NaN
    cfg.ensemble.n = 4
    with pytest.warns(UserWarning, match="malformed"):
        rows = run_ensemble(cfg, tmp_path, resume=True)
    assert len(rows) == 4
    lines = f.read_text().splitlines()
    assert len(lines) == 4 and all(json.loads(x) for x in lines)  # torn line and dup removed
    cfg.ensemble.n = 2  # a smaller n: the rows of members 2 and 3 are not part of it
    rows = run_ensemble(cfg, tmp_path, resume=True)
    assert [r["member"] for r in rows] == [0, 1]
    assert [r["member"] for r in load_rows(tmp_path)] == [0, 1]
    meta = yaml.safe_load((tmp_path / "ensemble.yaml").read_text())
    assert meta["swingersim"] and meta["python"] and meta["started"] and meta["resumed"]


def test_write_summary_falls_back_to_csv(tmp_path):
    pytest.importorskip("pyarrow")
    from swingersim.ensemble import write_summary

    (tmp_path / "summary.parquet").write_text("stale")
    rows = [{"member": 0, "x": "text"}, {"member": 1, "x": 1.5}, {"member": 2, "x": [1, 2]}]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        write_summary(rows, tmp_path)
    assert (tmp_path / "summary.csv").exists()
    assert not (tmp_path / "summary.parquet").exists() or (
        (tmp_path / "summary.parquet").read_bytes()[:4] == b"PAR1"
    )


def test_draw_members_shows_the_real_error():
    from swingersim.ensemble import draw_members

    cfg = _cfg({"n": 1, "max_tries": 3, "vary": {"orbits.A.e": {"uniform": [1.2, 1.5]}}})
    with pytest.raises(RuntimeError, match="no acceptable draw.*ValidationError") as ei:
        draw_members(cfg)
    assert isinstance(ei.value.__cause__, ValidationError)


def test_vectorised_draws():
    from swingersim.config import Distribution
    from swingersim.ensemble import draw_values

    rng = np.random.default_rng(1)
    v = draw_values(Distribution(normal=(1.0, 1.0), min=0.5, max=1.2), rng, 5000)
    assert v.min() >= 0.5 and v.max() <= 1.2
    iso = draw_values("isotropic", rng, 20000)
    assert abs(np.mean(np.cos(np.radians(iso)))) < 0.03
    m = draw_values(Distribution(uniform=(80.0, 85.0), mirror=True), rng, 4000)
    assert 0.4 < np.mean(m > 90) < 0.6 and np.all((m <= 85) | (m >= 95))
    c = draw_values(Distribution(choice=[1.0, 2.0]), rng, 1000)
    assert set(c) == {1.0, 2.0}


# --------------------------------------------------------------------------- configuration

FIT = {
    "minima": "data/x.csv",
    "stage2": {"bjd_t0": 2450000.0, "free": ["orbits.A.P", "masses.m1"]},
}


def test_config_names_in_keys_and_lists_are_checked():
    with pytest.raises(ValidationError, match="not free parameters"):
        _cfg(fit={**FIT, "stage2": {**FIT["stage2"], "bounds": {"orbits.B.P": [1, 2]}}})
    with pytest.raises(ValidationError, match="not free parameters"):
        _cfg(fit={**FIT, "stage2": {**FIT["stage2"], "fd_steps": {"masses.m2": 1e-4}}})
    with pytest.raises(ValidationError, match="physics.activity"):
        _cfg(physics={"gr": False, "activity": {"m9": {"p_mod_yr": 10, "amplitude_s": 5}}})
    with pytest.raises(ValidationError, match="magnetic_braking"):
        _cfg(
            physics={
                "evolution": True,
                "tides": True,
                "rotation": True,
                "magnetic_braking": {"stars": ["m7"]},
            },
            stellar_model={"kind": "sse", "age_myr": 100.0},
        )
    with pytest.raises(ValidationError, match="fit.rv.stars"):
        _cfg(fit={**FIT, "rv": {"file": "rv.csv", "stars": ["m5"]}})
    with pytest.raises(ValidationError, match="fit.ephemerides"):
        _cfg(fit={**FIT, "ephemerides": {"C": {"T0": 1.0, "P": 2.0}}})
    with pytest.raises(ValidationError, match="fit.stage1.pairs"):
        _cfg(fit={**FIT, "stage1": {"pairs": {"inner": ["lite"]}}})
    with pytest.raises(ValidationError, match="fit.stage2.pairs"):
        _cfg(fit={**FIT, "stage2": {**FIT["stage2"], "pairs": ["X"]}})
    with pytest.raises(ValidationError, match="predict.pairs"):
        _cfg(fit={**FIT, "output": {"predict": {"pairs": ["Z"]}}})
    ok = _cfg(fit={**FIT, "ephemerides": {"A": {"T0": 1.0, "P": 2.0}}}, fit_epoch_bjd_tdb=2.4e6)
    assert ok.fit_epoch_bjd_tdb == 2.4e6


def test_config_ensemble_and_fit_rules():
    with pytest.raises(ValidationError, match="directions"):
        _cfg({"n": 1, "directions": ["forward", "forward"]})
    with pytest.raises(ValidationError, match="ma01_window"):
        _cfg({"n": 1, "ma01_window": [2.0, 1.0]})
    with pytest.raises(ValidationError, match="should be stars"):
        _cfg({"n": 1, "vary": {"stars.m1.shape": 1.0}})
    assert _cfg({"n": 1, "vary": {"stars.m1.P_rot": 2.0, "orbits.A.contact_fillout": 0.2}})
    lc = [{"file": "lc.csv", "red_noise": True}]
    with pytest.raises(ValidationError, match="red_noise"):
        _cfg(fit={**FIT, "lc": lc})
    assert _cfg(fit={**FIT, "lc": lc, "error_scale": "global"})
    assert _cfg(fit={**FIT, "lc": lc, "stage2": {**FIT["stage2"], "clip_sigma": 5.0}})
    with pytest.warns(UserWarning, match="walkers"):
        _cfg(fit={**FIT, "stage2": {**FIT["stage2"], "walkers": 40}})
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        c = _cfg(fit=FIT)
        Config.model_validate(c.model_dump())  # a dump does not warn


def test_config_at_keeps_contact_fillout():
    from swingersim.backends.secular import run_secular

    cfg = _cfg()
    cfg.orbits["A"].contact_fillout = 0.3
    res = run_secular(cfg)
    assert config_at(cfg, res).orbits["A"].contact_fillout == 0.3


# --------------------------------------------------------------------------- command line


def _write(tmp_path, data, name="c.yaml"):
    p = tmp_path / name
    p.write_text(yaml.safe_dump(data))
    return p


def test_cli_overrides_are_validated(tmp_path):
    from swingersim.cli import main

    p = _write(tmp_path, {**BASE, "ensemble": {"n": 2}})
    with pytest.raises(ValidationError):
        main(["run", str(p), "--t-max-yr", "0", "--no-html"])
    with pytest.raises(ValidationError):
        main(["ensemble", str(p), "-o", str(tmp_path / "e"), "--n", "-5", "--no-html"])


def test_cli_rejects_options_that_do_not_apply(tmp_path, capsys):
    from swingersim.cli import main

    p = _write(tmp_path, {**BASE, "fit": {"minima": "m.csv"}})
    assert main(["fit", str(p), "--bjd-t0", "2450000", "--no-html"]) == 2
    assert main(["fit", str(p), "--processes", "2", "--no-html"]) == 2  # stage 1
    e = _write(tmp_path, {**BASE, "ensemble": {"n": 2}}, "e.yaml")
    assert main(["ensemble", str(e), "-o", str(tmp_path / "o"), "--retry-errors"]) == 2
    assert main(["ensemble-report", str(tmp_path / "o"), "--x", "orbits.A.i", "--no-html"]) == 2
    out = capsys.readouterr().out
    assert "--anchor" in out and "stage 2" in out and "--resume" in out and "--plot" in out


def test_cli_ensemble_survives_bad_observations(tmp_path, capsys):
    from swingersim.cli import main

    e = _write(
        tmp_path,
        {
            **BASE,
            "ensemble": {
                "n": 2,
                "vary": {"orbits.A.i": "isotropic"},
                "observed": {"no_such_column": {"value": 1.0, "sigma": 0.1}},
            },
        },
    )
    assert main(["ensemble", str(e), "-o", str(tmp_path / "o"), "--no-html"]) == 0
    assert "no conditioning on ensemble.observed" in capsys.readouterr().out
    bad = _write(
        tmp_path,
        {**BASE, "ensemble": {"n": 2, "observed": {"orbits.A.i": {"min": 0.0}}}},
        "bad.yaml",
    )
    with pytest.raises(ValueError, match="not an ensemble.vary parameter"):
        main(["ensemble", str(bad), "-o", str(tmp_path / "o2"), "--no-html"])


def test_mass_transfer_summary_with_underscore_orbit_names():
    from swingersim.cli import _mass_transfer_groups

    system = {
        "mass_transfer_A_B_t_onset_yr": 5.0,
        "mass_transfer_A_B_donor": "m1",
        "mass_transfer_A_t_onset_yr": None,
        "other": 1,
    }
    g = _mass_transfer_groups(system, {"A": 0, "A_B": 0, "out": 0})
    assert g == {"A": {"t_onset_yr": None}, "A_B": {"t_onset_yr": 5.0, "donor": "m1"}}


# --------------------------------------------------------------------------- worker


def test_second_worker_refuses_to_start(tmp_path):
    from swingersim.worker import Worker, WorkerLocked

    w = Worker(tmp_path / "jobs", tmp_path)
    with pytest.raises(WorkerLocked, match="already running"):
        Worker(tmp_path / "jobs", tmp_path)
    w.close()
    Worker(tmp_path / "jobs", tmp_path).close()


def test_recover_stops_a_paused_job_before_requeueing(tmp_path):
    from swingersim.worker import Worker

    jobs = tmp_path / "jobs"
    w = Worker(jobs, tmp_path)
    proc = subprocess.Popen(
        [sys.executable, "-c", "import time  # swingersim job\ntime.sleep(600)"],
        start_new_session=True,
    )
    threading.Thread(target=proc.wait, daemon=True).start()  # reap it like init would
    os.kill(proc.pid, signal.SIGSTOP)
    (jobs / "running" / "005_x.yaml").write_text(yaml.safe_dump({"command": "validate"}))
    (jobs / "running" / "005_x.status.json").write_text(json.dumps({"pid": proc.pid}))
    (jobs / "pause" / "005_x").write_text("")
    t0 = time.time()
    assert w.recover() == ["005_x"]
    assert time.time() - t0 < 20
    assert proc.poll() is not None or not Path(f"/proc/{proc.pid}").exists()
    assert (jobs / "inbox" / "005_x.yaml").exists()
    assert not (jobs / "pause" / "005_x").exists()  # a requeued job must not start paused
    w.close()


def test_job_cores_of_fit_and_validate(tmp_path):
    from swingersim.worker import job_cores

    cfg = {**BASE, "fit": {**FIT, "stage2": {**FIT["stage2"], "processes": 0}}}
    (tmp_path / "f.yaml").write_text(yaml.safe_dump(cfg))
    cfg["fit"]["stage2"]["processes"] = 3
    (tmp_path / "f3.yaml").write_text(yaml.safe_dump(cfg))
    (tmp_path / "s1.yaml").write_text(yaml.safe_dump({**BASE, "fit": {"minima": "m.csv"}}))
    fit = {"command": "fit", "config": "f.yaml"}
    assert job_cores(fit, 12, tmp_path) == 12  # fit.stage2.processes 0 = all cores
    assert job_cores({**fit, "config": "f3.yaml"}, 12, tmp_path) == 3
    assert job_cores({**fit, "args": {"processes": 2}}, 12, tmp_path) == 2
    assert job_cores({**fit, "args": {"stage": 1}}, 12, tmp_path) == 1
    assert job_cores({**fit, "args": {"anchor": True}}, 12, tmp_path) == 1
    assert job_cores({"command": "fit", "config": "s1.yaml"}, 12, tmp_path) == 1
    assert job_cores({"command": "validate", "args": {"all": True}}, 12) == 12
    assert job_cores({"command": "validate", "ids": ["x"]}, 12) == 1


def test_worker_exports_the_job_cores(tmp_path):
    from swingersim.worker import Worker

    root = tmp_path / "q"
    (root / "configs").mkdir(parents=True)
    cfg = {**BASE, "run": {"t_max_yr": 2.0e5, "n_output": 10, "backend": "nbody"}}
    cfg["ensemble"] = {"n": 4, "vary": {"orbits.A.i": "isotropic"}}
    (root / "configs" / "slow.yaml").write_text(yaml.safe_dump(cfg))
    w = Worker(root / "jobs", root, poll_s=0.1, cores=2)
    job = {"command": "ensemble", "config": "configs/slow.yaml", "output": "output/e"}
    (root / "jobs" / "inbox" / "001_e.yaml").write_text(
        yaml.safe_dump({**job, "args": {"processes": 0}})
    )
    w.schedule()
    a = w.active["001_e"]
    try:
        import time

        # the forked child shows the parent's environment until it has exec'd: wait for it
        for _ in range(100):
            env = Path(f"/proc/{a.proc.pid}/environ").read_bytes().split(b"\0")
            if b"SWINGERSIM_CORES=2" in env:
                break
            time.sleep(0.05)
        assert b"SWINGERSIM_CORES=2" in env and b"OMP_NUM_THREADS=1" in env
        st = json.loads((root / "jobs" / "running" / "001_e.status.json").read_text())
        assert st["cores"] == 2
        assert not list((root / "jobs" / "running").glob(".*tmp"))
    finally:
        w.stop_all("test")
        w.close()
    assert w.start_job(root / "jobs" / "inbox" / "missing.yaml") is None  # moved away: skip


def test_end_kills_orphaned_pool_children(tmp_path):
    from swingersim.worker import Worker, _Active

    w = Worker(tmp_path / "jobs", tmp_path)
    code = textwrap.dedent(
        """
        import subprocess, sys  # swingersim
        p = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(600)"])
        print(p.pid, flush=True)
        """
    )
    proc = subprocess.Popen(
        [sys.executable, "-c", code], start_new_session=True, stdout=subprocess.PIPE, text=True
    )
    child = int(proc.stdout.readline())
    proc.wait()
    logf = open(tmp_path / "log.txt", "w")  # noqa: SIM115
    a = _Active("009_x", {"command": "run"}, proc, logf, None, {"id": "009_x"}, 1)
    w.active["009_x"] = a
    w._end(a)

    def gone(pid):
        try:
            return Path(f"/proc/{pid}/stat").read_text().split(")")[-1].split()[0] in "ZX"
        except OSError:
            return True

    assert gone(child)
    w.close()


# --------------------------------------------------------------------------- chaos maps


CB = {
    "topology": "custom",
    "hierarchy": {"out": ["inner", "m3"], "inner": ["m1", "m2"]},
    "masses": {"m1": 1.0, "m2": 1.0, "m3": 1.0e-9},
    "orbits": {"inner": {"a": 1.0}, "out": {"a": 3.0}},
    "stars": {"m1": {"R": 0.001}, "m2": {"R": 0.001}, "m3": {"R": 0.001}},
    "physics": {"gr": False},
    "run": {"t_max_yr": 10, "collisions": False},
}


def test_chaosmap_resume_checks_config_and_classes(tmp_path, monkeypatch):
    from swingersim import massloss
    from swingersim.chaosmap import chaos_map, megno_array, summarize_map

    calls = []
    orig = massloss.resolve_masses
    monkeypatch.setattr(massloss, "resolve_masses", lambda c: calls.append(1) or orig(c))
    cfg = Config.model_validate(CB)
    xs = np.array([1.6, 3.0])
    rows = chaos_map(cfg, tmp_path, "orbits.out.a", xs, t_max_yr=30.0)
    assert calls  # megno_point resolves the masses as run_nbody does
    other = Config.model_validate({**CB, "masses": {**CB["masses"], "m2": 0.9}})
    with pytest.raises(ValueError, match="configuration changed"):
        chaos_map(other, tmp_path, "orbits.out.a", xs, t_max_yr=30.0, resume=True)
    with pytest.raises(ValueError, match="escape_factor"):
        chaos_map(cfg, tmp_path, "orbits.out.a", xs, t_max_yr=30.0, resume=True, escape_factor=5)
    f = tmp_path / "grid.jsonl"
    lines = f.read_text().splitlines()
    f.write_text(lines[0] + "\n" + lines[1][:10])  # torn last line: that point runs again
    with pytest.warns(UserWarning, match="malformed"):
        rows = chaos_map(cfg, tmp_path, "orbits.out.a", xs, t_max_yr=30.0, resume=True)
    assert len(rows) == 2 and all(json.loads(x) for x in f.read_text().splitlines())
    meta = json.loads((tmp_path / "chaosmap.json").read_text())
    assert meta["swingersim"] and meta["rebound"]
    fake = [
        {"i": 0, "j": 0, "megno": 2.0, "status": "ok"},
        {"i": 1, "j": 0, "megno": 2.1, "status": "wall_time"},
        {"i": 2, "j": 0, "megno": math.nan, "status": "error: boom"},
        {"i": 3, "j": 0, "megno": math.nan, "status": "invalid: x"},
        {"i": 4, "j": 0, "megno": math.nan, "status": "escape"},
        {"i": 5, "j": 0, "megno": 9.0, "status": "ok"},
    ]
    a = megno_array(fake)[:, 0]
    assert a[0] == 2.0 and np.isnan(a[1]) and np.isnan(a[2]) and np.isnan(a[3])
    assert np.isinf(a[4]) and a[5] == 9.0
    meta = {"x": "x", "y": None, "xvals": list(range(6)), "t_max_yr": 1.0}
    text = summarize_map(meta, fake)
    assert "regular (MEGNO < 4): 1" in text and "chaotic (MEGNO >= 4): 1" in text
    assert "escape or collision: 1" in text and "computing-time limit" in text
    assert "failed integration: 2" in text


def test_chaosmap_tasks_are_small():
    from swingersim.chaosmap import _init, _task

    _init(Config.model_validate(CB).model_dump(), "orbits.out.a", None, 1.0, 10.0, 20.0)
    r = _task((0, 0, 3.0, math.nan))
    assert r["status"] == "ok" and r["x"] == 3.0


# --------------------------------------------------------------------------- scattering


def test_scatter_summary_classes_and_linear_cross_section():
    from swingersim.scatter import ScatterConfig, summarize_scatter

    def cfg(dist):
        return ScatterConfig(
            name="t",
            scatter={
                "target": {"masses": [1.0, 1.0], "a_au": 1.0},
                "projectile": {"masses": [1.0]},
                "v_inf_kms": 10.0,
                "b_max_au": 2.0,
                "b_distribution": dist,
            },
        )

    rows = [
        {"id": 0, "b_au": 0.5, "outcome": "flyby", "status": "resolved", "wall_s": 1.0},
        {"id": 1, "b_au": 1.5, "outcome": "flyby", "status": "resolved", "wall_s": 1.0},
        {"id": 2, "b_au": 1.0, "outcome": "exchange", "status": "resolved", "wall_s": 1.0},
        {"id": 3, "b_au": 1.0, "outcome": "bound (1,2),3", "status": "wall_time", "wall_s": 1.0},
        {"id": 4, "b_au": 1.0, "outcome": "error", "status": "error", "error": "x"},
    ]
    text = summarize_scatter(cfg("linear"), rows)
    assert "1 encounters ended with an error" in text
    assert "unfinished (wall_time): bound (1,2),3" in text
    # linear in b: weights 2 b / b_max = 0.5, 1.5, 1, 1 -> flyby (0.5 + 1.5) / 4 of pi b_max^2
    flyby = next(x for x in text.splitlines() if x.strip().startswith("flyby"))
    assert f"{math.pi * 4 * 0.5:.4g}" in flyby and " 0.500 " in flyby
    area = summarize_scatter(cfg("area"), rows)
    assert f"{math.pi * 4 * 0.5:.4g}" in next(
        x for x in area.splitlines() if x.strip().startswith("flyby")
    )


# --------------------------------------------------------------------------- validation


def test_failing_benchmark_is_an_error_and_others_run(monkeypatch, tmp_path, capsys):
    from swingersim import validation as v
    from swingersim.cli import main

    b = next(x for x in v.BENCHMARKS if x.id == "polytrope-k2-n1.5")

    def boom():
        raise ZeroDivisionError("bad")

    broken = v.Benchmark(**{**b.__dict__, "id": "broken", "compute": boom})
    monkeypatch.setattr(v, "BENCHMARKS", [broken, b])
    out = v.run_benchmarks()
    assert [o.benchmark.id for o in out] == ["broken", "polytrope-k2-n1.5"]
    assert out[0].error.startswith("ZeroDivisionError") and v._verdict(out[0]) == "ERROR"
    assert out[1].passed
    text = v.format_report(out, command="swingersim validate")
    assert "Seconds" in text and "swingersim validate" in text and "ERROR" in text
    with pytest.raises(ValueError, match="unknown benchmark"):
        v.run_benchmarks(["nope"])
    monkeypatch.chdir(tmp_path)
    assert main(["validate", "nope"]) == 2
    assert main(["validate", "broken", b.id]) == 1  # an error is a failure
    assert (tmp_path / "output/validation/report-broken-polytrope-k2-n1.5.html").exists()
    assert not (tmp_path / "output/validation/report.html").exists()
    assert "FAILED: broken" in capsys.readouterr().out


def test_interrupted_validation_still_writes_the_report(monkeypatch, tmp_path):
    from swingersim import validation as v
    from swingersim.cli import main

    b = next(x for x in v.BENCHMARKS if x.id == "polytrope-k2-n1.5")

    def stop():
        raise KeyboardInterrupt

    monkeypatch.setattr(
        v, "BENCHMARKS", [b, v.Benchmark(**{**b.__dict__, "id": "stop", "compute": stop})]
    )
    monkeypatch.chdir(tmp_path)
    with pytest.raises(KeyboardInterrupt):
        main(["validate", b.id, "stop"])
    assert (tmp_path / f"output/validation/report-{b.id}-stop.html").exists()


def test_nodal_formula_and_mist_loo_are_shared():
    from swingersim import validation as v

    assert not hasattr(v, "_secular_nodal_formula") and not hasattr(v, "mist_loo_radius")
    exp = {b.id: b.expected for b in v.BENCHMARKS}
    assert exp["nodal-precession"] == exp["secular-nodal-precession"]
    assert exp["nodal-precession"] == pytest.approx(v._nodal_precession_formula())


# --------------------------------------------------------------------------- emulator, sbi


def test_emulator_residual_in_network_space_and_prior_rules():
    from swingersim import emulator as emu

    rng = np.random.default_rng(0)
    x = rng.uniform(0.0, 3.0, 300)
    y = 10**x * (1 + 0.1 * rng.normal(size=300))  # a relative error of 10 %
    em = emu.train(
        x[:, None], y[:, None], np.zeros((300, 0)), ["x"], ["y"], log_targets=["y"],
        hidden=(16,), n_models=2, max_iter=300,
    )  # fmt: skip
    assert em.resid_net is not None and em.resid_net[0] < 0.2  # dex
    X = np.array([[0.0], [3.0]])
    mu = em.predict(X)["mean"][:, 0]
    s = np.sqrt(emu._resid_var(em, X, 0, mu))
    assert s[1] / s[0] == pytest.approx(mu[1] / mu[0], rel=1e-6)  # relative, not absolute
    cfg = _cfg({"n": 2, "vary": {"orbits.A.i": "isotropic"}, "require": ["hierarchical"]})
    with pytest.raises(ValueError, match="require"):
        emu.prior_sample(cfg, ["orbits.A.i"], 10)
    cfg = _cfg({"n": 2, "vary": {"orbits.A.i": "isotropic", "orbits.A.e": 0.2}})
    X = emu.prior_sample(cfg, ["orbits.A.i", "orbits.A.e"], 1000)
    assert X.shape == (1000, 2) and np.all(X[:, 1] == 0.2) and 0 <= X[:, 0].min()


def test_emulator_dataset_uses_index_runs():
    from swingersim import emulator as emu

    rows = [
        {"member": 0, "direction": "forward", "a": 1.0, "y": 2.0, "fate": "intact"},
        {"member": 0, "direction": "backward", "a": 1.0, "y": 4.0, "fate": "intact"},
    ]
    _X, Y, _C, keep = emu.dataset(rows, ["a"], ["backward:y"])
    assert keep == [0] and Y[0, 0] == 4.0


def test_sbi_training_set_skips_non_positive_logs_and_fingerprint(tmp_path):
    from swingersim.sbi_fit import training_fingerprint, training_set

    rows = [
        {"member": k, "direction": "forward", "p": float(k), "y": float(k), "fate": "intact"}
        for k in range(4)
    ]
    obs = {"y": {"value": 1.0, "sigma": 0.1}}
    with pytest.warns(UserWarning, match="1 members left out"):
        th, *_ = training_set(rows, ["p"], obs, ["p"])
    assert len(th) == 3
    (tmp_path / "ensemble.yaml").write_text(yaml.safe_dump({"fingerprint": "abc"}))
    f1 = training_fingerprint(tmp_path, 4, obs)
    assert f1 == training_fingerprint(tmp_path, 4, obs)
    assert f1 != training_fingerprint(tmp_path, 5, obs)
    assert f1 != training_fingerprint(tmp_path, 4, {**obs, "fate": {"equals": "intact"}})


# --------------------------------------------------------------------------- stderr filter


def test_stderr_filter_survives_a_closed_reader_and_keeps_crash_tracebacks(tmp_path):
    script = tmp_path / "s.py"
    script.write_text(
        textwrap.dedent(
            """
            import os, sys
            from swingersim import cnoise
            cnoise.install()
            print("ready", flush=True)
            sys.stdin.readline()
            for k in range(20000):  # far more than a pipe holds
                try:
                    os.write(2, b"x" * 100 + b"\\n")
                except OSError:
                    pass
            print("finished", flush=True)
            """
        )
    )
    env = {k: v for k, v in os.environ.items() if k != "SWINGERSIM_NO_STDERR_FILTER"}
    p = subprocess.Popen(
        [sys.executable, str(script)],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
        text=True,
    )
    assert p.stdout.readline().strip() == "ready"
    p.stderr.close()  # the reader of the log goes away
    p.stdin.write("go\n")
    p.stdin.flush()
    try:
        assert p.stdout.readline().strip() == "finished"  # blocked forever before
    finally:
        p.kill()
        p.wait()
    crash = tmp_path / "c.py"
    crash.write_text("import faulthandler\nfrom swingersim import cnoise\ncnoise.install()\n"
                     "faulthandler._sigsegv()\n")  # fmt: skip
    out = subprocess.run(
        [sys.executable, str(crash)],
        capture_output=True,
        text=True,
        env=env,
        timeout=60,
        check=False,
    )
    assert out.returncode != 0 and "Segmentation fault" in out.stderr
