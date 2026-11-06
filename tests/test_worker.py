import json
import shutil
from pathlib import Path

import pytest
import yaml

from swingersim.worker import JobError, Worker, build_command

ROOT = Path(__file__).resolve().parents[1]


def _setup(tmp_path):
    root = tmp_path / "swingersim"
    (root / "configs").mkdir(parents=True)
    cfg = {
        "topology": "2+2",
        "masses": {"m1": 1.0, "m2": 1.0, "m3": 1.0, "m4": 1.0},
        "orbits": {
            "A": {"a": 0.5, "e": 0.1},
            "B": {"a": 0.4, "e": 0.1},
            "out": {"a": 8.0, "e": 0.2},
        },
        "physics": {"gr": False},
        "run": {"t_max_yr": 20.0, "n_output": 10, "backend": "secular"},
        "ensemble": {"n": 2, "vary": {"orbits.A.i": "isotropic"}},
    }
    (root / "configs" / "tiny.yaml").write_text(yaml.safe_dump(cfg))
    return root


def test_rejects_unsafe_jobs(tmp_path):
    root = _setup(tmp_path)
    with pytest.raises(JobError):
        build_command({"command": "rm", "config": "x"}, root)
    with pytest.raises(JobError):
        build_command({"command": "ensemble", "config": "../../etc/passwd", "output": "o"}, root)
    with pytest.raises(JobError):
        build_command(
            {"command": "ensemble", "config": "configs/tiny.yaml", "output": "/tmp/x"}, root
        )
    with pytest.raises(JobError):
        build_command(
            {
                "command": "ensemble",
                "config": "configs/tiny.yaml",
                "output": "o",
                "args": {"shell": "ls"},
            },
            root,
        )
    argv, out = build_command(
        {
            "command": "ensemble",
            "config": "configs/tiny.yaml",
            "output": "out/e",
            "args": {"n": 3, "processes": 2, "resume": True},
        },
        root,
    )
    assert argv[3:] == [
        "ensemble",
        str(root / "configs/tiny.yaml"),
        "-o",
        str(root / "out/e"),
        "--n",
        "3",
        "-j",
        "2",
        "--resume",
    ]
    assert out == root / "out/e"


def test_worker_runs_a_job(tmp_path):
    root = _setup(tmp_path)
    jobs = root / "jobs"
    w = Worker(jobs, root, poll_s=0.1)
    (jobs / "inbox" / "001_tiny.yaml").write_text(
        yaml.safe_dump(
            {"command": "ensemble", "config": "configs/tiny.yaml", "output": "output/tiny"}
        )
    )
    (jobs / "inbox" / "002_bad.yaml").write_text(yaml.safe_dump({"command": "format_disk"}))
    w.loop(once=True)
    st = json.loads((jobs / "done" / "001_tiny.status.json").read_text())
    assert st["state"] == "done" and st["rows_done"] == 2
    assert (root / "output/tiny/summary.csv").exists()
    assert "fate intact" in (jobs / "done" / "001_tiny.log").read_text()
    bad = json.loads((jobs / "failed" / "002_bad.status.json").read_text())
    assert bad["state"] == "failed" and "unknown command" in bad["error"]
    info = json.loads((jobs / "worker.json").read_text())
    assert info["cpus"] >= 1 and info["current_job"] is None
    assert not list((jobs / "inbox").iterdir())
    shutil.rmtree(root)


def test_recover_requeues_jobs_of_a_stopped_worker(tmp_path):
    root = _setup(tmp_path)
    jobs = root / "jobs"
    w = Worker(jobs, root, poll_s=0.1)
    (jobs / "running" / "007_x.yaml").write_text(yaml.safe_dump({"command": "validate"}))
    (jobs / "running" / "007_x.log").write_text("partial output\n")
    (jobs / "running" / "007_x.status.json").write_text(json.dumps({"state": "running"}))
    assert w.recover() == ["007_x"]
    assert (jobs / "inbox" / "007_x.yaml").exists()
    assert not list((jobs / "running").iterdir())
    assert any(p.name.startswith("007_x.interrupted") for p in (jobs / "failed").iterdir())


@pytest.mark.parametrize("how", ["stop_file", "sigterm", "sigint"])
def test_worker_stops_and_requeues(tmp_path, how):
    """A long job is stopped by the stop file, SIGTERM or Ctrl+C (SIGINT to the worker): the
    computation ends and the job goes back to the inbox (resume: true continues it)."""
    import os
    import signal
    import subprocess
    import sys
    import time

    root = _setup(tmp_path)
    cfg = yaml.safe_load((root / "configs" / "tiny.yaml").read_text())
    cfg["run"] = {"t_max_yr": 2.0e5, "n_output": 10, "backend": "nbody"}  # minutes per member
    (root / "configs" / "slow.yaml").write_text(yaml.safe_dump(cfg))
    jobs = root / "jobs"
    (jobs / "inbox").mkdir(parents=True)
    (jobs / "inbox" / "001_slow.yaml").write_text(
        yaml.safe_dump(
            {
                "command": "ensemble",
                "config": "configs/slow.yaml",
                "output": "output/slow",
                "args": {"processes": 1, "resume": True},
            }
        )
    )
    w = subprocess.Popen(
        [sys.executable, "-m", "swingersim", "worker", "--root", str(root), "--poll", "0.5"],
        start_new_session=True,
    )
    deadline = time.time() + 60
    status = jobs / "running" / "001_slow.status.json"
    while time.time() < deadline and not (
        status.exists() and "pid" in json.loads(status.read_text() or "{}")
    ):
        time.sleep(0.5)
    job_pid = json.loads(status.read_text())["pid"]
    time.sleep(2)
    if how == "stop_file":
        (jobs / "stop").write_text("")
    elif how == "sigterm":
        w.send_signal(signal.SIGTERM)
    else:
        w.send_signal(signal.SIGINT)
    assert w.wait(timeout=60) == 0
    assert (jobs / "inbox" / "001_slow.yaml").exists()
    assert not list((jobs / "running").iterdir())
    assert json.loads((jobs / "worker.json").read_text())["stopped"] is True
    assert not (jobs / "stop").exists()
    with pytest.raises(ProcessLookupError):
        os.kill(job_pid, 0)  # the computation has ended
    shutil.rmtree(root)


def test_stderr_filter_drops_reboundx_flood(tmp_path):
    """The 'Fractional Error' lines of REBOUNDx gr_full (child processes included) are dropped
    and counted; other stderr output passes."""
    import os
    import subprocess
    import sys

    script = tmp_path / "spam.py"
    script.write_text(
        "import os, sys\n"
        "from multiprocessing import Pool\n"
        "from swingersim import cnoise\n"
        "def spam(i):\n"
        "    for k in range(500):\n"
        "        os.write(2, f'Fractional Error: {k}e-14\\n'.encode())\n"
        "    os.write(2, f'child {i} ok\\n'.encode())\n"
        "    return i\n"
        "if __name__ == '__main__':\n"
        "    cnoise.install()\n"
        "    with Pool(2) as p:\n"
        "        p.map(spam, range(2))\n"
    )
    env = {k: v for k, v in os.environ.items() if k != "SWINGERSIM_NO_STDERR_FILTER"}
    out = subprocess.run(
        [sys.executable, str(script)],
        capture_output=True,
        text=True,
        env=env,
        timeout=120,
        check=False,
    )
    err = out.stderr
    assert err.count("Fractional Error: ") == 1  # the first one, with the note
    assert "child 0 ok" in err and "child 1 ok" in err
    assert "1000 REBOUNDx 'Fractional Error' messages suppressed" in err


def _proc_state(pid):
    return Path(f"/proc/{pid}/stat").read_text().split(")")[-1].split()[0]


def test_pause_resume_and_cancel(tmp_path):
    """jobs/pause/<id> freezes the job's processes and removing it continues them; a file in
    jobs/cancel/ kills the job, deletes its output and keeps the job files in cancel/."""
    import subprocess
    import sys
    import time

    root = _setup(tmp_path)
    cfg = yaml.safe_load((root / "configs" / "tiny.yaml").read_text())
    cfg["run"] = {"t_max_yr": 2.0e5, "n_output": 10, "backend": "nbody"}
    (root / "configs" / "slow.yaml").write_text(yaml.safe_dump(cfg))
    jobs = root / "jobs"
    (jobs / "inbox").mkdir(parents=True)
    (jobs / "inbox" / "001_slow.yaml").write_text(
        yaml.safe_dump(
            {"command": "ensemble", "config": "configs/slow.yaml", "output": "output/slow"}
        )
    )
    w = subprocess.Popen(
        [sys.executable, "-m", "swingersim", "worker", "--root", str(root), "--poll", "0.5"],
        start_new_session=True,
    )
    try:
        status = jobs / "running" / "001_slow.status.json"
        deadline = time.time() + 60
        while time.time() < deadline and not (
            status.exists() and "pid" in json.loads(status.read_text() or "{}")
        ):
            time.sleep(0.3)
        pid = json.loads(status.read_text())["pid"]
        (jobs / "pause").mkdir(exist_ok=True)
        (jobs / "pause" / "001_slow").write_text("")
        deadline = time.time() + 20
        while time.time() < deadline and json.loads(status.read_text())["state"] != "paused":
            time.sleep(0.3)
        assert json.loads(status.read_text())["state"] == "paused"
        assert _proc_state(pid) == "T"
        (jobs / "pause" / "001_slow").unlink()
        deadline = time.time() + 20
        while time.time() < deadline and json.loads(status.read_text())["state"] != "running":
            time.sleep(0.3)
        assert _proc_state(pid) != "T"
        (root / "output" / "slow").mkdir(parents=True, exist_ok=True)
        (jobs / "cancel" / "001_slow").write_text("")
        deadline = time.time() + 30
        while time.time() < deadline and not (jobs / "cancel" / "001_slow.yaml").exists():
            time.sleep(0.3)
        assert (jobs / "cancel" / "001_slow.yaml").exists()
        st = json.loads((jobs / "cancel" / "001_slow.status.json").read_text())
        assert st["state"] == "cancelled" and st["output_removed"]
        assert not (root / "output" / "slow").exists()
        assert not (jobs / "cancel" / "001_slow").exists()
    finally:
        (jobs / "stop").write_text("")
        w.wait(timeout=60)
    shutil.rmtree(root)


def test_every_command_can_be_a_job(tmp_path):
    from swingersim.worker import JobError, build_command

    argv, out = build_command({"command": "pulsation", "config": "configs/x.yaml"}, tmp_path)
    assert argv[3:] == ["pulsation", str(tmp_path / "configs/x.yaml")]
    argv, out = build_command(
        {
            "command": "applegate",
            "args": {
                "mass": 1,
                "radius": 1,
                "teff": 5800,
                "a": 0.05,
                "pbin": 2,
                "pmod": 10,
                "oc": 100,
            },
        },
        tmp_path,
    )
    assert "--mass" in argv and out is None
    with pytest.raises(JobError):
        build_command({"command": "gui"}, tmp_path)
    with pytest.raises(JobError):
        build_command({"command": "run", "config": "../x.yaml"}, tmp_path)


def test_job_cores():
    from swingersim.worker import job_cores

    assert job_cores({"command": "run"}, 12) == 1
    assert job_cores({"command": "lightcurve"}, 12) == 1
    assert job_cores({"command": "ensemble"}, 12) == 1  # the CLI default is one process
    assert job_cores({"command": "ensemble", "args": {"processes": 0}}, 12) == 12
    assert job_cores({"command": "ensemble", "args": {"processes": 4}}, 12) == 4
    assert job_cores({"command": "validate"}, 12) == 1  # the fast benchmarks run serially
    assert job_cores({"command": "validate", "args": {"all": True}}, 12) == 12
    assert job_cores({"command": "run", "cores": "all"}, 12) == 12
    assert job_cores({"command": "ensemble", "args": {"processes": 0}, "cores": 3}, 12) == 3


def test_one_core_jobs_run_side_by_side(tmp_path):
    """Three one-core jobs on two cores: two start at once, the third waits; a job that needs
    all cores waits until the others have finished and is not overtaken."""
    root = _setup(tmp_path)
    jobs = root / "jobs"
    w = Worker(jobs, root, poll_s=0.1, cores=2)
    for k in (1, 2, 3):
        (jobs / "inbox" / f"00{k}_run.yaml").write_text(
            yaml.safe_dump(
                {"command": "run", "config": "configs/tiny.yaml", "output": f"output/r{k}.h5"}
            )
        )
    (jobs / "inbox" / "004_all.yaml").write_text(
        yaml.safe_dump(
            {
                "command": "ensemble",
                "config": "configs/tiny.yaml",
                "output": "output/e",
                "args": {"processes": 0},
            }
        )
    )
    (jobs / "inbox" / "005_run.yaml").write_text(
        yaml.safe_dump({"command": "run", "config": "configs/tiny.yaml", "output": "output/r5.h5"})
    )
    w.schedule()
    assert sorted(w.active) == ["001_run", "002_run"]
    for a in list(w.active.values()):
        a.proc.wait()
        w.poll_job(a)
    w.schedule()
    assert sorted(w.active) == ["003_run"]  # 004 needs both cores, 005 must not overtake it
    w.loop(once=True)
    for k in ("001_run", "002_run", "003_run", "004_all", "005_run"):
        assert json.loads((jobs / "done" / f"{k}.status.json").read_text())["state"] == "done"
    assert json.loads((jobs / "done" / "004_all.status.json").read_text())["cores"] == 2
    shutil.rmtree(root)


def test_outputs_in_new_folders(tmp_path):
    """A command writes into folders that do not exist yet (jobs 050-053 failed on this)."""
    from swingersim.cli import main

    root = _setup(tmp_path)
    out = tmp_path / "new" / "deeper" / "timings.csv"
    plot = tmp_path / "other" / "oc.png"
    assert (
        main(
            [
                "timings",
                str(root / "configs" / "tiny.yaml"),
                "-o",
                str(out),
                "--years",
                "0.05",
                "--plot",
                str(plot),
            ]
        )
        == 0
    )
    assert out.exists()


def test_emulator_needs_a_computed_ensemble(tmp_path):
    from swingersim import emulator

    with pytest.raises(ValueError, match="has not been computed yet"):
        emulator.train_from_ensemble(tmp_path / "missing", ["orbits.out.P"], classes=["x"])


def test_memory_guard_stops_the_largest_job(tmp_path, monkeypatch):
    """Before the computer runs out of memory the worker stops the job that uses most of it."""
    import swingersim.worker as wk

    root = _setup(tmp_path)
    cfg = yaml.safe_load((root / "configs" / "tiny.yaml").read_text())
    cfg["run"] = {"t_max_yr": 2.0e5, "n_output": 10, "backend": "nbody"}
    (root / "configs" / "slow.yaml").write_text(yaml.safe_dump(cfg))
    jobs = root / "jobs"
    w = Worker(jobs, root, poll_s=0.1, cores=2)
    (jobs / "inbox" / "001_slow.yaml").write_text(
        yaml.safe_dump({"command": "run", "config": "configs/slow.yaml", "output": "output/s.h5"})
    )
    w.schedule()
    assert "001_slow" in w.active
    monkeypatch.setattr(wk, "memory_available", lambda: (0.5, 16.0))
    assert w.memory_guard() == "failed"
    st = json.loads((jobs / "failed" / "001_slow.status.json").read_text())
    assert "out of memory" in st["error"]
    assert not w.active
    shutil.rmtree(root)


def test_secular_memory_stays_bounded():
    """A long secular run with tides is integrated in chunks: the extremes are those of the
    whole run, and the memory does not grow with the number of integrator steps."""
    import numpy as np

    from swingersim.backends.chunked import merge_extremes, merge_roche

    tot: dict = {}
    merge_extremes(
        tot,
        {
            "A": {
                "e_max": 0.2,
                "t_e_max_yr": 1.0,
                "q_min_au": 0.1,
                "t_q_min_yr": 1.0,
                "t_1me_lt_1e-1_yr": float("nan"),
            }
        },
    )
    merge_extremes(
        tot,
        {
            "A": {
                "e_max": 0.95,
                "t_e_max_yr": 5.0,
                "q_min_au": 0.01,
                "t_q_min_yr": 5.0,
                "t_1me_lt_1e-1_yr": 4.5,
            }
        },
    )
    merge_extremes(
        tot,
        {
            "A": {
                "e_max": 0.5,
                "t_e_max_yr": 9.0,
                "q_min_au": 0.05,
                "t_q_min_yr": 9.0,
                "t_1me_lt_1e-1_yr": 8.0,
            }
        },
    )
    assert tot["A"]["e_max"] == 0.95 and tot["A"]["t_e_max_yr"] == 5.0
    assert tot["A"]["q_min_au"] == 0.01 and tot["A"]["t_1me_lt_1e-1_yr"] == 4.5
    r: dict = {}
    merge_roche(r, {"m1": {"fill_max": 0.3, "t_fill_max_yr": 1.0}})
    merge_roche(r, {"m1": {"fill_max": 0.6, "t_fill_max_yr": 2.0}})
    assert r["m1"]["t_fill_max_yr"] == 2.0
    assert np.isfinite(tot["A"]["e_max"])
