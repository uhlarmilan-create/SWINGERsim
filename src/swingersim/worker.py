"""Job worker: runs SWINGERsim jobs dropped as files into a folder.

Meant for a long-running computer (e.g. the user's PC) that shares a folder with whoever
prepares the jobs. Layout of the job folder (created on start):

    jobs/
      drafts/     configurations written with swingersim gui (not run by the worker)
      submitted/  frozen copies of the configurations of jobs made from drafts
      inbox/      <id>.yaml  queued jobs, picked up in order of their file names;
                  <id>.yaml.hold  jobs waiting until they are queued (renamed to .yaml)
      running/    the jobs being computed: <id>.yaml, <id>.log, <id>.status.json; jobs that
                  need one core (a run, a light curve, minima, ...) run side by side up to
                  the number of cores, a job that needs all cores (an ensemble with
                  processes: 0, validate --all, ...) runs alone (``cores:`` in a job file overrides)
      done/       finished jobs (same three files; the results are where `output` says)
      failed/     jobs that ended with an error
      cancel/     an empty file <id> here stops that job: its process is killed, its output
                  (inside output/) is deleted and the job files are kept here as <id>.*
      pause/      an empty file <id> here pauses the running job (SIGSTOP of its processes);
                  removing the file continues it exactly where it stopped (SIGCONT)
      stop        an empty file here stops the worker: running jobs go back to the inbox
                  (an ensemble with resume: true continues where it stopped) - the same as
                  Ctrl+C or SIGTERM, but possible from another computer
      worker.json heartbeat: host, CPUs, versions, current job, time of the last check
      worker.lock held by the running worker (a second worker on the same folder refuses to
                  start)

Every job gets SWINGERSIM_CORES = the cores it was booked for: ``processes: 0`` (all cores) of an
ensemble, scattering, chaos map or fit then means the job's cores, not the whole computer.

A job file:

    command: ensemble            # ensemble | run | lightcurve | ensemble-report | ...
    config: configs/squadra/czev4315_ensemble.yaml
    output: output/squadra/czev4315_20k
    args: {n: 20000, processes: 8, resume: true}
    note: free text

Only SWINGERsim commands are executed, never shell commands, and every path must lie inside the
worker's root directory (the SWINGERsim folder), so a job file cannot touch anything else.
"""

from __future__ import annotations

import fcntl
import json
import os
import platform
import shutil
import signal
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

import yaml

from . import __version__

SUBDIRS = ("drafts", "submitted", "inbox", "running", "done", "failed", "cancel", "pause")


class JobError(ValueError):
    pass


def _now() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")


def _inside(root: Path, value: str, what: str) -> Path:
    """Path of a job value under root (swingersim.rootpath.inside_root, the rule of the GUI)."""
    from .rootpath import inside_root

    try:
        return inside_root(root, value)
    except ValueError as exc:
        raise JobError(f"{what} '{value}' is outside the SWINGERsim folder ({exc})") from exc


NOT_JOBS = ("gui", "worker")  # interactive commands that cannot run as jobs
OUTPUT_POSITIONALS = ("directory", "h5")  # positional arguments filled from 'output'


def _subparser(cmd: str):
    import argparse

    from .cli import build_parser

    p = build_parser()
    sub = next(a for a in p._actions if isinstance(a, argparse._SubParsersAction))
    if cmd not in sub.choices or cmd in NOT_JOBS:
        allowed = sorted(set(sub.choices) - set(NOT_JOBS))
        raise JobError(f"unknown command '{cmd}' (allowed: {', '.join(allowed)})")
    return sub.choices[cmd]


def build_command(job: dict, root: Path) -> tuple[list[str], Path | None]:
    """Command line (python -m swingersim ...) of a job, after checking it; returns also the
    output path (for progress reports and the cancellation clean-up).

    Every SWINGERsim command except gui and worker can be a job; its arguments are read from the
    command-line parser of SWINGERsim: positional arguments from the top-level keys of the job
    (``config``, ``ensemble``, ``ids``, ...; ``output`` for ``directory`` and ``h5``), options
    from ``args`` by their name (``args: {n: 1000, processes: 0, resume: true}``), ``output``
    as ``-o``. Every path argument must lie inside the SWINGERsim folder."""
    import argparse

    cmd = job.get("command")
    sp = _subparser(str(cmd))
    argv = [sys.executable, "-m", "swingersim", cmd]
    args = dict(job.get("args") or {})
    output = None
    actions = {a.dest: a for a in sp._actions if not isinstance(a, argparse._HelpAction)}

    def path(value, what):
        return str(_inside(root, str(value), what))

    def plain(value, what):
        """A value as a command-line word: never something argparse reads as an option (a job
        value '--html=/elsewhere' would otherwise become an option and escape the folder)."""
        v = str(value)
        if v.startswith("-"):
            try:
                float(v)
            except ValueError:
                raise JobError(f"value '{v}' of '{what}' looks like an option") from None
        return v

    used_output = False
    for a in actions.values():
        if a.option_strings:
            continue
        value = job.get(a.dest, args.pop(a.dest, None))
        if value is None and a.dest in OUTPUT_POSITIONALS and "output" in job:
            value, used_output = job["output"], True
            output = _inside(root, str(value), "output")
        if value is None:
            if a.nargs in ("*", "?"):
                continue
            raise JobError(f"command '{cmd}' needs '{a.dest}'")
        values = value if isinstance(value, (list, tuple)) else [value]
        argv += [path(v, a.dest) if a.type is Path else plain(v, a.dest) for v in values]
    if "output" in job and not used_output:
        if "output" not in actions:
            raise JobError(f"command '{cmd}' has no output option")
        output = _inside(root, str(job["output"]), "output")
        argv += [actions["output"].option_strings[0], str(output)]
    if cmd in ("mesa-grid", "mesa-reduce") and "root" not in args:
        argv += ["--root", str(root)]
    for key, value in args.items():
        a = actions.get(key)
        if a is None or not a.option_strings or key == "output":
            allowed = [d for d, x in actions.items() if x.option_strings and d != "output"]
            raise JobError(f"argument '{key}' is not allowed for '{cmd}' ({', '.join(allowed)})")
        flag = a.option_strings[0]
        if isinstance(a, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
            if value:
                argv.append(flag)
            continue
        if value is None:
            continue
        values = value if isinstance(value, (list, tuple)) else [value]
        values = [path(v, key) if a.type is Path else plain(v, key) for v in values]
        if a.nargs in ("+", "*") or isinstance(a.nargs, int):
            argv += [flag, *values]
        else:
            if len(values) != 1:
                raise JobError(f"argument '{key}' takes one value")
            argv += [flag, values[0]]
    return argv, output


def _signal_group(pid: int, sig) -> None:
    try:
        os.killpg(pid, sig)
    except (ProcessLookupError, PermissionError):
        pass


def _group_alive(pgid: int) -> bool:
    try:
        os.killpg(pgid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _terminate_group(pid: int, timeout: float = 30.0, proc=None) -> None:
    """Stop a job's process group (the job and its pool workers): SIGCONT (a paused job must
    run to receive the signal), SIGTERM, then SIGKILL if anything is left after ``timeout``
    seconds. ``proc`` (the Popen of the job, when it is a child of this worker) is reaped."""
    _signal_group(pid, signal.SIGCONT)
    _signal_group(pid, signal.SIGTERM)
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if proc is not None:
            proc.poll()  # reap the job itself: a zombie still counts as a group member
        if not _group_alive(pid):
            break
        time.sleep(0.1)
    else:
        _signal_group(pid, signal.SIGKILL)
    if proc is not None:
        proc.wait()
    end = time.monotonic() + 5.0
    while _group_alive(pid) and time.monotonic() < end:
        time.sleep(0.05)


def remove_output(root: Path, output: Path | None) -> bool:
    """Delete the output of a cancelled job - only inside <root>/output, never anything else.
    A run also leaves <output>.partial.json."""
    if output is None:
        return False
    base = (root / "output").resolve()
    out = output.resolve()
    if out == base or base not in out.parents:
        return False
    done = False
    for p in (out, out.with_name(out.name + ".partial.json")):
        if p.is_dir():
            shutil.rmtree(p, ignore_errors=True)
            done = True
        elif p.exists():
            p.unlink()
            done = True
    # an emptied folder of its own (output/<id>/) goes too
    parent = out.parent
    if base in parent.parents and parent.is_dir() and not any(parent.iterdir()):
        parent.rmdir()
    return done


def _stop_job(proc) -> None:
    """Stop a job's process group (`_terminate_group`); further Ctrl+C presses during the
    cleanup are ignored so that the job is always requeued."""
    old = signal.signal(signal.SIGINT, signal.SIG_IGN)
    try:
        _terminate_group(proc.pid, 30.0, proc)
    finally:
        signal.signal(signal.SIGINT, old)


class StopRequested(KeyboardInterrupt):
    """The stop file appeared or SIGTERM arrived: stop like Ctrl+C."""


def _sigterm(_signum, _frame):
    raise StopRequested("SIGTERM")


def _job_processes(status: dict) -> list[int]:
    """Process ids of a job's computation that may still run: the recorded pid, or (for jobs
    started by an older worker) processes whose command line ends with the job's command."""
    pids = []
    if status.get("pid"):
        pids.append(int(status["pid"]))
    elif status.get("command_line"):
        target = status["command_line"].split()
        for d in Path("/proc").glob("[0-9]*"):
            try:
                args = (d / "cmdline").read_bytes().decode().split("\0")
            except OSError:
                continue
            args = [a for a in args if a]
            if len(args) >= len(target) and args[-len(target) :] == target and "-m" in args:
                pids.append(int(d.name))
    alive = []
    for pid in pids:
        try:
            cmd = Path(f"/proc/{pid}/cmdline").read_bytes()
        except OSError:
            continue
        if b"swingersim" in cmd:  # never touch anything that is not a SWINGERsim computation
            alive.append(pid)
    return alive


def _progress(cmd: str, output: Path | None) -> dict:
    """What the job has done so far: rows of an ensemble / scattering / chaos map, or the time
    reached by an N-body run (from its intermediate-results file)."""
    if output is None:
        return {}
    name = "grid.jsonl" if cmd == "chaosmap" else "rows.jsonl"
    if cmd in ("ensemble", "scatter", "chaosmap") and (output / name).exists():
        with open(output / name, encoding="utf-8") as f:
            return {"rows_done": sum(1 for line in f if line.strip())}
    partial = output.with_name(output.name + ".partial.json")
    if cmd == "run" and partial.exists():
        try:
            d = json.loads(partial.read_text())
            return {"t_yr": d["time_yr"][-1], "t_end_yr": d["t_end_yr"]}
        except (OSError, ValueError, KeyError, IndexError):
            return {}
    return {}


def _rows_total(job: dict, root: Path) -> int | None:
    """Rows an ensemble will write: members x directions (None when unknown)."""
    if job.get("command") != "ensemble" or "config" not in job:
        return None
    try:
        cfg = yaml.safe_load((root / job["config"]).read_text()) or {}
        ens = cfg.get("ensemble") or {}
        n = int((job.get("args") or {}).get("n") or ens.get("n") or 0)
        return n * len(ens.get("directions") or ["forward"]) or None
    except (OSError, ValueError, TypeError, yaml.YAMLError):
        return None


# commands that use several cores by themselves; every other command is one calculation on
# one core, and several of them can run at the same time
ALL_CORES = ("emulator-train", "sbi")  # parallel inside (numpy, torch threads)
THREADED = ALL_CORES  # their numerical libraries get the job's cores as threads
PROCESS_DEFAULT = {
    "ensemble": 1,
    "scatter": 1,
    "chaosmap": 1,
    "mesa-grid": 0,
    "mesa-reduce": 0,
}


def _fit_processes(job: dict, root: Path | None) -> int:
    """Processes of a fit job: --processes, else fit.stage2.processes of its configuration
    (0 = all cores) for a stage-2 fit; a stage-1 fit or an anchoring runs on one core."""
    args = job.get("args") or {}
    if args.get("processes") is not None:
        return int(args["processes"] or 0)
    if args.get("anchor"):
        return 1
    try:
        cfg = yaml.safe_load((Path(root or ".") / str(job["config"])).read_text()) or {}
        stage2 = (cfg.get("fit") or {}).get("stage2")
    except (OSError, KeyError, TypeError, AttributeError, yaml.YAMLError):
        stage2 = None
    stage = args.get("stage") or (2 if stage2 else 1)
    if int(stage) != 2 or not isinstance(stage2, dict):
        return 1
    return int(stage2.get("processes") or 0)


def job_cores(job: dict, total: int | None = None, root: Path | None = None) -> int:
    """Cores a job occupies while it runs: ``cores:`` of the job file if given (a number or
    "all"), else ``processes`` of a parallel command (0 = all cores; for a fit
    fit.stage2.processes of its configuration, read from ``root``), all cores for
    emulator-train / sbi and for validate --all (its slow benchmarks run ensembles), and one
    core for everything else (a run, a light curve, minima, timings, validate of the fast
    benchmarks, which run one after another, ...)."""
    total = total or os.cpu_count() or 1
    want = job.get("cores")
    if want is not None:
        return total if str(want) == "all" else max(1, min(int(want), total))
    cmd = job.get("command")
    args = job.get("args") or {}
    if cmd in ALL_CORES:
        return total
    if cmd == "validate":
        return total if args.get("all") else 1
    if cmd == "fit":
        n = _fit_processes(job, root)
        return total if n <= 0 else min(n, total)
    if cmd in PROCESS_DEFAULT:
        n = args.get("processes", PROCESS_DEFAULT[cmd])
        n = int(n or 0)
        return total if n <= 0 else min(n, total)
    return 1


class _Active:
    """A job that is being computed."""

    def __init__(self, job_id, job, proc, logf, output, status, cores):
        self.job_id, self.job, self.proc, self.logf = job_id, job, proc, logf
        self.output, self.status, self.cores = output, status, cores
        self.t0 = time.time()
        self.paused_since: float | None = None
        self.last_status = 0.0


def memory_available() -> tuple[float, float] | None:
    """(available, total) memory of the computer in GB, from /proc/meminfo (Linux)."""
    try:
        info = {}
        for line in Path("/proc/meminfo").read_text().splitlines():
            k, v = line.split(":", 1)
            info[k] = float(v.split()[0]) / 1024**2
        return info["MemAvailable"], info["MemTotal"]
    except (OSError, KeyError, ValueError):
        return None


def group_memory_gb(pgid: int) -> float:
    """Resident memory of all processes of a job's process group [GB]."""
    page = os.sysconf("SC_PAGE_SIZE") if hasattr(os, "sysconf") else 4096
    total = 0
    for d in Path("/proc").glob("[0-9]*"):
        try:
            stat = (d / "stat").read_text()
            fields = stat[stat.rindex(")") + 2 :].split()
            if int(fields[2]) != pgid:  # fields: state, ppid, pgrp, ...
                continue
            total += int((d / "statm").read_text().split()[1]) * page
        except (OSError, ValueError, IndexError):
            continue
    return total / 1024**3


# the worker stops the largest job before the computer runs out of memory (a desktop that
# swaps or whose session is killed by the kernel is worse than one failed job)
MEMORY_RESERVE_FRACTION = 0.08
MEMORY_RESERVE_GB = 1.5


class WorkerLocked(RuntimeError):
    """Another worker already runs on this job folder."""


class Worker:
    """Runs the queued jobs of the inbox in the order of their file names. Jobs that need one
    core run side by side up to ``cores`` (default: all cores of the computer); a job that
    needs all cores (an ensemble with processes: 0, validate --all, ...) waits until the others
    have finished and runs alone. Only one worker may run on a job folder (jobs/worker.lock)."""

    def __init__(self, jobdir: Path, root: Path, poll_s: float = 10.0, cores: int | None = None):
        self.jobdir = jobdir.resolve()
        self.root = root.resolve()
        self.poll_s = poll_s
        self.cores = max(1, int(cores or os.cpu_count() or 1))
        for d in SUBDIRS:
            (self.jobdir / d).mkdir(parents=True, exist_ok=True)
        self._lock = open(self.jobdir / "worker.lock", "a+", encoding="utf-8")  # noqa: SIM115
        try:
            fcntl.flock(self._lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            self._lock.seek(0)
            other = self._lock.read().strip()
            self._lock.close()
            raise WorkerLocked(
                f"another swingersim worker is already running on {self.jobdir}"
                + (f" ({other})" if other else "")
                + "; stop it first (an empty file 'stop' in the job folder) or use another folder"
            ) from None
        self._lock.seek(0)
        self._lock.truncate()
        self._lock.write(f"pid {os.getpid()} on {platform.node()} since {_now()}\n")
        self._lock.flush()
        self.active: dict[str, _Active] = {}

    def close(self) -> None:
        """Release the job-folder lock (also released when the process ends)."""
        if not self._lock.closed:
            fcntl.flock(self._lock, fcntl.LOCK_UN)
            self._lock.close()

    # ------------------------------------------------------------------ bookkeeping
    def heartbeat(self, extra: dict | None = None) -> None:
        import importlib.util

        info = {
            "host": platform.node(),
            "cpus": os.cpu_count(),
            "cores": self.cores,
            "cores_used": sum(a.cores for a in self.active.values()),
            "python": platform.python_version(),
            "swingersim": __version__,
            "numba": importlib.util.find_spec("numba") is not None,
            "root": str(self.root),
            "current_job": next(iter(self.active), None),
            "current_jobs": list(self.active),
            "paused": any(a.paused_since is not None for a in self.active.values()),
            "last_check": _now(),
            "pid": os.getpid(),
        }
        info.update(extra or {})
        tmp = self.jobdir / "worker.json.tmp"
        tmp.write_text(json.dumps(info, indent=1))
        tmp.replace(self.jobdir / "worker.json")

    def _status(self, where: Path, job_id: str, data: dict) -> None:
        """Write <id>.status.json atomically (the GUI reads it at any moment)."""
        data = {k: v for k, v in data.items() if not k.startswith("_")}
        path = where / f"{job_id}.status.json"
        tmp = where / f".{job_id}.status.json.tmp"
        tmp.write_text(json.dumps(data, indent=1))
        tmp.replace(path)

    # ------------------------------------------------------------------ jobs
    def queued(self) -> list[Path]:
        return sorted(p for p in (self.jobdir / "inbox").iterdir() if p.suffix in (".yaml", ".yml"))

    def _cores_of(self, path: Path) -> int:
        try:
            job = yaml.safe_load(path.read_text()) or {}
            return job_cores(job if isinstance(job, dict) else {}, self.cores, self.root)
        except (OSError, yaml.YAMLError, ValueError, TypeError):
            return 1  # a broken job fails at once when it starts

    def start_job(self, path: Path) -> str | None:
        """Move a job to running/ and start its computation; returns its final state when it
        could not start (a bad job file), else None."""
        job_id = path.stem
        running = self.jobdir / "running"
        job_file = running / path.name
        try:
            shutil.move(str(path), job_file)
        except FileNotFoundError:  # renamed (held) or removed since the inbox was listed
            return None
        status = {"id": job_id, "state": "running", "submitted_file": path.name, "started": _now()}
        try:
            job = yaml.safe_load(job_file.read_text()) or {}
            if not isinstance(job, dict):
                raise JobError("the job file must be a YAML mapping")
            argv, output = build_command(job, self.root)
            cores = job_cores(job, self.cores, self.root)
        except (JobError, yaml.YAMLError, ValueError, TypeError) as exc:
            status.update(state="failed", error=str(exc), finished=_now())
            return self._finish(job_id, status, "failed")
        status["command_line"] = " ".join(argv[2:])
        status["cores"] = cores
        total = _rows_total(job, self.root)
        if total:
            status["rows_total"] = total
        if output is not None:
            status["output"] = str(output.relative_to(self.root))
        status.setdefault("paused_s", 0)
        env = dict(os.environ)
        # processes: 0 inside the job = the cores it was booked for (ensemble.resolve_processes)
        env["SWINGERSIM_CORES"] = str(cores)
        # one calculation per process: no extra BLAS / OpenMP / numba threads, except for the
        # commands whose parallelism is threads (torch, scikit-learn)
        threads = str(cores) if job.get("command") in THREADED else "1"
        for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            env[k] = threads
        if cores == 1:
            env["NUMBA_NUM_THREADS"] = "1"
        logf = open(running / f"{job_id}.log", "w", encoding="utf-8")  # noqa: SIM115
        proc = subprocess.Popen(
            argv,
            cwd=self.root,
            stdout=logf,
            stderr=subprocess.STDOUT,
            start_new_session=True,
            env=env,
        )
        status["pid"] = proc.pid
        self.active[job_id] = _Active(job_id, job, proc, logf, output, status, cores)
        self._status(running, job_id, status)
        print(f"[{_now()}] job {job_id}: started ({cores} core{'s' * (cores > 1)})", flush=True)
        return None

    def poll_job(self, a: _Active) -> str | None:
        """One check of an active job: cancel, pause, progress; its final state when it ended."""
        running = self.jobdir / "running"
        pause = self.jobdir / "pause" / a.job_id
        if (self.jobdir / "cancel" / a.job_id).exists():
            _terminate_group(a.proc.pid, 30.0, a.proc)
            return self._end(a, cancelled=True)
        if a.proc.poll() is not None:
            return self._end(a)
        if pause.exists() and a.paused_since is None:
            _signal_group(a.proc.pid, signal.SIGSTOP)
            a.paused_since = time.time()
            a.status.update(state="paused", paused_at=_now())
        elif not pause.exists() and a.paused_since is not None:
            _signal_group(a.proc.pid, signal.SIGCONT)
            a.status["paused_s"] += round(time.time() - a.paused_since)
            a.paused_since = None
            a.status.update(state="running")
            a.status.pop("paused_at", None)
        now = time.time()
        if now - a.last_status >= min(self.poll_s, 5.0):
            a.status.update(_progress(a.job.get("command"), a.output), elapsed_s=round(now - a.t0))
            self._status(running, a.job_id, a.status)
            a.last_status = now
        return None

    def _end(self, a: _Active, cancelled: bool = False) -> str:
        # pool workers of a job whose main process has ended (crashed, killed) must not go on
        # computing as orphans on cores that the worker now gives to other jobs
        _terminate_group(a.proc.pid, 5.0, a.proc)
        a.logf.close()
        del self.active[a.job_id]
        status = a.status
        status.update(
            elapsed_s=round(time.time() - a.t0),
            returncode=a.proc.returncode,
            finished=_now(),
            **_progress(a.job.get("command"), a.output),
        )
        (self.jobdir / "pause" / a.job_id).unlink(missing_ok=True)
        if cancelled:
            (self.jobdir / "cancel" / a.job_id).unlink(missing_ok=True)
            status["state"] = "cancelled"
            status["output_removed"] = remove_output(self.root, a.output)
            state = self._finish(a.job_id, status, "cancel")
        else:
            status["state"] = "done" if a.proc.returncode == 0 else "failed"
            state = self._finish(a.job_id, status, status["state"])
        print(f"[{_now()}] job {a.job_id}: {state}", flush=True)
        return state

    def memory_guard(self) -> str | None:
        """When the free memory of the computer falls below the reserve, stop the job that
        uses the most memory: it ends as failed with an explanation (its results so far stay,
        so an ensemble with resume continues later, e.g. with fewer processes)."""
        mem = memory_available()
        if mem is None or not self.active:
            return None
        avail, total = mem
        if avail >= max(MEMORY_RESERVE_GB, MEMORY_RESERVE_FRACTION * total):
            return None
        use = {jid: group_memory_gb(a.proc.pid) for jid, a in self.active.items()}
        jid = max(use, key=use.get)
        a = self.active[jid]
        _terminate_group(a.proc.pid, 30.0, a.proc)
        a.status["error"] = (
            f"stopped by the worker: the computer was running out of memory "
            f"({avail:.1f} of {total:.1f} GB free; this job used {use[jid]:.1f} GB). "
            "Its results so far are kept; run it again with fewer processes "
            "(args: processes) or alone"
        )
        print(f"[{_now()}] job {jid}: stopped, out of memory ({use[jid]:.1f} GB)", flush=True)
        return self._end(a)

    def stop_all(self, reason: str) -> None:
        """Stop every active job and put it back into the inbox."""
        for a in list(self.active.values()):
            _stop_job(a.proc)
            a.logf.close()
            del self.active[a.job_id]
            self._requeue(a.job_id, reason)

    def _requeue(self, job_id: str, reason: str) -> None:
        """Put an interrupted job back into the inbox (an ensemble with resume: true continues
        where it stopped); its log is kept in failed/."""
        running = self.jobdir / "running"
        stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S")
        for f in running.glob(f"{job_id}.*"):
            if f.suffix in (".yaml", ".yml"):
                shutil.move(str(f), self.jobdir / "inbox" / f.name)
            elif f.name.endswith(".log"):
                shutil.move(str(f), self.jobdir / "failed" / f"{job_id}.interrupted-{stamp}.log")
            else:
                f.unlink()
        (self.jobdir / "failed" / f"{job_id}.interrupted-{stamp}.txt").write_text(reason + "\n")
        (self.jobdir / "pause" / job_id).unlink(missing_ok=True)  # requeued jobs start running

    def recover(self) -> list[str]:
        """Jobs left in running/ by a worker that stopped: stop their computation if it still
        runs (also when it was paused: SIGCONT, SIGTERM, SIGKILL after 30 s) and, once it has
        ended, put them back into the inbox."""
        recovered = []
        for st in (self.jobdir / "running").glob("*.status.json"):
            job_id = st.name.removesuffix(".status.json")
            try:
                status = json.loads(st.read_text())
            except ValueError:
                status = {}
            for pid in _job_processes(status):
                _terminate_group(pid, 30.0)
            self._requeue(job_id, "the worker stopped while the job was running; requeued")
            recovered.append(job_id)
        return recovered

    def _finish(self, job_id: str, status: dict, where: str) -> str:
        running, dest = self.jobdir / "running", self.jobdir / where
        for f in running.glob(f"{job_id}.*"):
            if not f.name.endswith(".status.json"):
                shutil.move(str(f), dest / f.name)
        (running / f"{job_id}.status.json").unlink(missing_ok=True)
        self._status(dest, job_id, status)
        self.heartbeat()
        return status["state"]

    def check_stop(self) -> None:
        stop = self.jobdir / "stop"
        if stop.exists():
            stop.unlink(missing_ok=True)
            raise StopRequested("stop file")

    # ------------------------------------------------------------------ loop
    def schedule(self) -> None:
        """Start queued jobs while cores are free, in the order of the inbox: the first job
        that does not fit waits (later jobs do not overtake it)."""
        for path in self.queued():
            need = self._cores_of(path)
            used = sum(a.cores for a in self.active.values())
            if self.active and used + need > self.cores:
                return
            if path.stem in self.active or not path.exists():
                continue
            state = self.start_job(path)
            if state is not None:
                print(f"[{_now()}] job {path.stem}: {state}", flush=True)

    def loop(self, once: bool = False) -> None:
        signal.signal(signal.SIGTERM, _sigterm)
        (self.jobdir / "stop").unlink(missing_ok=True)  # a stale stop file from before
        for job_id in self.recover():
            print(f"[{_now()}] job {job_id}: requeued after a previous stop", flush=True)
        last_beat = 0.0
        try:
            while True:
                self.check_stop()
                for a in list(self.active.values()):
                    self.poll_job(a)
                self.memory_guard()
                self.schedule()
                now = time.time()
                if now - last_beat >= min(self.poll_s, 5.0):
                    self.heartbeat()
                    last_beat = now
                if once and not self.active and not self.queued():
                    self.heartbeat()
                    return
                time.sleep(min(self.poll_s, 1.0) if self.active else self.poll_s)
        except KeyboardInterrupt as exc:
            self.stop_all(_reason(exc))
            raise


def _reason(exc: BaseException) -> str:
    return (
        "stopped by the stop file or SIGTERM"
        if isinstance(exc, StopRequested)
        else "interrupted by Ctrl+C"
    )
