"""Process pools for long parallel computations (ensembles, chaos maps, scattering).

``multiprocessing.Pool.imap_unordered`` waits forever when one of its worker processes dies
(a segmentation fault in a C library, the out-of-memory killer): the task is lost and the
result never arrives. `run_tasks` uses ``concurrent.futures.ProcessPoolExecutor`` instead, whose
pool breaks visibly (``BrokenProcessPool``) when a worker dies. At most ``processes`` tasks are
in flight at a time, so the tasks that may have killed the worker are known: they are run again
in a new pool, and a task that was in flight during two breaks is run alone in its own process;
if that process dies too, the task is reported as lost (``lost(task, message)`` gives the row
to record) and the computation goes on with the others.
"""

from __future__ import annotations

import multiprocessing as mp
import os
import signal
from collections import deque
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from concurrent.futures.process import BrokenProcessPool

LOST_MESSAGE = (
    "the worker process computing this task died (segmentation fault or killed, e.g. by the "
    "out-of-memory killer)"
)


def _context():
    """The fork start method where it exists (the swingersim default on Linux: the workers
    inherit the loaded modules, grids and the stderr filter), else the platform default."""
    try:
        return mp.get_context("fork")
    except ValueError:  # pragma: no cover - Windows
        return mp.get_context()


def _worker_init(initializer, initargs) -> None:
    """Start of a worker process: default signal handlers (a forked worker inherits Python
    handlers of its parent, e.g. the SIGTERM handler of swingersim worker, and could then ignore
    the termination by the pool and block its shutdown), then the caller's initializer."""
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    if initializer is not None:
        initializer(*initargs)


def _processes(ex: ProcessPoolExecutor) -> list:
    return list((getattr(ex, "_processes", None) or {}).values())


def _close(ex: ProcessPoolExecutor, wait_: bool, procs) -> None:
    """Shut a pool down; without waiting, its worker processes are killed first (Ctrl+C, an
    error in the consumer or a broken pool must not leave computing orphans behind, and the
    pool's own clean-up must not wait for a worker that does not end)."""
    if not wait_:
        for p in list(procs) + _processes(ex):
            try:
                if p.is_alive():
                    os.kill(p.pid, signal.SIGKILL)
            except (OSError, ValueError, AttributeError):
                pass
    ex.shutdown(wait=wait_, cancel_futures=True)
    if not wait_:
        for p in list(procs):
            try:
                p.join(timeout=5)
            except (OSError, ValueError, AssertionError):
                pass


def run_tasks(
    fn: Callable,
    tasks: Iterable,
    processes: int,
    lost: Callable[[object, str], object],
    initializer: Callable | None = None,
    initargs: tuple = (),
) -> Iterator:
    """Yield ``fn(task)`` for every task in the order the tasks finish, computed by up to
    ``processes`` worker processes (``initializer(*initargs)`` runs once in each). A task whose
    worker process died yields ``lost(task, message)`` instead. With ``processes <= 1`` the
    tasks run one after another in this process (``initializer`` first)."""
    tasks = list(tasks)
    if processes <= 1 or len(tasks) <= 1:
        if initializer is not None:
            initializer(*initargs)
        for t in tasks:
            yield fn(t)
        return
    pending = deque(enumerate(tasks))
    strikes: dict[int, int] = {}
    isolate: deque = deque()
    ctx = _context()
    init = (initializer, initargs)
    n = min(processes, len(tasks))
    while pending or isolate:
        if isolate:  # a suspect alone in its own process: a death now is its own
            k, t = isolate.popleft()
            ex = ProcessPoolExecutor(1, mp_context=ctx, initializer=_worker_init, initargs=init)
            ok = False
            procs: set = set()
            try:
                fut = ex.submit(fn, t)
                procs.update(_processes(ex))
                try:
                    result = fut.result()
                except BrokenProcessPool:
                    result = lost(t, LOST_MESSAGE)
                ok = True
            finally:
                _close(ex, ok, procs)
            yield result
            continue
        ex = ProcessPoolExecutor(n, mp_context=ctx, initializer=_worker_init, initargs=init)
        inflight: dict = {}
        procs = set()
        clean = False
        try:
            try:
                while pending or inflight:
                    while pending and len(inflight) < n:
                        k, t = pending.popleft()
                        inflight[ex.submit(fn, t)] = (k, t)
                    procs.update(_processes(ex))
                    done, _ = wait(list(inflight), return_when=FIRST_COMPLETED)
                    for f in done:
                        k, t = inflight[f]
                        result = f.result()  # BrokenProcessPool when a worker died
                        del inflight[f]
                        yield result
                clean = True
            except BrokenProcessPool:
                suspects = []
                for f, (k, t) in inflight.items():
                    if f.done() and not f.cancelled() and f.exception() is None:
                        yield f.result()  # finished before the pool broke
                    else:
                        suspects.append((k, t))
                if len(suspects) == 1:  # the only task in flight killed its worker
                    yield lost(suspects[0][1], LOST_MESSAGE)
                else:
                    for k, t in reversed(suspects):
                        strikes[k] = strikes.get(k, 0) + 1
                        if strikes[k] >= 2:
                            isolate.appendleft((k, t))
                        else:
                            pending.appendleft((k, t))
        finally:
            _close(ex, clean, procs)
