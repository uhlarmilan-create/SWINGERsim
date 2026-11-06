"""Filter of messages that compiled libraries print straight to stderr.

REBOUNDx (gr_full.c) prints "Fractional Error: ..." to the C stderr every time its iterative
1PN acceleration does not reach machine precision in 10 substitutions. Whether this happens
depends on the compiler flags and the CPU (on one machine never, on another at every force
evaluation), and the relative deviation it prints can be huge for an acceleration component
close to zero while the absolute error is negligible; the energy check of the run measures the
real accuracy. A worker job once wrote 2.9 GB of these lines in 20 minutes. The filter replaces
file descriptor 2 by a pipe, drops these lines (the first one is shown with a note, the rest
are counted and summarised at exit) and passes everything else through. Child processes started
afterwards (process pools) inherit the pipe, so one filter in the main process covers them.

If the filter cannot write to the original stderr any more (EPIPE: the reader of the log went
away) or fails otherwise, it puts the original stderr back as file descriptor 2 and goes on
reading (and dropping) what is still written into the pipe, so that no process blocks on a
full pipe. Python's fault handler writes to the original stderr, so the traceback of a crash
(segmentation fault) is not lost in the pipe of a filter that died with the process.
"""

from __future__ import annotations

import atexit
import faulthandler
import os
import threading

PATTERNS = (b"Fractional Error:",)
_state: dict = {}


def install() -> None:
    """Install the filter once per process tree (idempotent)."""
    if _state or os.environ.get("SWINGERSIM_NO_STDERR_FILTER"):
        return
    try:
        orig = os.dup(2)
        r, w = os.pipe()
        os.dup2(w, 2)
        os.close(w)
    except OSError:
        return
    _state.update(orig=orig, count=0, shown=False)
    try:  # crash tracebacks straight to the real stderr (all threads, child processes too)
        faulthandler.enable(file=orig, all_threads=True)
    except (OSError, ValueError, RuntimeError):
        pass

    def pump() -> None:
        buf = b""
        with os.fdopen(r, "rb", buffering=0) as fh:
            try:
                while True:
                    chunk = fh.read(65536)
                    if not chunk:
                        break
                    buf += chunk
                    *lines, buf = buf.split(b"\n")
                    out = []
                    for ln in lines:
                        if ln.startswith(PATTERNS):
                            _state["count"] += 1
                            if not _state["shown"]:
                                _state["shown"] = True
                                out.append(
                                    ln + b"  [REBOUNDx gr_full: 1PN iteration above machine "
                                    b"precision; further messages counted, not shown]"
                                )
                            continue
                        out.append(ln)
                    if out:
                        _write_all(orig, b"\n".join(out) + b"\n")
                if buf and not buf.startswith(PATTERNS):
                    _write_all(orig, buf)
            except Exception:  # noqa: BLE001 - EPIPE or anything else: never block the writers
                _state["failed"] = True
                try:
                    os.dup2(orig, 2)  # new output goes straight to the original stderr
                except OSError:
                    pass
                try:  # drain what processes (pool workers) still write into the pipe
                    while fh.read(65536):
                        pass
                except OSError:
                    pass

    t = threading.Thread(target=pump, name="swingersim-stderr-filter", daemon=True)
    t.start()
    _state["thread"] = t
    atexit.register(_finish)


def _write_all(fd: int, data: bytes) -> None:
    while data:
        n = os.write(fd, data)
        data = data[n:]


def _finish() -> None:
    orig = _state.get("orig")
    if orig is None:
        return
    try:
        os.dup2(orig, 2)  # closes this process's write end of the pipe
    except OSError:
        return
    t = _state.get("thread")
    if t is not None:
        t.join(timeout=2.0)
    n = _state.get("count", 0)
    if n > 1:
        try:
            os.write(
                orig, f"[swingersim] {n} REBOUNDx 'Fractional Error' messages suppressed\n".encode()
            )
        except OSError:
            pass
