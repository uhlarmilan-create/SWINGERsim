"""swingersim gui: one local web application for the inputs and the outputs of swingersim.

    swingersim gui                     # http://127.0.0.1:8765 in a browser
    swingersim gui --host 0.0.0.0      # reachable from the lab network

Five pages:

* Configuration: the drafts in ``jobs/drafts/`` - forms generated from the configuration
  schema (every key with its help text, default and immediate validation), new drafts, drafts
  opened from an example or another draft, saved under the same or a new name - and "What to
  compute" with the configuration: run, light curve, eclipse times (O - C) exact or as an
  observer, stability check, chaos map, pulsation periods, ensemble, scattering experiment;
  every ticked item becomes a job.
* Jobs: the job folder of the worker - drafts sent to the inbox, the inbox with the jobs to
  run (select and run), the running job with its progress, intermediate results and log
  (pause, continue, stop), and the finished, failed and cancelled jobs; the worker itself is
  started and stopped here.
* Outputs: every finished (or running) job as a result, grouped by its configuration: a card
  with facts read from the output and, for each kind of result, a page made for it (a run:
  system view, eclipses over the run, orbits over time, secular dynamics; an ensemble:
  observations and re-weighting, drawn vs. conditioned parameters, outcomes, members; a light
  curve, eclipse times, minima, a validation), the computations that can follow, and under
  Advanced the files, the log, the configuration, free plots and the 3D viewer
  (``gui/digest.py`` reads the outputs). Every file under ``output/`` can be opened too.
* O - C fit: the proposal of the fitting of measured minima.
* Tools: what belongs neither to a configuration nor to a result (validate, applegate,
  escape-theory, mesa-grid, mesa-reduce).

Only the Python standard library is used for the server; the plots use plotly.js, served from
the installed plotly package (extra ``viz``) or, without it, from the plotly CDN. Reading
parquet tables needs pandas + pyarrow (extra ``ensemble``); every ensemble also writes a CSV.
The server only reads and writes files inside the SWINGERsim folder.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import mimetypes
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import traceback
import webbrowser
from datetime import UTC, datetime
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, unquote, urlparse

import yaml

APP = Path(__file__).with_name("app.html")
PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.35.2.min.js"
MAX_ROWS = 200_000  # rows sent to the browser per column request
TEXT_LIMIT = 400_000  # characters of a text file shown


# --------------------------------------------------------------------------- helpers
class GuiError(ValueError):
    pass


def _inside(root: Path, rel: str) -> Path:
    """Path rel under root (swingersim.rootpath.inside_root); refused paths raise GuiError."""
    from ..rootpath import inside_root

    try:
        return inside_root(root, rel)
    except ValueError as exc:
        msg = str(exc)
        raise GuiError(
            f"path outside the SWINGERsim folder: {rel}" if msg.startswith("outside") else msg
        ) from exc


def _rel(root: Path, p: Path) -> str:
    return str(p.absolute().relative_to(root.absolute())).replace("\\", "/")


def _num(v):
    """JSON-safe number (NaN and infinities become None)."""
    if isinstance(v, float) and not math.isfinite(v):
        return None
    return v


def plotly_js_path() -> Path | None:
    try:
        import plotly
    except ImportError:
        return None
    p = Path(plotly.__file__).parent / "package_data" / "plotly.min.js"
    return p if p.exists() else None


# --------------------------------------------------------------------------- configuration
def config_schema() -> dict:
    from ..config import Config

    return Config.model_json_schema()


def list_configs(root: Path) -> list[str]:
    out = []
    for base in ("jobs/drafts", "configs", "jobs/submitted"):
        d = root / base
        if d.exists():
            out += [_rel(root, p) for p in sorted(d.rglob("*.yaml")) if "/." not in _rel(root, p)]
    return out


def load_config(root: Path, rel: str) -> dict:
    p = _inside(root, rel)
    text = p.read_text(encoding="utf-8")
    data = yaml.safe_load(text) or {}
    return {"path": rel, "text": text, "data": data}


def validate_config(data: dict) -> dict:
    """Validate a configuration dict; returns {ok, errors: [{loc, msg}], yaml}."""
    from pydantic import ValidationError

    from ..config import Config

    text = yaml.safe_dump(data, sort_keys=False, allow_unicode=True)
    if isinstance(data, dict) and "scatter" in data:  # a scattering experiment
        from ..scatter import ScatterConfig

        Model = ScatterConfig
    else:
        Model = Config
    try:
        Model.model_validate(data)
    except ValidationError as exc:
        errs = [{"loc": ".".join(str(x) for x in e["loc"]), "msg": e["msg"]} for e in exc.errors()]
        return {"ok": False, "errors": errs, "yaml": text}
    except (ValueError, TypeError) as exc:
        return {"ok": False, "errors": [{"loc": "", "msg": str(exc)}], "yaml": text}
    return {"ok": True, "errors": [], "yaml": text}


def parse_yaml(text: str) -> dict:
    try:
        data = yaml.safe_load(text) or {}
    except yaml.YAMLError as exc:
        raise GuiError(f"YAML error: {exc}") from exc
    if not isinstance(data, dict):
        raise GuiError("the configuration must be a mapping")
    return data


def save_config(root: Path, rel: str, text: str, overwrite: bool) -> dict:
    if not rel.endswith((".yaml", ".yml")):
        raise GuiError("the file name must end with .yaml")
    p = _inside(root, rel)
    if p.exists() and not overwrite:
        raise GuiError(f"{rel} exists; tick 'overwrite' to replace it")
    parse_yaml(text)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")
    return {"saved": rel}


# --------------------------------------------------------------------------- outputs
CSV_KINDS = (
    (("binary", "kind", "epoch", "t_min_d"), "minima"),
    (("binary", "kind", "t_mid_d"), "seasons"),
    (("binary", "kind", "epoch", "t_yr"), "timings"),
    (("t_day", "flux"), "lightcurve"),
)


def _csv_header(p: Path) -> list[str]:
    """Column names of a CSV file (comment lines starting with '#' before them are skipped)."""
    from ..tables import header

    return header(p)


_JSON_KIND: dict[str, tuple[int, int, str | None]] = {}  # path -> (mtime_ns, size, kind)


def _json_kind(p: Path) -> str | None:
    """emulator or json; parsed once per (path, mtime, size): an output folder can hold many
    large .json files and the Outputs page rescans it."""
    st = p.stat()
    hit = _JSON_KIND.get(str(p))
    if hit and hit[:2] == (st.st_mtime_ns, st.st_size):
        return hit[2]
    try:
        d = json.loads(p.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        kind = None
    else:
        kind = (
            "emulator"
            if isinstance(d, dict) and {"inputs", "targets", "metrics"} <= set(d)
            else "json"
        )
    if len(_JSON_KIND) > 5000:
        _JSON_KIND.clear()
    _JSON_KIND[str(p)] = (st.st_mtime_ns, st.st_size, kind)
    return kind


def classify(p: Path) -> str | None:
    """Kind of an output file or folder (None: not shown)."""
    if p.is_dir():
        names = {q.name for q in p.iterdir()}
        if "chaosmap.json" in names:
            return "chaosmap"
        if "scatter.yaml" in names:
            return "scatter"
        if "ensemble.yaml" in names or "rows.jsonl" in names or "summary.parquet" in names:
            return "ensemble"
        if "summary.csv" in names:
            return "ensemble"
        return None
    suf = p.suffix.lower()
    if suf == ".h5":
        return "run"
    if suf == ".csv":
        head = _csv_header(p)
        for cols, kind in CSV_KINDS:
            if tuple(head[: len(cols)]) == cols:
                return kind
        return "table"
    if suf == ".parquet":
        return "table"
    if suf == ".json":
        return _json_kind(p)
    if suf in (".png", ".jpg", ".jpeg", ".svg"):
        return "image"
    if suf in (".txt", ".md", ".log", ".yaml", ".yml"):
        return "text"
    if suf == ".html":
        return "html"
    return None


def scan_outputs(root: Path, base: str = "output", depth: int = 5) -> list[dict]:
    """All outputs under root/base: [{path, kind, size, mtime, parent}]. The files inside a
    recognised folder (ensemble, scattering, chaos map) are listed with that folder as parent."""
    top = root / base
    out: list[dict] = []
    if not top.exists():
        return out

    def walk(d: Path, level: int, parent: str | None):
        try:
            entries = sorted(d.iterdir(), key=lambda q: q.name)
        except OSError:
            return
        for q in entries:
            if q.name.startswith(".") or q.name == "__pycache__":
                continue
            try:
                kind = classify(q)
                st = q.stat()
            except OSError:
                continue
            rel = _rel(root, q)
            if q.is_dir():
                if kind:
                    out.append({"path": rel, "kind": kind, "mtime": st.st_mtime, "parent": parent})
                if level < depth:
                    walk(q, level + 1, rel if kind else parent)
            elif kind:
                out.append(
                    {
                        "path": rel,
                        "kind": kind,
                        "size": st.st_size,
                        "mtime": st.st_mtime,
                        "parent": parent,
                    }
                )

    walk(top, 0, None)
    return out


def _read_table(p: Path) -> dict:
    """Columns of a CSV, parquet or JSON-lines table as numpy arrays: float (NaN for an empty
    cell) when every value is a number, else an object array of the values."""
    if p.suffix == ".parquet":
        try:
            import pandas as pd
        except ImportError as exc:
            raise GuiError(
                "reading parquet needs pandas + pyarrow (pip install -e '.[ensemble]')"
            ) from exc
        df = pd.read_parquet(p)
        return {
            str(c): df[c].to_numpy(dtype=float)
            if df[c].dtype.kind in "iuf"
            else _as_array(df[c].tolist())
            for c in df.columns
        }
    if p.suffix == ".jsonl":
        rows = []
        with p.open(encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line:
                    try:
                        rows.append(json.loads(line))
                    except ValueError:
                        continue
        cols: dict[str, list] = {}
        for k, r in enumerate(rows):
            flat = _flatten(r)
            for c, v in flat.items():
                cols.setdefault(c, [None] * k).append(v)
            for col in cols.values():
                if len(col) < k + 1:
                    col.append(None)
        return {c: _as_array(v) for c, v in cols.items()}
    with p.open(newline="", encoding="utf-8", errors="replace") as f:
        from ..tables import data_lines

        reader = csv.reader(data_lines(f))
        head = [h.strip() for h in next(reader, [])]
        data = [[] for _ in head]
        for row in reader:
            if not row:
                continue
            for k in range(len(head)):
                data[k].append(_parse_cell(row[k]) if k < len(row) else None)
    return {h: _as_array(v) for h, v in zip(head, data, strict=True)}


def _is_number(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _as_array(values: list):
    """A column as a float array (None and NaN -> NaN) when all its values are numbers, else
    as an object array. A float array takes 8 bytes per value, a list of Python floats ~32."""
    import numpy as np

    if all(v is None or _is_number(v) for v in values):
        return np.array([math.nan if v is None else v for v in values], dtype=float)
    out = np.empty(len(values), dtype=object)
    out[:] = values
    return out


def _flatten(d: dict, prefix: str = "") -> dict:
    out = {}
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(_flatten(v, key + "."))
        elif isinstance(v, list):
            out[key] = json.dumps(v)
        else:
            out[key] = v
    return out


def _parse_cell(s: str):
    s = s.strip()
    if s == "":
        return None
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return float(s)
    except ValueError:
        pass
    if s in ("True", "true"):
        return True
    if s in ("False", "false"):
        return False
    return s


def _table_file(p: Path) -> Path:
    """The table of a folder output (ensemble, scattering, chaos map) or the file itself."""
    if p.is_file():
        return p
    for name in ("summary.parquet", "summary.csv", "rows.jsonl", "grid.jsonl"):
        q = p / name
        if q.exists():
            if name == "summary.parquet":
                try:
                    import pandas  # noqa: F401
                    import pyarrow  # noqa: F401
                except ImportError:
                    continue
            return q
    raise GuiError(f"no table found in {p}")


_TABLE_CACHE: dict[tuple, dict] = {}  # (file, mtime_ns, size) -> {"cols", "info"}
_LOCK = threading.Lock()
_LOAD_LOCKS: dict[str, threading.Lock] = {}  # one reader per file at a time
TABLES_CACHED = 3


def _table_entry(p: Path) -> dict:
    f = _table_file(p)
    st = f.stat()
    key = (str(f), st.st_mtime_ns, st.st_size)
    with _LOCK:
        hit = _TABLE_CACHE.get(key)
        if hit is not None:
            return hit
        load_lock = _LOAD_LOCKS.setdefault(str(f), threading.Lock())
    with load_lock:  # a second request for the same table waits for the first reader
        with _LOCK:
            hit = _TABLE_CACHE.get(key)
        if hit is not None:
            return hit
        cols = _read_table(f)
        entry = {"cols": cols, "info": None}
        with _LOCK:
            for k in [k for k in _TABLE_CACHE if k[0] == key[0]]:  # an older version of it
                del _TABLE_CACHE[k]
            while len(_TABLE_CACHE) >= TABLES_CACHED:  # the oldest entry goes
                del _TABLE_CACHE[next(iter(_TABLE_CACHE))]
            _TABLE_CACHE[key] = entry
    return entry


def table(p: Path) -> dict:
    """The columns of a table (numpy arrays), cached (TABLES_CACHED tables)."""
    return _table_entry(p)["cols"]


def _table_info(p: Path) -> list[dict]:
    entry = _table_entry(p)
    if entry["info"] is None:
        entry["info"] = [_col_info(c, v) for c, v in entry["cols"].items()]
    return entry["info"]


def _col_info(name: str, values) -> dict:
    import numpy as np

    values = np.asarray(values) if not isinstance(values, np.ndarray) else values
    info = {"name": name, "n": len(values)}
    if values.dtype != object:
        finite = values[np.isfinite(values)]
        if len(finite):
            info.update(type="number", min=float(finite.min()), max=float(finite.max()))
            return info
        info.update(type="category", n_cat=1 if len(values) else 0)
        if len(values):
            info["categories"] = {"None": len(values)}
        return info
    nums = [v for v in values if _is_number(v)]
    finite = [v for v in nums if math.isfinite(v)]
    if finite and len(nums) >= 0.5 * sum(v is not None for v in values):
        info.update(type="number", min=min(finite), max=max(finite))
    else:
        cats: dict[str, int] = {}
        for v in values:
            cats[str(v)] = cats.get(str(v), 0) + 1
        info.update(type="category", n_cat=len(cats))
        if len(cats) <= 30:
            info["categories"] = cats
    return info


PRESETS = {
    "ensemble": [
        {"title": "fates", "type": "bar", "x": "fate", "split": "direction"},
        {"title": "outer period", "type": "hist", "x": "orbits.out.P", "log": True},
        {"title": "eccentricity of A", "type": "hist", "x": "A_e_max"},
    ],
    "scatter": [{"title": "outcomes", "type": "bar", "x": "outcome"}],
    "chaosmap": [{"title": "MEGNO map", "type": "heat", "x": "x", "y": "y", "z": "megno"}],
    "minima": [{"title": "O − C", "type": "oc", "x": "t_min_d", "y": "o_c_d", "split": "binary"}],
    "seasons": [
        {"title": "depth", "type": "line", "x": "t_mid_d", "y": "depth_mag", "split": "binary"}
    ],
    "timings": [{"title": "O − C", "type": "oc", "x": "t_obs_yr", "y": "o_c_d", "split": "binary"}],
    "lightcurve": [{"title": "light curve", "type": "lc", "x": "t_day", "y": "flux"}],
}


def _dynamic_presets(kind: str, names: list[str]) -> list[dict]:
    """Ready-made plots for the columns an output actually has (an ensemble's columns depend
    on its topology and on what its configuration switched on)."""
    out: list[dict] = []
    has = set(names)
    if kind == "ensemble":
        orbits = [n[: -len("_e_max")] for n in names if n.endswith("_e_max")]
        for o in orbits:
            out.append({"title": f"max e of {o}", "type": "hist", "x": f"{o}_e_max"})
        for o in orbits:
            if f"{o}_log1me_min" in has:
                out.append(
                    {
                        "title": f"min log(1 − e) of {o}",
                        "type": "hist",
                        "x": f"{o}_log1me_min",
                        "split": "fate",
                    }
                )
        for n in names:
            if n.startswith("imut_") and n.endswith("_0") and n[:-2] + "_max" in has:
                pair = n[5:-2]
                out.append(
                    {
                        "title": f"i_mut {pair}: start vs max",
                        "type": "scatter",
                        "x": n,
                        "y": n[:-2] + "_max",
                        "split": "fate",
                    }
                )
        if "orbits.out.P" in has and orbits:
            out.append(
                {
                    "title": "fate vs outer period",
                    "type": "scatter",
                    "x": "orbits.out.P",
                    "y": f"{orbits[0]}_e_max",
                    "split": "fate",
                    "log": True,
                }
            )
        if "event_time_yr" in has:
            out.append(
                {
                    "title": "time of the event",
                    "type": "hist",
                    "x": "event_time_yr",
                    "split": "fate",
                    "log": True,
                }
            )
        for n in names:
            if n.endswith("_eclipsing0"):
                out.append({"title": f"{n[:-11]} eclipsing now", "type": "bar", "x": n})
        for n in ("eclipse_fraction_both",) + tuple(
            x for x in names if x.endswith("_eclipse_fraction")
        ):
            if n in has:
                out.append({"title": n.replace("_", " "), "type": "hist", "x": n})
        for n in names:
            if n.endswith("_eclipse_change_yr"):
                out.append(
                    {
                        "title": f"{n[:-18]}: eclipses end/begin after",
                        "type": "hist",
                        "x": n,
                        "log": True,
                    }
                )
        for n in names:
            if n.endswith("_roche_fill_max"):
                out.append({"title": f"Roche-lobe filling of {n[:-15]}", "type": "hist", "x": n})
        for n, title in (
            ("ma01_min", "stability margin (Mardling-Aarseth)"),
            ("v22_min", "stability (Vynatheya+22)"),
            ("max_energy_error", "energy error"),
            ("wall_s", "wall time per member"),
        ):
            if n in has:
                out.append(
                    {
                        "title": title,
                        "type": "hist",
                        "x": n,
                        "log": n in ("max_energy_error", "wall_s"),
                    }
                )
        if "final_hierarchy" in has:
            out.append({"title": "final hierarchy", "type": "bar", "x": "final_hierarchy"})
    elif kind == "scatter":
        if "b_au" in has:
            out.append(
                {
                    "title": "outcome vs impact parameter",
                    "type": "hist",
                    "x": "b_au",
                    "split": "outcome",
                }
            )
        if "t_end_yr" in has:
            out.append(
                {
                    "title": "duration of the encounter",
                    "type": "hist",
                    "x": "t_end_yr",
                    "split": "outcome",
                    "log": True,
                }
            )
    elif kind == "minima":
        for y, title in (
            ("o_c_geom_d", "O − C of the geometric conjunction"),
            ("depth_mag", "depth"),
            ("width_d", "width"),
        ):
            out.append({"title": title, "type": "oc", "x": "t_min_d", "y": y, "split": "kind"})
    elif kind == "seasons":
        for y in ("width_d", "gamma", "rms_mag"):
            out.append({"title": y, "type": "line", "x": "t_mid_d", "y": y, "split": "kind"})
    elif kind == "timings":
        out.append(
            {"title": "O − C by kind", "type": "oc", "x": "t_obs_yr", "y": "o_c_d", "split": "kind"}
        )
    return out


def output_info(root: Path, rel: str) -> dict:
    """Description of one output for the browser: kind, files, columns, presets, text."""
    p = _inside(root, rel)
    if not p.exists():
        raise GuiError(f"{rel} does not exist")
    kind = classify(p) or "file"
    info: dict = {"path": rel, "kind": kind}
    if p.is_dir():
        info["files"] = [
            {"path": _rel(root, q), "kind": classify(q) or "file", "size": q.stat().st_size}
            for q in sorted(p.iterdir())
            if q.is_file()
        ]
        for name in ("README.md", "report_conditioned.txt", "summary.txt"):
            if (p / name).exists():
                info.setdefault("text", (p / name).read_text(encoding="utf-8")[:TEXT_LIMIT])
                info["text_file"] = name
                break
    if kind in (
        "ensemble",
        "scatter",
        "chaosmap",
        "minima",
        "seasons",
        "timings",
        "lightcurve",
        "table",
    ):
        t = table(p)
        info["n_rows"] = len(next(iter(t.values()), []))
        info["columns"] = _table_info(p)
        names = set(t)
        info["presets"] = [
            pr
            for pr in PRESETS.get(kind, []) + _dynamic_presets(kind, list(t))
            if all(pr.get(k) is None or pr[k] in names for k in ("x", "y", "z", "split"))
        ]
        if kind == "chaosmap":
            meta = json.loads((p / "chaosmap.json").read_text(encoding="utf-8"))
            info["meta"] = {k: meta[k] for k in meta if k != "config"}
    elif kind == "emulator":
        info["emulator"] = json.loads(p.read_text(encoding="utf-8"))
    elif kind in ("text", "json"):
        info["text"] = p.read_text(encoding="utf-8", errors="replace")[:TEXT_LIMIT]
    elif kind == "run":
        info["run"] = run_summary(p)
    return info


def columns(root: Path, rel: str, cols: list[str]) -> dict:
    """The values of some columns of a table (the first MAX_ROWS rows; "_truncated" says
    whether there were more, "_n_rows" how many)."""
    t = table(_inside(root, rel))
    n = len(next(iter(t.values()), []))
    out: dict = {"_n_rows": n, "_truncated": n > MAX_ROWS}
    for c in cols:
        if c not in t:
            raise GuiError(f"no column {c}")
        out[c] = [_num(v) for v in t[c][:MAX_ROWS].tolist()]
    return out


def run_summary(p: Path) -> dict:
    import h5py

    with h5py.File(p, "r") as f:
        attrs = {k: (v.decode() if isinstance(v, bytes) else v) for k, v in f.attrs.items()}
        out = {"groups": list(f.keys())}
        if "swingersim_version" not in attrs and "quadsim_version" in attrs:
            attrs["swingersim_version"] = attrs[
                "quadsim_version"
            ]  # files written before the rename
        for k in ("event", "fate", "backend", "swingersim_version"):
            if k in attrs:
                out[k] = str(attrs[k])
        if "time_yr" in f:
            t = f["time_yr"][()]
            out["t_range_yr"] = [float(t[0]), float(t[-1])] if len(t) else None
            out["n_output"] = len(t)
    return out


def run_viewer(root: Path, rel: str) -> str:
    """Path (relative) of the interactive viewer page of a run: the <run>_viewer.html the run
    wrote when it is newer than the run file, else one built into output/.gui_cache (keyed by
    the run file and the SWINGERsim version)."""
    from .. import __version__
    from ..runproducts import product_paths

    p = _inside(root, rel)
    st = p.stat()
    own = product_paths(p)["viewer"]
    try:
        if own.is_file() and own.stat().st_mtime_ns >= st.st_mtime_ns:
            return _rel(root, own)
    except (OSError, ValueError):
        pass
    key = f"{p}:{st.st_mtime_ns}:{st.st_size}:{__version__}"
    h = hashlib.sha1(key.encode()).hexdigest()[:16]
    cache = root / "output" / ".gui_cache"
    html = cache / f"{p.stem}-{h}.html"
    if not html.exists():
        try:
            from ..viewer import write_html
        except ImportError as exc:
            raise GuiError("the run viewer needs plotly: pip install -e '.[viz]'") from exc
        write_html(p, html, offline=True)
    return _rel(root, html)


CACHE_DAYS = 30  # entries of output/.gui_cache not touched for this long are removed


def prune_cache(root: Path, days: float = CACHE_DAYS) -> int:
    """Remove the files of output/.gui_cache (viewer pages, model curves, cards) that were not
    written or used for ``days``. Returns how many were removed."""
    top = root / "output" / ".gui_cache"
    if not top.is_dir():
        return 0
    limit = time.time() - days * 86400.0
    n = 0
    for q in top.rglob("*"):
        try:
            st = q.stat()
            if q.is_file() and max(st.st_mtime, st.st_atime) < limit:
                q.unlink()
                n += 1
        except OSError:
            continue
    return n


# --------------------------------------------------------------------------- jobs
JOB_FOLDERS = ("inbox", "running", "done", "failed", "cancel")
NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.+-]*$")
# where a job made in the GUI writes: output/<job id>/<name>; "" = the folder itself
OUTPUT_NAMES = {
    "run": "run.h5",
    "ensemble": "",
    "scatter": "",
    "chaosmap": "",
    "timings": "timings.csv",
    "minima": "minima.csv",
    "lightcurve": "lightcurve.csv",
    "emulator-train": "model.pkl",
    "fit": "",
}
QUICK_TIMEOUT_S = 300  # commands run directly in the GUI


def _jobdir(root: Path) -> Path:
    jd = root / "jobs"
    for sub in ("drafts", "submitted", "inbox", "running", "done", "failed", "cancel", "pause"):
        (jd / sub).mkdir(parents=True, exist_ok=True)
    return jd


def _job_id(name: str) -> str:
    return name.split(".yaml")[0].split(".yml")[0]


def _pid_alive(pid) -> bool:
    try:
        os.kill(int(pid), 0)
    except (OSError, ValueError, TypeError):
        return False
    try:  # never mistake an unrelated process with a recycled pid for the worker
        return b"swingersim" in Path(f"/proc/{int(pid)}/cmdline").read_bytes()
    except OSError:
        return True  # no /proc (not Linux): trust the pid


def worker_state(root: Path) -> dict:
    jd = root / "jobs"
    w = {}
    if (jd / "worker.json").exists():
        try:
            w = json.loads((jd / "worker.json").read_text(encoding="utf-8"))
        except ValueError:
            w = {}
    alive = bool(w) and not w.get("stopped") and _pid_alive(w.get("pid"))
    w["alive"] = alive
    w["stop_requested"] = (jd / "stop").exists()
    return w


START_WAIT_S = 1.5  # how long a start request waits for the worker's first heartbeat


def start_worker(root: Path) -> dict:
    jd = _jobdir(root)
    st = worker_state(root)
    if st["alive"]:
        return {"started": False, "worker": st}
    (jd / "stop").unlink(missing_ok=True)
    with open(jd / "worker.out", "a", encoding="utf-8") as log:  # the child keeps its own copy
        proc = subprocess.Popen(
            [
                sys.executable,
                "-m",
                "swingersim",
                "worker",
                str(jd),
                "--root",
                str(root),
                "--poll",
                "3",
            ],
            cwd=root,
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,  # the worker outlives the GUI
        )
    # a short look for the first heartbeat; the Jobs page polls the state anyway, so the
    # request thread is not held for long
    deadline = time.monotonic() + START_WAIT_S
    st = worker_state(root)
    while not (st.get("pid") == proc.pid and st["alive"]) and time.monotonic() < deadline:
        if proc.poll() is not None:  # it ended at once (an error): worker.out says why
            break
        time.sleep(0.1)
        st = worker_state(root)
    return {
        "started": proc.poll() is None,
        "pid": proc.pid,
        "returncode": proc.poll(),
        "worker": st,
    }


def stop_worker(root: Path) -> dict:
    jd = _jobdir(root)
    (jd / "stop").write_text("")
    return {"stop_requested": True}


def _read_job(p: Path) -> dict:
    try:
        job = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
        return job if isinstance(job, dict) else {}
    except (OSError, yaml.YAMLError):
        return {}


def _status(folder: Path, job_id: str) -> dict:
    f = folder / f"{job_id}.status.json"
    if f.exists():
        try:
            return json.loads(f.read_text(encoding="utf-8"))
        except ValueError:
            return {}
    return {}


def job_report(root: Path, job: dict) -> str | None:
    """The HTML report a finished job wrote (relative path), if it exists."""
    args = job.get("args") or {}
    out = job.get("output")
    cands = []
    if args.get("html"):
        cands.append(args["html"])
    if job.get("command") == "validate":
        cands.append("output/validation/report.html")
    if out:
        q = root / out
        cands += (
            [str(Path(out) / "report.html"), str(Path(out).with_suffix(".html"))]
            if (q.is_dir() or not q.suffix)
            else [str(Path(out).with_suffix(".html"))]
        )
    elif args.get("plot"):
        cands.append(str(Path(args["plot"]).with_suffix(".html")))
    for c in cands:
        try:
            q = _inside(root, str(c))
        except GuiError:
            continue
        if q.is_file():
            return _rel(root, q)
    return None


def jobs(root: Path) -> dict:
    from ..worker import job_cores

    jd = _jobdir(root)
    out: dict = {"folders": {}, "worker": worker_state(root), "drafts": list_drafts(root)}
    total = int(out["worker"].get("cores") or os.cpu_count() or 1)
    out["cores"] = total
    for sub in JOB_FOLDERS:
        d = jd / sub
        items = []
        for q in sorted(d.iterdir()):
            if not (q.suffix in (".yaml", ".yml") or q.name.endswith((".yaml.hold", ".yml.hold"))):
                continue
            jid = _job_id(q.name)
            job = _read_job(q)
            item = {
                "id": jid,
                "file": q.name,
                "path": _rel(root, q),
                "held": q.name.endswith(".hold"),
                "mtime": q.stat().st_mtime,
                **{k: job.get(k) for k in ("command", "config", "output", "note", "args")},
                "status": _status(d, jid),
            }
            try:
                item["cores"] = job_cores(job, total, root)
            except (ValueError, TypeError):
                item["cores"] = 1
            if sub in ("done", "failed"):
                item["report"] = job_report(root, job)
            if sub == "running":
                item["paused_requested"] = (jd / "pause" / jid).exists()
                item["cancel_requested"] = (jd / "cancel" / jid).exists()
                log = d / f"{jid}.log"
                if log.exists():
                    item["log_tail"] = _tail(log, 25)
            items.append(item)
        out["folders"][sub] = items
    return out


def _tail(p: Path, n: int = 200, block: int = 65536) -> str:
    """The last n lines of a text file, read from its end (a log can be large)."""
    with p.open("rb") as f:
        f.seek(0, os.SEEK_END)
        end = f.tell()
        pos, data = end, b""
        while pos > 0 and data.count(b"\n") <= n:
            step = min(block, pos)
            pos -= step
            f.seek(pos)
            data = f.read(step) + data
    lines = data.decode("utf-8", errors="replace").splitlines()
    if pos > 0:  # the first line may be cut
        lines = lines[1:]
    return "\n".join(lines[-n:])


def log_tail(root: Path, rel: str, n: int = 200) -> str:
    return _tail(_inside(root, rel), n)


_NUMBER_LOCK = threading.RLock()  # job numbers: one submission at a time


def _next_number(jd: Path) -> int:
    best = 0
    for sub in JOB_FOLDERS + ("submitted",):
        for q in (jd / sub).glob("*"):
            m = re.match(r"^(\d+)_", q.name)
            if m:
                best = max(best, int(m.group(1)))
    return best + 1


def safe_name(name: str, what: str = "name") -> str:
    """A job id, draft or grid name: letters, digits and _.+- only (no path separators, no
    '..'), so that it cannot point outside its folder."""
    name = str(name)
    if not NAME_RE.match(name) or ".." in name:
        raise GuiError(f"invalid {what} '{name}'")
    return name


def _find_job(jd: Path, job_id: str) -> tuple[str, Path] | tuple[None, None]:
    safe_name(job_id, "job id")
    for sub in JOB_FOLDERS:
        for suffix in (".yaml", ".yaml.hold", ".yml", ".yml.hold"):
            q = jd / sub / f"{job_id}{suffix}"
            if q.exists():
                return sub, q
    return None, None


# drafts ----------------------------------------------------------------------
def list_drafts(root: Path) -> list[dict]:
    d = _jobdir(root) / "drafts"
    out = []
    for q in sorted(d.glob("*.yaml")):
        data = _read_job(q)
        out.append(
            {
                "name": q.stem,
                "path": _rel(root, q),
                "mtime": q.stat().st_mtime,
                "kind": data["command"]  # a complete job file (command, config, args)
                if isinstance(data.get("command"), str)
                else "scatter"
                if "scatter" in data
                else "ensemble"
                if "ensemble" in data
                else "run",
                "topology": data.get("topology"),
            }
        )
    return out


def _draft_path(root: Path, name: str) -> Path:
    name = name.strip().removesuffix(".yaml").removesuffix(".yml")
    if not NAME_RE.match(name):
        raise GuiError("a draft name may contain letters, digits and _ . + - (no spaces, no /)")
    return _jobdir(root) / "drafts" / f"{name}.yaml"


def save_draft(
    root: Path, name: str, text: str, overwrite: bool, source: str | None = None
) -> dict:
    """Save a draft. ``source``: the file the configuration was opened from; its relative paths
    (minima table, light curves, stellar table, ...) are rewritten relative to jobs/drafts/
    (the result lists them under "rebased" and gives the rewritten "data")."""
    p = _draft_path(root, name)
    if p.exists() and not overwrite:
        raise GuiError(f"draft {p.stem} exists; save it under another name or allow overwriting")
    data = parse_yaml(text)
    changed: list[str] = []
    if source and not isinstance(data.get("command"), str):  # a configuration, not a job file
        try:
            src_dir = _inside(root, source).parent
        except GuiError:
            src_dir = None
        if src_dir is not None:
            changed = rebase_paths(data, src_dir, p.parent)
            if changed:
                text = yaml.safe_dump(data, sort_keys=False, allow_unicode=True)
    p.write_text(text, encoding="utf-8")
    out = {"saved": _rel(root, p), "name": p.stem, "rebased": changed}
    if changed:
        out["data"] = data
    return out


def delete_draft(root: Path, name: str) -> dict:
    p = _draft_path(root, name)
    if not p.exists():
        raise GuiError(f"no draft {name}")
    trash = _jobdir(root) / "drafts" / ".deleted"
    trash.mkdir(exist_ok=True)
    shutil.move(str(p), trash / f"{p.stem}.{int(time.time())}.yaml")
    return {"deleted": p.stem}


# fields of a configuration that hold paths relative to the configuration file ("*": every
# item of a list); a copy of a configuration elsewhere rewrites them (rebase_paths)
PATH_FIELDS = (
    ("stellar_model", "file"),
    ("stellar_model", "entry"),  # 'file.py:function' (a module name is not a path)
    ("fit", "minima"),
    ("fit", "rv", "file"),
    ("fit", "lc", "*", "file"),
    ("fit", "inclination", "*", "file"),
    ("fit", "output", "dir"),
)


def _rebase_one(v, src_dir: Path, dest_dir: Path, entry: bool = False):
    if not isinstance(v, str) or not v.strip():
        return v
    func = ""
    if entry:
        target, _, func = v.rpartition(":")
        if not target.endswith(".py"):
            return v
        v, func = target, ":" + func
    if os.path.isabs(v):
        return v + func
    new = os.path.relpath(os.path.normpath(Path(src_dir) / v), Path(dest_dir))
    return new.replace("\\", "/") + func


def rebase_paths(data: dict, src_dir: Path, dest_dir: Path) -> list[str]:
    """Rewrite the relative paths of a configuration dict (PATH_FIELDS) written for a file in
    src_dir so that they hold for a file in dest_dir. Returns the fields that changed."""
    changed: list[str] = []
    if Path(src_dir).absolute() == Path(dest_dir).absolute() or not isinstance(data, dict):
        return changed

    def visit(node, keys, label):
        k = keys[0]
        if k == "*":
            if isinstance(node, list):
                for i, item in enumerate(node):
                    visit(item, keys[1:], f"{label}[{i}]")
            return
        if not isinstance(node, dict) or k not in node:
            return
        here = f"{label}.{k}" if label else k
        if len(keys) > 1:
            visit(node[k], keys[1:], here)
            return
        new = _rebase_one(node[k], src_dir, dest_dir, entry=k == "entry")
        if new != node[k]:
            node[k] = new
            changed.append(here)

    for keys in PATH_FIELDS:
        visit(data, keys, "")
    return changed


def _freeze_config(src: Path, dest: Path) -> None:
    """Copy a configuration; the paths inside it that are relative to its folder (minima table,
    radial velocities, light curves, inclinations, stellar table, fit output) are rewritten
    relative to the new place."""
    data = yaml.safe_load(src.read_text(encoding="utf-8")) or {}
    if isinstance(data, dict) and not isinstance(data.get("command"), str):
        rebase_paths(data, src.parent, dest.parent)
    dest.write_text(
        f"# frozen copy of {src.name} ({time.strftime('%Y-%m-%d %H:%M')})\n"
        + yaml.safe_dump(data, sort_keys=False, allow_unicode=True),
        encoding="utf-8",
    )


def _all_cores(cmd: str, args: dict) -> None:
    """A job sent from the GUI uses every core of the computer: commands with a
    --processes option get processes 0 (all cores) unless it was set."""
    from ..cli import command_catalog

    c = command_catalog().get(cmd, {})
    if "processes" not in args and any(a["dest"] == "processes" for a in c.get("arguments", [])):
        args["processes"] = 0


def hold_others(root: Path, keep: list[str]) -> list[str]:
    """Put every queued job of the inbox on hold except ``keep``: a "run" in the GUI runs
    only the jobs the user chose, never other jobs that happen to be queued."""
    jd = _jobdir(root)
    held = []
    for q in sorted((jd / "inbox").glob("*.yaml")):
        jid = _job_id(q.name)
        if jid not in keep:
            q.rename(q.with_name(q.name + ".hold"))
            held.append(jid)
    return held


def submit(root: Path, drafts: list[str], *args, **kwargs) -> dict:
    """Make inbox jobs of drafts (see _submit); one submission at a time, so that two requests
    never take the same job number."""
    with _NUMBER_LOCK:
        return _submit(root, drafts, *args, **kwargs)


def submit_job(root: Path, job: dict, queue: bool = False, name: str = "job") -> dict:
    """A job given directly (from the Tools page or a dialog): checked and written to the
    inbox (see _submit_job); one submission at a time."""
    with _NUMBER_LOCK:
        return _submit_job(root, job, queue, name)


def _submit(
    root: Path,
    drafts: list[str],
    command: str | None = None,
    args: dict | None = None,
    output: str | None = None,
    note: str = "",
    queue: bool = False,
    suffix: str = "",
) -> dict:
    """Make inbox jobs of drafts: a frozen copy of each draft in jobs/submitted/<id>.yaml and
    jobs/inbox/<id>.yaml.hold (or .yaml with queue=True). ``suffix`` is added to the job name
    (``012_vw_lmi_lightcurve``): the computations asked on the Configuration page."""
    from ..worker import JobError, build_command

    jd = _jobdir(root)
    made = []
    for name in drafts:
        src = _draft_path(root, name)
        if not src.exists():
            raise GuiError(f"no draft {name}")
        data = _read_job(src)
        if isinstance(data.get("command"), str) and command is None:
            # the draft is a complete job file (command, config, output, args): queue it as it
            # is, not as the configuration of a "run"
            made += submit_job(root, data, queue=queue, name=src.stem)["jobs"]
            continue
        cmd = command or (
            "scatter" if "scatter" in data else "ensemble" if "ensemble" in data else "run"
        )
        stem = re.sub(r"^\d+_", "", src.stem) or src.stem
        if suffix:
            stem += "_" + re.sub(r"[^A-Za-z0-9_.+-]", "_", suffix)
        jid = f"{_next_number(jd):03d}_{stem}"
        frozen = jd / "submitted" / f"{jid}.yaml"
        _freeze_config(src, frozen)
        job: dict = {"command": cmd, "config": _rel(root, frozen)}
        out = output
        if out is None and cmd in OUTPUT_NAMES:
            out = f"output/{jid}" + (f"/{OUTPUT_NAMES[cmd]}" if OUTPUT_NAMES[cmd] else "")
        if out:
            job["output"] = out.replace("{id}", jid)
        a = {
            k: (v.replace("{id}", jid) if isinstance(v, str) else v)
            for k, v in (args or {}).items()
        }
        if cmd in ("ensemble", "scatter") and "resume" not in a:
            a["resume"] = True
        _all_cores(cmd, a)
        if a:
            job["args"] = a
        job["note"] = note or f"from draft {src.stem}"
        try:
            build_command(job, root)
        except JobError as exc:
            frozen.unlink(missing_ok=True)
            raise GuiError(str(exc)) from exc
        target = jd / "inbox" / (f"{jid}.yaml" if queue else f"{jid}.yaml.hold")
        target.write_text(
            f"# made by swingersim gui from jobs/drafts/{src.name}\n"
            + yaml.safe_dump(job, sort_keys=False, allow_unicode=True),
            encoding="utf-8",
        )
        made.append(jid)
    return {"jobs": made}


def _submit_job(root: Path, job: dict, queue: bool = False, name: str = "job") -> dict:
    """A job given directly: checked and written to the inbox."""
    from ..worker import JobError, build_command

    jd = _jobdir(root)
    stem = re.sub(r"[^A-Za-z0-9_.+-]", "_", name) or "job"
    stem = re.sub(r"^\d+_", "", stem) or stem  # "060_validate" -> 061_validate, not 061_060_...
    jid = f"{_next_number(jd):03d}_{stem}"
    job = json.loads(json.dumps(job).replace("{id}", jid))
    args = dict(job.get("args") or {})
    _all_cores(job.get("command", ""), args)
    if args:
        job["args"] = args
    cfg = job.get("config")
    if isinstance(cfg, str) and cfg.startswith("jobs/drafts/"):  # a draft is frozen
        src = _inside(root, cfg)
        if src.exists():
            frozen = jd / "submitted" / f"{jid}.yaml"
            _freeze_config(src, frozen)
            job["config"] = _rel(root, frozen)
    try:
        build_command(job, root)
    except JobError as exc:
        raise GuiError(str(exc)) from exc
    target = jd / "inbox" / (f"{jid}.yaml" if queue else f"{jid}.yaml.hold")
    target.write_text(
        "# made by swingersim gui\n" + yaml.safe_dump(job, sort_keys=False), encoding="utf-8"
    )
    return {"jobs": [jid]}


def run_now(root: Path, job: dict) -> dict:
    """Run a (short) command at once and return what it printed."""
    from ..worker import JobError, build_command

    try:
        argv, _ = build_command(job, root)
    except JobError as exc:
        raise GuiError(str(exc)) from exc
    t0 = time.time()
    try:
        p = subprocess.run(
            argv, cwd=root, capture_output=True, text=True, timeout=QUICK_TIMEOUT_S, check=False
        )
    except subprocess.TimeoutExpired as exc:
        raise GuiError(
            f"the command did not finish in {QUICK_TIMEOUT_S} s - send it to the inbox instead"
        ) from exc
    return {
        "command_line": " ".join(argv[2:]),
        "returncode": p.returncode,
        "stdout": p.stdout[-TEXT_LIMIT:],
        "stderr": p.stderr[-20000:],
        "seconds": round(time.time() - t0, 1),
    }


def job_action(root: Path, ids: list[str], action: str) -> dict:
    """Apply an action to jobs. A job that is not where the action needs it (e.g. a selected
    inbox job that has meanwhile run) is skipped and reported; an error only when none fit."""
    if action not in ("queue", "hold", "remove", "pause", "resume", "cancel", "again"):
        raise GuiError(f"unknown action {action}")
    jd = _jobdir(root)
    done, skipped = [], []
    for jid in ids:
        try:
            _job_action_one(jd, jid, action)
            done.append(jid)
        except GuiError as exc:
            skipped.append({"id": jid, "reason": str(exc)})
    if skipped and not done:
        raise GuiError("; ".join(x["reason"] for x in skipped))
    return {"action": action, "jobs": done, "skipped": skipped}


def _job_action_one(jd: Path, jid: str, action: str) -> None:
    where, p = _find_job(jd, jid)
    if p is None:
        raise GuiError(f"no job {jid}")
    if action in ("queue", "hold"):
        if where != "inbox":
            raise GuiError(f"{jid} is not in the inbox")
        base = jd / "inbox" / f"{jid}.yaml"
        if action == "queue" and p.name.endswith(".hold"):
            p.rename(base)
        elif action == "hold" and not p.name.endswith(".hold"):
            p.rename(base.with_name(base.name + ".hold"))
    elif action == "remove":
        if where != "inbox":
            raise GuiError(f"{jid} is not in the inbox (stop a running job instead)")
        shutil.move(str(p), jd / "cancel" / f"{jid}.yaml")
        (jd / "cancel" / f"{jid}.status.json").write_text(
            json.dumps(
                {
                    "id": jid,
                    "state": "removed",
                    # UTC like the worker's times, so that the lists sort them together
                    "finished": datetime.now(UTC).isoformat(timespec="seconds"),
                }
            )
        )
    elif action in ("pause", "resume", "cancel"):
        if where != "running":
            raise GuiError(f"{jid} is not running")
        if action == "pause":
            (jd / "pause" / jid).write_text("")
        elif action == "resume":
            (jd / "pause" / jid).unlink(missing_ok=True)
        else:
            (jd / "cancel" / jid).write_text("")
    elif action == "again":
        if where not in ("done", "failed", "cancel"):
            raise GuiError(f"{jid} has not finished")
        target = jd / "inbox" / f"{jid}.yaml.hold"
        if target.exists() or (jd / "inbox" / f"{jid}.yaml").exists():
            raise GuiError(f"{jid} is already in the inbox")
        # the job leaves its folder; the files of the previous attempt are kept aside
        prev = p.parent / ".previous"
        prev.mkdir(exist_ok=True)
        for f in p.parent.glob(f"{jid}.*"):
            if f != p:
                shutil.move(str(f), prev / f.name)
        shutil.move(str(p), target)


def job_products(root: Path, job_id: str) -> dict:
    """A job and everything it wrote: its output (file or folder), the other files named in
    its arguments, the files of its output/<id>/ folder, its configuration and log."""
    jd = _jobdir(root)
    where, p = _find_job(jd, job_id)
    if p is None:
        raise GuiError(f"no job {job_id}")
    job = _read_job(p)
    st = _status(jd / where, job_id)
    files: list[str] = []

    def add(rel):
        if not rel or not isinstance(rel, str):
            return
        try:
            q = _inside(root, rel)
        except GuiError:
            return
        if q.exists() and _rel(root, q) not in files:
            files.append(_rel(root, q))

    out = job.get("output") or st.get("output")
    add(out)
    add(job_report(root, job))
    if isinstance(out, str) and out.endswith(".h5"):  # the standard products of a run
        from ..runproducts import product_paths

        for q in product_paths(Path(out)).values():
            add(q.as_posix())
    for v in (job.get("args") or {}).values():
        for x in v if isinstance(v, list) else [v]:
            if isinstance(x, str) and x.startswith("output/"):
                add(x)
    if out:
        q = _inside(root, out)
        folder = q if q.is_dir() else q.parent
        if folder.name == job_id or q.is_dir():
            for f in sorted(folder.iterdir()) if folder.exists() else []:
                if f.is_file() and not f.name.startswith("."):
                    add(_rel(root, f))
    partial = None
    if out and where == "running":
        pq = _inside(root, out)
        pq = pq.with_name(pq.name + ".partial.json")
        if pq.exists():
            partial = _rel(root, pq)
    products = []
    for f in files:
        q = root / f
        products.append(
            {
                "path": f,
                "kind": classify(q) or "file",
                "size": q.stat().st_size if q.is_file() else None,
                "main": f == out,
            }
        )
    log = jd / where / f"{job_id}.log"
    return {
        "id": job_id,
        "folder": where,
        "job": job,
        "job_file": _rel(root, p),
        "status": st,
        "products": products,
        "partial": partial,
        "log": _rel(root, log) if log.exists() else None,
        "config": job.get("config"),
    }


def finished_jobs(root: Path) -> list[dict]:
    """Done jobs (newest first), then the running one, for the Outputs page."""
    jd = _jobdir(root)
    out = []
    for sub in ("running", "done", "failed", "cancel"):
        for q in (jd / sub).glob("*.y*ml"):
            jid = _job_id(q.name)
            st = _status(jd / sub, jid)
            job = _read_job(q)
            out.append(
                {
                    "id": jid,
                    "folder": sub,
                    "command": job.get("command"),
                    "output": job.get("output"),
                    "note": job.get("note"),
                    "finished": st.get("finished") or st.get("started"),
                    "state": st.get("state", sub),
                }
            )
    # the running job first, then the others newest first
    out.sort(key=lambda d: (d["folder"] == "running", d["finished"] or ""), reverse=True)
    return out


# time series of a run --------------------------------------------------------
ANGLES = {"inc", "Omega", "omega", "M", "inc_inv", "Omega_inv"}
SERIES_LIMIT = 20_000


def _series_meta(name: str) -> dict:
    parts = name.split("/")
    leaf = parts[-1]
    unit = ""
    if parts[0] == "mutual_inclination" or (parts[0] == "orbits" and leaf in ANGLES):
        unit = "deg"
    elif parts[0] == "orbits" and leaf == "a":
        unit = "AU"
    return {
        "name": name,
        "group": parts[0],
        "item": "/".join(parts[1:-1]) if len(parts) > 2 else "",
        "quantity": leaf,
        "unit": unit,
    }


def h5_index(root: Path, rel: str) -> dict:
    """The time series of a run (every 1-D dataset as long as time_yr), or of the
    intermediate-results file of a running N-body run."""
    p = _inside(root, rel)
    if p.suffix == ".json":
        d = json.loads(p.read_text(encoding="utf-8"))
        names = [f"orbits/{o}/{k}" for o, el in d["orbits"].items() for k in el]
        names += [f"mutual_inclination/{k}" for k in d["mutual_inclination"]]
        names += ["energy_error"]
        return {
            "n": len(d["time_yr"]),
            "series": [_series_meta(n) for n in names],
            "partial": True,
            "t_end_yr": d.get("t_end_yr"),
            "t_yr": d["time_yr"][-1] if d["time_yr"] else 0.0,
            "written": d.get("written"),
        }
    import h5py

    out = []
    with h5py.File(p, "r") as f:
        n = len(f["time_yr"]) if "time_yr" in f else 0

        def visit(name, obj):
            if (
                isinstance(obj, h5py.Dataset)
                and obj.ndim == 1
                and len(obj) == n
                and name != "time_yr"
                and obj.dtype.kind in "fiub"
            ):
                out.append(_series_meta(name))

        f.visititems(visit)
        attrs = {
            k: (v.decode() if isinstance(v, bytes) else str(v))
            for k, v in f.attrs.items()
            if k not in ("config_yaml",)
        }
    return {"n": n, "series": out, "attrs": attrs, "partial": False}


def h5_series(root: Path, rel: str, names: list[str]) -> dict:
    p = _inside(root, rel)
    if p.suffix == ".json":
        d = json.loads(p.read_text(encoding="utf-8"))

        def get(name):
            parts = name.split("/")
            if parts[0] == "orbits":
                return d["orbits"][parts[1]][parts[2]]
            if parts[0] == "mutual_inclination":
                return d["mutual_inclination"][parts[1]]
            return d[name]

        t = d["time_yr"]
        out = {"time_yr": t}
        for nm in names:
            v = get(nm)
            meta = _series_meta(nm)
            out[nm] = [
                None if x is None else (math.degrees(x) if meta["unit"] == "deg" else x) for x in v
            ]
        return out
    import h5py
    import numpy as np

    from .digest import _decimate

    with h5py.File(p, "r") as f:
        t = f["time_yr"][()]
        vals = {}
        for nm in names:
            if nm not in f or not hasattr(f[nm], "shape"):
                raise GuiError(f"no series {nm}")
            v = np.asarray(f[nm][()], float)  # whole: the extremes must survive (below)
            if len(v) != len(t):
                raise GuiError(f"{nm} is not a time series of the run")
            vals[nm] = np.degrees(v) if _series_meta(nm)["unit"] == "deg" else v
    # at most SERIES_LIMIT points that keep the smallest and largest value of every series in
    # each slice of the run (a plain stride steps over short peaks)
    idx = _decimate(list(vals.values()), SERIES_LIMIT) if vals else _decimate([t], SERIES_LIMIT)
    out = {"time_yr": [_num(x) for x in t[idx].tolist()], "n_total": len(t)}
    for nm, v in vals.items():
        out[nm] = [_num(x) for x in v[idx].tolist()]
    return out


def ensemble_report(root: Path, rel: str) -> dict:
    from ..ensemble import load_rows, report

    d = _inside(root, rel)
    src = d / "rows.jsonl"
    if not src.exists():
        return {"text": "no members finished yet (rows.jsonl is missing)"}
    # a large ensemble takes a while to read: the report is kept next to its rows
    cache = d / ".report.txt"
    stamp = f"{src.stat().st_mtime_ns} {src.stat().st_size}"
    if cache.exists():
        head, _, text = cache.read_text(encoding="utf-8").partition("\n")
        if head == stamp:
            return {"text": text, "cached": True}
    rows = load_rows(d)
    if not rows:
        return {"text": "no members finished yet"}
    text = report(rows)
    try:
        cache.write_text(stamp + "\n" + text, encoding="utf-8")
    except OSError:
        pass
    return {"text": text, "rows": len(rows)}


# --------------------------------------------------------------------------- results
def _code_stamp() -> str:
    """Changes with the code that makes the cards and pages (digest.py, server.py): a SWINGERsim
    updated in place rebuilds its cached cards even without a new version number."""
    h = hashlib.sha1()
    for name in ("digest.py", "server.py"):
        try:
            h.update((Path(__file__).parent / name).read_bytes())
        except OSError:
            pass
    return h.hexdigest()[:12]


CARD_VERSION = _code_stamp()  # the cached cards are rebuilt when it changes
FULL_CACHED = 6  # full digests (result pages) kept in memory


def _job_entry(root: Path, jid: str) -> tuple[str, Path, dict, dict]:
    jd = _jobdir(root)
    where, p = _find_job(jd, jid)
    if p is None:
        raise GuiError(f"no job {jid}")
    return where, p, _read_job(p), _status(jd / where, jid)


def _group_name(root: Path, job: dict) -> str:
    """Name a result is grouped under: its configuration, or for a job on another result (an
    emulator of an ensemble, an SBI) the folder or file it works on, or the command."""
    from . import digest as D

    if job.get("config"):
        return D.config_name(root, job["config"])
    for k in ("ensemble", "directory", "model", "h5"):
        v = job.get(k) or (job.get("args") or {}).get(k)
        if isinstance(v, str) and v:
            return Path(v.rstrip("/")).stem
    return str(job.get("command") or "")


def _result_kind(job: dict) -> str:
    from .digest import ROUTE

    return ROUTE.get(str(job.get("command")), "other")


def _file_stamp(q: Path | None) -> str:
    try:
        st = q.stat()
    except (OSError, AttributeError):
        return "-"
    return f"{st.st_mtime_ns}:{st.st_size}"


def _out_stamp(root: Path, job: dict) -> str:
    """Changes when the output of a job changes: the file, or the newest file of a folder (a
    file rewritten inside a folder does not change the folder's own time)."""
    out = job.get("output")
    if not out:
        return ""
    try:
        q = _inside(root, str(out))
    except GuiError:
        return ""
    if not q.exists():
        return "missing"
    if q.is_dir():
        best = (0, 0)
        n = 0
        try:
            for c in q.iterdir():
                if c.is_file():
                    st = c.stat()
                    best = max(best, (st.st_mtime_ns, st.st_size))
                    n += 1
        except OSError:
            pass
        return f"{best[0]}:{best[1]}:{n}"
    return _file_stamp(q)


def _inputs_stamp(root: Path, job: dict) -> str:
    """Changes when a file a digest reads besides the output changes: the configuration, the
    measured minima of its system (timings) and the seasons table (minima)."""
    from . import digest as D

    parts = []
    cfg = job.get("config")
    if isinstance(cfg, str) and cfg:
        try:
            parts.append(_file_stamp(_inside(root, cfg)))
        except GuiError:
            pass
        if job.get("command") == "timings":
            parts.append(_file_stamp(D.measured_minima_for(root, cfg)))
    seasons = (job.get("args") or {}).get("seasons")
    if isinstance(seasons, str) and seasons:
        try:
            parts.append(_file_stamp(_inside(root, seasons)))
        except GuiError:
            pass
    return ",".join(parts)


def _stamp(root: Path, where: str, p: Path, job: dict, st: dict) -> str:
    from .. import __version__

    return "|".join(
        str(x)
        for x in (
            CARD_VERSION,
            __version__,
            where,
            p.stat().st_mtime_ns,
            st.get("finished") or st.get("state"),
            _out_stamp(root, job),
            _inputs_stamp(root, job),
        )
    )


TEXT_BYTES = 4 * TEXT_LIMIT  # bytes of a log read for a page that shows its end


def _digest(root: Path, jid: str, full: bool, use: list[str] | None = None) -> dict:
    """The digest of a job's result (see gui/digest.py), with the job's own facts."""
    from . import digest as D

    where, _, job, st = _job_entry(root, jid)
    kind = _result_kind(job)
    out = job.get("output")
    cfg_rel = job.get("config")
    log = _jobdir(root) / where / f"{jid}.log"

    def log_text(whole: bool = False) -> str:
        """The log: whole (it is parsed) or its end (it is shown)."""
        if not log.exists():
            return ""
        if whole:
            return log.read_text(encoding="utf-8", errors="replace")
        with log.open("rb") as f:
            f.seek(0, os.SEEK_END)
            f.seek(max(0, f.tell() - TEXT_BYTES))
            return f.read().decode("utf-8", errors="replace")[-TEXT_LIMIT:]

    base = {
        "id": jid,
        "folder": where,
        "state": st.get("state") or ("cancelled" if where == "cancel" else where),
        "command": job.get("command"),
        "kind": kind,
        "config": cfg_rel,
        "config_name": _group_name(root, job),
        "output": out,
        "note": job.get("note") or "",
        "date": st.get("finished") or st.get("started") or "",
        "elapsed_s": st.get("elapsed_s"),
        "error": st.get("error") or "",
        "facts": [],
        "quality": [],
        "setup": "",
        "spark": None,
    }
    if where == "failed" or where == "cancel":
        base["facts"] = [
            {"k": "error", "v": (st.get("error") or base["state"])[:160], "note": "", "tone": "bad"}
        ]
        base["setup"] = str(job.get("command") or "")
        return base
    q = None
    if out:
        try:
            q = _inside(root, str(out))
        except GuiError:
            q = None
    try:
        if kind == "run" and q is not None and q.exists():
            base.update(D.run_digest(q, full=full))
        elif kind == "ensemble" and q is not None and q.exists():
            base.update(D.ensemble_digest(q, use=use, full=full))
        elif kind == "lightcurve" and q is not None and q.exists():
            base.update(D.lightcurve_digest(q, D._load_config(root, cfg_rel), job, full=full))
        elif kind == "timings" and q is not None and q.exists():
            tab = D._read_csv(q)  # read once for the digest and the measured minima
            base.update(D.timings_digest(q, D._load_config(root, cfg_rel), full=full, table=tab))
            mf = D.measured_minima_for(root, cfg_rel)
            cmp_ = D.compare_minima(q, mf, root, table=tab) if mf is not None else None
            if cmp_ and cmp_["binaries"]:
                for b, ser in cmp_["binaries"].items():
                    n = sum(len(v["res"]) for v in ser.values())
                    rr = [x for v in ser.values() for x in v["res"] if x is not None]
                    rms = math.sqrt(sum(x * x for x in rr) / len(rr)) if rr else None
                    mp, ms_ = (
                        (ser.get("primary") or {}).get("mean"),
                        (ser.get("secondary") or {}).get("mean"),
                    )
                    base["facts"].append(
                        {
                            "k": f"{b}: measured − model",
                            "v": f"rms {D.fmt(rms, 3)} min",
                            "note": f"{n} measured minima"
                            + (
                                f" · mean p {mp:+.1f}, s {ms_:+.1f} min"
                                if mp is not None and ms_ is not None
                                else ""
                            ),
                            "tone": "",
                        }
                    )
                if full and isinstance(base.get("data"), dict):
                    base["data"]["measured"] = cmp_
        elif kind == "fit" and q is not None and q.exists() and q.is_dir():
            cfg_path = None
            if cfg_rel:
                try:
                    cfg_path = _inside(root, str(cfg_rel))
                except GuiError:
                    cfg_path = None
            base.update(D.fit_digest(q, full=full, root=root, cfg_path=cfg_path))
        elif kind == "minima" and q is not None and q.exists():
            seasons = (job.get("args") or {}).get("seasons")
            sq = None
            if seasons:
                try:
                    sq = _inside(root, str(seasons))
                except GuiError:
                    sq = None
            base.update(D.minima_digest(q, log_text(whole=True), sq, full=full))
        elif kind == "benchmarks":
            base.update(D.benchmarks_digest(log_text(whole=True), full=full))
        else:
            base["setup"] = " ".join(
                x for x in (str(job.get("command") or ""), str(cfg_rel or "")) if x
            )
            text = log_text() if full else (_tail(log, 5) if log.exists() else "")
            lines = [x for x in text.splitlines() if x.strip()]
            other = (
                D.other_digest(str(job.get("command") or ""), log_text(whole=True))
                if log.exists()
                else {}
            )
            base["facts"] = other.get("facts") or [
                {
                    "k": "last line of the log",
                    "v": (lines[-1] if lines else "–")[:120],
                    "note": "",
                    "tone": "",
                }
            ]
            if full:
                base["data"] = {"text": text}
                # figures and the HTML report the job wrote into its output folder
                if q is not None and q.is_dir():
                    base["data"]["figures"] = [_rel(root, f) for f in sorted(q.glob("*.png"))]
                    if (q / "report.html").is_file():
                        base["data"]["html"] = _rel(root, q / "report.html")
    except Exception as exc:  # noqa: BLE001 - the card says what went wrong
        traceback.print_exc(file=sys.stderr)
        base["facts"] = [
            {"k": "could not read the output", "v": str(exc)[:200], "note": "", "tone": "bad"}
        ]
        if full:
            base["data"] = {"text": log_text()}
    if base.get("state") == "running":
        base["facts"].insert(
            0, {"k": "state", "v": "running", "note": "intermediate", "tone": "warn"}
        )
    return base


def _card_cache(root: Path) -> Path:
    d = root / "output" / ".gui_cache" / "cards"
    d.mkdir(parents=True, exist_ok=True)
    return d


def card(root: Path, jid: str, cached_only: bool = False) -> dict | None:
    """The card of a job (cached in output/.gui_cache/cards/ until the job, its output, its
    configuration or swingersim changes); None with cached_only when it still has to be computed."""
    safe_name(jid, "job id")
    where, p, job, st = _job_entry(root, jid)
    stamp = _stamp(root, where, p, job, st)
    f = _card_cache(root) / f"{jid}.json"
    if f.exists():
        try:
            c = json.loads(f.read_text(encoding="utf-8"))
            if c.get("_stamp") == stamp:
                return c
        except ValueError:
            pass
    if cached_only:
        return None
    c = _digest(root, jid, full=False)
    c["_stamp"] = stamp
    if where != "running":
        try:
            f.write_text(json.dumps(c, default=_json_default), encoding="utf-8")
        except OSError:
            pass
    return c


def results(root: Path) -> dict:
    """Every finished, failed and running job as a result card (cards still to be computed
    are marked pending), grouped by the name of their configuration."""
    jd = _jobdir(root)
    cards = []
    for sub in ("running", "done", "failed", "cancel"):
        for q in (jd / sub).glob("*.y*ml"):
            jid = _job_id(q.name)
            st = _status(jd / sub, jid)
            if sub == "cancel" and st.get("state") == "removed":
                continue  # removed from the inbox before it ran: no result
            c = None
            try:
                c = card(root, jid, cached_only=True)
            except (GuiError, OSError):
                continue
            if c is None:
                job = _read_job(q)
                c = {
                    "id": jid,
                    "folder": sub,
                    "state": st.get("state") or sub,
                    "command": job.get("command"),
                    "kind": _result_kind(job),
                    "config": job.get("config"),
                    "config_name": _group_name(root, job),
                    "output": job.get("output"),
                    "date": st.get("finished") or st.get("started") or "",
                    "pending": True,
                    "facts": [],
                    "setup": "",
                }
            cards.append(c)
    # the running jobs first, then the others newest first
    cards.sort(key=lambda c: (c.get("folder") == "running", c.get("date") or ""), reverse=True)
    groups: dict[str, int] = {}
    for c in cards:
        groups[c["config_name"]] = groups.get(c["config_name"], 0) + 1
    return {"cards": cards, "groups": [{"name": k, "n": v} for k, v in sorted(groups.items())]}


_FULL: dict[tuple, dict] = {}  # (job, stamp, use) -> full digest, oldest first
_FULL_LOCK = threading.Lock()


def full_digest(root: Path, jid: str, use: list[str] | None = None) -> dict:
    """The full digest of a result (the page), cached in memory under the card stamp (and the
    observations used to weight an ensemble) - a page reopened or re-weighted back is not read
    again. Running jobs are always read anew."""
    safe_name(jid, "job id")
    where, p, job, st = _job_entry(root, jid)
    key = (str(root), jid, _stamp(root, where, p, job, st), None if use is None else tuple(use))
    if where != "running":
        with _FULL_LOCK:
            hit = _FULL.get(key)
        if hit is not None:
            return hit
    d = _digest(root, jid, full=True, use=use)
    if where != "running":
        with _FULL_LOCK:
            for k in [k for k in _FULL if k[:2] == key[:2] and k[2] != key[2]]:
                del _FULL[k]  # an older state of the same result
            while len(_FULL) >= FULL_CACHED:
                del _FULL[next(iter(_FULL))]
            _FULL[key] = d
    return d


def result(root: Path, jid: str, use: list[str] | None = None) -> dict:
    """The full digest of one result, its files, and the other results of its configuration."""
    d = dict(full_digest(root, jid, use=use))  # the page's own keys go into a copy
    prod = job_products(root, jid)
    d["products"] = prod["products"]
    d["log"] = prod["log"]
    d["job"] = prod["job"]
    d["status"] = prod["status"]
    rel = []
    for c in results(root)["cards"]:
        if c["config_name"] == d["config_name"] and c["id"] != jid:
            rel.append(
                {k: c.get(k) for k in ("id", "kind", "command", "setup", "date", "state", "folder")}
            )
    d["related"] = rel[:20]
    return d


def save_member_draft(root: Path, jid: str, member: int) -> dict:
    """Draft of one ensemble member (jobs/drafts/<configuration>_member<N>.yaml)."""
    from . import digest as D

    _, _, job, _ = _job_entry(root, jid)
    folder = _inside(root, str(job.get("output") or ""))
    data = D.member_draft(root, folder, int(member))
    name = f"{D.config_name(root, job.get('config'))}_member{int(member)}"
    try:  # the paths of the configuration are relative to the folder it was read from
        meta = yaml.safe_load((folder / "ensemble.yaml").read_text(encoding="utf-8")) or {}
        base = meta.get("base_dir") if isinstance(meta, dict) else None
    except (OSError, yaml.YAMLError):
        base = None
    if isinstance(base, str) and base:
        rebase_paths(data, Path(base), _draft_path(root, name).parent)
    text = (
        f"# member {int(member)} of the ensemble {job.get('output')} (job {jid})\n"
        + yaml.safe_dump(data, sort_keys=False, allow_unicode=True)
    )
    return save_draft(root, name, text, overwrite=True)


# --------------------------------------------------------------------------- HTTP
class Handler(BaseHTTPRequestHandler):
    root: Path = Path(".")
    server_version = "swingersim-gui"

    def log_message(self, fmt, *args):  # quiet
        pass

    @staticmethod
    def _run_only(root: Path, body: dict, ids: list[str], res: dict) -> dict:
        """With "exclusive": only these jobs stay queued; with "start_worker": start it."""
        if body.get("exclusive"):
            res["held"] = hold_others(root, ids)
        if body.get("start_worker"):
            res["worker"] = start_worker(root)
        return res

    # -- responses
    def _send(self, code: int, body: bytes, ctype: str, extra: dict | None = None):
        extra = dict(extra or {})
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", extra.pop("Cache-Control", "no-store"))
        for k, v in extra.items():
            self.send_header(k, v)
        self.end_headers()
        self.wfile.write(body)

    def _json(self, data, code: int = 200):
        body = json.dumps(data, default=_json_default, allow_nan=False).encode()
        self._send(code, body, "application/json; charset=utf-8")

    def _error(self, exc: Exception, code: int = 400):
        if code >= 500:  # an error of swingersim itself: the traceback goes to the terminal
            print(f"swingersim gui: {self.command} {self.path}", file=sys.stderr)
            traceback.print_exception(exc, file=sys.stderr)
        self._json({"error": str(exc)}, code)

    def _allowed_origin(self) -> bool:
        """Requests only from the GUI page itself: the Host must be this server under a local
        name (DNS rebinding: a foreign site resolving its name to 127.0.0.1 sends its own Host),
        and an Origin header, when the browser sends one, must be the same host. Bound to the
        network (--host 0.0.0.0) the Host is not checked, the Origin still is."""
        host = (self.headers.get("Host") or "").strip().lower()
        port = self.server.server_address[1]
        local = {f"127.0.0.1:{port}", f"localhost:{port}", f"[::1]:{port}"}
        bound = self.server.server_address[0]
        if bound not in ("0.0.0.0", "", "::") and host not in local:
            return False
        origin = self.headers.get("Origin")
        return origin is None or origin.strip().lower() == f"http://{host}"

    def _body(self) -> dict:
        ctype = (self.headers.get("Content-Type") or "").split(";")[0].strip().lower()
        if ctype != "application/json":  # a cross-site form or no-cors fetch cannot send this
            raise GuiError("POST requests need Content-Type: application/json")
        n = int(self.headers.get("Content-Length") or 0)
        if n > 20_000_000:
            raise GuiError("request too large")
        raw = self.rfile.read(n) if n else b"{}"
        return json.loads(raw or b"{}")

    # -- routes
    def do_GET(self):
        if not self._allowed_origin():
            return self._json({"error": "request from a foreign page refused"}, 403)
        u = urlparse(self.path)
        q = {k: v[0] for k, v in parse_qs(u.query).items()}
        root = self.root
        try:
            if u.path in ("/", "/index.html"):
                return self._send(200, APP.read_bytes(), "text/html; charset=utf-8")
            if u.path == "/plotly.js":
                p = plotly_js_path()
                if p is None:
                    return self._send(HTTPStatus.FOUND, b"", "text/plain", {"Location": PLOTLY_CDN})
                return self._send(
                    200,
                    p.read_bytes(),
                    "application/javascript",
                    {"Cache-Control": "max-age=86400"},
                )
            if u.path == "/api/schema":
                return self._json(config_schema())
            if u.path == "/api/help":
                from .assist import help_texts

                return self._json(help_texts())
            if u.path == "/api/mesa":
                from .assist import mesa_status

                return self._json(mesa_status(root))
            if u.path == "/api/mesa_grid":
                from .assist import mesa_grid

                return self._json(mesa_grid(root, q.get("name", "")))
            if u.path == "/api/escape":
                from .assist import escape_compare

                return self._json(escape_compare(root, _inside(root, q.get("path", ""))))
            if u.path == "/api/methods":
                from .assist import methods_text

                return self._json(methods_text(_inside(root, q.get("path", ""))))
            if u.path == "/api/presets":
                from ..config import PRESETS

                return self._json(PRESETS)
            if u.path == "/api/configs":
                return self._json(list_configs(root))
            if u.path == "/api/config":
                return self._json(load_config(root, q.get("path", "")))
            if u.path == "/api/outputs":
                return self._json(scan_outputs(root, q.get("base", "output")))
            if u.path == "/api/output":
                return self._json(output_info(root, q.get("path", "")))
            if u.path == "/api/columns":
                cols = [c for c in q.get("cols", "").split("\x1f") if c]
                return self._json(columns(root, q.get("path", ""), cols))
            if u.path == "/api/viewer":
                return self._json({"page": run_viewer(root, q.get("path", ""))})
            if u.path == "/api/jobs":
                return self._json(jobs(root))
            if u.path == "/api/drafts":
                return self._json(list_drafts(root))
            if u.path == "/api/commands":
                from ..cli import command_catalog

                cat = command_catalog()
                for k in ("gui", "worker"):
                    cat.pop(k, None)
                for c in cat.values():  # jobs from the GUI use all cores unless told otherwise
                    for a in c["arguments"]:
                        if a["dest"] == "processes":
                            a["default"] = 0
                            a["help"] = "parallel processes (0 = all cores; the GUI default)"
                return self._json({"commands": cat, "outputs": OUTPUT_NAMES})
            if u.path == "/api/results":
                return self._json(results(root))
            if u.path == "/api/card":
                return self._json(card(root, q.get("id", "")))
            if u.path == "/api/result":  # the full digest: cached (full_digest)
                use = q.get("use")
                return self._json(
                    result(
                        root,
                        q.get("id", ""),
                        None if use is None else [x for x in use.split("\x1f") if x],
                    )
                )
            if u.path == "/api/finished":
                return self._json(finished_jobs(root))
            if u.path == "/api/job":
                return self._json(job_products(root, q.get("id", "")))
            if u.path == "/api/series_index":
                return self._json(h5_index(root, q.get("path", "")))
            if u.path == "/api/series":
                names = [c for c in q.get("names", "").split("\x1f") if c]
                return self._json(h5_series(root, q.get("path", ""), names))
            if u.path == "/api/ensemble_report":
                return self._json(ensemble_report(root, q.get("path", "")))
            if u.path == "/api/log":
                return self._json({"text": log_tail(root, q.get("path", ""))})
            if u.path in ("/api/system", "/api/inc_curve", "/api/oc_curve"):
                from . import modelviews as MV

                p = _inside(root, q.get("path", ""))
                if not p.is_file():
                    raise GuiError("no such configuration: " + q.get("path", ""))
                if u.path == "/api/system":
                    return self._json(MV.system(root, p))
                try:
                    t0, y0, y1 = float(q["t0"]), float(q["y0"]), float(q["y1"])
                except (KeyError, ValueError) as exc:
                    raise GuiError("t0, y0 and y1 must be numbers") from exc
                if u.path == "/api/inc_curve":
                    return self._json(MV.inc_curve(root, p, q.get("pair", "inner"), t0, y0, y1))
                return self._json(MV.oc_curve(root, p, t0, y0, y1))
            if u.path == "/api/version":
                from .. import __version__

                return self._json({"version": __version__, "root": str(root)})
            if u.path == "/file" or u.path.startswith("/files/"):
                # /files/<path>: relative links of a page (its PNG figures) resolve next to it
                rel = unquote(u.path[len("/files/") :]) if u.path != "/file" else q.get("path", "")
                p = _inside(root, rel)
                if not p.is_file():
                    raise GuiError("not a file")
                ctype = mimetypes.guess_type(p.name)[0] or "application/octet-stream"
                if ctype.startswith("text/"):
                    ctype += "; charset=utf-8"
                return self._send(200, p.read_bytes(), ctype)
            return self._send(404, b"not found", "text/plain")
        except GuiError as exc:
            return self._error(exc)
        except (BrokenPipeError, ConnectionResetError):  # the browser went away (page left)
            return None
        except Exception as exc:  # noqa: BLE001 - shown to the user in the page
            return self._error(exc, 500)

    def do_POST(self):
        if not self._allowed_origin():
            return self._json({"error": "request from a foreign page refused"}, 403)
        u = urlparse(self.path)
        root = self.root
        try:
            body = self._body()
            if u.path == "/api/validate":
                data = body.get("data")
                if data is None:
                    data = parse_yaml(body.get("yaml", ""))
                res = validate_config(data)
                res["data"] = data
                return self._json(res)
            if u.path == "/api/quickcheck":
                from .assist import quickcheck

                # relative paths of the configuration (a stellar table) are relative to the
                # file it was opened from; a new configuration becomes a draft
                base = root / "jobs" / "drafts"
                if body.get("source"):
                    try:
                        base = _inside(root, str(body["source"])).parent
                    except GuiError:
                        pass
                return self._json(quickcheck(body.get("data") or {}, base))
            if u.path == "/api/applegate":
                from .assist import applegate

                return self._json(
                    applegate(
                        body.get("data"),
                        body.get("star"),
                        body.get("p_mod_yr"),
                        body.get("oc_s"),
                        body.get("values"),
                    )
                )
            if u.path == "/api/mesa/save":
                from .assist import mesa_save

                return self._json(
                    mesa_save(root, body.get("spec") or {}, bool(body.get("overwrite")))
                )
            if u.path == "/api/save":  # not used by the page (tests, scripts); drafts below
                return self._json(
                    save_config(
                        root,
                        body.get("path", ""),
                        body.get("yaml", ""),
                        bool(body.get("overwrite")),
                    )
                )
            if u.path == "/api/draft/save":
                return self._json(
                    save_draft(
                        root,
                        body.get("name", ""),
                        body.get("yaml", ""),
                        bool(body.get("overwrite")),
                        body.get("source") or None,
                    )
                )
            if u.path == "/api/draft/delete":
                return self._json(delete_draft(root, body.get("name", "")))
            if u.path == "/api/submit":
                res = submit(
                    root,
                    list(body.get("drafts") or []),
                    body.get("command") or None,
                    body.get("args") or None,
                    body.get("output") or None,
                    body.get("note") or "",
                    bool(body.get("queue")),
                    str(body.get("suffix") or ""),
                )
                return self._json(self._run_only(root, body, res["jobs"], res))
            if u.path == "/api/submit_job":
                res = submit_job(
                    root,
                    body.get("job") or {},
                    bool(body.get("queue")),
                    body.get("name") or "job",
                )
                return self._json(self._run_only(root, body, res["jobs"], res))
            if u.path == "/api/run_now":
                return self._json(run_now(root, body.get("job") or {}))
            if u.path == "/api/job/action":
                ids = list(body.get("ids") or [])
                res = job_action(root, ids, body.get("action", ""))
                return self._json(self._run_only(root, body, ids, res))
            if u.path == "/api/member_draft":
                return self._json(
                    save_member_draft(root, body.get("id", ""), body.get("member", 0))
                )
            if u.path == "/api/worker/start":
                return self._json(start_worker(root))
            if u.path == "/api/worker/stop":
                return self._json(stop_worker(root))
            return self._send(404, b"not found", "text/plain")
        except GuiError as exc:
            return self._error(exc)
        except (BrokenPipeError, ConnectionResetError):  # the browser went away (page left)
            return None
        except Exception as exc:  # noqa: BLE001
            return self._error(exc, 500)


def _json_default(o):
    try:
        import numpy as np

        if isinstance(o, np.generic):
            v = o.item()
            return _num(v)
        if isinstance(o, np.ndarray):
            return [_num(x) for x in o.tolist()]
    except ImportError:
        pass
    return str(o)


def make_server(root: Path, host: str = "127.0.0.1", port: int = 8765) -> ThreadingHTTPServer:
    handler = type("QuadsimHandler", (Handler,), {"root": Path(root).absolute()})
    srv = ThreadingHTTPServer((host, port), handler)
    # old entries of output/.gui_cache go (in the background: the cache can be large)
    threading.Thread(target=prune_cache, args=(Path(root).absolute(),), daemon=True).start()
    return srv


def serve(root: Path, host: str = "127.0.0.1", port: int = 8765, open_browser: bool = True):
    srv = make_server(root, host, port)
    url = f"http://{'127.0.0.1' if host in ('0.0.0.0', '') else host}:{srv.server_address[1]}/"
    print(f"swingersim gui on {url}  (folder {Path(root).absolute()}; Ctrl+C stops it)")
    if host in ("0.0.0.0", ""):
        print("  reachable from the network: anyone who can reach this port sees the folder")
    if open_browser:
        print("  if no browser window opens, open the address above yourself (the server runs)")
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        srv.server_close()


__all__ = [
    "classify",
    "columns",
    "finished_jobs",
    "h5_index",
    "h5_series",
    "job_action",
    "job_products",
    "jobs",
    "list_drafts",
    "make_server",
    "output_info",
    "run_now",
    "save_draft",
    "scan_outputs",
    "serve",
    "start_worker",
    "stop_worker",
    "submit",
    "validate_config",
]
