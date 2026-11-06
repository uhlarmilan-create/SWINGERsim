"""Model curves for the Outputs page of swingersim gui, computed from a configuration on request and
cached in output/.gui_cache/models/ (keyed by the file, its modification time, the files it
refers to, the SWINGERsim version and the request):

- system:    what the configuration looks like (systemview.system_model): stars, orbits, eclipse
             geometry, stability of the nested triples;
- inc_curve: the inclination of a pair against time (and its eclipse limit) over a span of years;
- oc_curve:  every eclipse of every pair of two stars over a span of years (BJD of the minima),
             for an O - C diagram with the model as a dense curve.

The pages draw these with plotly (interactive), not as images.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from ..units import DAYS_PER_YEAR as DPY

BJD_J2000 = 2451544.5  # BJD of the calendar year 2000.0 used by the pages (year = 2000 + ...)


def _referenced(path: Path) -> list[str]:
    """Stamps (path, mtime, size) of the files a configuration refers to: a changed stellar
    table or minima table changes the model curves too."""
    import yaml

    from .server import PATH_FIELDS

    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError):
        return []
    found: list[str] = []

    def visit(node, keys):
        k = keys[0]
        if k == "*":
            for item in node if isinstance(node, list) else []:
                visit(item, keys[1:])
            return
        if not isinstance(node, dict) or k not in node:
            return
        if len(keys) > 1:
            visit(node[k], keys[1:])
        elif isinstance(node[k], str):
            v = node[k].rpartition(":")[0] if k == "entry" else node[k]
            if v:
                found.append(v)

    for keys in PATH_FIELDS:
        if keys[-1] != "dir":  # an output folder is not an input
            visit(data, keys)
    out = []
    for v in found:
        q = Path(v) if Path(v).is_absolute() else path.parent / v
        try:
            st = q.stat()
            out.append(f"{q}:{st.st_mtime_ns}:{st.st_size}")
        except OSError:
            out.append(f"{q}:-")
    return out


def _cache(root: Path, path: Path, what: str, params: dict) -> Path:
    from .. import __version__

    st = path.stat()
    key = json.dumps(
        [str(path), st.st_mtime_ns, st.st_size, _referenced(path), __version__, what, params],
        sort_keys=True,
    )
    d = root / "output" / ".gui_cache" / "models"
    d.mkdir(parents=True, exist_ok=True)
    return d / (hashlib.sha1(key.encode()).hexdigest()[:20] + ".json")


def _cfg(path: Path):
    from ..config import load_config

    return load_config(path)


def _cached(root: Path, path: Path, what: str, params: dict, compute) -> dict:
    c = _cache(root, path, what, params)
    if c.exists():
        try:
            out = json.loads(c.read_text())
            c.touch()  # in use: kept by the pruning of the cache
            return out
        except (ValueError, OSError):
            pass
    out = compute()
    c.write_text(json.dumps(out))
    return out


def system(root: Path, path: Path) -> dict:
    from ..systemview import system_model

    return _cached(root, path, "system", {}, lambda: system_model(_cfg(path)))


def _span(y0: float, y1: float) -> None:
    if not (np.isfinite(y0) and np.isfinite(y1) and y0 < y1):
        raise ValueError(f"the span of years must have y0 < y1 (got {y0} and {y1})")


def inc_curve(
    root: Path, path: Path, pair: str, bjd_t0: float, y0: float, y1: float, n: int = 400
) -> dict:
    """Inclination of the pair [deg, 0-90] and its eclipse limit at n times between the calendar
    years y0 and y1 (BJD of t = 0 of the configuration: bjd_t0)."""
    _span(y0, y1)

    def compute():
        from ..incfit import model_inclination

        years = np.linspace(y0, y1, n)
        bjd = BJD_J2000 + (years - 2000.0) * DPY
        m = model_inclination(_cfg(path), pair, bjd, bjd_t0, with_limit=True)
        return {
            "pair": pair,
            "year": np.round(years, 4).tolist(),
            "i": [None if not np.isfinite(v) else round(float(v), 4) for v in m[0]],
            "limit": [None if not np.isfinite(v) else round(float(v), 4) for v in m[1]],
        }

    return _cached(
        root, path, "inc", {"pair": pair, "t0": bjd_t0, "y0": y0, "y1": y1, "n": n}, compute
    )


def oc_curve(
    root: Path, path: Path, bjd_t0: float, y0: float, y1: float, max_points: int = 6000
) -> dict:
    """BJD of every model eclipse between the calendar years y0 and y1, per pair and kind."""
    _span(y0, y1)

    def compute():
        from ..timing import eclipse_times

        cfg = _cfg(path)
        # the span in years from t = 0 of the configuration: [a, b]
        a = (BJD_J2000 + (y0 - 2000.0) * DPY - bjd_t0) / DPY
        b = (BJD_J2000 + (y1 - 2000.0) * DPY - bjd_t0) / DPY
        ecl = []
        if a < 0:  # backwards from t = 0 to a (or only to b when the span ends before t = 0)
            ecl += [e for e in eclipse_times(cfg, a)[0] if a <= e.t_obs_yr <= min(b, 0.0)]
        if b > 0:  # forwards from t = 0 to b (or from a when the span starts after t = 0)
            ecl += [e for e in eclipse_times(cfg, b)[0] if max(a, 0.0) <= e.t_obs_yr <= b]
        out: dict = {}
        for e in ecl:
            d = out.setdefault(e.binary, {}).setdefault(e.kind, [])
            d.append(bjd_t0 + e.t_obs_yr * DPY)
        for kinds in out.values():
            for k, v in kinds.items():
                v = np.unique(np.array(v))  # sorted; an eclipse at t = 0 only once
                if len(v) > max_points:
                    v = v[np.linspace(0, len(v) - 1, max_points).astype(int)]
                kinds[k] = np.round(v, 6).tolist()
        return {"bjd_t0": bjd_t0, "eclipses": out}

    return _cached(root, path, "oc", {"t0": bjd_t0, "y0": y0, "y1": y1}, compute)
