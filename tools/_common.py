"""Helpers shared by the report scripts in tools/ (squadra_report, squadra_figures, squadra_fits,
squadra3_report, squadra3_posterior, make_lc_windows, tess_minima, build_docs_pdf).

* calendar years <-> BJD, the matplotlib style and colours of the figures;
* figure output: a stale PNG is deleted before a new one is written, failures are collected and
  are fatal at the end unless the PNGs were declined (--no-png);
* caches that invalidate: the key holds the SWINGERsim version, the request parameters and the
  size and modification time of every input file (configurations, fitted configurations, data
  files), and every cache file is written atomically;
* the epoch (BJD of t = 0) of a configuration, ensemble weights from the ensemble configuration,
  and small numerical helpers used by several figures.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import tempfile
from collections.abc import Callable, Iterable
from pathlib import Path

import numpy as np

DPY = 365.25
BJD_2000 = 2451544.5  # BJD of 2000 January 1.0: the calendar-year axis of the figures


def year(bjd):
    """Calendar year of a BJD (2000.0 = BJD 2451544.5, years of 365.25 d)."""
    return 2000.0 + (np.asarray(bjd, float) - BJD_2000) / DPY


def bjd(yr):
    """BJD of a calendar year (inverse of year)."""
    return BJD_2000 + (np.asarray(yr, float) - 2000.0) * DPY


# ------------------------------------------------------------------ figures
C_DATA = "#3A3F45"

MPL_RC = {
    "font.size": 9,
    "axes.titlesize": 9.5,
    "axes.labelsize": 9,
    "legend.fontsize": 7.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": "#E4E4E0",
    "grid.linewidth": 0.6,
    "savefig.dpi": 160,
}


def mpl_pyplot():
    """matplotlib.pyplot with the Agg backend and the style of the report figures."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(MPL_RC)
    return plt


class FigureFailures(RuntimeError):
    pass


class FigureWriter:
    """Writes the figures of a report into a folder and remembers the PNGs that failed.

    A PNG is deleted before it is written, so a failed export cannot leave the PNG of an
    earlier run behind for the PDF. finish() raises FigureFailures when a PNG failed, unless
    png=False (the PNGs were declined, e.g. no kaleido on the machine: --no-png)."""

    def __init__(self, out: Path, png: bool = True):
        self.out = Path(out)
        self.out.mkdir(parents=True, exist_ok=True)
        self.png = png
        self.names: list[str] = []
        self.failed: list[tuple[str, str]] = []

    def _png(self, name: str, write: Callable[[Path], None]) -> None:
        p = self.out / f"{name}.png"
        p.unlink(missing_ok=True)
        if not self.png:
            return
        try:
            write(p)
        except Exception as exc:  # noqa: BLE001 - reported (and fatal) in finish()
            p.unlink(missing_ok=True)
            self.failed.append((name, f"{type(exc).__name__}: {exc}"))
            print(f"  {name}.png not written: {exc}", flush=True)

    def save_mpl(self, fig, name: str) -> str:
        import matplotlib.pyplot as plt

        try:
            self._png(name, lambda p: fig.savefig(p, bbox_inches="tight"))
        finally:
            plt.close(fig)
        self.names.append(name)
        return f"{name}.png"

    def save_plotly(self, fig, name: str, scale: float = 2.0) -> str:
        """The JSON (interactive HTML version) and the PNG (PDF) of a plotly figure."""
        atomic_write_text(self.out / f"{name}.json", fig.to_json())
        self._png(name, lambda p: fig.write_image(str(p), scale=scale))
        self.names.append(name)
        print("  " + name, flush=True)
        return f"{name}.png"

    def finish(self) -> None:
        if self.failed and self.png:
            lines = "\n".join(f"  {n}: {e}" for n, e in self.failed)
            raise FigureFailures(
                f"{len(self.failed)} PNG figure(s) could not be written (old ones were removed):\n"
                f"{lines}\nrun with --no-png to write the JSON figures only"
            )


# ------------------------------------------------------------------ caches
def atomic_write_text(path: Path, text: str) -> None:
    """Write a text file atomically (a temporary file in the same folder, then os.replace):
    an interrupted run never leaves a truncated cache or figure behind."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(text)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def file_stamp(path: Path) -> list:
    """(path, size, mtime) of an input file; a missing file is part of the key as well."""
    path = Path(path)
    try:
        st = path.stat()
    except FileNotFoundError:
        return [str(path), None, None]
    return [str(path.resolve()), st.st_size, st.st_mtime_ns]


def digest(obj) -> str:
    """Stable hash of a JSON-like object (configuration dicts, request parameters)."""
    text = json.dumps(obj, sort_keys=True, default=_jsonable)
    return hashlib.sha1(text.encode()).hexdigest()


def _jsonable(x):
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    if isinstance(x, Path):
        return str(x)
    return repr(x)


def cache_key(files: Iterable[Path] = (), params=None) -> str:
    """Key of a cached computation: the SWINGERsim version, the inputs (size and modification
    time of every file) and the request (parameters, e.g. t0 and the grid of years, or the
    configuration dict itself when it was modified in memory)."""
    import swingersim

    return digest(
        {
            "swingersim": swingersim.__version__,
            "files": [file_stamp(f) for f in files],
            "params": params,
        }
    )


class JsonCache:
    """JSON cache in a folder: cache.get(name, compute, files, params) returns the cached value
    when its key (cache_key) is unchanged, else computes and stores it (one file per name, so
    an outdated value is replaced, not kept beside the new one)."""

    def __init__(self, folder: Path):
        self.folder = Path(folder)

    def get(self, name: str, compute: Callable[[], object], files=(), params=None):
        key = cache_key(files, params)
        f = self.folder / f"{name}.json"
        if f.exists():
            try:
                d = json.loads(f.read_text(encoding="utf-8"))
                if isinstance(d, dict) and d.get("_key") == key:
                    return d["value"]
            except (ValueError, KeyError):
                pass
        value = compute()
        atomic_write_text(f, json.dumps({"_key": key, "value": value}, default=_jsonable))
        return value


def file_cached(path: Path, write: Callable[[Path], None], files=(), params=None) -> bool:
    """A cached output file (CSV of a model light curve or of eclipse times): write(tmp) makes
    it when the key (cache_key, kept in path + '.key') changed or the file is missing; the file
    is moved into place atomically. Returns True when it was (re)computed."""
    path = Path(path)
    key = cache_key(files, params)
    kf = path.with_name(path.name + ".key")
    if path.exists() and kf.exists() and kf.read_text(encoding="utf-8").strip() == key:
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp{path.suffix}")
    try:
        write(tmp)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)
    atomic_write_text(kf, key + "\n")
    return True


# ------------------------------------------------------------------ epochs and inputs
_EPOCH = re.compile(r"t = 0(?: of this configuration)? = BJD (\d+(?:\.\d+)?)")


def fit_epoch(path: Path, default: float | None = None) -> float:
    """BJD (TDB) of t = 0 of a configuration file: the key ``fit_epoch_bjd_tdb`` written by
    swingersim fit into stage2_fitted.yaml, else fit.stage2.bjd_t0, else the comment line
    '# t = 0 [of this configuration] = BJD ...'. Raises ValueError when none is found and no
    default is given."""
    import yaml

    text = Path(path).read_text(encoding="utf-8")
    d = yaml.safe_load(text) or {}
    if d.get("fit_epoch_bjd_tdb") is not None:
        return float(d["fit_epoch_bjd_tdb"])
    st2 = ((d.get("fit") or {}).get("stage2")) or {}
    if st2.get("bjd_t0") is not None:
        return float(st2["bjd_t0"])
    for line in text.splitlines():
        if line.lstrip().startswith("#"):
            m = _EPOCH.search(line)
            if m:
                return float(m.group(1))
    if default is not None:
        return float(default)
    raise ValueError(f"{path}: no epoch (fit_epoch_bjd_tdb, fit.stage2.bjd_t0 or '# t = 0 = BJD')")


def read_table(path: Path):
    """CSV with optional '#' comment lines (the files of SWINGERsim and of the tools)."""
    import pandas as pd

    return pd.read_csv(path, comment="#")


def ensemble_observed(config: Path, gaussian_only: bool = False) -> dict:
    """The ``ensemble.observed`` conditions of an ensemble configuration (column -> {value,
    sigma} and/or {min, max, equals}); gaussian_only: only the {value, sigma} entries."""
    import yaml

    d = yaml.safe_load(Path(config).read_text(encoding="utf-8")) or {}
    if "config" in d and "ensemble" not in d:  # ensemble.yaml of a run: the config under 'config'
        d = d["config"] or {}
    obs = (d.get("ensemble") or {}).get("observed") or {}
    out = {}
    for col, o in obs.items():
        o = {k: v for k, v in (o or {}).items() if v is not None}
        if gaussian_only:
            if o.get("value") is None or o.get("sigma") is None:
                continue
            o = {"value": o["value"], "sigma": o["sigma"]}
        out[col] = o
    return out


def ensemble_rows(summary_csv: Path) -> list[dict]:
    """Rows of an ensemble summary.csv as swingersim.ensemble expects them (NaN for empty)."""
    import pandas as pd

    d = pd.read_csv(summary_csv)
    rows = d.to_dict("records")
    for r in rows:
        r["member"] = int(r["member"])
    return rows


def ensemble_weights(rows: list[dict], observed: dict, members) -> np.ndarray:
    """Weights (swingersim.ensemble.observation_weights) of the given members, in their order."""
    from swingersim.ensemble import observation_weights

    w = observation_weights(rows, observed)
    return np.array([w.get(int(m), 0.0) for m in members], float)


# ------------------------------------------------------------------ numerical helpers
def nodal_period(years, h, L) -> float:
    """Period [yr] of the precession of the unit vector h (one per year) about L: the slope of
    its unwrapped azimuth around L."""
    h = np.asarray(h, float)
    L = np.asarray(L, float)
    e1 = np.cross(L, [1.0, 0.0, 0.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(L, e1)
    ang = np.unwrap(np.arctan2(h @ e2, h @ e1))
    return float(2.0 * np.pi / abs(np.polyfit(np.asarray(years, float), ang, 1)[0]))


def eclipse_switches(years, inc, limit) -> list[float]:
    """Calendar years at which the pair starts or stops eclipsing (i crosses its limit)."""
    years = np.asarray(years, float)
    ecl = np.asarray(inc, float) > np.asarray(limit, float)
    return [float(x) for x in years[1:][np.diff(ecl.astype(int)) != 0]]


def band_luminosity(radius_rsun: float, teff: float, band: str) -> float:
    """R^2 x passband surface brightness: the light of a star in a band (relative units)."""
    from swingersim.photometry import surface_brightness

    return float(radius_rsun) ** 2 * surface_brightness(float(teff), band)


def time_bins(t, f, per_day: float) -> tuple[np.ndarray, np.ndarray]:
    """Means of (t, f) in bins of 1/per_day days counted from the first time."""
    import pandas as pd

    t = np.asarray(t, float)
    b = np.floor((t - t[0]) * per_day)
    g = pd.DataFrame({"b": b, "t": t, "f": np.asarray(f, float)}).groupby("b").mean()
    return g.t.values, g.f.values


def wrap_deg(x):
    """Angle differences [deg] wrapped to (-180, 180]."""
    return -((180.0 - np.asarray(x, float)) % 360.0 - 180.0)


def nested_value(d: dict, dotted: str):
    """d['a']['b']['c'] for 'a.b.c'; e and omega of a fit in ecosw / esinw are derived from
    the orbit's e and omega."""
    parts = dotted.split(".")
    if parts[-1] in ("ecosw", "esinw") and len(parts) == 3:
        o = d[parts[0]][parts[1]]
        w = math.radians(float(o.get("omega", 0.0)))
        e = float(o.get("e", 0.0))
        return e * (math.cos(w) if parts[-1] == "ecosw" else math.sin(w))
    x = d
    for p in parts:
        x = x[p]
    return float(x)
