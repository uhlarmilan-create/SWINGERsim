"""Result digests for swingersim gui: what the Outputs page shows of a finished job.

Every finished (or running) job becomes a *result*: a card with a few facts read from its
output, and a page made for the kind of output - a run, an ensemble, a light curve, eclipse
timings, minima, a validation. The digests only read the files a job wrote (and the
configuration it used); they never interpret the system ("a 2+(2+1) quadruple"), they state
what the output contains (ranges, fractions, events, stability ratios, effective numbers of
members).

A digest is a plain JSON-able dict:

    {"id", "folder", "state", "command", "kind", "config_name", "config", "output",
     "date", "setup", "facts": [{"k", "v", "tone"}], "spark": {...} | None,
     "data": {...}}                      # the page (only in the full digest)

``tone`` is "good", "warn", "bad" or "" (neutral). The cards are cached in
``output/.gui_cache/cards/`` and rebuilt when the job or its output changes.
"""

from __future__ import annotations

import csv
import json
import math
import re
import threading
from pathlib import Path

import numpy as np
import yaml

# thresholds of the data-quality marks (the same for every result)
ENERGY_GOOD, ENERGY_WARN = 1e-6, 1e-3
NEFF_GOOD, NEFF_WARN = 100.0, 20.0
MAX_POINTS = 1500  # points per time series sent to the browser
MAX_FRAMES = 2400  # outputs of a run sent for the system view
LC_POINTS = 30000  # points of a light curve

ROUTE = {  # command -> page
    "run": "run",
    "ensemble": "ensemble",
    "lightcurve": "lightcurve",
    "timings": "timings",
    "minima": "minima",
    "validate": "benchmarks",
    "fit": "fit",
}


# --------------------------------------------------------------------------- helpers
_SUP = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")


def _f(x) -> float | None:
    """JSON-safe float (None for NaN, infinities and missing values)."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _tone_energy(err: float | None, dissipative: bool = False) -> str:
    if err is None:
        return ""
    if dissipative:
        return ""
    return "good" if err < ENERGY_GOOD else "warn" if err < ENERGY_WARN else "bad"


def _tone_ratio(r: float | None) -> str:
    if r is None:
        return ""
    return "good" if r >= 1.5 else "warn" if r >= 1.0 else "bad"


def _tone_neff(n: float | None) -> str:
    if n is None:
        return ""
    return "good" if n >= NEFF_GOOD else "warn" if n >= NEFF_WARN else "bad"


def fmt(v, digits: int = 4) -> str:
    """Compact number for facts and tables."""
    if v is None:
        return "–"
    if isinstance(v, str):
        return v
    x = float(v)
    if not math.isfinite(x):
        return "–"
    a = abs(x)
    if a == 0:
        return "0"
    if a >= 1e5 or a < 1e-3:
        m, e = f"{x:.{max(digits - 2, 1)}e}".split("e")
        return f"{m}·10{str(int(e)).translate(_SUP)}"
    out = f"{x:.{digits}g}"
    if "e" in out:  # 100 .. 99999 with fewer digits: rounded integer, not 1.3e+03
        return f"{float(out):.0f}"
    return out


def fmt_range(lo, hi, digits: int = 4, unit: str = "") -> str:
    if lo is None or hi is None:
        return "–"
    return f"{fmt(lo, digits)} – {fmt(hi, digits)}{unit}"


def _down(n: int, limit: int) -> np.ndarray:
    """Indices of at most ``limit`` evenly spaced samples (first and last kept)."""
    if n <= limit:
        return np.arange(n)
    return np.unique(np.linspace(0, n - 1, limit).round().astype(int))


def _envelope(y, buckets: int) -> np.ndarray:
    """Indices that keep the first and the last sample and, in each of at most ``buckets``
    equal slices of the series, its smallest and its largest value (and the first gap, NaN, of a
    slice that has one). A series decimated like this keeps its extremes - plain striding can
    step over a narrow eclipse or a short eccentricity peak. At most 3 * buckets + 2 indices."""
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n <= 2 * buckets + 2:
        return np.arange(n)
    k = -(-n // buckets)  # samples per slice
    m = -(-n // k)  # slices
    yy = np.full(m * k, np.nan)
    yy[:n] = y
    yy = yy.reshape(m, k)
    nan = np.isnan(yy)
    lo = np.where(nan, np.inf, yy).argmin(axis=1)
    hi = np.where(nan, -np.inf, yy).argmax(axis=1)
    base = np.arange(m) * k
    parts = [np.array([0, n - 1]), base + lo, base + hi]
    gap = nan.any(axis=1) & ~nan.all(axis=1)
    gap[-1] = gap[-1] and bool(np.isnan(y[(m - 1) * k :]).any())  # not the padding
    if gap.any():
        parts.append((base + nan.argmax(axis=1))[gap])
    keep = np.zeros(n, dtype=bool)
    for p in parts:
        keep[p[p < n]] = True
    return np.flatnonzero(keep)


def _envelope_union(series, buckets: int) -> np.ndarray:
    """Indices that keep the extremes (see _envelope) of every series of the same length."""
    series = [np.asarray(s, dtype=float) for s in series if s is not None]
    if not series:
        return np.arange(0)
    keep = np.zeros(len(series[0]), dtype=bool)
    for s in series:
        keep[_envelope(s, buckets)] = True
    return np.flatnonzero(keep)


def _decimate(series, limit: int) -> np.ndarray:
    """At most ``limit`` indices of series of the same length that keep the extremes of every
    series (_envelope) with as many slices as the limit allows: flat stretches need one point
    per slice, so the slices are made as fine as the limit permits (bisection)."""
    series = [np.asarray(s, dtype=float) for s in series if s is not None]
    n = len(series[0]) if series else 0
    if n <= limit:
        return np.arange(n)
    lo, hi = 1, max(limit // 2, 1)
    best = _envelope_union(series, lo)
    while lo <= hi:  # the largest number of slices whose indices fit into the limit
        mid = (lo + hi) // 2
        idx = _envelope_union(series, mid)
        if len(idx) <= limit:
            best, lo = idx, mid + 1
        else:
            hi = mid - 1
    return best


def _arr(x, idx=None, scale: float = 1.0, digits: int = 7) -> list:
    a = np.asarray(x, dtype=float)
    if idx is not None:
        a = a[idx]
    if scale != 1.0:
        a = a * scale
    spec = f".{digits}g"
    isfinite = math.isfinite
    return [float(format(v, spec)) if isfinite(v) else None for v in a.tolist()]


def _range(x) -> list:
    """[min, max] of the finite values of a series (the whole series, not a decimated one)."""
    a = np.asarray(x, dtype=float)
    a = a[np.isfinite(a)]
    return [float(a.min()), float(a.max())] if len(a) else [None, None]


def _nanmax(x) -> float | None:
    a = np.asarray(x, dtype=float)
    a = a[np.isfinite(a)]
    return float(a.max()) if len(a) else None


def _nanmin(x) -> float | None:
    a = np.asarray(x, dtype=float)
    a = a[np.isfinite(a)]
    return float(a.min()) if len(a) else None


def _series(t, y, limit: int, digits: int = 7, scale: float = 1.0) -> dict:
    """A time series for a plot: decimated with its extremes kept, and the min - max range of
    the whole series (for the caption)."""
    y = np.asarray(y, dtype=float) * scale
    idx = _decimate([y], limit)
    lo, hi = _range(y)
    return {
        "t": _arr(t, idx, digits=10),
        "y": _arr(y, idx, digits=digits),
        "range": [_f(lo), _f(hi)],
        "n_total": len(y),
    }


def _read_csv(path: Path) -> dict[str, np.ndarray]:
    """Columns of a CSV table: numbers as float arrays (an empty cell is NaN), text as object
    arrays of str. Lines starting with '#' (headers of the writers) are skipped. pandas parses
    when it is installed (extra 'ensemble'), else the csv module."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    if text.startswith("#") or "\n#" in text or "\n #" in text:
        text = "".join(
            ln for ln in text.splitlines(keepends=True) if not ln.lstrip().startswith("#")
        )
    if not text.strip():
        return {}
    try:
        import pandas as pd
    except ImportError:
        pd = None
    if pd is None:
        rows = [r for r in csv.reader(text.splitlines()) if r]
        head = [h.strip() for h in rows[0]]
        body = rows[1:]
        return {
            h: _column(np.array([r[k] if k < len(r) else "" for r in body], dtype=object))
            for k, h in enumerate(head)
        }
    import io

    df = pd.read_csv(
        io.StringIO(text),
        keep_default_na=False,
        na_values=["", "nan", "NaN"],
        skipinitialspace=True,
        skip_blank_lines=True,
    )
    out = {}
    for c in df.columns:
        s = df[c]
        if s.dtype.kind in "iufb":
            out[str(c).strip()] = s.to_numpy(dtype=float)
        else:
            s = s.astype(object).where(s.notna(), "").astype(str).str.strip()
            out[str(c).strip()] = _column(s.to_numpy(dtype=object))
    return out


def _column(v: np.ndarray) -> np.ndarray:
    """Floats if every non-empty cell is a number, else the stripped text."""
    s = np.array([x.strip() for x in v], dtype=object)
    try:
        return np.where(s == "", "nan", s).astype(float)
    except ValueError:
        return s


def _numbers(d: dict, key: str) -> np.ndarray:
    """A column as floats (text cells become NaN)."""
    v = d[key]
    if v.dtype != object:
        return v
    out = np.full(len(v), np.nan)
    for i, x in enumerate(v):
        try:
            out[i] = float(x)
        except ValueError:
            pass
    return out


def config_name(root: Path, rel: str | None) -> str:
    """Name of the configuration a job used: the file stem, or, for the frozen copy made when
    a draft was sent to the inbox, the name of the draft it was copied from."""
    if not rel:
        return ""
    p = root / rel
    name = Path(rel).stem
    try:
        first = p.open(encoding="utf-8").readline()
    except OSError:
        return name
    m = re.match(r"#\s*frozen copy of\s+(\S+?)\.ya?ml", first)
    return m.group(1) if m else name


def _load_config(root: Path, rel: str | None):
    """The Config of a job (None if it cannot be read: e.g. a scattering experiment)."""
    if not rel:
        return None
    from ..config import load_config

    try:
        return load_config(root / rel)
    except Exception:  # noqa: BLE001 - a page without the configuration still shows the data
        return None


def _tree(cfg) -> dict:
    """Hierarchy of a configuration: orbits outermost first, with their two children, and the
    two-star binaries (the orbits that can eclipse)."""
    if cfg is None:
        return {"orbits": [], "binaries": [], "stars": []}
    tree = cfg.tree
    children = {c for kids in tree.values() for c in kids}
    root = next(n for n in tree if n not in children)
    order: list[str] = []

    def walk(n):
        order.append(n)
        for c in tree[n]:
            if c in tree:
                walk(c)

    walk(root)
    orbits = [{"name": n, "children": list(tree[n])} for n in order]
    binaries = [n for n in order if all(c not in tree for c in tree[n])]
    return {"orbits": orbits, "binaries": binaries, "stars": list(cfg.star_names), "root": root}


# --------------------------------------------------------------------------- run (HDF5)
def _attrs(obj) -> dict:
    return {k: (v.decode() if isinstance(v, bytes) else v) for k, v in obj.attrs.items()}


class _RunFile:
    """The datasets of a run file, read when they are needed (a card reads a few of them; the
    files may be gzip-compressed and the hierarchy group may lack its detected dataset)."""

    def __init__(self, f):
        self.f = f
        self._memo: dict[str, np.ndarray | None] = {}

    def get(self, name: str) -> np.ndarray | None:
        if name not in self._memo:
            obj = self.f.get(name)
            ok = obj is not None and hasattr(obj, "shape") and obj.ndim == 1
            self._memo[name] = np.asarray(obj[()]) if ok else None
        return self._memo[name]

    def ends(self, name: str) -> tuple[float | None, float | None]:
        """First and last value of a dataset without reading it."""
        if name in self._memo and self._memo[name] is not None:
            v = self._memo[name]
            return (_f(v[0]), _f(v[-1])) if len(v) else (None, None)
        obj = self.f.get(name)
        if obj is None or not hasattr(obj, "shape") or obj.ndim != 1 or not len(obj):
            return None, None
        return _f(obj[0]), _f(obj[-1])

    def keys(self, group: str) -> list[str]:
        g = self.f.get(group)
        return list(g.keys()) if g is not None and hasattr(g, "keys") else []

    def length(self, name: str) -> int:
        obj = self.f.get(name)
        return len(obj) if obj is not None and hasattr(obj, "shape") and obj.ndim == 1 else 0


def run_digest(path: Path, full: bool = True) -> dict:
    import h5py

    from ..config import Config
    from ..units import RSUN_AU

    out: dict = {}
    with h5py.File(path, "r") as f:
        R = _RunFile(f)
        attrs = _attrs(f)
        try:
            cfg = Config.model_validate(yaml.safe_load(attrs.get("config_yaml", "")))
        except Exception:  # noqa: BLE001
            cfg = None
        tree = _tree(cfg)
        n = R.length("time_yr") if "time_yr" in f else 1
        t0, t_end = R.ends("time_yr") if "time_yr" in f else (0.0, 0.0)
        orbit_names = R.keys("orbits")
        mutual_names = R.keys("mutual_inclination")
        star_attrs, star_keys = {}, {}
        for name in R.keys("stars"):
            g = f["stars"][name]
            star_attrs[name] = _attrs(g)
            star_keys[name] = [
                k for k in g if isinstance(g[k], h5py.Dataset) and g[k].ndim == 1 and len(g[k])
            ]
        binaries = json.loads(attrs.get("binaries_json", "{}") or "{}")
        system = json.loads(attrs.get("system_json", "{}") or "{}")
        transitions = []
        configured = ""
        if "hierarchy" in f:  # attributes only: the detected dataset may be missing
            ha = _attrs(f["hierarchy"])
            transitions = json.loads(ha.get("transitions_json", "[]") or "[]")
            # a dynamically started mass transfer is a continuation, not a change of hierarchy
            transitions = [t for t in transitions if t.get("kind") != "mass_transfer_onset"]
            configured = str(ha.get("configured", ""))
        diags = [_attrs(f["diagnostics"][inner]) for inner in R.keys("diagnostics")]
        energy = R.get("energy_error")
        angmom = R.get("angmom_error")
        evolving = any("radius_rsun" in ks for ks in star_keys.values())

        def star_series(name: str, key: str) -> np.ndarray | None:
            return R.get(f"stars/{name}/{key}") if key in star_keys.get(name, []) else None

        def ends_or(name, key, attr):
            if key in star_keys.get(name, []):
                a, b = R.ends(f"stars/{name}/{key}")
                if a is not None or b is not None:
                    return a, b
            a = _f(star_attrs.get(name, {}).get(attr))
            return a, a

        backend = str(attrs.get("backend", ""))
        event = str(attrs.get("event", "") or "")
        event_t = _f(attrs.get("event_time_yr"))
        dissipative = bool(attrs.get("dissipative", False)) or bool(attrs.get("gas_phase", False))
        mass_loss = bool(attrs.get("mass_loss", False))
        e_max = _nanmax(np.abs(energy)) if energy is not None and len(energy) else None
        l_max = _nanmax(np.abs(angmom)) if angmom is not None and len(angmom) else None
        masses_cfg = dict(cfg.masses) if cfg is not None else {}

        # star table: start -> end (first and last values only)
        star_rows = []
        for name in tree["stars"] or sorted(star_attrs):
            sa = star_attrs.get(name, {})
            m0, m1 = ends_or(name, "mass_msun", "mass_msun")
            if m0 is None and name in masses_cfg:
                m0 = m1 = float(masses_cfg[name])
            phase = sa.get("phase")
            if "phase" in star_keys.get(name, []):
                from ..evolution import PHASES

                p0, p1 = R.ends(f"stars/{name}/phase")
                if p0 is not None and p1 is not None:
                    phase = f"{PHASES[int(p0)]} → {PHASES[int(p1)]}"
            star_rows.append(
                {
                    "name": name,
                    "mass": [m0, m1],
                    "radius": list(ends_or(name, "radius_rsun", "radius_rsun")),
                    "teff": list(ends_or(name, "teff_k", "teff_k")),
                    "roche_fill": _f(sa.get("roche_fill_max")),
                    "prot": list(ends_or(name, "P_rot_d", "P_rot_d")),
                    "obliquity": list(ends_or(name, "obliquity_deg", "obliquity_deg")),
                    "k2": _f(sa.get("apsidal_k2")),
                    "rg2": _f(sa.get("rg2")),
                    "lag_s": _f(sa.get("lag_s")),
                    "phase": str(phase) if phase else "",
                    "age": _f(sa.get("age_myr")),
                }
            )

        # mergers followed by the run (physics.mergers): the product takes the mass of the merger,
        # the other star ends there
        mergers = [m for m in system.get("_mergers", []) if isinstance(m, dict)]
        rows_by = {r["name"]: r for r in star_rows}
        for mg in mergers:
            a_, b_ = (mg.get("stars") or ["", ""])[:2]
            when = f"t = {fmt(mg.get('time_yr'))} yr"
            if a_ in rows_by:
                r = rows_by[a_]
                r["mass"] = [r["mass"][0], _f(mg.get("mass_msun"))]
                if mg.get("radius_rsun") is not None:  # the product on its new track
                    r["radius"] = [r["radius"][0], _f(mg.get("radius_rsun"))]
                r["phase"] = (
                    r["phase"] + "; " if r["phase"] else ""
                ) + f"merged with {b_} at {when}"
            if b_ in rows_by:
                r = rows_by[b_]
                r["mass"] = [r["mass"][0], None]
                r["phase"] = (
                    r["phase"] + "; " if r["phase"] else ""
                ) + f"merged into {a_} at {when}"

        # eclipses of the two-star binaries: |cos i| < (R1 + R2) / a
        eclipses = []
        t_full = R.get("time_yr") if full else None
        for b in tree["binaries"]:
            if b not in orbit_names:
                continue
            kids = next(x["children"] for x in tree["orbits"] if x["name"] == b)
            stored = "eclipse_fraction" in binaries.get(b, {})
            frac = _f(binaries[b]["eclipse_fraction"]) if stored else None
            need = full or not stored or b == tree["binaries"][-1]
            if not need:
                eclipses.append({"binary": b, "stars": kids, "fraction": frac})
                continue
            inc = R.get(f"orbits/{b}/inc")
            a = R.get(f"orbits/{b}/a")
            if inc is None or a is None:
                continue
            rs = 0.0
            for k in kids:
                r = star_series(k, "radius_rsun")
                rs = rs + (
                    np.asarray(r, float)
                    if r is not None and len(r) == len(a)
                    else float(star_attrs.get(k, {}).get("radius_rsun", 0.0) or 0.0)
                )
            rsum = np.asarray(rs, float) * RSUN_AU
            with np.errstate(invalid="ignore", divide="ignore"):
                ratio = np.clip(rsum / np.asarray(a, float), 0.0, 1.0)
                lim = np.degrees(np.arccos(ratio)) * np.ones(len(a))
                ecl = np.abs(np.cos(np.asarray(inc, float))) < ratio
            if not stored:
                frac = _f(float(np.mean(ecl))) if len(ecl) else None
            inc_deg = np.degrees(np.asarray(inc, float))
            e = {
                "binary": b,
                "stars": kids,
                "fraction": frac,
                "limit_now": _f(lim[0]) if len(lim) else None,
                "spark": _arr(inc_deg, _down(len(inc_deg), 120), digits=4),
            }
            if full:
                # windows as [t_start, t_end] intervals on the full time grid
                wins = []
                edges = np.flatnonzero(np.diff(np.concatenate(([0], ecl.astype(np.int8), [0]))))
                for s, q in zip(edges[::2], edges[1::2] - 1, strict=True):
                    wins.append([_f(t_full[s]), _f(t_full[q])])
                idx = np.unique(_decimate([inc_deg, lim], MAX_POINTS))
                e.update(
                    windows=wins[:400],
                    t=_arr(t_full, idx, digits=10),
                    inc=_arr(inc_deg, idx, digits=6),
                    limit=_arr(lim, idx, digits=5),
                    inc_range=[_f(x) for x in _range(inc_deg)],
                )
            eclipses.append(e)

        mutual0 = None
        if mutual_names:
            mutual0 = (mutual_names[0], np.degrees(R.get(f"mutual_inclination/{mutual_names[0]}")))

        if full:
            data = _run_page(
                R,
                tree,
                t_full,
                orbit_names,
                mutual_names,
                star_keys,
                star_attrs,
                masses_cfg,
                energy,
                n,
            )
    # secular table, one row per inner orbit in its parent
    secular = []
    for d in diags:
        b = binaries.get(d.get("inner"), {})
        secular.append(
            {
                "pair": f"{d.get('inner')} in {d.get('outer')}",
                "i_mut0": _f(d.get("i_mut0_deg")),
                "e_theory": [_f(d.get("quad_e_min")), _f(d.get("quad_e_max"))],
                "e_run": [_f(d.get("e_min")), _f(d.get("e_max_true", d.get("e_max")))],
                "prec_est": _f(d.get("t_prec_yr")),
                "prec_meas": _f(d.get("prec_period_measured_yr")),
                "apsidal_meas": _f(b.get("U_measured_yr")),
                "apsidal_own": _f(b.get("U_expected_yr")),
                "apsidal_tidal": _f(b.get("U_tidal_yr")),
                "apsidal_gr": _f(b.get("U_gr_yr")),
                "zlk_window": bool(d.get("in_zlk_window", False)),
                "t_zlk": _f(d.get("t_zlk_yr")),
                "ma01": _f(d.get("stability_ratio")),
                "v22": _f(d.get("stability_ratio_v22")),
                "eps_sa": _f(d.get("epsilon_sa")),
                "q_min_over_radii": _f(d.get("q_min_over_radii")),
            }
        )
    ratios = [r for r in (s["ma01"] for s in secular) if r is not None]
    worst = (
        min(secular, key=lambda s: s["ma01"] if s["ma01"] is not None else math.inf)
        if ratios
        else None
    )
    approx = [s["pair"] for s in secular if (s["eps_sa"] or 0) > 0.01]

    # facts
    facts = []
    if event:
        facts.append(
            {
                "k": "Stopping event",
                "v": event.replace("_", " "),
                "note": f"at t = {fmt(event_t)} yr" if event_t is not None else "",
                "tone": "warn",
            }
        )
    else:
        facts.append(
            {
                "k": "Stopping event",
                "v": "none",
                "note": f"reached t = {fmt(t_end)} yr",
                "tone": "",
            }
        )
    for mg in mergers:
        age = _f(mg.get("age_myr"))
        facts.append(
            {
                "k": "Merger"
                + (f" ({mg.get('kind', '').replace('_', ' ')})" if mg.get("kind") else ""),
                "v": f"{' + '.join(mg.get('stars', []))} → {fmt(mg.get('mass_msun'), 4)} M☉",
                "note": f"at t = {fmt(mg.get('time_yr'))} yr; {fmt(mg.get('mass_lost_msun'), 3)} M☉ "
                "ejected"
                + (f"; product rejuvenated to {fmt(age, 4)} Myr on its new track" if age else "")
                + (f"; the run continues as {mg['hierarchy']}" if mg.get("hierarchy") else "")
                + (f" ({mg['stop']})" if mg.get("stop") else ""),
                "tone": "warn",
            }
        )
    for og in system.get("_onsets", []) or []:  # physics.mass_transfer.dynamic_onset
        if not isinstance(og, dict):
            continue
        how = (
            "a → a (1 − e²), angular momentum kept"
            if og.get("circularisation") == "angular_momentum"
            else "a → a (1 − e), donor at its lobe"
        )
        facts.append(
            {
                "k": f"Mass transfer started by the dynamics (orbit {og.get('orbit', '?')})",
                "v": f"P {fmt(og.get('P_before_d'), 4)} → {fmt(og.get('P_after_d'), 4)} d",
                "note": f"at t = {fmt(og.get('time_yr'))} yr {'+'.join(og.get('stars', []))} touched "
                f"its Roche lobe at e = {fmt(og.get('e_onset'), 3)}; circularised ({how})"
                + (f"; not continued: {og['stop']}" if og.get("stop") else ""),
                "tone": "warn",
            }
        )
    facts.append(
        {
            "k": "Hierarchy changes",
            "v": str(len(transitions)),
            "note": f"configured {configured}" if configured else "",
            "tone": "warn" if transitions else "",
        }
    )
    if eclipses:
        facts.append(
            {
                "k": "Eclipsing (fraction of run)",
                "v": " · ".join(
                    f"{e['binary']} {fmt(100 * e['fraction'], 3)} %"
                    for e in eclipses
                    if e["fraction"] is not None
                ),
                "note": "observer along z",
                "tone": "",
            }
        )
    if mutual0 is not None:
        k0, d = mutual0
        facts.append(
            {
                "k": f"Mutual inclination {k0.replace('-', '–')}",
                "v": fmt_range(_nanmin(d), _nanmax(d), 3, "°"),
                "note": f"start {fmt(_f(d[0]) if len(d) else None, 3)}°",
                "tone": "",
            }
        )
    if worst is not None:
        facts.append(
            {
                "k": "Lowest stability ratio",
                "v": fmt(worst["ma01"], 3),
                "note": f"MA01, {worst['pair']}",
                "tone": _tone_ratio(worst["ma01"]),
            }
        )
    if e_max is not None:
        facts.append(
            {
                "k": "Max |ΔE/E|",
                "v": fmt(e_max, 2),
                "note": f"|ΔL/L| {fmt(l_max, 2)}" if l_max is not None else "",
                "tone": _tone_energy(e_max, dissipative or mass_loss),
            }
        )
    mt = [
        k
        for k in system
        if k.startswith("mass_transfer_") and k.endswith("_t_onset_yr") and system[k] is not None
    ]
    for k in mt:
        orb = k.split("_")[2]
        facts.append(
            {
                "k": f"Mass transfer in {orb}",
                "v": f"from {fmt(_f(system[k]))} yr",
                "note": f"donor {system.get(f'mass_transfer_{orb}_donor', '')}",
                "tone": "",
            }
        )
        e_tr = _f(system.get(f"mass_transfer_{orb}_e_max_transferring"))
        if e_tr is not None and e_tr > 0.1:
            facts.append(
                {
                    "k": f"Eccentric while transferring ({orb})",
                    "v": f"e up to {fmt(e_tr, 2)}",
                    "note": "the transfer is modelled for a circular orbit - physics.tides "
                    "would circularise it",
                    "tone": "warn",
                }
            )

    quality = []
    if e_max is not None:
        quality.append(
            {
                "k": "energy conservation",
                "v": fmt(e_max, 2)
                if not (dissipative or mass_loss)
                else "not conserved (dissipation)",
                "tone": _tone_energy(e_max, dissipative or mass_loss),
            }
        )
    if l_max is not None:
        quality.append(
            {
                "k": "angular momentum",
                "v": fmt(l_max, 2),
                "tone": "" if (dissipative or mass_loss) else _tone_energy(l_max),
            }
        )
    if approx:
        quality.append(
            {
                "k": "secular estimates",
                "v": "approximate (" + ", ".join(approx) + ")",
                "tone": "warn",
            }
        )
    quality.append({"k": "stars evolve", "v": "yes" if evolving else "no (fixed)", "tone": ""})
    if backend:
        quality.append({"k": "backend", "v": backend, "tone": ""})

    span = f"{fmt(t0)} → {fmt(t_end)} yr" if n else ""
    physics = []
    if cfg is not None:
        ph = cfg.physics
        for k in ("gr", "tides", "rotation"):
            if getattr(ph, k, False):
                physics.append({"gr": "GR"}.get(k, k))
    setup = " · ".join(
        x
        for x in (
            f"{len(tree['stars'])} stars" if tree["stars"] else "",
            "orbits " + ", ".join(o["name"] for o in tree["orbits"]) if tree["orbits"] else "",
            backend,
            ", ".join(physics),
            "evolving stars" if evolving else "",
            span,
            f"{n} outputs",
        )
        if x
    )
    out.update(
        {
            "setup": setup,
            "facts": facts,
            "quality": quality,
            "spark": None,
        }
    )
    with_inc = [e for e in eclipses if "spark" in e]
    if with_inc:
        e0 = with_inc[-1]
        out["spark"] = {"type": "line", "label": f"inclination of {e0['binary']}", "y": e0["spark"]}
    elif mutual0 is not None:
        k0, d = mutual0
        out["spark"] = {
            "type": "line",
            "label": f"i_mut {k0}",
            "y": _arr(d, _down(len(d), 120), digits=4),
        }
    for e in eclipses:
        e.pop("spark", None)
    if not full:
        return out

    out["data"] = {
        "tree": tree,
        "configured": configured,
        "transitions": transitions,
        "event": event,
        "event_time": event_t,
        "backend": backend,
        "stars": star_rows,
        "eclipses": eclipses,
        "secular": secular,
        "has_ensemble": bool(cfg is not None and getattr(cfg, "ensemble", None)),
        "has_pulsation": bool(
            cfg is not None
            and any(getattr(s, "pulsation", None) for s in (cfg.stars or {}).values())
        ),
        **data,
    }
    return out


def _run_page(R, tree, t, orbit_names, mutual_names, star_keys, star_attrs, masses_cfg, energy, n):
    """The series of a run page: elements for the system view (MAX_FRAMES samples) and the
    plotted series (MAX_POINTS each, with their extremes kept)."""
    orbits = {
        name: {k: R.get(f"orbits/{name}/{k}") for k in R.keys(f"orbits/{name}")}
        for name in orbit_names
    }
    fi = _down(n, MAX_FRAMES)
    view = {
        "t": _arr(t, fi, digits=10),
        "elements": {
            name: {
                k: _arr(o[k], fi, digits=9)
                for k in ("a", "e", "inc", "Omega", "omega", "M")
                if o.get(k) is not None
            }
            for name, o in orbits.items()
        },
        "masses": {},
        "radii": {},
    }
    for name in tree["stars"]:
        keys = star_keys.get(name, [])
        sa = star_attrs.get(name, {})
        ms = R.get(f"stars/{name}/mass_msun") if "mass_msun" in keys else None
        rs = R.get(f"stars/{name}/radius_rsun") if "radius_rsun" in keys else None
        view["masses"][name] = (
            _arr(ms, fi, digits=6)
            if ms is not None and len(ms) == n
            else (masses_cfg.get(name) or _f(sa.get("mass_msun")))
        )
        view["radii"][name] = (
            _arr(rs, fi, digits=5) if rs is not None and len(rs) == n else _f(sa.get("radius_rsun"))
        )
        view.setdefault("teff", {})[name] = _f(sa.get("teff_k"))

    def ser(y, digits, deg=False):
        if y is None or len(y) != n:
            return None
        y = np.degrees(np.asarray(y, float)) if deg else y
        return _series(t, y, MAX_POINTS, digits=digits)

    series = {
        "orbits": {
            name: {
                k: s
                for k, s in (
                    ("a", ser(o.get("a"), 7)),
                    ("e", ser(o.get("e"), 6)),
                    ("inc", ser(o.get("inc"), 6, deg=True)),
                )
                if s is not None
            }
            for name, o in orbits.items()
        },
        "mutual": {
            k: s
            for k in mutual_names
            if (s := ser(R.get(f"mutual_inclination/{k}"), 6, deg=True)) is not None
        },
        "energy": ser(np.abs(energy), 4) if energy is not None else None,
        "stars": {
            name: {
                k: s
                for k in ("radius_rsun", "mass_msun", "teff_k", "P_rot_d")
                if k in star_keys.get(name, [])
                and (s := ser(R.get(f"stars/{name}/{k}"), 6)) is not None
            }
            for name in tree["stars"]
        },
    }
    return {
        "orbits0": {
            name: {
                "P_d": _f(_period_days(o, name, tree, view["masses"])),
                "e": _f(o["e"][0]) if o.get("e") is not None and len(o["e"]) else None,
                "a_au": _f(o["a"][0]) if o.get("a") is not None and len(o["a"]) else None,
                "inc": _f(math.degrees(o["inc"][0]))
                if o.get("inc") is not None and len(o["inc"])
                else None,
            }
            for name, o in orbits.items()
        },
        "view": view,
        "series": series,
    }


def _period_days(o: dict, name: str, tree: dict, masses: dict) -> float:
    from .. import units

    def mass(x):
        if x in masses:
            v = masses[x]
            return float(v[0] if isinstance(v, list) else v or 0.0)
        kids = next((q["children"] for q in tree["orbits"] if q["name"] == x), [])
        return sum(mass(c) for c in kids)

    if o.get("a") is None or not len(o["a"]):
        return math.nan
    a = float(o["a"][0])
    m = mass(name)
    if not (a > 0 and m > 0):
        return math.nan
    return units.sma_to_period(a, m) * units.DAYS_PER_YEAR


# --------------------------------------------------------------------------- ensemble
_ROWS: dict = {}
_ROWS_LOCK = threading.Lock()


def _rows(folder: Path) -> list[dict]:
    from ..ensemble import load_rows

    src = folder / "rows.jsonl"
    if not src.exists():  # a copied output may carry only the table (summary.csv)
        src = folder / "summary.csv"
        if not src.exists():
            return []
    key = (str(src), src.stat().st_mtime_ns, src.stat().st_size)
    with _ROWS_LOCK:
        if key in _ROWS:
            return _ROWS[key]
    rows = load_rows(folder) if src.suffix == ".jsonl" else _summary_rows(src)
    with _ROWS_LOCK:
        _ROWS.clear()  # one ensemble at a time (a large one takes a lot of memory)
        _ROWS[key] = rows
    return rows


def _summary_rows(path: Path) -> list[dict]:
    """The rows of an ensemble from its summary.csv (numbers as floats, member as int)."""

    def val(v: str):
        if v == "":
            return math.nan
        try:
            return float(v)
        except ValueError:
            return v

    with path.open(newline="") as f:
        rows = [{k: val(v) for k, v in r.items()} for r in csv.DictReader(f)]
    for r in rows:
        r["member"] = int(r["member"])
    return sorted(rows, key=lambda r: (r["member"], r["direction"]))


def _ensemble_config(folder: Path):
    from ..config import Config

    meta = folder / "ensemble.yaml"
    if not meta.exists():
        return None
    try:
        d = yaml.safe_load(meta.read_text(encoding="utf-8")) or {}
        return Config.model_validate(d.get("config", d))
    except Exception:  # noqa: BLE001
        return None


def _hist(
    x: np.ndarray, w: np.ndarray, lo: float, hi: float, log: bool, bins: int = 24
) -> tuple[list, list, list]:
    ok = np.isfinite(x)
    if log:
        ok &= x > 0
    xs, ws = x[ok], w[ok]
    if log:
        xs = np.log10(xs)
        lo, hi = math.log10(lo), math.log10(hi)
    if not (hi > lo):
        hi = lo + 1.0
    edges = np.linspace(lo, hi, bins + 1)
    prior, _ = np.histogram(np.clip(xs, lo, hi), bins=edges)
    post, _ = np.histogram(np.clip(xs, lo, hi), bins=edges, weights=ws)
    p = prior / max(prior.sum(), 1)
    q = post / max(post.sum(), 1e-300)
    if log:
        edges = 10**edges
    return (
        [float(f"{e:.6g}") for e in edges],
        [float(f"{v:.5g}") for v in p],
        [float(f"{v:.5g}") for v in q],
    )


def _label(col: str) -> str:
    m = re.match(r"orbits\.(\w+)\.(\w+)$", col)
    if m:
        orb, k = m.groups()
        names = {
            "P": "period P [d]",
            "e": "eccentricity e",
            "i": "inclination i [°]",
            "Omega": "node Ω [°]",
            "omega": "periastron ω [°]",
            "a": "semi-major axis a [AU]",
            "M0": "mean anomaly M0 [°]",
        }
        return f"{orb}: {names.get(k, k)}"
    m = re.match(r"imut_(\w+)_(\w+)_0$", col)
    if m:
        return f"mutual inclination {m.group(1)}–{m.group(2)} at t = 0 [°]"
    m = re.match(r"masses\.(\w+)$", col)
    if m:
        return f"mass {m.group(1)} [Msun]"
    return col


def _observed_text(ob: dict) -> str:
    parts = []
    if ob.get("equals") is not None:
        parts.append(
            "= "
            + (
                ", ".join(map(str, ob["equals"]))
                if isinstance(ob["equals"], list)
                else str(ob["equals"])
            )
        )
    if ob.get("value") is not None:
        parts.append(f"{fmt(ob['value'])} ± {fmt(ob.get('sigma'))}")
    if ob.get("min") is not None:
        parts.append(f"≥ {fmt(ob['min'])}")
    if ob.get("max") is not None:
        parts.append(f"≤ {fmt(ob['max'])}")
    return ", ".join(parts)


def ensemble_digest(folder: Path, use: list[str] | None = None, full: bool = True) -> dict:
    """Digest of an ensemble folder. ``use``: the observed columns to condition on (default
    all of ``ensemble.observed`` of its configuration; [] = no conditioning)."""
    from ..ensemble import (
        FATES,
        _wquantile,
        effective_size,
        observation_weights,
        split_column,
        wilson,
    )

    rows = _rows(folder)
    cfg = _ensemble_config(folder)
    ens = getattr(cfg, "ensemble", None) if cfg is not None else None
    observed_cfg = {}
    if ens is not None and getattr(ens, "observed", None):
        observed_cfg = {
            k: (v.model_dump() if hasattr(v, "model_dump") else dict(v))
            for k, v in ens.observed.items()
        }
    n_target = getattr(ens, "n", None) if ens is not None else None
    if use is None:
        use = list(observed_cfg)
    obs_use = {k: observed_cfg[k] for k in use if k in observed_cfg}
    members = sorted({r["member"] for r in rows})
    dirs = [d for d in ("forward", "backward") if any(r["direction"] == d for r in rows)]
    weights: dict[int, float] = {m: 1.0 for m in members}
    weight_error = None
    if obs_use and rows:
        try:
            weights = observation_weights(rows, obs_use)
        except ValueError as exc:
            weight_error = str(exc)
    w_members = np.array([weights.get(m, 0.0) for m in members])
    n_eff = effective_size(w_members) if len(w_members) else 0.0
    conditioned = bool(obs_use) and weight_error is None

    def frac(mask: np.ndarray, w: np.ndarray) -> dict:
        total = float(np.sum(w))
        if total <= 0:
            return {"p": None, "lo": None, "hi": None}
        p = float(np.sum(w[mask]) / total)
        ne = effective_size(w) if conditioned else float(np.sum(w > 0))
        _, lo, hi = wilson(p * ne, ne)
        return {"p": p, "lo": lo, "hi": hi}

    by_dir = {}
    for d in dirs:
        rs = [r for r in rows if r["direction"] == d]
        w = np.array([weights.get(r["member"], 0.0) for r in rs])
        fate = np.array([str(r.get("fate", "")) for r in rs])
        fates = []
        for name in list(FATES) + sorted(set(fate) - set(FATES)):
            m = fate == name
            if m.any() and np.sum(w[m]) > 0 or (m.any() and not conditioned):
                fates.append({"fate": name, **frac(m, w), "n": int(m.sum())})
        ok = [i for i, r in enumerate(rs) if r.get("fate") != "error"]
        rs_ok = [rs[i] for i in ok]
        w_ok = w[ok] if len(ok) else np.array([])

        def col(name, rs_ok=rs_ok):
            return np.array(
                [r.get(name, np.nan) if r.get(name) not in (None, "") else np.nan for r in rs_ok],
                dtype=float,
            )

        orbit_names = sorted(
            {k[: -len("_e_max")] for r in rs_ok for k in r if k.endswith("_e_max")},
            key=lambda x: (
                ("A", "B", "inner", "mid", "out").index(x)
                if x in ("A", "B", "inner", "mid", "out")
                else 9
            ),
        )
        orbits = []
        for o in orbit_names:
            v = col(f"{o}_e_max")
            q = _wquantile(v, w_ok, [0.16, 0.5, 0.84])
            near = float(np.sum(w_ok[1 - v < 0.1]) / max(np.sum(w_ok), 1e-300)) if len(v) else None
            row = {"name": o, "emax": [_f(x) for x in q], "near": _f(near)}
            if any(f"{o}_didt0_deg_yr" in r for r in rs_ok):
                dv = col(f"{o}_didt0_deg_yr")
                row["didt"] = [_f(x) for x in _wquantile(dv, w_ok, [0.16, 0.5, 0.84])]
                ch = col(f"{o}_eclipse_change_yr")
                e0 = col(f"{o}_eclipsing0")
                now = (e0 == 1) & (w_ok > 0)
                if now.any():
                    lost = now & np.isfinite(ch)
                    row["change_frac"] = _f(np.sum(w_ok[lost]) / np.sum(w_ok[now]))
                    row["change_median"] = (
                        _f(_wquantile(np.abs(ch[lost]), w_ok[lost], [0.5])[0])
                        if lost.any()
                        else None
                    )
            if any(f"{o}_eclipse_fraction" in r for r in rs_ok):
                ef = col(f"{o}_eclipse_fraction")
                g = np.isfinite(ef) & (w_ok > 0)
                row["eclipse"] = _f(np.sum(w_ok[g] * ef[g]) / np.sum(w_ok[g])) if g.any() else None
            orbits.append(row)
        pairs = []
        for key in sorted(
            {k for r in rs_ok for k in r if k.startswith("imut_") and k.endswith("_0")}
        ):
            pair = key[5:-2]
            lo, hi, i0 = col(f"imut_{pair}_min"), col(f"imut_{pair}_max"), col(key)
            cross = (
                float(np.sum(w_ok[(lo < 90) & (hi > 90)]) / max(np.sum(w_ok), 1e-300))
                if len(lo)
                else None
            )
            pairs.append(
                {
                    "name": pair.replace("_", " – "),
                    "now": [_f(x) for x in _wquantile(i0, w_ok, [0.16, 0.5, 0.84])],
                    "cross": _f(cross),
                }
            )
        stab = []
        for c, label in (
            ("ma01_min", "Mardling–Aarseth"),
            ("v22_min", "Vynatheya+22"),
            ("mlp_stable", "MLP Vynatheya+23"),
            ("ghost_stable", "N-body ghost test"),
        ):
            v = col(c)
            g = np.isfinite(v)
            if g.any():
                st = v > 1.0 if c.endswith("_min") else v > 0.5
                stab.append({"name": label, **frac(st[g], w_ok[g])})
        both = None
        if any("eclipse_fraction_both" in r for r in rs_ok):
            v = col("eclipse_fraction_both")
            g = np.isfinite(v) & (w_ok > 0)
            both = _f(np.sum(w_ok[g] * v[g]) / np.sum(w_ok[g])) if g.any() else None
        by_dir[d] = {
            "n": len(rs),
            "errors": len(rs) - len(ok),
            "fates": fates,
            "orbits": orbits,
            "pairs": pairs,
            "stability": stab,
            "eclipse_both": both,
        }

    fw = by_dir.get("forward") or next(iter(by_dir.values()), None)
    facts = [
        {
            "k": "Members",
            "v": f"{len(members):,}".replace(",", " ")
            + (
                f" / {n_target:,}".replace(",", " ")
                if n_target and n_target != len(members)
                else ""
            ),
            "note": " and ".join(dirs) or "none yet",
            "tone": "",
        },
    ]
    if conditioned:
        facts.append(
            {
                "k": "Effective members (n_eff)",
                "v": fmt(n_eff, 4),
                "note": f"conditioned on {len(obs_use)} observations",
                "tone": _tone_neff(n_eff),
            }
        )
    if fw and fw["fates"]:
        top = max(fw["fates"], key=lambda x: x["p"] or 0)
        facts.append(
            {
                "k": "Most common fate",
                "v": f"{top['fate']} {fmt(100 * (top['p'] or 0), 3)} %",
                "note": f"{len(fw['fates'])} fate{'s' if len(fw['fates']) != 1 else ''} seen",
                "tone": "",
            }
        )
    if fw and fw["stability"]:
        s = fw["stability"][0]
        facts.append(
            {
                "k": f"Stable ({s['name']})",
                "v": f"{fmt(100 * (s['p'] or 0), 3)} %",
                "note": f"{fmt(100 * (s['lo'] or 0), 3)}–{fmt(100 * (s['hi'] or 0), 3)} %",
                "tone": "",
            }
        )
    ecl = [o for o in (fw or {}).get("orbits", []) if o.get("eclipse") is not None]
    if ecl:
        facts.append(
            {
                "k": "Eclipsing (time average)",
                "v": " · ".join(f"{o['name']} {fmt(100 * o['eclipse'], 2)} %" for o in ecl),
                "note": f"all at once {fmt(100 * fw['eclipse_both'], 2)} %"
                if fw.get("eclipse_both") is not None
                else "",
                "tone": "",
            }
        )
    quality = []
    if conditioned:
        quality.append(
            {"k": "n_eff with the observations", "v": fmt(n_eff, 4), "tone": _tone_neff(n_eff)}
        )
    if cfg is not None:
        quality.append({"k": "backend", "v": cfg.run.backend or "auto", "tone": ""})
    errs = sum(v["errors"] for v in by_dir.values())
    quality.append(
        {"k": "members with an error", "v": str(errs), "tone": "warn" if errs else "good"}
    )
    if weight_error:
        quality.append({"k": "observations", "v": weight_error, "tone": "bad"})
    setup_bits = []
    if cfg is not None:
        setup_bits += [
            f"{len(cfg.star_names)} stars",
            "orbits " + ", ".join(_tree(cfg)["binaries"] + [_tree(cfg)["root"]])
            if cfg.tree
            else "",
        ]
        setup_bits.append(str(cfg.run.backend or "auto backend"))
        setup_bits.append(f"{fmt(cfg.run.t_max_yr)} yr " + " and ".join(dirs))
    if ens is not None:
        setup_bits.append(f"{len(ens.vary or {})} varied parameters")
        setup_bits.append(f"{len(observed_cfg)} observations")
        if getattr(ens, "require", None):
            setup_bits.append("required: " + ", ".join(ens.require))
    out = {
        "setup": " · ".join(b for b in setup_bits if b),
        "facts": facts,
        "quality": quality,
        "n_eff": _f(n_eff),
        "conditioned": conditioned,
    }
    # spark: the most narrowed parameter, conditioned (or the first varied one)
    params = [k for k in (rows[0] if rows else {}) if "." in k]
    fwd = [r for r in rows if r["direction"] == (dirs[0] if dirs else "forward")]
    w_fwd = np.array([weights.get(r["member"], 0.0) for r in fwd])
    hists = []
    for p in params + sorted(
        {k for r in fwd[:50] for k in r if k.startswith("imut_") and k.endswith("_0")}
    ):
        x = np.array([r.get(p, np.nan) for r in fwd], dtype=float)
        if not np.isfinite(x).any():
            continue
        lo, hi = float(np.nanmin(x)), float(np.nanmax(x))
        spec = (getattr(ens, "vary", None) or {}).get(p) if ens is not None else None
        spec = spec.model_dump() if hasattr(spec, "model_dump") else spec
        log = bool(lo > 0 and isinstance(spec, dict) and spec.get("loguniform") is not None)
        e, pr, po = _hist(x, w_fwd, lo, hi, log)
        qa = _wquantile(x, np.ones(len(x)), [0.16, 0.5, 0.84])
        qc = _wquantile(x, w_fwd, [0.16, 0.5, 0.84])
        if log:
            wa, wc = (
                math.log(qa[2] / qa[0]) if qa[0] > 0 else math.nan,
                math.log(qc[2] / qc[0]) if qc[0] > 0 else math.nan,
            )
        else:
            wa, wc = qa[2] - qa[0], qc[2] - qc[0]
        narrow = wc / wa if wa and math.isfinite(wa) and wa > 0 and math.isfinite(wc) else 1.0
        hists.append(
            {
                "col": p,
                "label": _label(p),
                "log": log,
                "edges": e,
                "prior": pr,
                "post": po,
                "q_prior": [_f(v) for v in qa],
                "q_post": [_f(v) for v in qc],
                "narrow": _f(narrow),
                "drawn": "." in p,
            }
        )
    hists.sort(
        key=lambda h: (
            h["narrow"] if conditioned and h["narrow"] is not None else 1.0,
            not h["drawn"],
        )
    )
    if hists:
        h = hists[0]
        out["spark"] = {
            "type": "bars",
            "label": h["label"] + (", conditioned" if conditioned else ""),
            "y": h["post"] if conditioned else h["prior"],
        }
    else:
        out["spark"] = None
    if not full:
        return out
    outcome_cols = []
    names0 = list(rows[0]) if rows else []
    for o in [x["name"] for x in (fw or {}).get("orbits", [])]:
        outcome_cols.append((f"{o}_e_max", f"largest e of {o}", 0.0, 1.0))
    for c in names0:
        if c.endswith("_eclipse_fraction"):
            outcome_cols.append(
                (c, f"{c[: -len('_eclipse_fraction')]} eclipsing [fraction of run]", 0.0, 1.0)
            )
    if "ma01_min" in names0:
        outcome_cols.append(("ma01_min", "stability ratio MA01 (lowest)", None, None))
    outcomes = {}
    for d in dirs:
        rs = [r for r in rows if r["direction"] == d and r.get("fate") != "error"]
        w = np.array([weights.get(r["member"], 0.0) for r in rs])
        hs = []
        for c, label, lo, hi in outcome_cols:
            x = np.array(
                [r.get(c, np.nan) if r.get(c) not in (None, "") else np.nan for r in rs],
                dtype=float,
            )
            if not np.isfinite(x).any():
                continue
            l0 = float(np.nanmin(x)) if lo is None else lo
            h0 = float(np.nanpercentile(x, 99)) if hi is None else hi
            e, pr, po = _hist(x, w, l0, h0, False)
            above = int(np.sum(x > h0))
            if above:  # _hist puts them into the last bin
                label += f" · the last bin holds the {above} values above {fmt(h0, 3)}"
            hs.append({"col": c, "label": label, "edges": e, "prior": pr, "post": po})
        outcomes[d] = hs
    # members with the largest weight (forward rows: the drawn parameters)
    order = np.argsort(-w_fwd, kind="stable")[:15]
    cols = params[:8]
    w_sum = float(np.sum(w_fwd)) if len(w_fwd) else 0.0
    bw = {r["member"]: r for r in rows if r["direction"] == "backward"}
    top = []
    for k in order:
        r = fwd[int(k)]
        top.append(
            {
                "member": r["member"],
                "weight": _f(w_fwd[int(k)]),
                "values": [_f(r.get(c)) for c in cols],
                "fate": r.get("fate", ""),
                "fate_back": (bw.get(r["member"]) or {}).get("fate", ""),
                "ma01": _f(r.get("ma01_min")),
                "v22": _f(r.get("v22_min")),
                # w = exp(-chi^2/2): chi^2 of the member against the ticked observations
                "chi2": _f(-2.0 * math.log(w_fwd[int(k)])) if w_fwd[int(k)] > 0 else None,
                "share": _f(w_fwd[int(k)] / w_sum) if w_sum > 0 else None,
            }
        )
    observations = [
        {
            "col": k,
            "direction": split_column(k)[0],
            "label": split_column(k)[1],
            "constraint": _observed_text(v),
            "used": k in obs_use,
        }
        for k, v in observed_cfg.items()
    ]
    out["data"] = {
        "directions": dirs,
        "by_dir": by_dir,
        "observations": observations,
        "n_eff": _f(n_eff),
        "n_members": len(members),
        "conditioned": conditioned,
        "params": hists,
        "outcomes": outcomes,
        "top": {
            "columns": [_label(c) for c in cols],
            "cols": cols,
            "rows": top,
            "n_good": int(np.sum(w_fwd >= 0.5)) if conditioned else None,
        },
        "weight_error": weight_error,
    }
    return out


def member_draft(root: Path, folder: Path, member: int) -> dict:
    """A draft configuration for one ensemble member: the base configuration with the drawn
    values of that member (the 'vary' columns), without the ensemble section."""
    from ..ensemble import _set

    rows = _rows(folder)
    row = next(
        (r for r in rows if r["member"] == member and r["direction"] == "forward"), None
    ) or next((r for r in rows if r["member"] == member), None)
    if row is None:
        raise ValueError(f"no member {member}")
    meta = yaml.safe_load((folder / "ensemble.yaml").read_text(encoding="utf-8")) or {}
    data = meta.get("config", meta)
    data = json.loads(json.dumps(data))
    for k, v in row.items():
        if "." in k and isinstance(v, (int, float)):
            _set(data, k, float(v))
    data.pop("ensemble", None)
    return _drop_none(data)


def _drop_none(d):
    """Remove null values at every level (a null sub-object breaks the configuration form)."""
    if isinstance(d, dict):
        return {k: _drop_none(v) for k, v in d.items() if v is not None}
    if isinstance(d, list):
        return [_drop_none(v) for v in d]
    return d


# --------------------------------------------------------------------------- light curve
def lightcurve_digest(path: Path, cfg=None, job: dict | None = None, full: bool = True) -> dict:
    d = _read_csv(path)
    if "t_day" not in d or "flux" not in d:
        raise ValueError(f"{Path(path).name}: no light curve (columns t_day, flux)")
    t = _numbers(d, "t_day")
    flux = _numbers(d, "flux")
    comps = {k[5:]: _numbers(d, k) for k in d if k.startswith("flux_")}
    tree = _tree(cfg)
    n = len(t)
    cad = float(np.median(np.diff(t))) * 1440 if n > 1 else math.nan
    fmax = _nanmax(flux)
    fmin = _nanmin(flux)
    depth = -2.5 * math.log10(fmin / fmax) if fmax and fmin and fmax > 0 and fmin > 0 else math.nan
    total0 = fmax or 1.0

    def loss(c) -> float | None:
        hi, lo = _nanmax(c), _nanmin(c)
        return _f(1.0 - lo / hi) if hi and hi > 0 and lo is not None else None

    stars = []
    for name, c in comps.items():
        cm = _nanmax(c)
        stars.append(
            {
                "name": name,
                "share": _f(cm / total0) if cm is not None else None,
                "loss": loss(c),
            }
        )
    groups = []
    for b in tree["binaries"]:
        kids = next(x["children"] for x in tree["orbits"] if x["name"] == b)
        if all(k in comps for k in kids):
            light = sum(comps[k] for k in kids)
            P = None
            if cfg is not None and b in cfg.orbits and getattr(cfg.orbits[b], "P", None):
                P = float(cfg.orbits[b].P)
            lm = _nanmax(light)
            groups.append(
                {
                    "name": b,
                    "stars": kids,
                    "share": _f(lm / total0) if lm is not None else None,
                    "loss": loss(light),
                    "P_d": P,
                    "span_periods": _f((t[-1] - t[0]) / P) if P and n else None,
                }
            )
    args = (job or {}).get("args") or {}
    band = args.get("band", "TESS")
    deepest = int(np.nanargmin(flux)) if n and np.isfinite(flux).any() else None
    facts = [
        {
            "k": "Flux range",
            "v": fmt_range(_f(fmin), _f(fmax), 4),
            "note": "1 = all stars outside eclipse",
            "tone": "",
        },
        {
            "k": "Deepest minimum",
            "v": f"{fmt(_f(depth), 3)} mag",
            "note": f"at t = {fmt(_f(t[deepest]) if deepest is not None else None, 5)} d",
            "tone": "",
        },
    ]
    if groups:
        facts.append(
            {
                "k": "Light share per binary",
                "v": " · ".join(f"{g['name']} {fmt(100 * (g['share'] or 0), 3)} %" for g in groups),
                "note": " · ".join(f"{g['name']} = {' + '.join(g['stars'])}" for g in groups),
                "tone": "",
            }
        )
        facts.append(
            {
                "k": "Largest light loss per binary",
                "v": " · ".join(f"{g['name']} {fmt(100 * (g['loss'] or 0), 3)} %" for g in groups),
                "note": "eclipses inside that binary",
                "tone": "",
            }
        )
    length = _f(t[-1] - t[0]) if n else None
    out = {
        "setup": " · ".join(
            x
            for x in (
                f"{len(comps)} stars" if comps else "",
                f"band {band}",
                f"{fmt(length, 4)} d from t = {fmt(_f(t[0]) if n else None, 4)} d",
                f"cadence {fmt(_f(cad), 3)} min",
                f"{n} points",
            )
            if x
        ),
        "facts": facts,
        "quality": [],
        "spark": {
            "type": "line",
            "label": f"flux, {fmt(length, 3)} d",
            "y": _arr(flux, _decimate([flux], 160), digits=5),
        },
    }
    if not full:
        return out
    # the eclipses survive the decimation: the extremes of the total and of every star's light
    idx = _decimate([flux, *comps.values()], LC_POINTS)
    out["data"] = {
        "t": _arr(t, idx, digits=10),
        "flux": _arr(flux, idx, digits=8),
        "stars": {k: _arr(v, idx, digits=8) for k, v in comps.items()},
        "star_rows": stars,
        "binaries": groups,
        "band": band,
        "n": n,
        "n_shown": len(idx),
        "pulsation": bool(
            cfg is not None
            and any(getattr(s, "pulsation", None) for s in (cfg.stars or {}).values())
        ),
    }
    return out


# --------------------------------------------------------------------------- eclipse timings
TIMING_POINTS = 4000  # eclipses per binary and kind sent to the browser


def timings_digest(path: Path, cfg=None, full: bool = True, table: dict | None = None) -> dict:
    """Eclipse times of a run (swingersim timings). ``table``: the CSV already read (_read_csv)."""
    d = table if table is not None else _read_csv(path)
    tree = _tree(cfg)
    parent = {c: o["name"] for o in tree["orbits"] for c in o["children"]}
    binary = np.asarray(d.get("binary", np.array([], dtype=object))).astype(str)
    kind = np.asarray(d.get("kind", np.array([], dtype=object))).astype(str)
    t_all = _numbers(d, "t_obs_yr") if "t_obs_yr" in d else np.full(len(binary), np.nan)
    oc_all = _numbers(d, "o_c_d") * 1440.0 if "o_c_d" in d else np.full(len(binary), np.nan)
    imp_all_ = _numbers(d, "impact") if "impact" in d else np.full(len(binary), np.nan)
    names = sorted(set(binary.tolist()))
    per = []
    series = {}
    for b in names:
        sel = np.flatnonzero(binary == b)
        kinds = {}
        for k in ("primary", "secondary"):
            ii = sel[kind[sel] == k]
            if not len(ii):
                continue
            ii = ii[np.argsort(t_all[ii], kind="stable")]
            tt, oc, imp = t_all[ii], oc_all[ii], imp_all_[ii]
            idx = _decimate([oc], TIMING_POINTS)
            fin = oc[np.isfinite(oc)]
            kinds[k] = {
                "t": _arr(tt, idx, digits=10),
                "oc": _arr(oc, idx, digits=7),
                "impact": _arr(imp, idx, digits=4),
                "mean": _f(float(np.mean(fin))) if len(fin) else None,
                "n_total": len(ii),
                "range": [_f(x) for x in _range(oc)],
                "_impact": _range(imp),
            }
        P = _f(_numbers(d, "P_d")[sel[0]]) if len(sel) and "P_d" in d else None
        prim = kinds.get("primary")
        lo, hi = prim["range"] if prim else (None, None)
        amp = _f((hi - lo) / 2.0) if lo is not None and hi is not None else None
        p_out = None
        if (
            cfg is not None
            and parent.get(b)
            and parent[b] in cfg.orbits
            and getattr(cfg.orbits[parent[b]], "P", None)
        ):
            p_out = float(cfg.orbits[parent[b]].P) / 365.25
        kids = next((o["children"] for o in tree["orbits"] if o["name"] == b), [])
        imp_lo = [v["_impact"][0] for v in kinds.values() if v["_impact"][0] is not None]
        imp_hi = [v["_impact"][1] for v in kinds.values() if v["_impact"][1] is not None]
        for v in kinds.values():
            v.pop("_impact")
        per.append(
            {
                "binary": b,
                "n_primary": kinds.get("primary", {}).get("n_total", 0),
                "n_secondary": kinds.get("secondary", {}).get("n_total", 0),
                "P_d": P,
                "oc_range": [lo, hi],
                "semi_amplitude": amp,
                "outer": parent.get(b),
                "outer_P_yr": _f(p_out),
                "impact": [_f(min(imp_lo)), _f(max(imp_hi))] if imp_lo else [None, None],
                "T0_yr": _f(_numbers(d, "T0_yr")[sel[0]]) if len(sel) and "T0_yr" in d else None,
                "primary_star": kids[0] if kids else "",
                "secondary_star": kids[1] if len(kids) > 1 else "",
            }
        )
        series[b] = kinds
    t_lo, t_hi = _range(t_all)
    span = t_hi - t_lo if t_lo is not None else 0.0
    facts = [
        {
            "k": "Eclipses",
            "v": str(len(binary)),
            "note": " · ".join(f"{p['binary']} {p['n_primary']} + {p['n_secondary']}" for p in per),
            "tone": "",
        }
    ]
    for p in per:
        facts.append(
            {
                "k": f"{p['binary']}: O − C of primaries",
                "v": f"{fmt_range(p['oc_range'][0], p['oc_range'][1], 3)} min",
                "note": f"half range {fmt(p['semi_amplitude'], 4)} min"
                if p["semi_amplitude"] is not None
                else "",
                "tone": "",
            }
        )
    out = {
        "setup": f"eclipse times from the dynamics · {fmt(_f(span), 3)} yr · every eclipse of every eclipsing binary · times include the light travel time",
        "facts": facts,
        "quality": [],
        "spark": None,
    }
    if per and series[per[0]["binary"]].get("primary"):
        y = series[per[0]["binary"]]["primary"]["oc"]
        out["spark"] = {
            "type": "line",
            "label": f"O − C of {per[0]['binary']} primaries [min]",
            "y": y[:: max(1, len(y) // 160)],
        }
    if not full:
        return out
    out["data"] = {"binaries": per, "series": series}
    return out


# --------------------------------------------------------------------------- minima
_SHAPE = re.compile(
    r"^\s*(\S+) (primary|secondary)\s*: (\d+) minima; shape A = ([-\d.e+]+) mag, D = ([-\d.e+]+) h, "
    r"Gamma = ([-\d.e+]+), C = ([-\d.e+]+), rms ([-\d.e+]+) mag; minimum - conjunction ([-+\d.e]+) s "
    r"\(scatter ([-\d.e+]+) s\), median sigma ([-\d.e+]+) s"
)
_EPHEM = re.compile(r"^\s*(\S+): ephemeris of the minima T0 = ([-\d.e+]+) d, P = ([-\d.e+]+) d")
_SKIP = re.compile(r"skipped: (\d+) blended, (\d+) truncated")


def minima_digest(
    path: Path, log_text: str = "", seasons: Path | None = None, full: bool = True
) -> dict:
    """Minima measured on a synthetic light curve (swingersim minima). O - C in seconds."""
    d = _read_csv(path)
    shapes, ephem = {}, {}
    skipped = None
    for line in log_text.splitlines():
        m = _SHAPE.match(line)
        if m:
            b, k = m.group(1), m.group(2)
            shapes[f"{b}:{k}"] = {
                "n": int(m.group(3)),
                "A_mag": float(m.group(4)),
                "D_h": float(m.group(5)),
                "Gamma": float(m.group(6)),
                "C": float(m.group(7)),
                "rms_mag": float(m.group(8)),
                "offset_s": float(m.group(9)),
                "scatter_s": float(m.group(10)),
                "sigma_s": float(m.group(11)),
            }
        m = _EPHEM.match(line)
        if m:
            ephem[m.group(1)] = {"T0_d": float(m.group(2)), "P_d": float(m.group(3))}
        m = _SKIP.search(line)
        if m:
            skipped = {"blended": int(m.group(1)), "truncated": int(m.group(2))}
    binary = np.asarray(d.get("binary", np.array([], dtype=object))).astype(str)
    kinds = np.asarray(d.get("kind", np.array([], dtype=object))).astype(str)

    def col(name, scale=1.0):  # an empty or missing cell is None
        v = _numbers(d, name) if name in d else np.full(len(binary), np.nan)
        return [_f(x) for x in v * scale]

    cols = {
        "t": col("t_min_d"),
        "oc": col("o_c_d", 86400.0),
        "sigma": col("sigma_d", 86400.0),
        "geom": col("o_c_geom_d", 86400.0),
        "depth": col("depth_mag"),
        "width": col("width_d"),
    }
    series = {}
    for i, b in enumerate(binary.tolist()):
        s = series.setdefault(b, {}).setdefault(kinds[i], {k: [] for k in cols})
        for k, v in cols.items():
            s[k].append(v[i])
    rows = []
    if seasons is not None and seasons.exists():
        sd = _read_csv(seasons)
        n_s = len(next(iter(sd.values()), []))
        for i in range(n_s if "binary" in sd else 0):
            rows.append({k: (sd[k][i] if sd[k].dtype == object else _f(sd[k][i])) for k in sd})
    n = len(binary)
    facts = [
        {
            "k": "Minima",
            "v": str(n),
            "note": " · ".join(
                f"{b} {len(s.get('primary', {}).get('t', []))} + {len(s.get('secondary', {}).get('t', []))}"
                for b, s in series.items()
            ),
            "tone": "",
        }
    ]
    sig = [v for b in series.values() for k in b.values() for v in k["sigma"] if v is not None]
    if sig:
        facts.append(
            {
                "k": "Median σ of a minimum",
                "v": f"{fmt(float(np.median(sig)), 3)} s",
                "note": "from the fit of each minimum",
                "tone": "",
            }
        )
    if shapes:
        off = " · ".join(f"{k.replace(':', ' ')} {v['offset_s']:+.1f} s" for k, v in shapes.items())
        facts.append(
            {
                "k": "Minimum − conjunction",
                "v": off,
                "note": "systematic shift of the measured minima",
                "tone": "",
            }
        )
    for b, e in ephem.items():
        facts.append(
            {
                "k": f"Ephemeris of {b}",
                "v": f"P {e['P_d']:.8g} d",
                "note": f"T0 {e['T0_d']:.7g} d",
                "tone": "",
            }
        )
    if skipped is not None:
        facts.append(
            {
                "k": "Skipped eclipses",
                "v": str(skipped["blended"] + skipped["truncated"]),
                "note": f"{skipped['blended']} blended, {skipped['truncated']} truncated",
                "tone": "",
            }
        )
    t_lo, t_hi = _range(np.array([math.nan if x is None else x for x in cols["t"]], float))
    span = (t_hi - t_lo) / 365.25 if t_lo is not None else 0.0
    out = {
        "setup": f"minima measured on a synthetic light curve with the Mikulášek model · {fmt(span, 3)} yr",
        "facts": facts,
        "quality": [],
        "spark": None,
    }
    first = next(iter(series.values()), {}).get("primary")
    if first:
        out["spark"] = {
            "type": "dots",
            "label": f"O − C of {next(iter(series))} primaries [s]",
            "y": first["oc"][:: max(1, len(first["oc"]) // 160)],
        }
    if not full:
        return out
    out["data"] = {"series": series, "shapes": shapes, "ephemeris": ephem, "seasons": rows}
    return out


# --------------------------------------------------------------------------- validation
_BENCH = re.compile(r"\|\s*(PASS|EXPLAINED|OPEN|FAIL|SKIP)\s*\|")


_STAB_ROW = re.compile(
    r"^\s*(\S+) in (\S+) \(\+([^)]*)\).*?i_mut = ([\d.]+) deg.*?Mardling-Aarseth ([\d.]+), "
    r"Vynatheya\+22 ([\d.]+)"
)


def other_digest(command: str, text: str) -> dict:
    """Facts of the results without a page of their own, read from what the job printed:
    stability check, chaos map, scattering experiment, pulsation periods. Empty when the command
    is none of these (the page then shows the log)."""
    facts: list[dict] = []
    lines = text.splitlines()

    def tone_ratio(r: float) -> str:
        return "good" if r >= 1.5 else ("warn" if r >= 1.0 else "bad")

    if command == "stability":
        for ln in lines:
            m = _STAB_ROW.search(ln)
            if m:
                ma, v22 = float(m.group(5)), float(m.group(6))
                facts.append(
                    {
                        "k": f"{m.group(1)} in {m.group(2)}: stability ratio",
                        "v": f"MA01 {fmt(ma, 3)} · V22 {fmt(v22, 3)}",
                        "note": f"with {m.group(3)}, i_mut {m.group(4)}° · above 1: stable",
                        "tone": tone_ratio(min(ma, v22)),
                    }
                )
            elif ln.startswith("machine-learning classifier"):
                verdict = ln.split(":", 1)[1].strip()
                facts.append(
                    {
                        "k": "machine-learning classifier",
                        "v": verdict,
                        "note": "Vynatheya et al. (2023)",
                        "tone": "good" if verdict == "stable" else "warn",
                    }
                )
            elif ln.startswith("N-body ghost test"):
                m2 = re.search(
                    r"\((\d+) outer orbits.*?\): (\w+)(.*?)(?:, max relative divergence ([\d.e+-]+))?\s*(?:\[|$)",
                    ln,
                )
                ok = "stable" in ln and "unstable" not in ln
                facts.append(
                    {
                        "k": "N-body ghost test",
                        "v": "stable" if ok else "unstable",
                        "note": (
                            f"{m2.group(1)} outer orbits"
                            + (f", largest divergence {m2.group(4)}" if m2 and m2.group(4) else "")
                            if m2
                            else ln[:100]
                        ),
                        "tone": "good" if ok else "bad",
                    }
                )
    elif command == "chaosmap":
        head = next((x for x in lines if x.startswith("chaos map:")), None)
        if head:
            facts.append(
                {"k": "chaos map", "v": head.split(":", 1)[1].strip(), "note": "", "tone": ""}
            )
        for key, label, tone in (
            ("regular (MEGNO < 4)", "regular", "good"),
            ("chaotic (MEGNO >= 4)", "chaotic", "warn"),
            ("escape or collision", "escape or collision", "bad"),
        ):
            ln = next((x for x in lines if x.strip().startswith(key)), None)
            if ln:
                n = int(ln.split(":")[1])
                facts.append(
                    {
                        "k": label,
                        "v": str(n),
                        "note": key if key != label else "",
                        "tone": tone if n else "",
                    }
                )
        ln = next((x for x in lines if "with regular motion from there outwards" in x), None)
        if ln:
            k, v = ln.strip().split(":", 1)
            facts.append({"k": k, "v": v.strip(), "note": "", "tone": ""})
    elif command == "scatter":
        head = next((x for x in lines if re.search(r": \d+ encounters$", x)), None)
        if head:
            facts.append(
                {
                    "k": "encounters",
                    "v": head.rsplit(":", 1)[1].split()[0],
                    "note": head.split(":")[0],
                    "tone": "",
                }
            )
        rows, on = [], False
        for x in lines:
            if x.strip().startswith("outcome"):
                on = True
                continue
            if on:
                m = re.match(r"^\s+(.+?)\s+(\d+)\s+([\d.]+)\s+\[([\d.]+), ([\d.]+)\]\s+(\S+)", x)
                if not m:
                    if rows:
                        break
                    continue
                rows.append(m)
        for m in rows[:4]:
            facts.append(
                {
                    "k": m.group(1),
                    "v": f"{fmt(100 * float(m.group(3)), 3)} %",
                    "note": f"{m.group(2)} encounters · 68 %: {fmt(100 * float(m.group(4)), 3)}–{fmt(100 * float(m.group(5)), 3)} % · σ {m.group(6)} AU²",
                    "tone": "",
                }
            )
        ln = next((x for x in lines if x.strip().startswith("energy error:")), None)
        if ln:
            mx = re.search(r"max ([\d.e+-]+)", ln)
            facts.append(
                {
                    "k": "energy error",
                    "v": ln.split(":", 1)[1].strip(),
                    "note": "",
                    "tone": "good" if mx and float(mx.group(1)) < 1e-6 else "warn",
                }
            )
    elif command == "view":
        m = next(
            (
                re.search(r"(\d+) viewers written, (\d+) up to date, (\d+) failed", x)
                for x in reversed(lines)
                if "viewers written" in x
            ),
            None,
        )
        if m:
            facts.append(
                {
                    "k": "viewers written",
                    "v": m.group(1),
                    "note": f"{m.group(2)} were up to date",
                    "tone": "",
                }
            )
            facts.append(
                {
                    "k": "failed",
                    "v": m.group(3),
                    "note": "see the log",
                    "tone": "warn" if int(m.group(3)) else "",
                }
            )
        else:
            n = sum(x.startswith("saved ") for x in lines)
            facts.append({"k": "viewers written", "v": str(n), "note": "", "tone": ""})
    elif command == "pulsation":
        for x in lines:
            m = re.match(r"^(m\d+): (\S+) \((\w+)\)\s+P = (\S+) d\s+dP/dt = (\S+) s/yr", x)
            if m:
                facts.append(
                    {
                        "k": f"{m.group(1)}: {m.group(2).replace('_', ' ')} ({m.group(3)})",
                        "v": f"P {m.group(4)} d",
                        "note": f"dP/dt {m.group(5)} s/yr from the evolution of the star (t = 0)",
                        "tone": "",
                    }
                )
    return {"facts": facts} if facts else {}


def benchmarks_digest(log_text: str, full: bool = True) -> dict:
    counts = {"PASS": 0, "EXPLAINED": 0, "FAIL": 0, "SKIP": 0}
    rows = []
    for line in log_text.splitlines():
        m = _BENCH.search(line)
        if m:
            counts[{"OPEN": "EXPLAINED"}.get(m.group(1), m.group(1))] += 1  # OPEN: older logs
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            rows.append(cells)
    opens = [
        line.strip().split(" ", 1)[1]
        for line in log_text.splitlines()
        if line.startswith(("  EXPLAINED ", "  OPEN "))
    ]
    total = sum(counts.values())
    facts = [
        {
            "k": k,
            "v": str(v),
            "note": "",
            "tone": {"PASS": "good", "EXPLAINED": "warn", "FAIL": "bad", "SKIP": ""}[k]
            if v
            else "",
        }
        for k, v in counts.items()
        if v or k in ("PASS", "FAIL")
    ]
    out = {
        "setup": f"benchmarks against published results · {total}",
        "facts": facts,
        "quality": [],
        "spark": {
            "type": "stack",
            "label": "PASS · EXPLAINED · FAIL",
            "y": [counts["PASS"], counts["EXPLAINED"], counts["FAIL"]],
        }
        if total
        else None,
    }
    if full:
        rep = next(
            (ln.split()[1] for ln in reversed(log_text.splitlines()) if ln.startswith("report ")),
            None,
        )
        out["data"] = {"counts": counts, "open": opens, "rows": rows[:400], "report": rep}
    return out


# --------------------------------------------------------------------------- fit (swingersim fit)
_TAU = re.compile(r"autocorrelation time tau max ([0-9.]+|nan) steps \(steps / tau ([0-9.]+|nan)\)")
_CHI = re.compile(r"chi\^2 ([-\d.e+]+) -> ([-\d.e+]+) \(reduced ([-\d.e+]+)")
_COUNTS = re.compile(
    r"^\s+(\d+) minima(.*?), (\d+) free parameters, (\d+) model evaluations, (\d+) s", re.MULTILINE
)


def _read_rows(path: Path) -> list[dict]:
    from ..tables import read_rows

    return read_rows(path)


def _comments(path: Path) -> list[str]:
    with path.open(encoding="utf-8") as f:
        return [ln.lstrip()[1:].strip() for ln in f if ln.lstrip().startswith("#")]


def _num(s):
    try:
        v = float(s)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _year(bjd):
    return None if bjd is None else 2000.0 + (bjd - 2451544.5) / 365.25


def _rms(v) -> float | None:
    v = [x for x in v if x is not None]
    return _f(math.sqrt(sum(x * x for x in v) / len(v))) if v else None


def _fit_bounds(cfg_path: Path | None, rows: list[dict]) -> dict:
    """Bounds of the free parameters of a fit: the columns lower / upper of stage2_params.csv
    when it has them, else fit.stage2.bounds of the fit's configuration with the default bounds
    of swingersim fit (around the start values). {} when neither is known."""
    out = {}
    for r in rows:
        lo, hi = _num(r.get("lower")), _num(r.get("upper"))
        if lo is not None or hi is not None:
            out[r.get("parameter")] = (lo, hi)
    if out or cfg_path is None:
        return out
    try:
        data = yaml.safe_load(Path(cfg_path).read_text(encoding="utf-8")) or {}
        stage2 = (data.get("fit") or {}).get("stage2")
        if not isinstance(stage2, dict):  # not the configuration of the fit: bounds unknown
            return out
        given = stage2.get("bounds") or {}
        given = {str(k): (float(v[0]), float(v[1])) for k, v in given.items()}
    except (OSError, yaml.YAMLError, AttributeError, TypeError, ValueError, IndexError):
        return out
    from ..dynfit import default_bounds

    names = [r.get("parameter") for r in rows]
    starts = [_num(r.get("start")) for r in rows]
    for name, s0 in zip(names, starts, strict=True):
        if not name or s0 is None:
            continue
        try:
            lo, hi = default_bounds([name], [s0], given)
        except (ValueError, KeyError, IndexError):
            continue
        out[name] = (_f(lo[0]), _f(hi[0]))
    return out


def _on_bound(v, e, lo, hi) -> bool:
    """Whether a fitted value sits on a bound. Without known bounds only an error of exactly 0
    (the Jacobian column vanished) says so: a small error alone is a precise parameter."""
    if v is None:
        return False
    if lo is None and hi is None:
        return e is not None and e == 0
    width = (hi - lo) if lo is not None and hi is not None else None
    tol = 1e-6 * width if width else 1e-9 * max(abs(v), 1.0)
    return bool((lo is not None and v - lo <= tol) or (hi is not None and hi - v <= tol))


def fit_digest(
    folder: Path, full: bool = True, root: Path | None = None, cfg_path: Path | None = None
) -> dict:
    """The outputs of swingersim fit in a folder: stage 2 (joint N-body fit; stage2_*), stage 1
    (O - C models; stage1_*) and the anchoring (anchor_report.txt). ``cfg_path``: the
    configuration of the fit (for the bounds of the free parameters)."""
    rep = folder / "stage2_report.txt"
    text = rep.read_text(encoding="utf-8", errors="replace") if rep.exists() else ""
    s1 = folder / "stage1_summary.json"
    anchor = folder / "anchor_report.txt"
    facts, quality, warnings = [], [], []
    stage = (
        "stage 2"
        if text
        else ("stage 1" if s1.exists() else ("anchoring" if anchor.exists() else ""))
    )
    chi = _CHI.search(text)
    m = _COUNTS.search(text)
    setup = "swingersim fit"
    if m:
        setup = f"joint N-body fit · {m.group(1)} minima{m.group(2)} · {m.group(3)} free parameters · {m.group(4)} model evaluations"
    if chi:
        red = float(chi.group(3))
        # far below 1: the errors are larger than the scatter (e.g. red-noise factors of light
        # curves), so the parameter errors are conservative; far above: model or errors wrong
        tone = "warn" if red < 0.5 else "good" if red < 2 else ("warn" if red < 10 else "bad")
        low = " - below 0.5: the errors used are larger than the scatter, the parameter errors are conservative"
        facts.append(
            {
                "k": "χ² start → fitted",
                "v": f"{fmt(float(chi.group(1)), 4)} → {fmt(float(chi.group(2)), 4)}",
                "note": f"reduced {fmt(red, 4)}" + (low if red < 0.5 else ""),
                "tone": tone,
            }
        )
        quality.append({"k": "reduced χ²", "v": fmt(red, 3), "tone": tone})
    for ln in text.splitlines():
        s = ln.strip()
        if s.startswith("warning:"):
            warnings.append(s[len("warning:") :].strip())
        if s.startswith("strongly correlated parameters"):
            warnings.append(s)
        if s.startswith("note:") and "held at its configured value" in s:
            warnings.append(s[len("note:") :].strip())
    tau = _TAU.search(text)
    if tau:  # emcee convergence: the chain after burn-in (half the steps) should be >= 10 tau
        ratio = _num(tau.group(2))
        ok = ratio is not None and ratio >= 20.0
        facts.append(
            {
                "k": "MCMC convergence",
                "v": f"steps / τ = {fmt(ratio, 3)}",
                "note": f"autocorrelation time τ = {tau.group(1)} steps; "
                + (
                    "converged (≥ 10 τ after burn-in)"
                    if ok
                    else "too short: set fit.stage2.mcmc_tau_target: 50 (or more mcmc_steps)"
                ),
                "tone": "good" if ok else "warn",
            }
        )
        quality.append(
            {"k": "MCMC steps / τ", "v": fmt(ratio, 3), "tone": "good" if ok else "warn"}
        )
    data: dict = {"stage": stage, "report": text[-TEXT_LIMIT_FIT:], "warnings": warnings}
    # parameters
    params = []
    pp = folder / "stage2_params.csv"
    if pp.exists():
        rows_p = _read_rows(pp)
        bounds = _fit_bounds(cfg_path, rows_p)
        for r in rows_p:
            v, e, s0 = _num(r.get("value")), _num(r.get("error")), _num(r.get("start"))
            mu, sg = _num(r.get("prior_mean")), _num(r.get("prior_sigma"))
            name = r.get("parameter")
            lo, hi = bounds.get(name, (None, None))
            p = {
                "name": name,
                "start": s0,
                "value": v,
                "error": e,
                "unit": r.get("unit") or "",
                "prior": [mu, sg] if mu is not None and sg else None,
                "prior_dev": _f((v - mu) / sg) if mu is not None and sg and v is not None else None,
                "shift_sigma": _f((v - s0) / e) if e and v is not None and s0 is not None else None,
                "bounds": [lo, hi] if lo is not None or hi is not None else None,
                "tiny_error": _on_bound(v, e, lo, hi),
            }
            # columns of newer fits (emcee percentiles, the linear error next to them)
            for k in ("error_linear", "p16", "p50", "p84"):
                if _num(r.get(k)) is not None:
                    p[k] = _num(r.get(k))
            params.append(p)
        for p in params:
            if p["tiny_error"]:
                warnings.append(
                    f"{p['name']}: error {fmt(p['error'], 3)} - the parameter sits on the bound "
                    f"({fmt(p['bounds'][0]) if p['bounds'] else '?'} … "
                    f"{fmt(p['bounds'][1]) if p['bounds'] else '?'})"
                    if p["bounds"]
                    else f"{p['name']}: error 0 - the parameter probably sits on the bound"
                )
            if p["error"] is not None and p["error"] > 1e3 and p["unit"] == "deg":
                warnings.append(
                    f"{p['name']}: error {fmt(p['error'], 3)}° - not determined by the data (a degeneracy)"
                )
        devs = [p for p in params if p["prior_dev"] is not None and abs(p["prior_dev"]) > 2]
        if devs:
            facts.append(
                {
                    "k": "priors > 2σ off",
                    "v": str(len(devs)),
                    "note": ", ".join(f"{p['name']} {p['prior_dev']:+.1f}σ" for p in devs[:3]),
                    "tone": "warn",
                }
            )
    data["params"] = params
    # minima
    oc = {}
    rp = folder / "stage2_residuals.csv"
    if rp.exists():
        rows = _read_rows(rp)
        for r in rows:
            b = r.get("pair")
            t = _num(r.get("bjd_tdb"))
            d = oc.setdefault(
                b,
                {
                    "t": [],
                    "kind": [],
                    "res": [],
                    "res_start": [],
                    "sigma": [],
                    "source": [],
                    "status": [],
                },
            )
            d["t"].append(_f(_year(t)))
            d["kind"].append(r.get("type"))
            d["res"].append(_num(r.get("o_minus_model_min")))
            ms = _num(r.get("model_start_bjd"))
            d["res_start"].append(
                _f((t - ms) * 1440.0) if t is not None and ms is not None else None
            )
            d["sigma"].append(_f((_num(r.get("sigma_d")) or 0.0) * 1440.0))
            d["source"].append(r.get("source") or "")
            d["status"].append(r.get("status") or "fitted")
        per = []
        for b, d in oc.items():
            fit_ = [k for k, st in enumerate(d["status"]) if st == "fitted"]
            rp_ = [d["res"][k] for k in fit_ if d["kind"][k] == "p" and d["res"][k] is not None]
            rs_ = [d["res"][k] for k in fit_ if d["kind"][k] == "s" and d["res"][k] is not None]
            srcs = {}
            for k in fit_:
                srcs.setdefault(d["source"][k], []).append(d["res"][k])
            per.append(
                {
                    "pair": b,
                    "n": len(fit_),
                    "n_out": len(d["status"]) - len(fit_),
                    "rms": _rms([d["res"][k] for k in fit_]),
                    "rms_start": _rms([d["res_start"][k] for k in fit_]),
                    "mean_p": _f(float(np.mean(rp_))) if rp_ else None,
                    "mean_s": _f(float(np.mean(rs_))) if rs_ else None,
                    "sources": [
                        {"source": s, "n": len(v), "rms": _rms(v)} for s, v in sorted(srcs.items())
                    ],
                }
            )
            if rp_ and rs_ and len(rp_) >= 3 and len(rs_) >= 3:
                err = math.hypot(
                    np.std(rp_, ddof=1) / math.sqrt(len(rp_)),
                    np.std(rs_, ddof=1) / math.sqrt(len(rs_)),
                )
                diff = float(np.mean(rp_) - np.mean(rs_))
                if err > 0 and abs(diff / err) > 4 and abs(diff) > 1.0:
                    msg = (
                        f"{b}: primary and secondary minima are offset by {diff:+.1f} min ({diff / err:+.0f} sigma) - "
                        "e cos(omega) of the pair is probably wrong (a local minimum); restart near the other sign"
                    )
                    if msg not in warnings:
                        warnings.insert(0, msg)
            facts.append(
                {
                    "k": f"{b}: rms of the minima",
                    "v": f"{fmt(per[-1]['rms'], 3)} min",
                    "note": f"start {fmt(per[-1]['rms_start'], 3)} min · {per[-1]['n']} minima"
                    + (f" ({per[-1]['n_out']} left out)" if per[-1]["n_out"] else ""),
                    "tone": "",
                }
            )
        data["minima"] = {"series": oc, "pairs": per}
        if full is False:
            data.pop("minima", None)
    if any("primary and secondary minima" in w for w in warnings):
        quality.append({"k": "primary vs secondary", "v": "offset", "tone": "bad"})
    # radial velocities
    rvp = folder / "stage2_rv.csv"
    if rvp.exists() and full:
        stars = {}
        for r in _read_rows(rvp):
            d = stars.setdefault(
                r.get("star"),
                {
                    "t": [],
                    "bjd": [],
                    "rv": [],
                    "model": [],
                    "res": [],
                    "sigma": [],
                    "source": [],
                    "status": [],
                },
            )
            t = _num(r.get("bjd_tdb"))
            d["t"].append(_f(_year(t)))
            d["bjd"].append(t)
            d["rv"].append(_num(r.get("rv")))
            d["model"].append(_num(r.get("model")))
            d["res"].append(_num(r.get("o_minus_model")))
            d["sigma"].append(_num(r.get("sigma")))
            d["source"].append(r.get("source") or "")
            d["status"].append(r.get("status") or "fitted")
        data["rv"] = stars
    if rvp.exists():
        res_all = []
        for r in _read_rows(rvp):
            if (r.get("status") or "fitted") == "fitted":
                res_all.append(_num(r.get("o_minus_model")))
        facts.append(
            {
                "k": "radial velocities: rms",
                "v": f"{fmt(_rms(res_all), 3)} km/s",
                "note": f"{len(res_all)} velocities",
                "tone": "",
            }
        )
    # light curves
    lcp = folder / "stage2_lc.csv"
    if lcp.exists():
        segs = {}
        for r in _read_rows(lcp):
            key = f"{r.get('source', '')}:{r.get('segment', '')}"
            d = segs.setdefault(key, {"t": [], "flux": [], "model": [], "start": [], "res": []})
            d["t"].append(_num(r.get("bjd_tdb")))
            d["flux"].append(_num(r.get("flux")))
            d["model"].append(_num(r.get("model")))
            d["start"].append(_num(r.get("model_start")))
            d["res"].append(_num(r.get("o_minus_model")))
            d.setdefault("status", []).append(r.get("status") or "fitted")
        stats = []
        for k, d in segs.items():  # rms of the fitted points (outliers left out)
            rr = [
                x
                for x, st in zip(d["res"], d["status"], strict=True)
                if x is not None and st == "fitted"
            ]
            stats.append(
                {"segment": k, "n": len(rr), "rms_ppt": _f(1e3 * _rms(rr)) if rr else None}
            )
        # the rms of all fitted points (segments weighted by their number of points; the
        # outliers left out of the fit are left out here too, as in the report)
        all_res = [
            x
            for d in segs.values()
            for x, st in zip(d["res"], d["status"], strict=True)
            if x is not None and st == "fitted"
        ]
        n_out = sum(st != "fitted" for d in segs.values() for st in d["status"])
        facts.append(
            {
                "k": "light curves: rms",
                "v": f"{fmt(_f(1e3 * _rms(all_res)) if all_res else None, 3)} ppt",
                "note": f"{len(segs)} segment(s), {len(all_res)} points"
                + (f" ({n_out} outliers left out)" if n_out else ""),
                "tone": "",
            }
        )
        if full:
            for d in segs.values():
                if len(d["t"]) > LC_POINTS // 4:
                    nan = [math.nan if x is None else x for x in d["res"]]
                    fl = [math.nan if x is None else x for x in d["flux"]]
                    idx = _decimate([nan, fl], LC_POINTS // 4)
                    for c in d:
                        d[c] = [d[c][i] for i in idx]
                    d["n_total"] = len(nan)
            data["lc"] = {"segments": segs, "stats": stats}
    # inclinations
    ip = folder / "stage2_inclination.csv"
    if ip.exists():
        rows = _read_rows(ip)
        inc = [
            {
                "t": _f(_year(_num(r.get("time")))),
                "i_obs": _num(r.get("i_obs")),
                "sigma": _num(r.get("sigma")),
                "limit": r.get("limit") or "",
                "i_start": _num(r.get("i_start")),
                "i_model": _num(r.get("i_model")),
                "i_limit": _num(r.get("i_eclipse_limit")),
                "z": _num(r.get("z")),
                "source": r.get("source") or "",
                "label": r.get("label") or "",
                "pair": r.get("pair") or "",
            }
            for r in rows
        ]
        chi2 = sum((x["z"] or 0.0) ** 2 for x in inc if x["z"] is not None)
        if not any(x["z"] is not None for x in inc):  # older tables (no z column)
            chi2 = sum(
                ((x["i_obs"] - x["i_model"]) / x["sigma"]) ** 2
                for x in inc
                if None not in (x["i_obs"], x["i_model"], x["sigma"])
            )
        facts.append(
            {
                "k": "inclinations: χ²",
                "v": fmt(_f(chi2), 3),
                "note": f"{len(inc)} epochs",
                "tone": "",
            }
        )
        data["inclination"] = inc
        data["inc_pairs"] = sorted({x["pair"] for x in inc if x["pair"]})
    # predictions
    pr = folder / "stage2_predictions.csv"
    if pr.exists() and full:
        rows = _read_rows(pr)
        data["predictions"] = rows[:200]
        data["n_predictions"] = len(rows)
        note = _comments(pr)
        if note:
            data["predictions_note"] = " ".join(note)
    # stage 1
    if s1.exists():
        try:
            data["stage1"] = json.loads(s1.read_text(encoding="utf-8"))
        except ValueError:
            pass
        p1 = folder / "stage1_parameters.csv"
        if p1.exists():
            data["stage1_params"] = _read_rows(p1)
    if anchor.exists():
        data["anchor"] = anchor.read_text(encoding="utf-8", errors="replace")
    figs = sorted(folder.glob("*.png"))
    data["figures"] = [str(p.relative_to(root)) if root else p.name for p in figs]
    mt = re.search(r"t = 0 of (?:the|this) configuration[:=]? (?:= )?BJD ([\d.]+)", text) or (
        re.search(
            r"t = 0 of this configuration = BJD ([\d.]+)",
            (folder / "stage2_fitted.yaml").read_text(encoding="utf-8"),
        )
        if (folder / "stage2_fitted.yaml").exists()
        else None
    )
    data["bjd_t0"] = float(mt.group(1).rstrip(".")) if mt else None
    try:
        fy = yaml.safe_load((folder / "stage2_fitted.yaml").read_text(encoding="utf-8")) or {}
        data["periods"] = {
            k: float(v["P"])
            for k, v in (fy.get("orbits") or {}).items()
            if isinstance(v, dict) and v.get("P")
        }
    except (OSError, yaml.YAMLError, TypeError, ValueError):
        data["periods"] = {}
    # the binary each star belongs to (radial velocities are phased with it)
    data["star_orbit"] = {}
    try:
        from ..config import Config
        from ..hierarchy import skeleton

        fcfg = Config.model_validate(fy)
        for node in skeleton(fcfg.tree).nodes():
            for c in (node.left, node.right):
                if isinstance(c, int):
                    data["star_orbit"][f"m{c + 1}"] = node.name
    except Exception:  # noqa: BLE001, S110 - an old or partial file: phase with the first orbit
        pass
    if (folder / "stage2_fitted.yaml").exists():
        data["fitted_config"] = (
            str((folder / "stage2_fitted.yaml").relative_to(root)) if root else "stage2_fitted.yaml"
        )
    if (folder / "report.html").exists():
        data["html"] = str((folder / "report.html").relative_to(root)) if root else "report.html"
    if full:
        data["siblings"] = _sibling_fits(folder)
    for w in warnings[:3]:
        quality.append(
            {"k": "warning", "v": w if len(w) <= 90 else w[:88].rstrip() + " …", "tone": "warn"}
        )
    out = {
        "setup": setup if stage == "stage 2" else f"swingersim fit · {stage}",
        "facts": facts,
        "quality": quality,
        "spark": None,
    }
    if oc:
        b = min(oc)
        y = [v for v, st in zip(oc[b]["res"], oc[b]["status"], strict=True) if st == "fitted"]
        out["spark"] = {
            "type": "dots",
            "label": f"measured − model of {b} [min]",
            "y": y[:: max(1, len(y) // 200)],
        }
    if full:
        out["data"] = data
    return out


TEXT_LIMIT_FIT = 60000
_RMS_LINE = re.compile(r"^\s+(\S+)\s+(\d+) minima\s+([-\d.e+]+) ->\s+([-\d.e+]+)\s*$")


def _sibling_fits(folder: Path) -> list[dict]:
    """The other fits of the same system (stage2_report.txt in the folders next to this one):
    reduced chi^2, rms of the minima per pair, the data used - a quick comparison of fit variants."""
    out = []
    for rep in sorted(folder.parent.glob("*/stage2_report.txt")):
        t = rep.read_text(encoding="utf-8", errors="replace")
        chi = _CHI.search(t)
        m = _COUNTS.search(t)
        rms = {}
        if "rms of the residuals [min], start -> fitted:" in t:
            for ln in t.split("rms of the residuals [min], start -> fitted:")[1].splitlines()[1:]:
                mm = _RMS_LINE.match(ln)
                if not mm:
                    break
                rms[mm.group(1)] = _f(float(mm.group(4)))
        out.append(
            {
                "name": rep.parent.name,
                "this": rep.parent == folder,
                "chi2_red": _f(float(chi.group(3))) if chi else None,
                "data": (m.group(1) + " minima" + m.group(2)) if m else "",
                "free": int(m.group(3)) if m else None,
                "rms": rms,
                "date": rep.stat().st_mtime,
                "params": _param_table(rep.parent / "stage2_params.csv"),
                "bjd_t0": _f(m0.group(1))
                if (m0 := re.search(r"t = 0 of the configuration: BJD ([\d.]+)", t))
                else None,
            }
        )
    return out


def _param_table(path: Path) -> dict:
    """parameter -> [value, error] of a stage2_params.csv (empty if missing)."""
    if not path.is_file():
        return {}
    out = {}
    for r in _read_rows(path):
        name = r.get("parameter")
        if name:
            out[name] = [_f(r.get("value")), _f(r.get("error"))]
    return out


# --------------------------------------------------------------------------- measured minima
def measured_minima_for(root: Path, cfg_rel: str | None) -> Path | None:
    """The measured minima that belong to a configuration: its own fit.minima, or that of a
    configuration of the same system in the same folder (<system>_published.yaml ->
    <system>_fit*.yaml)."""
    if not root or not cfg_rel:
        return None
    p = root / cfg_rel
    cands = [p]
    stem = p.stem
    for suf in ("_published", "_anchored", "_fit", "_ensemble"):
        if suf in stem:
            stem = stem.split(suf)[0]
            break
    cands += sorted(p.parent.glob(f"{stem}_*.yaml"))
    for c in cands:
        try:
            d = yaml.safe_load(c.read_text(encoding="utf-8")) or {}
        except (OSError, yaml.YAMLError):
            continue
        f = ((d.get("fit") or {}).get("minima")) if isinstance(d, dict) else None
        if isinstance(f, str) and f:
            q = (c.parent / f).resolve()
            if q.exists():
                return q
    return None


def compare_minima(
    timings_csv: Path, minima_file: Path, root: Path | None = None, table: dict | None = None
) -> dict | None:
    """Measured minima against the model eclipses of a timings output (needs the column bjd_obs,
    i.e. the job was run with bjd_t0): nearest model eclipse of the same binary and kind.
    ``table``: the timings CSV already read (_read_csv)."""
    d = table if table is not None else _read_csv(timings_csv)
    if "bjd_obs" not in d:
        return None
    from ..ocfit import read_minima

    mins = read_minima(minima_file)
    # the model's time axis: years from t = 0 (t_obs_yr); BJD of t = 0 from the two columns
    bjd0 = float(d["bjd_obs"][0]) - float(d["t_obs_yr"][0]) * 365.25
    try:
        rel = str(minima_file.relative_to(root)) if root else minima_file.name
    except ValueError:
        rel = minima_file.name
    out = {"file": rel, "bjd_t0": bjd0, "binaries": {}}
    binary = np.asarray(d["binary"]).astype(str)
    kinds = np.asarray(d["kind"]).astype(str)
    bjd_obs, oc_d = _numbers(d, "bjd_obs"), _numbers(d, "o_c_d")
    for b in sorted(set(binary.tolist())):
        series = {}
        for kind, k1 in (("primary", 0), ("secondary", 1)):
            ii = np.flatnonzero((binary == b) & (kinds == kind))
            if not len(ii):
                continue
            tb = bjd_obs[ii]
            ocm = oc_d[ii]
            o = np.argsort(tb)
            tb, ocm = tb[o], ocm[o]
            sel = (mins.pair == b) & (mins.kind == k1)
            if not np.any(sel):
                continue
            tm = mins.time[sel]
            if len(tb) > 1:
                j = np.clip(np.searchsorted(tb, tm), 1, len(tb) - 1)
                near = np.where(np.abs(tb[j - 1] - tm) < np.abs(tb[j] - tm), j - 1, j)
            else:
                near = np.zeros(len(tm), dtype=int)
            half = 0.5 * float(np.median(np.diff(tb))) if len(tb) > 1 else np.inf
            ok = np.abs(tb[near] - tm) <= half
            res = (tm - tb[near]) * 1440.0
            series[kind] = {
                "t": _arr((tm[ok] - bjd0) / 365.25, digits=9),
                "oc": _arr((ocm[near][ok] * 1440.0 + res[ok]), digits=6),
                "res": _arr(res[ok], digits=6),
                "sigma": _arr(mins.sigma[sel][ok] * 1440.0, digits=4),
                "n_unmatched": int(np.sum(~ok)),
                "rms": _f(float(np.sqrt(np.mean(res[ok] ** 2)))) if np.any(ok) else None,
                "mean": _f(float(np.mean(res[ok]))) if np.any(ok) else None,
            }
        if series:
            out["binaries"][b] = series
    return out


__all__ = [
    "ROUTE",
    "benchmarks_digest",
    "compare_minima",
    "config_name",
    "ensemble_digest",
    "fit_digest",
    "fmt",
    "fmt_range",
    "lightcurve_digest",
    "measured_minima_for",
    "member_draft",
    "minima_digest",
    "run_digest",
    "timings_digest",
]
