"""Monte Carlo ensembles: many runs of one system with uncertain or unknown parameters.

The configuration file describes the nominal system; its ``ensemble`` section says which
parameters are drawn from which distribution (``vary``), how many members to run, and in
which direction(s) of time. Every member is run with the backend of the configuration and
reduced to one row of numbers (fate, extremes, mutual inclinations, eclipses, apsidal periods,
stability verdicts) in ``summary.csv`` (and ``summary.parquet`` if pandas + pyarrow are
installed). Rows are appended to ``rows.jsonl`` as members finish, so an interrupted ensemble
can be resumed.

The same machinery reproduces published population studies, e.g. Pejcha et al. (2013),
MNRAS 435, 943, Table 1 (see swingersim.validation).
"""

from __future__ import annotations

import copy
import csv
import hashlib
import json
import math
import os
import time
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml
from pydantic import ValidationError

from .config import Config, Distribution, EnsembleSpec
from .parallel import run_tasks

FATES = (
    "intact",
    "hierarchy_change",
    "ejection",
    "breakup",
    "exchange",
    "reconfiguration",
    "capture",
    "collision",
    "merger",
    "roche_overflow",
    "energy_error",
    "wall_time",
    "e_reached",
    "gw_unresolved",
    "error",
)


def resolve_processes(n: int | None) -> int:
    """Worker processes of a parallel command: ``n`` if positive; 0 (or less) = the cores of
    the job (environment variable SWINGERSIM_CORES, set by swingersim worker) or, outside the worker,
    all cores of the computer."""
    if n is not None and int(n) > 0:
        return int(n)
    try:
        k = int(os.environ.get("SWINGERSIM_CORES", "0"))
    except ValueError:
        k = 0
    return k if k > 0 else (os.cpu_count() or 1)


# --------------------------------------------------------------------------- drawing members


def draw_value(dist, rng: np.random.Generator) -> float:
    """One value of a distribution of the ``vary`` section."""
    if isinstance(dist, int | float):
        return float(dist)
    if dist == "isotropic":
        return math.degrees(math.acos(rng.uniform(-1.0, 1.0)))
    if dist == "angle":
        return float(rng.uniform(0.0, 360.0))
    d: Distribution = dist
    for _ in range(10000):
        if d.uniform is not None:
            x = rng.uniform(*d.uniform)
        elif d.loguniform is not None:
            lo, hi = np.log(d.loguniform)
            x = math.exp(rng.uniform(lo, hi))
        elif d.normal is not None:
            x = rng.normal(*d.normal)
        elif d.isotropic is not None:
            c1, c2 = (math.cos(math.radians(v)) for v in d.isotropic)
            x = math.degrees(math.acos(rng.uniform(min(c1, c2), max(c1, c2))))
        else:
            x = d.choice[rng.integers(len(d.choice))]
        if (d.min is None or x >= d.min) and (d.max is None or x <= d.max):
            if d.mirror and rng.uniform() < 0.5:
                x = 180.0 - x
            return float(x)
    raise ValueError(f"cannot draw a value of {d} within min/max")


def draw_values(dist, rng: np.random.Generator, n: int) -> np.ndarray:
    """n values of a distribution of the ``vary`` section at once (vectorised `draw_value`;
    the same distribution, not the same random sequence)."""
    if isinstance(dist, int | float):
        return np.full(n, float(dist))
    if dist == "isotropic":
        return np.degrees(np.arccos(rng.uniform(-1.0, 1.0, n)))
    if dist == "angle":
        return rng.uniform(0.0, 360.0, n)
    d: Distribution = dist
    out = np.empty(n)
    todo = np.arange(n)
    for _ in range(10000):
        m = len(todo)
        if d.uniform is not None:
            x = rng.uniform(*d.uniform, m)
        elif d.loguniform is not None:
            lo, hi = np.log(d.loguniform)
            x = np.exp(rng.uniform(lo, hi, m))
        elif d.normal is not None:
            x = rng.normal(*d.normal, m)
        elif d.isotropic is not None:
            c1, c2 = (math.cos(math.radians(v)) for v in d.isotropic)
            x = np.degrees(np.arccos(rng.uniform(min(c1, c2), max(c1, c2), m)))
        else:
            x = np.asarray(d.choice, float)[rng.integers(len(d.choice), size=m)]
        ok = np.ones(m, bool)
        if d.min is not None:
            ok &= x >= d.min
        if d.max is not None:
            ok &= x <= d.max
        out[todo[ok]] = x[ok]
        todo = todo[~ok]
        if not len(todo):
            break
    else:
        raise ValueError(f"cannot draw a value of {d} within min/max")
    if d.mirror:
        flip = rng.uniform(size=n) < 0.5
        out[flip] = 180.0 - out[flip]
    return out


def _set(data: dict, path: str, value: float) -> None:
    parts = path.split(".")
    if parts[0] == "masses":
        data["masses"][parts[1]] = value
        return
    group = data.setdefault(parts[0], {}).setdefault(parts[1], {})
    group[parts[2]] = value
    # an orbit is given by exactly one of P and a
    if parts[0] == "orbits" and parts[2] in ("P", "a"):
        group["a" if parts[2] == "P" else "P"] = None


@dataclass
class Member:
    index: int
    config: dict  # validated configuration (model_dump) without the ensemble section
    params: dict[str, float]
    rejected: int  # draws rejected before this one


def _criteria_ok(cfg: Config, require: list[str]) -> bool:
    if not require:
        return True
    from .hierarchy import build_tree, initialize
    from .stability import nested_triples

    root, masses = build_tree(cfg)
    pos, vel = initialize(root, masses)
    if "hierarchical" in require and not _hierarchical(root):
        return False
    if not ({"stable_ma01", "stable_v22"} & set(require)):
        return True
    trip = nested_triples(root, masses, pos, vel)
    if "stable_ma01" in require and min(t.ma01 for t in trip) <= 1.0:
        return False
    return not ("stable_v22" in require and min(t.v22 for t in trip) <= 1.0)


def _hierarchical(root) -> bool:
    """Every inner orbit's apocentre lies inside the pericentre of the orbit around it."""
    from .hierarchy import Node

    for n in root.nodes():
        el = n.elements
        for c in (n.left, n.right):
            if isinstance(c, Node) and c.elements.a * (1 + c.elements.e) >= el.a * (1 - el.e):
                return False
    return True


def _ma01_min(c: Config) -> float:
    from .hierarchy import build_tree, initialize
    from .stability import nested_triples

    root, masses = build_tree(c)
    pos, vel = initialize(root, masses)
    return min(t.ma01 for t in nested_triples(root, masses, pos, vel))


def draw_members(cfg: Config) -> list[Member]:
    """Draw the members of the ensemble (deterministic for a given seed)."""
    ens = cfg.ensemble
    if ens is None:
        raise ValueError("the configuration has no 'ensemble' section")
    rng = np.random.default_rng(ens.seed)
    base = cfg.model_dump(exclude={"ensemble"})
    members = []
    for k in range(ens.n):
        last_error: Exception | None = None
        for tries in range(ens.max_tries):
            data = copy.deepcopy(base)
            params = {p: draw_value(d, rng) for p, d in ens.vary.items()}
            for p, v in params.items():
                _set(data, p, v)
            enc = (data.get("physics") or {}).get("encounters")
            if enc is not None:  # every member meets its own sequence of passing stars
                enc["seed"] = int(base["physics"]["encounters"].get("seed", 1)) + k
            try:
                c = Config.model_validate(data)
                c._base_dir = cfg.base_dir
                from .hierarchy import build_tree, initialize

                root, masses = build_tree(c)
                initialize(root, masses)  # unbound / invalid orbits raise here
                if not _criteria_ok(c, ens.require):
                    continue
                win = getattr(ens, "ma01_window", None)
                if ens.ma01_unstable_keep < 1.0 or win is not None:
                    r = _ma01_min(c)
                    if win is not None and not (win[0] <= r <= win[1]):
                        continue
                    if (
                        ens.ma01_unstable_keep < 1.0
                        and r <= 1.0
                        and rng.uniform() >= ens.ma01_unstable_keep
                    ):
                        continue
            except (ValidationError, ValueError) as exc:
                last_error = exc
                continue
            members.append(Member(k, c.model_dump(), params, tries))
            break
        else:
            why = (
                f"; the last draw was invalid: {type(last_error).__name__}: {last_error}"
                if last_error is not None
                else ""
            )
            raise RuntimeError(
                f"member {k}: no acceptable draw in {ens.max_tries} tries; "
                f"check the distributions and 'require'{why}"
            ) from last_error
    return members


# --------------------------------------------------------------------------- one member


def _fate(res) -> str:
    """Fate of a member: the event that ended its run, else its first hierarchy change (a
    merger followed by physics.mergers counts as 'merger')."""
    if res.event in (
        "collision",
        "merger",
        "roche_overflow",
        "energy_error",
        "wall_time",
        "e_reached",
        "gw_unresolved",  # GW decay faster than an output step: the run stopped
    ):
        return res.event
    if res.transitions:
        return res.transitions[0]["kind"]
    if res.event == "hierarchy_change":
        return "hierarchy_change"
    return "intact"


def _parent_period(cfg: Config, res, name: str) -> float:
    """Period [yr] at t = 0 of the orbit that contains binary ``name`` (its perturber)."""
    from . import units

    tree = cfg.tree
    parent = next((n for n, kids in tree.items() if name in kids), None)
    if parent is None:
        return 0.0

    def mass(x):
        return cfg.masses[x] if x in cfg.masses else sum(mass(c) for c in tree[x])

    return units.sma_to_period(float(res.orbits[parent]["a"][0]), mass(parent))


def _observables(
    res, name: str, stars: list[str], p_parent: float, window_yr: float | None = None
) -> dict:
    """Quantities an observer of the binary measures now: the inclination to the sky and its
    rate of change at t = 0, and when the binary starts or stops eclipsing.

    ``tilt`` = |i - 90 deg| is the distance from edge-on; its rate is what changing eclipse
    depths measure (negative: the binary turns towards edge-on, the eclipses deepen), whatever
    the unknown sense of revolution. The rates are linear fits over two periods of the orbit
    that contains the binary (at least 6 outputs), so that the oscillations of the osculating
    N-body elements with that period average out; ``window_yr`` (run.rate_window_yr) sets the
    length of the fit instead. The fit is one-sided (forward or backward from t = 0); the
    conditioning on observations combines both directions into a central rate
    (`observation_weights`)."""
    from . import units
    from .diagnostics.binary import eclipsing

    o = res.orbits[name]
    t = np.asarray(res.times)
    inc = np.degrees(np.asarray(o["inc"]))
    if window_yr:
        k = max(min(len(t), 3), int(np.searchsorted(np.abs(t), window_yr, side="right")))
    else:
        k = max(min(len(t), 6), int(np.searchsorted(np.abs(t), 2.0 * p_parent)) + 1)
    ok = np.isfinite(inc[:k])
    rate = float(np.polyfit(t[:k][ok], inc[:k][ok], 1)[0]) if ok.sum() >= 2 else math.nan
    tilt = np.abs(inc - 90.0)
    rate_edge = float(np.polyfit(t[:k][ok], tilt[:k][ok], 1)[0]) if ok.sum() >= 2 else math.nan
    if res.stars:
        r_sum = (
            sum(
                np.asarray(
                    (res.stellar or {}).get(s, {}).get("radius_rsun", res.stars[s].radius_rsun)
                )
                for s in stars
            )
            * units.RSUN_AU
        )
    else:
        r_sum = 0.0
    ecl = eclipsing(o, r_sum)
    change = np.flatnonzero(ecl != ecl[0])
    return {
        f"{name}_i0": float(inc[0]),
        f"{name}_didt0_deg_yr": rate,
        f"{name}_dtiltdt0_deg_yr": rate_edge,
        f"{name}_eclipsing0": float(ecl[0]),
        f"{name}_eclipse_change_yr": float(t[change[0]]) if len(change) else math.nan,
    }


def summarize(cfg: Config, res) -> dict:
    """One row of numbers describing a finished run."""
    from .hierarchy import build_tree, initialize
    from .stability import nested_triples

    row: dict = {
        "fate": _fate(res),
        "event": res.event or "",
        "event_time_yr": res.event_time if res.event_time is not None else math.nan,
        "t_end_yr": float(res.times[-1]),
        "n_transitions": len(res.transitions),
        "first_transition_yr": res.transitions[0]["time_yr"] if res.transitions else math.nan,
        "final_hierarchy": res.hierarchy[-1] if res.hierarchy else "",
        "max_energy_error": float(np.nanmax(res.energy_error)),
    }
    for name, o in res.orbits.items():
        e_out = np.asarray(o["e"])
        ex = res.extremes.get(name, {})
        e_max = np.nanmax([np.nanmax(e_out), ex.get("e_max", math.nan)])
        row[f"{name}_e0"] = float(e_out[0])
        row[f"{name}_e_min"] = float(np.nanmin(e_out))
        row[f"{name}_e_max"] = float(e_max)
        row[f"{name}_log1me_min"] = float(np.log10(max(1.0 - e_max, 1e-16)))
        row[f"{name}_t_e_max_yr"] = ex.get("t_e_max_yr", math.nan)
        row[f"{name}_q_min_au"] = ex.get("q_min_au", math.nan)
        for key, val in ex.items():
            if key.startswith("t_1me_lt_"):
                row[f"{name}_{key}"] = val
        row[f"{name}_a_min"] = float(np.nanmin(o["a"]))
        row[f"{name}_a_max"] = float(np.nanmax(o["a"]))
        if name not in res.binaries:  # apsidal period of an outer orbit (binaries: below)
            from .diagnostics.binary import measured_apsidal_period

            row[f"{name}_U_measured_yr"] = measured_apsidal_period(
                res.times, o["Omega"], o["omega"]
            )
    for pair, s in res.mutual.items():
        d = np.degrees(np.asarray(s))
        key = pair.replace("-", "_")
        row[f"imut_{key}_0"] = float(d[0])
        row[f"imut_{key}_min"] = float(np.nanmin(d))
        row[f"imut_{key}_max"] = float(np.nanmax(d))
    for name, b in res.binaries.items():
        row[f"{name}_U_measured_yr"] = b.get("U_measured_yr", math.nan)
        row[f"{name}_eclipse_fraction"] = b.get("eclipse_fraction", math.nan)
        row.update(
            _observables(
                res,
                name,
                b["stars"].split("+"),
                _parent_period(cfg, res, name),
                cfg.run.rate_window_yr,
            )
        )
    for key in ("eclipse_fraction_both", "eclipse_fraction_all"):
        if key in res.system:
            row[key] = res.system[key]
    for star, r in res.roche.items():
        row[f"{star}_roche_fill_max"] = r["fill_max"]
    if res.megno is not None:
        row["megno"] = float(res.megno[-1])
    root, masses = build_tree(cfg)
    # the star that filled its Roche lobe and whether its companion is a star (binary) or a
    # subsystem (a tertiary overflowing onto an inner binary: triple common envelope)
    stars = (res.event_detail or {}).get("stars") if res.event == "roche_overflow" else None
    if stars:
        k = int(stars[0][1:]) - 1
        row["roche_first"] = stars[0]
        for n in root.nodes():
            if k in (n.left, n.right):
                other = n.right if n.left == k else n.left
                row["roche_first_kind"] = "binary" if isinstance(other, int) else "tertiary"
    for key, v in res.system.items():  # mass transfer and common envelope summaries
        if key.startswith("mass_transfer_") and (v is None or isinstance(v, (int, float, str))):
            row[key] = v
    pos, vel = initialize(root, masses)
    trip = nested_triples(root, masses, pos, vel)
    row["ma01_min"] = min(t.ma01 for t in trip)
    row["v22_min"] = min(t.v22 for t in trip)
    return row


def run_member(task: dict) -> dict:
    """Run one member in one direction (top-level function so that it can be pickled)."""
    from .backends.nbody import run_nbody
    from .backends.secular import run_secular

    t0 = time.perf_counter()
    row = {"member": task["index"], "direction": task["direction"], **task["params"]}
    try:
        cfg = Config.model_validate(task["config"])
        cfg._base_dir = task["base_dir"]
        sign = 1.0 if task["direction"] == "forward" else -1.0
        cfg.run.t_max_yr = sign * abs(cfg.run.t_max_yr)
        res = run_secular(cfg) if cfg.run.backend == "secular" else run_nbody(cfg)
        row.update(summarize(cfg, res))
        if task["mlp"]:
            from .stability import mlp_stable

            v = mlp_stable(cfg)
            row["mlp_stable"] = math.nan if v is None else float(v)
        if task["ghost"]:
            from .stability import ghost_test

            g = ghost_test(
                cfg, n_outer=task["ghost_orbits"], max_wall_s=task.get("ghost_max_wall_s")
            )
            row["ghost_stable"] = math.nan if g.reason == "timeout" else float(g.stable)
            row["ghost_reason"] = g.reason
            row["ghost_delta_max"] = g.delta_max
        if task["ghost"] or task["mlp"]:
            from .stability import mlp_parameters

            try:  # the 11 parameters of Vynatheya et al. (2023), for stability classifiers
                row.update({f"v23_{k}": float(v) for k, v in mlp_parameters(cfg).items()})
            except ValueError:
                pass
        if task["save_dir"]:
            from .io import save_run

            save_run(
                Path(task["save_dir"]) / f"member_{task['index']:05d}_{task['direction']}.h5",
                cfg,
                res,
            )
    except Exception as exc:  # noqa: BLE001 - one failing member must not stop the ensemble
        row.update({"fate": "error", "error": f"{type(exc).__name__}: {exc}"})
    row["wall_s"] = time.perf_counter() - t0
    return row


# --------------------------------------------------------------------------- the ensemble


def _fingerprint(cfg: Config) -> str:
    """Identity of an ensemble. ``n`` and ``observed`` are left out: an ensemble can be extended
    by more members (``--resume`` with a larger n; members depend only on the seed), and the
    observations only weight the finished members."""
    data = cfg.model_dump()
    data["ensemble"] = {k: v for k, v in data["ensemble"].items() if k not in ("n", "observed")}
    return hashlib.sha256(json.dumps(data, sort_keys=True, default=str).encode()).hexdigest()[:16]


def _config_differences(old: dict, new: dict, prefix: str = "") -> list[str]:
    """Settings that differ between a stored and the current configuration. Keys that the
    stored one does not have (added in a later SWINGERsim version) and ensemble.n / observed
    are ignored."""
    out = []
    for key, value in old.items():
        path = f"{prefix}{key}"
        if path in ("ensemble.n", "ensemble.observed"):
            continue
        if key not in new:
            out.append(path)
        elif isinstance(value, dict) and isinstance(new[key], dict):
            out += _config_differences(value, new[key], path + ".")
        elif value != new[key]:
            out.append(path)
    return out


def _json_safe(x):
    """NaN and infinities as null (bare NaN is not JSON; `_read_rows` reads both)."""
    if isinstance(x, float) and not math.isfinite(x):
        return None
    if isinstance(x, dict):
        return {k: _json_safe(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_json_safe(v) for v in x]
    if isinstance(x, np.floating):
        return _json_safe(float(x))
    if isinstance(x, np.integer):
        return int(x)
    return x


def row_line(row: dict) -> str:
    """One line of rows.jsonl (also of grid.jsonl and the scattering rows)."""
    return json.dumps(_json_safe(row), allow_nan=False) + "\n"


def read_jsonl(path: Path, null_as_nan: bool = True) -> tuple[list[dict], int]:
    """Rows of a JSON-lines file and the number of malformed lines that were skipped (a line
    torn by a crash or a full disk while it was written); top-level null values are read as
    NaN (``null_as_nan``; rows written by older versions (quadsim) have bare NaN, read as well)."""
    if not path.exists():
        return [], 0
    rows, bad = [], 0
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        try:
            r = json.loads(line)
        except ValueError:
            bad += 1
            continue
        if not isinstance(r, dict):
            bad += 1
            continue
        if null_as_nan:
            r = {k: (math.nan if v is None else v) for k, v in r.items()}
        rows.append(r)
    if bad:
        warnings.warn(
            f"{path}: {bad} malformed line(s) skipped (torn by an interrupted write)", stacklevel=2
        )
    return rows, bad


def rewrite_jsonl(path: Path, rows: list[dict]) -> None:
    """Replace a JSON-lines file atomically (resume after a torn last line, retry_errors)."""
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text("".join(row_line(r) for r in rows), encoding="utf-8")
    tmp.replace(path)


def _dedupe(rows: list[dict]) -> list[dict]:
    """One row per (member, direction): the last one written (a member run again)."""
    out: dict = {}
    for r in rows:
        out[(r["member"], r["direction"])] = r
    return list(out.values())


def _read_rows(path: Path) -> list[dict]:
    return _dedupe(read_jsonl(path)[0])


def _lost_row(task: dict, message: str) -> dict:
    return {
        "member": task["index"],
        "direction": task["direction"],
        **task["params"],
        "fate": "error",
        "error": message,
        "wall_s": math.nan,
    }


def run_ensemble(
    cfg: Config,
    out_dir: str | Path,
    processes: int = 1,
    resume: bool = False,
    progress=None,
    retry_errors: bool = False,
) -> list[dict]:
    """Run (or resume) the ensemble and write rows.jsonl, summary.csv (and .parquet).
    With ``retry_errors`` (needs ``resume``) the members that ended with an error are removed
    from rows.jsonl and run again (e.g. after a bug fix). A member whose worker process dies
    (segmentation fault, out-of-memory kill) is recorded with fate 'error'."""
    import platform
    from datetime import UTC, datetime

    from . import __version__

    ens: EnsembleSpec = cfg.ensemble
    processes = resolve_processes(processes)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows_file = out / "rows.jsonl"
    meta_file = out / "ensemble.yaml"
    fp = _fingerprint(cfg)
    if rows_file.exists() and not resume:
        raise FileExistsError(f"{rows_file} exists; use resume=True (--resume) or a new directory")
    meta = {}
    if resume and meta_file.exists():
        meta = yaml.safe_load(meta_file.read_text()) or {}
        if meta.get("fingerprint") != fp:
            current = json.loads(json.dumps(cfg.model_dump(), default=str))  # tuples -> lists
            diff = _config_differences(meta.get("config") or {}, current)
            if diff:
                raise ValueError(
                    "the configuration changed since the ensemble was started: " + ", ".join(diff)
                )
    now = datetime.now(UTC).isoformat(timespec="seconds")
    meta_file.write_text(
        yaml.safe_dump(
            {
                "fingerprint": fp,
                "base_dir": cfg.base_dir,
                "swingersim": __version__,
                "python": platform.python_version(),
                "started": meta.get("started", now),
                **({"resumed": now} if meta.get("started") else {}),
                "config": cfg.model_dump(),
            },
            sort_keys=False,
        )
    )
    raw, bad = read_jsonl(rows_file)
    old_rows = _dedupe(raw)
    if retry_errors and resume:
        old_rows = [r for r in old_rows if r.get("fate") != "error"]
    if bad or len(old_rows) < len(raw):  # torn lines, duplicates or retried errors
        rewrite_jsonl(rows_file, old_rows)
    done = {(r["member"], r["direction"]) for r in old_rows}
    members = draw_members(cfg)
    save_dir = str(out / "runs") if ens.save_runs else ""
    tasks = [
        {
            "index": m.index,
            "direction": d,
            "config": m.config,
            "params": m.params,
            "base_dir": cfg.base_dir,
            "mlp": ens.mlp,
            "ghost": ens.ghost_test,
            "ghost_max_wall_s": ens.ghost_max_wall_s,
            "ghost_orbits": ens.ghost_orbits,
            "save_dir": save_dir,
        }
        for m in members
        for d in ens.directions
        if (m.index, d) not in done
    ]
    total = len(members) * len(ens.directions)
    n_done = sum(1 for r in old_rows if r["member"] < ens.n)
    with open(rows_file, "a", encoding="utf-8") as f:
        for row in run_tasks(run_member, tasks, processes, _lost_row):
            f.write(row_line(row))
            f.flush()
            n_done += 1
            if progress:
                progress(n_done, total, row)
    rows = _members_of(_read_rows(rows_file), ens.n, ens.directions)
    write_summary(rows, out)
    return rows


def _members_of(rows: list[dict], n: int | None, directions=None) -> list[dict]:
    """Rows of members 0 .. n-1 (and of the given directions), sorted: rows left from a run
    with a larger n (resumed with a smaller one) stay in rows.jsonl but are not part of the
    ensemble."""
    keep = [
        r
        for r in rows
        if (n is None or r["member"] < n) and (directions is None or r["direction"] in directions)
    ]
    return sorted(keep, key=lambda r: (r["member"], r["direction"]))


def write_summary(rows: list[dict], out: Path) -> None:
    columns: list[str] = []
    for r in rows:
        for k in r:
            if k not in columns:
                columns.append(k)
    with open(out / "summary.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=columns, restval="")
        w.writeheader()
        w.writerows(rows)
    parquet = out / "summary.parquet"
    try:
        import pandas as pd

        pd.DataFrame(rows, columns=columns).to_parquet(parquet, index=False)
    except Exception as exc:  # noqa: BLE001 - Parquet is optional; summary.csv is always written
        parquet.unlink(missing_ok=True)  # never leave a stale Parquet next to a new CSV
        if not isinstance(exc, ImportError):
            warnings.warn(f"no summary.parquet ({type(exc).__name__}: {exc})", stacklevel=2)


def load_rows(out_dir: str | Path) -> list[dict]:
    """Rows of an ensemble folder, one per member and direction (the members 0 .. n-1 of the
    stored configuration when ensemble.yaml is there)."""
    out = Path(out_dir)
    n = directions = None
    meta = out / "ensemble.yaml"
    if meta.exists():
        try:
            ens = ((yaml.safe_load(meta.read_text()) or {}).get("config") or {}).get("ensemble")
            n = int(ens["n"])
            directions = list(ens.get("directions") or ["forward"])
        except (ValueError, TypeError, KeyError, AttributeError, yaml.YAMLError):
            n = directions = None
    return _members_of(_read_rows(out / "rows.jsonl"), n, directions)


def check_observed(cfg: Config, observed: dict | None = None) -> None:
    """Cheap checks of the observed columns before an ensemble runs for hours: a valid
    direction prefix, and drawn parameters (columns with a dot) that are varied."""
    observed = cfg.ensemble.observed if observed is None else observed
    for col in observed:
        _, c = split_column(col)
        if "." in c and c not in cfg.ensemble.vary:
            raise ValueError(f"observed column '{col}': '{c}' is not an ensemble.vary parameter")


# --------------------------------------------------------------------------- statistics


def wilson(k: float, n: float, z: float = 1.0) -> tuple[float, float, float]:
    """Fraction k/n with its Wilson score interval (z = 1: 68 %). With weighted members, k is
    the weighted count scaled to the effective number of members n."""
    if n <= 0:
        return math.nan, math.nan, math.nan
    p = k / n
    den = 1.0 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(max(p * (1 - p), 0.0) / n + z * z / (4 * n * n)) / den
    return p, max(0.0, c - h), min(1.0, c + h)


def _frac(k, n) -> str:
    p, lo, hi = wilson(k, n)
    return f"{100 * p:5.1f} % ({100 * lo:.1f}-{100 * hi:.1f}) [{k}/{n}]"


# --------------------------------------------------------------------------- conditioning


RATE_SUFFIXES = ("_didt0_deg_yr", "_dtiltdt0_deg_yr")


def is_rate_column(col: str) -> bool:
    """Present rates of change, measured by an observer over a baseline centred on now."""
    return col.endswith(RATE_SUFFIXES)


def split_column(col: str) -> tuple[str, str]:
    """'backward:B_e0' -> ('backward', 'B_e0'); without a prefix: 'central' for the present
    rates of change (mean of the forward and backward one-sided fits), else 'forward'."""
    if ":" in col:
        d, c = col.split(":", 1)
        if d not in ("forward", "backward", "central"):
            raise ValueError(f"unknown direction in '{col}' (forward:, backward: or central:)")
        return d, c
    return ("central" if is_rate_column(col) else "forward"), col


def member_value(runs: dict[str, dict[int, dict]], m: int, col: str) -> float:
    """Value of a (possibly prefixed) summary column for member m. 'central' is the mean of
    the forward and backward values; with only one direction run it falls back to that one.
    Other columns fall back to the other direction when the requested one is missing (the
    present state is shared by both)."""

    def num(r, c):
        x = r.get(c, math.nan) if r else math.nan
        return math.nan if x is None or x == "" else float(x)

    d, c = split_column(col)
    f, b = runs["forward"].get(m), runs["backward"].get(m)
    if d == "central":
        xs = [x for x in (num(f, c), num(b, c)) if math.isfinite(x)]
        if f and b:
            return 0.5 * (xs[0] + xs[1]) if len(xs) == 2 else math.nan
        return xs[0] if xs else math.nan
    r = runs[d].get(m) or (b if d == "forward" else f)
    return num(r, c)


def index_runs(rows: list[dict]) -> dict[str, dict[int, dict]]:
    """Rows of rows.jsonl by direction and member: runs['forward'][m]."""
    runs: dict[str, dict[int, dict]] = {"forward": {}, "backward": {}}
    for r in rows:
        runs[r["direction"]][r["member"]] = r
    return runs


def passes_limits(runs: dict[str, dict[int, dict]], m: int, observed: dict) -> bool:
    """Does member m pass every hard condition of ``observed`` - {equals: ...} on a text
    column, min / max limits, and a finite value of every limited or Gaussian column unless
    missing_ok? (The Gaussian part {value, sigma} is the likelihood, not a condition.)"""
    for col, ob in observed.items():
        ob = ob if isinstance(ob, dict) else ob.model_dump()
        if ob.get("equals") is not None:  # categorical column, e.g. fate
            allowed = ob["equals"] if isinstance(ob["equals"], list) else [ob["equals"]]
            d, c = split_column(col)
            r = runs[d].get(m) or runs["backward" if d == "forward" else "forward"].get(m)
            if r is None or r.get(c) not in allowed:
                return False
            continue
        x = member_value(runs, m, col)
        if not math.isfinite(x):
            if ob.get("missing_ok"):
                continue
            return False
        if ob.get("min") is not None and x < ob["min"]:
            return False
        if ob.get("max") is not None and x > ob["max"]:
            return False
    return True


def observation_weights(rows: list[dict], observed: dict) -> dict[int, float]:
    """Weight of every member given measured values: the Gaussian likelihood
    w = exp(-chi^2 / 2), chi^2 = sum ((x - value) / sigma)^2, times 0 if a hard limit
    (min / max) is violated.

    ``observed`` maps a summary column (e.g. ``B_dtiltdt0_deg_yr``) to {value, sigma} and/or
    {min, max}, or a text column (``fate``) to {equals: value or list} (dicts or
    ObservedValue). The column is taken from the forward run (the present
    state is shared by both directions) unless it is written ``backward:<column>``. The present
    rates of change (``*_didt0_deg_yr``, ``*_dtiltdt0_deg_yr``) are by default *central*: the
    mean of the forward and backward fits, i.e. a rate over a baseline centred on now, as an
    observer measures it (``forward:`` or ``backward:`` select one side). A member without a
    finite value of an observed column gets weight 0.
    """
    runs = index_runs(rows)
    members = sorted(set(runs["forward"]) | set(runs["backward"]))
    known = {k for r in rows for k in r}
    missing = [c for c in observed if split_column(c)[1] not in known]
    if missing:
        raise ValueError(f"observed columns not in the summary: {', '.join(missing)}")
    obs = {c: (o if isinstance(o, dict) else o.model_dump()) for c, o in observed.items()}
    gauss = {c: o for c, o in obs.items() if o.get("value") is not None}
    weights = {}
    for m in members:
        if not passes_limits(runs, m, obs):
            weights[m] = 0.0
            continue
        chi2 = 0.0
        for c, o in gauss.items():
            x = member_value(runs, m, c)
            if math.isfinite(x):  # a missing value got here only with missing_ok
                chi2 += ((x - o["value"]) / o["sigma"]) ** 2
        weights[m] = math.exp(-0.5 * chi2)
    return weights


def effective_size(w: np.ndarray) -> float:
    """Kish effective number of members, (sum w)^2 / sum w^2."""
    s2 = float(np.sum(w**2))
    return float(np.sum(w) ** 2 / s2) if s2 > 0 else 0.0


def evidence(rows: list[dict], weights: dict[int, float]) -> tuple[float, float, int]:
    """Bayesian evidence of the ensemble: the mean weight Z = (1/n) sum w over the n finished
    members (members removed by hard limits count with w = 0), its Monte Carlo error
    std(w)/sqrt(n) and n. Ratios of Z of two ensembles conditioned on the same observations are
    Bayes factors."""
    w = np.array([weights.get(r["member"], 0.0) for r in rows if r.get("fate") != "error"])
    if not len(w):
        return float("nan"), float("nan"), 0
    return float(w.mean()), float(w.std() / math.sqrt(len(w))), len(w)


def evidence_text(rows: list[dict], weights: dict[int, float]) -> str:
    z, dz, n = evidence(rows, weights)
    if n == 0:
        return "evidence Z: no finished member"
    if z < 1.0 / n and effective_size(np.array([weights.get(r["member"], 0.0) for r in rows])) < 2:
        return (
            f"evidence Z = {z:.2g} (fewer than two matching members: quote Z < 1/n = {1.0 / n:.2g})"
        )
    return f"evidence Z = {z:.3g} +- {dz:.2g}"


def _wquantile(x: np.ndarray, w: np.ndarray, qs) -> np.ndarray:
    ok = np.isfinite(x) & (w > 0)
    if not ok.any():
        return np.full(len(qs), np.nan)
    x, w = x[ok], w[ok]
    order = np.argsort(x)
    x, w = x[order], w[order]
    c = (np.cumsum(w) - 0.5 * w) / np.sum(w)
    return np.interp(qs, c, x)


def _wfrac(mask: np.ndarray, w: np.ndarray) -> str:
    total = float(np.sum(w))
    if total <= 0:
        return "  n/a"
    p = float(np.sum(w[mask]) / total)
    n_eff = effective_size(w)
    _, lo, hi = wilson(p * n_eff, n_eff)
    if np.all((w == 0) | (w == 1)):
        return _frac(int(np.sum(mask & (w > 0))), int(np.sum(w > 0)))
    return f"{100 * p:5.1f} % ({100 * lo:.1f}-{100 * hi:.1f}) [n_eff {n_eff:.1f}]"


def report(
    rows: list[dict],
    thresholds=(1e-1, 1e-2, 1e-3),
    weights: dict[int, float] | None = None,
) -> str:
    """Text summary of an ensemble: fates, high eccentricities, inclinations, eclipses.

    With ``weights`` (from observation_weights) every statistic is weighted: the result is the
    distribution conditioned on the observations, and the drawn parameters are summarised
    before (prior) and after (posterior) the conditioning."""
    lines = []
    for direction in ("forward", "backward"):
        rs = [r for r in rows if r["direction"] == direction]
        if not rs:
            continue
        n = len(rs)
        wall = np.array([r.get("wall_s", np.nan) for r in rs], dtype=float)
        errors = [r for r in rs if r["fate"] == "error"]
        rs_all = [r for r in rs if r["fate"] != "error"]
        if weights is not None:
            rs = [r for r in rs if weights.get(r["member"], 0.0) > 0]
        w_all = np.array([weights.get(r["member"], 0.0) if weights else 1.0 for r in rs])
        head = f"{direction}: {n} members"
        if weights is not None:
            head += (
                f", conditioned on the observations: effective number of members "
                f"{effective_size(w_all):.1f}, {evidence_text(rs_all, weights)}"
            )
        lines.append(head + " (fractions with 68 % Wilson intervals)")
        fate = np.array([r["fate"] for r in rs])
        for f in FATES:
            if np.any(fate == f):
                lines.append(f"  fate {f:16s} {_wfrac(fate == f, w_all)}")
        ok = [r for r in rs if r["fate"] != "error"]
        w = np.array([weights.get(r["member"], 0.0) if weights else 1.0 for r in ok])
        orbits = sorted(
            {k.removesuffix("_e_max") for r in ok for k in r if k.endswith("_e_max")},
            key=lambda x: (
                ("A", "B", "inner", "mid", "out").index(x)
                if x in ("A", "B", "inner", "mid", "out")
                else 9
            ),
        )
        for name in orbits:
            vals = np.array([r[f"{name}_e_max"] for r in ok], dtype=float)
            if not np.isfinite(vals).any():
                continue
            q = _wquantile(vals, w, [0.16, 0.5, 0.84])
            parts = [f"e_max median {q[1]:.3f} (16-84 %: {q[0]:.3f}-{q[2]:.3f})"]
            for t in thresholds:
                frac = np.sum(w[1.0 - vals < t]) / max(np.sum(w), 1e-300)
                parts.append(f"1-e<{t:g}: {100 * frac:.1f} %")
            lines.append(f"  {name:6s} " + ", ".join(parts))
        for key in sorted({k for r in ok for k in r if k.startswith("imut_") and k.endswith("_0")}):
            pair = key.removeprefix("imut_").removesuffix("_0")
            lo = np.array([r[f"imut_{pair}_min"] for r in ok], dtype=float)
            hi = np.array([r[f"imut_{pair}_max"] for r in ok], dtype=float)
            i0 = np.array([r[key] for r in ok], dtype=float)
            q0 = _wquantile(i0, w, [0.16, 0.5, 0.84])
            cross = np.sum(w[(lo < 90) & (hi > 90)]) / max(np.sum(w), 1e-300)
            lines.append(
                f"  i_mut {pair:7s} now {q0[1]:.1f} deg (16-84 %: {q0[0]:.1f}-{q0[2]:.1f}), "
                f"range over the run: median {_wquantile(lo, w, [0.5])[0]:.1f}..."
                f"{_wquantile(hi, w, [0.5])[0]:.1f} deg, members that cross 90 deg: "
                f"{100 * cross:.1f} %"
            )
        for name in [n for n in ("A", "B", "inner") if any(f"{n}_didt0_deg_yr" in r for r in ok)]:
            v = np.array([r[f"{name}_didt0_deg_yr"] for r in ok], dtype=float)
            q = _wquantile(v, w, [0.16, 0.5, 0.84])
            ch = np.array([r.get(f"{name}_eclipse_change_yr", np.nan) for r in ok], dtype=float)
            e0 = np.array([r.get(f"{name}_eclipsing0", np.nan) for r in ok], dtype=float)
            text = f"  {name:6s} di/dt now {q[1]:+.3g} deg/yr (16-84 %: {q[0]:+.3g}..{q[2]:+.3g})"
            ecl_now = (e0 == 1) & (w > 0)
            if ecl_now.any():
                lost = ecl_now & np.isfinite(ch)
                frac = np.sum(w[lost]) / np.sum(w[ecl_now])
                text += (
                    "; of the members eclipsing now, eclipses "
                    + ("end" if direction == "forward" else "began")
                    + f" within the run in {100 * frac:.0f} %"
                )
                if lost.any():
                    med = _wquantile(np.abs(ch[lost]), w[lost], [0.5])[0]
                    text += f" (median after {med:.4g} yr)"
            lines.append(text)
        if any("eclipse_fraction_both" in r for r in ok):
            parts = []
            for key, label in (
                ("A_eclipse_fraction", "A"),
                ("B_eclipse_fraction", "B"),
                ("eclipse_fraction_both", "both"),
            ):
                v = np.array([r.get(key, np.nan) for r in ok], dtype=float)
                good = np.isfinite(v) & (w > 0)
                mean = np.sum(w[good] * v[good]) / np.sum(w[good]) if good.any() else np.nan
                parts.append(f"{label} {100 * mean:.1f} %")
            lines.append("  eclipses (time averaged over members): " + ", ".join(parts))
        for col, label in (
            ("ma01_min", "Mardling-Aarseth"),
            ("v22_min", "Vynatheya+22"),
            ("mlp_stable", "MLP Vynatheya+23"),
            ("ghost_stable", "N-body ghost test"),
        ):
            v = np.array([r.get(col, np.nan) for r in ok], dtype=float)
            good = np.isfinite(v)
            if good.any():
                stable = v > 1.0 if col.endswith("_min") else v > 0.5
                lines.append(f"  stable by {label:18s} {_wfrac(stable[good], w[good])}")
        if weights is not None and direction == "forward":
            params = [k for k in rows[0] if "." in k]
            if params:
                lines.append("  drawn parameters, prior -> conditioned (median, 16-84 %):")
            every = [r for r in rows if r["direction"] == direction]
            for p in params:
                x_all = np.array([r[p] for r in every], dtype=float)
                qa = _wquantile(x_all, np.ones(len(x_all)), [0.16, 0.5, 0.84])
                x = np.array([r[p] for r in ok], dtype=float)
                qc = _wquantile(x, w, [0.16, 0.5, 0.84])
                lines.append(
                    f"    {p:20s} {qa[1]:.4g} ({qa[0]:.4g}-{qa[2]:.4g}) -> "
                    f"{qc[1]:.4g} ({qc[0]:.4g}-{qc[2]:.4g})"
                )
        for r in errors[:5]:
            lines.append(f"  error in member {r['member']}: {r.get('error', '')}")
        lines.append(f"  wall time per member: median {np.nanmedian(wall):.1f} s")
    return "\n".join(lines)


# --------------------------------------------------------------------------- figure


def _inner_orbits(rows) -> list[str]:
    names = {k.removesuffix("_e_max") for r in rows for k in r if k.endswith("_e_max")}
    return [n for n in ("A", "B", "inner", "mid") if n in names]


def plot_ensemble(
    rows: list[dict], path: str | Path, x: str | None = None, y: str | None = None
) -> Path:
    """Overview figure (PNG/PDF/SVG by extension): fates, high eccentricities, mutual
    inclinations or eclipses, and the dependence on one drawn parameter ``x``."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    params = [k for k in rows[0] if "." in k] if rows else []
    x = x or (params[0] if params else None)
    dirs = [d for d in ("forward", "backward") if any(r["direction"] == d for r in rows)]
    colors = {"forward": "tab:blue", "backward": "tab:orange"}
    inner = _inner_orbits(rows)
    e_all = np.array([r.get(f"{n}_e_max", np.nan) for r in rows for n in inner], dtype=float)
    high_e = bool(np.nanmax(e_all, initial=0.0) > 0.9)
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))

    # fates
    fates = [f for f in FATES if any(r["fate"] == f for r in rows)]
    w = 0.8 / max(len(dirs), 1)
    for k, d in enumerate(dirs):
        rs = [r for r in rows if r["direction"] == d]
        frac = [sum(r["fate"] == f for r in rs) / len(rs) for f in fates]
        ax[0, 0].bar(np.arange(len(fates)) + k * w, frac, w, label=d, color=colors[d])
    ax[0, 0].set_xticks(np.arange(len(fates)) + w * (len(dirs) - 1) / 2, fates, rotation=20)
    ax[0, 0].set_ylabel("fraction of members")
    ax[0, 0].set_title("fate")
    ax[0, 0].legend()

    # highest eccentricity
    if high_e:
        col, bins = "log1me_min", np.linspace(-7, 0, 36)
    else:
        col = "e_max"
        lo, hi = np.nanmin(e_all, initial=1.0), np.nanmax(e_all, initial=0.0)
        bins = np.linspace(lo, hi, 41) if hi > lo else 10
    for d in dirs:
        for name, ls in zip(inner, ("-", "--", ":")):
            v = np.array(
                [r.get(f"{name}_{col}", np.nan) for r in rows if r["direction"] == d],
                dtype=float,
            )
            v = v[np.isfinite(v)]
            if high_e:
                v = np.clip(v, -7, 0)
            ax[0, 1].hist(v, bins, histtype="step", ls=ls, color=colors[d], label=f"{name} ({d})")
    ax[0, 1].set_xlabel("log10(1 - e_max)" if high_e else "e_max")
    ax[0, 1].set_ylabel("members")
    ax[0, 1].set_title("highest eccentricity reached")
    ax[0, 1].legend(fontsize=8)

    # eclipses (2+2) or mutual inclination range
    if any("eclipse_fraction_both" in r for r in rows):
        for d in dirs:
            for key, ls in (
                ("A_eclipse_fraction", "-"),
                ("B_eclipse_fraction", "--"),
                ("eclipse_fraction_both", ":"),
            ):
                v = np.array([r.get(key, np.nan) for r in rows if r["direction"] == d], dtype=float)
                ax[1, 0].hist(
                    100 * v[np.isfinite(v)],
                    np.linspace(0, 100, 26),
                    histtype="step",
                    ls=ls,
                    color=colors[d],
                    label=f"{key.replace('_eclipse_fraction', '').replace('eclipse_fraction_', '')} ({d})",
                )
        ax[1, 0].set_xlabel("time with eclipses [% of the run]")
        ax[1, 0].set_ylabel("members")
        ax[1, 0].set_title("eclipses for the observer (A, B, both)")
        ax[1, 0].legend(fontsize=8)
    else:
        pair = next(
            (
                k.removeprefix("imut_").removesuffix("_0")
                for k in rows[0]
                if k.startswith("imut_") and k.endswith("_0")
            ),
            None,
        )
        if pair:
            for d in dirs:
                rs = [r for r in rows if r["direction"] == d and r["fate"] != "error"]
                i0 = np.array([r[f"imut_{pair}_0"] for r in rs])
                lo = np.array([r[f"imut_{pair}_min"] for r in rs])
                hi = np.array([r[f"imut_{pair}_max"] for r in rs])
                ax[1, 0].vlines(i0, lo, hi, color=colors[d], alpha=0.4, label=d)
            ax[1, 0].set_xlabel(f"initial mutual inclination {pair} [deg]")
            ax[1, 0].set_ylabel("range during the run [deg]")
            ax[1, 0].legend()

    # dependence on one drawn parameter
    if x:
        ys = [y] if y else [f"{n}_{'log1me_min' if high_e else 'e_max'}" for n in inner]
        for d in dirs:
            rs = [r for r in rows if r["direction"] == d and r["fate"] != "error"]
            for col, m in zip(ys, ("o", "s", "^")):
                ax[1, 1].scatter(
                    [r[x] for r in rs],
                    [r.get(col, np.nan) for r in rs],
                    s=8,
                    marker=m,
                    color=colors[d],
                    alpha=0.5,
                    label=f"{col} ({d})",
                )
        ax[1, 1].set_xlabel(x)
        ax[1, 1].set_ylabel(y or ("log10(1 - e_max)" if high_e else "e_max"))
        ax[1, 1].legend(fontsize=8)
    fig.suptitle(f"swingersim ensemble: {len({r['member'] for r in rows})} members")
    fig.tight_layout()
    path = Path(path)
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path
