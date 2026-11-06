"""Storage of single runs (HDF5).

Layout (units in the attribute ``unit`` of every dataset; angles are radians unless the name
ends in ``_deg``; time in Julian years of 365.25 d, lengths in AU, masses in Msun):

* ``time_yr``, ``energy_error``, ``angmom_error``, ``megno`` - one value per output;
* ``orbits/<orbit>/<a, e, inc, Omega, omega, M, inc_inv, Omega_inv>`` - osculating Jacobi
  elements of every orbit (sky frame; ``*_inv`` relative to the invariable plane at t = 0);
* ``mutual_inclination/<orbit1>-<orbit2>`` - angles between the orbit normals;
* ``stars/<mK>`` - attributes of the stellar model (radius, k2, rg2, lag ...), the masses
  ``mass_msun`` (dynamical, at t = 0) and ``mass_msun_config`` (as configured: the ZAMS labels
  of the tracks with track mass loss), and series (spins, evolving stars, masses with mass loss);
* ``hierarchy`` - attributes ``configured``, ``final`` (detected at the last output),
  ``transitions_json`` (confirmed changes) and ``detected_changes_json`` (every output at which
  the detected hierarchy differs from the previous one, before the persistence check);
* ``extremes``, ``diagnostics`` and the JSON attributes ``binaries_json``, ``system_json``.

Root attributes: ``swingersim_version``, ``rebound_version``, ``reboundx_version``,
``numpy_version``, ``h5py_version``, ``python_version``, ``created_utc``, ``wall_s`` (when
known), ``config_yaml``, the event and the flags of the run. The numerical series are stored
with gzip + shuffle compression (read transparently by h5py).
"""

from __future__ import annotations

import datetime as dt
import json
import platform
from dataclasses import asdict
from pathlib import Path

import h5py
import numpy as np
import yaml

from . import __version__
from .backends.nbody import RunResult
from .config import Config

COMPRESS_MIN = 64  # smaller datasets are stored uncompressed (compression costs more there)

ELEMENT_UNITS = {
    "a": "AU",
    "e": "1",
    "inc": "rad",
    "Omega": "rad",
    "omega": "rad",
    "M": "rad",
    "inc_inv": "rad",
    "Omega_inv": "rad",
}
CONVENTIONS = {
    "time": "yr (Julian year, 365.25 d), from the start of the run (t = 0)",
    "length": "AU",
    "mass": "Msun",
    "angle": "rad unless the name ends in _deg",
    "G": "Gaussian, (0.01720209895 * 365.25)^2 AU^3 Msun^-1 yr^-2",
    "frame": "x-y = reference (sky) plane, z its normal along the line of sight; Jacobi orbits",
}


def series_unit(name: str) -> str:
    """Unit of a stored series from its name (suffix convention of swingersim)."""
    if name in ELEMENT_UNITS:
        return ELEMENT_UNITS[name]
    for suffix, unit in (
        ("_yr", "yr"),
        ("_d", "d"),
        ("_deg", "deg"),
        ("_msun", "Msun"),
        ("_rsun", "Rsun"),
        ("_lsun", "Lsun"),
        ("_k", "K"),
        ("_au", "AU"),
        ("_myr", "Myr"),
    ):
        if name.endswith(suffix):
            return unit
    return "1"


def _put(group, name: str, data, unit: str | None = None) -> None:
    """A dataset, compressed (gzip + shuffle) when numerical and not tiny, with its unit."""
    arr = np.asarray(data)
    if arr.dtype.kind in "fiub" and arr.ndim >= 1 and arr.size >= COMPRESS_MIN:
        ds = group.create_dataset(
            name, data=arr, compression="gzip", compression_opts=4, shuffle=True
        )
    else:
        ds = group.create_dataset(name, data=arr)
    ds.attrs["unit"] = unit if unit is not None else series_unit(name)


def _versions() -> dict[str, str]:
    out = {"swingersim_version": __version__, "numpy_version": np.__version__}
    out["h5py_version"] = h5py.__version__
    out["python_version"] = platform.python_version()
    for mod in ("scipy", "rebound", "reboundx", "emcee"):
        try:
            out[f"{mod}_version"] = str(__import__(mod).__version__)
        except Exception:  # noqa: BLE001 - a missing package is recorded, not fatal
            out[f"{mod}_version"] = ""
    return out


def software_versions() -> str:
    """One line with the versions of SWINGERsim and the packages that compute its results (for
    the headers of the fit outputs; the run files store them as attributes)."""
    v = _versions()
    name = {"swingersim": "SWINGERsim"}
    return ", ".join(
        f"{name.get(k.removesuffix('_version'), k.removesuffix('_version'))} {x}"
        for k, x in v.items()
        if x
    )


def _star_masses(cfg: Config) -> tuple[dict[str, float], dict[str, float]]:
    """(configured, dynamical at t = 0) masses of the stars."""
    from . import massloss

    configured = {k: float(v) for k, v in cfg.masses.items()}
    try:
        rc = massloss.resolve_masses(cfg)
        names = list(rc.star_names)
        m0 = massloss.initial_masses(rc, np.array([rc.masses[k] for k in names], float))
        return configured, {k: float(m) for k, m in zip(names, m0, strict=True)}
    except Exception:  # noqa: BLE001 - e.g. tracks unavailable: the configured masses
        return configured, dict(configured)


def _detected_changes(times, detected: list[str]) -> list[dict]:
    """The outputs at which the detected hierarchy differs from the previous output."""
    out, last = [], None
    for t, h in zip(times, detected, strict=False):
        if h != last:
            out.append({"time_yr": float(t), "to": h})
            last = h
    return out


def save_run(path: str | Path, cfg: Config, res: RunResult) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        for k, v in _versions().items():
            f.attrs[k] = v
        f.attrs["created_utc"] = dt.datetime.now(dt.UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
        if getattr(res, "wall_s", None) is not None:
            f.attrs["wall_s"] = float(res.wall_s)
        f.attrs["units_json"] = json.dumps(CONVENTIONS)
        f.attrs["config_yaml"] = yaml.safe_dump(cfg.model_dump(), sort_keys=False)
        f.attrs["event"] = res.event or ""
        if res.event_time is not None:
            f.attrs["event_time_yr"] = res.event_time
        for k, v in res.event_detail.items():
            f.attrs[f"event_{k}"] = str(v)
        _put(f, "time_yr", res.times, "yr")
        _put(f, "energy_error", res.energy_error, "1")
        f.attrs["dissipative"] = res.dissipative
        f.attrs["gas_phase"] = res.gas_phase
        f.attrs["mass_loss"] = res.mass_loss
        f.attrs["backend"] = res.backend
        if res.angmom_error is not None:
            _put(f, "angmom_error", res.angmom_error, "1")
        st = f.create_group("stars")
        m_cfg, m_dyn = _star_masses(cfg)
        for name, pr in res.stars.items():
            g = st.create_group(name)
            for k, v in asdict(pr).items():
                if v is not None:
                    g.attrs[k] = v
            if name in m_dyn:
                g.attrs["mass_msun"] = m_dyn[name]
                g.attrs["mass_msun_config"] = m_cfg[name]
            if name in res.roche:
                for k, v in res.roche[name].items():
                    g.attrs[f"roche_{k}"] = v
            for k, v in res.spins.get(name, {}).items():
                _put(g, k, v)
            for k, v in res.stellar.get(name, {}).items():
                _put(g, k, v, "1" if k == "phase" else None)
        f.attrs["binaries_json"] = json.dumps(res.binaries)
        g = f.create_group("orbits")
        for name, series in res.orbits.items():
            sg = g.create_group(name)
            for k, s in series.items():
                _put(sg, k, s)
        m = f.create_group("mutual_inclination")
        for k, s in res.mutual.items():
            _put(m, k, s, "rad")
        h = f.create_group("hierarchy")
        h.attrs["configured"] = res.configured_hierarchy
        h.attrs["final"] = res.hierarchy[-1] if res.hierarchy else ""
        h.attrs["transitions_json"] = json.dumps(res.transitions)
        h.attrs["detected_changes_json"] = json.dumps(_detected_changes(res.times, res.hierarchy))
        if res.megno is not None:
            _put(f, "megno", res.megno, "1")
            f.attrs["lyapunov_time_yr"] = res.lyapunov_time_yr
        ex = f.create_group("extremes")
        for name, vals in res.extremes.items():
            g = ex.create_group(name)
            for k, v in vals.items():
                g.attrs[k] = v
        f.attrs["system_json"] = json.dumps(res.system)
        d = f.create_group("diagnostics")
        for diag in res.diagnostics:
            g = d.create_group(diag.inner)
            for k, v in diag.as_dict().items():
                g.attrs[k] = v
