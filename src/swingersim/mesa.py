"""Own stellar-model grids computed with MESA (Paxton et al. 2011-2019; Jermyn et al. 2023).

swingersim does not ship MESA; it drives a local installation (environment variables ``MESA_DIR``
and ``MESASDK_ROOT``, see the installation page of the manual):

* ``swingersim mesa-grid grid.yaml`` writes one MESA run directory per (mass, Z) of the grid,
  builds the standard ``star/work`` executable once, and runs the tracks in parallel (resumable;
  the worker can run it as a job). Every track starts on the pre-main sequence and writes
  structure profiles at regular model intervals.
* ``swingersim mesa-reduce grid.yaml`` reads the history and profiles and reduces every profile
  to the quantities swingersim needs, computed from the full radial structure:

  - radius, luminosity, effective temperature, age, mass, central hydrogen;
  - base of the convective envelope (radius R_c, mass M_c) and the eddy turnover time
    tau_c = integral dr / u_conv over the envelope;
  - apsidal-motion constant k2 from Radau's equation integrated over the profile
    (as in Claret 2023) and the moment of inertia I / (M R^2);
  - the frequency-averaged inertial-wave dissipation C = 1/(<Q'_IW> eps^2) of the convective
    envelope with the method of Barker (2020, 2022) (``swingersim.inertial_waves``).

  The result is ``mesa/<grid>/<grid>.npz`` (one record per profile).

Profile files follow the MESA format: two header lines (names, values), a blank line, column
numbers, column names, then one row per zone from the surface to the centre.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import yaml

from . import units
from .magnetic import MSUN_G, RSUN_CM

G_CGS = 6.67430e-8

HISTORY_COLUMNS = (
    "model_number",
    "star_age",
    "star_mass",
    "log_L",
    "log_Teff",
    "log_R",
    "log_g",
    "center_h1",
    "center_he4",
    "he_core_mass",
)
PROFILE_COLUMNS = (
    "zone",
    "radius",
    "mass",
    "logRho",
    "logP",
    "logT",
    "mixing_type",
    "gradr",
    "grada",
    "conv_vel",
    "brunt_N2",
)


# --------------------------------------------------------------------------- grid definition
@dataclass
class MesaGridSpec:
    name: str
    masses: list[float]
    z: list[float] = field(default_factory=lambda: [0.02])
    y: float | None = None  # None: MESA default helium (Y = 0.24 + 2 Z)
    alpha_mlt: float = 1.82  # solar calibration of MIST (Choi et al. 2016)
    stop: str = "tams"  # tams | rgb_tip | max_age
    max_age_yr: float = 1.4e10
    profile_interval: int = 10
    max_num_profile_models: int = 600
    overshoot_f_core: float = 0.016  # exponential core overshoot above overshoot_mass
    overshoot_mass: float = 1.1
    pre_ms_t_c: float = 5.0e5
    threads_per_run: int = 2
    extra_controls: dict = field(default_factory=dict)
    mesa_dir: str | None = None
    mesasdk_root: str | None = None

    @classmethod
    def load(cls, path: str | Path) -> MesaGridSpec:
        with open(path, encoding="utf-8") as f:
            d = yaml.safe_load(f)
        m = d.get("masses")
        if isinstance(m, dict):  # {from, to, step}
            n = round((m["to"] - m["from"]) / m["step"]) + 1
            d["masses"] = [round(m["from"] + k * m["step"], 4) for k in range(n)]
        for k in ("max_age_yr", "alpha_mlt", "overshoot_f_core", "overshoot_mass", "pre_ms_t_c"):
            if k in d:
                d[k] = float(d[k])  # YAML reads 1.4e10 (without a sign) as a string
        if "y" in d and d["y"] is not None:
            d["y"] = float(d["y"])
        d["z"] = [float(v) for v in d.get("z", [0.02])]
        spec = cls(**d)
        if spec.stop not in ("tams", "rgb_tip", "max_age"):
            raise ValueError("mesa grid: stop must be tams, rgb_tip or max_age")
        return spec

    def runs(self) -> list[tuple[float, float]]:
        return [(float(m), float(z)) for z in self.z for m in self.masses]


def run_name(mass: float, z: float) -> str:
    return f"M{mass:07.4f}_Z{z:.5f}"


def grid_dir(spec: MesaGridSpec, root: Path) -> Path:
    return Path(root).resolve() / "mesa" / spec.name  # absolute: MESA runs in the run folders


# --------------------------------------------------------------------------- inlists
def _fortran(v) -> str:
    if isinstance(v, bool):
        return ".true." if v else ".false."
    if isinstance(v, str):
        return f"'{v}'"
    if isinstance(v, int):
        return str(v)
    return f"{float(v):.10e}".replace("e", "d")


def inlist_text(spec: MesaGridSpec, mass: float, z: float, columns_dir: Path) -> str:
    """Complete MESA inlist of one track (namelists star_job, eos, kap, controls, pgstar)."""
    star_job = {
        "show_log_description_at_start": False,
        "create_pre_main_sequence_model": True,
        "pre_ms_T_c": spec.pre_ms_t_c,
        "history_columns_file": str(columns_dir / "history_columns.list"),
        "profile_columns_file": str(columns_dir / "profile_columns.list"),
        "save_model_when_terminate": False,
        "pgstar_flag": False,
    }
    kap = {"Zbase": z}
    controls = {
        "initial_mass": mass,
        "initial_z": z,
        "mixing_length_alpha": spec.alpha_mlt,
        "max_age": spec.max_age_yr,
        "history_interval": 1,
        "profile_interval": spec.profile_interval,
        "max_num_profile_models": spec.max_num_profile_models,
        "write_profiles_flag": True,
        "photo_interval": 100000,
        "log_directory": "LOGS",
    }
    if spec.y is not None:
        controls["initial_y"] = spec.y
    lines_extra = []
    if spec.stop == "tams":
        lines_extra += [
            "  xa_central_lower_limit_species(1) = 'h1'",
            "  xa_central_lower_limit(1) = 1d-4",
        ]
    elif spec.stop == "rgb_tip":
        controls["power_he_burn_upper_limit"] = 1.0
    if mass >= spec.overshoot_mass and spec.overshoot_f_core > 0:
        lines_extra += [
            "  overshoot_scheme(1) = 'exponential'",
            "  overshoot_zone_type(1) = 'burn_H'",
            "  overshoot_zone_loc(1) = 'core'",
            "  overshoot_bdy_loc(1) = 'top'",
            f"  overshoot_f(1) = {_fortran(spec.overshoot_f_core)}",
            f"  overshoot_f0(1) = {_fortran(0.5 * spec.overshoot_f_core)}",
        ]
    controls.update(spec.extra_controls)

    def block(name, d, extra=()):
        body = [f"  {k} = {_fortran(v)}" for k, v in d.items()] + list(extra)
        return f"&{name}\n" + "\n".join(body) + f"\n/ ! end of {name} namelist\n"

    return (
        f"! swingersim MESA grid '{spec.name}': M = {mass} Msun, Z = {z}\n"
        + block("star_job", star_job)
        + block("eos", {})
        + block("kap", kap)
        + block("controls", controls, lines_extra)
        + block("pgstar", {})
    )


def write_columns(dirpath: Path) -> None:
    dirpath.mkdir(parents=True, exist_ok=True)
    (dirpath / "history_columns.list").write_text("\n".join(HISTORY_COLUMNS) + "\n")
    (dirpath / "profile_columns.list").write_text("\n".join(PROFILE_COLUMNS) + "\n")


# --------------------------------------------------------------------------- running MESA
def _env(spec: MesaGridSpec) -> tuple[Path, Path]:
    mesa = spec.mesa_dir or os.environ.get("MESA_DIR")
    sdk = spec.mesasdk_root or os.environ.get("MESASDK_ROOT")
    if not mesa or not Path(mesa, "star", "work").is_dir():
        raise RuntimeError("MESA not found: set MESA_DIR (or mesa_dir in the grid file)")
    if not sdk or not Path(sdk, "bin", "mesasdk_init.sh").exists():
        raise RuntimeError(
            "MESA SDK not found: set MESASDK_ROOT (or mesasdk_root in the grid file)"
        )
    return Path(mesa), Path(sdk)


def _shell(cmd: str, cwd: Path, sdk: Path, mesa: Path, threads: int, log: Path) -> int:
    env = dict(os.environ, MESA_DIR=str(mesa), MESASDK_ROOT=str(sdk), OMP_NUM_THREADS=str(threads))
    full = f"source '{sdk}/bin/mesasdk_init.sh' && {cmd}"
    with open(log, "a", encoding="utf-8") as f:
        return subprocess.call(
            ["bash", "-c", full], cwd=cwd, env=env, stdout=f, stderr=subprocess.STDOUT
        )


def build_work(spec: MesaGridSpec, root: Path) -> Path:
    """Copy $MESA_DIR/star/work into the grid folder and compile it once (./mk)."""
    mesa, sdk = _env(spec)
    gdir = grid_dir(spec, root)
    work = gdir / "work"
    if not (work / "star").exists():
        if work.exists():
            shutil.rmtree(work)
        shutil.copytree(mesa / "star" / "work", work)
        code = _shell("./mk", work, sdk, mesa, spec.threads_per_run, gdir / "build.log")
        if code != 0 or not (work / "star").exists():
            raise RuntimeError(f"building the MESA work directory failed, see {gdir / 'build.log'}")
    return work


def _run_one(args) -> tuple[str, int]:
    spec_dict, root, mass, z = args
    spec = MesaGridSpec(**spec_dict)
    mesa, sdk = _env(spec)
    gdir = grid_dir(spec, Path(root))
    rdir = gdir / "runs" / run_name(mass, z)
    if (rdir / "done").exists():
        return run_name(mass, z), 0
    rdir.mkdir(parents=True, exist_ok=True)
    if (rdir / "LOGS").exists():
        shutil.rmtree(rdir / "LOGS")  # an interrupted run starts again from the pre-main sequence
    (rdir / "inlist").write_text(inlist_text(spec, mass, z, gdir))
    code = _shell(
        f"'{gdir / 'work' / 'star'}'", rdir, sdk, mesa, spec.threads_per_run, rdir / "run.log"
    )
    ok = finished(rdir, code)
    if ok:
        (rdir / "done").write_text(f"exit code {code}\n")
    return run_name(mass, z), 0 if ok else (code or 1)


def finished(rdir: Path, code: int) -> bool:
    """A MESA run that ended normally: its history and more than three profiles exist, and
    MESA exited with code 0 or reported its termination code ('termination code: ...') without
    an abnormal end ('terminated evolution: ...', printed e.g. for 'cannot find acceptable
    model' before the termination code min_timestep_limit)."""
    hist = rdir / "LOGS" / "history.data"
    if not hist.exists() or len(list((rdir / "LOGS").glob("profile*.data"))) <= 3:
        return False
    log = rdir / "run.log"
    text = log.read_text(errors="replace") if log.exists() else ""
    if "terminated evolution" in text:
        return False
    return code == 0 or "termination code:" in text


def run_grid(spec: MesaGridSpec, root: Path, processes: int = 0, test: bool = False) -> dict:
    """Build MESA once and compute all missing tracks; returns {run name: exit code}."""
    root = Path(root)
    gdir = grid_dir(spec, root)
    gdir.mkdir(parents=True, exist_ok=True)
    write_columns(gdir)
    build_work(spec, root)
    runs = spec.runs()
    if test:  # one short solar-mass track to check the installation and the inlists
        runs = [(1.0, spec.z[0])]
    from .ensemble import resolve_processes

    nproc = processes or max(1, resolve_processes(0) // max(spec.threads_per_run, 1))
    out = {}
    tasks = [(spec.__dict__, str(root), m, z) for m, z in runs]
    with ProcessPoolExecutor(max_workers=nproc) as ex:
        futs = [ex.submit(_run_one, t) for t in tasks]
        for f in as_completed(futs):
            name, code = f.result()
            out[name] = code
            print(f"{name}: {'ok' if code == 0 else f'FAILED (code {code})'}", flush=True)
    return out


# --------------------------------------------------------------------------- reading MESA output
def read_mesa_table(path: str | Path) -> tuple[dict, dict[str, np.ndarray]]:
    """(header dict, columns) of a MESA history or profile file."""
    lines = Path(path).read_text().splitlines()
    hnames, hvals = lines[1].split(), lines[2].split()
    header = {}
    for k, v in zip(hnames, hvals):
        try:
            header[k] = float(v.replace("D", "E"))
        except ValueError:
            header[k] = v.strip('"')
    names = lines[5].split()
    data = np.array(
        [[float(x.replace("D", "E")) for x in ln.split()] for ln in lines[6:] if ln.strip()]
    )
    return header, {n: data[:, i] for i, n in enumerate(names)}


# --------------------------------------------------------------------------- structure
def radau_k2(r: np.ndarray, rho: np.ndarray) -> float:
    """Apsidal-motion constant k2 of a spherical profile (centre to surface, any units):
    Radau's equation r eta' + 6 (rho/rho_mean)(eta + 1) + eta(eta - 1) = 6, eta(0) = 0,
    k2 = (3 - eta_s) / (2 (2 + eta_s))."""
    from scipy.integrate import solve_ivp

    m = (
        np.concatenate(
            [[0.0], np.cumsum(0.5 * (rho[1:] * r[1:] ** 2 + rho[:-1] * r[:-1] ** 2) * np.diff(r))]
        )
        * 4.0
        * np.pi
    )
    m += 4.0 / 3.0 * np.pi * rho[0] * r[0] ** 3
    ratio = rho / (3.0 * m / (4.0 * np.pi * r**3))
    edge = float(ratio[-1])

    def f(x, y):
        eta = y[0]
        rr = float(np.interp(x, r, ratio, left=1.0, right=edge))
        return [(6.0 - 6.0 * rr * (eta + 1.0) - eta * (eta - 1.0)) / x]

    sol = solve_ivp(f, (r[0], r[-1]), [0.0], rtol=1e-9, atol=1e-11, max_step=(r[-1] - r[0]) / 2000)
    eta = sol.y[0, -1]
    return float((3.0 - eta) / (2.0 * (2.0 + eta)))


def envelope(r, m, convective):
    """(index of the base, index of the top) of the convective envelope: the most massive
    convective zone whose top lies in the outer 10 % of the radius; None if there is none.
    Arrays from the centre to the surface."""
    zones, start = [], None
    for i, c in enumerate(convective):
        if c and start is None:
            start = i
        if (not c or i == len(convective) - 1) and start is not None:
            end = i if c else i - 1
            zones.append((start, end))
            start = None
    zones = [z for z in zones if r[z[1]] > 0.9 * r[-1]]
    if not zones:
        return None
    return max(zones, key=lambda z: m[z[1]] - (m[z[0] - 1] if z[0] > 0 else 0.0))


def reduce_profile(header: dict, cols: dict[str, np.ndarray], n_grid: int = 6000) -> dict:
    """Quantities swingersim needs from one MESA profile (see the module docstring)."""
    from scipy.interpolate import PchipInterpolator

    from .inertial_waves import barker_profile

    order = np.argsort(cols["radius"])
    r_sun = cols["radius"][order]
    m_sun = cols["mass"][order]
    rho = 10.0 ** cols["logRho"][order]
    p = 10.0 ** cols["logP"][order]
    conv = (
        cols["mixing_type"][order] == 1
        if "mixing_type" in cols
        else (cols["gradr"][order] > cols["grada"][order])
    )
    u_c = cols.get("conv_vel", np.zeros_like(r_sun))[order]
    radius, mass = float(r_sun[-1]), float(m_sun[-1])
    r_cm, m_g = radius * RSUN_CM, mass * MSUN_G
    x = r_sun / radius
    # dimensionless units G = M = R = 1
    rho_n = rho / (m_g / r_cm**3)
    p_n = p / (G_CGS * m_g**2 / r_cm**4)
    out = {
        "age_myr": float(header.get("star_age", np.nan)) / 1e6,
        "mass": mass,
        "radius_rsun": radius,
        "luminosity_lsun": float(header.get("photosphere_L", np.nan)),
        "teff_k": float(header.get("Teff", np.nan)),
        "center_h1": float(header.get("center_h1", np.nan)),
    }
    out["rg2"] = float(np.trapezoid(8.0 * np.pi / 3.0 * rho_n * x**4, x))
    out["apsidal_k2"] = radau_k2(x, rho_n)
    env = envelope(x, m_sun / mass, conv)
    if env is None:
        out.update(r_env=0.0, m_env=0.0, turnover_yr=0.0, iw_c=0.0, fully_convective=False)
        return out
    i0, i1 = env
    full = i0 == 0
    out["fully_convective"] = bool(full)
    out["r_env"] = radius * (1.0 - (0.0 if full else x[i0]))
    out["m_env"] = mass * (1.0 - (0.0 if full else m_sun[i0] / mass))
    seg = slice(i0, i1 + 1)
    with np.errstate(divide="ignore"):
        inv_u = np.where(u_c[seg] > 0, 1.0 / u_c[seg], 0.0)
    out["turnover_yr"] = float(np.trapezoid(inv_u, r_sun[seg] * RSUN_CM)) / units.SECONDS_PER_YEAR
    # inertial waves: smooth resampling of the profile on a fine grid
    keep = np.concatenate([[True], np.diff(x) > 0])
    xs, lr, lp = x[keep], np.log(rho_n[keep]), np.log(p_n[keep])
    grid = np.unique(
        np.concatenate(
            [
                np.linspace(xs[0], xs[-1], n_grid),
                np.clip(x[i0] + np.linspace(-0.02, 0.02, n_grid // 6), xs[0], xs[-1]),
            ]
        )
    )
    f_r, f_p = PchipInterpolator(xs, lr), PchipInterpolator(xs, lp)
    rho_g, p_g = np.exp(f_r(grid)), np.exp(f_p(grid))
    drho_dp = rho_g / p_g * f_r.derivative()(grid) / f_p.derivative()(grid)
    top = x[i1] if i1 < len(x) - 1 else grid[-1]
    cmask = (grid >= (x[i0] if not full else grid[0])) & (grid <= top)
    out["iw_c"] = barker_profile(grid, rho_g, drho_dp, cmask)
    return out


def reduce_grid(spec: MesaGridSpec, root: Path, processes: int = 0) -> Path:
    """Reduce every profile of every finished track; writes mesa/<grid>/<grid>.npz."""
    gdir = grid_dir(spec, Path(root))
    tasks = []
    skipped = []
    for mass, z in spec.runs():
        rdir = gdir / "runs" / run_name(mass, z)
        if not (rdir / "done").exists():  # missing, running, crashed or interrupted
            if rdir.exists():
                skipped.append(rdir.name)
            continue
        logs = rdir / "LOGS"
        for prof in sorted(logs.glob("profile*.data"), key=lambda p: int(p.stem[7:])):
            tasks.append((mass, z, str(prof)))
    if skipped:
        import warnings

        warnings.warn(
            f"MESA grid {spec.name}: runs without a normal end are left out: {skipped}",
            stacklevel=2,
        )
    if not tasks:
        raise RuntimeError(f"no finished MESA runs with profiles under {gdir / 'runs'}")
    rows = []
    from .ensemble import resolve_processes

    with ProcessPoolExecutor(max_workers=resolve_processes(processes or 0)) as ex:
        futs = {ex.submit(_reduce_file, t[2]): t for t in tasks}
        for f in as_completed(futs):
            mass, z, _ = futs[f]
            d = f.result()
            d.update(initial_mass=mass, z=z)
            rows.append(d)
    rows.sort(key=lambda d: (d["z"], d["initial_mass"], d["age_myr"]))
    keys = sorted({k for d in rows for k in d})
    out = gdir / f"{spec.name}.npz"
    np.savez_compressed(out, **{k: np.array([d.get(k, np.nan) for d in rows]) for k in keys})
    return out


def _reduce_file(path: str) -> dict:
    h, c = read_mesa_table(path)
    return reduce_profile(h, c)
