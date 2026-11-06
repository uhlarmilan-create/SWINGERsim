"""Stellar structure parameters needed by the dynamics: radius, apsidal constant, gyration.

A *stellar model* maps a mass to :class:`StarProperties`. Three kinds are available:

``zams`` (default)
    Radius of a zero-age main-sequence star of solar metallicity from the fit of
    Tout et al. (1996), MNRAS 281, 257; apsidal-motion constant and moment of inertia of a
    polytrope (n = 1.5 for fully convective stars below 0.35 Msun, n = 3 above), obtained by
    integrating the Lane-Emden and Radau equations.
``table``
    A CSV file supplied by the user with columns ``mass_msun, radius_rsun, apsidal_k2, rg2``
    (optionally ``lag_s``), interpolated linearly in mass. This is the way to use e.g. the
    tables of Claret (2019) or an own MESA grid.
``python``
    A function written by the user, ``"file.py:function"`` or ``"package.module:function"``,
    called as ``function(mass_msun)`` and returning a dict with the same keys.

Values given for an individual star in the ``stars`` section of the configuration always take
precedence over the model.

Conventions: ``apsidal_k2`` is the apsidal-motion constant k2 as tabulated by Claret; the
potential Love number used by REBOUNDx is twice that. ``rg2`` = I / (M R^2).
"""

from __future__ import annotations

import csv
import importlib
import importlib.util
import math
from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from .evolution import PHASES, SSE_Z_SUN, TEFF_SUN, grid_limits, track

LOGG_SUN = 4.438  # IAU 2015 nominal solar surface gravity [log cm s^-2]
FULLY_CONVECTIVE_MSUN = 0.35  # main-sequence stars below are fully convective (polytrope n = 1.5)


@dataclass(frozen=True)
class StarProperties:
    radius_rsun: float
    apsidal_k2: float  # apsidal-motion constant (Love number k2 = 2 * apsidal_k2)
    rg2: float  # I / (M R^2)
    lag_s: float = 0.0  # constant tidal time lag [s], Mignard (1979) convention
    # optional quantities of evolving stellar models (kind sse; ZAMS values for kind zams)
    luminosity_lsun: float | None = None
    teff_k: float | None = None
    age_myr: float | None = None
    phase: str | None = None  # "main sequence", "Hertzsprung gap", "giant branch"
    m_env: float | None = None  # convective envelope mass [Msun]
    r_env: float | None = None  # convective envelope depth [Rsun]
    turnover_yr: float | None = None  # convective eddy turnover time [yr]


# --- Tout et al. (1996): ZAMS radius, Z = 0.02 (Table 2 coefficients at log Z/Zsun = 0) ---
_TOUT_R = {
    "theta": 1.71535900,
    "iota": 6.59778800,
    "kappa": 10.08855000,
    "lam": 1.01249500,
    "mu": 0.07490166,
    "nu": 0.01077422,
    "xi": 3.08223400,
    "omicron": 17.84778000,
    "pi": 0.00022582,
}


def _clipf(x, lo, hi) -> float:
    """float(np.clip(x, lo, hi)) for scalars without the numpy call (hot paths of the tracks)."""
    return float(min(max(x, lo), hi))


def zams_radius_tout(m: float) -> float:
    """ZAMS radius [Rsun] for solar metallicity, Tout et al. (1996) eq. (2)."""
    c = _TOUT_R
    num = (
        c["theta"] * m**2.5
        + c["iota"] * m**6.5
        + c["kappa"] * m**11
        + c["lam"] * m**19
        + c["mu"] * m**19.5
    )
    den = c["nu"] + c["xi"] * m**2 + c["omicron"] * m**8.5 + m**18.5 + c["pi"] * m**19.5
    return float(num / den)


# --- polytropes: Lane-Emden + Radau equation ---
@lru_cache(maxsize=16)
def polytrope_constants(n: float) -> tuple[float, float]:
    """(apsidal-motion constant k2, I/(M R^2)) of a polytrope of index n.

    Lane-Emden: theta'' + 2 theta'/xi + theta^n = 0.  Radau's equation for the logarithmic
    derivative eta of the tidal distortion:  xi eta' + 6 (rho/rho_mean)(eta + 1) + eta(eta - 1)
    = 6, eta(0) = 0, and k2 = (3 - eta_s) / (2 (2 + eta_s)) at the surface.
    """

    def rhs(x, y):
        th, dth, eta, _inertia = y
        thn = max(th, 0.0) ** n
        ratio = x * thn / (-3.0 * dth)  # rho / mean density inside x
        return [
            dth,
            -thn - 2.0 * dth / x,
            (6.0 - 6.0 * ratio * (eta + 1.0) - eta * (eta - 1.0)) / x,
            thn * x**4,
        ]

    x0 = 1e-6
    y0 = [1.0 - x0**2 / 6.0, -x0 / 3.0, 0.0, 0.0]
    surface = lambda x, y: y[0]
    surface.terminal = True
    sol = solve_ivp(rhs, (x0, 50.0), y0, rtol=1e-11, atol=1e-13, events=surface)
    x1 = sol.t_events[0][0]
    _th, dth, eta, inertia = sol.y_events[0][0]
    k2 = (3.0 - eta) / (2.0 * (2.0 + eta))
    rg2 = (2.0 / 3.0) * inertia / (x1**4 * (-dth))
    return float(k2), float(rg2)


def zams_properties(m: float) -> StarProperties:
    n = 1.5 if m < FULLY_CONVECTIVE_MSUN else 3.0
    k2, rg2 = polytrope_constants(n)
    extra = {}
    lim = grid_limits()["mass"]
    if lim[0] <= m <= lim[1]:  # luminosity and Teff of the same ZAMS models (Tout et al. 1996)
        st = track(m, SSE_Z_SUN).at(0.0)
        extra = {
            "luminosity_lsun": st.luminosity_lsun,
            "teff_k": TEFF_SUN * st.luminosity_lsun**0.25 / math.sqrt(zams_radius_tout(m)),
            "age_myr": 0.0,
            "phase": "main sequence",
            "m_env": st.m_env,
            "r_env": st.r_env,
            "turnover_yr": st.convective_turnover_yr,
        }
    return StarProperties(radius_rsun=zams_radius_tout(m), apsidal_k2=k2, rg2=rg2, **extra)


@lru_cache(maxsize=1)
def _claret2023() -> dict[str, np.ndarray]:
    path = Path(__file__).with_name("data") / "claret2023_structure.npz"
    with np.load(path) as d:
        return {k: d[k] for k in d.files}


def claret_structure(
    mass: float, logg: float, feh: float = 0.0, segment: str = "main"
) -> tuple[float, float]:
    """(apsidal constant k2, I/(M R^2)) of a star of initial mass `mass` [Msun] at surface
    gravity `logg` [cgs] and [Fe/H] from the MESA models of Claret (2023), A&A 674, A67 (tables
    1-69, 0.2-8 Msun, [Fe/H] = -0.5, 0, +0.5), interpolated in log g (as Claret recommends)
    along one monotonic segment of each track, then in log M and [Fe/H]. Segments: 'pms'
    (pre-main sequence, start of the model to the ZAMS), 'main' (ZAMS to the tip of the giant
    branch), 'cheb' (core helium burning). A log g outside a segment takes its nearest end."""
    g = _claret2023()
    masses, fehs = g["mass"], g["feh"]
    if not masses[0] <= mass <= masses[-1]:
        raise ValueError(f"mass {mass} outside the Claret (2023) grid ({masses[0]}-{masses[-1]})")
    if segment not in ("pms", "main", "cheb"):
        raise ValueError(f"unknown Claret (2023) segment {segment}")
    lg, lk, bt = (g[f"{segment}_{q}"] for q in ("logg", "log_k2", "beta"))
    from .gridclip import clip_to_grid

    fe = clip_to_grid(feh, float(fehs[0]), float(fehs[-1]), "[Fe/H]", "Claret (2023) grid")
    iz = int(_clipf(np.searchsorted(fehs, fe) - 1, 0, len(fehs) - 2))
    wz = (fe - fehs[iz]) / (fehs[iz + 1] - fehs[iz])
    lm = np.log(masses)
    im = int(_clipf(np.searchsorted(lm, np.log(mass)) - 1, 0, len(masses) - 2))
    wm = (np.log(mass) - lm[im]) / (lm[im + 1] - lm[im])
    k_out = b_out = 0.0
    for a, wa in ((iz, 1 - wz), (iz + 1, wz)):
        for b, wb in ((im, 1 - wm), (im + 1, wm)):
            w = wa * wb
            if w == 0.0:
                continue
            gg = lg[a, b]  # increasing
            if not np.all(np.isfinite(gg)):
                raise ValueError(
                    f"Claret (2023) has no '{segment}' segment for {masses[b]} Msun, "
                    f"[Fe/H] = {fehs[a]}"
                )
            k_out += w * float(np.interp(logg, gg, lk[a, b]))
            b_out += w * float(np.interp(logg, gg, bt[a, b]))
    return 10.0**k_out, b_out**2


@lru_cache(maxsize=1)
def _claret2004() -> dict[str, np.ndarray]:
    path = Path(__file__).with_name("data") / "claret2004_structure.npz"
    with np.load(path) as d:
        return {k: d[k] for k in d.files}


def claret2004_structure(mass: float, logg: float, segment: str = "main") -> tuple[float, float]:
    """(apsidal constant k2, I/(M R^2)) of a massive star from the models of Claret (2004),
    A&A 424, 919 (0.8-125 Msun, Z = 0.02, ZAMS to carbon exhaustion), interpolated in log g along
    one monotonic segment of each track ('main': ZAMS to the largest radius; 'cheb': core helium
    burning, where the track has one), then linearly in log M. SWINGERsim uses it above the 8 Msun
    of Claret (2023). Above 125 Msun the 125 Msun track is used (no models exist); without a
    'cheb' segment for a mass the 'main' one is used; a log g outside a segment takes its
    nearest end."""
    g = _claret2004()
    masses = g["mass"]
    m = _clipf(mass, masses[0], masses[-1])
    lm = np.log(masses)
    im = int(_clipf(np.searchsorted(lm, math.log(m)) - 1, 0, len(masses) - 2))
    wm = (math.log(m) - lm[im]) / (lm[im + 1] - lm[im])
    k_out = b_out = 0.0
    for b, wb in ((im, 1 - wm), (im + 1, wm)):
        if wb == 0.0:
            continue
        seg = segment if np.all(np.isfinite(g[f"{segment}_logg"][b])) else "main"
        gg, lk, bt = (g[f"{seg}_{q}"][b] for q in ("logg", "log_k2", "beta"))
        k_out += wb * float(np.interp(logg, gg, lk))
        b_out += wb * float(np.interp(logg, gg, bt))
    return 10.0**k_out, b_out**2


def claret_covers(mass: float) -> bool:
    """True if the Claret (2023) grid covers the initial mass [Msun]."""
    masses = _claret2023()["mass"]
    return bool(masses[0] <= mass <= masses[-1])


def _structure(kind: str, mass: float, radius: float, feh: float, fallback):
    """k2 and rg2 from Claret (2023) or the polytrope fallback."""
    if kind == "claret2023":
        logg = LOGG_SUN + math.log10(mass) - 2.0 * math.log10(radius)
        return claret_structure(mass, logg, feh)
    return fallback


def sse_properties(m: float, age_myr: float, z: float) -> StarProperties:
    """Star of initial mass m [Msun] at age_myr after Hurley et al. (2000); moment of inertia
    after Hurley et al. (2002); apsidal constant of a polytrope (n = 1.5 below 0.35 Msun, 3 above)
    until stellar-structure tables are configured."""
    st = track(m, z).at(age_myr)
    n = 1.5 if m < FULLY_CONVECTIVE_MSUN else 3.0
    k2, _ = polytrope_constants(n)
    return StarProperties(
        radius_rsun=st.radius_rsun,
        apsidal_k2=k2,
        rg2=st.rg2,
        luminosity_lsun=st.luminosity_lsun,
        teff_k=st.teff_k,
        age_myr=st.age_myr,
        phase=PHASES[st.phase],
        m_env=st.m_env,
        r_env=st.r_env,
        turnover_yr=st.convective_turnover_yr,
    )


# --- user tables and user functions ---
@lru_cache(maxsize=8)
def _read_table(path: str) -> dict[str, np.ndarray]:
    with open(path, newline="", encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(f)]
    if not rows:
        raise ValueError(f"stellar table {path} is empty")
    cols = {k.strip(): np.array([float(r[k]) for r in rows]) for k in rows[0]}
    need = {"mass_msun", "radius_rsun", "apsidal_k2", "rg2"}
    if not need <= set(cols):
        raise ValueError(f"stellar table {path} needs columns {sorted(need)}")
    order = np.argsort(cols["mass_msun"])
    return {k: v[order] for k, v in cols.items()}


def table_properties(m: float, path: str) -> StarProperties:
    t = _read_table(path)
    if not t["mass_msun"][0] <= m <= t["mass_msun"][-1]:
        raise ValueError(f"mass {m} outside the stellar table {path}")
    v = {k: float(np.interp(m, t["mass_msun"], col)) for k, col in t.items()}
    return StarProperties(v["radius_rsun"], v["apsidal_k2"], v["rg2"], v.get("lag_s", 0.0))


def _load_function(entry: str, base: Path):
    target, _, func = entry.rpartition(":")
    if not target or not func:
        raise ValueError(
            "stellar_model.entry must look like 'file.py:function' or 'pkg.mod:function'"
        )
    if target.endswith(".py"):
        path = (base / target).resolve()
        spec = importlib.util.spec_from_file_location(path.stem, path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    else:
        module = importlib.import_module(target)
    return getattr(module, func)


def python_properties(m: float, entry: str, base: Path) -> StarProperties:
    out = _load_function(entry, base)(m)
    if isinstance(out, StarProperties):
        return out
    return StarProperties(
        float(out["radius_rsun"]),
        float(out["apsidal_k2"]),
        float(out["rg2"]),
        float(out.get("lag_s", 0.0)),
    )


def stellar_track(sm, mass: float):
    """Evolution track of a star for stellar_model sm (kind sse or mist)."""
    if sm.kind == "mist":
        from .mist import mist_track

        return mist_track(float(mass), sm.feh, sm.vvcrit)
    return track(float(mass), sm.z)


def mist_properties(m: float, age_myr: float, feh: float, vvcrit: float) -> StarProperties:
    """Star of initial mass m [Msun] at age_myr from the MIST v1.2 tracks (0.1-300 Msun); k2 and
    I/(MR^2) from the matching segment (pre-main sequence, main sequence/giant branch, core
    helium burning) of Claret (2023) for 0.2-8 Msun, of Claret (2004) above 8 Msun, and of a
    polytrope n = 1.5 below 0.2 Msun."""
    from .mist import MIST_PHASES, mist_track

    tr = mist_track(float(m), feh, vvcrit)
    st = tr.at(age_myr)
    # surface gravity with the current mass (massive stars lose a large part of their mass)
    logg = LOGG_SUN + math.log10(tr.current_mass(age_myr)) - 2.0 * math.log10(st.radius_rsun)
    # AGB stars (phases 5-6) take the structure of the tip of the giant branch (the end of the
    # 'main' segment, which log g below it selects): Claret's grids have no AGB models
    segment = {0: "pms", 4: "cheb"}.get(st.phase, "main")
    if claret_covers(m):
        try:
            k2, rg2 = claret_structure(m, logg, feh, segment)
        except ValueError:
            k2, rg2 = claret_structure(m, logg, feh, "main")
    elif m > 8.0:
        k2, rg2 = claret2004_structure(m, logg, "cheb" if segment == "cheb" else "main")
    else:
        k2, rg2 = polytrope_constants(1.5 if m < FULLY_CONVECTIVE_MSUN else 3.0)
    m_env, r_env = mist_envelope(tr, age_myr, feh)
    turnover = 0.0
    if m_env > 0.0 and r_env > 0.0:  # Hurley et al. (2002) eq. 31 (Rasio et al. 1996)
        x = m_env * r_env * (st.radius_rsun - 0.5 * r_env)
        turnover = 0.4311 * (x / (3.0 * st.luminosity_lsun)) ** (1.0 / 3.0)
    return StarProperties(
        radius_rsun=st.radius_rsun,
        apsidal_k2=k2,
        rg2=rg2,
        luminosity_lsun=st.luminosity_lsun,
        teff_k=st.teff_k,
        age_myr=st.age_myr,
        phase=MIST_PHASES[st.phase],
        m_env=m_env,
        r_env=r_env,
        turnover_yr=turnover,
    )


GIANT_TEFF_CONV = (5500.0, 7000.0)  # core-He-burning stars: convective below, radiative above


def mist_envelope(tr, age_myr: float, feh: float) -> tuple[float, float]:
    """Mass [Msun] and depth [Rsun] of the convective envelope of a MIST star, which the MIST
    tracks do not contain. Giants (giant branch, AGB; core helium burning while cool): the
    whole envelope above the helium core, M_env = M - M_c, depth ~ R (Hurley et al. 2002,
    Sect. 2.3). Main sequence and Hertzsprung gap: the envelope of the SSE star of the same
    mass at the same fraction of its main-sequence life (Hurley et al. 2002, fitted to the
    detailed models of Pols et al. 1998), with Z = 0.02 x 10^[Fe/H]."""
    from .mist import EEP_RGB_TIP, EEP_TAMS, EEP_ZAMS

    st = tr.at(age_myr)
    m_now = tr.current_mass(age_myr)
    if st.phase in (3, 5, 6) or (st.phase == 4 and st.teff_k < GIANT_TEFF_CONV[1]):
        f = 1.0
        if st.phase == 4:  # blue loop: the convective envelope shrinks as the star heats up
            lo, hi = GIANT_TEFF_CONV
            f = _clipf((hi - st.teff_k) / (hi - lo), 0.0, 1.0)
        return max(m_now - st.m_core, 0.0) * f, st.radius_rsun * f
    if st.phase == 0 and 0.1 <= tr.mass <= 2.5:
        return _pms_envelope(tr, age_myr, feh, st, m_now)
    if st.phase not in (1, 2) or not 0.1 <= tr.mass <= 50.0:
        return 0.0, 0.0
    t_zams = float(tr._age[EEP_ZAMS - 1])
    t_tams = float(tr._age[EEP_TAMS - 1]) if tr.n_eep >= EEP_TAMS else tr.t_end
    lim = grid_limits()
    z = _clipf(0.02 * 10.0**feh, *lim["z"])
    m = _clipf(tr.mass, *lim["mass"])
    try:
        sse = track(m, z)
        if st.phase == 1:
            frac = (age_myr - t_zams) / max(t_tams - t_zams, 1e-12)
            s = sse.at(_clipf(frac, 0.0, 1.0) * sse.t_ms * (1.0 - 1e-9))
        else:  # Hertzsprung gap: the SSE gap at the same fraction
            e_hg = getattr(tr, "eep_gb", 0.5 * (EEP_TAMS + EEP_RGB_TIP))  # end of phase 2
            t_rgb = float(np.interp(min(e_hg, tr.n_eep), tr._x, tr._age))
            frac = _clipf((age_myr - t_tams) / max(t_rgb - t_tams, 1e-12), 0.0, 1.0)
            s = sse.at(sse.t_ms + frac * (sse.t_bgb - sse.t_ms) * (1.0 - 1e-9))
    except ValueError:
        return 0.0, 0.0
    # the same fractions of mass and radius as the SSE star
    return s.m_env / m * m_now, s.r_env / s.radius_rsun * st.radius_rsun


def _pms_envelope(tr, age_myr: float, feh: float, st, m_now: float) -> tuple[float, float]:
    """Convective envelope of a pre-main-sequence MIST star: fully convective on the Hayashi
    track; a radiative core appears where the track turns from the Hayashi line (minimum Teff)
    and the envelope shrinks, linearly in the equivalent evolutionary point, to its ZAMS size
    (the SSE envelope at the ZAMS). Stars below 0.35 Msun stay fully convective."""
    from .mist import EEP_ZAMS

    if tr.mass <= FULLY_CONVECTIVE_MSUN:
        return m_now, st.radius_rsun
    n = min(EEP_ZAMS, tr.n_eep)
    cache = tr.__dict__.setdefault("_pms_env", {})
    if "onset" not in cache:
        teff = [tr.at_eep(e).teff_k for e in range(1, n + 1)]
        cache["onset"] = float(np.argmin(teff) + 1)
        lim = grid_limits()
        z = _clipf(0.02 * 10.0**feh, *lim["z"])
        m = _clipf(tr.mass, *lim["mass"])
        try:
            s0 = track(m, z).at(0.0)
            cache["zams"] = (s0.m_env / m, s0.r_env / s0.radius_rsun)
        except ValueError:
            cache["zams"] = (0.0, 0.0)
    eep = tr.eep_of(age_myr)
    onset, (fm, fr) = cache["onset"], cache["zams"]
    x = _clipf((eep - onset) / max(n - onset, 1e-9), 0.0, 1.0)
    return m_now * (1.0 - x * (1.0 - fm)), st.radius_rsun * (1.0 - x * (1.0 - fr))


Z_SUN_MESA = 0.0142  # present-day solar Z of Asplund et al. (2009), as in MIST


def find_mesa_grid(name: str, base: Path | None = None) -> Path:
    """mesa/<name>/<name>.npz of an own MESA grid (swingersim mesa-reduce), searched next to the
    configuration and in its parent folders, then in the working folder."""
    starts = []
    if base is not None:
        b = Path(base).resolve()
        starts += [b, *list(b.parents)[:5]]
    starts.append(Path.cwd())
    for d in starts:
        p = d / "mesa" / name / f"{name}.npz"
        if p.is_file():
            return p
    raise ValueError(
        f"stellar_model.structure mesa:{name}: mesa/{name}/{name}.npz not found "
        "(compute the grid with swingersim mesa-grid and mesa-reduce)"
    )


@lru_cache(maxsize=4)
def _mesa_grid(path: str) -> dict[str, np.ndarray]:
    with np.load(path) as d:
        return {k: d[k] for k in d.files}


def mesa_structure(
    path: str | Path, mass: float, age_myr: float | None, z: float
) -> tuple[float, float]:
    """(apsidal constant k2, I/(M R^2)) of a star of initial mass ``mass`` [Msun] at ``age_myr``
    (from the start of the track, as MESA and MIST count; None = the zero-age main sequence)
    from an own reduced MESA grid: the tracks of the metallicity nearest to ``z``, interpolated
    in age along the two tracks that bracket the mass, then linearly in log M. Ages beyond a
    track take its last profile."""
    g = _mesa_grid(str(path))
    zs = np.unique(g["z"])
    zz = zs[int(np.argmin(np.abs(np.log(zs) - math.log(z))))]
    sel = g["z"] == zz
    masses = np.unique(g["initial_mass"][sel])
    if not masses[0] <= mass <= masses[-1]:
        raise ValueError(
            f"mass {mass} outside the MESA grid {Path(path).stem} ({masses[0]}-{masses[-1]} Msun)"
        )
    lm = np.log(masses)
    im = int(_clipf(np.searchsorted(lm, math.log(mass)) - 1, 0, max(len(masses) - 2, 0)))
    pairs = [(im, 1.0)] if len(masses) == 1 else []
    if not pairs:
        wm = (math.log(mass) - lm[im]) / (lm[im + 1] - lm[im])
        pairs = [(im, 1.0 - wm), (im + 1, wm)]
    k2 = rg2 = 0.0
    for j, w in pairs:
        if w == 0.0:
            continue
        s = sel & (g["initial_mass"] == masses[j])
        o = np.argsort(g["age_myr"][s])
        age, kk, rr = g["age_myr"][s][o], g["apsidal_k2"][s][o], g["rg2"][s][o]
        if age_myr is None:  # zero-age main sequence: central hydrogen starts to fall
            h = g["center_h1"][s][o] if "center_h1" in g else None
            a0 = float(age[np.argmax(h < np.nanmax(h) - 0.003)]) if h is not None else float(age[0])
        else:
            a0 = float(age_myr)
        k2 += w * float(np.interp(a0, age, kk))
        rg2 += w * float(np.interp(a0, age, rr))
    return k2, rg2


def _mesa_zams_age(path: str, mass: float, z: float) -> float:
    """Age [Myr] of the zero-age main sequence of the nearest track of a MESA grid."""
    g = _mesa_grid(path)
    zs = np.unique(g["z"])
    zz = zs[int(np.argmin(np.abs(np.log(zs) - math.log(z))))]
    sel = g["z"] == zz
    masses = np.unique(g["initial_mass"][sel])
    mm = masses[int(np.argmin(np.abs(np.log(masses) - math.log(mass))))]
    s = sel & (g["initial_mass"] == mm)
    o = np.argsort(g["age_myr"][s])
    age = g["age_myr"][s][o]
    if "center_h1" not in g:
        return float(age[0])
    h = g["center_h1"][s][o]
    return float(age[np.argmax(h < np.nanmax(h) - 0.003)])


def structure_kind(sm) -> str:
    """Resolved stellar_model.structure: auto means Claret (2023) for sse, polytrope otherwise."""
    if sm.structure != "auto":
        return sm.structure
    return "claret2023" if sm.kind == "sse" else "polytrope"


def age_offset(cfg, name: str) -> float:
    """Effective age of star ``name`` minus stellar_model.age_myr [Myr] (stars.*.age_offset_myr)."""
    spec = cfg.stars.get(name)
    return float(spec.age_offset_myr) if spec is not None else 0.0


def star_properties(cfg, name: str, m: float) -> StarProperties:
    """Model values for star ``name`` with the per-star overrides of the configuration."""
    sm = cfg.stellar_model
    if sm.kind in ("sse", "mist") and age_offset(cfg, name):
        sm = sm.model_copy(update={"age_myr": sm.age_myr + age_offset(cfg, name)})

    def base() -> Path:  # only for files of the model (Path.cwd per call is costly)
        return Path(cfg.base_dir) if getattr(cfg, "base_dir", None) else Path.cwd()

    if sm.kind in ("zams", "sse"):
        p = zams_properties(m) if sm.kind == "zams" else sse_properties(m, sm.age_myr, sm.z)
        kind = structure_kind(sm)
        if sm.structure == "auto" and kind == "claret2023" and not claret_covers(m):
            kind = "polytrope"  # outside 0.2-8 Msun: polytrope, as for zams
        if kind == "claret2023":
            feh = math.log10(sm.z / SSE_Z_SUN) if sm.kind == "sse" else 0.0
            k2, rg2 = _structure("claret2023", m, p.radius_rsun, feh, None)
            p = replace(p, apsidal_k2=k2, rg2=rg2)
    elif sm.kind == "mist":
        p = mist_properties(m, sm.age_myr, sm.feh, sm.vvcrit)
    if sm.kind in ("zams", "sse", "mist") and str(sm.structure).startswith("mesa:"):
        z = (
            sm.z * Z_SUN_MESA / SSE_Z_SUN
            if sm.kind == "sse"
            else Z_SUN_MESA * 10.0**sm.feh
            if sm.kind == "mist"
            else Z_SUN_MESA
        )
        grid = find_mesa_grid(sm.structure[5:], base())
        if sm.kind == "zams":
            age = None
        elif sm.kind == "sse":  # SSE counts from the ZAMS, MESA from the pre-main sequence
            age = None if not sm.age_myr else sm.age_myr + _mesa_zams_age(str(grid), m, z)
        else:
            age = sm.age_myr
        k2, rg2 = mesa_structure(grid, m, age, z)
        p = replace(p, apsidal_k2=k2, rg2=rg2)
    if sm.kind == "table":
        p = table_properties(m, str((base() / sm.file).resolve()))
    elif sm.kind == "python":
        p = python_properties(m, sm.entry, base())
    spec = cfg.stars.get(name)
    if spec is None:
        return p
    radius = spec.R if spec.R is not None else p.radius_rsun * (1.0 + spec.radius_inflation)
    teff = p.teff_k
    lum = p.luminosity_lsun
    if spec.R is None and spec.radius_inflation and p.luminosity_lsun is not None:
        teff = TEFF_SUN * p.luminosity_lsun**0.25 / math.sqrt(radius)
    if spec.Teff is not None:
        # a given Teff (with the radius) fixes the luminosity, L = 4 pi R^2 sigma T^4
        teff = float(spec.Teff)
        lum = radius**2 * (teff / TEFF_SUN) ** 4
    elif spec.R is not None and p.luminosity_lsun is not None:
        # a given radius keeps the model luminosity: Teff follows from L = 4 pi R^2 sigma T^4
        teff = TEFF_SUN * p.luminosity_lsun**0.25 / math.sqrt(radius)
    return replace(
        p,
        radius_rsun=radius,
        luminosity_lsun=lum,
        apsidal_k2=spec.apsidal_k2 if spec.apsidal_k2 is not None else p.apsidal_k2,
        rg2=spec.rg2 if spec.rg2 is not None else p.rg2,
        lag_s=spec.lag_s if spec.lag_s is not None else p.lag_s,
        teff_k=teff,
    )


def roche_lobe_eggleton(q: float) -> float:
    """Roche-lobe radius over separation, Eggleton (1983) ApJ 268, 368; q = M_star / M_other."""
    q13 = q ** (1.0 / 3.0)
    return 0.49 * q13**2 / (0.6 * q13**2 + np.log1p(q13))
