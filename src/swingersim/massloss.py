"""Stellar mass loss in the orbital dynamics (physics.mass_loss).

Two sources of the stellar masses as functions of time:

``track``
    The masses follow the stellar tracks (``current_mass`` of MIST v1.2, which contains the
    wind prescriptions of the MESA models: Vink et al. 2001 for hot stars, Reimers 1975 and
    Bloecker 1995 for cool ones; Choi et al. 2016, sect. 3.6). Needs ``physics.evolution`` and
    ``stellar_model.kind: mist``. The SSE tracks of Hurley et al. (2000) in swingersim are computed
    without mass loss; for them swingersim adds the wind prescription of the SSE code (Hurley,
    Pols & Tout 2000, Sect. 7.1; routine mlwind of the SSE/BSE package): Nieuwenhuijzen & de
    Jager (1990) for L > 4000 Lsun, switched on linearly over 4000-4500 Lsun, and Reimers (1975)
    with eta = ``reimers_eta`` from the Hertzsprung gap on, the larger of the two; the mass of
    a star of initial mass M0 is the integral of that rate along its track.
``rates``
    Constant rates dM/dt [Msun/yr] given per star (negative = loss), for parametric studies and
    benchmarks; no stellar evolution is needed.

The wind is isotropic in the frame of the star and carries the star's momentum (the "Jeans
mode" of mass loss). The orbital response is

* N-body: the masses change during the integration (REBOUNDx operator ``modify_mass``,
  dm/dt = m / tau_mass), the velocities are unchanged. This is exact for any rate, from the
  adiabatic to the impulsive limit (Hadjidemetriou 1963; Veras et al. 2011).
* Secular: the orbit-averaged response of every orbit X with mass M_X (the total mass of the
  stars in its two components) to slow loss is da_X/dt = -a_X dM_X/dt / M_X with constant
  eccentricity vector and orientation (Hadjidemetriou 1963): a_X M_X is an adiabatic invariant.
  It holds when the mass-loss index Psi = |dM_X/dt| / (n_X M_X) << 1 (Veras et al. 2011,
  eq. 15); swingersim reports the largest Psi of every orbit.

Spins: the wind of a star of radius R spinning at Omega carries angular momentum at the rate
dJ/dt = (2/3) R^2 Omega dM/dt (a thin spherical shell; Hurley et al. 2002, eq. 11), so the spin
angular momentum J = r_g^2 M R^2 Omega scales as M^(2/(3 r_g^2)).

References: Jeans (1924), MNRAS 85, 2; Hadjidemetriou (1963), Icarus 2, 440; Blaauw (1961),
BAN 15, 265; Hills (1983), ApJ 267, 322; Veras et al. (2011), MNRAS 417, 2104; Perets & Kratter
(2012), ApJ 760, 99; Toonen et al. (2022), A&A 661, A61.
"""

from __future__ import annotations

import itertools
import math
from functools import lru_cache

import numpy as np

from . import units

YR_PER_MYR = 1.0e6


def enabled(cfg) -> bool:
    """The stellar masses change during the run: winds (physics.mass_loss) or mass transfer
    (physics.mass_transfer, swingersim.masstransfer)."""
    return (
        cfg.physics.mass_loss is not None or getattr(cfg.physics, "mass_transfer", None) is not None
    )


def wind(cfg) -> bool:
    return cfg.physics.mass_loss is not None


# ------------------------------------------------------------------ SSE winds
def sse_wind_rate(m: float, lum: float, r: float, z: float, phase: int, eta: float = 0.5) -> float:
    """Wind mass-loss rate [Msun/yr, negative] of an SSE star of mass m, luminosity lum [Lsun],
    radius r [Rsun] and phase (1 MS, 2 HG, 3 GB): Hurley et al. (2000), Sect. 7.1, as coded in
    the routine mlwind of SSE/BSE - the larger of Nieuwenhuijzen & de Jager (1990),
    9.6e-15 x (Z/0.02)^(1/2) R^0.81 L^1.24 M^0.16 for L > 4000 with the ramp
    x = min(1, (L - 4000)/500) (no jump of the rate at 4000 Lsun), and Reimers (1975),
    4e-13 eta L R / M, from the Hertzsprung gap on (kw >= 2; the paper says "on the GB and
    beyond", the code applies it to the gap as well)."""
    rate = 0.0
    if lum > 4000.0:
        x = min(1.0, (lum - 4000.0) / 500.0)
        rate = 9.6e-15 * x * math.sqrt(z / 0.02) * r**0.81 * lum**1.24 * m**0.16
    if phase >= 2:
        rate = max(rate, 4.0e-13 * eta * lum * r / m)
    return -rate


@lru_cache(maxsize=512)
def _sse_wind_history(m0: float, z: float, eta: float):
    """(ages [Myr], masses) of an SSE star of initial mass m0 with the SSE winds, from the
    ZAMS to helium ignition."""
    from scipy.integrate import solve_ivp

    from .evolution import track

    tr = track(m0, z)
    t_end = tr.t_hei * (1 - 1e-9)

    def rhs(age, y):
        st = tr.at(age)
        return [
            sse_wind_rate(max(y[0], 1e-3), st.luminosity_lsun, st.radius_rsun, z, st.phase, eta)
            * YR_PER_MYR
        ]

    edges = sorted({0.0, tr.t_ms, tr.t_bgb, t_end})
    ages, ms = [0.0], [m0]
    for a, b in itertools.pairwise(edges):
        if b <= a:
            continue
        s = solve_ivp(
            rhs, (a, b), [ms[-1]], method="LSODA", rtol=1e-8, atol=1e-10, max_step=(b - a) / 200.0
        )
        ages += list(s.t[1:])
        ms += list(s.y[0, 1:])
    return np.array(ages), np.array(ms)


class SSEWindTrack:
    """SSE track with the current mass of its winds (same interface as the MIST tracks)."""

    def __init__(self, m0: float, z: float, eta: float):
        from .evolution import track

        self.track = track(m0, z)
        self.m0, self.z, self.eta = m0, z, eta

    def current_mass(self, age_myr: float) -> float:
        self.track.phase_of(age_myr)  # raises beyond helium ignition
        a, m = _sse_wind_history(self.m0, self.z, self.eta)
        return float(np.interp(age_myr, a, m))

    def mdot(self, age_myr: float) -> float:
        st = self.track.at(age_myr)
        return sse_wind_rate(
            self.current_mass(age_myr),
            st.luminosity_lsun,
            st.radius_rsun,
            self.z,
            st.phase,
            self.eta,
        )


def wind_track(cfg, m0: float):
    """Track with ``current_mass`` and ``mdot`` (MIST, or SSE with the SSE winds)."""
    from .stellar import stellar_track

    sm = cfg.stellar_model
    if sm.kind == "sse":
        return SSEWindTrack(round(float(m0), 10), float(sm.z), _eta(cfg))
    return stellar_track(sm, float(m0))


def _eta(cfg) -> float:
    ml = getattr(cfg.physics, "mass_loss", None)
    return float(getattr(ml, "reimers_eta", 0.5)) if ml is not None else 0.5


def _track_masses(cfg, masses0, t_yr: float, skip=()) -> np.ndarray:
    from .stellar import age_offset

    age = float(cfg.stellar_model.age_myr) + t_yr / YR_PER_MYR
    out = []
    for k, m in enumerate(masses0):
        if k in skip:
            out.append(float(m))
            continue
        tr = wind_track(cfg, float(m))
        out.append(float(tr.current_mass(age + age_offset(cfg, f"m{k + 1}"))))
    return np.array(out)


def initial_masses(cfg, masses_cfg) -> np.ndarray:
    """Masses at t = 0. With source track the configured masses are the *initial* (ZAMS)
    masses of the tracks; the dynamics starts with the current masses at stellar_model.age_myr."""
    m = np.asarray(masses_cfg, float)
    if wind(cfg) and cfg.physics.mass_loss.source == "track":
        return _track_masses(cfg, m, 0.0)
    return m.copy()


def wind_masses_at(cfg, masses_cfg, t_yr: float, skip=()) -> np.ndarray:
    """Masses [Msun] with the winds alone (physics.mass_loss) at run time t_yr (stars in
    ``skip`` keep their configured masses)."""
    m = np.asarray(masses_cfg, float)
    ml = cfg.physics.mass_loss
    if ml is None:
        return m.copy()
    if ml.source == "track":
        return _track_masses(cfg, m, t_yr, skip)
    rates = np.array([ml.rates.get(f"m{k + 1}", 0.0) for k in range(len(m))])
    out = m + rates * t_yr
    if np.any(out <= 0.0):
        bad = [f"m{k + 1}" for k in np.flatnonzero(out <= 0.0)]
        raise ValueError(f"physics.mass_loss: stars {bad} have no mass left at t = {t_yr:g} yr")
    return out


def masses_at(cfg, masses_cfg, t_yr: float) -> np.ndarray:
    """Stellar masses [Msun] at run time t_yr; masses_cfg are the configured masses. Winds
    (physics.mass_loss) act on every star; the stars of binaries with mass transfer
    (physics.mass_transfer) follow their transfer histories, which include their winds."""
    m = np.asarray(masses_cfg, float)
    if not enabled(cfg):
        return m.copy()
    hist = {}
    if getattr(cfg.physics, "mass_transfer", None) is not None:
        from . import masstransfer

        hist = masstransfer.histories(cfg)
    skip = {k for h in hist.values() for k in h.pair}
    out = wind_masses_at(cfg, m, t_yr, skip)
    for h in hist.values():
        out[list(h.pair)] = h.masses_at(t_yr)
    return out


def mass_rates(cfg, masses_cfg, t_yr: float, dt_yr: float | None = None) -> np.ndarray:
    """dM/dt [Msun/yr] of every star at run time t_yr (central difference over dt_yr for
    tracks)."""
    m = np.asarray(masses_cfg, float)
    if not enabled(cfg):
        return np.zeros(len(m))
    ml = cfg.physics.mass_loss
    if ml is not None and ml.source == "rates":
        return np.array([ml.rates.get(f"m{k + 1}", 0.0) for k in range(len(m))])
    h = dt_yr or 1.0e3
    # central difference; one-sided at the ends of the tracks (birth, end of the grid)
    for a, b in ((t_yr - h, t_yr + h), (t_yr - h, t_yr), (t_yr, t_yr + h)):
        try:
            return (masses_at(cfg, m, b) - masses_at(cfg, m, a)) / (b - a)
        except ValueError:  # includes EvolvedBeyondGrid and ages before the track start
            continue
    return np.zeros(len(m))


def mass_series(cfg, masses_cfg, times) -> np.ndarray:
    """Masses of all stars at the run times (shape n_star x n_t)."""
    return np.array([masses_at(cfg, masses_cfg, float(t)) for t in np.asarray(times)]).T


def psi(mdot: float, m_total: float, a_au: float) -> float:
    """Mass-loss index Psi = |dM/dt| / (n M) of an orbit (Veras et al. 2011, eq. 15)."""
    n = math.sqrt(units.G * m_total / a_au**3)
    return abs(mdot) / (n * m_total)


def adiabatic_a(a0: float, m0: float, m: float) -> float:
    """Semi-major axis after slow isotropic mass loss, a M = const (Hadjidemetriou 1963)."""
    return a0 * m0 / m


def veras_adiabatic_e(e0: float, psi0: float, f: float) -> float:
    """Osculating eccentricity within one orbit during slow isotropic mass loss, to first order
    in Psi: the integral of de/dt = (alpha/mu)(e + cos f) (Veras et al. 2011, eq. 4, alpha =
    -dmu/dt > 0) over the true anomaly f from the pericentre,

        e = e0 + Psi0 (1 - e0^2)^(3/2) sin f / (1 + e0 cos f).

    (Veras et al. print the denominator of their eq. 17 as 1 - e0 cos f; the integral of their
    eq. 4 gives 1 + e0 cos f, and the N-body runs follow the latter - see validation.)
    Reference implementation (tests/validation)."""
    return e0 + psi0 * (1.0 - e0**2) ** 1.5 * math.sin(f) / (1.0 + e0 * math.cos(f))


def impulsive_circular(delta: float) -> tuple[float, float]:
    """Orbit after the instantaneous loss of a fraction delta of the total mass of a circular
    binary (Blaauw 1961; Hills 1983): a'/a = (1 - delta)/(1 - 2 delta), e' = delta/(1 - delta);
    unbound for delta >= 1/2. Reference implementation (tests/validation)."""
    if delta >= 0.5:
        return float("inf"), float("inf")
    return (1.0 - delta) / (1.0 - 2.0 * delta), delta / (1.0 - delta)


def wind_spin_factor(m_old: float, m_new: float, rg2: float) -> float:
    """Factor on the spin angular momentum after the wind took the mass m_old - m_new
    (Hurley et al. 2002, eq. 11: dJ = (2/3) R^2 Omega dM)."""
    if rg2 <= 0 or m_old <= 0:
        return 1.0
    return (m_new / m_old) ** (2.0 / (3.0 * rg2))


def orbit_masses(root, masses) -> dict[str, float]:
    from .hierarchy import mass_of

    return {n.name: mass_of(n.left, masses) + mass_of(n.right, masses) for n in root.nodes()}


def psi_summary(cfg, root, masses_cfg, a_of: dict[str, float], times) -> dict[str, float]:
    """Largest mass-loss index of every orbit over the run (orbit -> Psi_max), evaluated with
    the semi-major axes a_of (the initial ones scaled adiabatically)."""
    t = np.asarray(times, float)
    sample = t[:: max(1, len(t) // 50)]
    m_cfg = np.asarray(masses_cfg, float)
    m0 = orbit_masses(root, initial_masses(cfg, m_cfg))
    dt = abs(t[-1] - t[0]) / 1000.0 if len(t) > 1 else None
    best = dict.fromkeys(a_of, 0.0)
    for tt in sample:  # masses and rates once per sample, for all orbits
        M_all = orbit_masses(root, masses_at(cfg, m_cfg, float(tt)))
        Md_all = orbit_masses(root, mass_rates(cfg, m_cfg, float(tt), dt))
        for name, a0 in a_of.items():
            M, Md = M_all[name], Md_all[name]
            best[name] = max(best[name], psi(Md, M, adiabatic_a(a0, m0[name], M)))
    return {name: float(v) for name, v in best.items()}


def check_psi(cfg, root, masses_cfg, a_of: dict[str, float], times) -> dict[str, float]:
    """psi_summary, refusing orbit-averaged runs beyond physics.mass_loss.max_psi."""
    out = psi_summary(cfg, root, masses_cfg, a_of, times)
    if not wind(cfg):  # mass transfer: the binaries follow their histories, not a M = const
        return out
    lim = cfg.physics.mass_loss.max_psi
    bad = {n: v for n, v in out.items() if v > lim}
    if bad:
        txt = ", ".join(f"{n}: {v:.3g}" for n, v in bad.items())
        raise ValueError(
            f"physics.mass_loss: mass-loss index Psi = |dM/dt|/(n M) of orbits {txt} exceeds "
            f"max_psi = {lim}; the orbit-averaged (adiabatic) response is not valid there "
            "(Veras et al. 2011) - use run.backend: nbody"
        )
    return out


@lru_cache(maxsize=256)
def _initial_mass_cached(kind, feh, vvcrit, age, branch, m_now, z=0.02, eta=0.5):
    from types import SimpleNamespace

    sm = SimpleNamespace(kind=kind, feh=feh, vvcrit=vvcrit, age_myr=age, z=z)
    cfg = SimpleNamespace(
        stellar_model=sm,
        physics=SimpleNamespace(mass_loss=SimpleNamespace(branch=branch, reimers_eta=eta)),
    )
    return _initial_mass(cfg, m_now)


def initial_mass_for(cfg, m_now: float, offset_myr: float = 0.0) -> float:
    sm = cfg.stellar_model
    return _initial_mass_cached(
        sm.kind,
        float(sm.feh),
        float(sm.vvcrit),
        float(sm.age_myr) + float(offset_myr),
        cfg.physics.mass_loss.branch,
        round(float(m_now), 10),
        float(sm.z),
        _eta(cfg),
    )


def _initial_mass(cfg, m_now: float, n_grid: int = 400) -> float:
    """Initial mass of the MIST track whose current mass at stellar_model.age_myr is m_now.
    The current mass at a fixed age need not be monotonic in the initial mass (very massive
    stars lose most of their mass); all solutions are searched on a grid up to 300 Msun and an
    ambiguous case is refused."""
    from scipy.optimize import brentq

    from .evolution import EvolvedBeyondGrid

    sm = cfg.stellar_model
    age = float(sm.age_myr)

    def f(m0: float) -> float:
        try:
            return wind_track(cfg, float(m0)).current_mass(age) - m_now
        except EvolvedBeyondGrid:
            return float("nan")

    lo = max(0.1, m_now * 0.999)
    grid = np.geomspace(lo, 300.0 if sm.kind == "mist" else 8.0, n_grid)
    vals = np.array([f(m) for m in grid])
    roots = []
    for k in range(len(grid) - 1):
        a, b = vals[k], vals[k + 1]
        if np.isfinite(a) and np.isfinite(b) and a * b <= 0.0:
            if a == 0.0:
                roots.append(float(grid[k]))
            else:
                roots.append(float(brentq(f, grid[k], grid[k + 1], xtol=1e-8 * grid[k])))
    roots = sorted({round(r, 6) for r in roots})
    if not roots:
        raise ValueError(
            f"no MIST track ({sm.feh=}, vvcrit {sm.vvcrit}) has {m_now} Msun at {age} Myr "
            "(the star would be beyond the end of its track or the mass is out of range)"
        )
    if len(roots) > 1:
        msg = f"{m_now} Msun at {age} Myr is reached by tracks of initial masses {roots}"
        if cfg.physics.mass_loss.branch == "unique":
            raise ValueError(msg + "; give the initial masses (physics.mass_loss.masses: initial)")
        import warnings

        warnings.warn(msg + f"; using the least evolved star ({roots[0]} Msun)", stacklevel=2)
    return roots[0]


def resolve_masses(cfg):
    """Configuration with initial masses: with physics.mass_loss {source: track, masses:
    current} the configured present masses are replaced by the initial masses of their tracks
    (stored in cfg.masses, so the stars keep the track labels the rest of swingersim expects)."""
    ml = cfg.physics.mass_loss
    if ml is None or ml.source != "track" or ml.masses != "current":
        return cfg
    new = cfg.model_copy(deep=True)
    from .stellar import age_offset

    new.masses = {
        k: initial_mass_for(cfg, float(v), age_offset(cfg, k)) for k, v in cfg.masses.items()
    }
    new.physics.mass_loss = ml.model_copy(update={"masses": "initial"})
    if getattr(cfg, "base_dir", None):
        new._base_dir = cfg.base_dir
    return new
