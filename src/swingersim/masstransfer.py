"""Stable mass transfer by Roche-lobe overflow in two-star binaries (physics.mass_transfer).

Scope (F4 d): a star of a two-star binary that fills its Roche lobe transfers mass to its
companion. The binary is treated as circular (tides circularise close binaries before the
overflow; Hurley, Tout & Pols 2002, Sect. 2.6) and its evolution is computed before the
dynamics, from the stellar tracks and the initial separation. The orbits then follow that
history: the masses of the two stars, the separation of their orbit (scaled at every
operator-splitting step, keeping the eccentricity vector and the orientation), and every outer
orbit responds adiabatically to the mass that leaves the binary (a M = const, as for a wind).

Rates and stability follow the binary-star evolution algorithm (BSE) of Hurley et al. (2002):

* the donor of radius R overfilling its lobe R_L (Eggleton 1983) transfers
  3e-6 min(M, 5)^2 [ln(R/R_L)]^3 Msun/yr (eqs. 58-59), so it stays just above its lobe while
  its expansion or the shrinking lobe drives the transfer; the rate is limited by the thermal
  rate M / tau_KH of giant-like donors (Hertzsprung gap and later, eqs. 60-61) or by the
  dynamical rate (eqs. 62-63);
* transfer is dynamically unstable, ending in a common envelope (not followed by swingersim: the
  run stops), for q = M_donor / M_accretor above q_crit: giants eq. 57, Hertzsprung-gap donors
  4, main-sequence stars below 0.7 Msun 0.695;
* a main-sequence, Hertzsprung-gap or core-helium-burning accretor takes at most
  10 tau_Mdot / tau_KH,2 of the transferred mass (eqs. 64-65); the rest leaves the binary with
  the specific orbital angular momentum of the donor; giants take everything.

Stars that lose or gain mass are re-modelled as in BSE (Hurley et al. 2000, Sect. 7.1; 2002,
Sect. 2.6.6): a main-sequence star, or a Hertzsprung-gap star whose new mass would not shrink its
helium core, follows the track of its present mass at the same evolutionary state (equivalent
evolutionary point or fraction of its phase), so a donor shrinks with its mass and an accretor
brightens and is rejuvenated; giants keep the track of the mass they had when they left the gap,
and their radius follows the envelope mass with zeta_eq (eq. 56). The core mass never decreases.

Orbit response. Of the mass |dM_d| lost by the donor a fraction ``alpha`` leaves from the
vicinity of the donor (specific orbital angular momentum of the donor), a fraction ``beta``
from the vicinity of the accretor (isotropic re-emission) and 1 - alpha - beta is accreted.
With J = M_d M_a sqrt(G a / M),

    dln J / dM_d = (alpha M_a^2 + beta M_d^2) / (M M_d M_a),
    dln a / dM_d = 2 dln J / dM_d + (alpha + beta) / M - 2 / M_d + 2 (1 - alpha - beta) / M_a,

which gives a M_d^2 M_a^2 = const for conservative transfer and a M = const for alpha = 1
(Hurley et al. 2002, Sect. 2.6.1; Soberman, Phinney & van den Heuvel 1997).

Events that end the history (and the run): ``common_envelope``, ``contact`` (the accretor
fills its own lobe), ``donor_stripped`` (the donor is down to its helium core; the single-star
track cannot describe it), ``track_end``.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field, replace
from functools import lru_cache

import numpy as np

from . import units
from .commonenv import X_GIANT
from .evolution import EvolvedBeyondGrid
from .massloss import YR_PER_MYR, sse_wind_rate
from .stellar import roche_lobe_eggleton, stellar_track


def tau_kh_yr(m: float, r: float, lum: float, m_core: float = 0.0, giant: bool = False) -> float:
    """Kelvin-Helmholtz time scale [yr] of Hurley et al. (2002), eq. 61:
    1e7 M M' / (R L) with M' = M (main sequence) or M - M_c (giant-like stars)."""
    return 1.0e7 * m * ((m - m_core) if giant else m) / (r * lum)


def tau_dyn_yr(m: float, r: float) -> float:
    """Dynamical time scale [yr], Hurley et al. (2002) eq. 63."""
    return 5.05e-5 * math.sqrt(r**3 / m)


def bse_rate(m: float, r: float, rl: float) -> float:
    """Transfer rate [Msun/yr] of an overfilling donor, Hurley et al. (2002) eqs. 58-59:
    3e-6 min(M, 5)^2 [ln(R / R_L)]^3; 0 when R <= R_L."""
    if r <= rl:
        return 0.0
    return 3.0e-6 * min(m, 5.0) ** 2 * math.log(r / rl) ** 3


# critical mass ratio of main-sequence (above 0.7 Msun) and core-helium-burning donors: their
# adiabatic exponent zeta_ad = 4.72 equated with the lobe exponent zeta_L = 2.13 q - 1.67 of
# conservative transfer (Tout et al. 1997) gives q = 3.0 (Hurley et al. 2002, BSE; e.g. Iorio
# et al. 2023 / SEVN, "HW" criterion)
Q_CRIT_MS = 3.0


def q_crit(phase: int, m: float, m_core: float, q_ms: float = Q_CRIT_MS) -> float | None:
    """Critical mass ratio q = M_donor / M_accretor above which transfer is dynamical (Hurley et
    al. 2002, Sect. 2.6): giants (giant branch, AGB) eq. 57 with x = 0.3; Hertzsprung-gap stars 4;
    main-sequence stars below 0.7 Msun (deeply convective) 0.695; main sequence above 0.7 Msun and
    core helium burning 3.0 (zeta_ad = 4.72; q_ms overrides it: physics.mass_transfer.q_crit_ms,
    e.g. 1.6 of binary_c after de Mink et al. 2007, where case-A contact is taken as a merger)."""
    if phase in (3, 5, 6):
        return (1.67 - X_GIANT + 2.0 * (m_core / m) ** 5) / 2.13
    if phase == 2:
        return 4.0
    if phase in (0, 1) and m < 0.7:
        return 0.695
    if phase in (1, 4):
        return q_ms
    return None  # pre-main sequence, remnants: not limited here


def dln_a_dmd(md: float, ma: float, alpha: float, beta: float) -> float:
    """dln a / dM_d of a circular binary during mass transfer (see the module docstring)."""
    m = md + ma
    dlnj = (alpha * ma * ma + beta * md * md) / (m * md * ma)
    return 2.0 * dlnj + (alpha + beta) / m - 2.0 / md + 2.0 * (1.0 - alpha - beta) / ma


def zeta_lobe(md: float, ma: float, alpha: float, beta: float, h: float = 1e-6) -> float:
    """zeta_L = dln R_L / dln M_d of the donor's Roche lobe (Eggleton 1983) during transfer."""

    def ln_f(x):
        q = x / (ma - (1.0 - alpha - beta) * (x - md))
        return math.log(roche_lobe_eggleton(q))

    dlnf = (ln_f(md * (1 + h)) - ln_f(md * (1 - h))) / (2.0 * h * md)
    return md * (dln_a_dmd(md, ma, alpha, beta) + dlnf)


def orbit_path(md0, ma0, a0, alpha, beta, md_min, n=2000):
    """Separation and accretor mass as functions of the donor mass from md0 down to md_min:
    arrays (md, ma, a), md decreasing. Reference implementation (tests/validation)."""
    from scipy.integrate import solve_ivp

    def rhs(md, y):
        ma, _ = y
        return [-(1.0 - alpha - beta), dln_a_dmd(md, ma, alpha, beta)]

    grid = np.linspace(md0, md_min, n)
    sol = solve_ivp(rhs, (md0, md_min), [ma0, math.log(a0)], t_eval=grid, rtol=1e-11, atol=1e-13)
    return sol.t, sol.y[0], np.exp(sol.y[1])


class _Track:
    """Methods shared by the histories (History, TertiaryHistory): the arrays t, a, m, r_corr,
    age, m_ev, m_core (one row per star of ``pair``) and the common-envelope remnants."""

    def a_at(self, t: float) -> float:
        return float(np.interp(t, self.t, self.a))

    def state_at(self, sm, i: int, t: float):
        """Re-modelled state of star i (index within ``pair``) at run time t; the core left by
        a common envelope is a dark point mass of its core radius."""
        if i in self.remnants and self.ce is not None and t >= self.ce["t_yr"]:
            mc, rc = self.remnants[i]
            k = int(np.searchsorted(self.t, self.ce["t_yr"]))
            base = remodelled_state(sm, float(self.m_ev[i][k]), float(self.age[i][k]))
            return replace(
                base, radius_rsun=rc, luminosity_lsun=1e-4, m_core=mc, m_env=0.0, r_env=0.0
            )
        m_ev = float(np.interp(t, self.t, self.m_ev[i]))
        age = float(np.interp(t, self.t, self.age[i]))
        st = remodelled_state(sm, m_ev, age)
        corr = math.exp(float(np.interp(t, self.t, self.r_corr[i])))
        core = float(np.interp(t, self.t, self.m_core[i]))
        return replace(st, radius_rsun=st.radius_rsun * corr, m_core=max(st.m_core, core))

    def regime_at(self, t: float) -> str | None:
        """Regime of the history at run time t (None before its first recorded regime)."""
        out = None
        for t_p, r in getattr(self, "phases", ()):
            if t_p > t:
                break
            out = r
        return out

    def transferring(self, t: float) -> bool:
        if self.t_onset is None or t < self.t_onset:
            return False
        if self.regime_at(t) == "detached":  # the binary has detached after an episode
            return False
        return self.t_event is None or t <= self.t_event

    def fill_at(self, t: float) -> float:
        """Largest Roche-lobe filling R / R_L of the (potential) donors at run time t in the
        circular history (0 when not recorded)."""
        if self.fill is None:
            return 0.0
        return float(np.interp(t, self.t[: len(self.fill)], self.fill))

    def regimes_text(self) -> str:
        """Transfer regimes and the run times [yr] at which they begin."""
        return " ".join(f"{r}@{t:.6g}" for t, r in self.phases)

    def episode_summary(self) -> dict:
        """Regimes, peak rate and donor core mass at the onset (history summaries)."""
        return {
            "regimes": self.regimes_text(),
            "peak_rate_msun_yr": self.peak_rate,
            "t_peak_rate_yr": self.t_peak_rate,
            "peak_regime": self.peak_regime,
            "donor_core_onset_msun": self.donor_core_onset,
        }


@dataclass
class History(_Track):
    """Mass-transfer history of one binary (donor index d, accretor index a)."""

    pair: tuple[int, int]
    t: np.ndarray  # run time [yr]
    m: np.ndarray  # masses of the two stars of `pair`, shape (2, n_t)
    a: np.ndarray  # separation [AU]
    donor: int | None = None
    t_onset: float | None = None
    t_reversal: float | None = None  # q reversal (zeta_L = 0)
    event: str | None = None
    t_event: float | None = None
    phases: list[tuple[float, str]] = field(default_factory=list)  # (t, regime)
    r_corr: np.ndarray | None = None  # ln of the radius corrections of giants (eq. 56)
    age: np.ndarray | None = None  # ages [Myr] of the two stars on their tracks (m_ev)
    m_ev: np.ndarray | None = None  # masses of the tracks the two stars follow
    m_core: np.ndarray | None = None  # helium core masses [Msun] (never decrease)
    ce: dict | None = None  # outcome of a common envelope (physics.mass_transfer.common_envelope)
    remnants: dict = field(default_factory=dict)  # index in pair -> (mass, radius) after a CE
    fill: np.ndarray | None = None  # largest R / R_L at the history points (up to the event)
    peak_rate: float | None = None  # largest transfer rate [Msun/yr] at the history points
    t_peak_rate: float | None = None
    peak_regime: str | None = None
    donor_core_onset: float | None = None  # helium core mass of the donor at the onset [Msun]
    n_zero_rhs: int = 0  # evaluations beyond the track grid (derivatives set to zero)

    def masses_at(self, t: float) -> tuple[float, float]:
        return (
            float(np.interp(t, self.t, self.m[0])),
            float(np.interp(t, self.t, self.m[1])),
        )


def _limit_rate(st, m: float, rate: float) -> tuple[float, str]:
    """Rate of eq. 58 limited by the thermal rate (giant-like donors: Hertzsprung gap, giants,
    core helium burning, AGB) or the dynamical rate (other stars), Hurley et al. (2002)
    Sects. 2.6.2-2.6.3; returns (rate, regime)."""
    if st.phase >= 2:
        lim = m / tau_kh_yr(m, st.radius_rsun, st.luminosity_lsun, min(st.m_core, 0.99 * m), True)
        return (lim, "thermal") if rate > lim else (rate, "nuclear")
    lim = m / tau_dyn_yr(m, st.radius_rsun)
    return (lim, "dynamical") if rate > lim else (rate, "nuclear")


def _accreted_fraction(sa, ma: float, rate: float) -> float:
    """Fraction of the transferred mass the accretor takes, Hurley et al. (2002) eqs. 64-65:
    main-sequence, Hertzsprung-gap and core-helium-burning accretors take at most
    10 tau_Mdot / tau_KH,2 (tau_Mdot = M_2 / Mdot); giants take everything."""
    if sa.phase in (3, 5, 6) or rate <= 0.0:
        return 1.0
    tau_m = ma / rate
    tau2 = tau_kh_yr(
        ma, sa.radius_rsun, sa.luminosity_lsun, min(sa.m_core, 0.99 * ma), sa.phase >= 2
    )
    return min(10.0 * tau_m / tau2, 1.0)


# ------------------------------------------------------------------ evolutionary coordinate
CORE_RELAX_YR = 1.0e3  # yr: the running maximum of the core mass follows its growth
M_STEP = 0.005  # Msun: spacing of the tracks between which a re-modelled star is interpolated


def _nodes(sm, tr):
    """(coordinates, ages [Myr]) at which the age of the track is linear in the evolutionary
    coordinate: the equivalent evolutionary points (MIST) or the phase boundaries (SSE)."""
    nodes = tr.__dict__.get("_mt_nodes")
    if nodes is None:
        if sm.kind == "mist":
            nodes = (np.arange(1.0, tr.n_eep + 1.0), np.asarray(tr._age, float))
        else:
            nodes = (np.array([0.0, 1.0, 2.0, 3.0]), np.array([0.0, tr.t_ms, tr.t_bgb, tr.t_hei]))
        tr.__dict__["_mt_nodes"] = nodes
    return nodes


def age_at_coord(sm, tr, u: float) -> float:
    """Age [Myr] of the track tr at the evolutionary coordinate u: the equivalent evolutionary
    point (MIST) or phase - 1 + fraction of the phase (SSE: 0-1 main sequence, 1-2 Hertzsprung
    gap, 2-3 giant branch)."""
    us, a = _nodes(sm, tr)
    return float(np.interp(u, us, a))


def coord_end(sm, tr) -> float:
    """Last evolutionary coordinate the track describes."""
    return float(tr.n_eep) if sm.kind == "mist" else 3.0


def _state_on(sm, tr, u: float):
    if sm.kind == "mist":
        return tr.at_eep(u, age_at_coord(sm, tr, u))
    return tr.at(age_at_coord(sm, tr, u))


# grid tracks of the re-modelled stars (key: stellar model, index of the lower track). Held
# here rather than in the LRU cache of stellar_track, whose 256 entries other users of the
# tracks (fits, wind tracks of other masses) may evict while a history is integrated.
_GRID_TRACKS: dict = {}


def _grid_tracks(sm, m_ev: float):
    """The two grid tracks bracketing m_ev and the interpolation weight. At the top of the
    track grid (no upper neighbour) both tracks are the lower one, with weight 0."""
    m_ev = max(float(m_ev), 0.1)  # trial steps of the integrator may probe nonsense
    i = math.floor(m_ev / M_STEP + 1e-6)
    lo = i * M_STEP
    w = min(max((m_ev - lo) / M_STEP, 0.0), 1.0)
    key = (sm.kind, getattr(sm, "feh", None), getattr(sm, "vvcrit", None), getattr(sm, "z", None))
    key += (i,)
    pair = _GRID_TRACKS.get(key)
    if pair is None:
        t0 = stellar_track(sm, round(lo, 6))
        try:
            t1 = stellar_track(sm, round(lo + M_STEP, 6))
        except ValueError as exc:
            if w > 1e-9:  # the star has grown beyond the track grid: its evolution ends here
                raise EvolvedBeyondGrid(
                    f"mass {m_ev:.4f} Msun beyond the top of the stellar track grid"
                ) from exc
            t1 = t0  # m_ev at the top of the grid
        if len(_GRID_TRACKS) > 512:  # ~0.1 MB per MIST track
            _GRID_TRACKS.clear()
            _COORDS.clear()  # keyed by the ids of the tracks
            _FRACTIONS.clear()
        pair = _GRID_TRACKS[key] = (t0, t1)
    if pair[1] is pair[0]:
        w = 0.0
    return pair[0], pair[1], w


def track_end_age(sm, tr) -> float:
    return float(tr.t_end) if sm.kind == "mist" else float(tr.t_hei)


def virtual_end_age(sm, m_ev: float) -> float:
    """End age [Myr] of the track of mass m_ev (interpolated between the grid tracks): the end
    of the evolutionary coordinates both tracks describe (MIST tracks end at different EEPs;
    the interpolated star cannot go beyond the shorter one). 0 beyond the top of the grid."""
    try:
        t0, t1, w = _grid_tracks(sm, m_ev)
    except EvolvedBeyondGrid:
        return 0.0
    if sm.kind == "mist" and w >= 1e-9:
        u0, a0 = _nodes(sm, t0)
        u1, a1 = _nodes(sm, t1)
        n = min(len(u0), len(u1))
        return float((1 - w) * a0[n - 1] + w * a1[n - 1])
    return (1 - w) * track_end_age(sm, t0) + w * track_end_age(sm, t1)


_COORDS: dict = {}  # recent virtual coordinates (the RHS asks for each star's twice)


def _coord_between(sm, t0, t1, w: float, age_myr: float) -> float:
    key = (id(t0), id(t1), w, age_myr)
    u = _COORDS.get(key)
    if u is not None:
        return u
    u0, a0 = _nodes(sm, t0)
    if w < 1e-9:
        u = float(np.interp(age_myr, a0, u0))
    else:
        u1, a1 = _nodes(sm, t1)
        n = min(len(u0), len(u1))
        u = float(np.interp(age_myr, (1 - w) * a0[:n] + w * a1[:n], u0[:n]))
    if len(_COORDS) > 32:
        _COORDS.clear()
    _COORDS[key] = u
    return u


def remodelled_state(sm, m_ev: float, age_myr: float):
    """State of a star that follows the track of mass m_ev and is at age_myr on it; between
    the grid tracks (spacing M_STEP) the states at the same evolutionary coordinate are
    interpolated (log R, log L, core mass), so the state is continuous in m_ev."""
    t0, t1, w = _grid_tracks(sm, m_ev)
    u = _coord_between(sm, t0, t1, w, age_myr)
    s0 = _state_on(sm, t0, min(u, coord_end(sm, t0) - 1e-9))
    if w < 1e-9:
        return replace(s0, age_myr=float(age_myr))
    s1 = _state_on(sm, t1, min(u, coord_end(sm, t1) - 1e-9))
    return replace(
        s0,
        age_myr=float(age_myr),
        radius_rsun=s0.radius_rsun ** (1 - w) * s1.radius_rsun**w,
        luminosity_lsun=s0.luminosity_lsun ** (1 - w) * s1.luminosity_lsun**w,
        m_core=(1 - w) * s0.m_core + w * s1.m_core,
        m_env=(1 - w) * s0.m_env + w * s1.m_env,
        r_env=(1 - w) * s0.r_env + w * s1.r_env,
    )


def _core_at_gb(sm, tr) -> float:
    c = tr.__dict__.get("_mt_core_gb")
    if c is None:
        if sm.kind == "mist":
            from .mist import EEP_RGB_TIP, EEP_TAMS

            u = getattr(tr, "eep_gb", 0.5 * (EEP_TAMS + EEP_RGB_TIP))
        else:
            u = 2.0
        c = tr.__dict__["_mt_core_gb"] = _state_on(sm, tr, min(u, coord_end(sm, tr) - 1e-9)).m_core
    return c


def core_at_giant_branch(sm, m_ev: float) -> float:
    """Helium core mass [Msun] at the base of the giant branch of the track of mass m_ev."""
    t0, t1, w = _grid_tracks(sm, m_ev)
    return (1 - w) * _core_at_gb(sm, t0) + w * _core_at_gb(sm, t1)


def age_shift_per_mass(sm, m_ev: float, age_myr: float) -> float:
    """d(age)/d(m_ev) [Myr/Msun] at fixed evolutionary coordinate: how the age on the track
    must change when the star is re-modelled for a different mass (fractional age kept)."""
    t0, t1, w = _grid_tracks(sm, m_ev)
    u = min(
        _coord_between(sm, t0, t1, w, age_myr), coord_end(sm, t0) - 1e-9, coord_end(sm, t1) - 1e-9
    )
    return (age_at_coord(sm, t1, u) - age_at_coord(sm, t0, u)) / M_STEP


def _wind_columns(sm, tr):
    """(EEPs, current mass / initial mass, (dM/dt) / initial mass [1/yr]) of a MIST track,
    contiguous (cached on the track)."""
    cols = tr.__dict__.get("_mt_wind")
    if cols is None:
        x, _ = _nodes(sm, tr)
        cols = tr.__dict__["_mt_wind"] = (
            x,
            np.ascontiguousarray(tr._t[:, tr._cols["star_mass"]]) / tr.mass,
            np.ascontiguousarray(tr._t[:, tr._cols["star_mdot"]]) / tr.mass,
        )
    return cols


_FRACTIONS: dict = {}  # recent results (wind_rate and track_label ask for the same point)


def _mist_fractions(sm, m_ev: float, age_myr: float) -> tuple[float, float]:
    """(current mass / initial mass, (dM/dt) / initial mass [1/yr]) of the MIST wind of the
    track of mass m_ev at age_myr: the grid tracks at the same evolutionary coordinate,
    interpolated linearly in m_ev like the structure (remodelled_state). Raises
    EvolvedBeyondGrid beyond the end of that track."""
    t0, t1, w = _grid_tracks(sm, m_ev)
    key = (id(t0), id(t1), w, age_myr)
    out = _FRACTIONS.get(key)
    if out is not None:
        return out
    if age_myr > virtual_end_age(sm, m_ev):
        raise EvolvedBeyondGrid(f"track of {m_ev:.4f} Msun ends before {age_myr:.2f} Myr")
    u = _coord_between(sm, t0, t1, w, age_myr)
    x0, m0, d0 = _wind_columns(sm, t0)
    fm, fd = float(np.interp(u, x0, m0)), float(np.interp(u, x0, d0))
    if w >= 1e-9:
        x1, m1, d1 = _wind_columns(sm, t1)
        fm = (1 - w) * fm + w * float(np.interp(u, x1, m1))
        fd = (1 - w) * fd + w * float(np.interp(u, x1, d1))
    if len(_FRACTIONS) > 32:
        _FRACTIONS.clear()
    out = _FRACTIONS[key] = (fm, fd)
    return out


def wind_rate(cfg, k: int, m: float, st, age_myr: float, m_ev: float) -> float:
    """Wind rate dM/dt [Msun/yr, < 0 for a loss] of star k (index in the configuration) of mass m in
    the state st, on the track of mass m_ev at age_myr (physics.mass_loss; Jeans mode)."""
    ml = cfg.physics.mass_loss
    sm = cfg.stellar_model
    if ml is None:
        return 0.0
    if ml.source == "rates":
        return float(ml.rates.get(f"m{k + 1}", 0.0))
    if sm.kind == "sse":
        eta = float(getattr(ml, "reimers_eta", 0.5))
        return sse_wind_rate(m, st.luminosity_lsun, st.radius_rsun, float(sm.z), st.phase, eta)
    try:
        frac_m, frac_mdot = _mist_fractions(sm, m_ev, age_myr)
    except (EvolvedBeyondGrid, ValueError):
        return 0.0
    return m * frac_mdot / frac_m


def track_label(cfg, m: float, m_ev: float, age_myr: float) -> float:
    """Track mass of a star re-modelled at its present mass m: the mass itself, or with MIST
    winds the initial mass whose track has the present mass (the wind already lost)."""
    ml = cfg.physics.mass_loss
    if ml is None or ml.source != "track" or cfg.stellar_model.kind != "mist":
        return m
    try:
        frac_m, _ = _mist_fractions(cfg.stellar_model, m_ev, age_myr)
    except (EvolvedBeyondGrid, ValueError):
        return m
    return m / frac_m


# mass range of main-sequence stars with a convective core or fully convective (Hurley et al.
# 2000, Sect. 7.1): accreted hydrogen is mixed into the burning region and rejuvenates them
REJUV_LOW, REJUV_HIGH = 0.35, 1.25


def ms_fraction_age(sm, m_ev: float, age_myr: float) -> float:
    """Time spent on the main sequence [Myr] of the track of mass m_ev at age_myr: age - age at
    the ZAMS (MIST, equivalent evolutionary point 202) or the age itself (SSE)."""
    if sm.kind != "mist":
        return float(age_myr)
    from .mist import EEP_ZAMS

    t0, t1, w = _grid_tracks(sm, m_ev)
    z = (1 - w) * age_at_coord(sm, t0, EEP_ZAMS) + w * age_at_coord(sm, t1, EEP_ZAMS)
    return max(float(age_myr) - z, 0.0)


def star_response(cfg, st, dm: float, m: float, m_ev: float, age_myr: float, core: float, label):
    """Response of a star to a mass change dm [Msun/yr] (transfer, accretion and, for SSE
    tracks, its wind): (dln R_corr/dt, d age/dt [Myr/yr], d m_ev/dt).

    As in BSE (Hurley et al. 2000, Sect. 7.1; 2002, Sect. 2.6.6), a main-sequence star, or a
    Hertzsprung-gap star whose new mass would not shrink its helium core, follows the track of
    its present mass at the same fractional age (equivalent evolutionary point); a main-sequence
    accretor with a convective core (M > 1.25 Msun) or fully convective (M < 0.35 Msun) mixes the
    fresh hydrogen into its burning region and is rejuvenated further, t' = mu (t'_MS / t_MS) t
    with mu = M / M' (physics.mass_transfer.rejuvenation: hurley; fractional_age: mu = 1).
    Giants keep their track and respond with zeta_eq = -x + 2 (M_c / M)^5 (eq. 56)."""
    sm = cfg.stellar_model
    if dm == 0.0:
        return 0.0, 0.0, 0.0
    remodel = st.phase <= 1 or (
        st.phase == 2 and (dm >= 0.0 or core_at_giant_branch(sm, label) >= core)
    )
    if not remodel:
        zeta = -X_GIANT + 2.0 * (min(st.m_core, m) / m) ** 5
        return zeta * dm / m, 0.0, 0.0
    dage = age_shift_per_mass(sm, m_ev, age_myr) * dm
    mt = cfg.physics.mass_transfer
    rej = getattr(mt, "rejuvenation", "hurley") if mt is not None else "hurley"
    if rej == "hurley" and dm > 0.0 and st.phase == 1 and (m > REJUV_HIGH or m < REJUV_LOW):
        # fraction x of the main sequence -> x M / (M + dM): d age = -t_MS,done dM / M
        dage -= ms_fraction_age(sm, m_ev, age_myr) * dm / m
    # the track label follows the mass: m_ev = m / frac_m with MIST winds (track_label)
    return 0.0, dage, dm * track_label(cfg, m, m_ev, age_myr) / m


# ------------------------------------------------------------------ contact binaries
JDOT_GS08 = (
    4.9e41 / 1.41017e52
)  # Gazeas & Stepien (2008) eq. 19 in Msun AU^2 yr^-2 (R, M solar, P d)


def orbital_j(m1: float, m2: float, a: float) -> float:
    """Orbital angular momentum [Msun AU^2 / yr] of a circular binary."""
    return m1 * m2 * math.sqrt(units.G * a / (m1 + m2))


def contact_jdot(m1: float, m2: float, r1: float, r2: float, a: float) -> float:
    """Angular-momentum loss [Msun AU^2 / yr^2, < 0] of a contact binary by the magnetised wind
    of its common envelope, dJ/dt = -4.9e41 (R_1^2 M_1 + R_2^2 M_2) / P g cm^2 s^-1 per year
    (R, M solar, P in days; Gazeas & Stepien 2008, eq. 19; Stepien 2006)."""
    p_day = units.sma_to_period(a, m1 + m2) * units.DAYS_PER_YEAR
    return -JDOT_GS08 * (r1 * r1 * m1 + r2 * r2 * m2) / p_day


@lru_cache(maxsize=1)
def _outer_lobe_table() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """log q grid and the volume radii of both stars' parts of the outer critical surface."""
    from . import roche

    lq = np.linspace(math.log(0.01), math.log(100.0), 121)
    r = np.array(
        [
            [
                roche.binary_surface(
                    math.exp(x), roche.critical_potentials(math.exp(x))[1], k, 24, 48
                ).r_vol
                for x in lq
            ]
            for k in (0, 1)
        ]
    )
    return lq, r[0], r[1]


def outer_lobe_radius(q: float, k: int) -> float:
    """Volume-equivalent radius (units of a) of star k's part of the outer critical surface
    (L2/L3) of a binary with q = M_2 / M_1: a contact binary that overfills it loses mass through
    the outer Lagrangian point and merges (interpolated in a table of swingersim.roche volumes)."""
    lq, r0, r1 = _outer_lobe_table()
    return float(np.interp(math.log(min(max(q, 0.01), 100.0)), lq, r0 if k == 0 else r1))


# ------------------------------------------------------------------ shared machinery
@dataclass(frozen=True)
class _Layout:
    """Rows of the state vector of a history: the masses of the n stars come first, then
    (after the orbits) ln of their radius corrections, ages, track masses and core masses."""

    n: int
    rc: int
    age: int
    mev: int
    core: int

    @property
    def size(self) -> int:
        return self.core + self.n


_BINARY = _Layout(2, 3, 5, 7, 9)  # m1 m2 ln_a | rc age m_ev core (x2)
_TERTIARY = _Layout(3, 5, 8, 11, 14)  # m3 m1 m2 ln_a_out ln_a_in | rc age m_ev core (x3)


def _memoised(fn):
    """fn(t, y) memoised on the exact arguments (the RHS and every event function of one
    step evaluate the same point)."""
    memo: dict = {}

    def terms(t, y):
        key = (float(t), *(float(v) for v in y))
        out = memo.get(key)
        if out is None:
            if len(memo) > 64:
                memo.clear()
            out = memo[key] = fn(t, y)
        return out

    return terms


def _remodel(sm, y, lay: _Layout):
    """Re-modelled states of the stars of a history (core mass never below its running
    maximum, radius with the giant correction): (states, raw core masses, radii)."""
    st = [remodelled_state(sm, y[lay.mev + k], y[lay.age + k]) for k in range(lay.n)]
    raw_core = [s.m_core for s in st]
    rad = [s.radius_rsun * math.exp(y[lay.rc + k]) for k, s in enumerate(st)]
    st = [
        replace(s, m_core=max(s.m_core, y[lay.core + k]), radius_rsun=rad[k])
        for k, s in enumerate(st)
    ]
    return st, raw_core, rad


def _winds(cfg, stars, y, st, lay: _Layout) -> list[float]:
    """Wind rates [Msun/yr] of the stars of a history (physics.mass_loss; Jeans mode)."""
    return [
        wind_rate(cfg, stars[k], y[k], st[k], y[lay.age + k], y[lay.mev + k]) for k in range(lay.n)
    ]


def _star_derivatives(cfg, out, st, raw_core, dm, w, y, lay: _Layout) -> None:
    """Fill the derivatives of the radius corrections, ages, track masses and core masses in
    out. The MIST tracks contain their own winds; the SSE tracks do not, so an SSE star is
    re-modelled (or its giant envelope responds) for the wind as well, as in BSE. The core
    mass is a running maximum that relaxes upwards within CORE_RELAX_YR."""
    sse = cfg.stellar_model.kind == "sse"
    for k in range(lay.n):
        out[lay.core + k] = max(raw_core[k] - y[lay.core + k], 0.0) / CORE_RELAX_YR
    for k in range(lay.n):
        m, m_ev, age = y[k], y[lay.mev + k], y[lay.age + k]
        dm_star = dm[k] + (w[k] if sse else 0.0)
        out[lay.rc + k], dage, out[lay.mev + k] = star_response(
            cfg, st[k], dm_star, m, m_ev, age, y[lay.core + k], track_label(cfg, m, m_ev, age)
        )
        out[lay.age + k] = 1.0 / YR_PER_MYR + dage


def _start(cfg, stars, m_labels, lay: _Layout):
    """(masses at t = 0, labels, ages, initial core masses) of the stars of a history."""
    from .massloss import wind_masses_at
    from .stellar import age_offset

    sm = cfg.stellar_model
    labels = [float(m_labels[k]) for k in stars]
    ages = [float(sm.age_myr) + age_offset(cfg, f"m{k + 1}") for k in stars]
    c0 = [remodelled_state(sm, labels[k], ages[k]).m_core for k in range(lay.n)]
    m_start = list(labels)
    if cfg.physics.mass_loss is not None:  # winds before the run: the current masses
        m_now = wind_masses_at(cfg, np.asarray(m_labels, float), 0.0)
        m_start = [float(m_now[k]) for k in stars]
    return m_start, labels, ages, c0


def _integrate(rhs, t_end: float, y0, atol, events):
    from scipy.integrate import solve_ivp

    return solve_ivp(
        rhs,
        (0.0, t_end),
        y0,
        method="LSODA",
        rtol=1e-8,
        atol=atol,
        max_step=t_end / 2000.0,
        events=events,
    )


def _first_event(sol, names, ts, ys, skip=()):
    """Cut the solution at the earliest terminal event: (ts, ys, event, t_event)."""
    hits = [
        (float(te[0]), names[k], ye[0])
        for k, (te, ye) in enumerate(zip(sol.t_events, sol.y_events, strict=True))
        if len(te) and names[k] not in skip
    ]
    if not hits:
        return ts, ys, None, None
    te, ev, ye = min(hits, key=lambda x: x[0])
    keep = ts < te
    return np.append(ts[keep], te), np.column_stack([ys[:, keep], ye]), ev, te


def _check_solution(h, sol, where: str, n_zero: int) -> None:
    """Mark a history whose integration failed (event integration_failed: the run stops
    there instead of following a history padded to its end) and report RHS evaluations
    beyond the track grid (zero derivatives) that no track_end event explains."""
    h.n_zero_rhs = n_zero
    if n_zero and h.event != "track_end":
        warnings.warn(
            f"physics.mass_transfer: the history of {where} evaluated its derivatives beyond "
            f"the track grid {n_zero} times (set to zero there)",
            stacklevel=3,
        )
    if sol.status >= 0 or h.event is not None:
        return
    h.event, h.t_event = "integration_failed", float(sol.t[-1])
    warnings.warn(
        f"physics.mass_transfer: the history of {where} failed at t = {h.t_event:.6g} yr "
        f"({sol.message})",
        stacklevel=3,
    )


def _scan(h, ts, ys, point, reversed_, where: str) -> None:
    """Transfer episodes on the solution points: onset and donor, mass-ratio reversal,
    regimes, peak rate and the Roche-lobe filling. point(t, y) returns (donor star or None,
    regime, rate, filling, donor core mass, attributes set at the onset)."""
    prev = "detached"
    fills = []
    for t, y in zip(ts, ys.T, strict=True):
        try:
            donor, regime, rate, fill, core, extra = point(t, y)
        except EvolvedBeyondGrid:
            warnings.warn(
                f"physics.mass_transfer: the history of {where} left the track grid at "
                f"t = {t:.6g} yr without an event",
                stacklevel=3,
            )
            break
        fills.append(fill)
        if donor is not None and h.t_onset is None:
            h.t_onset, h.donor, h.donor_core_onset = float(t), donor, float(core)
            for k, v in (extra or {}).items():
                setattr(h, k, v)
        if donor is not None and h.t_reversal is None and reversed_(y):
            h.t_reversal = float(t)
        if donor is not None and (h.peak_rate is None or rate > h.peak_rate):
            h.peak_rate, h.t_peak_rate, h.peak_regime = float(rate), float(t), regime
        if regime != prev:
            h.phases.append((float(t), regime))
            prev = regime
    h.fill = np.array(fills)


def _exact_onset(h, sol, k_onset: int, point) -> None:
    """Onset of the overflow at the root of max(R / R_L) - 1 (non-terminal ODE event) instead
    of the first solution point in transfer (up to max_step = t_end / 2000 later)."""
    te = sol.t_events[k_onset]
    if not len(te):
        return
    t_on, y_on = float(te[0]), sol.y_events[k_onset][0]
    if h.t_event is not None and t_on > h.t_event:
        return
    if h.t_onset is not None and t_on >= h.t_onset:
        return
    h.t_onset = t_on
    if h.donor is None:
        donor, _, _, _, core, extra = point(t_on, y_on, at_root=True)
        h.donor, h.donor_core_onset = donor, float(core)
        for k, v in (extra or {}).items():
            setattr(h, k, v)


def _recorded(point):
    """(onset event, point): the event max(R / R_L) - 1 records point(t, y) wherever it is
    evaluated - solve_ivp evaluates the events at every accepted step - so the scan of the
    solution points after the integration reuses them."""
    seen: dict = {}

    def ev_onset(t, y):
        p = seen[(float(t), *(float(v) for v in y))] = point(t, y)
        return p[3] - 1.0

    def cached(t, y, at_root=False):
        p = None if at_root else seen.get((float(t), *(float(v) for v in y)))
        return point(t, y, at_root) if p is None else p

    return ev_onset, cached


def _pad(ts, ys, t_end: float, lay: _Layout | None = None):
    """Extend the history to t_end with its last state (lay: the stars keep ageing)."""
    if ts[-1] >= t_end:
        return ts, ys
    y_last = ys[:, -1].copy()
    if lay is not None:
        y_last[lay.age : lay.age + lay.n] += (t_end - ts[-1]) / YR_PER_MYR
    return np.append(ts, t_end), np.column_stack([ys, y_last])


def _store(h, ts, ys, lay: _Layout) -> None:
    h.t = ts
    h.m = ys[: lay.n]
    h.r_corr = ys[lay.rc : lay.rc + lay.n]
    h.age = ys[lay.age : lay.age + lay.n]
    h.m_ev = ys[lay.mev : lay.mev + lay.n]
    h.m_core = ys[lay.core : lay.core + lay.n]


def binary_history(cfg, pair, m_labels, a0: float, t_end: float) -> History:
    """Mass-transfer history of the two-star binary ``pair`` = (i, j) with initial separation
    a0 [AU] (circular) over run times 0..t_end [yr] (t_end > 0).

    Each star carries its mass m, an evolutionary coordinate u (equivalent evolutionary point
    or fraction of its phase) and the mass m_ev of the track it follows. As in BSE (Hurley et
    al. 2000, Sect. 7.1; 2002, Sect. 2.6.6), a main-sequence or Hertzsprung-gap star that loses
    or gains mass is re-modelled as a star of its new mass in the same evolutionary state
    (m_ev = m, u kept): a donor shrinks with its mass and crosses the gap on the longer time
    scale of its new mass, an accretor brightens and keeps its fractional age. Giants keep the
    track of the mass they had when they left the gap (their structure follows the core) and
    their radius responds to the envelope mass with zeta_eq = -x + 2 (M_c/M)^5 (eq. 56).

    The rate follows eqs. 58-63, so the donor overfills its lobe just enough to transfer at the
    rate its expansion (or the shrinking lobe) requires. Accretion: eqs. 64-65 (prescription
    bse; the rest leaves with the specific angular momentum of the donor) or the fixed
    fractions alpha, beta."""
    mt = cfg.physics.mass_transfer
    contact = getattr(mt, "contact", None)
    contact_on = contact is not None
    sm = cfg.stellar_model
    lay = _BINARY
    if t_end <= 0:
        raise ValueError("physics.mass_transfer needs a forward run (t_max_yr > 0)")
    n_zero = [0]

    def _terms(t, y):
        """(donor index or None, rate, regime, accreted fraction, alpha, beta, states, lobes,
        radii, raw core masses)."""
        m = (y[0], y[1])
        a = math.exp(y[2])
        st, raw_core, rad = _remodel(sm, y, lay)
        rl = [roche_lobe_eggleton(m[k] / m[1 - k]) * a / units.RSUN_AU for k in (0, 1)]
        over = [rad[k] / rl[k] for k in (0, 1)]
        d = int(np.argmax(over))
        if over[d] <= 1.0:
            return None, 0.0, "detached", 1.0, 0.0, 0.0, st, rl, rad, raw_core
        if contact_on and over[1 - d] > 1.0:
            # contact binary: the less massive star stays at its lobe and feeds the primary
            # through the common envelope (conservative; the mass ratio falls towards the
            # Darwin instability)
            d = int(np.argmin(m))
            rate, _ = _limit_rate(st[d], m[d], bse_rate(m[d], rad[d], rl[d]))
            return d, rate, "contact", 1.0, 0.0, 0.0, st, rl, rad, raw_core
        rate, regime = _limit_rate(st[d], m[d], bse_rate(m[d], rad[d], rl[d]))
        if mt.accretion == "bse":
            f = _accreted_fraction(st[1 - d], m[1 - d], rate)
            alpha, beta = 1.0 - f, 0.0
        else:
            alpha, beta = mt.alpha, mt.beta
            f = 1.0 - alpha - beta
        return d, rate, regime, f, alpha, beta, st, rl, rad, raw_core

    terms = _memoised(_terms)

    def rhs(t, y):
        try:
            d, rate, regime, f, alpha, beta, st, _, _, raw_core = terms(t, y)
        except EvolvedBeyondGrid:
            n_zero[0] += 1
            return [0.0] * lay.size
        dm = [0.0, 0.0]
        dlna = 0.0
        if d is not None:
            dm[d] = -rate
            dm[1 - d] = f * rate
            dlna = dln_a_dmd(y[d], y[1 - d], alpha, beta) * (-rate)
        w = _winds(cfg, pair, y, st, lay)
        # isotropic winds from both stars (Jeans mode): a M = const for the binary
        dlna -= (w[0] + w[1]) / (y[0] + y[1])
        if contact_on and contact.braking and d is not None and regime == "contact":
            dlna += (
                2.0
                * contact_jdot(y[0], y[1], st[0].radius_rsun, st[1].radius_rsun, math.exp(y[2]))
                / orbital_j(y[0], y[1], math.exp(y[2]))
            )
        out = [0.0] * lay.size
        out[0], out[1], out[2] = dm[0] + w[0], dm[1] + w[1], dlna
        _star_derivatives(cfg, out, st, raw_core, dm, w, y, lay)
        return out

    def ev_ce(t, y):
        d, _, regime, _, alpha, beta, st, _, _, _ = terms(t, y)
        if d is None or regime == "contact":  # a contact binary is limited by ev_darwin, ev_l2
            return -1.0
        if mt.stability == "zeta" and st[d].phase >= 3:
            return zeta_lobe(y[d], y[1 - d], alpha, beta) - zeta_ad_giant(y[d], st[d].m_core)
        qc = q_crit(st[d].phase, y[d], st[d].m_core, mt.q_crit_ms)
        return -1.0 if qc is None else y[d] / y[1 - d] - qc

    def ev_darwin(t, y):
        if not contact_on or terms(t, y)[2] != "contact":
            return 1.0
        return min(y[0], y[1]) / max(y[0], y[1]) - contact.q_min

    def ev_l2(t, y):
        _, _, regime, _, _, _, _, _, rad, _ = terms(t, y)
        if not contact_on or regime != "contact":
            return -1.0
        a = math.exp(y[2])
        q = y[1] / y[0]
        return max(rad[k] * units.RSUN_AU / (a * outer_lobe_radius(q, k)) for k in (0, 1)) - 1.0

    def ev_contact(t, y):
        d, _, _, _, _, _, _, rl, rad, _ = terms(t, y)
        if d is None:
            return -1.0
        return rad[1 - d] / rl[1 - d] - 1.0

    def ev_stripped(t, y):
        d, _, _, _, _, _, st, _, _, _ = terms(t, y)
        if d is None:
            return 1.0
        return y[d] - 1.001 * max(st[d].m_core, 1e-3)

    def ev_track_end(t, y):
        return min(
            virtual_end_age(sm, y[lay.mev + k]) * (1 - 1e-9) - y[lay.age + k] for k in (0, 1)
        )

    def point(t, y, at_root=False):
        d, rate, regime, _, _, _, st, rl, rad, _ = terms(t, y)
        fill = max(rad[0] / rl[0], rad[1] / rl[1])
        if at_root:  # R = R_L: the star that fills its lobe
            d = int(np.argmax([rad[0] / rl[0], rad[1] / rl[1]]))
        if d is None:
            return None, regime, 0.0, fill, None, None
        return pair[d], regime, rate, fill, st[d].m_core, None

    ev_onset, point = _recorded(point)

    for e_ in (ev_ce, ev_contact, ev_l2):
        e_.terminal, e_.direction = True, 1
    ev_darwin.terminal, ev_darwin.direction = True, -1
    if contact_on:
        ev_contact.terminal = False
    ev_stripped.terminal, ev_stripped.direction = True, -1
    ev_track_end.terminal, ev_track_end.direction = True, -1
    ev_onset.terminal, ev_onset.direction = False, 1
    m_start, labels, ages, c0 = _start(cfg, pair, m_labels, lay)
    y0 = [m_start[0], m_start[1], math.log(a0), 0.0, 0.0, *ages, *labels, *c0]
    h = History(pair=pair, t=None, m=None, a=None)
    names = (
        "common_envelope",
        "contact",
        "donor_stripped",
        "track_end",
        "contact_merger_darwin",
        "contact_merger_l2",
        "onset",
    )
    sol = _integrate(
        rhs,
        t_end,
        y0,
        [1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-7, 1e-7, 1e-10, 1e-10, 1e-8, 1e-8],
        (ev_ce, ev_contact, ev_stripped, ev_track_end, ev_darwin, ev_l2, ev_onset),
    )
    skip = ("onset", "contact") if contact_on else ("onset",)
    ts, ys, h.event, h.t_event = _first_event(sol, names, sol.t, sol.y, skip)
    where = f"binary m{pair[0] + 1} + m{pair[1] + 1}"
    _check_solution(h, sol, where, n_zero[0])

    def reversed_(y):
        dd = pair.index(h.donor)
        return y[dd] < y[1 - dd]

    _scan(h, ts, ys, point, reversed_, where)
    _exact_onset(h, sol, names.index("onset"), point)
    if h.event == "common_envelope" and h.donor is None:
        d = terms(h.t_event, ys[:, -1])[0]
        h.donor, h.t_onset = pair[d if d is not None else 0], h.t_event
    if h.event == "common_envelope" and mt.common_envelope is not None:
        ts, ys = _common_envelope(h, ts, ys, terms, mt.common_envelope, t_end)
    if h.event != "integration_failed":
        ts, ys = _pad(ts, ys, t_end, lay if h.event == "post_common_envelope" else None)
    _store(h, ts, ys, lay)
    h.a = np.exp(ys[2])
    return h


def _common_envelope(h: History, ts, ys, terms, spec, t_end: float):
    """Outcome of the common envelope at the end of the history (swingersim.commonenv): the
    history continues to the separation of the cores (event post_common_envelope) or ends in
    a merger (event common_envelope_merger)."""
    from . import commonenv

    y = ys[:, -1].copy()
    t_ce = float(ts[-1])
    d = h.pair.index(h.donor)
    acc = 1 - d
    _, _, _, _, _, _, st, _, rad, _ = terms(t_ce, y)
    sd, sa = st[d], st[acc]
    m1, m2 = float(y[d]), float(y[acc])
    mc1 = min(max(sd.m_core, float(y[9 + d])), m1)
    rc1 = commonenv.core_radius(mc1, float(y[7 + d]), sd.phase)
    giant2 = sa.phase >= 2 and sa.m_core > 0.0
    if giant2:  # a giant companion loses its envelope too
        mc2 = min(max(sa.m_core, float(y[9 + acc])), m2)
        rc2 = commonenv.core_radius(mc2, float(y[7 + acc]), sa.phase)
    else:  # a main-sequence companion counts as all core (Hurley et al. 2002, Sect. 2.7.1)
        mc2, rc2 = m2, rad[acc]
    a_i = math.exp(float(y[2]))
    o = commonenv.outcome(m1, mc1, rad[d], rc1, m2, mc2, rad[acc], rc2, a_i, spec.alpha, spec.lam)
    h.ce = {
        "outcome": o.outcome,
        "t_yr": t_ce,
        "a_i_au": a_i,
        "a_f_au": o.a_f,
        "donor_mass": m1,
        "donor_core": mc1,
        "companion_mass": m2,
        "core_radius_rsun": rc1,
        "merged_mass": o.m1_f if o.outcome == "merge" else None,
    }
    if o.outcome == "merge":
        h.event = "common_envelope_merger"
        return ts, ys
    t1 = min(t_ce + spec.duration_yr, t_end)
    y1 = y.copy()
    y1[d] = mc1
    y1[acc] = mc2 if giant2 else m2
    y1[2] = math.log(o.a_f)
    y1[5:7] += (t1 - t_ce) / YR_PER_MYR
    y1[9 + d] = mc1
    h.event = "post_common_envelope"
    h.remnants = {d: (mc1, rc1)}
    if giant2:
        h.remnants[acc] = (mc2, rc2)
    return np.append(ts, t1), np.column_stack([ys, y1])


# ------------------------------------------------------------------ tertiary mass transfer
def zeta_ad_giant(m: float, m_core: float) -> float:
    """Adiabatic mass-radius exponent of a giant (condensed polytrope n = 3/2 with a core of
    fractional mass m_c = M_c / M; Hjellming & Webbink 1987, fit of Soberman, Phinney & van den
    Heuvel 1997, eqs. 60-61)."""
    mc = min(max(m_core / m, 0.0), 0.95)
    z1 = 2.0 / 3.0 * mc / (1.0 - mc) - 1.0 / 3.0 * (1.0 - mc) / (1.0 + 2.0 * mc)
    return z1 - 0.03 * mc + 0.2 * mc / (1.0 + (1.0 - mc) ** -6)


def circularisation_distance(a_peri: float, q_donor: float) -> float:
    """Distance [same unit as a_peri] from the accreting subsystem within which the stream of a
    tertiary donor forms a circumbinary disc: a_c = a (1 - r_L / a)^2 with r_L the lobe radius
    of the donor of mass ratio q_donor = M_donor / M_inner (Leigh et al. 2020, eq. 5)."""
    return a_peri * (1.0 - roche_lobe_eggleton(q_donor)) ** 2


@dataclass
class TertiaryHistory(_Track):
    """Mass transfer from a star (the tertiary) onto an inner two-star binary it orbits
    (swingersim.masstransfer.tertiary_history). Stars: (donor, inner primary, inner secondary)."""

    pair: tuple[int, int, int]  # donor, inner stars (named 'pair' like History: star indices)
    outer: str
    inner: str
    t: np.ndarray
    m: np.ndarray  # (3, n_t)
    a: np.ndarray  # outer separation [AU]
    a_in: np.ndarray  # inner separation [AU]
    donor: int | None = None
    t_onset: float | None = None
    t_reversal: float | None = None
    event: str | None = None
    t_event: float | None = None
    phases: list[tuple[float, str]] = field(default_factory=list)
    r_corr: np.ndarray | None = None
    age: np.ndarray | None = None
    m_ev: np.ndarray | None = None
    m_core: np.ndarray | None = None
    disc: bool | None = None  # the stream forms a circumbinary disc around the inner binary
    accreted_fraction: float | None = None
    ce: dict | None = None
    remnants: dict = field(default_factory=dict)
    e_out: float = 0.0
    fill: np.ndarray | None = None  # R / R_L of the donor at the history points
    peak_rate: float | None = None
    t_peak_rate: float | None = None
    peak_regime: str | None = None
    donor_core_onset: float | None = None
    n_zero_rhs: int = 0

    def masses_at(self, t: float) -> np.ndarray:
        return np.array([float(np.interp(t, self.t, self.m[k])) for k in range(3)])

    def a_in_at(self, t: float) -> float:
        return float(np.interp(t, self.t, self.a_in))

    def orbit_a(self, t: float) -> dict[str, float]:
        return {self.outer: self.a_at(t), self.inner: self.a_in_at(t)}

    def a_factors(self, t0: float, t1: float) -> dict[str, float]:
        """Factors on a beyond a M = const (applied by the backends themselves)."""
        m0, m1 = self.masses_at(t0), self.masses_at(t1)
        f_out = self.a_at(t1) / self.a_at(t0) * (m1.sum() / m0.sum())
        f_in = self.a_in_at(t1) / self.a_in_at(t0) * ((m1[1] + m1[2]) / (m0[1] + m0[2]))
        return {self.outer: f_out, self.inner: f_in}

    def donor_lobe_rsun(self, t: float) -> float:
        m = self.masses_at(t)
        peri = self.a_at(t) * (1.0 - self.e_out)
        return roche_lobe_eggleton(m[0] / (m[1] + m[2])) * peri / units.RSUN_AU

    def describe(self) -> str:
        return (
            f"{self.event} (tertiary m{self.pair[0] + 1} onto the binary {self.inner} "
            f"in orbit {self.outer})"
        )

    def summary(self) -> dict:
        return {
            "donor": None if self.donor is None else f"m{self.donor + 1}",
            "t_onset_yr": self.t_onset,
            "t_reversal_yr": self.t_reversal,
            "event": self.event,
            "t_event_yr": self.t_event,
            "m_final": [float(x) for x in self.m[:, -1]],
            "a_initial_au": float(self.a[0]),
            "a_final_au": float(self.a[-1]),
            "a_inner_initial_au": float(self.a_in[0]),
            "a_inner_final_au": float(self.a_in[-1]),
            "circumbinary_disc": self.disc,
            "accreted_fraction": self.accreted_fraction,
            **self.episode_summary(),
        }


def tertiary_history(cfg, stars, names, m_labels, orbit, t_end: float) -> TertiaryHistory:
    """History of a tertiary (stars[0]) orbiting the two-star binary of stars[1], stars[2]:
    its evolution, its Roche-lobe overflow onto the inner binary and the response of both
    orbits (physics.mass_transfer.tertiary).

    The inner binary is a point mass for the donor: its lobe is that of a binary of mass ratio
    M_3 / (M_1 + M_2) at the periastron of the outer orbit (Eggleton 1983), the rate follows
    BSE (Hurley et al. 2002, eqs. 58-63; Hamers et al. 2022). The stream forms a circumbinary
    disc when the inner separation is below a_c = a_peri (1 - r_L/a)^2 (Leigh et al. 2020,
    eq. 5); then 90 % of the transferred mass is accreted (60 % by the primary, 30 % by the
    secondary), otherwise 20 % (10 % each; Hamers et al. 2022, Sect. 2.2). The rest leaves the
    inner binary isotropically: for the outer orbit it is re-emitted by the accretor
    (Soberman et al. 1997), and its ejection costs inner orbital energy with efficiency
    alpha lambda (de Vries et al. 2014, eqs. 2-3):

        dln a_in = dM_1/M_1 + dM_2/M_2 - 2 M_in dM_ej / (alpha lambda M_1 M_2),

    or, with inner_orbit: angular_momentum, the accreted gas brings the specific angular
    momentum h sqrt(G M_in a_in): dln a_in = 2 h dM_acc / mu - 2 (dM_1/M_1 + dM_2/M_2) +
    dM_acc / M_in. Transfer is dynamically unstable (triple common envelope, the run stops) for
    q = M_3 / M_in above q_crit (Hurley et al. 2002, eq. 57); it also stops when an inner star
    fills its own lobe or the triple becomes unstable (Mardling & Aarseth 2001; a triple that
    is unstable already at t = 0 is reported with a warning)."""
    from .stability import ma01_ratio

    mt = cfg.physics.mass_transfer
    spec = mt.tertiary
    sm = cfg.stellar_model
    lay = _TERTIARY
    a_out0, e_out, a_in0, e_in, i_mut = orbit
    n_zero = [0]

    def _terms(t, y):
        m = y[0:3]
        a_out, a_in = math.exp(y[3]), math.exp(y[4])
        st, raw_core, rad = _remodel(sm, y, lay)
        m_in = m[1] + m[2]
        peri = a_out * (1.0 - e_out)
        rl0 = roche_lobe_eggleton(m[0] / m_in) * peri / units.RSUN_AU
        rl_in = [
            roche_lobe_eggleton(m[1] / m[2]) * a_in * (1.0 - e_in) / units.RSUN_AU,
            roche_lobe_eggleton(m[2] / m[1]) * a_in * (1.0 - e_in) / units.RSUN_AU,
        ]
        if rad[0] <= rl0:
            return 0.0, "detached", 1.0, (0.5, 0.5), st, rl0, rl_in, rad, raw_core, None
        rate, regime = _limit_rate(st[0], m[0], bse_rate(m[0], rad[0], rl0))
        disc = a_in < circularisation_distance(peri, m[0] / m_in)
        if spec.accretion == "auto":
            f = 0.9 if disc else 0.2
        else:
            f = float(spec.accretion)
        if spec.split == "auto":
            s1 = 2.0 / 3.0 if disc else 0.5
        else:
            s1 = float(spec.split)
        split = (s1, 1.0 - s1) if m[1] >= m[2] else (1.0 - s1, s1)
        return rate, regime, f, split, st, rl0, rl_in, rad, raw_core, disc

    terms = _memoised(_terms)

    def rhs(t, y):
        try:
            rate, _, f, split, st, _, _, _, raw_core, _ = terms(t, y)
        except EvolvedBeyondGrid:
            n_zero[0] += 1
            return [0.0] * lay.size
        m = y[0:3]
        m_in = m[1] + m[2]
        dm = [-rate, f * split[0] * rate, f * split[1] * rate]
        w = _winds(cfg, stars, y, st, lay)
        dlna_out = dln_a_dmd(m[0], m_in, 0.0, 1.0 - f) * (-rate) if rate > 0 else 0.0
        dlna_out -= (w[0] + w[1] + w[2]) / (m[0] + m_in)
        if rate > 0:
            if spec.inner_orbit == "energy":
                dlna_in = dm[1] / m[1] + dm[2] / m[2]
                dlna_in -= 2.0 * m_in * (1.0 - f) * rate / (spec.alpha_lambda * m[1] * m[2])
            else:
                mu = m[1] * m[2] / m_in
                dacc = dm[1] + dm[2]
                dlna_in = 2.0 * spec.h * dacc / mu - 2.0 * (dm[1] / m[1] + dm[2] / m[2])
                dlna_in += dacc / m_in
        else:
            dlna_in = 0.0
        dlna_in -= (w[1] + w[2]) / m_in
        out = [0.0] * lay.size
        out[0:5] = [dm[k] + w[k] for k in range(3)] + [dlna_out, dlna_in]
        _star_derivatives(cfg, out, st, raw_core, dm, w, y, lay)
        return out

    def ev_tce(t, y):
        rate, _, f, _, st, *_ = terms(t, y)
        if rate <= 0:
            return -1.0
        q = y[0] / (y[1] + y[2])
        if mt.stability == "zeta" and st[0].phase >= 3:  # the accreted fraction f of the orbit
            return zeta_lobe(y[0], y[1] + y[2], 0.0, 1.0 - f) - zeta_ad_giant(y[0], st[0].m_core)
        qc = q_crit(st[0].phase, y[0], st[0].m_core, mt.q_crit_ms)
        return -1.0 if qc is None else q - qc

    def ev_inner(t, y):
        _, _, _, _, _, _, rl_in, rad, _, _ = terms(t, y)
        return max(rad[1] / rl_in[0], rad[2] / rl_in[1]) - 1.0

    def ev_unstable(t, y):
        r = ma01_ratio(math.exp(y[4]), math.exp(y[3]), e_out, y[0] / (y[1] + y[2]), i_mut)
        return r - 1.0

    def ev_stripped(t, y):
        rate, _, _, _, st, *_ = terms(t, y)
        if rate <= 0:
            return 1.0
        return y[0] - 1.001 * max(st[0].m_core, 1e-3)

    def ev_track_end(t, y):
        return min(
            virtual_end_age(sm, y[lay.mev + k]) * (1 - 1e-9) - y[lay.age + k] for k in range(3)
        )

    def point(t, y, at_root=False):
        rate, regime, f, _, st, rl0, _, rad, _, disc = terms(t, y)
        fill = rad[0] / rl0
        if at_root:  # R = R_L: the stream starts
            m_in = y[1] + y[2]
            peri = math.exp(y[3]) * (1.0 - e_out)
            disc = math.exp(y[4]) < circularisation_distance(peri, y[0] / m_in)
            f = (0.9 if disc else 0.2) if spec.accretion == "auto" else float(spec.accretion)
        elif rate <= 0:
            return None, regime, 0.0, fill, None, None
        extra = {"disc": bool(disc), "accreted_fraction": float(f)}
        return stars[0], regime, rate, fill, st[0].m_core, extra

    ev_onset, point = _recorded(point)

    for e_ in (ev_tce, ev_inner):
        e_.terminal, e_.direction = True, 1
    for e_ in (ev_unstable, ev_stripped, ev_track_end):
        e_.terminal, e_.direction = True, -1
    ev_onset.terminal, ev_onset.direction = False, 1
    m_start, labels, ages, c0 = _start(cfg, stars, m_labels, lay)
    y0 = [*m_start, math.log(a_out0), math.log(a_in0), 0.0, 0.0, 0.0, *ages, *labels, *c0]
    h = TertiaryHistory(
        pair=tuple(stars),
        outer=names[0],
        inner=names[1],
        t=None,
        m=None,
        a=None,
        a_in=None,
        e_out=e_out,
    )
    ev_names = (
        "triple_common_envelope",
        "inner_roche_overflow",
        "dynamical_instability",
        "donor_stripped",
        "track_end",
        "onset",
    )
    where = f"tertiary m{stars[0] + 1} onto {names[1]}"
    if ev_unstable(0.0, y0) <= 0.0:
        # the event (direction -1) catches a triple that becomes unstable during the transfer;
        # one that is unstable already at the start is the user's choice of a secular run
        warnings.warn(
            f"physics.mass_transfer.tertiary: the triple of {where} is dynamically unstable "
            "already at t = 0 (Mardling & Aarseth 2001); the secular treatment is not valid",
            stacklevel=2,
        )
    sol = _integrate(
        rhs,
        t_end,
        y0,
        [1e-10] * 8 + [1e-7] * 3 + [1e-10] * 3 + [1e-8] * 3,
        (ev_tce, ev_inner, ev_unstable, ev_stripped, ev_track_end, ev_onset),
    )
    ts, ys, h.event, h.t_event = _first_event(sol, ev_names, sol.t, sol.y, ("onset",))
    _check_solution(h, sol, where, n_zero[0])

    y_first = ys[:, 0].copy()

    def reversed_(y):
        return y_first[0] > y_first[1] + y_first[2] and y[0] < y[1] + y[2]

    _scan(h, ts, ys, point, reversed_, where)
    _exact_onset(h, sol, ev_names.index("onset"), point)
    if h.event == "triple_common_envelope" and h.t_onset is None:
        h.t_onset, h.donor = h.t_event, stars[0]  # unstable from the first overflow
    if h.event == "track_end" and h.t_onset is None:
        h.event, h.t_event = None, None  # ended by the track like any other run
    if h.event != "integration_failed":
        ts, ys = _pad(ts, ys, t_end)
    _store(h, ts, ys, lay)
    h.a, h.a_in = np.exp(ys[3]), np.exp(ys[4])
    return h


# ------------------------------------------------------------------ coupling to the dynamics
_CACHE: dict[str, dict[str, History]] = {}


def enabled(cfg) -> bool:
    return getattr(cfg.physics, "mass_transfer", None) is not None


def dynamic_onset(cfg) -> str | None:
    """physics.mass_transfer.dynamic_onset (None without mass transfer)."""
    spec = getattr(cfg.physics, "mass_transfer", None)
    return None if spec is None else spec.dynamic_onset


def stop_on_overflow(cfg) -> bool:
    """Does an unexpected Roche-lobe overflow end the run (run.stop_on_roche_overflow, or
    physics.mass_transfer.dynamic_onset, which then continues it: swingersim.merger)?"""
    return bool(cfg.run.stop_on_roche_overflow) or dynamic_onset(cfg) is not None


def histories(cfg) -> dict:
    """Mass-transfer histories (orbit name -> History or TertiaryHistory), cached per
    configuration: every two-star binary, and (secular backend, physics.mass_transfer.tertiary)
    every star orbiting a two-star binary, whose history then also carries that binary."""
    if not enabled(cfg):
        return {}
    key = cfg.model_dump_json()
    if key in _CACHE:
        return _CACHE[key]
    import warnings

    from .hierarchy import build_tree, initialize, node_h_unit

    root, m_cfg = build_tree(cfg)
    out = {}
    for n in root.nodes():
        if isinstance(n.left, int) and isinstance(n.right, int):
            if n.elements.e > 0.1:
                warnings.warn(
                    f"physics.mass_transfer: binary {n.name} has e = {n.elements.e:.2f}; its "
                    "mass transfer is computed for a circular orbit of the same a",
                    stacklevel=2,
                )
            out[n.name] = binary_history(
                cfg, (n.left, n.right), m_cfg, float(n.elements.a), float(cfg.run.t_max_yr)
            )
    spec = cfg.physics.mass_transfer
    if spec.tertiary is not None and cfg.run.backend == "secular":
        pos, vel = initialize(root, m_cfg)
        for n in root.nodes():
            kids = (n.left, n.right)
            star = [c for c in kids if isinstance(c, int)]
            sub = [c for c in kids if not isinstance(c, int)]
            if len(star) != 1 or len(sub) != 1:
                continue
            inner = sub[0]
            if not (isinstance(inner.left, int) and isinstance(inner.right, int)):
                continue
            hin = out.get(inner.name)
            if hin is not None and hin.t_onset is not None:
                warnings.warn(
                    f"physics.mass_transfer.tertiary: the inner binary {inner.name} transfers "
                    f"mass itself; a Roche-lobe overflow of m{star[0] + 1} ends the run",
                    stacklevel=2,
                )
                continue
            h_out = node_h_unit(n, pos, vel, m_cfg)
            h_in = node_h_unit(inner, pos, vel, m_cfg)
            i_mut = float(np.arccos(np.clip(h_out @ h_in, -1.0, 1.0)))
            th = tertiary_history(
                cfg,
                (star[0], inner.left, inner.right),
                (n.name, inner.name),
                m_cfg,
                (
                    float(n.elements.a),
                    float(n.elements.e),
                    float(inner.elements.a),
                    float(inner.elements.e),
                    i_mut,
                ),
                float(cfg.run.t_max_yr),
            )
            if th.t_onset is None and th.event is None:
                continue  # never overflows: the two-star history describes the inner binary
            out.pop(inner.name, None)
            out[n.name] = th
    if len(_CACHE) > 64:
        _CACHE.clear()
    _CACHE[key] = out
    return out


def orbit_a_at(cfg, t_yr: float) -> dict[str, float]:
    """Semi-major axes [AU] of the orbits whose separation a mass-transfer history prescribes."""
    out = {}
    for name, h in histories(cfg).items():
        if isinstance(h, TertiaryHistory):
            out.update(h.orbit_a(t_yr))
        else:
            out[name] = h.a_at(t_yr)
    return out


def a_factors(cfg, t0: float, t1: float) -> dict[str, float]:
    """Factors on the semi-major axes from t0 to t1 beyond a M = const, for every orbit with a
    mass-transfer history."""
    out = {}
    for name, h in histories(cfg).items():
        if isinstance(h, TertiaryHistory):
            out.update(h.a_factors(t0, t1))
        else:
            out[name] = a_ratio(cfg, name, t0, t1)
    return out


def _thin(h, thr: float) -> np.ndarray:
    """History times at which ln m, ln a, ln m_ev or the radius correction of any star has
    moved by more than thr since the previous kept time (bands of width thr; between two kept
    times every quantity stays within one band), plus the first and last points, the onset,
    the reversal, the event, the common envelope and the start of every regime."""
    rows = [np.log(np.maximum(h.m, 1e-300)), np.log(np.maximum(h.m_ev, 1e-300)), h.r_corr]
    rows.append(np.log(np.atleast_2d(h.a)))
    if getattr(h, "a_in", None) is not None:
        rows.append(np.log(np.atleast_2d(h.a_in)))
    q = np.floor(np.vstack(rows) / thr)
    keep = np.zeros(len(h.t), bool)
    keep[[0, -1]] = True
    keep[1:] |= np.any(q[:, 1:] != q[:, :-1], axis=0)
    extra = [h.t_onset, h.t_reversal, h.t_event, h.t_peak_rate, *(t for t, _ in h.phases)]
    if h.ce is not None:
        extra.append(h.ce["t_yr"])
    extra = np.array([t for t in extra if t is not None], float)
    return np.concatenate([h.t[keep], extra])


def time_grid(cfg) -> np.ndarray:
    """Times that resolve the transfer in the operator splitting: the history points at which
    a mass, separation or radius correction has changed by more than run.evolution_max_dr / 4
    since the previous one, and the onsets, reversals, events and regime changes (the
    segmentation itself refines the grid wherever the structure changes faster)."""
    hs = histories(cfg).values()
    if not hs:
        return np.array([])
    thr = 0.25 * float(cfg.run.evolution_max_dr)
    return np.unique(np.concatenate([_thin(h, thr) for h in hs]))


def a_ratio(cfg, name: str, t0: float, t1: float) -> float:
    """Factor on the semi-major axis of orbit ``name`` from t0 to t1 beyond the adiabatic
    response to its mass change (a M = const), which the backends apply themselves."""
    h = histories(cfg).get(name)
    if h is None:
        return 1.0
    m0, m1 = sum(h.masses_at(t0)), sum(h.masses_at(t1))
    return h.a_at(t1) / h.a_at(t0) * (m1 / m0)


def stop(cfg) -> tuple[float, str] | None:
    """Earliest event that ends a transfer (and the run): (time, description)."""
    best = None
    for name, h in histories(cfg).items():
        if h.event is None or h.event == "track_end":
            continue
        if h.event == "post_common_envelope" and cfg.run.backend == "nbody":
            continue  # the N-body backend follows the impulsive loss of the envelope
        if isinstance(h, TertiaryHistory):
            what = h.describe()
        else:
            d = h.donor
            if d is None:  # integration_failed before any transfer
                what = f"{h.event} (orbit {name})"
            else:
                acc = h.pair[1] if d == h.pair[0] else h.pair[0]
                what = f"{h.event} (orbit {name}, donor m{d + 1}, accretor m{acc + 1})"
            if h.event == "post_common_envelope":
                what += "; the envelope leaves impulsively: use backend nbody to follow it"
        if best is None or h.t_event < best[0]:
            best = (h.t_event, what)
    return best


def transferring_stars(cfg) -> set[int]:
    """Stars of binaries that go through a transfer episode (their Roche-lobe filling is
    expected and does not end the run)."""
    return {k for h in histories(cfg).values() if h.t_onset is not None for k in h.pair}


def star_state(cfg, k: int, t_yr: float):
    """Re-modelled state of star k after the onset of its binary's mass transfer (its mass,
    luminosity and radius differ from the track of its initial mass); None otherwise."""
    for h in histories(cfg).values():
        if k in h.pair and h.t_onset is not None and t_yr >= h.t_onset and h.m_ev is not None:
            try:
                return h.state_at(cfg.stellar_model, h.pair.index(k), t_yr)
            except ValueError:
                return None
    return None


def track_point(cfg, k: int, t_yr: float) -> tuple[float, float] | None:
    """(m_ev, age on that track [Myr]) of star k at t_yr after the onset of its binary's mass
    transfer: the track the re-modelled star follows; None before the onset or without
    transfer."""
    for h in histories(cfg).values():
        if k in h.pair and h.t_onset is not None and t_yr >= h.t_onset and h.m_ev is not None:
            i = h.pair.index(k)
            m_ev = float(np.interp(t_yr, h.t, h.m_ev[i]))
            age = float(np.interp(t_yr, h.t, h.age[i]))
            if i in h.remnants and h.ce is not None and t_yr >= h.ce["t_yr"]:
                j = int(np.searchsorted(h.t, h.ce["t_yr"]))
                m_ev, age = float(h.m_ev[i][j]), float(h.age[i][j])
            return m_ev, age
    return None


def track_limit(cfg, k: int) -> float | None:
    """Run time [yr] at which star k leaves the track it follows after a mass transfer (inf for
    the core left by a common envelope); None when star k never transfers mass."""
    for h in histories(cfg).values():
        if (
            isinstance(h, TertiaryHistory)
            and h.event == "track_end"
            and h.t_event is not None
            and k == h.pair[0]
        ):
            return float(h.t_event)  # the tertiary donor reached the end of its track
        if k in h.pair and h.t_onset is not None and h.m_ev is not None:
            i = h.pair.index(k)
            if i in h.remnants:
                return math.inf
            t_last = float(h.t[-1])
            left = virtual_end_age(cfg.stellar_model, float(h.m_ev[i][-1])) - float(h.age[i][-1])
            return t_last + left * YR_PER_MYR
    return None


FILL_TOL = 0.02  # relative margin of the Roche-lobe filling (2 x run.evolution_max_dr default)


def roche_excluded(cfg, k: int, t_yr: float, e: float = 0.0) -> bool:
    """True when the Roche-lobe filling of star k at t_yr belongs to its modelled mass
    transfer (it does not end the run): from the onset to the end of the transfer, or to the
    end of the common envelope.

    Shortly before the onset the dynamics may already see an overflow that belongs to it: the
    history is circular, while the dynamical orbit (eccentricity e, e.g. ~0.01 forced by a
    tertiary) reaches the lobe at periastron, a (1 - e), and the segments hold the stellar
    structure (radius) at their middle, within run.evolution_max_dr. Such an overflow is
    excluded when the circular history fills the lobe to at least (1 - e)(1 - FILL_TOL)."""
    for h in histories(cfg).values():
        if k not in h.pair or h.t_onset is None:
            continue
        if t_yr >= h.t_onset:
            if h.event == "post_common_envelope":
                return t_yr <= h.t[np.searchsorted(h.t, h.ce["t_yr"]) + 1]
            if h.regime_at(t_yr) != "detached":
                return True  # during the transfer (after it: as before the onset)
        if h.fill_at(t_yr) >= (1.0 - min(max(e, 0.0), 0.5)) * (1.0 - FILL_TOL):
            return True
    return False


def donor_radius(cfg, k: int, t_yr: float) -> float | None:
    """Radius [Rsun] of star k while it transfers mass: its Roche lobe (the donor is held at
    its lobe); None otherwise."""
    for h in histories(cfg).values():
        if h.donor == k and h.transferring(t_yr):
            if isinstance(h, TertiaryHistory):
                return h.donor_lobe_rsun(t_yr)
            d = h.pair.index(k)
            m = h.masses_at(t_yr)
            return roche_lobe_eggleton(m[d] / m[1 - d]) * h.a_at(t_yr) / units.RSUN_AU
    return None


def summary(cfg) -> dict[str, dict]:
    """Per binary: onset, q reversal, event, final masses and separation."""
    out = {}
    for name, h in histories(cfg).items():
        if isinstance(h, TertiaryHistory):
            out[name] = h.summary()
            continue
        out[name] = {
            "donor": None if h.donor is None else f"m{h.donor + 1}",
            "t_onset_yr": h.t_onset,
            "t_reversal_yr": h.t_reversal,
            "event": h.event,
            "t_event_yr": h.t_event,
            "m_final": [float(x) for x in h.m[:, -1]],
            "a_initial_au": float(h.a[0]),
            "a_final_au": float(h.a[-1]),
            **h.episode_summary(),
        }
        if h.ce is not None:
            out[name]["common_envelope"] = h.ce
    return out


E_TRANSFER_WARN = 0.1  # eccentricity of a transferring binary above which the circular
# history of the transfer no longer describes the dynamics


def eccentric_transfer(cfg, res) -> dict:
    """The transfer histories are circular (as BSE): when the dynamics of the run drives the
    transferring binary eccentric (e.g. Kozai-Lidov cycles once the binary widened, without
    tides that would damp them), the history is inconsistent with the orbit. Returns
    mass_transfer_<orbit>_e_max_transferring and warns above E_TRANSFER_WARN."""
    out = {}
    times = np.asarray(getattr(res, "times", []), float)
    for name, h in histories(cfg).items():
        if isinstance(h, TertiaryHistory) or h.t_onset is None or name not in res.orbits:
            continue
        e = np.asarray(res.orbits[name].get("e", []), float)
        if len(e) != len(times) or not len(e):
            continue
        on = np.array([h.transferring(float(t)) for t in times])
        if not on.any():
            continue
        e_max = float(np.nanmax(np.where(on, e, np.nan)))
        out[f"mass_transfer_{name}_e_max_transferring"] = e_max
        if e_max > E_TRANSFER_WARN:
            warnings.warn(
                f"physics.mass_transfer: binary {name} reaches e = {e_max:.2f} while it "
                "transfers mass; the transfer is modelled for a circular orbit (physics.tides "
                "would circularise it) - the history does not describe this orbit",
                stacklevel=2,
            )
    return out


def system_entries(cfg) -> dict:
    """Flat summary for RunResult.system (mass_transfer_<orbit>_<key>)."""
    out = {}
    for name, d in summary(cfg).items():
        for k, v in d.items():
            if isinstance(v, dict):
                for kk, vv in v.items():
                    out[f"mass_transfer_{name}_ce_{kk}"] = vv
            else:
                out[f"mass_transfer_{name}_{k}"] = v
    return out
