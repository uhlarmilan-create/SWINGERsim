"""Analytic tidal and rotational formulae for a binary of two extended stars.

* Apsidal motion from the tidal and rotational distortion of both stars, Sterne (1939),
  MNRAS 99, 451 (classical form used for eclipsing binaries, e.g. Claret & Gimenez 1993):

      d(omega)/dt / n = sum_i k2_i (R_i/a)^5 [ 15 (m_j/m_i) f(e)
                                              + (1 + m_j/m_i) (Omega_i/n)^2 g(e) ]
      f(e) = (1 + 3/2 e^2 + 1/8 e^4) / (1 - e^2)^5,   g(e) = 1 / (1 - e^2)^2

  with k2_i the apsidal-motion constant and spins aligned with the orbit.

* Secular evolution under the constant-time-lag model (Mignard 1979; Hut 1981), in the form of
  Leconte et al. (2010), A&A 516, A64, eqs. (1)-(3) with aligned spins, potential Love number
  k2 = 2 x apsidal constant and time lag dt_i:

      K_i = 3/2 k2_i dt_i (G m_j^2 / R_i) (R_i/a)^6 n^2
      da/dt = 4 a^2 / (G m_1 m_2) sum_i K_i [ N(e) Omega_i/n - Na1(e) ]
      de/dt = 11 a e / (G m_1 m_2) sum_i K_i [ Om_e(e) Omega_i/n - 18/11 Ne1(e) ]
      dOmega_i/dt = -2 K_i / (C_i n) [ Om(e) Omega_i/n - N(e) ],   C_i = rg2_i m_i R_i^2

  The spin equation is written from the balance of the tidal torque with the change of the
  orbital angular momentum (it reduces to the circular-orbit result of Murray & Dermott 1999,
  eq. 4.159); the same functions of e appear in Hut (1981), eqs. (9)-(11).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import units


def _f_apsidal(e: float) -> float:
    e2 = e * e
    return (1.0 + 1.5 * e2 + e2 * e2 / 8.0) / (1.0 - e2) ** 5


def _g_apsidal(e: float) -> float:
    return 1.0 / (1.0 - e * e) ** 2


def apsidal_rate_sterne(
    m: tuple[float, float],
    r_au: tuple[float, float],
    apsidal_k2: tuple[float, float],
    a: float,
    e: float,
    spin: tuple[float, float] = (0.0, 0.0),
) -> float:
    """Classical apsidal advance d(omega)/dt [rad/yr] from tides and rotation (Sterne 1939).

    ``spin`` are the angular velocities of the stars [rad/yr], aligned with the orbit.
    """
    n = np.sqrt(units.G * (m[0] + m[1]) / a**3)
    rate = 0.0
    for i in (0, 1):
        j = 1 - i
        q = m[j] / m[i]
        x = (r_au[i] / a) ** 5
        rate += (
            apsidal_k2[i]
            * x
            * (15.0 * q * _f_apsidal(e) + (1.0 + q) * (spin[i] / n) ** 2 * _g_apsidal(e))
        )
    return float(rate * n)


def pseudo_sync_factor(e: float) -> float:
    """Pseudo-synchronous spin over the mean motion, Hut (1981) eq. 42."""
    e2 = e * e
    return (1.0 + 7.5 * e2 + 45.0 / 8.0 * e2**2 + 5.0 / 16.0 * e2**3) / (
        (1.0 + 3.0 * e2 + 3.0 / 8.0 * e2**2) * (1.0 - e2) ** 1.5
    )


# backward runs: a spin whose synchronisation time is shorter than |t_max| / SLAVE_FACTOR is
# held pseudo-synchronous (see slaved_spin)
SLAVE_FACTOR = 3.0


def sync_time_yr(m: float, m_other: float, r_au: float, k2: float, rg2: float, lag_yr: float,
                 a: float) -> float:  # fmt: skip
    """e-folding time [yr] of the spin of a star towards synchronism in a circular orbit,
    t_sync = C n^2 / (2 K) (constant time lag, Leconte et al. 2010; the rates of
    swingersim.secular.kernel.leconte_pair); inf without dissipation or moment of inertia."""
    if lag_yr <= 0 or rg2 <= 0 or k2 <= 0:
        return float("inf")
    n2 = units.G * (m + m_other) / a**3
    K = 1.5 * (2.0 * k2) * lag_yr * (units.G * m_other**2 / r_au) * (r_au / a) ** 6 * n2
    return float(rg2 * m * r_au**2 * n2 / (2.0 * K))


def slaved_spin(backward_span_yr: float, t_sync_yr: float) -> bool:
    """Is the spin held pseudo-synchronous? Only in a backward run (span > 0): integrated into
    the past, tidal dissipation drives the spin away from synchronism, a deviation growing as
    exp(|t| / t_sync); over |t| >> t_sync the spin equation is unstable (the integration ends
    in nonsense or fails). In the forward direction the same spin is attracted to
    pseudo-synchronism on t_sync, so holding it there is the solution the backward run
    should reproduce (the adiabatic limit used for the past tidal evolution of orbits,
    e.g. Jackson, Greenberg & Barnes 2008)."""
    return backward_span_yr > 0 and SLAVE_FACTOR * t_sync_yr < backward_span_yr


def apsidal_rate_gr(m_total: float, a: float, e: float) -> float:
    """Relativistic apsidal advance [rad/yr] (Einstein)."""
    n = np.sqrt(units.G * m_total / a**3)
    return float(3.0 * units.G * m_total * n / (units.C**2 * a * (1.0 - e * e)))


@dataclass
class TidalRates:
    da_dt: float  # [AU/yr]
    de_dt: float  # [1/yr]
    dspin_dt: tuple[float, float]  # [rad/yr^2]
    t_sync: tuple[float, float] = (np.inf, np.inf)  # e-folding time of the spin [yr]
    spin_pseudo_sync: float = np.nan  # pseudo-synchronous angular velocity [rad/yr]

    def timescales(self, a: float, e: float) -> dict:
        """Characteristic times [yr]: a/|da/dt|, e/|de/dt| and the e-folding time of the
        approach of each spin to pseudo-synchronism. (Further positional arguments - spins and
        mean motion, passed by older callers - are not needed and ignored.)"""
        return {
            "t_a_yr": abs(a / self.da_dt) if self.da_dt else np.inf,
            "t_circ_yr": abs(e / self.de_dt) if (self.de_dt and e > 0) else np.inf,
            "t_sync1_yr": self.t_sync[0],
            "t_sync2_yr": self.t_sync[1],
            "P_pseudo_sync_d": 2.0 * np.pi / self.spin_pseudo_sync * units.DAYS_PER_YEAR,
        }


def leconte_rates(
    m: tuple[float, float],
    r_au: tuple[float, float],
    apsidal_k2: tuple[float, float],
    rg2: tuple[float, float],
    lag_yr: tuple[float, float],
    a: float,
    e: float,
    spin: tuple[float, float],
) -> TidalRates:
    """Orbit-averaged tidal evolution, constant time lag, aligned spins (Leconte et al. 2010).
    A star with lag <= 0 does not dissipate (no contribution, spin constant); a dissipating
    star with rg2 <= 0 is taken as pseudo-synchronous with its spin held (t_sync = 0)."""
    G = units.G
    n = np.sqrt(G * (m[0] + m[1]) / a**3)
    e2 = e * e
    b = 1.0 - e2
    N = (1 + 15 / 2 * e2 + 45 / 8 * e2**2 + 5 / 16 * e2**3) / b**6
    Na1 = (1 + 31 / 2 * e2 + 255 / 8 * e2**2 + 185 / 16 * e2**3 + 25 / 64 * e2**4) / b**7.5
    Om_e = (1 + 3 / 2 * e2 + 1 / 8 * e2**2) / b**5
    Ne1 = (1 + 15 / 4 * e2 + 15 / 8 * e2**2 + 5 / 64 * e2**3) / b**6.5
    Om = (1 + 3 * e2 + 3 / 8 * e2**2) / b**4.5
    da = de = 0.0
    dw, tsync = [], []
    for i in (0, 1):
        j = 1 - i
        if lag_yr[i] <= 0.0:  # no dissipation: no contribution, spin constant
            dw.append(0.0)
            tsync.append(np.inf)
            continue
        love = 2.0 * apsidal_k2[i]
        K = 1.5 * love * lag_yr[i] * (G * m[j] ** 2 / r_au[i]) * (r_au[i] / a) ** 6 * n**2
        C = rg2[i] * m[i] * r_au[i] ** 2
        # no moment of inertia (rg2 <= 0): the limit rg2 -> 0, a pseudo-synchronous star
        # (spin/n = N/Om, Hut 1981) whose spin is held (same rule as secular.kernel.leconte_pair)
        x = spin[i] / n if C > 0.0 else N / Om
        da += 4.0 * a**2 / (G * m[0] * m[1]) * K * (N * x - Na1)
        de += 11.0 * a * e / (G * m[0] * m[1]) * K * (Om_e * x - 18.0 / 11.0 * Ne1)
        if C > 0.0:
            dw.append(-2.0 * K / (C * n) * (Om * x - N))
            tsync.append(C * n**2 / (2.0 * K * Om) if K > 0 else np.inf)  # e-folding of the spin
        else:
            dw.append(0.0)
            tsync.append(0.0)
    return TidalRates(
        float(da),
        float(de),
        (float(dw[0]), float(dw[1])),
        (float(tsync[0]), float(tsync[1])),
        float(n * N / Om),  # Hut (1981): dOmega/dt = 0
    )


# --------------------------------------------------------------------------- tidal lag from
# the stellar structure (physics.tidal_lag: hurley2002)

E2_COEFF = 1.592e-9  # Hurley et al. (2002) eq. (43): E2 = 1.592e-9 M^2.84 (Zahn 1977 fit)


def hurley_k_over_t(
    mass: float,
    radius_rsun: float,
    m_env: float,
    turnover_yr: float,
    m_other: float,
    a_au: float,
    p_orb_yr: float,
    spin_rad_yr: float,
    reduction: str = "quadratic",
) -> float:
    """Tidal dissipation rate (k/T) [1/yr] of a star in a binary, Hurley, Tout & Pols (2002)
    MNRAS 329, 897, Sect. 2.3: equilibrium tide with turbulent convective damping in a
    convective envelope (eq. 30, after Rasio et al. 1996 and Zahn 1977) or the dynamical tide
    with radiative damping (eqs. 42-43, after Zahn 1975, 1977) for a radiative envelope.

    The turbulent viscosity is reduced when the tidal pumping period P_tid, 1/P_tid =
    |1/P_orb - 1/P_spin| (eq. 33), is shorter than twice the convective turnover time (the factor
    two: the tide comes around twice per period): quadratic, f = min(1, (P_tid / 2 tau)^2) (Goldreich & Nicholson
    1977; Hurley et al. 2002), or linear, f = min(1, P_tid / 2 tau) (Zahn 1966, 1989; used by
    Zahn & Bouchet 1989 for pre-main-sequence binaries)."""
    convective = m_env > 1e-3 * mass and turnover_yr > 0.0
    if convective:
        spin_freq = spin_rad_yr / (2.0 * np.pi)
        beat = abs(1.0 / p_orb_yr - spin_freq)
        p_tid = 1.0 / beat if beat > 0 else np.inf  # Hurley et al. (2002) eq. 33
        x = p_tid / (2.0 * turnover_yr)
        f_conv = min(1.0, x if reduction == "linear" else x * x)
        return 2.0 / 21.0 * f_conv / turnover_yr * m_env / mass
    a_rsun = a_au / units.RSUN_AU
    q = m_other / mass
    e2 = E2_COEFF * mass**2.84
    return 1.9782e4 * mass * radius_rsun**2 / a_rsun**5 * (1.0 + q) ** (5.0 / 6.0) * e2


def inertial_wave_q(
    mass: float, radius_rsun: float, m_env: float, r_env: float, spin_rad_yr: float
) -> float:
    """Modified tidal quality factor <Q'_IW> of inertial waves in the convective envelope
    (swingersim.inertial_waves: two-layer model of Ogilvie 2013 with a radiative core, n = 1.5
    polytrope when fully convective; Barker 2020, eq. 39); inf without an envelope."""
    from .inertial_waves import structure_c
    from .stellar import FULLY_CONVECTIVE_MSUN

    if m_env <= 1e-4 * mass or r_env <= 0.0:
        return np.inf
    alpha = max(1.0 - r_env / radius_rsun, 0.0)
    beta = max(1.0 - m_env / mass, 0.0)
    if mass <= FULLY_CONVECTIVE_MSUN:  # main-sequence stars below 0.35 Msun are fully
        alpha = beta = 0.0  # convective (the BSE envelope fit gives half the mass there)
    r_au = radius_rsun * units.RSUN_AU
    eps2 = spin_rad_yr**2 * r_au**3 / (units.G * mass)
    c = structure_c(alpha, beta)
    return 1.0 / (c * eps2) if c * eps2 > 0.0 else np.inf


def inertial_wave_k_over_t(
    mass: float,
    radius_rsun: float,
    m_env: float,
    r_env: float,
    m_other: float,
    a_au: float,
    spin_rad_yr: float,
) -> float:
    """(k/T) [1/yr] of the constant-time-lag model that damps the eccentricity at the rate of
    inertial waves, d ln e/dt = -(225/16) (M2/M) (R/a)^5 n / Q'_IW (Barker 2022, eq. 1);
    the constant-lag rate of a synchronous star is -(21/2) (k/T) q (1+q) (R/a)^8 (Hut 1981).
    Inertial waves exist only for tidal frequencies below 2 Omega, i.e. n <= 2 Omega."""
    n = np.sqrt(units.G * (mass + m_other) / a_au**3)
    if spin_rad_yr <= 0.0 or n > 2.0 * spin_rad_yr:
        return 0.0
    qp = inertial_wave_q(mass, radius_rsun, m_env, r_env, spin_rad_yr)
    if not np.isfinite(qp):
        return 0.0
    x = radius_rsun * units.RSUN_AU / a_au
    q = m_other / mass
    rate = 225.0 / 16.0 * q * x**5 * n / qp
    return rate / (10.5 * q * (1.0 + q) * x**8)


def lag_from_k_over_t(k_over_t_yr: float, mass: float, radius_rsun: float, apsidal_k2: float):
    """Constant time lag [s] (swingersim / Mignard convention, Love number 2 k2) equivalent to
    Hut's (1981) (k/T): T = R^3 / (G M tau_Hut) and tau_Hut = 2 lag."""
    r_au = radius_rsun * units.RSUN_AU
    tau_hut_yr = k_over_t_yr * r_au**3 / (units.G * mass * apsidal_k2)
    return 0.5 * tau_hut_yr * units.SECONDS_PER_YEAR


def with_structure_lags(
    cfg, props: list, pairs: list[tuple[int, int, float]], spins, masses=None
) -> list:
    """props with lag_s from the stellar structure (physics.tidal_lag: hurley2002) for the stars
    of the two-star binaries ``pairs`` = [(i, j, a_au), ...]; spins [rad/yr] per star; masses:
    the current masses [Msun] (mass transfer, winds), default the configured ones. A lag given
    in stars.<name>.lag_s is kept."""
    from dataclasses import replace

    if cfg.physics.tidal_lag != "hurley2002":
        return props
    out = list(props)
    eff = cfg.physics.tidal_efficiency

    def mass(k):
        return float(masses[k]) if masses is not None else cfg_mass(cfg, k)

    for i, j, a in pairs:
        p_orb = 2.0 * np.pi / np.sqrt(units.G * (mass(i) + mass(j)) / a**3)
        for s, o in ((i, j), (j, i)):
            spec = cfg.stars.get(f"m{s + 1}")
            if spec is not None and spec.lag_s is not None:
                continue
            p = props[s]
            if p.m_env is None or p.turnover_yr is None:
                raise ValueError(
                    "physics.tidal_lag: hurley2002 needs the envelope of the stars "
                    "(stellar_model.kind zams or sse)"
                )
            kt = hurley_k_over_t(
                mass(s),
                p.radius_rsun,
                p.m_env,
                p.turnover_yr,
                mass(o),
                a,
                p_orb,
                float(spins[s]),
                getattr(cfg.physics, "tidal_convective_reduction", "quadratic"),
            )
            if cfg.physics.inertial_waves and p.r_env is not None:
                kt += inertial_wave_k_over_t(
                    mass(s),
                    p.radius_rsun,
                    p.m_env,
                    p.r_env,
                    mass(o),
                    a,
                    float(spins[s]),
                )
            lag = eff * lag_from_k_over_t(kt, mass(s), p.radius_rsun, p.apsidal_k2)
            out[s] = replace(p, lag_s=lag)
    return out


def cfg_mass(cfg, k: int) -> float:
    return float(cfg.masses[f"m{k + 1}"])


def k2_rotating(cfg, k: int, k2: float, spin_rad_yr: float, mass: float, radius_rsun: float):
    """Apsidal-motion constant of star k corrected for its rotation (physics.k2_rotation):
    log k2 = log k2_static - 0.87 lambda, lambda = 2 V^2 / (3 g R) = (2/3) Omega^2 R^3 / (G M)
    (Claret 1999; Claret & Gimenez 2010, Sect. 3): a rotating star is more centrally condensed. A k2
    given in stars.<name>.apsidal_k2 is taken as it is."""
    if not getattr(cfg.physics, "k2_rotation", False) or not spin_rad_yr:
        return k2
    spec = cfg.stars.get(f"m{k + 1}")
    if spec is not None and spec.apsidal_k2 is not None:
        return k2
    r = radius_rsun * units.RSUN_AU
    lam = 2.0 / 3.0 * spin_rad_yr**2 * r**3 / (units.G * mass)
    return k2 * 10.0 ** (-0.87 * lam)
