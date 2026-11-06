"""Magnetism of low-mass stars: wind braking of the spins and the Applegate mechanism.

Magnetic braking (physics.magnetic_braking). The magnetised wind of a star with a convective
envelope removes spin angular momentum with the torque of Matt et al. (2015):

    T = -T0 (tau_cz / tau_cz,sun)^p (Omega / Omega_sun)^(p+1)   (unsaturated, Ro > Ro_sun / chi)
    T = -T0 chi^p (Omega / Omega_sun)                          (saturated),
    T0 = 9.5e30 erg (R / Rsun)^3.1 (M / Msun)^0.5,

with the Rossby number Ro = 1 / (Omega tau_cz), chi = 10, p = 2, Omega_sun = 2.6e-6 s^-1,
tau_cz,sun = 12.9 d, and the convective turnover time of Cranmer & Saar (2011, eq. 36):

    tau_cz = 314.24 exp[-(Teff / 1952.5 K) - (Teff / 6250 K)^18] + 0.002 d,

which vanishes for hot stars without convective envelopes. Over one operator-splitting segment
the moment of inertia is constant and dOmega/dt = T / I is integrated exactly (power law in the
unsaturated, exponential in the saturated regime). With tides, the braked spin of a star in a
close binary takes angular momentum from the orbit: magnetic braking shrinks close binaries.

Applegate mechanism. A magnetically active star whose quadrupole moment changes through its
activity cycle modulates the orbital period of its binary (Applegate 1992). The energy that
the cycle needs, relative to the energy the star radiates over one modulation period, decides
whether activity can explain an observed eclipse-timing variation: ``applegate_energy_ratio``
(thin-shell approximation of Tian et al. 2009, as used by Volschow et al. 2016) and its
inverse ``applegate_max_oc``, the largest O - C semi-amplitude the star can power.
"""

from __future__ import annotations

import math

import numpy as np

from . import units

OMEGA_SUN = 2.6e-6  # s^-1 (Matt et al. 2015, Table 1)
TAU_CZ_SUN_D = 12.9  # d
T0_SUN = 9.5e30  # erg
MSUN_G = units.MSUN_G  # IAU 2015 nominal solar mass [g]
RSUN_CM = 6.957e10  # IAU 2015 nominal solar radius [cm]
SEC_PER_YR = units.SECONDS_PER_YEAR  # Julian year [s]
DAY_S = 86400.0


def tau_cz_days(teff: float) -> float:
    """Convective turnover time [d], Cranmer & Saar (2011) eq. 36 (fit for 3300-7000 K)."""
    return 314.24 * math.exp(-(teff / 1952.5) - (teff / 6250.0) ** 18) + 0.002


def torque0(m: float, r: float, t0: float = T0_SUN) -> float:
    """T0 [erg] of Matt et al. (2015) eq. 8 for mass m [Msun] and radius r [Rsun]."""
    return t0 * r**3.1 * m**0.5


def omega_sat(teff: float, chi: float = 10.0) -> float:
    """Spin rate [s^-1] below which the torque is unsaturated (Ro = Ro_sun / chi)."""
    return chi * OMEGA_SUN * TAU_CZ_SUN_D / tau_cz_days(teff)


def torque(omega: float, m: float, r: float, teff: float, chi=10.0, p=2.0, t0=T0_SUN) -> float:
    """Spin-down torque [erg] (positive number) on a star spinning at omega [s^-1]."""
    t0m = torque0(m, r, t0)
    if omega >= omega_sat(teff, chi):
        return t0m * chi**p * (omega / OMEGA_SUN)
    return t0m * (tau_cz_days(teff) / TAU_CZ_SUN_D) ** p * (omega / OMEGA_SUN) ** (p + 1)


def spin_down(
    omega: float,
    m: float,
    r: float,
    teff: float,
    rg2: float,
    dt_yr: float,
    chi=10.0,
    p=2.0,
    t0=T0_SUN,
) -> float:
    """Spin rate after dt_yr of braking at constant moment of inertia I = rg2 M R^2; omega in
    s^-1 (any sign; the magnitude is braked). Exact in each regime, crossing at Omega_sat.
    dt_yr < 0: the spin rate dt_yr earlier (a backward run), the exact inverse."""
    if dt_yr == 0 or omega == 0.0:
        return omega
    sign = 1.0 if omega > 0 else -1.0
    w = abs(omega)
    inertia = rg2 * m * MSUN_G * (r * RSUN_CM) ** 2
    t0m = torque0(m, r, t0)
    w_sat = omega_sat(teff, chi)
    t = dt_yr * SEC_PER_YR
    if t < 0:
        return sign * _spin_earlier(w, -t, t0m, inertia, w_sat, teff, chi, p)
    if w > w_sat:  # saturated: exponential decay down to Omega_sat
        rate = t0m * chi**p / (inertia * OMEGA_SUN)
        t_cross = math.log(w / w_sat) / rate
        if t <= t_cross:
            return sign * w * math.exp(-rate * t)
        w, t = w_sat, t - t_cross
    k = t0m * (tau_cz_days(teff) / TAU_CZ_SUN_D) ** p / (inertia * OMEGA_SUN ** (p + 1))
    return sign * (w**-p + p * k * t) ** (-1.0 / p)


def _spin_earlier(w, t, t0m, inertia, w_sat, teff, chi, p) -> float:
    """Spin rate a time t [s] > 0 before the rate w: the braking laws of spin_down run back
    (faster in the past), unsaturated up to Omega_sat, then saturated."""
    k = t0m * (tau_cz_days(teff) / TAU_CZ_SUN_D) ** p / (inertia * OMEGA_SUN ** (p + 1))
    if w < w_sat:
        t_cross = (w**-p - w_sat**-p) / (p * k)
        if t <= t_cross:
            return (w**-p - p * k * t) ** (-1.0 / p)
        w, t = w_sat, t - t_cross
    rate = t0m * chi**p / (inertia * OMEGA_SUN)
    return w * math.exp(rate * t)


# ------------------------------------------------------------------ Applegate mechanism
def applegate_energy_ratio(
    m: float,
    r: float,
    teff: float,
    a_rsun: float,
    p_bin_d: float,
    dp_over_p: float,
    p_mod_yr: float,
) -> float:
    """Energy needed by the Applegate mechanism over the energy the active star radiates in one
    modulation period, thin-shell approximation of Tian et al. (2009):
    dE/E = 0.233 (M/Msun)^3 (R/Rsun)^-10 (T/6000 K)^-4 (a/Rsun)^4 (dP/s)^2 (P_mod/yr)^-1."""
    dp_s = dp_over_p * p_bin_d * DAY_S
    return 0.233 * m**3 * r**-10 * (teff / 6000.0) ** -4 * a_rsun**4 * dp_s**2 / p_mod_yr


def dp_over_p_from_oc(oc_semi_amplitude_s: float, p_mod_yr: float) -> float:
    """Relative period modulation from the O - C semi-amplitude K: dP/P = 2 pi (O - C) / P_mod
    with the full amplitude O - C = 2 K (Volschow et al. 2016, eq. 1 and Table 2)."""
    return 4.0 * math.pi * oc_semi_amplitude_s / (p_mod_yr * SEC_PER_YR)


def applegate_max_oc(
    m: float, r: float, teff: float, a_rsun: float, p_bin_d: float, p_mod_yr: float
) -> float:
    """Largest O - C semi-amplitude [s] the Applegate mechanism of the star can produce with
    a modulation period p_mod_yr: the one that needs all the energy it radiates (dE/E = 1)."""
    e1 = applegate_energy_ratio(m, r, teff, a_rsun, p_bin_d, 1.0, p_mod_yr)
    dpp = 1.0 / math.sqrt(e1)
    return dpp * p_mod_yr * SEC_PER_YR / (4.0 * math.pi)


# ------------------------------------------------------------------ activity in the timings
def activity_shift_yr(cfg, binary_stars: tuple[int, int], t_yr) -> np.ndarray | float:
    """O - C shift [yr] of a binary's eclipses by the activity cycles of its stars
    (physics.activity): sum of K sin(2 pi t / P_mod + phase). The Applegate mechanism changes
    the quadrupole moment of the active star and so the orbital period; the eclipse times
    follow the integral of the period change, a sinusoid of the cycle period."""
    act = getattr(cfg.physics, "activity", None) or {}
    out = 0.0 * np.asarray(t_yr, float)
    for k in binary_stars:
        spec = act.get(f"m{k + 1}")
        if spec is None:
            continue
        ph = 2 * math.pi * np.asarray(t_yr, float) / spec.p_mod_yr + math.radians(spec.phase_deg)
        out = out + spec.amplitude_s / SEC_PER_YR * np.sin(ph)
    return out


def activity_report(cfg) -> list[dict]:
    """Energy check of every activity cycle: the energy the Applegate mechanism needs over the
    energy the star radiates in one cycle (> 1: the cycle cannot be magnetic activity)."""
    from . import units
    from .hierarchy import build_tree
    from .stellar import star_properties

    act = getattr(cfg.physics, "activity", None) or {}
    if not act:
        return []
    root, masses = build_tree(cfg)
    out = []
    for n in root.nodes():
        if not (isinstance(n.left, int) and isinstance(n.right, int)):
            continue
        for k in (n.left, n.right):
            spec = act.get(f"m{k + 1}")
            if spec is None:
                continue
            p = star_properties(cfg, f"m{k + 1}", masses[k])
            teff = p.teff_k
            if teff is None:
                lum = p.luminosity_lsun or masses[k] ** 4
                teff = 5772.0 * lum**0.25 / math.sqrt(p.radius_rsun)
            a_rsun = n.elements.a / units.RSUN_AU
            p_bin_d = units.sma_to_period(n.elements.a, masses[n.left] + masses[n.right]) * 365.25
            dpp = dp_over_p_from_oc(spec.amplitude_s, spec.p_mod_yr)
            ratio = applegate_energy_ratio(
                masses[k], p.radius_rsun, teff, a_rsun, p_bin_d, dpp, spec.p_mod_yr
            )
            out.append(
                {
                    "star": f"m{k + 1}",
                    "binary": n.name,
                    "dE_over_E": ratio,
                    "max_amplitude_s": applegate_max_oc(
                        masses[k], p.radius_rsun, teff, a_rsun, p_bin_d, spec.p_mod_yr
                    ),
                    "possible": ratio < 1.0,
                }
            )
    return out


__all__ = [
    "activity_report",
    "activity_shift_yr",
    "applegate_energy_ratio",
    "applegate_max_oc",
    "dp_over_p_from_oc",
    "omega_sat",
    "spin_down",
    "tau_cz_days",
    "torque",
]
