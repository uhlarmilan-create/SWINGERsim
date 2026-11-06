"""Internal unit system and constants.

Internally swingersim works in solar masses, astronomical units and Julian years.
In these units the gravitational constant is the Gaussian value
G = (k * 365.25)^2 ~ 39.4769 AU^3 Msun^-1 yr^-2.

Angles are radians internally; the configuration file uses degrees.
"""

import math

GAUSS_K = 0.01720209895  # Gaussian gravitational constant [AU^1.5 Msun^-0.5 day^-1]
DAYS_PER_YEAR = 365.25
G = (GAUSS_K * DAYS_PER_YEAR) ** 2  # AU^3 / (Msun yr^2)

AU_KM = 149_597_870.7
C_KM_S = 299_792.458
SECONDS_PER_YEAR = DAYS_PER_YEAR * 86_400.0
C = C_KM_S * SECONDS_PER_YEAR / AU_KM  # speed of light [AU/yr] ~ 63241.08

RSUN_AU = 695_700.0 / AU_KM  # nominal solar radius (IAU 2015) [AU]
# solar mass [g]: IAU 2015 Resolution B3 nominal GM_sun = 1.3271244e20 m^3 s^-2 divided by
# CODATA 2018 G = 6.67430e-11 m^3 kg^-1 s^-2
MSUN_G = 1.98841e33
KMS_TO_AUYR = SECONDS_PER_YEAR / AU_KM  # 1 km/s in AU/yr (~0.2109)

DEG = math.pi / 180.0


def days_to_years(p_days: float) -> float:
    return p_days / DAYS_PER_YEAR


def period_to_sma(p_years: float, m_total: float) -> float:
    """Kepler's third law: semi-major axis [AU] from period [yr] and total mass [Msun]."""
    return (G * m_total * p_years**2 / (4.0 * math.pi**2)) ** (1.0 / 3.0)


def sma_to_period(a_au: float, m_total: float) -> float:
    """Period [yr] from semi-major axis [AU] and total mass [Msun]."""
    return 2.0 * math.pi * math.sqrt(a_au**3 / (G * m_total))
