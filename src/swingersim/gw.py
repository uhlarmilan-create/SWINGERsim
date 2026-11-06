"""Gravitational-wave emission of binaries (physics.gw), orbit-averaged (Peters 1964).

    da/dt = -(64/5) G^3 m1 m2 M / (c^5 a^3 (1 - e^2)^(7/2)) (1 + 73/24 e^2 + 37/96 e^4)
    de/dt = -(304/15) e G^3 m1 m2 M / (c^5 a^4 (1 - e^2)^(5/2)) (1 + 121/304 e^2)

(Peters 1964, Phys. Rev. 136, B1224, eqs. 5.6-5.7). The orbital period then changes at
dP/dt = (3/2) (P/a) da/dt; for the Hulse-Taylor pulsar PSR B1913+16 this is the general-relativity
prediction of -2.40e-12 (Weisberg & Huang 2016). The secular backends integrate the rates with
the orbits; the N-body backend applies them to every orbit of two single stars once per output
step (operator splitting; the step must be short compared to a / |da/dt|), keeping the
orientation and the mean anomaly. Relevant only for compact remnants and very close binaries.
"""

from __future__ import annotations

import math

from . import units


def peters_rates(m1: float, m2: float, a: float, e: float) -> tuple[float, float]:
    """(da/dt [AU/yr], de/dt [1/yr]) of gravitational-wave emission."""
    G, c = units.G, units.C
    M = m1 + m2
    k = G**3 * m1 * m2 * M / c**5
    e2 = e * e
    b = 1.0 - e2
    da = -64.0 / 5.0 * k / (a**3 * b**3.5) * (1 + 73.0 / 24.0 * e2 + 37.0 / 96.0 * e2 * e2)
    de = -304.0 / 15.0 * e * k / (a**4 * b**2.5) * (1 + 121.0 / 304.0 * e2)
    return da, de


def merger_time_circular_yr(m1: float, m2: float, a: float) -> float:
    """Time to coalescence of a circular binary [yr] (Peters 1964, eq. 5.10). Reference
    implementation (tests/validation)."""
    G, c = units.G, units.C
    return 5.0 / 256.0 * c**5 * a**4 / (G**3 * m1 * m2 * (m1 + m2))


def pdot(m1: float, m2: float, p_days: float, e: float) -> float:
    """dP/dt (dimensionless) of a binary of period p_days by gravitational waves. Reference
    implementation (tests/validation)."""
    M = m1 + m2
    p_yr = p_days / units.DAYS_PER_YEAR
    a = (units.G * M * (p_yr / (2 * math.pi)) ** 2) ** (1 / 3)
    da, _ = peters_rates(m1, m2, a, e)
    return 1.5 * p_yr / a * da


__all__ = ["merger_time_circular_yr", "pdot", "peters_rates"]
