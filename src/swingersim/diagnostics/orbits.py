"""Orbit-level diagnostics."""

from __future__ import annotations

import numpy as np

from ..tides import apsidal_rate_gr


def mutual_inclination(h1: np.ndarray, h2: np.ndarray) -> np.ndarray:
    """Angle [rad] between orbital angular momentum unit vectors (works on arrays (..., 3))."""
    c = np.clip(np.sum(h1 * h2, axis=-1), -1.0, 1.0)
    return np.arccos(c)


def gr_periastron_advance_rate(m_total: float, a: float, e: float) -> float:
    """1PN apsidal advance d(omega)/dt [rad/yr] of an isolated binary (tides.apsidal_rate_gr,
    kept under this name for the diagnostics)."""
    return apsidal_rate_gr(m_total, a, e)
