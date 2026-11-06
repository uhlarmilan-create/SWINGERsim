"""One policy for quantities outside the grids of tabulated models.

* Model atmospheres and structure tables (limb and gravity darkening of Claret & Bloemen 2011,
  TLUSTY limb darkening of Reeve & Howarth 2016, apsidal constants of Claret 2023): a value
  outside the grid is clipped to the nearest edge before interpolating, with one warning per
  quantity and value that names the quantity, the value and the range of the grid.
* Stellar evolution tracks (MIST), which define the star itself: the nearest installed grid is
  used up to ``MIST_FEH_TOL`` dex beyond the installed range, with a warning; farther out it is
  an error.

Values inside a grid are never changed.
"""

from __future__ import annotations

import warnings

MIST_FEH_TOL = 0.1  # dex beyond the installed MIST grids that is still clipped
_WARNED: set = set()


def clip_to_grid(
    value: float,
    lo: float,
    hi: float,
    what: str,
    grid: str,
    tol: float | None = None,
    eps: float = 1e-9,
) -> float:
    """``value`` clipped to [lo, hi] (a warning once per (what, grid, value)); with ``tol`` a
    value farther than tol outside the range raises ValueError."""
    v = float(value)
    if lo - eps <= v <= hi + eps:
        return min(max(v, lo), hi)
    c = min(max(v, lo), hi)
    if tol is not None and abs(v - c) > tol + eps:
        raise ValueError(
            f"{what} = {v:g} is outside the {grid} ({lo:g} .. {hi:g}) by more than {tol:g}"
        )
    key = (what, grid, round(v, 6))
    if key not in _WARNED:
        _WARNED.add(key)
        warnings.warn(
            f"{what} = {v:g} is outside the {grid} ({lo:g} .. {hi:g}); using {c:g}",
            stacklevel=3,
        )
    return c


__all__ = ["MIST_FEH_TOL", "clip_to_grid"]
