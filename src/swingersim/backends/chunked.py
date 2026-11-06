"""Secular integration in chunks of output intervals.

A dense solution of a whole run keeps an interpolant for every integrator step. A tight binary
with tides over 10^6 yr takes millions of steps, so one run grew to gigabytes and an ensemble
on all cores exhausted the memory of the computer. Here the run is integrated chunk by chunk:
each chunk's dense solution gives the outputs and the extremes between them (every integrator
step plus a uniform grid, as before) and is then dropped, so the memory stays bounded
(secular_common.integrate_until_overflow). The extremes of the chunks are merged here."""

from __future__ import annotations

import math


def merge_extremes(total: dict, part: dict) -> None:
    """Fold the extremes of one chunk into those of the run so far."""
    for n, ex in part.items():
        cur = total.get(n)
        if cur is None:
            total[n] = dict(ex)
            continue
        if ex["e_max"] > cur["e_max"]:
            cur["e_max"], cur["t_e_max_yr"] = ex["e_max"], ex["t_e_max_yr"]
        if ex["q_min_au"] < cur["q_min_au"]:
            cur["q_min_au"], cur["t_q_min_yr"] = ex["q_min_au"], ex["t_q_min_yr"]
        for k, v in ex.items():  # first crossings of 1 - e thresholds: the earliest chunk
            if k.startswith("t_1me") and _nan(cur.get(k)) and not _nan(v):
                cur[k] = v


def merge_roche(total: dict, part: dict) -> None:
    for s, r in part.items():
        if s not in total or r["fill_max"] > total[s]["fill_max"]:
            total[s] = dict(r)


def _nan(v) -> bool:
    return v is None or (isinstance(v, float) and math.isnan(v))
