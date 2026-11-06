"""Orbital inclinations against time in the joint fit of stage 2 (swingersim fit; fit.inclination).

A table of inclinations of an eclipsing pair measured at several epochs (light-curve solutions
of separate data sets, as Zasche et al. 2023 for V907 Sco) is compared with the inclination of
the pair's orbit in the N-body model at the same times (forward and backward from t = 0 of the
configuration). Only |cos i| is measured by eclipses, so the model inclination is folded to
0-90 deg, and so are measured values above 90 deg (i -> 180 - i). Columns: mid_hjd (or time), i_deg, sigma_deg[, limit][, source].

One-sided data, column ``limit``:
  (empty)     a measurement: residual (i - i_model) / sigma;
  upper       i_model <= i_deg (e.g. "i ~ 75 deg or less" from a light curve without eclipses);
  lower       i_model >= i_deg;
  no_eclipse  the pair did not eclipse at that epoch: i_model below the eclipse limit of the
              model, arccos((R1 + R2) / r_conj), with the radii of the model stars and the
              separation at the nearer conjunction; i_deg is ignored (may be empty);
  eclipse     the pair did eclipse: i_model above the eclipse limit.
A one-sided row contributes (bound - i_model) / sigma only when the model violates it (a hinge),
and nothing otherwise; sigma sets how sharply (a soft margin of a few tenths of a degree).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import units
from .config import Config

DPY = units.DAYS_PER_YEAR


@dataclass
class ObservedInclination:
    time: np.ndarray  # BJD
    inc: np.ndarray  # deg
    sigma: np.ndarray  # deg
    source: np.ndarray
    pair: str
    label: str
    limit: np.ndarray | None = None  # "", "upper", "lower", "no_eclipse", "eclipse"

    def __len__(self) -> int:
        return len(self.time)


def load(spec, base_dir) -> ObservedInclination:
    path = Path(spec.file)
    if not path.is_absolute():
        path = Path(base_dir or ".") / path
    with open(path, encoding="utf-8") as f:
        lines = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    t, i, s, src, lim = [], [], [], [], []
    for r in csv.DictReader(lines):
        r = {k.strip().lower(): (v or "").strip() for k, v in r.items() if k}
        t.append(float(r.get("time") or r.get("mid_hjd") or r.get("bjd") or r["hjd"]))
        kind = (r.get("limit") or "").lower().replace("-", "_").replace(" ", "_")
        if kind in ("none", "value", "measured"):
            kind = ""
        if kind not in LIMITS:
            raise ValueError(f"fit.inclination: unknown limit {kind!r} in {path} (use {LIMITS})")
        lim.append(kind)
        i.append(float(r["i_deg"]) if r.get("i_deg") else np.nan)
        if kind in ("", "upper", "lower") and not np.isfinite(i[-1]):
            raise ValueError(f"fit.inclination: i_deg missing in {path} (limit {kind!r})")
        if np.isfinite(i[-1]) and not 0.0 <= i[-1] <= 180.0:
            raise ValueError(f"fit.inclination: i_deg {i[-1]} outside 0-180 deg in {path}")
        if i[-1] > 90.0:  # eclipses measure |cos i|: the model is folded to 0-90, so is the data
            if kind in ("upper", "lower"):
                raise ValueError(
                    f"fit.inclination: a limit i {kind} {i[-1]} deg > 90 is ambiguous in {path} "
                    "(eclipses give |cos i|); give it as 180 - i with the opposite limit"
                )
            i[-1] = 180.0 - i[-1]
        s.append(max(float(r.get("sigma_deg") or 0.0), spec.min_sigma_deg))
        if s[-1] <= 0:
            raise ValueError(f"fit.inclination: sigma_deg must be > 0 in {path}")
        src.append(r.get("source") or "")
    if not t:
        raise ValueError(f"fit.inclination: no rows in {path}")
    return ObservedInclination(
        np.array(t),
        np.array(i),
        np.array(s),
        np.array(src, dtype=object),
        spec.pair,
        spec.label or path.stem,
        np.array(lim, dtype=object),
    )


LIMITS = ("", "upper", "lower", "no_eclipse", "eclipse")


def residuals(q: ObservedInclination, m) -> np.ndarray:
    """(bound - model) / sigma per row; m = model_inclination(..., with_limit=True), i.e. the
    array [i_model, i_eclipse_limit] of shape (2, n), or i_model alone. One-sided rows are zero
    when satisfied."""
    m = np.asarray(m, float)
    i_mod, i_lim = (m[0], m[1]) if m.ndim == 2 else (m, np.full(len(m), np.nan))
    lim = q.limit if q.limit is not None else np.full(len(q), "", dtype=object)
    bound = np.where(np.isin(lim, ("no_eclipse", "eclipse")), i_lim, q.inc)
    r = (bound - i_mod) / q.sigma
    upper = np.isin(lim, ("upper", "no_eclipse"))
    lower = np.isin(lim, ("lower", "eclipse"))
    r = np.where(upper, np.minimum(r, 0.0), r)
    r = np.where(lower, np.maximum(r, 0.0), r)
    return r


def model_i(m) -> np.ndarray:
    """The model inclination from what model_for returns (with or without the limit)."""
    m = np.asarray(m, float)
    return m[0] if m.ndim == 2 else m


def model_limit(m) -> np.ndarray:
    """The eclipse limit from what model_for returns (NaN where it was not computed)."""
    m = np.asarray(m, float)
    return m[1] if m.ndim == 2 else np.full(m.shape, np.nan)


def model_inclination(
    cfg: Config, pair: str, time_bjd, bjd_t0: float, with_limit: bool = False
) -> np.ndarray:
    """Inclination [deg, folded to 0-90] of the orbit of the pair at the given times; with
    with_limit also the eclipse limit arccos((R1 + R2) / r_conj) [deg] (array of shape (2, n))."""
    from .obsmodel import observables

    return observables(cfg, bjd_t0, incs=[(pair, time_bjd, with_limit)])["inc"][0]


def needs_limit(q: ObservedInclination) -> bool:
    """Whether the rows of q need the eclipse limit of the model (no_eclipse / eclipse)."""
    return bool(q.limit is not None and np.any(np.isin(q.limit, ("no_eclipse", "eclipse"))))


def model_for(d: dict, base_dir, incs, bjd_t0: float) -> list:
    """Model inclinations of every table of incs (one integration per direction for all)."""
    from .obsmodel import observables

    cfg = Config.model_validate(d)
    cfg._base_dir = base_dir
    return observables(cfg, bjd_t0, incs=[(q.pair, q.time, needs_limit(q)) for q in incs])["inc"]
