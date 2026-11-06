"""Light curves in the joint fit of stage 2 (swingersim fit; fit.lc).

Segments of a measured light curve (e.g. the windows around the eclipses of a third body, or
whole TESS sectors binned) are compared with the synthetic light curve of swingersim (N-body
positions with light time, limb-darkened spheres with ellipsoidal, reflection and beaming
terms; photometry.light_curve) at the measured times. The calibration of every segment is
solved linearly at each evaluation and does not enter the free parameters:

    measured(t) = a * model(t) + b + c (t - t_mid) + d (t - t_mid)^2

``a`` scales the eclipse depths (dilution by third light or a different normalisation), ``b``
the zero point (with ``a``: any extra light), ``c``, ``d`` slow trends of the photometry
(baseline_degree 0, 1 or 2); with scale: false only b (+ trends) are free and the depths must
come out of the stars as given. The weighted residuals of all points join those of the minima
and velocities; per-file error corrections (jitter) work as for the other data.

Table (CSV, # header lines): columns ``time`` (BJD_TDB), ``flux`` (normalised), ``sigma`` and
optionally ``segment`` (a label; otherwise segments are split at gaps longer than gap_d).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class ObservedLC:
    time: np.ndarray  # BJD_TDB [d]
    flux: np.ndarray
    sigma: np.ndarray
    segment: np.ndarray  # str label per point
    source: str  # label of the file (error corrections, report)
    band: str = "TESS"
    exposure_min: float = 0.0
    n_sub: int = 1
    baseline_degree: int = 1
    scale: bool = True
    shapes: str = "auto"  # auto, as_configured, spheres (swingersim.shapes)
    roche_threshold: float = 0.5
    mesh_res: tuple = (64, 128)
    fidelity: str = "mapped"
    coarse_res: tuple = (24, 48)
    delta: np.ndarray | None = None  # mapped: accurate - coarse model at the current solution
    red_noise: bool = False  # inflate the errors by the red-noise factor beta per segment

    def __len__(self) -> int:
        return len(self.time)

    def subset(self, keep, sigma=None) -> ObservedLC:
        s = self.sigma if sigma is None else sigma
        return ObservedLC(
            self.time[keep],
            self.flux[keep],
            s[keep],
            self.segment[keep],
            self.source,
            self.band,
            self.exposure_min,
            self.n_sub,
            self.baseline_degree,
            self.scale,
            self.shapes,
            self.roche_threshold,
            self.mesh_res,
            self.fidelity,
            self.coarse_res,
            None if self.delta is None else self.delta[keep],
            self.red_noise,
        )


def read_lc(path: str | Path, gap_d: float = 0.5) -> tuple:
    """(time, flux, sigma, segment) of a light-curve table."""
    with open(path, encoding="utf-8") as f:
        lines = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    t, fl, s, seg = [], [], [], []
    for r in csv.DictReader(lines):
        r = {k.strip().lower(): (v or "").strip() for k, v in r.items() if k}
        t.append(float(r["time"]))
        fl.append(float(r["flux"]))
        s.append(float(r.get("sigma") or "nan"))
        seg.append(r.get("segment") or "")
    t, fl, s = np.array(t), np.array(fl), np.array(s)
    seg = np.array(seg, dtype=object)
    order = np.argsort(t)
    t, fl, s, seg = t[order], fl[order], s[order], seg[order]
    if not np.all(seg != ""):
        k = np.concatenate([[0], np.cumsum(np.diff(t) > gap_d)])
        seg = np.array([f"s{j}" for j in k], dtype=object)
    if not len(t):
        raise ValueError(f"{path}: no light-curve points (columns time, flux, sigma)")
    return t, fl, s, seg.astype(str)


def load(spec, base_dir) -> ObservedLC:
    path = Path(spec.file)
    if not path.is_absolute():
        path = Path(base_dir or ".") / path
    t, fl, s, seg = read_lc(path, spec.gap_d)
    ok = np.isfinite(t) & np.isfinite(fl) & np.isfinite(s) & (s > 0)
    if spec.segments:
        ok &= np.isin(seg, spec.segments)
    return ObservedLC(
        t[ok],
        fl[ok],
        s[ok],
        seg[ok],
        spec.label or path.stem,
        spec.band,
        spec.exposure_min,
        spec.n_sub,
        spec.baseline_degree,
        spec.scale,
        spec.shapes,
        spec.roche_threshold,
        tuple(spec.mesh_resolution),
        spec.fidelity,
        tuple(spec.coarse_resolution),
        red_noise=spec.red_noise,
    )


def red_noise_beta(time, resid, windows_min=(20.0, 40.0, 60.0, 120.0), gap_factor: float = 3.0):
    """Red-noise factor beta >= 1 of a set of residuals by the time-averaging method (Pont,
    Zucker & Queloz 2006; Winn et al. 2008): the residuals are averaged in bins of N points
    spanning each window; beta_N = sigma_N / (sigma_1 / sqrt(N) * sqrt(M / (M - 1))) with M
    bins; beta = median over the windows (1 if there are too few points)."""
    t = np.asarray(time, float)
    r = np.asarray(resid, float)
    ok = np.isfinite(t) & np.isfinite(r)
    t, r = t[ok], r[ok]
    if len(t) < 20:
        return 1.0
    o = np.argsort(t)
    t, r = t[o], r[o]
    cad = float(np.median(np.diff(t)))
    if not cad > 0:
        return 1.0
    chunks = np.split(np.arange(len(t)), np.nonzero(np.diff(t) > gap_factor * cad)[0] + 1)
    s1 = float(np.std(r))
    if not s1 > 0:
        return 1.0
    betas = []
    for w in windows_min:
        n = round(w / 1440.0 / cad)
        if n < 2:
            continue
        means = [r[c[k : k + n]].mean() for c in chunks for k in range(0, len(c) - n + 1, n)]
        m = len(means)
        if m < 4:
            continue
        expected = s1 / np.sqrt(n) * np.sqrt(m / (m - 1.0))
        betas.append(float(np.std(means)) / expected)
    return max(1.0, float(np.median(betas))) if betas else 1.0


def model_flux(
    cfg, lc: ObservedLC, bjd_t0: float, resolution=None, table_phases: int | None = None
) -> np.ndarray:
    """Synthetic normalised flux of the configuration at the measured times (meshes with the
    given (theta, phi) resolution, default lc.mesh_res)."""
    from . import mesh
    from .photometry import light_curve

    t = lc.time - bjd_t0
    if np.any(t < 0):
        raise ValueError(
            "fit.lc: light-curve points before t = 0 of the configuration (fit.stage2.bjd_t0); "
            "choose bjd_t0 before the first point"
        )
    res = tuple(resolution or lc.mesh_res)
    old = (mesh.N_THETA, mesh.N_PHI, mesh.TABLE_PHASES)
    mesh.N_THETA, mesh.N_PHI = int(res[0]), int(res[1])
    if table_phases:
        mesh.TABLE_PHASES = int(table_phases)
    try:
        return light_curve(cfg, t, lc.band, exposure_min=lc.exposure_min, n_sub=lc.n_sub).flux
    finally:
        mesh.N_THETA, mesh.N_PHI, mesh.TABLE_PHASES = old


def lc_config(d: dict, lc: ObservedLC) -> tuple[dict, list[str]]:
    """The configuration dict of the light-curve model of lc (star shapes, swingersim.shapes)."""
    from .shapes import apply_policy

    return apply_policy(d, lc.shapes, lc.roche_threshold)


def needs_mesh(d_lc: dict) -> bool:
    from . import mesh
    from .config import Config

    return mesh.needed(Config.model_validate(d_lc))


def calibrate(lc: ObservedLC, model: np.ndarray):
    """Per segment the weighted linear least-squares calibration of the model to the data.
    Returns (calibrated model at every point, {segment: coefficients})."""
    out = np.full(len(lc), np.nan)
    coef = {}
    for s in np.unique(lc.segment):
        m = lc.segment == s
        t = lc.time[m]
        x = t - 0.5 * (t.min() + t.max())
        cols = [model[m]] if lc.scale else []
        cols.append(np.ones(m.sum()))
        for p in range(1, lc.baseline_degree + 1):
            cols.append(x**p)
        A = np.column_stack(cols)
        w = 1.0 / lc.sigma[m]
        y = lc.flux[m] - (0.0 if lc.scale else model[m])
        c, *_ = np.linalg.lstsq(A * w[:, None], y * w, rcond=None)
        out[m] = A @ c + (0.0 if lc.scale else model[m])
        coef[str(s)] = c.tolist()
    return out, coef


def apply_calibration(lc: ObservedLC, model: np.ndarray, coef: dict, t_mid=None) -> np.ndarray:
    """The calibrated model a model + b + c x + d x^2 at the points of lc with the coefficients
    of calibrate ({segment: coefficients}); x = t - t_mid of the segment (default: the centre
    of the points of lc in that segment, as in calibrate; a float, or {segment: t_mid}). NaN for
    segments without coefficients."""
    out = np.full(len(lc), np.nan)
    model = np.asarray(model, float)
    for s in np.unique(lc.segment):
        c = coef.get(str(s))
        if c is None:
            continue
        m = lc.segment == s
        t = lc.time[m]
        tm = t_mid.get(str(s)) if isinstance(t_mid, dict) else t_mid
        x = t - (0.5 * (t.min() + t.max()) if tm is None else tm)
        k = 0
        val = np.zeros(int(m.sum()))
        if lc.scale:
            val += c[0] * model[m]
            k = 1
        else:
            val += model[m]
        val += c[k]
        for p in range(1, lc.baseline_degree + 1):
            val += c[k + p] * x**p
        out[m] = val
    return out


def residuals(lc: ObservedLC, model: np.ndarray) -> np.ndarray:
    cal, _ = calibrate(lc, model)
    return (lc.flux - cal) / lc.sigma
