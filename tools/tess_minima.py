"""Times of minima from TESS light curves (template fit per sector).

For every sector and eclipse type (primary, secondary) the eclipses are folded with a local
linear ephemeris, a mean eclipse profile is built and made symmetric about its own centre
(mirror symmetry), and every single eclipse is then fitted with that profile:

    flux(t) = a * T(t - t_min) + b + c (t - t_n)   (t_min, a, b, c free: depth, level, local trend)

The error of t_min comes from the covariance of the fit scaled by the residual rms. Eclipses
with less than 70 % of the window covered, or with the centre not covered on both sides, are
skipped. Times are BJD_TDB (TESS).

Usage:  python tools/tess_minima.py <csv files> --period P --t0 T0 [--sec-phase 0.5]
        [--half-width 0.08] [--clip 5] [--label TESS] [--out minima.csv]
The input files are written by tools/fetch_tess.py (bjd_tdb, flux, flux_err).
"""

from __future__ import annotations

import argparse
import re
from itertools import pairwise
from pathlib import Path

import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.optimize import least_squares, minimize_scalar


def load(path):
    """Times and flux of a light-curve CSV (columns bjd_tdb, flux; '#' comment lines allowed, e.g.
    the *_noosc.csv and *_sap.csv files), sorted, with isolated outliers removed."""
    import pandas as pd

    d = pd.read_csv(path, comment="#")
    d.columns = [c.strip().lower() for c in d.columns]
    t, f = d["bjd_tdb"].to_numpy(float), d["flux"].to_numpy(float)
    ok = np.isfinite(t) & np.isfinite(f)
    t, f = t[ok], f[ok]
    o = np.argsort(t)
    t, f = t[o], f[o]
    from scipy.ndimage import median_filter

    r = f - median_filter(f, size=15, mode="nearest")
    mad = 1.4826 * np.median(np.abs(r))
    keep = np.abs(r) < 10 * max(mad, 1e-5)  # isolated outliers (QLP), not the eclipses
    return t[keep], f[keep]


def template(t, f, tc, P, hw):
    """Symmetric mean profile around the eclipses at tc + n P (x in days)."""
    x = ((t - tc) / P + 0.5) % 1.0 - 0.5
    x *= P
    m = np.abs(x) < 1.6 * hw
    x, y = x[m], f[m]
    nb = 60
    edges = np.linspace(-1.6 * hw, 1.6 * hw, nb + 1)
    xc = 0.5 * (edges[1:] + edges[:-1])
    yb = np.array(
        [
            np.median(y[(x >= a) & (x < b)]) if np.any((x >= a) & (x < b)) else np.nan
            for a, b in pairwise(edges)
        ]
    )
    ok = np.isfinite(yb)
    sp = UnivariateSpline(xc[ok], yb[ok], s=0, k=3)

    def asym(c):
        z = np.linspace(0, hw, 80)
        return np.sum((sp(c + z) - sp(c - z)) ** 2)

    c = minimize_scalar(asym, bounds=(-0.3 * hw, 0.3 * hw), method="bounded").x
    z = np.linspace(-1.5 * hw, 1.5 * hw, 301)
    sym = 0.5 * (sp(c + z) + sp(c - z))
    return UnivariateSpline(z, sym, s=0, k=3), c


def locate(t, f, tc, P, hw, search=0.08):
    """Centre of the folded eclipse near tc (searched within +-search in phase)."""
    x = (((t - tc) / P + 0.5) % 1.0 - 0.5) * P
    m = np.abs(x) < search * P + hw
    if m.sum() < 20:
        return tc
    edges = np.linspace(-search * P - hw, search * P + hw, 121)
    xc = 0.5 * (edges[1:] + edges[:-1])
    yb = np.array(
        [
            np.median(f[m][(x[m] >= a) & (x[m] < b)])
            if np.any((x[m] >= a) & (x[m] < b))
            else np.nan
            for a, b in pairwise(edges)
        ]
    )
    ok = np.isfinite(yb) & (np.abs(xc) <= search * P)
    if not ok.any():
        return tc
    return tc + xc[ok][np.argmin(yb[ok])]


def fit_eclipses(t, f, tc, P, hw, clip):
    tmpl, c = template(t, f, tc, P, hw)
    tc = tc + c
    out = []
    n0, n1 = int(np.floor((t.min() - tc) / P)) - 1, int(np.ceil((t.max() - tc) / P)) + 1
    cad = np.median(np.diff(t))
    for n in range(n0, n1 + 1):
        tn = tc + n * P
        m = np.abs(t - tn) < hw
        if m.sum() < 0.7 * (2 * hw / cad):
            continue
        if (t[m] < tn - 0.3 * hw).sum() < 5 or (t[m] > tn + 0.3 * hw).sum() < 5:
            continue
        tt, ff = t[m], f[m]

        def model_res(p, tt, ff, tn=tn):
            return p[1] * tmpl(tt - p[0]) + p[2] + p[3] * (tt - tn) - ff

        r = least_squares(model_res, [tn, 1.0, 0.0, 0.0], args=(tt, ff))
        for _ in range(2):  # outlier clipping
            rr = model_res(r.x, tt, ff)
            s = 1.4826 * np.median(np.abs(rr))
            keep = np.abs(rr) < clip * max(s, 1e-6)
            tt, ff = tt[keep], ff[keep]
            r = least_squares(model_res, r.x, args=(tt, ff))
        J = r.jac
        rms = np.sqrt(np.sum(r.fun**2) / max(len(tt) - 4, 1))
        try:
            cov = np.linalg.inv(J.T @ J) * rms**2
            err = float(np.sqrt(cov[0, 0]))
        except np.linalg.LinAlgError:
            continue
        if abs(r.x[0] - tn) > 0.3 * hw or not np.isfinite(err):
            continue
        out.append((float(r.x[0]), err, len(tt), rms))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--period", type=float, required=True)
    ap.add_argument("--t0", type=float, required=True, help="a primary minimum (BJD)")
    ap.add_argument("--sec-phase", type=float, default=0.5)
    ap.add_argument("--half-width", type=float, default=0.08, help="window half-width [phase]")
    ap.add_argument("--clip", type=float, default=5.0)
    ap.add_argument("--binary", default="inner")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    P = a.period
    hw = a.half_width * P
    rows = []
    for fn in a.files:
        t, f = load(fn)
        sec = re.search(r"_s(\d+)_", Path(fn).name)
        sec = sec.group(1) if sec else "?"
        for kind, ph in (("primary", 0.0), ("secondary", a.sec_phase)):
            tc = a.t0 + ph * P
            # local ephemeris: shift tc by whole periods into the sector
            tc += np.round((t.mean() - tc) / P) * P
            tc = locate(t, f, tc, P, hw)  # eclipse-timing shifts (ETV, apsidal motion)
            try:
                res = fit_eclipses(t, f, tc, P, hw, a.clip)
            except Exception as e:  # noqa: BLE001
                print(f"  {fn} {kind}: {e}")
                continue
            for tm, err, npt, rms in res:
                rows.append((a.binary, kind, tm, err, sec, npt, rms))
            print(f"{Path(fn).name} {kind}: {len(res)} minima")
    rows.sort(key=lambda r: r[2])
    lines = ["binary,kind,bjd_tdb,sigma_d,cycle,ref"]
    for b, k, tm, err, sec, npt, _rms in rows:
        lines.append(f"{b},{k},{tm:.6f},{err:.6f},,TESS S{sec} template fit ({npt} pts)")
    text = "\n".join(lines) + "\n"
    if a.out:
        Path(a.out).write_text(text)
    else:
        print(text)


if __name__ == "__main__":
    main()
