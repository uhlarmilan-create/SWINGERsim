"""Remove a strictly periodic oscillation from a TESS light curve before an eclipse fit.

    python tools/remove_oscillation.py DATA.csv MODEL.csv BJD_T0 PERIOD OUT.csv [--harmonics 2]

DATA: bjd_tdb, flux, flux_err (tools/fetch_tess.py). MODEL: the synthetic light curve of a
configuration that already follows the eclipses (columns t_day = BJD - BJD_T0, flux; e.g. the
light curves written by tools/squadra_figures.py). The model is calibrated linearly to the data
(depth scale, zero point, linear trend), a Fourier series of the given period and its harmonics
is fitted to the residuals, and this series is subtracted from the measured flux. The eclipses
themselves are not changed; the series is fitted to the whole sector, so its amplitude and phase
are set by the out-of-eclipse light and are not biased by the eclipse residuals.

Used for V994 Her: Zasche et al. (2023, Sect. 3) found an oscillation of about 2 % with the
period P_A - P_B = 0.663 d in the TESS residuals and removed it before their photodynamical fit.
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("data")
    ap.add_argument("model")
    ap.add_argument("bjd_t0", type=float)
    ap.add_argument("period", type=float)
    ap.add_argument("output")
    ap.add_argument("--harmonics", type=int, default=2)
    ap.add_argument(
        "--extra",
        type=int,
        default=0,
        help="further frequencies found in the residuals (pre-whitening, 0.3-3 d)",
    )
    a = ap.parse_args(argv)
    d = pd.read_csv(a.data, comment="#")
    m = pd.read_csv(a.model, comment="#")
    t = d["bjd_tdb"].to_numpy()
    f = d["flux"].to_numpy() / np.nanmedian(d["flux"])
    mod = np.interp(t - a.bjd_t0, m["t_day"], m["flux"])
    x = (t - t.mean()) / (t.max() - t.min())
    X = np.column_stack([mod, np.ones_like(t), x])
    c = np.linalg.lstsq(X, f, rcond=None)[0]
    r = f - X @ c
    oot = mod > np.percentile(mod, 60)  # out of eclipse: the series is fitted there only
    freqs = [k / a.period for k in range(1, a.harmonics + 1)]

    def series(fr):
        w = 2 * np.pi * (t - a.bjd_t0)
        return np.column_stack([g(f_ * w) for f_ in fr for g in (np.sin, np.cos)])

    for _ in range(a.extra):
        from astropy.timeseries import LombScargle

        F = series(freqs)
        q = np.linalg.lstsq(F[oot], r[oot], rcond=None)[0]
        fr, pw = LombScargle(t[oot], (r - F @ q)[oot]).autopower(
            minimum_frequency=1 / 3.0, maximum_frequency=1 / 0.3, samples_per_peak=20
        )
        freqs.append(float(fr[np.argmax(pw)]))
    F = series(freqs)
    q = np.linalg.lstsq(F[oot], r[oot], rcond=None)[0]
    osc = F @ q
    amp = [float(np.hypot(q[2 * k], q[2 * k + 1])) for k in range(len(freqs))]
    out = pd.DataFrame(
        {
            "bjd_tdb": t,
            "flux": f - osc,
            "flux_err": d["flux_err"].to_numpy() / np.nanmedian(d["flux"]),
        }
    )
    with open(a.output, "w") as fh:
        fh.write(
            f"# {a.data}: oscillation of period {a.period} d ({a.harmonics} harmonics) and {a.extra} further "
            f"periods {[round(1 / f_, 4) for f_ in freqs[a.harmonics :]]} d, semi-amplitudes "
            + ", ".join(f"{1e3 * x_:.2f}" for x_ in amp)
            + f" ppt) removed (tools/remove_oscillation.py); model {a.model}, depth scale {c[0]:.3f}\n"
            f"# rms of the residuals {1e3 * r.std():.2f} -> {1e3 * (r - osc).std():.2f} ppt "
            f"(out of eclipse {1e3 * r[oot].std():.2f} -> {1e3 * (r - osc)[oot].std():.2f})\n"
        )
        out.to_csv(fh, index=False, float_format="%.7f")
    print(
        f"{a.output}: periods {[round(1 / f_, 4) for f_ in freqs]} d, out of eclipse "
        f"{1e3 * r[oot].std():.2f} -> {1e3 * (r - osc)[oot].std():.2f} ppt; scale {c[0]:.3f}, semi-amplitudes {[round(1e3 * x_, 2) for x_ in amp]} ppt, "
        f"rms {1e3 * r.std():.2f} -> {1e3 * (r - osc).std():.2f} ppt"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
