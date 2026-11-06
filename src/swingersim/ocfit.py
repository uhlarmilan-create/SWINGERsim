"""Fitting of observed times of minima, stage 1 (swingersim fit).

Stage 1 is the analytic O - C model that the SQUADRA group fits today with OCFit
(Gajdos & Parimucha 2019), one pair at a time or both pairs with one outer orbit:

    T(n, j) = T0 + P (n + j/2) + Q n^2 + LiTE(t) + D_j(omega(n))

with the cycle number n, the minimum type j (0 primary, 1 secondary) and

* the light-travel-time effect of an outer orbit (Irwin 1952),

      LiTE = (a12 sin i3 / c) [ (1 - e3^2) / (1 + e3 cos nu3) sin(nu3 + w3) + e3 sin w3 ],

  where nu3 is the true anomaly of the outer orbit at the time of the minimum (periastron time
  T3, period P3); its semi-amplitude is K = (a12 sin i3 / c) sqrt(1 - e3^2 cos^2 w3);

* the apsidal motion of an eccentric pair: omega(n) = w0 + dw n, and the minimum occurs where
  the projected separation of the two (spherical) stars is smallest. At i = 90 deg that is the
  conjunction nu + omega = pi/2 (primary) or 3 pi/2 (secondary), and

      D_j = (P / 2 pi) [ M(nu_j) + omega - (pi/2 + j pi) ],

  the time from the mean to the true conjunction (M = mean anomaly). This is the exact form of
  the series of Gimenez & Bastero (1995) used by OCFit; for i < 90 deg the minimum of the
  projected separation is found numerically (the e cot^2 i terms of Gimenez & Garcia-Pelayo
  1983). The apsidal period is U = 2 pi P / dw.

With ``shared_outer`` both pairs have one outer orbit: P3, T3 and e3 in common, the pericentre
of the second pair opposite to the first (anti-phase), and the ratio of the a sin i values
gives the mass ratio of the pairs, m_A / m_B = (a sin i)_B / (a sin i)_A.

The fit is a weighted least-squares fit (scipy) started from a grid over the outer period;
errors come from the covariance matrix (scaled by the reduced chi^2 when it exceeds 1), or
from an MCMC run with emcee. The joint dynamical model of all pairs is stage 2 (dynfit.py).
"""

from __future__ import annotations

import contextlib
import csv
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

C_AU_PER_DAY = 173.1446326742403  # speed of light [AU/d]
DAY_S = 86400.0
TWO_PI = 2.0 * math.pi


# ------------------------------------------------------------------ data
@dataclass
class Minima:
    pair: np.ndarray  # str
    kind: np.ndarray  # 0 primary, 1 secondary
    time: np.ndarray  # BJD_TDB [d]
    sigma: np.ndarray  # [d]
    source: np.ndarray  # str
    method: np.ndarray  # str
    depth: np.ndarray
    width: np.ndarray
    header: dict = field(default_factory=dict)

    def select(self, mask) -> Minima:
        return Minima(
            *(
                getattr(self, k)[mask]
                for k in ("pair", "kind", "time", "sigma", "source", "method", "depth", "width")
            ),
            header=dict(self.header),
        )

    def __len__(self) -> int:
        return len(self.time)


def _float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


COLUMN_ALIASES = {
    "binary": "pair",
    "kind": "type",
    "bjd_tdb": "time",
    "bjd": "time",
    "sigma_d": "sigma",
    "ref": "source",
}


def read_minima(path: str | Path) -> Minima:
    """CSV table of minima. Header lines start with '#' (``# key: value``); columns pair, type
    (p/s or 1/2 or primary/secondary), time, sigma; optional source, method, depth, width.
    The columns of the measured-minima tables of swingersim gui (binary, kind, bjd_tdb, sigma_d,
    ref) are accepted as well."""
    header, rows = {}, []
    with open(path, encoding="utf-8") as f:
        lines = []
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                if ":" in s:
                    k, v = s[1:].split(":", 1)
                    header[k.strip()] = v.split("(")[0].strip()
                continue
            lines.append(s)
    reader = csv.DictReader(lines)
    for r in reader:
        row = {k.strip().lower(): (v or "").strip() for k, v in r.items() if k}
        for alias, name in COLUMN_ALIASES.items():  # the table of swingersim gui / data/observed
            if alias in row and name not in row:
                row[name] = row.pop(alias)
        rows.append(row)
    need = {"pair", "type", "time", "sigma"}
    if not rows or not need <= set(rows[0]):
        raise ValueError(f"{path}: the table needs the columns {sorted(need)}")
    kinds = {"p": 0, "1": 0, "primary": 0, "pri": 0, "s": 1, "2": 1, "secondary": 1, "sec": 1}
    bad = sorted({r["type"] for r in rows if r["type"].lower() not in kinds})
    if bad:
        raise ValueError(f"{path}: unknown minimum types {bad} (p or s)")
    col = lambda k, d="": np.array([r.get(k, d) for r in rows])
    return Minima(
        pair=col("pair"),
        kind=np.array([kinds[r["type"].lower()] for r in rows]),
        time=np.array([float(r["time"]) for r in rows]),
        sigma=np.array([float(r["sigma"]) for r in rows]),
        source=col("source", "unknown"),
        method=col("method", "individual"),
        depth=np.array([_float(r.get("depth")) for r in rows]),
        width=np.array([_float(r.get("width")) for r in rows]),
        header=header,
    )


def _light_time(t, star, kind: str):
    """Light travel time [d] of the geocentre (kind "barycentric") or towards the Sun's centre
    ("heliocentric"): the projection of the Earth's barycentric (heliocentric) position on the
    direction to the star / c, as astropy's Time.light_travel_time for an observer at the
    geocentre - but without its Earth-rotation frames, which need the IERS table of UT1 - UTC
    and fail for dates beyond its predictions (minima predicted decades ahead) or offline."""
    from astropy import constants as const
    from astropy import units as u
    from astropy.coordinates import get_body_barycentric

    pos = get_body_barycentric("earth", t)
    if kind == "heliocentric":
        pos = pos - get_body_barycentric("sun", t)
    n = star.icrs.cartesian.xyz.to_value()
    r = pos.xyz.to_value(u.m)
    return (np.tensordot(n, r, axes=1) / const.c.to_value(u.m / u.s)) / 86400.0


JD_UTC_1960 = 2436934.5  # 1960 Jan 1.0: UTC (and astropy's table of TAI - UTC) starts here


def delta_t_s(jd) -> np.ndarray:
    """Delta T = TT - UT [s] for Julian dates before 1960 (1600-1961): the polynomial expressions
    of Espenak & Meeus (2006, Five Millennium Canon of Solar Eclipses, NASA/TP-2006-214141;
    eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html). y is the decimal (Julian) year."""
    y = 2000.0 + (np.asarray(jd, float) - 2451545.0) / 365.25
    if np.any(y < 1600.0) or np.any(y > 1961.0):
        raise ValueError("delta_t_s: Delta T polynomials are used for 1600-1961 only")
    out = np.empty_like(y)
    pieces = [
        (1600, 1700, lambda t: 120 - 0.9808 * t - 0.01532 * t**2 + t**3 / 7129, 1600),
        (1700, 1800, lambda t: 8.83 + 0.1603 * t - 0.0059285 * t**2 + 0.00013336 * t**3
         - t**4 / 1174000, 1700),  # fmt: skip
        (1800, 1860, lambda t: 13.72 - 0.332447 * t + 0.0068612 * t**2 + 0.0041116 * t**3
         - 0.00037436 * t**4 + 0.0000121272 * t**5 - 0.0000001699 * t**6
         + 0.000000000875 * t**7, 1800),  # fmt: skip
        (1860, 1900, lambda t: 7.62 + 0.5737 * t - 0.251754 * t**2 + 0.01680668 * t**3
         - 0.0004473624 * t**4 + t**5 / 233174, 1860),  # fmt: skip
        (1900, 1920, lambda t: -2.79 + 1.494119 * t - 0.0598939 * t**2 + 0.0061966 * t**3
         - 0.000197 * t**4, 1900),  # fmt: skip
        (1920, 1941, lambda t: 21.20 + 0.84493 * t - 0.076100 * t**2 + 0.0020936 * t**3, 1920),
        (1941, 1961.0001, lambda t: 29.07 + 0.407 * t - t**2 / 233 + t**3 / 2547, 1950),
    ]
    for lo, hi, f, y0 in pieces:
        m = (y >= lo) & (y < hi)
        out[m] = f(y[m] - y0)
    return out


@contextlib.contextmanager
def _quiet_erfa():
    """TT -> TDB in astropy evaluates a microsecond-level diurnal term with an approximate UT and
    warns 'dubious year' for dates before 1960; the warning concerns that term only."""
    import warnings

    from erfa import ErfaWarning

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ErfaWarning)
        yield


def _utc_time(jd):
    """astropy Time of geocentric 'UTC' Julian dates. From 1960 on: UTC (astropy/ERFA, including
    the rubber seconds of 1960-1972). Before 1960 UTC did not exist and the times are Universal
    Time; ERFA would then assume TT - UTC = 32.184 s, wrong by up to ~35 s around 1900. Those
    times are converted with TT = UT + Delta T (delta_t_s)."""
    from astropy.time import Time

    jd = np.atleast_1d(np.asarray(jd, float))
    old = jd < JD_UTC_1960
    if not old.any():
        return Time(jd, format="jd", scale="utc")
    tt = jd.copy()
    tt[old] = jd[old] + delta_t_s(jd[old]) / 86400.0
    if not (~old).any():
        return Time(tt, format="jd", scale="tt")
    tt[~old] = Time(jd[~old], format="jd", scale="utc").tt.jd
    return Time(tt, format="jd", scale="tt")


def to_bjd_tdb(times, system: str, ra_deg: float, dec_deg: float):
    """Convert HJD_UTC or JD_UTC (geocentric) to BJD_TDB (astropy ephemeris and time scales).
    Times before 1960 are Universal Time and get TT = UT + Delta T (see _utc_time)."""
    if system == "BJD_TDB":
        return np.asarray(times, float)
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    star = SkyCoord(ra_deg * u.deg, dec_deg * u.deg)
    jd = np.asarray(times, float)
    scalar = jd.ndim == 0
    with _quiet_erfa():
        t = _utc_time(jd)
        if system == "HJD_UTC":  # remove the heliocentric correction, then add the barycentric one
            t = _utc_time(np.atleast_1d(jd) - _light_time(t.tdb, star, "heliocentric"))
        out = t.tdb.jd + _light_time(t.tdb, star, "barycentric")
    return out[0] if scalar else out


def from_bjd_tdb(bjd, ra_deg: float, dec_deg: float):
    """Inverse of to_bjd_tdb for JD_UTC: the geocentric UTC Julian date at which an event with
    the given BJD_TDB is seen (three fixed-point iterations; error far below 1 ms); before
    1960 the result is Universal Time (TT - Delta T)."""
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    star = SkyCoord(ra_deg * u.deg, dec_deg * u.deg)
    target = np.asarray(bjd, float)
    scalar = target.ndim == 0
    target = np.atleast_1d(target)
    jd = target.copy()
    with _quiet_erfa():
        for _ in range(3):
            t = _utc_time(jd)
            b = t.tdb.jd + _light_time(t.tdb, star, "barycentric")
            jd = jd - (b - target)
    return jd[0] if scalar else jd


# ------------------------------------------------------------------ model pieces
def kepler_E(M, e: float) -> np.ndarray:
    """Eccentric anomaly of the mean anomalies M [rad] (Newton; the one solver of stage 1 and
    stage 2): start M + e sin M (pi for e >= 0.8), until the step is below 1e-14."""
    M = np.asarray(M, float)
    E = M + e * np.sin(M) if e < 0.8 else np.full_like(M, math.pi)
    for _ in range(60):
        dE = (E - e * np.sin(E) - M) / (1.0 - e * np.cos(E))
        E = E - dE
        if np.max(np.abs(dE), initial=0.0) < 1e-14:
            break
    return E


def true_anomaly(t, T3, P3, e3):
    M = TWO_PI * (np.asarray(t, float) - T3) / P3
    E = kepler_E(np.mod(M, TWO_PI), e3)
    return 2.0 * np.arctan2(math.sqrt(1 + e3) * np.sin(E / 2), math.sqrt(1 - e3) * np.cos(E / 2))


def lite(t, asini_au, e3, w3, T3, P3):
    """Light-time effect [d] of Irwin (1952) (same definition as OCFit)."""
    nu = true_anomaly(t, T3, P3, e3)
    return (
        asini_au
        / C_AU_PER_DAY
        * ((1 - e3**2) / (1 + e3 * np.cos(nu)) * np.sin(nu + w3) + e3 * math.sin(w3))
    )


def lite_semi_amplitude_d(asini_au, e3, w3):
    return asini_au / C_AU_PER_DAY * math.sqrt(1 - e3**2 * math.cos(w3) ** 2)


def mass_function(asini_au, P3_d):
    """f(m3) = (m3 sin i3)^3 / (m12 + m3)^2 = (a12 sin i3)^3 / P3^2 [Msun], AU and yr."""
    return asini_au**3 / (P3_d / 365.25) ** 2


def mean_from_true(nu, e):
    E = 2.0 * np.arctan2(math.sqrt(1 - e) * np.sin(nu / 2), math.sqrt(1 + e) * np.cos(nu / 2))
    return E - e * np.sin(E)


def _wrap(x):
    return (x + math.pi) % TWO_PI - math.pi


def minimum_anomaly(omega, e, kind, inc_rad=math.pi / 2):
    """True anomaly of the minimum of the projected separation near conjunction j."""
    nu = math.pi / 2 - omega + kind * math.pi
    if abs(inc_rad - math.pi / 2) < 1e-9 or e == 0.0:
        return nu
    s2 = math.sin(inc_rad) ** 2

    def d2(v):  # projected separation squared / a^2
        r = (1 - e * e) / (1 + e * np.cos(v))
        return r * r * (1 - s2 * np.sin(v + omega) ** 2)

    h = 1e-4
    for _ in range(20):
        f1 = (d2(nu + h) - d2(nu - h)) / (2 * h)
        f2 = (d2(nu + h) - 2 * d2(nu) + d2(nu - h)) / (h * h)
        step = np.where(f2 > 0, f1 / np.where(f2 > 0, f2, 1.0), 0.0)
        nu = nu - np.clip(step, -0.2, 0.2)
        if np.all(np.abs(step) < 1e-10):
            break
    return nu


def apsidal_shift(n, kind, P, e, w0, dw, inc_rad=math.pi / 2):
    """D_j [d]: time of the minimum after the mean conjunction (apsidal motion)."""
    n = np.asarray(n, float)
    kind = np.asarray(kind)
    omega = w0 + dw * (n + kind / 2.0)
    nu = minimum_anomaly(omega, e, kind, inc_rad)
    M = mean_from_true(nu, e)
    return P / TWO_PI * _wrap(M + omega - (math.pi / 2 + kind * math.pi))


# ------------------------------------------------------------------ parameters
LITE_KEYS = ("P3", "T3", "e3", "w3")


@dataclass
class PairModel:
    name: str
    terms: list[str]
    inc_rad: float = math.pi / 2

    def keys(self, shared: bool) -> list[str]:
        k = ["T0", "P"]
        if "quadratic" in self.terms:
            k.append("Q")
        if "apsidal" in self.terms:
            k += ["e", "w0", "dw"]
        if "lite" in self.terms:
            k.append("asini")
            if not shared:
                k += list(LITE_KEYS)
        return k


@dataclass
class Model:
    pairs: list[PairModel]
    shared: bool
    names: list[str] = field(default_factory=list)  # flat parameter names "A.T0", "out.P3"
    t_ref: float = 0.0  # times inside the fit are relative to t_ref (numerical precision)

    def __post_init__(self):
        self.names = []
        if self.shared:
            self.names += [f"out.{k}" for k in LITE_KEYS]
        for p in self.pairs:
            self.names += [f"{p.name}.{k}" for k in p.keys(self.shared)]

    def unpack(self, x) -> dict[str, dict[str, float]]:
        d = dict(zip(self.names, x, strict=True))
        out = {}
        for p in self.pairs:
            q = {k.split(".", 1)[1]: v for k, v in d.items() if k.startswith(p.name + ".")}
            if "lite" in p.terms and self.shared:
                q.update({k: d[f"out.{k}"] for k in LITE_KEYS})
                if p is not self.pairs[0]:
                    q["w3"] = d["out.w3"] + math.pi
            out[p.name] = q
        return out

    def predict(self, x, pair, n, kind, t_obs=None):
        q = self.unpack(x)[pair.name]
        t = q["T0"] + q["P"] * (n + kind / 2.0)
        if "quadratic" in pair.terms:
            t = t + q["Q"] * n * n
        if "apsidal" in pair.terms:
            t = t + apsidal_shift(n, kind, q["P"], q["e"], q["w0"], q["dw"], pair.inc_rad)
        if "lite" in pair.terms:
            tt = t if t_obs is None else t_obs
            t = t + lite(tt, q["asini"], q["e3"], q["w3"], q["T3"], q["P3"])
        return t

    def components(self, x, pair, n, kind, t_obs):
        q = self.unpack(x)[pair.name]
        c = {}
        if "quadratic" in pair.terms:
            c["quadratic"] = q["Q"] * n * n
        if "apsidal" in pair.terms:
            c["apsidal"] = apsidal_shift(n, kind, q["P"], q["e"], q["w0"], q["dw"], pair.inc_rad)
        if "lite" in pair.terms:
            c["lite"] = lite(t_obs, q["asini"], q["e3"], q["w3"], q["T3"], q["P3"])
        return c


def bounds_for(name: str, span: float):
    k = name.split(".", 1)[1]
    return {
        "e3": (0.0, 0.95),
        "e": (0.0, 0.95),
        "asini": (0.0, 1e4),
        "P3": (1.0, 1e3 * span),
        "dw": (-0.5, 0.5),
    }.get(k, (-np.inf, np.inf))


# ------------------------------------------------------------------ fit
@dataclass
class PairData:
    model: PairModel
    n: np.ndarray
    kind: np.ndarray
    t: np.ndarray
    sigma: np.ndarray
    source: np.ndarray
    T0: float
    P: float


@dataclass
class FitResult:
    names: list[str]
    x: np.ndarray
    err: np.ndarray
    cov: np.ndarray
    chi2: float
    dof: int
    model: Model
    data: list[PairData]
    samples: np.ndarray | None = None
    scale: dict = field(default_factory=dict)
    outliers: list = field(default_factory=list)  # minima left out by stage1.clip_sigma
    best_fit: dict | None = None  # emcee: the least-squares solution (x, chi2, err)
    tau_max: float = float("nan")  # emcee: largest autocorrelation time [steps]

    def params(self) -> dict[str, float]:
        return {n: self.absolute(n, v) for n, v in zip(self.names, self.x, strict=True)}

    def absolute(self, name: str, v: float) -> float:
        """Parameter value with times back in the time system of the data."""
        return v + self.model.t_ref if name.endswith((".T0", ".T3")) else v

    def errors(self) -> dict[str, float]:
        return dict(zip(self.names, self.err, strict=True))


def assign_epochs(t, kind, T0, P):
    n = np.round((t - T0) / P - kind / 2.0)
    return n


def _residuals(x, model: Model, data: list[PairData]):
    out = []
    for d in data:
        out.append((d.t - model.predict(x, d.model, d.n, d.kind, d.t)) / d.sigma)
    return np.concatenate(out)


def _linear_start(d: PairData):
    A = np.vstack([np.ones_like(d.n), d.n + d.kind / 2.0]).T / d.sigma[:, None]
    sol, *_ = np.linalg.lstsq(A, d.t / d.sigma, rcond=None)
    return sol


def _sign(kind):
    return 2.0 * np.asarray(kind, float) - 1.0  # -1 primary, +1 secondary


def _apsidal_grid(d: PairData, span_cycles: float):
    """First guess of the apsidal motion from the anti-phase of primary and secondary minima:
    D_s - D_p ~ (2 e P / pi) cos(omega) (first order in e). Returns (e, w0, dw)."""
    T0, P = _linear_start(d)
    r = d.t - (T0 + P * (d.n + d.kind / 2.0))
    s = _sign(d.kind)
    best = (np.inf, None)
    for U_cyc in np.exp(
        np.linspace(
            math.log(max(span_cycles / 50, 10.0)), math.log(1e4 * max(span_cycles, 1.0)), 1500
        )
    ):
        dw = TWO_PI / U_cyc
        ph = dw * d.n
        A = (
            np.vstack([np.ones_like(d.n), d.n, d.kind, d.n**2, s * np.cos(ph), s * np.sin(ph)]).T
            / d.sigma[:, None]
        )
        sol, *_ = np.linalg.lstsq(A, r / d.sigma, rcond=None)
        chi = float(np.sum((A @ sol - r / d.sigma) ** 2))
        if chi < best[0]:
            best = (chi, (dw, sol))
    dw, sol = best[1]
    amp = math.hypot(sol[4], sol[5])
    e = min(max(math.pi * amp / P, 0.01), 0.8)
    # D_s - D_p = (2 e P / pi) cos(w0 + dw n)  =>  2 * amp cos(dw n + phi)
    w0 = math.atan2(-sol[5], sol[4]) % TWO_PI
    return e, w0, dw, best[0]


def _lite_grid(data: list[PairData], prange, n_best=3):
    """Best outer periods from a grid of circular light-time fits on the residuals of the
    linear ephemeris; an anti-phase (apsidal) term is fitted along for apsidal pairs."""
    lo, hi = prange
    res = []
    for d in data:
        T0, P = _linear_start(d)
        res.append((d, d.t - (T0 + P * (d.n + d.kind / 2.0))))
    periods = np.exp(np.linspace(math.log(lo), math.log(hi), 3000))
    chis = np.zeros(len(periods))
    for i, P3 in enumerate(periods):
        chi = 0.0
        for d, r in res:
            ph = TWO_PI * d.t / P3
            cols = [np.ones_like(d.t), d.n, np.sin(ph), np.cos(ph)]
            if "apsidal" in d.model.terms:
                s = _sign(d.kind)
                cols += [d.kind, s * d.n, s * d.n**2, s * d.n**3]
            A = np.vstack(cols).T / d.sigma[:, None]
            sol, *_ = np.linalg.lstsq(A, r / d.sigma, rcond=None)
            chi += float(np.sum((A @ sol - r / d.sigma) ** 2))
        chis[i] = chi
    # local minima, best first
    idx = [i for i in range(1, len(chis) - 1) if chis[i] <= chis[i - 1] and chis[i] <= chis[i + 1]]
    idx = sorted(idx, key=lambda i: chis[i])[:n_best] or [int(np.argmin(chis))]
    return [float(periods[i]) for i in idx]


def _diff_steps(names, x):
    """Relative finite-difference steps for scipy (absolute steps for the times: a relative step
    of a Julian date of 2.46e6 would be days)."""
    absolute = {"T0": 1e-6, "T3": 1e-3, "P": 1e-10, "P3": 1e-4}
    out = []
    for n, v in zip(names, x, strict=True):
        h = absolute.get(n.split(".", 1)[1])
        out.append(h / max(1.0, abs(v)) if h else 1.5e-8)
    return np.array(out)


def fit_stage1(model: Model, data: list[PairData], prange=None):
    from scipy.optimize import least_squares

    tmin = min(d.t.min() for d in data)
    tmax = max(d.t.max() for d in data)
    span = tmax - tmin
    has_lite = any("lite" in d.model.terms for d in data)
    if prange is None:
        gaps = np.diff(np.sort(np.concatenate([d.t for d in data])))
        gaps = gaps[gaps > 0]
        lo = max(2.0 * float(np.median(gaps)) if len(gaps) else 1.0, 4.0 * max(d.P for d in data))
        prange = (lo, 3.0 * span)
    P3_list = _lite_grid(data, prange) if has_lite else [None]
    aps = {d.model.name: _apsidal_grid(d, span / d.P) for d in data if "apsidal" in d.model.terms}
    lin = {d.model.name: _linear_start(d) for d in data}

    lb = np.array([bounds_for(n, span)[0] for n in model.names])
    ub = np.array([bounds_for(n, span)[1] for n in model.names])

    def start(P3, e3, w3, flip):
        x = []
        for nm in model.names:
            pre, key = nm.split(".", 1)
            if key in LITE_KEYS:
                v = {"P3": P3, "T3": tmin + 0.25 * P3, "e3": e3, "w3": w3}[key]
            elif key == "T0":
                v = lin[pre][0]
            elif key == "P":
                v = lin[pre][1]
            elif key == "Q":
                v = 0.0
            elif key == "asini":
                d = next(dd for dd in data if dd.model.name == pre)
                r = d.t - (lin[pre][0] + lin[pre][1] * (d.n + d.kind / 2.0))
                v = max(0.5 * (np.percentile(r, 95) - np.percentile(r, 5)) * C_AU_PER_DAY, 1e-4)
            elif key == "e":
                v = aps[pre][0]
            elif key == "w0":
                v = (aps[pre][1] + (math.pi if flip else 0.0)) % TWO_PI
            elif key == "dw":
                v = aps[pre][2]
            x.append(v)
        return np.clip(np.array(x, float), lb + 1e-9, ub - 1e-9)

    starts = []
    flips = (False, True) if aps else (False,)
    for P3 in P3_list:
        for e3 in (0.1, 0.4, 0.7) if has_lite else (0.0,):
            for w3 in (0.0, math.pi / 2, math.pi, 1.5 * math.pi) if has_lite else (0.0,):
                for fl in flips:
                    starts.append(start(P3, e3, w3, fl))
    first = []
    for x0 in starts:  # short runs from every start, then the best ones to convergence
        try:
            r = least_squares(
                _residuals,
                x0,
                args=(model, data),
                bounds=(lb, ub),
                x_scale="jac",
                max_nfev=60,
                diff_step=_diff_steps(model.names, x0),
            )
        except (ValueError, FloatingPointError):
            continue
        first.append(r)
    if not first:
        raise RuntimeError("the fit did not converge from any start")
    first.sort(key=lambda r: r.cost)
    best = None
    for r0 in first[:4]:
        r = least_squares(
            _residuals,
            r0.x,
            args=(model, data),
            bounds=(lb, ub),
            x_scale="jac",
            max_nfev=5000,
            ftol=1e-12,
            xtol=1e-12,
            gtol=1e-12,
            diff_step=_diff_steps(model.names, r0.x),
        )
        if best is None or r.cost < best.cost:
            best = r
    x = best.x.copy()
    for i, nm in enumerate(model.names):  # angles into [0, 2 pi)
        if nm.endswith((".w3", ".w0")):
            x[i] = x[i] % TWO_PI
    chi2 = float(np.sum(_residuals(x, model, data) ** 2))
    dof = sum(len(d.t) for d in data) - len(x)
    from .dynfit import covariance

    # (J^T J)^-1 from the singular values of the column-scaled Jacobian (as stage 2): inverting
    # J^T J squares the condition number (T0 of 2.4e6 d next to e of 0.1)
    cov = covariance(np.asarray(best.jac, float), names=list(model.names))
    red = chi2 / max(dof, 1)
    if red > 1:
        cov = cov * red
    err = np.sqrt(np.clip(np.diag(cov), 0, None))
    return FitResult(list(model.names), x, err, cov, chi2, dof, model, data)


def run_emcee(res: FitResult, walkers=48, steps=4000, seed=1):
    import emcee

    span = max(d.t.max() for d in res.data) - min(d.t.min() for d in res.data)
    lb = np.array([bounds_for(n, span)[0] for n in res.names])
    ub = np.array([bounds_for(n, span)[1] for n in res.names])
    red = max(res.chi2 / max(res.dof, 1), 1.0)

    def lnp(x):
        if np.any(x < lb) or np.any(x > ub):
            return -np.inf
        return -0.5 * float(np.sum(_residuals(x, res.model, res.data) ** 2)) / red

    rng = np.random.default_rng(seed)
    ndim = len(res.x)
    scale = np.where(res.err > 0, res.err, 1e-8)
    p0 = res.x + 1e-2 * scale * rng.normal(size=(walkers, ndim))
    p0 = np.clip(p0, lb + 1e-12, ub - 1e-12)
    s = emcee.EnsembleSampler(walkers, ndim, lnp)
    s.random_state = np.random.RandomState(seed).get_state()  # reproducible chain
    s.run_mcmc(p0, steps, progress=False)
    chain = s.get_chain(discard=steps // 2, thin=max(steps // 400, 1), flat=True)
    # the result describes the samples: x their median, err half the 16-84 % range, cov their
    # covariance and chi^2 that of the median; the least-squares solution is kept in best_fit
    res.best_fit = {"x": res.x.copy(), "chi2": res.chi2, "err": res.err.copy()}
    res.samples = chain
    res.x = np.median(chain, axis=0)
    res.err = 0.5 * (np.percentile(chain, 84, axis=0) - np.percentile(chain, 16, axis=0))
    res.cov = np.atleast_2d(np.cov(chain, rowvar=False))
    res.chi2 = float(np.sum(_residuals(res.x, res.model, res.data) ** 2))
    try:
        tau = s.get_autocorr_time(discard=steps // 2, quiet=True)
        res.tau_max = float(np.nanmax(tau))
    except (ValueError, FloatingPointError):
        res.tau_max = float("nan")
    return res


def derived(res: FitResult, masses: dict | None = None) -> dict[str, tuple[float, float]]:
    """Derived quantities with errors (from the samples or by linear propagation)."""
    out = {}
    samples = res.samples
    if samples is None:
        rng = np.random.default_rng(0)
        try:
            samples = rng.multivariate_normal(res.x, res.cov, size=4000)
        except (np.linalg.LinAlgError, ValueError):
            samples = np.tile(res.x, (2, 1))
        span = max(d.t.max() for d in res.data) - min(d.t.min() for d in res.data)
        lb = np.array([bounds_for(n, span)[0] for n in res.names])
        ub = np.array([bounds_for(n, span)[1] for n in res.names])
        samples = samples[np.all((samples >= lb) & (samples <= ub), axis=1)]
    vals: dict[str, list] = {}
    for s in samples:
        par = res.model.unpack(s)
        for name, q in par.items():
            if "asini" in q:
                vals.setdefault(f"{name}.K_s", []).append(
                    lite_semi_amplitude_d(q["asini"], q["e3"], q["w3"]) * DAY_S
                )
                vals.setdefault(f"{name}.f_m3", []).append(mass_function(q["asini"], q["P3"]))
                vals.setdefault(f"{name}.P3_yr", []).append(q["P3"] / 365.25)
            if "dw" in q and q["dw"] != 0:
                vals.setdefault(f"{name}.U_yr", []).append(TWO_PI * q["P"] / q["dw"] / 365.25)
                vals.setdefault(f"{name}.Pa_d", []).append(q["P"] / (1 - q["dw"] / TWO_PI))
        if res.model.shared and len(par) >= 2:
            a = [q["asini"] for q in par.values()]
            if a[0] > 0:
                vals.setdefault("mass_ratio_A_over_B", []).append(a[1] / a[0])
    for k, v in vals.items():
        v = np.asarray(v)
        v = v[np.isfinite(v)]
        if len(v):
            out[k] = (
                float(np.median(v)),
                float(0.5 * (np.percentile(v, 84) - np.percentile(v, 16))),
            )
    return out


# ------------------------------------------------------------------ driver
def prepare_minima(cfg) -> Minima:
    """The minima of fit.minima in BJD_TDB, without those removed by fit.exclude."""
    fit = cfg.fit
    base = Path(cfg.base_dir or ".")
    path = Path(fit.minima)
    if not path.is_absolute():
        path = base / path
    mins = read_minima(path)
    if fit.time_system != "BJD_TDB":
        mins.time = to_bjd_tdb(
            mins.time, fit.time_system, fit.coordinates.ra_deg, fit.coordinates.dec_deg
        )
    ex = fit.exclude
    keep = np.ones(len(mins), bool)
    if ex.before is not None:
        keep &= mins.time >= ex.before
    if ex.after is not None:
        keep &= mins.time <= ex.after
    if ex.sources:
        keep &= ~np.isin(mins.source, ex.sources)
    if ex.pairs:
        keep &= ~np.isin(mins.pair, ex.pairs)
    return mins.select(keep)


def prepare(cfg) -> tuple[Model, list[PairData], Minima]:
    fit = cfg.fit
    mins = prepare_minima(cfg)
    pairs = sorted(set(mins.pair.tolist()))
    if not pairs:
        raise ValueError("no minima left after fit.exclude")
    models, data = [], []
    for name in pairs:
        terms = list(fit.stage1.pairs.get(name, fit.stage1.terms))
        orb = cfg.orbits.get(name)
        inc = math.radians(orb.i) if orb is not None and orb.i is not None else math.pi / 2
        pm = PairModel(name, terms, inc)
        sel = mins.pair == name
        t, kind = mins.time[sel], mins.kind[sel]
        if name in fit.ephemerides:
            T0, P = fit.ephemerides[name].T0, fit.ephemerides[name].P
        else:
            if orb is None or orb.P is None:
                raise ValueError(f"fit: give fit.ephemerides.{name} or orbits.{name}.P")
            # orbits.<pair>.P is the osculating period at t = 0; the mean period of the
            # eclipses can differ by 1e-4 and more (a third body), which over thousands of
            # cycles numbers the minima wrongly: the epochs are counted from one primary to the
            # next and the linear ephemeris of the primaries is fitted first
            from .timing import linear_ephemeris

            prim = np.sort(t[kind == 0])
            if len(prim) >= 3:
                T0, P, _ = linear_ephemeris(prim, float(orb.P))
            else:
                P = float(orb.P)
                T0 = float(prim[0] if len(prim) else t[0] - P / 2)
        n = assign_epochs(t, kind, T0, P)
        models.append(pm)
        data.append(
            PairData(pm, n, kind, t, mins.sigma[sel].astype(float), mins.source[sel], T0, P)
        )
    model = Model(models, fit.stage1.shared_outer and sum("lite" in m.terms for m in models) >= 2)
    model.t_ref = float(np.round(np.median(mins.time)))
    for d in data:
        d.t = d.t - model.t_ref
        d.T0 = d.T0 - model.t_ref
    return model, data, mins


def rescale_errors(res: FitResult, how: str) -> dict:
    """Error factors that make the reduced chi^2 one (global or per source)."""
    if how == "none":
        return {}
    r = _residuals(res.x, res.model, res.data)
    src = np.concatenate([d.source for d in res.data])
    scale = {}
    groups = (
        {"all": np.ones(len(r), bool)}
        if how == "global"
        else {s: src == s for s in sorted(set(src.tolist()))}
    )
    npar = len(res.x)
    for g, m in groups.items():
        k = int(m.sum())
        dof = max(k - npar * k / len(r), 1.0)
        f = math.sqrt(float(np.sum(r[m] ** 2)) / dof)
        scale[g] = max(f, 1.0) if how == "per_source" else f
    for d in res.data:
        if how == "global":
            d.sigma = d.sigma * scale["all"]
        else:
            d.sigma = d.sigma * np.array([scale[s] for s in d.source])
    return scale


def clip_stage1(res: FitResult, model, data: list[PairData], st, log=print, rounds: int = 3):
    """Leave out the minima that deviate by more than st.clip_sigma robust standard deviations
    of their source (sources with fewer than 5 minima together) and fit again; repeated until
    nothing more is left out (at most `rounds` times, at most st.max_clip_fraction of all)."""
    from dataclasses import replace

    from .dynfit import clip_outliers, source_groups

    n_all = sum(len(d.t) for d in data)
    outliers: list = []
    for _ in range(rounds):
        z = _residuals(res.x, model, data)
        src = np.concatenate([d.source for d in data])
        groups = source_groups(src, 5)
        frac = max(st.max_clip_fraction - len(outliers) / n_all, 0.0) * n_all / len(z)
        keep, ratio = clip_outliers(z, groups, np.ones(len(z), bool), st.clip_sigma, frac)
        if np.all(keep):
            break
        new, k0 = [], 0
        for d in data:
            m = keep[k0 : k0 + len(d.t)]
            rr = ratio[k0 : k0 + len(d.t)]
            for j in np.where(~m)[0]:
                outliers.append(
                    (d.model, d.n[j], d.kind[j], d.t[j], d.sigma[j], d.source[j], rr[j])
                )
            new.append(
                replace(d, n=d.n[m], kind=d.kind[m], t=d.t[m], sigma=d.sigma[m], source=d.source[m])
            )
            k0 += len(d.t)
        data = new
        log(f"stage 1: {int(np.sum(~keep))} outlying minima left out, fitting again")
        res = fit_stage1(model, data, st.P3_range_d)
    rows = []
    for pm, n, kind, t, sig, src, rat in outliers:  # residuals against the final fit
        pred = model.predict(res.x, pm, np.array([n]), np.array([kind]), np.array([t]))[0]
        rows.append(
            {
                "pair": pm.name,
                "type": "ps"[int(kind)],
                "time": float(t + model.t_ref),
                "sigma": float(sig),
                "source": str(src),
                "residual_min": float((t - pred) * 1440.0),
                "robust_sigma": float(rat),
            }
        )
    res.outliers = rows
    return res, data


def run_fit(cfg, out_dir: str | Path | None = None, plots: bool = True, sampler=None) -> FitResult:
    fit = cfg.fit
    if fit is None:
        raise ValueError("the configuration has no 'fit' section")
    model, data, _ = prepare(cfg)
    st = fit.stage1
    res = fit_stage1(model, data, st.P3_range_d)
    if st.clip_sigma is not None:
        res, data = clip_stage1(res, model, data, st, log=lambda *a: None)
    if fit.error_scale != "none":
        res.scale = rescale_errors(res, fit.error_scale)
        scale = res.scale
        outl = res.outliers
        res = fit_stage1(model, data, st.P3_range_d)
        res.scale = scale
        res.outliers = outl
    if (sampler or st.sampler) == "emcee":
        res = run_emcee(res, st.walkers, st.steps, st.seed)
    out = Path(out_dir or fit.output.dir)
    if not out.is_absolute() and cfg.base_dir and out_dir is None:
        out = Path(cfg.base_dir) / out
    write_outputs(res, cfg, out, plots)
    return res


def _fmt(v, e):
    if not np.isfinite(e) or e <= 0:
        return f"{v:.6g}", ""
    digits = max(0, -math.floor(math.log10(e)) + 1)
    return f"{v:.{digits}f}", f"{e:.{digits}f}"


UNITS = {
    "T0": "d",
    "P": "d",
    "Q": "d",
    "e": "",
    "w0": "rad",
    "dw": "rad/cycle",
    "asini": "AU",
    "P3": "d",
    "T3": "d",
    "e3": "",
    "w3": "rad",
    "K_s": "s",
    "f_m3": "Msun",
    "P3_yr": "yr",
    "U_yr": "yr",
    "Pa_d": "d",
    "mass_ratio_A_over_B": "",
}


def table_rows(res: FitResult):
    rows = []
    for n, v, e in zip(res.names, res.x, res.err, strict=True):
        rows.append((n, res.absolute(n, v), e, UNITS.get(n.split(".", 1)[1], "")))
    for n, (v, e) in derived(res).items():
        rows.append((n, v, e, UNITS.get(n.split(".", 1)[-1], "")))
    return rows


def write_outputs(res: FitResult, cfg, out: Path, plots: bool = True):
    out.mkdir(parents=True, exist_ok=True)
    fit = cfg.fit
    rows = table_rows(res)
    if "csv" in fit.output.tables:
        with open(out / "stage1_parameters.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["parameter", "value", "error", "unit"])
            for r in rows:
                w.writerow([r[0], f"{r[1]:.10g}", f"{r[2]:.3g}", r[3]])
    if "latex" in fit.output.tables:
        with open(out / "stage1_parameters.tex", "w", encoding="utf-8") as f:
            f.write("\\begin{tabular}{lll}\n\\hline\nParameter & Value & Unit \\\\\n\\hline\n")
            for n, v, e, u in rows:
                s, se = _fmt(v, e)
                val = f"${s} \\pm {se}$" if se else f"${s}$"
                f.write(f"{n.replace('_', chr(92) + '_')} & {val} & {u} \\\\\n")
            f.write("\\hline\n\\end{tabular}\n")
    # residual table
    with open(out / "stage1_oc.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(
            ["pair", "type", "epoch", "time", "sigma", "source", "o_c_d", "model_d", "residual_d"]
        )
        for d in res.data:
            lin = d.T0 + d.P * (d.n + d.kind / 2.0)
            model_t = res.model.predict(res.x, d.model, d.n, d.kind, d.t)
            for k in range(len(d.t)):
                w.writerow(
                    [
                        d.model.name,
                        "ps"[d.kind[k]],
                        f"{d.n[k] + d.kind[k] / 2:g}",
                        f"{d.t[k] + res.model.t_ref:.6f}",
                        f"{d.sigma[k]:.6f}",
                        d.source[k],
                        f"{d.t[k] - lin[k]:.6f}",
                        f"{model_t[k] - lin[k]:.6f}",
                        f"{d.t[k] - model_t[k]:.6f}",
                    ]
                )
    summary = {
        "chi2": res.chi2,
        "dof": res.dof,
        "reduced_chi2": res.chi2 / max(res.dof, 1),
        "parameters": {
            n: [float(res.absolute(n, v)), float(e)]
            for n, v, e in zip(res.names, res.x, res.err, strict=True)
        },
        "derived": {k: list(v) for k, v in derived(res).items()},
        "error_scale": res.scale,
        "outliers": res.outliers,
        "sampler": "emcee" if res.samples is not None else "none",
        "values": (
            "median of the emcee samples, errors half the 16-84 % range, chi2 at the median"
            if res.samples is not None
            else "least-squares solution, errors of the covariance matrix"
        ),
        "ephemerides": {d.model.name: {"T0": d.T0 + res.model.t_ref, "P": d.P} for d in res.data},
        "t_ref": res.model.t_ref,
    }
    if fit.stage2 is not None:
        summary["stage2"] = "configured: swingersim fit --stage 2 fits the joint N-body model"
    if fit.output.predict:
        summary["predictions_file"] = "stage1_predictions.csv"
        pr = fit.output.predict
        write_predictions(res, out / "stage1_predictions.csv", pr.years, pr.pairs or None, pr.start)
    if res.best_fit is not None:
        summary["best_fit"] = {
            "chi2": res.best_fit["chi2"],
            "parameters": {
                n: [float(res.absolute(n, v)), float(e)]
                for n, v, e in zip(res.names, res.best_fit["x"], res.best_fit["err"], strict=True)
            },
        }
        summary["tau_max_steps"] = res.tau_max
    from .io import software_versions

    summary["software"] = software_versions()
    (out / "stage1_summary.json").write_text(json.dumps(summary, indent=2))
    if plots:
        plot_oc(res, out, components="oc_components" in fit.output.plots)
        if "corner" in fit.output.plots and res.samples is not None:
            try:
                import corner

                fig = corner.corner(res.samples, labels=res.names)
                fig.savefig(out / "stage1_corner.png", dpi=90)
            except ImportError:
                pass


def write_predictions(res: FitResult, path, years=5, pairs=None, start=None):
    """Future minima with the 1-sigma spread of the model (half the 16-84 % range over 300
    parameter sets from the samples or the covariance), vectorised over the epochs."""
    rng = np.random.default_rng(2)
    if res.samples is not None:
        samp = res.samples[rng.integers(0, len(res.samples), 300)]
    else:
        try:
            samp = rng.multivariate_normal(res.x, res.cov, size=300)
        except (np.linalg.LinAlgError, ValueError):
            samp = np.tile(res.x, (2, 1))
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pair", "type", "epoch", "time_bjd_tdb", "sigma_d"])
        for d in res.data:
            if pairs and d.model.name not in pairs:
                continue
            t_start = (start - res.model.t_ref) if start is not None else float(d.t.max())
            n0 = math.ceil((t_start - d.T0) / d.P)
            n1 = n0 + int(years * 365.25 / d.P)
            step = max(1, (n1 - n0) // 2000)
            n = np.arange(n0, n1 + 1, step, dtype=float)
            nn = np.repeat(n, 2)  # epochs in the order n p, n s, n+1 p, ...
            kk = np.tile([0, 1], len(n))
            t = res.model.predict(res.x, d.model, nn, kk)
            t = res.model.predict(res.x, d.model, nn, kk, t)
            ts = np.array([res.model.predict(x, d.model, nn, kk, t) for x in samp])
            p16, p84 = np.percentile(ts, [16, 84], axis=0)
            sig = 0.5 * (p84 - p16)
            for j in range(len(nn)):
                w.writerow(
                    [
                        d.model.name,
                        "ps"[kk[j]],
                        nn[j] + kk[j] / 2,
                        f"{t[j] + res.model.t_ref:.6f}",
                        f"{sig[j]:.6f}",
                    ]
                )


def plot_oc(res: FitResult, out: Path, components=True):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    k = len(res.data)
    x0 = 2450000.0 - res.model.t_ref  # plot BJD - x0
    fig, axes = plt.subplots(
        2 + int(components),
        k,
        figsize=(5.5 * k, 7 + 3 * int(components)),
        squeeze=False,
        sharex="col",
    )
    for c, d in enumerate(res.data):
        lin = d.T0 + d.P * (d.n + d.kind / 2.0)
        tt = np.linspace(d.t.min(), d.t.max(), 1500)
        nn = np.round((tt - d.T0) / d.P)
        ax = axes[0][c]
        for j, mk in ((0, "o"), (1, "s")):
            m = d.kind == j
            ax.errorbar(
                d.t[m] - x0,
                (d.t[m] - lin[m]) * 1440,
                d.sigma[m] * 1440,
                fmt=mk,
                ms=3,
                lw=0.6,
                label="primary" if j == 0 else "secondary",
            )
            tm = res.model.predict(res.x, d.model, nn, np.full(len(nn), j), tt)
            ax.plot(tt - x0, (tm - (d.T0 + d.P * (nn + j / 2))) * 1440, lw=1)
        ax.set_ylabel("O−C [min]")
        ax.set_title(f"pair {d.model.name}")
        ax.legend(fontsize=8)
        r = d.t - res.model.predict(res.x, d.model, d.n, d.kind, d.t)
        axes[1][c].errorbar(d.t - x0, r * 1440, d.sigma * 1440, fmt="o", ms=3, lw=0.6)
        axes[1][c].axhline(0, color="0.5", lw=0.8)
        axes[1][c].set_ylabel("residual [min]")
        if components:
            comps = res.model.components(res.x, d.model, nn, np.zeros(len(nn), int), tt)
            for name, v in comps.items():
                axes[2][c].plot(tt - x0, v * 1440, label=name)
            if "apsidal" in comps:
                v2 = res.model.components(res.x, d.model, nn, np.ones(len(nn), int), tt)["apsidal"]
                axes[2][c].plot(tt - x0, v2 * 1440, label="apsidal (secondary)", ls="--")
            axes[2][c].legend(fontsize=8)
            axes[2][c].set_ylabel("component [min]")
        axes[-1][c].set_xlabel("BJD − 2450000")
    fig.tight_layout()
    fig.savefig(out / "stage1_oc.png", dpi=110)
    plt.close(fig)


def report(res: FitResult) -> str:
    lines = [
        (
            f"stage 1: chi^2 = {res.chi2:.1f} for {res.dof} degrees of freedom "
            f"(reduced {res.chi2 / max(res.dof, 1):.2f})"
        )
    ]
    for n, v, e, u in table_rows(res):
        s, se = _fmt(v, e)
        lines.append(f"  {n:24s} {s:>18s} ± {se:<12s} {u}")
    if res.outliers:
        lines.append(f"  {len(res.outliers)} outlying minima left out:")
        for o in res.outliers:
            lines.append(
                f"    {o['pair']} {o['type']} {o['time']:.5f} {o['residual_min']:+9.2f} min "
                f"({o['robust_sigma']:.1f} robust sigma) {o['source']}"
            )
    return "\n".join(lines)


# ------------------------------------------------------------------ synthetic data
def synthetic_minima(
    path, pairs: dict, outer: dict | None, t0, t1, sigma=1e-4, seed=1, every=1, shared_out=True
):
    """Write a minima table from known parameters (tests and validation)."""
    rng = np.random.default_rng(seed)
    rows = []
    for k, (name, p) in enumerate(pairs.items()):
        pm = PairModel(name, p["terms"], math.radians(p.get("i", 90.0)))
        names = (
            ["T0", "P"]
            + (["Q"] if "quadratic" in pm.terms else [])
            + (["e", "w0", "dw"] if "apsidal" in pm.terms else [])
            + (["asini", *LITE_KEYS] if "lite" in pm.terms else [])
        )
        q = dict(p)
        if "lite" in pm.terms and outer is not None:
            q.update(outer)
            if k == 1 and shared_out:
                q["w3"] = outer["w3"] + math.pi
        m = Model([pm], False)
        x = np.array([q[nn] for nn in names])
        n0 = math.ceil((t0 - p["T0"]) / p["P"])
        n1 = math.floor((t1 - p["T0"]) / p["P"])
        for n in range(n0, n1 + 1, every):
            for j in (0, 1):
                nn, kk = np.array([n], float), np.array([j])
                t = float(m.predict(x, pm, nn, kk)[0])
                t = float(m.predict(x, pm, nn, kk, np.array([t]))[0])
                rows.append((name, "ps"[j], t + sigma * rng.normal(), sigma))
    with open(path, "w", encoding="utf-8") as f:
        f.write("# time_system: BJD_TDB\npair,type,time,sigma,source,method\n")
        f.writelines(f"{r[0]},{r[1]},{r[2]:.7f},{r[3]:.6f},synthetic,individual\n" for r in rows)
    return len(rows)


__all__ = [
    "apsidal_shift",
    "fit_stage1",
    "lite",
    "lite_semi_amplitude_d",
    "mass_function",
    "read_minima",
    "run_fit",
    "synthetic_minima",
    "to_bjd_tdb",
]
