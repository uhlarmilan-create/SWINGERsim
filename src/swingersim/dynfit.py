"""Fitting of observed times of minima, stage 2 (swingersim fit, ``fit.stage2``): the joint
dynamical model of all eclipsing pairs.

Stage 1 (ocfit.py) describes the O - C of every pair with analytic terms (light-time effect,
apsidal motion). Stage 2 computes the times of minima of all pairs from one N-body model of the
whole system - the same model as ``swingersim run`` and ``swingersim timings`` (Newton, 1PN, tides as
configured) - and fits the osculating elements and masses at t = 0 to the measured minima.
Light-time effect, apsidal motion, the dynamical perturbations of a compact system and the
precession of the orbits therefore come out of the same fit, and the fitted configuration can
be run directly (stability, evolution forward and backward, light curves, predictions).

Model time of a measured minimum (pair, type, time t): the N-body system is integrated to
shortly before the expected conjunction; the two-body motion of the pair from that state gives
the nearest conjunction of the right type (the minimum of the projected separation with the
first star farther for a primary minimum, as in ``swingersim timings``); the N-body system is then
integrated to that conjunction and the time refined once more. The observed time adds the
light-time of the pair's centre of mass, t_obs = t + z_cm / c. Minima before t = 0 are computed
by a second integration backward in time. One evaluation costs one N-body integration over
the span of the data in each direction.

The fit is a weighted least-squares fit (scipy.optimize.least_squares, trust-region reflective
with bounds) of the free parameters listed in ``fit.stage2.free``; the Jacobian is computed by
central differences (one-sided at a bound), two N-body integrations per free parameter, in
parallel on several cores.
Errors come from the covariance matrix, scaled by sqrt(chi^2_red) when it exceeds 1.

What O - C alone can and cannot give: periods, phases, eccentricities and arguments of
pericentre of all orbits, and - through the light-time amplitude and the dynamical
perturbations - the mass ratios of the hierarchy. The inclinations of the eclipsing pairs and
the absolute masses need light curves and radial velocities; keep them fixed at the values of
a light-curve / spectroscopic solution. Mutual inclinations and nodes are constrained only in
compact systems with a measurable dynamical (not light-time) signal.
"""

from __future__ import annotations

import copy
import csv
import json
import math
import time
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import yaml
from scipy.optimize import brentq, least_squares

from . import units
from .config import Config
from .kepler import rotation_matrix, state_to_elements
from .ocfit import kepler_E
from .timing import _binaries, linear_ephemeris

DPY = units.DAYS_PER_YEAR
J2000_BJD = 2451544.5  # year 2000.0 of the plots (JD of 2000 January 1, 0 h)

# default finite-difference steps: absolute, or relative to the value (rel_)
E_MAX = 0.95  # largest eccentricity of a fitted orbit
REINCLUDE = 0.8  # outliers come back below this fraction of the clipping threshold

DEFAULT_STEPS = {
    "P": ("rel", 1e-6),
    "M0": ("abs", 0.01),
    "e": ("abs", 1e-4),
    "ecosw": ("abs", 1e-4),
    "esinw": ("abs", 1e-4),
    "omega": ("abs", 0.01),
    "Omega": ("abs", 0.01),
    "i": ("abs", 0.01),
    "mass": ("rel", 1e-4),
    "rv": ("abs", 0.01),
    "R": ("rel", 1e-4),
    "Teff": ("abs", 1.0),
    "apsidal_k2": ("rel", 1e-3),
    "rg2": ("rel", 1e-3),
    "Pdot": ("abs", 1e-8),  # [d/yr]: the O - C grows as (1/2) (Pdot / P) t^2
}
STAR_KEYS = {"R": "radius_rsun", "Teff": "teff_k", "apsidal_k2": "apsidal_k2", "rg2": "rg2"}
UNITS = {
    "P": "d",
    "M0": "deg",
    "omega": "deg",
    "Omega": "deg",
    "i": "deg",
    "mass": "Msun",
    "rv": "km/s",
    "R": "Rsun",
    "Teff": "K",
    "Pdot": "d/yr",
}


# ------------------------------------------------------------------ data
@dataclass
class Observed:
    pair: np.ndarray  # str
    kind: np.ndarray  # 0 primary, 1 secondary
    time: np.ndarray  # BJD_TDB [d]
    sigma: np.ndarray  # [d]
    source: np.ndarray  # str

    def __len__(self) -> int:
        return len(self.time)


def observed_from_minima(mins, pairs: list[str] | None, max_sigma_d: float | None) -> Observed:
    """The minima of ocfit.read_minima, restricted to the selected pairs and errors."""
    keep = np.ones(len(mins), bool)
    if pairs:
        keep &= np.isin(mins.pair, pairs)
    if max_sigma_d is not None:
        keep &= mins.sigma <= max_sigma_d
    keep &= np.isfinite(mins.time) & np.isfinite(mins.sigma) & (mins.sigma > 0)
    return Observed(
        mins.pair[keep].astype(str),
        mins.kind[keep].astype(int),
        mins.time[keep].astype(float),
        mins.sigma[keep].astype(float),
        mins.source[keep].astype(str),
    )


# ------------------------------------------------------------------ two-body helpers
def _relative_orbit(mu: float, r: np.ndarray, v: np.ndarray):
    """Function tau -> (x, y, z, vx, vy) of the two-body relative orbit from the state (r, v),
    vectorised in tau [yr]."""
    el = state_to_elements(mu, r, v)
    if not 0.0 <= el.e < 1.0:
        return None
    a, e = el.a, el.e
    n = math.sqrt(mu / a**3)
    R = rotation_matrix(el.inc, el.Omega, el.omega)
    b = math.sqrt(1.0 - e * e)

    def at(tau):
        tau = np.atleast_1d(np.asarray(tau, float))
        E = kepler_E(np.mod(el.M + n * tau, 2.0 * math.pi), e)
        cE, sE = np.cos(E), np.sin(E)
        rpf = np.vstack([a * (cE - e), a * b * sE, np.zeros_like(E)])
        k = n * a / (1.0 - e * cE)
        vpf = np.vstack([-k * sE, k * b * cE, np.zeros_like(E)])
        x, y, z = R @ rpf
        vx, vy, _ = R @ vpf
        return x, y, z, vx, vy

    return at, 2.0 * math.pi / n


def _conjunction(mu, r, v, lo: float, hi: float, target: float, kind: int):
    """Time (relative to the state) of the conjunction of the given type nearest to `target`
    within [lo, hi]: a minimum of the projected separation (x vx + y vy from - to +) with the
    first star farther (z of the relative vector < 0) for a primary minimum."""
    orb = _relative_orbit(mu, r, v)
    if orb is None:
        return None
    at, period = orb
    n = max(math.ceil((hi - lo) / period * 48.0) + 1, 8)
    tau = np.linspace(lo, hi, n)
    x, y, _, vx, vy = at(tau)
    g = x * vx + y * vy
    best = None
    for k in np.nonzero((g[:-1] < 0.0) & (g[1:] >= 0.0))[0]:

        def gf(s):
            xx, yy, _, vxx, vyy = at(s)
            return float(xx[0] * vxx[0] + yy[0] * vyy[0])

        try:
            s = brentq(gf, tau[k], tau[k + 1], xtol=1e-13)
        except ValueError:
            continue
        zz = at(s)[2][0]
        if (kind == 0) != (zz < 0.0):
            continue
        if best is None or abs(s - target) < abs(best - target):
            best = s
    return best


# ------------------------------------------------------------------ model times
def build_tree_root(cfg):
    from .hierarchy import build_tree

    return build_tree(cfg)[0]


def _leaves(c) -> list[int]:
    return c.stars() if hasattr(c, "stars") else [c]


def model_times(
    cfg: Config, obs: Observed, bjd_t0: float, separation: np.ndarray | None = None
) -> np.ndarray:
    """Model BJD of every measured minimum (NaN where no conjunction of that type was found).

    ``obs.pair`` names an orbit of the configuration: a binary of two stars (eclipses), or a
    wider orbit (e.g. the 51-d orbit of a third star around an inner binary), whose
    conjunctions - the centres of mass of its two components aligned along the line of sight,
    type p when the first component (left in the hierarchy) is behind - are the centres of
    third-body eclipse events. With ``separation`` (an array of the length of obs) the
    projected separation [AU] of the two components at each conjunction is stored there.
    The integration is obsmodel.observables (one per direction; the module docstring of
    dynfit describes the search of a conjunction)."""
    from .obsmodel import observables

    return observables(cfg, bjd_t0, minima=obs, separation=separation)["t"]


# ------------------------------------------------------------------ parameters
def _split(name: str) -> tuple[str, str, str]:
    parts = name.split(".")
    if parts[0] == "masses" and len(parts) == 2:
        return "masses", parts[1], "mass"
    if name == "rv.gamma":
        return "rv", "gamma", "rv"
    if parts[0] == "stars" and len(parts) == 3:
        if parts[2] not in STAR_KEYS:
            raise ValueError(
                f"fit.stage2.free: unknown star parameter '{name}' (R, Teff, apsidal_k2, rg2)"
            )
        return "stars", parts[1], parts[2]
    if parts[0] == "rv" and len(parts) >= 3 and parts[1] == "offset":
        return "rv", "offset." + ".".join(parts[2:]), "rv"
    if parts[0] == "orbits" and len(parts) == 3:
        if parts[2] not in ("P", "M0", "e", "omega", "Omega", "i", "ecosw", "esinw", "Pdot"):
            raise ValueError(f"fit.stage2.free: unknown orbit parameter '{name}'")
        return "orbits", parts[1], parts[2]
    raise ValueError(
        f"fit.stage2.free: '{name}' - use orbits.<orbit>.<P|M0|e|omega|Omega|i|ecosw|esinw|Pdot>, "
        "masses.<star>, stars.<star>.<R|Teff|apsidal_k2|rg2>, rv.gamma or rv.offset.<source>"
    )


def start_values(raw: dict, names: list[str]) -> np.ndarray:
    x = []
    for name in names:
        grp, key, p = _split(name)
        if grp == "masses":
            x.append(float(raw["masses"][key]))
            continue
        if grp == "stars":
            v = ((raw.get("stars") or {}).get(key) or {}).get(p)
            if v is None:  # the value of the stellar model
                from .stellar import star_properties

                cfg = Config.model_validate(raw)
                v = getattr(star_properties(cfg, key, float(raw["masses"][key])), STAR_KEYS[p])
                if v is None:
                    raise ValueError(f"fit.stage2: no {p} of star {key} in the stellar model")
            x.append(float(v))
            continue
        if grp == "rv":
            spec = (raw.get("fit") or {}).get("rv") or {}
            if key == "gamma":
                x.append(float(spec.get("gamma", 0.0)))
            else:
                x.append(float((spec.get("offsets") or {}).get(key[7:], 0.0)))
            continue
        o = raw["orbits"][key]
        if p in ("ecosw", "esinw"):
            e, w = float(o.get("e", 0.0)), math.radians(float(o.get("omega", 0.0)))
            x.append(e * (math.cos(w) if p == "ecosw" else math.sin(w)))
        else:
            if p == "P" and o.get("P") is None:
                raise ValueError(f"fit.stage2: orbit {key} is given by a; give its period P")
            x.append(float(o.get(p, 0.0)))
    return np.array(x)


def apply(raw: dict, names: list[str], x) -> dict:
    """Copy of the raw configuration with the free parameters set to x."""
    d = copy.deepcopy(raw)
    ew: dict[str, dict[str, float]] = {}
    for name, val in zip(names, x):
        grp, key, p = _split(name)
        if grp == "masses":
            d["masses"][key] = float(val)
        elif grp == "stars":
            stars = d.get("stars") or {}
            d["stars"] = stars
            st = stars.get(key) or {}
            stars[key] = st
            st[p] = float(val)
        elif grp == "rv":
            spec = d.setdefault("fit", {}).setdefault("rv", {})
            if key == "gamma":
                spec["gamma"] = float(val)
            else:
                spec.setdefault("offsets", {})[key[7:]] = float(val)
        elif p in ("ecosw", "esinw"):
            ew.setdefault(key, {})[p] = float(val)
        else:
            d["orbits"][key][p] = float(val)
    for key, q in ew.items():
        o = d["orbits"][key]
        e0, w0 = float(o.get("e", 0.0)), math.radians(float(o.get("omega", 0.0)))
        ec = q.get("ecosw", e0 * math.cos(w0))
        es = q.get("esinw", e0 * math.sin(w0))
        # the box of ecosw, esinw (each +-0.95) holds e up to 1.34: the corners are cut at the
        # bound of e (0.95), so that the model exists everywhere in the box (no cliff of
        # failed models next to the solution for the finite differences)
        o["e"] = float(min(math.hypot(ec, es), E_MAX))
        o["omega"] = float(math.degrees(math.atan2(es, ec)) % 360.0)
    return d


ANGLES = ("M0", "omega", "Omega")


def prior_residuals(raw: dict, names: list[str], x, priors: dict) -> np.ndarray:
    """(value - mean) / sigma of every prior at the free parameters x; angles M0, omega and
    Omega wrapped to +-180 deg."""
    if not priors:
        return np.zeros(0)
    d = apply(raw, names, x)
    keys = list(priors)
    v = start_values(d, keys)
    out = np.empty(len(keys))
    for k, key in enumerate(keys):
        mu, sig = priors[key]
        dv = v[k] - mu
        if _split(key)[2] in ANGLES:
            dv = (dv + 180.0) % 360.0 - 180.0
        out[k] = dv / sig
    return out


def default_steps(names, x0, steps: dict) -> np.ndarray:
    out = []
    for name, v in zip(names, x0):
        if name in steps:
            out.append(float(steps[name]))
            continue
        kind, s = DEFAULT_STEPS[_split(name)[2]]
        out.append(s * abs(v) if kind == "rel" and v != 0 else s)
    return np.array(out)


def default_bounds(names, x0, bounds: dict) -> tuple[np.ndarray, np.ndarray]:
    lo, hi = [], []
    for name, v in zip(names, x0):
        if name in bounds:
            a, b = bounds[name]
        else:
            p = _split(name)[2]
            a, b = {
                "P": (0.5 * v, 2.0 * v),
                "mass": (0.05 * v, 20.0 * v),
                "e": (0.0, E_MAX),
                "ecosw": (-E_MAX, E_MAX),
                "esinw": (-E_MAX, E_MAX),
                "i": (0.0, 180.0),
                "R": (0.2 * v, 5.0 * v),
                "Teff": (0.5 * v, 2.0 * v),
                "apsidal_k2": (1e-5, 0.74),
                "rg2": (0.005, 0.4),
            }.get(p, (-np.inf, np.inf))
        lo.append(a)
        hi.append(b)
    return np.array(lo, float), np.array(hi, float)


# ------------------------------------------------------------------ fit
@dataclass
class Stage2Result:
    names: list[str]
    x0: np.ndarray
    x: np.ndarray
    err: np.ndarray
    cov: np.ndarray
    obs: Observed
    t_start: np.ndarray  # model BJD with the start values
    t_model: np.ndarray  # model BJD with the fitted values
    chi2_start: float
    chi2: float
    chi2_red: float
    n_eval: int
    wall_s: float
    raw: dict
    bjd_t0: float
    converged: str = ""
    history: list = field(default_factory=list)
    err_linear: np.ndarray | None = None  # errors of the covariance matrix alone
    samples: np.ndarray | None = None  # emcee samples (second half of the chain)
    priors: dict = field(default_factory=dict)  # name -> (mean, sigma)
    chi2_priors: float = 0.0  # contribution of the priors to chi2
    sigma_in: np.ndarray | None = None  # quoted errors of the fitted minima [d]
    error_scale: str = "none"
    error_model: str = "scale"
    error_factors: dict = field(default_factory=dict)  # group -> factor or jitter [d]
    outliers: Observed | None = None  # minima left out by clip_sigma
    t_outliers: np.ndarray | None = None  # model BJD of the left-out minima
    z_outliers: np.ndarray | None = None  # their deviation [robust sigma of the source]
    passes: list = field(default_factory=list)  # summary of the fits in a row
    note: str = ""  # a warning of the fits in a row
    prediction: object = None  # Prediction of fit.output.predict
    rv: object = None  # radvel.ObservedRV fitted (errors as used)
    rv_start: np.ndarray | None = None  # model velocities with the start values [km/s]
    rv_model: np.ndarray | None = None  # with the fitted values
    rv_sigma_in: np.ndarray | None = None  # quoted errors of the fitted velocities
    rv_outliers: object = None  # velocities left out by clip_sigma
    rv_model_outliers: np.ndarray | None = None
    rv_outliers_z: np.ndarray | None = None
    inc: list | None = None  # incfit.ObservedInclination fitted
    inc_start: list | None = None  # model inclinations with the start values [deg]
    inc_model: list | None = None  # with the fitted values
    lc: list = field(default_factory=list)  # lcfit.ObservedLC fitted (errors as used)
    lc_start: list = field(default_factory=list)  # calibrated model, start values
    lc_model: list = field(default_factory=list)  # calibrated model, fitted values
    lc_coef: list = field(default_factory=list)  # per light curve {segment: coefficients}
    lc_sigma_in: list = field(default_factory=list)
    fidelity: list = field(default_factory=list)  # passes of the mapped mesh fit
    shape_notes: list = field(default_factory=list)  # star shapes of the light-curve models
    lc_beta: dict = field(default_factory=dict)  # red-noise factor per light-curve segment
    base_dir: str | None = None  # folder of the configuration (relative data paths)
    dof: int = 0  # degrees of freedom of chi2_red (kept data rows - free parameters)
    loss: str = "linear"
    lc_outliers: list = field(default_factory=list)  # per light curve (data, model, z) or None
    samples_logp: np.ndarray | None = None  # log-probability of the samples
    edges: dict = field(default_factory=dict)  # name -> side where the model fails ("+", "-")
    notes: list = field(default_factory=list)  # remarks on the setup of the fit
    mcmc: dict | None = None  # emcee: walkers, steps, discard, seed, acceptance, tau
    fidelity_note: str = ""  # mapped meshes: what chi^2 and the errors refer to
    dense: dict | None = None  # cache of _dense_model (plots)

    def rms_lc(self, which="model") -> dict:
        """{(source, segment): rms [ppt]} of measured - calibrated model."""
        out = {}
        models = self.lc_model if which == "model" else self.lc_start
        for lc, m in zip(self.lc, models, strict=False):
            for s in np.unique(lc.segment):
                k = lc.segment == s
                out[(lc.source, str(s))] = float(
                    np.sqrt(np.nanmean((lc.flux[k] - m[k]) ** 2)) * 1e3
                )
        return out

    def rms_rv(self, model) -> dict[str, float]:
        if self.rv is None or model is None:
            return {}
        r = self.rv.rv - model
        return {
            s: float(np.sqrt(np.nanmean(r[self.rv.star == s] ** 2)))
            for s in sorted(set(self.rv.star.tolist()))
        }

    def rms_min(self, times) -> dict[str, float]:
        res = (self.obs.time - times) * 1440.0
        return {
            p: float(np.sqrt(np.nanmean(res[self.obs.pair == p] ** 2)))
            for p in sorted(set(self.obs.pair.tolist()))
        }


# ------------------------------------------------------------------ model of the data
_CTX: dict = {}


def _init(ctx: dict) -> None:
    _CTX.clear()
    _CTX.update(ctx)
    from . import cnoise

    cnoise.install()


def _cfg_of(d: dict, base_dir) -> Config:
    cfg = Config.model_validate(d)
    cfg._base_dir = base_dir
    return cfg


def _times_for(x) -> np.ndarray:
    d = apply(_CTX["raw"], _CTX["names"], x)
    return model_times(_cfg_of(d, _CTX["base_dir"]), _CTX["obs"], _CTX["bjd_t0"])


def rv_shift(d: dict, rv) -> np.ndarray:
    """gamma + the offset of the source of every measured velocity [km/s] (from fit.rv of the
    configuration dict d)."""
    spec = (d.get("fit") or {}).get("rv") or {}
    off = spec.get("offsets") or {}
    return float(spec.get("gamma", 0.0)) + np.array([float(off.get(s, 0.0)) for s in rv.source])


def rv_light_offsets(d: dict, base_dir, rv, bjd_t0: float) -> np.ndarray | float:
    """fit.rv.centre light: centre-of-light minus centre-of-mass velocity of every measured
    velocity [km/s] (swingersim.rvlight); centre mass: 0."""
    spec = (d.get("fit") or {}).get("rv") or {}
    if spec.get("centre", "mass") != "light":
        return 0.0
    from .rvlight import light_offset

    t_day = np.asarray(rv.time, float) - bjd_t0
    return light_offset(_cfg_of(d, base_dir), rv.star, t_day, spec.get("band", "V"))


def _mapped(lc) -> bool:
    """Whether the fit model of lc is the coarse mesh plus the stored difference to the
    accurate one (fidelity mapped of a light curve with meshes)."""
    return lc.fidelity == "mapped" and lc.delta is not None


def lc_models_for(d: dict, base_dir, lcs, bjd_t0: float, mode: str = "fit") -> list:
    """Synthetic (uncalibrated) flux of the configuration dict d at the times of every light
    curve. mode fit: what the fit uses (with fidelity mapped the coarse mesh plus the stored
    difference to the accurate one); raw: the same without the stored difference (the part
    that depends on the parameters); accurate / coarse: that model alone."""
    from .lcfit import lc_config, model_flux

    out = []
    for lc in lcs:
        d_lc, _ = lc_config(d, lc)
        cfg = _cfg_of(d_lc, base_dir)
        if mode == "coarse" or (mode in ("fit", "raw") and _mapped(lc)):
            m = model_flux(cfg, lc, bjd_t0, lc.coarse_res, table_phases=64)
            out.append(m + lc.delta if mode == "fit" else m)
        else:
            out.append(model_flux(cfg, lc, bjd_t0, lc.mesh_res))
    return out


# extra data sets set by run_stage2 for the whole fit (fit.inclination), passed to the workers
_EXTRA: dict = {}


def model_at(d: dict, ctx: dict) -> dict:
    """Model of all data of ctx (the full sets, outliers included) at the configuration dict d:
    {"t": BJD of the minima, "rv": velocities with gamma and offsets [km/s], "inc": [model of
    every inclination table], "lc": [raw flux of every light curve (lc_models_for mode raw)]};
    None for a data type whose model failed. Minima, velocities and inclinations come from one
    N-body integration per direction (obsmodel.observables); if it fails, every type is tried
    alone so that the others still count (as with separate integrations)."""
    from .incfit import needs_limit
    from .obsmodel import observables

    obs, rv = ctx["obs"], ctx.get("rv")
    incs, lcs = ctx.get("inc") or [], ctx.get("lc") or []
    base_dir, t0 = ctx["base_dir"], ctx["bjd_t0"]
    has_rv = rv is not None and len(rv) > 0
    req = {
        "minima": obs if len(obs) else None,
        "rv": (rv.star, rv.time) if has_rv else None,
        "incs": [(q.pair, q.time, needs_limit(q)) for q in incs],
    }
    out = {"t": None, "rv": None, "inc": None, "lc": None}

    def store(o):
        if req["minima"] is not None and o.get("t") is not None:
            out["t"] = o["t"]
        if has_rv and o.get("rv") is not None:
            out["rv"] = o["rv"] + rv_shift(d, rv) + rv_light_offsets(d, base_dir, rv, t0)
        if incs and o.get("inc"):
            out["inc"] = o["inc"]

    try:
        store(observables(_cfg_of(d, base_dir), t0, **req))
    except (ValueError, RuntimeError):
        for key in ("minima", "rv", "incs"):
            if not req[key]:
                continue
            try:
                store(observables(_cfg_of(d, base_dir), t0, **{key: req[key]}))
            except (ValueError, RuntimeError):
                pass
    if req["minima"] is None:
        out["t"] = np.zeros(0)
    if lcs:
        try:
            out["lc"] = lc_models_for(d, base_dir, lcs, t0, mode="raw")
        except (ValueError, RuntimeError):
            pass
    return out


def _model_worker(x) -> dict:
    return model_at(apply(_CTX["raw"], _CTX["names"], x), _CTX)


def _model_failed(m: dict, ctx: dict) -> bool:
    """True when a data type of ctx has no model in m (the integration failed)."""
    if m is None:
        return True
    obs, rv = ctx["obs"], ctx.get("rv")
    return bool(
        (len(obs) and m.get("t") is None)
        or (rv is not None and len(rv) and m.get("rv") is None)
        or (ctx.get("inc") and m.get("inc") is None)
        or (ctx.get("lc") and m.get("lc") is None)
    )


def _fail(n: int) -> np.ndarray:
    return np.full(n, 1e3)


def residuals_of(m: dict, x, ctx: dict, w: dict) -> np.ndarray:
    """(observed - model) / sigma of the kept minima, radial velocities, light-curve points
    (calibrated per segment), the inclination rows, then the prior terms. w holds (keep, sigma)
    of every data set; a point without a model counts 1e3 sigma."""
    from .incfit import residuals as inc_residuals
    from .lcfit import residuals as lc_residuals

    obs = ctx["obs"]
    k, s = w["obs"]
    r = (obs.time[k] - m["t"][k]) / s[k] if m["t"] is not None else _fail(int(k.sum()))
    parts = [r]
    rv = ctx.get("rv")
    if rv is not None and len(rv):
        k, s = w["rv"]
        parts.append((rv.rv[k] - m["rv"][k]) / s[k] if m["rv"] is not None else _fail(k.sum()))
    for i, q in enumerate(ctx.get("lc") or []):
        k, s = w["lc"][i]
        if m["lc"] is None:
            parts.append(_fail(int(k.sum())))
            continue
        mm = m["lc"][i] + q.delta if _mapped(q) else m["lc"][i]
        try:
            parts.append(lc_residuals(q.subset(k, s), mm[k]))
        except (ValueError, RuntimeError):
            parts.append(_fail(int(k.sum())))
    for i, q in enumerate(ctx.get("inc") or []):
        try:
            parts.append(inc_residuals(q, m["inc"][i]) if m["inc"] is not None else _fail(len(q)))
        except (ValueError, RuntimeError):
            parts.append(_fail(len(q)))
    out = np.concatenate([np.asarray(p, float) for p in parts])
    out[~np.isfinite(out)] = 1e3
    return np.concatenate([out, prior_residuals(ctx["raw"], ctx["names"], x, ctx.get("priors"))])


class _Problem:
    """Data (full sets, outliers included), current weights (kept points and errors) and the
    model cache of one stage-2 fit. The cache holds the raw model of the last evaluations keyed
    by the parameter vector: start values, solution and the Jacobian columns are computed once
    and reused by the next fit in a row (new errors / outliers) and by the outputs."""

    def __init__(self, raw, names, base_dir, bjd_t0, priors, obs, rv, lcs, incs, cache_size=64):
        self.names = list(names)
        self.ctx = {
            "raw": raw,
            "names": self.names,
            "obs": obs,
            "rv": rv if rv is not None and len(rv) else None,
            "lc": list(lcs or []),
            "inc": list(incs or []),
            "bjd_t0": bjd_t0,
            "base_dir": base_dir,
            "priors": dict(priors or {}),
        }
        self.w = {
            "obs": (np.ones(len(obs), bool), obs.sigma.copy()),
            "rv": (np.ones(len(rv), bool), rv.sigma.copy()) if self.ctx["rv"] else None,
            "lc": [(np.ones(len(q), bool), q.sigma.copy()) for q in self.ctx["lc"]],
        }
        self.cache: OrderedDict = OrderedDict()
        self.cache_size = max(int(cache_size), 8)
        self.n_model = 0

    @property
    def obs(self):
        return self.ctx["obs"]

    @property
    def rv(self):
        return self.ctx["rv"]

    @property
    def lcs(self):
        return self.ctx["lc"]

    @property
    def incs(self):
        return self.ctx["inc"]

    def _store(self, key, m):
        self.cache[key] = m
        self.cache.move_to_end(key)
        while len(self.cache) > self.cache_size:
            self.cache.popitem(last=False)

    def model(self, x) -> dict:
        x = np.asarray(x, float)
        key = x.tobytes()
        if key in self.cache:
            self.cache.move_to_end(key)
            return self.cache[key]
        m = model_at(apply(self.ctx["raw"], self.names, x), self.ctx)
        self.n_model += 1
        self._store(key, m)
        return m

    def models(self, xs, pool, store: bool = True) -> list:
        """Models at all xs; those not in the cache in parallel."""
        xs = [np.asarray(x, float) for x in xs]
        keys = [x.tobytes() for x in xs]
        todo = [i for i, k in enumerate(keys) if k not in self.cache]
        todo = list({keys[i]: i for i in todo}.values())  # each point once
        if pool is not None and len(todo) > 1:
            new = list(pool.map(_model_worker, [xs[i] for i in todo]))
        else:
            new = [model_at(apply(self.ctx["raw"], self.names, xs[i]), self.ctx) for i in todo]
        self.n_model += len(todo)
        got = {keys[i]: m for i, m in zip(todo, new, strict=True)}
        out = []
        for k in keys:
            if k in got:
                out.append(got[k])
            else:
                self.cache.move_to_end(k)
                out.append(self.cache[k])
        if store:
            for k, m in got.items():
                self._store(k, m)
        return out

    def resid(self, m: dict, x) -> np.ndarray:
        return residuals_of(m, x, self.ctx, self.w)

    def n_data(self) -> int:
        """Rows of data in the residual vector (kept points and inclination rows)."""
        n = int(self.w["obs"][0].sum())
        if self.rv is not None:
            n += int(self.w["rv"][0].sum())
        n += sum(int(k.sum()) for k, _ in self.w["lc"])
        return n + sum(len(q) for q in self.incs)

    def n_calib(self) -> int:
        """Linear calibration coefficients of the light curves (per segment with kept points:
        zero point, scale and baseline polynomial, lcfit.calibrate): fitted, so they use up
        degrees of freedom like the free parameters."""
        n = 0
        for (k, _), q in zip(self.w["lc"], self.ctx["lc"], strict=True):
            segs = np.unique(np.asarray(q.segment)[k])
            n += len(segs) * lc_calib_terms(q)
        return n

    def n_inactive(self, r) -> int:
        """Inclination rows of one-sided limits that the model satisfies (residual 0): they
        carry no information at this solution and do not count as degrees of freedom."""
        n_inc = sum(len(q) for q in self.incs)
        if not n_inc:
            return 0
        lim = np.concatenate(
            [q.limit if q.limit is not None else np.full(len(q), "", object) for q in self.incs]
        )
        rows = np.asarray(r[self.n_data() - n_inc : self.n_data()])
        return int(np.sum((lim != "") & (rows == 0.0)))

    def counts(self) -> dict:
        return {
            "minima": int(self.w["obs"][0].sum()),
            "rv": int(self.w["rv"][0].sum()) if self.rv is not None else 0,
            "lc": sum(int(k.sum()) for k, _ in self.w["lc"]),
            "inc": sum(len(q) for q in self.incs),
        }


def _pool(processes: int, ctx: dict):
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor

    from .ensemble import resolve_processes

    n = resolve_processes(processes)  # <= 0: all cores, or SWINGERSIM_CORES of the job worker
    if n <= 1:
        return None, 1
    method = "fork" if "fork" in mp.get_all_start_methods() else "spawn"
    return (
        ProcessPoolExecutor(
            n, mp_context=mp.get_context(method), initializer=_init, initargs=(ctx,)
        ),
        n,
    )


def _fit_once(pb: _Problem, spec, x0, steps, pool, ncore, log, rec: dict) -> dict:
    """One least-squares fit with the current weights of pb, from x0."""
    names = pb.names
    lo, hi = default_bounds(names, x0, spec.bounds)
    x0 = np.clip(x0, lo, hi)
    last = {"key": None}

    def fun(x):
        r = pb.resid(pb.model(x), x)
        key = x.tobytes()
        if key != last["key"]:
            last["key"] = key
            rec["n_fun"] += 1
            chi2 = float(np.sum(r**2))
            rec["hist"].append({"eval": pb.n_model, "chi2": chi2, "x": x.tolist()})
            log(
                f"  evaluation {pb.n_model}: chi^2 = {chi2:.6g} (reduced "
                f"{chi2 / max(len(r) - len(x), 1):.4g})  "
                f"[{time.perf_counter() - rec['t_wall']:.0f} s]",
            )
        return r

    def jac(x):
        """Central differences (forward ones are biased by the curvature of the residuals,
        e.g. the cumulative effect of a period change); one-sided at a bound."""
        n = len(x)
        up = x + steps <= hi
        dn = x - steps >= lo
        xs = [x + np.eye(n)[k] * steps[k] for k in range(n) if up[k]]
        xs += [x - np.eye(n)[k] * steps[k] for k in range(n) if dn[k]]
        ms = pb.models(xs, pool)
        failed = [_model_failed(m, pb.ctx) for m in ms]
        cols = [pb.resid(m, xk) for m, xk in zip(ms, xs, strict=True)]
        it, itf = iter(cols), iter(failed)
        fp = [next(it) if up[k] else None for k in range(n)]
        fm = [next(it) if dn[k] else None for k in range(n)]
        bad_p = [next(itf) if up[k] else False for k in range(n)]
        bad_m = [next(itf) if dn[k] else False for k in range(n)]
        # a step on which the model fails (a star beyond its track, an unstable orbit, ...):
        # the edge of the model's validity acts as a bound - one-sided difference from the
        # other side, and the parameter is reported
        edge = {}
        for k in range(n):
            if up[k] and bad_p[k]:
                up[k], edge[names[k]] = False, "+"
            if dn[k] and bad_m[k]:
                dn[k] = False
                edge[names[k]] = "both" if names[k] in edge else "-"
        rec.setdefault("edges", {})[x.tobytes()] = edge
        f0 = fun(x) if not np.all(up & dn) else None
        J = []
        for k in range(n):
            if up[k] and dn[k]:
                J.append((fp[k] - fm[k]) / (2.0 * steps[k]))
            elif up[k]:
                J.append((fp[k] - f0) / steps[k])
            elif dn[k]:
                J.append((f0 - fm[k]) / steps[k])
            else:  # the model fails on both sides: no derivative
                J.append(np.zeros(len(f0)))
        return np.column_stack(J)

    c = pb.counts()
    log(
        f"stage 2: {len(names)} free parameters, {c['minima']} minima of "
        f"{', '.join(sorted(set(pb.obs.pair[pb.w['obs'][0]].tolist())))}"
        + (f", {c['rv']} radial velocities" if pb.rv is not None else "")
        + (f", {c['lc']} light-curve points" if pb.lcs else "")
        + f", {ncore} core(s)"
    )
    r0 = fun(x0)
    res = least_squares(
        fun,
        x0,
        jac=jac,
        bounds=(lo, hi),
        method="trf",
        x_scale="jac",
        loss=spec.loss,
        f_scale=spec.f_scale,
        max_nfev=spec.max_iterations,
        xtol=1e-8,
        # stop when chi^2 changes by less than 1e-5 relative or 0.05 absolute (far below
        # the change of 1 that moves a parameter by one sigma); the larger of the two
        ftol=max(1e-5, 0.05 / max(float(np.sum(r0**2)), 1.0)),
        gtol=1e-8,
    )
    x = res.x
    r = fun(x)
    # trf evaluates the Jacobian at every accepted point: with the linear loss res.jac is the
    # Jacobian at res.x; robust losses scale it by the loss, then the cached columns give it
    J = np.asarray(res.jac, float).copy() if spec.loss == "linear" else jac(x)
    edges = rec.get("edges", {}).get(x.tobytes())
    if edges is None:  # the Jacobian of res.x was not computed by jac (cached by trf)
        J = jac(x)
        edges = rec["edges"][x.tobytes()]
    for name, side in edges.items():
        log(
            f"  warning: {name}: the model fails "
            + {"+": "just above", "-": "just below", "both": "on both sides of"}[side]
            + " the fitted value - the fit stopped at the edge of the model's validity "
            "(e.g. a star beyond the end of its stellar track); its error is not valid"
        )
    chi2, chi2_red, dof = _chi2(pb, r, len(x))
    return {
        "x0": x0,
        "r0": r0,
        "x": x,
        "r": r,
        "J": J,
        "lo": lo,
        "hi": hi,
        "message": res.message,
        "edges": edges,
        "chi2": chi2,
        "chi2_red": chi2_red,
        "dof": dof,
    }


def _chi2(pb: _Problem, r, npar: int) -> tuple[float, float, int]:
    """(chi^2 of all rows, reduced chi^2 of the data, degrees of freedom): the prior rows and
    the satisfied one-sided inclination rows are not data, n_eff = kept data rows minus them."""
    n = pb.n_data()
    dof = max(n - pb.n_inactive(r) - npar - pb.n_calib(), 1)
    return float(np.sum(r**2)), float(np.sum(r[:n] ** 2)) / dof, dof


def _errors(pb: _Problem, f: dict, log) -> tuple[np.ndarray, np.ndarray]:
    """Covariance and errors of a fit: the data rows of the Jacobian scaled by
    sqrt(max(chi^2_red, 1)) (errors of all data inflated alike); the prior terms keep theirs."""
    Js = f["J"].copy()
    Js[: pb.n_data()] /= math.sqrt(max(f["chi2_red"], 1.0))
    cov = covariance(Js, log=log, names=pb.names)
    return cov, np.sqrt(np.clip(np.diag(cov), 0.0, None))


def covariance(J: np.ndarray, log=None, names=None) -> np.ndarray:
    """(J^T J)^-1 from the singular values of the column-scaled Jacobian. Forming J^T J squares
    the condition number; with parameters of very different scale (a period known to 1e-5 d,
    angles to 0.1 deg) and correlations of 0.999 a pseudo-inverse of J^T J silently drops the
    weakest direction and the errors come out far too small. Directions with a singular value
    below 1e-12 of the largest are not constrained by the data; they are left out (their errors
    come out too small) and a warning names the parameters that dominate them."""
    d = np.linalg.norm(J, axis=0)
    d[d == 0] = 1.0
    _, sv, vt = np.linalg.svd(J / d, full_matrices=False)
    keep = sv > sv[0] * 1e-12
    if not np.all(keep):
        worst = [int(np.argmax(np.abs(v))) for v in vt[~keep]]
        who = ", ".join(names[k] if names else str(k) for k in worst)
        msg = (
            f"covariance: {int(np.sum(~keep))} direction(s) of the parameters not constrained "
            f"by the data (singular Jacobian; mostly {who}) were dropped - their errors are "
            "too small"
        )
        if log is not None:
            log("  warning: " + msg)
        else:
            import warnings

            warnings.warn(msg, RuntimeWarning, stacklevel=2)
    inv = (vt[keep].T / sv[keep] ** 2) @ vt[keep]
    return inv / np.outer(d, d)


# ------------------------------------------------------------------ errors and outliers
OTHER = "other sources"


def source_groups(source, min_n: int) -> np.ndarray:
    """Group of every minimum for the error factors: its source, or OTHER for the sources with
    fewer than min_n minima (a factor of two or three minima is not determined)."""
    src = np.asarray(source).astype(str)
    names, counts = np.unique(src, return_counts=True)
    small = set(names[counts < min_n].tolist())
    return np.array([OTHER if x in small else x for x in src], dtype=object).astype(str)


def clip_outliers(z, groups, keep, clip_sigma: float, max_fraction: float):
    """Minima to leave out: |z| larger than clip_sigma robust standard deviations of their group
    (1.4826 x the median |z| of the kept minima of the group, at least 1, so that a source that
    fits better than its quoted errors does not lose good minima). The worst first; together
    with those already left out at most max_fraction of all. Returns (new keep, ratio)."""
    a = np.abs(np.asarray(z, float))
    a[~np.isfinite(a)] = np.inf
    ratio = np.zeros(len(a))
    for g in np.unique(groups):
        m = groups == g
        mk = m & keep
        if not np.any(mk):
            continue
        s = max(1.4826 * float(np.median(a[mk])) if np.isfinite(np.median(a[mk])) else 1.0, 1.0)
        ratio[m] = a[m] / s
    # a minimum left out against an earlier solution comes back when the present model fits it
    # (below REINCLUDE x clip_sigma: the margin keeps a point from flipping between the fits)
    new = keep.copy()
    new[~keep & (ratio < REINCLUDE * clip_sigma)] = True
    cand = np.where(new & (ratio > clip_sigma))[0]
    room = int(np.floor(max_fraction * len(a))) - int(np.sum(~new))
    cand = cand[np.argsort(-ratio[cand])][: max(room, 0)]
    new[cand] = False
    return new, ratio


def lc_calib_terms(q) -> int:
    """Calibration coefficients per light-curve segment (lcfit.calibrate)."""
    return int(bool(q.scale)) + 1 + int(q.baseline_degree)


def error_factors(
    resid_d,
    sigma_d,
    groups,
    how: str,
    model: str,
    npar: int,
    n_total: int | None = None,
    extra_dof: dict | None = None,
) -> dict:
    """Per group the factor (model scale) or the jitter [d] (model jitter) that makes the
    reduced chi^2 of the group one; the npar degrees of freedom of the fit are shared among the
    groups in proportion to their points, out of n_total points of all data sets fitted
    together (default: the points given here). per_source never goes below the quoted errors
    (factor >= 1, jitter >= 0); global may also shrink them (factor < 1)."""
    from scipy.optimize import brentq

    r = np.asarray(resid_d, float)
    s = np.asarray(sigma_d, float)
    n = len(r)
    n_all = max(int(n_total) if n_total is not None else n, n, 1)
    g = np.full(n, "all") if how == "global" else np.asarray(groups)
    out = {}
    for name in sorted(set(g.tolist())):
        m = (g == name) & np.isfinite(r)
        k = int(m.sum())
        if k == 0:
            continue
        dof = max(k - npar * k / n_all - (extra_dof or {}).get(name, 0), 1.0)
        rm, sm = r[m], s[m]
        if model == "scale":
            f = math.sqrt(float(np.sum((rm / sm) ** 2)) / dof)
            out[name] = max(f, 1.0) if how == "per_source" else f
        else:

            def excess(j, rm=rm, sm=sm, dof=dof):
                return float(np.sum(rm**2 / (sm**2 + j**2))) - dof

            if excess(0.0) <= 0.0:
                out[name] = 0.0
            else:
                hi = float(np.sqrt(np.sum(rm**2) / dof)) + 1e-12
                out[name] = float(brentq(excess, 0.0, hi * 1.01 + 1e-9))
    return out


def apply_errors(sigma_d, groups, how: str, model: str, factors: dict) -> np.ndarray:
    s = np.asarray(sigma_d, float)
    g = np.full(len(s), "all") if how == "global" else np.asarray(groups)
    v = np.array([factors.get(x, 1.0 if model == "scale" else 0.0) for x in g], float)
    return s * v if model == "scale" else np.sqrt(s**2 + v**2)


def _subset(obs: Observed, keep, sigma) -> Observed:
    return Observed(obs.pair[keep], obs.kind[keep], obs.time[keep], sigma[keep], obs.source[keep])


def _times_at(raw: dict, base_dir, obs: Observed, bjd_t0: float) -> np.ndarray:
    try:
        return model_times(_cfg_of(raw, base_dir), obs, bjd_t0)
    except (ValueError, RuntimeError):
        return np.full(len(obs), np.nan)


# ------------------------------------------------------------------ stage 2: fits in a row
def _new_state(pb: _Problem, spec, x0) -> dict:
    """State of the fits in a row (carried through the passes of the mapped-mesh fit): kept
    points and errors of every data set, error factors, start values and their model."""
    groups = source_groups(pb.obs.source, spec.min_source_minima)
    sets = [
        {
            "kind": "minima",
            "groups": groups,
            "keep": np.ones(len(pb.obs), bool),
            "sigma": pb.obs.sigma.copy(),
            "ratio": np.zeros(len(pb.obs)),
        }
    ]
    if pb.rv is not None:
        sets.append(
            {
                "kind": "rv",
                "groups": np.char.add("rv:", source_groups(pb.rv.source, spec.min_source_minima)),
                "keep": np.ones(len(pb.rv), bool),
                "sigma": pb.rv.sigma.copy(),
                "ratio": np.zeros(len(pb.rv)),
            }
        )
    for i, q in enumerate(pb.lcs):
        sets.append(
            {
                "kind": "lc",
                "i": i,
                "groups": np.array([f"lc:{q.source}:{g}" for g in q.segment]),
                "keep": np.ones(len(q), bool),
                "sigma": q.sigma.copy(),
                "ratio": np.zeros(len(q)),
            }
        )
    x0 = np.asarray(x0, float)
    lo, hi = default_bounds(pb.names, x0, spec.bounds)
    x0 = np.clip(x0, lo, hi)
    return {
        "sets": sets,
        "factors": {},
        "beta": {},
        "x0": x0,
        "x": x0.copy(),
        "steps": default_steps(pb.names, x0, spec.fd_steps),
        "start": None,  # model, calibrated light curves and chi^2 parts at x0 (first fit)
        "passes": [],
        "hist": [],
        "note": "",
        "n_fun": 0,
    }


def _set_weights(pb: _Problem, st: dict) -> None:
    for x in st["sets"]:
        if x["kind"] == "minima":
            pb.w["obs"] = (x["keep"], x["sigma"])
        elif x["kind"] == "rv":
            pb.w["rv"] = (x["keep"], x["sigma"])
        else:
            pb.w["lc"][x["i"]] = (x["keep"], x["sigma"])


def _data_resid(pb: _Problem, x: dict, m: dict) -> np.ndarray:
    """measured - model of every point of one data set (outliers included); light curves
    calibrated per segment with all points and their quoted errors."""
    if x["kind"] == "minima":
        return pb.obs.time - (m["t"] if m["t"] is not None else np.nan)
    if x["kind"] == "rv":
        return pb.rv.rv - (m["rv"] if m["rv"] is not None else np.nan)
    from .lcfit import calibrate

    q = pb.lcs[x["i"]]
    if m["lc"] is None:
        return np.full(len(q), np.nan)
    mm = m["lc"][x["i"]] + q.delta if _mapped(q) else m["lc"][x["i"]]
    try:
        return q.flux - calibrate(q, mm)[0]
    except (ValueError, RuntimeError):
        return np.full(len(q), np.nan)


def _quoted(pb: _Problem, x: dict) -> np.ndarray:
    if x["kind"] == "minima":
        return pb.obs.sigma
    if x["kind"] == "rv":
        return pb.rv.sigma
    return pb.lcs[x["i"]].sigma


def _run_fits(pb: _Problem, spec, how: str, robust: bool, st: dict, pool, ncore, log) -> dict:
    """Fit, then (robust) rescale the errors and leave out outliers, fit again from the
    solution, ... until neither the data left out nor the factors change (by more than 2 %) or
    fit.stage2.refits fits were done. Minima, radial velocities (groups rv:<source>) and light
    curves (groups lc:<file>:<segment>) are treated alike. Returns the last fit."""
    from .lcfit import calibrate

    clip = spec.clip_sigma
    names = pb.names
    rec = {"t_wall": time.perf_counter(), "hist": st["hist"], "n_fun": 0}
    f = None
    for n in range(spec.refits if robust else 1):
        _set_weights(pb, st)
        c = pb.counts()
        if c["minima"] + c["rv"] + c["lc"] <= len(names) and robust:
            raise ValueError(f"fit.stage2: {c['minima']} minima left for {len(names)} parameters")
        if robust:
            log(
                f"stage 2, fit {len(st['passes']) + 1}: "
                + ", ".join(
                    f"{int(x['keep'].sum())} "
                    + {"minima": "minima", "rv": "velocities", "lc": "LC points"}[x["kind"]]
                    + f" ({int(np.sum(~x['keep']))} left out)"
                    for x in st["sets"]
                )
            )
        f = _fit_once(pb, spec, st["x"], st["steps"], pool, ncore, log, rec)
        if st["start"] is None:  # the start values of the very first fit
            m0 = pb.model(f["x0"])
            lc0 = []
            for i, q in enumerate(pb.lcs):
                k, s = pb.w["lc"][i]
                full = np.full(len(q), np.nan)
                if m0["lc"] is not None:
                    mm = m0["lc"][i] + q.delta if _mapped(q) else m0["lc"][i]
                    full[k] = calibrate(q.subset(k, s), mm[k])[0]
                lc0.append(full)
            st["start"] = {"x0": f["x0"], "model": m0, "lc": lc0}
            st["x0"] = f["x0"]
        if robust:
            st["passes"].append(
                {
                    "fit": len(st["passes"]) + 1,
                    "minima": c["minima"],
                    "other_data": c["rv"] + c["lc"],
                    "left_out": sum(int(np.sum(~x["keep"])) for x in st["sets"]),
                    "chi2_red": f["chi2_red"],
                    "factors": dict(st["factors"]),
                }
            )
        st["x"] = f["x"]
        if not robust:
            break
        m = pb.model(f["x"])
        # 1. outliers of every data set
        new = []
        for x in st["sets"]:
            r = _data_resid(pb, x, m)
            nk = x["keep"]
            ratio = None
            if clip is not None:
                nk, ratio = clip_outliers(
                    r / x["sigma"], x["groups"], x["keep"], clip, spec.max_clip_fraction
                )
            new.append((x, r, nk, ratio))
        # 2. error factors: the degrees of freedom shared among all points fitted together
        n_total = sum(int(nk.sum()) for _, _, nk, _ in new) + sum(len(q) for q in pb.incs)
        new_f = dict(st["factors"])
        changed = False
        for x, r, nk, ratio in new:
            q0 = _quoted(pb, x)
            ns = q0.copy()
            if how != "none":
                # minima: global or per source; the others always per group (global would mix
                # days, km/s and flux)
                h = how if x["kind"] == "minima" else "per_source"
                # a light-curve segment also fits its own calibration coefficients
                extra = (
                    {g: lc_calib_terms(pb.lcs[x["i"]]) for g in np.unique(x["groups"])}
                    if x["kind"] == "lc"
                    else None
                )
                fr = error_factors(
                    r[nk], q0[nk], x["groups"][nk], h, spec.error_model, len(names), n_total,
                    extra,
                )  # fmt: skip
                new_f.update(fr)
                ns = apply_errors(q0, x["groups"], h, spec.error_model, fr)
            if x["kind"] == "lc" and getattr(pb.lcs[x["i"]], "red_noise", False):
                from .lcfit import red_noise_beta

                t = pb.lcs[x["i"]].time
                for g in np.unique(x["groups"]):
                    sel = (x["groups"] == g) & nk
                    b = red_noise_beta(t[sel], r[sel])
                    if abs(st["beta"].get(g, 0.0) - b) > 0.01:
                        log(f"  red-noise factor of {g}: errors x {b:.2f}")
                    st["beta"][g] = b
                    ns = np.where(x["groups"] == g, ns * b, ns)
            changed = (
                changed
                or bool(np.any(nk != x["keep"]))
                or bool(np.any(np.abs(ns[nk] / x["sigma"][nk] - 1.0) > 0.02))
            )
            x["new"] = (nk, ns, ratio)
        if not changed:
            break
        if n == spec.refits - 1:
            st["note"] = (
                f"the errors / outliers still changed after {spec.refits} fits "
                "(fit.stage2.refits); the last fit used the factors listed"
            )
            log("  warning: " + st["note"])
            break
        for g in sorted(new_f):
            v = new_f[g]
            if spec.error_model == "scale":
                txt = f"x {v:.3g}"
            elif g.startswith("rv:"):
                txt = f"jitter {v:.3g} km/s"
            elif g.startswith("lc:"):
                txt = f"jitter {v * 1e3:.3g} ppt"
            else:
                txt = f"jitter {v * 1440:.3g} min"
            log(f"  errors of {g}: {txt}")
        st["factors"] = new_f
        for x in st["sets"]:
            nk, ns, ratio = x.pop("new")
            if np.any(nk != x["keep"]):
                what = {"minima": "minima", "rv": "rv points", "lc": "lc points"}[x["kind"]]
                log(f"  {int(np.sum(x['keep'] & ~nk))} outlying {what} left out")
                x["ratio"][x["keep"] & ~nk] = ratio[x["keep"] & ~nk]
            x["keep"], x["sigma"] = nk, ns
    for x in st["sets"]:
        x.pop("new", None)
    st["n_fun"] += rec["n_fun"]
    return f


def _mcmc_info(pb: _Problem, f: dict, cov, err, spec, pool, ncore, log) -> dict | None:
    """emcee around the solution f with its covariance (fit.stage2.sampler: emcee)."""
    if spec.sampler != "emcee" or not np.all(np.isfinite(err)):
        return None
    return run_mcmc(
        pb,
        f["x"],
        err,
        cov,
        f["lo"],
        f["hi"],
        max(f["chi2_red"], 1.0),
        spec.walkers,
        spec.mcmc_steps,
        spec.seed,
        pool,
        ncore,
        log,
        tau_target=getattr(spec, "mcmc_tau_target", 0.0),
        max_steps=getattr(spec, "mcmc_max_steps", None),
    )


def _result(pb: _Problem, st: dict, f: dict, cov, err, mc, spec, how, t_wall) -> Stage2Result:
    """Stage2Result of the last fit f: data, errors and models of the kept points (start values
    of the first fit, fitted values of the last), the outliers with the fitted model."""
    from .lcfit import apply_calibration, calibrate

    sets = {x["kind"] if x["kind"] != "lc" else ("lc", x["i"]): x for x in st["sets"]}
    _set_weights(pb, st)
    x = f["x"]
    m1, m0 = pb.model(x), st["start"]["model"]
    keep, sig = pb.w["obs"]
    obs = pb.obs
    tm1 = m1["t"] if m1["t"] is not None else np.full(len(obs), np.nan)
    tm0 = m0["t"] if m0["t"] is not None else np.full(len(obs), np.nan)
    err_lin = err.copy()
    samples = None
    message = f["message"]
    if mc is not None:
        samples = mc["samples"]
        p16, p84 = np.percentile(samples, [16, 84], axis=0)
        err = 0.5 * (p84 - p16)
        message = f"{message}; emcee acceptance {mc['acceptance']:.2f}"
    # chi^2 of the start values with the final errors (light curves calibrated as in the first
    # fit): measured - start model of the kept points
    r0 = [(obs.time[keep] - tm0[keep]) / sig[keep]]
    res = Stage2Result(
        pb.names,
        st["x0"],
        x,
        err,
        cov,
        _subset(obs, keep, sig),
        tm0[keep],
        tm1[keep],
        0.0,
        f["chi2"],
        f["chi2_red"],
        pb.n_model,
        time.perf_counter() - t_wall,
        pb.ctx["raw"],
        pb.ctx["bjd_t0"],
        message,
        st["hist"],
        err_lin,
        samples,
        dict(pb.ctx["priors"]),
        float(np.sum(prior_residuals(pb.ctx["raw"], pb.names, x, pb.ctx["priors"]) ** 2)),
    )
    res.base_dir = pb.ctx["base_dir"]
    res.edges = dict(f.get("edges") or {})
    # angles of the outputs in 0-360 deg (the fit moves freely across the wrap); the samples
    # shifted alike
    res.x = np.array(res.x, float)
    for k, name in enumerate(res.names):
        if _split(name)[2] in ("Omega", "omega", "M0"):
            shift = 360.0 * math.floor(res.x[k] / 360.0)
            if shift:
                res.x[k] -= shift
                if res.samples is not None:
                    res.samples = np.array(res.samples, float)
                    res.samples[:, k] -= shift
    res.dof = f["dof"]
    res.loss = spec.loss
    res.sigma_in = obs.sigma[keep]
    res.error_scale, res.error_model, res.error_factors = how, spec.error_model, st["factors"]
    res.passes = st["passes"]
    res.note = st["note"]
    res.lc_beta = dict(st["beta"])
    res.mcmc = mc
    if mc is not None:
        res.samples_logp = mc["log_prob"]
    if np.any(~keep):
        res.outliers = _subset(obs, ~keep, obs.sigma)
        res.t_outliers = tm1[~keep]
        res.z_outliers = sets["minima"]["ratio"][~keep]
    if pb.rv is not None:
        k, s = pb.w["rv"]
        rv1 = m1["rv"] if m1["rv"] is not None else np.full(len(pb.rv), np.nan)
        rv0 = m0["rv"] if m0["rv"] is not None else np.full(len(pb.rv), np.nan)
        res.rv = pb.rv.subset(k, s)
        res.rv_start, res.rv_model = rv0[k], rv1[k]
        res.rv_sigma_in = pb.rv.sigma[k]
        r0.append((res.rv.rv - res.rv_start) / res.rv.sigma)
        if np.any(~k):
            res.rv_outliers = pb.rv.subset(~k)
            res.rv_model_outliers = rv1[~k]
            res.rv_outliers_z = sets["rv"]["ratio"][~k]
    if pb.lcs:
        res.lc, res.lc_start, res.lc_model, res.lc_coef, res.lc_sigma_in = [], [], [], [], []
        res.lc_outliers = []
        for i, q in enumerate(pb.lcs):
            k, s = pb.w["lc"][i]
            sub = q.subset(k, s)
            if m1["lc"] is not None:
                mm = m1["lc"][i] + q.delta if _mapped(q) else m1["lc"][i]
                cal, coef = calibrate(sub, mm[k])
            else:
                mm, cal, coef = np.full(len(q), np.nan), np.full(int(k.sum()), np.nan), {}
            res.lc.append(sub)
            res.lc_start.append(st["start"]["lc"][i][k])
            res.lc_model.append(cal)
            res.lc_coef.append(coef)
            res.lc_sigma_in.append(q.sigma[k])
            r0.append((sub.flux - res.lc_start[-1]) / sub.sigma)
            if np.any(~k):
                o = q.subset(~k)
                mids = {
                    str(g): 0.5
                    * (sub.time[sub.segment == g].min() + sub.time[sub.segment == g].max())
                    for g in np.unique(sub.segment)
                }
                res.lc_outliers.append(
                    (o, apply_calibration(o, mm[~k], coef, mids), sets[("lc", i)]["ratio"][~k])
                )
            else:
                res.lc_outliers.append(None)
    if pb.incs:
        from .incfit import residuals as inc_residuals

        res.inc = pb.incs
        nan = [np.full(len(q), np.nan) for q in pb.incs]
        res.inc_start = m0["inc"] if m0["inc"] is not None else nan
        res.inc_model = m1["inc"] if m1["inc"] is not None else nan
        r0 += [inc_residuals(q, m) for q, m in zip(pb.incs, res.inc_start, strict=True)]
    r0 = np.concatenate([np.asarray(v, float) for v in r0])
    r0[~np.isfinite(r0)] = 1e3
    pr0 = prior_residuals(pb.ctx["raw"], pb.names, st["x0"], pb.ctx["priors"])
    res.chi2_start = float(np.sum(r0**2) + np.sum(pr0**2))
    return res


def _problem(raw, base_dir, obs, spec, rv, lc) -> _Problem:
    n = len(spec.free)
    return _Problem(
        raw,
        spec.free,
        base_dir,
        spec.bjd_t0,
        dict(spec.priors),
        obs,
        rv,
        lc,
        list(_EXTRA.get("inc") or []),
        cache_size=4 * n + 12,
    )


def _stage2(pb, spec, how, robust, mcmc, processes, log, x0=None) -> Stage2Result:
    t_wall = time.perf_counter()
    _init(pb.ctx)
    pool, ncore = _pool(spec.processes if processes is None else processes, pb.ctx)
    try:
        st = _new_state(pb, spec, start_values(pb.ctx["raw"], pb.names) if x0 is None else x0)
        f = _run_fits(pb, spec, how, robust, st, pool, ncore, log)
        cov, err = _errors(pb, f, log)
        mc = _mcmc_info(pb, f, cov, err, spec, pool, ncore, log) if mcmc else None
    finally:
        if pool is not None:
            pool.shutdown()
    return _result(pb, st, f, cov, err, mc, spec, how, t_wall)


def fit_stage2(
    raw: dict,
    base_dir: str | None,
    obs: Observed,
    spec,
    processes: int | None = None,
    log=print,
    mcmc: bool = True,
    rv=None,
    lc=None,
) -> Stage2Result:
    """One weighted least-squares fit of all given data (no error rescaling or outlier
    rejection); emcee afterwards with fit.stage2.sampler: emcee and mcmc."""
    pb = _problem(raw, base_dir, obs, spec, rv, lc)
    return _stage2(pb, spec, "none", False, mcmc, processes, log)


def fit_stage2_robust(
    raw: dict,
    base_dir: str | None,
    obs: Observed,
    spec,
    how: str = "none",
    processes: int | None = None,
    log=print,
    rv=None,
    lc=None,
) -> Stage2Result:
    """Stage 2 with error factors per source and outlier rejection: fit, then rescale the
    errors and leave out outliers, fit again from the solution, ... until neither the data
    left out nor the factors change (by more than 2 %) or fit.stage2.refits fits were done.
    Minima, radial velocities (groups rv:<source>) and light curves (groups
    lc:<file>:<segment>) are treated alike. Without error_scale and clip_sigma it is one plain
    fit. The model of every parameter set is computed once (start, solution and Jacobian are
    reused by the next fit in a row); emcee (fit.stage2.sampler) runs once around the last
    fit."""
    robust = not (how == "none" and spec.clip_sigma is None)
    lc = list(lc or [])
    if not robust and any(getattr(q, "red_noise", False) for q in lc):
        log(
            "  warning: fit.lc red_noise needs error_scale or clip_sigma (the fits in a row that "
            "estimate it); it is ignored"
        )
    pb = _problem(raw, base_dir, obs, spec, rv, lc)
    return _stage2(pb, spec, how, robust, True, processes, log)


def run_mcmc(
    pb,
    x,
    err,
    cov,
    lo,
    hi,
    scale,
    walkers,
    steps,
    seed,
    pool,
    ncore,
    log,
    tau_target: float = 0.0,
    max_steps: int | None = None,
) -> dict:
    """MCMC (emcee) around the least-squares solution for errors that do not rely on the
    linear model: in a long, curved valley of chi^2 - typical for an outer orbit covered by
    less than about two periods, where M0, omega and P trade against each other - the errors
    of the covariance matrix come out too small. chi^2 of the data is divided by
    max(chi^2_red, 1) (the priors keep theirs); the first half of the chain is discarded.

    The walkers are at least 2 x parameters + 2 and even (emcee moves half of them at a time);
    their number does not depend on the computer, so the same seed gives the same chain on any
    number of cores (cores only share the models of a half-step). The log-probability
    of all walkers of a half-step is computed at once (emcee vectorize), their models in
    parallel. With tau_target > 0 the chain is extended after `steps` (by a quarter of its
    length at a time) until it is at least tau_target integrated autocorrelation times long
    (the emcee documentation advises at least 50) or reaches max_steps. Returns {samples, log_prob, acceptance, tau
    (integrated autocorrelation time per parameter, steps), walkers, steps, discard, seed,
    warning}."""
    import emcee

    rng = np.random.default_rng(seed)
    ndim = len(x)
    walkers = max(int(walkers), 2 * ndim + 2)
    walkers += walkers % 2
    try:  # a ball of half the covariance: the walkers start spread along the valley
        L = np.linalg.cholesky(cov + np.diag((1e-9 * np.where(err > 0, err, 1.0)) ** 2))
        p0 = x + 0.5 * rng.normal(size=(walkers, ndim)) @ L.T
    except np.linalg.LinAlgError:
        p0 = x + 0.1 * np.where(err > 0, err, 1e-9) * rng.normal(size=(walkers, ndim))
    inner_lo, inner_hi = lo.copy(), hi.copy()
    fl, fh = np.isfinite(lo), np.isfinite(hi)
    inner_lo[fl] += 1e-9 * np.abs(lo[fl])
    inner_hi[fh] -= 1e-9 * np.abs(hi[fh])
    p0 = np.clip(p0, inner_lo, inner_hi)
    n_data = pb.n_data()

    # orbits fitted with ecosw and esinw: a flat prior in (e cos w, e sin w) is a prior on e
    # proportional to e (area element e de dw), which pushes nearly circular orbits away from
    # e = 0 (Lucy & Sweeney 1971; Ford 2006; Eastman et al. 2013). The Jacobian 1/e makes the
    # prior flat in e and w, as for the orbits fitted with e and omega.
    ew_pairs = []
    names = list(getattr(pb, "names", ()))
    for k, nm in enumerate(names):
        if nm.endswith(".ecosw") and nm[: -len("ecosw")] + "esinw" in names:
            ew_pairs.append((k, names.index(nm[: -len("ecosw")] + "esinw")))

    def lnprob(X):
        X = np.atleast_2d(X)
        out = np.full(len(X), -np.inf)
        ok = np.all((X >= lo) & (X <= hi), axis=1)
        for a, b in ew_pairs:  # e beyond its bound: outside the prior
            ok &= np.hypot(X[:, a], X[:, b]) < E_MAX
        idx = np.nonzero(ok)[0]
        for i, m in zip(idx, pb.models(X[idx], pool, store=False), strict=True):
            r = pb.resid(m, X[i])
            out[i] = -0.5 * (
                float(np.sum(r[:n_data] ** 2)) / scale + float(np.sum(r[n_data:] ** 2))
            )
            for a, b in ew_pairs:
                out[i] -= math.log(max(math.hypot(X[i, a], X[i, b]), 1e-300))
        return out

    log(
        f"emcee: {walkers} walkers x {steps} steps = {walkers * steps} model evaluations "
        f"({walkers // 2} per half-step on {max(int(ncore), 1)} core(s))"
    )
    sampler = emcee.EnsembleSampler(walkers, ndim, lnprob, vectorize=True)
    # the moves of emcee draw from the sampler's own generator (by default a copy of the global
    # numpy state): seeded, the chain is reproducible
    sampler.random_state = np.random.RandomState(seed).get_state()
    t0 = time.perf_counter()
    max_steps = max(int(max_steps or steps), int(steps)) if tau_target > 0 else int(steps)

    def tau_of(n: int) -> tuple[np.ndarray, float]:
        t = np.asarray(sampler.get_autocorr_time(discard=n // 2, quiet=True), float)
        return t, float(np.nanmax(t)) if np.any(np.isfinite(t)) else float("nan")

    state, done, chunk = p0, 0, int(steps)
    while True:
        start = state
        for k, state in enumerate(  # the last state continues the chain
            sampler.sample(start, iterations=chunk, skip_initial_state_check=True),
            start=done + 1,
        ):
            if k % max(steps // 20, 1) == 0 or k == done + chunk:
                log(
                    f"  emcee step {k}/{max(done + chunk, steps)}, acceptance "
                    f"{np.mean(sampler.acceptance_fraction):.2f}  [{time.perf_counter() - t0:.0f} s]"
                )
        done += chunk
        tau, tau_max = tau_of(done)
        if tau_target <= 0 or done >= max_steps:
            break
        if np.isfinite(tau_max) and done >= tau_target * tau_max:
            break
        chunk = min(max(done // 4, 20), max_steps - done)
        log(
            f"  emcee: chain {done} steps = {done / tau_max if tau_max > 0 else float('nan'):.1f} "
            f"tau (target {tau_target:g}): {chunk} more steps"
        )
    steps = done
    discard = steps // 2
    warning = ""
    if not np.isfinite(tau_max) or (steps - discard) < 10.0 * tau_max:
        warning = (
            f"emcee: the chain after burn-in ({steps - discard} steps) is shorter than 10 "
            f"autocorrelation times (tau max {tau_max:.1f} steps) - the errors are not "
            "converged; increase fit.stage2.mcmc_steps or set fit.stage2.mcmc_tau_target: 50"
        )
        log("  warning: " + warning)
    log(
        f"  emcee: autocorrelation time tau max {tau_max:.1f} steps, "
        f"steps / tau = {steps / tau_max if tau_max > 0 else float('nan'):.1f}"
    )
    return {
        "samples": sampler.get_chain(discard=discard, flat=True),
        "log_prob": sampler.get_log_prob(discard=discard, flat=True),
        "acceptance": float(np.mean(sampler.acceptance_fraction)),
        "tau": tau.tolist(),
        "tau_max": tau_max,
        "walkers": walkers,
        "steps": int(steps),
        "discard": discard,
        "seed": seed,
        "warning": warning,
    }


# ------------------------------------------------------------------ outputs
def _unit(name: str) -> str:
    return UNITS.get(_split(name)[2], "")


def metadata(res: Stage2Result) -> dict:
    """How the fit was made (report, header of stage2_fitted.yaml)."""
    from . import __version__

    mc = res.mcmc or {}
    return {
        "swingersim": __version__,
        "time_system": "BJD_TDB",
        "bjd_t0": res.bjd_t0,
        "loss": res.loss,
        "error_scale": res.error_scale,
        "error_model": res.error_model,
        "sampler": "emcee" if mc else "none",
        "walkers": mc.get("walkers"),
        "steps": mc.get("steps"),
        "discard": mc.get("discard"),
        "seed": mc.get("seed"),
        "tau_max": mc.get("tau_max"),
    }


def metadata_lines(res: Stage2Result) -> list[str]:
    md = metadata(res)
    out = [
        (
            f"SWINGERsim {md['swingersim']}; times BJD_TDB, t = 0 at BJD_TDB {md['bjd_t0']}; loss "
            f"{md['loss']}; errors {md['error_scale']} ({md['error_model']})"
        )
    ]
    if md["sampler"] == "emcee":
        tau = md["tau_max"]
        out.append(
            f"sampler emcee: {md['walkers']} walkers x {md['steps']} steps, first "
            f"{md['discard']} discarded, seed {md['seed']}, autocorrelation time tau max "
            f"{tau:.1f} steps (steps / tau {md['steps'] / tau if tau and tau > 0 else float('nan'):.1f})"
        )
    else:
        out.append("sampler none (linear errors of the covariance matrix)")
    return out


def report(res: Stage2Result) -> str:
    lines = [
        "swingersim fit, stage 2 - joint dynamical model (N-body) of all eclipsing pairs",
        f"  t = 0 of the configuration: BJD {res.bjd_t0}",
        (
            f"  {len(res.obs)} minima, "
            + (f"{len(res.rv)} radial velocities, " if res.rv is not None else "")
            + (f"{sum(len(q) for q in res.lc)} light-curve points, " if res.lc else "")
            + f"{len(res.priors)} priors, {len(res.names)} free "
            f"parameters, {res.n_eval} model evaluations, {res.wall_s:.0f} s"
        ),
        (
            f"  chi^2 {res.chi2_start:.6g} -> {res.chi2:.6g} (reduced {res.chi2_red:.4g}"
            + (f" for {res.dof} degrees of freedom" if res.dof else "")
            + (f"; priors {res.chi2_priors:.4g}" if res.priors else "")
            + "); errors of all data scaled by sqrt(max(chi^2_red, 1)) = "
            f"{math.sqrt(max(res.chi2_red, 1.0)):.3g}"
        ),
        (
            "  (chi^2_red: chi^2 of the data without the priors / (kept data rows - satisfied "
            "one-sided inclination rows - free parameters))"
        ),
        f"  optimiser: {res.converged}",
    ]
    lines += ["  " + x for x in metadata_lines(res)]
    if res.fidelity_note:
        lines.append("  mapped meshes: " + res.fidelity_note)
    if res.mcmc and res.mcmc.get("warning"):
        lines.append("  warning: " + res.mcmc["warning"])
    if res.priors:
        lines.append("  priors (mean +- sigma -> fitted value, deviation in sigma):")
        pr = prior_residuals(res.raw, res.names, res.x, res.priors)
        vals = start_values(apply(res.raw, res.names, res.x), list(res.priors))
        for k, (name, (mu, sig)) in enumerate(res.priors.items()):
            lines.append(
                f"    {name:<22} {mu:.8g} +- {sig:.3g} -> {vals[k]:.8g} {_unit(name)}"
                f"  ({pr[k]:+.2f})"
            )
    if res.error_scale != "none" and res.error_factors:
        if res.error_model == "scale":
            what = "factors of the quoted errors"
            vals = [f"{g} x {v:.3g}" for g, v in res.error_factors.items()]
        else:
            what = "jitter added in quadrature to the quoted errors"
            vals = [
                f"{g} {v:.3g} km/s"
                if g.startswith("rv:")
                else f"{g} {v * 1e3:.3g} ppt"
                if g.startswith("lc:")
                else f"{g} {v * 1440.0:.3g} min"
                for g, v in res.error_factors.items()
            ]
        lines.append(f"  errors ({res.error_scale}), {what}: " + "; ".join(vals))
    if res.lc_beta:
        lines.append(
            "  light curves, red-noise factor beta (errors multiplied by it): "
            + "; ".join(f"{g} {b:.2f}" for g, b in res.lc_beta.items())
        )
    if len(res.passes) > 1:
        lines.append(
            "  fits in a row: "
            + "; ".join(
                f"{p['fit']}: {p['minima']} minima, chi^2_red {p['chi2_red']:.4g}"
                for p in res.passes
            )
        )
    if res.note:
        lines.append("  warning: " + res.note)
    for x in getattr(res, "notes", []):
        lines.append("  note: " + x)
    if res.outliers is not None and len(res.outliers):
        lines.append(
            f"  {len(res.outliers)} outlying minima left out (deviation in robust sigma of "
            "the source):"
        )
        for j in range(len(res.outliers)):
            d = (res.outliers.time[j] - res.t_outliers[j]) * 1440.0
            lines.append(
                f"    {res.outliers.pair[j]} {'p' if res.outliers.kind[j] == 0 else 's'} "
                f"BJD {res.outliers.time[j]:.5f}  {d:+9.2f} min  {res.z_outliers[j]:6.1f}  "
                f"{res.outliers.source[j]}"
            )
    lines.append("  rms of the residuals [min], start -> fitted:")
    r0, r1 = res.rms_min(res.t_start), res.rms_min(res.t_model)
    for p in r0:
        n = int(np.sum(res.obs.pair == p))
        lines.append(f"    {p:<8} {n:4d} minima  {r0[p]:10.3f} -> {r1[p]:8.3f}")
    split = primary_secondary_split(res.obs, res.t_model)
    if split:
        lines.append(
            "  primary vs secondary minima, mean residual [min] (+- error of the mean), fitted:"
        )
        for p, d in split.items():
            lines.append(
                f"    {p:<8} primary {d['mean_p']:+8.2f} +- {d['err_p']:.2f} ({d['n_p']}), "
                f"secondary {d['mean_s']:+8.2f} +- {d['err_s']:.2f} ({d['n_s']}): "
                f"difference {d['diff']:+.2f} min ({d['z']:+.1f} sigma)"
            )
            if d["warn"]:
                lines.append(
                    f"    warning: the primary and secondary minima of {p} are offset against "
                    "each other - e cos(omega) of the pair is probably wrong (a local minimum: "
                    "restart the fit with e cos(omega) of the other sign or larger e); the error "
                    "floors per source may hide it"
                )
    for k, q in enumerate(res.inc or []):
        from .incfit import model_i, model_limit
        from .incfit import residuals as inc_residuals

        m0, m1 = res.inc_start[k], res.inc_model[k]
        z0, z1 = inc_residuals(q, m0), inc_residuals(q, m1)
        lines.append(
            f"  inclinations {q.label} (pair {q.pair}), {len(q)} epochs: chi^2 "
            f"{float(np.nansum(z0**2)):.1f} -> {float(np.nansum(z1**2)):.1f}"
        )
        i1, l1 = model_i(m1), model_limit(m1)
        lim = q.limit if q.limit is not None else [""] * len(q)
        for j in range(len(q)):
            kind = lim[j]
            if kind in ("no_eclipse", "eclipse"):
                obs = f"{'no eclipse' if kind == 'no_eclipse' else 'eclipse':>10}"
                bound = f"limit {l1[j]:6.2f}"
            else:
                sym = {"upper": "<=", "lower": ">="}.get(kind, "  ")
                obs = f"i {sym}{q.inc[j]:6.2f} +- {q.sigma[j]:.2f}"
                bound = ""
            ok = " (satisfied)" if kind and z1[j] == 0 else ""
            lines.append(
                f"    BJD {q.time[j]:.1f}  {obs}  model {i1[j]:6.2f} {bound} "
                f"({z1[j]:+.1f} sigma){ok}  {q.source[j]}"
            )
    if res.rv is not None:
        lines.append("  radial velocities, rms of the residuals [km/s], start -> fitted:")
        v0, v1 = res.rms_rv(res.rv_start), res.rms_rv(res.rv_model)
        for s in v1:
            n = int(np.sum(res.rv.star == s))
            lines.append(f"    {s:<8} {n:4d} velocities {v0.get(s, np.nan):8.2f} -> {v1[s]:8.2f}")
        r = res.rv.rv - res.rv_model
        lines.append(
            "    mean residual per source [km/s]: "
            + "; ".join(
                f"{src} {np.nanmean(r[res.rv.source == src]):+.2f}"
                for src in sorted(set(res.rv.source.tolist()))
            )
        )
        if res.rv_outliers is not None and len(res.rv_outliers):
            o = res.rv_outliers
            lines.append(
                f"    {len(o)} outlying velocities left out: "
                + "; ".join(
                    f"{o.star[k]} BJD {o.time[k]:.4f} {o.rv[k] - res.rv_model_outliers[k]:+.1f}"
                    for k in range(len(o))
                )
            )
    for x in res.shape_notes:
        lines.append(f"  star shapes: {x}")
    if res.fidelity:
        lines.append(
            "  mapped meshes (coarse + accurate - coarse), passes: "
            + "; ".join(
                f"{p['pass']}: chi^2_red {p['chi2_red']:.4g}, difference changed by "
                f"{p['change_sigma']:.2f} sigma"
                for p in res.fidelity
            )
        )
    if res.lc:
        lines.append(
            "  light curves, rms of measured - calibrated model [ppt], start -> fitted, scale a:"
        )
        r0, r1 = res.rms_lc("start"), res.rms_lc("model")
        for q, coef in zip(res.lc, res.lc_coef, strict=False):
            for seg in np.unique(q.segment):
                key = (q.source, str(seg))
                n = int(np.sum(q.segment == seg))
                a = coef.get(str(seg), [np.nan])[0] if q.scale else 1.0
                lines.append(
                    f"    {q.source}:{seg:<6} {n:5d} points {r0.get(key, np.nan):8.2f} -> "
                    f"{r1[key]:7.2f}   a = {a:.3f}"
                )
        n_out = sum(len(o[0]) for o in res.lc_outliers or [] if o is not None)
        if n_out:
            lines.append(
                f"    {n_out} outlying light-curve points left out (status outlier in stage2_lc.csv)"
            )
    mc = res.samples is not None
    lines.append(
        "  parameters (start -> fitted +- 1 sigma"
        + (
            "; errors: half the 16-84 % range of the emcee samples, in brackets the linear "
            "error of the covariance matrix):"
            if mc
            else "; linear errors of the covariance matrix):"
        )
    )
    for k, name in enumerate(res.names):
        lin = f" (linear {res.err_linear[k]:.3g})" if mc else ""
        lines.append(
            f"    {name:<22} {res.x0[k]:.10g} -> {res.x[k]:.10g} +- {res.err[k]:.3g}"
            f" {_unit(name)}{lin}"
        )
    for name, side in getattr(res, "edges", {}).items():
        where = {"+": "just above", "-": "just below", "both": "on both sides of"}[side]
        lines.append(
            f"  warning: {name}: the model fails {where} the fitted value - the fit stopped at "
            "the edge of the model's validity (e.g. a star beyond the end of its stellar track at "
            "the configured age); the data pull the parameter further, its error is not valid"
        )
    el = res.err_linear if res.err_linear is not None else res.err
    corr = res.cov / np.outer(el, el) if np.all(el > 0) else None
    if corr is not None:
        high = [
            f"{res.names[i]} / {res.names[j]} {corr[i, j]:+.3f}"
            for i in range(len(res.names))
            for j in range(i + 1, len(res.names))
            if abs(corr[i, j]) > 0.95
        ]
        if high:
            lines.append("  strongly correlated parameters (|r| > 0.95): " + "; ".join(high))
            if not mc:
                lines.append(
                    "  -> their errors assume a straight valley of chi^2 (linear model); "
                    "fit.stage2.sampler: emcee checks them"
                )
    return "\n".join(lines)


EPOCH_KEY = "fit_epoch_bjd_tdb"  # top-level key of the configuration: BJD_TDB of t = 0


def fitted_config(res: Stage2Result) -> dict:
    d = apply(res.raw, res.names, res.x)
    d.pop("fit", None)
    return d


def primary_secondary_split(obs: Observed, t_model, min_n: int = 3) -> dict:
    """Mean measured - model residual of the primary and the secondary minima of every pair
    [min]. A significant difference (> 4 sigma and > 1 min) means that the model puts the
    secondary eclipses at the wrong phase, i.e. e cos(omega) of the pair is off - easy to miss
    when per-source error floors absorb it (V907 Sco, SQUADRA test part III)."""
    out = {}
    r = (obs.time - np.asarray(t_model, float)) * 1440.0
    for p in sorted(set(obs.pair.tolist())):
        sel = (obs.pair == p) & np.isfinite(r)
        rp, rs = r[sel & (obs.kind == 0)], r[sel & (obs.kind == 1)]
        if len(rp) < min_n or len(rs) < min_n:
            continue
        ep = float(np.std(rp, ddof=1) / math.sqrt(len(rp)))
        es = float(np.std(rs, ddof=1) / math.sqrt(len(rs)))
        diff = float(np.mean(rp) - np.mean(rs))
        err = math.hypot(ep, es)
        z = diff / err if err > 0 else 0.0
        out[p] = {
            "mean_p": float(np.mean(rp)),
            "err_p": ep,
            "n_p": len(rp),
            "mean_s": float(np.mean(rs)),
            "err_s": es,
            "n_s": len(rs),
            "diff": diff,
            "z": float(z),
            "warn": bool(abs(z) > 4.0 and abs(diff) > 1.0),
        }
    return out


def write_inc_table(res: Stage2Result, path: Path) -> Path:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        from .incfit import model_i, model_limit
        from .incfit import residuals as inc_residuals

        w.writerow(
            [
                "label",
                "pair",
                "time",
                "i_obs",
                "sigma",
                "limit",
                "i_start",
                "i_model",
                "i_eclipse_limit",
                "z",
                "source",
            ]
        )
        for k, q in enumerate(res.inc):
            i0, i1 = model_i(res.inc_start[k]), model_i(res.inc_model[k])
            l1, z1 = model_limit(res.inc_model[k]), inc_residuals(q, res.inc_model[k])
            lim = q.limit if q.limit is not None else [""] * len(q)
            for j in range(len(q)):
                w.writerow(
                    [
                        q.label,
                        q.pair,
                        f"{q.time[j]:.3f}",
                        "" if not np.isfinite(q.inc[j]) else q.inc[j],
                        q.sigma[j],
                        lim[j],
                        f"{i0[j]:.4f}",
                        f"{i1[j]:.4f}",
                        "" if not np.isfinite(l1[j]) else f"{l1[j]:.4f}",
                        f"{z1[j]:.3f}",
                        q.source[j],
                    ]
                )
    return path


def plot_inclination(res: Stage2Result, path: Path) -> Path:
    """Measured inclinations and the model i(t) of the pair over the data span (+ 30 yr)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from .incfit import model_inclination

    cfg = _cfg_of(apply(res.raw, res.names, res.x), res.base_dir)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    for k, q in enumerate(res.inc):
        t0, t1 = float(q.time.min()) - 3650.0, float(q.time.max()) + 30 * 365.25
        tg = np.linspace(t0, t1, 600)
        try:
            ig, lg = model_inclination(cfg, q.pair, tg, res.bjd_t0, with_limit=True)
        except (ValueError, RuntimeError):
            ig = lg = np.full(len(tg), np.nan)
        yr = 2000.0 + (tg - J2000_BJD) / DPY
        ax.plot(yr, ig, "-", color="C0", lw=1.5, label=f"model {q.pair}" if k == 0 else None)
        if k == 0:
            ax.fill_between(yr, 0.0, lg, color="0.85", lw=0, label="no eclipses (model limit)")
        lim = q.limit if q.limit is not None else np.full(len(q), "", dtype=object)
        ty = 2000.0 + (q.time - J2000_BJD) / DPY
        m = lim == ""
        ax.errorbar(ty[m], q.inc[m], q.sigma[m], fmt="o", color="C3", ms=4, label=q.label)
        for kind, mk, txt in (("upper", "v", "upper limit"), ("lower", "^", "lower limit")):
            mm = lim == kind
            if np.any(mm):
                ax.plot(ty[mm], q.inc[mm], mk, color="C3", ms=7, mfc="none", label=txt)
        for kind, mk, txt in (("no_eclipse", "x", "no eclipse"), ("eclipse", "+", "eclipse")):
            mm = lim == kind
            if np.any(mm):
                yy = np.interp(q.time[mm], tg, lg)
                ax.plot(ty[mm], yy, mk, color="C2", ms=8, mew=1.5, label=txt)
    ax.set_xlabel("year")
    ax.set_ylabel("inclination [deg] (folded to 0-90)")
    lo = np.nanmin([np.nanmin(q.inc) if np.any(np.isfinite(q.inc)) else 90.0 for q in res.inc])
    ax.set_ylim(max(0.0, lo - 25.0), 91.0)
    ax.legend(frameon=False)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return path


def write_outputs(res: Stage2Result, out: Path, source_name: str, plots: bool = True) -> list[Path]:
    out.mkdir(parents=True, exist_ok=True)
    files = []
    p = out / "stage2_params.csv"
    lin = res.err_linear if res.err_linear is not None else res.err
    pct = None
    if res.samples is not None and len(res.samples):
        pct = np.percentile(res.samples, [16, 50, 84], axis=0)
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "parameter",
                "start",
                "value",
                "error",
                "unit",
                "prior_mean",
                "prior_sigma",
                "error_linear",
                "p16",
                "p50",
                "p84",
            ]
        )
        for k, name in enumerate(res.names):
            mu, sig = res.priors.get(name, ("", ""))
            w.writerow(
                [
                    name,
                    repr(float(res.x0[k])),
                    repr(float(res.x[k])),
                    f"{res.err[k]:.6g}",
                    _unit(name),
                    mu,
                    sig,
                    f"{lin[k]:.6g}",
                    *([repr(float(v)) for v in pct[:, k]] if pct is not None else ["", "", ""]),
                ]
            )
    files.append(p)
    p = out / "stage2_covariance.csv"  # of the linear model (covariance matrix of the fit)
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["parameter", *res.names])
        for k, name in enumerate(res.names):
            w.writerow([name, *(repr(float(v)) for v in res.cov[k])])
    files.append(p)
    p = out / "stage2_residuals.csv"
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "pair",
                "type",
                "bjd_tdb",
                "sigma_d",
                "model_start_bjd",
                "model_bjd",
                "o_minus_model_min",
                "source",
                "sigma_quoted_d",
                "status",
            ]
        )
        sig_in = res.sigma_in if res.sigma_in is not None else res.obs.sigma
        for j in range(len(res.obs)):
            w.writerow(
                [
                    res.obs.pair[j],
                    "p" if res.obs.kind[j] == 0 else "s",
                    f"{res.obs.time[j]:.6f}",
                    f"{res.obs.sigma[j]:.6g}",
                    f"{res.t_start[j]:.6f}",
                    f"{res.t_model[j]:.6f}",
                    f"{(res.obs.time[j] - res.t_model[j]) * 1440.0:.4f}",
                    res.obs.source[j],
                    f"{sig_in[j]:.6g}",
                    "fitted",
                ]
            )
        if res.outliers is not None:
            o = res.outliers
            for j in range(len(o)):
                w.writerow(
                    [
                        o.pair[j],
                        "p" if o.kind[j] == 0 else "s",
                        f"{o.time[j]:.6f}",
                        f"{o.sigma[j]:.6g}",
                        "",
                        f"{res.t_outliers[j]:.6f}",
                        f"{(o.time[j] - res.t_outliers[j]) * 1440.0:.4f}",
                        o.source[j],
                        f"{o.sigma[j]:.6g}",
                        "outlier",
                    ]
                )
    files.append(p)
    from .io import software_versions

    rep = report(res)
    rms = res.rms_min(res.t_model)
    head = [
        f"# Fitted by SWINGERsim (swingersim fit, stage 2, joint N-body model) from {source_name}:",
        f"# {len(res.obs)} minima, chi^2_red {res.chi2_red:.4g}; rms "
        + ", ".join(f"{k} {v:.2f} min" for k, v in rms.items()),
        "# free parameters: " + ", ".join(res.names),
        f"# t = 0 of this configuration = BJD {res.bjd_t0}.",
    ]
    if res.inc:
        files.append(write_inc_table(res, out / "stage2_inclination.csv"))
        if plots:
            files.append(plot_inclination(res, out / "stage2_inclination.png"))
    if res.lc:
        files.append(write_lc_table(res, out / "stage2_lc.csv"))
    if res.rv is not None:
        files.append(write_rv_table(res, out / "stage2_rv.csv"))
        fr = (apply(res.raw, res.names, res.x).get("fit") or {}).get("rv") or {}
        head.append(
            f"# radial velocities: gamma {float(fr.get('gamma', 0.0)):.3f} km/s"
            + "".join(f", offset {k} {v:+.3f}" for k, v in (fr.get("offsets") or {}).items())
            + "; rms "
            + ", ".join(f"{k} {v:.2f}" for k, v in res.rms_rv(res.rv_model).items())
            + " km/s"
        )
    head += ["# " + x for x in metadata_lines(res)]
    head.append(f"# software: {software_versions()}")
    fitted = fitted_config(res)
    if EPOCH_KEY in Config.model_fields:  # the epoch as a key of the configuration
        fitted = {EPOCH_KEY: float(res.bjd_t0), **fitted}
    p = out / "stage2_fitted.yaml"
    p.write_text("\n".join(head) + "\n" + yaml.safe_dump(fitted, sort_keys=False))
    files.append(p)
    p = out / "stage2_report.txt"
    p.write_text(rep + "\n")
    files.append(p)
    if res.samples is not None:
        p = out / "stage2_samples.csv"
        thin = max(len(res.samples) // 4000, 1)
        cols = [res.samples[::thin]]
        head_s = list(res.names)
        if res.samples_logp is not None and len(res.samples_logp) == len(res.samples):
            cols.append(np.asarray(res.samples_logp)[::thin, None])
            head_s.append("log_prob")
        np.savetxt(p, np.hstack(cols), delimiter=",", header=",".join(head_s), comments="")
        files.append(p)
    p = out / "stage2_history.json"
    p.write_text(json.dumps({"names": res.names, "history": res.history}, indent=1))
    files.append(p)
    if plots:
        files.append(plot_oc(res, out / "stage2_oc.png"))
        if res.rv is not None:
            files.append(plot_rv(res, out / "stage2_rv.png"))
        if res.lc:
            files.append(plot_lc(res, out / "stage2_lc.png"))
    return files


def write_lc_table(res: Stage2Result, path: Path) -> Path:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "source",
                "segment",
                "bjd_tdb",
                "flux",
                "sigma",
                "sigma_quoted",
                "model_start",
                "model",
                "o_minus_model",
                "status",
            ]
        )
        outl = res.lc_outliers or [None] * len(res.lc)
        for q, m0, m1, si, o in zip(
            res.lc, res.lc_start, res.lc_model, res.lc_sigma_in, outl, strict=False
        ):
            for k in range(len(q)):
                w.writerow(
                    [
                        q.source,
                        q.segment[k],
                        f"{q.time[k]:.6f}",
                        f"{q.flux[k]:.6f}",
                        f"{q.sigma[k]:.6f}",
                        f"{si[k]:.6f}",
                        f"{m0[k]:.6f}",
                        f"{m1[k]:.6f}",
                        f"{q.flux[k] - m1[k]:.6f}",
                        "fitted",
                    ]
                )
            if o is not None:
                d, mo, _ = o
                for k in range(len(d)):
                    w.writerow(
                        [
                            d.source,
                            d.segment[k],
                            f"{d.time[k]:.6f}",
                            f"{d.flux[k]:.6f}",
                            f"{d.sigma[k]:.6f}",
                            f"{d.sigma[k]:.6f}",
                            "",
                            f"{mo[k]:.6f}",
                            f"{d.flux[k] - mo[k]:.6f}",
                            "outlier",
                        ]
                    )
    return path


def plot_lc(res: Stage2Result, path: Path) -> Path:
    """Every light-curve segment (at most 16) with the calibrated model of the start values
    and of the fit, and the residuals."""
    import matplotlib

    matplotlib.use("Agg")
    from dataclasses import replace as _replace

    import matplotlib.pyplot as plt

    from .lcfit import apply_calibration

    segs = []
    for qi, q in enumerate(res.lc):
        for s in np.unique(q.segment):
            segs.append((qi, str(s)))
    segs = segs[:16]
    ncol = min(3, len(segs))
    nrow = -(-len(segs) // ncol)
    fig = plt.figure(figsize=(5.0 * ncol, 3.6 * nrow))
    gs = fig.add_gridspec(nrow * 3, ncol)
    # dense model curves of the fit: one accurate model per light curve on the dense grids of
    # all its plotted segments, calibrated with the coefficients of the fit (res.lc_coef)
    fitted = apply(res.raw, res.names, res.x)
    curves: dict = {}
    for qi in sorted({qi for qi, _ in segs}):
        q = res.lc[qi]
        mine = [s for i, s in segs if i == qi]
        grids = [
            np.linspace(q.time[q.segment == s].min(), q.time[q.segment == s].max(), 600)
            for s in mine
        ]
        tt = np.concatenate(grids)
        seg = np.concatenate(
            [np.full(len(g), s, dtype=object) for g, s in zip(grids, mine, strict=True)]
        )
        order = np.argsort(tt, kind="stable")
        dense = _replace(
            q,
            time=tt[order],
            flux=np.zeros(len(tt)),
            sigma=np.ones(len(tt)),
            segment=seg[order].astype(str),
            delta=None,
        )
        try:
            fm = lc_models_for(fitted, res.base_dir, [dense], res.bjd_t0, mode="accurate")[0]
            # the calibration of a segment uses its time centre: that of the measured points
            for s in mine:
                k = dense.segment == s
                ts = q.time[q.segment == s]
                cal = apply_calibration(
                    dense.subset(k), fm[k], res.lc_coef[qi], t_mid=0.5 * (ts.min() + ts.max())
                )
                curves[(qi, s)] = (dense.time[k], cal)
        except (ValueError, RuntimeError, IndexError, KeyError):
            pass
    for k, (qi, s) in enumerate(segs):
        q = res.lc[qi]
        m = q.segment == s
        r, c = divmod(k, ncol)
        ax = fig.add_subplot(gs[3 * r : 3 * r + 2, c])
        axr = fig.add_subplot(gs[3 * r + 2, c], sharex=ax)
        t = q.time[m]
        ax.plot(t - t[0], q.flux[m], ".", ms=2.5, color="#1f5fa8", label="measured")
        ax.plot(t - t[0], res.lc_start[qi][m], "-", lw=0.7, color="0.6", label="start")
        if (qi, s) in curves:
            tt, curve = curves[(qi, s)]
            ax.plot(tt - t[0], curve, "-", lw=0.9, color="#c0392b", label="fitted")
        else:
            ax.plot(t - t[0], res.lc_model[qi][m], "-", lw=0.9, color="#c0392b", label="fitted")
        resid = (q.flux[m] - res.lc_model[qi][m]) * 1e3
        axr.plot(t - t[0], resid, ".", ms=2, color="0.25")
        axr.axhline(0, color="k", lw=0.7)
        ax.set_title(
            f"{q.source}: {s}  (BJD {t[0]:.2f}), rms {np.sqrt(np.mean(resid**2)):.2f} ppt",
            fontsize=8,
        )
        axr.set_xlabel("days")
        if c == 0:
            ax.set_ylabel("flux")
            axr.set_ylabel("O - C [ppt]")
        if k == 0:
            ax.legend(fontsize=6, frameon=False)
    fig.suptitle("swingersim fit, stage 2: light curves (calibrated per segment)", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path


def write_rv_table(res: Stage2Result, path: Path) -> Path:
    rv = res.rv
    sig_in = res.rv_sigma_in if res.rv_sigma_in is not None else rv.sigma
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "star",
                "bjd_tdb",
                "rv",
                "sigma",
                "sigma_quoted",
                "model_start",
                "model",
                "o_minus_model",
                "source",
                "status",
            ]
        )
        for k in range(len(rv)):
            w.writerow(
                [
                    rv.star[k],
                    f"{rv.time[k]:.5f}",
                    f"{rv.rv[k]:.3f}",
                    f"{rv.sigma[k]:.3f}",
                    f"{sig_in[k]:.3f}",
                    f"{res.rv_start[k]:.3f}",
                    f"{res.rv_model[k]:.3f}",
                    f"{rv.rv[k] - res.rv_model[k]:.3f}",
                    rv.source[k],
                    "fitted",
                ]
            )
        o = res.rv_outliers
        if o is not None:
            for k in range(len(o)):
                w.writerow(
                    [
                        o.star[k],
                        f"{o.time[k]:.5f}",
                        f"{o.rv[k]:.3f}",
                        f"{o.sigma[k]:.3f}",
                        f"{o.sigma[k]:.3f}",
                        "",
                        f"{res.rv_model_outliers[k]:.3f}",
                        f"{o.rv[k] - res.rv_model_outliers[k]:.3f}",
                        o.source[k],
                        "outlier",
                    ]
                )
    return path


def plot_rv(res: Stage2Result, path: Path) -> Path:
    """Radial velocities: top, the densest stretch of the data (at most 30 days) with the model
    curves of all stars; bottom, the residuals of all velocities against time."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from .radvel import model_rv

    rv = res.rv
    stars = sorted(set(rv.star.tolist()), key=lambda s: int(s[1:]))
    palette = ["#c0392b", "#e67e22", "#1f5fa8", "#7f8c8d", "#27ae60", "#8e44ad"]
    cols = {s: palette[k % len(palette)] for k, s in enumerate(stars)}
    # densest 30-day window
    t = np.sort(rv.time)
    best = (t[0], 0)
    for a in t:
        nwin = int(np.sum((t >= a) & (t <= a + 30.0)))
        if nwin > best[1]:
            best = (a, nwin)
    w0 = best[0] - 1.0
    w1 = min(best[0] + 31.0, float(t[-1]) + 1.0)
    fitted = apply(res.raw, res.names, res.x)
    cfg = _cfg_of(fitted, res.base_dir)
    shift = rv_shift(fitted, rv)
    gam = float(((fitted.get("fit") or {}).get("rv") or {}).get("gamma", 0.0))
    grid = np.linspace(w0, w1, 1500)
    try:  # all stars in one integration
        all_s, all_t = np.repeat(stars, len(grid)), np.tile(grid, len(stars))
        v = model_rv(cfg, all_s, all_t, res.bjd_t0)
        spec_rv = (fitted.get("fit") or {}).get("rv") or {}
        if spec_rv.get("centre", "mass") == "light":  # the model curve of the light, as fitted
            from .rvlight import light_offset

            v = v + light_offset(cfg, all_s, all_t - res.bjd_t0, spec_rv.get("band", "V"))
        curves = {s: v[k * len(grid) : (k + 1) * len(grid)] + gam for k, s in enumerate(stars)}
    except (ValueError, RuntimeError):
        curves = {}
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(9.0, 6.6), gridspec_kw={"height_ratios": [3, 1.6]})
    for s in stars:
        m = (rv.star == s) & (rv.time >= w0) & (rv.time <= w1)
        if s in curves:
            ax.plot(grid, curves[s], "-", color=cols[s], lw=0.8)
        ax.errorbar(
            rv.time[m],
            rv.rv[m] - (shift[m] - gam),
            yerr=rv.sigma[m],
            fmt="o",
            ms=3.5,
            lw=0.6,
            color=cols[s],
            label=s,
        )
        k = rv.star == s
        axr.errorbar(
            rv.time[k],
            rv.rv[k] - res.rv_model[k],
            yerr=rv.sigma[k],
            fmt="o",
            ms=2.5,
            lw=0.5,
            color=cols[s],
        )
    if res.rv_outliers is not None and len(res.rv_outliers):
        o = res.rv_outliers
        axr.plot(o.time, o.rv - res.rv_model_outliers, "x", color="#b0182c", ms=5, label="left out")
    ax.set_xlim(w0, w1)
    ax.set_xlabel("BJD")
    ax.set_ylabel("radial velocity [km/s]")
    ax.set_title(
        f"densest {w1 - w0:.0f} days of the data with the fitted model (offsets of the "
        "sources removed)",
        fontsize=9,
        pad=22,
    )
    ax.legend(
        fontsize=7, ncol=len(stars), loc="lower center", bbox_to_anchor=(0.5, 1.07), frameon=False
    )
    axr.axhline(0.0, color="k", lw=0.8)
    axr.set_xlabel("BJD")
    axr.set_ylabel("measured - model [km/s]")
    rms = res.rms_rv(res.rv_model)
    axr.set_title("rms " + ", ".join(f"{s} {v:.2f}" for s, v in rms.items()) + " km/s", fontsize=8)
    fig.suptitle("swingersim fit, stage 2: radial velocities", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def _ephemeris(t, kind, P0):
    """Linear ephemeris fitted to the primaries of one pair: epochs, T0, P."""
    pri = kind == 0
    T0 = float(t[pri][0] if np.any(pri) else t[0])
    P = P0
    for _ in range(3):
        E = np.round((t - T0) / P - 0.5 * kind) + 0.5 * kind
        if np.sum(pri) >= 2:
            P, T0 = np.polyfit(E[pri], t[pri], 1)
    E = np.round((t - T0) / P - 0.5 * kind) + 0.5 * kind
    return E, T0, P


def _dense_model(res: Stage2Result, fitted: dict) -> dict:
    """All model eclipses of the fitted configuration over the span of the data (swingersim
    timings), and per pair the mean linear ephemeris of its model primaries: {pair: (times,
    kinds, T0, P)}; {} if it cannot be computed. The ephemeris of the model, not of the data,
    numbers the cycles: data with gaps of years and an apsidal motion of hours cannot be
    counted reliably by themselves."""
    from .timing import eclipse_times

    if res.dense is not None:
        return res.dense
    try:
        cfg = _cfg_of(fitted, res.base_dir)
        y0 = (float(np.min(res.obs.time)) - res.bjd_t0) / DPY - 0.01
        y1 = (float(np.max(res.obs.time)) - res.bjd_t0) / DPY + 0.01
        bins = {n.name for n in _binaries(build_tree_root(cfg))}
        pairs = sorted(p for p in set(res.obs.pair.tolist()) if p in bins)
        ecl = []
        if y0 < 0:
            ecl += eclipse_times(cfg, y0, pairs)[0]
        if y1 > 0:
            ecl += eclipse_times(cfg, y1, pairs)[0]
    except (ValueError, RuntimeError):
        res.dense = {}
        return res.dense
    out = {}
    for p in pairs:
        ee = sorted((e for e in ecl if e.binary == p), key=lambda e: e.t_obs_yr)
        t = np.array([res.bjd_t0 + e.t_obs_yr * DPY for e in ee])
        k = np.array([0 if e.kind == "primary" else 1 for e in ee])
        tp = t[k == 0]
        if len(tp) < 2:
            continue
        T0, P, _ = linear_ephemeris(tp, float(fitted["orbits"][p]["P"]))
        out[p] = (t, k, T0, P)
    res.dense = out
    return out


def plot_oc(res: Stage2Result, path: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pairs = sorted(set(res.obs.pair.tolist()))
    fig, axes = plt.subplots(
        2,
        len(pairs),
        figsize=(5.2 * len(pairs), 6.4),
        squeeze=False,
        gridspec_kw={"height_ratios": [3, 1.3]},
        sharex="col",
    )
    fitted = fitted_config(res)
    dense = _dense_model(res, fitted)
    for k, p in enumerate(pairs):
        m = res.obs.pair == p
        t, kind, sig = res.obs.time[m], res.obs.kind[m], res.obs.sigma[m]
        if p in dense:
            T0, P = dense[p][2], dense[p][3]
            E = np.round((res.t_model[m] - T0) / P - 0.5 * kind) + 0.5 * kind
        else:
            E, T0, P = _ephemeris(t, kind, float(fitted["orbits"][p]["P"]))
        oc = (t - (T0 + P * E)) * 1440.0
        ocm = (res.t_model[m] - (T0 + P * E)) * 1440.0
        oc0 = (res.t_start[m] - (T0 + P * E)) * 1440.0
        yr = 2000.0 + (t - J2000_BJD) / DPY
        ax, axr = axes[0, k], axes[1, k]
        for kk, mk, lab in ((0, "o", "primary"), (1, "s", "secondary")):
            s = kind == kk
            ax.errorbar(
                yr[s],
                oc[s],
                yerr=sig[s] * 1440.0,
                fmt=mk,
                ms=3.5,
                lw=0.6,
                alpha=0.75,
                label=f"measured {lab}",
                color="#1f5fa8" if kk == 0 else "#c2571a",
                mfc="none" if kk else None,
            )
        ax.plot(yr, oc0, ".", color="0.6", ms=2.5, label="model, start values", zorder=1)
        if p in dense:
            tt, kd = dense[p][0], dense[p][1]
            ee = np.round((tt - T0) / P - 0.5 * kd) + 0.5 * kd
            oo = (tt - (T0 + P * ee)) * 1440.0
            for kk, ls in ((0, "-"), (1, ":")):
                s = kd == kk
                ax.plot(
                    2000.0 + (tt[s] - J2000_BJD) / DPY,
                    oo[s],
                    ls,
                    color="k",
                    lw=0.9,
                    label="model, fitted" + (" (secondary)" if kk else ""),
                    zorder=3,
                )
        else:
            order = np.argsort(t)
            ax.plot(yr[order], ocm[order], "-", color="k", lw=1.1, label="model, fitted")
        ax.set_title(
            f"pair {p}: O - C against the mean ephemeris of the model\nT = {T0:.5f} + {P:.7f} E",
            fontsize=9,
        )
        ax.set_ylabel("O - C [min]")
        rr = oc - ocm
        axr.errorbar(yr, rr, yerr=sig * 1440.0, fmt=".", ms=3, lw=0.5, color="0.25")
        if res.outliers is not None and np.any(res.outliers.pair == p):
            mo = res.outliers.pair == p
            to, ko = res.outliers.time[mo], res.outliers.kind[mo]
            Eo = np.round((res.t_outliers[mo] - T0) / P - 0.5 * ko) + 0.5 * ko
            yo = 2000.0 + (to - J2000_BJD) / DPY
            oco = (to - (T0 + P * Eo)) * 1440.0
            ro = (to - res.t_outliers[mo]) * 1440.0
            # outliers far away would squash the panels: they are drawn at the edge (arrows)
            for a, yv, ref in ((ax, oco, np.concatenate([oc, ocm])), (axr, ro, rr)):
                ref = ref[np.isfinite(ref)]
                lo_, hi_ = (float(np.min(ref)), float(np.max(ref))) if len(ref) else (-1.0, 1.0)
                pad = 0.25 * max(hi_ - lo_, 1e-3)
                lo_, hi_ = lo_ - pad, hi_ + pad
                inside = (yv >= lo_) & (yv <= hi_)
                yc = np.clip(yv, lo_, hi_)
                for sel, mk, lab in (
                    (inside, "x", "left out (outlier)"),
                    (~inside & (yc >= hi_), "^", "left out, off scale"),
                    (~inside & (yc <= lo_), "v", "left out, off scale"),
                ):
                    if np.any(sel):
                        a.plot(yo[sel], yc[sel], mk, color="#b0182c", ms=5, mew=1.1, label=lab)
                a.set_ylim(lo_, hi_)
        axr.axhline(0.0, color="k", lw=0.8)
        axr.set_ylabel("measured - fitted [min]")
        axr.set_xlabel("year")
        rms = float(np.sqrt(np.nanmean(rr**2)))
        axr.set_title(f"rms {rms:.2f} min", fontsize=8)
    fig.suptitle("swingersim fit, stage 2: joint N-body model of the eclipse times", fontsize=10)
    handles, labels = [], []
    for a in axes.flat:
        for h, lab in zip(*a.get_legend_handles_labels()):
            if lab not in labels:
                handles.append(h)
                labels.append(lab)
    ncol = min(len(labels), 3 if len(pairs) == 1 else 5)
    nrow = -(-len(labels) // max(ncol, 1))
    bottom = 0.03 * nrow + 0.01  # room for the legend rows
    fig.legend(handles, labels, loc="lower center", ncol=ncol, fontsize=7, frameon=False)
    fig.tight_layout(rect=(0, bottom, 1, 1))
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


# ------------------------------------------------------------------ predictions
def _times_or_nan(x) -> np.ndarray:
    try:
        return _times_for(x)
    except (ValueError, RuntimeError):
        return np.full(len(_CTX["obs"]), np.nan)


@dataclass
class Prediction:
    pair: np.ndarray
    kind: np.ndarray  # 0 primary, 1 secondary
    t: np.ndarray  # BJD_TDB of the fitted model
    sigma: np.ndarray  # [d] half the 16-84 % range of the samples
    lo: np.ndarray  # 16 % [BJD]
    hi: np.ndarray  # 84 % [BJD]
    duration: np.ndarray  # first to last contact [d] (fitted model)
    impact: np.ndarray  # projected separation / (R1 + R2) at minimum
    n_samples: int
    source: str  # "emcee" or "covariance"
    jd_utc: np.ndarray | None = None  # geocentric UTC (with fit.coordinates)


def _reach_au(cfg, node) -> float:
    """How far the stars of the two components of an orbit reach from their centres of mass
    (sum over both components of the largest star distance plus radius) [AU]: a conjunction
    of the orbit with a projected separation below it can give an eclipse."""
    from .backends.nbody import build_simulation
    from .stellar import star_properties

    _, _, _, root, masses = build_simulation(cfg)

    def extent(c) -> float:
        if not hasattr(c, "stars"):
            return star_properties(cfg, f"m{c + 1}", float(masses[c])).radius_rsun * units.RSUN_AU
        el = c.elements
        mt = float(np.sum(masses[c.stars()]))
        best = 0.0
        for side, other in ((c.left, c.right), (c.right, c.left)):
            m_other = float(np.sum(masses[_leaves(other)]))
            best = max(best, el.a * (1.0 + el.e) * m_other / mt + extent(side))
        return best

    nd = {n.name: n for n in root.nodes()}[node.name]
    return extent(nd.left) + extent(nd.right)


def predict_minima(
    res: Stage2Result,
    base_dir,
    years: float = 5.0,
    pairs: list[str] | None = None,
    start: float | None = None,
    samples: int = 100,
    processes: int | None = None,
    seed: int = 7,
    log=print,
) -> Prediction:
    """Eclipses of the fitted model from `start` (default: now, or the last measured minimum
    if it is later) for
    `years`, each with the spread of its time over parameter sets drawn from the emcee
    samples (when the fit ran emcee) or from the covariance matrix: every set is one N-body
    integration that times the same eclipses (the conjunction nearest the central time), in
    parallel. Only eclipses that occur in the fitted model are listed (impact < 1)."""
    from .timing import eclipse_times

    if start is None:  # from now, or from the last minimum if the data reach the future
        from astropy.time import Time

        t_start = max(float(np.max(res.obs.time)), float(Time.now().jd))
    else:
        t_start = float(start)
    t_end = t_start + years * DPY
    if t_end <= res.bjd_t0:
        raise ValueError("predictions: the interval ends before t = 0 of the configuration")
    names = sorted(set(res.obs.pair.tolist())) if not pairs else list(pairs)
    fitted = apply(res.raw, res.names, res.x)
    cfg = Config.model_validate(fitted)
    cfg._base_dir = base_dir
    from .backends.nbody import build_simulation

    root = build_simulation(cfg)[3]
    nodes = {n.name: n for n in root.nodes()}
    unknown = [p for p in names if p not in nodes]
    if unknown:
        raise ValueError(f"predictions: {unknown} are not orbits of the configuration")
    binaries = [p for p in names if not hasattr(nodes[p].left, "stars")]
    binaries = [p for p in binaries if not hasattr(nodes[p].right, "stars")]
    outer = [p for p in names if p not in binaries]
    ecl = eclipse_times(cfg, (t_end - res.bjd_t0) / DPY, binaries)[0] if binaries else []
    ecl = [e for e in ecl if res.bjd_t0 + e.t_obs_yr * DPY >= t_start]
    rows = [
        (
            e.binary,
            0 if e.kind == "primary" else 1,
            res.bjd_t0 + e.t_obs_yr * DPY,
            e.duration_yr * DPY,
            e.impact,
        )
        for e in ecl
    ]
    for p in outer:  # conjunctions of a wider orbit: centres of third-body eclipse events
        P = float(fitted["orbits"][p]["P"])
        k = np.arange(int(np.ceil(2.0 * (t_end - t_start) / P)) + 2)
        seeds = Observed(
            np.array([p] * len(k)),
            k % 2,
            t_start + 0.5 * P * k,
            np.ones(len(k)),
            np.array(["seed"] * len(k)),
        )
        sep = np.full(len(k), np.nan)
        tc = model_times(cfg, seeds, res.bjd_t0, sep)
        reach = _reach_au(cfg, nodes[p])
        for j in range(len(k)):
            if np.isfinite(tc[j]) and t_start <= tc[j] <= t_end:
                rows.append((p, int(k[j] % 2), float(tc[j]), np.nan, float(sep[j] / reach)))
    rows.sort(key=lambda r: r[2])
    if not rows:
        raise ValueError("predictions: no eclipses of the selected pairs in the interval")
    pred = Observed(
        np.array([r[0] for r in rows]),
        np.array([r[1] for r in rows]),
        np.array([r[2] for r in rows]),
        np.ones(len(rows)),
        np.array(["prediction"] * len(rows)),
    )
    ecl_duration = np.array([r[3] for r in rows])
    ecl_impact = np.array([r[4] for r in rows])
    rng = np.random.default_rng(seed)
    if res.samples is not None and len(res.samples) >= 2:
        xs = res.samples[rng.integers(0, len(res.samples), samples)]
        src = "emcee"
    else:
        try:
            xs = rng.multivariate_normal(res.x, res.cov, size=samples, method="svd")
        except (np.linalg.LinAlgError, ValueError):
            xs = np.tile(res.x, (samples, 1))
        src = "covariance"
    ctx = {
        "raw": res.raw,
        "names": res.names,
        "obs": pred,
        "bjd_t0": res.bjd_t0,
        "base_dir": base_dir,
        "priors": {},
    }
    _init(ctx)
    pool, ncore = _pool(0 if processes is None else processes, ctx)
    log(
        f"predictions: {len(rows)} eclipses of {', '.join(names)} in {years:g} yr, "
        f"{samples} parameter sets ({src}) on {ncore} core(s)"
    )
    try:
        if pool is not None:
            T = np.array(list(pool.map(_times_or_nan, list(xs))))
        else:
            T = np.array([_times_or_nan(x) for x in xs])
    finally:
        if pool is not None:
            pool.shutdown()
    # a draw that lost an eclipse to the neighbouring cycle would widen the spread artificially
    period = {p: float(fitted["orbits"][p]["P"]) for p in names}
    half = np.array([0.5 * period[p] for p in pred.pair])
    T[np.abs(T - pred.time) > half] = np.nan
    lo, hi = np.nanpercentile(T, [16, 84], axis=0)
    return Prediction(
        pred.pair,
        pred.kind,
        pred.time,
        0.5 * (hi - lo),
        lo,
        hi,
        ecl_duration,
        ecl_impact,
        len(xs),
        src,
    )


PREDICTION_HEADER = (
    "pair,type,bjd_tdb,sigma_min,bjd_16,bjd_84,duration_h,impact,window_start_bjd,window_end_bjd,"
    "utc,window_start_utc,window_end_utc"
)


def write_predictions(pr: Prediction, path: Path, coordinates=None) -> Path:
    """stage2_predictions.csv; with fit.coordinates also the geocentric UTC of every minimum
    and of its observing window (minimum +- half the duration +- 3 sigma)."""
    from astropy.time import Time

    win = 0.5 * np.nan_to_num(pr.duration, nan=0.0) + 3.0 * np.nan_to_num(pr.sigma, nan=0.0)
    utc = w0 = w1 = None
    if coordinates is not None:
        from .ocfit import from_bjd_tdb

        jd = from_bjd_tdb(pr.t, coordinates.ra_deg, coordinates.dec_deg)
        pr.jd_utc = jd
        iso = lambda v: Time(v, format="jd", scale="utc").iso
        utc, w0, w1 = iso(jd), iso(jd - win), iso(jd + win)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "pair",
                "type",
                "bjd_tdb",
                "sigma_min",
                "bjd_16",
                "bjd_84",
                "duration_h",
                "impact",
                "window_start_bjd",
                "window_end_bjd",
                "utc",
                "window_start_utc",
                "window_end_utc",
            ]
        )
        for j in range(len(pr.t)):
            w.writerow(
                [
                    pr.pair[j],
                    "p" if pr.kind[j] == 0 else "s",
                    f"{pr.t[j]:.6f}",
                    f"{pr.sigma[j] * 1440.0:.3f}",
                    f"{pr.lo[j]:.6f}",
                    f"{pr.hi[j]:.6f}",
                    f"{pr.duration[j] * 24.0:.3f}",
                    f"{pr.impact[j]:.3f}",
                    f"{pr.t[j] - win[j]:.6f}",
                    f"{pr.t[j] + win[j]:.6f}",
                    utc[j][:19] if utc is not None else "",
                    w0[j][:16] if w0 is not None else "",
                    w1[j][:16] if w1 is not None else "",
                ]
            )
    return path


def plot_predictions(res: Stage2Result, pr: Prediction, path: Path) -> Path:
    """O - C of every pair against the mean ephemeris of the fitted model: the measured minima,
    the fitted model, and the predicted minima with their 1-sigma band."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pairs = sorted(set(pr.pair.tolist()))
    fitted = fitted_config(res)
    dense = _dense_model(res, fitted)
    fig, axes = plt.subplots(len(pairs), 1, figsize=(8.0, 2.9 * len(pairs)), squeeze=False)
    for k, p in enumerate(pairs):
        ax = axes[k, 0]
        mp = pr.pair == p
        if p in dense:
            T0, P = dense[p][2], dense[p][3]
        else:
            tp = pr.t[mp & (pr.kind == 0)]
            P = float(fitted["orbits"][p]["P"])
            T0 = float(tp[0]) if len(tp) else float(pr.t[mp][0])
        yr = lambda t: 2000.0 + (np.asarray(t) - J2000_BJD) / DPY

        def oc(t, kind, ref=None, T0=T0, P=P):
            ref = t if ref is None else ref
            E = np.round((ref - T0) / P - 0.5 * kind) + 0.5 * kind
            return (t - (T0 + P * E)) * 1440.0

        mo = res.obs.pair == p
        for kk, mk, col, lab in ((0, "o", "#1f5fa8", "primary"), (1, "s", "#c2571a", "secondary")):
            s = mo & (res.obs.kind == kk)
            ax.errorbar(
                yr(res.obs.time[s]),
                oc(res.obs.time[s], kk, res.t_model[s]),
                yerr=res.obs.sigma[s] * 1440.0,
                fmt=mk,
                ms=3,
                lw=0.5,
                alpha=0.6,
                color=col,
                mfc="none" if kk else None,
                label=f"measured {lab}",
            )
            q = mp & (pr.kind == kk)
            if np.any(q):
                c = oc(pr.t[q], kk)
                ax.fill_between(
                    yr(pr.t[q]),
                    c - pr.sigma[q] * 1440.0,
                    c + pr.sigma[q] * 1440.0,
                    color=col,
                    alpha=0.25,
                    lw=0,
                    label=f"predicted {lab}, 1 sigma",
                )
                ax.plot(yr(pr.t[q]), c, "-" if kk == 0 else ":", color=col, lw=0.9)
        if p in dense:
            tt, kd = dense[p][0], dense[p][1]
            for kk, ls in ((0, "-"), (1, ":")):
                s = kd == kk
                ax.plot(yr(tt[s]), oc(tt[s], kk), ls, color="k", lw=0.7, label=None)
        ax.axvline(yr(float(np.min(pr.t))), color="0.5", lw=0.6, ls="--")
        ax.set_ylabel("O - C [min]")
        ax.set_title(
            f"pair {p}: predicted minima, typical 1 sigma {np.nanmedian(pr.sigma[mp]) * 1440:.2f}"
            f" min, largest {np.nanmax(pr.sigma[mp]) * 1440:.2f} min ({pr.source}, "
            f"{pr.n_samples} sets)",
            fontsize=8,
        )
    axes[-1, 0].set_xlabel("year")
    handles, labels = [], []
    for a in axes.flat:
        for h, lab in zip(*a.get_legend_handles_labels()):
            if lab not in labels:
                handles.append(h)
                labels.append(lab)
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=7, frameon=False)
    fig.tight_layout(rect=(0, 0.07 if len(pairs) == 1 else 0.04, 1, 1))
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


# ------------------------------------------------------------------ anchoring
def data_ephemeris(t, kind, P0: float, t_ref: float):
    """Linear ephemeris (T0 at the minimum nearest t_ref, P) of one pair's measured minima,
    numbered outwards from t_ref: a window of a few cycles around t_ref is fitted first and
    doubled until it covers all data, so that gaps of years are bridged with the period of the
    nearer data (no cycle slips as long as the O - C changes by less than a quarter period
    between a window and the next). Primaries only when there are at least two; returns
    (T0, P, epochs of all minima)."""
    t = np.asarray(t, float)
    kind = np.asarray(kind, int)
    use = kind == 0 if np.sum(kind == 0) >= 2 else np.ones(len(t), bool)
    order = np.argsort(np.abs(t - t_ref))
    T0 = float(t[order[0]] - 0.5 * P0 * kind[order[0]])
    P = float(P0)
    span = float(np.max(np.abs(t - T0))) + P
    w = 10.0 * P
    while True:
        m = use & (np.abs(t - T0) <= w)
        E = np.round((t - T0) / P - 0.5 * kind) + 0.5 * kind
        if np.sum(m) >= 2 and np.ptp(E[m]) > 0:
            P, T0 = (float(v) for v in np.polyfit(E[m], t[m], 1))
        if w >= span:
            break
        w *= 2.0
    E = np.round((t - T0) / P - 0.5 * kind) + 0.5 * kind
    return T0, P, E


def _model_ephemeris(raw: dict, base_dir, pairs: list[str], y0: float, y1: float, bjd_t0: float):
    """Mean linear ephemeris {pair: (T0, P)} of the model's primary eclipses between y0 and y1
    [yr from t = 0], from swingersim timings (one integration for all pairs)."""
    from .timing import eclipse_times

    cfg = _cfg_of(raw, base_dir)
    ecl = []
    if y0 < 0:
        ecl += eclipse_times(cfg, y0, pairs)[0]
    if y1 > 0:
        ecl += eclipse_times(cfg, y1, pairs)[0]
    out = {}
    for p in pairs:
        tp = np.array(
            sorted(bjd_t0 + e.t_obs_yr * DPY for e in ecl if e.binary == p and e.kind == "primary")
        )
        if len(tp) < 2:
            raise ValueError(f"anchoring: the model has fewer than two primary eclipses of {p}")
        T0, P, _ = linear_ephemeris(tp, float(raw["orbits"][p]["P"]))
        out[p] = (T0, P)
    return out


@dataclass
class AnchorResult:
    raw: dict  # the anchored configuration
    pairs: dict  # pair -> {P0, P, M0_0, M0, dt0_min, drift_min, P_data, P_model}
    rms_min: dict  # pair -> rms of measured - model [min] after anchoring
    iterations: int


def anchor(
    raw: dict,
    base_dir,
    obs: Observed,
    bjd_t0: float,
    pairs: list[str] | None = None,
    max_iter: int = 10,
    log=print,
) -> AnchorResult:
    """Adjust the period P and the mean anomaly M0 of every eclipsing pair so that the mean
    linear ephemeris of the model's eclipses equals that of the measured minima: P is scaled by
    P_data / P_model and M0 moves by 360 dt / P, where dt is the delay of the model eclipse at
    the minimum nearest t = 0 (coarse step, robust against drifts of many periods). Then a
    straight line through measured - model at the epochs of the minima of each pair gives the
    remaining phase and period correction (fine step, unbiased by uneven sampling of a curved
    O - C), repeated until no minimum moves by more than 0.1 min. Everything else stays as in
    the configuration. This removes the drift of published
    osculating elements (Lightcurvefactory and other photodynamical codes) in the N-body model
    of swingersim, the first step before any comparison with the minima."""
    d = copy.deepcopy(raw)
    d.pop("fit", None)
    if not pairs:
        bins = {n.name for n in _binaries(build_tree_root(Config.model_validate(d)))}
        names = sorted(p for p in set(obs.pair.tolist()) if p in bins)
        if not names:
            raise ValueError("anchoring: no minima of binaries of two stars")
    else:
        names = list(pairs)
    y0 = (float(np.min(obs.time)) - bjd_t0) / DPY - 0.02
    y1 = (float(np.max(obs.time)) - bjd_t0) / DPY + 0.02
    info = {}
    data = {}
    for p in names:
        o = d["orbits"].get(p)
        if o is None or o.get("P") is None:
            raise ValueError(f"anchoring: orbit {p} needs a period P in the configuration")
        m = obs.pair == p
        if np.sum(m) < 2:
            raise ValueError(f"anchoring: pair {p} has fewer than two minima")
        T0d, Pd, _ = data_ephemeris(obs.time[m], obs.kind[m], float(o["P"]), bjd_t0)
        data[p] = (T0d, Pd)
        info[p] = {"P0": float(o["P"]), "M0_0": float(o.get("M0", 0.0)), "P_data": Pd}
    # 1. coarse: the mean ephemeris of the model eclipses against that of the minima (robust
    #    against drifts of many periods, but biased by the uneven sampling of the data when the
    #    O - C is curved, e.g. by apsidal motion)
    for it in range(1, 4):
        worst = 0.0
        meph = _model_ephemeris(d, base_dir, names, y0, y1, bjd_t0)
        for p in names:
            T0d, Pd = data[p]
            T0m, Pm = meph[p]
            Em = np.round((T0d - T0m) / Pm)
            dt = (T0m + Pm * Em) - T0d  # > 0: the model eclipse comes late
            o = d["orbits"][p]
            o["P"] = float(o["P"]) * Pd / Pm
            o["M0"] = float((float(o.get("M0", 0.0)) + 360.0 * dt / float(o["P"])) % 360.0)
            if it == 1:
                span = float(np.max(np.abs(obs.time[obs.pair == p] - T0d)))
                info[p]["dt0_min"] = dt * 1440.0
                info[p]["drift_min"] = (Pm - Pd) / Pd * span * 1440.0
            info[p]["P_model"] = Pm
            worst = max(worst, abs((Pm - Pd) / Pd * (y1 - y0) * DPY))
            log(f"  anchoring, mean ephemeris {it}: {p} delay {dt * 1440:+.2f} min")
        if worst < 0.05 * min(float(d["orbits"][p]["P"]) for p in names):
            break
    # 2. fine: the model at the epochs of the measured minima; a straight line through
    #    measured - model of every pair (weights 1/sigma^2, outliers beyond 5 robust sigma left
    #    out) gives the phase at t = 0 and the period correction; repeated until both vanish
    it = 0
    best = (np.inf, None)  # (median |O - C| / sigma, configuration) - the fine step can
    # oscillate when the minima scatter by hours or primaries and secondaries disagree (a wrong
    # e cos(omega)); then the best iterate is kept instead of the last one
    worse = 0
    for it in range(1, max_iter + 1):
        tm = _times_at(d, base_dir, obs, bjd_t0)
        cur = _anchor_metric(obs, tm)
        if cur < best[0]:
            best, worse = (cur, copy.deepcopy(d)), 0
        else:
            worse += 1
            if worse >= 3:
                log(
                    "  anchoring: the corrections do not converge (minima scattered by hours, or "
                    "primary and secondary minima disagree - check e cos(omega) of the pair); the "
                    "best iterate is kept"
                )
                break
        done = True
        for p in names:
            m = (obs.pair == p) & np.isfinite(tm)
            o = d["orbits"][p]
            P = float(o["P"])
            r = obs.time[m] - tm[m]
            x = (obs.time[m] - bjd_t0) / P
            z = r / obs.sigma[m]
            mad = 1.4826 * float(np.median(np.abs(z - np.median(z)))) or 1.0
            k = np.abs(z - np.median(z)) <= 5.0 * mad
            if np.sum(k) < 2:
                raise ValueError(f"anchoring: too few usable minima of {p}")
            c1, c0 = np.polyfit(x[k], r[k], 1, w=1.0 / obs.sigma[m][k])
            o["P"] = P + float(c1)
            o["M0"] = float((float(o.get("M0", 0.0)) - 360.0 * float(c0) / P) % 360.0)
            log(f"  anchoring, minima {it}: {p} phase {-c0 * 1440:+.3f} min, period {c1 / P:+.3e}")
            # converged when the correction moves no minimum by more than 0.1 min
            if (abs(c0) + abs(c1) * float(np.max(np.abs(x)))) * 1440.0 > 0.1:
                done = False
        if done:
            break
    tm = _times_at(d, base_dir, obs, bjd_t0)
    if best[1] is not None and _anchor_metric(obs, tm) > 1.2 * best[0]:
        d = best[1]
        tm = _times_at(d, base_dir, obs, bjd_t0)
    rms = {}
    for p in names:
        m = obs.pair == p
        info[p]["P"] = float(d["orbits"][p]["P"])
        info[p]["M0"] = float(d["orbits"][p]["M0"])
        rms[p] = float(np.sqrt(np.nanmean(((obs.time[m] - tm[m]) * 1440.0) ** 2)))
    return AnchorResult(d, info, rms, it)


def _anchor_metric(obs: Observed, tm) -> float:
    """Median |measured - model| / sigma: how well the anchored eclipses match (robust)."""
    z = np.abs(obs.time - np.asarray(tm, float)) / obs.sigma
    z = z[np.isfinite(z)]
    return float(np.median(z)) if len(z) else np.inf


def anchor_report(a: AnchorResult, bjd_t0: float) -> str:
    lines = [
        "swingersim fit --anchor: period and phase of every eclipsing pair anchored to the minima",
        f"  t = 0 of the configuration: BJD {bjd_t0}; {a.iterations} iteration(s)",
        (
            "  pair     P before -> after [d]          M0 before -> after [deg]   delay at t=0   "
            "drift over the data   rms after"
        ),
    ]
    for p, i in a.pairs.items():
        lines.append(
            f"  {p:<8} {i['P0']:.9g} -> {i['P']:.9g} ({i['P'] / i['P0'] - 1:+.2e})   "
            f"{i['M0_0']:.4f} -> {i['M0']:.4f}   {i['dt0_min']:+9.2f} min   "
            f"{i['drift_min']:+10.1f} min   {a.rms_min[p]:8.2f} min"
        )
    lines.append(
        "  (drift: how far the eclipses of the configuration as given would move against the "
        "measured ones over the span of the data)"
    )
    return "\n".join(lines)


def write_anchored(a: AnchorResult, out: Path, source_name: str, bjd_t0: float) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    head = [
        f"# Anchored by SWINGERsim (swingersim fit --anchor) from {source_name}: the period P and the mean anomaly",
        "# M0 of every eclipsing pair adjusted so that the model eclipses follow the minima;",
        "# everything else as in the source. "
        + "; ".join(
            f"{p}: P {i['P0']:.9g} -> {i['P']:.9g}, M0 {i['M0_0']:.4f} -> {i['M0']:.4f}, rms "
            f"{a.rms_min[p]:.2f} min"
            for p, i in a.pairs.items()
        ),
        f"# t = 0 of this configuration = BJD {bjd_t0}.",
    ]
    p = out / "anchored.yaml"
    p.write_text("\n".join(head) + "\n" + yaml.safe_dump(a.raw, sort_keys=False))
    (out / "anchor_report.txt").write_text(anchor_report(a, bjd_t0) + "\n")
    return p


def load_rv(cfg):
    """The measured radial velocities of fit.rv (stars, errors filtered), or None."""
    spec = cfg.fit.rv if cfg.fit is not None else None
    if spec is None:
        return None
    from .radvel import read_rv

    path = Path(spec.file)
    if not path.is_absolute():
        path = Path(cfg.base_dir or ".") / path
    rv = read_rv(path)
    keep = np.isfinite(rv.rv) & np.isfinite(rv.sigma) & (rv.sigma > 0)
    if spec.stars:
        keep &= np.isin(rv.star, spec.stars)
    if spec.max_sigma is not None:
        keep &= rv.sigma <= spec.max_sigma
    rv = rv.subset(keep)
    if not len(rv):
        raise ValueError(f"fit.rv: no velocities left in {path}")
    return rv


# ------------------------------------------------------------------ driver
def hold_one_node(cfg, spec, log=print):
    """A rotation of the whole system about the line of sight changes none of the fitted data
    (minima, velocities, light curves, inclinations): with the node Omega of every orbit free,
    only their differences are measured and the fit is degenerate (the covariance matrix
    singular, errors of 1e5 deg and spurious errors of the other parameters). The node of the
    outermost orbit is then held at its configured value. Returns (spec, note or None)."""
    from .hierarchy import build_tree

    nodes = [f"orbits.{o}.Omega" for o in cfg.orbits]
    if not nodes or not all(n in spec.free for n in nodes):
        return spec, None
    root = build_tree(cfg)[0]
    held = f"orbits.{root.name}.Omega"
    note = (
        f"{held} held at its configured value - with the node of every orbit free only their "
        "differences are measured (a rotation of the whole system on the sky changes none of "
        "the fitted data)"
    )
    log("  note: " + note)
    return spec.model_copy(update={"free": [n for n in spec.free if n != held]}), note


def run_stage2(
    config_path: str | Path, out_dir=None, plots=True, processes=None, log=print, sampler=None
) -> Stage2Result:
    """swingersim fit --stage 2; sampler ("none" / "emcee", swingersim fit --sampler) overrides
    fit.stage2.sampler."""
    from .config import load_config
    from .ocfit import prepare_minima

    config_path = Path(config_path)
    cfg = load_config(config_path)
    if cfg.fit is None or cfg.fit.stage2 is None:
        raise ValueError(f"{config_path}: no 'fit.stage2' section")
    spec = cfg.fit.stage2
    if sampler is not None:
        if sampler not in ("none", "emcee"):
            raise ValueError(f"--sampler {sampler!r}: none or emcee")
        spec = spec.model_copy(update={"sampler": sampler})
    if spec.sampler == "emcee":  # fail now, not after an hour of least squares
        try:
            import emcee  # noqa: F401
        except ImportError:
            raise ValueError(
                "fit.stage2.sampler: emcee needs the package emcee: pip install emcee "
                '(or pip install -e ".[fit]")'
            ) from None
    raw = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    mins = prepare_minima(cfg)
    obs = observed_from_minima(mins, spec.pairs or None, spec.max_sigma_d)
    rv = load_rv(cfg)
    from .lcfit import load as load_lc

    lcs = [load_lc(q, cfg.base_dir) for q in cfg.fit.lc]
    for q in lcs:
        if len(q) and float(q.time.min()) < spec.bjd_t0:
            raise ValueError(
                f"fit.lc {q.source}: points before fit.stage2.bjd_t0 (light curves are computed "
                "forward from t = 0)"
            )
    from .incfit import load as load_inc

    _EXTRA["inc"] = [load_inc(q, cfg.base_dir) for q in cfg.fit.inclination]
    if any(n.startswith("rv.") for n in spec.free) and rv is None:
        raise ValueError("fit.stage2.free: rv.gamma / rv.offset.<source> need fit.rv")
    if rv is not None:
        unknown = {n[10:] for n in spec.free if n.startswith("rv.offset.")} - set(rv.source)
        if unknown:
            raise ValueError(f"fit.stage2.free: no velocities of the source(s) {sorted(unknown)}")
    if len(obs) + (len(rv) if rv is not None else 0) + sum(len(q) for q in _EXTRA["inc"]) <= len(
        spec.free
    ):
        raise ValueError(f"fit.stage2: {len(obs)} minima for {len(spec.free)} free parameters")
    spec, held_note = hold_one_node(cfg, spec, log)
    out = Path(out_dir) if out_dir is not None else Path(cfg.base_dir or ".") / cfg.fit.output.dir
    if spec.anchor:
        log("stage 2: anchoring the eclipsing pairs to the minima first")
        a = anchor(raw, cfg.base_dir, obs, spec.bjd_t0, log=log)
        write_anchored(a, out, config_path.name, spec.bjd_t0)
        log(anchor_report(a, spec.bjd_t0))
        raw = dict(raw, orbits=a.raw["orbits"])
    how = spec.error_scale if spec.error_scale is not None else cfg.fit.error_scale
    res = _fit_with_fidelity(raw, cfg.base_dir, obs, spec, how, processes, log, rv, lcs)
    if held_note:
        res.notes.append(held_note)
    write_outputs(res, out, config_path.name, plots)
    pspec = cfg.fit.output.predict
    if pspec is not None:
        try:
            pr = predict_minima(
                res,
                cfg.base_dir,
                pspec.years,
                pspec.pairs or None,
                pspec.start,
                pspec.samples,
                processes if processes is not None else spec.processes,
                log=log,
            )
        except ValueError as exc:
            # e.g. the fitted pair does not eclipse in the interval (QX Cas): a result, not an error
            msg = f"{exc} (the next {pspec.years:g} yr)"
            if not msg.startswith("predictions"):
                msg = "predictions: " + msg
            log("  " + msg)
            (out / "stage2_predictions.csv").write_text(f"# {msg}\n" + PREDICTION_HEADER + "\n")
            with open(out / "stage2_report.txt", "a", encoding="utf-8") as f:
                f.write(f"  {msg}\n")
        else:
            write_predictions(pr, out / "stage2_predictions.csv", cfg.fit.coordinates)
            if plots:
                plot_predictions(res, pr, out / "stage2_predictions.png")
            res.prediction = pr
    return res


def _set_deltas(pb: _Problem, x, idx: list[int]) -> float:
    """Store accurate - coarse model in every mapped light curve (indices idx of pb.lcs) at the
    parameters x (the coarse model from the cache of the fit); returns the largest change
    against the previous difference [median sigma]."""
    d = apply(pb.ctx["raw"], pb.names, x)
    m = pb.model(x)
    t0, base_dir = pb.ctx["bjd_t0"], pb.ctx["base_dir"]
    acc = lc_models_for(d, base_dir, [pb.lcs[i] for i in idx], t0, mode="accurate")
    change = 0.0
    for i, a in zip(idx, acc, strict=True):
        q = pb.lcs[i]
        if m["lc"] is not None:
            crs = m["lc"][i]
        else:
            crs = lc_models_for(d, base_dir, [q], t0, mode="coarse")[0]
        new = a - crs
        old = q.delta if q.delta is not None else np.zeros(len(q))
        change = max(change, float(np.max(np.abs(new - old)) / np.median(q.sigma)))
        q.delta = new
    return change


FIDELITY_PASSES = 5


def _fit_with_fidelity(raw, base_dir, obs, spec, how, processes, log, rv, lcs):
    """fit_stage2_robust; with light curves of deformed stars (meshes) and fidelity mapped the
    fit runs on the coarse mesh plus the difference accurate - coarse, which is recomputed at
    each solution until it changes by less than 0.3 sigma (at most 5 passes): output space
    mapping, converging to the solution of the accurate model with one accurate evaluation
    per pass. The passes carry the kept points, errors and factors from one to the next and
    sample nothing; emcee (fit.stage2.sampler) runs once at the end with the final difference.
    The start values and their chi^2 are those of the first pass."""
    from .lcfit import lc_config, needs_mesh

    notes = []
    idx = []
    lcs = list(lcs or [])
    for i, q in enumerate(lcs):
        d_lc, nt = lc_config(raw, q)
        notes += [f"{q.source}: {x}" for x in nt]
        if needs_mesh(d_lc) and q.fidelity == "mapped":
            idx.append(i)
    for x in notes:
        log("  shapes " + x)
    if not idx:
        res = fit_stage2_robust(raw, base_dir, obs, spec, how, processes, log, rv=rv, lc=lcs)
        res.shape_notes = notes
        return res
    t_wall = time.perf_counter()
    robust = not (how == "none" and spec.clip_sigma is None)
    if not robust and any(getattr(q, "red_noise", False) for q in lcs):
        log(
            "  warning: fit.lc red_noise needs error_scale or clip_sigma (the fits in a row that "
            "estimate it); it is ignored"
        )
    for i in idx:  # the fit model of these is coarse + delta from now on
        lcs[i].delta = np.zeros(len(lcs[i]))
    pb = _problem(raw, base_dir, obs, spec, rv, lcs)
    _init(pb.ctx)
    x0 = start_values(raw, pb.names)
    change = _set_deltas(pb, x0, idx)
    log(
        f"stage 2, fidelity: coarse mesh + (accurate - coarse), initial difference {change:.2f} sigma"
    )
    passes = []
    note = ""
    pool, ncore = _pool(spec.processes if processes is None else processes, pb.ctx)
    try:
        st = _new_state(pb, spec, x0)
        for k in range(FIDELITY_PASSES):
            f = _run_fits(pb, spec, how, robust, st, pool, ncore, log)
            change = _set_deltas(pb, f["x"], idx)
            passes.append({"pass": k + 1, "chi2_red": f["chi2_red"], "change_sigma": change})
            log(f"stage 2, fidelity pass {k + 1}: difference changed by {change:.3f} sigma")
            if change < 0.3:
                break
        else:
            note = (
                f"the difference accurate - coarse still changed by {change:.2f} sigma after "
                f"{FIDELITY_PASSES} passes (not converged): the solution belongs to the coarse "
                "mesh plus the previous difference"
            )
            log("  warning: " + note)
        # chi^2 with the final difference: the accurate model at the solution
        f["r"] = pb.resid(pb.model(f["x"]), f["x"])
        f["chi2"], f["chi2_red"], f["dof"] = _chi2(pb, f["r"], len(f["x"]))
        cov, err = _errors(pb, f, log)
        mc = _mcmc_info(pb, f, cov, err, spec, pool, ncore, log)
    finally:
        if pool is not None:
            pool.shutdown()
    res = _result(pb, st, f, cov, err, mc, spec, how, t_wall)
    res.fidelity = passes
    res.shape_notes = notes
    res.fidelity_note = (
        "chi^2 and the light-curve models at the solution are those of the accurate meshes; the "
        "linear errors (covariance) come from the Jacobian of the last pass (coarse mesh + "
        "difference)" + ("; emcee sampled with the final difference" if mc is not None else "")
    )
    if note:
        res.note = (res.note + "; " if res.note else "") + note
    return res


def run_anchor(config_path, out_dir=None, bjd_t0: float | None = None, log=print):
    """swingersim fit --anchor: only the anchoring of the eclipsing pairs to the minima of
    fit.minima; writes anchored.yaml and anchor_report.txt to the output folder."""
    from .config import load_config
    from .ocfit import prepare_minima

    config_path = Path(config_path)
    cfg = load_config(config_path)
    if cfg.fit is None:
        raise ValueError(f"{config_path}: no 'fit' section (fit.minima is needed)")
    st = cfg.fit.stage2
    if bjd_t0 is None:
        if st is None:
            raise ValueError("--anchor needs fit.stage2.bjd_t0 or --bjd-t0")
        bjd_t0 = st.bjd_t0
    raw = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    obs = observed_from_minima(
        prepare_minima(cfg),
        (st.pairs or None) if st else None,
        st.max_sigma_d if st else None,
    )
    a = anchor(raw, cfg.base_dir, obs, bjd_t0, log=log)
    out = Path(out_dir) if out_dir is not None else Path(cfg.base_dir or ".") / cfg.fit.output.dir
    write_anchored(a, out, config_path.name, bjd_t0)
    return a, bjd_t0, out
