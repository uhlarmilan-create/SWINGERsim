"""Phenomenological models of synthetic light curves (F4, phase 1).

swingersim computes light curves from physics (orbits, evolving stars, atmospheres). Observers do
not measure physics, they measure minima: times, depths, widths and shapes of eclipses, which
they describe with a phenomenological model and collect in O-C diagrams. This module reduces a
synthetic light curve in the same way, so that the forward model ends in the quantities that
are fitted to observations in phase 2.

The eclipse model is that of Mikulášek (2015, A&A 584, A8, eq. 13): an eclipse of central depth
A (magnitudes), width D (days), kurtosis Gamma and correction C,

    dm(t) = A (1 + C x^2) {1 - {1 - exp[1 - cosh(x)]}^Gamma},   x = (t - t_min) / D,

plus a quadratic baseline (light of the other stars, a slowly changing level and the
ellipsoidal and reflection variation, which is symmetric about the minimum, to second order). Every eclipse is fitted in two steps, as observers do with a template:

1. a *shape* (A, D, Gamma, C) is fitted for every binary and kind of eclipse (primary,
   secondary) to all its eclipses together, aligned at their geometric times (or per season,
   which gives the slow changes of depth and width, e.g. from the precession of the orbit);
2. the *time of minimum* of every single eclipse is fitted with that shape, together with the
   depth scale and the quadratic baseline; its standard error comes from the covariance of the least
   squares fit with the scatter of the residuals.

The geometric times of conjunction (swingersim timings, with the light-time effect) locate the
eclipses; the phenomenological minima are then independent measurements on the light curve and
include everything that shifts a minimum in real data: asymmetric eclipses of eccentric orbits,
light of the other stars, simultaneous variability, noise and sampling. Their difference to the
geometric times is the error an observer makes by taking a minimum for a conjunction.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

from . import units

DAY = 1.0 / units.DAYS_PER_YEAR
MAX_SHAPE_ECLIPSES = 60  # eclipses (evenly spread) used for one common shape


def profile(dt, depth: float, width: float, gamma: float, corr: float = 0.0):
    """Mikulášek (2015) eq. 13 for one eclipse: dimming [mag] at time dt [d] from minimum."""
    x = np.asarray(dt, float) / width
    core = 1.0 - np.exp(1.0 - np.cosh(np.clip(x, -50.0, 50.0)))
    return depth * (1.0 + corr * x * x) * (1.0 - np.power(np.clip(core, 0.0, 1.0), gamma))


def _profile_jac(dt, depth: float, width: float, gamma: float, corr: float):
    """Profile and its derivatives with respect to (depth, log width, Gamma, C, t_min)."""
    x = np.clip(np.asarray(dt, float) / width, -50.0, 50.0)
    ex = np.exp(1.0 - np.cosh(x))
    core = np.clip(1.0 - ex, 0.0, 1.0)
    pos = core > 0.0
    lc = np.where(pos, np.log(np.where(pos, core, 1.0)), 0.0)
    cg = np.where(pos, np.exp(gamma * lc), 0.0)  # core^Gamma
    poly = 1.0 + corr * x * x
    val = depth * poly * (1.0 - cg)
    # d core^Gamma / dx = Gamma core^(Gamma - 1) e^(1 - cosh x) sinh x
    dcg = np.where(pos, gamma * np.exp((gamma - 1.0) * lc) * ex * np.sinh(x), 0.0)
    dval_dx = depth * (2.0 * corr * x * (1.0 - cg) - poly * dcg)
    return val, np.column_stack(
        [
            poly * (1.0 - cg),  # d/d depth
            -x * dval_dx,  # d/d log width (x = dt / width)
            -depth * poly * cg * lc,  # d/d Gamma
            depth * x * x * (1.0 - cg),  # d/d C
            -dval_dx / width,  # d/d t_min
        ]
    )


@dataclass
class Shape:
    depth: float  # A [mag]
    width: float  # D [d]
    gamma: float
    corr: float
    rms: float  # of the shape fit [mag]
    n_eclipses: int


@dataclass
class Minimum:
    binary: str
    kind: str
    epoch: float
    t_geom_d: float  # geometric conjunction incl. light time [d]
    t_min_d: float  # phenomenological time of minimum [d]
    sigma_d: float
    depth: float  # A of this eclipse [mag]
    rms: float
    o_c_d: float = math.nan  # against the linear ephemeris of the phenomenological minima
    o_c_geom_d: float = math.nan  # the same for the geometric times
    width_d: float = math.nan  # D of this eclipse [d]


@dataclass
class MinimaResult:
    minima: list[Minimum]
    shapes: dict[tuple[str, str], Shape]
    seasons: list[dict] = field(default_factory=list)
    ephemeris: dict[str, tuple[float, float]] = field(default_factory=dict)  # T0 [d], P [d]
    n_blended: int = 0  # eclipses skipped: overlap with an eclipse of another binary
    n_truncated: int = 0  # eclipses skipped: window cut by the start or end of the run
    bjd_t0: float = math.nan  # BJD_TDB of t = 0 of the configuration (time column of the CSV)
    meta: dict = field(default_factory=dict)  # seed, band, cadence, noise: header of the CSV


# BJD_TDB of t = 0 when neither the call nor the configuration (fit.stage2.bjd_t0) gives one:
# J2000.0, so that the times of minima are absolute times on a plain, documented scale
DEFAULT_BJD_T0 = 2451545.0


def _bjd_t0_of(cfg, bjd_t0: float | None) -> float:
    """t = 0 of the run in BJD_TDB: the argument, else fit.stage2.bjd_t0 of the
    configuration, else DEFAULT_BJD_T0."""
    if bjd_t0 is not None:
        return float(bjd_t0)
    stage2 = getattr(getattr(cfg, "fit", None), "stage2", None)
    if stage2 is not None and getattr(stage2, "bjd_t0", None) is not None:
        return float(stage2.bjd_t0)
    return DEFAULT_BJD_T0


def _duration_guess(radii_rsun, p_d: float, a_au: float) -> float:
    """Total duration [d] of a central eclipse of two spheres on a circular orbit."""
    rsum = (radii_rsun[0] + radii_rsun[1]) * units.RSUN_AU
    return p_d / math.pi * min(math.asin(min(rsum / a_au, 1.0)), math.pi / 2)


def _fit_shape_from(
    segments: list[tuple[np.ndarray, np.ndarray]],
    d_guess: float,
    start: tuple,
    tol: float = 1e-10,
    max_nfev: int = 4000,
) -> tuple[np.ndarray, float]:
    """Least-squares fit of a common shape (depth, log width, Gamma, C) and per-eclipse
    quadratic baselines and shifts, starting from `start`; returns (parameters, rms)."""
    n = len(segments)

    sizes = [len(dt) for dt, _ in segments]
    rows = np.concatenate([[0], np.cumsum(sizes)])
    n_par = 4 + 4 * n

    def resid(p):
        depth, lw, gamma, corr = p[:4]
        out = []
        for k, (dt, mag) in enumerate(segments):
            b0, b1, b2, sh = p[4 + 4 * k : 8 + 4 * k]
            base = b0 + dt * (b1 + b2 * dt)
            out.append(mag - (base + profile(dt - sh, depth, math.exp(lw), gamma, corr)))
        return np.concatenate(out)

    def jac(p):
        from scipy.sparse import csr_matrix

        depth, lw, gamma, corr = p[:4]
        data, ri, ci = [], [], []
        for k, (dt, _) in enumerate(segments):
            sh = p[7 + 4 * k]
            _, d = _profile_jac(dt - sh, depth, math.exp(lw), gamma, corr)
            m = len(dt)
            r_idx = np.arange(rows[k], rows[k] + m)
            block = np.column_stack([d[:, :4], np.ones(m), dt, dt * dt, d[:, 4]])
            cols = np.r_[0:4, 4 + 4 * k : 8 + 4 * k]
            data.append(-block.ravel())
            ri.append(np.repeat(r_idx, 8))
            ci.append(np.tile(cols, m))
        return csr_matrix(
            (np.concatenate(data), (np.concatenate(ri), np.concatenate(ci))),
            shape=(rows[-1], n_par),
        )

    p0 = list(start)
    for dt, mag in segments:
        p0 += [float(np.median(mag[np.abs(dt) > 0.8 * np.abs(dt).max()])), 0.0, 0.0, 0.0]
    lo = [0.0, math.log(d_guess / 100.0), 0.05, -0.3] + [-np.inf, -np.inf, -np.inf, -d_guess] * n
    hi = [np.inf, math.log(d_guess * 2.0), 20.0, 5.0] + [np.inf, np.inf, np.inf, d_guess] * n
    p0 = np.clip(p0, np.array(lo) + 1e-12, np.array(hi) - 1e-12)
    r = least_squares(
        resid,
        p0,
        jac=jac,
        bounds=(lo, hi),
        x_scale="jac",
        tr_solver="lsmr",
        max_nfev=max_nfev,
        ftol=tol,
        xtol=tol,
        gtol=tol,
    )
    return r.x, float(np.sqrt(np.mean(r.fun**2)))


def _fit_shape(segments: list[tuple[np.ndarray, np.ndarray]], d_guess: float) -> Shape:
    """Common shape of eclipses aligned at their geometric minima (each with own baseline).

    The shape parameters are correlated (Gamma, C and D all change the wings), so the fit is
    started from several shapes on a few eclipses and the best start is then used for all."""
    n_all = len(segments)
    if n_all > MAX_SHAPE_ECLIPSES:  # a common shape needs no more eclipses than this
        pick = np.linspace(0, n_all - 1, MAX_SHAPE_ECLIPSES).round().astype(int)
        segments = [segments[k] for k in pick]
    depth0 = max(float(np.median([m.max() - m.min() for _, m in segments])), 1e-6)
    probe = segments[:: max(len(segments) // 3, 1)][:3]
    best = None
    for gamma in (0.8, 1.5, 2.5):
        for frac in (4.0, 8.0, 12.0):  # total duration / D
            x, rms = _fit_shape_from(
                probe, d_guess, (depth0, math.log(d_guess / frac), gamma, 0.0), 1e-8, 200
            )
            if best is None or rms < best[1]:
                best = (x, rms)
    x, rms = _fit_shape_from(segments, d_guess, tuple(best[0][:4]))
    return Shape(x[0], math.exp(x[1]), x[2], x[3], rms, n_all)


def _fit_minimum(dt: np.ndarray, mag: np.ndarray, shape: Shape) -> tuple[float, ...]:
    """Time of minimum (offset from dt = 0) of one eclipse with the shape (Gamma, C) of its
    group; depth, width and a quadratic baseline are free (eclipses of precessing orbits change
    both). Returns (offset, sigma, depth, width, rms)."""

    def model(p):
        sh, lscale, lstretch, b0, b1, b2 = p
        val, d = _profile_jac(
            dt - sh,
            shape.depth * math.exp(lscale),
            shape.width * math.exp(lstretch),
            shape.gamma,
            shape.corr,
        )
        return b0 + dt * (b1 + b2 * dt) + val, val, d

    def resid(p):
        return mag - model(p)[0]

    def jac(p):
        _, val, d = model(p)
        return -np.column_stack([d[:, 4], val, d[:, 1], np.ones_like(dt), dt, dt * dt])

    base = float(np.median(mag[np.abs(dt) > 0.8 * np.abs(dt).max()]))
    w = shape.width
    ln2 = math.log(2.0)
    r = least_squares(
        resid,
        [0.0, 0.0, 0.0, base, 0.0, 0.0],
        jac=jac,
        bounds=([-w, -5.0, -ln2, -np.inf, -np.inf, -np.inf], [w, 5.0, ln2, np.inf, np.inf, np.inf]),
        x_scale="jac",
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
    )
    n, k = len(dt), 6
    rms = float(np.sqrt(np.sum(r.fun**2) / max(n - k, 1)))
    try:
        cov = np.linalg.inv(r.jac.T @ r.jac) * rms**2
        sigma = float(math.sqrt(max(cov[0, 0], 0.0)))
    except np.linalg.LinAlgError:
        sigma = math.nan
    depth = shape.depth * math.exp(r.x[1])
    width = shape.width * math.exp(r.x[2])
    return float(r.x[0]), sigma, depth, width, rms


def synthetic_minima(
    cfg,
    years: float,
    binaries: list[str] | None = None,
    band: str = "TESS",
    cadence_min: float = 2.0,
    noise_ppm: float = 0.0,
    seed: int = 1,
    window: float = 2.5,
    season_days: float | None = None,
    samples_per_period: int = 40,
    proximity: bool = True,
    light_time: bool = True,
    isolate: bool = True,
    engine: str = "swingersim",
    phoebe_options: dict | None = None,
    bjd_t0: float | None = None,
) -> MinimaResult:
    """Synthetic minima of the inner binaries over `years`: light curve around every eclipse
    (cadence `cadence_min`, optional white noise), Mikulášek (2015) shape per binary and kind
    (and per season of `season_days` if given), times of minima and O-C.

    isolate=True measures every binary on its own light plus the mean light of the other
    stars (ideal disentangling, as after subtracting a model of the other binary);
    isolate=False uses the total light curve as it is observed: eclipses of other binaries
    in a window are masked or, when they reach the central half of the eclipse, the minimum
    is skipped, and the proximity variation of the other binaries stays in the data.

    engine "phoebe" computes the light curve with PHOEBE for every inner binary
    (swingersim.phoebe_lc; phoebe_options are passed on, e.g. atm, distortion).

    bjd_t0: BJD_TDB of t = 0 of the configuration, for the absolute times of the minima (the
    ``time`` column of write_minima_csv); default fit.stage2.bjd_t0 of the configuration,
    else DEFAULT_BJD_T0 (J2000.0 = BJD_TDB 2451545.0)."""
    from .hierarchy import build_tree
    from .photometry import light_curve
    from .stellar import star_properties
    from .timing import eclipse_times, linear_ephemeris

    # eclipses of all inner binaries (those not followed are needed for masking)
    ecl, p_all = eclipse_times(cfg, years, None, samples_per_period)
    wanted = set(p_all) if binaries is None else set(binaries)
    if not wanted <= set(p_all):
        raise ValueError(f"no inner binary of two stars named {sorted(wanted - set(p_all))}")
    p_nom = {b: p for b, p in p_all.items() if b in wanted}
    if not any(e.binary in wanted for e in ecl):
        raise ValueError("no eclipses of the selected binaries in the run")
    root, masses = build_tree(cfg)
    nodes = {n.name: n for n in root.nodes()}
    radii = [star_properties(cfg, f"m{k + 1}", m).radius_rsun for k, m in enumerate(masses)]
    circ = {}
    for name, p_yr in p_all.items():
        n = nodes[name]
        circ[name] = _duration_guess((radii[n.left], radii[n.right]), p_yr / DAY, n.elements.a)
    # eclipse durations from the sky-plane velocity at minimum (eccentric orbits: primary and
    # secondary differ), the circular estimate if not available
    dur = np.array(
        [e.duration_yr / DAY if math.isfinite(e.duration_yr) else circ[e.binary] for e in ecl]
    )
    step = cadence_min / 1440.0
    centres = np.array([e.t_obs_yr / DAY for e in ecl])
    # the window must contain enough light outside the eclipse to fix the baseline (a linear
    # baseline and a shift of a wide eclipse are otherwise degenerate); at most a quarter period
    p_of = np.array([p_all[e.binary] / DAY for e in ecl])
    halves = np.minimum(window * dur / 2.0, 0.24 * p_of) + step
    t_end = years * units.DAYS_PER_YEAR
    wins, n_blend, n_trunc = [], 0, 0
    for k, e in enumerate(ecl):
        if e.binary not in wanted:
            continue
        tc, h = centres[k], halves[k]
        if tc - h < 0.0 or tc + h > t_end:
            n_trunc += 1  # incomplete window at the ends of the run
            continue
        t = np.arange(tc - h, tc + h + 0.5 * step, step)
        keep = np.ones(len(t), bool)
        blended = False
        for j in np.nonzero(np.abs(centres - tc) < h + 0.5 * dur + step)[0]:
            if j == k or (isolate and ecl[j].binary != e.binary):
                continue
            # another eclipse (other binary, or the other kind of a short-period one) in this
            # window: if it reaches the central half of this eclipse, the minimum is blended
            # and not measured (observers discard such minima too), otherwise its points
            # are masked
            if abs(centres[j] - tc) < 0.25 * dur[k] + 0.5 * dur[j]:
                blended = True
                break
            keep &= np.abs(t - centres[j]) > 0.55 * dur[j] + step
        if blended or keep.sum() < 0.5 * len(t):
            n_blend += 1
            continue
        wins.append((k, tc, t[keep]))
    if not wins:
        raise ValueError("no complete unblended eclipse windows in the run")
    t_all = np.concatenate([w[2] for w in wins])
    order = np.argsort(t_all)
    if engine == "phoebe":
        from .phoebe_lc import phoebe_light_curve

        lc = phoebe_light_curve(
            cfg, t_all[order], band=band, light_time=light_time, **(phoebe_options or {})
        )
    elif engine == "swingersim":
        lc = light_curve(cfg, t_all[order], band=band, proximity=proximity, light_time=light_time)
    else:
        raise ValueError(f"unknown light-curve engine {engine!r} (swingersim or phoebe)")
    flux = np.empty_like(t_all)
    flux[order] = lc.flux
    comps = {}
    for key, val in lc.components.items():
        comps[key] = np.empty_like(t_all)
        comps[key][order] = val
    rng = np.random.default_rng(seed)
    segs, i0 = {}, 0
    for k, tc, t in wins:
        sl = slice(i0, i0 + len(t))
        i0 += len(t)
        e = ecl[k]
        f = flux[sl]
        if isolate:
            # ideal disentangling: the light of the other stars is replaced by its mean over
            # the window (third light), as after subtracting a model of the other binary
            n = nodes[e.binary]
            if e.binary in comps:  # PHOEBE: light of the whole binary
                own = comps[e.binary][sl]
            else:
                own = comps[f"m{n.left + 1}"][sl] + comps[f"m{n.right + 1}"][sl]
            f = own + float(np.mean(f - own))
        if noise_ppm > 0:
            f = f + rng.normal(0.0, noise_ppm * 1e-6, len(f))
        m = -2.5 * np.log10(np.clip(f, 1e-12, None))
        segs.setdefault((e.binary, e.kind), []).append((k, tc, t - tc, m))

    def shape_of(items):
        d_guess = float(np.median([dur[it[0]] for it in items]))
        return _fit_shape([(dt, m) for _, _, dt, m in items], d_guess)

    shapes = {key: shape_of(items) for key, items in segs.items() if items}
    seasons = []
    if season_days:
        for key, items in segs.items():
            t_first = min(tc for _, tc, _, _ in items)
            groups: dict[int, list] = {}
            for it in items:
                groups.setdefault(int((it[1] - t_first) // season_days), []).append(it)
            for g, its in sorted(groups.items()):
                if len(its) < 2:
                    continue
                s = shape_of(its)
                seasons.append(
                    {
                        "binary": key[0],
                        "kind": key[1],
                        "t_mid_d": float(np.mean([it[1] for it in its])),
                        "n": len(its),
                        "depth_mag": s.depth,
                        "width_d": s.width,
                        "gamma": s.gamma,
                        "corr": s.corr,
                        "rms_mag": s.rms,
                    }
                )
    minima = []
    for key, items in segs.items():
        shape = shapes[key]
        for k, tc, dt, m in items:
            off, sig, depth, width, rms = _fit_minimum(dt, m, shape)
            minima.append(
                Minimum(key[0], key[1], math.nan, tc, tc + off, sig, depth, rms, width_d=width)
            )
    if getattr(cfg.physics, "activity", None):  # magnetic activity cycles (Applegate 1992)
        from .magnetic import activity_shift_yr

        pair = {n.name: (n.left, n.right) for n in build_tree(cfg)[0].nodes()}
        for x in minima:
            sh = float(activity_shift_yr(cfg, pair[x.binary], x.t_geom_d * DAY)) / DAY
            x.t_min_d += sh
            x.t_geom_d += sh
    ephem = {}
    for name, p_yr in p_nom.items():
        prim = [x for x in minima if x.binary == name and x.kind == "primary"]
        if len(prim) < 2:
            continue
        p_d = p_yr / DAY
        T0, P, _ = linear_ephemeris(np.array([x.t_min_d for x in prim]), p_d)
        Tg, Pg, _ = linear_ephemeris(np.array([x.t_geom_d for x in prim]), p_d)
        ephem[name] = (T0, P)
        for x in minima:
            if x.binary != name:
                continue
            shift = 0.0 if x.kind == "primary" else 0.5
            E = np.round((x.t_min_d - T0) / P - shift) + shift
            x.epoch = float(E)
            x.o_c_d = x.t_min_d - (T0 + P * E)
            x.o_c_geom_d = x.t_geom_d - (Tg + Pg * E)
    minima.sort(key=lambda x: x.t_min_d)
    from . import __version__

    meta = {
        "band": band,
        "cadence_min": cadence_min,
        "noise_ppm": noise_ppm,
        "seed": seed,
        "engine": engine,
        "isolate": isolate,
        "swingersim_version": __version__,
    }
    return MinimaResult(
        minima,
        shapes,
        seasons,
        ephem,
        n_blended=n_blend,
        n_truncated=n_trunc,
        bjd_t0=_bjd_t0_of(cfg, bjd_t0),
        meta=meta,
    )


def write_minima_csv(
    res: MinimaResult, path, config: str | None = None, bjd_t0: float | None = None
) -> None:
    """CSV of the synthetic minima, readable by ``swingersim fit`` (swingersim.ocfit.read_minima):
    '#' header lines (time_system, bjd_t0, config, SWINGERsim version, seed, band, cadence, noise),
    then the columns binary, kind (aliases of pair, type), epoch, time [BJD_TDB], t_min_d [d
    from t = 0], sigma_d, t_geom_d, o_c_d, o_c_geom_d, depth_mag, width_d, rms_mag.

    bjd_t0: BJD_TDB of t = 0 (default: that of the result, see synthetic_minima); config: the
    name of the configuration file, written to the header."""
    import csv

    t0 = float(bjd_t0) if bjd_t0 is not None else res.bjd_t0
    if not math.isfinite(t0):
        t0 = DEFAULT_BJD_T0
    meta = dict(res.meta)
    head = [
        "# SWINGERsim synthetic minima (swingersim minima, phenomenological fits of the light curve)",
        "# time_system: BJD_TDB",
        f"# bjd_t0: {t0:.6f} (BJD_TDB of t = 0 of the configuration; time = bjd_t0 + t_min_d)",
        f"# config: {config if config is not None else 'unknown'}",
        f"# swingersim_version: {meta.pop('swingersim_version', 'unknown')}",
    ]
    head += [f"# {k}: {v}" for k, v in meta.items()]
    cols = [
        "binary",
        "kind",
        "epoch",
        "time",
        "t_min_d",
        "sigma_d",
        "t_geom_d",
        "o_c_d",
        "o_c_geom_d",
        "depth_mag",
        "width_d",
        "rms_mag",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        f.write("\n".join(head) + "\n")
        w = csv.writer(f)
        w.writerow(cols)
        for x in res.minima:
            w.writerow(
                [
                    x.binary,
                    x.kind,
                    f"{x.epoch:g}",
                    f"{t0 + x.t_min_d:.8f}",
                    f"{x.t_min_d:.8f}",
                    f"{x.sigma_d:.3e}",
                    f"{x.t_geom_d:.8f}",
                    f"{x.o_c_d:.3e}",
                    f"{x.o_c_geom_d:.3e}",
                    f"{x.depth:.6f}",
                    f"{x.width_d:.6f}",
                    f"{x.rms:.2e}",
                ]
            )


def write_seasons_csv(res: MinimaResult, path) -> None:
    import csv

    cols = ["binary", "kind", "t_mid_d", "n", "depth_mag", "width_d", "gamma", "corr", "rms_mag"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for row in res.seasons:
            w.writerow(row)


def plot_minima(res: MinimaResult, path) -> None:
    """O-C of the phenomenological minima (points with error bars) and of the geometric
    conjunctions (line) for every binary and kind; secondary minima shifted by the median
    O-C of their conjunctions (given in the legend)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = sorted(res.ephemeris)
    fig, axes = plt.subplots(len(names), 1, figsize=(10, 3.2 * len(names)), squeeze=False)
    for ax, b in zip(axes[:, 0], names):
        for kind, col in (("primary", "C0"), ("secondary", "C3")):
            m = [x for x in res.minima if x.binary == b and x.kind == kind]
            if not m:
                continue
            t = np.array([x.t_min_d for x in m])
            oc = np.array([x.o_c_d for x in m]) * 86400.0
            ocg = np.array([x.o_c_geom_d for x in m]) * 86400.0
            # displaced secondary minima of eccentric orbits: shown shifted by their median
            shift = float(np.median(ocg)) if kind == "secondary" else 0.0
            label = f"{kind} minima" + (f" ({-shift:+.0f} s)" if abs(shift) >= 1.0 else "")
            ax.errorbar(
                t,
                oc - shift,
                yerr=np.array([x.sigma_d for x in m]) * 86400.0,
                fmt="o",
                ms=2.5,
                color=col,
                label=label,
            )
            ax.plot(t, ocg - shift, "-", color=col, lw=0.8)
        ax.set_ylabel(f"O-C {b} [s]")
        ax.legend(fontsize=8)
    axes[-1, 0].set_xlabel("time [d]")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def summarize(res: MinimaResult) -> str:
    lines = []
    for (b, k), s in sorted(res.shapes.items()):
        mins = [x for x in res.minima if x.binary == b and x.kind == k]
        d = np.array([x.t_min_d - x.t_geom_d for x in mins]) * 86400.0
        sig = np.array([x.sigma_d for x in mins]) * 86400.0
        lines.append(
            f"  {b} {k:9s}: {len(mins)} minima; shape A = {s.depth:.4f} mag, D = "
            f"{s.width * 24:.3f} h, Gamma = {s.gamma:.3f}, C = {s.corr:+.3f}, rms {s.rms:.1e} mag; "
            f"minimum - conjunction {np.median(d):+.1f} s (scatter {np.std(d):.1f} s), "
            f"median sigma {np.median(sig):.1f} s"
        )
    if res.n_blended or res.n_truncated:
        lines.append(f"  skipped: {res.n_blended} blended, {res.n_truncated} truncated eclipses")
    for b, (T0, P) in res.ephemeris.items():
        lines.append(f"  {b}: ephemeris of the minima T0 = {T0:.6f} d, P = {P:.8f} d")
    return "\n".join(lines)
