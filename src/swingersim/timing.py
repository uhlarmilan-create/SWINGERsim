"""Synthetic eclipse timings and O-C diagrams of the inner binaries.

The configuration is integrated with the N-body backend physics (Newton, 1PN, tides as
configured). Every binary whose two components are single stars is followed: its projected
separation on the sky has a minimum at each conjunction, found by sampling the integration
(`samples_per_period` per period of the binary) and refining the minimum with the two-body
motion of the pair between samples. A conjunction is an eclipse when the projected separation
at minimum is smaller than R1 + R2.

Conventions: the line of sight is the z axis and points away from the observer (radial
velocity = v_z); a *primary* eclipse is an eclipse of the binary's first star (m1 of A, m3 of B)
by its companion, i.e. the first star is the farther one. Observed times include the light-time
effect of the binary's centre of mass relative to the barycentre of the whole system,
t_obs = t + z_com / c (Roemer delay; Irwin 1952). The dynamical (physical) delays come from
the integration itself.

O-C: a linear ephemeris T = T0 + P E is fitted (least squares) to the primary eclipses of each
binary; secondary eclipses use the same ephemeris with half-integer epochs.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from . import units
from .config import Config
from .hierarchy import Node
from .kepler import Elements, elements_to_state, state_to_elements

DAY = 1.0 / units.DAYS_PER_YEAR


@dataclass
class Eclipse:
    binary: str
    kind: str  # "primary" or "secondary"
    t_yr: float  # dynamical time of minimum projected separation
    t_obs_yr: float  # including the light-time effect
    impact: float  # projected separation / (R1 + R2)
    duration_yr: float = math.nan  # first to last contact (straight-line chord at minimum)


def _binaries(root: Node) -> list[Node]:
    return [n for n in root.nodes() if isinstance(n.left, int) and isinstance(n.right, int)]


def _advance(mu: float, r: np.ndarray, v: np.ndarray, tau: float):
    """Two-body propagation of a bound relative orbit by tau."""
    el = state_to_elements(mu, r, v)
    n = math.sqrt(mu / el.a**3)
    M = (el.M + n * tau) % (2.0 * math.pi)
    return elements_to_state(mu, Elements(el.a, el.e, el.inc, el.Omega, el.omega, M))


def _g(r, v) -> float:
    """d/dt of half the squared projected separation (x, y on the sky)."""
    return r[0] * v[0] + r[1] * v[1]


def eclipse_times(
    cfg: Config,
    years: float,
    binaries: list[str] | None = None,
    samples_per_period: int = 40,
) -> tuple[list[Eclipse], dict[str, float]]:
    """All eclipses of the selected inner binaries within `years` from t = 0; a negative
    `years` integrates into the past (t from 0 down to `years`), for minima observed before
    the epoch of the configuration's elements.

    Returns the eclipses (in time order) and the nominal period [yr] of each followed binary."""
    from .backends.nbody import build_simulation
    from .stellar import star_properties

    sim, _, _, root, masses = build_simulation(cfg)
    from .pdot import integrator

    advance = integrator(cfg, sim, root, masses)  # orbits.<pair>.Pdot
    sim.collision = "none"
    props = [star_properties(cfg, f"m{k + 1}", m) for k, m in enumerate(masses)]
    nodes = [n for n in _binaries(root) if binaries is None or n.name in binaries]
    if not nodes:
        raise ValueError("no inner binary of two single stars selected")
    p_nom = {
        n.name: units.sma_to_period(n.elements.a, masses[n.left] + masses[n.right]) for n in nodes
    }
    dt = min(p_nom.values()) / samples_per_period
    c = units.C
    ps = sim.particles

    def state(n: Node):
        a, b = ps[n.left], ps[n.right]
        r = np.array([b.x - a.x, b.y - a.y, b.z - a.z])
        v = np.array([b.vx - a.vx, b.vy - a.vy, b.vz - a.vz])
        ma, mb = masses[n.left], masses[n.right]
        zc = (ma * a.z + mb * b.z) / (ma + mb)
        vzc = (ma * a.vz + mb * b.vz) / (ma + mb)
        return r, v, zc, vzc

    prev = {n.name: (0.0, *state(n)) for n in nodes}
    out: list[Eclipse] = []
    t = 0.0
    sign = 1.0 if years >= 0 else -1.0
    while sign * t < sign * years:
        t = t + sign * dt
        t = min(t, years) if sign > 0 else max(t, years)
        advance(t)
        for n in nodes:
            ta, ra, va, zca, vzca = prev[n.name]
            rb, vb, zcb, vzcb = state(n)
            prev[n.name] = (t, rb, vb, zcb, vzcb)
            # the earlier and the later of the two samples (backward runs step into the past)
            if sign > 0:
                t0, r0, v0, zc0, vzc0, t1, r1, v1 = ta, ra, va, zca, vzca, t, rb, vb
            else:
                t0, r0, v0, zc0, vzc0, t1, r1, v1 = t, rb, vb, zcb, vzcb, ta, ra, va
            g0, g1 = _g(r0, v0), _g(r1, v1)
            if not (g0 < 0.0 <= g1):
                continue
            mu = units.G * (masses[n.left] + masses[n.right])
            try:
                tau = brentq(
                    lambda s, mu=mu, r0=r0, v0=v0: _g(*_advance(mu, r0, v0, s)),
                    0.0,
                    t1 - t0,
                    xtol=1e-12,
                )
            except ValueError:
                continue
            r, v = _advance(mu, r0, v0, tau)
            rsum = (props[n.left].radius_rsun + props[n.right].radius_rsun) * units.RSUN_AU
            impact = math.hypot(r[0], r[1]) / rsum
            if impact >= 1.0:
                continue
            v_sky = math.hypot(v[0], v[1])
            dur = 2.0 * rsum * math.sqrt(1.0 - impact**2) / v_sky if v_sky > 0 else math.nan
            # the second star is farther (larger z) -> it is eclipsed: secondary eclipse
            kind = "secondary" if r[2] > 0 else "primary"
            zc = zc0 + vzc0 * tau  # centre of mass moves slowly: linear in tau
            t_ecl = t0 + tau
            ma, mb = masses[n.left], masses[n.right]
            t_obs = t_ecl + zc / c + roemer_in_binary(r[2], ma, mb)
            out.append(Eclipse(n.name, kind, t_ecl, t_obs, impact, dur))
    out.sort(key=lambda e: e.t_yr)
    if getattr(cfg.physics, "activity", None):  # magnetic activity cycles (Applegate 1992)
        from .magnetic import activity_shift_yr

        pair = {n.name: (n.left, n.right) for n in nodes}
        for e in out:
            e.t_obs_yr += float(activity_shift_yr(cfg, pair[e.binary], e.t_yr))
    return out, p_nom


def roemer_in_binary(z_rel: float, m_left: float, m_right: float) -> float:
    """Light-travel time across a binary at an eclipse [yr] (Kaplan 2010, Fabrycky 2010): the
    eclipse is seen when the retarded positions of the two stars line up, not at the barycentric
    conjunction. With z_rel = z_right - z_left (z away from the observer) the minimum is shifted by
    z_rel (m_left - m_right) / (M c) after the light time of the centre of mass: for unequal
    masses the eclipses of the lighter (farther-orbiting) star come later, of the heavier one
    earlier, by up to (a_2 - a_1) / c."""
    return z_rel * (m_left - m_right) / ((m_left + m_right) * units.C)


def linear_ephemeris(times: np.ndarray, p_guess: float) -> tuple[float, float, np.ndarray]:
    """Least-squares T = T0 + P E. The integer epochs are counted from one eclipse to the
    next (robust when the mean period differs from p_guess), then T0 and P are fitted."""
    times = np.sort(np.asarray(times, float))
    # over a gap (eclipses stop while the orbit precesses) the epochs are counted with the
    # mean period of the eclipses before it, not with p_guess: a difference of 1e-3 of the
    # period over a gap of 500 cycles would put the eclipses after it half a period off
    E = np.zeros(len(times))
    p_run = p_guess
    for k in range(1, len(times)):
        E[k] = E[k - 1] + np.round((times[k] - times[k - 1]) / p_run)
        if k >= 2 and E[k] > E[0]:
            p_run = (times[k] - times[0]) / (E[k] - E[0])
    A = np.vstack([np.ones_like(E), E]).T
    (T0, P), *_ = np.linalg.lstsq(A, times, rcond=None)
    return float(T0), float(P), E


def o_minus_c(
    eclipses: list[Eclipse], p_nom: dict[str, float], observed: bool = True
) -> list[dict]:
    """O-C of every eclipse against the linear ephemeris of its binary's primary eclipses."""
    rows = []
    for name, p0 in p_nom.items():
        ecl = [e for e in eclipses if e.binary == name]
        prim = np.array([(e.t_obs_yr if observed else e.t_yr) for e in ecl if e.kind == "primary"])
        if len(prim) < 2:
            continue
        T0, P, _ = linear_ephemeris(prim, p0)
        for e in ecl:
            t = e.t_obs_yr if observed else e.t_yr
            shift = 0.0 if e.kind == "primary" else 0.5
            # nearest epoch of this kind on the fitted ephemeris (O-C stays within P/2)
            E = np.round((t - T0) / P - shift) + shift
            rows.append(
                {
                    "binary": name,
                    "kind": e.kind,
                    "epoch": float(E),
                    "t_yr": e.t_yr,
                    "t_obs_yr": e.t_obs_yr,
                    "o_c_d": (t - (T0 + P * E)) / DAY,
                    "impact": e.impact,
                    "T0_yr": T0,
                    "P_d": P / DAY,
                }
            )
    return rows


def sinusoid_amplitude(t: np.ndarray, oc: np.ndarray, period: float) -> float:
    """Semi-amplitude of a sinusoid of the given period fitted together with a straight line
    (the light-time effect of a circular outer orbit, Irwin 1952)."""
    w = 2.0 * math.pi / period
    A = np.vstack([np.ones_like(t), t, np.sin(w * t), np.cos(w * t)]).T
    coef, *_ = np.linalg.lstsq(A, oc, rcond=None)
    return float(math.hypot(coef[2], coef[3]))


def summarize_timings(rows: list[dict], years: float) -> str:
    lines = [f"eclipse timings over {years:g} yr"]
    for name in sorted({r["binary"] for r in rows}):
        sel = [r for r in rows if r["binary"] == name]
        prim = [r for r in sel if r["kind"] == "primary"]
        sec = [r for r in sel if r["kind"] == "secondary"]
        oc = np.array([r["o_c_d"] for r in prim])
        lines.append(
            f"  {name}: {len(prim)} primary and {len(sec)} secondary eclipses, P = "
            f"{sel[0]['P_d']:.8f} d, O-C of primaries {oc.min() * 1440:+.2f} .. "
            f"{oc.max() * 1440:+.2f} min (semi-amplitude {np.ptp(oc) / 2 * 1440:.2f} min)"
        )
        t = np.array([r["t_yr"] for r in prim])
        if len(t) > 2:
            gaps = np.diff(t)
            big = gaps > 3 * np.median(gaps)
            if big.any():
                lines.append(
                    f"    eclipses missing between {t[:-1][big][0]:.3f} and {t[1:][big][0]:.3f} yr"
                )
            last = t[-1]
            if years - last > 3 * np.median(gaps):
                lines.append(f"    no primary eclipses after t = {last:.3f} yr")
    return "\n".join(lines)


def write_csv(rows: list[dict], path, bjd_t0: float | None = None) -> None:
    """The eclipses as CSV; with `bjd_t0` (BJD of t = 0) also the observed BJD of each one."""
    import csv

    keys = ["binary", "kind", "epoch", "t_yr", "t_obs_yr", "o_c_d", "impact", "T0_yr", "P_d"]
    if bjd_t0 is not None:
        keys.append("bjd_obs")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            row = {k: r[k] for k in keys if k != "bjd_obs"}
            if bjd_t0 is not None:
                row["bjd_obs"] = f"{bjd_t0 + r['t_obs_yr'] * units.DAYS_PER_YEAR:.6f}"
            w.writerow(row)


def plot_oc(rows: list[dict], path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = sorted({r["binary"] for r in rows})
    if not names:  # no eclipses in the time span: an empty figure that says so
        fig, ax = plt.subplots(figsize=(7, 2.6))
        ax.text(0.5, 0.5, "no eclipses in this time span", ha="center", va="center")
        ax.set_axis_off()
        fig.savefig(path, dpi=130)
        plt.close(fig)
        return
    fig, axes = plt.subplots(
        len(names), 1, figsize=(7, 2.6 * len(names)), sharex=True, squeeze=False
    )
    for ax, name in zip(axes[:, 0], names):
        for kind, color in (("primary", "#3b6ea5"), ("secondary", "#d08c2f")):
            sel = [r for r in rows if r["binary"] == name and r["kind"] == kind]
            ax.plot(
                [r["t_obs_yr"] for r in sel],
                [r["o_c_d"] * 1440 for r in sel],
                ".",
                ms=2.5,
                color=color,
                label=kind,
            )
        ax.set_ylabel(f"O-C {name} [min]")
        ax.legend(frameon=False, fontsize=8)
    axes[-1, 0].set_xlabel("t [yr]")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
