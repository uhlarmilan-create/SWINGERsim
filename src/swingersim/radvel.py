"""Radial velocities of the stars of a configuration (N-body model) and tables of measured ones.

The radial velocity of star k is its velocity along the line of sight relative to the centre of
mass of the whole system plus the systemic velocity gamma: RV_k = gamma + v_z,k (positive away
from the observer, the convention of the eclipse timings, where t_obs = t + z / c). The measured
velocity at the observed time t_obs belongs to the emission time t_e with t_e + z_k(t_e) / c =
t_obs; the light time of the star (up to tens of minutes in a wide outer orbit) is therefore
included, which matters for the steep RV curves of close binaries (a K of 100 km/s in a 3-day
orbit changes by 0.15 km/s per minute).

Tables of measured velocities (CSV, header lines start with #): columns ``star`` (the star of the
configuration, m1, m2, ...), ``time`` (BJD_TDB [d]), ``rv`` [km/s], ``sigma`` [km/s] and
optionally ``source`` (instrument or paper; per-source zero-point offsets and error corrections
use it).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import units
from .config import Config

DPY = units.DAYS_PER_YEAR
KMS_PER_AUYR = units.AU_KM / units.SECONDS_PER_YEAR  # 1 AU/yr = 4.74 km/s


@dataclass
class ObservedRV:
    star: np.ndarray  # str, m1, m2, ...
    time: np.ndarray  # BJD_TDB [d]
    rv: np.ndarray  # [km/s]
    sigma: np.ndarray  # [km/s]
    source: np.ndarray  # str

    def __len__(self) -> int:
        return len(self.time)

    def subset(self, keep, sigma=None) -> ObservedRV:
        s = self.sigma if sigma is None else sigma
        return ObservedRV(
            self.star[keep], self.time[keep], self.rv[keep], s[keep], self.source[keep]
        )


def read_rv(path: str | Path) -> ObservedRV:
    rows = []
    with open(path, encoding="utf-8") as f:
        lines = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    for r in csv.DictReader(lines):
        r = {k.strip().lower(): (v or "").strip() for k, v in r.items() if k}
        if not r.get("rv") or r["rv"] in ("...", "nan"):
            continue
        rows.append(
            (
                r["star"],
                float(r["time"]),
                float(r["rv"]),
                float(r.get("sigma") or "nan"),
                r.get("source") or "rv",
            )
        )
    if not rows:
        raise ValueError(f"{path}: no radial velocities (columns star, time, rv, sigma, source)")
    star, t, v, s, src = zip(*rows, strict=True)
    return ObservedRV(np.array(star), np.array(t), np.array(v), np.array(s), np.array(src))


def model_rv(cfg: Config, star, time_bjd, bjd_t0: float) -> np.ndarray:
    """Model radial velocity [km/s, without gamma] of the given stars at the observed times
    (BJD_TDB); forward and backward from t = 0 of the configuration as the minima.

    The velocity belongs to the emission time t_e with t_e + z(t_e)/c = t_obs: every point is
    approached from before its emission time (by more than any change of the light time) and
    t_e found by Newton steps with the acceleration of the star until the mismatch is below
    1 ms (obsmodel.RV_FINE_YR); the velocity is then corrected to first order with the acceleration
    (obsmodel.observables, one integration per direction for all points)."""
    from .obsmodel import observables

    return observables(cfg, bjd_t0, rv=(np.asarray(star).astype(str), time_bjd))["rv"]


def rv_curve(cfg: Config, bjd_t0: float, t0: float, t1: float, step_d: float, stars=None):
    """Model radial velocities of all (or the given) stars on a grid of observed times:
    (times [BJD], {star: rv [km/s] without gamma})."""
    n = len(cfg.masses)
    names = list(stars) if stars else [f"m{k + 1}" for k in range(n)]
    t = np.arange(t0, t1 + 0.5 * step_d, step_d)
    v = model_rv(cfg, np.repeat(names, len(t)), np.tile(t, len(names)), bjd_t0)  # one integration
    return t, {s: v[k * len(t) : (k + 1) * len(t)] for k, s in enumerate(names)}


def write_curve(path: Path, t, curves: dict, gamma: float = 0.0) -> Path:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["bjd_tdb", *[f"rv_{s}" for s in curves]])
        for k in range(len(t)):
            w.writerow([f"{t[k]:.6f}", *[f"{curves[s][k] + gamma:.4f}" for s in curves]])
    return Path(path)


def plot_curve(path: Path, t, curves: dict, gamma: float = 0.0, measured=None) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9.0, 4.2))
    cols = ["#c0392b", "#e67e22", "#1f5fa8", "#7f8c8d", "#27ae60", "#8e44ad"]
    for k, (s, v) in enumerate(curves.items()):
        ax.plot(t, v + gamma, "-", lw=0.8, color=cols[k % len(cols)], label=s)
        if measured is not None:
            m = (measured.star == s) & (measured.time >= t[0]) & (measured.time <= t[-1])
            ax.errorbar(
                measured.time[m],
                measured.rv[m],
                yerr=measured.sigma[m],
                fmt="o",
                ms=3,
                lw=0.5,
                color=cols[k % len(cols)],
            )
    ax.set_xlabel("BJD")
    ax.set_ylabel("radial velocity [km/s]")
    ax.set_title(f"swingersim radial velocities (gamma {gamma:g} km/s)", fontsize=9)
    ax.legend(fontsize=7, ncol=len(curves), frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return Path(path)
