"""Stellar evolution: radius, luminosity and envelope structure as functions of mass and age.

Built-in model: the single-star evolution (SSE) formulae of Hurley, Pols & Tout (2000), MNRAS
315, 543, without wind mass loss, for the main sequence (MS), the Hertzsprung gap (HG) and the
first giant branch (GB), 0.1-8 Msun and Z = 0.004-0.03. The envelope quantities (mass and
radius of the convective envelope, gyration radius of the envelope) are those of the BSE code,
Hurley, Tout & Pols (2002), MNRAS 329, 897, Sect. 7.2 (routine ``mrenv``; calibrated for
Z = 0.02 only).

The formulae are evaluated once, offline, with the published SSE Fortran routines
(tools/sse_grid/ssegrid.f) on a grid of 320 masses x 5 metallicities and stored in
``data/sse_tracks.npz``. Every phase is sampled at the same set of 'equivalent evolutionary
points': a fraction x in [0, 1] of the phase (cosine spacing, dense at both ends). A star of
any mass and metallicity is interpolated at fixed (phase, x) between the neighbouring grid
masses (linearly in log M) and metallicities (linearly in log Z), as in METISSE (Agrawal et al.
2020) or the EEP interpolation of MIST (Dotter 2016). The phase boundaries (end of the MS
t_MS, base of the giant branch t_BGB, helium ignition t_HeI) are interpolated in log t.

Beyond helium ignition (core helium burning, AGB, remnants) swingersim does not follow the star:
``StellarTrack.at`` raises ``EvolvedBeyondGrid``; the dynamics stops the run there.

Units: masses Msun, radii Rsun, luminosities Lsun, times Myr, temperatures K.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

DATA = Path(__file__).with_name("data")
TEFF_SUN = 5772.0  # IAU 2015 nominal solar effective temperature [K]
SSE_Z_SUN = 0.02  # solar metallicity of the SSE calibration (Pols et al. 1998)

PHASES = {
    0: "pre-main sequence",
    1: "main sequence",
    2: "Hertzsprung gap",
    3: "giant branch",
    4: "core helium burning",
    5: "early AGB (massive stars: core carbon burning)",
    6: "thermally pulsing AGB",
}  # phase codes shared by the SSE (1-3) and MIST (0-6) tracks


class EvolvedBeyondGrid(ValueError):
    """The star has evolved past the phases covered by the track grid."""


@dataclass(frozen=True)
class StellarState:
    """State of one star at a given age."""

    age_myr: float
    phase: int  # 1 MS, 2 HG, 3 GB (SSE types 0/1, 2, 3)
    sse_type: int  # SSE stellar type kw (0 = deeply convective MS, 1 = MS, 2 = HG, 3 = GB)
    radius_rsun: float
    luminosity_lsun: float
    m_core: float  # helium core mass [Msun] (0 on the MS)
    m_env: float  # mass of the convective envelope [Msun]
    r_env: float  # depth of the convective envelope [Rsun]
    k2_env: float  # squared gyration radius of the envelope (Hurley et al. 2002)

    @property
    def teff_k(self) -> float:
        return TEFF_SUN * self.luminosity_lsun**0.25 / math.sqrt(self.radius_rsun)

    @property
    def rg2(self) -> float:
        """Moment of inertia I/(M R^2) on the MS and HG, where the SSE core is negligible
        (Hurley et al. 2002, eq. 109 with M_c = 0)."""
        return self.k2_env

    @property
    def convective_turnover_yr(self) -> float:
        """Eddy turnover time of the convective envelope, Hurley et al. (2002) eq. (31)
        (after Rasio et al. 1996); 0 without a convective envelope."""
        if self.m_env <= 0.0 or self.r_env <= 0.0:
            return 0.0
        x = self.m_env * self.r_env * (self.radius_rsun - 0.5 * self.r_env)
        return 0.4311 * (x / (3.0 * self.luminosity_lsun)) ** (1.0 / 3.0)


@lru_cache(maxsize=1)
def _grid() -> dict[str, np.ndarray]:
    with np.load(DATA / "sse_tracks.npz") as d:
        return {k: d[k] for k in d.files}


def grid_limits() -> dict[str, tuple[float, float]]:
    g = _grid()
    return {
        "mass": (float(g["mass"][0]), float(g["mass"][-1])),
        "z": (float(g["z"][0]), float(g["z"][-1])),
    }


def _bracket(values: np.ndarray, v: float, log: bool, what: str) -> tuple[int, float]:
    lo, hi = float(values[0]), float(values[-1])
    if not lo * (1 - 1e-9) <= v <= hi * (1 + 1e-9):
        raise ValueError(f"{what} {v} outside the stellar track grid ({lo}-{hi})")
    xs = np.log(values) if log else values
    xv = math.log(v) if log else v
    i = int(min(max(np.searchsorted(xs, xv) - 1, 0), len(xs) - 2))
    w = (xv - xs[i]) / (xs[i + 1] - xs[i])
    return i, float(min(max(w, 0.0), 1.0))


class StellarTrack:
    """Evolution of one star of given initial mass and metallicity (no mass loss)."""

    def __init__(self, mass: float, z: float = SSE_Z_SUN):
        g = _grid()
        self.mass = float(mass)
        self.z = float(z)
        im, wm = _bracket(g["mass"], self.mass, True, "mass")
        iz, wz = _bracket(g["z"], self.z, True, "metallicity Z")
        corners = [
            (iz, im, (1 - wz) * (1 - wm)),
            (iz, im + 1, (1 - wz) * wm),
            (iz + 1, im, wz * (1 - wm)),
            (iz + 1, im + 1, wz * wm),
        ]

        def mix(key, log=False, scale_mass=False):
            out = 0.0
            for a, b, w in corners:
                if w == 0.0:
                    continue
                v = np.asarray(g[key][a, b], dtype=float)
                if scale_mass:
                    v = v / g["mass"][b]
                out = out + w * (np.log(v) if log else v)
            out = np.exp(out) if log else out
            return out * self.mass if scale_mass else out

        self.t_ms = float(mix("t_ms", log=True))
        self.t_bgb = float(mix("t_bgb", log=True))
        self.t_hei = float(mix("t_hei", log=True))
        self._phase = g["phase"]
        self._x = g["x"]
        self._log_l = mix("log_l")
        self._log_r = mix("log_r")
        self._m_core = mix("m_core", scale_mass=True)
        self._m_env = mix("m_env", scale_mass=True)
        # depth of the envelope relative to the radius, then back to Rsun
        rel = 0.0
        for a, b, w in corners:
            if w:
                rel = rel + w * g["r_env"][a, b] / 10.0 ** g["log_r"][a, b]
        self._r_env_rel = rel
        self._k2_env = mix("k2_env")
        # stellar type (kw 0 deeply convective / 1) of the nearest grid track
        a, b, _ = max(corners, key=lambda c: c[2])
        kw = g["kw"][a, b].astype(int)
        self._kw = kw

    # ---------------------------------------------------------------- phases
    def phase_of(self, age_myr: float) -> tuple[int, float]:
        """(phase, fraction x of that phase) at the given age."""
        t = float(age_myr)
        if t < 0:
            raise ValueError("age must not be negative")
        if t <= self.t_ms:
            return 1, t / self.t_ms
        if t <= self.t_bgb:
            return 2, (t - self.t_ms) / (self.t_bgb - self.t_ms)
        if t <= self.t_hei:
            return 3, (t - self.t_bgb) / (self.t_hei - self.t_bgb)
        raise EvolvedBeyondGrid(
            f"star of {self.mass:.3f} Msun at {t:.1f} Myr is past helium ignition "
            f"({self.t_hei:.1f} Myr); swingersim follows stars only to the tip of the giant branch"
        )

    def _phase_sel(self, ph: int):
        cache = self.__dict__.setdefault("_sel_cache", {})
        if ph not in cache:
            sel = np.flatnonzero(self._phase == ph)
            cache[ph] = (sel, self._x[sel])
        return cache[ph]

    def at(self, age_myr: float) -> StellarState:
        ph, x = self.phase_of(age_myr)
        sel, xs = self._phase_sel(ph)

        def f(arr):
            return float(np.interp(x, xs, arr[sel]))

        log_r = f(self._log_r)
        r = 10.0**log_r
        if ph == 1:
            kw = 0 if int(self._kw[sel][min(int(np.searchsorted(xs, x)), len(sel) - 1)]) == 0 else 1
        else:
            kw = ph
        return StellarState(
            age_myr=float(age_myr),
            phase=ph,
            sse_type=kw,
            radius_rsun=r,
            luminosity_lsun=10.0 ** f(self._log_l),
            m_core=max(f(self._m_core), 0.0),
            m_env=max(f(self._m_env), 0.0),
            r_env=max(f(self._r_env_rel), 0.0) * r,
            k2_env=f(self._k2_env),
        )


@lru_cache(maxsize=256)
def track(mass: float, z: float = SSE_Z_SUN) -> StellarTrack:
    return StellarTrack(round(float(mass), 12), round(float(z), 12))
