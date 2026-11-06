"""MIST v1.2 stellar tracks (MESA Isochrones and Stellar Tracks; Choi et al. 2016, ApJ 823, 102;
Dotter 2016, ApJS 222, 8) as an alternative to the analytic SSE tracks.

The grids are built from the published EEP archives with ``tools/mist/build_grid.py`` and stored
in ``data/mist/feh_<feh>_vvcrit<rot>.npz`` (one file per [Fe/H] and initial rotation
v/v_crit = 0 or 0.4). Each track is given at the equivalent evolutionary points (EEPs) 1-1409:
pre-main sequence from the birthline-like start of the MIST models (EEP 1-201), main sequence
(202-453), subgiant and red giant branch (454-604), core helium burning (605-707), the early
AGB (708-807) and the thermally pulsing AGB to its departure (808-1409; massive stars: to the end
of carbon burning).
A star of any mass and [Fe/H] is interpolated at fixed EEP between the neighbouring masses
(linearly in log M) and metallicities (linearly in [Fe/H]) – the method of the MIST isochrones
(Dotter 2016) – and in age between EEPs.

Ages are counted from the start of the MIST pre-main-sequence models (MIST convention), not
from the zero-age main sequence as in the SSE tracks. Masses are the initial masses, 0.1-300 Msun.
The MIST wind mass loss is recorded (``current_mass``, ``mdot``); the
dynamics follows it with physics.mass_loss: {source: track}. The rotating grids (v/v_crit = 0.4)
lose much more mass above ~50 Msun than the non-rotating ones (rotationally boosted winds and
rotational mixing; Choi et al. 2016, eq. 21): a 100 Msun star leaves the main sequence with
~25 Msun (rotating) instead of ~60 Msun.
The MIST tracks contain no convective-envelope structure; for the tidal lag from the envelope
(physics.tidal_lag: hurley2002) the envelope is estimated from the track (swingersim.stellar).
"""

from __future__ import annotations

import math
from functools import lru_cache
from pathlib import Path

import numpy as np

from .evolution import PHASES, EvolvedBeyondGrid, StellarState

MIST_DIR = Path(__file__).with_name("data") / "mist"
EEP_ZAMS, EEP_TAMS, EEP_RGB_TIP, EEP_ZACHEB, EEP_TACHEB = 202, 454, 605, 631, 707
EEP_TPAGB, EEP_POSTAGB = 808, 1409  # first thermal pulse (massive: end of C burning); end
MASS_COLS = ("star_mass", "he_core_mass", "star_mdot")
MIST_PHASES = PHASES  # 2 = subgiant (Hertzsprung gap), first half of the TAMS - RGB tip span


def available() -> list[tuple[float, float]]:
    """[(feh, vvcrit)] of the installed MIST grids."""
    out = []
    for p in sorted(MIST_DIR.glob("feh_*_vvcrit*.npz")):
        s = p.stem
        feh = float(s[4:9].replace("p", "+").replace("m", "-"))
        out.append((feh, float(s.split("vvcrit")[1])))
    return sorted(out)


@lru_cache(maxsize=32)
def _grid(feh_tag: str, vvcrit: float) -> dict:
    path = MIST_DIR / f"feh_{feh_tag}_vvcrit{vvcrit:.1f}.npz"
    if not path.exists():
        raise FileNotFoundError(
            f"MIST grid {path.name} is not installed; build it with tools/mist/build_grid.py "
            "from the MIST v1.2 EEP archive (https://mist.science)"
        )
    with np.load(path) as d:
        return {k: d[k] for k in d.files}


def _tag(feh: float) -> str:
    return ("p" if feh >= 0 else "m") + f"{abs(feh):.2f}"


HAYASHI_TEFF = 5500.0  # [K] a giant-branch star is cooler
HAYASHI_SLOPE = 0.25  # |d log Teff / d log L| on the giant branch is below (Hertzsprung gap: L
# about constant while Teff falls)


def giant_branch_start(log_teff, log_l, window: int = 10, half: int = 5) -> int | None:
    """Index of the base of the giant branch in a track section from the end of the main
    sequence to the helium ignition: the first point from which, for `window` points, the star
    is on the Hayashi line - cooler than HAYASHI_TEFF and climbing it (L rising, Teff only slowly
    falling) - (Hurley et al. 2000, Sect. 5.2: the Hertzsprung gap ends where the star reaches the
    Hayashi track). None: the star ignites helium before (massive stars cross the gap then)."""
    lt, ll = np.asarray(log_teff, float), np.asarray(log_l, float)
    n = len(lt)
    bad = np.ones(n, bool)
    for i in range(n):
        i0, i1 = max(i - half, 0), min(i + half, n - 1)
        dl = ll[i1] - ll[i0]
        if dl > 0:
            bad[i] = abs((lt[i1] - lt[i0]) / dl) > HAYASHI_SLOPE
    bad |= lt > np.log10(HAYASHI_TEFF)
    for i in range(n - 2):
        if not bad[i : i + window].any():
            return i
    return None


def _phase(eep: float, eep_gb: float | None = None) -> int:
    """Phase of a (1-based, fractional) EEP: the primary EEPs of Dotter (2016, Table 2) begin
    their phases (ZAMS = EEP 202 is the first main-sequence point, as t_ms = age[EEP_TAMS - 1]
    is the end of the main sequence); eep_gb: the base of the giant branch of the track
    (MistTrack.eep_gb; default the middle of the TAMS - RGB tip interval)."""
    if eep < EEP_ZAMS:
        return 0
    if eep < EEP_TAMS:
        return 1
    if eep < (0.5 * (EEP_TAMS + EEP_RGB_TIP) if eep_gb is None else eep_gb):
        return 2  # Hertzsprung gap / subgiant
    if eep < EEP_RGB_TIP:
        return 3
    if eep < EEP_TACHEB:
        return 4
    if eep < EEP_TPAGB:
        return 5  # early AGB (massive stars: core carbon burning, red supergiant)
    return 6


class MistTrack:
    """Evolution of one star of initial mass ``mass`` and metallicity ``feh`` from MIST."""

    def __init__(self, mass: float, feh: float = 0.0, vvcrit: float = 0.4):
        grids = [g for g in available() if abs(g[1] - vvcrit) < 1e-6]
        if not grids:
            raise FileNotFoundError(f"no MIST grid with v/vcrit = {vvcrit} installed in {MIST_DIR}")
        fehs = np.array([g[0] for g in grids])
        from .gridclip import MIST_FEH_TOL, clip_to_grid

        # gridclip policy: the nearest installed grid up to MIST_FEH_TOL dex beyond the range
        feh = clip_to_grid(
            feh, float(fehs[0]), float(fehs[-1]), "[Fe/H]", "installed MIST grids", MIST_FEH_TOL
        )
        self.mass, self.feh, self.vvcrit = float(mass), float(feh), float(vvcrit)
        iz = int(np.clip(np.searchsorted(fehs, feh) - 1, 0, max(len(fehs) - 2, 0)))
        if len(fehs) == 1:
            pairs = [(fehs[0], 1.0)]
        else:
            wz = float(np.clip((feh - fehs[iz]) / (fehs[iz + 1] - fehs[iz]), 0.0, 1.0))
            pairs = [(fehs[iz], 1.0 - wz), (fehs[iz + 1], wz)]
        cols = None
        acc = None
        for fz, wzz in pairs:
            if wzz == 0.0:
                continue
            g = _grid(_tag(fz), vvcrit)
            cols = list(g["columns"])
            masses = g["mass"]
            if not masses[0] <= mass <= masses[-1]:
                raise ValueError(f"mass {mass} outside the MIST grid ({masses[0]}-{masses[-1]})")
            lm = np.log(masses)
            im = int(np.clip(np.searchsorted(lm, math.log(mass)) - 1, 0, len(lm) - 2))
            wm = float((math.log(mass) - lm[im]) / (lm[im + 1] - lm[im]))
            for b, wb in ((im, 1.0 - wm), (im + 1, wm)):
                w = wzz * wb
                if w == 0.0:
                    continue
                t = g["eep"][b].astype(float)
                if str(g.get("age_mode", "")) == "steps":  # increments, summed in float64
                    t[:, 0] = np.cumsum(t[:, 0])
                t = np.column_stack([np.log10(t[:, 0]), t[:, 1:]])  # log age
                # masses and mass rates are interpolated as fractions of the initial mass
                # (exact at the grid masses, no bias between them; as for the SSE tracks)
                for c in MASS_COLS:
                    if c in cols:
                        t[:, cols.index(c)] /= masses[b]
                acc = w * t if acc is None else acc + w * t  # NaN where any track ends
        for c in MASS_COLS:
            if c in cols:
                acc[:, cols.index(c)] *= self.mass
        self._cols = {c: i for i, c in enumerate(cols)}
        ok = np.flatnonzero(np.all(np.isfinite(acc), axis=1))
        n = int(ok[-1]) + 1 if len(ok) else 0
        if n < 2:
            raise ValueError("MIST track has no valid points")
        self._t = acc[:n]
        self._age = 10.0 ** self._t[:, 0] / 1e6  # Myr
        # enforce monotonic age (interpolation between masses keeps it, check anyway)
        self._age = np.maximum.accumulate(self._age)
        self.t_start = float(self._age[0])
        self.t_end = float(self._age[-1])
        self.n_eep = n
        # EEP axis and contiguous columns for np.interp (the strided self._t[:, i] is copied by
        # every call; same values, so the interpolation is bit-identical)
        self._x = np.arange(1, n + 1)
        self._c = {c: np.ascontiguousarray(self._t[:, i]) for c, i in self._cols.items()}
        self.t_ms = float(self._age[EEP_TAMS - 1]) if n >= EEP_TAMS else self.t_end
        self.t_hei = float(self._age[EEP_RGB_TIP - 1]) if n >= EEP_RGB_TIP else self.t_end
        # base of the giant branch (EEP); = EEP_RGB_TIP: no giant branch (massive stars)
        self.eep_gb = float(EEP_RGB_TIP)
        if n >= EEP_RGB_TIP:
            sl = slice(EEP_TAMS - 1, EEP_RGB_TIP - 1)
            g = giant_branch_start(self._c["log_Teff"][sl], self._c["log_L"][sl])
            if g is not None:
                self.eep_gb = float(EEP_TAMS + g)
        elif n > EEP_TAMS:
            self.eep_gb = float(min(0.5 * (EEP_TAMS + EEP_RGB_TIP), n))

    def eep_of(self, age_myr: float) -> float:
        t = float(age_myr)
        if t < self.t_start:
            raise ValueError(
                f"age {t:.3f} Myr precedes the start of the MIST pre-main-sequence model "
                f"({self.t_start:.3f} Myr)"
            )
        if t > self.t_end:
            raise EvolvedBeyondGrid(
                f"star of {self.mass:.3f} Msun at {t:.1f} Myr is beyond the end of its MIST track "
                f"({self.t_end:.1f} Myr, end of core helium burning or of the model age range)"
            )
        return float(np.interp(t, self._age, self._x))

    def at(self, age_myr: float) -> StellarState:
        return self.at_eep(self.eep_of(age_myr), age_myr)

    def at_eep(self, eep: float, age_myr: float | None = None) -> StellarState:
        """State at an equivalent evolutionary point."""
        eep = float(min(max(eep, 1.0), self.n_eep))
        if age_myr is None:
            age_myr = float(np.interp(eep, self._x, self._age))
        x, cols = self._x, self._c

        def f(col):
            return float(np.interp(eep, x, cols[col]))

        ph = _phase(eep, self.eep_gb)
        return StellarState(
            age_myr=float(age_myr),
            phase=ph,
            sse_type=-1,
            radius_rsun=10.0 ** f("log_R"),
            luminosity_lsun=10.0 ** f("log_L"),
            m_core=max(f("he_core_mass"), 0.0),
            m_env=0.0,
            r_env=0.0,
            k2_env=0.0,
        )

    def _col(self, age_myr: float, col: str) -> float:
        eep = self.eep_of(age_myr)
        return float(np.interp(eep, self._x, self._c[col]))

    def current_mass(self, age_myr: float) -> float:
        return self._col(age_myr, "star_mass")

    def mdot(self, age_myr: float) -> float:
        """Wind mass-loss rate of the MIST model [Msun/yr] (negative = loss)."""
        return self._col(age_myr, "star_mdot")

    def surface_h1(self, age_myr: float) -> float:
        """Surface hydrogen mass fraction."""
        return self._col(age_myr, "surface_h1")

    def is_wolf_rayet(self, age_myr: float) -> bool:
        """Wolf-Rayet star in the sense of the MIST wind scheme: X_surf < 0.4 and
        Teff > 1e4 K (Choi et al. 2016, sect. 3.6)."""
        st = self.at(age_myr)
        return self.surface_h1(age_myr) < 0.4 and st.teff_k > 1.0e4


@lru_cache(maxsize=2048)  # massloss._initial_mass scans ~400 masses
def mist_track(mass: float, feh: float = 0.0, vvcrit: float = 0.4) -> MistTrack:
    return MistTrack(round(float(mass), 12), round(float(feh), 12), round(float(vvcrit), 3))


def evolutionary_state(
    radius_rsun: float,
    teff_k: float,
    feh: float = 0.0,
    vvcrit: float = 0.0,
    m_range: tuple[float, float] = (2.0, 150.0),
    n_mass: int = 120,
) -> dict:
    """Single-star MIST star with the given radius and effective temperature: initial mass,
    age and present mass, on the first passage of the radius (main sequence or the start of
    the Hertzsprung gap). The 'evolutionary mass' of a star placed in the HR diagram; compared
    with a dynamical mass it shows the mass discrepancy of products of binary interaction."""
    rows = []
    for mi in np.geomspace(m_range[0], m_range[1], n_mass):
        tr = mist_track(float(mi), feh, vvcrit)
        ages = np.linspace(tr._age[0], tr.t_end * (1.0 - 1e-9), 1500)
        eeps = np.clip(np.interp(ages, tr._age, tr._x), 1.0, tr.n_eep)  # as MistTrack.at
        r = 10.0 ** np.interp(eeps, tr._x, tr._c["log_R"])
        # first upward crossing after the ZAMS minimum (the pre-main sequence contracts)
        k0 = int(np.argmin(r[: max(len(r) // 3, 2)]))
        idx = np.nonzero((r[k0:-1] < radius_rsun) & (r[k0 + 1 :] >= radius_rsun))[0]
        if len(idx) == 0:
            continue
        k = k0 + idx[0]
        a = ages[k] + (radius_rsun - r[k]) / (r[k + 1] - r[k]) * (ages[k + 1] - ages[k])
        rows.append((mi, a, tr.at(a).teff_k, tr.current_mass(a)))
    if len(rows) < 2:
        raise ValueError(f"no MIST track reaches R = {radius_rsun} Rsun in {m_range} Msun")
    rows = np.array(rows)
    d = rows[:, 2] - teff_k
    j = np.nonzero(d[:-1] * d[1:] <= 0.0)[0]
    if len(j) == 0:
        raise ValueError(f"no MIST track has R = {radius_rsun} Rsun at Teff = {teff_k} K")
    j = int(j[0])
    w = -d[j] / (d[j + 1] - d[j])
    m_init, age, _, m_now = rows[j] + (rows[j + 1] - rows[j]) * w
    return {"m_init": float(m_init), "age_myr": float(age), "m_now": float(m_now)}
