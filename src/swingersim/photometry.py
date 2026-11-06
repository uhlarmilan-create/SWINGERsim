"""Synthetic photometry of multiple stars (F4 b, tier A): light curves from the N-body positions.

Model
-----
* Stars are spheres with a radially symmetric limb darkening, the four-coefficient law of
  Claret (2000), I(mu)/I(1) = 1 - sum_k a_k (1 - mu^(k/2)), coefficients from the tables of
  Claret & Bloemen (2011, A&A 529, A75, LTE ATLAS/PHOENIX) interpolated in (Teff, log g,
  [M/H]) below 15 000 K and from the non-LTE TLUSTY atmospheres of Reeve & Howarth (2016,
  MNRAS 456, 1294; OSTAR2002 and BSTAR2006, 15 000-55 000 K) above; any law with fewer terms
  (linear, quadratic in sqrt(mu)) is a special case. Coefficients can be given per star.
  Gravity darkening y from Claret & Bloemen; beyond their grid (O stars above ~39 000 K) von
  Zeipel, T^4 ~ g, in the band y = (1/4) x e^x / (e^x - 1), x = h c / (lambda_eff k Teff).
  Values beyond the edges of the tables: swingersim.gridclip (clipped, with a warning).
* Surface brightness in a passband: black body, photon-counting detector,
  w = integral B_lambda(Teff) S(lambda) lambda d lambda (passband curves in data/passbands:
  Bessell UBVRI, TESS, Gaia G/BP/RP). The flux of a star outside eclipse is pi R^2 w.
* Eclipses: every star k is occulted by all stars nearer to the observer whose disks overlap it.
  The blocked flux is integral_0^R I(r) alpha(r) 2 pi r dr, where alpha(r) is the fraction of
  the circle of radius r around the centre of k covered by the union of the occulting disks;
  the radial integral is done with Gauss-Legendre quadrature between the radii at which a disk
  edge crosses (as in ellc, Maxted 2016). Any number of simultaneous occulters is exact for
  spheres, so mutual eclipses of the two binaries of a 2+2 system and triple overlaps need no
  special treatment.
* Light-travel time: each star is seen at its own emission time t_obs - z/c (z away from the
  observer, barycentre at z = 0), to first order: x(t_obs) - v z / c.
* Proximity effects of the stars of two-star binaries (small-distortion approximations):
  ellipsoidal variation -alpha_ell q (R/d)^3 (2 cos^2 psi - sin^2 i) (Morris 1985; Morris &
  Naftilan 1993) with alpha_ell = 0.15 (15 + u)(1 + tau)/(3 - u); Doppler beaming
  -(3 - alpha_nu) v_z / c (Loeb & Gaudi 2003); reflection of the companion's light by a Lambert
  sphere of Bond albedo A (geometric albedo 2A/3), A = 0.5 for convective and 1 for radiative
  envelopes (Rucinski 1969).

Validity: detached systems (the stars well inside their Roche lobes). The validation compares
the eclipses with batman (Kreidberg 2015) and ellc (Maxted 2016); tier B (ellc with Roche
geometry) and tier C (PHOEBE) are for close binaries.

Units: times in days in the user interface (yr internally), fluxes normalised to the total
flux of all stars outside eclipse.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import lru_cache
from itertools import pairwise
from pathlib import Path

import numpy as np

from . import units

DATA = Path(__file__).with_name("data")
H_PLANCK = 6.62607015e-34
C_LIGHT = 2.99792458e8
K_BOLTZ = 1.380649e-23
LD_PROXY = {"TESS": "I", "Gaia_G": "R", "Gaia_RP": "I", "Gaia_BP": "V"}
GL_NODES = 24


# --------------------------------------------------------------------------- passbands


def bands() -> list[str]:
    return sorted(p.stem for p in (DATA / "passbands").glob("*.csv"))


@lru_cache(maxsize=32)
def passband(name: str) -> tuple[np.ndarray, np.ndarray]:
    """(wavelength [m], response) of a passband."""
    path = DATA / "passbands" / f"{name}.csv"
    if not path.exists():
        raise ValueError(f"unknown passband {name}; available: {', '.join(bands())}")
    d = np.loadtxt(path, delimiter=",", comments="#")
    return d[:, 0] * 1e-10, d[:, 1]


def _planck(lam: np.ndarray, teff: float) -> np.ndarray:
    x = H_PLANCK * C_LIGHT / (lam * K_BOLTZ * teff)
    with np.errstate(over="ignore"):
        return 2 * H_PLANCK * C_LIGHT**2 / lam**5 / np.expm1(x)


@lru_cache(maxsize=4096)
def surface_brightness(teff: float, band: str) -> float:
    """Photon-counting passband intensity of a black body (arbitrary but common units)."""
    lam, s = passband(band)
    return float(np.trapezoid(_planck(lam, teff) * s * lam, lam))


def beaming_alpha(teff: float, band: str) -> float:
    """alpha_nu = d ln F_nu / d ln nu of a black body averaged over the passband (photon
    weights); the beaming factor is 3 - alpha_nu (Loeb & Gaudi 2003)."""
    lam, s = passband(band)
    nu = C_LIGHT / lam
    x = H_PLANCK * nu / (K_BOLTZ * teff)
    alpha = 3.0 - x * np.exp(x) / np.expm1(x)
    w = _planck(lam, teff) * s * lam
    return float(np.trapezoid(alpha * w, lam) / np.trapezoid(w, lam))


# --------------------------------------------------------------------------- limb darkening


@lru_cache(maxsize=1)
def _claret2011():
    with np.load(DATA / "claret2011_ld.npz") as d:
        return {k: d[k] for k in d.files}


CLARET2011_GRID = "Claret & Bloemen (2011) limb-darkening grid"


def _claret_interp(teff: float, logg: float, m_h: float, b: str) -> list[float] | None:
    """[y, a1..a4] of Claret & Bloemen (2011) interpolated at (Teff, log g, [M/H]) (log g and
    [M/H] already clipped to the grid), None where the grid has no models."""
    from scipy.interpolate import RegularGridInterpolator

    g = _claret2011()
    grid = g["coeff"][list(g["bands"]).index(b)]  # (5, m_h, logg, log_teff): y, a1..a4
    pts = (g["m_h"], g["logg"], g["log_teff"])
    x = np.array([[m_h, logg, math.log10(teff)]])
    vals = [
        float(RegularGridInterpolator(pts, grid[c], bounds_error=False)(x)[0]) for c in range(5)
    ]
    return None if any(math.isnan(v) for v in vals) else vals


def _claret_clip(logg: float, m_h: float) -> tuple[float, float]:
    """log g and [M/H] clipped to the Claret & Bloemen (2011) grid (gridclip policy)."""
    from .gridclip import clip_to_grid

    g = _claret2011()
    m_h = clip_to_grid(m_h, float(g["m_h"][0]), float(g["m_h"][-1]), "[M/H]", CLARET2011_GRID)
    logg = clip_to_grid(logg, float(g["logg"][0]), float(g["logg"][-1]), "log g", CLARET2011_GRID)
    return logg, m_h


def _claret_band(band: str) -> str:
    b = LD_PROXY.get(band, band)
    if b not in list(_claret2011()["bands"]):
        raise ValueError(f"no limb-darkening table for band {band}")
    return b


def _claret_teff_max(g, grid, iz: int, logg: float) -> float:
    """Hottest model of the Claret & Bloemen grid at the log g node nearest to logg."""
    j = int(np.argmin(np.abs(g["logg"] - logg)))
    ok = np.isfinite(grid[1, iz, j])
    return float(10 ** g["log_teff"][ok].max()) if ok.any() else math.nan


def claret_coefficients(
    teff: float, logg: float, m_h: float, band: str
) -> tuple[np.ndarray, float]:
    """Four-coefficient limb darkening (a1..a4) and gravity-darkening exponent y of Claret &
    Bloemen (2011) for a band (TESS, Gaia: nearest Johnson-Cousins band, see LD_PROXY).

    Outside the grid (swingersim.gridclip policy): log g and [M/H] are clipped to the grid with a
    warning; where the grid has no model at (Teff, log g) - stars hotter than the hottest model
    atmosphere of their log g - the nearest log g with a model (0.5 or 1 dex away) is used, and
    beyond that (or below 3500 K, e.g. AGB stars) the nearest node; each with a warning that
    names the cause."""
    g = _claret2011()
    b = _claret_band(band)
    logg, m_h = _claret_clip(logg, m_h)
    vals = _claret_interp(teff, logg, m_h, b)
    if vals is not None:
        return np.array(vals[1:]), vals[0]
    grid = g["coeff"][list(g["bands"]).index(b)]
    iz = int(np.argmin(np.abs(g["m_h"] - m_h)))
    t_lo, t_hi = 10 ** float(g["log_teff"][0]), 10 ** float(g["log_teff"][-1])
    if teff < t_lo:
        cause = f"Teff {teff:.0f} K is below the coolest model ({t_lo:.0f} K)"
    elif teff > t_hi:
        cause = f"Teff {teff:.0f} K is above the hottest model ({t_hi:.0f} K)"
    else:
        cause = (
            f"Teff {teff:.0f} K is hotter than the hottest model atmosphere at log g "
            f"{logg:.2f} (~{_claret_teff_max(g, grid, iz, logg):.0f} K)"
        )
    if t_lo <= teff <= t_hi:
        for dg in (0.5, -0.5, 1.0, -1.0):
            v2 = _claret_interp(teff, logg + dg, m_h, b)
            if v2 is not None:
                _warn_ld(
                    (round(teff, -1), round(logg, 1), "dg"),
                    f"{cause} in the {CLARET2011_GRID}; using log g {logg + dg:.2f}",
                )
                return np.array(v2[1:]), v2[0]
    vals = _nearest_ld_node(g, grid, iz, logg, teff, cause)
    if vals is None:
        raise ValueError(
            f"no limb darkening for Teff {teff:.0f} K, log g {logg:.2f} in band {band}"
        )
    return np.array(vals[1:]), vals[0]


_LD_WARNED: set = set()


def _warn_ld(key, msg: str) -> None:
    import warnings

    if key not in _LD_WARNED:
        _LD_WARNED.add(key)
        warnings.warn(msg, stacklevel=4)


def _nearest_ld_node(g, grid, iz: int, logg: float, teff: float, cause: str) -> list[float] | None:
    """Coefficients of the nearest tabulated (Teff, log g) node at the [M/H] node iz for stars
    outside the model-atmosphere grid, e.g. AGB stars (Teff ~ 3000 K, log g ~ 0); a warning
    naming the cause is issued once per star."""
    lt, lg = np.meshgrid(g["log_teff"], g["logg"])
    ok = np.isfinite(grid[1, iz])
    if not ok.any():
        return None
    d = ((lt - math.log10(teff)) / 0.05) ** 2 + ((lg - logg) / 0.5) ** 2
    d[~ok] = np.inf
    j, k = np.unravel_index(int(np.argmin(d)), d.shape)
    _warn_ld(
        (round(teff, -1), round(logg, 1)),
        f"{cause} in the {CLARET2011_GRID}; using the nearest node (Teff "
        f"{10 ** g['log_teff'][k]:.0f} K, log g {g['logg'][j]:.1f})",
    )
    return [float(grid[c, iz, j, k]) for c in range(5)]


TLUSTY_TEFF_MIN = 15000.0  # hotter stars use the TLUSTY tables (Reeve & Howarth 2016)


@lru_cache(maxsize=1)
def _tlusty():
    with np.load(DATA / "tlusty_ld.npz") as d:
        return {k: d[k] for k in d.files}


@lru_cache(maxsize=64)
def _tlusty_interp(m_h: float, band: str):
    from scipy.interpolate import LinearNDInterpolator

    g = _tlusty()
    sel = np.isclose(g["m_h"], m_h)
    b = list(g["bands"]).index(band)
    pts = np.column_stack([g["teff"][sel] / 1e4, g["logg"][sel]])
    return LinearNDInterpolator(pts, g["coeff"][sel][:, b, :]), pts


def _lambda_eff(band: str) -> float:
    lam, s = passband(band)
    return float(np.trapezoid(lam * s, lam) / np.trapezoid(s, lam))


def von_zeipel_y(teff: float, band: str) -> float:
    """Gravity-darkening exponent of a radiative star in a band: T^4 ~ g (von Zeipel 1924)
    and black-body intensity, y = d ln I / d ln g = (1/4) x e^x / (e^x - 1)."""
    x = H_PLANCK * C_LIGHT / (_lambda_eff(band) * K_BOLTZ * teff)
    return 0.25 * x / (-math.expm1(-x))


def tlusty_coefficients(teff: float, logg: float, m_h: float, band: str) -> np.ndarray:
    """Four-coefficient limb darkening (a1..a4) of the TLUSTY non-LTE atmospheres (Reeve &
    Howarth 2016), interpolated linearly in (Teff, log g) at the two nearest metallicities
    (C, G, L, S, T = 2, 1, 1/2, 1/5, 1/10 of solar) and linearly in [M/H] between them. A point
    just outside the model grid (up to 0.25 in log g or 1500 K) takes the nearest model; farther
    out it is an error."""
    g = _tlusty()
    b = LD_PROXY.get(band, band)
    if b not in g["bands"]:
        raise ValueError(f"no TLUSTY limb-darkening table for band {band}")
    from .gridclip import clip_to_grid

    mhs = np.unique(g["m_h"])
    mh = clip_to_grid(
        m_h, float(mhs[0]), float(mhs[-1]), "[M/H]", "TLUSTY limb-darkening grid (OSTAR/BSTAR)"
    )
    k = int(np.clip(np.searchsorted(mhs, mh) - 1, 0, len(mhs) - 2))
    w = (mh - mhs[k]) / (mhs[k + 1] - mhs[k])
    out = np.zeros(4)
    for z, wz in ((mhs[k], 1 - w), (mhs[k + 1], w)):
        if wz == 0.0:
            continue
        f, pts = _tlusty_interp(float(z), b)
        x = np.array([teff / 1e4, logg])
        v = f(x[None, :])[0]
        if np.any(np.isnan(v)):
            d = np.hypot((pts[:, 0] - x[0]) / 0.15, (pts[:, 1] - x[1]) / 0.25)
            j = int(np.argmin(d))
            if d[j] > 1.0:
                raise ValueError(
                    f"no TLUSTY limb darkening for Teff {teff:.0f} K, log g {logg:.2f} "
                    f"([M/H] = {z:+.2f}; grid 15 000-55 000 K, log g limited by the Eddington limit)"
                )
            f2 = g["coeff"][np.isclose(g["m_h"], z)][:, list(g["bands"]).index(b), :]
            v = f2[j]
        out += wz * v
    return out


def limb_darkening(teff: float, logg: float, m_h: float, band: str) -> tuple[np.ndarray, float]:
    """Limb darkening (a1..a4) and gravity-darkening exponent y of a star: Claret & Bloemen
    (2011) below 15 000 K, the TLUSTY non-LTE tables (Reeve & Howarth 2016) above; y from Claret
    & Bloemen where their grid has a model at the star's Teff and log g, von Zeipel beyond (the
    hottest O stars). Grid edges: swingersim.gridclip."""
    if teff >= TLUSTY_TEFF_MIN:
        a = tlusty_coefficients(teff, logg, m_h, band)
        # gravity darkening of Claret & Bloemen where their grid has a model at this Teff and
        # log g, von Zeipel beyond it (O stars hotter than the hottest model of their log g)
        lg, mh = _claret_clip(logg, m_h)
        vals = _claret_interp(teff, lg, mh, _claret_band(band))
        return a, (vals[0] if vals is not None else von_zeipel_y(teff, band))
    return claret_coefficients(teff, logg, m_h, band)


def ld_profile(a: np.ndarray, mu: np.ndarray) -> np.ndarray:
    """I(mu)/I(1) of the four-coefficient law (a may have fewer than four terms)."""
    out = np.ones_like(mu)
    for k, ak in enumerate(a, start=1):
        out -= ak * (1.0 - mu ** (k / 2.0))
    return out


def ld_mean(a: np.ndarray) -> float:
    """Disk average of I(mu)/I(1): 1 - sum a_k k / (k + 4)."""
    return 1.0 - sum(ak * k / (k + 4.0) for k, ak in enumerate(a, start=1))


# --------------------------------------------------------------------------- occultation


_GL_X, _GL_W = np.polynomial.legendre.leggauss(GL_NODES)


def covered_fraction(r: np.ndarray, occ: list[tuple[float, float, float]]) -> np.ndarray:
    """Fraction of the circle of radius r (around the origin) inside the union of the disks
    occ = [(x, y, R), ...]."""
    r = np.atleast_1d(r)
    out = np.zeros_like(r)
    arcs_all = []
    for x, y, rad in occ:
        d = math.hypot(x, y)
        phi = math.atan2(y, x)
        with np.errstate(divide="ignore", invalid="ignore"):
            c = (r * r + d * d - rad * rad) / (2.0 * r * d)
        half = np.where(c <= -1.0, np.pi, np.where(c >= 1.0, 0.0, np.arccos(np.clip(c, -1, 1))))
        if d == 0.0:
            half = np.where(r <= rad, np.pi, 0.0)
        arcs_all.append((phi, half))
    if len(arcs_all) == 1:
        return arcs_all[0][1] / np.pi
    for i in range(len(r)):
        segs = []
        for phi, half in arcs_all:
            h = half[i]
            if h <= 0.0:
                continue
            if h >= np.pi:
                segs = [(0.0, 2 * np.pi)]
                break
            lo = (phi - h) % (2 * np.pi)
            hi = lo + 2 * h
            if hi > 2 * np.pi:
                segs += [(lo, 2 * np.pi), (0.0, hi - 2 * np.pi)]
            else:
                segs.append((lo, hi))
        if not segs:
            continue
        segs.sort()
        total, cur_lo, cur_hi = 0.0, segs[0][0], segs[0][1]
        for lo, hi in segs[1:]:
            if lo > cur_hi:
                total += cur_hi - cur_lo
                cur_lo, cur_hi = lo, hi
            else:
                cur_hi = max(cur_hi, hi)
        total += cur_hi - cur_lo
        out[i] = min(total / (2 * np.pi), 1.0)
    return out


def blocked_fraction(radius: float, a: np.ndarray, occ: list[tuple[float, float, float]]) -> float:
    """Fraction of the flux of a limb-darkened disk of radius `radius` (centre at the origin)
    hidden by the occulting disks occ = [(x, y, R), ...] (same length units)."""
    breaks = {0.0, 1.0}
    for x, y, rad in occ:
        d = math.hypot(x, y)
        for b in (abs(d - rad) / radius, (d + rad) / radius):
            if 0.0 < b < 1.0:
                breaks.add(b)
    edges = sorted(breaks)
    total = 0.0
    for s0, s1 in pairwise(edges):
        # s = s0 + (s1 - s0)(1 - cos th)/2 clusters the nodes at both ends of the interval, where
        # the covered arc (disk edges) and mu (limb) have square-root singularities
        th = 0.5 * np.pi * (_GL_X + 1.0)
        s = s0 + 0.5 * (s1 - s0) * (1.0 - np.cos(th))
        ds = 0.5 * (s1 - s0) * np.sin(th) * 0.5 * np.pi
        mu = np.sqrt(np.clip(1.0 - s * s, 0.0, 1.0))
        f = ld_profile(a, mu) * covered_fraction(s * radius, occ) * 2.0 * s * ds
        total += float(np.dot(_GL_W, f))
    return total / ld_mean(a)  # normalised by the integral over the whole disk


# --------------------------------------------------------------------------- light curve


@dataclass
class StarLight:
    name: str
    radius_au: float
    weight: float  # passband intensity (disk average)
    ld: np.ndarray
    gd_y: float
    albedo: float
    beaming: float  # 3 - alpha_nu


@dataclass
class LightCurve:
    t_day: np.ndarray
    flux: np.ndarray  # normalised: 1 = all stars outside eclipse
    band: str
    components: dict[str, np.ndarray] = field(default_factory=dict)  # per-star flux
    stars: list[StarLight] = field(default_factory=list)
    pulsation: dict[str, np.ndarray] = field(default_factory=dict)  # flux factor of pulsators
    meta: dict = field(default_factory=dict)  # how the curve was computed (write_csv header)

    def header_lines(self) -> list[str]:
        """The '#' comment lines of write_csv: band, time system, options, limb darkening,
        SWINGERsim version, noise."""
        from . import __version__

        m = dict(self.meta)
        lines = [f"SWINGERsim {__version__} light curve ({m.pop('engine', 'swingersim')})"]
        lines.append(f"band: {self.band}")
        bjd = m.pop("bjd_t0", None)
        ts = m.pop("time_system", "BJD_TDB")
        lines.append(
            "time: t_day = days after t = 0 of the configuration (observer's time; "
            "light-travel time counted from the plane through the barycentre of the system); "
            + (f"{ts} = {bjd:.6f} + t_day" if bjd is not None else f"{ts} of t = 0: not given")
        )
        for key, label in (
            ("light_time", "light-travel time"),
            ("proximity", "proximity effects (ellipsoidal, reflection, beaming)"),
        ):
            if key in m:
                lines.append(f"{label}: {'on' if m.pop(key) else 'off'}")
        if "exposure_min" in m:
            ex, ns = m.pop("exposure_min"), m.pop("n_sub", 1)
            lines.append(
                f"exposure: {ex:g} min, {ns} sub-exposures"
                if ex > 0 and ns > 1
                else "exposure: instantaneous (no integration)"
            )
        lines.append(
            "normalisation: "
            + m.pop("normalisation", "1 = summed flux of all stars out of eclipse")
        )
        law = m.pop("ld_law", None)
        if law:
            lines.append(f"limb darkening: {law}")
        for s in self.stars:
            lines.append(
                f"  {s.name}: R = {s.radius_au / units.RSUN_AU:.5f} Rsun, a_k = "
                + ", ".join(f"{x:.5f}" for x in np.atleast_1d(s.ld))
                + f", gravity darkening y = {s.gd_y:.4f}"
            )
        if "noise_ppm" in m:
            lines.append(
                f"noise: {m.pop('noise_ppm'):g} ppm (Gaussian), seed {m.pop('noise_seed')}"
            )
        else:
            lines.append("noise: none")
        lines += [f"{k}: {v}" for k, v in m.items()]
        return lines

    def write_csv(self, path, header: bool = True) -> None:
        """CSV with columns t_day, flux, flux_<component>; ``header``: '#' comment lines with
        the metadata (header_lines) before the column names."""
        cols = ["t_day", "flux"] + [f"flux_{k}" for k in self.components]
        data = np.column_stack([self.t_day, self.flux, *self.components.values()])
        meta = "".join(f"# {ln}\n" for ln in self.header_lines()) if header else ""
        np.savetxt(
            path, data, delimiter=",", header=meta + ",".join(cols), comments="", fmt="%.10g"
        )


def star_lights(cfg, masses, props, band: str, ld: dict | None = None) -> list[StarLight]:
    from .stellar import LOGG_SUN

    out = []
    sm = cfg.stellar_model
    m_h = {"sse": math.log10(sm.z / 0.02), "mist": sm.feh}.get(sm.kind, 0.0)
    for k, (m, p) in enumerate(zip(masses, props, strict=True)):
        name = f"m{k + 1}"
        if p.teff_k is None:
            raise ValueError(
                f"star {name}: no effective temperature (use stellar_model zams or sse)"
            )
        logg = LOGG_SUN + math.log10(m) - 2.0 * math.log10(p.radius_rsun)
        a, y = limb_darkening(p.teff_k, logg, m_h, band)
        if ld and name in ld:  # given limb darkening; the gravity darkening of the tables stays
            a = np.asarray(ld[name], float)
        conv = (p.m_env or 0.0) > 1e-3 * m
        out.append(
            StarLight(
                name=name,
                radius_au=p.radius_rsun * units.RSUN_AU,
                weight=surface_brightness(round(p.teff_k, 1), band),
                ld=a,
                gd_y=y,
                albedo=0.5 if conv else 1.0,
                beaming=3.0 - beaming_alpha(p.teff_k, band),
            )
        )
    return out


def proximity_terms(
    pos: np.ndarray,
    vel: np.ndarray,
    masses: np.ndarray,
    lights: list[StarLight],
    pairs: list[tuple[int, int]],
    base: np.ndarray,
    ellipsoidal: set | None = None,
    reflection: set | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Proximity effects of the stars of two-star binaries, per star: (relative ellipsoidal
    change of the flux, reflected flux [units of base], Doppler-beaming factor).

    ``ellipsoidal`` / ``reflection``: the pairs (i, j) that get the term (None = all); the mesh
    photometry leaves out the ellipsoidal term of Roche bodies (their shape carries it) and the
    reflection inside a contact envelope."""
    n = len(lights)
    ell = np.zeros(n)
    refl = np.zeros(n)
    for i, j in pairs:
        do_ell = ellipsoidal is None or (i, j) in ellipsoidal
        do_refl = reflection is None or (i, j) in reflection
        if not (do_ell or do_refl):
            continue
        r = pos[j] - pos[i]
        d = float(np.linalg.norm(r))
        cos_psi = r[2] / d  # direction to the companion against the line of sight
        if do_ell:
            h = np.cross(r, vel[j] - vel[i])
            hn = float(np.linalg.norm(h))
            sin2i = 1.0 - (h[2] / hn) ** 2 if hn > 0 else 1.0
        for s, o, sign in ((i, j, 1.0), (j, i, -1.0)):
            L = lights[s]
            if do_ell:
                u_lin = 3.0 * (1.0 - ld_mean(L.ld))
                tau = L.gd_y  # gravity-darkening coefficient of the intensity, I ~ g^y
                alpha_ell = 0.15 * (15.0 + u_lin) * (1.0 + tau) / (3.0 - u_lin)
                q = masses[o] / masses[s]
                ell[s] -= alpha_ell * q * (L.radius_au / d) ** 3 * (2.0 * cos_psi**2 - sin2i)
            if do_refl:
                # reflection: star s lit by o; phase angle between o->s and s->observer
                cos_beta = -sign * cos_psi  # +z points away from the observer
                beta = math.acos(max(-1.0, min(1.0, cos_beta)))
                phase = (math.sin(beta) + (math.pi - beta) * math.cos(beta)) / math.pi
                refl[s] += (2.0 / 3.0) * L.albedo * base[o] * (L.radius_au / d) ** 2 * phase
    beam = np.array([1.0 - lights[s].beaming * vel[s][2] / units.C for s in range(n)])
    return ell, refl, beam


def visible_fractions(pos: np.ndarray, lights: list[StarLight]) -> np.ndarray:
    """Visible fraction of each star's disk (spheres of their radii, limb-darkened), the stars
    nearer to the observer in front."""
    vis = np.ones(len(lights))
    order = np.argsort(pos[:, 2])  # nearest to the observer first
    for idx_b, k in enumerate(order):
        occ = []
        for j in order[:idx_b]:
            dx, dy = pos[j, 0] - pos[k, 0], pos[j, 1] - pos[k, 1]
            if math.hypot(dx, dy) < lights[j].radius_au + lights[k].radius_au:
                occ.append((dx, dy, lights[j].radius_au))
        if occ:
            vis[k] = 1.0 - blocked_fraction(lights[k].radius_au, lights[k].ld, occ)
    return vis


def flux_at(
    pos: np.ndarray,
    vel: np.ndarray,
    masses: np.ndarray,
    lights: list[StarLight],
    pairs: list[tuple[int, int]],
    proximity: bool = True,
) -> tuple[float, np.ndarray]:
    """Total flux and per-star fluxes for one instant (positions already light-time
    corrected), in units of the summed unocculted fluxes."""
    base = np.array([math.pi * s.radius_au**2 * s.weight for s in lights])
    f = base.copy()
    if proximity:
        ell, refl, beam = proximity_terms(pos, vel, masses, lights, pairs, base)
        f = (base + base * ell + refl) * beam
    f = f * visible_fractions(pos, lights)  # eclipses
    return float(f.sum() / base.sum()), f / base.sum()


# dense output of the integrations of SeenStates: states at grid nodes spaced by DENSE_STEPS per
# periastron time scale of the fastest orbit of the stars read, cubic Hermite interpolation in
# between (position error ~ (w h)^4 / 384 a: < 1e-9 a for 256 steps; velocity ~ (w h)^3)
DENSE_STEPS = 256


class _Dense:
    """States of a REBOUND simulation at arbitrary times from nodes of a grid of spacing h (cubic
    Hermite interpolation): one exact-finish integration per node instead of one per observing
    time, which costs IAS15 a shortened step each (a light curve has many points per orbit).
    The last MAX_NODES nodes are kept, so the groups of SeenStates (each read at its own emission
    time, light-minutes to hours apart) share one integration."""

    MAX_NODES = 64

    def __init__(self, sim, h: float, advance=None):
        self.sim, self.h = sim, float(h)
        self.advance = advance or (lambda t, exact=1: sim.integrate(t, exact_finish_time=exact))
        self.nodes: dict[int, tuple] = {}  # grid index -> (t, pos, vel, m)

    def _node(self, k: int):
        node = self.nodes.get(k)
        if node is None:
            t = k * self.h
            self.advance(t)
            n = self.sim.N
            xyz, vxyz, m = np.empty((n, 3)), np.empty((n, 3)), np.empty(n)
            self.sim.serialize_particle_data(xyz=xyz, vxvyvz=vxyz, m=m)
            node = self.nodes[k] = (t, xyz, vxyz, m)
            if len(self.nodes) > self.MAX_NODES:  # drop the node farthest from this one
                del self.nodes[max(self.nodes, key=lambda q: abs(q - k))]
        return node

    def state(self, t: float, idx):
        k = math.floor(t / self.h)
        # the node nearer to the simulation's present time first (it moves on monotonically)
        if abs((k + 1) * self.h - self.sim.t) < abs(k * self.h - self.sim.t):
            b = self._node(k + 1)
            a = self._node(k)
        else:
            a = self._node(k)
            b = self._node(k + 1)
        ta, pa, va, ma = a
        tb, pb, vb, mb = b
        dt = tb - ta
        u = (t - ta) / dt
        u2, u3 = u * u, u * u * u
        h00, h10, h01, h11 = 2 * u3 - 3 * u2 + 1, u3 - 2 * u2 + u, -2 * u3 + 3 * u2, u3 - u2
        d00, d10, d01, d11 = 6 * u2 - 6 * u, 3 * u2 - 4 * u + 1, -6 * u2 + 6 * u, 3 * u2 - 2 * u
        pa_, pb_, va_, vb_ = pa[idx], pb[idx], va[idx], vb[idx]
        pos = h00 * pa_ + h10 * dt * va_ + h01 * pb_ + h11 * dt * vb_
        vel = (d00 * pa_ + d10 * dt * va_ + d01 * pb_ + d11 * dt * vb_) / dt
        m = ma[idx] + u * (mb[idx] - ma[idx])
        return pos, vel, m


def _dense_step(root, masses, stars) -> float:
    """Node spacing [yr] for reading the stars `stars`: DENSE_STEPS per periastron time scale
    P (1 - e)^1.5 / (1 + e)^0.5 of the fastest orbit with one of them as a direct child."""
    best = math.inf
    for n in root.nodes():
        kids = [c for c in (n.left, n.right) if isinstance(c, int)]
        if not any(k in stars for k in kids):
            continue
        from .hierarchy import mass_of

        el = n.elements
        P = units.sma_to_period(el.a, mass_of(n.left, masses) + mass_of(n.right, masses))
        e = min(max(float(el.e), 0.0), 0.95)
        best = min(best, P * (1.0 - e) ** 1.5 / math.sqrt(1.0 + e))
    return best / DENSE_STEPS if math.isfinite(best) else math.inf


class SeenStates:
    """Positions and velocities of all stars as seen at observing times (monotonic calls).

    With the light-travel time every star is seen at its own emission time, x_i(t - z_i/c).
    A first-order correction x - v z/c is exact enough for slowly moving bodies, but not for an
    inner binary whose centre of mass is far from the system's centre (QZ Car: z/c ~ 0.2 d for a
    6-d binary; the third-order term shifts eclipses by ~200 s). Every two-star binary therefore
    has its own integration, evaluated at the emission time of its centre of mass,
    t_e = t - z_com(t_e)/c (solved iteratively), and only the small offsets of its two stars from
    that centre are corrected to first order (the Romer delay inside the binary). Stars outside
    two-star binaries use the first-order correction."""

    def __init__(self, cfg, light_time: bool = True, dense: bool = True):
        from .backends.nbody import build_simulation

        self.light_time = light_time
        from .pdot import integrator

        sim, _, _, root, masses = build_simulation(cfg)
        sim.collision = "none"
        self.n = len(masses)
        # integration of every simulation with the kicks of orbits.<pair>.Pdot
        self._advance = {id(sim): integrator(cfg, sim, root, masses)}
        self.root = root
        pairs = [
            (n.left, n.right)
            for n in root.nodes()
            if isinstance(n.left, int) and isinstance(n.right, int)
        ]
        self.pairs = pairs
        self.emission = {}  # (i, j) -> (t_e [yr], positions, velocities, masses) unretarded
        in_pairs = {k for p in pairs for k in p}
        rest = [k for k in range(self.n) if k not in in_pairs]
        self.groups = []
        if not light_time:
            self.groups.append(("none", list(range(self.n)), sim, None))
        if rest and light_time:
            self.groups.append(("rest", rest, sim, None))
        h = _dense_step(root, masses, list(range(self.n))) if dense else math.inf
        for i, j in pairs if light_time else []:
            # dense output: every group reads the one integration at its own times; otherwise
            # every two-star binary has its own integration (one exact step per observing time)
            s = (
                sim
                if math.isfinite(h) or not rest and not any(g[2] is sim for g in self.groups)
                else None
            )
            if s is None:
                s = build_simulation(cfg)[0]
                s.collision = "none"
                self._advance[id(s)] = integrator(cfg, s, root, masses)
            self.groups.append(("binary", [i, j], s, [None]))
        self.dense = {id(sim): _Dense(sim, h, self._advance[id(sim)])} if math.isfinite(h) else {}

    @staticmethod
    def _state(sim, idx):
        n = sim.N
        xyz, vxyz, m = np.empty((n, 3)), np.empty((n, 3)), np.empty(n)
        sim.serialize_particle_data(xyz=xyz, vxvyvz=vxyz, m=m)
        return xyz[idx], vxyz[idx], m[idx]

    def _read(self, sim, idx, t: float):
        """(positions, velocities, masses) of the stars idx of sim at time t."""
        d = self.dense.get(id(sim))
        if d is not None:
            return d.state(t, idx)
        self._advance[id(sim)](t)
        return self._state(sim, idx)

    def at(self, t_yr: float) -> tuple[np.ndarray, np.ndarray]:
        pos = np.empty((self.n, 3))
        vel = np.empty((self.n, 3))
        t_emit = np.full(self.n, t_yr)
        c = units.C
        for kind, idx, sim, last in self.groups:
            if kind == "binary":
                # t_e = t - z_com(t_e)/c: predicted from the last state (linear motion of the
                # centre of mass), then corrected until it is consistent (usually once)
                if last[0] is None:
                    p, v, m = self._state(sim, idx)
                    last[0] = (sim.t, float(m @ p[:, 2] / m.sum()), float(m @ v[:, 2] / m.sum()))
                te0, z0, vz0 = last[0]
                t_e = (t_yr - (z0 - vz0 * te0) / c) / (1.0 + vz0 / c)
                for _ in range(10):
                    p, v, m = self._read(sim, idx, t_e)
                    zb = float(m @ p[:, 2] / m.sum())
                    vzb = float(m @ v[:, 2] / m.sum())
                    new = (t_yr - (zb - vzb * t_e) / c) / (1.0 + vzb / c)
                    if abs(new - t_e) < 1e-12:
                        break
                    t_e = new
                last[0] = (t_e, zb, vzb)
                self.emission[tuple(idx)] = (t_e, p.copy(), v.copy(), m)
                t_emit[idx] = t_e - (p[:, 2] - zb) / c
                p = p - v * ((p[:, 2:3] - zb) / c)  # Romer delay inside the binary
            else:
                p, v, _ = self._read(sim, idx, t_yr)
                if kind == "rest":
                    t_emit[idx] = t_yr - p[:, 2] / c
                    p = p - v * (p[:, 2:3] / c)
                else:
                    ms = np.array([sim.particles[k].m for k in idx])
                    for i, j in self.pairs:
                        a, b = idx.index(i), idx.index(j)
                        self.emission[(i, j)] = (t_yr, p[[a, b]], v[[a, b]], ms[[a, b]])
            pos[idx], vel[idx] = p, v
        self.t_emit = t_emit  # emission time of every star [yr]
        return pos, vel


def light_curve(
    cfg,
    t_day: np.ndarray,
    band: str = "TESS",
    ld: dict | None = None,
    exposure_min: float = 0.0,
    n_sub: int = 1,
    proximity: bool = True,
    light_time: bool = True,
) -> LightCurve:
    """Synthetic light curve at the observing times t_day [d after t = 0] (N-body positions)."""
    from . import massloss
    from .backends.nbody import model_props
    from .hierarchy import build_tree, initialize

    cfg = massloss.resolve_masses(cfg)
    seen = SeenStates(cfg, light_time)
    m_lab = np.array(build_tree(cfg)[1], float)  # track labels (initial masses)
    root, masses = seen.root, massloss.initial_masses(cfg, m_lab)
    props = model_props(cfg, root, m_lab, *initialize(root, masses), masses)
    lights = star_lights(cfg, masses, props, band, ld)
    from .pulsation import pulsators

    puls = pulsators(cfg, masses, props)
    from .discs import build_discs, transmission

    discs = build_discs(cfg, root)
    disc_cache: dict = {}
    from .discs import disc_light

    shining = any(d.spec.emission or d.spec.albedo > 0 for d in discs)
    teffs = [p.teff_k for p in props]
    pairs = seen.pairs
    t_day = np.asarray(t_day, float)
    if exposure_min > 0 and n_sub > 1:
        offs = (np.arange(n_sub) + 0.5) / n_sub - 0.5
        t_eval = (t_day[:, None] + offs[None, :] * exposure_min / 1440.0).ravel()
    else:
        n_sub, t_eval = 1, t_day
    order = np.argsort(t_eval)
    flux = np.empty(len(t_eval))
    comp = np.empty((len(t_eval), len(masses)))
    nstar = len(masses)
    pfac = np.ones((len(t_eval), len(puls)))
    dl = np.zeros((len(t_eval), len(discs)))
    scene = None
    from . import mesh

    if mesh.needed(cfg):
        p0, v0 = initialize(root, masses)
        scene = mesh.build_scene(cfg, root, masses, props, lights, p0, v0, band)
    base_sum = sum(math.pi * s_.radius_au**2 * s_.weight for s_ in lights)
    for k in order:
        pos, vel = seen.at(t_eval[k] / units.DAYS_PER_YEAR)
        if scene is not None:
            flux[k], comp[k] = mesh.scene_flux(*scene, pos, vel, masses, lights, pairs, proximity)
        else:
            flux[k], comp[k] = flux_at(pos, vel, masses, lights, pairs, proximity)
        if discs:
            comp[k] *= transmission(
                discs, t_eval[k] / units.DAYS_PER_YEAR, pos, masses, lights, disc_cache
            )
            flux[k] = comp[k].sum()
            if shining:
                dl[k] = (
                    disc_light(
                        discs, t_eval[k] / units.DAYS_PER_YEAR, pos, masses, lights, teffs, band
                    )
                    / base_sum
                )
                flux[k] += dl[k].sum()
        if puls:
            t_e = seen.t_emit * units.DAYS_PER_YEAR
            for q, pz in enumerate(puls):
                pfac[k, q] = pz.flux_factor(t_e[pz.star])
                comp[k, pz.star] *= pfac[k, q]
            flux[k] = comp[k].sum() + dl[k].sum()
    flux = flux.reshape(-1, n_sub).mean(axis=1)
    comp = comp.reshape(-1, n_sub, nstar).mean(axis=1)
    components = {f"m{k + 1}": comp[:, k] for k in range(nstar)}
    if shining:
        dl = dl.reshape(-1, n_sub, len(discs)).mean(axis=1)
        for q, d in enumerate(discs):
            components[f"disc_{d.name}"] = dl[:, q]
    meta = {
        "engine": "tier A+, surface meshes" if scene is not None else "tier A, spheres",
        **time_meta(cfg),
        "light_time": light_time,
        "positions": (
            f"N-body states at nodes every {seen.dense[next(iter(seen.dense))].h * units.DAYS_PER_YEAR:.5g} d, "
            "cubic Hermite interpolation in between"
            if seen.dense
            else "N-body integration to every observing time"
        ),
        "proximity": proximity,
        "exposure_min": float(exposure_min),
        "n_sub": int(n_sub),
        "normalisation": "1 = summed flux of all stars out of eclipse, without proximity "
        "effects (pi R^2 times the passband surface brightness)",
        "ld_law": "Claret (2000) four-coefficient law, I(mu)/I(1) = 1 - sum_k a_k "
        "(1 - mu^(k/2))" + (" (given per star)" if ld else " (Claret & Bloemen 2011 / TLUSTY)"),
    }
    return LightCurve(
        t_day=t_day,
        flux=flux,
        band=band,
        components=components,
        stars=lights,
        pulsation={
            f"m{pz.star + 1}": pfac[:, q].reshape(-1, n_sub).mean(axis=1)
            for q, pz in enumerate(puls)
        },
        meta=meta,
    )


def time_meta(cfg) -> dict:
    """BJD of t = 0 and the time system of the configuration (fit.stage2.bjd_t0,
    fit.time_system), as far as given."""
    fit = getattr(cfg, "fit", None)
    st = getattr(fit, "stage2", None) if fit is not None else None
    out = {"time_system": getattr(fit, "time_system", "BJD_TDB") if fit else "BJD_TDB"}
    if st is not None and getattr(st, "bjd_t0", None) is not None:
        out["bjd_t0"] = float(st.bjd_t0)
    return out


def add_noise(lc: LightCurve, ppm: float, seed: int = 1) -> LightCurve:
    rng = np.random.default_rng(seed)
    return LightCurve(
        lc.t_day,
        lc.flux + rng.normal(0.0, ppm * 1e-6, len(lc.flux)),
        lc.band,
        lc.components,
        lc.stars,
        lc.pulsation,
        {**lc.meta, "noise_ppm": float(ppm), "noise_seed": seed},
    )
