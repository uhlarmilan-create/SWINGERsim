"""Light curves of multiple stars with PHOEBE 2 for every inner binary (F4 b, tier C).

PHOEBE (Prša et al. 2016; Horvat et al. 2018; Conroy et al. 2020) models one binary: Roche
geometry, gravity darkening, irradiation, model atmospheres, light-travel time inside the
binary. swingersim models the multiple system. The hybrid joins them:

1. the N-body integration of swingersim gives, in short chunks (default one binary period), the
   osculating Keplerian elements of every inner binary at the emission time of its centre of
   mass for the chunk centre, t_e = t_c - z_com(t_e)/c (photometry.SeenStates), and the
   velocity of the centre of mass along the line of sight;
2. PHOEBE computes the light of the binary through the chunk with these elements, at the
   emission times t_e + (t - t_c)/(1 + v_z/c) (light-time effect of the binary in the system);
3. stars outside two-star binaries shine with their passband luminosity / 4 pi;
4. eclipses between stars of different binaries are taken from the geometry of swingersim
   (limb-darkened spheres, tier A) and subtracted.

Conventions: swingersim's z axis points away from the observer and its elements describe the
orbit of the second star about the first; PHOEBE's argument of periastron is that of the
primary, so per0 = omega + 180 deg, incl = i. Fluxes are normalised to the summed passband
luminosities of all stars / 4 pi (1 = all stars outside eclipse, as in tier A); PHOEBE's mesh
makes the absolute level a few 0.1 % lower (1500 triangles per star).

PHOEBE is an optional dependency (``pip install phoebe``). Passbands other than Johnson:V
are downloaded on first use (``phoebe.download_passband``).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from . import units

PASSBANDS = {
    "U": "Johnson:U",
    "B": "Johnson:B",
    "V": "Johnson:V",
    "R": "Cousins:Rc",
    "I": "Cousins:Ic",
    "TESS": "TESS:T",
    "Kepler_mean": "Kepler:mean",
    "Gaia_G": "Gaia:G",
    "Gaia_BP": "Gaia:BP",
    "Gaia_RP": "Gaia:RP",
}


def _phoebe():
    try:
        import phoebe
    except ImportError as exc:  # pragma: no cover - depends on the installation
        raise ImportError(
            "PHOEBE light curves need the phoebe package: pip install phoebe"
        ) from exc
    phoebe.logger(clevel="ERROR")
    return phoebe


def passband(band: str) -> str:
    """PHOEBE passband of a swingersim band (installed or downloaded)."""
    phoebe = _phoebe()
    pb = PASSBANDS.get(band, band)
    if pb not in phoebe.list_installed_passbands():
        try:
            phoebe.download_passband(pb)
        except Exception as exc:
            raise RuntimeError(
                f"PHOEBE passband {pb} is not installed and could not be downloaded ({exc}); "
                f"run phoebe.download_passband('{pb}') on a computer with internet access"
            ) from exc
    return pb


@dataclass
class ChunkElements:
    a_au: float
    e: float
    inc_deg: float
    Omega_deg: float
    omega_deg: float  # of the second star about the first (swingersim)
    period_d: float
    q: float  # m2 / m1
    t_e_d: float  # emission (dynamical) time of the chunk centre [d]
    vz_au_d: float  # velocity of the centre of mass along the line of sight [AU/d]
    t_peri_d: float = math.nan  # time of periastron passage (dynamical) [d]


def _settings(
    b,
    ds: str,
    pb: str,
    star_teff: dict,
    atm: str,
    distortion: str,
    irradiation: bool,
    ld: dict | None,
    ntriangles: int,
) -> None:
    b.set_value(f"passband@{ds}", pb)
    for comp, teff in star_teff.items():
        b.set_value(f"teff@{comp}@component", teff)
        radiative = teff >= 7000.0  # same rule as tier A (Rucinski 1969)
        b.set_value(f"gravb_bol@{comp}", 1.0 if radiative else 0.32)
        b.set_value(f"irrad_frac_refl_bol@{comp}", 1.0 if radiative else 0.5)
        b.set_value(f"distortion_method@{comp}@compute", distortion)
        b.set_value(f"ntriangles@{comp}@compute", ntriangles)
        b.set_value(f"atm@{comp}@compute", atm)
        if atm == "blackbody" or ld is not None:
            b.set_value(f"ld_mode@{comp}@{ds}", "manual")
            b.set_value(f"ld_func@{comp}@{ds}", "power")
            b.set_value(f"ld_coeffs@{comp}@{ds}", list(ld[comp]))
            b.set_value(f"ld_mode_bol@{comp}", "manual")
            b.set_value(f"ld_func_bol@{comp}", "linear")
            b.set_value(f"ld_coeffs_bol@{comp}", [0.5])
    b.set_value("irrad_method@compute", "horvat" if irradiation else "none")
    b.set_value("pblum_mode@" + ds, "absolute")
    b.set_value("distance@system", 1.0)


def _binary_bundle(pb, teffs, atm, distortion, irradiation, ld, ntriangles, light_time):
    phoebe = _phoebe()
    b = phoebe.default_binary()
    b.add_dataset("lc", compute_times=[0.0], dataset="lc01")
    b.flip_constraint("t0_perpass@binary", solve_for="t0_supconj")
    _settings(
        b,
        "lc01",
        pb,
        {"primary": teffs[0], "secondary": teffs[1]},
        atm,
        distortion,
        irradiation,
        None if ld is None else {"primary": ld[0], "secondary": ld[1]},
        ntriangles,
    )
    b.set_value_all("ltte", bool(light_time))
    return b


def _single_pblum(pb, teff, radius_rsun, atm, ld) -> float:
    """Passband luminosity / 4 pi of a single spherical star (flux at 1 m)."""
    phoebe = _phoebe()
    b = phoebe.default_star()
    b.add_dataset("lc", compute_times=[0.0], dataset="lc01")
    b.set_value("requiv@starA", radius_rsun)
    _settings(
        b,
        "lc01",
        pb,
        {"starA": teff},
        atm,
        "sphere",
        False,
        None if ld is None else {"starA": ld},
        500,
    )
    pbl = b.compute_pblums(pblum_abs=True)
    return float(next(v.value for k, v in pbl.items() if k.startswith("pblum_abs"))) / (4 * math.pi)


def _mesh_correction(b, pbl_sum: float) -> float:
    """Ratio of the exact flux of two uneclipsed spheres (passband luminosities / 4 pi) to
    PHOEBE's flux of the same spheres on its mesh; restores the bundle's settings."""
    saved = {}
    for comp in ("primary", "secondary"):
        saved[comp] = b.get_value(f"distortion_method@{comp}@compute")
        b.set_value(f"distortion_method@{comp}@compute", "sphere")
    irr = b.get_value("irrad_method@compute")
    ecl = b.get_value("eclipse_method@compute")
    times = b.get_value("compute_times@lc01")
    b.set_value("irrad_method@compute", "none")
    b.set_value("eclipse_method@compute", "only_horizon")
    b.set_value("compute_times@lc01", times[:1])
    b.run_compute(progressbar=False)
    f0 = float(b["value@fluxes@lc01@model"][0])
    for comp, v in saved.items():
        b.set_value(f"distortion_method@{comp}@compute", v)
    b.set_value("irrad_method@compute", irr)
    b.set_value("eclipse_method@compute", ecl)
    b.set_value("compute_times@lc01", times)
    return pbl_sum / (4.0 * math.pi) / f0


def phoebe_light_curve(
    cfg,
    t_day: np.ndarray,
    band: str = "V",
    atm: str = "auto",
    distortion: str = "roche",
    irradiation: bool = True,
    light_time: bool = True,
    chunk_days: float | None = None,
    ntriangles: int = 1500,
    swingersim_ld: bool = False,
    ld: dict | None = None,
    exposure_min: float = 0.0,
    n_sub: int = 1,
    proximity: bool = True,
):
    """Light curve at the observing times t_day [d] with PHOEBE for every inner binary.

    atm "ck2004" (Castelli & Kurucz atmospheres, PHOEBE's limb darkening), "blackbody" or
    "auto" (ck2004, black bodies for binaries outside its grid, e.g. O giants);
    distortion "roche" or "sphere"; swingersim_ld=True (always with blackbody) uses the limb
    darkening of tier A (Claret four-coefficient law = PHOEBE's power law), which makes the
    two tiers directly comparable; ``ld`` (per-star coefficients, as for photometry.light_curve)
    are used in the same way. exposure_min / n_sub integrate each exposure over n_sub
    sub-exposures (as photometry.light_curve). proximity=False cannot switch off the Roche
    distortion and reflection of PHOEBE separately: it is rejected unless distortion='sphere'
    and irradiation=False (which is what it means)."""
    if not proximity and (distortion != "sphere" or irradiation):
        raise ValueError(
            "the PHOEBE engine has no proximity switch: use distortion 'sphere' and no "
            "irradiation (--phoebe-distortion sphere --no-irradiation) for a light curve "
            "without proximity effects"
        )
    t_day = np.asarray(t_day, float)
    if exposure_min > 0 and n_sub > 1:  # integration over the exposure: sub-exposures
        offs = (np.arange(n_sub) + 0.5) / n_sub - 0.5
        t_eval = (t_day[:, None] + offs[None, :] * exposure_min / 1440.0).ravel()
        lc = phoebe_light_curve(
            cfg,
            t_eval,
            band,
            atm,
            distortion,
            irradiation,
            light_time,
            chunk_days,
            ntriangles,
            swingersim_ld,
            ld,
        )
        lc.t_day = t_day
        lc.flux = lc.flux.reshape(-1, n_sub).mean(axis=1)
        lc.components = {k: v.reshape(-1, n_sub).mean(axis=1) for k, v in lc.components.items()}
        lc.meta.update(exposure_min=float(exposure_min), n_sub=int(n_sub))
        return lc
    if any(s.pulsation is not None for s in cfg.stars.values()):
        import warnings

        warnings.warn(
            "pulsations are ignored by the PHOEBE engine (use engine 'swingersim')", stacklevel=2
        )
    if cfg.discs:
        import warnings

        warnings.warn(
            "dust discs are ignored by the PHOEBE engine (use engine 'swingersim')", stacklevel=2
        )
    from . import massloss
    from .backends.nbody import build_simulation, model_props
    from .hierarchy import build_tree, initialize
    from .kepler import state_to_elements
    from .photometry import LightCurve, SeenStates, blocked_fraction, star_lights

    pb = passband(band)
    cfg = massloss.resolve_masses(cfg)
    sim, _, _, root, masses = build_simulation(cfg)
    sim.collision = "none"
    m_lab = np.array(build_tree(cfg)[1], float)
    props = model_props(cfg, root, m_lab, *initialize(root, masses), masses)
    lights = star_lights(cfg, masses, props, band if band in PASSBANDS else "V", ld)
    for s_ in lights:  # PHOEBE's power law takes four coefficients (missing terms are zero)
        s_.ld = np.pad(np.asarray(s_.ld, float), (0, max(4 - len(s_.ld), 0)))
    n_star = len(masses)
    pairs = [
        (n.name, n.left, n.right)
        for n in root.nodes()
        if isinstance(n.left, int) and isinstance(n.right, int)
    ]
    partner = {}
    for _, i, j in pairs:
        partner[i], partner[j] = j, i
    singles = [k for k in range(n_star) if k not in partner]
    use_ld = swingersim_ld or ld is not None

    order = np.argsort(t_day)
    ts = t_day[order]
    # chunks per binary: centres at the middle of each occupied interval
    chunks = {}
    p_nom = {}
    for name, i, j in pairs:
        p_d = (
            units.sma_to_period(
                next(n for n in root.nodes() if n.name == name).elements.a, masses[i] + masses[j]
            )
            * units.DAYS_PER_YEAR
        )
        width = chunk_days or p_d
        p_nom[name] = p_d
        idx = np.floor(ts / width).astype(int)
        _, start = np.unique(idx, return_index=True)
        stops = list(start[1:]) + [len(ts)]
        chunks[name] = [
            (0.5 * (ts[s] + ts[e - 1]), s, e) for s, e in zip(start, stops, strict=True)
        ]
    events = [(t, 1, k) for k, t in enumerate(ts)]
    for name, lst in chunks.items():
        events += [(tc, 0, (name, c)) for c, (tc, _, _) in enumerate(lst)]
        # mean motion over one orbit around every chunk centre (osculating elements of a close,
        # tidally distorted binary give a Keplerian period off by minutes per orbit)
        h = 0.5 * p_nom[name]
        events += [
            (tc + sg * h, 2, (name, c, sg)) for c, (tc, _, _) in enumerate(lst) for sg in (-1, 1)
        ]
    events.sort(key=lambda x: (x[0], x[1]))

    seen = SeenStates(cfg, light_time)
    elems: dict[str, list] = {name: [None] * len(lst) for name, lst in chunks.items()}
    lam: dict[tuple, tuple[float, float]] = {}  # (name, c, sign) -> (t_e [d], mean longitude)
    # eclipses between different binaries: the blocked fraction of each star's own light
    # (limb-darkened spheres, tier A), applied to the star's PHOEBE flux below
    blocked = np.zeros((len(ts), n_star))
    for t, kind, key in events:
        pos, _ = seen.at(t / units.DAYS_PER_YEAR)
        if kind == 2:
            name = key[0]
            _, i, j = next(p for p in pairs if p[0] == name)
            t_e, p2, v2, m2 = seen.emission[(i, j)]
            el = state_to_elements(units.G * m2.sum(), p2[1] - p2[0], v2[1] - v2[0])
            lam[key] = (t_e * units.DAYS_PER_YEAR, el.Omega + el.omega + el.M)
            continue
        if kind == 0:
            name, c = key
            _, i, j = next(p for p in pairs if p[0] == name)
            # the binary at the emission time of its centre of mass (dynamical time t_e)
            t_e, p2, v2, m2 = seen.emission[(i, j)]
            el = state_to_elements(units.G * m2.sum(), p2[1] - p2[0], v2[1] - v2[0])
            vzc = float(m2 @ v2[:, 2] / m2.sum()) if light_time else 0.0
            elems[name][c] = ChunkElements(
                el.a,
                el.e,
                math.degrees(el.inc),
                math.degrees(el.Omega),
                math.degrees(el.omega),
                units.sma_to_period(el.a, m2.sum()) * units.DAYS_PER_YEAR,
                m2[1] / m2[0],
                t_e * units.DAYS_PER_YEAR,  # emission (dynamical) time of the chunk centre [d]
                vzc / units.DAYS_PER_YEAR,
            )
            continue
        z_order = np.argsort(pos[:, 2])
        for ib, k in enumerate(z_order):
            occ = []
            for j in z_order[:ib]:
                if partner.get(k) == j:
                    continue  # eclipses inside a binary come from PHOEBE
                dx, dy = pos[j, 0] - pos[k, 0], pos[j, 1] - pos[k, 1]
                if math.hypot(dx, dy) < lights[j].radius_au + lights[k].radius_au:
                    occ.append((dx, dy, lights[j].radius_au))
            if occ:
                blocked[key, k] = blocked_fraction(lights[k].radius_au, lights[k].ld, occ)

    for name, lst in elems.items():
        for c, el in enumerate(lst):
            (t0, l0), (t1, l1) = lam[(name, c, -1)], lam[(name, c, 1)]
            n_kep = 2.0 * math.pi / el.period_d
            dl = l1 - l0
            dl += 2.0 * math.pi * round((n_kep * (t1 - t0) - dl) / (2.0 * math.pi))
            n_eff = dl / (t1 - t0)  # mean motion over one orbit [rad/d]
            el.period_d = 2.0 * math.pi / n_eff
            # phase from the mean longitude averaged over that orbit (the osculating mean
            # anomaly of a distorted close binary oscillates within the orbit)
            lam_mid = l0 + 0.5 * dl
            m_mid = lam_mid - math.radians(el.Omega_deg + el.omega_deg)
            el.t_peri_d = 0.5 * (t0 + t1) - m_mid / n_eff

    c_au_d = units.C / units.DAYS_PER_YEAR
    comps = {}
    star_flux = np.zeros(n_star)  # unocculted flux of every star in PHOEBE's units
    norm = 0.0

    def binary_flux(name, i, j, atm_b):
        use = use_ld or atm_b == "blackbody"
        b = _binary_bundle(
            pb,
            (props[i].teff_k, props[j].teff_k),
            atm_b,
            distortion,
            irradiation,
            (lights[i].ld, lights[j].ld) if use else None,
            ntriangles,
            light_time,
        )
        b.set_value("requiv@primary", props[i].radius_rsun)
        b.set_value("requiv@secondary", props[j].radius_rsun)
        flux = np.empty(len(ts))
        pbl_sum, corr = None, 1.0
        for (tc, s, e), el in zip(chunks[name], elems[name], strict=True):
            b.set_value("period@binary", el.period_d)
            b.set_value("q@binary", el.q)
            b.set_value("sma@binary", el.a_au / units.RSUN_AU)
            b.set_value("ecc@binary", el.e)
            b.set_value("incl@binary", el.inc_deg)
            b.set_value("long_an@binary", el.Omega_deg % 360.0)
            b.set_value("per0@binary", (el.omega_deg + 180.0) % 360.0)
            t_obs = ts[s:e]
            # emission times: t_e(tc) + (t - tc) / (1 + v_z / c) (centre of mass moving along z)
            t_em = el.t_e_d + (t_obs - tc) / (1.0 + el.vz_au_d / c_au_d)
            b.set_value("t0_perpass@binary", el.t_peri_d)
            b.set_value("compute_times@lc01", t_em)
            b.run_compute(progressbar=False)
            flux[s:e] = b["value@fluxes@lc01@model"]
            if pbl_sum is None:
                pbl = b.compute_pblums(pblum_abs=True)
                pbl_sum = sum(v.value for k, v in pbl.items() if k.startswith("pblum_abs"))
                for k, v in pbl.items():
                    if k.startswith("pblum_abs"):
                        star_flux[i if "primary" in k else j] = v.value / (4.0 * math.pi)
                corr = _mesh_correction(b, pbl_sum)
            flux[s:e] *= corr
        return flux, pbl_sum

    for name, i, j in pairs:
        if atm == "auto":
            try:
                flux, pbl_sum = binary_flux(name, i, j, "ck2004")
            except ValueError as exc:
                if "extrapolation" not in str(exc):
                    raise
                import warnings

                warnings.warn(
                    f"binary {name}: outside the ck2004 atmospheres of PHOEBE (hot or low-gravity "
                    "stars); black bodies with the limb darkening of tier A (Claret / TLUSTY)",
                    stacklevel=2,
                )
                flux, pbl_sum = binary_flux(name, i, j, "blackbody")
        else:
            flux, pbl_sum = binary_flux(name, i, j, atm)
        comps[name] = flux
        norm += pbl_sum / (4.0 * math.pi)
    for k in singles:
        args = (pb, props[k].teff_k, props[k].radius_rsun)
        try:
            a_k = "ck2004" if atm == "auto" else atm
            f = _single_pblum(*args, a_k, lights[k].ld if (use_ld or a_k == "blackbody") else None)
        except ValueError:
            if atm != "auto":
                raise
            f = _single_pblum(*args, "blackbody", lights[k].ld)
        comps[f"m{k + 1}"] = np.full(len(ts), f)
        star_flux[k] = f
        norm += f
    # the light lost in eclipses between binaries, in the same (PHOEBE) units as the fluxes:
    # each star's blocked fraction times its own share of the light
    inter = blocked @ star_flux / norm
    total = sum(comps.values()) / norm - inter
    out = np.empty(len(ts))
    out[order] = total
    comp_out = {}
    for k, v in comps.items():
        arr = np.empty(len(ts))
        arr[order] = v / norm
        comp_out[k] = arr
    if np.any(inter > 0):
        arr = np.empty(len(ts))
        arr[order] = -inter
        comp_out["mutual"] = arr
    from .photometry import time_meta

    meta = {
        "engine": f"PHOEBE {getattr(_phoebe(), '__version__', '')} for the inner binaries "
        f"(tier C; atm {atm}, distortion {distortion}, irradiation "
        f"{'on' if irradiation else 'off'}, {ntriangles} triangles)",
        **time_meta(cfg),
        "light_time": light_time,
        "exposure_min": 0.0,
        "n_sub": 1,
        "normalisation": "1 = summed passband luminosities / 4 pi of all stars (PHOEBE)",
        "ld_law": "PHOEBE (ck2004 atmospheres)"
        if not use_ld and atm != "blackbody"
        else "power law = Claret (2000) four-coefficient law, coefficients below",
    }
    return LightCurve(
        t_day=t_day, flux=out, band=band, components=comp_out, stars=lights, meta=meta
    )
