"""Benchmarks: Roche geometry and surface meshes (rotating and contact stars),
contact-binary evolution, rejuvenation of accretors, rotation-corrected k2, pre-main-sequence
tides, rigid discs around binaries and stars (appended to swingersim.validation.BENCHMARKS)."""

from __future__ import annotations

import math
from functools import cache

import numpy as np

from . import units
from ._benchmark import Benchmark

# --------------------------------------------------------------------------- Roche geometry


def roche_lobe_vs_eggleton(q: float = 0.3) -> float:
    """Volume-equivalent radius of the Roche lobe of the primary (q = M2/M1) from the surface
    mesh of swingersim.roche, divided by Eggleton's fit."""
    from .roche import roche_lobe_volume_radius
    from .stellar import roche_lobe_eggleton

    return roche_lobe_volume_radius(q, 0) / roche_lobe_eggleton(1.0 / q)


def rotating_critical_flattening() -> float:
    from .roche import rotating_radius

    return float(rotating_radius(np.array([math.pi / 2]), 1.0)[0])


def rotating_equator(w: float = 0.6) -> float:
    """R_eq/R_pol - (1 + w_K^2 / 2), w_K = Omega / Omega_K(R_eq): zero for the Roche model."""
    from .roche import rotating_radius

    x = float(rotating_radius(np.array([math.pi / 2]), w)[0])
    wk2 = w * w * 8.0 / 27.0 * x**3
    return x - (1.0 + 0.5 * wk2)


def _two_sphere_cfg():
    from .config import Config

    return Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"A": ["m1", "m2"]},
            "masses": {"m1": 1.0, "m2": 0.8},
            "orbits": {"A": {"P": 5.0, "i": 89.0}},
            "stars": {"m1": {"R": 1.0, "Teff": 5800}, "m2": {"R": 0.8, "Teff": 5000}},
            "physics": {"gr": False},
            "run": {"t_max_yr": 0.1},
        }
    )


def mesh_sphere_eclipse() -> float:
    """Largest difference of the light curves of two spheres through both eclipses: surface
    meshes (swingersim.mesh) against the exact sphere eclipses of tier A."""
    from . import mesh, photometry

    cfg = _two_sphere_cfg()
    t = np.concatenate([np.linspace(1.1, 1.4, 61), np.linspace(3.6, 3.9, 61)])
    a = photometry.light_curve(cfg, t, "V", proximity=False, light_time=False)
    orig = mesh.needed
    mesh.needed = lambda c: True
    try:
        b = photometry.light_curve(cfg, t, "V", proximity=False, light_time=False)
    finally:
        mesh.needed = orig
    d = b.flux - a.flux
    return float(np.max(np.abs(d - np.median(d))))


def _contact_cfg(q=0.42, f=0.3, incl=80.0):
    from .config import Config

    return Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"A": ["m1", "m2"]},
            "masses": {"m1": 1.2, "m2": 1.2 * q},
            "orbits": {"A": {"P": 0.4775, "i": incl, "contact_fillout": f}},
            "stars": {"m1": {"R": 1.0, "Teff": 6200}, "m2": {"R": 0.7, "Teff": 6200}},
            "physics": {"gr": False},
            "run": {"t_max_yr": 0.01},
        }
    )


def contact_amplitude_difference() -> float:
    """V amplitude of the contact binary, swingersim minus PHOEBE [mag]."""
    return contact_amplitude_swingersim() - contact_amplitude_phoebe()


def contact_amplitude_swingersim() -> float:
    from .photometry import light_curve

    cfg = _contact_cfg()
    t = np.linspace(0.0, 0.4775, 193)
    lc = light_curve(cfg, t, "V", proximity=False, light_time=False)
    return float(-2.5 * math.log10(lc.flux.min() / lc.flux.max()))


@cache
def contact_amplitude_phoebe() -> float:
    """Amplitude [mag] of the V light curve of the same contact binary with PHOEBE 2 (black
    bodies, no irradiation, the same four-coefficient limb darkening and gravity darkening)."""
    import phoebe

    from .photometry import limb_darkening

    q, f, incl = 0.42, 0.3, 80.0
    b = phoebe.default_binary(contact_binary=True)
    b.set_value("q", context="component", value=q)
    b.set_value("period", component="binary", context="component", value=0.4775)
    b.set_value("incl", component="binary", context="component", value=incl)
    b.flip_constraint("pot", component="contact_envelope", solve_for="requiv@primary")
    b.flip_constraint("fillout_factor", component="contact_envelope", solve_for="pot")
    b.set_value("fillout_factor", component="contact_envelope", context="component", value=f)
    b.set_value(
        "sma",
        component="binary",
        value=(units.G * 1.2 * (1 + q) * (0.4775 / 365.25) ** 2 / (4 * math.pi**2)) ** (1 / 3)
        / units.RSUN_AU,
    )
    a, y = limb_darkening(6200.0, 4.2, 0.0, "V")
    for c in ("primary", "secondary"):
        b.set_value("teff", component=c, value=6200.0)
        b.set_value("gravb_bol", component=c, value=0.32)
    b.add_dataset("lc", compute_phases=np.linspace(0, 1, 101), passband="Johnson:V")
    b.set_value_all("atm", "blackbody")
    b.set_value_all("ld_mode", "manual")
    b.set_value_all("ld_func", "power")
    b.set_value_all("ld_coeffs", list(a))
    b.set_value_all("ld_mode_bol", "manual")
    b.set_value_all("ld_func_bol", "power")
    b.set_value_all("ld_coeffs_bol", list(a))
    b.set_value_all("irrad_method", "none")
    b.run_compute(model="qsim_check")
    fl = np.asarray(b.get_value("fluxes", model="qsim_check"))
    del y
    return float(-2.5 * math.log10(fl.min() / fl.max()))


# --------------------------------------------------------------------------- rigid discs


def ring_rate_circular() -> float:
    """Disc-averaged precession rate of a narrow ring at 4 a_b around a circular binary,
    divided by the test-particle rate (3/4)(M1 M2/M^2)(a_b/r)^2 n(r)."""
    from .discs import ring_average_rate

    m1, m2, a, r = 1.0, 0.5, 1.0, 4.0
    w = ring_average_rate(m1, m2, a, r * 0.9999, r, 1.0)
    n = math.sqrt(units.G * (m1 + m2) / r**3)
    return w / (0.75 * m1 * m2 / (m1 + m2) ** 2 * (a / r) ** 2 * n)


def polar_critical_inclination(e: float = 0.5) -> float:
    """Smallest inclination [deg] at which a test-particle disc around an eccentric binary
    librates about the binary's eccentricity vector (node at 90 deg from the pericentre),
    found by integrating swingersim's rigid-disc equation (swingersim.discs.rigid_rhs)."""
    from scipy.integrate import solve_ivp

    from .discs import rigid_rhs

    rg = {"lb": np.array([0.0, 0.0, 1.0]), "eb": np.array([e, 0.0, 0.0]), "w": 1.0, "e2": e * e}

    def librates(inc_deg):
        i = math.radians(inc_deg)
        l0 = np.array([math.sin(i), 0.0, math.cos(i)])  # node perpendicular to e_b
        sol = solve_ivp(
            lambda t, l: rigid_rhs(l, rg), (0, 200.0), l0, rtol=1e-9, atol=1e-12, max_step=0.05
        )
        # librating about e_b: l.e_b never changes sign
        return bool(np.all(sol.y[0] > 0))

    lo, hi = 1.0, 89.0
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        if librates(mid):
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def polar_critical_expected(e: float = 0.5) -> float:
    return math.degrees(math.acos(math.sqrt(5 * e * e / (1 + 4 * e * e))))


def ee_cep_rigid_period() -> float:
    """Precession period [yr] of the illustrative EE Cep disc (configs/catalog/ee_cep.yaml)
    from the torque of the B star on a rigid disc (dynamics: rigid)."""
    import yaml

    from .config import Config
    from .discs import build_discs
    from .hierarchy import build_tree
    from .validation_discs import ROOT

    d = yaml.safe_load((ROOT / "configs/catalog/ee_cep.yaml").read_text())
    d["discs"]["disc"].pop("precession_period_yr", None)
    d["discs"]["disc"]["dynamics"] = "rigid"
    cfg = Config.model_validate(d)
    root, _ = build_tree(cfg)
    return abs(2.0 * math.pi / build_discs(cfg, root)[0].rate)


# --------------------------------------------------------------------------- stars


def rejuvenation_ratio() -> float:
    """Ratio of the main-sequence age after accreting 1 % of the mass to the age before, for a
    2 Msun SSE star half-way through its main sequence (fractional age kept, then the Hurley et
    al. 2000 factor M / M') - divided by the analytic value (t'_MS / t_MS)(M / M')."""
    from dataclasses import replace

    from .config import Config
    from .evolution import track
    from .masstransfer import star_response

    cfg = Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"A": ["m1", "m2"]},
            "masses": {"m1": 2.0, "m2": 1.0},
            "orbits": {"A": {"P": 10.0}},
            "stellar_model": {"kind": "sse", "age_myr": 100.0, "z": 0.02},
            "physics": {"evolution": True, "mass_transfer": {}},
            "run": {"backend": "secular", "t_max_yr": 1.0},
        }
    )
    m, dm = 2.0, 0.02
    tr = track(m, 0.02)
    age = 0.5 * tr.t_ms
    st = replace(tr.at(age), phase=1)
    _, dage, _ = star_response(cfg, st, dm, m, m, age, 0.0, m)
    age_new = age + dage
    t_ms_new = track(m + dm, 0.02).t_ms
    expected = age * (t_ms_new / tr.t_ms) * (m / (m + dm))
    return age_new / expected


def k2_rotation_factor() -> float:
    """k2 of a star rotating at V = 200 km/s (2 Msun, 2 Rsun) relative to the static value,
    against 10^(-0.87 * 2 V^2 / (3 g R)) (Claret 1999)."""
    from .config import Config
    from .tides import k2_rotating

    cfg = Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"A": ["m1", "m2"]},
            "masses": {"m1": 2.0, "m2": 1.0},
            "orbits": {"A": {"P": 10.0}},
            "physics": {"tides": True, "rotation": True, "k2_rotation": True},
            "run": {"t_max_yr": 1.0},
        }
    )
    v = 200.0 * units.SECONDS_PER_YEAR / units.AU_KM  # AU/yr
    r = 2.0 * units.RSUN_AU
    omega = v / r
    got = k2_rotating(cfg, 0, 0.005, omega, 2.0, 2.0) / 0.005
    g = units.G * 2.0 / r**2
    return got / 10 ** (-0.87 * 2.0 * v * v / (3.0 * g * r))


def pms_population_period() -> float:
    from .validation import circularization_period_mist

    return circularization_period_mist(3.0, 1.0, True, 1.0, "linear", final=True)


def m67_final_period() -> float:
    from .validation import circularization_period_mist

    return circularization_period_mist(4000.0, 1.0, True, 1.0, "linear", final=True)


def contact_jdot_time() -> float:
    """Angular-momentum loss time J / |dJ/dt| [Gyr] of a 1.0 + 0.5 Msun contact binary with
    P = 0.35 d and radii 1.0 and 0.7 Rsun (Gazeas & Stepien 2008, eq. 19)."""
    from .masstransfer import contact_jdot, orbital_j

    a = units.period_to_sma(0.35 / units.DAYS_PER_YEAR, 1.5)
    return orbital_j(1.0, 0.5, a) / abs(contact_jdot(1.0, 0.5, 1.0, 0.7, a)) / 1e9


def contact_jdot_expected() -> float:
    """The same from eq. 19 in cgs units: J = M1 M2 (G a / M)^(1/2)."""
    g, msun, rsun_cm, day = 6.674e-8, 1.98847e33, 6.957e10, 86400.0
    m1, m2, p = 1.0 * msun, 0.5 * msun, 0.35 * day
    a = (g * (m1 + m2) * p * p / (4 * math.pi**2)) ** (1 / 3)
    j = m1 * m2 * math.sqrt(g * a / (m1 + m2))
    jdot = 4.9e41 * (1.0 * 1.0 + 0.7**2 * 0.5) / 0.35  # per year
    del rsun_cm
    return j / jdot / 1e9


def algol_giant_phase_period() -> float:
    """Period [d] at the end of the BSE Algol example's transfer from its state at the end of
    the Hertzsprung gap (0.53 + 2.74 Msun, 17 d) to the final donor mass 0.42 Msun, conservative
    transfer (a M_d^2 M_a^2 = const)."""
    from .masstransfer import orbit_path

    a0 = units.period_to_sma(17.0 / units.DAYS_PER_YEAR, 0.53 + 2.74)
    md, ma, a = orbit_path(0.53, 2.74, a0, 0.0, 0.0, 0.42, n=400)
    return units.sma_to_period(float(a[-1]), float(md[-1] + ma[-1])) * units.DAYS_PER_YEAR


BENCHMARKS_V027: list[Benchmark] = [
    Benchmark(
        id="mt-bse-algol-giant-phase",
        tier="B",
        quantity="period at the end of the BSE Algol example's transfer, computed from its state "
        "at the end of the Hertzsprung gap (0.53 + 2.74 Msun, 17 d; conservative)",
        unit="d",
        expected=30.0,
        tolerance=0.03,
        relative=True,
        source="Hurley, Tout & Pols (2002) MNRAS 329, 897, Sect. 3.1",
        link="https://arxiv.org/abs/astro-ph/0201220",
        rationale="the published period is given to 1 d",
        compute=algol_giant_phase_period,
    ),
    Benchmark(
        id="roche-lobe-volume",
        tier="A",
        quantity="volume-equivalent radius of the Roche lobe (q = M2/M1 = 0.3) from the surface "
        "mesh of swingersim.roche / Eggleton's formula",
        unit="",
        expected=1.0,
        tolerance=0.01,
        relative=True,
        source="Eggleton (1983) ApJ 268, 368: the fit is accurate to 1 % for all q",
        link="https://doi.org/10.1086/160960",
        rationale="the fit's own accuracy; the mesh integration is better than 1e-3",
        compute=roche_lobe_vs_eggleton,
    ),
    Benchmark(
        id="rotating-star-critical",
        tier="A",
        quantity="equatorial / polar radius of a Roche-model star at critical rotation",
        unit="",
        expected=1.5,
        tolerance=1e-6,
        relative=True,
        source="Maeder (2009) Physics, Formation and Evolution of Rotating Stars, Sect. 4.2",
        link="https://doi.org/10.1007/978-3-540-76949-1",
        rationale="exact (bisection to 1e-12)",
        compute=rotating_critical_flattening,
    ),
    Benchmark(
        id="rotating-star-equator",
        tier="A",
        quantity="R_eq/R_pol - (1 + w^2/2) with w = Omega / Omega_K(R_eq), Roche model, "
        "Omega / Omega_crit = 0.6",
        unit="",
        expected=0.0,
        tolerance=1e-9,
        relative=False,
        source="Roche model: GM/R_eq + Omega^2 R_eq^2 / 2 = GM/R_pol (Espinosa Lara & "
        "Rieutord 2011, eq. 29)",
        link="https://doi.org/10.1051/0004-6361/201117252",
        rationale="exact identity",
        compute=rotating_equator,
    ),
    Benchmark(
        id="mesh-sphere-eclipse",
        tier="A",
        quantity="largest difference of the light curves of two spheres through both eclipses: "
        "surface meshes (swingersim.mesh) against the exact tier-A eclipses",
        unit="",
        expected=0.0,
        tolerance=1e-3,
        relative=False,
        source="exact sphere eclipses of swingersim.photometry (validated against batman, ellc)",
        link="https://doi.org/10.1051/0004-6361/201628079",
        rationale="discretisation of the mesh (64 x 128 directions; ~3e-4 obtained)",
        compute=mesh_sphere_eclipse,
    ),
    Benchmark(
        id="contact-vs-phoebe",
        tier="C",
        quantity="V amplitude of a contact binary (q = 0.42, fill-out 0.3, i = 80 deg): "
        "swingersim mesh minus PHOEBE 2.5 with the same physics",
        unit="mag",
        expected=0.0,
        tolerance=0.03,
        relative=False,
        source="PHOEBE 2.5 (Prsa et al. 2016; Conroy et al. 2020), contact_envelope, black "
        "bodies, power limb darkening, no irradiation",
        link="https://doi.org/10.3847/1538-4365/ab9ee4",
        rationale="different meshes and gravity-darkening laws (PHOEBE: bolometric beta = 0.32; "
        "swingersim: band exponent y)",
        compute=contact_amplitude_difference,
        slow=True,
    ),
    Benchmark(
        id="disc-ring-precession",
        tier="A",
        quantity="precession rate of a narrow rigid ring at 4 a_b around a circular binary / "
        "test-particle rate",
        unit="",
        expected=1.0,
        tolerance=1e-3,
        relative=True,
        source="Zanazzi & Lai (2018) MNRAS 473, 603, eq. 2",
        link="https://doi.org/10.1093/mnras/stx2375",
        rationale="the disc average of a narrow ring is the test-particle rate",
        compute=ring_rate_circular,
    ),
    Benchmark(
        id="disc-polar-critical",
        tier="A",
        quantity="critical inclination [deg] for polar libration of a disc around a binary with "
        "e = 0.5, from integrations of swingersim's rigid-disc equation",
        unit="deg",
        expected=polar_critical_expected(0.5),
        tolerance=0.2,
        relative=False,
        source="Farago & Laskar (2010) MNRAS 401, 1189; Martin & Lubow (2017) ApJL 835, L28: "
        "cos I_crit = sqrt(5 e^2 / (1 + 4 e^2))",
        link="https://doi.org/10.1111/j.1365-2966.2009.15711.x",
        rationale="bisection on finite integrations near the separatrix",
        compute=polar_critical_inclination,
    ),
    Benchmark(
        id="eecep-disc-precession",
        tier="D",
        quantity="precession period [yr] of the EE Cep disc from the torque of the B star on a "
        "rigid disc (illustrative disc of configs/catalog/ee_cep.yaml, r_out 0.8 AU)",
        unit="yr",
        expected=65.0,
        tolerance=0.15,
        relative=True,
        source="Pienkowski et al. (2020) A&A 639, A23: precession in 11-12 orbital periods "
        "(62-67 yr)",
        link="https://doi.org/10.1051/0004-6361/201937181",
        rationale="the disc size is illustrative; the rate scales as r_out^(3/2)",
        compute=ee_cep_rigid_period,
        explanation="a rigid disc of 0.8 AU precesses in ~180 yr; the observed 65 yr needs "
        "r_out ~ 1.6 AU (the rate scales as r_out^1.5) or a heavier B star - a constraint on the "
        "size of the disc",
    ),
    Benchmark(
        id="rejuvenation-hurley",
        tier="A",
        quantity="main-sequence age of a 2 Msun accretor after gaining 1 % of its mass / "
        "(t'_MS / t_MS)(M / M') t",
        unit="",
        expected=1.0,
        tolerance=2e-3,
        relative=True,
        source="Hurley, Pols & Tout (2000) MNRAS 315, 543, Sect. 7.1",
        link="https://doi.org/10.1046/j.1365-8711.2000.03426.x",
        rationale="first-order step (dM/M = 0.01)",
        compute=rejuvenation_ratio,
    ),
    Benchmark(
        id="k2-rotation",
        tier="A",
        quantity="rotational correction of k2 (2 Msun, 2 Rsun, V = 200 km/s) / 10^(-0.87 * 2 V^2 / 3 g R)",
        unit="",
        expected=1.0,
        tolerance=1e-9,
        relative=True,
        source="Claret (1999) A&A 350, 56; Claret & Gimenez (2010) A&A 519, A57, Sect. 3",
        link="https://doi.org/10.1051/0004-6361/201014008",
        rationale="identity",
        compute=k2_rotation_factor,
    ),
    Benchmark(
        id="contact-jdot",
        tier="A",
        quantity="angular-momentum loss time [Gyr] of a 1.0 + 0.5 Msun contact binary, "
        "P = 0.35 d (swingersim units vs cgs)",
        unit="Gyr",
        expected=contact_jdot_expected(),
        tolerance=1e-3,
        relative=True,
        source="Gazeas & Stepien (2008) MNRAS 390, 1577, eq. 19",
        link="https://doi.org/10.1111/j.1365-2966.2008.13844.x",
        rationale="unit conversion",
        compute=contact_jdot_time,
    ),
    Benchmark(
        id="tides-circularization-pms-population",
        tier="D",
        quantity="circularization period (final period, e < 0.01) of twin 1 Msun binaries at "
        "3 Myr (MIST pre-main sequence from 1 Myr, equilibrium tide with linear damping and "
        "inertial waves)",
        unit="d",
        expected=7.1,
        tolerance=1.2,
        relative=False,
        source="Meibom & Mathieu (2005) ApJ 620, 970, Table 1 (pre-main-sequence binaries: "
        "7.1 +- 1.2 d; quoted by Mazeh 2008)",
        link="https://arxiv.org/abs/0801.0134",
        rationale="the published uncertainty",
        compute=pms_population_period,
        slow=True,
    ),
    Benchmark(
        id="tides-circularization-m67-final",
        tier="D",
        quantity="circularization period (final period) of twin 1 Msun binaries at 4 Gyr "
        "(MIST from 1 Myr, equilibrium tide with linear damping and inertial waves)",
        unit="d",
        expected=12.1,
        tolerance=1.2,
        relative=False,
        source="Meibom & Mathieu (2005) ApJ 620, 970, Table 3 (M67: 12.1 +1.0/-1.5 d; quoted by "
        "Mazeh 2008)",
        link="https://arxiv.org/abs/0801.0134",
        rationale="the published uncertainty",
        compute=m67_final_period,
        slow=True,
        explanation="the pre-main sequence circularises up to ~6.6 d (as observed, 7.1 +- 1.2 d) "
        "but the main sequence adds only ~1 d. Barker (2022) finds that inertial waves in "
        "evolving stars with spin-down can explain the growth to ~12 d; the frequency-averaged "
        "two-layer estimate used here with MIST structure does not reproduce that (an open "
        "difference in the stellar structure and spin history, not a numerical error)",
    ),
]
