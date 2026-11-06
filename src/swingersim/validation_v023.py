"""Benchmarks added in version 0.23 (appended to swingersim.validation.BENCHMARKS)."""

from __future__ import annotations

import math
import tempfile
from functools import cache
from pathlib import Path

import numpy as np

from ._benchmark import Benchmark

# ------------------------------------------------------------------ O - C fitting (stage 1)


def gimenez_bastero_shift(n, kind, P, e, w0, dw, terms=10):
    """Apsidal-motion term of the times of minima at i = 90 deg as the series of Gimenez &
    Bastero (1995), in the form used by OCFit (Gajdos & Parimucha 2019)."""
    w = w0 + dw * (np.asarray(n, float) + np.asarray(kind) / 2.0)
    nu = math.pi / 2 - w
    b = e / (1 + math.sqrt(1 - e * e))
    s1 = np.zeros_like(w)
    s2 = np.zeros_like(w)
    for k in range(1, terms):
        term = (-b) ** k * (1 / k + math.sqrt(1 - e * e)) * np.sin(k * nu)
        s1 += term
        s2 += -term if k % 2 else term
    return np.where(np.asarray(kind) == 0, P / math.pi * s1, P / math.pi * s2)


def ocfit_gb95_max_diff_s(e: float = 0.3) -> float:
    """Largest difference [s] between swingersim's exact apsidal term and the Gimenez & Bastero
    (1995) series, P = 2.5 d, over a full apsidal cycle (constant offset removed)."""
    from .ocfit import apsidal_shift

    n = np.arange(0, 3200, 3.0)
    worst = 0.0
    for j in (0, 1):
        k = np.full(len(n), j)
        d = apsidal_shift(n, k, 2.5, e, 1.1, 2e-3) - gimenez_bastero_shift(n, k, 2.5, e, 1.1, 2e-3)
        worst = max(worst, float(np.abs(d - d.mean()).max()) * 86400.0)
    return worst


def ocfit_lite_amplitude_ratio() -> float:
    """Half the peak-to-peak light-time effect over one outer orbit, divided by the
    semi-amplitude K = (a sin i / c) sqrt(1 - e^2 cos^2 w) of Irwin (1952)."""
    from .ocfit import lite, lite_semi_amplitude_d

    t = np.linspace(0, 1000.0, 200001)
    v = lite(t, 1.7, 0.55, 2.2, 100.0, 1000.0)
    return 0.5 * float(v.max() - v.min()) / lite_semi_amplitude_d(1.7, 0.55, 2.2)


@cache
def ocfit_nbody_recovery() -> dict:
    """swingersim timings of an N-body 2+2 (outer P = 1500 d, e = 0.3, i = 60 deg), fitted with
    stage 1 (shared outer orbit): fitted / Keplerian a sin i of both pairs and fitted P3."""
    from .config import Config
    from .ocfit import fit_stage1, prepare
    from .timing import eclipse_times
    from .units import DAYS_PER_YEAR

    mA, mB = (1.2, 1.0), (1.1, 0.9)
    cfg = Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": mA[0], "m2": mA[1], "m3": mB[0], "m4": mB[1]},
            "orbits": {
                "A": {"P": 1.2, "e": 0.0, "i": 89.0, "Omega": 10.0},
                "B": {"P": 1.5, "e": 0.0, "i": 88.0, "Omega": 40.0},
                "out": {"P": 1500.0, "e": 0.3, "i": 60.0, "omega": 70.0, "Omega": 20.0},
            },
            "physics": {"gr": False},
            "run": {"t_max_yr": 13.0},
        }
    )
    ecl, _ = eclipse_times(cfg, 13.0)
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "minima.csv"
        with open(path, "w", encoding="utf-8") as f:
            f.write("pair,type,time,sigma,source,method\n")
            for e in ecl:
                t = 2459000.0 + e.t_obs_yr * DAYS_PER_YEAR
                f.write(f"{e.binary},{e.kind[0]},{t:.8f},0.00001,nbody,individual\n")
        cfg2 = Config.model_validate(
            {
                **cfg.model_dump(exclude_none=True, exclude={"fit"}),
                "fit": {"minima": str(path), "stage1": {"terms": ["lite"], "shared_outer": True}},
            }
        )
        model, data, _ = prepare(cfg2)
        res = fit_stage1(model, data)
    par = res.params()
    from .units import G

    M = sum(mA) + sum(mB)
    a_out = (G * M * (1500.0 / DAYS_PER_YEAR / (2 * math.pi)) ** 2) ** (1 / 3)
    s = math.sin(math.radians(60.0))
    return {
        "A": par["A.asini"] / (a_out * sum(mB) / M * s),
        "B": par["B.asini"] / (a_out * sum(mA) / M * s),
        "P3": par["out.P3"] / 1500.0,
        "e3": par["out.e3"],
    }


# ------------------------------------------------------------------ winds
def sse_wind_rgb_loss_1msun(eta: float = 0.1) -> float:
    """Mass [Msun] a 1 Msun SSE star (Z = 0.02) loses by the SSE winds up to helium ignition
    with the Reimers coefficient eta (MIST uses eta = 0.1 on the red giant branch)."""
    from .massloss import _sse_wind_history

    _, ms = _sse_wind_history(1.0, 0.02, eta)
    return 1.0 - float(ms[-1])


def wind_transfer_jeans_invariant() -> float:
    """Largest relative change of a (M1 + M2) of a MIST binary (25 + 18 Msun, P = 6 d) with
    winds before its mass transfer starts (Jeans mode: a M = const)."""
    from .config import Config
    from .masstransfer import histories

    cfg = Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"out": ["A", "m3"], "A": ["m1", "m2"]},
            "masses": {"m1": 25.0, "m2": 18.0, "m3": 10.0},
            "orbits": {"A": {"P": 6.0}, "out": {"P": 3000.0, "e": 0.2}},
            "stellar_model": {"kind": "mist", "age_myr": 1.0, "vvcrit": 0.4},
            "physics": {
                "gr": False,
                "evolution": True,
                "mass_transfer": {},
                "mass_loss": {"source": "track"},
            },
            "run": {"backend": "secular", "t_max_yr": 8.0e6},
        }
    )
    h = histories(cfg)["A"]
    i = int(np.searchsorted(h.t, h.t_onset)) if h.t_onset is not None else len(h.t)
    am = h.a[:i] * (h.m[0, :i] + h.m[1, :i])
    return float(np.max(np.abs(am / am[0] - 1.0)))


# ------------------------------------------------------------------ tides in nested secular
def nested_tides_hut81_ratio() -> float:
    """de/dt of the inner binary of a 3+1 (1 + 0.8 Msun, P = 4 d, e = 0.05, synchronous spins,
    lag 1000 s) in the nested secular backend over the rate of Hut (1981), eq. 10, with
    T = R^3 / (G M tau_Hut), tau_Hut = 2 lag."""
    from .backends.secular_nested import build_nested
    from .config import Config
    from .stellar import star_properties
    from .units import RSUN_AU, SECONDS_PER_YEAR, G

    cfg = Config.model_validate(
        {
            "topology": "3+1",
            "masses": {"m1": 1.0, "m2": 0.8, "m3": 0.001, "m4": 0.001},
            "orbits": {
                "inner": {"P": 4.0, "e": 0.05},
                "mid": {"P": 4.0e4, "e": 0.1, "i": 10},
                "out": {"P": 2.0e6, "e": 0.1, "i": 5},
            },
            "stars": {
                "m1": {"lag_s": 1000.0, "P_rot": "sync"},
                "m2": {"lag_s": 1000.0, "P_rot": "sync"},
            },
            "physics": {"gr": False, "tides": True, "rotation": True},
            "run": {"backend": "secular", "t_max_yr": 1.0},
        }
    )
    sys_, y0, _root, *_ = build_nested(cfg)
    q = sys_.diss[0]
    p = y0[6 * q + 3 : 6 * q + 6]
    e = float(np.linalg.norm(p))
    dedt = float(sys_.rhs(0.0, y0)[6 * q + 3 : 6 * q + 6] @ p / e)
    a = float(sys_.a0[q])
    m = (1.0, 0.8)
    n = math.sqrt(G * sum(m) / a**3)
    e2 = e * e
    f3 = 1 + 15 / 4 * e2 + 15 / 8 * e2**2 + 5 / 64 * e2**3
    f4 = 1 + 1.5 * e2 + e2**2 / 8
    hut = 0.0
    for i in (0, 1):
        pr = star_properties(cfg, f"m{i + 1}", m[i])
        R = pr.radius_rsun * RSUN_AU
        tau = 2.0 * 1000.0 / SECONDS_PER_YEAR
        T = R**3 / (G * m[i] * tau)
        qq = m[1 - i] / m[i]
        w = sys_.spins_state(y0)[i] / n
        hut += (
            -27
            * (pr.apsidal_k2 / T)
            * qq
            * (1 + qq)
            * (R / a) ** 8
            * e
            / (1 - e2) ** 6.5
            * (f3 - 11 / 18 * (1 - e2) ** 1.5 * f4 * w)
        )
    return dedt / hut


# ------------------------------------------------------------------ environment in N-body
def galactic_tide_emax_nbody(inc_deg: float = 60.0) -> float:
    """Largest eccentricity of a 2 x 10^4 AU binary (1 + 1 Msun, e0 = 0.01) inclined by inc_deg
    to the Galactic disc, N-body backend with physics.galactic_tide, over 3 Gyr."""
    import warnings

    from .backends.nbody import run_nbody
    from .config import Config

    cfg = Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"out": ["m1", "m2"]},
            "masses": {"m1": 1.0, "m2": 1.0},
            "orbits": {"out": {"a": 2.0e4, "e": 0.01, "i": inc_deg}},
            "physics": {"gr": False, "galactic_tide": {}},
            "run": {"backend": "nbody", "t_max_yr": 3.0e9, "n_output": 3001},
        }
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return float(run_nbody(cfg).extremes["out"]["e_max"])


def _e_and_j(e, inc, Om, om):
    ci, si, cO, sO, co, so = (
        math.cos(inc),
        math.sin(inc),
        math.cos(Om),
        math.sin(Om),
        math.cos(om),
        math.sin(om),
    )
    P = np.array([cO * co - sO * so * ci, sO * co + cO * so * ci, so * si])
    k = np.array([sO * si, -cO * si, ci])
    return e * P, math.sqrt(1 - e * e) * k


def flyby_nbody_backend_error(q_au: float = 20.0) -> float:
    """|De(swingersim N-body run with one passing star) - De(first order)| / |De(first order)|:
    a parabolic passage of 1 Msun at Q = 20 AU past a 1 + 1 Msun binary (a = 1 AU, e = 0.5),
    the passing star integrated as a particle by the N-body backend (physics.encounters)."""
    import warnings

    from . import environment as envmod
    from .backends.nbody import run_nbody
    from .config import Config

    enc = envmod.Encounter(t_yr=400.0, q_au=q_au, ecc=1.0, rotation=np.eye(3), m_per=1.0)
    cfg = Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"out": ["m1", "m2"]},
            "masses": {"m1": 1.0, "m2": 1.0},
            "orbits": {"out": {"a": 1.0, "e": 0.5, "i": 50.0, "omega": 30.0, "Omega": 40.0}},
            "physics": {"gr": False, "encounters": {"n_star_pc3": 0.1}},
            "run": {"backend": "nbody", "t_max_yr": 800.0, "n_output": 3},
        }
    )
    orig = envmod.encounter_sequence
    envmod.encounter_sequence = lambda *a, **k: [enc]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = run_nbody(cfg)
    finally:
        envmod.encounter_sequence = orig
    o = r.orbits["out"]
    e0, j0 = _e_and_j(o["e"][0], o["inc"][0], o["Omega"][0], o["omega"][0])
    e1, _ = _e_and_j(o["e"][-1], o["inc"][-1], o["Omega"][-1], o["omega"][-1])
    de_fo, _ = envmod.first_order_kick(e0, j0, float(o["a"][0]), 2.0, 1.0, q_au, 1.0)
    return float(np.linalg.norm((e1 - e0) - de_fo) / np.linalg.norm(de_fo))


# ------------------------------------------------------------------ gravitational waves
def gw_hulse_taylor_pdot() -> float:
    """dP/dt of PSR B1913+16 from the Peters (1964) rates with the masses and orbit of
    Weisberg & Huang (2016)."""
    from .gw import pdot

    return pdot(1.438, 1.390, 0.322997448918, 0.6171340)


def _gw_triple(backend: str, t_max: float, a: float = 0.003, e: float = 0.3, n_out: int = 11):
    from .config import Config

    return Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"out": ["A", "m3"], "A": ["m1", "m2"]},
            "masses": {"m1": 1.4, "m2": 1.4, "m3": 1.0},
            "orbits": {"A": {"a": a, "e": e}, "out": {"a": 100.0, "e": 0.1, "i": 10}},
            "stars": {"m1": {"R": 1e-5}, "m2": {"R": 1e-5}},
            "physics": {"gr": False, "gw": True},
            "run": {"backend": backend, "t_max_yr": t_max, "n_output": n_out},
        }
    )


def gw_secular_circular_decay() -> float:
    """a(t) / a0 (1 - t/tau)^(1/4) at t = 0.9 tau of a circular 1.4 + 1.4 Msun binary in a
    triple, secular backend with physics.gw (tau = Peters' merger time)."""
    import warnings

    from .backends.secular import run_secular
    from .gw import merger_time_circular_yr

    tau = merger_time_circular_yr(1.4, 1.4, 0.003)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = run_secular(_gw_triple("secular", 0.9 * tau, e=0.0, n_out=21))
    return float(r.orbits["A"]["a"][-1] / (0.003 * 0.1**0.25))


def gw_nbody_da_ratio() -> float:
    """Change of a over 10 yr of a 1.4 + 1.4 Msun binary (a = 0.003 AU, e = 0.3) in a triple,
    N-body backend with physics.gw, over the Peters (1964) rate times 10 yr."""
    import warnings

    from .backends.nbody import run_nbody
    from .gw import peters_rates

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = run_nbody(_gw_triple("nbody", 10.0))
    da, _ = peters_rates(1.4, 1.4, 0.003, 0.3)
    a = r.orbits["A"]["a"]
    return float((a[-1] - a[0]) / (da * 10.0))


# ------------------------------------------------------------------ tidal capture
def tidal_capture_radius(n: float) -> float:
    """Largest periastron distance / R for tidal capture of two equal stars (0.6 Msun,
    0.65 Rsun, polytrope n) at a relative velocity of 10 km/s at infinity, with the fits of
    Portegies Zwart & Meinen (1993)."""
    from .tidalcapture import capture_radius
    from .units import AU_KM, RSUN_AU, SECONDS_PER_YEAR

    r = 0.65 * RSUN_AU
    v = 10.0 * SECONDS_PER_YEAR / AU_KM
    return capture_radius(0.6, r, n, 0.6, r, n, v) / r


OCFIT = "Gajdos & Parimucha (2019) OEJV 197, 71 (OCFit); Gimenez & Bastero (1995) Ap&SS 226, 99"
OCFIT_LINK = "https://github.com/pavolgaj/OCFit"
IRWIN = "Irwin (1952) ApJ 116, 211, eqs. 1-3"
IRWIN_LINK = "https://doi.org/10.1086/145604"

BENCHMARKS_V023: list[Benchmark] = [
    Benchmark(
        id="ocfit-gimenez-bastero",
        tier="A",
        quantity="apsidal term of the times of minima, e = 0.3: largest difference between the "
        "exact conjunction times and the series of Gimenez & Bastero (1995) used by OCFit",
        unit="s",
        expected=0.0,
        tolerance=0.01,
        relative=False,
        source=OCFIT,
        link=OCFIT_LINK,
        rationale="the series is truncated after the ninth power of e/(1+sqrt(1-e^2)); its "
        "remainder at e = 0.3 is below 10^-3 s, so 0.01 s checks the sign conventions and the "
        "phase of primary and secondary minima",
        compute=ocfit_gb95_max_diff_s,
    ),
    Benchmark(
        id="ocfit-lite-amplitude",
        tier="A",
        quantity="half peak-to-peak light-time effect over one orbit / K of Irwin (1952), "
        "e3 = 0.55",
        unit="",
        expected=1.0,
        tolerance=1e-6,
        relative=True,
        source=IRWIN,
        link=IRWIN_LINK,
        rationale="analytic identity; the sampling of the orbit limits the agreement to 10^-7",
        compute=ocfit_lite_amplitude_ratio,
    ),
    Benchmark(
        id="ocfit-nbody-asini-A",
        tier="A",
        quantity="stage-1 fit of N-body eclipse timings of a 2+2 (outer P = 1500 d, e = 0.3, "
        "i = 60 deg): fitted / Keplerian a sin i of pair A",
        unit="",
        expected=1.0,
        tolerance=0.01,
        relative=True,
        source="Kepler's third law and the configured outer orbit; Irwin (1952)",
        link=IRWIN_LINK,
        rationale="the dynamical perturbations of the pairs by each other change the timings by "
        "about 0.5 % of the light-time amplitude for P_out/P_in ~ 1000 (dynamical delay ~ (3/8 pi) "
        "P_in^2/P_out)",
        compute=lambda: ocfit_nbody_recovery()["A"],
        slow=True,
    ),
    Benchmark(
        id="ocfit-nbody-asini-B",
        tier="A",
        quantity="the same for pair B (anti-phase, shared outer orbit)",
        unit="",
        expected=1.0,
        tolerance=0.01,
        relative=True,
        source="Kepler's third law and the configured outer orbit; Irwin (1952)",
        link=IRWIN_LINK,
        rationale="as for pair A; the ratio of the two gives the mass ratio of the pairs",
        compute=lambda: ocfit_nbody_recovery()["B"],
        slow=True,
    ),
    Benchmark(
        id="ocfit-nbody-P3",
        tier="A",
        quantity="the same fit: fitted / configured outer period",
        unit="",
        expected=1.0,
        tolerance=0.005,
        relative=True,
        source="configured outer orbit",
        link=IRWIN_LINK,
        rationale="three outer orbits are covered; the osculating period differs from the mean "
        "one by the perturbations of the pairs, < 0.5 %",
        compute=lambda: ocfit_nbody_recovery()["P3"],
        slow=True,
    ),
    Benchmark(
        id="sse-wind-rgb-1msun",
        tier="C",
        quantity="mass lost by a 1 Msun star (Z = 0.02) on the red giant branch with the Reimers "
        "wind, eta = 0.1: SSE track with the SSE winds against the MIST v1.2 (MESA) model at "
        "the tip of the giant branch",
        unit="Msun",
        expected=1.0 - 0.95842,
        tolerance=0.1,
        relative=True,
        source="MIST v1.2 track of 1 Msun, [Fe/H] = 0, vvcrit 0 (Choi et al. 2016, Reimers "
        "eta = 0.1 on the RGB); Hurley et al. (2000) Sect. 7.1",
        link="https://doi.org/10.3847/0004-637X/823/2/102",
        rationale="the same Reimers law integrated along two different tracks (SSE fits against "
        "MESA): luminosity and radius on the upper giant branch differ by several per cent",
        compute=sse_wind_rgb_loss_1msun,
    ),
    Benchmark(
        id="wind-transfer-jeans",
        tier="A",
        quantity="25 + 18 Msun MIST binary with winds and mass transfer: largest relative change "
        "of a (M1 + M2) before the transfer starts",
        unit="",
        expected=0.0,
        tolerance=1e-8,
        relative=False,
        source="Jeans (1924); Hadjidemetriou (1963): isotropic slow mass loss keeps a M constant",
        link="https://doi.org/10.1016/0019-1035(63)90072-6",
        rationale="analytic invariant; limited by the ODE tolerance (1e-8)",
        compute=wind_transfer_jeans_invariant,
    ),
    Benchmark(
        id="secular-nested-tides-hut81",
        tier="A",
        quantity="3+1 in the nested secular backend: tidal de/dt of the inner binary (P = 4 d, "
        "e = 0.05, synchronous) / Hut (1981) eq. 10",
        unit="",
        expected=1.0,
        tolerance=1e-6,
        relative=True,
        source="Hut (1981) A&A 99, 126, eqs. 9-11 (tau_Hut = 2 x the lag of Leconte et al. 2010)",
        link="https://ui.adsabs.harvard.edu/abs/1981A%26A....99..126H",
        rationale="analytic: the equations are the same constant-time-lag model; the agreement "
        "checks the wiring of the dissipation into the nested hierarchy",
        compute=nested_tides_hut81_ratio,
    ),
    Benchmark(
        id="galactic-tide-emax-60-nbody",
        tier="A",
        quantity="largest eccentricity of a 2e4 AU binary inclined by 60 deg to the disc, N-body "
        "backend with the Galactic tide (from e = 0.01)",
        unit="",
        expected=math.sqrt(1.0 - 1.25 * math.cos(math.radians(60.0)) ** 2),
        tolerance=0.01,
        relative=True,
        source="Heisler & Tremaine (1986); Grishin & Perets (2022) MNRAS 512, 4993, eq. 6",
        link="https://doi.org/10.1093/mnras/stac706",
        rationale="the orbit-averaged e_max holds to O(P/tau_GT) ~ 10^-3; the tidal kicks are "
        "applied every P/400",
        compute=galactic_tide_emax_nbody,
        slow=True,
    ),
    Benchmark(
        id="flyby-nbody-backend",
        tier="B",
        quantity="parabolic passage at Q = 20 a integrated by the N-body backend "
        "(physics.encounters): relative difference of the eccentricity kick from the first-order "
        "result of Hamers & Samsing (2019)",
        unit="",
        expected=0.0,
        tolerance=0.05,
        relative=False,
        source="Hamers & Samsing (2019) MNRAS 487, 5630, eq. 9 (first order in the perturbation)",
        link="https://doi.org/10.1093/mnras/stz1631",
        rationale="second-order terms are of relative size (a/Q)^(3/2) ~ 1 % at Q = 20 a; the "
        "passage starts and ends at 20 Q",
        compute=flyby_nbody_backend_error,
    ),
    Benchmark(
        id="gw-hulse-taylor-pdot",
        tier="D",
        quantity="orbital decay dP/dt of PSR B1913+16 by gravitational waves (masses 1.438 and "
        "1.390 Msun, P = 0.323 d, e = 0.617)",
        unit="",
        expected=-2.40263e-12,
        tolerance=1e-3,
        relative=True,
        source="Weisberg & Huang (2016) ApJ 829, 55, Table 3 (general-relativity prediction); "
        "Peters (1964)",
        link="https://doi.org/10.3847/0004-637X/829/1/55",
        rationale="the masses are rounded to 1e-3 Msun; 0.1 % covers that and the constants",
        compute=gw_hulse_taylor_pdot,
    ),
    Benchmark(
        id="gw-secular-circular-decay",
        tier="A",
        quantity="circular 1.4 + 1.4 Msun binary in a triple, secular backend with physics.gw: "
        "a / a0 (1 - t/tau)^(1/4) at 90 % of the merger time",
        unit="",
        expected=1.0,
        tolerance=1e-4,
        relative=True,
        source="Peters (1964) Phys. Rev. 136, B1224, eqs. 5.6 and 5.10",
        link="https://doi.org/10.1103/PhysRev.136.B1224",
        rationale="analytic solution; the integration tolerance and the third star (100 AU, "
        "negligible) limit the agreement",
        compute=gw_secular_circular_decay,
    ),
    Benchmark(
        id="gw-nbody-decay",
        tier="A",
        quantity="the same binary with e = 0.3, N-body backend: change of a over 10 yr / Peters "
        "rate x 10 yr",
        unit="",
        expected=1.0,
        tolerance=0.01,
        relative=True,
        source="Peters (1964) Phys. Rev. 136, B1224, eq. 5.6",
        link="https://doi.org/10.1103/PhysRev.136.B1224",
        rationale="the emission is applied once per output step; the osculating a of the N-body "
        "binary scatters by the tertiary's perturbation at the 1e-4 level of the change",
        compute=gw_nbody_da_ratio,
        slow=True,
    ),
    Benchmark(
        id="tidal-capture-radius-n3",
        tier="B",
        quantity="largest periastron distance for tidal capture / R, two n = 3 stars (0.6 Msun, "
        "0.65 Rsun) at 10 km/s",
        unit="",
        expected=2.72,
        tolerance=0.05,
        relative=True,
        source="Ray, Kembhavi & Antia (1987) A&A 184, 164, as quoted by Portegies Zwart & Meinen "
        "(1993) A&A 280, 174, Sect. 2",
        link="https://ui.adsabs.harvard.edu/abs/1993A%26A...280..174P",
        rationale="the fits reproduce the T_l of Lee & Ostriker (1986) to a few per cent in "
        "log T; 5 % in the radius",
        compute=lambda: tidal_capture_radius(3.0),
        explanation="the fits of Portegies Zwart & Meinen give 2.37 R (n = 3) and 3.33 R (n = 1.5) "
        "instead of the quoted 2.72 and 3.84 R; the energy at the quoted radii is 5 and 13 times "
        "below the orbital energy, so the quoted radii belong to other conventions of Ray et "
        "al. (velocity, structure); the same fits agree with McMillan et al. (1987) "
        "(tidal-capture-radius-n15-mcmillan)",
    ),
    Benchmark(
        id="tidal-capture-radius-n15",
        tier="B",
        quantity="the same for two n = 1.5 stars",
        unit="",
        expected=3.84,
        tolerance=0.05,
        relative=True,
        source="Ray, Kembhavi & Antia (1987) A&A 184, 164, as quoted by Portegies Zwart & Meinen "
        "(1993) A&A 280, 174, Sect. 2",
        link="https://ui.adsabs.harvard.edu/abs/1993A%26A...280..174P",
        rationale="as for n = 3",
        compute=lambda: tidal_capture_radius(1.5),
        explanation="the fits of Portegies Zwart & Meinen (1993) in swingersim are those of NBODY6 "
        "(efac2.f, efac3.f); for n = 1.5 they give 3.33 R, as McMillan et al. (1987) quoted by "
        "Mardling (1995): ~3.4 R (tidal-capture-radius-n15-mcmillan); the 3.84 R attributed to "
        "Ray et al. (1987) belong to their different conventions (tidal-capture-radius-n3)",
    ),
    Benchmark(
        id="tidal-capture-radius-n15-mcmillan",
        tier="B",
        quantity="largest periastron distance for tidal capture / R, two equal n = 1.5 stars at "
        "10 km/s",
        unit="",
        expected=3.4,
        tolerance=0.05,
        relative=True,
        source="McMillan, McDermott & Taam (1987) ApJ 318, 261, as quoted by Mardling (1995) "
        "ApJ 450, 722 ('capture up to ~3.4 R')",
        link="https://arxiv.org/abs/astro-ph/9403047",
        rationale="the quoted value is rounded to 0.1 R",
        compute=lambda: tidal_capture_radius(1.5),
    ),
]
