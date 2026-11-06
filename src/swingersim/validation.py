"""Validation against published results.

Every benchmark reproduces a result taken from the literature and compares it with an
independent expected value. Internal consistency checks (round trips, energy conservation)
live in the normal test suite; they are necessary but not sufficient.

Run all benchmarks:  swingersim validate
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass
from functools import lru_cache, partial
from pathlib import Path

import numpy as np
import rebound
import reboundx

from . import units
from ._benchmark import ROOT, Benchmark
from .backends.nbody import run_nbody
from .config import Config
from .hierarchy import build_tree
from .kepler import state_to_elements


@dataclass
class Outcome:
    benchmark: Benchmark
    obtained: float
    deviation: float
    passed: bool
    seconds: float
    skipped: str = ""  # reason: an optional package is not installed
    error: str = ""  # the computation raised: counts as a failure


# --------------------------------------------------------------------------- helpers


def _hierarchical_triple(i_in_deg: float, e_in: float, t_max: float, n_output: int) -> Config:
    """Test-particle ZLK setup as a 3+1 configuration with a negligible, very distant 4th body.

    inner: m1 = 1 Msun with a test particle (1e-6 Msun) at a = 1 AU
    mid:   perturber m3 = 1 Msun on a circular orbit at a = 20 AU (no octupole term)
    """
    return Config.model_validate(
        {
            "topology": "3+1",
            "masses": {"m1": 1.0, "m2": 1e-6, "m3": 1.0, "m4": 1e-9},
            "orbits": {
                "inner": {"a": 1.0, "e": e_in, "i": i_in_deg, "Omega": 0.0, "omega": 0.0},
                "mid": {"a": 20.0, "e": 0.0, "i": 0.0},
                "out": {"a": 2000.0, "e": 0.0, "i": 0.0},
            },
            "stars": {k: {"R": 0.001} for k in ("m1", "m2", "m3", "m4")},
            "physics": {"gr": False},
            "run": {"t_max_yr": t_max, "n_output": n_output},
        }
    )


@lru_cache(maxsize=8)
def _zlk_result(i_in_deg: float, e_in: float, t_max: float):
    return run_nbody(_hierarchical_triple(i_in_deg, e_in, t_max, 20_000))


def _zlk_run(i_in_deg: float, e_in: float, t_max: float):
    res = _zlk_result(i_in_deg, e_in, t_max)
    return res.orbits["inner"]["e"], np.degrees(res.mutual["mid-inner"])


def _zlk_emax_formula(i0_deg: float) -> float:
    return float(np.sqrt(1.0 - 5.0 / 3.0 * np.cos(np.radians(i0_deg)) ** 2))


# --------------------------------------------------------------------------- computations


def mercury_perihelion(gr: bool = True) -> float:
    """1PN perihelion advance of Mercury around the Sun [arcsec per century]."""
    sim = rebound.Simulation()
    sim.G = units.G
    sim.integrator = "ias15"
    m_merc = 1.6601e-7
    sim.add(m=1.0)
    sim.add(m=m_merc, a=0.38709893, e=0.20563069)
    sim.move_to_com()
    rebx = reboundx.Extras(sim)  # kept alive for the whole integration
    if gr:
        force = rebx.load_force("gr_full")
        rebx.add_force(force)
        force.params["c"] = units.C
    mu = units.G * (1.0 + m_merc)
    times = np.linspace(0.0, 1000.0, 201)
    omega = []
    for t in times:
        sim.integrate(t)
        p0, p1 = sim.particles[0], sim.particles[1]
        r = np.array([p1.x - p0.x, p1.y - p0.y, p1.z - p0.z])
        v = np.array([p1.vx - p0.vx, p1.vy - p0.vy, p1.vz - p0.vz])
        omega.append(state_to_elements(mu, r, v).omega)
    rate = np.polyfit(times, np.unwrap(omega), 1)[0]  # rad / yr
    return float(np.degrees(rate) * 3600.0 * 100.0)


def zlk_emax_65() -> float:
    e, _ = _zlk_run(65.0, 0.01, 20_000.0)
    return float(e.max())


def zlk_emax_80() -> float:
    e, _ = _zlk_run(80.0, 0.01, 8_000.0)
    return float(e.max())


def zlk_imin_65() -> float:
    _, i = _zlk_run(65.0, 0.01, 20_000.0)
    return float(i.min())


def zlk_below_critical() -> float:
    e, _ = _zlk_run(30.0, 0.001, 6_000.0)
    return float(e.max())


def zlk_kozai_constant() -> float:
    """Relative secular drift of sqrt(1 - e^2) cos(i) over a ZLK run (i0 = 65 deg)."""
    e, i = _zlk_run(65.0, 0.01, 20_000.0)
    jz = np.sqrt(1.0 - e**2) * np.cos(np.radians(i))
    n = len(jz) // 4
    return float(abs(jz[-n:].mean() - jz[:n].mean()) / abs(jz[0]))


def _ma_critical(e_out: float, q_out: float) -> float:
    return 2.8 * (1.0 + q_out) ** 0.4 * (1.0 + e_out) ** 0.4 / (1.0 - e_out) ** 1.2


def _stability_run(ratio: float) -> float:
    """Coplanar equal-mass 2+2 with a_out / a_in = ratio x Mardling-Aarseth critical value.

    Returns 1.0 if the hierarchy changed or stars collided within 100 outer periods, else 0.0.
    """
    e_out, p_in_days = 0.3, 50.0
    a_in = units.period_to_sma(units.days_to_years(p_in_days), 2.0)
    a_out = ratio * _ma_critical(e_out, 1.0) * a_in
    t_max = 100.0 * units.sma_to_period(a_out, 4.0)
    cfg = Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 1.0, "m2": 1.0, "m3": 1.0, "m4": 1.0},
            "orbits": {
                "A": {"P": p_in_days, "e": 0.1},
                "B": {"P": p_in_days, "e": 0.1, "omega": 90.0, "M0": 120.0},
                "out": {"a": a_out, "e": e_out, "M0": 180.0},
            },
            "physics": {"gr": False},
            "run": {"t_max_yr": t_max, "n_output": 2000},
        }
    )
    res = run_nbody(cfg)
    return 1.0 if (res.transitions or res.event == "collision") else 0.0


def stability_inside() -> float:
    return _stability_run(0.5)


def stability_outside() -> float:
    return _stability_run(2.5)


# --- Naoz et al. (2013) MNRAS 431, 2155 -------------------------------------------------
MJUP = 9.5458e-4  # Jupiter mass [Msun]


def _naoz_config(m, a, e, i_tot, g1, t_max, n_out) -> Config:
    return Config.model_validate(
        {
            "topology": "3+1",
            "masses": {"m1": m[0], "m2": m[1], "m3": m[2], "m4": 1e-9},
            "orbits": {
                "inner": {"a": a[0], "e": e[0], "i": i_tot, "Omega": 0.0, "omega": g1},
                "mid": {"a": a[1], "e": e[1], "omega": 0.0},
                "out": {"a": 1e6},
            },
            "stars": {k: {"R": 0.001} for k in ("m1", "m2", "m3", "m4")},
            "physics": {"gr": False},
            "run": {"t_max_yr": t_max, "n_output": n_out},
        }
    )


def naoz_fig4_a2() -> float:
    """Outer semi-major axis that gives the paper's angular-momentum ratio.

    The caption states i1 = 6.75 deg and i2 = 63.25 deg with respect to the total angular
    momentum, i.e. G1/G2 = sin(i2)/sin(i1) = 7.60. With the quoted a2 = 50 AU the standard
    definitions give G1/G2 = 6.79 (i1 = 7.5 deg); the published curve corresponds to 7.60.
    """
    m1, m2, m3, a1, e1 = 1.4, 0.3, 0.01, 5.0, 0.5
    ratio = np.sin(np.radians(63.25)) / np.sin(np.radians(6.75))
    g1 = m1 * m2 * np.sqrt(units.G * a1 / (m1 + m2)) * np.sqrt(1 - e1**2)
    g2 = g1 / ratio
    return float((g2 / ((m1 + m2) * m3)) ** 2 * (m1 + m2 + m3) / units.G)


@lru_cache(maxsize=2)
def _naoz_fig4_run():
    cfg = _naoz_config((1.4, 0.3, 0.01), (5.0, naoz_fig4_a2()), (0.5, 0.0), 70.0, 120.0, 5e5, 4000)
    return np.degrees(run_nbody(cfg).mutual["mid-inner"])


def naoz_fig4_imin() -> float:
    return float(_naoz_fig4_run().min())


def naoz_fig4_imax() -> float:
    return float(_naoz_fig4_run().max())


@lru_cache(maxsize=1)
def _naoz_fig3_run():
    cfg = _naoz_config((1.0, MJUP, 40 * MJUP), (6.0, 100.0), (0.001, 0.6), 65.0, 45.0, 5e6, 5000)
    res = run_nbody(cfg)
    return res.times, np.degrees(res.mutual["mid-inner"]), res.orbits["inner"]["e"]


def naoz_fig3_flip_time() -> float:
    t, i, _ = _naoz_fig3_run()
    up = np.flatnonzero((i[:-1] < 90.0) & (i[1:] >= 90.0))
    return float(t[up[0]]) if len(up) else float("nan")


def naoz_fig3_imax() -> float:
    _, i, _ = _naoz_fig3_run()
    return float(i.max())


# --- quadrupole secular theory, MEGNO, Brouwer's law ---------------------------------------


def quad_tp_emax() -> float:
    from .secular.quadrupole import predict

    return predict(1.0, 1e-6, 1.0, 1.0, 20.0, 0.01, 0.0, 65.0, 0.0).e_max


def quad_fig4_imin() -> float:
    from .secular.quadrupole import predict

    return predict(1.4, 0.3, 0.01, 5.0, naoz_fig4_a2(), 0.5, 0.0, 70.0, 120.0).i_min_deg


def quad_fig4_imax() -> float:
    from .secular.quadrupole import predict

    return predict(1.4, 0.3, 0.01, 5.0, naoz_fig4_a2(), 0.5, 0.0, 70.0, 120.0).i_max_deg


def _zlk_run_full():
    return _zlk_result(65.0, 0.01, 20_000.0)


def zlk_period_nbody_vs_quadrupole() -> float:
    """Relative difference of the N-body ZLK period from the quadrupole secular period."""
    res = _zlk_run_full()
    d = next(x for x in res.diagnostics if x.inner == "inner")
    return float(d.zlk_period_yr / d.quad_period_yr - 1.0)


def _megno_config(stable: bool) -> Config:
    if stable:
        cfg = _hierarchical_triple(20.0, 0.01, 2_000.0, 200)
    else:
        cfg = Config.model_validate(
            {
                "topology": "2+2",
                "masses": {"m1": 1.0, "m2": 1.0, "m3": 1.0, "m4": 1.0},
                "orbits": {
                    "A": {"P": 50.0, "e": 0.1, "i": 30.0},
                    "B": {"P": 50.0, "e": 0.1, "i": 60.0, "Omega": 90.0, "M0": 50.0},
                    "out": {"P": 400.0, "e": 0.4},
                },
                "physics": {"gr": False},
                "run": {"t_max_yr": 30.0, "n_output": 300, "collisions": False},
            }
        )
    cfg.run.megno = True
    return cfg


def megno_regular() -> float:
    return float(run_nbody(_megno_config(True)).megno[-1])


def megno_chaotic() -> float:
    m = run_nbody(_megno_config(False)).megno[-1]
    return 1.0 if m > 4.0 else 0.0


def brouwer_exponent() -> float:
    """Growth exponent of the IAS15 energy error, averaged over 6 binaries (3e4 orbits)."""
    errs = []
    for m0 in np.linspace(0.0, 5.0, 6):
        sim = rebound.Simulation()
        sim.G = units.G
        sim.integrator = "ias15"
        sim.add(m=1.0)
        sim.add(m=0.5, a=1.0, e=0.5, M=m0)
        sim.move_to_com()
        period, e0 = sim.particles[1].P, sim.energy()
        norb = np.logspace(1, np.log10(3e4), 30)
        err = []
        for n in norb:
            sim.integrate(n * period)
            err.append(abs(sim.energy() / e0 - 1.0))
        errs.append(err)
    rms = np.sqrt(np.mean(np.square(errs), axis=0))
    return float(np.polyfit(np.log10(norb[5:]), np.log10(rms[5:]), 1)[0])


# --- Vokrouhlicky (2016) MNRAS 461, 3964: V994 Her-like 2+2 system ---------------------------


def v994_config(t_max: float = 1.0, n_out: int = 2) -> Config:
    """Section 3.1 / Fig. 2 setup: masses and periods of V994 Her, e_A = e_B = 0, e = 0.76."""
    return Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 3.01, "m2": 2.58, "m3": 1.84, "m4": 1.93},
            "orbits": {
                "A": {"P": 2.08327, "e": 0.0, "i": 84.0, "Omega": 10.0},
                "B": {"P": 1.42004, "e": 0.0, "i": 85.0, "Omega": 15.0},
                "out": {"P": 1062.8775, "e": 0.76, "i": 90.0, "Omega": 0.0, "omega": 236.0},
            },
            "physics": {"gr": False},
            "run": {"t_max_yr": t_max, "n_output": n_out, "extremes_every_orbits": 0},
        }
    )


@lru_cache(maxsize=1)
def _v994_initial():
    res = run_nbody(v994_config())
    return {k: float(np.degrees(v[0])) for k, v in res.mutual.items()}


def v994_ja() -> float:
    return _v994_initial()["out-A"]


def v994_jb() -> float:
    return _v994_initial()["out-B"]


def nodal_precession_period() -> float:
    """Period [yr] of the nodal precession of a circular inner orbit tilted 20 deg."""
    res = run_nbody(_hierarchical_triple(20.0, 0.001, 3_000.0, 3_000))
    rate = np.polyfit(res.times, np.unwrap(res.orbits["inner"]["Omega"]), 1)[0]
    return float(2.0 * np.pi / abs(rate))


def _nodal_precession_formula() -> float:
    """Quadrupole nodal precession period of a circular inner orbit (a = 1 AU, 1 Msun + test
    mass) tilted 20 deg to a 1-Msun perturber at 20 AU (N-body and secular benchmarks)."""
    p_in = units.sma_to_period(1.0, 1.0 + 1e-6)
    p_out = units.sma_to_period(20.0, 2.0 + 1e-6)
    m_tot, m_pert = 2.0 + 1e-6, 1.0
    return 4.0 / 3.0 * p_out**2 / p_in * m_tot / m_pert / np.cos(np.radians(20.0))


# --- F2: stellar structure, tides, rotation, Roche lobes -------------------------------------


def roche_volume_radius(q: float, n_theta: int = 40, n_phi: int = 40) -> float:
    """Exact volume-equivalent Roche-lobe radius / separation of the star with mass ratio q.

    Direct integration of the volume inside the critical equipotential through L1 in the
    corotating frame (circular orbit, synchronous rotation, point-mass potential); the lobe is
    symmetric about the orbital plane, so half of it is integrated.
    """
    from scipy.optimize import brentq

    m1, m2 = q / (1.0 + q), 1.0 / (1.0 + q)
    xcm = m2

    def phi(x, y, z):
        r1 = np.sqrt(x * x + y * y + z * z)
        r2 = np.sqrt((x - 1.0) ** 2 + y * y + z * z)
        return -m1 / r1 - m2 / r2 - 0.5 * ((x - xcm) ** 2 + y * y)

    xl1 = brentq(lambda x: m1 / x**2 - m2 / (1.0 - x) ** 2 - (x - xcm), 1e-6, 1.0 - 1e-6)
    phi_l1 = phi(xl1, 0.0, 0.0)
    mu, w_mu = np.polynomial.legendre.leggauss(n_theta)
    ph, w_ph = np.polynomial.legendre.leggauss(n_phi)
    ph, w_ph = (ph + 1.0) * np.pi / 2.0, w_ph * np.pi / 2.0
    rs = np.linspace(1e-4, 1.5, 400)
    volume = 0.0
    for c, wc in zip(mu, w_mu):
        sn = np.sqrt(1.0 - c * c)
        for p_, wp in zip(ph, w_ph):
            d = np.array([c, sn * np.cos(p_), sn * np.sin(p_)])
            vals = phi(*(rs[:, None] * d).T) - phi_l1
            k = int(np.argmax(vals > 0))
            r = brentq(lambda r, d=d: phi(*(r * d)) - phi_l1, rs[k - 1], rs[k]) if k else xl1
            volume += wc * wp * r**3 / 3.0
    return float((3.0 * 2.0 * volume / (4.0 * np.pi)) ** (1.0 / 3.0))


def eggleton_q01() -> float:
    from .stellar import roche_lobe_eggleton

    return roche_lobe_eggleton(0.1)


def eggleton_q1() -> float:
    from .stellar import roche_lobe_eggleton

    return roche_lobe_eggleton(1.0)


def polytrope_k2(n: float) -> float:
    from .stellar import polytrope_constants

    return polytrope_constants(n)[0]


def _tidal_binary_config(
    e: float,
    t_max: float,
    n_out: int,
    rotation: bool = False,
    lag_s: float = 0.0,
    spin_factor: float | str = "sync",
) -> Config:
    """2+2 with one tidally active binary A (1.0 + 0.8 Msun, a = 0.05 AU) and a negligible,
    distant pair B, so that A is effectively isolated."""
    a, m_a = 0.05, 1.8
    if spin_factor == "sync":
        p_rot = "sync"
    else:
        p_rot = units.sma_to_period(a, m_a) * units.DAYS_PER_YEAR / spin_factor
    star = {"apsidal_k2": 0.015, "rg2": 0.08, "lag_s": lag_s, "P_rot": p_rot}
    return Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 1.0, "m2": 0.8, "m3": 1e-9, "m4": 1e-9},
            "orbits": {
                "A": {"a": a, "e": e, "omega": 30.0},
                "B": {"a": 1.0, "e": 0.0},
                "out": {"a": 1.0e4, "e": 0.0},
            },
            "stars": {
                "m1": {"R": 1.5, **star},
                "m2": {"R": 1.2, **star},
                "m3": {"R": 1e-4},
                "m4": {"R": 1e-4},
            },
            "physics": {"gr": False, "tides": True, "rotation": rotation},
            "run": {"t_max_yr": t_max, "n_output": n_out, "extremes_every_orbits": 0},
        }
    )


def _sterne_expected(e: float, rotation: bool) -> float:
    from .tides import apsidal_rate_sterne

    m, a = (1.0, 0.8), 0.05
    n = np.sqrt(units.G * 1.8 / a**3)
    w = (n, n) if rotation else (0.0, 0.0)
    r = (1.5 * units.RSUN_AU, 1.2 * units.RSUN_AU)
    return 2.0 * np.pi / apsidal_rate_sterne(m, r, (0.015, 0.015), a, e, w)


def apsidal_period_static() -> float:
    res = run_nbody(_tidal_binary_config(0.3, 60.0, 600))
    return res.binaries["A"]["U_measured_yr"]


def apsidal_period_rotating() -> float:
    res = run_nbody(_tidal_binary_config(0.3, 60.0, 600, rotation=True))
    return res.binaries["A"]["U_measured_yr"]


@lru_cache(maxsize=4)
def _tidal_rates(e: float, spin_factor: float) -> tuple[float, float, float, tuple]:
    """Tidal da/dt, de/dt, dOmega_2/dt from N-body: run with and without dissipation and fit
    the difference, which removes the conservative short-period terms."""
    from .backends.nbody import build_simulation

    lag_s = 63.0
    n = np.sqrt(units.G * 1.8 / 0.05**3)
    ts = np.linspace(0.0, 60 * 2 * np.pi / n, 600)
    series = []
    for lag in (lag_s, 0.0):
        cfg = _tidal_binary_config(
            e, 1.0, 2, rotation=spin_factor > 0, lag_s=lag, spin_factor=spin_factor or 1.0
        )
        sim, *_ = build_simulation(cfg)
        a_, e_, w_ = [], [], []
        for t in ts:
            sim.integrate(t)
            o = sim.particles[1].orbit(primary=sim.particles[0])
            a_.append(o.a)
            e_.append(o.e)
            w_.append(np.linalg.norm(sim.particles[1].params["Omega"]))
        series.append(np.array([a_, e_, w_]))
    diff = series[0] - series[1]
    rates = tuple(np.polyfit(ts, d, 1)[0] for d in diff)
    return (*rates, (lag_s, n))


def _leconte_expected(e: float, spin_factor: float) -> tuple[float, float, float]:
    from .tides import leconte_rates

    n = np.sqrt(units.G * 1.8 / 0.05**3)
    lag = 63.0 / units.SECONDS_PER_YEAR
    r = (1.5 * units.RSUN_AU, 1.2 * units.RSUN_AU)
    w = (spin_factor * n, spin_factor * n)
    t = leconte_rates((1.0, 0.8), r, (0.015, 0.015), (0.08, 0.08), (lag, lag), 0.05, e, w)
    return t.da_dt, t.de_dt, t.dspin_dt[1]


def tidal_ratio(e: float, spin_factor: float, k: int) -> float:
    """N-body tidal rate / Leconte et al. (2010) rate for quantity k (0 a, 1 e, 2 spin)."""
    return float(_tidal_rates(e, spin_factor)[k] / _leconte_expected(e, spin_factor)[k])


# --- F2: secular 2+2 backend ----------------------------------------------------------------


def _run_secular(cfg: Config):
    from .backends.secular import run_secular

    cfg.run.backend = "secular"
    return run_secular(cfg)


def _secular_kozai_config(i0: float) -> Config:
    """Test-particle ZLK as a 2+2: A = 1 Msun + test particle at 1 AU, B = the perturber
    (1 Msun) split into a negligible 0.001-AU pair, circular outer orbit at 20 AU."""
    return Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 1.0, "m2": 1e-6, "m3": 0.5, "m4": 0.5},
            "orbits": {
                "A": {"a": 1.0, "e": 0.01, "i": i0},
                "B": {"a": 0.001, "e": 0.0},
                "out": {"a": 20.0, "e": 0.0},
            },
            "stars": {k: {"R": 1e-4} for k in ("m1", "m2", "m3", "m4")},
            "physics": {"gr": False},
            "run": {"t_max_yr": 20_000.0, "n_output": 400},
        }
    )


def secular_kozai_emax() -> float:
    return float(_run_secular(_secular_kozai_config(65.0)).extremes["A"]["e_max"])


def secular_nodal_precession() -> float:
    cfg = _secular_kozai_config(20.0)
    cfg.orbits["A"].e = 0.001
    res = _run_secular(cfg)
    rate = np.polyfit(res.times, np.unwrap(res.orbits["A"]["Omega"]), 1)[0]
    return float(2.0 * np.pi / abs(rate))


@lru_cache(maxsize=4)
def _secular_v994(fig: int, t_max: float):
    cfg = v994_config(t_max, 4000)
    if fig == 3:
        cfg.orbits["A"].e, cfg.orbits["A"].omega = 0.031, 10.0
        cfg.orbits["B"].e, cfg.orbits["B"].omega = 0.125, 301.0
    res = _run_secular(cfg)
    jab = np.degrees(res.mutual["A-B"])
    return float(jab.min()), float(jab.max()), res.system["J_AB_max_analytic_deg"]


def secular_gr_apsidal() -> float:
    """Secular GR apsidal period of an isolated binary (the A pair of the Kozai setup at
    a = 0.05 AU, e = 0.3, coplanar with the far perturber)."""
    cfg = _secular_kozai_config(0.0)
    cfg.masses["m2"] = 0.8
    cfg.orbits["A"].a, cfg.orbits["A"].e = 0.05, 0.3
    cfg.orbits["out"].a = 2000.0
    cfg.physics.gr = True
    cfg.run.t_max_yr, cfg.run.n_output = 30_000.0, 300
    res = _run_secular(cfg)
    return res.binaries["A"]["U_measured_yr"]


def _gr_apsidal_formula() -> float:
    from .tides import apsidal_rate_gr

    return float(2.0 * np.pi / apsidal_rate_gr(1.8, 0.05, 0.3))


# --------------------------------------------------------------------------- ensembles


def _pejcha_ensemble_config(a_out: float, t_max: float, n: int, seed: int, backend: str) -> Config:
    """Population of equal-mass 2+2 quadruples of Pejcha et al. (2013), Fig. 1 / Table 1:
    inclinations of A and B to the outer orbit isotropic, all other angles random."""
    vary = {"orbits.A.i": "isotropic", "orbits.B.i": "isotropic"}
    for k in ("Omega", "omega"):
        vary.update({f"orbits.A.{k}": "angle", f"orbits.B.{k}": "angle"})
    vary["orbits.out.omega"] = "angle"
    if a_out == 6.0:  # as configs/validation/pejcha13_F6*.yaml
        vary.update({f"orbits.{o}.M0": "angle" for o in ("A", "B", "out")})
    run = {"t_max_yr": t_max, "n_output": 500, "collisions": False, "backend": backend}
    if backend == "nbody":
        run.update(
            extremes_every_orbits=1.0,
            stop_when_e_reached={"A": 0.999, "B": 0.999},
            max_energy_error=1e-3,
            max_wall_s=1500.0,
        )
    return Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 1.0, "m2": 1.0, "m3": 1.0, "m4": 1.0},
            "orbits": {
                "A": {"a": 0.573 if a_out == 6.0 else 0.57, "e": 0.3},
                "B": {"a": 0.351 if a_out == 6.0 else 0.35, "e": 0.1},
                "out": {"a": a_out, "e": 0.3},
            },
            "physics": {"gr": False},
            "run": run,
            "ensemble": {
                "n": n,
                "seed": seed,
                "vary": vary,
            },
        }
    )


@lru_cache(maxsize=4)
def _high_e_fractions(a_out: float, t_max: float, n: int, seed: int, backend: str):
    """Fractions of members in which A and B reach 1 - e < 1e-3."""
    import os
    import tempfile

    from .ensemble import run_ensemble

    cfg = _pejcha_ensemble_config(a_out, t_max, n, seed, backend)
    with tempfile.TemporaryDirectory() as d:
        rows = run_ensemble(cfg, d, processes=os.cpu_count() or 1)
    ok = [r for r in rows if r["fate"] != "error"]
    return tuple(float(np.mean([r[f"{b}_log1me_min"] < -3.0 for r in ok])) for b in ("A", "B"))


def vok16_pejcha_fraction(binary: str) -> float:
    return _high_e_fractions(5.0, 11_000.0, 400, 16, "secular")["AB".index(binary)]


def pejcha13_f6_fraction(binary: str, backend: str = "nbody") -> float:
    n = 100 if backend == "nbody" else 400
    return _high_e_fractions(6.0, 19_470.0, n, 13, backend)["AB".index(binary)]


# --------------------------------------------------------------------------- fifth body


def quintuple_out_zlk() -> tuple[float, float]:
    """2+2+1: ZLK cycle of the orbit of the two binaries ("out") driven by a fifth star.

    The binaries A and B have negligible secondaries (1e-6 Msun), so that they act as point
    masses; the out orbit is then the inner orbit of a hierarchical triple (1 + 1 Msun, perturber
    4 Msun at 60 AU), whose quadrupole theory (Naoz et al. 2013) gives e_max. Returns
    (N-body e_max, theory e_max).
    """
    from .secular.quadrupole import predict

    cfg = Config.model_validate(
        {
            "topology": "2+2+1",
            "masses": {"m1": 1.0, "m2": 1e-6, "m3": 1.0, "m4": 1e-6, "m5": 4.0},
            "orbits": {
                "A": {"a": 0.3, "e": 0.0},
                "B": {"a": 0.3, "e": 0.0, "i": 30.0},
                "out": {"a": 8.0, "e": 0.01, "i": 65.0},
                "out5": {"a": 60.0, "e": 0.0, "i": 0.0},
            },
            "stars": {f"m{k}": {"R": 0.001} for k in range(1, 6)},
            "physics": {"gr": False},
            "run": {"t_max_yr": 6000.0, "n_output": 3000, "extremes_every_orbits": 20.0},
        }
    )
    res = run_nbody(cfg)
    theory = predict(1.0, 1.0, 4.0, 8.0, 60.0, 0.0, 0.01, 65.0, 0.0).e_max
    return float(res.extremes["out"]["e_max"]), float(theory)


@lru_cache(maxsize=1)
def _quintuple_out_zlk_cached():
    return quintuple_out_zlk()


# --------------------------------------------------------------------------- CzeV343


CZEV343 = {  # Pejcha et al. (2022) A&A 667, A53, Tables 1-3 (fiducial model)
    "P_A": 1.20388,
    "e_A": 0.109,
    "m": (1.27, 1.76),
    "R": (1.26, 1.75),
    "log_k2": (-2.07, -2.47),
}


def _czev343_spin_ratio(e: float) -> float:
    """Pseudo-synchronous spin / mean motion (Hut 1981): f2(e) / (f5(e) (1-e^2)^1.5)."""
    e2 = e * e
    n = (1 + 15 / 2 * e2 + 45 / 8 * e2**2 + 5 / 16 * e2**3) / (1 - e2) ** 6
    om = (1 + 3 * e2 + 3 / 8 * e2**2) / (1 - e2) ** 4.5
    return n / om


def czev343_apsidal_rate(part: str) -> float:
    """Apsidal rate of binary A of CzeV343 [rad/day]: tidal, rotational or GR part."""
    from .tides import apsidal_rate_gr, apsidal_rate_sterne

    c = CZEV343
    m, e = c["m"], c["e_A"]
    a = units.period_to_sma(units.days_to_years(c["P_A"]), sum(m))
    r = tuple(x * units.RSUN_AU for x in c["R"])
    k2 = tuple(10**x for x in c["log_k2"])
    n = 2.0 * np.pi / units.days_to_years(c["P_A"])
    w = n * _czev343_spin_ratio(e)
    tidal = apsidal_rate_sterne(m, r, k2, a, e, (0.0, 0.0))
    rates = {
        "tidal": tidal,
        "rotational": apsidal_rate_sterne(m, r, k2, a, e, (w, w)) - tidal,
        "gr": apsidal_rate_gr(sum(m), a, e),
    }
    return float(rates[part] / units.DAYS_PER_YEAR)


def czev343_apsidal_period_nbody() -> float:
    """Apsidal period [yr] of binary A of CzeV343 measured in a full N-body run of the
    quadruple (tides, pseudo-synchronous rotation, GR, outer orbit from the LTTE)."""
    c = CZEV343
    p_rot = c["P_A"] / _czev343_spin_ratio(c["e_A"])
    stars = {
        "m1": {"R": 1.26, "apsidal_k2": 10**-2.07, "P_rot": p_rot},
        "m2": {"R": 1.75, "apsidal_k2": 10**-2.47, "P_rot": p_rot},
        "m3": {"R": 1.37, "P_rot": "sync"},
        "m4": {"R": 1.17, "P_rot": "sync"},
    }
    cfg = Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 1.27, "m2": 1.76, "m3": 1.36, "m4": 1.20},
            "orbits": {
                "A": {"P": 1.20388, "e": 0.109, "i": 83.54, "Omega": 0.0, "omega": 30.0},
                "B": {"P": 0.80259, "e": 0.004, "i": 83.20, "Omega": 0.0},
                "out": {"P": 1453.9, "e": 0.703, "i": 11.8, "Omega": 90.0, "omega": 12.9},
            },
            "stars": stars,
            "physics": {"gr": True, "tides": True, "rotation": True},
            "run": {"t_max_yr": 40.0, "n_output": 400, "max_energy_error": 1e-4},
        }
    )
    res = run_nbody(cfg)
    return float(res.binaries["A"]["U_measured_yr"])


# --------------------------------------------------------------------------- secular 3+1


def secular_3p1_kozai_emax() -> float:
    """3+1 secular backend: test-particle ZLK of the inner orbit, i0 = 65 deg."""
    from .backends.secular import run_secular

    cfg = _hierarchical_triple(65.0, 0.01, 20_000.0, 400)
    cfg.run.backend = "secular"
    return float(run_secular(cfg).extremes["inner"]["e_max"])


@lru_cache(maxsize=1)
def _secular_naoz_fig3():
    from .backends.secular import run_secular

    cfg = _naoz_config((1.0, MJUP, 40 * MJUP), (6.0, 100.0), (0.001, 0.6), 65.0, 45.0, 5e6, 5000)
    cfg.run.backend = "secular"
    res = run_secular(cfg)
    return res.times, np.degrees(res.mutual["mid-inner"])


def secular_naoz_fig3_flip_time() -> float:
    t, i = _secular_naoz_fig3()
    up = np.flatnonzero((i[:-1] < 90.0) & (i[1:] >= 90.0))
    return float(t[up[0]]) if len(up) else float("nan")


def secular_naoz_fig3_imax() -> float:
    return float(_secular_naoz_fig3()[1].max())


# --------------------------------------------------------------------------- gas phase


def _scatter_spec(**kw):
    from .scatter import ScatterSpec

    return ScatterSpec(**kw)


def scatter_deflection() -> tuple[float, float]:
    """Single star (1 Msun) passing a single star (0.5 Msun), v_inf = 5 km/s, b = 3 AU: the
    deflection angle 2 arcsin(1/e) from the eccentricity of the relative orbit after the
    passage (Laplace-Runge-Lenz vector of the final N-body state) and the analytic value."""
    from .scatter import KMS_TO_AUYR, initial_state

    spec = _scatter_spec(
        target={"masses": [1.0]},
        projectile={"masses": [0.5]},
        v_inf_kms=5.0,
        b_au=3.0,
        r_start_factor=2.0,
    )
    m, _, pos, vel, r0, _h = initial_state(
        spec, {"id": 0, "b_au": 3.0, "q_au": None, "target": None, "projectile": None}
    )
    sim = rebound.Simulation()
    sim.G = units.G
    sim.integrator = "ias15"
    for mk, p, v in zip(m, pos, vel):
        sim.add(m=mk, x=p[0], y=p[1], z=p[2], vx=v[0], vy=v[1], vz=v[2])
    mu = units.G * m.sum()
    v_inf = 5.0 * KMS_TO_AUYR
    sim.integrate(4.0 * r0 / v_inf)
    p0, p1 = sim.particles[0], sim.particles[1]
    r = np.array([p1.x - p0.x, p1.y - p0.y, p1.z - p0.z])
    v = np.array([p1.vx - p0.vx, p1.vy - p0.vy, p1.vz - p0.vz])
    lrl = np.cross(v, np.cross(r, v)) / mu - r / np.linalg.norm(r)
    e_obt = np.linalg.norm(lrl)
    e_exp = np.sqrt(1.0 + (v_inf**2 * 3.0 / mu) ** 2)
    return float(np.degrees(2 * np.arcsin(1 / e_obt))), float(np.degrees(2 * np.arcsin(1 / e_exp)))


@lru_cache(maxsize=1)
def _at16_run() -> dict:
    """Antognini & Thompson (2016) binary-binary set-up: equal masses (1 Msun), a = 1 AU,
    circular, v_inf = v_crit, b uniform in area up to b_max = (4 v_c/v_inf + 3) a = 7 AU,
    12 000 encounters. Returns sigma / (pi a^2) per outcome class."""
    import tempfile

    from .scatter import ScatterConfig, run_scatter

    cfg = ScatterConfig(
        name="at16",
        scatter={
            "target": {"masses": [1.0, 1.0], "a_au": 1.0},
            "projectile": {"masses": [1.0, 1.0], "a_au": 1.0},
            "v_inf_vc": 1.0,
            "b_max_au": 7.0,
            "r_start_factor": 46.0,
            "n": 12000,
            "seed": 2016,
            "max_wall_s": 120.0,
        },
    )
    with tempfile.TemporaryDirectory() as d:
        rows = run_scatter(cfg, d, processes=0)
    # encounters that failed or were not finished (wall-time limit, t_max) have no outcome:
    # left out of the denominator, so that the result does not depend on the machine's load
    done = [r for r in rows if r.get("status") not in ("error", "wall_time", "unresolved")]
    print(
        f"at16: {len(rows) - len(done)} of {len(rows)} encounters without a final outcome "
        f"({100 * (len(rows) - len(done)) / len(rows):.2f} %, left out of the cross-sections)"
    )
    out = {}
    for r in done:
        out[r["outcome"]] = out.get(r["outcome"], 0) + 49.0 / max(len(done), 1)
    return out


MSE_CASES = {  # same as tools/compare_mse.py
    "3p1": {"config": "3+1", "m": [1.0, 0.8, 0.7, 0.6], "a": [1.0, 15.0, 300.0],
            "e": [0.1, 0.2, 0.3], "i": [2.0, 85.0, 10.0], "w": [30.0, 60.0, 0.0],
            "O": [0.0, 30.0, 70.0], "t_end": 1e5},
    "2p2": {"config": "2+2", "m": [1.0, 0.9, 0.8, 0.7], "a": [1.0, 1.5, 40.0],
            "e": [0.1, 0.1, 0.2], "i": [2.0, 60.0, 80.0], "w": [30.0, 60.0, 0.0],
            "O": [0.0, 100.0, 40.0], "t_end": 8e4},
}  # fmt: skip


@lru_cache(maxsize=2)
def mse_case(key: str) -> dict:
    """swingersim secular run of an MSE comparison case: {orbit: (e_max, t_e_max)}."""
    from .backends.secular import run_secular

    c = MSE_CASES[key]
    names = ["inner", "mid", "out"] if c["config"] == "3+1" else ["A", "B", "out"]
    orbits = {
        n: {"a": c["a"][k], "e": c["e"][k], "i": c["i"][k], "Omega": c["O"][k], "omega": c["w"][k]}
        for k, n in enumerate(names)
    }
    cfg = Config.model_validate(
        {
            "topology": c["config"],
            "masses": {f"m{k + 1}": m for k, m in enumerate(c["m"])},
            "orbits": orbits,
            "physics": {"gr": True, "tides": False},
            "run": {"backend": "secular", "t_max_yr": c["t_end"], "n_output": 8001},
        }
    )
    res = run_secular(cfg)
    return {n: (res.extremes[n]["e_max"], res.extremes[n]["t_e_max_yr"]) for n in names}


def mse_rlof_time() -> float:
    """3+1 quadruple whose inner binary (0.3 AU) is driven by ZLK cycles into Roche-lobe
    overflow: time when star m1 first fills its Roche lobe at pericentre (Eggleton 1983),
    with the ZAMS radii that MSE (SSE) gives the stars."""
    from .backends.secular import run_secular
    from .stellar import roche_lobe_eggleton

    radii = [0.0041295875, 0.0033901858, 0.0030302416, 0.0025896044]  # AU, from MSE
    m = [1.0, 0.8, 0.7, 0.6]
    orbits = {
        "inner": {"a": 0.3, "e": 0.1, "i": 0.5, "Omega": 0.0, "omega": 30.0},
        "mid": {"a": 5.0, "e": 0.2, "i": 89.0, "Omega": 30.0, "omega": 60.0},
        "out": {"a": 200.0, "e": 0.3, "i": 89.0, "Omega": 30.0, "omega": 0.0},
    }
    cfg = Config.model_validate(
        {
            "topology": "3+1",
            "masses": {f"m{k + 1}": x for k, x in enumerate(m)},
            "orbits": orbits,
            "stars": {f"m{k + 1}": {"R": radii[k] / units.RSUN_AU} for k in range(4)},
            "physics": {"gr": True, "tides": False},
            "run": {"backend": "secular", "t_max_yr": 600.0, "n_output": 6001},
        }
    )
    res = run_secular(cfg)
    o = res.orbits["inner"]
    fill = radii[0] / (roche_lobe_eggleton(m[0] / m[1]) * o["a"] * (1.0 - o["e"]))
    return float(res.times[np.argmax(fill >= 1.0)]) if fill.max() >= 1.0 else float("nan")


@lru_cache(maxsize=1)
def chaosmap_circumbinary() -> float:
    """Inner stability limit of a circumbinary test particle (equal-mass circular binary,
    a = 1 AU) from a MEGNO map over a (2.0-2.7, step 0.05) and 8 initial phases, 1e4 binary
    orbits: the smallest a beyond which every phase is regular."""
    import tempfile

    from .chaosmap import chaos_map, megno_array, stable_from
    from .config import load_config

    cfg = load_config(ROOT / "configs/validation/quarles18_circumbinary.yaml")
    xs = np.round(np.arange(2.0, 2.701, 0.05), 3)
    ys = np.arange(0.0, 360.0, 45.0)
    with tempfile.TemporaryDirectory() as d:
        rows = chaos_map(
            cfg, d, "orbits.out.a", xs, "orbits.out.M0", ys, t_max_yr=7071.0, processes=0
        )
    return stable_from(xs, megno_array(rows))


@lru_cache(maxsize=4)
def _escape(l_ratio: float):
    from .escape_theory import angular_momentum_scale, escape_outcomes

    m = (1.0, 1.0, 1.0)
    E0 = -units.G / 2.0  # equal masses, binary a = 1 AU
    return escape_outcomes(m, E0, l_ratio * angular_momentum_scale(m, E0))


def escape_energy_slope(l_ratio: float, u_lo: float, u_hi: float) -> float:
    """Power-law slope of dP/dE_B of the statistical escape theory between u_lo and u_hi
    (u = E_B / E0)."""
    r = _escape(l_ratio)
    sel = (r.u_grid > u_lo) & (r.u_grid < u_hi) & (r.u_density_all > 0)
    return float(np.polyfit(np.log(r.u_grid[sel]), np.log(r.u_density_all[sel]), 1)[0])


def timing_ltte_amplitude() -> float:
    """Semi-amplitude [min] of the light-time effect in the synthetic O-C of binary A of a
    2+2 with a circular, edge-on outer orbit (a = 10 AU, equal masses): A's centre of mass
    moves on a 5 AU circle, so the expected semi-amplitude is 5 AU / c."""
    from .timing import eclipse_times, o_minus_c, sinusoid_amplitude

    orb = {"e": 0.0, "i": 90.0, "Omega": 0.0}
    cfg = Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 1.0, "m2": 1.0, "m3": 1.0, "m4": 1.0},
            "orbits": {"A": {"P": 10.0, **orb}, "B": {"P": 13.0, **orb}, "out": {"a": 10.0, **orb}},
            "stars": {f"m{k}": {"R": 1.0} for k in range(1, 5)},
            "physics": {"gr": False},
            "run": {"t_max_yr": 1.0},
        }
    )
    p_out = float(np.sqrt(1000.0 / 4.0))
    ecl, p_nom = eclipse_times(cfg, p_out, ["A"])
    rows = [r for r in o_minus_c(ecl, p_nom) if r["kind"] == "primary"]
    t = np.array([r["t_obs_yr"] for r in rows])
    oc = np.array([r["o_c_d"] for r in rows])
    return sinusoid_amplitude(t, oc, p_out) * 1440.0


def _romer_cfg() -> Config:
    orb = {"e": 0.0, "Omega": 0.0}
    return Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 2.0, "m2": 1.0, "m3": 0.6, "m4": 0.6},
            "orbits": {
                "A": {"P": 5.0, "i": 89.0, **orb},
                "B": {"P": 30.0, "i": 0.0, **orb},
                "out": {"a": 500.0, "i": 0.0, **orb},
            },
            "stellar_model": {"kind": "zams"},
            "physics": {"gr": False},
            "run": {"t_max_yr": 1.0},
        }
    )


def phenom_romer_delay_expected() -> float:
    """Kaplan (2010) delay [s] of the secondary minimum of a circular binary:
    (P / pi c)(K2 - K1) = 2 a sin i (m1 - m2) / (M c), for the setup of _romer_cfg."""
    m1, m2, p_yr = 2.0, 1.0, 5.0 / units.DAYS_PER_YEAR
    a = (units.G * (m1 + m2) * p_yr**2 / (4.0 * math.pi**2)) ** (1.0 / 3.0)
    return (
        2.0
        * a
        * math.sin(math.radians(89.0))
        * (m1 - m2)
        / (m1 + m2)
        / units.C
        * (units.SECONDS_PER_YEAR)
    )


def phenom_romer_delay() -> float:
    """Delay [s] of the phenomenological secondary minima (Mikulasek 2015 model fitted to the
    synthetic light curve) after the half period of the primary ephemeris; 2 Msun + 1 Msun on
    the ZAMS, P = 5 d, circular. Proximity effects are off: Doppler beaming shifts the minima
    by a further +0.4 s (a real effect, not part of the Romer delay)."""
    from .phenom import synthetic_minima

    res = synthetic_minima(_romer_cfg(), 0.15, ["A"], cadence_min=1.0, proximity=False)
    sec = [x.o_c_d for x in res.minima if x.kind == "secondary"]
    return float(np.mean(sec)) * 86400.0


def _phoebe_test_cfg() -> Config:
    d = _romer_cfg().model_dump(exclude_none=True)
    d["orbits"]["A"].update({"e": 0.1, "omega": 40.0, "i": 88.0, "Omega": 30.0})
    return Config.model_validate(d)


def phoebe_sphere_vs_tier_a() -> float:
    """Largest difference of the normalised light curves of tier A and of PHOEBE with the same
    physics (spheres, black bodies, the same four-coefficient limb darkening, light time, no
    proximity effects) over one orbit of an eccentric eclipsing binary (needs phoebe)."""
    from .phoebe_lc import phoebe_light_curve
    from .photometry import light_curve

    cfg = _phoebe_test_cfg()
    t = np.linspace(0.0, 5.0, 601)
    a = light_curve(cfg, t, band="V", proximity=False)
    p = phoebe_light_curve(
        cfg,
        t,
        band="V",
        atm="blackbody",
        distortion="sphere",
        irradiation=False,
        swingersim_ld=True,
    )
    d = a.flux - p.flux
    return float(np.max(np.abs(d - np.median(d))))


def phoebe_romer_delay() -> float:
    """Kaplan (2010) delay [s] of the secondary minimum measured on PHOEBE light curves of the
    hybrid (with minus without light time, which removes the mesh bias of PHOEBE's minima)."""
    from .phenom import synthetic_minima

    opts = {"atm": "blackbody", "distortion": "sphere", "irradiation": False, "swingersim_ld": True}
    out = []
    for lt in (True, False):
        res = synthetic_minima(
            _romer_cfg(),
            0.08,
            ["A"],
            band="V",
            engine="phoebe",
            light_time=lt,
            phoebe_options=opts,
        )
        sec = [x.o_c_d for x in res.minima if x.kind == "secondary"]
        out.append(float(np.mean(sec)) * 86400.0)
    return out[0] - out[1]


# QZ Car = HD 93206, Brož et al. (2022), A&A 666, A24: Table 2 (fitted) and Table 3 (derived)
QZCAR_TABLE3 = {
    "a1": 54.032,  # Ac1+Ac2 [Rsun]
    "a2": 136.132,  # Aa1+Aa2
    "a3": 13053.0,  # Ac+Aa
    "R": {"m1": 16.626, "m2": 11.975, "m3": 28.286, "m4": 4.340},
    "L_V": {"m1": 0.220, "m2": 0.133, "m3": 0.635, "m4": 0.010},
}
QZCAR_TEFF = {"m1": 29687.0, "m2": 32979.0, "m3": 29564.0, "m4": 23446.0}  # Table 2


def _qzcar_cfg(with_table3_stars: bool = False) -> Config:
    import yaml

    path = ROOT / "configs/benchmarks/qzcar_nominal.yaml"
    d = yaml.safe_load(path.read_text())
    if with_table3_stars:
        for k, r in QZCAR_TABLE3["R"].items():
            d["stars"][k].update({"R": r, "Teff": QZCAR_TEFF[k]})
    return Config.model_validate(d)


def qzcar_semi_major_axis(orbit: str) -> float:
    """Semi-major axis [Rsun] of an orbit of QZ Car from the Table 2 periods and masses of
    Brož et al. (2022) (swingersim's hierarchical initialisation)."""
    from .hierarchy import build_tree

    root, _ = build_tree(_qzcar_cfg())
    node = next(n for n in root.nodes() if n.name == orbit)
    return float(node.elements.a / units.RSUN_AU)


def qzcar_light_fraction(star: str, band: str = "V") -> float:
    """Share of star in the light of QZ Car in a band, from the radii and Teff of Brož et al.
    (2022) (black body with Claret / TLUSTY limb darkening, swingersim tier A)."""
    from .photometry import star_lights
    from .stellar import star_properties

    cfg = _qzcar_cfg(True)
    m = np.array([cfg.masses[k] for k in cfg.star_names])
    props = [star_properties(cfg, f"m{k + 1}", mm) for k, mm in enumerate(m)]
    lights = star_lights(cfg, m, props, band)
    f = np.array([math.pi * s.radius_au**2 * s.weight for s in lights])
    return float(f[int(star[1:]) - 1] / f.sum())


@lru_cache(maxsize=8)
def qzcar_evolutionary_mass(star: str) -> float:
    """Present mass [Msun] of the single-star MIST model (v/vcrit = 0) with the radius and
    Teff of a component of QZ Car (Brož et al. 2022, Tables 2 and 3)."""
    from .mist import evolutionary_state

    return evolutionary_state(QZCAR_TABLE3["R"][star], QZCAR_TEFF[star])["m_now"]


def cummings2018_ifmr(m_init: float) -> float:
    """White-dwarf mass of the MIST-based semi-empirical initial-final mass relation of
    Cummings et al. (2018), ApJ 866, 21, eqs. 4-6 (0.83-7.20 Msun)."""
    if m_init < 2.85:
        return 0.080 * m_init + 0.489
    if m_init < 3.60:
        return 0.187 * m_init + 0.184
    return 0.107 * m_init + 0.471


def mist_final_mass(m_init: float) -> float:
    """Mass [Msun] of the white dwarf a MIST star (solar [Fe/H], non-rotating) leaves: its
    helium-core mass at the departure from the AGB (EEP 1409, the end of the swingersim tracks,
    defined by an envelope of 20 % of the mass left; Dotter 2016), the envelope being lost on
    the post-AGB."""
    from .mist import mist_track

    tr = mist_track(m_init, 0.0, 0.0)
    return float(tr.at(tr.t_end * (1.0 - 1e-9)).m_core)


# --------------------------------------------------------------------------- pulsations
MURPHY14_KIC9651065 = {"P": 272.70, "e": 0.468, "varpi": 2.01, "a1sini": 0.37, "m1": 1.70}


def _pulsator_binary_cfg(m1, m2, p_day, e, omega_deg, i_deg, p_puls=0.1) -> Config:
    return Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"A": ["m1", "m2"]},
            "masses": {"m1": m1, "m2": m2},
            "orbits": {"A": {"P": p_day, "e": e, "omega": omega_deg, "i": i_deg}},
            "stars": {
                "m1": {"pulsation": {"type": "manual", "period_d": p_puls, "amplitude_mag": 0.01}}
            },
            "stellar_model": {"kind": "zams"},
            "physics": {"gr": False},
            "run": {"t_max_yr": 1.0},
        }
    )


def pulsation_time_delays(cfg, span_d, p_puls=0.1, cadence_d=1.0 / 48.0, segment_d=2.0):
    """Time delays [s] of a pulsating star m1 measured as in Murphy et al. (2014): the light
    curve of the pulsator (Kepler long cadence) is cut into segments, the pulsation phase is
    fitted in each one and converted to a delay; returns (segment mid-times [d], delays)."""
    from .photometry import light_curve

    t = np.arange(0.0, span_d, cadence_d)
    lc = light_curve(cfg, t, band="V", proximity=False)
    y = -2.5 * np.log10(lc.components["m1"] / np.median(lc.components["m1"]))
    mids, tau = [], []
    for s in np.arange(0.0, span_d - segment_d / 2, segment_d):
        k = (t >= s) & (t < s + segment_d)
        w = 2.0 * math.pi * t[k] / p_puls
        x = np.column_stack([np.sin(w), np.cos(w), np.ones(k.sum())])
        a, b, _ = np.linalg.lstsq(x, y[k], rcond=None)[0]
        mids.append(t[k].mean())
        tau.append(-math.atan2(b, a))  # phase of A sin(w + phi): a later phase = negative delay
    tau = np.unwrap(np.array(tau)) / (2.0 * math.pi) * p_puls * 86400.0
    return np.array(mids), tau - tau.mean()


def time_delay_harmonics(mids, tau, p_orb, n=3, segment_d=2.0) -> list[float]:
    """Amplitudes [s] of the first n orbital harmonics of the time delays, corrected for the
    averaging over a segment (sinc factor)."""
    cols = [np.ones_like(mids)]
    for k in range(1, n + 1):
        w = 2.0 * math.pi * k * mids / p_orb
        cols += [np.sin(w), np.cos(w)]
    co = np.linalg.lstsq(np.column_stack(cols), tau, rcond=None)[0]
    out = []
    for k in range(1, n + 1):
        x = math.pi * k * segment_d / p_orb
        out.append(math.hypot(co[2 * k - 1], co[2 * k]) / (math.sin(x) / x))
    return out


def pulsation_circular_delay_expected() -> float:
    """a1 sin i / c [s] of a 1.7 + 0.5 Msun binary, P = 100 d, i = 60 deg."""
    m1, m2 = 1.7, 0.5
    a = (units.G * (m1 + m2) * (100.0 / units.DAYS_PER_YEAR) ** 2 / (4.0 * math.pi**2)) ** (1 / 3)
    return a * m2 / (m1 + m2) * math.sin(math.radians(60.0)) / units.C * units.SECONDS_PER_YEAR


def pulsation_circular_delay() -> float:
    cfg = _pulsator_binary_cfg(1.7, 0.5, 100.0, 0.0, 0.0, 60.0)
    mids, tau = pulsation_time_delays(cfg, 200.0)
    return time_delay_harmonics(mids, tau, 100.0, 2)[0]


@lru_cache(maxsize=1)
def _kic9651065_harmonics() -> tuple[float, ...]:
    """Time-delay harmonics of a model of KIC 9651065 with the orbit of Murphy et al. (2014),
    Table A1 (a1 sin i = 0.37 au, e = 0.468, varpi = 2.01 rad, P = 272.70 d; m1 = 1.70 Msun,
    i = 45 deg and m2 chosen to give that a1 sin i)."""
    from scipy.optimize import brentq

    d = MURPHY14_KIC9651065
    s_i = math.sin(math.radians(45.0))

    def a1sini(m2):
        a = units.G * (d["m1"] + m2) * (d["P"] / units.DAYS_PER_YEAR) ** 2 / (4 * math.pi**2)
        return a ** (1 / 3) * m2 / (d["m1"] + m2) * s_i

    m2 = brentq(lambda m: a1sini(m) - d["a1sini"], 0.05, 10.0)
    cfg = _pulsator_binary_cfg(d["m1"], m2, d["P"], d["e"], math.degrees(d["varpi"]), 45.0)
    mids, tau = pulsation_time_delays(cfg, 3.0 * d["P"])
    return tuple(time_delay_harmonics(mids, tau, d["P"], 4))


def delta_scuti_pl_offset(m: float) -> float:
    """M_V of a MIST star of mass m [Msun] (solar [Fe/H]) on the main sequence where it has cooled
    to Teff = 7500 K, minus M_V of the fundamental-mode P-L relation of Poro et al. (2024, eq. 10)
    at its period from Q = 0.0329 d; BC_V = 0 (|BC_V| < 0.04 mag at 6500-8500 K)."""
    from scipy.optimize import brentq

    from .pulsation import Q_RADIAL, q_period
    from .stellar import mist_properties

    def teff(age):
        return mist_properties(m, age, 0.0, 0.4).teff_k - 7500.0

    from .mist import EEP_ZAMS, mist_track

    tr = mist_track(m, 0.0, 0.4)
    ages = np.linspace(float(tr._age[EEP_ZAMS - 1]), tr.t_ms, 400)
    vals = [teff(a) for a in ages]
    k = next(j for j in range(1, len(ages)) if vals[j - 1] > 0 >= vals[j])
    age = brentq(teff, ages[k - 1], ages[k], xtol=1e-6)
    p = mist_properties(m, age, 0.0, 0.4)
    period = q_period(Q_RADIAL["F"], tr.current_mass(age), p.radius_rsun)
    mv = 4.74 - 2.5 * math.log10(p.luminosity_lsun)
    return mv - (-3.236 * math.log10(period) - 1.636)


def cepheid_radius_offset(m: float) -> float:
    """log R of a MIST Cepheid (solar [Fe/H], v/vcrit = 0.4) minus log R of the period-radius
    relation of Anderson et al. (2016), Table 5 (Z = 0.014, second crossing, average over the
    strip): the star is taken on its core-helium-burning blue loop where it is closest to the
    middle of the fundamental instability strip of Bono et al. (2000, Table 5), its period from
    their Table 6."""
    from .mist import EEP_TACHEB, EEP_ZACHEB, mist_track
    from .pulsation import cepheid_period, cepheid_strip
    from .stellar import mist_properties

    z = 0.0142
    tr = mist_track(m, 0.0, 0.4)
    best = None
    for age in np.linspace(tr._age[EEP_ZACHEB - 1], tr._age[EEP_TACHEB - 1], 2000):
        p = mist_properties(m, age, 0.0, 0.4)
        red, blue = cepheid_strip(p.luminosity_lsun, z)
        d = abs(math.log10(p.teff_k) - 0.5 * math.log10(red * blue))
        if red <= p.teff_k <= blue and (best is None or d < best[0]):
            best = (d, age, p)
    if best is None:
        raise ValueError(f"the MIST star of {m} Msun does not cross the Cepheid strip")
    _, age, p = best
    period = cepheid_period(p.luminosity_lsun, p.teff_k, tr.current_mass(age), z, "F")
    return math.log10(p.radius_rsun) - (0.652 * math.log10(period) + 1.180)


# --------------------------------------------------------------------------- mass transfer
def mt_orbit_ratio(md0, ma0, md1, alpha, beta) -> float:
    """a_f / a_i of a circular binary integrated by swingersim.masstransfer.orbit_path when the
    donor goes from md0 to md1 with constant fractions alpha (donor wind), beta (isotropic
    re-emission) and 1 - alpha - beta accreted."""
    from .masstransfer import orbit_path

    _, _, a = orbit_path(md0, ma0, 1.0, alpha, beta, md1, n=50)
    return float(a[-1])


def soberman1997_a_ratio(md0, ma0, md1, alpha, beta) -> float:
    """Soberman, Phinney & van den Heuvel (1997), eqs. 20-28 (A = 1): a/a0 as a function of
    q = m_donor / m_accretor for constant fractions alpha, beta and accreted eps."""
    eps = 1.0 - alpha - beta
    ma1 = ma0 + eps * (md0 - md1)
    q0, q = md0 / ma0, md1 / ma1
    aw = alpha
    bw = (alpha + beta) / (1.0 - eps)
    cw = alpha * eps / (1.0 - eps) + beta / (eps * (1.0 - eps))
    return (
        (q / q0) ** (2 * aw - 2)
        * ((1 + q) / (1 + q0)) ** (1 - 2 * bw)
        * ((1 + eps * q) / (1 + eps * q0)) ** (3 + 2 * cw)
    )


def _zeta_lobe(q: float) -> float:
    from .masstransfer import zeta_lobe

    return zeta_lobe(q, 1.0, 0.0, 0.0)


def postnov_conservative(md0, ma0, md1) -> float:
    """Postnov & Yungelson (2014), eq. 37: a_f/a_i = (M1i M2i / M1f M2f)^2."""
    ma1 = ma0 + md0 - md1
    return (md0 * ma0 / (md1 * ma1)) ** 2


def postnov_isotropic(md0, ma0, md1) -> float:
    """Postnov & Yungelson (2014), eq. 43 (isotropic re-emission, beta = 1, accretor mass
    constant): a_f/a_i = (M_i/M_f) (M2i/M2f)^2 exp(-2 (M2i - M2f)/M1)."""
    return (md0 + ma0) / (md1 + ma0) * (md0 / md1) ** 2 * math.exp(-2.0 * (md0 - md1) / ma0)


@lru_cache(maxsize=1)
def bse_algol_example() -> dict:
    """swingersim's mass-transfer history of the Algol example of Hurley et al. (2002), Sect. 3.1:
    2.9 + 0.9 Msun, Z = 0.02 (SSE tracks), circular orbit of P = 2.7 d (the period after the
    tidal circularisation of their e = 0.7, P = 8 d orbit)."""
    from .masstransfer import binary_history

    cfg = Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"A": ["m1", "m2"]},
            "masses": {"m1": 2.9, "m2": 0.9},
            "orbits": {"A": {"P": 2.7}},
            "stellar_model": {"kind": "sse", "age_myr": 0.0, "z": 0.02},
            "physics": {"gr": False, "evolution": True, "mass_transfer": {}},
            "run": {"t_max_yr": 1.5e9, "backend": "secular"},
        }
    )
    m = 3.8
    a0 = (units.G * m * (2.7 / units.DAYS_PER_YEAR) ** 2 / (4 * math.pi**2)) ** (1 / 3)
    h = binary_history(cfg, (0, 1), [2.9, 0.9], a0, 1.5e9)
    mf = h.m[:, -1]
    p = 2 * math.pi * math.sqrt(h.a[-1] ** 3 / (units.G * mf.sum())) * units.DAYS_PER_YEAR
    return {"m1": float(mf[0]), "m2": float(mf[1]), "P": float(p), "event": h.event}


# --------------------------------------------------------------------------- common envelope
def _ce_period_d(m1, mc1, m2, mc2, rc2, p_i_d, alpha, lam=0.5, r2=None) -> float:
    """Period [d] after a common envelope whose giant donor (m1, core mc1) fills its Roche lobe
    in a circular orbit of period p_i_d (swingersim.commonenv; the donor radius is its Roche
    lobe, Eggleton 1983)."""
    from . import commonenv
    from .stellar import roche_lobe_eggleton

    m = m1 + m2
    a_i = (units.G * m * (p_i_d / units.DAYS_PER_YEAR) ** 2 / (4 * math.pi**2)) ** (1 / 3)
    r1 = roche_lobe_eggleton(m1 / m2) * a_i / units.RSUN_AU
    rc1 = commonenv.r_white_dwarf(mc1)
    o = commonenv.outcome(m1, mc1, r1, rc1, m2, mc2, r2 or rc2, rc2, a_i, alpha, lam)
    mf = o.m1_f + o.m2_f
    return 2 * math.pi * math.sqrt(o.a_f**3 / (units.G * mf)) * units.DAYS_PER_YEAR


def ce_tout97_period_hr() -> float:
    """Tout et al. (1997) CV progenitor as quoted by Hurley et al. (2002), Sect. 3.3: 7 + 1 Msun,
    P = 1200 d; the AGB primary fills its lobe at 5.36 Msun (orbit widened by the wind,
    a M = const) and the envelope ejection (alpha_CE = 1, lambda = 0.5) leaves a 1.014 Msun
    white dwarf with the 1 Msun main-sequence star."""
    p_i = 1200.0 * (8.0 / 6.36) ** 2  # P ~ a^1.5 M^-0.5 with a M const: P ~ M^-2
    return 24.0 * _ce_period_d(5.36, 1.014, 1.0, 1.0, 0.9, p_i, 1.0)


def ce_bse_algol_second_period_d() -> float:
    """Second common envelope of the BSE Algol example (Hurley et al. 2002, Sect. 3.1): the
    2.85 Msun giant (q = 6.83, He core 0.416 Msun) fills its lobe at P = 18 d next to the
    0.417 Msun CO white dwarf; alpha_CE = 3."""
    from .commonenv import r_white_dwarf

    rwd = r_white_dwarf(0.417)
    return _ce_period_d(6.83 * 0.417, 0.416, 0.417, 0.417, rwd, 18.0, 3.0)


def _ce_triple_config(age: float = 1652.75, t_max: float = 3836.0) -> Config:
    return Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"out": ["A", "m3"], "A": ["m1", "m2"]},
            "masses": {"m1": 1.8, "m2": 0.5, "m3": 1.0},
            "orbits": {"A": {"P": 200.0}, "out": {"P": 3650.0, "e": 0.0, "i": 10.0}},
            "stellar_model": {"kind": "sse", "age_myr": age, "z": 0.02},
            "physics": {
                "gr": False,
                "evolution": True,
                "mass_transfer": {"common_envelope": {"alpha": 1.0, "lam": 0.5}},
            },
            "run": {"backend": "nbody", "t_max_yr": t_max, "n_output": 401},
        }
    )


@lru_cache(maxsize=1)
def ce_impulsive_triple() -> dict:
    """A 1.8 + 0.5 Msun binary whose giant goes through a common envelope inside a circular
    10-yr orbit of a 1 Msun tertiary (N-body): the outer semi-major axis and eccentricity
    after the impulsive loss of the envelope, against Hadjidemetriou (1963) for a circular
    orbit: a'/a = (1 - f)/(1 - 2 f), e' = f/(1 - f), f = lost fraction of the total mass."""
    import warnings

    from .backends.nbody import run_nbody
    from .masstransfer import histories

    cfg = _ce_triple_config()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        h = histories(cfg)["A"]
        res = run_nbody(cfg)
    t_ce = h.ce["t_yr"]
    k0 = int(np.searchsorted(res.times, t_ce - 30.0))
    k1 = int(np.searchsorted(res.times, t_ce + 30.0))
    a, e = res.orbits["out"]["a"], res.orbits["out"]["e"]
    m0 = 1.8 + 0.5 + 1.0
    f = (m0 - (h.ce["donor_core"] + 0.5 + 1.0)) / m0
    return {
        "a_ratio": float(np.median(a[k1:]) / a[k0]),
        "a_ratio_expected": (1 - f) / (1 - 2 * f),
        "e_after": float(np.median(e[k1:])),
        "e_expected": f / (1 - f),
        "inner_a_after": float(np.median(res.orbits["A"]["a"][k1:])),
        "a_f": float(h.ce["a_f_au"]),
        "event": res.event,
    }


# --------------------------------------------------------------------------- magnetism
SUN_RG2 = 1.05e54 / (1.989e33 * 6.957e10**2)  # Matt et al. (2015), Table 1: I_sun


def matt15_sun_period_d() -> float:
    """Rotation period [d] of a solar model (1 Msun, 1 Rsun, 5772 K, I_sun) braked with the
    torque of Matt et al. (2015) from P = 1 d to 4.55 Gyr (constant structure)."""
    from . import magnetic

    w0 = 2 * math.pi / 86400.0
    w = magnetic.spin_down(w0, 1.0, 1.0, 5772.0, SUN_RG2, 4.55e9)
    return 2 * math.pi / w / 86400.0


def matt15_unsat_ratio() -> float:
    """Omega(4.55 Gyr) / Omega(1.1375 Gyr) of the same model: the asymptotic unsaturated
    spin-down Omega ~ t^(-1/p) of Matt et al. (2015, eq. 11) gives 4^(-1/2) = 0.5."""
    from . import magnetic

    w0 = 2 * math.pi / 86400.0
    w1 = magnetic.spin_down(w0, 1.0, 1.0, 5772.0, SUN_RG2, 1.1375e9)
    w4 = magnetic.spin_down(w0, 1.0, 1.0, 5772.0, SUN_RG2, 4.55e9)
    return w4 / w1


VOLSCHOW16 = {  # Volschow et al. (2016), Tables 1-2: M, R, Teff, a [Rsun], P_bin [d], K [s], P_mod [yr]
    "HU Aqr": (0.18, 0.22, 3400, 0.69, 0.0868, 87.7, 19.4),
    "NN Ser": (0.111, 0.149, 2920, 0.934, 0.130, 27.65, 15.482),
    "V471 Tau": (0.93, 0.96, 5040, 3.3, 0.522, 137.2, 30.5),
}


def applegate_ratio(system: str) -> float:
    """Energy of the Applegate mechanism over the energy radiated in one modulation period,
    from the O - C semi-amplitude of the system (swingersim.magnetic, thin-shell model)."""
    from . import magnetic

    m, r, teff, a, p_bin, k, p_mod = VOLSCHOW16[system]
    dpp = magnetic.dp_over_p_from_oc(k, p_mod)
    return magnetic.applegate_energy_ratio(m, r, teff, a, p_bin, dpp, p_mod)


def cmdra_radius_excess(star: str, feh: float, age_myr: float = 4000.0) -> float:
    """Observed radius of CM Dra A or B over the MIST radius of its mass, minus 1."""
    from .mist import mist_track

    m, r = {"A": (0.23102, 0.2534), "B": (0.21409, 0.2398)}[star]
    return r / mist_track(m, feh, 0.0).at(age_myr).radius_rsun - 1.0


# --------------------------------------------------------------------------- environment
@lru_cache(maxsize=4)
def galactic_tide_emax(inc_deg: float = 60.0) -> float:
    """Largest eccentricity of a 2 x 10^4 AU orbit (third star of a triple) inclined by inc_deg
    to the Galactic disc, secular backend with physics.galactic_tide, over 3 Gyr."""
    import warnings

    from .backends.secular import run_secular

    cfg = Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"out": ["A", "m3"], "A": ["m1", "m2"]},
            "masses": {"m1": 1.0, "m2": 1.0, "m3": 1.0},
            "orbits": {"A": {"a": 1.0}, "out": {"a": 2.0e4, "e": 0.01, "i": inc_deg}},
            "physics": {"gr": False, "galactic_tide": {}},
            "run": {"backend": "secular", "t_max_yr": 3.0e9, "n_output": 3001},
        }
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return float(np.nanmax(run_secular(cfg).orbits["out"]["e"]))


@lru_cache(maxsize=4)
def flyby_kick_error(q_over_a: float = 20.0) -> float:
    """|De(N-body) - De(first order)| / |De(first order)| for a parabolic passage of a 1 Msun
    star at Q = q_over_a x a past an equal-mass binary (a = 1 AU, e = 0.5, inclined): the
    first-order kick of swingersim.environment against a direct IAS15 integration."""
    import rebound

    from .environment import first_order_kick

    sim = rebound.Simulation()
    sim.G = units.G
    sim.integrator = "ias15"
    sim.add(m=1.0)
    sim.add(
        m=1.0,
        a=1.0,
        e=0.5,
        inc=math.radians(50),
        omega=math.radians(30),
        Omega=math.radians(40),
        M=0.3,
    )
    mp, q = 1.0, q_over_a
    mu, p = sim.G * 3.0, 2.0 * q
    f0 = math.radians(170.0)
    r = p / (1 + math.cos(-f0))
    vf = math.sqrt(mu / p)
    com = sim.com()
    sim.add(
        m=mp,
        x=com.x + r * math.cos(-f0),
        y=com.y + r * math.sin(-f0),
        z=com.z,
        vx=com.vx + vf * math.sin(f0),
        vy=com.vy + vf * (1 + math.cos(f0)),
        vz=com.vz,
    )
    sim.move_to_com()

    def vectors():
        p0, p1 = sim.particles[0], sim.particles[1]
        rr = np.array([p1.x - p0.x, p1.y - p0.y, p1.z - p0.z])
        vv = np.array([p1.vx - p0.vx, p1.vy - p0.vy, p1.vz - p0.vz])
        mu2 = sim.G * 2.0
        hh = np.cross(rr, vv)
        a = 1 / (2 / np.linalg.norm(rr) - vv @ vv / mu2)
        return np.cross(vv, hh) / mu2 - rr / np.linalg.norm(rr), hh / math.sqrt(mu2 * a), a

    e0, j0, a0 = vectors()
    d = math.tan(f0 / 2)
    sim.integrate(math.sqrt(p**3 / mu) * (d + d**3 / 3))  # Barker's equation, -f0 -> +f0
    e1, _, _ = vectors()
    de_fo, _ = first_order_kick(e0, j0, a0, 2.0, mp, q, 1.0)
    return float(np.linalg.norm((e1 - e0) - de_fo) / np.linalg.norm(de_fo))


# --------------------------------------------------------------------------- emulators
VYN23_RUNS = {
    "2p2": ("output/validation/vyn23_2p2", "output/validation/vyn23_2p2_mlp.pkl"),
    "3p1": ("output/validation/vyn23_3p1", "output/validation/vyn23_3p1_mlp.pkl"),
}


def _vyn23_paths(kind: str):

    root = ROOT
    ens, model = (root / p for p in VYN23_RUNS[kind])
    if not ens.exists():
        raise FileNotFoundError(
            f"{ens} missing: run configs/validation/vynatheya23_{kind}_ghost.yaml with the worker "
            "(jobs 024-027)"
        )
    return ens, model


def vyn23_ma01_score(kind: str) -> float:
    """Overall score of the MA01 criterion (both nested triples) against the N-body ghost test
    on swingersim's training set of Vynatheya et al. (2023) type."""
    from .emulator import classification_scores
    from .ensemble import load_rows

    ens, _ = _vyn23_paths(kind)
    rows = [r for r in load_rows(ens) if r.get("direction") == "forward"]
    t, p = [], []
    for r in rows:
        g, m = r.get("ghost_stable"), r.get("ma01_min")
        if g in (None, "") or m in (None, ""):
            continue
        g, m = float(g), float(m)
        if math.isfinite(g) and math.isfinite(m):
            t.append(g > 0.5)
            p.append(m > 1.0)
    return classification_scores(t, p)["score"]


def vyn23_mlp_score(kind: str) -> float:
    """Held-out score of swingersim's MLP emulator trained on that set (11 parameters of
    Vynatheya et al. 2023 -> ghost-test stability)."""
    from .emulator import load

    _, model = _vyn23_paths(kind)
    if not model.exists():
        raise FileNotFoundError(f"{model} missing: run the emulator-train job (025/027)")
    return float(load(model).metrics["ghost_stable:score"])


def sbi_linear_gaussian_std() -> float:
    """Posterior standard deviation of the first parameter of a linear Gaussian model
    x = A theta + noise (sigma 0.3), recovered by neural posterior estimation (needs sbi)."""
    from .sbi_fit import train

    rng = np.random.default_rng(0)
    th = rng.uniform(-3, 3, size=(2500, 2))
    A = np.array([[1.0, 0.5], [0.0, 2.0]])
    sig = np.array([0.3, 0.3])
    model = train(th, th @ A.T, sig, ["a", "b"], ["x1", "x2"], max_epochs=300)
    return float(model.sample(np.array([0.5, 1.0]), 20000).std(axis=0)[0])


@lru_cache(maxsize=1)
def escape_vs_nbody() -> dict:
    """swingersim N-body binary-single scatterings (configs/validation/
    stone_leigh19_binary_single.yaml, 4000 encounters) compared with the statistical escape
    theory for their resonant part."""
    import tempfile

    from .escape_theory import compare_with_scatter
    from .scatter import load_scatter_config, run_scatter

    cfg = load_scatter_config(ROOT / "configs/validation/stone_leigh19_binary_single.yaml")
    with tempfile.TemporaryDirectory() as d:
        rows = run_scatter(cfg, d, processes=0)
    return compare_with_scatter(rows, [1.0, 0.6, 0.8])


@lru_cache(maxsize=1)
def _gas_run():
    """Outer orbit of a 2+2 (binaries tight enough not to interfere) in a gas disc: full gas
    density for 3000 yr with tau_a = 5000 yr, tau_e = 2000 yr, tau_i = 1500 yr."""
    from .backends.secular import run_secular

    cfg = Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 1.0, "m2": 1.0, "m3": 1.0, "m4": 1.0},
            "orbits": {
                "A": {"a": 0.01},
                "B": {"a": 0.01},
                "out": {"a": 10.0, "e": 0.5, "i": 40.0, "Omega": 30.0},
            },
            "physics": {
                "gr": False,
                "gas": {
                    "t_end_yr": 3000.0,
                    "t_decay_yr": 1000.0,
                    "plane": {"i": 0.0, "Omega": 0.0},
                    "orbits": {"out": {"tau_a_yr": 5000.0, "tau_e_yr": 2000.0, "tau_i_yr": 1500.0}},
                },
            },
            "run": {"t_max_yr": 3000.0, "n_output": 301, "backend": "secular"},
        }
    )
    o = run_secular(cfg).orbits["out"]
    return float(o["a"][-1]), float(o["e"][-1]), float(np.degrees(o["inc"][-1]))


def v994_eclipse_probability(binary: str) -> float:
    """Fraction of isotropic orientations in which binary A or B of the V994 Her setup eclipses
    (swingersim's eclipse criterion, radii R ~ M^0.8 of Vokrouhlicky 2016, Sect. 3.3)."""
    from .diagnostics.binary import eclipsing

    rng = np.random.default_rng(3)
    inc = np.arccos(rng.uniform(-1.0, 1.0, 400_000))
    m, p_d, radii = (
        ((3.01, 2.58), 2.08327, (2.4, 2.1))
        if binary == "A"
        else ((1.84, 1.93), 1.42004, (1.6, 1.7))
    )
    a = units.period_to_sma(p_d / units.DAYS_PER_YEAR, sum(m))
    orbit = {"inc": inc, "a": np.full(inc.shape, a)}
    return float(np.mean(eclipsing(orbit, sum(radii) * units.RSUN_AU)))


# --------------------------------------------------------------------------- registry

KOZAI = "Kozai (1962) AJ 67, 591; Lidov (1962) PSS 9, 719; Naoz (2016) ARA&A 54, 441"
KOZAI_LINK = "https://doi.org/10.1146/annurev-astro-081915-023315"


# --------------------------------------------------------------------------- sextuples (v0.11)

SEXTUPLE_NBODY_EMAX_A = 0.75207  # swingersim N-body (IAS15, 1PN), 400 yr, run of sextuple_config()
SEXTUPLE_NBODY_ZLK_A_YR = 70.2


def sextuple_config(backend: str = "secular") -> Config:
    """A 2+2+2 sextuple: quadruple AB (binaries A, B on a 400-d orbit) and binary C on a
    33-yr orbit; A is inclined by 60 deg to the AB orbit (ZLK cycles), the AB orbit by 50 deg to
    the outer orbit."""
    return Config.model_validate(
        {
            "topology": "2+2+2",
            "masses": {"m1": 1.0, "m2": 0.8, "m3": 0.9, "m4": 0.7, "m5": 1.1, "m6": 0.6},
            "orbits": {
                "A": {"P": 8.0, "e": 0.05, "i": 60.0, "Omega": 0.0},
                "B": {"P": 6.0, "e": 0.02, "i": 10.0, "Omega": 40.0},
                "AB": {"P": 400.0, "e": 0.2, "i": 0.0},
                "C": {"P": 12.0, "e": 0.0, "i": 30.0, "Omega": 100.0},
                "out": {"P": 12000.0, "e": 0.3, "i": 50.0, "Omega": 20.0},
            },
            "physics": {"gr": True, "tides": False},
            "run": {"t_max_yr": 400.0, "n_output": 2001, "backend": backend},
        }
    )


@lru_cache(maxsize=2)
def _sextuple_diag(backend: str = "secular"):
    from .backends.secular import run_secular

    cfg = sextuple_config(backend)
    res = run_secular(cfg) if backend == "secular" else run_nbody(cfg)
    return {d.inner: d for d in res.diagnostics}["A"]


def sextuple_emax_secular() -> float:
    return float(_sextuple_diag().e_max_true)


def sextuple_zlk_secular() -> float:
    return float(_sextuple_diag().zlk_period_yr)


def sextuple_emax_nbody() -> float:
    return float(_sextuple_diag("nbody").e_max_true)


# V994 Her as a sextuple: Zasche et al. (2023) MNRAS 520, 3127, Table 3 and Sect. 4.1
ZASCHE23 = "Zasche et al. (2023) MNRAS 520, 3127, Table 3 and Sect. 4.1"
ZASCHE23_LINK = "https://doi.org/10.1093/mnras/stad328"
V994S_OMEGADOT_OBS = {"A": 1.86, "B": 3.66, "C": 1.68}  # deg/yr, observed (Table 3)


def v994_sextuple_config(t_max: float = 30.0, gr: bool = False) -> Config:
    """The V994 Her sextuple of Zasche et al. (2023), Table 3: masses, radii, periods and
    eccentricities of A, B, C and of the A-B orbit; all orbits coplanar (the authors find the
    system plausibly close to coplanar). The C orbit (unconstrained; median period ~3000 yr in
    their Monte Carlo) is set to 3000 yr, e = 0.3."""
    inc = 84.66
    return Config.model_validate(
        {
            "topology": "2+2+2",
            "masses": {
                "m1": 2.929,
                "m2": 2.216,
                "m3": 1.913,
                "m4": 1.889,
                "m5": 1.810,
                "m6": 1.077,
            },
            "orbits": {
                "A": {"P": 2.0832039, "e": 0.0276, "i": inc},
                "B": {"P": 1.4200981, "e": 0.1186, "i": inc},
                "AB": {"P": 1062.3, "e": 0.687, "i": inc},
                "C": {"P": 1.9601064, "e": 0.1893, "i": inc},
                "out": {"P": 3000.0 * 365.25, "e": 0.3, "i": inc},
            },
            "stars": {
                "m1": {"R": 2.118},
                "m2": {"R": 1.721},
                "m3": {"R": 1.572},
                "m4": {"R": 1.579},
                "m5": {"R": 1.544},
                "m6": {"R": 0.961},
            },
            "physics": {"gr": gr, "tides": False},
            "run": {"t_max_yr": t_max, "n_output": 3001},
        }
    )


@lru_cache(maxsize=1)
def _v994_sextuple_rates() -> dict[str, float]:
    """Apsidal advance [deg/yr] of A, B, C caused by the other stars alone (Newtonian point
    masses, no GR, no tides): linear fit of the longitude of periastron over 30 yr."""
    res = run_nbody(v994_sextuple_config())
    t = np.asarray(res.times)
    out = {}
    for b in ("A", "B", "C"):
        o = res.orbits[b]
        w = np.unwrap(np.asarray(o["omega"]) + np.asarray(o["Omega"]))
        out[b] = float(np.degrees(np.polyfit(t, w, 1)[0]))
    return out


def v994s_dyn_fraction(b: str) -> float:
    return _v994_sextuple_rates()[b] / V994S_OMEGADOT_OBS[b]


def v994s_gr_fraction(b: str) -> float:
    cfg = v994_sextuple_config()
    root, masses = build_tree(cfg)
    node = {n.name: n for n in root.nodes()}[b]
    m = masses[node.left] + masses[node.right]
    a, e = node.elements.a, node.elements.e
    rate = 3.0 * (units.G * m) ** 1.5 / (units.C**2 * a**2.5 * (1.0 - e**2))  # rad/yr
    return float(np.degrees(rate)) / V994S_OMEGADOT_OBS[b]


# --------------------------------------------------------------------------- stellar evolution (F4)

HURLEY2000 = "Hurley, Pols & Tout (2000) MNRAS 315, 543"
HURLEY2000_LINK = "https://doi.org/10.1046/j.1365-8711.2000.03426.x"


def sse_sun() -> dict:
    """1 Msun, Z = 0.02 at the solar age 4.57 Gyr from the tabulated SSE tracks."""
    from .evolution import track

    s = track(1.0, 0.02).at(4570.0)
    return {
        "R": s.radius_rsun,
        "L": s.luminosity_lsun,
        "bcz": 1.0 - s.r_env / s.radius_rsun,
    }


def sse_zams_radius(m: float) -> float:
    from .evolution import track

    return track(m, 0.02).at(0.0).radius_rsun


def sse_ms_lifetime(m: float) -> float:
    from .evolution import track

    return track(m, 0.02).t_ms


def sse_tams_radius(m: float) -> float:
    from .evolution import track

    t = track(m, 0.02)
    return t.at(t.t_ms * (1.0 - 1e-6)).radius_rsun


MIST = "MIST v1.2, [Fe/H] = 0, v/vcrit = 0 (Choi et al. 2016, ApJ 823, 102), EEP tracks"
MIST_LINK = "https://doi.org/10.3847/0004-637X/823/2/102"
MIST_WHY = (
    "independent stellar-evolution code (MESA) for the SSE tracks; differences of ~10 % come "
    "from convective overshooting and the solar abundance scale (Z = 0.02 vs 0.0142)"
)


def mist_loo(mass: float, eep: int, col: int, feh_tag: str = "p0.00", vvcrit: float = 0.0) -> float:
    """Column ``col`` of a MIST grid track at an EEP interpolated from the neighbouring grid
    masses with the track itself left out (twice the grid spacing), divided by the tabulated
    value; log columns (log_L, log_Teff, log_R) are compared as linear quantities."""
    from .mist import _grid

    g = _grid(feh_tag, vvcrit)
    masses = list(np.round(g["mass"], 4))
    k = masses.index(mass)
    lm = np.log(g["mass"])
    w = (lm[k] - lm[k - 1]) / (lm[k + 1] - lm[k - 1])
    e = g["eep"].astype(float)
    if col in (1, 5, 6):  # masses: interpolated as fractions of the initial mass (as in mist.py)
        m = g["mass"]
        v = m[k] * (
            (1 - w) * e[k - 1, eep - 1, col] / m[k - 1] + w * e[k + 1, eep - 1, col] / m[k + 1]
        )
        return float(v / e[k, eep - 1, col])
    v = (1 - w) * e[k - 1, eep - 1, col] + w * e[k + 1, eep - 1, col]
    if col in (2, 3, 4):
        return float(10.0 ** (v - e[k, eep - 1, col]))
    return float(v / e[k, eep - 1, col])


def mist_tau_h(mass: float, vvcrit: float) -> float:
    """Main-sequence (hydrogen-burning) lifetime [Myr] of a MIST track, ZAMS (EEP 202) to TAMS
    (EEP 454), [Fe/H] = 0."""
    from .mist import mist_track

    tr = mist_track(mass, 0.0, vvcrit)
    return float(tr._age[453] - tr._age[201])


def vink2001_hot(lum: float, mass: float, teff: float, z: float) -> float:
    """log10 of the mass-loss rate [Msun/yr] of Vink et al. (2001), eq. 24 (27.5-50 kK,
    v_inf/v_esc = 2.6, Z_sun = 0.019)."""
    x = math.log10(teff / 40000.0)
    return (
        -6.697
        + 2.194 * math.log10(lum / 1e5)
        - 1.313 * math.log10(mass / 30.0)
        - 1.226 * math.log10(2.6 / 2.0)
        + 0.933 * x
        - 10.92 * x**2
        + 0.85 * math.log10(z / 0.019)
    )


def mist_mdot_over_vink(mass: float = 40.0, frac: float = 0.3) -> float:
    """Wind mass-loss rate of the non-rotating solar MIST track at the fraction ``frac`` of its
    main-sequence life divided by the Vink et al. (2001) rate for its L, M, Teff and Z."""
    from .mist import mist_track

    tr = mist_track(mass, 0.0, 0.0)
    age = tr._age[201] + frac * (tr._age[453] - tr._age[201])
    st = tr.at(age)
    z = 0.0142  # MIST protosolar Z for [Fe/H] = 0 (Choi et al. 2016)
    log_v = vink2001_hot(st.luminosity_lsun, tr.current_mass(age), st.teff_k, z)
    return float(-tr.mdot(age) / 10.0**log_v)


def claret2004_row(quantity: str) -> float:
    """Claret (2004) table of 39.81 Msun (log M = 1.6), model at 3.654 Myr (log g = 3.416212,
    log k2 = -2.815231, beta = 0.161706): log k2 or beta interpolated by swingersim at that log g."""
    from .stellar import claret2004_structure

    k2, rg2 = claret2004_structure(10.0**1.6, 3.416212, "main")
    return math.log10(k2) if quantity == "log_k2" else math.sqrt(rg2)


def tlusty_row_a1() -> float:
    """a1 of the TLUSTY model OG40000g400v10 (solar, 40 000 K, log g 4.0) in Bessell V as
    interpolated by swingersim at that grid point (Reeve & Howarth 2016 table: +1.18262)."""
    from .photometry import tlusty_coefficients

    return float(tlusty_coefficients(40000.0, 4.0, 0.0, "V")[0])


def tlusty_vs_atlas_mean(teff: float = 25000.0) -> float:
    """Disk-averaged intensity <I>/I(1) in V of a 25 000 K, log g 4.0 star: TLUSTY non-LTE
    (Reeve & Howarth 2016) / ATLAS LTE (Claret & Bloemen 2011)."""
    from .photometry import claret_coefficients, ld_mean, tlusty_coefficients

    return float(
        ld_mean(tlusty_coefficients(teff, 4.0, 0.0, "V"))
        / ld_mean(claret_coefficients(teff, 4.0, 0.0, "V")[0])
    )


def mist_tams_age(mass: float = 1.0) -> float:
    """Age [Myr] at the terminal-age main sequence (EEP 454) of the swingersim MIST track."""
    from .mist import mist_track

    return mist_track(mass, 0.0, 0.0).t_ms


def iw_two_layer_ratio(alpha: float = 0.6, beta: float = 0.95) -> float:
    """Inertial-wave dissipation of a (smoothed) two-layer star from the impulsive formalism of
    Barker (2020, eqs. 12-38) divided by the closed two-layer formula of Ogilvie (2013) in
    the form of Barker (2020, eq. 39)."""
    from .inertial_waves import barker_profile, two_layer, two_layer_profile

    return barker_profile(*two_layer_profile(alpha, beta)) / two_layer(alpha, beta)


def iw_log_q(mass: float, age_frac: float = 0.3, p_rot_d: float = 10.0) -> float:
    """log10 <Q'_IW> of an SSE star at a fraction of its main-sequence life, rotating with
    p_rot_d days."""
    from .evolution import track
    from .tides import inertial_wave_q

    t = track(mass, 0.02)
    st = t.at(age_frac * t.t_ms)
    w = 2.0 * np.pi / (p_rot_d / units.DAYS_PER_YEAR)
    return float(np.log10(inertial_wave_q(mass, st.radius_rsun, st.m_env, st.r_env, w)))


def claret_k2_beta(mass: float, logg: float, feh: float = 0.0) -> tuple[float, float]:
    """(k2, beta) of the Claret (2023) grid as swingersim interpolates it."""
    from .stellar import claret_structure

    k2, rg2 = claret_structure(mass, logg, feh)
    return k2, math.sqrt(rg2)


def circularization_period(
    age_myr: float = 4000.0, efficiency: float = 1.0, inertial_waves: bool = False
) -> float:
    """Initial period [d] below which twin 1 Msun binaries (e0 = 0.3, near-synchronous spins)
    are circular (e < 0.01) at age_myr: SSE stars, k2 and I of Claret (2023), constant-time-lag
    tides with the lag of Hurley et al. (2002), orbit-averaged (Leconte et al. 2010)."""
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq

    from .evolution import track
    from .stellar import LOGG_SUN, claret_structure
    from .tides import hurley_k_over_t, inertial_wave_k_over_t, lag_from_k_over_t, leconte_rates

    m = (1.0, 1.0)
    tr = track(1.0, 0.02)

    def rhs(t, y):
        a, e, w1, w2 = y
        e = min(max(e, 1e-9), 0.99)
        st = tr.at(1.0 + t / 1e6)
        k2, rg2 = claret_structure(1.0, LOGG_SUN - 2.0 * np.log10(st.radius_rsun))
        r_au = st.radius_rsun * units.RSUN_AU
        p_orb = 2 * np.pi / np.sqrt(units.G * 2.0 / a**3)

        def k_over_t(w):
            kt = hurley_k_over_t(
                1.0, st.radius_rsun, st.m_env, st.convective_turnover_yr, 1.0, a, p_orb, w
            )
            if inertial_waves:
                kt += inertial_wave_k_over_t(1.0, st.radius_rsun, st.m_env, st.r_env, 1.0, a, w)
            return kt

        lags = tuple(
            efficiency
            * lag_from_k_over_t(k_over_t(w), 1.0, st.radius_rsun, k2)
            / units.SECONDS_PER_YEAR
            for w in (w1, w2)
        )
        r = leconte_rates(m, (r_au, r_au), (k2, k2), (rg2, rg2), lags, a, e, (w1, w2))
        return [r.da_dt, r.de_dt, r.dspin_dt[0], r.dspin_dt[1]]

    def log_e_final(p_day: float) -> float:
        a = (units.G * 2.0 * (p_day / units.DAYS_PER_YEAR) ** 2 / (4 * np.pi**2)) ** (1 / 3)
        w = 0.9 * np.sqrt(units.G * 2.0 / a**3)
        sol = solve_ivp(
            rhs,
            (0.0, age_myr * 1e6),
            [a, 0.3, w, w],
            method="LSODA",
            rtol=1e-6,
            atol=[1e-12, 1e-9, 1e-6, 1e-6],
        )
        return float(np.log10(max(sol.y[1, -1], 1e-9) / 0.01))

    return float(brentq(log_e_final, 1.0, 40.0, xtol=0.02))


def circularization_period_mist(
    age_myr: float = 4000.0,
    start_myr: float = 1.0,
    inertial_waves: bool = True,
    mass=1.0,
    reduction: str = "quadratic",
    final: bool = False,
) -> float:
    """As circularization_period, but with MIST stars (solar metallicity, no rotation) from the
    pre-main sequence (start_myr): large radii and convective envelopes of young stars
    (swingersim.stellar._pms_envelope), then the main sequence."""
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq

    from .mist import mist_track
    from .stellar import LOGG_SUN, claret_structure, mist_envelope
    from .tides import hurley_k_over_t, inertial_wave_k_over_t, lag_from_k_over_t, leconte_rates

    m = (mass, mass)
    tr = mist_track(mass, 0.0, 0.0)
    cache: dict = {}

    def star(age):
        key = round(age, 3)
        if key not in cache:
            st = tr.at(age)
            m_env, r_env = mist_envelope(tr, age, 0.0)
            seg = "pms" if st.phase == 0 else "main"
            logg = LOGG_SUN + math.log10(mass) - 2.0 * math.log10(st.radius_rsun)
            try:
                k2, rg2 = claret_structure(mass, logg, 0.0, seg)
            except ValueError:
                k2, rg2 = claret_structure(mass, logg, 0.0, "main")
            x = m_env * r_env * (st.radius_rsun - 0.5 * r_env)
            turn = 0.4311 * (x / (3.0 * st.luminosity_lsun)) ** (1 / 3) if x > 0 else 0.0
            cache[key] = (st.radius_rsun, m_env, r_env, k2, rg2, turn)
        return cache[key]

    def rhs(t, y):
        a, e, w1, w2 = y
        if a <= 0.0:
            return [0.0, 0.0, 0.0, 0.0]
        e = min(max(e, 1e-9), 0.99)
        r, m_env, r_env, k2, rg2, turn = star(start_myr + t / 1e6)
        r_au = r * units.RSUN_AU
        p_orb = 2 * np.pi / np.sqrt(units.G * 2.0 * mass / a**3)

        def k_over_t(w):
            kt = hurley_k_over_t(mass, r, m_env, turn, mass, a, p_orb, w, reduction)
            if inertial_waves:
                kt += inertial_wave_k_over_t(mass, r, m_env, r_env, mass, a, w)
            return kt

        lags = tuple(
            lag_from_k_over_t(k_over_t(w), mass, r, k2) / units.SECONDS_PER_YEAR for w in (w1, w2)
        )
        rr = leconte_rates(m, (r_au, r_au), (k2, k2), (rg2, rg2), lags, a, e, (w1, w2))
        return [rr.da_dt, rr.de_dt, rr.dspin_dt[0], rr.dspin_dt[1]]

    def log_e_final(p_day: float) -> float:
        a = (units.G * 2.0 * mass * (p_day / units.DAYS_PER_YEAR) ** 2 / (4 * np.pi**2)) ** (1 / 3)
        w = 0.9 * np.sqrt(units.G * 2.0 * mass / a**3)
        sol = solve_ivp(
            rhs,
            (0.0, (age_myr - start_myr) * 1e6),
            [a, 0.3, w, w],
            method="LSODA",
            rtol=1e-6,
            atol=[1e-12, 1e-9, 1e-6, 1e-6],
        )
        last["a"] = float(sol.y[0, -1])
        return float(np.log10(max(sol.y[1, -1], 1e-9) / 0.01))

    last: dict = {}
    p0 = float(brentq(log_e_final, 2.0, 40.0, xtol=0.02))
    if not final:
        return p0
    log_e_final(p0)  # the final period of the critical binary (observers' definition)
    return float(units.sma_to_period(last["a"], 2.0 * mass) * units.DAYS_PER_YEAR)


# --------------------------------------------------------------------------- photometry (F4 b)


def lc_occultation(rp: float, z: float) -> float:
    """Flux of a quadratically limb-darkened star (u = 0.4, 0.25) occulted by a dark disk of
    radius rp at the projected distance z (both in stellar radii)."""
    from .photometry import blocked_fraction

    u1, u2 = 0.4, 0.25
    a = np.array([0.0, u1 + 2 * u2, 0.0, -u2])
    return 1.0 - blocked_fraction(1.0, a, [(z, 0.0, rp)])


def _two_star_curve(r1, r2, w2, inc_deg, a1, a2, phases, gd=0.0, proximity=False, q=0.8):
    from .photometry import StarLight, flux_at

    inc = np.radians(inc_deg)
    m = np.array([1.0, q])
    lights = [
        StarLight("m1", r1, 1.0, np.asarray(a1, float), gd, 0.0, 0.0),
        StarLight("m2", r2, w2, np.asarray(a2, float), gd, 0.0, 0.0),
    ]
    out = []
    for p in phases:
        th = 2 * np.pi * p
        r = np.array([np.sin(th), -np.cos(th) * np.cos(inc), -np.cos(th) * np.sin(inc)])
        v = np.array([np.cos(th), np.sin(th) * np.cos(inc), np.sin(th) * np.sin(inc)]) * 1e-9
        pos = np.array([-m[1] / m.sum() * r, m[0] / m.sum() * r])
        vel = np.array([-m[1] / m.sum() * v, m[0] / m.sum() * v])
        out.append(flux_at(pos, vel, m, lights, [(0, 1)], proximity)[0])
    return np.array(out)


def lc_detached_depths() -> dict[str, float]:
    """Primary and secondary eclipse depths of a detached binary of spheres (R/a = 0.12, 0.08,
    i = 87 deg, surface-brightness ratio 0.4, four-coefficient limb darkening)."""
    a1 = [0.5, 0.2, 0.1, -0.05]
    a2 = [0.6, 0.1, 0.15, -0.1]
    prim = _two_star_curve(0.12, 0.08, 0.4, 87.0, a1, a2, np.linspace(-0.01, 0.01, 401))
    sec = _two_star_curve(0.12, 0.08, 0.4, 87.0, a1, a2, np.linspace(0.49, 0.51, 401))
    return {"primary": 1.0 - prim.min(), "secondary": 1.0 - sec.min()}


def lc_ellipsoidal_amplitude(r1: float = 0.15) -> float:
    """Semi-amplitude of the ellipsoidal variation (R1/a = 0.15, q = 0.8, i = 60 deg, quadratic
    limb darkening 0.5, 0.2, gravity darkening y = 0.32; companion dark)."""
    u1, u2 = 0.5, 0.2
    a = [0.0, u1 + 2 * u2, 0.0, -u2]
    f = _two_star_curve(r1, 0.8 * r1, 1e-12, 60.0, a, a, np.linspace(0, 1, 401), 0.32, True)
    return float((f.max() - f.min()) / 2.0 / f.mean())


# --------------------------------------------------------------------------- mass loss (F4 d)


def _ml_binary(rate: float, t_max: float, e: float = 0.0, n_output: int = 101) -> Config:
    """Equal-mass binary, 1 + 1 Msun at a = 1 AU, star m1 losing mass at a constant rate."""
    return Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"A": ["m1", "m2"]},
            "masses": {"m1": 1.0, "m2": 1.0},
            "orbits": {"A": {"a": 1.0, "e": e, "i": 0.0, "Omega": 0.0, "omega": 0.0, "M0": 0.0}},
            "physics": {"gr": False, "mass_loss": {"source": "rates", "rates": {"m1": rate}}},
            "run": {"t_max_yr": t_max, "n_output": n_output, "collisions": False},
        }
    )


def ml_adiabatic_invariant() -> float:
    """Largest relative deviation of a M from its initial value while m1 loses 0.3 Msun
    slowly (Psi = 5.6e-6)."""
    res = run_nbody(_ml_binary(-1e-4, 3000.0))
    m = res.stellar["m1"]["mass_msun"] + 1.0
    am = res.orbits["A"]["a"] * m
    return float(np.max(np.abs(am / am[0] - 1.0)))


def _ml_true_anomaly(o) -> np.ndarray:
    from .kepler import solve_kepler

    ecc = np.array([solve_kepler(m, e) for m, e in zip(o["M"], o["e"])])
    return 2.0 * np.arctan2(
        np.sqrt(1.0 + o["e"]) * np.sin(ecc / 2.0), np.sqrt(1.0 - o["e"]) * np.cos(ecc / 2.0)
    )


def ml_osculating_e(e0: float) -> float:
    """Largest |e_Nbody - e_theory| / Psi over two orbits, e_theory the integral of Veras et
    al. (2011) eq. 4 at first order in Psi (massloss.veras_adiabatic_e)."""
    from . import massloss

    rate = -1e-3
    res = run_nbody(_ml_binary(rate, 2.0, e=e0, n_output=801))
    o = res.orbits["A"]
    psi = massloss.psi(rate, 2.0, 1.0)
    f = _ml_true_anomaly(o)
    pred = np.array([massloss.veras_adiabatic_e(e0, psi, x) for x in f])
    return float(np.max(np.abs(o["e"] - pred)) / psi)


def ml_circular_emax() -> float:
    """Largest eccentricity / Psi of an initially circular binary under slow mass loss."""
    from . import massloss

    rate = -1e-3
    res = run_nbody(_ml_binary(rate, 2.0, e=0.0, n_output=801))
    return float(np.max(res.orbits["A"]["e"]) / massloss.psi(rate, 2.0, 1.0))


@lru_cache(maxsize=1)
def _ml_impulsive() -> tuple[float, float]:
    """A circular binary loses 15 % of its total mass within 1e-3 of its period."""
    p = units.sma_to_period(1.0, 2.0)
    dt = 1e-3 * p
    res = run_nbody(_ml_binary(-0.3 / dt, dt, n_output=11))
    return float(res.orbits["A"]["a"][-1]), float(res.orbits["A"]["e"][-1])


def _ml_zlk_triple(loss: float, backend: str) -> Config:
    """Test-particle ZLK triple (i0 = 65 deg, e0 = 0.2) whose perturber m3 loses the
    fraction ``loss`` of its mass in 30 000 yr."""
    d = _hierarchical_triple(65.0, 0.2, 30000.0, 6001).model_dump(exclude_none=True)
    d["orbits"]["out"]["a"] = 200.0  # keep the distant dummy orbit adiabatic (Psi << 1)
    d["run"]["backend"] = backend
    if loss:
        d["physics"]["mass_loss"] = {"source": "rates", "rates": {"m3": -loss / 30000.0}}
    return Config.model_validate(d)


def _ml_up_crossings(t, e, level: float = 0.6) -> np.ndarray:
    k = np.flatnonzero((e[:-1] < level) & (e[1:] >= level))
    tc = t[k] + (level - e[k]) / (e[k + 1] - e[k]) * (t[k + 1] - t[k])
    out: list[float] = []
    for x in tc:
        if not out or x - out[-1] > 1000.0:
            out.append(float(x))
    return np.array(out)


def ml_zlk_nbody_vs_quadrupole() -> float:
    """Time of the 4th rise of e above 0.6 while the perturber loses half its mass: N-body
    relative to the quadrupole secular theory (backend secular, whose ZLK flow is the
    constant-mass flow re-parametrised in time by the coefficient c2(t), see the tests)."""
    from .backends.secular import run_secular

    rn = run_nbody(_ml_zlk_triple(0.5, "nbody"))
    rs = run_secular(_ml_zlk_triple(0.5, "secular"))
    tn = _ml_up_crossings(rn.times, rn.orbits["inner"]["e"])
    ts = _ml_up_crossings(rs.times, rs.orbits["inner"]["e"])
    return float(tn[3] / ts[3])


def _ml_tedi_config(t_loss: float):
    """Coplanar triple with the masses of Perets & Kratter (2012, sect. 2: 5 + 2 Msun inner
    binary, 2.25 Msun tertiary); m1 loses 3 Msun in t_loss; a_in = 1 AU, e_out = 0.1 and a_out
    at 1.15 times the Mardling & Aarseth (2001) limit."""
    from .stability import ma01_ratio

    m1, m2, m3, e_out = 5.0, 2.0, 2.25, 0.1
    a_out0 = 1.15 / ma01_ratio(1.0, 1.0, e_out, m3 / (m1 + m2), 0.0)
    rate = -3.0 / t_loss
    cfg = Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"out": ["A", "m3"], "A": ["m1", "m2"]},
            "masses": {"m1": m1, "m2": m2, "m3": m3},
            "orbits": {
                "A": {"a": 1.0, "e": 0.0, "i": 0.0, "M0": 0.0},
                "out": {"a": a_out0, "e": e_out, "i": 0.5, "M0": 90.0},
            },
            "physics": {"gr": False, "mass_loss": {"source": "rates", "rates": {"m1": rate}}},
            "run": {"t_max_yr": t_loss, "n_output": 2001},
        }
    )

    def ratio(t):
        m_in0, m_in = m1 + m2, m1 + rate * t + m2
        return ma01_ratio(m_in0 / m_in, a_out0 * (m_in0 + m3) / (m_in + m3), e_out, m3 / m_in, 0.0)

    return cfg, ratio


def ml_tedi_onset() -> float:
    """Triple evolution dynamical instability (Perets & Kratter 2012): time at which the inner
    eccentricity first exceeds 0.1, relative to the time at which the adiabatically expanding
    orbits cross the Mardling & Aarseth (2001) stability limit."""
    from scipy.optimize import brentq

    cfg, ratio = _ml_tedi_config(8000.0)
    t_ma = brentq(lambda t: ratio(t) - 1.0, 0.0, 8000.0)
    res = run_nbody(cfg)
    e = res.orbits["A"]["e"]
    k = int(np.argmax(e > 0.1))
    if e[k] <= 0.1:
        return float("inf")
    return float(res.times[k] / t_ma)


BENCHMARKS: list[Benchmark] = [
    Benchmark(
        id="v994s-dyn-A",
        tier="D",
        quantity="V994 Her sextuple: share of the apsidal motion of A caused by B (and C)",
        unit="",
        expected=0.12,
        tolerance=0.25,
        relative=True,
        source=ZASCHE23,
        link=ZASCHE23_LINK,
        rationale="the paper quotes a rounded 12 %; the rate scales as (1 - e_AB^2)^-1.5 and the "
        "error of e_AB (0.687 +0.050/-0.037) alone changes it by up to ~25 %; coplanarity assumed",
        compute=partial(v994s_dyn_fraction, "A"),
    ),
    Benchmark(
        id="v994s-dyn-B",
        tier="D",
        quantity="V994 Her sextuple: share of the apsidal motion of B caused by A (and C)",
        unit="",
        expected=0.05,
        tolerance=0.25,
        relative=True,
        source=ZASCHE23,
        link=ZASCHE23_LINK,
        rationale="as for A (rounded 5 %, error of e_AB)",
        compute=partial(v994s_dyn_fraction, "B"),
    ),
    Benchmark(
        id="v994s-dyn-C",
        tier="D",
        quantity="V994 Her sextuple: share of the apsidal motion of C caused by A and B",
        unit="",
        expected=0.0,
        tolerance=0.01,
        relative=False,
        source=ZASCHE23,
        link=ZASCHE23_LINK,
        rationale="the paper finds the classical perturbations of the distant C negligible; "
        "with a 3000-yr C orbit the share must stay below 1 %",
        compute=partial(v994s_dyn_fraction, "C"),
    ),
    Benchmark(
        id="v994s-gr-A",
        tier="D",
        quantity="V994 Her: relativistic share of the apsidal motion of A",
        unit="",
        expected=0.10,
        tolerance=0.10,
        relative=True,
        source=ZASCHE23,
        link=ZASCHE23_LINK,
        rationale="1PN rate of A divided by the observed 1.86 deg/yr; the paper quotes a rounded "
        "10 %",
        compute=partial(v994s_gr_fraction, "A"),
    ),
    Benchmark(
        id="v994s-gr-B",
        tier="D",
        quantity="V994 Her: relativistic share of the apsidal motion of B",
        unit="",
        expected=0.08,
        tolerance=0.10,
        relative=True,
        source=ZASCHE23,
        link=ZASCHE23_LINK,
        rationale="1PN rate of B divided by the observed 3.66 deg/yr; the paper quotes 8 %",
        compute=partial(v994s_gr_fraction, "B"),
    ),
    Benchmark(
        id="v994s-gr-C",
        tier="D",
        quantity="V994 Her: relativistic share of the apsidal motion of C",
        unit="",
        expected=0.06,
        tolerance=0.10,
        relative=True,
        source=ZASCHE23,
        link=ZASCHE23_LINK,
        rationale="1PN rate of C divided by the observed 1.68 deg/yr",
        compute=partial(v994s_gr_fraction, "C"),
        explanation="the GR rate from the paper's own Table 3 (1.81 + 1.077 Msun, P = 1.960 d, "
        "e = 0.189; Gimenez 1985) is 0.136 deg/yr, 8.1 % of the observed 1.68 deg/yr - "
        "swingersim's value; the ~10 % and ~8 % quoted for A and B are reproduced the same way, so "
        "the 6 % quoted for C is probably a slip of the paper (the primary's mass alone gives "
        "6 %)",
    ),
    Benchmark(
        id="sextuple-secular-emax",
        tier="C",
        quantity="2+2+2 sextuple: maximum eccentricity of A, secular backend vs N-body",
        unit="",
        expected=SEXTUPLE_NBODY_EMAX_A,
        tolerance=0.02,
        relative=True,
        source="swingersim N-body backend (REBOUND IAS15 + 1PN), validated separately",
        link="https://rebound.readthedocs.io",
        rationale="the secular backend averages over all orbits and neglects the binary-triplet "
        "terms; with a_A/a_AB = 0.08 and a_AB/a_out = 0.10 the expected agreement is ~1 %",
        compute=sextuple_emax_secular,
    ),
    Benchmark(
        id="sextuple-secular-zlk",
        tier="C",
        quantity="2+2+2 sextuple: ZLK period of A, secular backend vs N-body",
        unit="yr",
        expected=SEXTUPLE_NBODY_ZLK_A_YR,
        tolerance=0.05,
        relative=True,
        source="swingersim N-body backend (REBOUND IAS15 + 1PN), validated separately",
        link="https://rebound.readthedocs.io",
        rationale="periods are measured between eccentricity maxima on the 0.2-yr output grid "
        "(~5 cycles in 400 yr); 5 % covers the grid and the averaging",
        compute=sextuple_zlk_secular,
    ),
    Benchmark(
        id="sextuple-nbody-emax",
        tier="C",
        quantity="2+2+2 sextuple: maximum eccentricity of A, N-body (reproducibility)",
        unit="",
        expected=SEXTUPLE_NBODY_EMAX_A,
        tolerance=0.002,
        relative=True,
        source="SWINGERsim reference run (regression value)",
        link="https://rebound.readthedocs.io",
        rationale="guards the general-hierarchy N-body path (six stars, two nested quadruples) "
        "against regressions",
        compute=sextuple_emax_nbody,
        slow=True,
    ),
    Benchmark(
        id="gr-mercury",
        tier="A",
        quantity="1PN perihelion advance of Mercury",
        unit="arcsec/century",
        expected=42.98,
        tolerance=0.005,
        relative=True,
        source="Will (2014) Living Rev. Relativ. 17, 4",
        link="https://doi.org/10.12942/lrr-2014-4",
        rationale="classic value 42.98 arcsec/century; only Sun + Mercury, so planetary "
        "perturbations (the other ~532 arcsec/century) are absent by construction",
        compute=mercury_perihelion,
    ),
    Benchmark(
        id="zlk-emax-65",
        tier="A",
        quantity="ZLK maximum eccentricity, i0 = 65 deg, test particle",
        unit="",
        expected=_zlk_emax_formula(65.0),
        tolerance=0.01,
        relative=True,
        source=KOZAI,
        link=KOZAI_LINK,
        rationale="formula is the quadrupole test-particle limit; a_out/a_in = 20 and "
        "e0 = 0.01 add corrections of a few 0.1 %",
        compute=zlk_emax_65,
    ),
    Benchmark(
        id="zlk-emax-80",
        tier="A",
        quantity="ZLK maximum eccentricity, i0 = 80 deg, test particle",
        unit="",
        expected=_zlk_emax_formula(80.0),
        tolerance=0.01,
        relative=True,
        source=KOZAI,
        link=KOZAI_LINK,
        rationale="as above",
        compute=zlk_emax_80,
    ),
    Benchmark(
        id="zlk-imin-65",
        tier="A",
        quantity="mutual inclination at maximum eccentricity, i0 = 65 deg",
        unit="deg",
        expected=float(np.degrees(np.arccos(np.sqrt(3.0 / 5.0)))),
        tolerance=0.5,
        relative=False,
        source=KOZAI,
        link=KOZAI_LINK,
        rationale="cos^2 i = 3/5 at e_max for an initially circular orbit; 0.5 deg covers "
        "sampling of the output and the finite hierarchy",
        compute=zlk_imin_65,
    ),
    Benchmark(
        id="zlk-below-critical",
        tier="A",
        quantity="maximum eccentricity below the critical angle, i0 = 30 deg",
        unit="",
        expected=0.0,
        tolerance=0.01,
        relative=False,
        source=KOZAI,
        link=KOZAI_LINK,
        rationale="no ZLK growth below 39.2 deg; e stays at its initial 0.001 plus "
        "short-period terms",
        compute=zlk_below_critical,
    ),
    Benchmark(
        id="zlk-kozai-constant",
        tier="A",
        quantity="secular drift of sqrt(1-e^2) cos i during ZLK cycles, i0 = 65 deg",
        unit="",
        expected=0.0,
        tolerance=1e-3,
        relative=False,
        source=KOZAI,
        link=KOZAI_LINK,
        rationale="exactly conserved in the averaged quadrupole problem; osculating "
        "fluctuations (~0.5 %) average out in the means over the first and last quarter",
        compute=zlk_kozai_constant,
    ),
    Benchmark(
        id="stability-unstable",
        tier="A",
        quantity="2+2 at 0.5 x Mardling-Aarseth limit changes hierarchy (1 = yes)",
        unit="",
        expected=1.0,
        tolerance=0.0,
        relative=False,
        source="Mardling & Aarseth (2001) MNRAS 321, 398",
        link="https://doi.org/10.1046/j.1365-8711.2001.03974.x",
        rationale="far inside the unstable region the hierarchy must break within 100 outer "
        "periods; applied to each binary of a 2+2 the criterion is approximate, hence the "
        "factor-2 margin",
        compute=stability_inside,
    ),
    Benchmark(
        id="stability-stable",
        tier="A",
        quantity="2+2 at 2.5 x Mardling-Aarseth limit keeps hierarchy (0 = no change)",
        unit="",
        expected=0.0,
        tolerance=0.0,
        relative=False,
        source="Mardling & Aarseth (2001) MNRAS 321, 398",
        link="https://doi.org/10.1046/j.1365-8711.2001.03974.x",
        rationale="far outside the unstable region the hierarchy must survive 100 outer "
        "periods; same margin as above",
        compute=stability_outside,
    ),
    Benchmark(
        id="naoz13-fig4-imin",
        tier="B",
        quantity="Naoz et al. 2013 Fig. 4 (PSR B1620-26 setup): minimum mutual inclination",
        unit="deg",
        expected=57.5,
        tolerance=0.5,
        relative=False,
        source="Naoz et al. (2013) MNRAS 431, 2155, Fig. 4 and text of Sect. 5.1",
        link="https://arxiv.org/abs/1107.2414",
        rationale="published curve is quadrupole order; neglected higher orders at "
        "a1/a2 ~ 0.13 are ~1-2 % of the quadrupole torque, i.e. a few 0.1 deg on a 50 deg swing",
        compute=naoz_fig4_imin,
    ),
    Benchmark(
        id="naoz13-fig4-imax",
        tier="B",
        quantity="Naoz et al. 2013 Fig. 4: maximum mutual inclination (orbit flip to retrograde)",
        unit="deg",
        expected=106.7,
        tolerance=0.5,
        relative=False,
        source="Naoz et al. (2013) MNRAS 431, 2155, Fig. 4 and text of Sect. 5.1",
        link="https://arxiv.org/abs/1107.2414",
        rationale="as above",
        compute=naoz_fig4_imax,
    ),
    Benchmark(
        id="naoz13-fig3-flip-time",
        tier="B",
        quantity="Naoz et al. 2013 Fig. 3 (eccentric ZLK): time of the first flip through 90 deg",
        unit="yr",
        expected=3.8e6,
        tolerance=0.10,
        relative=True,
        source="Naoz et al. (2013) MNRAS 431, 2155, Fig. 3 (read from the figure)",
        link="https://arxiv.org/abs/1107.2414",
        rationale="value read from the plot to about 3 %; the paper's own direct integration "
        "and octupole curves differ by several % in flip times",
        compute=naoz_fig3_flip_time,
        slow=True,
    ),
    Benchmark(
        id="naoz13-fig3-imax",
        tier="B",
        quantity="Naoz et al. 2013 Fig. 3: maximum inclination after the flip",
        unit="deg",
        expected=142.6,
        tolerance=2.0,
        relative=False,
        source="Naoz et al. (2013) MNRAS 431, 2155, Fig. 3 (read from the figure)",
        link="https://arxiv.org/abs/1107.2414",
        rationale="read from the plot to about 1 deg; spikes are sampled at 1000-yr steps",
        compute=naoz_fig3_imax,
        slow=True,
    ),
    Benchmark(
        id="quad-tp-emax",
        tier="A",
        quantity="secular quadrupole code: ZLK e_max, test particle, i0 = 65 deg",
        unit="",
        expected=_zlk_emax_formula(65.0),
        tolerance=1e-3,
        relative=True,
        source=KOZAI,
        link=KOZAI_LINK,
        rationale="the secular equations reduce exactly to the Kozai limit; e0 = 0.01 adds 4e-5",
        compute=quad_tp_emax,
    ),
    Benchmark(
        id="quad-naoz13-fig4-imin",
        tier="B",
        quantity="secular quadrupole code: Naoz et al. 2013 Fig. 4, minimum inclination",
        unit="deg",
        expected=57.5,
        tolerance=0.1,
        relative=False,
        source="Naoz et al. (2013) MNRAS 431, 2155, Fig. 4 (quadrupole order)",
        link="https://arxiv.org/abs/1107.2414",
        rationale="same equations as the paper (A27, A29, A36); the paper quotes 0.1 deg",
        compute=quad_fig4_imin,
    ),
    Benchmark(
        id="quad-naoz13-fig4-imax",
        tier="B",
        quantity="secular quadrupole code: Naoz et al. 2013 Fig. 4, maximum inclination",
        unit="deg",
        expected=106.7,
        tolerance=0.1,
        relative=False,
        source="Naoz et al. (2013) MNRAS 431, 2155, Fig. 4 (quadrupole order)",
        link="https://arxiv.org/abs/1107.2414",
        rationale="as above",
        compute=quad_fig4_imax,
    ),
    Benchmark(
        id="zlk-period-nbody",
        tier="A",
        quantity="ZLK period, N-body relative to the quadrupole secular period (i0 = 65 deg)",
        unit="",
        expected=0.0,
        tolerance=0.02,
        relative=False,
        source="Naoz et al. (2013) MNRAS 431, 2155, eqs. (A27), (A29)",
        link="https://arxiv.org/abs/1107.2414",
        rationale="a_out/a_in = 20, eps_SA ~ 0.01: non-secular corrections of order 1 %",
        compute=zlk_period_nbody_vs_quadrupole,
    ),
    Benchmark(
        id="megno-regular",
        tier="A",
        quantity="MEGNO of a regular (quasi-periodic) hierarchical triple",
        unit="",
        expected=2.0,
        tolerance=0.05,
        relative=False,
        source="Cincotta & Simo (2000) A&AS 147, 205",
        link="https://doi.org/10.1051/aas:2000108",
        rationale="MEGNO tends to 2 for quasi-periodic orbits; 2000 yr = 2000 inner orbits",
        compute=megno_regular,
    ),
    Benchmark(
        id="megno-chaotic",
        tier="A",
        quantity="MEGNO > 4 for a 2+2 inside the unstable region (1 = yes)",
        unit="",
        expected=1.0,
        tolerance=0.0,
        relative=False,
        source="Cincotta & Simo (2000) A&AS 147, 205",
        link="https://doi.org/10.1051/aas:2000108",
        rationale="MEGNO grows linearly in time for chaotic orbits",
        compute=megno_chaotic,
    ),
    Benchmark(
        id="brouwer-law",
        tier="A",
        quantity="growth exponent of the IAS15 energy error (Brouwer's law)",
        unit="",
        expected=0.5,
        tolerance=0.15,
        relative=False,
        source="Rein & Spiegel (2015) MNRAS 446, 1424",
        link="https://doi.org/10.1093/mnras/stu2164",
        rationale="round-off errors random-walk as t^0.5; 6 realisations leave ~0.1 scatter",
        compute=brouwer_exponent,
        slow=True,
    ),
    Benchmark(
        id="vok16-v994-JA",
        tier="B",
        quantity="Vokrouhlicky 2016, V994 Her setup: initial tilt J_A of binary A to the outer orbit",
        unit="deg",
        expected=11.6,
        tolerance=0.1,
        relative=False,
        source="Vokrouhlicky (2016) MNRAS 461, 3964, Sect. 2.3 (V994 Her, Zasche & Uhlar 2016)",
        link="https://doi.org/10.1093/mnras/stw1596",
        rationale="pure geometry of the sky-frame elements; the paper quotes 0.1 deg",
        compute=v994_ja,
    ),
    Benchmark(
        id="vok16-v994-JB",
        tier="B",
        quantity="Vokrouhlicky 2016, V994 Her setup: initial tilt J_B of binary B to the outer orbit",
        unit="deg",
        expected=15.8,
        tolerance=0.1,
        relative=False,
        source="Vokrouhlicky (2016) MNRAS 461, 3964, Sect. 2.3 (V994 Her, Zasche & Uhlar 2016)",
        link="https://doi.org/10.1093/mnras/stw1596",
        rationale="as above",
        compute=v994_jb,
    ),
    Benchmark(
        id="nodal-precession",
        tier="A",
        quantity="nodal precession period, circular inner orbit, i = 20 deg",
        unit="yr",
        expected=float(_nodal_precession_formula()),
        tolerance=0.01,
        relative=True,
        source="Söderhjelm (1975) A&A 42, 229 (quadrupole order)",
        link="https://ui.adsabs.harvard.edu/abs/1975A&A....42..229S",
        rationale="quadrupole formula; higher orders at a_out/a_in = 20 are below 1 %",
        compute=nodal_precession_period,
    ),
    Benchmark(
        id="polytrope-k2-n1.5",
        tier="A",
        quantity="apsidal-motion constant of a polytrope n = 1.5",
        unit="",
        expected=0.1433,
        tolerance=0.002,
        relative=True,
        source="Brooker & Olle (1955) MNRAS 115, 101",
        link="https://doi.org/10.1093/mnras/115.1.101",
        rationale="tabulated to 4 digits; the n = 1 case also agrees with the analytic "
        "(15 - pi^2) / (2 pi^2) to 1e-9 (test suite)",
        compute=lambda: polytrope_k2(1.5),
    ),
    Benchmark(
        id="polytrope-k2-n3",
        tier="A",
        quantity="apsidal-motion constant of a polytrope n = 3",
        unit="",
        expected=0.0144,
        tolerance=0.005,
        relative=True,
        source="Brooker & Olle (1955) MNRAS 115, 101",
        link="https://doi.org/10.1093/mnras/115.1.101",
        rationale="tabulated value 0.0144; 0.5 % covers the rounding to three digits",
        compute=lambda: polytrope_k2(3.0),
    ),
    Benchmark(
        id="roche-eggleton-q1",
        tier="A",
        quantity="Roche-lobe radius / separation, q = 1 (Eggleton formula vs exact)",
        unit="",
        expected=roche_volume_radius(1.0),
        tolerance=0.01,
        relative=True,
        source="Eggleton (1983) ApJ 268, 368",
        link="https://doi.org/10.1086/160960",
        rationale="the formula is a fit to the exact volume-equivalent radius, stated to be "
        "better than 1 % for all q; the expected value is integrated here directly",
        compute=eggleton_q1,
    ),
    Benchmark(
        id="roche-eggleton-q0.1",
        tier="A",
        quantity="Roche-lobe radius / separation, q = 0.1 (Eggleton formula vs exact)",
        unit="",
        expected=roche_volume_radius(0.1),
        tolerance=0.01,
        relative=True,
        source="Eggleton (1983) ApJ 268, 368",
        link="https://doi.org/10.1086/160960",
        rationale="as above",
        compute=eggleton_q01,
    ),
    Benchmark(
        id="apsidal-sterne-static",
        tier="A",
        quantity="apsidal period from the tidal distortion, e = 0.3, no rotation",
        unit="yr",
        expected=float(_sterne_expected(0.3, False)),
        tolerance=0.005,
        relative=True,
        source="Sterne (1939) MNRAS 99, 451",
        link="https://doi.org/10.1093/mnras/99.5.451",
        rationale="first-order theory in (R/a)^5 with R/a ~ 0.14; the N-body value also "
        "contains short-period terms, which average out over the run",
        compute=apsidal_period_static,
    ),
    Benchmark(
        id="apsidal-sterne-rotating",
        tier="A",
        quantity="apsidal period from tidal + rotational distortion, e = 0.3, Omega = n",
        unit="yr",
        expected=float(_sterne_expected(0.3, True)),
        tolerance=0.005,
        relative=True,
        source="Sterne (1939) MNRAS 99, 451",
        link="https://doi.org/10.1093/mnras/99.5.451",
        rationale="as above; spins aligned with the orbit",
        compute=apsidal_period_rotating,
    ),
    Benchmark(
        id="tides-leconte-a",
        tier="A",
        quantity="tidal da/dt, N-body / theory, e = 0.3, non-rotating stars",
        unit="",
        expected=1.0,
        tolerance=0.01,
        relative=True,
        source="Leconte et al. (2010) A&A 516, A64; Hut (1981) A&A 99, 126",
        link="https://doi.org/10.1051/0004-6361/201014337",
        rationale="orbit-averaged constant-time-lag theory; the N-body rate is the difference "
        "of runs with and without dissipation over 60 orbits",
        compute=lambda: tidal_ratio(0.3, 0.0, 0),
    ),
    Benchmark(
        id="tides-leconte-e",
        tier="A",
        quantity="tidal de/dt, N-body / theory, e = 0.3, non-rotating stars",
        unit="",
        expected=1.0,
        tolerance=0.01,
        relative=True,
        source="Leconte et al. (2010) A&A 516, A64; Hut (1981) A&A 99, 126",
        link="https://doi.org/10.1051/0004-6361/201014337",
        rationale="as above",
        compute=lambda: tidal_ratio(0.3, 0.0, 1),
    ),
    Benchmark(
        id="tides-leconte-spin",
        tier="A",
        quantity="tidal spin-up dOmega/dt, N-body / theory, e = 0.6, Omega = 1.2 n",
        unit="",
        expected=1.0,
        tolerance=0.01,
        relative=True,
        source="Leconte et al. (2010) A&A 516, A64; Hut (1981) A&A 99, 126",
        link="https://doi.org/10.1051/0004-6361/201014337",
        rationale="as above, with a fast, eccentric case",
        compute=lambda: tidal_ratio(0.6, 1.2, 2),
    ),
    Benchmark(
        id="tides-leconte-e-eccentric",
        tier="A",
        quantity="tidal de/dt, N-body / theory, e = 0.6, Omega = 1.2 n",
        unit="",
        expected=1.0,
        tolerance=0.01,
        relative=True,
        source="Leconte et al. (2010) A&A 516, A64; Hut (1981) A&A 99, 126",
        link="https://doi.org/10.1051/0004-6361/201014337",
        rationale="as above",
        compute=lambda: tidal_ratio(0.6, 1.2, 1),
    ),
    Benchmark(
        id="secular-kozai-emax",
        tier="A",
        quantity="secular backend: ZLK maximum eccentricity, i0 = 65 deg, test particle",
        unit="",
        expected=_zlk_emax_formula(65.0),
        tolerance=0.005,
        relative=True,
        source=KOZAI,
        link=KOZAI_LINK,
        rationale="the secular equations contain no short-period terms, so the quadrupole "
        "test-particle limit is reproduced more closely than by N-body",
        compute=secular_kozai_emax,
    ),
    Benchmark(
        id="secular-nodal-precession",
        tier="A",
        quantity="secular backend: nodal precession period, i = 20 deg [yr]",
        unit="yr",
        expected=float(_nodal_precession_formula()),
        tolerance=0.005,
        relative=True,
        source="Söderhjelm (1975) A&A 42, 229 (quadrupole order)",
        link="https://ui.adsabs.harvard.edu/abs/1975A&A....42..229S",
        rationale="same quadrupole formula as nodal-precession; the secular model contains "
        "exactly this order",
        compute=secular_nodal_precession,
    ),
    Benchmark(
        id="secular-gr-apsidal",
        tier="A",
        quantity="secular backend: relativistic apsidal period, a = 0.05 AU, e = 0.3",
        unit="yr",
        expected=_gr_apsidal_formula(),
        tolerance=0.005,
        relative=True,
        source="Will (2014) Living Rev. Relativ. 17, 4",
        link="https://doi.org/10.12942/lrr-2014-4",
        rationale="1PN formula; the far perturber adds less than 0.1 %",
        compute=secular_gr_apsidal,
    ),
    Benchmark(
        id="vok16-jab-max-analytic",
        tier="B",
        quantity="Vokrouhlický 2016, analytic limit J_AB^max for V994 Her (Fig. 2 setup)",
        unit="deg",
        expected=5.2,
        tolerance=0.05,
        relative=False,
        source="Vokrouhlicky (2016) MNRAS 461, 3964, Sect. 2.3, eqs. (33)-(36)",
        link="https://doi.org/10.1093/mnras/stw1596",
        rationale="the paper quotes about 5.2 deg for this setup",
        compute=lambda: _secular_v994(2, 1.0)[2],
    ),
    Benchmark(
        id="vok16-fig2-jab-min",
        tier="B",
        quantity="secular backend: Vokrouhlický 2016 Fig. 2, minimum J_AB over 100 kyr",
        unit="deg",
        expected=3.5,
        tolerance=0.15,
        relative=False,
        source="Vokrouhlicky (2016) MNRAS 461, 3964, Fig. 2",
        link="https://doi.org/10.1093/mnras/stw1596",
        rationale="value quoted in the text to 0.1 deg; swingersim N-body gives 3.41 deg",
        compute=lambda: _secular_v994(2, 100_000.0)[0],
    ),
    Benchmark(
        id="vok16-fig2-jab-max",
        tier="B",
        quantity="secular backend: Vokrouhlický 2016 Fig. 2, maximum J_AB over 100 kyr",
        unit="deg",
        expected=5.2,
        tolerance=0.15,
        relative=False,
        source="Vokrouhlicky (2016) MNRAS 461, 3964, Fig. 2",
        link="https://doi.org/10.1093/mnras/stw1596",
        rationale="value quoted in the text to 0.1 deg; swingersim N-body gives 5.19 deg",
        compute=lambda: _secular_v994(2, 100_000.0)[1],
    ),
    Benchmark(
        id="vok16-fig2-jab-1myr",
        tier="B",
        quantity="secular backend: Fig. 2 setup, maximum J_AB over 1 Myr",
        unit="deg",
        expected=5.2,
        tolerance=0.15,
        relative=False,
        source="Vokrouhlicky (2016) MNRAS 461, 3964, Sect. 2.3 (stable over 100 Myr)",
        link="https://doi.org/10.1093/mnras/stw1596",
        rationale="the paper finds no long-term growth of J_AB beyond the analytic limit "
        "over 100 Myr; 1 Myr keeps the benchmark at a few minutes",
        compute=lambda: _secular_v994(2, 1_000_000.0)[1],
        slow=True,
    ),
    Benchmark(
        id="vok16-fig3-jab-max",
        tier="B",
        quantity="secular backend: Vokrouhlický 2016 Fig. 3, maximum J_AB over 20 kyr",
        unit="deg",
        expected=6.0,
        tolerance=0.5,
        relative=False,
        source="Vokrouhlicky (2016) MNRAS 461, 3964, Fig. 3",
        link="https://doi.org/10.1093/mnras/stw1596",
        rationale="the text says 'between about 1 and about 6 deg'; swingersim N-body gives 6.30",
        compute=lambda: _secular_v994(3, 20_000.0)[1],
    ),
    Benchmark(
        id="vok16-eclipse-prob-A",
        tier="B",
        quantity="V994 Her: chance that binary A eclipses, isotropic orientation",
        unit="",
        expected=0.366,
        tolerance=0.02,
        relative=True,
        source="Vokrouhlicky (2016) MNRAS 461, 3964, Sect. 3.3",
        link="https://doi.org/10.1093/mnras/stw1596",
        rationale="the paper quotes 36.6 % with radii 'approximately' 2.4 and 2.1 Rsun; the "
        "rounding of the radii alone changes the value by about 1 %",
        compute=lambda: v994_eclipse_probability("A"),
    ),
    Benchmark(
        id="vok16-eclipse-prob-B",
        tier="B",
        quantity="V994 Her: chance that binary B eclipses, isotropic orientation",
        unit="",
        expected=0.397,
        tolerance=0.02,
        relative=True,
        source="Vokrouhlicky (2016) MNRAS 461, 3964, Sect. 3.3",
        link="https://doi.org/10.1093/mnras/stw1596",
        rationale="as above, radii 1.6 and 1.7 Rsun",
        compute=lambda: v994_eclipse_probability("B"),
    ),
    Benchmark(
        id="secular-3p1-kozai-emax",
        tier="A",
        quantity="3+1 secular backend: ZLK e_max of the inner orbit, test particle, i0 = 65 deg",
        unit="",
        expected=_zlk_emax_formula(65.0),
        tolerance=0.005,
        relative=True,
        source=KOZAI,
        link=KOZAI_LINK,
        rationale="the nested secular equations reduce to the Kozai limit; e0 = 0.01",
        compute=secular_3p1_kozai_emax,
    ),
    Benchmark(
        id="secular-naoz13-fig3-flip-time",
        tier="B",
        quantity="3+1 secular backend: Naoz et al. 2013 Fig. 3, time of the first flip [yr]",
        unit="yr",
        expected=3.8e6,
        tolerance=0.10,
        relative=True,
        source="Naoz et al. (2013) MNRAS 431, 2155, Fig. 3 (read from the figure)",
        link="https://arxiv.org/abs/1107.2414",
        rationale="as for the N-body benchmark naoz13-fig3-flip-time; tests the octupole terms",
        compute=secular_naoz_fig3_flip_time,
    ),
    Benchmark(
        id="secular-naoz13-fig3-imax",
        tier="B",
        quantity="3+1 secular backend: Naoz et al. 2013 Fig. 3, maximum inclination after the flip",
        unit="deg",
        expected=142.6,
        tolerance=2.0,
        relative=False,
        source="Naoz et al. (2013) MNRAS 431, 2155, Fig. 3 (read from the figure)",
        link="https://arxiv.org/abs/1107.2414",
        rationale="as for naoz13-fig3-imax",
        compute=secular_naoz_fig3_imax,
    ),
    Benchmark(
        id="gas-migration-a",
        tier="A",
        quantity="gas phase: semi-major axis after 3000 yr with tau_a = 5000 yr [AU]",
        unit="AU",
        expected=float(10.0 * np.exp(-0.6)),
        tolerance=1e-6,
        relative=True,
        source="definition of the parametric gas phase: da/dt = -a/tau_a (exponential migration, "
        "e.g. Papaloizou & Larwood 2000)",
        link="https://doi.org/10.1046/j.1365-8711.2000.03466.x",
        rationale="analytic solution; only the integration error remains",
        compute=lambda: _gas_run()[0],
    ),
    Benchmark(
        id="gas-damping-e",
        tier="A",
        quantity="gas phase: eccentricity after 3000 yr with tau_e = 2000 yr",
        unit="",
        expected=float(0.5 * np.exp(-1.5)),
        tolerance=1e-6,
        relative=True,
        source="definition of the parametric gas phase: de/dt = -e/tau_e",
        link="https://doi.org/10.1046/j.1365-8711.2000.03466.x",
        rationale="analytic solution; coplanar-free outer orbit of equal-mass binaries has no "
        "secular eccentricity change of its own",
        compute=lambda: _gas_run()[1],
    ),
    Benchmark(
        id="gas-alignment-i",
        tier="A",
        quantity="gas phase: inclination to the disc after 3000 yr, i0 = 40 deg, tau_i = 1500 yr",
        unit="deg",
        expected=float(np.degrees(2.0 * np.arctan(np.tan(np.radians(20.0)) * np.exp(-2.0)))),
        tolerance=0.01,
        relative=False,
        source="definition of the parametric gas phase: dk/dt = (n - (n.k) k)/tau_i, i.e. "
        "tan(i/2) = tan(i0/2) exp(-t/tau_i)",
        link="https://doi.org/10.1046/j.1365-8711.2000.03466.x",
        rationale="analytic solution; the tight binaries make the outer orbit precess by less "
        "than 0.01 deg",
        compute=lambda: _gas_run()[2],
    ),
    Benchmark(
        id="scatter-hyperbola-deflection",
        tier="A",
        quantity="scattering set-up: deflection angle of a 1 + 0.5 Msun hyperbolic passage, "
        "v_inf = 5 km/s, b = 3 AU",
        unit="deg",
        expected=float(
            np.degrees(
                2
                * np.arcsin(
                    1
                    / np.sqrt(
                        1
                        + (
                            (5.0 * units.SECONDS_PER_YEAR / units.AU_KM) ** 2
                            * 3.0
                            / (units.G * 1.5)
                        )
                        ** 2
                    )
                )
            )
        ),
        tolerance=1e-6,
        relative=False,
        source="two-body problem: deflection 2 arcsin(1/e), e = sqrt(1 + b^2 v_inf^4 / (G M)^2)",
        link="https://doi.org/10.1093/mnras/stv2938",
        rationale="analytic; IAS15 conserves the Laplace-Runge-Lenz vector to ~1e-12",
        compute=lambda: scatter_deflection()[0],
    ),
    Benchmark(
        id="at16-binbin-single-ionization",
        tier="B",
        quantity="binary-binary scattering, equal masses, a1 = a2, e = 0, v_inf = v_crit: "
        "cross-section of 'single ionization' / (pi a^2)",
        unit="",
        expected=2.645,
        tolerance=0.1,
        relative=True,
        source="Antognini & Thompson (2016) MNRAS 456, 4219, Table 3 (sum of the rows of the class)",
        link="https://doi.org/10.1093/mnras/stv2938",
        rationale="1-sigma statistics 4 % (12 000 encounters vs 10^6); the paper writes the normalisation as "
        "pi (a1^2 + a2^2), its numbers agree with pi a^2 of one binary (factor-2 ambiguity noted)",
        compute=lambda: _at16_run().get("single ionization", 0.0),
        slow=True,
    ),
    Benchmark(
        id="at16-binbin-exchange-ionization",
        tier="B",
        quantity="binary-binary scattering, equal masses, a1 = a2, e = 0, v_inf = v_crit: "
        "cross-section of 'exchange + ionization' / (pi a^2)",
        unit="",
        expected=3.259,
        tolerance=0.1,
        relative=True,
        source="Antognini & Thompson (2016) MNRAS 456, 4219, Table 3 (sum of the rows of the class)",
        link="https://doi.org/10.1093/mnras/stv2938",
        rationale="1-sigma statistics 4 % (12 000 encounters vs 10^6); the paper writes the normalisation as "
        "pi (a1^2 + a2^2), its numbers agree with pi a^2 of one binary (factor-2 ambiguity noted)",
        compute=lambda: _at16_run().get("exchange + ionization", 0.0),
        slow=True,
    ),
    Benchmark(
        id="at16-binbin-exchange",
        tier="B",
        quantity="binary-binary scattering, equal masses, a1 = a2, e = 0, v_inf = v_crit: "
        "cross-section of 'exchange' / (pi a^2)",
        unit="",
        expected=0.302,
        tolerance=0.3,
        relative=True,
        source="Antognini & Thompson (2016) MNRAS 456, 4219, Table 3 (sum of the rows of the class)",
        link="https://doi.org/10.1093/mnras/stv2938",
        rationale="1-sigma statistics 10 % here and 3 % there (12 000 encounters vs 10^6); the paper writes the normalisation as "
        "pi (a1^2 + a2^2), its numbers agree with pi a^2 of one binary (factor-2 ambiguity noted)",
        compute=lambda: _at16_run().get("exchange", 0.0),
        slow=True,
    ),
    Benchmark(
        id="mse-3p1-inner-emax",
        tier="C",
        quantity="3+1 (a = 1, 15, 300 AU, i_mid = 85 deg): maximum e of the inner binary in 1e5 yr",
        unit="",
        expected=0.9707,
        tolerance=0.002,
        relative=False,
        source="MSE (Hamers et al. 2021, MNRAS 502, 4479), github.com/hamers/mse commit 5095596, run by tools/compare_mse.py",
        link="https://doi.org/10.1093/mnras/stab287",
        rationale="same Hamiltonian (binary-pair quadrupole + octupole, 1PN); MSE sampled every "
        "12.5 / 10 yr; MSE's extra terms (binary-triplet, hexadecupole) change e_max by < 0.001",
        compute=lambda: mse_case("3p1")["inner"][0],
    ),
    Benchmark(
        id="mse-3p1-inner-t-emax",
        tier="C",
        quantity="3+1: time of the maximum e of the inner binary",
        unit="yr",
        expected=49325.0,
        tolerance=0.002,
        relative=True,
        source="MSE (Hamers et al. 2021, MNRAS 502, 4479), github.com/hamers/mse commit 5095596, run by tools/compare_mse.py",
        link="https://doi.org/10.1093/mnras/stab287",
        rationale="same Hamiltonian (binary-pair quadrupole + octupole, 1PN); MSE sampled every "
        "12.5 / 10 yr; MSE's extra terms (binary-triplet, hexadecupole) change e_max by < 0.001",
        compute=lambda: mse_case("3p1")["inner"][1],
    ),
    Benchmark(
        id="mse-3p1-mid-emax",
        tier="C",
        quantity="3+1: maximum e of the middle orbit",
        unit="",
        expected=0.2245,
        tolerance=0.002,
        relative=False,
        source="MSE (Hamers et al. 2021, MNRAS 502, 4479), github.com/hamers/mse commit 5095596, run by tools/compare_mse.py",
        link="https://doi.org/10.1093/mnras/stab287",
        rationale="same Hamiltonian (binary-pair quadrupole + octupole, 1PN); MSE sampled every "
        "12.5 / 10 yr; MSE's extra terms (binary-triplet, hexadecupole) change e_max by < 0.001",
        compute=lambda: mse_case("3p1")["mid"][0],
    ),
    Benchmark(
        id="mse-2p2-A-emax",
        tier="C",
        quantity="2+2 (a = 1, 1.5, 40 AU): maximum e of binary A in 8e4 yr",
        unit="",
        expected=0.9494,
        tolerance=0.002,
        relative=False,
        source="MSE (Hamers et al. 2021, MNRAS 502, 4479), github.com/hamers/mse commit 5095596, run by tools/compare_mse.py",
        link="https://doi.org/10.1093/mnras/stab287",
        rationale="same Hamiltonian (binary-pair quadrupole + octupole, 1PN); MSE sampled every "
        "12.5 / 10 yr; MSE's extra terms (binary-triplet, hexadecupole) change e_max by < 0.001",
        compute=lambda: mse_case("2p2")["A"][0],
    ),
    Benchmark(
        id="mse-2p2-B-emax",
        tier="C",
        quantity="2+2: maximum e of binary B in 8e4 yr",
        unit="",
        expected=0.7758,
        tolerance=0.002,
        relative=False,
        source="MSE (Hamers et al. 2021, MNRAS 502, 4479), github.com/hamers/mse commit 5095596, run by tools/compare_mse.py",
        link="https://doi.org/10.1093/mnras/stab287",
        rationale="same Hamiltonian (binary-pair quadrupole + octupole, 1PN); MSE sampled every "
        "12.5 / 10 yr; MSE's extra terms (binary-triplet, hexadecupole) change e_max by < 0.001",
        compute=lambda: mse_case("2p2")["B"][0],
    ),
    Benchmark(
        id="mse-3p1-rlof-time",
        tier="C",
        quantity="3+1 (a = 0.3, 5, 200 AU, i = 89 deg): onset of Roche-lobe overflow of the "
        "inner binary driven by ZLK cycles",
        unit="yr",
        expected=380.03,
        tolerance=0.002,
        relative=True,
        source="MSE (Hamers et al. 2021, MNRAS 502, 4479), commit 5095596: time of the "
        "'RLOF start' event; tools/compare_mse.py",
        link="https://doi.org/10.1093/mnras/stab287",
        rationale="same radii (MSE's SSE ZAMS radii), Eggleton Roche lobe at pericentre; swingersim "
        "sampled every 0.1 yr (380.10 yr; 380.04 yr with 0.01 yr)",
        compute=mse_rlof_time,
    ),
    Benchmark(
        id="chaosmap-circumbinary-limit",
        tier="B",
        quantity="MEGNO map: inner stability limit of a circumbinary orbit, equal-mass "
        "circular binary, all of 8 phases regular [a_bin]",
        unit="",
        expected=2.27,
        tolerance=0.06,
        relative=False,
        source="Quarles et al. (2018) ApJ 856, 150, published grid a_crit.txt (mu = 0.5, e = 0); "
        "Holman & Wiegert (1999) fit gives 2.39",
        link="https://github.com/saturnaxis/CBP_stability",
        rationale="grid step 0.05, 8 phases and 1e4 binary orbits instead of finer grids and "
        "1e5 orbits",
        compute=chaosmap_circumbinary,
        slow=True,
    ),
    Benchmark(
        id="escape-theory-energy-slope-L0",
        tier="A",
        quantity="statistical escape theory, equal masses, L0 = 0: slope of dP/dE_B for "
        "|E_B| >> |E0|",
        unit="",
        expected=-3.5,
        tolerance=0.05,
        relative=False,
        source="Stone & Leigh (2019) Nature 576, 406" + ": dsigma/dE_B ~ |E_B|^(-7/2) for L0 = 0",
        link="https://www.nature.com/articles/s41586-019-1833-8",
        rationale="asymptotic limit of their eq. 2; fitted between 20 and 300 |E0|",
        compute=lambda: escape_energy_slope(0.0, 20.0, 300.0),
    ),
    Benchmark(
        id="escape-theory-energy-slope-highL",
        tier="A",
        quantity="statistical escape theory, L0 = L_c: slope of dP/dE_B",
        unit="",
        expected=-4.0,
        tolerance=0.1,
        relative=False,
        source="Stone & Leigh (2019) Nature 576, 406" + ": dsigma/dE_B ~ |E_B|^(-4) for L0 ~ L_c",
        link="https://www.nature.com/articles/s41586-019-1833-8",
        rationale="limit of their eq. 2; the range ends at the cut-off |E_max| ~ L0^-2",
        compute=lambda: escape_energy_slope(1.0, 1.5, 3.0),
    ),
    Benchmark(
        id="escape-theory-mean-e-L0",
        tier="A",
        quantity="statistical escape theory, L0 = 0: mean eccentricity of the binary",
        unit="",
        expected=0.8133,
        tolerance=0.01,
        relative=False,
        source="Stone & Leigh (2019) Nature 576, 406"
        + ", eq. 24: dsigma/de ~ e (1 + e) / sqrt(1 - e^2)",
        link="https://www.nature.com/articles/s41586-019-1833-8",
        rationale="mean of eq. 24; grid of 200 eccentricities",
        compute=lambda: _escape(0.0).mean_e,
    ),
    Benchmark(
        id="escape-theory-mean-e-highL",
        tier="A",
        quantity="statistical escape theory, L0 = L_c: mean eccentricity of the binary",
        unit="",
        expected=0.70,
        tolerance=0.02,
        relative=False,
        source="Stone & Leigh (2019) Nature 576, 406" + ", eq. 23: dsigma/de = 6/5 e (1 + e)",
        link="https://www.nature.com/articles/s41586-019-1833-8",
        rationale="eq. 23 is the high-angular-momentum approximation",
        compute=lambda: _escape(1.0).mean_e,
    ),
    Benchmark(
        id="timing-ltte-amplitude",
        tier="A",
        quantity="synthetic O-C: light-time semi-amplitude of binary A, circular edge-on outer "
        "orbit, a = 10 AU, equal masses",
        unit="min",
        expected=float(5.0 / units.C * units.DAYS_PER_YEAR * 1440.0),
        tolerance=0.001,
        relative=True,
        source="Irwin (1952) ApJ 116, 211: semi-amplitude a_12 sin i / c (circular orbit)",
        link="https://doi.org/10.1086/145604",
        rationale="analytic; the dynamical delays of this wide system are below 0.01 min",
        compute=timing_ltte_amplitude,
    ),
    Benchmark(
        id="phenom-romer-delay",
        tier="A",
        quantity="phenomenological minima: delay of the secondary minimum after half the "
        "period, 2 + 1 Msun ZAMS, P = 5 d, circular, i = 89 deg",
        unit="s",
        expected=phenom_romer_delay_expected(),
        tolerance=0.015,
        relative=True,
        source="Kaplan (2010) ApJ 717, L108: Romer delay (P / pi c)(K2 - K1) of the secondary "
        "eclipse",
        link="https://doi.org/10.1088/2041-8205/717/2/L108",
        rationale="analytic for point masses; the minima are measured on the synthetic light "
        "curve with the Mikulasek (2015) profile, which leaves a -0.5 % bias of the finite "
        "eclipse and the fitted profile",
        compute=phenom_romer_delay,
    ),
    Benchmark(
        id="phoebe-sphere-vs-tier-a",
        tier="C",
        quantity="largest difference of the normalised light curves, tier A vs PHOEBE 2.5 "
        "with the same physics (spheres, black bodies, Claret limb darkening), e = 0.1",
        unit="",
        expected=0.0,
        tolerance=1e-3,
        relative=False,
        source="PHOEBE 2 (Prsa et al. 2016, ApJS 227, 29; Conroy et al. 2020, ApJS 250, 34) "
        "run by the swingersim-PHOEBE hybrid on the same orbit",
        link="https://doi.org/10.3847/0067-0049/227/2/29",
        rationale="PHOEBE integrates over a triangle mesh (1500 per star), tier A exactly; the "
        "mesh changes the flux by a few 1e-4 in ingress and egress. Checks the mapping of the "
        "elements (incl = i, per0 = omega + 180 deg, time of periastron)",
        compute=phoebe_sphere_vs_tier_a,
        manual=True,
    ),
    Benchmark(
        id="phoebe-romer-delay",
        tier="A",
        quantity="PHOEBE hybrid: delay of the secondary minimum after half the period, 2 + 1 "
        "Msun, P = 5 d, circular (with minus without light time)",
        unit="s",
        expected=phenom_romer_delay_expected(),
        tolerance=0.02,
        relative=True,
        source="Kaplan (2010) ApJ 717, L108: Romer delay (P / pi c)(K2 - K1)",
        link="https://doi.org/10.1088/2041-8205/717/2/L108",
        rationale="checks the light-time effect of PHOEBE together with the barycentric delay "
        "of the hybrid; PHOEBE's mesh shifts single minima by a few s at 1500 triangles, the "
        "difference with and without light time removes this bias",
        compute=phoebe_romer_delay,
        manual=True,
    ),
    Benchmark(
        id="qzcar-broz-a1",
        tier="B",
        quantity="QZ Car: semi-major axis of Ac1+Ac2 (P = 6.0 d) from the Table 2 elements and masses",
        unit="Rsun",
        expected=QZCAR_TABLE3["a1"],
        tolerance=0.002,
        relative=True,
        source="Brož et al. (2022) A&A 666, A24, Tables 2 and 3",
        link="https://doi.org/10.1051/0004-6361/202243596",
        rationale="checks the hierarchical initialisation with a published four-body model; the "
        "masses of Table 2 are rounded to three digits (0.1 % in a)",
        compute=lambda: qzcar_semi_major_axis("A"),
    ),
    Benchmark(
        id="qzcar-broz-a2",
        tier="B",
        quantity="QZ Car: semi-major axis of Aa1+Aa2 (P = 20.7 d) from the Table 2 elements and masses",
        unit="Rsun",
        expected=QZCAR_TABLE3["a2"],
        tolerance=0.002,
        relative=True,
        source="Brož et al. (2022) A&A 666, A24, Tables 2 and 3",
        link="https://doi.org/10.1051/0004-6361/202243596",
        rationale="checks the hierarchical initialisation with a published four-body model; the "
        "masses of Table 2 are rounded to three digits (0.1 % in a)",
        compute=lambda: qzcar_semi_major_axis("B"),
    ),
    Benchmark(
        id="qzcar-broz-a3",
        tier="B",
        quantity="QZ Car: semi-major axis of Ac+Aa (P = 14 728 d) from the Table 2 elements and masses",
        unit="Rsun",
        expected=QZCAR_TABLE3["a3"],
        tolerance=0.002,
        relative=True,
        source="Brož et al. (2022) A&A 666, A24, Tables 2 and 3",
        link="https://doi.org/10.1051/0004-6361/202243596",
        rationale="checks the hierarchical initialisation with a published four-body model; the "
        "masses of Table 2 are rounded to three digits (0.1 % in a)",
        compute=lambda: qzcar_semi_major_axis("out"),
    ),
    Benchmark(
        id="qzcar-broz-lv-ac1",
        tier="C",
        quantity="QZ Car: share of Ac1 in the V-band light from its radius and Teff",
        unit="",
        expected=QZCAR_TABLE3["L_V"]["m1"],
        tolerance=0.005,
        relative=False,
        source="Brož et al. (2022) A&A 666, A24, Table 3 (synthetic spectra of the Xitau model)",
        link="https://doi.org/10.1051/0004-6361/202243596",
        rationale="black bodies with Claret / TLUSTY limb darkening against their synthetic "
        "spectra; the four stars have similar Teff, so the ratios are robust",
        compute=lambda: qzcar_light_fraction("m1"),
    ),
    Benchmark(
        id="qzcar-broz-lv-ac2",
        tier="C",
        quantity="QZ Car: share of Ac2 in the V-band light from its radius and Teff",
        unit="",
        expected=QZCAR_TABLE3["L_V"]["m2"],
        tolerance=0.005,
        relative=False,
        source="Brož et al. (2022) A&A 666, A24, Table 3 (synthetic spectra of the Xitau model)",
        link="https://doi.org/10.1051/0004-6361/202243596",
        rationale="black bodies with Claret / TLUSTY limb darkening against their synthetic "
        "spectra; the four stars have similar Teff, so the ratios are robust",
        compute=lambda: qzcar_light_fraction("m2"),
    ),
    Benchmark(
        id="qzcar-broz-lv-aa1",
        tier="C",
        quantity="QZ Car: share of Aa1 in the V-band light from its radius and Teff",
        unit="",
        expected=QZCAR_TABLE3["L_V"]["m3"],
        tolerance=0.005,
        relative=False,
        source="Brož et al. (2022) A&A 666, A24, Table 3 (synthetic spectra of the Xitau model)",
        link="https://doi.org/10.1051/0004-6361/202243596",
        rationale="black bodies with Claret / TLUSTY limb darkening against their synthetic "
        "spectra; the four stars have similar Teff, so the ratios are robust",
        compute=lambda: qzcar_light_fraction("m3"),
    ),
    Benchmark(
        id="qzcar-broz-lv-aa2",
        tier="C",
        quantity="QZ Car: share of Aa2 in the V-band light from its radius and Teff",
        unit="",
        expected=QZCAR_TABLE3["L_V"]["m4"],
        tolerance=0.005,
        relative=False,
        source="Brož et al. (2022) A&A 666, A24, Table 3 (synthetic spectra of the Xitau model)",
        link="https://doi.org/10.1051/0004-6361/202243596",
        rationale="black bodies with Claret / TLUSTY limb darkening against their synthetic "
        "spectra; the four stars have similar Teff, so the ratios are robust",
        compute=lambda: qzcar_light_fraction("m4"),
    ),
    Benchmark(
        id="qzcar-evolutionary-mass-ac1",
        tier="D",
        quantity="QZ Car: present mass of the single-star MIST model with the radius and Teff "
        "of Ac1, against its dynamical mass",
        unit="Msun",
        expected=25.5,
        tolerance=0.15,
        relative=True,
        source="Brož et al. (2022) A&A 666, A24, Tables 2 and 3 (dynamical mass 25.5 +- 2)",
        link="https://doi.org/10.1051/0004-6361/202243596",
        rationale="evolutionary and dynamical masses of O stars agree to 10-20 % for single "
        "stars; MIST v1.2, v/vcrit = 0; the age of this model is 5.6 Myr",
        compute=lambda: qzcar_evolutionary_mass("m1"),
    ),
    Benchmark(
        id="qzcar-evolutionary-mass-ac2",
        tier="D",
        quantity="the same for Ac2",
        unit="Msun",
        expected=33.2,
        tolerance=0.15,
        relative=True,
        source="Brož et al. (2022) A&A 666, A24, Tables 2 and 3 (dynamical mass 33.2 +- 2)",
        link="https://doi.org/10.1051/0004-6361/202243596",
        rationale="as for Ac1",
        compute=lambda: qzcar_evolutionary_mass("m2"),
        explanation="Ac2 rotates at 410 km/s (fitted, Broz et al. 2022, Table 2) – the likely "
        "mass gainer of the Ac1 -> Ac2 transfer, which looks younger and lighter on single-star "
        "tracks; its MIST mass (26 Msun) and age (5.1 Myr) are then only indicative",
    ),
    Benchmark(
        id="qzcar-evolutionary-mass-aa1",
        tier="D",
        quantity="the same for Aa1",
        unit="Msun",
        expected=69.8,
        tolerance=0.15,
        relative=True,
        source="Brož et al. (2022) A&A 666, A24, Tables 2 and 3 (dynamical mass 69.8 +- 5)",
        link="https://doi.org/10.1051/0004-6361/202243596",
        rationale="as for Ac1",
        compute=lambda: qzcar_evolutionary_mass("m3"),
        explanation="the luminosity of Aa1 fits a 43 Msun single star of 3.9 Myr, younger than "
        "Ac1 (5.6 Myr); Brož et al. note that Aa1 is extremely massive compared with Aa2 and "
        "suggest past mass transfer Aa2 -> Aa1 – a rejuvenated, overmassive gainer, which single-"
        "star tracks cannot describe",
    ),
    Benchmark(
        id="agb-final-mass-1p5",
        tier="D",
        quantity="white-dwarf mass (core at the end of the AGB) of a 1.5 Msun MIST star (solar [Fe/H]) against "
        "the white-dwarf initial-final mass relation",
        unit="Msun",
        expected=cummings2018_ifmr(1.5),
        tolerance=0.1,
        relative=False,
        source="Cummings et al. (2018) ApJ 866, 21, eqs. 4-6 (MIST-based IFMR from 80 white "
        "dwarfs in 13 clusters)",
        link="https://arxiv.org/abs/1809.01673",
        rationale="0.1 Msun: the observed scatter (0.06 Msun) plus the known offset of the "
        "models, which predict white dwarfs ~0.1 Msun lighter than observed for 3-6 Msun "
        "(Cummings et al. 2018, Sect. 5); checks the AGB extension of the MIST grids",
        compute=lambda: mist_final_mass(1.5),
    ),
    Benchmark(
        id="agb-final-mass-2p5",
        tier="D",
        quantity="white-dwarf mass (core at the end of the AGB) of a 2.5 Msun MIST star (solar [Fe/H]) against "
        "the white-dwarf initial-final mass relation",
        unit="Msun",
        expected=cummings2018_ifmr(2.5),
        tolerance=0.1,
        relative=False,
        source="Cummings et al. (2018) ApJ 866, 21, eqs. 4-6 (MIST-based IFMR from 80 white "
        "dwarfs in 13 clusters)",
        link="https://arxiv.org/abs/1809.01673",
        rationale="0.1 Msun: the observed scatter (0.06 Msun) plus the known offset of the "
        "models, which predict white dwarfs ~0.1 Msun lighter than observed for 3-6 Msun "
        "(Cummings et al. 2018, Sect. 5); checks the AGB extension of the MIST grids",
        compute=lambda: mist_final_mass(2.5),
    ),
    Benchmark(
        id="agb-final-mass-3p0",
        tier="D",
        quantity="white-dwarf mass (core at the end of the AGB) of a 3.0 Msun MIST star (solar [Fe/H]) against "
        "the white-dwarf initial-final mass relation",
        unit="Msun",
        expected=cummings2018_ifmr(3.0),
        tolerance=0.1,
        relative=False,
        source="Cummings et al. (2018) ApJ 866, 21, eqs. 4-6 (MIST-based IFMR from 80 white "
        "dwarfs in 13 clusters)",
        link="https://arxiv.org/abs/1809.01673",
        rationale="0.1 Msun: the observed scatter (0.06 Msun) plus the known offset of the "
        "models, which predict white dwarfs ~0.1 Msun lighter than observed for 3-6 Msun "
        "(Cummings et al. 2018, Sect. 5); checks the AGB extension of the MIST grids",
        compute=lambda: mist_final_mass(3.0),
    ),
    Benchmark(
        id="agb-final-mass-4p0",
        tier="D",
        quantity="white-dwarf mass (core at the end of the AGB) of a 4.0 Msun MIST star (solar [Fe/H]) against "
        "the white-dwarf initial-final mass relation",
        unit="Msun",
        expected=cummings2018_ifmr(4.0),
        tolerance=0.1,
        relative=False,
        source="Cummings et al. (2018) ApJ 866, 21, eqs. 4-6 (MIST-based IFMR from 80 white "
        "dwarfs in 13 clusters)",
        link="https://arxiv.org/abs/1809.01673",
        rationale="0.1 Msun: the observed scatter (0.06 Msun) plus the known offset of the "
        "models, which predict white dwarfs ~0.1 Msun lighter than observed for 3-6 Msun "
        "(Cummings et al. 2018, Sect. 5); checks the AGB extension of the MIST grids",
        compute=lambda: mist_final_mass(4.0),
    ),
    Benchmark(
        id="pulsation-time-delay-circular",
        tier="A",
        quantity="amplitude of the light-travel-time delays of a pulsating star in a circular "
        "binary (1.7 + 0.5 Msun, P = 100 d, i = 60 deg), measured from the synthetic light curve "
        "like Kepler data (30-min cadence, pulsation phases in 2-d segments)",
        unit="s",
        expected=pulsation_circular_delay_expected(),
        tolerance=1e-4,
        relative=True,
        source="Murphy et al. (2014) MNRAS 441, 2515, Sect. 3: a1 sin i = A c for a circular "
        "orbit (A = amplitude of the time delays)",
        link="https://arxiv.org/abs/1404.5649",
        rationale="analytic; the segment averaging is corrected by its sinc factor, the "
        "pulsation phase fit leaves ~1e-6",
        compute=pulsation_circular_delay,
    ),
    Benchmark(
        id="murphy2014-kic9651065-a1",
        tier="D",
        quantity="first orbital harmonic A1 of the time delays of KIC 9651065 (delta Sct in a "
        "272.7-d eccentric orbit) from its orbit (a1 sin i = 0.37 au, e = 0.468, varpi = 2.01)",
        unit="s",
        expected=165.57,
        tolerance=0.06,
        relative=True,
        source="Murphy et al. (2014) MNRAS 441, 2515, Table A1 (A1 = 165.57 +- 2.03 s)",
        link="https://arxiv.org/abs/1404.5649",
        rationale="a1 sin i is known to 0.02 au (5.4 per cent) and A1 to 1.2 per cent",
        compute=lambda: _kic9651065_harmonics()[0],
    ),
    Benchmark(
        id="murphy2014-kic9651065-a2-a1",
        tier="D",
        quantity="ratio of the second to the first orbital harmonic of the time delays of "
        "KIC 9651065 (sets the eccentricity; e = 0.468, varpi = 2.01 rad)",
        unit="",
        expected=38.41 / 165.57,
        tolerance=0.025,
        relative=False,
        source="Murphy et al. (2014) MNRAS 441, 2515, Table A1 (A2 = 38.41 +- 2.03 s) and "
        "Sect. A5 (e = 0.468 +- 0.025)",
        link="https://arxiv.org/abs/1404.5649",
        rationale="twice the uncertainty of the ratio (0.013) propagated from A1 and A2",
        compute=lambda: _kic9651065_harmonics()[1] / _kic9651065_harmonics()[0],
    ),
    Benchmark(
        id="delta-scuti-pl-1p8",
        tier="D",
        quantity="M_V of a 1.8 Msun MIST star at Teff = 7500 K on the main sequence minus M_V "
        "of the observed fundamental-mode delta Scuti P-L relation at its period from Q = 0.0329 d",
        unit="mag",
        expected=0.0,
        tolerance=0.25,
        relative=False,
        source="Poro et al. (2024), eq. 10, 2375 Galactic delta Scuti stars (Gaia DR3, TESS); "
        "Q of the fundamental mode, Breger (1979) PASP 91, 5, Table II",
        link="https://arxiv.org/abs/2401.01091",
        rationale="Q is uncertain by up to 18 per cent (Breger 1990, quoted by Poro et al.): 0.07 "
        "dex in P, 0.23 mag on the relation",
        compute=lambda: delta_scuti_pl_offset(1.8),
    ),
    Benchmark(
        id="delta-scuti-pl-2p0",
        tier="D",
        quantity="as delta-scuti-pl-1p8 for a 2.0 Msun star",
        unit="mag",
        expected=0.0,
        tolerance=0.25,
        relative=False,
        source="Poro et al. (2024), eq. 10; Breger (1979) PASP 91, 5",
        link="https://arxiv.org/abs/2401.01091",
        rationale="Q is uncertain by up to 18 per cent (Breger 1990, quoted by Poro et al.): 0.07 "
        "dex in P, 0.23 mag on the relation",
        compute=lambda: delta_scuti_pl_offset(2.0),
    ),
    Benchmark(
        id="cepheid-period-radius-6",
        tier="C",
        quantity="log R of a 6 Msun MIST Cepheid (blue loop, middle of the strip, period from "
        "Bono et al. 2000) minus the period-radius relation of rotating Geneva Cepheid models",
        unit="dex",
        expected=0.0,
        tolerance=0.05,
        relative=False,
        source="Anderson et al. (2016) A&A 591, A8, Table 5 (Z = 0.014, 2nd crossing, strip "
        "average); periods: Bono, Castellani & Marconi (2000) ApJ 529, 293, Tables 5-6",
        link="https://doi.org/10.1051/0004-6361/201528031",
        rationale="Anderson et al. note that their predicted radii differ from empirical ones by "
        "up to ~0.05 dex in log R; two independent model chains (MIST + Bono nonlinear "
        "pulsation, Geneva + linear pulsation)",
        compute=lambda: cepheid_radius_offset(6.0),
    ),
    Benchmark(
        id="cepheid-period-radius-7",
        tier="C",
        quantity="as cepheid-period-radius-6 for a 7 Msun Cepheid",
        unit="dex",
        expected=0.0,
        tolerance=0.05,
        relative=False,
        source="Anderson et al. (2016) A&A 591, A8, Table 5; Bono, Castellani & Marconi (2000) "
        "ApJ 529, 293, Tables 5-6",
        link="https://doi.org/10.1051/0004-6361/201528031",
        rationale="Anderson et al. note that their predicted radii differ from empirical ones by "
        "up to ~0.05 dex in log R",
        compute=lambda: cepheid_radius_offset(7.0),
    ),
    Benchmark(
        id="mt-soberman-general",
        tier="A",
        quantity="separation after transfer with constant fractions (donor 3 -> 1.5 Msun, "
        "accretor 2 Msun; 30 % donor wind, 40 % isotropic re-emission, 30 % accreted): "
        "numerical orbit response against the closed form",
        unit="a/a0",
        expected=soberman1997_a_ratio(3.0, 2.0, 1.5, 0.3, 0.4),
        tolerance=1e-7,
        relative=True,
        source="Soberman, Phinney & van den Heuvel (1997) A&A 327, 620, eqs. 20-28",
        link="https://arxiv.org/abs/astro-ph/9703016",
        rationale="analytic; integration tolerance 1e-11",
        compute=lambda: mt_orbit_ratio(3.0, 2.0, 1.5, 0.3, 0.4),
    ),
    Benchmark(
        id="mt-conservative",
        tier="A",
        quantity="separation after conservative transfer (donor 3 -> 1 Msun, accretor 2 -> 4)",
        unit="a/a0",
        expected=postnov_conservative(3.0, 2.0, 1.0),
        tolerance=1e-7,
        relative=True,
        source="Postnov & Yungelson (2014) Living Rev. Relativ. 17, 3, eq. 37",
        link="https://arxiv.org/abs/1403.4754",
        rationale="analytic",
        compute=lambda: mt_orbit_ratio(3.0, 2.0, 1.0, 0.0, 0.0),
    ),
    Benchmark(
        id="mt-isotropic-reemission",
        tier="A",
        quantity="separation after transfer with complete isotropic re-emission (donor "
        "3 -> 1 Msun, accretor 2 Msun)",
        unit="a/a0",
        expected=postnov_isotropic(3.0, 2.0, 1.0),
        tolerance=1e-7,
        relative=True,
        source="Postnov & Yungelson (2014) Living Rev. Relativ. 17, 3, eq. 43",
        link="https://arxiv.org/abs/1403.4754",
        rationale="analytic",
        compute=lambda: mt_orbit_ratio(3.0, 2.0, 1.0, 0.0, 1.0),
    ),
    Benchmark(
        id="mt-zeta-lobe-q2",
        tier="B",
        quantity="response zeta_L = dln R_L / dln M_d of the donor's Roche lobe (Eggleton 1983) "
        "to conservative transfer at q = 2",
        unit="",
        expected=2.13 * 2.0 - 1.67,
        tolerance=0.05,
        relative=False,
        source="Hurley, Tout & Pols (2002) MNRAS 329, 897, Sect. 2.6.1: zeta_L ~ 2.13 q - 1.67 "
        "(fit of Tout et al. 1997)",
        link="https://arxiv.org/abs/astro-ph/0201220",
        rationale="the published relation is a linear fit; it agrees with the exact zeta_L to "
        "~0.03 between q = 0.5 and 3",
        compute=lambda: _zeta_lobe(2.0),
    ),
    Benchmark(
        id="mt-bse-algol-donor",
        tier="C",
        quantity="donor mass at the end of the first mass transfer of the BSE Algol example "
        "(2.9 + 0.9 Msun, P = 2.7 d, Z = 0.02)",
        unit="Msun",
        expected=0.42,
        tolerance=0.05,
        relative=False,
        source="Hurley, Tout & Pols (2002) MNRAS 329, 897, Sect. 3.1 (M1 = 0.42 Msun when RLOF "
        "ends)",
        link="https://arxiv.org/abs/astro-ph/0201220",
        rationale="the donor ends near its helium core; the core masses of the SSE tracks in "
        "swingersim and in BSE agree to a few hundredths",
        compute=lambda: bse_algol_example()["m1"],
    ),
    Benchmark(
        id="mt-bse-algol-accretor",
        tier="C",
        quantity="accretor mass at the end of the first mass transfer of the BSE Algol example",
        unit="Msun",
        expected=2.85,
        tolerance=0.15,
        relative=False,
        source="Hurley, Tout & Pols (2002) MNRAS 329, 897, Sect. 3.1 (M2 = 2.85 Msun)",
        link="https://arxiv.org/abs/astro-ph/0201220",
        rationale="the accreted mass depends on the accretion limit of eq. 64; 0.15 Msun allows "
        "for the different stellar models",
        compute=lambda: bse_algol_example()["m2"],
        explanation="as mt-bse-algol-period: the published history needs angular momentum "
        "to be gained in the Hertzsprung gap (spin-orbit exchange in the 2002 BSE code); the "
        "BSE code in COSMIC itself does not reproduce it",
    ),
    Benchmark(
        id="mt-bse-algol-period",
        tier="C",
        quantity="orbital period at the end of the first mass transfer of the BSE Algol example",
        unit="d",
        expected=30.0,
        tolerance=0.15,
        relative=True,
        source="Hurley, Tout & Pols (2002) MNRAS 329, 897, Sect. 3.1 (P = 30 d when RLOF ends)",
        link="https://arxiv.org/abs/astro-ph/0201220",
        rationale="the period follows from the masses and the angular momentum lost with the "
        "non-accreted mass",
        compute=lambda: bse_algol_example()["P"],
        explanation="the published numbers are not self-consistent: from 2.7 d to 17 d at the "
        "end of the Hertzsprung gap the orbital angular momentum would have to grow by >= 8 % "
        "(an upper bound without any loss gives 13.5 d), while the giant phase from 17 d to "
        "30 d is exactly conservative (mt-bse-algol-giant-phase); COSMIC 4.2.1 ends at 3.2 d, "
        "swingersim at 22 d, below the 35 d upper bound",
    ),
    Benchmark(
        id="mt-cosmic-algol-donor",
        tier="C",
        quantity="donor mass when the first mass transfer of the BSE Algol example ends, against "
        "the BSE code in COSMIC 4.2.1 (2.9 + 0.9 Msun, circular 2.7 d, Z = 0.02, acc_lim = -1)",
        unit="Msun",
        expected=0.401,
        tolerance=0.03,
        relative=False,
        source="COSMIC 4.2.1 (Breivik et al. 2020, ApJ 898, 71; BSE of Hurley et al. 2002), run "
        "for this comparison: the donor leaves the Roche lobe at 0.401 Msun, 2.6 Myr after the "
        "end of its main sequence",
        link="https://cosmic-popsynth.github.io/",
        rationale="the donor is stripped to its helium core; core masses of the SSE tracks in "
        "swingersim and BSE agree to a few hundredths",
        compute=lambda: bse_algol_example()["m1"],
    ),
    Benchmark(
        id="ce-tout97-period",
        tier="B",
        quantity="orbital period after the common envelope of the CV progenitor of Tout et al. "
        "(1997): 5.36 Msun AGB star (1.014 Msun core) with a 1 Msun companion, alpha_CE = 1",
        unit="h",
        expected=35.0,
        tolerance=0.1,
        relative=True,
        source="Hurley, Tout & Pols (2002) MNRAS 329, 897, Sect. 3.3 (Tout et al. 1997: a "
        "1.014 Msun white dwarf in a 35-hr orbit)",
        link="https://arxiv.org/abs/astro-ph/0201220",
        rationale="checks the energy formalism (eqs. 69-72) with lambda = 0.5; the donor fills "
        "its Roche lobe in the orbit widened by the wind (a M = const)",
        compute=ce_tout97_period_hr,
    ),
    Benchmark(
        id="ce-bse-algol-second-period",
        tier="B",
        quantity="orbital period after the second common envelope of the BSE Algol example "
        "(2.85 Msun giant with a 0.416 Msun core, 0.417 Msun white dwarf, P = 18 d, alpha_CE = 3)",
        unit="d",
        expected=0.04,
        tolerance=0.15,
        relative=True,
        source="Hurley, Tout & Pols (2002) MNRAS 329, 897, Sect. 3.1 (orbit of 0.04 d)",
        link="https://arxiv.org/abs/astro-ph/0201220",
        rationale="the published period has one significant digit; 15 % covers its rounding",
        compute=ce_bse_algol_second_period_d,
    ),
    Benchmark(
        id="ce-impulsive-outer-a",
        tier="A",
        quantity="outer semi-major axis after a common envelope inside a triple (N-body) "
        "relative to before, divided by the impulsive-mass-loss prediction for a circular orbit",
        unit="",
        expected=1.0,
        tolerance=0.03,
        relative=True,
        source="Hadjidemetriou (1963) Icarus 2, 440; Veras et al. (2011) MNRAS 417, 2104: "
        "a'/a = (1 - f)/(1 - 2 f) for an instantaneous loss of the fraction f",
        link="https://arxiv.org/abs/1107.1239",
        rationale="the envelope (here 60 % of the inner binary) leaves within a year, much "
        "faster than the 10-yr outer orbit; 3 % allows for the initial outer eccentricity "
        "(1e-3) and the phase",
        compute=lambda: (
            ce_impulsive_triple()["a_ratio"] / ce_impulsive_triple()["a_ratio_expected"]
        ),
    ),
    Benchmark(
        id="ce-impulsive-outer-e",
        tier="A",
        quantity="outer eccentricity after the common envelope of the same triple minus the "
        "impulsive-mass-loss prediction",
        unit="",
        expected=0.0,
        tolerance=0.02,
        relative=False,
        source="Hadjidemetriou (1963): e' = f/(1 - f) for a circular orbit",
        link="https://arxiv.org/abs/1107.1239",
        rationale="as for ce-impulsive-outer-a",
        compute=lambda: ce_impulsive_triple()["e_after"] - ce_impulsive_triple()["e_expected"],
    ),
    Benchmark(
        id="matt15-sun-period",
        tier="B",
        quantity="rotation period of a solar model at 4.55 Gyr with the wind torque of Matt et "
        "al. (2015), started at P = 1 d",
        unit="d",
        expected=2 * math.pi / 2.6e-6 / 86400.0,
        tolerance=0.02,
        relative=True,
        source="Matt et al. (2015) ApJ 799, L23, eqs. 6-8 and Table 1 (the torque is normalised "
        "so that the Sun spins at Omega_sun = 2.6e-6 s^-1 at its age)",
        link="https://arxiv.org/abs/1412.4786",
        rationale="checks the torque, the Cranmer & Saar (2011) turnover time (12.9 d for the "
        "Sun) and the integration across saturation; 2 % allows for the initial spin",
        compute=matt15_sun_period_d,
    ),
    Benchmark(
        id="matt15-unsaturated-powerlaw",
        tier="A",
        quantity="Omega(4.55 Gyr) / Omega(1.14 Gyr) of the braked solar model",
        unit="",
        expected=0.5,
        tolerance=0.02,
        relative=False,
        source="Matt et al. (2015) eq. 11: Omega ~ t^(-1/p) with p = 2 (the Skumanich law)",
        link="https://arxiv.org/abs/1412.4786",
        rationale="asymptotic power law; the memory of the initial spin is below 2 % after 1 Gyr",
        compute=matt15_unsat_ratio,
    ),
    Benchmark(
        id="applegate-huaqr",
        tier="B",
        quantity="energy needed by the Applegate mechanism over the energy radiated in one "
        "modulation period, HU Aqr (thin-shell model)",
        unit="",
        expected=0.1,
        tolerance=0.1,
        relative=True,
        source="Volschow et al. (2016) A&A 587, A34, Tables 1, 2 and 4 (thin-shell model, eq. 7 "
        "after Tian et al. 2009: column Tian et al. 2009 of Table 4)",
        link="https://www.aanda.org/articles/aa/pdf/2016/03/aa27333-15.pdf",
        rationale="the published ratios have two significant digits; checks the conversion of "
        "the O - C semi-amplitude to dP/P (eq. 1) and the energy formula",
        compute=lambda: applegate_ratio("HU Aqr"),
    ),
    Benchmark(
        id="applegate-nnser",
        tier="B",
        quantity="energy needed by the Applegate mechanism over the energy radiated in one "
        "modulation period, NN Ser (thin-shell model)",
        unit="",
        expected=3.3,
        tolerance=0.1,
        relative=True,
        source="Volschow et al. (2016) A&A 587, A34, Tables 1, 2 and 4 (thin-shell model, eq. 7 "
        "after Tian et al. 2009: column Tian et al. 2009 of Table 4)",
        link="https://www.aanda.org/articles/aa/pdf/2016/03/aa27333-15.pdf",
        rationale="the published ratios have two significant digits; checks the conversion of "
        "the O - C semi-amplitude to dP/P (eq. 1) and the energy formula",
        compute=lambda: applegate_ratio("NN Ser"),
    ),
    Benchmark(
        id="applegate-v471tau",
        tier="B",
        quantity="energy needed by the Applegate mechanism over the energy radiated in one "
        "modulation period, V471 Tau (thin-shell model)",
        unit="",
        expected=0.014,
        tolerance=0.1,
        relative=True,
        source="Volschow et al. (2016) A&A 587, A34, Tables 1, 2 and 4 (thin-shell model, eq. 7 "
        "after Tian et al. 2009: column Tian et al. 2009 of Table 4)",
        link="https://www.aanda.org/articles/aa/pdf/2016/03/aa27333-15.pdf",
        rationale="the published ratios have two significant digits; checks the conversion of "
        "the O - C semi-amplitude to dP/P (eq. 1) and the energy formula",
        compute=lambda: applegate_ratio("V471 Tau"),
    ),
    Benchmark(
        id="cmdra-a-radius-excess-solar",
        tier="D",
        quantity="observed radius of CM Dra A over the MIST radius of its mass (4 Gyr, "
        "[Fe/H] = 0), minus 1",
        unit="",
        expected=0.035,
        tolerance=0.015,
        relative=False,
        source="Feiden & Chaboyer (2014) ApJ 789, 53, Sect. 3.2: standard models are 3-4 % "
        "smaller than CM Dra at solar metallicity",
        link="https://arxiv.org/abs/1405.1767",
        rationale="the radius inflation of active fully convective stars that "
        "stars.<name>.radius_inflation represents; MIST should show the same deficit as the "
        "Dartmouth models",
        compute=lambda: cmdra_radius_excess("A", 0.0),
    ),
    Benchmark(
        id="cmdra-a-radius-excess-lowz",
        tier="D",
        quantity="the same for [Fe/H] = -0.3 (the metallicity of Terrien et al. 2012)",
        unit="",
        expected=0.06,
        tolerance=0.015,
        relative=False,
        source="Feiden & Chaboyer (2014), Sect. 3.2: with the metallicity of Terrien et al. "
        "(2012) the disagreement doubles to 5-7 %",
        link="https://arxiv.org/abs/1405.1767",
        rationale="as above",
        compute=lambda: cmdra_radius_excess("A", -0.3),
    ),
    Benchmark(
        id="galactic-tide-emax-60",
        tier="A",
        quantity="largest eccentricity of a 2e4 AU orbit inclined by 60 deg to the Galactic "
        "disc (secular backend, physics.galactic_tide)",
        unit="",
        expected=math.sqrt(1 - 1.25 * 0.25),
        tolerance=0.01,
        relative=False,
        source="Grishin & Perets (2022) MNRAS 512, 4993, eq. 6 with gamma = 0 (Heisler & "
        "Tremaine 1986): e_max = sqrt(1 - (5/4) cos^2 i)",
        link="https://doi.org/10.1093/mnras/stac706",
        rationale="analytic result of the averaged vertical tide; the inner binary barely "
        "couples at this separation",
        compute=lambda: galactic_tide_emax(60.0),
    ),
    Benchmark(
        id="flyby-first-order-kick",
        tier="A",
        quantity="relative difference of the eccentricity-vector change of a binary in a "
        "parabolic passage at Q = 20 a: first-order formula against direct N-body",
        unit="",
        expected=0.0,
        tolerance=0.03,
        relative=False,
        source="Hamers & Samsing (2019) MNRAS 487, 5630, eq. 9 (first-order secular kick, "
        "valid for distant passages)",
        link="https://arxiv.org/abs/1904.09624",
        rationale="equal binary masses (no octupole term); second-order terms are of relative "
        "size ~ (a/Q)^(3/2)",
        compute=lambda: flyby_kick_error(20.0),
    ),
    Benchmark(
        id="vyn23-2p2-ma01-score",
        tier="B",
        quantity="overall score of the MA01 criterion against the N-body ghost test on swingersim's "
        "2+2 training set (sampled as in Vynatheya et al. 2023, Sect. 3.1)",
        unit="",
        expected=0.83,
        tolerance=0.05,
        relative=False,
        source="Vynatheya, Mardling & Hamers (2023) MNRAS 525, 2388, Table 3",
        link="https://arxiv.org/abs/2301.09930",
        rationale="checks that the training set resembles theirs (the score depends on the "
        "sampling, which the paper describes only in part)",
        compute=lambda: vyn23_ma01_score("2p2"),
        manual=True,
    ),
    Benchmark(
        id="vyn23-2p2-mlp-score",
        tier="B",
        quantity="held-out score of swingersim's MLP stability emulator (2+2; 11 parameters, "
        "4 x 50 neurons) against the N-body ghost test",
        unit="",
        expected=0.94,
        tolerance=0.03,
        relative=False,
        source="Vynatheya, Mardling & Hamers (2023) MNRAS 525, 2388, Table 3 ('2+2 MLP')",
        link="https://arxiv.org/abs/2301.09930",
        rationale="their network was trained on 4 x 10^5 systems, swingersim's on 8000: a "
        "slightly lower score is expected",
        compute=lambda: vyn23_mlp_score("2p2"),
        manual=True,
    ),
    Benchmark(
        id="vyn23-3p1-ma01-score",
        tier="B",
        quantity="overall score of the MA01 criterion against the N-body ghost test on swingersim's "
        "3+1 training set (sampled as in Vynatheya et al. 2023, Sect. 3.1)",
        unit="",
        expected=0.56,
        tolerance=0.05,
        relative=False,
        source="Vynatheya, Mardling & Hamers (2023) MNRAS 525, 2388, Table 4",
        link="https://arxiv.org/abs/2301.09930",
        rationale="checks that the training set resembles theirs (the score depends on the "
        "sampling, which the paper describes only in part)",
        compute=lambda: vyn23_ma01_score("3p1"),
        manual=True,
        explanation="two causes: (1) the published 0.56 is inconsistent with the precision and "
        "recall of the same table (P_s 0.54, P_u 0.95, R_s 0.95, R_u 0.55 imply a stable fraction "
        "of ~0.36 and an accuracy of ~0.69); (2) swingersim's set draws a_inner independently "
        "(0.5-30 AU) without the window near the MA01 boundary used by the paper (Sect. 3.1), so "
        "it has only 9.5 % stable systems and MA01 (calling 28.5 % stable) is right about most "
        "unstable ones; with members drawn in the window near the MA01 boundary (job 085, 10 000 "
        "members) the score is 0.51, close to the published 0.56",
    ),
    Benchmark(
        id="vyn23-3p1-mlp-score",
        tier="B",
        quantity="held-out score of swingersim's MLP stability emulator (3+1; 11 parameters, "
        "4 x 50 neurons) against the N-body ghost test",
        unit="",
        expected=0.93,
        tolerance=0.03,
        relative=False,
        source="Vynatheya, Mardling & Hamers (2023) MNRAS 525, 2388, Table 4 ('3+1 MLP')",
        link="https://arxiv.org/abs/2301.09930",
        rationale="their network was trained on 4 x 10^5 systems, swingersim's on 8000: a "
        "slightly lower score is expected",
        compute=lambda: vyn23_mlp_score("3p1"),
        manual=True,
    ),
    Benchmark(
        id="sbi-linear-gaussian-std",
        tier="A",
        quantity="SBI (neural posterior estimation): posterior width of a linear Gaussian model",
        unit="",
        expected=0.30923,
        tolerance=0.1,
        relative=True,
        source="analytic posterior of x = A theta + N(0, 0.3^2), A = [[1, 0.5], [0, 2]]",
        link="https://doi.org/10.21105/joss.02505",
        rationale="2500 training simulations; needs the optional sbi package",
        compute=sbi_linear_gaussian_std,
        manual=True,
    ),
    Benchmark(
        id="escape-theory-vs-nbody-light-escaper",
        tier="B",
        quantity="resonant binary-single encounters (1 + 0.6 Msun binary, 0.8 Msun star): "
        "N-body fraction in which the 0.6 Msun star escapes vs the statistical theory",
        unit="",
        expected=0.769,
        tolerance=0.07,
        relative=False,
        source="Stone & Leigh (2019) Nature 576, 406, eq. 2 (theory value for these encounters); "
        "they report the same kind of agreement with FEWBODY for resonant encounters",
        link="https://www.nature.com/articles/s41586-019-1833-8",
        rationale="about 200 resonant encounters: 1-sigma 0.03; swingersim N-body gives 0.754 and "
        "a mean binary eccentricity of 0.70 vs 0.72",
        compute=lambda: escape_vs_nbody()["observed"][1],
        slow=True,
    ),
    Benchmark(
        id="czev343-apsidal-tidal",
        tier="B",
        quantity="CzeV343 binary A: tidal apsidal rate [rad/day]",
        unit="rad/d",
        expected=4.14e-4,
        tolerance=0.03,
        relative=True,
        source="Pejcha et al. (2022) A&A 667, A53, Table 3 (CzeV343)",
        link="https://doi.org/10.1051/0004-6361/202244335",
        rationale="same masses, radii and log k2; the radii are given to 0.01 Rsun, which "
        "changes R^5 by up to 2 %",
        compute=lambda: czev343_apsidal_rate("tidal"),
    ),
    Benchmark(
        id="czev343-apsidal-rotational",
        tier="B",
        quantity="CzeV343 binary A: rotational apsidal rate, pseudo-synchronous spins [rad/day]",
        unit="rad/d",
        expected=6.40e-5,
        tolerance=0.05,
        relative=True,
        source="Pejcha et al. (2022) A&A 667, A53, Table 3 (CzeV343)",
        link="https://doi.org/10.1051/0004-6361/202244335",
        rationale="as above; rotation pseudo-synchronous (Hut 1981), as the 6.40e-5 of the paper "
        "implies (synchronous rotation would give 5.5e-5)",
        compute=lambda: czev343_apsidal_rate("rotational"),
    ),
    Benchmark(
        id="czev343-apsidal-gr",
        tier="B",
        quantity="CzeV343 binary A: relativistic apsidal rate [rad/day]",
        unit="rad/d",
        expected=1.47e-5,
        tolerance=0.02,
        relative=True,
        source="Pejcha et al. (2022) A&A 667, A53, Table 3 (CzeV343)",
        link="https://doi.org/10.1051/0004-6361/202244335",
        rationale="depends only on the masses, period and eccentricity",
        compute=lambda: czev343_apsidal_rate("gr"),
    ),
    Benchmark(
        id="czev343-apsidal-period",
        tier="D",
        quantity="CzeV343: observed apsidal period of binary A, full N-body quadruple [yr]",
        unit="yr",
        expected=33.4,
        tolerance=0.06,
        relative=True,
        source="Pejcha et al. (2022) A&A 667, A53: U = 33.4 +- 0.2 yr from 11 years of timings",
        link="https://doi.org/10.1051/0004-6361/202244335",
        rationale="the published theory itself explains the observed rate to 4 %; the masses "
        "are uncertain by ~10 % and the orientation of the outer orbit only partly known",
        compute=czev343_apsidal_period_nbody,
        slow=True,
    ),
    Benchmark(
        id="quint-out-zlk-emax",
        tier="A",
        quantity="2+2+1: ZLK e_max of the orbit of the two binaries driven by a fifth star, i0 = 65 deg",
        unit="",
        expected=0.81272,
        tolerance=0.01,
        relative=True,
        source="Naoz et al. (2013) MNRAS 431, 2155, quadrupole equations (A27)-(A36); "
        "Kozai (1962); Lidov (1962)",
        link="https://arxiv.org/abs/1107.2414",
        rationale="binaries with negligible secondaries act as point masses; expected value "
        "from the quadrupole theory for these masses (the test-particle limit would give 0.838); "
        "1 % as for the other ZLK benchmarks",
        compute=lambda: _quintuple_out_zlk_cached()[0],
        slow=True,
    ),
    Benchmark(
        id="vok16-pejcha-fA",
        tier="B",
        quantity="Pejcha-like 2+2 population, secular: binary A reaches 1 - e < 1e-3",
        unit="",
        expected=0.310,
        tolerance=0.07,
        relative=False,
        source="Vokrouhlicky (2016) MNRAS 461, 3964, Sect. 4 (his secular model)",
        link="https://doi.org/10.1093/mnras/stw1596",
        rationale="400 members: statistical error 0.023, tolerance 3 sigma; SWINGERsim uses the same "
        "secular equations (configs/validation/vokrouhlicky16_pejcha.yaml)",
        compute=lambda: vok16_pejcha_fraction("A"),
        slow=True,
        explanation="swingersim gives about 0.20; swingersim N-body and the independent code "
        "SecularMultiple agree with swingersim; the fraction grows logarithmically with time and "
        "reaches the published values at 4-5 times the stated length: at 50 kyr (job 080) the "
        "fractions are 0.36 (A) and 0.57 (B) and bracket the published ones",
    ),
    Benchmark(
        id="vok16-pejcha-fB",
        tier="B",
        quantity="Pejcha-like 2+2 population, secular: binary B reaches 1 - e < 1e-3",
        unit="",
        expected=0.481,
        tolerance=0.07,
        relative=False,
        source="Vokrouhlicky (2016) MNRAS 461, 3964, Sect. 4 (his secular model)",
        link="https://doi.org/10.1093/mnras/stw1596",
        rationale="as for A",
        compute=lambda: vok16_pejcha_fraction("B"),
        slow=True,
        explanation="swingersim gives about 0.38; see vok16-pejcha-fA",
    ),
    Benchmark(
        id="pejcha13-F6-fA",
        tier="B",
        quantity="Pejcha et al. F6 population, N-body: binary A reaches 1 - e < 1e-3",
        unit="",
        expected=0.358,
        tolerance=0.14,
        relative=False,
        source="Pejcha et al. (2013) MNRAS 435, 943, Table 1, entry F6",
        link="https://doi.org/10.1093/mnras/stt1281",
        rationale="100 N-body members (hours of computing): statistical error 0.048, tolerance "
        "3 sigma (configs/validation/pejcha13_F6_nbody.yaml)",
        compute=lambda: pejcha13_f6_fraction("A"),
        manual=True,
    ),
    Benchmark(
        id="pejcha13-F6-fB",
        tier="B",
        quantity="Pejcha et al. F6 population, N-body: binary B reaches 1 - e < 1e-3",
        unit="",
        expected=0.530,
        tolerance=0.15,
        relative=False,
        source="Pejcha et al. (2013) MNRAS 435, 943, Table 1, entry F6",
        link="https://doi.org/10.1093/mnras/stt1281",
        rationale="as for A; statistical error 0.050",
        compute=lambda: pejcha13_f6_fraction("B"),
        manual=True,
    ),
    Benchmark(
        id="mist-grid-tams-age-1",
        tier="A",
        quantity="MIST grid in swingersim: age of 1 Msun at the end of the main sequence (EEP 454), "
        "[Fe/H] = 0, v/vcrit = 0",
        unit="Myr",
        expected=9919.423,
        tolerance=1e-5,
        relative=True,
        source=MIST,
        link=MIST_LINK,
        rationale="checks the conversion of the MIST EEP archive and the track interpolation",
        compute=lambda: mist_tams_age(1.0),
    ),
    Benchmark(
        id="mist-interpolation",
        tier="A",
        quantity="MIST: radius of 1.14 Msun at the TAMS interpolated from 1.12 and 1.16 Msun / "
        "tabulated",
        unit="",
        expected=1.0,
        tolerance=0.02,
        relative=True,
        source=MIST,
        link=MIST_LINK,
        rationale="interpolation at fixed EEP (Dotter 2016) with twice the grid spacing; the "
        "error at the real spacing is smaller",
        compute=lambda: mist_loo(1.14, 454, 4),  # column 4: log R
    ),
    Benchmark(
        id="mist-ms-lifetime-1",
        tier="B",
        quantity="main-sequence lifetime of 1 Msun (ZAMS to TAMS), SSE vs MIST",
        unit="Myr",
        expected=9878.0,
        tolerance=0.15,
        relative=True,
        source=MIST,
        link=MIST_LINK,
        rationale=MIST_WHY,
        compute=lambda: sse_ms_lifetime(1.0),
    ),
    Benchmark(
        id="mist-ms-lifetime-2",
        tier="B",
        quantity="main-sequence lifetime of 2 Msun (ZAMS to TAMS), SSE vs MIST",
        unit="Myr",
        expected=1054.0,
        tolerance=0.15,
        relative=True,
        source=MIST,
        link=MIST_LINK,
        rationale=MIST_WHY,
        compute=lambda: sse_ms_lifetime(2.0),
    ),
    Benchmark(
        id="mist-tams-radius-2",
        tier="B",
        quantity="radius of 2 Msun at the end of the main sequence, SSE vs MIST",
        unit="Rsun",
        expected=3.913,
        tolerance=0.15,
        relative=True,
        source=MIST,
        link=MIST_LINK,
        rationale=MIST_WHY + "; the radius at the TAMS sets when a binary reaches its Roche lobe",
        compute=lambda: sse_tams_radius(2.0),
    ),
    Benchmark(
        id="claret2023-k2-2msun",
        tier="A",
        quantity="Claret (2023) grid: log k2 of 2.0 Msun, [Fe/H] = 0, log g = 4.069",
        unit="",
        expected=-2.478,
        tolerance=2e-3,
        relative=False,
        source="Claret (2023) A&A 674, A67, table of 2.0 Msun, Z = 0.0134 (CDS J/A+A/674/A67)",
        link="https://doi.org/10.1051/0004-6361/202346250",
        rationale="checks the transcription and the log g interpolation of the apsidal constants "
        "used by stellar_model.structure: claret2023",
        compute=lambda: math.log10(claret_k2_beta(2.0, 4.069)[0]),
    ),
    Benchmark(
        id="claret2023-beta-1msun",
        tier="A",
        quantity="Claret (2023) grid: radius of gyration beta of 1.0 Msun, [Fe/H] = 0, "
        "log g = 4.216",
        unit="",
        expected=0.2265,
        tolerance=5e-4,
        relative=False,
        source="Claret (2023) A&A 674, A67, table of 1.0 Msun, Z = 0.0134 (CDS J/A+A/674/A67)",
        link="https://doi.org/10.1051/0004-6361/202346250",
        rationale="checks the moment of inertia I/(M R^2) = beta^2 used for spins and tides",
        compute=lambda: claret_k2_beta(1.0, 4.216)[1],
    ),
    Benchmark(
        id="sse-zams-radius-1",
        tier="A",
        quantity="SSE track grid: ZAMS radius of 1 Msun, Z = 0.02",
        unit="Rsun",
        expected=0.88825,
        tolerance=1e-3,
        relative=True,
        source="Tout et al. (1996) MNRAS 281, 257, eq. (2) and Table 2",
        link="https://doi.org/10.1093/mnras/281.1.257",
        rationale="the SSE ZAMS is the Tout et al. fit; checks tabulation and interpolation",
        compute=lambda: sse_zams_radius(1.0),
    ),
    Benchmark(
        id="sse-zams-radius-2",
        tier="A",
        quantity="SSE track grid: ZAMS radius of 2.5 Msun, Z = 0.02",
        unit="Rsun",
        expected=1.79434,
        tolerance=1e-3,
        relative=True,
        source="Tout et al. (1996) MNRAS 281, 257, eq. (2) and Table 2",
        link="https://doi.org/10.1093/mnras/281.1.257",
        rationale="as above, a mass between grid points (expected value: eq. 2 evaluated with the Table 2 coefficients)",
        compute=lambda: sse_zams_radius(2.5),
    ),
    Benchmark(
        id="sse-sun-radius",
        tier="D",
        quantity="Sun: radius of a 1 Msun, Z = 0.02 SSE star at 4.57 Gyr",
        unit="Rsun",
        expected=1.0,
        tolerance=0.05,
        relative=True,
        source=HURLEY2000 + "; the Sun",
        link=HURLEY2000_LINK,
        rationale="the SSE formulae reproduce detailed models to about 5 per cent (Hurley et al. "
        "2000); the solar calibration of the underlying models is not exact",
        compute=lambda: sse_sun()["R"],
    ),
    Benchmark(
        id="sse-sun-luminosity",
        tier="D",
        quantity="Sun: luminosity of a 1 Msun, Z = 0.02 SSE star at 4.57 Gyr",
        unit="Lsun",
        expected=1.0,
        tolerance=0.05,
        relative=True,
        source=HURLEY2000 + "; the Sun",
        link=HURLEY2000_LINK,
        rationale="as for the radius",
        compute=lambda: sse_sun()["L"],
    ),
    Benchmark(
        id="sse-sun-convection-zone",
        tier="D",
        quantity="Sun: base of the convective envelope r_bcz/R (BSE envelope model)",
        unit="",
        expected=0.713,
        tolerance=0.03,
        relative=False,
        source="Christensen-Dalsgaard, Gough & Thompson (1991) ApJ 378, 413 (helioseismology); "
        "envelope model of Hurley, Tout & Pols (2002) MNRAS 329, 897",
        link="https://doi.org/10.1086/170441",
        rationale="the BSE envelope formulae are approximate fits; 0.03 in r/R is a 10 per cent "
        "error in the depth of the envelope, which enters the tidal dissipation and the "
        "convective turnover time",
        compute=lambda: sse_sun()["bcz"],
    ),
    Benchmark(
        id="iw-two-layer",
        tier="B",
        quantity="inertial-wave dissipation of a two-layer star (alpha 0.6, beta 0.95): "
        "impulsive formalism / closed formula",
        unit="",
        expected=1.0,
        tolerance=0.005,
        relative=True,
        source="Ogilvie (2013) MNRAS 429, 613 (closed formula; Barker 2020 MNRAS 498, 2270, "
        "eq. 39) against the numerical method of Barker (2020, eqs. 12-38; 2022, App. A)",
        link="https://doi.org/10.1093/mnras/sts362",
        rationale="two independent published methods must agree for the structure where both "
        "apply; checks the solver used for realistic profiles and its normalisation",
        compute=lambda: iw_two_layer_ratio(0.6, 0.95),
    ),
    Benchmark(
        id="iw-q-solar",
        tier="D",
        quantity="log10 Q' of inertial waves, 1 Msun main-sequence star, P_rot = 10 d",
        unit="",
        expected=7.0,
        tolerance=0.3,
        relative=False,
        source="Barker (2022) ApJL 927, L36: <Q'_IW> ~ 10^7 (P_rot/10 d)^2 for 0.2-1.2 Msun "
        "(MESA profiles)",
        link="https://doi.org/10.3847/2041-8213/ac5b63",
        rationale="the two-layer model with the SSE envelope against the full-profile result",
        compute=lambda: iw_log_q(1.0),
    ),
    Benchmark(
        id="iw-q-m-dwarf",
        tier="D",
        quantity="log10 Q' of inertial waves, fully convective 0.3 Msun star, P_rot = 10 d",
        unit="",
        expected=7.0,
        tolerance=0.3,
        relative=False,
        source="Barker (2022) ApJL 927, L36: <Q'_IW> ~ 10^7 (P_rot/10 d)^2 for 0.2-1.2 Msun",
        link="https://doi.org/10.3847/2041-8213/ac5b63",
        rationale="n = 1.5 polytrope (computed with the impulsive formalism) for fully "
        "convective stars",
        compute=lambda: iw_log_q(0.3),
    ),
    Benchmark(
        id="tides-circularization-m67",
        tier="D",
        quantity="circularization period of solar-type binaries at 4 Gyr (M67), Hurley et al. "
        "(2002) tides",
        unit="d",
        expected=12.1,
        tolerance=0.3,
        relative=True,
        source="Meibom & Mathieu (2005) ApJ 620, 970, Table 3 (M67: 12.1 +1.0/-1.5 d)",
        link="https://doi.org/10.1086/427082",
        rationale="tests the tidal-lag recipe (physics.tidal_lag: hurley2002) against the "
        "observed transition to circular orbits",
        compute=circularization_period,
        slow=True,
        explanation="equilibrium tides with the convective damping of Zahn/Hurley circularise "
        "only up to ~3.8 d in 4 Gyr on the main sequence; the pre-main sequence and inertial "
        "waves give 7.5 d (final period, tides-circularization-m67-final)",
    ),
    Benchmark(
        id="tides-circularization-m67-iw",
        tier="D",
        quantity="circularization period at 4 Gyr (M67) with inertial waves, main sequence only",
        unit="d",
        expected=12.1,
        tolerance=0.3,
        relative=True,
        source="Meibom & Mathieu (2005) ApJ 620, 970, Table 3 (12.1 +1.0/-1.5 d); mechanism: Barker (2022) ApJL 927, L36",
        link="https://doi.org/10.3847/2041-8213/ac5b63",
        rationale="physics.inertial_waves: frequency-averaged dynamical tide in the convective "
        "envelope added to the equilibrium tide",
        compute=lambda: circularization_period(4000.0, inertial_waves=True),
        explanation="7.5 d from the main sequence alone (equilibrium tide only: 3.8 d); with "
        "the pre-main sequence (MIST) 7.5 d in the final period (tides-circularization-m67-final)",
    ),
    Benchmark(
        id="tides-circularization-pms",
        tier="B",
        quantity="period below which twin 1 Msun binaries (e0 = 0.3) are circular when they "
        "reach the main sequence (MIST pre-main sequence from 1 Myr to 50 Myr; equilibrium tide "
        "and inertial waves)",
        unit="d",
        expected=7.85,
        tolerance=0.1,
        relative=True,
        source="Zahn & Bouchet (1989) A&A 223, 112, Summary: the cutoff period from the "
        "pre-main-sequence evolution lies between 7.2 and 8.5 d (0.5-1.25 Msun)",
        link="https://ui.adsabs.harvard.edu/abs/1989A%26A...223..112Z",
        rationale="10 % spans the published range; Zahn & Bouchet damp the equilibrium tide by "
        "turbulent convection, swingersim adds the inertial waves that dominate on the pre-main "
        "sequence (Barker 2020); the Hurley et al. (2002) equilibrium tide alone gives < 3 d",
        compute=lambda: circularization_period_mist(50.0),
    ),
    Benchmark(
        id="lc-occultation-small",
        tier="C",
        quantity="flux of a limb-darkened star behind a dark disk, rp = 0.1, central",
        unit="",
        expected=0.9878909547186898,
        tolerance=1e-7,
        relative=False,
        source="batman (Kreidberg 2015, PASP 127, 1161), analytic Mandel & Agol (2002) model",
        link="https://doi.org/10.1086/683602",
        rationale="both codes are exact for spheres; the tolerance is the quadrature accuracy",
        compute=lambda: lc_occultation(0.1, 0.0),
    ),
    Benchmark(
        id="lc-occultation-limb",
        tier="C",
        quantity="flux of a limb-darkened star behind a dark disk, rp = 0.3 at the limb (z = 0.9)",
        unit="",
        expected=0.9452362815022254,
        tolerance=1e-7,
        relative=False,
        source="batman (Kreidberg 2015, PASP 127, 1161)",
        link="https://doi.org/10.1086/683602",
        rationale="as above, a partial occultation crossing the limb",
        compute=lambda: lc_occultation(0.3, 0.9),
    ),
    Benchmark(
        id="lc-occultation-large",
        tier="C",
        quantity="flux of a limb-darkened star behind a larger dark disk, rp = 1.2, z = 0.5",
        unit="",
        expected=0.1289794706606836,
        tolerance=1e-7,
        relative=False,
        source="batman (Kreidberg 2015, PASP 127, 1161)",
        link="https://doi.org/10.1086/683602",
        rationale="occulter larger than the star (eclipse of a smaller companion)",
        compute=lambda: lc_occultation(1.2, 0.5),
    ),
    Benchmark(
        id="lc-detached-primary",
        tier="C",
        quantity="primary eclipse depth of a detached binary of spheres",
        unit="",
        expected=0.3808208,
        tolerance=1e-4,
        relative=False,
        source="ellc (Maxted 2016, A&A 591, A111), shape sphere, grid very_fine",
        link="https://doi.org/10.1051/0004-6361/201628579",
        rationale="ellc integrates on a finite grid; swingersim agrees with the exact batman model "
        "to 1e-8, the difference to ellc (6e-5) is ellc's grid",
        compute=lambda: lc_detached_depths()["primary"],
    ),
    Benchmark(
        id="lc-detached-secondary",
        tier="C",
        quantity="secondary eclipse depth of a detached binary of spheres",
        unit="",
        expected=0.1446576,
        tolerance=1e-4,
        relative=False,
        source="ellc (Maxted 2016, A&A 591, A111), shape sphere, grid very_fine",
        link="https://doi.org/10.1051/0004-6361/201628579",
        rationale="as above",
        compute=lambda: lc_detached_depths()["secondary"],
    ),
    Benchmark(
        id="lc-ellipsoidal",
        tier="C",
        quantity="semi-amplitude of the ellipsoidal variation, R/a = 0.15, q = 0.8, i = 60 deg",
        unit="",
        expected=0.0027949,
        tolerance=0.1,
        relative=True,
        source="ellc (Maxted 2016, A&A 591, A111) with Roche geometry",
        link="https://doi.org/10.1051/0004-6361/201628579",
        rationale="SWINGERsim uses the leading term of Morris (1985); it underestimates the "
        "variation by 3, 7, 11 and 17 per cent at R/a = 0.10, 0.15, 0.20, 0.25 (Gomel et al. "
        "2021 describe the same trend) - closer pairs need tier B (ellc)",
        compute=lc_ellipsoidal_amplitude,
    ),
    Benchmark(
        id="ml-adiabatic-invariant",
        tier="A",
        quantity="max abs(a M / (a M)_0 - 1) while m1 loses 0.3 Msun slowly (Psi = 5.6e-6), N-body",
        unit="",
        expected=0.0,
        tolerance=1e-6,
        relative=False,
        source="Jeans (1924) MNRAS 85, 2; Hadjidemetriou (1963) Icarus 2, 440: a M is an "
        "adiabatic invariant of isotropic mass loss",
        link="https://ui.adsabs.harvard.edu/abs/1963Icar....2..440H",
        rationale="the invariant holds to O(Psi) within one orbit and exactly on average; "
        "IAS15 and the mass operator keep it to ~4e-9",
        compute=ml_adiabatic_invariant,
    ),
    Benchmark(
        id="ml-osculating-e-03",
        tier="A",
        quantity="max abs(e - e_theory) / Psi over two orbits, e0 = 0.3, Psi = 5.6e-5, N-body",
        unit="",
        expected=0.0,
        tolerance=0.05,
        relative=False,
        source="Veras et al. (2011) MNRAS 417, 2104, eq. 4 integrated over the true anomaly "
        "(first order in Psi); their printed eq. 17 has 1 - e0 cos f where the integral gives "
        "1 + e0 cos f",
        link="https://arxiv.org/abs/1107.1239",
        rationale="the first-order formula neglects O(Psi^2) and the drift of the pericentre; "
        "measured residual 1.3 % of Psi",
        compute=lambda: ml_osculating_e(0.3),
    ),
    Benchmark(
        id="ml-osculating-e-07",
        tier="A",
        quantity="max abs(e - e_theory) / Psi over two orbits, e0 = 0.7, Psi = 5.6e-5, N-body",
        unit="",
        expected=0.0,
        tolerance=0.05,
        relative=False,
        source="Veras et al. (2011) MNRAS 417, 2104, eq. 4 (as above)",
        link="https://arxiv.org/abs/1107.1239",
        rationale="as above; measured residual 1.2 % of Psi",
        compute=lambda: ml_osculating_e(0.7),
    ),
    Benchmark(
        id="ml-circular-emax",
        tier="A",
        quantity="largest e / Psi of an initially circular binary under slow mass loss",
        unit="",
        expected=2.0,
        tolerance=0.02,
        relative=True,
        source="Veras et al. (2011) eq. 4 in vector form for e0 = 0: the eccentricity vector "
        "grows as Psi times the integral of the radial unit vector, |e| = 2 Psi |sin(f/2)|",
        link="https://arxiv.org/abs/1107.1239",
        rationale="first order in Psi = 5.6e-5; measured 2.003",
        compute=ml_circular_emax,
    ),
    Benchmark(
        id="ml-impulsive-a",
        tier="A",
        quantity="a'/a after a circular binary loses 15 % of its mass within 1e-3 of a period",
        unit="",
        expected=0.85 / 0.70,
        tolerance=0.01,
        relative=True,
        source="Blaauw (1961) BAN 15, 265; Hills (1983) ApJ 267, 322: instantaneous loss of a "
        "fraction d of the mass, a'/a = (1 - d)/(1 - 2 d)",
        link="https://ui.adsabs.harvard.edu/abs/1961BAN....15..265B",
        rationale="the loss takes 1e-3 of an orbit, not zero time: deviation 0.05 %",
        compute=lambda: _ml_impulsive()[0] / 1.0,
    ),
    Benchmark(
        id="ml-impulsive-e",
        tier="A",
        quantity="eccentricity after a circular binary loses 15 % of its mass within 1e-3 P",
        unit="",
        expected=0.15 / 0.85,
        tolerance=0.01,
        relative=True,
        source="Blaauw (1961); Hills (1983): e' = d / (1 - d)",
        link="https://ui.adsabs.harvard.edu/abs/1983ApJ...267..322H",
        rationale="as above: deviation 0.24 %",
        compute=lambda: _ml_impulsive()[1],
    ),
    Benchmark(
        id="ml-zlk-nbody-secular",
        tier="A",
        quantity="time of the 4th ZLK rise (e > 0.6) while the perturber loses half its mass, "
        "N-body / secular quadrupole",
        unit="",
        expected=1.0,
        tolerance=0.01,
        relative=True,
        source="ZLK quadrupole theory (Kozai 1962; Lidov 1962) with time-dependent coefficient "
        "and adiabatic orbit expansion (Hadjidemetriou 1963)",
        link="https://doi.org/10.1086/108790",
        rationale="the constant-mass N-body / secular agreement of this triple is 0.3-0.8 %; "
        "measured 0.2 %. Initial e0 = 0.2: starting at e ~ 0 (near the separatrix) the cycle "
        "timing is ill-conditioned and no secular theory reproduces it",
        compute=ml_zlk_nbody_vs_quadrupole,
        slow=True,
    ),
    Benchmark(
        id="ml-tedi-onset",
        tier="B",
        quantity="onset of the inner instability (e_in > 0.1) / crossing of the MA01 limit, "
        "inner primary losing 3 of 5 Msun in 8000 yr",
        unit="",
        expected=1.0,
        tolerance=0.15,
        relative=True,
        source="Perets & Kratter (2012) ApJ 760, 99 (TEDI: mass loss widens the inner orbit "
        "more than the outer); Mardling & Aarseth (2001) MNRAS 321, 398 stability limit",
        link="https://arxiv.org/abs/1203.2914",
        rationale="MA01 is a fit to the boundary of a chaotic region, the instability grows "
        "over many outer orbits after the crossing (cf. Toonen et al. 2022); measured 1.065",
        compute=ml_tedi_onset,
    ),
    Benchmark(
        id="mist-interpolation-massive-r",
        tier="A",
        quantity="MIST: radius of 60 Msun in the middle of the main sequence (EEP 330) "
        "interpolated from 55 and 65 Msun / tabulated",
        unit="",
        expected=1.0,
        tolerance=0.01,
        relative=True,
        source=MIST,
        link=MIST_LINK,
        rationale="interpolation at fixed EEP with twice the grid spacing; on the main sequence "
        "of massive stars 0.1-0.6 %, but from the late main sequence on the envelope inflation of "
        "40-60 Msun stars is not resolved by the grid (radius errors up to 20-40 %, see "
        "Validation)",
        compute=lambda: mist_loo(60.0, 330, 4),
    ),
    Benchmark(
        id="mist-interpolation-massive-m",
        tier="A",
        quantity="MIST: current mass of 100 Msun at the TAMS interpolated from 95 and 105 Msun / "
        "tabulated",
        unit="",
        expected=1.0,
        tolerance=0.03,
        relative=True,
        source=MIST,
        link=MIST_LINK,
        rationale="as above; the mass at the TAMS is what physics.mass_loss uses",
        compute=lambda: mist_loo(100.0, 454, 1),
    ),
    Benchmark(
        id="mist-vink2001",
        tier="A",
        quantity="MIST wind of 40 Msun at 30 % of its main sequence / Vink et al. (2001) eq. 24",
        unit="",
        expected=1.0,
        tolerance=0.1,
        relative=True,
        source="Vink, de Koter & Lamers (2001), A&A 369, 574, eq. 24 (v_inf/v_esc = 2.6, "
        "Z_sun = 0.019), the hot-star wind of MIST (Choi et al. 2016, sect. 3.6)",
        link="https://arxiv.org/abs/astro-ph/0101509",
        rationale="MIST uses the recipe with eta = 1; the 5 % (0.023 dex) offset, constant from "
        "20 to 40 Msun, is the metallicity MESA passes to the recipe (surface Z vs. 0.0142)",
        compute=mist_mdot_over_vink,
    ),
    Benchmark(
        id="claret2004-log-k2",
        tier="A",
        quantity="Claret (2004): log k2 of 39.8 Msun at log g = 3.416 (3.65 Myr)",
        unit="",
        expected=-2.815231,
        tolerance=0.002,
        relative=False,
        source="Claret (2004), A&A 424, 919, CDS J/A+A/424/919, table for log M = 1.6",
        link="https://doi.org/10.1051/0004-6361:20040470",
        rationale="resampling of the track on 400 points in log g; measured 5e-4",
        compute=lambda: claret2004_row("log_k2"),
    ),
    Benchmark(
        id="claret2004-beta",
        tier="A",
        quantity="Claret (2004): gyration radius beta of 39.8 Msun at log g = 3.416",
        unit="",
        expected=0.161706,
        tolerance=0.001,
        relative=False,
        source="Claret (2004), A&A 424, 919, CDS J/A+A/424/919, table for log M = 1.6",
        link="https://doi.org/10.1051/0004-6361:20040470",
        rationale="as above; measured 7e-6",
        compute=lambda: claret2004_row("beta"),
    ),
    Benchmark(
        id="mist-geneva-tauh-15-0",
        tier="B",
        quantity="hydrogen-burning lifetime of 15 Msun, non-rotating, MIST vs Geneva",
        unit="Myr",
        expected=11.015,
        tolerance=0.15,
        relative=True,
        source="Ekstrom et al. (2012), A&A 537, A146, Table 2 (Z = 0.014, v/v_crit = 0.0)",
        link="https://doi.org/10.1051/0004-6361/201117751",
        rationale="two model families (MESA vs Geneva): overshooting, rotational mixing and "
        "winds differ; measured +12 %",
        compute=lambda: mist_tau_h(15.0, 0.0),
    ),
    Benchmark(
        id="mist-geneva-tauh-15-4",
        tier="B",
        quantity="hydrogen-burning lifetime of 15 Msun, rotating, MIST vs Geneva",
        unit="Myr",
        expected=13.447,
        tolerance=0.15,
        relative=True,
        source="Ekstrom et al. (2012), A&A 537, A146, Table 2 (Z = 0.014, v/v_crit = 0.4)",
        link="https://doi.org/10.1051/0004-6361/201117751",
        rationale="two model families (MESA vs Geneva): overshooting, rotational mixing and "
        "winds differ; measured -6 %",
        compute=lambda: mist_tau_h(15.0, 0.4),
    ),
    Benchmark(
        id="mist-geneva-tauh-60-0",
        tier="B",
        quantity="hydrogen-burning lifetime of 60 Msun, non-rotating, MIST vs Geneva",
        unit="Myr",
        expected=3.53,
        tolerance=0.15,
        relative=True,
        source="Ekstrom et al. (2012), A&A 537, A146, Table 2 (Z = 0.014, v/v_crit = 0.0)",
        link="https://doi.org/10.1051/0004-6361/201117751",
        rationale="two model families (MESA vs Geneva): overshooting, rotational mixing and "
        "winds differ; measured +5 %",
        compute=lambda: mist_tau_h(60.0, 0.0),
    ),
    Benchmark(
        id="mist-geneva-tauh-60-4",
        tier="B",
        quantity="hydrogen-burning lifetime of 60 Msun, rotating, MIST vs Geneva",
        unit="Myr",
        expected=4.465,
        tolerance=0.15,
        relative=True,
        source="Ekstrom et al. (2012), A&A 537, A146, Table 2 (Z = 0.014, v/v_crit = 0.4)",
        link="https://doi.org/10.1051/0004-6361/201117751",
        rationale="two model families (MESA vs Geneva): overshooting, rotational mixing and "
        "winds differ; measured -12 %",
        compute=lambda: mist_tau_h(60.0, 0.4),
    ),
    Benchmark(
        id="mist-geneva-tauh-120-0",
        tier="B",
        quantity="hydrogen-burning lifetime of 120 Msun, non-rotating, MIST vs Geneva",
        unit="Myr",
        expected=2.672,
        tolerance=0.15,
        relative=True,
        source="Ekstrom et al. (2012), A&A 537, A146, Table 2 (Z = 0.014, v/v_crit = 0.0)",
        link="https://doi.org/10.1051/0004-6361/201117751",
        rationale="two model families (MESA vs Geneva): overshooting, rotational mixing and "
        "winds differ; measured +4 %",
        compute=lambda: mist_tau_h(120.0, 0.0),
    ),
    Benchmark(
        id="mist-geneva-tauh-120-4",
        tier="B",
        quantity="hydrogen-burning lifetime of 120 Msun, rotating, MIST vs Geneva",
        unit="Myr",
        expected=3.178,
        tolerance=0.15,
        relative=True,
        source="Ekstrom et al. (2012), A&A 537, A146, Table 2 (Z = 0.014, v/v_crit = 0.4)",
        link="https://doi.org/10.1051/0004-6361/201117751",
        rationale="two model families (MESA vs Geneva): overshooting, rotational mixing and "
        "winds differ; measured +5 %",
        compute=lambda: mist_tau_h(120.0, 0.4),
    ),
    Benchmark(
        id="ld-tlusty-row",
        tier="A",
        quantity="TLUSTY limb darkening: a1 of OG40000g400v10 in Bessell V",
        unit="",
        expected=1.18262,
        tolerance=1e-4,
        relative=False,
        source="Reeve & Howarth (2016), MNRAS 456, 1294, CDS J/MNRAS/456/1294 (summary2, "
        "4-coefficient law)",
        link="https://ui.adsabs.harvard.edu/abs/2016MNRAS.456.1294R",
        rationale="a grid point is reproduced exactly by the (Teff, log g) interpolation",
        compute=tlusty_row_a1,
    ),
    Benchmark(
        id="ld-tlusty-vs-atlas",
        tier="C",
        quantity="disk-averaged intensity <I>/I(1) in V at 25 000 K, log g 4: TLUSTY non-LTE / "
        "ATLAS LTE",
        unit="",
        expected=1.0,
        tolerance=0.02,
        relative=True,
        source="Reeve & Howarth (2016) vs Claret & Bloemen (2011), two model-atmosphere codes "
        "where both apply (15 000-30 000 K)",
        link="https://doi.org/10.1051/0004-6361/201116451",
        rationale="non-LTE and line blanketing change the limb darkening of B stars by a few "
        "0.1 %; the switch from ATLAS to TLUSTY at 15 000 K is therefore smooth (0.7 % there)",
        compute=tlusty_vs_atlas_mean,
    ),
]


def run_benchmarks(
    ids: list[str] | None = None, include_slow: bool = False, out: list | None = None
) -> list[Outcome]:
    """Run the benchmarks (the given ids, else the fast ones, with ``include_slow`` also the
    slow ones). A benchmark whose computation raises is recorded as ERROR (a failure) and the
    others still run. ``out`` (a list) is filled as the benchmarks finish, so that an
    interrupted run keeps what it has done. Unknown ids raise ValueError."""
    import traceback

    if ids:
        known = {b.id for b in BENCHMARKS}
        unknown = [i for i in ids if i not in known]
        if unknown:
            raise ValueError(f"unknown benchmark ids: {', '.join(unknown)}")
    out = [] if out is None else out
    for b in BENCHMARKS:
        if ids and b.id not in ids:
            continue
        if b.slow and not include_slow and not ids:
            continue
        if b.manual and not ids:
            continue
        t0 = time.perf_counter()
        try:
            value = b.compute()
        except ImportError as exc:  # an optional package (torch, phoebe, ...) is missing
            name = getattr(exc, "name", None) or str(exc)
            out.append(
                Outcome(
                    b,
                    math.nan,
                    math.nan,
                    False,
                    time.perf_counter() - t0,
                    f"needs the optional package {name}",
                )
            )
            continue
        except Exception as exc:  # noqa: BLE001 - one broken benchmark must not stop the others
            where = traceback.extract_tb(exc.__traceback__)[-1]
            out.append(
                Outcome(
                    b,
                    math.nan,
                    math.nan,
                    False,
                    time.perf_counter() - t0,
                    error=f"{type(exc).__name__}: {exc} ({Path(where.filename).name}:"
                    f"{where.lineno})",
                )
            )
            continue
        try:
            dev, ok = b.check(float(value))
        except (TypeError, ValueError) as exc:
            out.append(
                Outcome(b, math.nan, math.nan, False, time.perf_counter() - t0, error=repr(exc))
            )
            continue
        out.append(Outcome(b, value, dev, ok, time.perf_counter() - t0))
    return out


def _verdict(o: Outcome) -> str:
    if o.error:
        return "ERROR"
    if o.skipped:
        return "SKIP"
    if o.passed:
        return "PASS"
    return "EXPLAINED" if o.benchmark.explanation else "FAIL"


def format_report(outcomes: list[Outcome], command: str | None = None) -> str:
    """Markdown table of the outcomes, headed by the SWINGERsim version, the date and the
    command."""
    import platform
    from datetime import UTC, datetime

    from . import __version__

    head = (
        f"SWINGERsim {__version__} (Python {platform.python_version()}, REBOUND "
        f"{rebound.__version__}), {datetime.now(UTC).strftime('%Y-%m-%d %H:%M UTC')}"
        + (f", command: {command.strip()}" if command else "")
    )
    lines = [
        head,
        "",
        "| ID | Quantity | Expected | Obtained | Deviation | Tolerance | Seconds | Result |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for o in outcomes:
        b = o.benchmark
        dev = f"{100 * o.deviation:+.2f} %" if b.relative else f"{o.deviation:+.3g} {b.unit}"
        tol = f"{100 * b.tolerance:.1f} %" if b.relative else f"{b.tolerance:g} {b.unit}"
        lines.append(
            f"| {b.id} | {b.quantity} | {b.expected:.5g} {b.unit} | {o.obtained:.5g} {b.unit} "
            f"| {dev} | {tol} | {o.seconds:.1f} | {_verdict(o)} |"
        )
    return "\n".join(lines)


from .validation_v023 import BENCHMARKS_V023

BENCHMARKS += BENCHMARKS_V023

from .validation_discs import BENCHMARKS_DISCS
from .validation_v027 import BENCHMARKS_V027

BENCHMARKS += BENCHMARKS_DISCS
BENCHMARKS += BENCHMARKS_V027

from .validation_v037 import BENCHMARKS_V037

BENCHMARKS += BENCHMARKS_V037
