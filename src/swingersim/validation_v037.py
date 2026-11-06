"""Benchmarks of semi-detached and contact pairs in multiple systems: secular period change
(swingersim.pdot), centre-of-light radial velocities (swingersim.rvlight) and mass transfer that the
dynamics starts (physics.mass_transfer.dynamic_onset)."""

from __future__ import annotations

import math
from functools import cache

import numpy as np

from .validation import Benchmark

P_CONTACT = 0.4775  # the contact binary of the benchmark contact-vs-phoebe


def pdot_quadratic_ratio() -> float:
    """Quadratic coefficient of the primary eclipse times of a binary with dP/dt = 4e-6 d/yr
    over 60 yr (N-body with the kicks of swingersim.pdot), divided by the analytic (1/2) P dP/dt."""
    import yaml

    from . import units
    from .config import Config
    from .timing import eclipse_times

    d = yaml.safe_load(
        """
topology: custom
hierarchy: {A: [m1, m2]}
masses: {m1: 1.5, m2: 0.8}
orbits: {A: {P: 2.5, i: 89.0, Pdot: 4.0e-6}}
physics: {gr: false}
run: {t_max_yr: 60, n_output: 3}
"""
    )
    ecl, _ = eclipse_times(Config.model_validate(d), 60.0)
    t = np.array(sorted(e.t_obs_yr for e in ecl if e.kind == "primary")) * units.DAYS_PER_YEAR
    E = np.round((t - t[0]) / 2.5)
    c2 = np.polyfit(E, t, 2)[0]
    return float(c2 / (0.5 * 2.5 * 4.0e-6 / units.DAYS_PER_YEAR))


@cache
def _phoebe_rv_offsets() -> dict:
    """PHOEBE 2: flux-weighted minus dynamical radial velocities of the contact binary of
    contact-vs-phoebe (black bodies, the same power limb darkening, no irradiation) over one
    period, by component."""
    import phoebe

    from . import units
    from .photometry import limb_darkening

    q, f, incl = 0.42, 0.3, 80.0
    b = phoebe.default_binary(contact_binary=True)
    b.set_value("q", context="component", value=q)
    b.set_value("period", component="binary", context="component", value=P_CONTACT)
    b.set_value("incl", component="binary", context="component", value=incl)
    b.flip_constraint("pot", component="contact_envelope", solve_for="requiv@primary")
    b.flip_constraint("fillout_factor", component="contact_envelope", solve_for="pot")
    b.set_value("fillout_factor", component="contact_envelope", context="component", value=f)
    sma = (units.G * 1.2 * (1 + q) * (P_CONTACT / 365.25) ** 2 / (4 * math.pi**2)) ** (1 / 3)
    b.set_value("sma", component="binary", value=sma / units.RSUN_AU)
    a, _ = limb_darkening(6200.0, 4.2, 0.0, "V")
    for c in ("primary", "secondary"):
        b.set_value("teff", component=c, value=6200.0)
        b.set_value("gravb_bol", component=c, value=0.32)
    ph = np.linspace(0, 1, 81)
    b.add_dataset("rv", compute_phases=ph, passband="Johnson:V", dataset="rvfw")
    b.add_dataset("rv", compute_phases=ph, passband="Johnson:V", dataset="rvdyn")
    b.set_value_all("rv_method", dataset="rvdyn", value="dynamical")
    b.set_value_all("rv_method", dataset="rvfw", value="flux-weighted")
    b.set_value_all("atm", "blackbody")
    for k in ("ld_mode", "ld_mode_bol"):
        b.set_value_all(k, "manual")
    for k in ("ld_func", "ld_func_bol"):
        b.set_value_all(k, "power")
    for k in ("ld_coeffs", "ld_coeffs_bol"):
        b.set_value_all(k, list(a))
    b.set_value_all("irrad_method", "none")
    b.run_compute(model="qsimrv")
    out = {"phase": ph}
    for comp in ("primary", "secondary"):
        fw = np.asarray(b.get_value("rvs", dataset="rvfw", component=comp, model="qsimrv"))
        dy = np.asarray(b.get_value("rvs", dataset="rvdyn", component=comp, model="qsimrv"))
        out[comp] = (fw - dy, dy)
    return out


def rv_light_vs_phoebe() -> float:
    """Largest difference of the centre-of-light offsets (light - mass) of the secondary of the
    contact binary, swingersim minus PHOEBE, divided by the amplitude of the offset."""
    from .radvel import model_rv
    from .rvlight import light_offset
    from .validation_v027 import _contact_cfg

    ref = _phoebe_rv_offsets()
    cfg = _contact_cfg()
    t = np.linspace(0.0, P_CONTACT, 401)
    dyn = model_rv(cfg, ["m1"] * len(t), 2459000.0 + t, 2459000.0)
    ph = ref["phase"]
    # the phase of swingersim's t = 0 in PHOEBE's convention: match the dynamical velocities
    shifts = np.linspace(0.0, 1.0, 2001)
    err = [
        np.max(
            np.abs(np.interp((ph + s) % 1.0, t / P_CONTACT, dyn, period=1.0) - ref["primary"][1])
        )
        for s in shifts
    ]
    tt = ((ph + shifts[int(np.argmin(err))]) % 1.0) * P_CONTACT
    off = light_offset(cfg, ["m2"] * len(tt), tt, light_time=False)
    want = ref["secondary"][0]
    return float(np.max(np.abs(off - want)) / np.max(np.abs(want)))


ZLK_ONSET = """
topology: custom
hierarchy: {out: [A, m3], A: [m1, m2]}
masses: {m1: 1.0, m2: 0.8, m3: 1.0}
orbits:
  A: {a: 0.1, e: 0.01}
  out: {a: 1.5, e: 0.2, i: 87.0, omega: 30}
stellar_model: {kind: sse, age_myr: 1000, z: 0.02}
physics: {gr: true, evolution: true, mass_transfer: {dynamic_onset: angular_momentum}}
run: {t_max_yr: 300, n_output: 61, backend: secular}
"""


@cache
def _zlk_onset(backend: str) -> dict:
    """The first dynamically started transfer of the Kozai-Lidov triple ZLK_ONSET."""
    import yaml

    from .backends.nbody import run_nbody
    from .backends.secular import run_secular
    from .config import Config

    d = yaml.safe_load(ZLK_ONSET)
    d["run"]["backend"] = backend
    cfg = Config.model_validate(d)
    res = (run_nbody if backend == "nbody" else run_secular)(cfg)
    return res.system["_onsets"][0]


def onset_circular_fill() -> float:
    """R / R_L of the donor after the circularisation of the Kozai-Lidov triple, times
    (1 + e_c): 1 when the orbital angular momentum is kept (R = r_L a (1 - e_c) at the contact,
    a_c = a (1 - e_c^2))."""
    og = _zlk_onset("secular")
    return float(og["fill_after"] * (1.0 + og["e_onset"]))


def onset_time_nbody_over_secular() -> float:
    og_n, og_s = _zlk_onset("nbody"), _zlk_onset("secular")
    return float(og_n["time_yr"] / og_s["time_yr"])


BENCHMARKS_V037: list[Benchmark] = [
    Benchmark(
        id="pdot-parabola",
        tier="A",
        quantity="quadratic coefficient of the eclipse times of a binary with dP/dt = 4e-6 d/yr "
        "over 60 yr (N-body with period kicks) / analytic (1/2) P dP/dt",
        unit="",
        expected=1.0,
        tolerance=1e-3,
        relative=True,
        source="O - C = (1/2) (dP/dt) P E^2 (Kwee 1958; Hilditch 2001, Sect. 4.7)",
        link="https://ui.adsabs.harvard.edu/abs/1958BAN....14..131K",
        rationale="piecewise-constant period between kicks: < 0.05 s per eclipse",
        compute=pdot_quadratic_ratio,
    ),
    Benchmark(
        id="rv-light-vs-phoebe",
        tier="C",
        quantity="centre-of-light minus centre-of-mass velocity of the secondary of a contact "
        "binary (q = 0.42, fill-out 0.3, i = 80 deg): largest difference swingersim - PHOEBE / its "
        "amplitude (52 km/s)",
        unit="",
        expected=0.0,
        tolerance=0.05,
        relative=False,
        source="PHOEBE 2.5 flux-weighted radial velocities (Prsa et al. 2016; Horvat et al. "
        "2018), the same physics as contact-vs-phoebe",
        link="https://doi.org/10.3847/1538-4365/ab9ee4",
        rationale="different meshes (largest differences inside the eclipses: ~2 %)",
        compute=rv_light_vs_phoebe,
        slow=True,
    ),
    Benchmark(
        id="dynamic-onset-circularisation",
        tier="A",
        quantity="Roche-lobe filling of the donor after the circularisation of a Kozai-Lidov "
        "triple at the contact (e_c = 0.89), times (1 + e_c)",
        unit="",
        expected=1.0,
        tolerance=1e-4,
        relative=True,
        source="tidal circularisation at constant orbital angular momentum, a_c = a (1 - e^2) "
        "(Hut 1981; Fabrycky & Tremaine 2007, eq. 2)",
        link="https://ui.adsabs.harvard.edu/abs/1981A%26A....99..126H",
        rationale="the radius of the dynamics and of the new mass-transfer history come from "
        "the same tracks (agreement ~1e-9)",
        compute=onset_circular_fill,
    ),
    Benchmark(
        id="dynamic-onset-nbody-vs-secular",
        tier="C",
        quantity="time of the Roche-lobe contact driven by Kozai-Lidov cycles (first maximum, "
        "~61 yr): N-body / secular (double-averaged, octupole)",
        unit="",
        expected=1.0,
        tolerance=0.03,
        relative=True,
        source="Kozai-Lidov cycles of the secular equations (Naoz 2016) against the direct "
        "integration of the same triple",
        link="https://doi.org/10.1146/annurev-astro-081915-023315",
        rationale="the double averaging misses the short-period terms of the outer orbit "
        "(a_out / a_in = 15): ~2 % in the time",
        compute=onset_time_nbody_over_secular,
        slow=True,
    ),
]
