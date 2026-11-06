"""Configuration schema (YAML -> validated pydantic model).

Units in the file: masses [Msun], periods [day], semi-major axes [AU],
angles [deg], radii [Rsun], times [yr].
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Literal, get_args

import yaml
from pydantic import BaseModel, ConfigDict, Field, PrivateAttr, field_validator, model_validator

# Hierarchies: orbit name -> [first child, second child]; a child is a star (m1, m2, ...) or
# another orbit. The orbit of a node is the motion of the second child around the first.
PRESETS: dict[str, dict[str, list[str]]] = {
    "2+2": {"out": ["A", "B"], "A": ["m1", "m2"], "B": ["m3", "m4"]},
    "3+1": {"out": ["mid", "m4"], "mid": ["inner", "m3"], "inner": ["m1", "m2"]},
    # a fifth star on a wide orbit around a quadruple
    "2+2+1": {"out5": ["out", "m5"], "out": ["A", "B"], "A": ["m1", "m2"], "B": ["m3", "m4"]},
    "3+1+1": {
        "out5": ["out", "m5"],
        "out": ["mid", "m4"],
        "mid": ["inner", "m3"],
        "inner": ["m1", "m2"],
    },
    # binary A and a triple: binary B with a close third star m5
    "2+(2+1)": {"out": ["A", "Bt"], "A": ["m1", "m2"], "Bt": ["B", "m5"], "B": ["m3", "m4"]},
    # sextuples: a 2+2 quadruple AB with a third binary C (e.g. V994 Her), two triples
    "2+2+2": {
        "out": ["AB", "C"],
        "AB": ["A", "B"],
        "A": ["m1", "m2"],
        "B": ["m3", "m4"],
        "C": ["m5", "m6"],
    },
    "3+3": {
        "out": ["T1", "T2"],
        "T1": ["A", "m3"],
        "A": ["m1", "m2"],
        "T2": ["B", "m6"],
        "B": ["m4", "m5"],
    },
}


def check_hierarchy(spec: dict[str, list[str]]) -> tuple[str, list[str]]:
    """Check a hierarchy and return (root orbit, star names m1..mN)."""
    children = [c for pair in spec.values() for c in pair]
    for name, pair in spec.items():
        if len(pair) != 2:
            raise ValueError(f"hierarchy: orbit '{name}' needs exactly two children")
        if name.startswith("m") and name[1:].isdigit():
            raise ValueError(f"hierarchy: orbit name '{name}' looks like a star")
    stars = [c for c in children if c not in spec]
    for c in stars:
        if not (c.startswith("m") and c[1:].isdigit()):
            raise ValueError(f"hierarchy: unknown child '{c}' (stars are m1, m2, ...)")
    if len(set(children)) != len(children):
        raise ValueError("hierarchy: every star and orbit may appear only once as a child")
    roots = [n for n in spec if n not in children]
    if len(roots) != 1:
        raise ValueError(f"hierarchy: needs exactly one outermost orbit, found {roots}")
    n = len(stars)
    if sorted(stars, key=lambda x: int(x[1:])) != [f"m{k + 1}" for k in range(n)]:
        raise ValueError(f"hierarchy: stars must be m1 ... m{n}")
    if len(spec) != n - 1:
        raise ValueError("hierarchy: the orbits do not form one tree")
    return roots[0], [f"m{k + 1}" for k in range(n)]


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid")


class OrbitSpec(_Strict):
    P: float | None = Field(None, gt=0, description="period [day]")
    a: float | None = Field(None, gt=0, description="semi-major axis [AU]")
    e: float = Field(0.0, ge=0, lt=1)
    i: float = Field(0.0, ge=0, le=180, description="inclination [deg]")
    Omega: float = Field(0.0, description="longitude of ascending node [deg]")
    omega: float = Field(0.0, description="argument of pericentre [deg]")
    M0: float = Field(0.0, description="mean anomaly at t = 0 [deg]")
    Pdot: float = Field(
        0.0,
        description="secular period change dP/dt of an orbit of two stars [day/yr] (mass "
        "transfer, angular-momentum loss of semi-detached and contact binaries; P(t) = P + Pdot "
        "t from t = 0; N-body backend; swingersim.pdot); a free fit parameter as orbits.<name>.Pdot",
    )
    Pdot_mode: Literal["angular_momentum", "conservative"] = Field(
        "angular_momentum",
        description="angular_momentum: the masses stay (angular-momentum loss); conservative: "
        "conservative transfer between the two stars, dM1/dt = (dP/dt / P) M1 M2 / (3 (M1 - M2))",
    )
    contact_fillout: float | None = Field(
        None,
        ge=0,
        le=1,
        description="light curves: the two stars of this binary are a contact binary in a "
        "common envelope with fill-out factor f (Lucy 1968): their radii follow from f and "
        "they share one effective temperature (swingersim.mesh)",
    )

    @model_validator(mode="after")
    def _exactly_one_size(self) -> OrbitSpec:
        if (self.P is None) == (self.a is None):
            raise ValueError("give exactly one of P (days) or a (AU)")
        return self


class ExtraModeSpec(_Strict):
    """A further pulsation mode of the same star."""

    period_d: float | None = Field(None, gt=0, description="period [day]")
    frequency_ratio: float | None = Field(
        None, gt=0, description="frequency of this mode / frequency of the main mode"
    )
    amplitude_mag: float = Field(..., ge=0, description="semi-amplitude [mag]")
    phase: float = Field(0.0, description="phase at emission time t = 0 [cycles]")

    @model_validator(mode="after")
    def _one_period(self) -> ExtraModeSpec:
        if (self.period_d is None) == (self.frequency_ratio is None):
            raise ValueError("extra mode: give exactly one of period_d or frequency_ratio")
        return self


class PulsationSpec(_Strict):
    """Radial or non-radial pulsations of a star (light curves; swingersim engine)."""

    type: Literal[
        "manual", "lpv", "delta_scuti", "gamma_dor", "beta_cep", "spb", "rr_lyrae", "cepheid"
    ] = Field(
        "manual",
        description="lpv: period from Trabucchi et al. (2021); delta_scuti, beta_cep: from Q; "
        "cepheid: Bono et al. (2000); rr_lyrae: Bono et al. (1996); gamma_dor, spb, manual need "
        "period_d (or Q)",
    )
    mode: Literal["F", "1O", "2O", "3O"] = Field(
        "F",
        description="radial mode for delta_scuti (Q), beta_cep (F), cepheid and rr_lyrae (F, 1O)",
    )
    period_d: float | None = Field(None, gt=0, description="period [day]; overrides the relation")
    Q: float | None = Field(
        None, gt=0, description="pulsation constant [day]: P = Q / sqrt(rho / rho_sun)"
    )
    amplitude_mag: float = Field(0.01, ge=0, description="semi-amplitude [mag] in the band")
    phase: float = Field(0.0, description="phase at emission time t = 0 [cycles]")
    period_spacing_s: float | None = Field(
        None,
        gt=0,
        description="gravity modes (spb, gamma_dor): period spacing Delta P [s] of a series of "
        "n_modes modes P_0 + k Delta P",
    )
    n_modes: int = Field(1, ge=1, le=100, description="number of modes of the series")
    extra_modes: list[ExtraModeSpec] = Field(default_factory=list)


class StarSpec(_Strict):
    """Per-star values; anything left out comes from the stellar model."""

    R: float | None = Field(None, gt=0, description="radius [Rsun]")
    Teff: float | None = Field(
        None, gt=0, description="effective temperature [K] (light curves; with R it fixes L)"
    )
    apsidal_k2: float | None = Field(
        None, gt=0, lt=0.75, description="apsidal-motion constant k2 (Love number = 2 k2)"
    )
    rg2: float | None = Field(None, gt=0, le=0.4, description="moment of inertia I / (M R^2)")
    lag_s: float | None = Field(
        None, ge=0, description="constant tidal time lag [s] (Mignard 1979); 0 = no dissipation"
    )
    P_rot: float | Literal["sync", "pseudosync"] = Field(
        "sync",
        description="rotation period [day], or 'sync' = mean motion of its binary, or "
        "'pseudosync' = pseudo-synchronous rotation of an eccentric binary (Hut 1981, eq. 42; "
        "faster than the mean motion, 1.24x at e = 0.2)",
    )
    obliquity: float = Field(0.0, ge=0, le=180, description="spin-orbit angle [deg]")
    spin_azimuth: float = Field(
        0.0, description="azimuth of the spin around the orbit normal, from the node [deg]"
    )
    age_offset_myr: float = Field(
        0.0,
        description="sse, mist: effective age of this star minus stellar_model.age_myr [Myr]; "
        "an apparent age of a product of binary interaction (a stripped donor looks older, a "
        "rejuvenated gainer younger than the system)",
    )
    pulsation: PulsationSpec | None = Field(
        None, description="pulsations (light curves: flux modulated at the emission time)"
    )
    shape: Literal["sphere", "rotating", "roche"] = Field(
        "sphere",
        description="light curves: sphere; rotating: oblate Roche model of a rigidly rotating "
        "star (P_rot, obliquity), gravity-darkened; roche: equipotential surface in its "
        "circular synchronised two-star binary (both stars); see swingersim.mesh",
    )
    radius_inflation: float = Field(
        0.0,
        ge=-0.5,
        le=1.0,
        description="fractional radius change against the stellar model at fixed luminosity "
        "(e.g. 0.05: magnetically active M dwarfs, observed 3-7 % larger than models; Feiden & "
        "Chaboyer 2014); the effective temperature follows",
    )

    @model_validator(mode="after")
    def _positive_period(self) -> StarSpec:
        if self.P_rot not in ("sync", "pseudosync") and self.P_rot <= 0:
            raise ValueError("P_rot must be positive or 'sync'")
        return self


class StellarModelSpec(_Strict):
    kind: Literal["zams", "table", "python", "sse", "mist"] = Field(
        "zams",
        description="zams: zero-age main sequence (Tout et al. 1996, polytropes); sse: evolved "
        "stars of age 'age_myr' after Hurley et al. (2000); mist: MIST v1.2 tracks (Choi et "
        "al. 2016) from the pre-main sequence to the departure from the AGB; table / python: own model",
    )
    file: str | None = Field(None, description="table: CSV file (relative to the config file)")
    entry: str | None = Field(None, description="python: 'file.py:function' or 'pkg.mod:function'")
    age_myr: float | None = Field(
        None,
        ge=0,
        description="sse, mist: age of the stars at t = 0 [Myr] (coeval system; "
        "sse counts from the ZAMS, mist from the start of the pre-main sequence)",
    )
    feh: float = Field(0.0, description="mist: metallicity [Fe/H] (installed grids)")
    vvcrit: Literal[0.0, 0.4] = Field(
        0.4, description="mist: initial rotation v/v_crit of the MIST models (0 or 0.4)"
    )
    structure: Literal["auto", "polytrope", "claret2023"] | str = Field(
        "auto",
        description="apsidal constant k2 and I/(M R^2): claret2023 = MESA models of Claret "
        "(2023) at the star's mass and log g; polytrope = n 1.5 / 3; auto = claret2023 for sse, "
        "polytrope for zams; mesa:<grid> = an own MESA grid (swingersim mesa-grid / mesa-reduce, "
        "mesa/<grid>/<grid>.npz) at the star's mass, age and metallicity",
    )
    z: float = Field(
        0.02,
        gt=0,
        description="sse: metallicity Z (mass fraction; 0.02 is solar in the SSE calibration, "
        "grid 0.004-0.03)",
    )

    @model_validator(mode="after")
    def _needs_source(self) -> StellarModelSpec:
        if self.kind == "table" and not self.file:
            raise ValueError("stellar_model.kind = table needs 'file'")
        if self.kind == "python" and not self.entry:
            raise ValueError("stellar_model.kind = python needs 'entry'")
        if self.kind in ("sse", "mist") and self.age_myr is None:
            raise ValueError(f"stellar_model.kind = {self.kind} needs 'age_myr'")
        if self.structure not in ("auto", "polytrope", "claret2023") and not (
            self.structure.startswith("mesa:") and len(self.structure) > 5
        ):
            raise ValueError(
                "stellar_model.structure must be auto, polytrope, claret2023 or mesa:<grid>"
            )
        if self.kind in ("table", "python") and self.structure != "auto":
            raise ValueError("stellar_model.structure applies to kind zams, sse and mist only")
        if self.kind == "mist" and self.structure == "polytrope":
            raise ValueError(
                "stellar_model.kind = mist uses the Claret (2023) structure (Claret 2004 above 8 Msun)"
            )
        if self.kind not in ("sse", "mist") and self.age_myr is not None:
            raise ValueError("stellar_model.age_myr is used only with kind = sse or mist")
        if not 0.004 <= self.z <= 0.03:
            raise ValueError("stellar_model.z must lie in the track grid, 0.004-0.03")
        return self


class GasDamping(_Strict):
    """Damping of one orbit by the surrounding gas (time scales at full gas density)."""

    tau_a_yr: float | None = Field(
        None, description="a changes as exp(-t/tau_a): migration inwards (> 0) or outwards (< 0)"
    )
    tau_e_yr: float | None = Field(None, gt=0, description="e decays as exp(-t/tau_e)")
    tau_i_yr: float | None = Field(
        None, gt=0, description="the orbit normal turns towards the disc normal on this time scale"
    )


class GasPlane(_Strict):
    i: float = Field(0.0, ge=0, le=180)
    Omega: float = 0.0


class GasSpec(_Strict):
    """Parametric gas phase of formation (secular backend): orbits migrate, circularise and
    align with the plane of the gas while the gas is present."""

    t_end_yr: float = Field(..., ge=0, description="full gas density until this time [yr]")
    t_decay_yr: float = Field(
        ..., gt=0, description="afterwards the density falls as exp(-(t - t_end)/t_decay)"
    )
    plane: Literal["invariable"] | GasPlane = Field(
        "invariable", description="disc plane: the initial invariable plane, or {i, Omega}"
    )
    orbits: dict[str, GasDamping] = Field(
        ..., description="orbit name -> damping time scales; orbits not listed feel no gas"
    )


class MassLossSpec(_Strict):
    """Stellar mass loss in the orbital dynamics (swingersim.massloss)."""

    source: Literal["track", "rates"] = Field(
        "track",
        description="track: masses follow the stellar tracks (needs physics.evolution; MIST: "
        "current_mass of the models, SSE: the SSE wind prescription of Hurley et al. 2000); "
        "rates: constant dM/dt per star",
    )
    reimers_eta: float = Field(
        0.5,
        gt=0,
        description="SSE stars: coefficient of the Reimers (1975) giant-branch wind (Hurley et al. "
        "2000 use 0.5)",
    )
    rates: dict[str, float] = Field(
        default_factory=dict,
        description="source rates: star -> dM/dt [Msun/yr] (negative = loss); stars not listed "
        "keep their mass",
    )
    masses: Literal["initial", "current"] = Field(
        "initial",
        description="source track: the configured masses are the initial masses of the tracks "
        "(initial) or the present masses at stellar_model.age_myr (current; swingersim finds the "
        "initial masses on the tracks)",
    )
    branch: Literal["lowest", "unique"] = Field(
        "lowest",
        description="masses: current - when several tracks reach the present mass at the age "
        "(a very massive star stripped to a Wolf-Rayet star can have the mass of a less massive "
        "main-sequence star), take the lowest initial mass (the least evolved star; a warning "
        "lists the others) or refuse (unique)",
    )
    max_psi: float = Field(
        0.05,
        gt=0,
        description="secular backend: refuse runs in which an orbit's mass-loss index Psi = "
        "|dM/dt|/(n M) exceeds this (the orbit-averaged adiabatic response needs Psi << 1, "
        "Veras et al. 2011); use the N-body backend there",
    )

    @model_validator(mode="after")
    def _rates(self) -> MassLossSpec:
        if self.source == "rates" and not self.rates:
            raise ValueError("physics.mass_loss: source rates needs rates: {star: dM/dt}")
        if self.source == "track" and self.rates:
            raise ValueError("physics.mass_loss: rates are used with source: rates only")
        if self.source == "rates" and self.masses != "initial":
            raise ValueError("physics.mass_loss.masses: current needs source: track")
        return self


class GalacticTideSpec(_Strict):
    """Vertical tide of the Galactic disc on wide orbits (Heisler & Tremaine 1986; Grishin &
    Perets 2022, eqs. 3-4); secular backends."""

    rho0_msun_pc3: float = Field(
        0.185, ge=0, description="local mass density of the disc [Msun/pc^3] (solar neighbourhood)"
    )
    pole_i_deg: float = Field(
        0.0,
        ge=0,
        le=180,
        description="angle between the Galactic pole and the z axis of swingersim",
    )
    pole_Omega_deg: float = Field(
        0.0, description="azimuth of the Galactic pole (same convention as the node of an orbit)"
    )


class EncountersSpec(_Strict):
    """Distant encounters with passing stars, first-order secular kicks of every orbit
    (Hamers & Samsing 2019, eq. 9) in the secular backends; the N-body backend integrates
    every passage."""

    n_star_pc3: float = Field(0.1, gt=0, description="number density of perturbers [pc^-3]")
    sigma_kms: float = Field(
        40.0, gt=0, description="one-dimensional dispersion of the relative velocities [km/s]"
    )
    m_perturber: float = Field(0.5, gt=0, description="mass of a perturber [Msun]")
    b_max_au: float | None = Field(
        None,
        gt=0,
        description="largest impact parameter [AU]; default 30 times the outermost semi-major axis",
    )
    seed: int = Field(1, description="random seed of the encounter sequence")


class ActivitySpec(_Strict):
    """Magnetic activity cycle of a star (Applegate 1992) in the synthetic eclipse timings."""

    p_mod_yr: float = Field(..., gt=0, description="period of the activity cycle [yr]")
    amplitude_s: float = Field(
        ..., ge=0, description="semi-amplitude K of the O - C modulation of its binary [s]"
    )
    phase_deg: float = Field(0.0, description="phase of the cycle at t = 0 [deg]")


class MagneticBrakingSpec(_Strict):
    """Wind braking of stellar spins, Matt et al. (2015) (swingersim.magnetic)."""

    chi: float = Field(10.0, gt=0, description="inverse critical Rossby number of saturation")
    p: float = Field(2.0, gt=0, description="exponent of the Rossby-number dependence")
    t0_erg: float = Field(
        9.5e30, gt=0, description="torque scale T0 for the Sun [erg] (eq. 8: x R^3.1 M^0.5)"
    )
    stars: list[str] | None = Field(None, description="stars braked (default: all)")
    max_step_yr: float = Field(
        1.0e6, gt=0, description="longest operator-splitting segment while braking [yr]"
    )


class CommonEnvelopeSpec(_Strict):
    """Outcome of a common envelope with the energy formalism (swingersim.commonenv)."""

    alpha: float = Field(
        1.0,
        gt=0,
        description="efficiency alpha_CE with which orbital energy unbinds the envelope "
        "(Hurley et al. 2002, eq. 71; values above 1 stand for other energy sources)",
    )
    lam: float = Field(
        0.5,
        gt=0,
        description="structure parameter lambda of the envelope binding energy (eq. 69; BSE: 0.5)",
    )
    duration_yr: float = Field(
        1.0,
        gt=0,
        description="time over which the envelope leaves the binary; much shorter than an "
        "outer period = impulsive mass loss for the outer orbits (backend nbody)",
    )


class TertiaryTransferSpec(_Strict):
    """Mass transfer from a star onto the two-star binary it orbits (tertiary mass transfer;
    swingersim.masstransfer.tertiary_history; secular backend)."""

    accretion: float | Literal["auto"] = Field(
        "auto",
        description="fraction of the transferred mass the inner binary accretes; auto: 0.9 when "
        "the stream forms a circumbinary disc (inner separation below a_c of Leigh et al. 2020), "
        "otherwise 0.2 (Hamers et al. 2022)",
    )
    split: float | Literal["auto"] = Field(
        "auto",
        description="share of the accreted mass taken by the more massive inner star; auto: 2/3 "
        "with a disc, 1/2 without (Hamers et al. 2022)",
    )
    inner_orbit: Literal["energy", "angular_momentum"] = Field(
        "energy",
        description="energy: the ejected mass costs inner orbital energy with efficiency "
        "alpha_lambda (de Vries et al. 2014); angular_momentum: the accreted gas brings the "
        "specific angular momentum h sqrt(G M_in a_in)",
    )
    alpha_lambda: float = Field(
        5.0, gt=0, description="energy efficiency of the ejection (de Vries et al. 2014: 3-7)"
    )
    h: float = Field(
        0.0, description="angular_momentum: specific angular momentum of the accreted gas"
    )

    @model_validator(mode="after")
    def _ranges(self) -> TertiaryTransferSpec:
        for k in ("accretion", "split"):
            v = getattr(self, k)
            if v != "auto" and not 0.0 <= float(v) <= 1.0:
                raise ValueError(f"physics.mass_transfer.tertiary.{k} must lie in [0, 1]")
        return self


class ContactSpec(_Strict):
    """Contact binaries (W UMa): a binary whose two stars both fill their Roche lobes shares a
    common envelope (Lucy 1968) instead of ending the mass-transfer history."""

    braking: bool = Field(
        True,
        description="angular-momentum loss by the magnetised wind of the common envelope "
        "(Gazeas & Stepien 2008, eq. 19)",
    )
    q_min: float = Field(
        0.09,
        gt=0,
        lt=1,
        description="mass ratio below which the contact binary becomes Darwin-unstable and merges "
        "(Rasio 1995: 0.07-0.09)",
    )


class MassTransferSpec(_Strict):
    """Stable Roche-lobe overflow in two-star binaries (swingersim.masstransfer)."""

    rejuvenation: Literal["hurley", "fractional_age"] = Field(
        "hurley",
        description="hurley: a main-sequence accretor with a convective core (M > 1.25 Msun) or "
        "fully convective (M < 0.35 Msun) is rejuvenated by t' = (M / M') (t'_MS / t_MS) t "
        "(Hurley et al. 2000, Sect. 7.1); fractional_age: every star keeps its fractional age",
    )
    stability: Literal["bse", "zeta"] = Field(
        "bse",
        description="bse: critical mass ratios of Hurley et al. (2002, eq. 57); zeta: giant donors "
        "are unstable when the adiabatic exponent of a condensed polytrope (Soberman et al. "
        "1997, eq. 61) falls below that of the Roche lobe",
    )
    q_crit_ms: float = Field(
        3.0,
        gt=0,
        description="critical mass ratio M_donor / M_accretor of main-sequence (M > 0.7 Msun) and "
        "core-helium-burning donors above which transfer is dynamical (common envelope): 3.0 of "
        "BSE (Hurley et al. 2002; adiabatic q_ad 2.6-4.6, Ge et al. 2015, 2020); 1.6 of "
        "binary_c (de Mink et al. 2007: case-A contact counted as a merger)",
    )
    dynamic_onset: Literal["angular_momentum", "periastron"] | None = Field(
        None,
        description="continue when the dynamics (e.g. Kozai-Lidov cycles with tidal friction, a "
        "perturbation by another star) drives a star of a two-star binary over its Roche lobe "
        "that the mass-transfer history did not expect: the binary is circularised, conserving "
        "its orbital angular momentum, a -> a (1 - e^2) (angular_momentum; tides, Hut 1981) or "
        "at its periastron, a -> a (1 - e) (periastron: the donor keeps filling its lobe), and "
        "the run continues with the mass transfer computed from that moment. null: such an "
        "overflow ends the run (secular backends; N-body with run.stop_on_roche_overflow)",
    )
    tertiary: TertiaryTransferSpec | None = Field(
        default_factory=TertiaryTransferSpec,
        description="follow the overflow of a star onto the two-star binary it orbits (secular "
        "backend); null: the overflow ends the run",
    )
    contact: ContactSpec | None = Field(
        None,
        description="follow a contact binary (both stars fill their lobes) instead of ending "
        "the history with the event contact",
    )

    common_envelope: CommonEnvelopeSpec | None = Field(
        None,
        description="follow a common envelope with the energy formalism instead of stopping: "
        "the cores spiral in to a closer orbit or merge, the envelope leaves the system",
    )

    accretion: Literal["bse", "fractions"] = Field(
        "bse",
        description="bse: the accretor takes what its thermal time scale allows (Hurley et al. "
        "2002, eqs. 64-65), the rest leaves with the donor's specific angular momentum; "
        "fractions: fixed alpha and beta",
    )

    alpha: float = Field(
        0.0,
        ge=0,
        le=1,
        description="fraction of the transferred mass lost from the donor's "
        "vicinity (fast wind, donor's specific orbital angular momentum)",
    )
    beta: float = Field(
        0.0,
        ge=0,
        le=1,
        description="fraction lost from the accretor's vicinity (isotropic "
        "re-emission); 1 - alpha - beta is accreted (0, 0 = conservative transfer)",
    )

    @model_validator(mode="after")
    def _fractions(self) -> MassTransferSpec:
        if self.alpha + self.beta > 1.0:
            raise ValueError("physics.mass_transfer: alpha + beta must not exceed 1")
        return self


class MergerSpec(_Strict):
    """Continue a run after two stars merge (swingersim.merger): a collision of the N-body backend
    or the merger ending a common envelope or a contact binary (physics.mass_transfer)."""

    mass_loss: float | Literal["glebbeek2013"] = Field(
        "glebbeek2013",
        description="fraction of M1 + M2 ejected by the merger; glebbeek2013: 0.3 q / (1 + q)^2 "
        "with q = M2 / M1 <= 1 (head-on collisions of main-sequence stars, Glebbeek et al. 2013; "
        "2-7.5 %; Lombardi et al. 2002: at most ~8 %); a common-envelope merger keeps the mass "
        "of the energy formalism (swingersim.commonenv)",
    )
    rejuvenation: float | Literal["glebbeek", "hurley"] = Field(
        "glebbeek",
        description="fractional main-sequence age of the product = factor x the mass-weighted "
        "fractional ages of the two stars: glebbeek: 0.81 (M <= 6 Msun) to 0.64 (M >= 60 Msun), "
        "log-interpolated (Glebbeek & Pols 2008; Schneider et al. 2016); hurley: 0.1 (BSE, Hurley "
        "et al. 2002); or a number. A star past the main sequence keeps its evolutionary point "
        "on the track of the new mass",
    )
    spin_fraction: float = Field(
        0.3,
        gt=0,
        le=1,
        description="rotation of the product as a fraction of break-up (tides/rotation): born "
        "near critical, merger products spin down quickly (Schneider et al. 2019)",
    )


class Physics(_Strict):
    gr: bool = True
    gw: bool = Field(
        False,
        description="gravitational-wave emission (2.5PN) of every orbit of two single stars, "
        "orbit-averaged (Peters 1964): shrinks and circularises compact binaries",
    )
    tides: bool = False
    rotation: bool = False
    evolution: bool = False
    tidal_lag: Literal["fixed", "hurley2002"] = Field(
        "fixed",
        description="fixed: stars.<name>.lag_s (default 0); hurley2002: from the envelope "
        "(convective or radiative damping, Hurley et al. 2002)",
    )
    tidal_efficiency: float = Field(
        1.0, gt=0, description="hurley2002: factor on the dissipation rate (calibration)"
    )
    tidal_convective_reduction: Literal["quadratic", "linear"] = Field(
        "quadratic",
        description="hurley2002: reduction of the turbulent viscosity when the tidal period is "
        "shorter than the convective turnover time: quadratic (Goldreich & Nicholson 1977; "
        "Hurley et al. 2002) or linear (Zahn 1966, 1989; Zahn & Bouchet 1989)",
    )
    k2_rotation: bool = Field(
        False,
        description="correct the model apsidal-motion constants for the stars' rotation, "
        "log k2 = log k2_static - 2 V^2 / (3 g R) (Claret 1999; Claret & Gimenez 2010)",
    )
    inertial_waves: bool = Field(
        False,
        description="hurley2002: add the frequency-averaged dissipation of inertial waves in "
        "convective envelopes (Ogilvie 2013 two-layer model, as in Barker 2022)",
    )
    gas: GasSpec | None = Field(None, description="parametric gas phase (secular backend)")
    mass_loss: MassLossSpec | None = Field(
        None, description="stellar mass loss (winds) in the dynamics; see swingersim.massloss"
    )
    mass_transfer: MassTransferSpec | None = Field(
        None,
        description="stable mass transfer by Roche-lobe overflow in two-star binaries (needs "
        "physics.evolution with stellar_model sse or mist); see swingersim.masstransfer",
    )
    mergers: MergerSpec | None = Field(
        None,
        description="continue after two stars merge (collision, common-envelope or contact "
        "merger) with the merged star and the remaining orbits (swingersim.merger); null: the "
        "merger ends the run",
    )
    galactic_tide: GalacticTideSpec | None = Field(
        None,
        description="tide of the Galactic disc on the orbits (all backends); matters for "
        "orbits wider than ~1000 AU",
    )
    encounters: EncountersSpec | None = Field(
        None,
        description="encounters with passing stars (secular backends: first order in the "
        "perturbation, closer encounters are counted but not applied; N-body: every passage is "
        "integrated)",
    )
    activity: dict[str, ActivitySpec] = Field(
        default_factory=dict,
        description="star -> magnetic activity cycle (Applegate 1992): its binary's eclipses "
        "are shifted by K sin(2 pi t / P_mod + phase) in swingersim timings and minima; the energy "
        "the cycle needs is checked (thin-shell model of Tian et al. 2009)",
    )
    magnetic_braking: MagneticBrakingSpec | None = Field(
        None,
        description="spin-down of stars with convective envelopes by their magnetised winds "
        "(Matt et al. 2015; needs physics.evolution and physics.rotation); see swingersim.magnetic",
    )

    @model_validator(mode="after")
    def _implemented(self) -> Physics:
        if self.tidal_lag != "fixed" and not self.tides:
            raise ValueError("physics.tidal_lag needs physics.tides: true")
        if self.inertial_waves and self.tidal_lag != "hurley2002":
            raise ValueError("physics.inertial_waves needs physics.tidal_lag: hurley2002")
        if self.magnetic_braking is not None and not (self.evolution and self.rotation):
            raise ValueError(
                "physics.magnetic_braking needs physics.evolution and physics.rotation (with "
                "physics.tides: the spins and their coupling to the orbits)"
            )
        if self.rotation and not self.tides:
            raise ValueError(
                "physics.rotation needs physics.tides (both use the same stellar distortion model)"
            )
        return self


class RunSpec(_Strict):
    backend: Literal["nbody", "secular"] = Field(
        "nbody", description="nbody (full) or secular (orbit-averaged)"
    )
    integrator: Literal["ias15"] = "ias15"
    t_max_yr: float = Field(
        ..., description="end time [yr]; negative = integrate into the past (backward)"
    )
    n_output: int = Field(1000, ge=2)
    collisions: bool = True
    stop_on_hierarchy_change: bool = Field(
        False, description="stop at the first confirmed change of the hierarchy"
    )
    extremes_every_orbits: float = Field(
        10.0,
        ge=0,
        description="internal sampling for e_max and closest periastron, in periods of the "
        "shortest inner orbit; 0 = only at outputs",
    )
    rate_window_yr: float | None = Field(
        None,
        gt=0,
        description="length of the linear fit for the present rates of change of the "
        "inclinations (*_didt0_deg_yr, *_dtiltdt0_deg_yr), on each side of t = 0; e.g. half "
        "of the observing baseline. Default: two periods of the perturbing orbit, at least "
        "6 outputs",
    )
    megno: bool = Field(False, description="compute the MEGNO chaos indicator (Newtonian)")
    secular_octupole: bool = Field(True, description="secular backend: include octupole terms")
    secular_rtol: float = Field(
        1e-10, gt=0, lt=1e-3, description="secular backend: relative tolerance of the integrator"
    )
    stop_on_roche_overflow: bool = Field(
        False,
        description="stop when a star fills its Roche lobe (in a two-star binary or orbiting a "
        "subsystem)",
    )
    hierarchy_persistence: int = Field(
        3, ge=1, description="outputs a new hierarchy must persist before it is confirmed"
    )
    stop_when_e_reached: dict[str, float] = Field(
        default_factory=dict,
        description="N-body: stop once every listed orbit has reached its eccentricity, e.g. "
        "{A: 0.999, B: 0.999} (event e_reached)",
    )
    max_wall_s: float | None = Field(
        None,
        gt=0,
        description="stop the run after this many seconds of computing (event wall_time)",
    )
    evolution_max_dr: float = Field(
        0.01,
        gt=0,
        le=0.2,
        description="physics.evolution: largest relative change of a stellar radius within one "
        "operator-splitting segment",
    )
    max_energy_error: float = Field(
        1e-5,
        gt=0,
        description="abort if |dE/E| exceeds this; with GR the 1PN Hamiltonian is checked (~1e-11)",
    )


class Distribution(_Strict):
    """A random distribution of one ensemble parameter; give exactly one kind."""

    uniform: tuple[float, float] | None = None
    loguniform: tuple[float, float] | None = None
    normal: tuple[float, float] | None = Field(None, description="(mean, sigma)")
    isotropic: tuple[float, float] | None = Field(
        None, description="inclination [deg] with cos(i) uniform between the two limits"
    )
    choice: list[float] | None = None
    min: float | None = Field(None, description="truncation (the value is drawn again)")
    max: float | None = None
    mirror: bool = Field(
        False,
        description="inclinations: replace i by 180 - i with probability 1/2 (photometry does "
        "not tell the sense of revolution of an eclipsing binary)",
    )

    @model_validator(mode="after")
    def _one_kind(self) -> Distribution:
        kinds = [
            k
            for k in ("uniform", "loguniform", "normal", "isotropic", "choice")
            if getattr(self, k) is not None
        ]
        if len(kinds) != 1:
            raise ValueError("give exactly one of uniform, loguniform, normal, isotropic, choice")
        if self.loguniform is not None and min(self.loguniform) <= 0:
            raise ValueError("loguniform limits must be positive")
        if self.normal is not None and self.normal[1] < 0:
            raise ValueError("normal: sigma must not be negative")
        if self.choice is not None and not self.choice:
            raise ValueError("choice needs at least one value")
        return self


class ObservedValue(_Strict):
    """A measured quantity or a hard limit on a column of the ensemble summary.

    ``value`` and ``sigma`` give a Gaussian likelihood; ``min`` and ``max`` give hard limits
    (weight 0 outside), e.g. that a star of the observed system has not filled its Roche lobe
    in the past. The column name may start with ``forward:``, ``backward:`` or ``central:``
    (mean of both) to take the value from that run. Default: central for the present rates of
    change (``*_didt0_deg_yr``, ``*_dtiltdt0_deg_yr``), forward for everything else."""

    value: float | None = None
    sigma: float | None = Field(None, gt=0)
    min: float | None = None
    max: float | None = None
    missing_ok: bool = Field(
        False,
        description="limits only: a missing value (NaN, e.g. no change during the run) passes",
    )
    equals: str | list[str] | None = Field(
        None,
        description="text columns (e.g. fate): the member passes only with this value (or one of "
        "these values); e.g. 'backward:fate: {equals: intact}' requires that the system "
        "survived the run into the past",
    )

    @model_validator(mode="after")
    def _kind(self) -> ObservedValue:
        gauss = self.value is not None or self.sigma is not None
        if self.equals is not None:
            if gauss or self.min is not None or self.max is not None:
                raise ValueError("equals cannot be combined with value, sigma, min or max")
            return self
        if gauss and (self.value is None or self.sigma is None):
            raise ValueError("give both value and sigma")
        if not gauss and self.min is None and self.max is None:
            raise ValueError("give value and sigma, or min and/or max (or equals)")
        if self.missing_ok and gauss:
            raise ValueError("missing_ok applies to limits (min / max) only")
        return self


VARY_PREFIXES = {"masses": None, "orbits": OrbitSpec, "stars": StarSpec}


def _is_numeric(annotation) -> bool:
    """Can a field hold a drawn number (float or int, possibly among other choices)?"""
    if annotation in (float, int):
        return True
    return any(a in (float, int) for a in get_args(annotation))


class EphemerisSpec(_Strict):
    T0: float = Field(
        ..., description="reference time of a primary minimum [d], in the time system of the table"
    )
    P: float = Field(..., gt=0, description="period of the linear ephemeris [d]")


class FitExcludeSpec(_Strict):
    before: float | None = Field(None, description="drop minima before this time [d]")
    after: float | None = Field(None, description="drop minima after this time [d]")
    sources: list[str] = Field(default_factory=list, description="drop minima of these sources")
    pairs: list[str] = Field(default_factory=list, description="drop these pairs")


class CoordinatesSpec(_Strict):
    ra_deg: float = Field(..., description="right ascension [deg] (ICRS)")
    dec_deg: float = Field(..., description="declination [deg] (ICRS)")


class RVSpec(_Strict):
    """Measured radial velocities fitted together with the minima in stage 2 (radvel.py)."""

    file: str = Field(
        ...,
        description="CSV table: star (m1, m2, ...), time (BJD_TDB), rv [km/s], sigma [km/s], "
        "source; path relative to the configuration",
    )
    gamma: float = Field(
        0.0, description="systemic velocity [km/s] (start value; free as rv.gamma)"
    )
    offsets: dict[str, float] = Field(
        default_factory=dict,
        description="zero-point offset [km/s] of a source added to its velocities (start "
        "values; free as rv.offset.<source>)",
    )
    stars: list[str] = Field(default_factory=list, description="stars used (default: all)")
    max_sigma: float | None = Field(
        None, gt=0, description="leave out velocities with a larger error [km/s]"
    )
    centre: Literal["mass", "light"] = Field(
        "mass",
        description="mass: the velocity of the star's centre of mass; light: the flux-weighted "
        "velocity of its visible surface (swingersim.rvlight: tidal distortion, synchronous "
        "rotation, the Rossiter-McLaughlin effect in eclipses) - for close, semi-detached and "
        "contact binaries",
    )
    band: str = Field("V", description="centre light: passband of the spectra (weights)")


class InclinationSpec(_Strict):
    """Measured orbital inclinations of an eclipsing pair against time (e.g. from light-curve
    solutions of separate epochs or eclipse depths), fitted together with the minima in stage 2
    (incfit.py). Only |cos i| is compared (the sense of the orbit is not measured)."""

    file: str = Field(
        ...,
        description="CSV table: mid_hjd or time (BJD/HJD), i_deg, sigma_deg[, source]; path "
        "relative to the configuration",
    )
    pair: str = Field(..., description="orbit name of the eclipsing pair (e.g. inner)")
    label: str | None = Field(None, description="name in the report (default: the file name)")
    min_sigma_deg: float = Field(
        0.0, ge=0, description="lower limit of the errors [deg] (model systematics)"
    )


class LCSpec(_Strict):
    """A measured light curve fitted together with the minima in stage 2 (lcfit.py)."""

    file: str = Field(
        ...,
        description="CSV table: time (BJD_TDB), flux (normalised), sigma[, segment]; path "
        "relative to the configuration",
    )
    label: str | None = Field(None, description="name in the report (default: the file name)")
    band: str = Field(
        "TESS", description="passband: U, B, V, R, I, TESS, Gaia_G, Gaia_BP or Gaia_RP"
    )
    exposure_min: float = Field(
        0.0, ge=0, description="exposure or binning time [min] (smearing of the model)"
    )
    n_sub: int = Field(1, ge=1, le=15, description="model points per exposure")
    baseline_degree: int = Field(
        1, ge=0, le=2, description="polynomial trend per segment solved with the scale (0-2)"
    )
    scale: bool = Field(
        True,
        description="scale the model per segment (dilution by extra light); false: the "
        "depths follow from the stars as given",
    )
    gap_d: float = Field(
        0.5, gt=0, description="without a segment column: split segments at gaps > gap_d [d]"
    )
    segments: list[str] = Field(default_factory=list, description="use only these segments")
    shapes: Literal["auto", "as_configured", "spheres"] = Field(
        "auto",
        description="star shapes of the light-curve model: auto = Roche surfaces (mesh) for "
        "the stars of circular two-star binaries that fill their Roche lobe beyond "
        "roche_threshold (and contact binaries as configured), spheres otherwise; "
        "as_configured = stars.<m>.shape and contact_fillout of the configuration; spheres = "
        "spheres everywhere (fastest)",
    )
    roche_threshold: float = Field(
        0.5,
        gt=0,
        le=1.0,
        description="auto: volume radius / Roche-lobe radius (Eggleton 1983) above which a "
        "star gets its Roche surface (spheres err by ~2 ppt at 0.5, ~5 ppt at 0.7, ~17 ppt "
        "at 0.95)",
    )
    mesh_resolution: tuple[int, int] = Field(
        (64, 128), description="mesh directions (theta, phi) of the accurate model"
    )
    fidelity: Literal["mapped", "direct"] = Field(
        "mapped",
        description="with meshes: mapped = the fit runs on a coarse mesh plus the difference "
        "accurate - coarse at the current solution, refreshed until it settles (output "
        "space mapping; converges to the accurate solution at a fraction of the cost); "
        "direct = the accurate mesh in every evaluation",
    )
    coarse_resolution: tuple[int, int] = Field(
        (24, 48), description="mapped: mesh directions (theta, phi) of the coarse model"
    )
    red_noise: bool = Field(
        False,
        description="inflate the errors of each segment by the red-noise factor beta of the "
        "time-averaging method (Pont et al. 2006, Winn et al. 2008): how much more the "
        "residuals averaged over 20-120 min scatter than white noise would; applied in the "
        "fits in a row (needs fit.stage2.error_scale or clip_sigma), so that correlated "
        "residuals (instrumental trends, oscillations, spots) do not give too small errors",
    )


class Stage1Spec(_Strict):
    terms: list[Literal["lite", "apsidal", "quadratic"]] = Field(
        default_factory=lambda: ["lite"],
        description="terms of the analytic O - C model of every pair (as OCFit): lite = light-time "
        "effect of an outer orbit (Irwin 1952), apsidal = apsidal motion of an eccentric pair, "
        "quadratic = quadratic ephemeris term",
    )
    pairs: dict[str, list[Literal["lite", "apsidal", "quadratic"]]] = Field(
        default_factory=dict,
        description="terms per pair, overriding terms (e.g. {A: [lite, apsidal], B: [lite]})",
    )
    shared_outer: bool = Field(
        False,
        description="fit the pairs together with one outer orbit: common P3, T3, e3 and the "
        "pericentre of B opposite to A (anti-phase); separate a sin i of each pair",
    )
    P3_range_d: tuple[float, float] | None = Field(
        None,
        description="search range of the outer period [d] (default: 2 x the shortest "
        "spacing of the data to 3 x the time span)",
    )
    clip_sigma: float | None = Field(
        None,
        gt=0,
        description="leave out minima that deviate by more than this many robust standard "
        "deviations of their source and fit again (e.g. 5)",
    )
    max_clip_fraction: float = Field(
        0.1, gt=0, le=0.5, description="at most this fraction of the minima is left out"
    )
    sampler: Literal["none", "emcee"] = Field(
        "none", description="none: least squares with covariance errors; emcee: MCMC afterwards"
    )
    walkers: int = Field(48, ge=8, description="emcee walkers")
    steps: int = Field(4000, ge=100, description="emcee steps per walker")
    seed: int = 1


class Stage2Spec(_Strict):
    """Joint dynamical model of all eclipsing pairs (swingersim fit, stage 2): the N-body times of
    minima of the whole system are fitted to the measured minima (dynfit.py)."""

    bjd_t0: float = Field(
        ..., description="BJD_TDB of t = 0 of the configuration (the epoch of its elements)"
    )
    free: list[str] = Field(
        ...,
        min_length=1,
        description="free parameters: orbits.<orbit>.<P|M0|e|omega|Omega|i|ecosw|esinw|Pdot>, "
        "masses.<star>, stars.<star>.<R|Teff|apsidal_k2|rg2>, rv.gamma and rv.offset.<source>; "
        "everything else stays as in the configuration. ecosw / esinw replace e and omega of an "
        "orbit with a small eccentricity",
    )
    pairs: list[str] = Field(
        default_factory=list,
        description="pairs whose minima are fitted (default: all in the table)",
    )
    max_sigma_d: float | None = Field(
        None, gt=0, description="leave out minima with a larger error [d] (e.g. visual timings)"
    )
    anchor: bool = Field(
        False,
        description="first anchor the period P and the mean anomaly M0 of every eclipsing pair "
        "to the minima (as swingersim fit --anchor) and start the fit from there; for published "
        "photodynamical elements, whose eclipses drift away in the N-body model",
    )
    priors: dict[str, tuple[float, float]] = Field(
        default_factory=dict,
        description="Gaussian priors [mean, sigma] (or {normal: [mean, sigma]} as in "
        "ensemble.vary) on free parameters or on e / omega of an orbit fitted with ecosw / "
        "esinw (and the reverse), e.g. from published light-curve and radial-velocity "
        "solutions: {orbits.AB.e: [0.687, 0.044], masses.m1: [2.929, 0.11]}; each adds "
        "((value - mean) / sigma)^2 to chi^2",
    )
    fd_steps: dict[str, float] = Field(
        default_factory=dict,
        description="finite-difference steps of the Jacobian (default: P 1e-6 relative, angles "
        "0.01 deg, e 1e-4, masses 1e-4 relative)",
    )
    bounds: dict[str, tuple[float, float]] = Field(
        default_factory=dict,
        description="bounds of free parameters (default: P within 0.5-2x the start, e 0-0.95, "
        "i 0-180 deg, masses 0.05-20x)",
    )
    loss: Literal["linear", "soft_l1", "huber", "cauchy"] = Field(
        "linear",
        description="linear: ordinary least squares; soft_l1 / huber / cauchy: robust against "
        "outlying minima beyond f_scale sigma",
    )
    f_scale: float = Field(3.0, gt=0, description="residual [sigma] where a robust loss starts")
    error_scale: Literal["none", "global", "per_source"] | None = Field(
        None,
        description="errors of the minima after a first fit: per_source = every source (column "
        "source) gets its own factor so that its reduced chi^2 is 1 (never below its quoted "
        "errors), global = one factor for all, none = quoted errors; the fit is then repeated "
        "(default: fit.error_scale)",
    )
    error_model: Literal["scale", "jitter"] = Field(
        "scale",
        description="scale: sigma x factor; jitter: sigma^2 + jitter^2 (an error floor of the "
        "source, e.g. for timings with unrealistically small quoted errors)",
    )
    min_source_minima: int = Field(
        5,
        ge=1,
        description="sources with fewer minima share one factor ('other sources')",
    )
    clip_sigma: float | None = Field(
        None,
        gt=0,
        description="after a fit, leave out minima that deviate by more than this many robust "
        "standard deviations of their source, then repeat the fit (e.g. 5)",
    )
    max_clip_fraction: float = Field(
        0.1, gt=0, le=0.5, description="at most this fraction of the minima is left out"
    )
    refits: int = Field(
        4,
        ge=1,
        le=10,
        description="at most this many fits in a row when error_scale or clip_sigma is used "
        "(fit, rescale / clip, fit again, ... until nothing changes)",
    )
    max_iterations: int = Field(
        30, ge=1, description="iterations of the optimiser (one Jacobian each)"
    )
    processes: int = Field(
        0,
        ge=0,
        description="cores for the Jacobian (0 = all; swingersim fit --processes overrides)",
    )
    sampler: Literal["none", "emcee"] = Field(
        "none",
        description="none: linear errors of the covariance matrix; emcee: MCMC around the "
        "solution afterwards (walkers x steps model evaluations, in parallel) for errors that "
        "hold in a curved valley of chi^2 (strongly correlated parameters)",
    )
    walkers: int = Field(
        0, ge=0, description="emcee walkers (0: 2 x free parameters + 2; at least that many, even)"
    )
    mcmc_steps: int = Field(
        300, ge=20, description="emcee steps per walker; the first half is discarded"
    )
    mcmc_tau_target: float = Field(
        0.0,
        ge=0,
        description="0: exactly mcmc_steps; e.g. 50: after mcmc_steps the chain is extended (in "
        "steps of a quarter of its length) until it is at least this many integrated "
        "autocorrelation times tau long (the emcee documentation, 'Autocorrelation analysis & "
        "convergence', advises at least 50 tau) "
        "or reaches mcmc_max_steps; the first half is discarded",
    )
    mcmc_max_steps: int = Field(
        5000, ge=20, description="mcmc_tau_target: the longest chain [steps per walker]"
    )
    seed: int = 1

    @field_validator("priors", mode="before")
    @classmethod
    def _priors_form(cls, v):
        if not isinstance(v, dict):
            return v
        out = {}
        for k, p in v.items():
            if isinstance(p, dict):
                if set(p) != {"normal"}:
                    raise ValueError(
                        f"fit.stage2.priors.{k}: [mean, sigma] or {{normal: [mean, sigma]}}"
                    )
                p = p["normal"]
            out[k] = p
        return out

    @model_validator(mode="after")
    def _check_priors(self) -> Stage2Spec:
        free = set(self.free)
        for what in ("bounds", "fd_steps"):
            bad = sorted(set(getattr(self, what)) - free)
            if bad:
                raise ValueError(
                    f"fit.stage2.{what}: {bad} are not free parameters (fit.stage2.free)"
                )
        for k, (lo, hi) in self.bounds.items():
            if not lo < hi:
                raise ValueError(f"fit.stage2.bounds.{k}: the lower bound must be below the upper")
        if self.sampler == "none":
            given = sorted(
                k
                for k in ("walkers", "mcmc_steps", "mcmc_tau_target", "mcmc_max_steps")
                if getattr(self, k) != type(self).model_fields[k].default
            )
            if given:
                warnings.warn(
                    f"fit.stage2: {', '.join(given)} have no effect with sampler: none",
                    UserWarning,
                    stacklevel=2,
                )
        for k, (_, s) in self.priors.items():
            if not s > 0:
                raise ValueError(f"fit.stage2.priors.{k}: sigma must be > 0")
            parts = k.split(".")
            ok = k in free
            if len(parts) == 3 and parts[0] == "orbits":
                o, el = parts[1], parts[2]
                if el in ("e", "omega"):
                    ok |= {f"orbits.{o}.ecosw", f"orbits.{o}.esinw"} <= free
                if el in ("ecosw", "esinw"):
                    ok |= {f"orbits.{o}.e", f"orbits.{o}.omega"} <= free
            if not ok:
                raise ValueError(
                    f"fit.stage2.priors.{k}: not a free parameter (a prior on a fixed one has "
                    "no effect)"
                )
        return self


class PredictSpec(_Strict):
    """Predicted minima with their uncertainty (fit.output.predict)."""

    years: float = Field(5.0, gt=0, le=200, description="how far ahead [yr]")
    pairs: list[str] = Field(default_factory=list, description="pairs (default: all fitted)")
    start: float | None = Field(
        None,
        description="BJD_TDB from which to predict (default: now, or the last measured "
        "minimum if it is later)",
    )
    samples: int = Field(
        100,
        ge=10,
        le=2000,
        description="stage 2: parameter sets drawn from the emcee samples or the covariance "
        "matrix, one N-body integration each (in parallel), for the uncertainty of every minimum",
    )


class FitOutputSpec(_Strict):
    dir: str = Field("output/fit", description="output folder (relative to the configuration)")
    tables: list[Literal["csv", "latex"]] = Field(default_factory=lambda: ["csv"])
    plots: list[Literal["oc", "oc_components", "corner"]] = Field(
        default_factory=lambda: ["oc", "oc_components"]
    )
    predict: PredictSpec | None = Field(
        None,
        description="{years: 5, pairs: [A, B]}: predicted minima with their uncertainty "
        "(stage 1: stage1_predictions.csv; stage 2: stage2_predictions.csv and "
        "stage2_predictions.png)",
    )


class FitSpec(_Strict):
    """Fitting of observed times of minima (swingersim fit). Stage 1 is the analytic O - C model
    of every pair, equivalent to OCFit; stage 2 is the joint dynamical (N-body) model of all
    pairs."""

    minima: str = Field(
        ...,
        description="CSV table of minima: pair, type (p/s), time, sigma, "
        "source, method[, depth, width]; path relative to the configuration",
    )
    time_system: Literal["BJD_TDB", "HJD_UTC", "JD_UTC"] = Field(
        "BJD_TDB",
        description="time system of the table; HJD_UTC and JD_UTC are converted to "
        "BJD_TDB (needs coordinates)",
    )
    coordinates: CoordinatesSpec | None = None
    ephemerides: dict[str, EphemerisSpec] = Field(
        default_factory=dict,
        description="linear ephemeris of every pair for the epochs and the O - C (default: "
        "fitted to the primary minima, period from orbits.<pair>.P)",
    )
    exclude: FitExcludeSpec = Field(default_factory=FitExcludeSpec)
    lc: list[LCSpec] = Field(
        default_factory=list,
        description="light curves (segments) fitted together with the minima (stage 2), e.g. "
        "the windows around eclipses of a third body",
    )
    inclination: list[InclinationSpec] = Field(
        default_factory=list,
        description="orbital inclinations of a pair against time fitted together with the "
        "minima (stage 2), e.g. [{file: data/inclination.csv, pair: inner}]",
    )
    rv: RVSpec | None = Field(
        None,
        description="radial velocities fitted together with the minima (stage 2); free "
        "parameters rv.gamma and rv.offset.<source>",
    )
    error_scale: Literal["none", "global", "per_source"] = Field(
        "none",
        description="rescale the errors so that the reduced chi^2 is 1: of all minima "
        "(global) or of each source separately (per_source)",
    )
    stage1: Stage1Spec = Field(default_factory=Stage1Spec)
    stage2: Stage2Spec | None = Field(
        None,
        description="joint dynamical model of all pairs: N-body times of minima fitted to the "
        "measured minima (swingersim fit --stage 2)",
    )
    output: FitOutputSpec = Field(default_factory=FitOutputSpec)

    @model_validator(mode="after")
    def _check(self) -> FitSpec:
        if self.time_system != "BJD_TDB" and self.coordinates is None:
            raise ValueError(f"fit.time_system {self.time_system} needs fit.coordinates")
        st = self.stage2
        if st is not None and any(lc.red_noise for lc in self.lc):
            scale = st.error_scale if st.error_scale is not None else self.error_scale
            if scale == "none" and st.clip_sigma is None:
                raise ValueError(
                    "fit.lc red_noise: true is applied between the fits in a row and needs "
                    "fit.stage2.error_scale (or fit.error_scale) global / per_source or "
                    "fit.stage2.clip_sigma; otherwise it would be ignored"
                )
        return self


class EnsembleSpec(_Strict):
    """Monte Carlo ensemble: the configuration is the nominal system, ``vary`` draws the
    uncertain or unknown parameters for every member."""

    n: int = Field(..., ge=1, description="number of members (drawn systems)")
    seed: int = 1
    vary: dict[str, float | Literal["isotropic", "angle"] | Distribution] = Field(
        default_factory=dict,
        description="'orbits.out.e', 'masses.m1', 'stars.m2.R', ... -> value or distribution",
    )
    directions: list[Literal["forward", "backward"]] = Field(
        default_factory=lambda: ["forward"],
        description="forward: to +|t_max_yr|, backward: to -|t_max_yr| (each member both ways "
        "if both are given)",
    )
    require: list[Literal["stable_ma01", "stable_v22", "hierarchical"]] = Field(
        default_factory=list,
        description="redraw members that the criterion classifies as unstable; hierarchical: "
        "every inner orbit's apocentre inside the pericentre of the orbit around it",
    )
    ma01_unstable_keep: float = Field(
        1.0,
        ge=0,
        le=1,
        description="keep only this fraction of the draws that the MA01 criterion classifies "
        "as unstable (balanced training sets for stability classifiers; Vynatheya et al. 2023 "
        "kept 0.3 for 2+2 and 0.2 for 3+1)",
    )
    ma01_window: tuple[float, float] | None = Field(
        None,
        description="keep only draws whose smallest MA01 ratio (R_p / R_p,crit over the nested "
        "triples) lies in [low, high]: training sets close to the stability boundary "
        "(Vynatheya et al. 2023, Sect. 3.1)",
    )
    max_tries: int = Field(1000, ge=1, description="draws allowed per member")
    mlp: bool = Field(False, description="add the Vynatheya+23 MLP verdict to every member")
    ghost_test: bool = Field(False, description="add the N-body ghost test to every member")
    ghost_orbits: float = Field(100.0, gt=0)
    ghost_max_wall_s: float | None = Field(
        None,
        gt=0,
        description="abandon a ghost test after this wall time (ghost_stable = NaN; Vynatheya "
        "et al. 2023 left out systems that took longer than 5 h)",
    )
    save_runs: bool = Field(False, description="also store the HDF5 file of every member")
    observed: dict[str, ObservedValue] = Field(
        default_factory=dict,
        description="summary column -> {value, sigma}: weight the members by these "
        "measurements (swingersim ensemble-report)",
    )

    @model_validator(mode="after")
    def _paths(self) -> EnsembleSpec:
        if len(set(self.directions)) != len(self.directions):
            raise ValueError("ensemble.directions: every direction only once")
        if self.ma01_window is not None and not self.ma01_window[0] < self.ma01_window[1]:
            raise ValueError("ensemble.ma01_window: [low, high] with low < high")
        for path in self.vary:
            parts = path.split(".")
            model = VARY_PREFIXES.get(parts[0], False)
            if model is False:
                raise ValueError(f"vary: unknown parameter '{path}' (masses., orbits., stars.)")
            if parts[0] == "masses":
                if len(parts) != 2 or not (parts[1].startswith("m") and parts[1][1:].isdigit()):
                    raise ValueError(f"vary: '{path}' should be masses.m1, masses.m2, ...")
                continue
            numeric = [k for k, f in model.model_fields.items() if _is_numeric(f.annotation)]
            if len(parts) != 3 or parts[2] not in numeric:
                raise ValueError(
                    f"vary: '{path}' should be {parts[0]}.<name>.<one of {', '.join(numeric)}>"
                )
            if parts[0] == "stars" and not (parts[1].startswith("m") and parts[1][1:].isdigit()):
                raise ValueError(f"vary: unknown star in '{path}'")
        return self


class DiscSpec(_Strict):
    """A dust disc or ring around a star or a binary that occults stars behind it (tier A
    photometry: light curves and minima). The disc is a slab of radius r_out (with an optional
    inner hole r_in) and full thickness ``thickness_au``; its light, mass and scattering are
    neglected. A star seen through it is dimmed by exp(-tau L / thickness) along each line of
    sight, L being the path length inside the slab (opaque if tau is null)."""

    host: str = Field(
        ...,
        description="star (m2) or orbit (A): the disc is centred on the star or on the centre of "
        "mass of the stars of that orbit (a circumbinary disc)",
    )
    r_out_au: float = Field(..., gt=0, description="outer radius [AU]")
    r_in_au: float = Field(0.0, ge=0, description="inner radius (hole) [AU]")
    thickness_au: float = Field(..., gt=0, description="full thickness of the slab [AU]")
    tau: float | None = Field(
        None,
        gt=0,
        description="optical depth through the full thickness, perpendicular to the disc "
        "(null = opaque)",
    )
    i: float = Field(
        90.0, ge=0, le=180, description="inclination of the disc plane [deg], as for orbits"
    )
    Omega: float = Field(0.0, description="longitude of the ascending node [deg], as for orbits")
    relative_to: str | None = Field(
        None,
        description="orbit whose plane the disc is tilted from (i and Omega are then ignored; "
        "the axis of the precession)",
    )
    tilt: float = Field(
        0.0, ge=0, le=180, description="relative_to: angle between the disc and the orbit [deg]"
    )
    tilt_node: float = Field(
        0.0,
        description="relative_to: node of the disc on the orbit plane, from the orbit's "
        "ascending node [deg] (the phase of the precession at t = 0)",
    )
    precession_period_yr: float | None = Field(
        None,
        description="relative_to: the disc normal precesses about the orbit normal with this "
        "period (negative = retrograde, the usual nodal regression of a ring around a binary)",
    )

    dynamics: Literal["prescribed", "rigid"] = Field(
        "prescribed",
        description="prescribed: precession_period_yr (or none); rigid: the disc precesses as a "
        "rigid body under the torque of the binary it surrounds (Zanazzi & Lai 2018; polar "
        "libration around eccentric binaries, Farago & Laskar 2010) or of the companion of its "
        "host star (Larwood et al. 1996); needs relative_to",
    )
    sigma_index: float = Field(
        1.0,
        ge=0,
        lt=2.5,
        description="rigid: surface density Sigma ~ r^-p of the disc (weights of its torque)",
    )
    damping_yr: float | None = Field(
        None,
        gt=0,
        description="rigid: time scale on which viscous dissipation aligns the disc with the "
        "binary (coplanar) or with its eccentricity vector (polar), whichever it librates around "
        "(Martin & Lubow 2017; Lubow & Martin 2018)",
    )
    edge_au: float | None = Field(
        None,
        gt=0,
        description="soft outer edge: the opacity falls as 1 / (1 + exp((r - r_out) / edge_au)) "
        "(epsilon Aur has a slow ingress; null = sharp edge)",
    )
    emission: bool = Field(
        False,
        description="add the thermal light of the disc: a flat reprocessing disc, "
        "T = (2/3 pi)^(1/4) (R_*/r)^(3/4) T_* (Chiang & Goldreich 1997) heated by its host "
        "star(s), emitting as a black body where it is optically thick",
    )
    albedo: float = Field(
        0.0,
        ge=0,
        le=1,
        description="single-scattering albedo of the grains: starlight scattered by the "
        "illuminated disc surface (Henyey-Greenstein phase function with asymmetry g)",
    )
    g: float = Field(0.5, gt=-1, lt=1, description="scattering asymmetry parameter")

    @model_validator(mode="after")
    def _check(self) -> DiscSpec:
        if self.dynamics == "rigid" and self.relative_to is None:
            raise ValueError("disc: dynamics rigid needs relative_to (the perturbing orbit)")
        if self.dynamics == "rigid" and self.precession_period_yr is not None:
            raise ValueError(
                "disc: dynamics rigid computes the precession; drop precession_period_yr"
            )
        if self.r_in_au >= self.r_out_au:
            raise ValueError("disc: r_in_au must be smaller than r_out_au")
        if self.edge_au is not None and self.tau is None:
            raise ValueError("disc: a soft edge (edge_au) needs a finite optical depth tau")
        if self.precession_period_yr is not None:
            if self.relative_to is None:
                raise ValueError("disc: precession_period_yr needs relative_to (the axis)")
            if self.precession_period_yr == 0:
                raise ValueError("disc: precession_period_yr must not be 0")
        return self


class Config(_Strict):
    topology: str = Field(
        ...,
        description="2+2, 3+1, 2+2+1, 3+1+1, 2+(2+1), 2+2+2, 3+3, or custom (with 'hierarchy'; up to 10 stars tested)",
    )
    hierarchy: dict[str, list[str]] | None = Field(
        None, description="custom topology: orbit -> [child 1, child 2]"
    )
    masses: dict[str, float]
    orbits: dict[str, OrbitSpec]
    stars: dict[str, StarSpec] = Field(default_factory=dict)
    discs: dict[str, DiscSpec] = Field(
        default_factory=dict,
        description="dust discs and rings that occult stars in light curves (name -> disc)",
    )
    stellar_model: StellarModelSpec = Field(default_factory=StellarModelSpec)
    physics: Physics = Field(default_factory=Physics)
    run: RunSpec
    ensemble: EnsembleSpec | None = Field(
        None, description="Monte Carlo ensemble (swingersim ensemble); ignored by swingersim run"
    )
    fit: FitSpec | None = Field(
        None,
        description="fit of the measurements (swingersim fit): times of minima, radial velocities, light curves and inclinations; ignored by swingersim run",
    )
    fit_epoch_bjd_tdb: float | None = Field(
        None,
        description="written by swingersim fit into stage2_fitted.yaml: BJD_TDB of t = 0 of the "
        "fitted elements; ignored by runs",
    )
    _base_dir: str | None = PrivateAttr(None)

    @property
    def tree(self) -> dict[str, list[str]]:
        """The hierarchy (orbit -> children) of this configuration."""
        return self.hierarchy if self.topology == "custom" else PRESETS[self.topology]

    @property
    def star_names(self) -> list[str]:
        return [f"m{k + 1}" for k in range(len(self.tree) + 1)]

    @property
    def base_dir(self) -> str | None:
        """Directory of the configuration file (relative paths are resolved against it)."""
        return self._base_dir

    @model_validator(mode="before")
    @classmethod
    def _default_run_for_fit(cls, data):
        if isinstance(data, dict) and "fit" in data and "run" not in data:
            data = {**data, "run": {"t_max_yr": 100.0}}
        return data

    @model_validator(mode="after")
    def _check_structure(self) -> Config:
        if self.topology == "custom":
            if not self.hierarchy:
                raise ValueError("topology custom needs 'hierarchy'")
        elif self.topology not in PRESETS:
            raise ValueError(
                f"unknown topology {self.topology}; use one of {list(PRESETS)} or custom"
            )
        elif self.hierarchy is not None:
            raise ValueError("'hierarchy' is only used with topology custom")
        _root, star_names = check_hierarchy(self.tree)
        if set(self.masses) != set(star_names):
            raise ValueError(f"masses must be exactly {star_names}")
        for k, m in self.masses.items():
            if m <= 0:
                raise ValueError(f"mass {k} must be positive")
        expected = set(self.tree)
        if set(self.orbits) != expected:
            raise ValueError(f"topology {self.topology} needs orbits {sorted(expected)}")
        if self.run.t_max_yr == 0:
            raise ValueError("run.t_max_yr must not be 0 (negative = backward in time)")
        for name, o in self.orbits.items():
            if o.Pdot == 0.0:
                continue
            kids = self.tree.get(name, [])
            if any(k in self.tree for k in kids):
                raise ValueError(
                    f"orbits.{name}.Pdot: only an orbit of two stars can change its period"
                )
            if o.Pdot_mode == "conservative" and (
                self.physics.mass_loss is not None or self.physics.mass_transfer is not None
            ):
                raise ValueError(
                    f"orbits.{name}.Pdot_mode conservative: the masses would change twice "
                    "(physics.mass_loss / mass_transfer); use angular_momentum"
                )
        if self.physics.gas is not None:
            bad = set(self.physics.gas.orbits) - expected
            if bad:
                raise ValueError(f"physics.gas.orbits: no orbit {sorted(bad)} in {self.topology}")
            if self.run.backend != "secular":
                raise ValueError("physics.gas needs run.backend: secular")
            if self.run.t_max_yr < 0:
                raise ValueError(
                    "physics.gas: integrate forward in time (the gas phase is the past)"
                )
        for name, e in self.run.stop_when_e_reached.items():
            if name not in expected or not 0 < e < 1:
                raise ValueError(f"run.stop_when_e_reached: bad entry {name}: {e}")
        if self.ensemble is not None:
            for path in self.ensemble.vary:
                parts = path.split(".")
                if parts[0] == "orbits" and parts[1] not in expected:
                    raise ValueError(f"ensemble.vary: no orbit '{parts[1]}' in {self.topology}")
                if parts[0] in ("masses", "stars") and parts[1] not in star_names:
                    raise ValueError(f"ensemble.vary: no star '{parts[1]}' in {self.topology}")
        unknown = set(self.stars) - set(star_names)
        if unknown:
            raise ValueError(f"unknown stars {sorted(unknown)}")
        self._check_names(star_names, expected)
        for name, d in self.discs.items():
            if d.host not in star_names and d.host not in expected:
                raise ValueError(f"discs.{name}.host: no star or orbit '{d.host}'")
            if d.relative_to is not None and d.relative_to not in expected:
                raise ValueError(f"discs.{name}.relative_to: no orbit '{d.relative_to}'")
        if self.physics.mass_transfer is not None:
            if not (self.physics.evolution and self.stellar_model.kind in ("sse", "mist")):
                raise ValueError(
                    "physics.mass_transfer needs physics.evolution: true and stellar_model.kind "
                    "sse or mist"
                )
            if self.physics.gas is not None:
                raise ValueError("physics.mass_transfer and physics.gas cannot be combined")
            if self.run.megno:
                raise ValueError("run.megno assumes constant masses; switch it off")
        ml = self.physics.mass_loss
        if ml is not None:
            bad = sorted(set(ml.rates) - set(star_names))
            if bad:
                raise ValueError(f"physics.mass_loss.rates: unknown stars {bad}")
            if ml.source == "track" and not (
                self.physics.evolution and self.stellar_model.kind in ("sse", "mist")
            ):
                raise ValueError(
                    "physics.mass_loss source track needs physics.evolution: true and "
                    "stellar_model.kind: sse or mist; or give source: rates"
                )
            if self.physics.gas is not None:
                raise ValueError("physics.mass_loss and physics.gas cannot be combined")
            if self.run.megno:
                raise ValueError(
                    "run.megno assumes constant masses; switch it off with physics.mass_loss"
                )
        offsets = {n: s.age_offset_myr for n, s in self.stars.items() if s.age_offset_myr}
        if offsets:
            if self.stellar_model.kind not in ("sse", "mist"):
                raise ValueError("stars.*.age_offset_myr needs stellar_model.kind: sse or mist")
            young = [n for n, o in offsets.items() if self.stellar_model.age_myr + o < 0]
            if young:
                raise ValueError(f"stars {young}: stellar_model.age_myr + age_offset_myr < 0")
        if self.physics.evolution:
            if self.stellar_model.kind not in ("sse", "mist"):
                raise ValueError(
                    "physics.evolution needs stellar_model.kind: sse or mist (with age_myr)"
                )
            fixed = sorted(
                n for n, s in self.stars.items() if s.R is not None or s.Teff is not None
            )
            if fixed:
                raise ValueError(
                    f"physics.evolution: stars {fixed} have a fixed radius R or Teff; these follow "
                    "from mass and age"
                )
        return self

    def _check_names(self, star_names: list[str], orbits: set[str]) -> None:
        """Star and orbit names used as keys or list entries (pydantic's extra='forbid' checks
        only the fields of the models)."""

        def check(what: str, names, allowed, kind: str) -> None:
            bad = sorted(set(names) - set(allowed))
            if bad:
                raise ValueError(f"{what}: no {kind} {bad} in {self.topology}")

        ph = self.physics
        check("physics.activity", ph.activity, star_names, "star")
        if ph.magnetic_braking is not None and ph.magnetic_braking.stars:
            check("physics.magnetic_braking.stars", ph.magnetic_braking.stars, star_names, "star")
        f = self.fit
        if f is None:
            return
        if f.rv is not None:
            check("fit.rv.stars", f.rv.stars, star_names, "star")
        check("fit.ephemerides", f.ephemerides, orbits, "orbit")
        check("fit.stage1.pairs", f.stage1.pairs, orbits, "orbit")
        if f.stage2 is not None:
            check("fit.stage2.pairs", f.stage2.pairs, orbits, "orbit")
        if f.output.predict is not None:
            check("fit.output.predict.pairs", f.output.predict.pairs, orbits, "orbit")


def load_config(path: str | Path) -> Config:
    with open(path, encoding="utf-8") as f:
        cfg = Config.model_validate(yaml.safe_load(f))
    cfg._base_dir = str(Path(path).resolve().parent)
    return cfg


def config_at(cfg: Config, res, j: int = -1) -> Config:
    """A copy of ``cfg`` whose orbits are the osculating elements of output ``j`` of an N-body
    run, e.g. to continue a run or to integrate back from its end state."""
    import math

    new = cfg.model_copy(deep=True)
    for name in list(new.orbits):
        o = res.orbits[name]
        new.orbits[name] = OrbitSpec(
            a=float(o["a"][j]),
            e=float(o["e"][j]),
            i=math.degrees(float(o["inc"][j])),
            Omega=math.degrees(float(o["Omega"][j])) % 360.0,
            omega=math.degrees(float(o["omega"][j])) % 360.0,
            M0=math.degrees(float(o["M"][j])) % 360.0,
            contact_fillout=cfg.orbits[name].contact_fillout,
        )
    new._base_dir = cfg.base_dir
    return new
