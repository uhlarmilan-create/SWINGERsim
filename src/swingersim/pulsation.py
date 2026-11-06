"""Stellar pulsations: periods from published relations and the light modulation.

A pulsating star varies its flux as m(t) = sum_k A_k sin(2 pi (t_e / P_k + phi_k)) in magnitudes,
where t_e is the star's own emission time. Using t_e (not the observing time) is what makes an
orbiting pulsator a clock: its pulsation phase lags by the light-travel time across the orbit
(phase modulation, Murphy et al. 2014, MNRAS 441, 2515).

Periods (when not given explicitly):

* ``lpv`` (Miras, semiregulars): non-linear fundamental-mode period-mass-radius relation of
  Trabucchi et al. (2021), MNRAS 500, 1575, eqs. 1-2 and Table 2 (O-rich AGB models, Z = 0.006;
  solar-metallicity periods are 10-15 per cent shorter, their Sect. 5.4). The fundamental mode
  dominates above R_dom,0 (their eq. 4); below it an overtone dominates and the period is only
  indicative.
* ``delta_scuti``: pulsation constant Q = P sqrt(rho / rho_sun) [d] of the radial modes of a
  typical delta Scuti model (fundamental 0.0329, overtones 0.0251, 0.0203, 0.0172; Breger 1979,
  PASP 91, 5, Table II, computed by A. N. Cox); the mode ranges of Poro et al. (2024,
  arXiv:2401.01091) contain these values.
* ``beta_cep``: radial fundamental mode with Q = 0.033 d; over the whole HR diagram Q of the radial
  fundamental mode lies between 0.03 and 0.04 d (Handler 2013, arXiv:1205.6407, after eq. 11).
* ``cepheid``: non-linear convective pulsation relations of Bono, Castellani & Marconi (2000),
  ApJ 529, 293: fundamental log P = a + b log Teff + c log M + d log L (Table 6, Z = 0.004, 0.008,
  0.02, interpolated in log Z; models of 5-11 Msun), first overtone log P = 10.763 + 0.672 log L
  - 3.305 log Teff (their Sect. 3); the edges of the fundamental instability strip from Table 5.
* ``spb``, ``gamma_dor`` (high-order gravity modes): a series of dipole modes equally spaced in
  period, P_k = P_0 + k Delta P (first-order asymptotic theory of Tassoul 1980, as in Wu, Li &
  Deng 2018, arXiv:1809.02620, eq. 2); P_0 and Delta P are given by the user (Delta P of SPB models is
  ~6400 s, Wu et al. 2018, Fig. 1).
* ``rr_lyrae``: pulsation equation of Bono et al. (1996), as quoted by Alcock et al. (1997, ApJ 482,
  89, eq. 2): log P0 = 11.627 + 0.823 log L - 3.506 log Teff - 0.582 log M (0.65-0.75 Msun); the
  first overtone from the period ratio P1 / P0 = 0.744 of double-mode RR Lyrae stars of the LMC
  (median of Alcock et al. 1997, Table 1).
* ``manual`` needs ``period_d`` or ``Q``.

Teff warnings: the Cepheid strip from Bono et al. (2000, Table 5), SPB stars 10000-21000 K (the
catalogue of Shi et al. 2024, arXiv:2412.03855); for the other types rough plausibility ranges.
"""

from __future__ import annotations

import itertools
import math
import warnings
from dataclasses import dataclass

import numpy as np

# Trabucchi et al. (2021), Table 2: K -> (k_b, M_bk [Msun], gamma_1, gamma_2), K = log k_b + ...
TRABUCCHI_2021 = {
    "Rb": (421.0, 2.6, 0.952, 0.114),
    "Pb": (440.0, 2.6, 0.976, -0.264),
    "Rs": (311.0, 1.0, 1.590, 0.654),
    "Ps": (388.0, 1.0, 1.808, 0.502),
    "alpha": (49.7, 2.6, -0.279, 0.544),
}
Q_RADIAL = {"F": 0.0329, "1O": 0.0251, "2O": 0.0203, "3O": 0.0172}  # Breger (1979), Table II
Q_BETA_CEP_F = 0.033  # Handler (2013): 0.03-0.04 d for the radial fundamental mode
# Bono, Castellani & Marconi (2000), Table 6: Z -> (a, b [log Teff], c [log M], d [log L])
BONO_2000_F = {
    0.004: (10.971, -3.387, -0.813, 0.929),
    0.008: (10.557, -3.279, -0.795, 0.931),
    0.02: (9.874, -3.108, -0.767, 0.942),
}
BONO_2000_1O = (10.763, 0.672, -3.305)  # log P = a + b log L + c log Teff
# Table 5: Z -> (blue edge a, b; red edge a, b), log Teff = a + b log L
BONO_2000_EDGES = {
    0.004: (3.905, -0.035, 3.995, -0.080),
    0.008: (3.919, -0.043, 4.030, -0.094),
    0.02: (3.959, -0.061, 4.125, -0.124),
}
POROS_2024_Q_RANGES = {
    "F": (0.027, 1.0),
    "1O": (0.021, 0.027),
    "2O": (0.018, 0.021),
    "3O": (0.0, 0.018),
}
RR_P1_P0 = 0.744  # Alcock et al. (1997), double-mode RR Lyrae in the LMC

# rough Teff windows [K] for plausibility warnings only
TEFF_WINDOWS = {
    "delta_scuti": (6300.0, 8600.0),
    "gamma_dor": (6700.0, 7600.0),
    "beta_cep": (20000.0, 33000.0),
    "spb": (10000.0, 21000.0),  # Shi et al. (2024)
    "rr_lyrae": (5900.0, 7600.0),
    "lpv": (0.0, 4300.0),
}


def _trabucchi_param(key: str, m: float) -> float:
    kb, mbk, g1, g2 = TRABUCCHI_2021[key]
    g = g1 if m < mbk else g2
    return math.log10(kb) + g * math.log10(m / mbk)


def lpv_fundamental_period(m: float, r: float) -> float:
    """Non-linear fundamental-mode period [d] of an O-rich AGB star of current mass m [Msun]
    and radius r [Rsun] (Trabucchi et al. 2021, eqs. 1-2, Table 2)."""
    log_rb, log_pb = _trabucchi_param("Rb", m), _trabucchi_param("Pb", m)
    log_rs, log_ps = _trabucchi_param("Rs", m), _trabucchi_param("Ps", m)
    alpha = _trabucchi_param("alpha", m)
    lr = math.log10(r)
    if lr < log_rb:
        return 10.0 ** (log_pb + alpha * (lr - log_rb))
    if lr < log_rs:
        beta = (log_pb - log_ps) / (log_rb - log_rs)
        return 10.0 ** (log_pb + beta * (lr - log_rb))
    return 10.0**log_ps


def lpv_dominant_fm_radius(m: float) -> float:
    """Radius [Rsun] above which the fundamental mode dominates (Trabucchi et al. 2021, eq. 4)."""
    lm = math.log10(m)
    return 10.0 ** (2.130 + 1.150 * lm - 0.496 * lm * lm)


def q_period(q: float, m: float, r: float) -> float:
    """Period [d] of a mode with pulsation constant Q [d]: P = Q / sqrt(rho / rho_sun)."""
    return q / math.sqrt(m / r**3)


def _bono_interp(table: dict, z: float) -> tuple[float, ...]:
    """Coefficients of a Bono et al. (2000) table interpolated linearly in log Z (clamped)."""
    zs = sorted(table)
    lz = min(max(math.log10(z), math.log10(zs[0])), math.log10(zs[-1]))
    for z0, z1 in itertools.pairwise(zs):
        if lz <= math.log10(z1):
            w = (lz - math.log10(z0)) / (math.log10(z1) - math.log10(z0))
            return tuple((1 - w) * a + w * b for a, b in zip(table[z0], table[z1], strict=True))
    return table[zs[-1]]


def cepheid_period(l: float, teff: float, m: float, z: float = 0.02, mode: str = "F") -> float:
    """Period [d] of a classical Cepheid (Bono, Castellani & Marconi 2000): fundamental mode from
    Table 6 (interpolated in log Z), first overtone from their Sect. 3 relation."""
    if mode == "F":
        a, b, c, d = _bono_interp(BONO_2000_F, z)
        return 10.0 ** (a + b * math.log10(teff) + c * math.log10(m) + d * math.log10(l))
    if mode == "1O":
        a, b, c = BONO_2000_1O
        return 10.0 ** (a + b * math.log10(l) + c * math.log10(teff))
    raise ValueError("cepheid: mode must be F or 1O")


def cepheid_strip(l: float, z: float = 0.02) -> tuple[float, float]:
    """Teff [K] of the red and blue edges of the fundamental Cepheid instability strip at
    luminosity l [Lsun] (Bono et al. 2000, Table 5)."""
    ab, bb, ar, br = _bono_interp(BONO_2000_EDGES, z)
    return 10.0 ** (ar + br * math.log10(l)), 10.0 ** (ab + bb * math.log10(l))


def rr_lyrae_fundamental_period(l: float, teff: float, m: float) -> float:
    """Fundamental period [d] of an RR Lyrae star (Bono et al. 1996 via Alcock et al. 1997)."""
    return 10.0 ** (
        11.627 + 0.823 * math.log10(l) - 3.506 * math.log10(teff) - 0.582 * math.log10(m)
    )


@dataclass
class Mode:
    period_d: float
    amplitude_mag: float
    phase: float  # cycles at t_e = 0


@dataclass
class Pulsator:
    star: int  # index 0..N-1
    kind: str
    modes: list[Mode]
    notes: list[str]

    def delta_mag(self, t_e_day: np.ndarray | float) -> np.ndarray | float:
        out = 0.0
        for md in self.modes:
            out = out + md.amplitude_mag * np.sin(
                2.0 * math.pi * (t_e_day / md.period_d + md.phase)
            )
        return out

    def flux_factor(self, t_e_day):
        return 10.0 ** (-0.4 * self.delta_mag(t_e_day))


def metallicity(cfg) -> float:
    """Metal mass fraction Z of the stellar model (MIST: protosolar Z = 0.0142 scaled by [Fe/H],
    Choi et al. 2016; SSE: its z; otherwise 0.02)."""
    sm = cfg.stellar_model
    if sm.kind == "mist":
        return 0.0142 * 10.0**sm.feh
    if sm.kind == "sse":
        return sm.z
    return 0.02


def main_period(spec, m: float, props, z: float = 0.02) -> tuple[float, list[str]]:
    """Period [d] of the main mode of a pulsation spec for a star of mass m, properties
    ``props`` (StarProperties) and metallicity z, with notes (warnings) about its validity."""
    notes: list[str] = []
    r, teff, lum = props.radius_rsun, props.teff_k, props.luminosity_lsun
    kind = spec.type
    if kind == "cepheid" and lum is not None and teff is not None:
        red, blue = cepheid_strip(lum, z)
        if not red <= teff <= blue:
            notes.append(
                f"Teff {teff:.0f} K is outside the fundamental Cepheid instability strip at this "
                f"luminosity ({red:.0f}-{blue:.0f} K, Bono et al. 2000)"
                + (
                    "; first-overtone Cepheids lie near and beyond its blue edge"
                    if spec.mode == "1O" and teff > blue
                    else ""
                )
            )
    window = TEFF_WINDOWS.get(kind)
    if window and teff is not None and not window[0] <= teff <= window[1]:
        notes.append(
            f"Teff {teff:.0f} K is outside the usual range of {kind} stars "
            f"({window[0]:.0f}-{window[1]:.0f} K)"
        )
    if spec.period_d is not None:
        if kind == "spb" and not 0.14 <= spec.period_d <= 6.5:
            notes.append("SPB periods of Shi et al. (2024) lie between 0.14 and 6.5 d")
        return spec.period_d, notes
    if spec.Q is not None:
        return q_period(spec.Q, m, r), notes
    if kind == "lpv":
        if not 0.6 <= m <= 7.0:
            notes.append(f"mass {m:.2f} Msun outside the models of Trabucchi et al. 2021 (0.6-7)")
        r_dom = lpv_dominant_fm_radius(m)
        if r < r_dom:
            notes.append(
                f"R = {r:.0f} Rsun < {r_dom:.0f} Rsun: an overtone, not the fundamental mode, "
                "dominates (Trabucchi et al. 2021, eq. 4)"
            )
        return lpv_fundamental_period(m, r), notes
    if kind == "delta_scuti":
        return q_period(Q_RADIAL[spec.mode], m, r), notes
    if kind == "beta_cep":
        if spec.mode != "F":
            raise ValueError("beta_cep: only the radial fundamental mode has a default Q; give Q")
        return q_period(Q_BETA_CEP_F, m, r), notes
    if kind == "cepheid":
        if lum is None or teff is None:
            raise ValueError("cepheid: the stellar model gives no L and Teff")
        if not 4.5 <= m <= 12.0:
            notes.append(f"mass {m:.2f} Msun outside the Cepheid models of Bono et al. (5-11)")
        if spec.mode not in ("F", "1O"):
            raise ValueError("cepheid: mode must be F or 1O")
        return cepheid_period(lum, teff, m, z, spec.mode), notes
    if kind == "rr_lyrae":
        if lum is None or teff is None:
            raise ValueError("rr_lyrae: the stellar model gives no L and Teff")
        if not 0.55 <= m <= 0.8:
            notes.append(f"mass {m:.2f} Msun outside the RR Lyrae range of Bono et al. (0.65-0.75)")
        p0 = rr_lyrae_fundamental_period(lum, teff, m)
        if spec.mode == "F":
            return p0, notes
        if spec.mode == "1O":
            return p0 * RR_P1_P0, notes
        raise ValueError("rr_lyrae: mode must be F or 1O")
    raise ValueError(f"pulsation type {kind}: give period_d (or Q for radial modes)")


def pulsators(cfg, masses, props) -> list[Pulsator]:
    """Pulsating stars of a configuration (current masses and stellar properties at t = 0)."""
    out = []
    for k, m in enumerate(masses):
        name = f"m{k + 1}"
        st = cfg.stars.get(name)
        spec = st.pulsation if st is not None else None
        if spec is None:
            continue
        p0, notes = main_period(spec, float(m), props[k], metallicity(cfg))
        modes = [Mode(p0, spec.amplitude_mag, spec.phase)]
        if spec.period_spacing_s is not None:  # g-mode series, P_k = P_0 + k Delta P
            dp = spec.period_spacing_s / 86400.0
            for j in range(1, spec.n_modes):
                modes.append(Mode(p0 + j * dp, spec.amplitude_mag, spec.phase + 0.618 * j))
        for md in spec.extra_modes:
            per = md.period_d if md.period_d is not None else p0 / md.frequency_ratio
            modes.append(Mode(per, md.amplitude_mag, md.phase))
        for n in notes:
            warnings.warn(f"star {name}: {n}", stacklevel=2)
        out.append(Pulsator(k, spec.type, modes, notes))
    return out


def period_table(cfg) -> list[dict]:
    """Main periods of all pulsating stars at t = 0, and their change dP/dt from the stellar
    model (evolving models only), for the 'pulsation' command."""
    import copy

    from . import massloss
    from .hierarchy import build_tree
    from .stellar import star_properties

    cfg = massloss.resolve_masses(cfg)
    m_lab = build_tree(cfg)[1]
    rows = []
    for k, m in enumerate(m_lab):
        name = f"m{k + 1}"
        st = cfg.stars.get(name)
        if st is None or st.pulsation is None:
            continue
        props = star_properties(cfg, name, float(m))
        m_now = massloss.initial_masses(cfg, m_lab)[k]
        z = metallicity(cfg)
        p, notes = main_period(st.pulsation, float(m_now), props, z)
        row = {
            "star": name,
            "type": st.pulsation.type,
            "mode": st.pulsation.mode,
            "period_d": p,
            "R_rsun": props.radius_rsun,
            "Teff_K": props.teff_k,
            "L_lsun": props.luminosity_lsun,
            "dP_dt_s_per_yr": None,
            "notes": notes,
        }
        sm = cfg.stellar_model
        if sm.kind in ("sse", "mist") and st.pulsation.period_d is None:
            dt = 1e-3 if sm.kind == "mist" else 1e-2  # Myr
            c2 = copy.deepcopy(cfg)
            c2.stellar_model.age_myr = sm.age_myr + dt
            pr2 = star_properties(c2, name, float(m))
            m2 = massloss.initial_masses(c2, m_lab)[k]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                p2, _ = main_period(st.pulsation, float(m2), pr2, z)
            row["dP_dt_s_per_yr"] = (p2 - p) * 86400.0 / (dt * 1e6)
        rows.append(row)
    return rows
