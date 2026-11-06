"""Secular (orbit-averaged) backend for 2+2 quadruples, Vokrouhlicky (2016) MNRAS 461, 3964.

The state of each orbit X in {A, B, out} is described by the non-singular vectorial elements
of Vokrouhlicky (2016), Sect. 2.2: the normalised angular momentum K = sqrt(1 - e^2) k and the
Laplace vector pi = e p (k orbit normal, p unit vector to the periastron). The semi-major axes
are constant in the point-mass problem. The equations of motion (his eqs. 12-17)

    dK_X/dt  = -(1/L_X) [ pi_X x dH/dpi_X + K_X x dH/dK_X ]
    dpi_X/dt = -(1/L_X) [ K_X x dH/dpi_X + pi_X x dH/dK_X ],   L_X = mu_X n_X a_X^2

are driven by the averaged interaction energy H = H2 + H3: the quadrupole terms (his eq. 18)
and octupole terms (eqs. 21-25) of two triples, A perturbed by the centre of mass of B and B
perturbed by the centre of mass of A. Added to H are, for every orbit, terms that depend on
its eccentricity only and give the apsidal advance of

* general relativity (1PN, for A, B and out), and
* the tidal and rotational distortion of the stars of A and B (Sterne 1939), with spins
  aligned with the orbit of their binary,

so that H stays conserved and checks the integration. With tidal dissipation (time lag > 0),
a_A, a_B, e_A, e_B and the stellar spin rates evolve in addition according to the
constant-time-lag theory (Hut 1981; Leconte et al. 2010, see swingersim.tides).

Limitations: first-order secular theory, so the short-period terms and the "Hint'"
coupling of the two binaries (order epsilon^4, his eq. 10) are neglected; mean-motion
resonances are absent; stellar spins stay aligned with their orbits.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .. import units
from ..tides import slaved_spin, sync_time_yr
from .kernel import _lock_spins, dissipative_rhs, secular_rhs


def _cross(a, b):
    """Cross product of two 3-vectors (np.cross has a large per-call overhead)."""
    return np.array(
        (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
    )


def _unit_normal(inc: float, Omega: float) -> np.ndarray:
    return np.array([np.sin(inc) * np.sin(Omega), -np.sin(inc) * np.cos(Omega), np.cos(inc)])


def _periastron_dir(inc: float, Omega: float, omega: float) -> np.ndarray:
    cO, sO, cw, sw, ci, si = (
        np.cos(Omega),
        np.sin(Omega),
        np.cos(omega),
        np.sin(omega),
        np.cos(inc),
        np.sin(inc),
    )
    return np.array([cO * cw - sO * sw * ci, sO * cw + cO * sw * ci, sw * si])


# e = 0 exactly is a fixed point of the secular equations (the eccentricity vector stays zero
# even where Kozai-Lidov cycles are unstable); an N-body run grows them from round-off.
# The secular backends start such an orbit at this eccentricity.
E_SEED = 1e-6


def vectors_from_elements(e: float, inc: float, Omega: float, omega: float):
    """(K, pi) from eccentricity and angles [rad]."""
    eta = np.sqrt(1.0 - e * e)
    return eta * _unit_normal(inc, Omega), e * _periastron_dir(inc, Omega, omega)


def elements_from_vectors(K: np.ndarray, pi: np.ndarray) -> tuple[float, float, float, float]:
    """(e, inc, Omega, omega) [rad] from the vectorial elements."""
    k = K / np.linalg.norm(K)
    e = float(np.linalg.norm(pi))
    inc = float(np.arccos(np.clip(k[2], -1.0, 1.0)))
    Omega = float(np.mod(np.arctan2(k[0], -k[1]), 2.0 * np.pi))
    node = np.array([np.cos(Omega), np.sin(Omega), 0.0])
    if e > 0:
        omega = float(np.mod(np.arctan2(pi @ np.cross(k, node), pi @ node), 2.0 * np.pi))
    else:
        omega = 0.0
    return e, inc, Omega, omega


@dataclass
class Binary:
    """Physical data of an inner binary for the precession and tidal terms."""

    m: tuple[float, float]
    r_au: tuple[float, float]
    apsidal_k2: tuple[float, float]
    rg2: tuple[float, float]
    lag_yr: tuple[float, float]
    tides: bool  # include the conservative tidal / rotational distortion


class SecularTwoTwo:
    """Right-hand side, energy and angular momentum of the 2+2 secular system."""

    def __init__(
        self,
        masses: np.ndarray,
        a: dict[str, float],
        binaries: dict[str, Binary],
        spins: dict[str, tuple[float, float]],
        gr: bool,
        octupole: bool = True,
    ):
        self._set_mass_terms(masses)
        self.a_out = a["out"]
        self.a0 = {"A": a["A"], "B": a["B"]}
        self.binaries = binaries
        self.spin0 = spins  # initial spin rates [rad/yr] of the two stars of A and of B
        self.gr = gr
        self.octupole = octupole
        self.dissipative = any(b.tides and max(b.lag_yr) > 0 for b in binaries.values())
        self._cached = None  # kernel coefficients (constant without dissipation)
        self._diss_args = None  # constant arrays of the compiled dissipative right-hand side

    def _set_mass_terms(self, masses) -> None:
        masses = np.asarray(masses, float)
        self.m = masses
        self.MA, self.MB = masses[0] + masses[1], masses[2] + masses[3]
        self.MAB = self.MA + self.MB
        self.muA = masses[0] * masses[1] / self.MA
        self.muB = masses[2] * masses[3] / self.MB
        self.mu = self.MA * self.MB / self.MAB
        self.XAa, self.XAb = masses[0] / self.MA, masses[1] / self.MA
        self.XBa, self.XBb = masses[2] / self.MB, masses[3] / self.MB

    def set_masses(self, masses, y=None) -> None:
        """New stellar masses (slow isotropic wind, operator splitting): every orbit keeps its
        eccentricity vector and a M (Hadjidemetriou 1963); the semi-major axes are rescaled in
        the system and, for the dissipative model, in the state y (in place)."""
        old = (self.MA, self.MB, self.MAB)
        self._set_mass_terms(masses)
        fA, fB, fout = old[0] / self.MA, old[1] / self.MB, old[2] / self.MAB
        self.a0 = {"A": self.a0["A"] * fA, "B": self.a0["B"] * fB}
        self.a_out *= fout
        if self.dissipative and y is not None:
            y[18] *= fA
            y[19] *= fB
        self._cached = None
        self._diss_args = None

    def scale_orbits(self, y, factors: dict[str, float]) -> None:
        """Multiply the semi-major axes of the named inner orbits (mass transfer), in the
        system and, for the dissipative model, in the state y (in place)."""
        if not factors:
            return
        a0 = dict(self.a0)
        for n, f in factors.items():
            if n == "out":
                self.a_out *= f
                continue
            a0[n] *= f
            if self.dissipative and y is not None:
                y[{"A": 18, "B": 19}[n]] *= f
        self.a0 = a0
        self._cached = None
        self._diss_args = None

    def mass_state(self) -> dict:
        """Masses and semi-major axes in force (snapshot for the outputs)."""
        return {"m": self.m.copy(), "a0": dict(self.a0), "a_out": self.a_out}

    def restore_mass_state(self, st: dict) -> None:
        self._set_mass_terms(st["m"])
        self.a0, self.a_out = dict(st["a0"]), st["a_out"]
        self._cached = None
        self._diss_args = None

    def set_structure(self, binaries: dict[str, Binary], spins=None) -> None:
        """Replace the stellar structure (evolving stars, operator splitting); spins: new
        constant spin rates of the conservative model."""
        self.binaries = binaries
        if spins is not None:
            self.spin0 = spins
        self._cached = None
        self._diss_args = None

    # ------------------------------------------------------------------ helpers
    def _L(self, mu, M, a):
        return mu * np.sqrt(units.G * M * a)  # mu n a^2

    def _coeffs(self, aA, aB):
        G, a = units.G, self.a_out
        c2A = G * self.muA * self.MB / (8.0 * a) * (aA / a) ** 2
        c2B = G * self.MA * self.muB / (8.0 * a) * (aB / a) ** 2
        c3A = 15.0 / 64.0 * G * self.muA * self.MB / a * (aA / a) ** 3 * (self.XAa - self.XAb)
        c3B = 15.0 / 64.0 * G * self.MA * self.muB / a * (aB / a) ** 3 * (self.XBb - self.XBa)
        return c2A, c2B, c3A, c3B

    def _precession_terms(self, name, a, spin):
        """Coefficients (alpha, beta, gamma, delta) of the e-only energy of orbit `name`:

        H_x = -alpha/eta + L n [ beta * I1(eta) + gamma * I2(eta) ],
        I1 = -2.625/9 eta^-9 + 1.75/7 eta^-7 - 0.125/5 eta^-5,  I2 = -eta^-3 / 3,
        whose derivative dH/deta = L * domega/dt gives the GR and Sterne apsidal rates.
        """
        G = units.G
        if name == "out":
            M, mu = self.MAB, self.mu
        elif name == "A":
            M, mu = self.MA, self.muA
        else:
            M, mu = self.MB, self.muB
        L = self._L(mu, M, a)
        n = np.sqrt(G * M / a**3)
        alpha = L * 3.0 * (G * M) ** 1.5 / (units.C**2 * a**2.5) if self.gr else 0.0
        beta = gamma = 0.0
        b = self.binaries.get(name)
        if b is not None and b.tides:
            for i in (0, 1):
                j = 1 - i
                q = b.m[j] / b.m[i]
                x5 = b.apsidal_k2[i] * (b.r_au[i] / a) ** 5
                beta += 15.0 * q * x5
                gamma += (1.0 + q) * (spin[i] / n) ** 2 * x5
        return L, n, alpha, beta, gamma

    @staticmethod
    def _I(eta):
        I1 = -2.625 / 9.0 * eta**-9 + 1.75 / 7.0 * eta**-7 - 0.125 / 5.0 * eta**-5
        dI1 = (2.625 - 1.75 * eta**2 + 0.125 * eta**4) * eta**-10
        I2 = -(eta**-3) / 3.0
        dI2 = eta**-4
        return I1, dI1, I2, dI2

    def unpack(self, y):
        v = {
            "KA": y[0:3],
            "pA": y[3:6],
            "KB": y[6:9],
            "pB": y[9:12],
            "K": y[12:15],
            "p": y[15:18],
        }
        if self.dissipative:
            v["aA"], v["aB"] = y[18], y[19]
            v["sA"] = (y[20], y[21])
            v["sB"] = (y[22], y[23])
        else:
            v["aA"], v["aB"] = self.a0["A"], self.a0["B"]
            v["sA"], v["sB"] = self.spin0["A"], self.spin0["B"]
        return v

    # ------------------------------------------------------------------ energy and gradients
    def energy_and_gradients(self, v):
        KA, pA, KB, pB, K, p = v["KA"], v["pA"], v["KB"], v["pB"], v["K"], v["p"]
        c2A, c2B, c3A, c3B = self._coeffs(v["aA"], v["aB"])
        e2 = p @ p
        eta2 = 1.0 - e2
        eta = np.sqrt(eta2)
        g = {k: np.zeros(3) for k in ("KA", "pA", "KB", "pB", "K", "p")}
        H = 0.0
        for Ki, pi_, c2, c3, sk, sp in (
            (KA, pA, c2A, c3A, "KA", "pA"),
            (KB, pB, c2B, c3B, "KB", "pB"),
        ):
            ei2 = pi_ @ pi_
            a_ = Ki @ K
            b_ = pi_ @ K
            # quadrupole: H2 = -c2 [eta^-3 (6 ei2 - 1) + eta^-5 (3 a^2 - 15 b^2)]
            em3, em5, em7, em9 = eta**-3, eta**-5, eta**-7, eta**-9
            H += -c2 * (em3 * (6 * ei2 - 1) + em5 * (3 * a_**2 - 15 * b_**2))
            g[sp] += -c2 * (12.0 * em3 * pi_ - 30.0 * em5 * b_ * K)
            g[sk] += -c2 * (6.0 * em5 * a_ * K)
            g["K"] += -c2 * em5 * (6.0 * a_ * Ki - 30.0 * b_ * pi_)
            g["p"] += -c2 * (3.0 * em5 * (6 * ei2 - 1) + 5.0 * em7 * (3 * a_**2 - 15 * b_**2)) * p
            if self.octupole and c3 != 0.0:
                # H3 = -c3 [ eta^-5 (1 - 8 ei2) d + eta^-7 d (35 b^2 - 5 a^2) - 10 eta^-7 a b c ]
                d = pi_ @ p
                c = p @ Ki
                w = 35 * b_**2 - 5 * a_**2
                H += -c3 * (em5 * (1 - 8 * ei2) * d + em7 * d * w - 10.0 * em7 * a_ * b_ * c)
                g[sp] += -c3 * (
                    em5 * (-16.0 * d * pi_ + (1 - 8 * ei2) * p)
                    + em7 * (w * p + 70.0 * d * b_ * K)
                    - 10.0 * em7 * a_ * c * K
                )
                g[sk] += -c3 * (-10.0 * em7 * d * a_ * K - 10.0 * em7 * (b_ * c * K + a_ * b_ * p))
                g["K"] += -c3 * (
                    em7 * d * (-10.0 * a_ * Ki + 70.0 * b_ * pi_)
                    - 10.0 * em7 * (b_ * c * Ki + a_ * c * pi_)
                )
                g["p"] += -c3 * (
                    (1 - 8 * ei2) * (5.0 * em7 * d * p + em5 * pi_)
                    + w * (7.0 * em9 * d * p + em7 * pi_)
                    - 10.0 * a_ * b_ * (7.0 * em9 * c * p + em7 * Ki)
                )
        # e-only precession terms (GR, stellar distortion), written with eta = |K|
        for name, sk, a_x, spin in (
            ("A", "KA", v["aA"], v["sA"]),
            ("B", "KB", v["aB"], v["sB"]),
            ("out", "K", self.a_out, (0.0, 0.0)),
        ):
            L, n, alpha, beta, gamma = self._precession_terms(name, a_x, spin)
            Kx = v[sk]
            ex = np.linalg.norm(Kx)
            I1, dI1, I2, dI2 = self._I(ex)
            H += -alpha / ex + L * n * (beta * I1 + gamma * I2)
            dH = alpha / ex**2 + L * n * (beta * dI1 + gamma * dI2)
            g[sk] += dH * Kx / ex
        return H, g

    def _dissipative_arrays(self):
        cst = np.array([units.G, units.C, *self.m], dtype=float)
        names = ("A", "B")
        r = np.array([self.binaries[n].r_au for n in names], dtype=float)
        k2 = np.array([self.binaries[n].apsidal_k2 for n in names], dtype=float)
        rg2 = np.array([self.binaries[n].rg2 for n in names], dtype=float)
        lag = np.array([self.binaries[n].lag_yr for n in names], dtype=float)
        tides = np.array([self.binaries[n].tides for n in names], dtype=np.bool_)
        # backward runs: spins much faster to synchronise than |t_max| held pseudo-synchronous
        span = float(getattr(self, "backward_span", 0.0))
        slaved = np.zeros((2, 2))
        self.slaved = set()
        for b, n in enumerate(names):
            bn = self.binaries[n]
            if not bn.tides:
                continue
            for i in (0, 1):
                ts = sync_time_yr(bn.m[i], bn.m[1 - i], bn.r_au[i], bn.apsidal_k2[i], bn.rg2[i],
                                  bn.lag_yr[i], self.a0[n])  # fmt: skip
                if slaved_spin(span, ts):
                    slaved[b, i] = 1.0
                    self.slaved.add(2 * b + i)
        return (cst, r, k2, rg2, lag, tides, float(self.a_out), bool(self.gr), bool(self.octupole),
                slaved)  # fmt: skip

    # ------------------------------------------------------------------ equations of motion
    def _kernel_coefficients(self, v) -> tuple[np.ndarray, np.ndarray]:
        c = np.array(self._coeffs(v["aA"], v["aB"]))
        prec = np.array(
            [
                self._precession_terms(name, a_x, spin)
                for name, a_x, spin in (
                    ("A", v["aA"], v["sA"]),
                    ("B", v["aB"], v["sB"]),
                    ("out", self.a_out, (0.0, 0.0)),
                )
            ]
        )
        return c, prec

    def rhs(self, _t, y):
        """Equations of motion, evaluated by the compiled kernels (swingersim.secular.kernel):
        with tidal dissipation the full state of 24 (dissipative_rhs), otherwise the 18
        vectorial elements with coefficients that stay constant until the masses, the orbits
        or the structure change (secular_rhs)."""
        if self.dissipative:
            if self._diss_args is None:
                self._diss_args = self._dissipative_arrays()
            return dissipative_rhs(np.ascontiguousarray(y, dtype=float), *self._diss_args)
        if self._cached is None:
            v = self._kernel_coefficients(self.unpack(y))
            S = np.array([self._spin_momentum(n, self.spin0[n]) for n in ("A", "B")])
            self._cached = (*v, S)
        c, prec, S = self._cached
        yy = np.ascontiguousarray(y[:18], dtype=float)
        out = secular_rhs(yy, c, prec, self.octupole)[1]
        if S.any():  # the constant spins locked to their binary (kernel._lock_spins)
            out = np.array(out, dtype=float)
            for b in range(2):
                if S[b]:
                    eta = float(np.linalg.norm(yy[6 * b : 6 * b + 3]))
                    _lock_spins(yy, out, 6 * b, float(prec[b, 0]) * eta, float(S[b]))
        return out

    def _spin_momentum(self, name: str, spins) -> float:
        """Spin angular momentum of the two stars of binary name (aligned spins)."""
        b = self.binaries[name]
        return float(sum(b.rg2[i] * b.m[i] * b.r_au[i] ** 2 * spins[i] for i in (0, 1)))

    def angular_momentum(self, y) -> np.ndarray:
        """Total angular momentum (orbits + aligned spins) [Msun AU^2 / yr]."""
        v = self.unpack(y)
        L = (
            self._L(self.muA, self.MA, v["aA"]) * v["KA"]
            + self._L(self.muB, self.MB, v["aB"]) * v["KB"]
            + self._L(self.mu, self.MAB, self.a_out) * v["K"]
        )
        for name, sk, ss in (("A", "KA", "sA"), ("B", "KB", "sB")):
            b = self.binaries[name]
            k = v[sk] / np.linalg.norm(v[sk])
            for i in (0, 1):
                L = L + b.rg2[i] * b.m[i] * b.r_au[i] ** 2 * v[ss][i] * k
        return L


def integrate(system, y0: np.ndarray, times: np.ndarray, rtol: float = 1e-10):
    """Integrate a secular system (SecularTwoTwo or swingersim.secular.nested.SecularNested: any
    object with rhs(t, y)) from times[0] to times[-1]; returns the dense solution object of
    solve_ivp."""
    return solve_ivp(
        system.rhs,
        (times[0], times[-1]),
        y0,
        method="DOP853",
        rtol=rtol,
        atol=rtol * 1e-2,
        dense_output=True,
    )
