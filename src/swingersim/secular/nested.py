"""Secular (orbit-averaged) dynamics of any hierarchy of nested binaries.

Used for 3+1 quadruples (Hamers et al. 2015, MNRAS 449, 4221) and for hierarchies with a
fifth star; 2+2 quadruples keep the dedicated swingersim.secular.twotwo, of which this is the
generalisation. Each orbit X has the vectorial elements K_X = sqrt(1 - e^2) k and
pi_X = e p (Vokrouhlicky 2016) and the equations of motion

    dK_X/dt  = -(1/L_X) [ pi_X x dH/dpi_X + K_X x dH/dK_X ]
    dpi_X/dt = -(1/L_X) [ K_X x dH/dpi_X + pi_X x dH/dK_X ],   L_X = mu_X n_X a_X^2.

The averaged interaction energy H is a sum of binary-pair terms, one for every orbit k and
every orbit l that contains it (the "binary pairs" of Hamers et al. 2015 and Hamers &
Portegies Zwart 2016):

* quadrupole, c2 = G mu_k M_o a_k^2 / (8 a_l^3), with M_o the mass of the other component of
  orbit l (the one that does not contain k);
* octupole, c3 = (15/64) G mu_k M_o a_k^3 / a_l^4 (X_k1 - X_k2) s, with X the mass fractions of
  the components of k and s = +1 (-1) if k lies in the first (second) component of l.

For a 2+2 quadruple these are exactly the terms of swingersim.secular.twotwo. Neglected: the
"cross" (binary-triplet) terms that depend on three orbits at once, of octupole order
(a_k/a_m)^2 (a_m/a_l)^3 for k inside m inside l (Hamers et al. 2015, eq. 2), the terms of
hexadecupole and higher order.

Tidal dissipation (time lag > 0) acts in every orbit of two single stars, as in the 2+2
backend: the semi-major axis, the eccentricity and the aligned spin rates of its stars evolve by
the constant-time-lag equations (Hut 1981; Leconte et al. 2010; swingersim.tides). The state then
holds, after the vectors of all orbits, the semi-major axes of these binaries and the spins of
their stars.

Optionally a parametric gas phase (physics.gas): while gas is present each listed orbit
migrates (da/dt = -a/tau_a), circularises (de/dt = -e/tau_e) and its normal turns towards the
disc normal n (dk/dt = (n - (n.k) k)/tau_i), with the gas density constant until t_end and
then decaying as exp(-(t - t_end)/t_decay). The semi-major axes are then state variables. Added to H, as in the 2+2 code, are the
e-only terms of the relativistic and stellar-distortion apsidal motion of every orbit.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .. import units
from ..hierarchy import Node, mass_of
from ..tides import slaved_spin, sync_time_yr
from .kernel import _lock_spins, nested_dissipative_rhs, secular_rhs_pairs
from .twotwo import integrate

__all__ = ["SecularNested", "StarStructure", "integrate"]


def _stars(child) -> list[int]:
    return child.stars() if isinstance(child, Node) else [child]


@dataclass
class StarStructure:
    r_au: float
    apsidal_k2: float
    spin: float  # rad/yr
    rg2: float = 0.0  # I / (M R^2)
    lag_yr: float = 0.0  # constant tidal time lag [yr] (> 0: dissipation)


class SecularNested:
    def __init__(
        self,
        root: Node,
        masses: np.ndarray,
        structure: dict[int, StarStructure] | None = None,
        gr: bool = True,
        octupole: bool = True,
        gas=None,
        gw: bool = False,
    ):
        self.root, self.masses = root, masses
        self.nodes = root.nodes()
        self.names = [n.name for n in self.nodes]
        self.index = {n.name: k for k, n in enumerate(self.nodes)}
        self.gr, self.octupole = gr, octupole
        self.structure = structure or {}
        self.a0 = np.array([n.elements.a for n in self.nodes])
        pk, pl, others, signs = [], [], [], []
        for k, node in enumerate(self.nodes):
            stars_k = set(node.stars())
            for l, outer in enumerate(self.nodes):
                if outer is node or not stars_k < set(outer.stars()):
                    continue
                first = stars_k <= set(_stars(outer.left))
                pk.append(k)
                pl.append(l)
                others.append(outer.right if first else outer.left)
                signs.append(1.0 if first else -1.0)
        self.pk, self.pl = np.array(pk, dtype=np.int64), np.array(pl, dtype=np.int64)
        self._others, self._signs = others, np.array(signs)
        self._set_mass_terms(masses)
        self.prec, self.c2, self.c3 = self.coefficients(self.a0)
        self.L = list(self.prec[:, 0])
        # gas phase: damping time scales per orbit (inf = none) and the disc normal
        self.gas = gas
        self.n_orb = len(self.nodes)
        # tidal dissipation: two-star orbits with a time lag; their a and spins join the state
        self.gw = gw
        self.diss = [
            q
            for q, nd in enumerate(self.nodes)
            if not isinstance(nd.left, Node)
            and not isinstance(nd.right, Node)
            and (
                gw
                or any(
                    self.structure.get(s) is not None and self.structure[s].lag_yr > 0
                    for s in (nd.left, nd.right)
                )
            )
        ]
        if self.diss and gas is not None:
            raise ValueError("tidal dissipation and physics.gas cannot be combined (secular)")
        n6, nd = 6 * self.n_orb, len(self.diss)
        self.a_index = {q: n6 + i for i, q in enumerate(self.diss)}
        self.spin_index = {q: n6 + nd + 2 * i for i, q in enumerate(self.diss)}
        self.tide = None  # (rho0 [Msun/AU^3], pole unit vector): Galactic tide
        self._kargs = None  # arguments of the compiled dissipative right-hand side
        if gas is not None:
            self.tau = np.full((self.n_orb, 3), np.inf)
            for name, d in gas.orbits.items():
                k = self.index[name]
                for j, v in enumerate((d.tau_a_yr, d.tau_e_yr, d.tau_i_yr)):
                    if v is not None:
                        self.tau[k, j] = v
            self.disc_normal = None  # set by the backend (invariable plane or given plane)

    def _set_mass_terms(self, masses) -> None:
        masses = np.asarray(masses, float)
        self.masses = masses
        self.mk = np.array(
            [(mass_of(n.left, masses), mass_of(n.right, masses)) for n in self.nodes]
        )
        mk1, mk2 = self.mk[self.pk, 0], self.mk[self.pk, 1]
        self._mo = np.array([mass_of(o, masses) for o in self._others])
        self._xs = (mk1 - mk2) / (mk1 + mk2) * self._signs

    def orbit_masses(self) -> np.ndarray:
        return self.mk[:, 0] + self.mk[:, 1]

    @property
    def dissipative(self) -> bool:
        return bool(self.diss)

    def extra_state(self) -> list[float]:
        """Initial values of the dissipative state: a of the binaries, then their spins."""
        out = [float(self.a0[q]) for q in self.diss]
        for q in self.diss:
            nd = self.nodes[q]
            out += [
                float(self.structure[s].spin) if s in self.structure else 0.0
                for s in (nd.left, nd.right)
            ]
        return out

    def a_state(self, y) -> np.ndarray:
        """Semi-major axes of all orbits for the state y."""
        n6 = 6 * self.n_orb
        if self.gas is not None and len(y) > n6:
            return np.asarray(y[n6 : n6 + self.n_orb], float)
        a = self.a0.copy()
        for q, i in self.a_index.items():
            if len(y) > i:
                a[q] = y[i]
        return a

    def spins_state(self, y) -> dict[int, float]:
        """Spin rate [rad/yr] of every star: from the state for dissipative binaries."""
        out = {k: st.spin for k, st in self.structure.items()}
        for q, i in self.spin_index.items():
            if len(y) > i + 1:
                nd = self.nodes[q]
                out[nd.left], out[nd.right] = float(y[i]), float(y[i + 1])
        return out

    def set_masses(self, masses, y=None) -> None:
        """New stellar masses (slow isotropic wind): a_X M_X = const for every orbit X with
        constant eccentricity vector (Hadjidemetriou 1963)."""
        old = self.orbit_masses()
        self._set_mass_terms(masses)
        f = old / self.orbit_masses()
        self.a0 = self.a0 * f
        if y is not None:
            for q, i in self.a_index.items():
                y[i] *= f[q]
        self.prec, self.c2, self.c3 = self.coefficients(self.a0)
        self.L = list(self.prec[:, 0])

    def scale_orbits(self, y, factors: dict[str, float]) -> None:
        """Multiply the semi-major axes of the named orbits (mass transfer)."""
        if not factors:
            return
        for n, f in factors.items():
            q = self.names.index(n)
            self.a0[q] *= f
            if y is not None and q in self.a_index:
                y[self.a_index[q]] *= f
        self._kargs = None  # a0 changed in place
        self.prec, self.c2, self.c3 = self.coefficients(self.a0)
        self.L = list(self.prec[:, 0])

    def mass_state(self) -> dict:
        return {"m": self.masses.copy(), "a0": self.a0.copy()}

    def restore_mass_state(self, st: dict) -> None:
        self._set_mass_terms(st["m"])
        self.a0 = st["a0"].copy()
        self.prec, self.c2, self.c3 = self.coefficients(self.a0)
        self.L = list(self.prec[:, 0])

    def _two_body_arrays(self):
        """(orbit, star, companion, R, k2, spin) of every star of a two-star orbit with a
        structure (cached until the structure changes)."""
        c = self.__dict__.get("_tb_cache")
        if c is not None and c[0] is self.structure:
            return c[1]
        rows = []
        for k, node in enumerate(self.nodes):
            if isinstance(node.left, Node) or isinstance(node.right, Node):
                continue
            for i, j in ((node.left, node.right), (node.right, node.left)):
                st = self.structure.get(i)
                if st is not None:
                    rows.append((k, i, j, st.r_au, st.apsidal_k2, st.spin))
        out = None
        if rows:
            cols = list(zip(*rows, strict=True))
            out = (
                np.array(cols[0], dtype=np.int64),
                np.array(cols[1], dtype=np.int64),
                np.array(cols[2], dtype=np.int64),
                np.array(cols[3], float),
                np.array(cols[4], float),
                np.array(cols[5], float),
            )
        self._tb_cache = (self.structure, out)
        return out

    def set_structure(self, structure: dict[int, StarStructure] | None) -> None:
        """Replace the stellar structure (evolving stars, operator splitting)."""
        self.structure = structure or {}
        self.prec, self.c2, self.c3 = self.coefficients(self.a0)
        self.L = list(self.prec[:, 0])

    def coefficients(self, a: np.ndarray, spins: dict | None = None):
        """(prec, c2, c3) for the semi-major axes a (one per orbit) and the spin rates."""
        G = units.G
        m1, m2 = self.mk[:, 0], self.mk[:, 1]
        M = m1 + m2
        mu = m1 * m2 / M
        L = mu * np.sqrt(G * M * a)
        n = np.sqrt(G * M / a**3)
        alpha = L * 3.0 * (G * M) ** 1.5 / (units.C**2 * a**2.5) if self.gr else 0.0 * a
        beta = np.zeros(len(a))
        gamma = np.zeros(len(a))
        tb = self._two_body_arrays()
        if tb is not None:
            qk, si, sj, r_au, k2, w0 = tb
            w = w0 if spins is None else np.array([spins[s] for s in si])
            q = self.masses[sj] / self.masses[si]
            x5 = k2 * (r_au / a[qk]) ** 5
            np.add.at(beta, qk, 15.0 * q * x5)
            np.add.at(gamma, qk, (1.0 + q) * (w / n[qk]) ** 2 * x5)
        prec = np.column_stack([L, n, alpha, beta, gamma])
        ak, al = a[self.pk], a[self.pl]
        mu_k = mu[self.pk]
        c2 = G * mu_k * self._mo * ak**2 / (8.0 * al**3)
        c3 = 15.0 / 64.0 * G * mu_k * self._mo * ak**3 / al**4 * self._xs
        return prec, c2, c3

    def gas_factor(self, t: float) -> float:
        g = self.gas
        if t <= g.t_end_yr:
            return 1.0
        return float(np.exp(-(t - g.t_end_yr) / g.t_decay_yr))

    def tide_coefficients(self, a=None) -> np.ndarray:
        """c_X = pi G rho0 mu_X a_X^2 / L_X of every orbit (Galactic tide, swingersim.environment)."""
        a = self.a0 if a is None else a
        m1, m2 = self.mk[:, 0], self.mk[:, 1]
        M = m1 + m2
        mu = m1 * m2 / M
        L = mu * np.sqrt(units.G * M * a)
        return np.pi * units.G * self.tide[0] * mu * a**2 / L

    def _add_tide(self, y, out, a=None):
        from ..environment import tide_derivatives

        c = self.tide_coefficients(a)
        z = self.tide[1]
        for x in range(self.n_orb):
            sk, sp = slice(6 * x, 6 * x + 3), slice(6 * x + 3, 6 * x + 6)
            dK, dp = tide_derivatives(y[sk], y[sp], c[x], z)
            out[sk] += dK
            out[sp] += dp
        return out

    def _kernel_args(self):
        """Arguments of swingersim.secular.kernel.nested_dissipative_rhs for the masses,
        semi-major axes and structure in force (cached until one of them changes)."""
        key = (self.structure, self.masses, self.a0)
        c = self._kargs
        if c is not None and all(u is v for u, v in zip(c[0], key, strict=True)):
            return c[1]
        i64 = np.int64
        a_idx = np.full(self.n_orb, -1, dtype=i64)
        for q, i in self.a_index.items():
            a_idx[q] = i
        spin_idx = {}
        for q, i in self.spin_index.items():
            nd = self.nodes[q]
            spin_idx[nd.left], spin_idx[nd.right] = i, i + 1
        tb = self._two_body_arrays()
        tb_i, tb_f = np.zeros((0, 2), dtype=i64), np.zeros((0, 5))
        if tb is not None:
            qk, si, sj, r_au, k2, w0 = tb
            tb_i = np.column_stack([qk, [spin_idx.get(int(s), -1) for s in si]]).astype(i64)
            tb_f = np.column_stack([self.masses[si], self.masses[sj], r_au, k2, w0])
        d_i = np.zeros((len(self.diss), 6), dtype=i64)
        d_f = np.zeros((len(self.diss), 10))
        span = float(getattr(self, "backward_span", 0.0))
        slaved = set()
        for k, q in enumerate(self.diss):
            nd = self.nodes[q]
            pair = (nd.left, nd.right)
            d_i[k, :3] = q, self.a_index[q], self.spin_index[q]
            d_f[k, :2] = [self.masses[s] for s in pair]
            if all(s in self.structure for s in pair):  # tides need both stars' structure
                d_i[k, 3] = 1
                sts = [self.structure[s] for s in pair]
                d_f[k, 2:] = [
                    v for f in ("r_au", "apsidal_k2", "rg2", "lag_yr") for v in
                    (getattr(sts[0], f), getattr(sts[1], f))
                ]  # fmt: skip
                for j, (s, st) in enumerate(zip(pair, sts, strict=True)):
                    other = pair[1 - j]
                    ts = sync_time_yr(self.masses[s], self.masses[other], st.r_au,
                                      st.apsidal_k2, st.rg2, st.lag_yr, self.a0[q])  # fmt: skip
                    if slaved_spin(span, ts):
                        d_i[k, 4 + j] = 1
                        slaved.add(s)
        self.slaved = slaved
        args = (
            np.array([units.G, units.C, self.gr, self.octupole, self.gw], dtype=float),
            np.column_stack([self.mk, self.a0]).astype(float),
            a_idx,
            np.column_stack([self.pk, self.pl]).astype(i64),
            np.column_stack([self._mo, self._xs]).astype(float),
            np.ascontiguousarray(tb_i),
            np.ascontiguousarray(tb_f, dtype=float),
            d_i,
            d_f,
        )
        self._kargs = (key, args)
        return args

    def _spin_momenta(self) -> np.ndarray:
        """Spin angular momentum of the rotating stars of every two-star orbit (aligned spins,
        constant without dissipation; cached until the structure changes)."""
        c = self.__dict__.get("_sm_cache")
        if c is not None and c[0] is self.structure:
            return c[1]
        S = np.zeros(self.n_orb)
        for k, node in enumerate(self.nodes):
            if isinstance(node.left, Node) or isinstance(node.right, Node):
                continue
            for s in (node.left, node.right):
                st = self.structure.get(s)
                if st is not None:
                    S[k] += st.rg2 * self.masses[s] * st.r_au**2 * st.spin
        self._sm_cache = (self.structure, S)
        return S

    def _lock(self, y, out, L) -> np.ndarray:
        """The spins locked to their orbit (swingersim.secular.kernel._lock_spins): without
        dissipation, with the constant spins of the structure."""
        S = self._spin_momenta()
        if not S.any():
            return out
        out = np.array(out, dtype=float)
        for k in np.nonzero(S)[0]:
            sk = 6 * int(k)
            eta = float(np.linalg.norm(y[sk : sk + 3]))
            _lock_spins(
                np.ascontiguousarray(y, dtype=float), out, sk, float(L[k]) * eta, float(S[k])
            )
        return out

    def _rhs_dissipative(self, t, y):
        """Compiled right-hand side (swingersim.secular.kernel.nested_dissipative_rhs) plus the
        Galactic tide."""
        out = nested_dissipative_rhs(np.ascontiguousarray(y, dtype=float), *self._kernel_args())
        if self.tide is not None:
            self._add_tide(y, out, self.a_state(y))
        return out

    def rhs(self, t, y):
        if self.diss:
            return self._rhs_dissipative(t, y)
        if self.gas is None and self.tide is not None:
            out = secular_rhs_pairs(
                np.ascontiguousarray(y),
                self.pk,
                self.pl,
                self.c2,
                self.c3,
                self.prec,
                self.octupole,
            )[1].copy()
            return self._add_tide(y, self._lock(y, out, self.L))
        if self.gas is None:
            out = secular_rhs_pairs(
                np.ascontiguousarray(y),
                self.pk,
                self.pl,
                self.c2,
                self.c3,
                self.prec,
                self.octupole,
            )[1]
            return self._lock(y, out, self.L)
        n6 = 6 * self.n_orb
        a = y[n6:]
        prec, c2, c3 = self.coefficients(a)
        out = np.empty_like(y)
        out[:n6] = secular_rhs_pairs(
            np.ascontiguousarray(y[:n6]), self.pk, self.pl, c2, c3, prec, self.octupole
        )[1]
        out[:n6] = self._lock(y[:n6], out[:n6], prec[:, 0])
        f = self.gas_factor(t)
        out[n6:] = -f * a / self.tau[:, 0]
        nrm = self.disc_normal
        for x in range(self.n_orb):
            sk, sp = slice(6 * x, 6 * x + 3), slice(6 * x + 3, 6 * x + 6)
            K, p = y[sk], y[sp]
            if np.isfinite(self.tau[x, 1]):
                eta2 = K @ K
                out[sp] += -f * p / self.tau[x, 1]
                out[sk] += f * (p @ p) / eta2 / self.tau[x, 1] * K
            if np.isfinite(self.tau[x, 2]):
                k = K / np.sqrt(K @ K)
                w = f / self.tau[x, 2] * np.cross(k, nrm)
                out[sk] += np.cross(w, K)
                out[sp] += np.cross(w, p)
        if self.tide is not None:
            self._add_tide(y, out, a)
        return out

    def energy(self, y) -> float:
        h = float(
            secular_rhs_pairs(
                np.ascontiguousarray(y[: 6 * self.n_orb]),
                self.pk,
                self.pl,
                self.c2,
                self.c3,
                self.prec,
                self.octupole,
            )[0]
        )
        if self.tide is not None:
            from ..environment import tide_energy

            c = self.tide_coefficients() * np.array(self.L)  # energies, not rates
            for x in range(self.n_orb):
                K, p = y[6 * x : 6 * x + 3], y[6 * x + 3 : 6 * x + 6]
                h += tide_energy(K, p, c[x], self.tide[1])
        return h

    def angular_momentum(self, y) -> np.ndarray:
        """Orbits, plus the aligned spins of the dissipative binaries."""
        if not self.diss:  # the constant spins of the structure, locked to their orbit
            tot = sum(L * y[6 * k : 6 * k + 3] for k, L in enumerate(self.L))
            for k, Sk in enumerate(self._spin_momenta()):
                if Sk:
                    K = y[6 * k : 6 * k + 3]
                    tot = tot + Sk * K / np.linalg.norm(K)
            return tot
        a = self.a_state(y)
        prec, _, _ = self.coefficients(a, self.spins_state(y))
        tot = sum(prec[k, 0] * y[6 * k : 6 * k + 3] for k in range(self.n_orb))
        spins = self.spins_state(y)
        for q in self.diss:
            nd = self.nodes[q]
            K = y[6 * q : 6 * q + 3]
            khat = K / np.linalg.norm(K)
            for s in (nd.left, nd.right):
                st = self.structure.get(s)
                if st is not None:
                    tot = tot + st.rg2 * self.masses[s] * st.r_au**2 * spins[s] * khat
        return tot
