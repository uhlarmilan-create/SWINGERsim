"""Extremes between outputs: maximum eccentricity and closest periastron of every orbit.

Short eccentricity peaks (e.g. during eccentric ZLK cycles) can be much shorter than the
output spacing. The N-body backend therefore samples the state at a finer internal cadence
(a few orbits of the shortest inner binary) and keeps, for every orbit of the hierarchy, the
largest eccentricity and the smallest periastron distance q = a(1 - e) with their times.
Sampling a few inner orbits apart is enough because these peaks are secular: they last many
inner orbits. Collisions are detected separately at every integrator step by REBOUND.

For the binaries made of two stars the tracker also follows the Roche-lobe filling factor
R / R_L at periastron, with R_L from Eggleton (1983) and the periastron distance a(1 - e) of
the osculating orbit (the usual approximation for eccentric orbits, cf. Sepinsky et al. 2007).
"""

from __future__ import annotations

import numpy as np

from ..hierarchy import Node, mass_of
from ..kepler import cross3
from ..stellar import roche_lobe_eggleton


def _weights(node: Node, masses: np.ndarray) -> np.ndarray:
    """Row vector w with r_rel = w @ pos (position of the right COM minus the left COM)."""
    w = np.zeros(len(masses))
    for child, sign in ((node.left, -1.0), (node.right, 1.0)):
        stars = child.stars() if isinstance(child, Node) else [child]
        m = masses[stars]
        w[stars] += sign * m / m.sum()
    return w


THRESHOLDS_1ME = (1e-1, 1e-2, 1e-3, 1e-4)  # first times with 1 - e below these are recorded


def threshold_key(thr: float) -> str:
    return f"t_1me_lt_{thr:.0e}_yr".replace("e-0", "e-")


class ExtremeTracker:
    def __init__(
        self, root: Node, masses: np.ndarray, G: float, radii_au: np.ndarray | None = None
    ):
        self.nodes = root.nodes()
        self.names = [n.name for n in self.nodes]
        self.G = G
        self._set_masses(masses)
        k = len(self.nodes)
        self.e_max = np.full(k, -np.inf)
        self.t_e_max = np.full(k, np.nan)
        self.q_min = np.full(k, np.inf)
        self.t_q_min = np.full(k, np.nan)
        self.t_first = np.full((k, len(THRESHOLDS_1ME)), np.nan)
        self.e_last = np.full(k, np.nan)
        # Roche lobes: (node index, star index, R / (R_L / separation)) for star-star binaries
        self.roche: list[tuple[int, int, float]] = []
        self.radii_au = None if radii_au is None else np.asarray(radii_au, float).copy()
        self._set_roche(masses, self.radii_au)
        self.fill_max = {s: 0.0 for _, s, _ in self.roche}
        self.t_fill_max = {s: float("nan") for _, s, _ in self.roche}

    def _set_masses(self, masses) -> None:
        masses = np.asarray(masses, float)
        self.W = np.array([_weights(n, masses) for n in self.nodes])
        self.mu = np.array(
            [self.G * (mass_of(n.left, masses) + mass_of(n.right, masses)) for n in self.nodes]
        )

    def _set_roche(self, masses, radii_au) -> None:
        self.roche = []
        if radii_au is None:
            return
        for k, n in enumerate(self.nodes):
            for s, o in ((n.left, n.right), (n.right, n.left)):
                if isinstance(s, Node):
                    continue  # lobes of single stars: against a star or a whole subsystem
                scale = radii_au[s] / roche_lobe_eggleton(masses[s] / mass_of(o, masses))
                self.roche.append((k, s, scale))

    def set_stars(self, masses, radii_au=None) -> None:
        """New masses and radii (evolving stars, mass loss); the extremes are kept. Without
        ``radii_au`` the radii stay, but the Roche lobes still follow the new mass ratios."""
        self._set_masses(masses)
        if radii_au is not None:
            self.radii_au = np.asarray(radii_au, float).copy()
        if self.radii_au is not None:
            self._set_roche(np.asarray(masses, float), self.radii_au)

    def update(self, t: float, pos: np.ndarray, vel: np.ndarray) -> list[int]:
        """Record a sample; returns the stars that fill their Roche lobe at periastron."""
        r = self.W @ pos
        v = self.W @ vel
        rn = np.sqrt(np.einsum("ij,ij->i", r, r))
        v2 = np.einsum("ij,ij->i", v, v)
        energy = 0.5 * v2 - self.mu / rn
        h = cross3(r, v)
        evec = cross3(v, h) / self.mu[:, None] - r / rn[:, None]
        e = np.sqrt(np.einsum("ij,ij->i", evec, evec))
        bound = energy < 0
        self.e_last = np.where(bound, e, np.nan)  # osculating e of every node (this sample)
        a = np.where(bound, -self.mu / (2.0 * np.where(bound, energy, -1.0)), np.nan)
        q = np.where(bound, a * (1.0 - e), rn)  # unbound: use the current separation
        better_e = bound & (e > self.e_max)
        self.e_max[better_e] = e[better_e]
        self.t_e_max[better_e] = t
        for j, thr in enumerate(THRESHOLDS_1ME):
            first = bound & (1.0 - e < thr) & np.isnan(self.t_first[:, j])
            self.t_first[first, j] = t
        better_q = q < self.q_min
        self.q_min[better_q] = q[better_q]
        self.t_q_min[better_q] = t
        overflow = []
        for k, s, scale in self.roche:
            fill = scale / q[k]
            if fill > self.fill_max[s]:
                self.fill_max[s], self.t_fill_max[s] = float(fill), t
            if fill >= 1.0:
                overflow.append(s)
        return overflow

    def roche_result(self) -> dict[str, dict[str, float]]:
        """Largest Roche-lobe filling factor at periastron for every star of a two-star binary
        and every star orbiting a subsystem."""
        return {
            f"m{s + 1}": {"fill_max": self.fill_max[s], "t_fill_max_yr": self.t_fill_max[s]}
            for s in self.fill_max
        }

    def result(self) -> dict[str, dict[str, float]]:
        return {
            name: {
                "e_max": float(self.e_max[i]) if np.isfinite(self.e_max[i]) else float("nan"),
                "t_e_max_yr": float(self.t_e_max[i]),
                "q_min_au": float(self.q_min[i]),
                "t_q_min_yr": float(self.t_q_min[i]),
                **{
                    threshold_key(thr): float(self.t_first[i, j])
                    for j, thr in enumerate(THRESHOLDS_1ME)
                },
            }
            for i, name in enumerate(self.names)
        }
