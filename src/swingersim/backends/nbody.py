"""N-body backend: REBOUND (IAS15) + REBOUNDx gr_full (1PN for all bodies) + tides_spin.

Stellar structure (physics.tides) uses the REBOUNDx effect ``tides_spin`` (Lu et al. 2023),
the equilibrium-tide model of Eggleton, Kiseleva & Hut (1998): every star with a radius,
Love number and spin is distorted by its own rotation and by the tides of all other stars,
and with a non-zero time lag the tides dissipate energy and the spins evolve.
"""

from __future__ import annotations

import time
import warnings
from dataclasses import dataclass, field
from itertools import combinations

import numpy as np
import rebound
import reboundx

from .. import massloss, masstransfer, units
from .. import masstransfer as _mt
from ..config import Config
from ..diagnostics import hierarchy as hdet
from ..diagnostics.binary import binary_diagnostics, eclipse_summary
from ..diagnostics.extremes import ExtremeTracker, _weights
from ..diagnostics.zlk import PairDiagnostics, analyse, system_diagnostics
from ..hierarchy import Node, build_tree, initialize, mass_of, node_elements, node_h_unit
from ..kepler import cross3, states_to_elements
from ..stellar import StarProperties, star_properties, zams_radius_tout
from . import evolving as ev

# mass loss uses the REBOUNDx operator modify_mass with IAS15; its accuracy is verified by the
# ml-* benchmarks (swingersim validate), so the generic REBOUNDx caution is not repeated per run
warnings.filterwarnings("ignore", message="REBOUNDx: Operators that affect particle")

ELEMENT_KEYS = ("a", "e", "inc", "Omega", "omega", "M")


MEGNO_SEED = 1  # seed of the random initial tangent vector of REBOUND's MEGNO: reproducible values


def lyapunov_time_yr(sim: rebound.Simulation) -> float:
    """Lyapunov time [yr] (e-folding time of the separation of neighbouring trajectories) from
    the MEGNO of a REBOUND simulation (time unit yr).

    REBOUND's ``sim.lyapunov()`` is the least-squares slope of the mean MEGNO <Y>(t) against t
    (Cincotta & Simo 2000, eq. 24). With REBOUND's definition Y(t) = (2/t) int s (d|delta|/ds)/|delta|
    ds (the factor 2 is explicit in WHFast and in the Gauss-Radau weights of IAS15, which sum to 2),
    a chaotic orbit with Lyapunov exponent lambda has Y -> lambda t and <Y> -> lambda t / 2, so the
    slope is lambda / 2 and the Lyapunov time is 1 / lambda = 1 / (2 sim.lyapunov())."""
    slope = sim.lyapunov()
    return float(1.0 / (2.0 * slope)) if slope > 0 else float("inf")


def default_radius_rsun(m: float) -> float:
    """Zero-age main-sequence radius [Rsun] of the default stellar model (Tout et al. 1996)."""
    return zams_radius_tout(m)


# REBOUNDx tides_spin: its constant time lag tau gives tidal rates twice those of the
# constant-time-lag model of Mignard (1979) / Hut (1981) / Leconte et al. (2010) with the same
# lag (found from the source code and checked by the benchmark tides-leconte-*), so the lag
# of the configuration is halved before it is handed to REBOUNDx.
REBX_TAU_PER_LAG = 0.5

# Intermediate results of a long run (swingersim run -o X writes X.partial.json while it runs,
# for swingersim gui): the file, the seconds between writes, the most points per series.
PARTIAL = {"path": None, "every_s": 20.0, "max_points": 2000}


class _PartialWriter:
    """Writes the elements, mutual inclinations and energy error computed so far, every
    PARTIAL['every_s'] seconds, as JSON (atomically: a reader never sees half a file)."""

    def __init__(self, t_end: float):
        self.path = PARTIAL["path"]
        self.t_end = float(t_end)
        self.last = time.perf_counter()

    def __call__(self, j, times, orbits, mutual, d_energy):
        if self.path is None:
            return
        now = time.perf_counter()
        if now - self.last < PARTIAL["every_s"]:
            return
        self.last = now
        import json
        from pathlib import Path

        n = j + 1
        step = max(1, n // PARTIAL["max_points"])
        sl = slice(0, n, step)

        def clean(a):
            return [None if not np.isfinite(v) else float(v) for v in np.asarray(a)[sl]]

        data = {
            "t_end_yr": self.t_end,
            "time_yr": clean(times),
            "orbits": {o: {k: clean(v) for k, v in el.items()} for o, el in orbits.items()},
            "mutual_inclination": {k: clean(v) for k, v in mutual.items()},
            "energy_error": clean(d_energy),
            "written": time.strftime("%Y-%m-%dT%H:%M:%S"),
        }
        path = Path(self.path)
        tmp = path.with_name(path.name + ".tmp")
        try:  # intermediate results are a convenience: never let them stop the run
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp.write_text(json.dumps(data))
            tmp.replace(path)
        except OSError:
            pass


def parent_of(root: Node, star: int) -> Node:
    """The binary in which ``star`` is a direct component."""
    for n in root.nodes():
        if n.left == star or n.right == star:
            return n
    raise ValueError(f"star {star} not in the tree")


def spin_vector(cfg: Config, name: str, node: Node, pos, vel, masses) -> np.ndarray:
    """Initial spin angular velocity [rad/yr] of a star; zero without physics.rotation."""
    if not cfg.physics.rotation:
        return np.zeros(3)
    spec = cfg.stars.get(name)
    p_rot = spec.P_rot if spec else "sync"
    if p_rot in ("sync", "pseudosync"):
        m_tot = mass_of(node.left, masses) + mass_of(node.right, masses)
        omega = 2.0 * np.pi / units.sma_to_period(node.elements.a, m_tot)
        if p_rot == "pseudosync":
            from ..tides import pseudo_sync_factor

            omega *= pseudo_sync_factor(float(node.elements.e))
    else:
        omega = 2.0 * np.pi / units.days_to_years(p_rot)
    h = node_h_unit(node, pos, vel, masses)
    nvec = np.cross([0.0, 0.0, 1.0], h)
    nvec = nvec / np.linalg.norm(nvec) if np.linalg.norm(nvec) > 1e-12 else np.array([1.0, 0, 0])
    ob = np.radians(spec.obliquity if spec else 0.0)
    az = np.radians(spec.spin_azimuth if spec else 0.0)
    axis = np.cos(ob) * h + np.sin(ob) * (np.cos(az) * nvec + np.sin(az) * np.cross(h, nvec))
    return omega * axis


@dataclass
class RunResult:
    times: np.ndarray  # [yr]
    orbits: dict[str, dict[str, np.ndarray]]  # node name -> element -> series
    mutual: dict[str, np.ndarray]  # "A-B" -> mutual inclination [rad]
    energy_error: np.ndarray  # relative energy error
    event: str | None = (
        None  # None, collision, energy_error, hierarchy_change, roche_overflow, wall_time, e_reached
    )
    event_time: float | None = None
    event_detail: dict = field(default_factory=dict)
    hierarchy: list[str] = field(default_factory=list)  # detected hierarchy at each output
    configured_hierarchy: str = ""  # hierarchy given in the configuration
    transitions: list[dict] = field(default_factory=list)  # confirmed hierarchy changes
    diagnostics: list[PairDiagnostics] = field(default_factory=list)
    extremes: dict[str, dict[str, float]] = field(default_factory=dict)  # per orbit
    megno: np.ndarray | None = None  # MEGNO at each output (Newtonian), if requested
    lyapunov_time_yr: float | None = None
    system: dict = field(default_factory=dict)  # 2+2: resonances between the binaries
    stars: dict[str, StarProperties] = field(default_factory=dict)  # stellar model used
    spins: dict[str, dict[str, np.ndarray]] = field(default_factory=dict)  # P_rot, obliquity
    roche: dict[str, dict[str, float]] = field(default_factory=dict)  # max Roche filling
    angmom_error: np.ndarray | None = None  # relative error of total angular momentum
    dissipative: bool = False  # tidal dissipation on: energy is not conserved
    gas_phase: bool = False  # physics.gas: energy and angular momentum are not conserved
    mass_loss: bool = False  # physics.mass_loss: energy and angular momentum are not conserved
    backend: str = "nbody"
    binaries: dict[str, dict] = field(default_factory=dict)  # apsidal motion, tidal times
    # evolving stars (physics.evolution): star -> {radius_rsun, luminosity_lsun, teff_k,
    # phase} at every output
    stellar: dict[str, dict[str, np.ndarray]] = field(default_factory=dict)
    wall_s: float | None = None  # wall-clock time of the run [s]
    # N-body: state of the stars at the end (pos, vel [AU, AU/yr], masses, radii [AU], spin
    # vectors [rad/yr] or None) - the start of a continuation after a merger (swingersim.merger)
    final_state: dict | None = None


def model_props(cfg: Config, root, masses, pos, vel, m_dyn=None) -> list[StarProperties]:
    """Stellar properties of all stars, with tidal lags from the structure when
    physics.tidal_lag is not 'fixed' (Hurley et al. 2002). ``masses`` label the stellar
    models; ``m_dyn`` are the masses of the dynamics (differ with mass loss from the tracks)."""
    props = [star_properties(cfg, f"m{k + 1}", m) for k, m in enumerate(masses)]
    masses = masses if m_dyn is None else m_dyn
    if cfg.physics.tides and cfg.physics.tidal_lag != "fixed":
        from ..hierarchy import node_elements
        from ..tides import with_structure_lags

        pairs = [
            (n.left, n.right, node_elements(n, pos, vel, masses).a)
            for n in root.nodes()
            if isinstance(n.left, int) and isinstance(n.right, int)
        ]
        spins = [
            float(
                np.linalg.norm(spin_vector(cfg, f"m{k + 1}", parent_of(root, k), pos, vel, masses))
            )
            for k in range(len(masses))
        ]
        props = with_structure_lags(cfg, props, pairs, spins, masses)
    return props


def build_simulation(
    cfg: Config, gr: bool | None = None, structure: bool | None = None
) -> tuple[rebound.Simulation, object, dict, Node, np.ndarray]:
    """REBOUND simulation of the configuration.

    Returns (sim, rebx, forces, root, masses); ``forces`` maps "gr" / "tides" to the loaded
    REBOUNDx effects. ``gr`` and ``structure`` override physics.gr and physics.tides (the
    MEGNO shadow simulation is purely Newtonian). build_setup also returns the configured
    masses (labels of the stellar models) and the stellar properties.
    """
    return build_setup(cfg, gr, structure)[:5]


def build_setup(
    cfg: Config, gr: bool | None = None, structure: bool | None = None
) -> tuple[rebound.Simulation, object, dict, Node, np.ndarray, np.ndarray, list[StarProperties]]:
    """build_simulation with (sim, rebx, forces, root, masses, m_lab, props): ``m_lab`` the
    configured masses (labels of the tracks), ``props`` the stellar models at t = 0."""
    root, m_lab = build_tree(cfg)
    masses = massloss.initial_masses(cfg, m_lab)  # dynamics at t = 0 (current masses)
    if not np.array_equal(masses, m_lab):  # periods with the masses of the dynamics
        root = build_tree(cfg, kepler_masses=masses)[0]
    pos, vel = initialize(root, masses)
    props = model_props(cfg, root, m_lab, pos, vel, masses)

    sim = rebound.Simulation()
    sim.G = units.G
    sim.integrator = cfg.run.integrator
    for m, x, v, pr in zip(masses, pos, vel, props):
        r = pr.radius_rsun * units.RSUN_AU
        sim.add(m=m, x=x[0], y=x[1], z=x[2], vx=v[0], vy=v[1], vz=v[2], r=r)
    sim.move_to_com()
    if cfg.run.collisions:
        sim.collision = "direct"  # REBOUND raises rebound.Collision (default resolve = halt)

    use_gr = cfg.physics.gr if gr is None else gr
    use_tides = cfg.physics.tides if structure is None else structure
    rebx, forces = None, {}
    if use_gr or use_tides:
        rebx = reboundx.Extras(sim)
    if use_gr:
        forces["gr"] = rebx.load_force("gr_full")
        rebx.add_force(forces["gr"])
        forces["gr"].params["c"] = units.C
    if use_tides:
        with warnings.catch_warnings():  # notice about the v4.5 fix of the EKH typo
            warnings.simplefilter("ignore", RuntimeWarning)
            forces["tides"] = rebx.load_force("tides_spin")
        rebx.add_force(forces["tides"])
        from ..tides import k2_rotating

        for k, (p, pr) in enumerate(zip(sim.particles, props)):
            w = spin_vector(cfg, f"m{k + 1}", parent_of(root, k), pos, vel, masses)
            k2 = k2_rotating(cfg, k, pr.apsidal_k2, float(np.linalg.norm(w)), p.m, pr.radius_rsun)
            p.params["k2"] = 2.0 * k2
            p.params["I"] = pr.rg2 * p.m * p.r**2
            p.params["tau"] = REBX_TAU_PER_LAG * pr.lag_s / units.SECONDS_PER_YEAR
            p.params["Omega"] = rebound.Vec3d(*w)
        rebx.initialize_spin_ode(forces["tides"])
    if massloss.enabled(cfg) and (gr is None and structure is None):
        if rebx is None:
            rebx = reboundx.Extras(sim)
        with warnings.catch_warnings():  # IAS15 + operator: accuracy checked by benchmarks
            warnings.simplefilter("ignore", RuntimeWarning)
            forces["mass_loss"] = rebx.load_operator("modify_mass")
            rebx.add_operator(forces["mass_loss"])
    return sim, rebx, forces, root, masses, np.asarray(m_lab, float), props


def _state(sim: rebound.Simulation, n_star: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Positions and velocities of the stars (the first n_star particles; passing stars of
    physics.encounters are appended after them)."""
    n = len(sim.particles)
    pos, vel = np.zeros((n, 3)), np.zeros((n, 3))
    sim.serialize_particle_data(xyz=pos, vxvyvz=vel)
    if n_star is not None and n_star < n:
        return pos[:n_star].copy(), vel[:n_star].copy()
    return pos, vel


def _node_geometry(nodes: list[Node], masses: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(W, mu): relative positions of all nodes r = W @ pos (right minus left centre of mass)
    and mu = G (M_left + M_right)."""
    W = np.array([_weights(n, masses) for n in nodes])
    mu = np.array([units.G * (mass_of(n.left, masses) + mass_of(n.right, masses)) for n in nodes])
    return W, mu


def _e_of(tracker: ExtremeTracker, node: int) -> float:
    """Osculating eccentricity of a node at the tracker's last sample (0 if unbound/unknown)."""
    e = float(tracker.e_last[node])
    return e if np.isfinite(e) else 0.0


def _invariable_basis(pos, vel, masses):
    """Unit vectors (x, y, z) of the invariable plane; z along the total angular momentum."""
    L = (masses[:, None] * np.cross(pos, vel)).sum(axis=0)
    z = L / np.linalg.norm(L)
    x = np.cross([0.0, 0.0, 1.0], z)
    if np.linalg.norm(x) < 1e-12:
        x = np.array([1.0, 0.0, 0.0])
    x /= np.linalg.norm(x)
    return x, np.cross(z, x), z


class _EccentricityReached(Exception):
    pass


class _WallTime(Exception):
    """run.max_wall_s reached inside an integration step (heartbeat, see run_nbody)."""


class _RocheOverflow(Exception):
    def __init__(self, stars: list[int]):
        super().__init__(stars)
        self.stars = stars


def _touching(sim: rebound.Simulation) -> list[str]:
    """Names of stars whose surfaces overlap (closest pair if none overlaps exactly)."""
    pos, _ = _state(sim)
    radii = np.array([p.r for p in sim.particles])
    best, hit = np.inf, []
    for i, j in combinations(range(len(radii)), 2):
        gap = np.linalg.norm(pos[i] - pos[j]) - radii[i] - radii[j]
        if gap < best:
            best, hit = gap, [i, j]
    return [f"m{k + 1}" for k in hit]


def _pair_frame(sim, i: int, j: int):
    """(m_i, m_j, centre-of-mass position, its velocity, relative position r_j - r_i, relative
    velocity) of stars i and j."""
    pi, pj = sim.particles[i], sim.particles[j]
    mi, mj = pi.m, pj.m
    M = mi + mj
    ri, rj = np.array([pi.x, pi.y, pi.z]), np.array([pj.x, pj.y, pj.z])
    vi, vj = np.array([pi.vx, pi.vy, pi.vz]), np.array([pj.vx, pj.vy, pj.vz])
    return mi, mj, (mi * ri + mj * rj) / M, (mi * vi + mj * vj) / M, rj - ri, vj - vi


def _set_pair(sim, i: int, j: int, rc, vc, r, v) -> None:
    """Place stars i and j (their current masses) at relative position r and velocity v about
    the centre of mass rc moving with vc."""
    pi, pj = sim.particles[i], sim.particles[j]
    mi, mj = pi.m, pj.m
    M = mi + mj
    pi.x, pi.y, pi.z = rc - mj / M * r
    pj.x, pj.y, pj.z = rc + mi / M * r
    pi.vx, pi.vy, pi.vz = vc - mj / M * v
    pj.vx, pj.vy, pj.vz = vc + mi / M * v


def _scale_pair(sim, i: int, j: int, g: float) -> None:
    """Scale the relative orbit of stars i and j by g in semi-major axis at fixed eccentricity
    vector and phase (r -> g r, v -> v / sqrt(g)), keeping their centre of mass and its
    velocity (mass transfer, operator splitting)."""
    _, _, rc, vc, r, v = _pair_frame(sim, i, j)
    _set_pair(sim, i, j, rc, vc, r * g, v / np.sqrt(g))


GW_MAX_SHRINK = 0.1  # largest relative shrinking of a by GW emission in one output step
GW_SUBSTEP = 0.01  # relative change of a (and absolute of e) per Runge-Kutta sub-step


def peters_evolve(m1: float, m2: float, a: float, e: float, dt: float) -> tuple[float, float]:
    """a and e after dt [yr] of gravitational-wave emission: Peters' (1964) orbit-averaged
    equations integrated with classical Runge-Kutta sub-steps, each changing a by at most
    GW_SUBSTEP a and e by GW_SUBSTEP (dt < 0: backward in time). Returns (0, e) when the
    binary merges within dt."""
    from ..gw import peters_rates

    def rates(a_, e_):
        if a_ <= 0.0:
            return 0.0, 0.0
        return peters_rates(m1, m2, a_, min(max(e_, 0.0), 0.999999))

    sgn = 1.0 if dt >= 0 else -1.0
    left = abs(dt)
    a0 = a
    for _ in range(100_000):
        if left <= 0.0:
            break
        da, de = rates(a, e)
        if da == 0.0 and de == 0.0:
            break
        h = min(left, GW_SUBSTEP * a / max(abs(da), 1e-300), GW_SUBSTEP / max(abs(de), 1e-300))
        hs = sgn * h
        k1 = (da, de)
        k2 = rates(a + 0.5 * hs * k1[0], e + 0.5 * hs * k1[1])
        k3 = rates(a + 0.5 * hs * k2[0], e + 0.5 * hs * k2[1])
        k4 = rates(a + hs * k3[0], e + hs * k3[1])
        a = a + hs / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0])
        e = e + hs / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1])
        left -= h
        if a <= 1e-3 * a0:  # merged (Peters' equations end at a = 0)
            return 0.0, max(e, 0.0)
    return a, min(max(e, 0.0), 0.999999)


def _gw_pair(sim, i: int, j: int, dt: float) -> dict | None:
    """Gravitational-wave emission of the pair i, j over dt [yr] (orbit-averaged, Peters 1964):
    new a and e at the same orientation and mean anomaly, centre of mass kept. If a would
    shrink by more than GW_MAX_SHRINK in this step (or the pair merge) the step does not
    resolve the inspiral: nothing is applied and the details are returned (event)."""
    from ..kepler import Elements, elements_to_state, state_to_elements

    mi, mj, rc, vc, r, v = _pair_frame(sim, i, j)
    mu = units.G * (mi + mj)
    try:
        el = state_to_elements(mu, r, v)
    except ValueError:
        return None
    a_new, e_new = peters_evolve(mi, mj, el.a, el.e, dt)
    if a_new < (1.0 - GW_MAX_SHRINK) * el.a:
        return {
            "stars": [f"m{i + 1}", f"m{j + 1}"],
            "a_au": float(el.a),
            "a_after_step_au": float(a_new),
            "note": "GW inspiral faster than the output step: more outputs (run.n_output) needed",
        }
    r, v = elements_to_state(mu, Elements(a_new, e_new, el.inc, el.Omega, el.omega, el.M))
    _set_pair(sim, i, j, rc, vc, r, v)
    return None


def _common_envelope_jump(sim, i: int, j: int, mi: float, mj: float, a_f: float) -> None:
    """End of a common envelope (physics.mass_transfer.common_envelope): the envelope leaves at
    once, so the centre of mass of stars i and j keeps its position and velocity (the outer
    orbits feel an impulsive mass loss, Hadjidemetriou 1963), the stars take their new masses
    (core and companion) and a circular relative orbit of the separation a_f in the same plane
    and at the same phase."""
    _, _, rc, vc, r, v = _pair_frame(sim, i, j)
    rhat = r / np.linalg.norm(r)
    n = np.cross(r, v)
    n = n / np.linalg.norm(n)
    M = mi + mj
    r_new = a_f * rhat
    v_new = np.sqrt(sim.G * M / a_f) * np.cross(n, rhat)
    sim.particles[i].m, sim.particles[j].m = mi, mj
    _set_pair(sim, i, j, rc, vc, r_new, v_new)


def run_nbody(cfg: Config) -> RunResult:
    """N-body run of the configuration; with physics.mergers a merger of two stars is followed
    by a run of the merged system (swingersim.merger)."""
    if cfg.physics.mergers is not None or _mt.dynamic_onset(cfg) is not None:
        from ..merger import run_with_mergers

        return run_with_mergers(cfg, _run_nbody_once)
    return _run_nbody_once(cfg)


def _run_nbody_once(cfg: Config) -> RunResult:
    t_start = time.perf_counter()
    cfg = massloss.resolve_masses(cfg)
    if cfg.run.t_max_yr == 0:
        raise ValueError("run.t_max_yr must not be 0 (negative values integrate backward)")
    # m_lab: configured masses (labels of the tracks); props: stellar models at t = 0
    sim, rebx, forces, root, masses, m_lab, props = build_setup(cfg)
    tides = "tides" in forces
    evolving = bool(cfg.physics.evolution)
    mass_loss = massloss.enabled(cfg)
    varying = evolving or mass_loss
    dissipative = tides and any(pr.lag_s > 0 for pr in props)
    n_star = len(masses)
    env = None
    if cfg.physics.galactic_tide is not None or cfg.physics.encounters is not None:
        from ..environment import NBodyEnvironment

        top = root
        m_sys = float(np.sum(masses))
        env = NBodyEnvironment(
            cfg,
            n_star,
            m_sys,
            float(top.elements.a),
            units.sma_to_period(top.elements.a, m_sys),
            float(cfg.run.t_max_yr),
        )

    def energy() -> float:
        # gr_full_hamiltonian is the whole 1PN Hamiltonian (kinetic + potential + 1PN terms,
        # REBOUNDx gr_full.c), not a correction to add to sim.energy()
        e = rebx.gr_full_hamiltonian(forces["gr"]) if "gr" in forces else sim.energy()
        if tides:
            e += rebx.tides_spin_energy()
        if env is not None:
            e += env.tide_energy(sim)
        return e

    def angmom() -> np.ndarray:
        L = np.array(sim.angular_momentum())
        if tides:
            L = L + np.array(rebx.spin_angular_momentum())
        return L

    star_parent = [parent_of(root, k) for k in range(len(masses))]
    mt_hist = masstransfer.histories(cfg)
    mt_stars = masstransfer.transferring_stars(cfg) if mt_hist else set()
    # common envelopes: discrete jumps at the end of the envelope ejection
    ce_jumps = []
    for name, hist in mt_hist.items():
        if hist.event == "post_common_envelope" and hist.ce is not None:
            k = int(np.searchsorted(hist.t, hist.ce["t_yr"])) + 1
            m_post = hist.m[:, k]
            ce_jumps.append((float(hist.t[k]), name, hist, float(m_post[0]), float(m_post[1])))

    nodes = root.nodes()
    pairs = list(combinations([n.name for n in nodes], 2))
    node_idx = {n.name: i for i, n in enumerate(nodes)}
    pair_idx = np.array([(node_idx[a], node_idx[b]) for a, b in pairs], int).reshape(-1, 2)
    star_parent_idx = [node_idx[star_parent[k].name] for k in range(len(masses))]
    geo = None  # (weights W of the relative positions, G M) of the nodes for these masses
    t_end, stop_evolution = ev.time_limit(cfg, m_lab, cfg.run.t_max_yr)
    times = np.linspace(0.0, t_end, cfg.run.n_output)
    keys = (*ELEMENT_KEYS, "inc_inv", "Omega_inv")
    orbits = {n.name: {k: np.full(len(times), np.nan) for k in keys} for n in nodes}
    pos0, vel0 = _state(sim)
    inv_x, inv_y, inv_z = _invariable_basis(pos0, vel0, masses)

    radii_au = np.array([p.r for p in sim.particles])
    tracker = ExtremeTracker(root, masses, units.G, radii_au)
    p_min = min(
        units.sma_to_period(n.elements.a, mass_of(n.left, masses) + mass_of(n.right, masses))
        for n in nodes
    )
    sample_dt = cfg.run.extremes_every_orbits * p_min

    megno_sim = None
    megno = np.full(len(times), np.nan) if cfg.run.megno else None
    if cfg.run.megno:
        megno_sim, _, _, _, _ = build_simulation(cfg, gr=False, structure=False)
        megno_sim.collision = "none"
        megno_sim.init_megno(seed=MEGNO_SEED)
    mutual = {f"{a}-{b}": np.full(len(times), np.nan) for a, b in pairs}
    d_energy = np.full(len(times), np.nan)
    d_angmom = np.full(len(times), np.nan)
    spins = (
        {
            f"m{k + 1}": {
                "P_rot_d": np.full(len(times), np.nan),
                "obliquity_deg": np.full(len(times), np.nan),
                # unit spin vector in the frame of the orbital elements (z: the line of sight
                # of the inclinations): sky-projected spin-orbit angles (Rossiter-McLaughlin)
                "spin_x": np.full(len(times), np.nan),
                "spin_y": np.full(len(times), np.nan),
                "spin_z": np.full(len(times), np.nan),
            }
            for k in range(len(masses))
        }
        if tides
        else {}
    )
    e0 = energy()
    L0 = angmom()
    e_ref, L_ref = e0, L0
    # the stellar models and the masses of the moments of inertia of the last stellar step
    cur = {"props": props, "m_mid": np.array([p.m for p in sim.particles[:n_star]])}

    def step_stars(t_a: float, t_b: float) -> None:
        """Evolving stars / mass loss over the output step t_a -> t_b: structure at the middle
        of the step, masses reaching their values at t_b (REBOUNDx modify_mass, dm/dt =
        m/tau_mass with tau_mass fitted to the step)."""
        nonlocal e_ref, L_ref
        m_old = np.array([p.m for p in sim.particles[:n_star]])
        jumping = {name: tx for tx, name, *_ in ce_jumps if min(t_a, t_b) < tx <= max(t_a, t_b)}
        if mt_hist:  # mass transfer: the change of a beyond a M = const, applied at once
            for name, hist in mt_hist.items():
                if name in jumping:  # the common envelope is applied as a jump
                    g = masstransfer.a_ratio(cfg, name, t_a, hist.ce["t_yr"])
                else:
                    g = masstransfer.a_ratio(cfg, name, t_a, t_b)
                if g != 1.0:
                    _scale_pair(sim, *hist.pair, g)
        if mass_loss:
            target = massloss.masses_at(cfg, m_lab, t_b)
            for name in jumping:  # masses jump at the end of the envelope
                pre = massloss.masses_at(cfg, m_lab, mt_hist[name].ce["t_yr"])
                for k in mt_hist[name].pair:
                    target[k] = pre[k]
            for p, mt in zip(sim.particles[:n_star], target):
                lr = np.log(mt / p.m)
                p.params["tau_mass"] = (t_b - sim.t) / lr if abs(lr) > 1e-15 else np.inf
        if evolving:
            new = ev.props_at(cfg, m_lab, 0.5 * (t_a + t_b))
            m_mid = massloss.masses_at(cfg, m_lab, 0.5 * (t_a + t_b)) if mass_loss else m_old
            if tides and cfg.physics.tidal_lag != "fixed":
                # lags from the current structure, separation and spins (Hurley et al. 2002),
                # as model_props at t = 0 and the secular backends at every segment
                from ..tides import with_structure_lags

                ps_ = sim.particles
                pos_ = np.array([[ps_[k].x, ps_[k].y, ps_[k].z] for k in range(n_star)])
                vel_ = np.array([[ps_[k].vx, ps_[k].vy, ps_[k].vz] for k in range(n_star)])
                m_now = np.array([ps_[k].m for k in range(n_star)])
                pairs_ = [
                    (n.left, n.right, node_elements(n, pos_, vel_, m_now).a)
                    for n in root.nodes()
                    if isinstance(n.left, int) and isinstance(n.right, int)
                ]
                spins_ = [
                    float(np.linalg.norm(np.array(ps_[k].params["Omega"]))) for k in range(n_star)
                ]
                new = with_structure_lags(cfg, new, pairs_, spins_, m_now)
            for k, (p, pn, po) in enumerate(zip(sim.particles[:n_star], new, cur["props"])):
                p.r = pn.radius_rsun * units.RSUN_AU
                if tides:
                    w = np.array(p.params["Omega"])
                    # angular momentum I w with I = rg2 m R^2: from the mass the moment of
                    # inertia was set with in the last step to the mass of this step
                    w0 = float(np.linalg.norm(w))
                    if w0 > 0:  # the accreted angular momentum and the break-up limit
                        w = w * ev.scale_spin(w0, po, pn, cur["m_mid"][k], m_mid[k]) / w0
                    wn = float(np.linalg.norm(w))
                    if wn > 0:  # magnetic braking over this output step
                        w = w * ev.brake_spin(cfg, k, wn, pn, m_mid[k], abs(t_b - t_a)) / wn
                    p.params["Omega"] = rebound.Vec3d(*w)
                    from ..tides import k2_rotating

                    k2 = k2_rotating(
                        cfg, k, pn.apsidal_k2, float(np.linalg.norm(w)), m_mid[k], pn.radius_rsun
                    )
                    p.params["k2"] = 2.0 * k2
                    p.params["I"] = pn.rg2 * m_mid[k] * p.r**2
                    p.params["tau"] = REBX_TAU_PER_LAG * pn.lag_s / units.SECONDS_PER_YEAR
            cur["props"] = new
            cur["m_mid"] = np.asarray(m_mid, float).copy()
            tracker.set_stars(m_old, np.array([p.r for p in sim.particles[:n_star]]))
        e_ref, L_ref = energy(), angmom()

    event, event_time, detail = None, None, {}

    configured = (hdet.from_node(root),)
    confirmed = configured
    candidate, candidate_since, candidate_count = None, None, 0
    hier_strings: list[str] = []
    transitions: list[dict] = []

    s = 1.0 if cfg.run.t_max_yr > 0 else -1.0  # direction of time (backward runs: -1)
    if s < 0:
        sim.dt = -abs(sim.dt)
        if megno_sim is not None:
            megno_sim.dt = -abs(megno_sim.dt)
        if cfg.run.collisions:
            # REBOUND's "direct" search counts only pairs that approach each other in forward
            # time (dv.dr < 0); integrating backward, two stars plunging into each other are
            # skipped, the step size collapses at the point-mass singularity and the run hangs.
            # The "line" search tests the segment of the whole step in either direction.
            sim.collision = "line"
    next_sample = s * sample_dt if sample_dt > 0 else s * np.inf
    wall0 = time.perf_counter()
    wall_hit = [False]
    if cfg.run.max_wall_s:
        # the limit is checked after every step (REBOUND heartbeat, ~0.5 us per step): a step
        # size that collapses inside one sim.integrate() call would never return to the checks
        # between output steps
        deadline = wall0 + cfg.run.max_wall_s

        def _heartbeat(sim_pointer) -> None:
            if time.perf_counter() > deadline:
                wall_hit[0] = True
                sim_pointer.contents.stop()

        # set the C field directly: some REBOUND releases cannot store the callback through
        # the sim.heartbeat setter ('_hb' missing from Simulation.__slots__); the reference
        # kept in this frame keeps the callback alive for the whole run
        heartbeat_type = dict(rebound.Simulation._fields_)["_heartbeat"]
        heartbeat_cb = heartbeat_type(_heartbeat)
        sim._heartbeat = heartbeat_cb

    from ..pdot import PeriodDrift

    drift = PeriodDrift.of(cfg, root, masses)  # orbits.<pair>.Pdot

    def integrate(t_to: float, exact: int) -> None:
        if drift is not None:
            drift.integrate(sim, t_to, exact)
        else:
            sim.integrate(t_to, exact_finish_time=exact)
        if wall_hit[0]:
            raise _WallTime

    e_targets = [(tracker.names.index(n), e) for n, e in cfg.run.stop_when_e_reached.items()]
    gw_pairs = (
        [(n.left, n.right) for n in nodes if isinstance(n.left, int) and isinstance(n.right, int)]
        if cfg.physics.gw
        else []
    )

    def excluded(k: int) -> bool:
        """Roche overflow of star k expected by the mass-transfer history (not an event), with
        the osculating eccentricity of its binary at the last sample."""
        e_k = _e_of(tracker, star_parent_idx[k])
        return masstransfer.roche_excluded(cfg, k, sim.t, e=e_k)

    partial = _PartialWriter(t_end)
    for j, t in enumerate(times):
        if varying and j > 0:
            step_stars(times[j - 1], t)
        if gw_pairs and j > 0:  # orbit-averaged emission over the step (operator splitting)
            for i_, j_ in gw_pairs:
                gw_detail = _gw_pair(sim, i_, j_, t - times[j - 1])
                if gw_detail is not None:
                    # the step does not resolve the inspiral: stop at the last output
                    event, event_time, detail = "gw_unresolved", sim.t, gw_detail
                    break
            if event:
                j -= 1
                break
            e_ref, L_ref = energy(), angmom()
        try:
            if j > 0:
                for tx, name, hist, mi, mj in ce_jumps:
                    if min(times[j - 1], t) < tx <= max(times[j - 1], t) and not event:
                        integrate(tx, 1)
                        i_, j_ = hist.pair
                        _common_envelope_jump(sim, i_, j_, mi, mj, float(hist.ce["a_f_au"]))
                        for k in (i_, j_):
                            sim.particles[k].params["tau_mass"] = np.inf
                        new = ev.props_at(cfg, m_lab, tx)
                        for k in (i_, j_):
                            sim.particles[k].r = new[k].radius_rsun * units.RSUN_AU
                        tracker.set_stars(
                            np.array([p.m for p in sim.particles[:n_star]]),
                            np.array([p.r for p in sim.particles[:n_star]]),
                        )
                        e_ref, L_ref = energy(), angmom()
            chunks = env.boundaries(sim.t, t) if env is not None and env.enabled else [t]
            for tc in chunks:
                while s * next_sample < s * tc:
                    # a sample within a few steps of tc goes exactly to tc: integrating with
                    # exact_finish_time=0 would overshoot tc and REBOUND would then integrate
                    # backward to it
                    near = s * next_sample + 4.0 * abs(sim.dt) >= s * tc
                    integrate(tc if near else next_sample, 1 if near else 0)
                    over = [
                        k
                        for k in tracker.update(sim.t, *_state(sim, n_star))
                        if not (k in mt_stars and excluded(k))
                    ]
                    if over and _mt.stop_on_overflow(cfg):
                        raise _RocheOverflow(over)
                    if e_targets and all(tracker.e_max[k] >= e for k, e in e_targets):
                        raise _EccentricityReached
                    next_sample = s * max(s * next_sample + sample_dt, s * sim.t + 0.5 * sample_dt)
                integrate(tc, 1)
                if env is not None and env.at(sim, tc):
                    e_ref, L_ref = energy(), angmom()
        except rebound.Collision:
            event, event_time = "collision", sim.t
            detail = {"stars": _touching(sim)}
        except _RocheOverflow as ov:
            event, event_time = "roche_overflow", sim.t
            detail = {"stars": [f"m{k + 1}" for k in ov.stars]}
        except _EccentricityReached:
            event, event_time = "e_reached", sim.t
        except _WallTime:
            event, event_time = "wall_time", sim.t
        pos, vel = _state(sim, n_star)
        if mass_loss or (drift is not None and drift.changes_masses):
            masses = np.array([p.m for p in sim.particles[:n_star]])
            tracker.set_stars(masses)
        over = [k for k in tracker.update(sim.t, pos, vel) if not (k in mt_stars and excluded(k))]
        if over and _mt.stop_on_overflow(cfg) and not event:
            event, event_time = "roche_overflow", sim.t
            detail = {"stars": [f"m{k + 1}" for k in over]}
        if megno_sim is not None:
            megno_sim.integrate(sim.t, exact_finish_time=1)
            megno[j] = megno_sim.megno()
        # relative orbits of all nodes at once (r = W @ pos, as ExtremeTracker.update): the
        # elements, orbit normals, invariable-plane angles, mutual inclinations and spin
        # obliquities all come from these (unbound nodes: NaN elements)
        if mass_loss or geo is None or (drift is not None and drift.changes_masses):
            geo = _node_geometry(nodes, masses)
        r_rel, v_rel = geo[0] @ pos, geo[0] @ vel
        els = states_to_elements(geo[1], r_rel, v_rel)
        hv = cross3(r_rel, v_rel)
        hu = hv / np.sqrt(np.einsum("ij,ij->i", hv, hv))[:, None]
        inc_inv = np.arccos(np.clip(hu @ inv_z, -1.0, 1.0))
        om_inv = np.mod(np.arctan2(hu @ inv_x, -(hu @ inv_y)), 2.0 * np.pi)
        for i, n in enumerate(nodes):
            o = orbits[n.name]
            for k in ELEMENT_KEYS:
                o[k][j] = els[k][i]
            o["inc_inv"][j] = inc_inv[i]
            o["Omega_inv"][j] = om_inv[i]
        if len(pair_idx):
            hl, hr = hu[pair_idx[:, 0]], hu[pair_idx[:, 1]]
            imut = np.arccos(np.clip(np.einsum("ij,ij->i", hl, hr), -1.0, 1.0))
            for (pa, pb), c in zip(pairs, imut):
                mutual[f"{pa}-{pb}"][j] = c
        # evolving stars: error within the output step (structure fixed in it); mass loss:
        # energy and angular momentum change physically (the wind), reported, not checked
        d_energy[j] = abs(energy() / e_ref - 1.0)
        d_angmom[j] = np.linalg.norm(angmom() - L_ref) / np.linalg.norm(L_ref)
        for k, p in enumerate(sim.particles[:n_star] if tides else []):
            w = np.array(p.params["Omega"])
            wn = np.linalg.norm(w)
            spins[f"m{k + 1}"]["P_rot_d"][j] = (
                2.0 * np.pi / wn * units.DAYS_PER_YEAR if wn > 0 else np.inf
            )
            if wn > 0:
                for c, wc in zip(("spin_x", "spin_y", "spin_z"), w / wn):
                    spins[f"m{k + 1}"][c][j] = wc
                hk = hu[star_parent_idx[k]]
                spins[f"m{k + 1}"]["obliquity_deg"][j] = np.degrees(
                    np.arccos(np.clip(w @ hk / wn, -1.0, 1.0))
                )

        found = hdet.detect(pos, vel, masses, units.G)
        hier_strings.append(hdet.to_string(found))
        if found == confirmed:
            candidate, candidate_count = None, 0
        elif found == candidate:
            candidate_count += 1
        else:
            candidate, candidate_since, candidate_count = found, sim.t, 1
        if candidate is not None and candidate_count >= cfg.run.hierarchy_persistence:
            transitions.append(
                {
                    "time_yr": float(candidate_since),
                    "kind": hdet.classify_change(confirmed, candidate),
                    "from": hdet.to_string(confirmed),
                    "to": hdet.to_string(candidate),
                    "label": hdet.label(candidate),
                }
            )
            confirmed, candidate, candidate_count = candidate, None, 0
            if cfg.run.stop_on_hierarchy_change and not event:
                event, event_time = "hierarchy_change", transitions[-1]["time_yr"]
                detail = {k: v for k, v in transitions[-1].items() if k != "time_yr"}

        partial(j, times, orbits, mutual, d_energy)
        if event:
            break
        conserved_error = d_angmom[j] if dissipative else d_energy[j]
        # mass loss and orbits.<pair>.Pdot change energy and angular momentum physically
        if not mass_loss and drift is None and conserved_error > cfg.run.max_energy_error:
            event, event_time = "energy_error", sim.t
            break
        if not event and e_targets and all(tracker.e_max[k] >= e for k, e in e_targets):
            event, event_time = "e_reached", sim.t
            break
        if cfg.run.max_wall_s and time.perf_counter() - wall0 > cfg.run.max_wall_s:
            event, event_time = "wall_time", sim.t
            break

    n = j + 1
    if event is not None:
        # the last row holds the state at the moment the run stopped (collision, Roche
        # overflow, e reached, wall time ...), not at the nominal output time
        times[j] = sim.t
    res = RunResult(
        times=times[:n],
        orbits={k: {e: s[:n] for e, s in v.items()} for k, v in orbits.items()},
        mutual={k: s[:n] for k, s in mutual.items()},
        energy_error=d_energy[:n],
        event=event,
        event_time=event_time,
        event_detail=detail,
        hierarchy=hier_strings,
        configured_hierarchy=hdet.to_string(configured),
        transitions=transitions,
        extremes=tracker.result(),
        stars={f"m{k + 1}": pr for k, pr in enumerate(props)},
        spins={k: {e: v[:n] for e, v in d.items()} for k, d in spins.items()},
        roche=tracker.roche_result(),
        angmom_error=d_angmom[:n],
        dissipative=dissipative,
        mass_loss=mass_loss,
    )
    if varying:
        res.stellar = ev.stellar_series(cfg, m_lab, res.times) if evolving else {}
        if mass_loss:
            ms = massloss.mass_series(cfg, m_lab, res.times)
            for k in range(len(masses)):
                res.stellar.setdefault(f"m{k + 1}", {})["mass_msun"] = ms[k]
        if stop_evolution is not None and not event and n == len(times):
            res.event, res.event_time = stop_evolution, float(times[-1])
            res.event_detail = {"age_myr": cfg.stellar_model.age_myr + times[-1] / 1e6}
    masses = massloss.initial_masses(cfg, m_lab)
    res.final_state = {
        "pos": pos,
        "vel": vel,
        "masses": np.array([p.m for p in sim.particles[:n_star]]),
        "radii": np.array([p.r for p in sim.particles[:n_star]]),
        "spins": (
            np.array([np.array(p.params["Omega"]) for p in sim.particles[:n_star]])
            if tides
            else None
        ),
    }
    if megno is not None:
        res.megno = megno[:n]
        res.lyapunov_time_yr = lyapunov_time_yr(megno_sim)
    res.wall_s = time.perf_counter() - t_start  # whole run (diagnostics below are fast)
    res.diagnostics = analyse(root, masses, cfg.physics.gr, res, (pos0, vel0), radii_au)
    res.system = system_diagnostics(root, masses, res.diagnostics, (pos0, vel0))
    spins0 = [
        spin_vector(cfg, f"m{k + 1}", star_parent[k], pos0, vel0, masses) if tides else np.zeros(3)
        for k in range(len(masses))
    ]
    res.binaries = binary_diagnostics(
        root, masses, props, spins0, pos0, vel0, res, cfg.physics.gr, tides, cfg=cfg
    )
    res.system.update(eclipse_summary(res, root, props))
    res.system.update(masstransfer.system_entries(cfg))
    res.system.update(masstransfer.eccentric_transfer(cfg, res))
    if env is not None and cfg.physics.encounters is not None:
        res.system["encounters_n"] = env.counts["n"]
        res.system["encounters_close"] = 0  # every passage is integrated as a particle
    return res
