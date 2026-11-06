"""Few-body scattering: a binary (or a single star) meets another binary (or single star).

Used to study how a quadruple can form dynamically: binary-binary encounters that end as a
bound 2+2 or 3+1 system, or as a triple plus an escaping star ("zachytenie dvojhviezdy
dvojhviezdou"). Newtonian point masses (optionally with 1PN relativity) integrated with
REBOUND IAS15; collisions are detected when stellar radii are given.

Set-up of one encounter (as in Hut & Bahcall 1983 and FEWBODY, Fregeau et al. 2004):

* each binary has a random orientation (isotropic), a random mean anomaly and a fixed or
  thermal (f(e) = 2e) eccentricity;
* the centres of mass approach each other on a two-body hyperbola with velocity at infinity
  v_inf and impact parameter b; b is fixed or drawn uniformly in area up to b_max
  (r_peri_max gives b_max through gravitational focusing, b^2 = q^2 + 2 G M q / v_inf^2);
* the integration starts where the separation is r_start_factor times the largest
  semi-major axis (or just outside the pericentre of the hyperbola if that is larger).

The encounter ends when the top-level objects found by the hierarchy detector
(swingersim.diagnostics.hierarchy) are pairwise unbound, receding and farther apart than the
starting distance, or when the system stays bound until t_max (a long-lived resonance or a
captured quadruple / triple; its dynamical stability is then judged with the Mardling &
Aarseth 2001 criterion applied to every nested triple).
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Literal

import numpy as np
import rebound
import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator

from . import units
from .diagnostics import hierarchy as hdet
from .kepler import Elements, elements_to_state, state_to_elements
from .stability import ma01_ratio
from .stellar import roche_lobe_eggleton

KMS_TO_AUYR = units.KMS_TO_AUYR


# --------------------------------------------------------------------------- configuration


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid")


class ScatterBody(_Strict):
    """A single star (one mass) or a binary (two masses, a_au, e)."""

    masses: list[float] = Field(min_length=1, max_length=2)
    a_au: float | None = Field(default=None, gt=0)
    e: float = Field(default=0.0, ge=0, lt=1)
    e_thermal: bool = False  # draw e from f(e) = 2e (redrawn if a star overfills its Roche lobe)
    radii_rsun: list[float] | None = None  # enables collision detection
    # polytropic index of every star for tidal capture (1.5, 2 or 3; default 1.5 below
    # 0.4 Msun, 3 above)
    polytrope_n: list[Literal[1.5, 2.0, 3.0]] | None = None

    @model_validator(mode="after")
    def _check(self):
        if len(self.masses) == 2 and self.a_au is None:
            raise ValueError("a binary needs a_au")
        if self.radii_rsun is not None and len(self.radii_rsun) != len(self.masses):
            raise ValueError("radii_rsun needs one radius per mass")
        if self.polytrope_n is not None and len(self.polytrope_n) != len(self.masses):
            raise ValueError("polytrope_n needs one index per mass")
        if len(self.masses) == 2 and self.radii_rsun:
            # e = 0 must be possible: otherwise e_thermal would redraw forever (and any e
            # means a binary that is already interacting)
            m1, m2 = self.masses
            r1, r2 = (r * units.RSUN_AU for r in self.radii_rsun)
            room = self.a_au * min(
                roche_lobe_eggleton(m1 / m2) / r1, roche_lobe_eggleton(m2 / m1) / r2
            )
            if room < 1.0:
                raise ValueError(
                    f"a star of the binary (a = {self.a_au} AU) overfills its Roche lobe even on a "
                    "circular orbit: increase a_au or reduce radii_rsun"
                )
        return self


class ScatterSpec(_Strict):
    target: ScatterBody
    projectile: ScatterBody
    v_inf_kms: float | None = Field(default=None, ge=0)
    v_inf_vc: float | None = Field(default=None, ge=0)  # in units of the critical velocity
    b_au: float | None = Field(default=None, ge=0)
    b_max_au: float | None = Field(default=None, gt=0)
    r_peri_max_au: float | None = Field(default=None, gt=0)
    b_distribution: Literal["area", "linear"] = "area"  # uniform in b^2 (default) or in b
    r_start_factor: float = Field(default=50.0, gt=1)
    t_max_yr: float | None = Field(default=None, gt=0)
    max_wall_s: float = Field(default=600.0, gt=0)
    check_periods: float = Field(default=5.0, gt=0)  # hierarchy checks every this many periods
    gr: bool = False
    tidal_capture: bool = False  # energy lost to tides at close periastron passages (PZM93)
    n: int = Field(default=100, ge=1)
    seed: int = 1

    @model_validator(mode="after")
    def _check(self):
        if (self.v_inf_kms is None) == (self.v_inf_vc is None):
            raise ValueError("give exactly one of v_inf_kms and v_inf_vc")
        if sum(x is not None for x in (self.b_au, self.b_max_au, self.r_peri_max_au)) != 1:
            raise ValueError("give exactly one of b_au, b_max_au and r_peri_max_au")
        if (self.v_inf_kms == 0 or self.v_inf_vc == 0) and self.r_peri_max_au is None:
            # a parabolic encounter has no impact parameter (b = q (1 + 2 G M / (q v^2))^1/2
            # diverges): with b the angular momentum b v_inf would be 0, a head-on collision
            raise ValueError(
                "v_inf = 0 (parabolic encounters) needs r_peri_max_au, not b_au / b_max_au"
            )
        if len(self.target.masses) + len(self.projectile.masses) < 2:
            raise ValueError("need at least two stars")
        if self.tidal_capture and not (self.target.radii_rsun and self.projectile.radii_rsun):
            raise ValueError("tidal_capture needs radii_rsun of all stars")
        return self


class ScatterConfig(_Strict):
    name: str = "scatter"
    description: str = ""
    scatter: ScatterSpec


def load_scatter_config(path) -> ScatterConfig:
    return ScatterConfig.model_validate(yaml.safe_load(Path(path).read_text()))


# --------------------------------------------------------------------------- set-up


def _binding(body: ScatterBody) -> float:
    if len(body.masses) == 1:
        return 0.0
    m1, m2 = body.masses
    return units.G * m1 * m2 / (2.0 * body.a_au)


def critical_velocity_kms(spec: ScatterSpec) -> float:
    """v_c: total energy zero (Hut & Bahcall 1983): mu_inf v_c^2 / 2 = sum of binding energies."""
    Mt, Mp = sum(spec.target.masses), sum(spec.projectile.masses)
    mu = Mt * Mp / (Mt + Mp)
    eb = _binding(spec.target) + _binding(spec.projectile)
    return math.sqrt(2.0 * eb / mu) / KMS_TO_AUYR


def v_inf_auyr(spec: ScatterSpec) -> float:
    if spec.v_inf_kms is not None:
        return spec.v_inf_kms * KMS_TO_AUYR
    return spec.v_inf_vc * critical_velocity_kms(spec) * KMS_TO_AUYR


def b_max(spec: ScatterSpec) -> float | None:
    """Largest impact parameter [AU] (None when v_inf = 0 and only r_peri_max is given)."""
    if spec.b_au is not None:
        return spec.b_au
    if spec.b_max_au is not None:
        return spec.b_max_au
    v = v_inf_auyr(spec)
    if v == 0:
        return None
    mu = units.G * (sum(spec.target.masses) + sum(spec.projectile.masses))
    q = spec.r_peri_max_au
    return math.sqrt(q * q + 2.0 * mu * q / v**2)


def _random_binary(body: ScatterBody, rng) -> dict | None:
    if len(body.masses) == 1:
        return None
    e = body.e
    if body.e_thermal:
        e_max = 1.0
        if body.radii_rsun:
            # no star may overfill its Roche lobe at pericentre (Eggleton 1983): f(e) = 2e
            # truncated at the e where a lobe touches its star (ScatterBody checked e = 0 fits)
            m1, m2 = body.masses
            r1, r2 = (r * units.RSUN_AU for r in body.radii_rsun)
            room = body.a_au * min(
                roche_lobe_eggleton(m1 / m2) / r1, roche_lobe_eggleton(m2 / m1) / r2
            )
            e_max = 1.0 - 1.0 / room
        e = e_max * math.sqrt(rng.uniform())
    return {
        "a": body.a_au,
        "e": e,
        "inc": math.acos(rng.uniform(-1, 1)),
        "Omega": rng.uniform(0, 2 * math.pi),
        "omega": rng.uniform(0, 2 * math.pi),
        "M": rng.uniform(0, 2 * math.pi),
    }


def draw_encounters(spec: ScatterSpec) -> list[dict]:
    """Parameters of all n encounters (reproducible from the seed)."""
    rng = np.random.default_rng(spec.seed)
    out = []
    for k in range(spec.n):
        u = rng.uniform()
        bm = b_max(spec)
        if spec.b_au is not None:
            b, q = spec.b_au, None
        elif bm is not None:
            b, q = bm * (math.sqrt(u) if spec.b_distribution == "area" else u), None
        else:  # parabolic approach: uniform in pericentre = uniform in area (focusing limit)
            b, q = None, spec.r_peri_max_au * u
        out.append(
            {
                "id": k,
                "b_au": b,
                "q_au": q,
                "target": _random_binary(spec.target, rng),
                "projectile": _random_binary(spec.projectile, rng),
            }
        )
    return out


def _hyperbola_state(mu, v_inf, h, r0):
    """Relative position and velocity at separation r0 before pericentre, orbit in the x-y
    plane, pericentre on the +x axis."""
    if h == 0.0:  # head-on
        return np.array([r0, 0.0, 0.0]), np.array([-math.sqrt(v_inf**2 + 2 * mu / r0), 0, 0])
    p = h * h / mu
    e = math.sqrt(1.0 + v_inf**2 * h * h / mu**2)
    cf = (p / r0 - 1.0) / e
    f = -math.acos(max(-1.0, min(1.0, cf)))
    r = r0 * np.array([math.cos(f), math.sin(f), 0.0])
    v = mu / h * np.array([-math.sin(f), e + math.cos(f), 0.0])
    return r, v


def _peri(mu, v_inf, h):
    if h == 0.0:
        return 0.0
    e = math.sqrt(1.0 + v_inf**2 * h * h / mu**2)
    return h * h / mu / (1.0 + e)


def initial_state(spec: ScatterSpec, enc: dict):
    """Masses, radii [AU], positions and velocities (centre-of-mass frame) and the start
    distance of one encounter."""
    bodies = (spec.target, spec.projectile)
    masses = np.array(spec.target.masses + spec.projectile.masses, float)
    radii = []
    for body in bodies:
        radii += (
            [r * units.RSUN_AU for r in body.radii_rsun]
            if body.radii_rsun
            else [0.0] * len(body.masses)
        )
    Mt, Mp = sum(spec.target.masses), sum(spec.projectile.masses)
    mu = units.G * (Mt + Mp)
    v = v_inf_auyr(spec)
    if enc["b_au"] is not None:
        h = enc["b_au"] * v
    else:
        q = enc["q_au"]
        h = math.sqrt(q * q * v * v + 2.0 * mu * q)
    a_max = max([b.a_au or 0.0 for b in bodies])
    r0 = spec.r_start_factor * a_max if a_max > 0 else 0.0
    r0 = max(r0, 1.5 * _peri(mu, v, h), 1e-3)
    rel_r, rel_v = _hyperbola_state(mu, v, h, r0)
    pos, vel = [], []
    for body, key, sgn, M in (
        (spec.target, "target", -1, Mt),
        (spec.projectile, "projectile", 1, Mp),
    ):
        frac = (Mp if key == "target" else Mt) / (Mt + Mp)
        cr, cv = sgn * frac * rel_r, sgn * frac * rel_v
        if len(body.masses) == 1:
            pos.append(cr)
            vel.append(cv)
            continue
        m1, m2 = body.masses
        el = Elements(**enc[key])
        r, w = elements_to_state(units.G * (m1 + m2), el)
        pos += [cr - m2 / (m1 + m2) * r, cr + m1 / (m1 + m2) * r]
        vel += [cv - m2 / (m1 + m2) * w, cv + m1 / (m1 + m2) * w]
    return masses, np.array(radii), np.array(pos), np.array(vel), r0, h


# --------------------------------------------------------------------------- analysis


def _com(tree, masses, pos, vel):
    stars = _star_list(tree)
    m = masses[stars].sum()
    return (
        m,
        (masses[stars, None] * pos[stars]).sum(0) / m,
        (masses[stars, None] * vel[stars]).sum(0) / m,
    )


def _star_list(tree) -> list[int]:
    return [tree] if isinstance(tree, int) else _star_list(tree[0]) + _star_list(tree[1])


def _orbit(tree, masses, pos, vel):
    """Elements and angular-momentum unit vector of the relative orbit of a tree node."""
    ma, ra, va = _com(tree[0], masses, pos, vel)
    mb, rb, vb = _com(tree[1], masses, pos, vel)
    r, v = rb - ra, vb - va
    el = state_to_elements(units.G * (ma + mb), r, v)
    hv = np.cross(r, v)
    return el, hv / np.linalg.norm(hv)


def stability_ratios(tree, masses, pos, vel) -> list[float]:
    """MA01 ratio (> 1 = stable) for every nested triple in a bound tree."""
    out = []
    if isinstance(tree, int):
        return out
    el_out, h_out = _orbit(tree, masses, pos, vel)
    for k, child in enumerate(tree):
        if isinstance(child, int):
            continue
        other = tree[1 - k]
        el_in, h_in = _orbit(child, masses, pos, vel)
        q_out = masses[_star_list(other)].sum() / masses[_star_list(child)].sum()
        imut = math.acos(max(-1.0, min(1.0, float(h_in @ h_out))))
        out.append(ma01_ratio(el_in.a, el_out.a, el_out.e, q_out, imut))
        out += stability_ratios(child, masses, pos, vel)
    return out


def _binaries(trees) -> set[frozenset]:
    out = set()

    def walk(t):
        if isinstance(t, int):
            return
        if isinstance(t[0], int) and isinstance(t[1], int):
            out.add(frozenset(t))
        walk(t[0])
        walk(t[1])

    for t in trees:
        walk(t)
    return out


def _separated(trees, masses, pos, vel, r_min) -> bool:
    objs = [_com(t, masses, pos, vel) for t in trees]
    for i in range(len(objs)):
        for j in range(i + 1, len(objs)):
            dr, dv = objs[j][1] - objs[i][1], objs[j][2] - objs[i][2]
            if np.linalg.norm(dr) < r_min or dr @ dv <= 0.0:
                return False
    return True


def _bound_total(masses, pos, vel) -> bool:
    """Is the whole system bound (negative total energy)?"""
    kin = 0.5 * float((masses * (vel**2).sum(1)).sum())
    pot = 0.0
    for i in range(len(masses)):
        for j in range(i + 1, len(masses)):
            pot -= units.G * masses[i] * masses[j] / np.linalg.norm(pos[i] - pos[j])
    return kin + pot < 0.0


def _check_interval(dt0, trees, masses, pos, vel) -> float:
    """Time to the next check of the hierarchy: a tenth of the outer period of a fully bound
    system, or a twentieth of the shortest separation / relative speed of unbound objects,
    never less than dt0."""
    if len(trees) == 1:
        if isinstance(trees[0], int):
            return dt0
        el, _ = _orbit(trees[0], masses, pos, vel)
        return max(dt0, 0.1 * units.sma_to_period(el.a, float(masses.sum())))
    objs = [_com(t, masses, pos, vel) for t in trees]
    times = []
    for i in range(len(objs)):
        for j in range(i + 1, len(objs)):
            dr, dv = objs[j][1] - objs[i][1], objs[j][2] - objs[i][2]
            times.append(np.linalg.norm(dr) / max(np.linalg.norm(dv), 1e-30))
    return max(dt0, 0.05 * min(times))


def classify(trees, initial_trees, collided) -> str:
    if collided:
        return "collision"
    lab = hdet.label(trees)
    if len(trees) == 1:
        return "bound " + lab
    if lab == hdet.label(initial_trees):
        return "flyby" if _binaries(trees) == _binaries(initial_trees) else "exchange"
    if all(isinstance(t, int) for t in trees):
        return "full ionization"
    if lab == "binary + single + single":  # outcome classes of Antognini & Thompson (2016)
        return (
            "single ionization"
            if _binaries(trees) <= _binaries(initial_trees)
            else ("exchange + ionization")
        )
    return lab


def run_encounter(spec: ScatterSpec, enc: dict) -> dict:
    masses, radii, pos, vel, r0, _h = initial_state(spec, enc)
    n_t = len(spec.target.masses)
    initial = (
        (0, 1) if n_t == 2 else 0,
        (n_t, n_t + 1) if len(spec.projectile.masses) == 2 else n_t,
    )
    sim = rebound.Simulation()
    sim.G = units.G
    sim.integrator = "ias15"
    for m, r, p, v in zip(masses, radii, pos, vel):
        sim.add(m=m, r=r, x=p[0], y=p[1], z=p[2], vx=v[0], vy=v[1], vz=v[2])
    collisions = radii.sum() > 0
    if collisions:
        sim.collision = "direct"
        sim.collision_resolve = "halt"
    rebx = None
    if spec.gr:
        import reboundx

        rebx = reboundx.Extras(sim)
        gr = rebx.load_force("gr_full")
        gr.params["c"] = units.C
        rebx.add_force(gr)
    tides = None
    if spec.tidal_capture:
        from .tidalcapture import PeriastronDissipation, default_polytrope

        poly = []
        for body in (spec.target, spec.projectile):
            for k, m in enumerate(body.masses):
                poly.append(body.polytrope_n[k] if body.polytrope_n else default_polytrope(m))
        tides = PeriastronDissipation(sim, masses, radii, poly)
        sim.heartbeat = tides
    E0 = sim.energy()
    E_total = E0
    L_total = float(np.linalg.norm(sum(m * np.cross(p, v) for m, p, v in zip(masses, pos, vel))))
    # energy scale: the total energy can be ~0 (v_inf = v_c), so errors are measured
    # relative to the kinetic energy plus the absolute potential energy at the start
    kin = 0.5 * float((masses * (vel**2).sum(1)).sum())
    E_scale = kin + abs(E0 - kin)
    v = v_inf_auyr(spec)
    Mt = sum(spec.target.masses) + sum(spec.projectile.masses)
    periods = [
        units.sma_to_period(b.a_au, sum(b.masses))
        for b in (spec.target, spec.projectile)
        if len(b.masses) == 2
    ]
    t_cross = r0 / max(v, math.sqrt(units.G * Mt / r0))  # time to cover r0
    dt0 = min(periods) * spec.check_periods if periods else t_cross / 20.0
    dt0 = min(dt0, t_cross / 5.0)
    dt = dt0  # interval between checks of the hierarchy (IAS15 chooses its own steps)
    t_max = spec.t_max_yr or max(20.0 * t_cross, 1e5 * (max(periods) if periods else 0.0))
    wall0 = time.time()
    collided, pair = False, None
    stable_since = None
    scrambles, in_scramble = 0, False  # phases with no bound pair (Stone & Leigh 2019)
    trees, status = initial, "unresolved"
    while sim.t < t_max:
        try:
            sim.integrate(min(sim.t + dt, t_max))
        except rebound.Collision:
            collided = True
            P = np.array([[p.x, p.y, p.z] for p in sim.particles])
            d = np.linalg.norm(P[:, None] - P[None], axis=2) - (radii[:, None] + radii[None])
            np.fill_diagonal(d, np.inf)
            pair = list(np.unravel_index(int(np.argmin(d)), d.shape))
            break
        P = np.array([[p.x, p.y, p.z] for p in sim.particles])
        V = np.array([[p.vx, p.vy, p.vz] for p in sim.particles])
        trees = hdet.detect(P, V, masses, units.G)
        all_single = all(isinstance(t, int) for t in trees)
        if all_single and not in_scramble and _bound_total(masses, P, V):
            scrambles += 1
        in_scramble = all_single
        dt = _check_interval(dt0, trees, masses, P, V)
        if periods:  # keep chunks short enough for the wall-time limit
            dt = min(dt, max(dt0, 50.0 * min(periods)))
        # a bound subsystem of three or more stars must also be dynamically stable
        if (
            len(trees) > 1
            and sim.t > t_cross
            and _separated(trees, masses, P, V, r0)
            and all(
                min(stability_ratios(t, masses, P, V), default=2.0) > 1.0
                for t in trees
                if not isinstance(t, int) and len(_star_list(t)) >= 3
            )
        ):
            status = "resolved"
            break
        # a captured (fully bound) system counts as resolved once it has stayed dynamically
        # stable (MA01 > 1 for every nested triple) for 10 of its outermost orbits
        if len(trees) == 1 and sim.t > t_cross and len(masses) >= 3:
            ratios = stability_ratios(trees[0], masses, P, V)
            if ratios and min(ratios) > 1.0:
                if stable_since is None:
                    stable_since = sim.t
                el, _ = _orbit(trees[0], masses, P, V)
                p_out = units.sma_to_period(el.a, float(masses.sum()))
                if sim.t - stable_since >= 10.0 * p_out:
                    status = "resolved"
                    break
            else:
                stable_since = None
        if time.time() - wall0 > spec.max_wall_s:
            status = "wall_time"
            break
    P = np.array([[p.x, p.y, p.z] for p in sim.particles])
    V = np.array([[p.vx, p.vy, p.vz] for p in sim.particles])
    trees = hdet.detect(P, V, masses, units.G)
    outcome = classify(trees, initial, collided)
    row = {
        "id": enc["id"],
        "b_au": enc["b_au"],
        "q_au": enc["q_au"],
        "outcome": outcome,
        "status": "collision" if collided else status,
        "hierarchy": hdet.to_string(trees),
        "t_end_yr": sim.t,
        "energy_error": (
            abs(sim.energy() - E0 + (tides.dissipated if tides else 0.0)) / E_scale
            if not collided
            else None
        ),
        "wall_s": round(time.time() - wall0, 2),
        "scrambles": scrambles,
        "E0": float(E_total),
        "L0": float(L_total),
        "initial": enc,
    }
    if collided:
        row["collision_stars"] = sorted(int(k) + 1 for k in pair)
    if tides is not None:
        row["tidal_events"] = tides.events
        row["tidal_dE"] = float(tides.dissipated)
    orbits = {}
    for t in trees:
        stack = [t]
        while stack:
            node = stack.pop()
            if isinstance(node, int):
                continue
            try:
                el, _ = _orbit(node, masses, P, V)
            except ValueError:  # pragma: no cover - detect() only nests bound pairs
                continue
            orbits[hdet.to_string((node,))] = {"a_au": el.a, "e": el.e}
            stack += list(node)
    row["orbits"] = orbits
    bound = [t for t in trees if not isinstance(t, int) and len(_star_list(t)) >= 3]
    if bound:
        ratios = [r for t in bound for r in stability_ratios(t, masses, P, V)]
        row["ma01_min"] = float(min(ratios))
        row["stable"] = bool(min(ratios) > 1.0)
    return row


# --------------------------------------------------------------------------- ensembles


def _error_row(enc: dict, message: str) -> dict:
    return {
        "id": enc["id"],
        "b_au": enc.get("b_au"),
        "q_au": enc.get("q_au"),
        "outcome": "error",
        "status": "error",
        "error": message,
    }


_SPEC: dict = {}  # the scattering experiment, set once per worker process by _init


def _init(spec: ScatterSpec) -> None:
    _SPEC["spec"] = spec


def _task(enc):
    try:
        return run_encounter(_SPEC["spec"], enc)
    except Exception as exc:  # noqa: BLE001  (keep the ensemble going)
        return _error_row(enc, repr(exc))


def run_scatter(
    cfg: ScatterConfig, out: Path, processes: int = 1, resume: bool = False, progress=None
) -> list[dict]:
    """Run (or resume) the encounters; rows are appended to out/rows.jsonl as they finish. An
    encounter whose worker process dies is recorded as an error."""
    import platform

    from . import __version__
    from .ensemble import read_jsonl, resolve_processes, rewrite_jsonl, row_line
    from .parallel import run_tasks

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    spec = cfg.scatter
    meta = out / "scatter.yaml"
    rows_file = out / "rows.jsonl"
    stored = json.loads(json.dumps(cfg.model_dump()))
    done: dict[int, dict] = {}
    if rows_file.exists():
        if not resume:
            raise FileExistsError(f"{rows_file} exists; use --resume or another folder")
        old = yaml.safe_load(meta.read_text()) or {}
        o, c = dict(old.get("config", {}).get("scatter", {})), dict(stored["scatter"])
        for key in ("n", "max_wall_s"):  # more encounters / more computing time are allowed
            o.pop(key, None)
            c.pop(key, None)
        if o != c:
            raise ValueError("the scatter configuration changed since the run was started")
        raw, bad = read_jsonl(rows_file, null_as_nan=False)
        for r in raw:
            done[r["id"]] = r
        if bad or len(done) < len(raw):  # a torn last line (or duplicates): rewrite
            rewrite_jsonl(rows_file, list(done.values()))
    meta.write_text(
        yaml.safe_dump(
            {
                "config": stored,
                "v_c_kms": critical_velocity_kms(spec),
                "swingersim": __version__,
                "rebound": rebound.__version__,
                "python": platform.python_version(),
            }
        )
    )
    todo = [e for e in draw_encounters(spec) if e["id"] not in done]
    rows = [r for r in done.values() if r["id"] < spec.n]
    with rows_file.open("a") as fh:
        for row in run_tasks(_task, todo, resolve_processes(processes), _error_row, _init, (spec,)):
            fh.write(row_line(row))
            fh.flush()
            rows.append(row)
            if progress:
                progress(f"  [{len(rows)}/{spec.n}] encounter {row['id']}: {row['outcome']}")
    rows.sort(key=lambda r: r["id"])
    return rows


UNFINISHED = ("wall_time", "unresolved")  # stopped by the wall-time limit or at t_max


def outcome_key(r: dict) -> str:
    """Outcome class of a finished encounter, with the stability of a bound subsystem."""
    key = r["outcome"]
    if r.get("stable") is not None and (key.startswith("bound") or "triple" in key):
        key += " (stable)" if r["stable"] else " (unstable)"
    return key


def cross_section_weights(spec: ScatterSpec, rows: list[dict]) -> np.ndarray | None:
    """Weights that turn the drawn encounters into a cross-section sigma = pi b_max^2 x
    (weighted fraction): 1 for b drawn uniformly in area; 2 b / b_max for b drawn uniformly
    in b (b_distribution linear: the pdf 1/b_max instead of 2 b / b_max^2). None when there is
    no b_max (fixed b, or parabolic encounters drawn in pericentre)."""
    bm = b_max(spec)
    if bm is None or spec.b_au is not None:
        return None
    if spec.b_distribution == "area":
        return np.ones(len(rows))
    b = np.array([r.get("b_au") if r.get("b_au") is not None else np.nan for r in rows], float)
    return 2.0 * b / bm


def summarize_scatter(cfg: ScatterConfig, rows: list[dict]) -> str:
    """Outcome fractions and cross-sections of the finished encounters. Encounters that ended
    with an error are left out of the fractions (their number is given); encounters stopped
    by the wall-time limit or at t_max before the outcome was settled are listed separately
    with their provisional outcome and count in the denominator."""
    from .ensemble import effective_size, wilson

    spec = cfg.scatter
    errors = [r for r in rows if r.get("status") == "error" or r.get("outcome") == "error"]
    valid = [r for r in rows if not (r.get("status") == "error" or r.get("outcome") == "error")]
    n = len(valid)
    bm = b_max(spec)
    w = cross_section_weights(spec, valid)
    lines = [
        f"{cfg.name}: {len(rows)} encounters",
        (
            f"  v_inf = {v_inf_auyr(spec) / KMS_TO_AUYR:.4g} km/s, v_c = "
            f"{critical_velocity_kms(spec):.4g} km/s"
        ),
    ]
    if w is not None:
        how = (
            "weighted by 2 b / b_max (b uniform in b)"
            if spec.b_distribution == "linear"
            else "b uniform in area"
        )
        lines.append(f"  b_max = {bm:.4g} AU (cross-sections: sigma = pi b_max^2 fraction, {how})")
    if errors:
        lines.append(
            f"  {len(errors)} encounters ended with an error ({100 * len(errors) / len(rows):.1f} "
            "%): left out of the fractions"
        )
    lines.append("")
    lines.append(
        f"  {'outcome':44s} {'n':>6s} {'fraction':>9s}  {'1-sigma interval':>18s}"
        + ("  sigma [AU^2]" if w is not None else "")
    )
    keys = []
    for r in valid:
        k = outcome_key(r)
        if r.get("status") in UNFINISHED:
            k = f"unfinished ({r['status']}): {k}"
        keys.append(k)
    keys_arr = np.array(keys, dtype=object)
    counts: dict[str, int] = {}
    for k in keys:
        counts[k] = counts.get(k, 0) + 1
    for key, k in sorted(counts.items(), key=lambda x: (x[0].startswith("unfinished"), -x[1])):
        f, lo, hi = wilson(k, n)
        s = f"  {key:44s} {k:6d} {f:9.3f}  [{lo:.3f}, {hi:.3f}]"
        if w is not None:
            ok = np.isfinite(w)
            sel = (keys_arr == key) & ok
            fw = float(np.sum(w[sel]) / np.sum(w[ok])) if np.sum(w[ok]) > 0 else math.nan
            if spec.b_distribution == "linear":  # interval of the weighted fraction
                n_eff = effective_size(w[ok])
                _, lo_w, hi_w = wilson(fw * n_eff, n_eff)
                s += f"  {math.pi * bm * bm * fw:.4g} ({math.pi * bm * bm * lo_w:.3g}-"
                s += f"{math.pi * bm * bm * hi_w:.3g})"
            else:
                s += f"  {math.pi * bm * bm * fw:.4g}"
        lines.append(s)
    err = [
        r["energy_error"]
        for r in valid
        if r.get("energy_error") is not None and math.isfinite(r["energy_error"])
    ]
    if err:
        lines.append("")
        lines.append(f"  energy error: median {np.median(err):.1e}, max {np.max(err):.1e}")
    wall = [r["wall_s"] for r in rows if r.get("wall_s") is not None]
    if wall:
        lines.append(f"  wall time per encounter: median {np.median(wall):.2f} s")
    return "\n".join(lines)
