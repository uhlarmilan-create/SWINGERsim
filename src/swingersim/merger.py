"""Continuation of a run after two stars merge (physics.mergers).

A merger ends the run of a configuration: a collision of the N-body backend (two stars touch;
REBOUND collision search), or a common envelope or a contact binary that ends in a merger
(physics.mass_transfer: events common_envelope_merger, contact_merger_darwin,
contact_merger_l2). With physics.mergers SWINGERsim replaces the two stars by one and continues
with a run of the merged system (a quadruple becomes a triple) to the end of the original run.

Physics of the merger (see the review of Schneider et al. 2025, arXiv:2509.18421):

* collisions in bound few-body systems are slow (relative velocity ~ escape velocity): the two
  stars merge, they do not destroy each other. The merger ejects a small fraction of the mass,
  phi = 0.3 q / (1 + q)^2 (q = M2 / M1 <= 1; head-on collisions of main-sequence stars,
  Glebbeek et al. 2013: 2-7.5 %; Lombardi et al. 2002: at most ~8 %). A common-envelope merger
  keeps the mass of the energy formalism (swingersim.commonenv, Hurley et al. 2002).
* the product sits at the centre of mass of the two stars with their total momentum; the mass
  leaves at once, so the outer orbits feel an impulsive mass loss (Hadjidemetriou 1963) - the
  new elements follow from the unchanged positions and velocities with the new mass.
* rejuvenation: the product mixes fresh hydrogen into its core. Its fractional main-sequence age
  is f = r (M1 f1 + M2 f2) / (M1 + M2) with r = 0.81 for M <= 6 Msun to 0.64 for M >= 60 Msun
  (Glebbeek & Pols 2008; Schneider et al. 2016) or r = 0.1 of BSE (Hurley et al. 2002); its
  age on the track of the new mass is t_ZAMS + f (t_TAMS - t_ZAMS). A star past the main
  sequence keeps its equivalent evolutionary point on the track of the new mass.
* rotation: born near break-up, the product spins down quickly; it starts at spin_fraction of
  break-up (Omega_crit = sqrt(G M / (1.5 R)^3)).

Dynamically triggered mass transfer (physics.mass_transfer.dynamic_onset). The mass-transfer
histories (swingersim.masstransfer) are computed before the dynamics, for circular binaries of the
configured separations. The dynamics can drive a star over its lobe that no history expects: the
eccentricity cycles of a third body (Kozai-Lidov cycles with tidal friction, Eggleton &
Kisseleva-Eggleton 2006; Fabrycky & Tremaine 2007), a close passage. Tides then circularise the
binary quickly (Hurley et al. 2002, Sect. 2.6; triple codes such as TRES, Toonen et al. 2016,
assume circularisation at the onset of mass transfer). SWINGERsim ends the run there and continues
with the circularised binary: tidal dissipation conserves the orbital angular momentum, so
a_c = a (1 - e^2) (Hut 1981; the final separation of tidal friction, Fabrycky & Tremaine 2007,
eq. 2), or, option periastron, a_c = a (1 - e) (the donor still fills its lobe; the transfer
starts at once). The pair keeps its centre of mass, its orbital plane and its phase, its stars
rotate synchronously; the other orbits and stars are unchanged. The mass-transfer histories of
the new configuration start at that moment from the present masses and ages, so the binary may
become semi-detached, a contact binary (W UMa type, physics.mass_transfer.contact) or a common
envelope, while the outer orbits follow the dynamics of the whole system.

The hierarchy of the merged system is found from the positions and velocities (as
swingersim.diagnostics.hierarchy.detect); orbits whose stars are those of a configured orbit keep
its name. If the remaining stars are no longer one bound system, or a single star remains, the
run ends with the merger. Not carried over: discs (light curves), the fit and ensemble
sections; a mass-transfer episode of another binary restarts from its state at the merger.
"""

from __future__ import annotations

import math
import re
import warnings

import numpy as np

from . import units
from .config import Config

MT_MERGER_EVENTS = ("common_envelope_merger", "contact_merger_darwin", "contact_merger_l2")
MAX_RESTARTS = 16  # continuations of one run (mergers and dynamically started transfers)


# ------------------------------------------------------------------ merger physics
def mass_lost_fraction(spec, m1: float, m2: float) -> float:
    """Fraction of M1 + M2 the merger ejects (physics.mergers.mass_loss)."""
    if spec.mass_loss == "glebbeek2013":
        q = min(m1, m2) / max(m1, m2)
        return 0.3 * q / (1.0 + q) ** 2
    return float(spec.mass_loss)


def rejuvenation_factor(spec, m_new: float) -> float:
    """Ratio of the fractional main-sequence age of the product to the mass-weighted fractional
    age of the two stars (physics.mergers.rejuvenation)."""
    r = spec.rejuvenation
    if r == "hurley":
        return 0.1
    if r == "glebbeek":
        x = (math.log(max(m_new, 6.0)) - math.log(6.0)) / (math.log(60.0) - math.log(6.0))
        return 0.81 + (0.64 - 0.81) * min(x, 1.0)
    return float(r)


def _ms_window(sm, m: float) -> tuple[float, float] | None:
    """(age at the ZAMS, age at the TAMS) [Myr] of the track of mass m, or None (no tracks)."""
    if sm.kind == "mist":
        from .mist import EEP_ZAMS, mist_track

        tr = mist_track(float(m), sm.feh, sm.vvcrit)
        return float(np.interp(EEP_ZAMS, tr._x, tr._age)), float(tr.t_ms)
    if sm.kind == "sse":
        from .evolution import track

        return 0.0, float(track(float(m), sm.z).t_ms)
    return None


def merged_age(sm, ages, labels, masses, m_new: float, spec) -> float | None:
    """Age [Myr] of the merged star on the track of mass m_new (None: stellar models without
    ages). ages: the ages of the two stars on their tracks, labels: the track masses, masses:
    their present masses (weights)."""
    wins = [_ms_window(sm, lab) for lab in labels]
    new = _ms_window(sm, m_new)
    if new is None or any(w is None for w in wins):
        return None
    t0n, t1n = new
    past = [a > w[1] for a, w in zip(ages, wins, strict=True)]
    if any(past):
        # the most evolved star keeps its evolutionary point (equivalent evolutionary point of
        # MIST; SSE: the time after the TAMS scaled with the main-sequence lifetime)
        k = max(range(2), key=lambda q: (ages[q] - wins[q][1]) / max(wins[q][1] - wins[q][0], 1e-9))
        if sm.kind == "mist":
            from .mist import EvolvedBeyondGrid, mist_track

            try:
                eep = mist_track(float(labels[k]), sm.feh, sm.vvcrit).eep_of(ages[k])
            except (EvolvedBeyondGrid, ValueError):
                return None
            tr = mist_track(float(m_new), sm.feh, sm.vvcrit)
            return float(np.interp(eep, tr._x, tr._age))
        from .evolution import EvolvedBeyondGrid, track

        try:  # SSE: same phase and fraction of it (Hertzsprung gap, giant branch)
            ph, x = track(float(labels[k]), sm.z).phase_of(ages[k])
        except (EvolvedBeyondGrid, ValueError):
            return None
        tn = track(float(m_new), sm.z)
        if ph == 2:
            return tn.t_ms + x * (tn.t_bgb - tn.t_ms)
        return tn.t_bgb + min(x, 0.999) * (tn.t_hei - tn.t_bgb)
    f = [
        min(max((a - w[0]) / max(w[1] - w[0], 1e-9), 0.0), 1.0)
        for a, w in zip(ages, wins, strict=True)
    ]
    f_new = (
        rejuvenation_factor(spec, m_new)
        * (masses[0] * f[0] + masses[1] * f[1])
        / (masses[0] + masses[1])
    )
    return t0n + f_new * (t1n - t0n)


def _star_age(cfg: Config, k: int, t_yr: float) -> tuple[float, float]:
    """(age [Myr] on its track, track mass) of star k at run time t_yr."""
    from . import masstransfer
    from .hierarchy import build_tree
    from .stellar import age_offset

    m_lab = build_tree(cfg)[1]
    age = cfg.stellar_model.age_myr + age_offset(cfg, f"m{k + 1}")
    if cfg.physics.evolution:
        age += t_yr / 1e6
    if masstransfer.enabled(cfg):
        for h in masstransfer.histories(cfg).values():
            if k in h.pair and getattr(h, "age", None) is not None and len(h.t):
                q = h.pair.index(k)
                return float(np.interp(t_yr, h.t, h.age[q])), float(np.interp(t_yr, h.t, h.m_ev[q]))
    return float(age), float(m_lab[k])


# ------------------------------------------------------------------ the merged system
def _orbit_name(stars_old: frozenset, old_nodes: dict, used: set) -> str:
    name = old_nodes.get(stars_old)
    if name is None or name in used:
        n = 1
        while f"post{n}" in used or f"post{n}" in old_nodes.values():
            n += 1
        name = f"post{n}"
    used.add(name)
    return name


def post_merger_config(
    cfg: Config,
    state: dict,
    i: int,
    j: int,
    m_new: float,
    t_yr: float,
    t_end: float,
    kind: str,
) -> tuple[Config | None, dict]:
    """Configuration of the system after stars i and j merged at run time t_yr (state: pos,
    vel, masses, spins of all stars at that moment) and the details of the merger. Returns
    (None, details with 'stop' = reason) when the run cannot continue."""
    from .diagnostics import hierarchy as hdet
    from .hierarchy import Node, node_elements, skeleton
    from .stellar import star_properties

    spec = cfg.physics.mergers
    i, j = min(i, j), max(i, j)
    pos, vel = np.asarray(state["pos"], float), np.asarray(state["vel"], float)
    masses = np.asarray(state["masses"], float)
    spins = state.get("spins")
    m1, m2 = float(masses[i]), float(masses[j])
    a1, l1 = _star_age(cfg, i, t_yr)
    a2, l2 = _star_age(cfg, j, t_yr)
    age_new = merged_age(cfg.stellar_model, (a1, a2), (l1, l2), (m1, m2), m_new, spec)
    info = {
        "time_yr": float(t_yr),
        "kind": kind,
        "stars": [f"m{i + 1}", f"m{j + 1}"],
        "masses_msun": [m1, m2],
        "mass_msun": float(m_new),
        "mass_lost_msun": float(m1 + m2 - m_new),
        "ages_myr": [a1, a2],
        "age_myr": age_new,
    }
    # the product: centre of mass and total momentum of the two stars
    keep = [k for k in range(len(masses)) if k != j]
    p_new, v_new = pos[keep].copy(), vel[keep].copy()
    m_post = masses[keep].copy()
    ii = keep.index(i)
    p_new[ii] = (m1 * pos[i] + m2 * pos[j]) / (m1 + m2)
    v_new[ii] = (m1 * vel[i] + m2 * vel[j]) / (m1 + m2)
    m_post[ii] = m_new
    if len(keep) < 2:
        info["stop"] = "a single star remains"
        return None, info
    trees = hdet.detect(p_new, v_new, m_post, units.G)
    if len(trees) != 1:
        info["stop"] = "the remaining stars are not one bound system: " + hdet.to_string(trees)
        return None, info

    # orbit names: a configured orbit with the same stars keeps its name
    old_root = skeleton(cfg.tree)
    old_nodes = {frozenset(n.stars()): n.name for n in old_root.nodes()}

    def old_stars(t) -> frozenset:
        if isinstance(t, int):
            k = keep[t]
            return frozenset({i, j}) if k == i else frozenset({k})
        return old_stars(t[0]) | old_stars(t[1])

    used: set = set()

    by_name = {n.name: n for n in old_root.nodes()}

    def stars_of(c) -> frozenset:
        return frozenset({c}) if isinstance(c, int) else frozenset(c.stars())

    def build(t):
        if isinstance(t, int):
            return t
        name = _orbit_name(old_stars(t), old_nodes, used)
        old = by_name.get(name)
        if old is not None and (
            stars_of(old.left) == old_stars(t[1]) or stars_of(old.right) == old_stars(t[0])
        ):
            t = (t[1], t[0])  # the configured order of the children: the same relative vector
        return Node(name, build(t[0]), build(t[1]))

    root = build(trees[0])
    hierarchy, orbits = {}, {}
    for node in root.nodes():
        el = node_elements(node, p_new, v_new, m_post)
        if not (np.isfinite(el.a) and el.a > 0 and el.e < 1):
            info["stop"] = f"orbit {node.name} is not bound after the merger"
            return None, info
        hierarchy[node.name] = [
            c.name if isinstance(c, Node) else f"m{c + 1}" for c in (node.left, node.right)
        ]
        orbits[node.name] = {
            "a": float(el.a),
            "e": float(el.e),
            "i": math.degrees(el.inc),
            "Omega": math.degrees(el.Omega) % 360.0,
            "omega": math.degrees(el.omega) % 360.0,
            "M0": math.degrees(el.M) % 360.0,
        }

    d = cfg.model_dump()
    old_name = {f"m{n + 1}": f"m{k + 1}" for n, k in enumerate(keep)}  # new name -> old name
    stars_old = d.get("stars") or {}
    stars = {}
    for new, old in old_name.items():
        s = dict(stars_old.get(old) or {})
        if new == f"m{ii + 1}":  # the product: none of the old star's structure overrides
            s = {"age_offset_myr": 0.0}
        stars[new] = s
    d.update(
        topology="custom",
        hierarchy=hierarchy,
        masses={f"m{n + 1}": float(m) for n, m in enumerate(m_post)},
        orbits=orbits,
        stars=stars,
        discs={},
        ensemble=None,
        fit=None,
        fit_epoch_bjd_tdb=None,
    )
    sm = d["stellar_model"]
    if cfg.physics.evolution and sm.get("age_myr") is not None:
        sm["age_myr"] = float(cfg.stellar_model.age_myr + t_yr / 1e6)
    if age_new is not None and sm.get("age_myr") is not None:
        stars[f"m{ii + 1}"]["age_offset_myr"] = float(age_new - sm["age_myr"])
    for new, old in old_name.items():  # ages of the other stars: their offsets stay
        if new != f"m{ii + 1}":
            stars[new]["age_offset_myr"] = (stars_old.get(old) or {}).get("age_offset_myr", 0.0)
    ph = d["physics"]
    ml = ph.get("mass_loss")
    if ml:
        if ml.get("source") == "track":
            ml["masses"] = "current"  # the masses at the merger are present masses
        if ml.get("rates"):
            inv = {old: new for new, old in old_name.items()}
            ml["rates"] = {inv[k]: v for k, v in ml["rates"].items() if k in inv}
    if ph.get("gas"):
        ph["gas"]["orbits"] = [o for o in ph["gas"]["orbits"] if o in orbits]
        if not ph["gas"]["orbits"]:
            ph["gas"] = None
    run = d["run"]
    total = float(cfg.run.t_max_yr)
    remaining = float(t_end - t_yr)
    if remaining * total <= 0 or abs(remaining) < 1e-9 * abs(total):
        info["stop"] = "the merger happened at the end of the run"
        return None, info
    run["t_max_yr"] = remaining
    run["n_output"] = max(3, round((cfg.run.n_output - 1) * remaining / total) + 1)
    run["stop_when_e_reached"] = {
        k: v for k, v in (run.get("stop_when_e_reached") or {}).items() if k in orbits
    }
    new_cfg = Config.model_validate(d)
    new_cfg._base_dir = cfg.base_dir
    # spins: the other stars keep theirs (N-body), the product rotates at spin_fraction of
    # break-up about the orbit of its former binary
    if spins is not None or cfg.physics.tides or cfg.physics.rotation:
        assigned = _set_spins(
            new_cfg, root, p_new, v_new, m_post, spins, keep, i, j, ii, pos, vel, masses
        )
        new_cfg = Config.model_validate(new_cfg.model_dump())
        new_cfg._base_dir = cfg.base_dir
        new_cfg = _rotation_for_spins(cfg, new_cfg, assigned)
    pr = star_properties(new_cfg, f"m{ii + 1}", m_new)
    info["radius_rsun"] = float(pr.radius_rsun)
    info["names"] = old_name
    info["hierarchy"] = _remap_hierarchy(hdet.to_string(trees), old_name)  # original numbers
    info["product"] = f"m{i + 1}"  # the product keeps the number of the first star
    # the orbits of the merged system with the original star names (for the 3D viewer)
    info["orbits_tree"] = {o: [old_name.get(c, c) for c in kids] for o, kids in hierarchy.items()}
    return new_cfg, info


def _set_spins(new_cfg, root, p_new, v_new, m_post, spins, keep, i, j, ii, pos, vel, masses):
    """P_rot, obliquity and spin azimuth of the stars of the merged system (in place); returns
    the names of the stars whose spin was set."""
    from .backends.nbody import parent_of
    from .hierarchy import node_h_unit
    from .stellar import star_properties

    spec = new_cfg.physics.mergers
    assigned = set()
    for n in range(len(keep)):
        name = f"m{n + 1}"
        st = new_cfg.stars.get(name)
        if st is None:
            from .config import StarSpec

            st = StarSpec()
            new_cfg.stars[name] = st
        parent = parent_of(root, n)
        h = node_h_unit(parent, p_new, v_new, m_post)
        if n == ii:  # the product: about the orbit of the merged pair
            L = np.cross(pos[j] - pos[i], vel[j] - vel[i])
            w_hat = L / np.linalg.norm(L) if np.linalg.norm(L) > 0 else h
            pr = star_properties(new_cfg, name, float(m_post[n]))
            r_eq = 1.5 * pr.radius_rsun * units.RSUN_AU
            w_mag = spec.spin_fraction * math.sqrt(units.G * m_post[n] / r_eq**3)
            w = w_mag * w_hat
        elif spins is not None:
            w = np.asarray(spins[keep[n]], float)
        else:
            continue
        wn = float(np.linalg.norm(w))
        if wn <= 0:
            continue
        u = w / wn
        nvec = np.cross([0.0, 0.0, 1.0], h)
        nvec = (
            nvec / np.linalg.norm(nvec) if np.linalg.norm(nvec) > 1e-12 else np.array([1.0, 0, 0])
        )
        cos_ob = float(np.clip(u @ h, -1.0, 1.0))
        perp = u - cos_ob * h
        az = math.degrees(math.atan2(perp @ np.cross(h, nvec), perp @ nvec)) % 360.0
        st.P_rot = 2.0 * math.pi / wn * units.DAYS_PER_YEAR
        st.obliquity = math.degrees(math.acos(cos_ob))
        st.spin_azimuth = az
        assigned.add(name)
    return assigned


NON_ROTATING_D = 1e30  # P_rot [d] of a star that does not rotate


def _rotation_for_spins(cfg: Config, new_cfg: Config, assigned: set) -> Config:
    """Tides without physics.rotation start from non-rotating stars, so the P_rot of the
    continued run would be ignored: switch physics.rotation on, with the stars whose spin was
    not set (assigned) not rotating, as at the start of the first run."""
    if not cfg.physics.tides or cfg.physics.rotation:
        return new_cfg
    d = new_cfg.model_dump()
    d["physics"]["rotation"] = True
    stars = d.get("stars") or {}
    for name in d["masses"]:
        if name not in assigned:
            st = dict(stars.get(name) or {})
            st.update(P_rot=NON_ROTATING_D, obliquity=0.0)
            stars[name] = st
    d["stars"] = stars
    out = Config.model_validate(d)
    out._base_dir = cfg.base_dir
    return out


def _state_from_elements(cfg: Config, res, masses: np.ndarray) -> dict:
    """Positions and velocities at the last output of a run without them (secular backends):
    from the elements of the orbits (mean anomaly 0 where it is not followed)."""
    from .hierarchy import build_tree, initialize
    from .kepler import Elements

    root, _ = build_tree(cfg)
    for node in root.nodes():
        o = res.orbits[node.name]

        def last(k, default=0.0, o=o):
            v = o.get(k)
            x = float(v[-1]) if v is not None and len(v) else float("nan")
            return x if np.isfinite(x) else default

        node.elements = Elements(
            last("a", node.elements.a),
            last("e", node.elements.e),
            last("inc", node.elements.inc),
            last("Omega", node.elements.Omega),
            last("omega", node.elements.omega),
            last("M", 0.0),
        )
    pos, vel = initialize(root, masses)
    return {"pos": pos, "vel": vel, "masses": masses, "spins": None}


def _masses_at(cfg: Config, t_yr: float) -> np.ndarray:
    """Masses of the stars at run time t_yr (mass loss, mass transfer)."""
    from . import massloss
    from .hierarchy import build_tree

    m_lab = build_tree(cfg)[1]
    if massloss.enabled(cfg):
        return np.asarray(massloss.masses_at(cfg, m_lab, t_yr), float)
    return np.asarray(massloss.initial_masses(cfg, m_lab), float)


def merger_of(cfg: Config, res) -> tuple[int, int, float, str, dict] | None:
    """(i, j, merged mass, kind, state) of the merger that ended the run res, or None."""
    spec = cfg.physics.mergers
    if res.event == "collision" and res.final_state is not None:
        names = res.event_detail.get("stars") or []
        if len(names) != 2:
            return None
        i, j = (int(n[1:]) - 1 for n in names)
        st = res.final_state
        m1, m2 = float(st["masses"][i]), float(st["masses"][j])
        m_new = (m1 + m2) * (1.0 - mass_lost_fraction(spec, m1, m2))
        return i, j, m_new, "collision", st
    if not (isinstance(res.event, str) and res.event.startswith("mass_transfer:")):
        return None
    from . import masstransfer

    hits = [
        h
        for h in masstransfer.histories(cfg).values()
        if getattr(h, "event", None) in MT_MERGER_EVENTS and h.t_event is not None
    ]
    if not hits:
        return None
    h = min(hits, key=lambda x: x.t_event)
    i, j = h.pair
    t = float(res.times[-1])
    if res.final_state is not None:
        st = dict(res.final_state)
    else:
        st = _state_from_elements(cfg, res, _masses_at(cfg, t))
    masses = np.asarray(st["masses"], float).copy()
    m1, m2 = h.masses_at(h.t_event)
    masses[i], masses[j] = m1, m2
    st["masses"] = masses
    if h.event == "common_envelope_merger" and h.ce and h.ce.get("merged_mass"):
        m_new = float(h.ce["merged_mass"])
    else:
        m_new = (m1 + m2) * (1.0 - mass_lost_fraction(spec, m1, m2))
    return i, j, m_new, h.event, st


# ------------------------------------------------------------------ dynamic onset of mass transfer
def _nested(t):
    """Structure of a tree (detected tuple or Node) as nested frozensets of star indices."""
    if isinstance(t, (int, np.integer)):
        return int(t)
    kids = (t.left, t.right) if hasattr(t, "left") else (t[0], t[1])
    return frozenset(_nested(k) for k in kids)


def onset_of(cfg: Config, res) -> tuple[str, int, int, dict] | None:
    """(orbit, i, j, state) of an unexpected Roche-lobe overflow in a two-star binary that ended
    the run res and that physics.mass_transfer.dynamic_onset continues; None otherwise."""
    from . import masstransfer
    from .hierarchy import skeleton

    if masstransfer.dynamic_onset(cfg) is None or res.event != "roche_overflow":
        return None
    names = (res.event_detail or {}).get("stars") or []
    stars = {int(str(n)[1:]) - 1 for n in names}
    if not stars:
        return None
    pair = None
    for node in skeleton(cfg.tree).nodes():
        two = isinstance(node.left, int) and isinstance(node.right, int)
        if two and stars <= {node.left, node.right}:
            pair = node
    if pair is None:  # a star orbiting a subsystem, or stars of two binaries at once
        return None
    t = float(res.times[-1])
    if res.final_state is not None:
        st = dict(res.final_state)
    else:
        st = _state_from_elements(cfg, res, _masses_at(cfg, t))
    return pair.name, pair.left, pair.right, st


def post_onset_config(
    cfg: Config, state: dict, name: str, i: int, j: int, t_yr: float, t_end: float
) -> tuple[Config | None, dict]:
    """Configuration of the system after binary `name` (stars i, j) started its mass transfer
    at run time t_yr: the binary circularised (physics.mass_transfer.dynamic_onset), the rest as
    in the state. Returns (None, details with 'stop') when the run cannot continue."""
    from . import masstransfer
    from .diagnostics import hierarchy as hdet
    from .hierarchy import node_elements, skeleton
    from .kepler import state_to_elements
    from .stellar import roche_lobe_eggleton

    mode = masstransfer.dynamic_onset(cfg)
    pos = np.array(state["pos"], float)
    vel = np.array(state["vel"], float)
    masses = np.asarray(state["masses"], float)
    mi, mj = float(masses[i]), float(masses[j])
    M = mi + mj
    rc, vc = (mi * pos[i] + mj * pos[j]) / M, (mi * vel[i] + mj * vel[j]) / M
    r, v = pos[j] - pos[i], vel[j] - vel[i]
    el = state_to_elements(units.G * M, r, v)
    a, e = float(el.a), float(el.e)
    info = {
        "time_yr": float(t_yr),
        "kind": "dynamic_onset",
        "orbit": name,
        "stars": [f"m{i + 1}", f"m{j + 1}"],
        "circularisation": mode,
        "e_before": e,
        "a_before_au": a,
        "P_before_d": units.sma_to_period(a, M) * units.DAYS_PER_YEAR if a > 0 else float("nan"),
        "names": {f"m{k + 1}": f"m{k + 1}" for k in range(len(masses))},
    }
    if not (np.isfinite(a) and a > 0 and e < 1):
        info["stop"] = f"orbit {name} is not bound at the overflow"
        return None, info
    if t_end - t_yr <= 1e-9 * abs(cfg.run.t_max_yr) or cfg.run.t_max_yr < 0:
        info["stop"] = (
            "a backward run"
            if cfg.run.t_max_yr < 0
            else "the overflow happened at the end of the run"
        )
        return None, info
    # the eccentricity at which the first star touches its lobe at periastron (the outputs of a
    # run sample the growth of e: the state may lie a little past the contact)
    radii = state.get("radii")
    if radii is None:
        from .backends.evolving import props_at
        from .hierarchy import build_tree

        props = props_at(cfg, build_tree(cfg)[1], t_yr)
        radii = [p.radius_rsun * units.RSUN_AU for p in props]
    touch = [
        1.0 - float(radii[k]) / (roche_lobe_eggleton(masses[k] / masses[o]) * a)
        for k, o in ((i, j), (j, i))
    ]
    e_on = min(e, max(min(touch), 0.0))
    info["e_onset"] = float(e_on)
    root = skeleton(cfg.tree)
    rhat = r / np.linalg.norm(r)
    nrm = np.cross(r, v)
    nrm = nrm / np.linalg.norm(nrm)

    def make(a_c: float) -> Config | None:
        """The configuration with the binary circular at separation a_c (None: info['stop'])."""
        p_new, v_new_all = pos.copy(), vel.copy()
        r_new = a_c * rhat
        v_new = math.sqrt(units.G * M / a_c) * np.cross(nrm, rhat)
        p_new[i], p_new[j] = rc - mj / M * r_new, rc + mi / M * r_new
        v_new_all[i], v_new_all[j] = vc - mj / M * v_new, vc + mi / M * v_new
        p_c = units.sma_to_period(a_c, M) * units.DAYS_PER_YEAR
        info.update(a_after_au=float(a_c), P_after_d=float(p_c))
        info["label"] = (
            f"orbit {name} circularised at the overflow of {'+'.join(info['stars'])}: "
            f"e {e_on:.3f} -> 0, P {info['P_before_d']:.4g} -> {p_c:.4g} d"
        )
        trees = hdet.detect(p_new, v_new_all, masses, units.G)
        info["hierarchy"] = hdet.to_string(trees)
        if len(trees) != 1 or _nested(trees[0]) != _nested(root):
            info["stop"] = "the hierarchy is not the configured one: " + info["hierarchy"]
            return None
        d = cfg.model_dump()
        for node in root.nodes():
            o = d["orbits"][node.name]
            en = node_elements(node, p_new, v_new_all, masses)
            if not (np.isfinite(en.a) and en.a > 0 and en.e < 1):
                info["stop"] = f"orbit {node.name} is not bound at the overflow"
                return None
            circ = node.name == name
            o.update(
                P=None,
                a=float(a_c if circ else en.a),
                e=0.0 if circ else float(en.e),
                i=math.degrees(en.inc),
                Omega=math.degrees(en.Omega) % 360.0,
                omega=math.degrees(en.omega) % 360.0,
                M0=math.degrees(en.M) % 360.0,
            )
        d.update(
            masses={f"m{k + 1}": float(m) for k, m in enumerate(masses)},
            discs={},
            ensemble=None,
            fit=None,
            fit_epoch_bjd_tdb=None,
        )
        sm = d["stellar_model"]
        if cfg.physics.evolution and sm.get("age_myr") is not None:
            sm["age_myr"] = float(cfg.stellar_model.age_myr + t_yr / 1e6)
        ml = d["physics"].get("mass_loss")
        if ml and ml.get("source") == "track":
            ml["masses"] = "current"
        run = d["run"]
        total = float(cfg.run.t_max_yr)
        remaining = float(t_end - t_yr)
        run["t_max_yr"] = remaining
        run["n_output"] = max(3, round((cfg.run.n_output - 1) * remaining / total) + 1)
        stars = d.get("stars") or {}
        d["stars"] = stars
        for k in (i, j):  # tides synchronise and align the spins of the circularised binary
            st = dict(stars.get(f"m{k + 1}") or {})
            st.update(P_rot="sync", obliquity=0.0)
            stars[f"m{k + 1}"] = st
        new_cfg = Config.model_validate(d)
        new_cfg._base_dir = cfg.base_dir
        spins = state.get("spins")
        assigned = {f"m{i + 1}", f"m{j + 1}"}
        if spins is not None:  # the other stars keep their spins (N-body tides)
            keep = list(range(len(masses)))
            assigned |= _set_spins(new_cfg, skeleton(new_cfg.tree), p_new, v_new_all, masses,
                                   spins, keep, i, j, -1, pos, vel, masses)  # fmt: skip
            for k in (i, j):
                new_cfg.stars[f"m{k + 1}"].P_rot = "sync"
                new_cfg.stars[f"m{k + 1}"].obliquity = 0.0
            new_cfg = Config.model_validate(new_cfg.model_dump())
            new_cfg._base_dir = cfg.base_dir
        return _rotation_for_spins(cfg, new_cfg, assigned)

    a_c = a * (1.0 - e_on * e_on) if mode == "angular_momentum" else a * (1.0 - e_on)
    new_cfg = make(a_c)
    if new_cfg is None:
        return None, info
    h = masstransfer.histories(new_cfg).get(name)
    fill0 = h.fill_at(0.0) if h is not None else 0.0
    if mode == "periastron" and 0.0 < fill0 < 1.0 + 1e-7:
        # the donor at its lobe in the circular history of the new configuration (whose stellar
        # model may differ slightly from the radius the dynamics used)
        new_cfg = make(a_c * fill0 / (1.0 + 1e-6))
        if new_cfg is None:
            return None, info
        h = masstransfer.histories(new_cfg).get(name)
        fill0 = h.fill_at(0.0) if h is not None else 0.0
    # the circular history of the new binary must see the overflow the dynamics saw: else the
    # same overflow would end the next run at once
    info["fill_after"] = float(fill0)
    if h is not None and h.t_onset is None and fill0 >= 1.0 - masstransfer.FILL_TOL:
        info["stop"] = "the circular binary fills its lobe, but its mass transfer is not followed"
        return None, info
    return new_cfg, info


# ------------------------------------------------------------------ joining the runs
def _remap_hierarchy(text: str, names: dict) -> str:
    """Star numbers of a hierarchy string of the merged system -> those of the original."""
    num = {int(new[1:]): int(old[1:]) for new, old in names.items()}
    return re.sub(r"\d+", lambda mo: str(num.get(int(mo.group()), mo.group())), text)


def _nan(v) -> bool:
    return v is None or (isinstance(v, float) and math.isnan(v))


def _finite(v) -> bool:
    return isinstance(v, (int, float)) and math.isfinite(v)


def _join(a, b, info: dict, t0: float):
    """Run a followed by run b of the merged system (b starts at run time t0 of a)."""
    names = info["names"]  # new star name -> old star name
    n_a, n_b = len(a.times), len(b.times) - 1  # b's first row is the state at the merger

    def cat(x, y):
        x = np.asarray(x, float) if x is not None else np.full(n_a, np.nan)
        y = np.asarray(y, float)[1:] if y is not None else np.full(n_b, np.nan)
        return np.concatenate([x, y])

    def series(da: dict, db: dict, key_map=None):
        out = {}
        keys = list(da) + [k for k in db if (key_map or {}).get(k, k) not in da]
        for k in keys:
            if key_map and k in da:  # a star of run a: its name in run b, none if absorbed
                kb = next((n for n, o in key_map.items() if o == k), None)
            else:
                kb = k
            sa, sb = da.get(k), (db.get(kb) if kb is not None else None)
            if isinstance(sa, dict) or isinstance(sb, dict):
                sa, sb = sa or {}, sb or {}
                out[k] = {
                    e: cat(sa.get(e, np.full(n_a, np.nan)), sb.get(e, np.full(n_b + 1, np.nan)))
                    for e in dict.fromkeys([*sa, *sb])
                }
            else:
                out[k] = cat(sa, sb if sb is not None else np.full(n_b + 1, np.nan))
        return out

    a.times = np.concatenate([a.times, t0 + np.asarray(b.times, float)[1:]])
    a.orbits = series(a.orbits, b.orbits)
    a.mutual = series(a.mutual, b.mutual)
    a.energy_error = cat(a.energy_error, b.energy_error)
    if a.angmom_error is not None or b.angmom_error is not None:
        a.angmom_error = cat(a.angmom_error, b.angmom_error)
    a.spins = series(a.spins, b.spins, names)
    a.stellar = series(a.stellar, b.stellar, names)
    if a.megno is not None or b.megno is not None:
        a.megno = cat(a.megno, b.megno)
    a.hierarchy = list(a.hierarchy) + [_remap_hierarchy(h, names) for h in b.hierarchy[1:]]
    onset = info["kind"] == "dynamic_onset"
    a.transitions = list(a.transitions) + [
        {
            "time_yr": float(t0),
            "kind": "mass_transfer_onset" if onset else "merger",
            "from": a.hierarchy[n_a - 1] if n_a else "",
            "to": info["hierarchy"],
            "label": info["label"]
            if onset
            else f"{info['stars'][0]} + {info['stars'][1]} -> {info['mass_msun']:.3f} Msun",
        }
    ]
    for tr in b.transitions:
        tr = dict(tr)
        tr["time_yr"] = float(t0 + tr["time_yr"])
        for k in ("from", "to"):
            tr[k] = _remap_hierarchy(str(tr.get(k, "")), names)
        a.transitions.append(tr)
    for name, ex in b.extremes.items():  # extremes of the whole run
        ea = a.extremes.setdefault(name, {})
        for k, v in ex.items():
            if k == "e_max":
                if not (v <= ea.get(k, -np.inf)):
                    ea[k], ea["t_e_max_yr"] = v, t0 + ex.get("t_e_max_yr", 0.0)
            elif k == "q_min_au":
                if not (v >= ea.get(k, np.inf)):
                    ea[k], ea["t_q_min_yr"] = v, t0 + ex.get("t_q_min_yr", 0.0)
            elif k not in ("t_e_max_yr", "t_q_min_yr") and _nan(ea.get(k)) and _finite(v):
                ea[k] = t0 + v  # a threshold first crossed in the later run
    for new, r in b.roche.items():
        old = names.get(new, new)
        ra = a.roche.setdefault(old, {})
        v = r.get("fill_max")
        if _finite(v) and not (v <= ra.get("fill_max", -np.inf)):
            ra["fill_max"] = float(v)
            if _finite(r.get("t_fill_max_yr")):
                ra["t_fill_max_yr"] = float(t0 + r["t_fill_max_yr"])
        for k, v in r.items():  # other fill keys: the larger one
            if k.startswith("fill") and k != "fill_max" and isinstance(v, float):
                ra[k] = max(ra.get(k, -np.inf), v)

    def shifted(lst):  # events of the later run: its times from its own start, its stars
        out = []  # with the names of this run
        for m in lst:
            m = {**m, "time_yr": float(t0 + m["time_yr"])}
            if "stars" in m:
                m["stars"] = [names.get(x, x) for x in m["stars"]]
            if "product" in m:
                m["product"] = names.get(m["product"], m["product"])
            if isinstance(m.get("orbits_tree"), dict):
                m["orbits_tree"] = {
                    o: [names.get(c, c) for c in kids] for o, kids in m["orbits_tree"].items()
                }
            out.append(m)
        return out

    own = {k: v for k, v in info.items() if k not in ("names",)}
    mergers = list(a.system.get("_mergers", [])) + ([] if onset else [own])
    mergers += shifted(b.system.get("_mergers", []))
    onsets = list(a.system.get("_onsets", [])) + ([own] if onset else [])
    onsets += shifted(b.system.get("_onsets", []))
    for k in [k for k in a.system if k.startswith(("merger_", "onset_"))]:
        del a.system[k]
    if mergers:
        a.system["_mergers"] = mergers
        a.system["mergers_n"] = len(mergers)
    for n, mg in enumerate(mergers, 1):
        a.system[f"merger_{n}_time_yr"] = mg["time_yr"]
        a.system[f"merger_{n}_kind"] = mg["kind"]
        a.system[f"merger_{n}_stars"] = "+".join(mg["stars"])
        a.system[f"merger_{n}_mass_msun"] = mg["mass_msun"]
    if onsets:
        a.system["_onsets"] = onsets
        a.system["onsets_n"] = len(onsets)
    for n, og in enumerate(onsets, 1):
        a.system[f"onset_{n}_time_yr"] = og["time_yr"]
        a.system[f"onset_{n}_orbit"] = og["orbit"]
        a.system[f"onset_{n}_e_before"] = og["e_before"]
        a.system[f"onset_{n}_P_before_d"] = og["P_before_d"]
        a.system[f"onset_{n}_P_after_d"] = og["P_after_d"]
    if onset:  # the mass transfer of the circularised system replaces that of the first run
        for k in [k for k in a.system if k.startswith("mass_transfer_")]:
            del a.system[k]
        for k, v in b.system.items():
            if k.startswith("mass_transfer_"):
                shift = k.endswith("_yr") and isinstance(v, (int, float)) and v is not None
                a.system[k] = float(t0 + v) if shift else v
    prefix = "post_onset_" if onset else "post_merger_"
    skip = ("_mergers", "mergers_n", "merger_", "_onsets", "onsets_n", "onset_", "post_")
    for k, v in b.system.items():  # diagnostics of the continued system
        if onset and k.startswith("mass_transfer_"):
            continue
        if not k.startswith(skip):
            a.system[f"{prefix}{k}"] = v
        elif k.startswith("post_"):
            a.system[k] = v
    a.event, a.event_time, a.event_detail = (
        b.event,
        (None if b.event_time is None else t0 + b.event_time),
        b.event_detail,
    )
    a.final_state = None
    a.wall_s = (a.wall_s or 0.0) + (b.wall_s or 0.0)
    a.mass_loss = a.mass_loss or b.mass_loss
    a.dissipative = a.dissipative or b.dissipative
    return a


def run_with_mergers(cfg: Config, run_once, depth: int = 0):
    """Run the configuration with run_once (a backend); after a merger continue with the merged
    system (recursively) and join the runs."""
    from . import massloss
    from .backends.evolving import time_limit
    from .hierarchy import build_tree

    cfg = massloss.resolve_masses(cfg)
    res = run_once(cfg)
    found = merger_of(cfg, res) if cfg.physics.mergers is not None else None
    if found is None:
        onset = onset_of(cfg, res)
        if onset is None:
            return res
        return _continue_onset(cfg, res, onset, run_once, depth)
    i, j, m_new, kind, state = found
    t_merge = float(res.times[-1])
    t_end = float(time_limit(cfg, build_tree(cfg)[1], cfg.run.t_max_yr)[0])
    if kind != "collision":  # the mass-transfer event cut the run: its own end is t_max
        t_end = float(cfg.run.t_max_yr)
    new_cfg, info = post_merger_config(cfg, state, i, j, m_new, t_merge, t_end, kind)
    res.event_detail = {
        **(res.event_detail or {}),
        "merger": {k: v for k, v in info.items() if k not in ("names",)},
    }
    if new_cfg is None or depth >= MAX_RESTARTS:
        res.event = "merger"
        res.system["_mergers"] = [res.event_detail["merger"]]
        res.system["mergers_n"] = 1
        return res
    try:
        post = run_with_mergers(new_cfg, run_once, depth + 1)
    except (ValueError, RuntimeError, ArithmeticError) as exc:  # cannot run: end with the merger
        warnings.warn(f"run after the merger failed: {exc}", stacklevel=2)
        res.event = "merger"
        res.event_detail["merger"]["stop"] = f"run of the merged system failed: {exc}"
        res.system["_mergers"] = [res.event_detail["merger"]]
        res.system["mergers_n"] = 1
        return res
    return _join(res, post, info, t_merge)


def _continue_onset(cfg: Config, res, onset, run_once, depth: int):
    """The run res ended at a dynamically started Roche-lobe overflow: continue it with the
    circularised binary (physics.mass_transfer.dynamic_onset) and join the runs."""
    from .backends.evolving import time_limit
    from .hierarchy import build_tree

    name, i, j, state = onset
    t_on = float(res.times[-1])
    t_end = float(time_limit(cfg, build_tree(cfg)[1], cfg.run.t_max_yr)[0])
    new_cfg, info = post_onset_config(cfg, state, name, i, j, t_on, t_end)
    own = {k: v for k, v in info.items() if k not in ("names",)}
    res.event_detail = {**(res.event_detail or {}), "dynamic_onset": own}
    if new_cfg is None or depth >= MAX_RESTARTS:
        if new_cfg is not None:
            own["stop"] = f"more than {MAX_RESTARTS} continuations of the run"
        warnings.warn(f"dynamic onset of orbit {name} not continued: {own['stop']}", stacklevel=2)
        res.system["_onsets"] = [*res.system.get("_onsets", []), own]
        return res
    try:
        post = run_with_mergers(new_cfg, run_once, depth + 1)
    except (ValueError, RuntimeError, ArithmeticError) as exc:
        warnings.warn(f"run after the onset of mass transfer failed: {exc}", stacklevel=2)
        own["stop"] = f"run of the circularised system failed: {exc}"
        res.system["_onsets"] = [*res.system.get("_onsets", []), own]
        return res
    return _join(res, post, info, t_on)


__all__ = [
    "mass_lost_fraction",
    "merged_age",
    "merger_of",
    "onset_of",
    "post_merger_config",
    "post_onset_config",
    "rejuvenation_factor",
    "run_with_mergers",
]
