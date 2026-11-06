"""Secular backend for 3+1 quadruples and other nested hierarchies (swingersim.secular.nested).

Same output as the N-body and the 2+2 secular backend: mean elements, mean anomaly 0, no
hierarchy changes or collisions. Relativity and the tidal and rotational distortion of the
stars as apsidal motion; tidal dissipation (constant time lag, aligned spins) in every orbit of
two single stars.
"""

from __future__ import annotations

from itertools import combinations, pairwise

import numpy as np

from .. import massloss, masstransfer, units
from .. import masstransfer as _mt
from ..config import Config
from ..diagnostics import hierarchy as hdet
from ..diagnostics.binary import eclipse_summary
from ..diagnostics.extremes import THRESHOLDS_1ME, threshold_key
from ..diagnostics.orbits import mutual_inclination
from ..hierarchy import Node, build_tree, initialize, mass_of
from ..secular.nested import SecularNested, StarStructure, integrate
from ..secular.twotwo import E_SEED, elements_from_vectors, vectors_from_elements
from ..stellar import roche_lobe_eggleton, star_properties
from ..tides import with_structure_lags
from . import evolving as ev
from . import secular_common as sc
from .nbody import RunResult, parent_of, spin_vector


def build_nested(cfg: Config):
    root, masses = build_tree(cfg)
    props = [star_properties(cfg, f"m{k + 1}", m) for k, m in enumerate(masses)]
    m_dyn = massloss.initial_masses(cfg, masses)  # current masses at t = 0
    if not np.array_equal(m_dyn, masses):  # periods with the masses of the dynamics
        root = build_tree(cfg, kepler_masses=m_dyn)[0]
    pos, vel = initialize(root, m_dyn)
    structure = None
    if cfg.physics.tides:
        spins = [
            float(
                np.linalg.norm(spin_vector(cfg, f"m{k + 1}", parent_of(root, k), pos, vel, m_dyn))
            )
            for k in range(len(masses))
        ]
        if cfg.physics.tidal_lag != "fixed":
            props = with_structure_lags(cfg, props, _pairs(root), spins, m_dyn)
        structure = _structure(props, spins, cfg, m_dyn)
    gas = cfg.physics.gas
    if cfg.physics.evolution and gas is not None:
        raise ValueError("physics.evolution and physics.gas cannot be combined")
    sys_ = SecularNested(
        root,
        m_dyn,
        structure,
        gr=cfg.physics.gr,
        octupole=cfg.run.secular_octupole,
        gas=gas,
        gw=cfg.physics.gw,
    )
    # backward runs: spins that synchronise much faster than |t_max| held pseudo-synchronous
    sys_.backward_span = -float(cfg.run.t_max_yr) if cfg.run.t_max_yr < 0 else 0.0
    sys_._kargs = None
    y0 = []
    for node in sys_.nodes:
        el = node.elements
        K, p = vectors_from_elements(max(el.e, E_SEED), el.inc, el.Omega, el.omega)
        y0.extend([*K, *p])
    y0 = np.array(y0)
    if sys_.dissipative:
        y0 = np.concatenate([y0, sys_.extra_state()])
    gt = cfg.physics.galactic_tide
    if gt is not None:
        from ..environment import pole_vector, rho_au

        sys_.tide = (rho_au(gt.rho0_msun_pc3), pole_vector(gt.pole_i_deg, gt.pole_Omega_deg))
    if gas is not None:
        if gas.plane == "invariable":
            L0 = sys_.angular_momentum(y0)
            sys_.disc_normal = L0 / np.linalg.norm(L0)
        else:
            i, O = np.radians(gas.plane.i), np.radians(gas.plane.Omega)
            sys_.disc_normal = np.array([np.sin(i) * np.sin(O), -np.sin(i) * np.cos(O), np.cos(i)])
        y0 = np.concatenate([y0, sys_.a0])
    return sys_, y0, root, masses, props, (pos, vel)


def _pairs(root, a=None) -> list[tuple[int, int, float]]:
    """(i, j, a) of every orbit of two single stars."""
    out = []
    for q, nd in enumerate(root.nodes()):
        if not isinstance(nd.left, Node) and not isinstance(nd.right, Node):
            out.append((nd.left, nd.right, float(nd.elements.a if a is None else a[q])))
    return out


def _structure(props, spins, cfg=None, masses=None) -> dict[int, StarStructure]:
    from ..tides import k2_rotating

    def k2(k):
        if cfg is None:
            return props[k].apsidal_k2
        m = float(masses[k])
        return k2_rotating(cfg, k, props[k].apsidal_k2, float(spins[k]), m, props[k].radius_rsun)

    return {
        k: StarStructure(
            r_au=props[k].radius_rsun * units.RSUN_AU,
            apsidal_k2=k2(k),
            spin=float(spins[k]),
            rg2=props[k].rg2,
            lag_yr=props[k].lag_s / units.SECONDS_PER_YEAR,
        )
        for k in range(len(props))
    }


def _extremes_of(sys_, masses, props):
    """extremes_of(sol, t, t_stop, overflow) of a nested run (secular_common.dense_extremes);
    the semi-major axes as SecularNested.a_state: the gas phase carries all of them in the
    state, a dissipative orbit its own, the others keep a0."""
    n6 = 6 * sys_.n_orb
    a0 = sys_.a0.copy()

    def a_of(Y, q):
        if sys_.gas is not None and Y.shape[0] > n6:
            return Y[n6 + q]
        if q in sys_.a_index:
            return Y[sys_.a_index[q]]
        return a0[q]

    lobes = {}
    for node in sys_.nodes:
        lobes[node.name] = [
            (
                f"m{s + 1}",
                roche_lobe_eggleton(masses[s] / mass_of(o, masses)),
                props[s].radius_rsun * units.RSUN_AU,
            )
            for s, o in ((node.left, node.right), (node.right, node.left))
            if not isinstance(s, Node)  # single stars: lobe against the companion (subsystem)
        ]
    orbits = [(n, 6 * q + 3) for q, n in enumerate(sys_.names)]

    def extremes_of(sol, t, t_stop=None, overflow=None, n_fine=None):
        if n_fine is None:
            n_fine = 20 * (len(t) - 1) + 1
        return sc.dense_extremes(sol, t, orbits, a_of, lobes, n_fine, t_stop, overflow)

    return extremes_of


def _resync_slaved(sys_, y, spins) -> None:
    """Backward runs: the spins held pseudo-synchronous (swingersim.tides.slaved_spin) are set to
    the pseudo-synchronous rate of their orbit at the start of every segment; the spin angular
    momentum they gained or lost since (structure change, magnetic braking) is exchanged with
    their orbit, as the tides do in the forward direction (total angular momentum kept)."""
    from ..tides import pseudo_sync_factor

    sys_._kernel_args()  # decides which spins are held (sys_.slaved) for this structure
    for q, i in sys_.spin_index.items():
        nd = sys_.nodes[q]
        for j, s in enumerate((nd.left, nd.right)):
            if s not in sys_.slaved:
                continue
            st = sys_.structure[s]
            ia = sys_.a_index[q]
            m0, m1 = sys_.masses[nd.left], sys_.masses[nd.right]
            M, mu = m0 + m1, m0 * m1 / (m0 + m1)
            e = float(np.linalg.norm(y[6 * q + 3 : 6 * q + 6]))
            C = st.rg2 * sys_.masses[s] * st.r_au**2
            for _ in range(2):  # the orbit changes a little: once more at the new a
                a = y[ia]
                w_ps = np.sqrt(units.G * M / a**3) * pseudo_sync_factor(e)
                dJ = C * (spins[s] - w_ps)
                L = mu * np.sqrt(units.G * M * a * (1.0 - e * e))
                y[ia] = a * ((L + dJ) / L) ** 2
                spins[s] = w_ps
            y[i + j] = spins[s]
    sys_._kargs = None


def _stream_nested(cfg, sys_, y0, masses, props, times, encs=(), enc_counts=None):
    """Evolving stars (backends/evolving.py): operator splitting, streamed integration.
    Encounters (physics.encounters) are applied at the start of the segment that begins at
    their time: the segment edges are cut at the encounter times."""
    spins = {k: st.spin for k, st in sys_.structure.items()}
    sign = 1.0 if times[-1] >= 0 else -1.0
    enc_next = [0]

    def apply(pm, prev, y, mm, m_prev):
        t0 = ev.SEGMENT.get("t0")
        if encs and t0 is not None:  # the encounters due at the start of this segment
            k = enc_next[0]
            due = []
            while k < len(encs) and sign * encs[k].t_yr <= sign * t0 + 1e-9 * max(abs(t0), 1):
                due.append(encs[k])
                k += 1
            enc_next[0] = k
            if due:
                _kick(sys_, y, due, enc_counts)
        if np.any(mm != sys_.masses):
            sys_.set_masses(mm, y)
        if not sys_.structure:
            return
        now = sys_.spins_state(y)  # dissipative binaries carry their spins in the state
        for k in list(spins):
            w = ev.scale_spin(now[k], prev[k], pm[k], m_prev[k], mm[k])
            spins[k] = ev.brake_spin(cfg, k, w, pm[k], mm[k], ev.SEGMENT["dt"])
        for q, i in sys_.spin_index.items():
            nd = sys_.nodes[q]
            y[i], y[i + 1] = spins[nd.left], spins[nd.right]
        if cfg.physics.tidal_lag != "fixed":
            pm[:] = with_structure_lags(cfg, pm, _pairs(sys_.root, sys_.a_state(y)), spins, mm)
        sys_.set_structure(_structure(pm, [spins[k] for k in range(len(pm))], cfg, mm))
        if getattr(sys_, "backward_span", 0.0) > 0:
            _resync_slaved(sys_, y, spins)
            sys_.set_structure(_structure(pm, [spins[k] for k in range(len(pm))], cfg, mm))

    specs = []
    for q, node in enumerate(sys_.nodes):
        two = not isinstance(node.left, Node) and not isinstance(node.right, Node)
        lobes = None
        if not two and not (isinstance(node.left, Node) and isinstance(node.right, Node)):
            star, sub = (
                (node.left, node.right)
                if isinstance(node.right, Node)
                else (
                    node.right,
                    node.left,
                )
            )
            lobes = ((star, tuple(sub.stars())),)
        specs.append(
            ev.OrbitSpec(
                node.name,
                6 * q + 3,
                sys_.a_index.get(q),
                node.elements.a,
                (node.left, node.right) if two else None,
                lobes,
                a_of=lambda q=q: sys_.a0[q],
            )
        )
    return ev.stream_evolving(
        cfg,
        masses,
        props,
        times,
        y0,
        sys_.rhs,
        sys_.energy,
        apply,
        lambda: (dict(sys_.structure), sys_.mass_state()),
        specs,
        THRESHOLDS_1ME,
        threshold_key,
        sc.refine,
        roche_lobe_eggleton,
        angmom=sys_.angular_momentum,
        scale_orbits=sys_.scale_orbits,
        extra_edges=[e.t_yr for e in encs],
    )


def _encounters(cfg, sys_):
    """Encounter sequence of the run (physics.encounters), or []."""
    spec = cfg.physics.encounters
    if spec is None:
        return []
    from ..environment import encounter_sequence

    q_out = int(np.argmax(sys_.a0))
    m_sys = float(np.sum(sys_.masses))
    return encounter_sequence(spec, m_sys, float(sys_.a0[q_out]), float(cfg.run.t_max_yr))


def _kick(sys_, y, encs, counts) -> None:
    """Apply the encounters encs to the state y in place (first-order secular kicks)."""
    from ..environment import apply_encounter

    m_orb = sys_.orbit_masses()
    for enc in encs:
        counts["n"] += 1
        close = False
        for q in range(sys_.n_orb):
            sk, sp = slice(6 * q, 6 * q + 3), slice(6 * q + 3, 6 * q + 6)
            r = apply_encounter(enc, y[sk], y[sp], float(sys_.a0[q]), float(m_orb[q]))
            if r is None:
                close = True
                continue
            y[sk], y[sp] = r
        counts["close"] += int(close)


def _integrate_with_encounters(cfg, sys_, y0, times, encs, counts):
    """Secular integration interrupted by the encounters (operator splitting)."""
    from .evolving import Segment, SegmentedSolution

    sign = 1.0 if times[-1] >= 0 else -1.0
    cuts = [0.0] + [e.t_yr for e in encs] + [float(times[-1])]
    segs = []
    y = np.array(y0, float)
    k = 0
    for ta, tb in pairwise(cuts):
        if sign * (tb - ta) > 0:
            sol = integrate(sys_, y, np.array([ta, tb]), rtol=cfg.run.secular_rtol)
            if sol.status < 0:
                raise RuntimeError(f"secular integration failed: {sol.message}")
            segs.append(Segment(ta, tb, sol, None, None))
            y = sol.y[:, -1].copy()
        if k < len(encs):
            _kick(sys_, y, [encs[k]], counts)
            k += 1
    return SegmentedSolution(segs)


def run_secular_nested(cfg: Config) -> RunResult:
    sys_, y0, root, masses, props, (pos0, vel0) = build_nested(cfg)
    names = sys_.names
    encs = _encounters(cfg, sys_)
    enc_counts = {"n": 0, "close": 0}
    evolving = ev.time_varying(cfg)
    mass_loss = massloss.enabled(cfg)
    m0 = massloss.initial_masses(cfg, masses)
    stop_event = None
    psi_max = None
    overflow = None
    if evolving:
        t_end, stop_event = ev.time_limit(cfg, masses, cfg.run.t_max_yr)
        times = np.linspace(0.0, t_end, cfg.run.n_output)
        if mass_loss:
            psi_max = massloss.check_psi(
                cfg, root, masses, dict(zip(names, sys_.a0, strict=True)), times
            )
        stream = _stream_nested(cfg, sys_, y0, masses, props, times, encs, enc_counts)
        Y = stream.Y
        if stream.t_stop is not None:  # Roche-lobe overflow: the run ends there
            times, Y = sc.truncate_at_stream_stop(times, Y, stream)
            stop_event = "roche_overflow"
    else:
        times = np.linspace(0.0, cfg.run.t_max_yr, cfg.run.n_output)
        extremes_of = _extremes_of(sys_, masses, props)
        stop = _mt.stop_on_overflow(cfg)
        if encs:
            sol = _integrate_with_encounters(cfg, sys_, y0, times, encs, enc_counts)
            times, Y, extremes, roche, overflow = sc.extremes_until_overflow(
                sol,
                times,
                sol.sol(times),
                lambda s_, t, t_stop=None, over=None: extremes_of(
                    s_, t, t_stop, over, max(20 * len(t), 2000)
                ),
                stop,
            )
        else:  # in chunks: bounded memory (see backends/chunked.py)
            times, Y, extremes, roche, overflow = sc.integrate_until_overflow(
                lambda y, t: integrate(sys_, y, t, rtol=cfg.run.secular_rtol),
                y0,
                times,
                extremes_of,
                stop,
            )

    def use(j):
        st, ms = stream.snapshots[j]
        sys_.restore_mass_state(ms)
        sys_.set_structure(st)

    if evolving:  # the reference values refer to the structure at t = 0
        use(0)
    L0 = sys_.angular_momentum(y0)
    basis = sc.invariable_basis(L0)
    H0 = sys_.energy(y0)

    keys = ("a", "e", "inc", "Omega", "omega", "M", "inc_inv", "Omega_inv")
    orbits = {n: {k: np.zeros(len(times)) for k in keys} for n in names}
    k_hat = {n: np.zeros((len(times), 3)) for n in names}
    energy_error = np.zeros(len(times))
    angmom_error = np.zeros(len(times))
    n6 = 6 * len(names)
    for j in range(len(times)):
        y = Y[:, j]
        if evolving:
            use(j)
        a_rep = sys_.a0.copy()
        if mass_loss:  # report a at the masses of time t_j (the segment holds its masses)
            m_seg = sys_.orbit_masses()
            sys_mt = massloss.orbit_masses(root, massloss.masses_at(cfg, masses, float(times[j])))
            a_rep = a_rep * m_seg / np.array([sys_mt[n] for n in names])
            for n, a_h in masstransfer.orbit_a_at(cfg, float(times[j])).items():  # transfer
                a_rep[names.index(n)] = a_h
        for q, n in enumerate(names):
            K, p = y[6 * q : 6 * q + 3], y[6 * q + 3 : 6 * q + 6]
            e, inc, Om, om = elements_from_vectors(K, p)
            o = orbits[n]
            if q in sys_.a_index:
                a_rep[q] = y[sys_.a_index[q]]
                if mass_loss:
                    a_rep[q] *= m_seg[q] / sys_mt[n]
            o["a"][j], o["e"][j], o["inc"][j], o["Omega"][j], o["omega"][j] = (
                y[n6 + q] if (len(y) > n6 and sys_.gas is not None) else a_rep[q],
                e,
                inc,
                Om,
                om,
            )
            k = K / np.linalg.norm(K)
            k_hat[n][j] = k
            o["inc_inv"][j], o["Omega_inv"][j] = sc.invariable_angles(k, basis)
        h_ref = stream.h_ref[j] if evolving else H0
        energy_error[j] = abs(sys_.energy(y) / h_ref - 1.0) if h_ref != 0 else 0.0
        # the wind carries angular momentum (mass loss; magnetic braking: the spins): the
        # reference is that of the segment
        braked = evolving and getattr(cfg.physics, "magnetic_braking", None) is not None
        L_ref = stream.l_ref[j] if (mass_loss or braked) else L0
        angmom_error[j] = np.linalg.norm(sys_.angular_momentum(y) - L_ref) / np.linalg.norm(L_ref)

    spin_out = {}
    if sys_.dissipative:
        for q in sys_.diss:
            nd = sys_.nodes[q]
            for s, off in ((nd.left, 0), (nd.right, 1)):
                w = Y[sys_.spin_index[q] + off]
                spin_out[f"m{s + 1}"] = {
                    "P_rot_d": np.where(w > 0, 2 * np.pi / np.where(w > 0, w, 1.0), np.inf)
                    * units.DAYS_PER_YEAR,
                    "obliquity_deg": np.zeros(len(times)),
                }
    mutual = {}
    for a, b in combinations(names, 2):  # same order as the N-body backend (root.nodes())
        mutual[f"{a}-{b}"] = np.array(
            [mutual_inclination(k_hat[a][j], k_hat[b][j]) for j in range(len(times))]
        )

    if evolving:
        extremes, roche = stream.extremes, stream.roche

    configured = hdet.to_string((hdet.from_node(root),))
    res = RunResult(
        times=times,
        orbits=orbits,
        mutual=mutual,
        energy_error=energy_error,
        hierarchy=[configured] * len(times),
        configured_hierarchy=configured,
        extremes=extremes,
        stars={f"m{k + 1}": p for k, p in enumerate(props)},
        spins=spin_out,
        roche=roche,
        angmom_error=angmom_error,
        dissipative=sys_.gas is not None or sys_.dissipative,
        gas_phase=sys_.gas is not None,
    )
    if evolving:
        use(0)
    sc.finish_result(
        cfg,
        res,
        root=root,
        masses=masses,
        m0=m0,
        props=props,
        pos0=pos0,
        vel0=vel0,
        evolving=evolving,
        mass_loss=mass_loss,
        stream=stream if evolving else None,
        stop_event=stop_event,
        psi_max=psi_max,
        overflow=overflow,
    )
    if cfg.physics.encounters is not None:
        res.system["encounters_n"] = enc_counts["n"]
        res.system["encounters_close"] = enc_counts["close"]
    res.system.update(eclipse_summary(res, root, props))
    return res
