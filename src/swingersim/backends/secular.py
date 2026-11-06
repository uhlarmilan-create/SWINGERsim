"""Secular backend: orbit-averaged 2+2 dynamics (swingersim.secular.twotwo) with the same output
as the N-body backend, so that diagnostics, files and the viewer work unchanged.

The elements in the output are mean (orbit-averaged) elements; the mean anomaly is not
followed and is stored as 0. Hierarchy changes, collisions and close encounters cannot occur in
a secular model; they are reported by the N-body backend only.
"""

from __future__ import annotations

import numpy as np

from .. import massloss, masstransfer, units
from .. import masstransfer as _mt
from ..config import Config
from ..diagnostics import hierarchy as hdet
from ..diagnostics.binary import both_eclipsing_fraction
from ..diagnostics.extremes import THRESHOLDS_1ME, threshold_key
from ..diagnostics.orbits import mutual_inclination
from ..hierarchy import build_tree, initialize
from ..secular.twotwo import (
    E_SEED,
    Binary,
    SecularTwoTwo,
    elements_from_vectors,
    integrate,
    vectors_from_elements,
)
from ..stellar import roche_lobe_eggleton, star_properties
from ..tides import with_structure_lags
from . import evolving as ev
from . import secular_common as sc
from .nbody import RunResult, parent_of, spin_vector

ORDER = ("A", "B", "out")


_refine = sc.refine  # the refinement of the maxima (also used by secular_nested)


def _binaries(cfg: Config, masses, props, spins=None) -> dict[str, Binary]:
    from ..tides import k2_rotating

    binaries = {}
    w = list(spins) if spins is not None else [0.0] * 4
    for name, (i, j) in (("A", (0, 1)), ("B", (2, 3))):
        pi, pj = props[i], props[j]
        binaries[name] = Binary(
            m=(masses[i], masses[j]),
            r_au=(pi.radius_rsun * units.RSUN_AU, pj.radius_rsun * units.RSUN_AU),
            apsidal_k2=(
                k2_rotating(cfg, i, pi.apsidal_k2, w[i], masses[i], pi.radius_rsun),
                k2_rotating(cfg, j, pj.apsidal_k2, w[j], masses[j], pj.radius_rsun),
            ),
            rg2=(pi.rg2, pj.rg2),
            lag_yr=(pi.lag_s / units.SECONDS_PER_YEAR, pj.lag_s / units.SECONDS_PER_YEAR),
            tides=cfg.physics.tides,
        )
    return binaries


def _scaled_spins(sys_, y, prev, pm, m_old=None, m_new=None, cfg=None):
    """Spins after a change of the stellar structure at constant spin angular momentum (less
    the angular momentum carried by the wind when the masses change), then braked by the
    magnetised winds over the coming segment (physics.magnetic_braking)."""

    def sc(w, k):
        mo = None if m_old is None else m_old[k]
        mn = None if m_new is None else m_new[k]
        w = ev.scale_spin(w, prev[k], pm[k], mo, mn)
        if cfg is not None:
            m = mn if mn is not None else sys_.m[k]
            w = ev.brake_spin(cfg, k, w, pm[k], m, ev.SEGMENT["dt"])
        return w

    if sys_.dissipative:
        for k in range(4):
            y[20 + k] = sc(y[20 + k], k)
        return None
    s0 = sys_.spin0
    return {
        "A": (sc(s0["A"][0], 0), sc(s0["A"][1], 1)),
        "B": (sc(s0["B"][0], 2), sc(s0["B"][1], 3)),
    }


def _resync_slaved(sys_, y) -> None:
    """Backward runs (2+2): the spins held pseudo-synchronous (swingersim.tides.slaved_spin) are
    set to the pseudo-synchronous rate of their binary at the start of every segment; the spin
    angular momentum gained or lost since is exchanged with the binary's orbit (as
    secular_nested._resync_slaved)."""
    from ..tides import pseudo_sync_factor

    sys_._diss_args = sys_._dissipative_arrays()  # decides sys_.slaved for this structure
    for b, name in enumerate(("A", "B")):
        bn = sys_.binaries[name]
        M = bn.m[0] + bn.m[1]
        mu = bn.m[0] * bn.m[1] / M
        e = float(np.linalg.norm(y[6 * b + 3 : 6 * b + 6]))
        for i in (0, 1):
            if 2 * b + i not in sys_.slaved:
                continue
            C = bn.rg2[i] * bn.m[i] * bn.r_au[i] ** 2
            for _ in range(2):
                a = y[18 + b]
                w_ps = np.sqrt(units.G * M / a**3) * pseudo_sync_factor(e)
                dJ = C * (y[20 + 2 * b + i] - w_ps)
                L = mu * np.sqrt(units.G * M * a * (1.0 - e * e))
                y[18 + b] = a * ((L + dJ) / L) ** 2
                y[20 + 2 * b + i] = w_ps
    sys_._diss_args = None


def _integrate_evolving(cfg: Config, sys_, y0, masses, props, times) -> ev.Stream:
    """Evolving stars: operator splitting and streamed integration (backends/evolving.py)."""

    def apply(pm, prev, y, mm, m_prev):
        spins = _scaled_spins(sys_, y, prev, pm, m_prev, mm, cfg)
        if np.any(mm != sys_.m):
            sys_.set_masses(mm, y)
        if cfg.physics.tidal_lag != "fixed":
            v = sys_.unpack(y)
            w = list(v["sA"]) + list(v["sB"]) if spins is None else list(spins["A"] + spins["B"])
            pm[:] = with_structure_lags(cfg, pm, [(0, 1, v["aA"]), (2, 3, v["aB"])], w, mm)
        v = sys_.unpack(y)
        w = list(v["sA"]) + list(v["sB"]) if spins is None else list(spins["A"] + spins["B"])
        sys_.set_structure(_binaries(cfg, mm, pm, w), spins)
        if sys_.dissipative and getattr(sys_, "backward_span", 0.0) > 0:
            _resync_slaved(sys_, y)
            v = sys_.unpack(y)
            sys_.set_structure(_binaries(cfg, mm, pm, list(v["sA"]) + list(v["sB"])), spins)

    orbits = [
        ev.OrbitSpec(
            "A",
            3,
            18 if sys_.dissipative else None,
            sys_.a0["A"],
            (0, 1),
            a_of=lambda: sys_.a0["A"],
        ),
        ev.OrbitSpec(
            "B",
            9,
            19 if sys_.dissipative else None,
            sys_.a0["B"],
            (2, 3),
            a_of=lambda: sys_.a0["B"],
        ),
        ev.OrbitSpec("out", 15, None, sys_.a_out, None, a_of=lambda: sys_.a_out),
    ]
    return ev.stream_evolving(
        cfg,
        masses,
        props,
        times,
        y0,
        sys_.rhs,
        lambda y: sys_.energy_and_gradients(sys_.unpack(y))[0],
        apply,
        lambda: (sys_.binaries, sys_.spin0, sys_.mass_state()),
        orbits,
        THRESHOLDS_1ME,
        threshold_key,
        _refine,
        roche_lobe_eggleton,
        angmom=sys_.angular_momentum,
        scale_orbits=sys_.scale_orbits,
    )


def _extremes_of(sys_, masses, props):
    """extremes_of(sol, t, t_stop, overflow) of a 2+2 run (secular_common.dense_extremes)."""
    a0 = {"A": sys_.a0["A"], "B": sys_.a0["B"], "out": sys_.a_out}
    a_idx = {"A": 18, "B": 19} if sys_.dissipative else {}

    def a_of(Y, q):
        n = ORDER[q]
        return Y[a_idx[n]] if n in a_idx else a0[n]

    lobes = {}
    for n, (i, j) in (("A", (0, 1)), ("B", (2, 3))):
        lobes[n] = [
            (
                f"m{s + 1}",
                roche_lobe_eggleton(masses[s] / masses[o]),
                props[s].radius_rsun * units.RSUN_AU,
            )
            for s, o in ((i, j), (j, i))
        ]
    orbits = list(zip(ORDER, (3, 9, 15), strict=True))

    def extremes_of(sol, t, t_stop=None, overflow=None):
        return sc.dense_extremes(
            sol, t, orbits, a_of, lobes, 20 * (len(t) - 1) + 1, t_stop, overflow
        )

    return extremes_of


def build_secular(cfg: Config):
    if cfg.topology != "2+2":
        raise ValueError("build_secular handles 2+2 only; other hierarchies: secular_nested")
    root, masses = build_tree(cfg)
    nodes = {n.name: n for n in root.nodes()}
    props = [star_properties(cfg, f"m{k + 1}", m) for k, m in enumerate(masses)]
    # masses of the dynamics at t = 0 (with mass loss from the tracks: the current masses at
    # stellar_model.age_myr; the configured masses label the tracks)
    m_dyn = massloss.initial_masses(cfg, masses)
    if not np.array_equal(m_dyn, masses):  # periods with the masses of the dynamics
        root = build_tree(cfg, kepler_masses=m_dyn)[0]
        nodes = {n.name: n for n in root.nodes()}
    pos, vel = initialize(root, m_dyn)
    spins = [
        np.linalg.norm(spin_vector(cfg, f"m{k + 1}", parent_of(root, k), pos, vel, m_dyn))
        if cfg.physics.tides
        else 0.0
        for k in range(4)
    ]
    props = with_structure_lags(
        cfg, props, [(0, 1, nodes["A"].elements.a), (2, 3, nodes["B"].elements.a)], spins, m_dyn
    )
    binaries = _binaries(cfg, m_dyn, props, spins)
    sys_ = SecularTwoTwo(
        m_dyn,
        {n: nodes[n].elements.a for n in ORDER},
        binaries,
        {"A": (spins[0], spins[1]), "B": (spins[2], spins[3])},
        gr=cfg.physics.gr,
        octupole=cfg.run.secular_octupole,
    )
    # backward runs: spins that synchronise much faster than |t_max| held pseudo-synchronous
    sys_.backward_span = -float(cfg.run.t_max_yr) if cfg.run.t_max_yr < 0 else 0.0
    y0 = []
    for n in ORDER:
        el = nodes[n].elements
        K, p = vectors_from_elements(max(el.e, E_SEED), el.inc, el.Omega, el.omega)
        y0.extend([*K, *p])
    if sys_.dissipative:
        y0.extend([nodes["A"].elements.a, nodes["B"].elements.a, *spins])
    return sys_, np.array(y0), root, masses, props, (pos, vel)


def run_secular(cfg: Config) -> RunResult:
    """Secular run of the configuration; with physics.mergers a merger of two stars (common
    envelope, contact) is followed by a run of the merged system (swingersim.merger)."""
    if cfg.physics.mergers is not None or _mt.dynamic_onset(cfg) is not None:
        from ..merger import run_with_mergers

        return run_with_mergers(cfg, _run_secular_once)
    return _run_secular_once(cfg)


def _run_secular_once(cfg: Config) -> RunResult:
    drifting = [n for n, o in cfg.orbits.items() if o.Pdot]
    if drifting:
        raise ValueError(
            f"orbits.{drifting[0]}.Pdot: a period change of a binary needs run.backend: nbody "
            "(the secular backends average over its orbit)"
        )
    cfg = massloss.resolve_masses(cfg)
    env = cfg.physics.galactic_tide is not None or cfg.physics.encounters is not None
    if cfg.topology != "2+2" or cfg.physics.gas is not None or env or cfg.physics.gw:
        from .secular_nested import run_secular_nested

        return run_secular_nested(cfg)
    sys_, y0, root, masses, props, (pos0, vel0) = build_secular(cfg)
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
            a_init = {"A": sys_.a0["A"], "B": sys_.a0["B"], "out": sys_.a_out}
            psi_max = massloss.check_psi(cfg, root, masses, a_init, times)
        stream = _integrate_evolving(cfg, sys_, y0, masses, props, times)
        Y = stream.Y
        if stream.t_stop is not None:  # Roche-lobe overflow: the run ends there
            times, Y = sc.truncate_at_stream_stop(times, Y, stream)
            stop_event = "roche_overflow"
    else:
        times = np.linspace(0.0, cfg.run.t_max_yr, cfg.run.n_output)
        times, Y, extremes, roche, overflow = sc.integrate_until_overflow(
            lambda y, t: integrate(sys_, y, t, rtol=cfg.run.secular_rtol),
            y0,
            times,
            _extremes_of(sys_, masses, props),
            _mt.stop_on_overflow(cfg),
        )

    def use_structure(j):
        """Stellar structure at output j (evolving stars) and its energy reference."""
        if not evolving:
            return None
        b, sp, ms = stream.snapshots[j]
        sys_.restore_mass_state(ms)
        sys_.set_structure(b, sp)
        return stream.h_ref[j]

    use_structure(0)  # evolving stars: the reference values refer to the structure at t = 0
    L0 = sys_.angular_momentum(Y[:, 0])
    basis = sc.invariable_basis(L0)
    H0, _ = sys_.energy_and_gradients(sys_.unpack(y0))

    keys = ("a", "e", "inc", "Omega", "omega", "M", "inc_inv", "Omega_inv")
    orbits = {n: {k: np.zeros(len(times)) for k in keys} for n in ORDER}
    k_hat = {n: np.zeros((len(times), 3)) for n in ORDER}
    energy_error = np.zeros(len(times))
    angmom_error = np.zeros(len(times))
    spins = (
        {
            f"m{k + 1}": {
                "P_rot_d": np.full(len(times), np.nan),
                "obliquity_deg": np.zeros(len(times)),
            }
            for k in range(4)
        }
        if cfg.physics.tides
        else {}
    )
    for j in range(len(times)):
        H_ref = use_structure(j)
        v = sys_.unpack(Y[:, j])
        a_now = {"A": v["aA"], "B": v["aB"], "out": sys_.a_out}
        if mass_loss:  # the segment holds its masses; report a at the masses of time t_j
            mt = massloss.masses_at(cfg, masses, float(times[j]))
            for n, (i0, i1) in (("A", (0, 2)), ("B", (2, 4)), ("out", (0, 4))):
                a_now[n] *= sys_.m[i0:i1].sum() / mt[i0:i1].sum()
            for n, a_h in masstransfer.orbit_a_at(cfg, float(times[j])).items():  # transfer
                a_now[n] = a_h
        for n, (sk, sp) in zip(ORDER, (("KA", "pA"), ("KB", "pB"), ("K", "p"))):
            e, inc, Om, om = elements_from_vectors(v[sk], v[sp])
            o = orbits[n]
            o["a"][j], o["e"][j], o["inc"][j], o["Omega"][j], o["omega"][j] = (
                a_now[n],
                e,
                inc,
                Om,
                om,
            )
            o["M"][j] = 0.0
            k = v[sk] / np.linalg.norm(v[sk])
            k_hat[n][j] = k
            o["inc_inv"][j], o["Omega_inv"][j] = sc.invariable_angles(k, basis)
        H, _ = sys_.energy_and_gradients(v)
        energy_error[j] = abs(H / (H_ref if H_ref is not None else H0) - 1.0)
        # the wind carries angular momentum (mass loss; magnetic braking: the spins): the
        # reference is that of the segment
        braked = evolving and getattr(cfg.physics, "magnetic_braking", None) is not None
        L_ref = stream.l_ref[j] if (mass_loss or braked) else L0
        angmom_error[j] = np.linalg.norm(sys_.angular_momentum(Y[:, j]) - L_ref) / np.linalg.norm(
            L_ref
        )
        for k in range(len(spins)):
            w = (v["sA"] + v["sB"])[k]
            spins[f"m{k + 1}"]["P_rot_d"][j] = (
                2.0 * np.pi / w * units.DAYS_PER_YEAR if w > 0 else np.inf
            )

    mutual = {}
    for a, b in (("out", "A"), ("out", "B"), ("A", "B")):
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
        spins=spins,
        roche=roche,
        angmom_error=angmom_error,
        dissipative=sys_.dissipative,
    )
    if evolving:
        sys_.restore_mass_state(stream.snapshots[0][2])
        sys_.set_structure(_binaries(cfg, sys_.m, props), stream.snapshots[0][1])
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
    if evolving:
        res.system["integration_steps"] = int(stream.n_steps)
    res.system["eclipse_fraction_both"] = both_eclipsing_fraction(res, root, props)
    return res
