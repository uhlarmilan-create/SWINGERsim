"""Command line: swingersim run config.yaml [-o out.h5]"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

from .backends.nbody import run_nbody
from .config import load_config
from .io import save_run


def _mass_transfer_groups(system: dict, orbits) -> dict[str, dict]:
    """{orbit: {key: value}} of the mass_transfer_<orbit>_<key> entries of res.system; the
    orbit is the longest orbit name that fits (orbit names may contain '_')."""
    names = sorted(orbits, key=len, reverse=True)
    out: dict[str, dict] = {}
    for k, v in system.items():
        if not k.startswith("mass_transfer_"):
            continue
        rest = k.removeprefix("mass_transfer_")
        orb = next((o for o in names if rest.startswith(o + "_")), None)
        if orb is None:  # an orbit the run does not report: the first word
            orb = rest.split("_", 1)[0]
        out.setdefault(orb, {})[rest[len(orb) + 1 :]] = v
    return dict(sorted(out.items()))


def _with(cfg, **changes):
    """A copy of a configuration with command-line overrides ('run.backend': 'nbody', ...),
    validated again like the file (an override must not bypass the checks)."""
    data = cfg.model_dump()
    for path, value in changes.items():
        if value is None:
            continue
        node = data
        *head, last = path.split(".")
        for h in head:
            node = node[h]
        node[last] = value
    new = type(cfg).model_validate(data)
    if hasattr(cfg, "_base_dir"):
        new._base_dir = cfg._base_dir
    return new


def _summary(res) -> str:
    kind = "secular (orbit-averaged) backend, " if res.backend == "secular" else ""
    elements = "mean" if res.backend == "secular" else "osculating"
    lines = [
        (
            f"integrated to t = {res.times[-1]:.4g} yr "
            f"({len(res.times)} outputs, {kind}{elements} elements)"
        )
    ]
    if res.event:
        lines.append(f"EVENT: {res.event} at t = {res.event_time:.6g} yr {res.event_detail}")
    for mg in res.system.get("_mergers", []):  # physics.mergers: the run continued
        age = mg.get("age_myr")
        lines.append(
            f"MERGER ({mg['kind']}) at t = {mg['time_yr']:.6g} yr: "
            f"{' + '.join(mg['stars'])} ({mg['masses_msun'][0]:.3f} + {mg['masses_msun'][1]:.3f} "
            f"Msun) -> {mg['mass_msun']:.3f} Msun (lost {mg['mass_lost_msun']:.3f})"
            + (f", age on the new track {age:.4g} Myr" if age is not None else "")
            + f"; the run continues as {mg.get('hierarchy', '?')}"
        )
    for og in res.system.get("_onsets", []):  # physics.mass_transfer.dynamic_onset
        lines.append(
            f"MASS TRANSFER STARTED BY THE DYNAMICS at t = {og['time_yr']:.6g} yr: orbit "
            f"{og['orbit']} ({'+'.join(og['stars'])} at their lobes, e = {og.get('e_onset', 0):.4f})"
            + (
                f" circularised ({og['circularisation']}), P {og['P_before_d']:.5g} -> "
                f"{og['P_after_d']:.5g} d"
                if og.get("P_after_d") is not None
                else ""
            )
            + (f"; not continued: {og['stop']}" if og.get("stop") else "")
        )
    for orb, g in _mass_transfer_groups(res.system, res.orbits).items():
        if g.get("t_onset_yr") is None:
            lines.append(f"mass transfer in {orb}: none (detached)")
            continue
        m = ", ".join(f"{x:.3f}" for x in g["m_final"])
        lines.append(
            f"mass transfer in {orb}: donor {g['donor']} from t = {g['t_onset_yr']:.6g} yr"
            + (
                f", mass ratio reversed at {g['t_reversal_yr']:.6g} yr"
                if g["t_reversal_yr"]
                else ""
            )
            + f"; masses at the end {m} Msun, a {g['a_initial_au']:.4g} -> {g['a_final_au']:.4g} AU"
            + (f"; ends with {g['event']} at {g['t_event_yr']:.6g} yr" if g["event"] else "")
        )
        if g.get("ce_outcome") == "survive":
            lines.append(
                f"  common envelope: core {g['ce_donor_core']:.3f} Msun (radius "
                f"{g['ce_core_radius_rsun']:.3g} Rsun) and companion spiral in from "
                f"{g['ce_a_i_au']:.4g} to {g['ce_a_f_au']:.4g} AU; the envelope "
                f"({g['ce_donor_mass'] - g['ce_donor_core']:.3f} Msun) leaves the system"
            )
        elif g.get("ce_outcome") == "merge":
            lines.append(
                f"  common envelope: the cores merge into a star of {g['ce_merged_mass']:.3f} Msun"
            )
    for name, s in res.orbits.items():
        lines.append(
            f"  {name:6s} a = {np.nanmin(s['a']):.5g}..{np.nanmax(s['a']):.5g} AU"
            f"  e = {np.nanmin(s['e']):.4f}..{np.nanmax(s['e']):.4f}"
        )
    for k, s in res.mutual.items():
        d = np.degrees(s)
        lines.append(f"  i_mut {k:10s} {np.nanmin(d):7.2f}..{np.nanmax(d):7.2f} deg")
    if res.gas_phase:
        e_note = "  (gas phase: energy and angular momentum change)"
    elif getattr(res, "mass_loss", False):
        e_note = "  (mass loss: checked within each segment / per output step)"
    elif res.dissipative:
        e_note = "  (tidal dissipation: energy decreases)"
    else:
        e_note = ""
    lines.append(f"  max |dE/E| = {np.nanmax(res.energy_error):.2e}{e_note}")
    if res.angmom_error is not None:
        lines.append(f"  max |dL/L| = {np.nanmax(res.angmom_error):.2e}")
    lines.append(f"hierarchy: configured {res.configured_hierarchy}")
    if not res.transitions:
        lines.append("  no change detected")
    for tr in res.transitions:
        lines.append(
            f"  t = {tr['time_yr']:.6g} yr  {tr['kind']}: {tr['from']} -> {tr['to']} ({tr['label']})"
        )
    lines.append("orbit pairs:")
    for d in res.diagnostics:
        lines.append(
            f"  {d.inner} in {d.outer}: i_mut0 = {d.i_mut0_deg:.1f} deg, g1_0 = {d.g1_0_deg:.0f} deg"
        )
        win = "inside" if d.in_zlk_window else "outside"
        lines.append(
            f"    quadrupole theory : e = {d.quad_e_min:.3f}..{d.quad_e_max:.3f}, "
            f"i_mut = {d.quad_i_min_deg:.1f}..{d.quad_i_max_deg:.1f} deg, "
            + (f"e-cycle {d.quad_period_yr:.3g} yr " if np.isfinite(d.quad_period_yr) else "")
            + f"({win} ZLK window)"
        )
        cyc = (
            f", {d.n_zlk_cycles} ZLK cycles of ~{d.zlk_period_yr:.3g} yr" if d.n_zlk_cycles else ""
        )
        lines.append(
            f"    simulation        : e = {d.e_min:.3f}..{d.e_max_true:.3f}, "
            f"i_mut = {d.i_mut_min_deg:.1f}..{d.i_mut_max_deg:.1f} deg{cyc}"
        )
        close = (
            f", closest periastron = {d.q_min_over_radii:.2f} x (R1+R2)"
            if np.isfinite(d.q_min_over_radii)
            else ""
        )
        meas = (
            f"{d.prec_period_measured_yr:.3g} yr measured"
            if np.isfinite(d.prec_period_measured_yr)
            else "not completed in the run"
        )
        lines.append(f"    precession period : {d.t_prec_yr:.3g} yr estimated, {meas}{close}")
        warn = "  (mild hierarchy: secular estimates approximate)" if d.epsilon_sa > 0.01 else ""
        lines.append(
            f"    stability ratio Mardling-Aarseth = {d.stability_ratio:.2f}"
            f" ({'stable' if d.stability_ratio > 1 else 'UNSTABLE'}), Vynatheya+22 ="
            f" {d.stability_ratio_v22:.2f} ({'stable' if d.stability_ratio_v22 > 1 else 'UNSTABLE'}),"
            f" eps_oct = {d.epsilon_oct:.2g},"
            f" eps_SA = {d.epsilon_sa:.2g}{warn}"
        )
        gr_note = " -> GR faster than ZLK" if d.gr_suppresses_zlk else ""
        lines.append(f"    GR apsidal period ~ {d.t_gr_yr:.3g} yr{gr_note}")
        b = res.binaries.get(d.inner)
        if b and np.isfinite(b["U_tidal_yr"]):
            note = " -> faster than ZLK, ZLK suppressed" if b["U_tidal_yr"] < d.t_zlk_yr else ""
            lines.append(
                f"    tidal+rotational apsidal period {b['U_tidal_yr']:.3g} yr, "
                f"ZLK time {d.t_zlk_yr:.3g} yr{note}"
            )
    s = res.system
    quads = [q for q in s.get("quadruples", []) if "period_ratio_A_B" in q]
    if not quads and "period_ratio_A_B" in s:  # results saved by older versions
        quads = [dict(s, binaries=["A", "B"])]
    for q in quads:
        a, b = q.get("binaries", ["A", "B"])
        where = f" (quadruple {q['core']})" if len(s.get("quadruples", [])) > 1 else ""
        lines.append(
            f"binaries {a} and {b}{where}: P_{a}/P_{b} = {q['period_ratio_A_B']:.4f} "
            f"(nearest {q['nearest_period_commensurability']}, "
            f"{100 * q['period_commensurability_deviation']:+.2f} %), "
            f"precession ratio {a}/{b} = {q['precession_ratio_analytic']:.2f} estimated"
            + (
                f", {q['precession_ratio_measured']:.2f} measured"
                if np.isfinite(q["precession_ratio_measured"])
                else ""
            )
        )
        if np.isfinite(q.get("J_AB_max_analytic_deg", np.nan)) and f"{a}-{b}" in res.mutual:
            jab = np.degrees(res.mutual[f"{a}-{b}"])
            lines.append(
                f"  mutual inclination of {a} and {b}: {np.nanmin(jab):.2f}..{np.nanmax(jab):.2f}"
                f" deg, analytic limit {q['J_AB_max_analytic_deg']:.2f} deg (Vokrouhlicky 2016, "
                f"circular {a} and {b})"
            )
        if "eclipse_fraction_both" in q:
            fa = res.binaries.get(a, {}).get("eclipse_fraction", np.nan)
            fb = res.binaries.get(b, {}).get("eclipse_fraction", np.nan)
            lines.append(
                f"  eclipses (observer along z): {a} {100 * fa:.0f} %, {b} {100 * fb:.0f} %, "
                f"both at once {100 * q['eclipse_fraction_both']:.0f} % of the run"
            )
    if "eclipse_fraction_all" in s:
        lines.append(
            "eclipses (observer along z): "
            + ", ".join(
                f"{n} {100 * b.get('eclipse_fraction', np.nan):.0f} %"
                for n, b in res.binaries.items()
            )
            + f"; all at once {100 * s['eclipse_fraction_all']:.0f} % of the run"
        )
    betas = {k.removeprefix("beta_hamers_lai_"): v for k, v in s.items() if k.startswith("beta_")}
    if betas:
        tilt = {d.inner: d.i_mut0_deg for d in res.diagnostics}
        flag = any(
            0.3 <= abs(b) <= 3.0 and 40.0 <= tilt.get(k, 0.0) <= 140.0 for k, b in betas.items()
        )
        lines.append(
            "secular coupling beta (Hamers & Lai 2017): "
            + ", ".join(f"{k} {v:.3g}" for k, v in betas.items())
            + (
                "  -> beta ~ 0.3-3 with a tilt above 40 deg: strong, possibly chaotic eccentricity "
                "excitation possible"
                if flag
                else ""
            )
        )
    if res.stars:
        lines.append("stars:")
        for name, pr in res.stars.items():
            row = f"  {name}  R = {pr.radius_rsun:.3g} Rsun"
            sp = res.spins.get(name) if res.spins else None
            if sp is not None:
                p0, p1 = sp["P_rot_d"][0], sp["P_rot_d"][-1]
                rot = (
                    f", P_rot = {p0:.4g} -> {p1:.4g} d, obliquity "
                    f"{sp['obliquity_deg'][0]:.1f} -> {sp['obliquity_deg'][-1]:.1f} deg"
                    if np.isfinite(p0)
                    else ", not rotating"
                )
                row += (
                    f", k2 = {pr.apsidal_k2:.3g}, rg2 = {pr.rg2:.3g}, lag = {pr.lag_s:.3g} s{rot}"
                )
            if name in res.roche:
                row += f", Roche lobe filled up to {res.roche[name]['fill_max']:.2f}"
            if pr.teff_k is not None:
                row += f", Teff = {pr.teff_k:.0f} K"
            if pr.age_myr is not None and not res.stellar:
                row += f", {pr.phase}, age {pr.age_myr:.4g} Myr"
            lines.append(row)
            st = res.stellar.get(name) if res.stellar else None
            if st is not None and "phase" not in st and "mass_msun" in st:
                m = st["mass_msun"]
                lines.append(f"      mass {m[0]:.4g} -> {m[-1]:.4g} Msun")
            elif st is not None and "phase" in st:
                from .evolution import PHASES

                wr = st.get("wolf_rayet")

                def phases_at(j, wr=wr, ph=st["phase"]):
                    if wr is not None and len(wr) and bool(wr[j]):
                        return "Wolf-Rayet star"
                    return PHASES[int(ph[j])]

                lines.append(
                    f"      evolving: age {st['age_myr'][0]:.5g} -> {st['age_myr'][-1]:.5g} Myr, "
                    f"R {st['radius_rsun'][0]:.3g} -> {st['radius_rsun'][-1]:.3g} Rsun, "
                    f"Teff {st['teff_k'][0]:.0f} -> {st['teff_k'][-1]:.0f} K, "
                    f"{phases_at(0)} -> {phases_at(-1)}"
                )
    for name, b in res.binaries.items():
        own = b["U_expected_yr"]
        parts = [
            f"{own:.4g} yr from the binary itself"
            if np.isfinite(own)
            else "none from the binary itself"
        ]
        if np.isfinite(b["U_tidal_yr"]):
            parts.append(f"tides+rotation {b['U_tidal_yr']:.4g} yr")
        parts.append(f"GR {b['U_gr_yr']:.4g} yr" + ("" if b["gr_included"] else " (GR off)"))
        meas = b.get("U_measured_yr", float("nan"))
        lines.append(
            f"  binary {name} ({b['stars']}): apsidal period {', '.join(parts)}"
            + (
                f"; measured {meas:.4g} yr (with the other stars"
                + (", regressing)" if meas < 0 else ")")
                if np.isfinite(meas)
                else "; not measurable in the run"
            )
        )
        if "t_circ_yr" in b:
            lines.append(
                f"    tidal time scales: circularisation {b['t_circ_yr']:.3g} yr, "
                f"orbital decay {b['t_a_yr']:.3g} yr, synchronisation "
                f"{b['t_sync1_yr']:.3g} / {b['t_sync2_yr']:.3g} yr "
                f"(pseudo-synchronous P_rot = {b['P_pseudo_sync_d']:.4g} d)"
            )
    if res.megno is not None:
        state = "regular" if res.megno[-1] < 2.5 else "CHAOTIC"
        lines.append(
            f"MEGNO = {res.megno[-1]:.2f} ({state}; 2 = regular), "
            f"Lyapunov time = {res.lyapunov_time_yr:.3g} yr"
        )
    return "\n".join(lines)


def _pulsation(args) -> int:
    """Periods of the pulsating stars of a configuration at t = 0."""
    import warnings

    from .pulsation import period_table

    cfg = load_config(args.config)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rows = period_table(cfg)
    if not rows:
        print("no star has a 'pulsation' section")
        return 0
    for r in rows:
        dp = r["dP_dt_s_per_yr"]
        teff = r["Teff_K"]
        print(
            f"{r['star']}: {r['type']} ({r['mode']})  P = {r['period_d']:.6g} d"
            + (f"  dP/dt = {dp:.3g} s/yr" if dp is not None else "")
            + f"  R = {r['R_rsun']:.3g} Rsun"
            + (f"  Teff = {teff:.0f} K" if teff is not None else "")
        )
        for n in r["notes"]:
            print(f"    note: {n}")
    return 0


def _verdict_text(ok: bool | None) -> str:
    """Verdict of a stability classifier; None: no trained model for this kind of system."""
    return "no model" if ok is None else ("stable" if ok else "UNSTABLE")


def _stability(args) -> int:
    from .hierarchy import build_tree, initialize
    from .stability import ghost_test, mlp_verdicts, nested_triples

    cfg = load_config(args.config)
    root, masses = build_tree(cfg)
    pos, vel = initialize(root, masses)
    print(f"stability of {args.config} ({cfg.topology})")
    print("nested triples (ratio > 1 = stable):")
    for t in nested_triples(root, masses, pos, vel):
        print(
            f"  {t.label:22s} a_in/a_out = {t.a_in / t.a_out:.4f}, e_in = {t.e_in:.3f},"
            f" e_out = {t.e_out:.3f}, i_mut = {t.i_mut_deg:.1f} deg, q_out = {t.q_out:.3f}:"
            f" Mardling-Aarseth {t.ma01:.2f}, Vynatheya+22 {t.v22:.2f}"
        )
    verdicts = mlp_verdicts(cfg)
    if verdicts is None:
        if len(masses) >= 4:
            print("machine-learning classifier (Vynatheya+23): not installed (see the manual)")
    elif len(verdicts) == 1:
        print(f"machine-learning classifier (Vynatheya+23): {_verdict_text(verdicts[0][1])}")
    else:
        print("machine-learning classifier (Vynatheya+23), every sub-quadruple (approximation):")
        for label, ok in verdicts:
            print(f"  {label:34s} {_verdict_text(ok)}")
    if not args.no_nbody:
        t0 = time.perf_counter()
        g = ghost_test(cfg, n_outer=args.orbits)
        verdict = "stable" if g.stable else f"UNSTABLE ({g.reason})"
        print(
            f"N-body ghost test ({args.orbits:g} outer orbits, Newtonian): {verdict} after"
            f" {g.outer_orbits:.1f} outer orbits, max relative divergence {g.delta_max:.2e}"
            f"  [{time.perf_counter() - t0:.1f} s]"
        )
    return 0


def _weights_for(directory: Path, observed_file: Path | None, rows: list[dict]):
    """Observation weights from --observed or from ensemble.observed of the stored config."""
    import yaml

    from .config import ObservedValue
    from .ensemble import observation_weights

    if observed_file:
        raw = yaml.safe_load(observed_file.read_text())
    else:
        meta = directory / "ensemble.yaml"
        cfg = yaml.safe_load(meta.read_text()).get("config", {}) if meta.exists() else {}
        raw = (cfg.get("ensemble") or {}).get("observed") or {}
    if not raw:
        return None
    observed = {k: ObservedValue.model_validate(v) for k, v in raw.items()}
    for k, v in observed.items():
        parts = [] if v.value is None else [f"{v.value:g} +- {v.sigma:g}"]
        parts += [f"{s} {x:g}" for s, x in ((">=", v.min), ("<=", v.max)) if x is not None]
        print(f"observed: {k} " + ", ".join(parts))
    return observation_weights(rows, observed)


def _ensemble(args) -> int:
    from .ensemble import check_observed, load_rows, report, run_ensemble

    if args.cmd == "ensemble-report":
        if (args.x or args.y) and not args.plot:
            print("--x and --y choose the last panel of the figure: they need --plot")
            return 2
        rows = load_rows(args.directory)
        weights = None if args.prior else _weights_for(args.directory, args.observed, rows)
        print(report(rows, tuple(args.thresholds), weights))
        if args.plot:
            from .ensemble import plot_ensemble

            print(f"saved {plot_ensemble(rows, args.plot, args.x, args.y)}")
        return 0
    if args.retry_errors and not args.resume:
        print("--retry-errors runs the failed members of an existing ensemble: it needs --resume")
        return 2
    cfg = load_config(args.config)
    if cfg.ensemble is None:
        print(f"{args.config} has no 'ensemble' section")
        return 1
    cfg = _with(cfg, **{"ensemble.n": args.n, "run.backend": args.backend})
    check_observed(cfg)
    t0 = time.perf_counter()

    def progress(k, total, row):
        print(
            f"  [{k}/{total}] member {row['member']} {row['direction']}: {row['fate']}"
            f" ({row['wall_s']:.1f} s)",
            flush=True,
        )

    rows = run_ensemble(
        cfg, args.output, args.processes, args.resume, progress, retry_errors=args.retry_errors
    )
    print(report(rows))
    if cfg.ensemble.observed:
        from .ensemble import observation_weights

        try:  # the members are computed: a bad observation must not fail the job
            weights = observation_weights(rows, cfg.ensemble.observed)
        except ValueError as exc:
            print(f"\nno conditioning on ensemble.observed: {exc}")
        else:
            print("\nconditioned on ensemble.observed:")
            print(report(rows, weights=weights))
    print(f"  total wall time {time.perf_counter() - t0:.1f} s; results in {args.output}")
    return 0


def _target_scales(specs) -> dict:
    """['B_didt0_deg_yr=orbits.out.P^-2'] -> {'B_didt0_deg_yr': ('orbits.out.P', -2.0)}."""
    out = {}
    for spec in specs:
        try:
            t, rest = spec.split("=", 1)
            inp, power = rest.rsplit("^", 1)
            out[t.strip()] = (inp.strip(), float(power))
        except ValueError as err:
            raise SystemExit(f"--scale-targets {spec}: expected target=input^power") from err
    return out


def _sbi(args) -> int:
    import yaml

    from .ensemble import load_rows
    from .sbi_fit import (
        calibration,
        coverage,
        coverage_curve,
        load_model,
        plot_diagnostics,
        plot_posterior,
        summarize_posterior,
        train,
        training_fingerprint,
        training_set,
    )

    rows = load_rows(args.directory)
    observed = yaml.safe_load(args.observed.read_text())
    observed = observed.get("observed", observed)
    theta, x, cols, sig, val = training_set(rows, args.params, observed, args.log)
    print(f"{len(theta)} members pass the limits and have all values")
    fp = training_fingerprint(args.directory, len(rows), observed)
    model = None
    u = None
    if args.model and args.model.exists() and not args.retrain:
        model = load_model(args.model)
        if (model.params, model.log_params, model.observables) != (args.params, args.log, cols) or (
            not np.allclose(model.sigmas, sig)
        ):
            print("the saved model does not match (parameters, observables or sigmas): retraining")
            model = None
        elif getattr(model, "fingerprint", "") != fp:
            print(
                "the saved model was trained on another training set (ensemble, number of "
                "members or limits): retraining"
            )
            model = None
        else:
            print(f"using the trained model {args.model}")
            if args.coverage:
                print(
                    "--coverage is computed after a training: skipped for the saved model "
                    "(use --retrain)"
                )
    if model is None:
        t0 = time.perf_counter()
        n_train = len(theta)
        if args.coverage:
            rng = np.random.default_rng(0)
            perm = rng.permutation(len(theta))
            theta, x = theta[perm], x[perm]
            n_train = int(0.8 * len(theta))
        model = train(
            theta[:n_train],
            x[:n_train],
            sig,
            args.params,
            cols,
            args.log,
            max_epochs=args.epochs,
            fingerprint=fp,
        )
        print(f"trained in {time.perf_counter() - t0:.0f} s")
        if args.model:
            model.save(args.model)
            print(f"saved {args.model}")
        h = model.history
        if h.get("training_loss"):
            print(
                f"  {len(h['training_loss'])} epochs, final loss: training "
                f"{h['training_loss'][-1]:.3f}, validation {h.get('validation_loss', [float('nan')])[-1]:.3f}"
            )
        if args.coverage:
            c = coverage(model, theta[n_train:], x[n_train:], n_test=100)
            print(f"coverage of the 68 % intervals on held-out members: {c:.2f} (ideal 0.68)")
            u = calibration(model, theta[n_train:], x[n_train:], n_test=200)
            lv, cov = coverage_curve(u, [0.5, 0.68, 0.9, 0.95])
            for k, p in enumerate(model.params):
                print(
                    f"  {p}: coverage at 50/68/90/95 %: "
                    + "/".join(f"{cov[i, k] * 100:.0f}" for i in range(len(lv)))
                )
    if args.diagnostics:
        print(f"saved {plot_diagnostics(model, u, args.diagnostics)}")
    samples = model.sample(val, args.samples)
    print(summarize_posterior(model, samples, val))
    if args.save_samples:
        header = ",".join(f"log10({p})" if p in model.log_params else p for p in model.params)
        np.savetxt(args.save_samples, samples, delimiter=",", header=header, comments="")
        print(f"saved {args.save_samples}")
    if args.plot:
        print(f"saved {plot_posterior(model, samples, args.plot)}")
    return 0


def _escape_theory(args) -> int:
    from .escape_theory import (
        angular_momentum_scale,
        binary_single_integrals,
        compare_with_scatter,
        escape_outcomes,
    )

    m = args.masses
    if args.compare:
        from .ensemble import read_jsonl

        rows, _ = read_jsonl(args.compare / "rows.jsonl", null_as_nan=False)
        c = compare_with_scatter(rows, m)
        print(f"resonant binary + single outcomes: {c['n']}")
        for s in range(3):
            print(
                f"  star {s + 1} ({m[s]:g} Msun) escapes: N-body {c['observed'][s]:.3f}, "
                f"theory {c['predicted'][s]:.3f}"
            )
        print(
            f"  mean eccentricity of the binary: N-body {c['observed_mean_e']:.3f}, "
            f"theory {c['predicted_mean_e']:.3f}"
        )
        return 0
    if args.a is None:
        print("--a (binary semi-major axis) is needed")
        return 1
    E0, L0 = binary_single_integrals(m[:2], m[2], args.a, args.e, args.v_inf, args.b)
    if args.L0_ratio is not None:
        L0 = args.L0_ratio * angular_momentum_scale(m, E0)
    r = escape_outcomes(m, E0, L0, args.alpha, args.criterion)
    print(
        f"Stone & Leigh (2019), {args.criterion} criterion, alpha = {args.alpha:g}: "
        f"E0 = {E0:.5g}, L0 / L_c = {L0 / angular_momentum_scale(m, E0):.3f}"
    )
    for s in range(3):
        print(f"  star {s + 1} ({m[s]:g} Msun) escapes with probability {r.probabilities[s]:.3f}")
    print(f"  binary left behind: mean e = {r.mean_e:.3f}, mean E_B / E0 = {r.mean_u:.2f}")
    print(
        "  dP/de: " + " ".join(f"{v:.2f}" for v in r.e_density_all[::2]) + "  (e = 0.025 .. 0.925)"
    )
    return 0


def _phoebe_args(p) -> None:
    p.add_argument(
        "--engine",
        choices=["swingersim", "phoebe"],
        default="swingersim",
        help="light of the inner binaries: swingersim (tier A) or PHOEBE (tier C, needs phoebe)",
    )
    p.add_argument(
        "--phoebe-atm",
        choices=["auto", "ck2004", "blackbody"],
        default="auto",
        help="auto: ck2004, black bodies with tier-A limb darkening outside its grid (O stars)",
    )
    p.add_argument("--phoebe-distortion", choices=["roche", "sphere"], default="roche")
    p.add_argument("--phoebe-ntriangles", type=int, default=1500, help="mesh of each star")
    p.add_argument("--no-irradiation", action="store_true", help="PHOEBE without reflection")


def _phoebe_options(args) -> dict:
    return {
        "atm": args.phoebe_atm,
        "distortion": args.phoebe_distortion,
        "ntriangles": args.phoebe_ntriangles,
        "irradiation": not args.no_irradiation,
    }


def _lightcurve(args) -> int:
    import numpy as np
    import yaml

    from . import units
    from .photometry import add_noise, light_curve

    cfg = load_config(args.config)
    ld = yaml.safe_load(args.ld.read_text()) if args.ld else None
    n = max(round(args.days * 1440.0 / args.cadence_min), 2)
    t = args.start_day + np.arange(n) * args.cadence_min / 1440.0
    t0 = time.perf_counter()
    if args.engine == "phoebe":
        from .phoebe_lc import phoebe_light_curve

        lc = phoebe_light_curve(
            cfg,
            t,
            band=args.band,
            light_time=not args.no_light_time,
            ld=ld,
            exposure_min=args.exposure_min,
            n_sub=args.supersample,
            proximity=not args.no_proximity,
            **_phoebe_options(args),
        )
    else:
        lc = light_curve(
            cfg,
            t,
            band=args.band,
            ld=ld,
            exposure_min=args.exposure_min,
            n_sub=args.supersample,
            proximity=not args.no_proximity,
            light_time=not args.no_light_time,
        )
    if args.noise_ppm > 0:
        lc = add_noise(lc, args.noise_ppm, args.seed)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    lc.write_csv(args.output)
    print(
        f"light curve {args.band} ({args.engine}): {n} points, {args.days:g} d, cadence {args.cadence_min:g} min"
    )
    # share of each star out of eclipse (its brightest moment), not at the first time
    bright = {k: float(np.max(v)) for k, v in lc.components.items() if k.startswith("m")}
    total = sum(bright.values()) or 1.0
    for s in lc.stars if args.engine == "swingersim" else []:
        print(
            f"  {s.name}: R = {s.radius_au / units.RSUN_AU:.3f} Rsun, limb darkening "
            f"{', '.join(f'{x:.3f}' for x in s.ld)}, share of the light "
            f"{bright.get(s.name, 0.0) / total * 100:.1f} % (out of eclipse)"
        )
    for name, comp in lc.components.items():
        if name in lc.pulsation:  # pulsations are not eclipses
            comp = comp / lc.pulsation[name]
            f = lc.pulsation[name]
            print(f"  {name} pulsates: flux factor {f.min():.4f} .. {f.max():.4f}")
        drop = 1.0 - comp.min() / comp.max() if comp.max() > 0 else 0.0
        if drop > 1e-4:
            print(f"  {name} is eclipsed: deepest loss of its light {100 * drop:.2f} %")
    print(
        f"  minimum flux {lc.flux.min():.5f}; saved {args.output} ({time.perf_counter() - t0:.1f} s)"
    )
    if args.plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 4))
        ax.plot(lc.t_day, lc.flux, ".", ms=1.5, color="k")
        ax.set_xlabel("time [d]")
        ax.set_ylabel(f"normalised flux ({args.band})")
        fig.tight_layout()
        fig.savefig(args.plot, dpi=130)
        print(f"saved {args.plot}")
    return 0


def _print_activity(cfg) -> None:
    from .magnetic import activity_report

    for r in activity_report(cfg):
        verdict = "possible" if r["possible"] else "NOT possible (more energy than the star has)"
        print(
            f"  activity of {r['star']} (binary {r['binary']}): dE/E = {r['dE_over_E']:.3g} -> "
            f"{verdict}; largest semi-amplitude {r['max_amplitude_s']:.3g} s"
        )


def build_parser() -> argparse.ArgumentParser:
    """The command-line parser of every swingersim command (also read by swingersim gui and the
    worker, which build their forms and job command lines from it)."""
    p = argparse.ArgumentParser(prog="swingersim")
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run one simulation from a YAML config")
    r.add_argument("config", type=Path)
    r.add_argument("-o", "--output", type=Path, help="HDF5 output file")
    r.add_argument(
        "--backend", choices=("nbody", "secular"), help="override run.backend of the config"
    )
    r.add_argument(
        "--t-max-yr",
        type=float,
        help="override run.t_max_yr of the config (negative = integrate into the past)",
    )
    r.add_argument(
        "--no-viewer", action="store_true", help="do not write the 3D viewer (<output>_viewer.html)"
    )
    r.add_argument(
        "--no-plots", action="store_true", help="do not write the PNG figures (<output>_*.png)"
    )
    v = sub.add_parser("validate", help="reproduce published results (benchmarks)")
    v.add_argument("ids", nargs="*", help="benchmark ids (default: all fast ones)")
    v.add_argument("--all", action="store_true", help="include slow benchmarks (minutes)")
    st = sub.add_parser(
        "stability", help="stability criteria and the N-body ghost test of Vynatheya et al. 2023"
    )
    st.add_argument("config", type=Path)
    st.add_argument("--no-nbody", action="store_true", help="skip the N-body ghost test")
    st.add_argument("--orbits", type=float, default=100.0, help="outer orbits for the test")
    et_ = sub.add_parser(
        "emulator-train", help="train a neural emulator on a finished ensemble (F4 e)"
    )
    et_.add_argument("ensemble", type=Path, help="ensemble output folder")
    et_.add_argument("-o", "--output", type=Path, required=True, help="model file (.pkl)")
    et_.add_argument("--inputs", nargs="+", required=True, help="input columns (parameters)")
    et_.add_argument("--targets", nargs="*", default=[], help="regression columns")
    et_.add_argument(
        "--classes", nargs="*", default=[], help="binary columns or column==value (fate==intact)"
    )
    et_.add_argument("--log-inputs", nargs="*", default=[], help="inputs used as log10")
    et_.add_argument("--log-targets", nargs="*", default=[], help="targets used as log10")
    et_.add_argument(
        "--angle-inputs", nargs="*", default=[], help="inputs in degrees, used as (cos, sin)"
    )
    et_.add_argument(
        "--scale-targets",
        nargs="*",
        default=[],
        help="target=input^power: fit target / input^power (e.g. a rate times P_out^2)",
    )
    et_.add_argument("--hidden", nargs="+", type=int, default=[50, 50, 50, 50])
    et_.add_argument("--models", type=int, default=5, help="networks in the deep ensemble")
    et_.add_argument("--test-fraction", type=float, default=0.2)
    et_.add_argument("--activation", default="relu", choices=["relu", "tanh", "logistic"])
    et_.add_argument("--alpha", type=float, default=1e-4, help="L2 regularisation")
    et_.add_argument("--learning-rate", type=float, default=1e-3)
    et_.add_argument("--batch-size", default="auto")
    et_.add_argument("--max-iter", type=int, default=2000)
    et_.add_argument("--seed", type=int, default=0)
    et_.add_argument(
        "--plot", type=Path, help="PNG of the training and held-out diagnostics (matplotlib)"
    )
    er_ = sub.add_parser(
        "emulator-report",
        help="posterior of the varied parameters from the emulator and the observations",
    )
    er_.add_argument("model", type=Path)
    er_.add_argument("config", type=Path, help="ensemble configuration (vary, observed)")
    er_.add_argument("--n", type=int, default=1_000_000, help="prior draws")
    er_.add_argument("--seed", type=int, default=0)
    pl = sub.add_parser("pulsation", help="periods of the pulsating stars (stars.*.pulsation)")
    pl.add_argument("config", type=Path)
    en = sub.add_parser("ensemble", help="Monte Carlo ensemble (the 'ensemble' section)")
    en.add_argument("config", type=Path)
    en.add_argument("-o", "--output", type=Path, required=True, help="output directory")
    en.add_argument(
        "-j", "--processes", type=int, default=1, help="parallel processes (0 = all cores)"
    )
    en.add_argument("--n", type=int, help="override ensemble.n")
    en.add_argument("--backend", choices=("nbody", "secular"), help="override run.backend")
    en.add_argument("--resume", action="store_true", help="continue an interrupted ensemble")
    en.add_argument(
        "--retry-errors",
        action="store_true",
        help="with --resume: run again the members that ended with an error",
    )
    sc = sub.add_parser("scatter", help="few-body scattering experiments (binary-binary etc.)")
    sc.add_argument("config", type=Path)
    sc.add_argument("-o", "--output", type=Path, required=True, help="output directory")
    sc.add_argument(
        "-j", "--processes", type=int, default=1, help="parallel processes (0 = all cores)"
    )
    sc.add_argument("--n", type=int, help="override scatter.n")
    sc.add_argument("--resume", action="store_true", help="continue an interrupted run")
    cm = sub.add_parser("chaosmap", help="MEGNO chaos map over one or two parameters")
    cm.add_argument("config", type=Path)
    cm.add_argument("-o", "--output", type=Path, required=True, help="output directory")
    cm.add_argument("--x", required=True, help="parameter path, e.g. orbits.out.P")
    cm.add_argument("--x-range", type=float, nargs=2, required=True, metavar=("MIN", "MAX"))
    cm.add_argument("--nx", type=int, default=30)
    cm.add_argument("--log-x", action="store_true", help="logarithmic spacing in x")
    cm.add_argument("--y", help="second parameter path (optional)")
    cm.add_argument("--y-range", type=float, nargs=2, metavar=("MIN", "MAX"))
    cm.add_argument("--ny", type=int, default=30)
    cm.add_argument("--log-y", action="store_true", help="logarithmic spacing in y")
    cm.add_argument(
        "--orbits",
        type=float,
        default=1000.0,
        help="length in periods of the widest orbit (default 1000)",
    )
    cm.add_argument("--t-max", type=float, help="length in years (overrides --orbits)")
    cm.add_argument("--max-wall", type=float, default=300.0, help="seconds per grid point")
    cm.add_argument("-j", "--processes", type=int, default=1, help="0 = all cores")
    cm.add_argument("--resume", action="store_true")
    tm = sub.add_parser("timings", help="synthetic eclipse timings and O-C of the inner binaries")
    tm.add_argument("config", type=Path)
    tm.add_argument("--years", type=float, required=True, help="end of the series [yr after t = 0]")
    tm.add_argument(
        "--from-years",
        type=float,
        default=0.0,
        help="start of the series [yr, <= 0]: negative values add the eclipses before t = 0 "
        "(integrated into the past), e.g. for minima observed before the epoch of the elements",
    )
    tm.add_argument(
        "--bjd-t0",
        type=float,
        help="BJD of the model time t = 0: adds the column bjd_obs (observed BJD of each minimum)",
    )
    tm.add_argument("-o", "--output", type=Path, required=True, help="CSV file")
    tm.add_argument("--binaries", nargs="+", help="orbit names (default: all inner binaries)")
    tm.add_argument("--samples", type=int, default=40, help="samples per binary period")
    tm.add_argument("--plot", type=Path, help="PNG of the O-C diagram")
    mn = sub.add_parser(
        "minima", help="phenomenological minima (Mikulasek 2015) of the synthetic light curve"
    )
    mn.add_argument("config", type=Path)
    mn.add_argument("--years", type=float, required=True, help="length of the series [yr]")
    mn.add_argument("-o", "--output", type=Path, required=True, help="CSV file of the minima")
    mn.add_argument("--binaries", nargs="+", help="orbit names (default: all inner binaries)")
    mn.add_argument(
        "--band", default="TESS", help="passband: U, B, V, R, I, TESS, Gaia_G, Gaia_BP or Gaia_RP"
    )
    mn.add_argument("--cadence-min", type=float, default=2.0, help="sampling [min]")
    mn.add_argument("--noise-ppm", type=float, default=0.0, help="white noise [ppm]")
    mn.add_argument("--seed", type=int, default=1)
    mn.add_argument("--window", type=float, default=2.5, help="window / eclipse duration")
    mn.add_argument("--season-days", type=float, help="shape per season of this length [d]")
    mn.add_argument("--seasons", type=Path, help="CSV of the per-season shapes")
    mn.add_argument("--samples", type=int, default=40, help="samples per binary period")
    mn.add_argument("--no-proximity", action="store_true", help="no ellipsoidal/beaming/reflection")
    mn.add_argument("--no-light-time", action="store_true", help="ignore the light-travel time")
    _phoebe_args(mn)
    mn.add_argument(
        "--blended",
        action="store_true",
        help="measure on the total light (default: each binary on its own light + third light)",
    )
    mn.add_argument("--plot", type=Path, help="PNG of the O-C of minima and conjunctions")
    mn.add_argument(
        "--bjd-t0",
        type=float,
        help="BJD_TDB of t = 0 for the time column of the CSV (default fit.stage2.bjd_t0, "
        "else J2000)",
    )
    lc = sub.add_parser("lightcurve", help="synthetic light curve from the N-body positions (F4)")
    lc.add_argument("config", type=Path)
    lc.add_argument("--days", type=float, required=True, help="length of the series [d]")
    lc.add_argument("--start-day", type=float, default=0.0, help="first time [d after t = 0]")
    lc.add_argument("--cadence-min", type=float, default=2.0, help="sampling [min]")
    lc.add_argument(
        "--band", default="TESS", help="passband: U, B, V, R, I, TESS, Gaia_G, Gaia_BP or Gaia_RP"
    )
    lc.add_argument("--exposure-min", type=float, default=0.0, help="exposure for smearing [min]")
    lc.add_argument("--supersample", type=int, default=1, help="sub-exposures per point")
    lc.add_argument("--noise-ppm", type=float, default=0.0, help="white noise [ppm]")
    lc.add_argument("--seed", type=int, default=1)
    lc.add_argument("--ld", type=Path, help="YAML: star -> limb-darkening coefficients a1..a4")
    lc.add_argument("--no-proximity", action="store_true", help="no ellipsoidal/beaming/reflection")
    lc.add_argument("--no-light-time", action="store_true", help="ignore the light-travel time")
    _phoebe_args(lc)
    lc.add_argument("-o", "--output", type=Path, required=True, help="CSV file")
    lc.add_argument("--plot", type=Path, help="PNG of the light curve")
    sb = sub.add_parser("sbi", help="simulation-based inference (neural posterior) on an ensemble")
    sb.add_argument("directory", type=Path, help="finished ensemble folder (training set)")
    sb.add_argument(
        "--observed",
        type=Path,
        required=True,
        help="YAML: column -> {value, sigma} and/or {min, max}",
    )
    sb.add_argument("--params", nargs="+", required=True, help="parameters to infer")
    sb.add_argument("--log", nargs="*", default=[], help="parameters inferred in log10")
    sb.add_argument("--model", type=Path, help="trained model file (reused if it matches)")
    sb.add_argument("--retrain", action="store_true", help="train even if --model exists")
    sb.add_argument("--epochs", type=int, default=2000, help="maximum training epochs")
    sb.add_argument("-n", "--samples", type=int, default=20000, help="posterior samples")
    sb.add_argument("--save-samples", type=Path, help="CSV of the posterior samples")
    sb.add_argument("--plot", type=Path, help="PNG corner plot")
    sb.add_argument(
        "--coverage",
        action="store_true",
        help="calibration test on 20 %% of the members held out from training",
    )
    sb.add_argument(
        "--diagnostics",
        type=Path,
        help="PNG: loss per epoch and, with --coverage, the coverage and calibration plots",
    )
    et = sub.add_parser(
        "escape-theory", help="statistical break-up of a three-body system (Stone & Leigh 2019)"
    )
    et.add_argument(
        "--masses",
        type=float,
        nargs=3,
        required=True,
        help="binary star 1, binary star 2, single star [Msun]",
    )
    et.add_argument("--a", type=float, help="binary semi-major axis [AU]")
    et.add_argument("--e", type=float, default=0.0, help="binary eccentricity")
    et.add_argument("--v-inf", type=float, default=0.0, help="velocity at infinity [km/s]")
    et.add_argument("--b", type=float, default=0.0, help="impact parameter [AU]")
    et.add_argument(
        "--L0-ratio", type=float, help="total angular momentum in units of L_c(E0) (instead of --b)"
    )
    et.add_argument("--alpha", type=float, default=2.0)
    et.add_argument("--criterion", choices=("AE", "SE"), default="AE")
    et.add_argument(
        "--compare",
        type=Path,
        help="a swingersim scatter folder of binary-single runs to compare with",
    )
    er = sub.add_parser("ensemble-report", help="statistics of a finished ensemble")
    er.add_argument("directory", type=Path)
    er.add_argument(
        "--thresholds",
        type=float,
        nargs="+",
        default=[1e-1, 1e-2, 1e-3],
        help="report the fraction of members with 1 - e_max below these values",
    )
    er.add_argument(
        "--observed",
        type=Path,
        help="YAML file: summary column -> {value, sigma}; overrides ensemble.observed",
    )
    er.add_argument(
        "--prior", action="store_true", help="ignore the observations (unconditioned report)"
    )
    er.add_argument("--plot", type=Path, help="also save an overview figure (needs matplotlib)")
    er.add_argument("--x", help="drawn parameter for the last panel, e.g. orbits.A.i")
    er.add_argument("--y", help="summary column for the last panel, e.g. B_didt0_deg_yr")
    wk = sub.add_parser(
        "worker", help="run jobs dropped as YAML files into a folder (e.g. on a lab computer)"
    )
    wk.add_argument("jobdir", type=Path, nargs="?", help="job folder (default: <root>/jobs)")
    wk.add_argument("--root", type=Path, default=Path("."), help="SWINGERsim folder (default: .)")
    wk.add_argument("--poll", type=float, default=10.0, help="seconds between checks")
    wk.add_argument("--once", action="store_true", help="run the waiting jobs and stop")
    wk.add_argument(
        "--cores",
        type=int,
        default=None,
        help="cores the jobs may use together (default: all); one-core jobs run side by side",
    )
    mg = sub.add_parser("mesa-grid", help="compute an own stellar-model grid with MESA (F4)")
    mg.add_argument("config", type=Path, help="grid file (YAML), e.g. configs/mesa/*.yaml")
    mg.add_argument("--root", type=Path, default=Path("."), help="SWINGERsim folder (default: .)")
    mg.add_argument("-j", "--processes", type=int, default=0, help="parallel MESA runs")
    mg.add_argument(
        "--test", action="store_true", help="only one 1 Msun track (checks the installation)"
    )
    mr = sub.add_parser("mesa-reduce", help="reduce MESA profiles to SWINGERsim structure grids")
    mr.add_argument("config", type=Path, help="grid file (YAML)")
    mr.add_argument("--root", type=Path, default=Path("."), help="SWINGERsim folder (default: .)")
    mr.add_argument("-j", "--processes", type=int, default=0, help="parallel processes")
    w = sub.add_parser(
        "view",
        help="interactive HTML viewer of saved runs (needs plotly): one run, several, or every "
        "run in folders",
    )
    w.add_argument(
        "h5",
        type=Path,
        nargs="+",
        help="HDF5 files written by 'swingersim run -o', or folders searched for them (e.g. output)",
    )
    w.add_argument(
        "-o",
        "--output",
        type=Path,
        help="HTML file of a single run (default: <run>_viewer.html next to the HDF5)",
    )
    w.add_argument(
        "--force",
        action="store_true",
        help="folders: rebuild every viewer (default: only missing ones and those older than "
        "their run)",
    )
    w.add_argument("--frames", type=int, default=150, help="number of animation frames")
    w.add_argument(
        "--cdn", action="store_true", help="load plotly.js from the internet (smaller file)"
    )
    ap = sub.add_parser(
        "applegate",
        help="can magnetic activity (Applegate 1992) explain an eclipse-timing variation?",
    )
    ap.add_argument("--mass", type=float, required=True, help="active star [Msun]")
    ap.add_argument("--radius", type=float, required=True, help="active star [Rsun]")
    ap.add_argument("--teff", type=float, required=True, help="active star [K]")
    ap.add_argument("--a", type=float, required=True, help="binary separation [Rsun]")
    ap.add_argument("--pbin", type=float, required=True, help="binary period [d]")
    ap.add_argument("--pmod", type=float, required=True, help="modulation period [yr]")
    ap.add_argument("--oc", type=float, help="observed O - C semi-amplitude [s]")
    ft = sub.add_parser(
        "fit",
        help="fit observed times of minima (stage 1: LiTE, apsidal motion, quadratic; "
        "stage 2: joint N-body model of all pairs)",
    )
    ft.add_argument("config", type=Path, help="configuration with a 'fit' section")
    ft.add_argument("-o", "--output", type=Path, help="output folder (default: fit.output.dir)")
    ft.add_argument(
        "--stage",
        type=int,
        choices=[1, 2],
        help="1: analytic O - C (default without fit.stage2); 2: joint N-body model "
        "(default when the configuration has fit.stage2)",
    )
    ft.add_argument(
        "--processes",
        type=int,
        help="stage 2: cores for the Jacobian (0 = all; default fit.stage2.processes)",
    )
    ft.add_argument(
        "--sampler",
        choices=["none", "emcee"],
        help="override the sampler of the stage that runs (fit.stage1.sampler or "
        "fit.stage2.sampler)",
    )
    ft.add_argument(
        "--anchor",
        action="store_true",
        help="only anchor the period and phase of every eclipsing pair to the minima and write "
        "anchored.yaml (start values for stage 2 from published elements)",
    )
    ft.add_argument(
        "--bjd-t0",
        type=float,
        help="--anchor: BJD_TDB of t = 0 of the configuration (default fit.stage2.bjd_t0)",
    )
    ft.add_argument("--no-plots", action="store_true", help="do not draw the figures")
    rvp = sub.add_parser(
        "rv",
        help="radial velocities of the stars from the N-body model (curves; with measured "
        "velocities for comparison)",
    )
    rvp.add_argument("config", type=Path, help="configuration")
    rvp.add_argument("--bjd-t0", type=float, help="BJD_TDB of t = 0 (default fit.stage2.bjd_t0)")
    rvp.add_argument("--start", type=float, required=True, help="first time [BJD_TDB]")
    rvp.add_argument("--end", type=float, required=True, help="last time [BJD_TDB]")
    rvp.add_argument("--step", type=float, default=0.01, help="step of the curves [d]")
    rvp.add_argument("--stars", nargs="+", help="stars (default: all), e.g. m1 m2")
    rvp.add_argument("--gamma", type=float, help="systemic velocity [km/s] (default fit.rv.gamma)")
    rvp.add_argument("--measured", type=Path, help="CSV of measured velocities (default fit.rv)")
    rvp.add_argument("-o", "--output", type=Path, help="output folder (default output/rv)")
    ap_ = sub.add_parser(
        "apsidal",
        help="apsidal motion of the eclipsing pairs split into tidal, rotational, relativistic "
        "(and dynamical) parts; k2 factor needed for a target apsidal period",
    )
    ap_.add_argument("config", type=Path, help="configuration")
    ap_.add_argument(
        "--target",
        nargs="+",
        default=[],
        help="target apsidal periods [yr] as pair=value, e.g. A=194 B=98 C=214",
    )
    ap_.add_argument(
        "--nbody-years",
        type=float,
        help="also run the N-body model for so many years and measure the total advance "
        "(adds the perturbations by the other stars)",
    )
    ap_.add_argument("-o", "--output", type=Path, help="CSV file of the table")
    g = sub.add_parser("gui", help="web application for inputs and outputs (opens a browser)")
    g.add_argument(
        "--root",
        type=Path,
        default=None,
        help="SWINGERsim folder (default: the current folder or the nearest folder above it "
        "that has configs/)",
    )
    g.add_argument("--host", default="127.0.0.1", help="0.0.0.0 = reachable from the network")
    g.add_argument("--port", type=int, default=8765)
    g.add_argument("--no-browser", action="store_true", help="do not open a browser")
    from .report import NO_REPORT

    for name, sp in sub.choices.items():
        if name in NO_REPORT:
            continue
        sp.add_argument(
            "--html",
            type=Path,
            help="HTML report of the results (one file, opens offline in a browser; default "
            "next to the output)",
        )
        sp.add_argument("--no-html", action="store_true", help="do not write the HTML report")
    return p


def find_root(start: Path) -> Path:
    """The SWINGERsim folder for the GUI: ``start`` or the nearest folder above it with configs/
    (started from jobs/inbox/, the GUI must not create a second jobs/ tree there)."""
    start = start.resolve()
    for d in (start, *start.parents):
        if (d / "configs").is_dir():
            if d != start:
                print(f"swingersim gui: using the SWINGERsim folder {d}")
            return d
    return start


def command_catalog() -> dict:
    """Every command with its arguments, for forms: {name: {help, arguments: [...]}}."""
    p = build_parser()
    sub = next(a for a in p._actions if isinstance(a, argparse._SubParsersAction))
    helps = {a.dest: a.help for a in sub._choices_actions}
    out = {}
    for name, sp in sub.choices.items():
        args = []
        for a in sp._actions:
            if isinstance(a, argparse._HelpAction):
                continue
            if isinstance(a, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
                kind = "flag"
            elif a.type is Path:
                kind = "path"
            elif a.type is int:
                kind = "int"
            elif a.type is float:
                kind = "float"
            else:
                kind = "str"
            default = a.default if not isinstance(a.default, Path) else str(a.default)
            args.append(
                {
                    "dest": a.dest,
                    "flags": list(a.option_strings),
                    "positional": not a.option_strings,
                    "kind": kind,
                    "nargs": a.nargs if not isinstance(a.nargs, int) else str(a.nargs),
                    "choices": list(a.choices) if a.choices else None,
                    "default": default if default is not argparse.SUPPRESS else None,
                    "required": bool(a.required)
                    or (not a.option_strings and a.nargs not in ("*", "?")),
                    "help": (a.help or "").replace("%%", "%"),
                }
            )
        out[name] = {"help": helps.get(name, ""), "arguments": args}
    return out


OUTPUT_FILE_ARGS = ("output", "plot", "seasons", "diagnostics")


def _make_output_folders(args) -> None:
    """Create the folders of the files a command will write (-o, --plot, --seasons,
    --diagnostics), so that e.g. 'timings -o output/x/new/t.csv' works in a new folder."""
    if args.cmd in ("gui", "worker", "mesa-grid", "mesa-reduce"):
        return
    for name in OUTPUT_FILE_ARGS:
        value = getattr(args, name, None)
        if isinstance(value, (str, Path)) and str(value):
            Path(value).parent.mkdir(parents=True, exist_ok=True)


REPORT_CONTEXT: dict = {}  # filled by a command for its HTML report (e.g. the outcomes)


def main(argv: list[str] | None = None) -> int:
    from . import report

    p = build_parser()
    argv = list(sys.argv[1:] if argv is None else argv)
    args = p.parse_args(argv)
    args.argv_rest = argv[1:]  # the command line after the command (for reports)
    _make_output_folders(args)
    if args.cmd in report.NO_REPORT or getattr(args, "no_html", False):
        return _main(args)
    REPORT_CONTEXT.clear()
    tee = report.capture()

    def write_report() -> None:
        path = REPORT_CONTEXT.get("path") or report.default_path(args)
        if path is None:
            return
        try:
            out = report.write(args, argv, tee.getvalue(), path, REPORT_CONTEXT)
            print(f"report {out}  (open it in a web browser)")
        except Exception as exc:  # noqa: BLE001 - the results are saved; the report is extra
            print(f"(no HTML report: {type(exc).__name__}: {exc})")

    try:
        rc = _main(args)
    except BaseException:
        report.release(tee)
        if args.cmd == "validate" and REPORT_CONTEXT.get("outcomes"):
            write_report()  # hours of benchmarks are not lost to an interruption
        raise
    report.release(tee)
    if rc == 0 or (args.cmd == "validate" and "outcomes" in REPORT_CONTEXT):
        write_report()
    return rc


def _main(args) -> int:
    if args.cmd == "applegate":
        from . import magnetic as mg

        kmax = mg.applegate_max_oc(args.mass, args.radius, args.teff, args.a, args.pbin, args.pmod)
        print(
            f"largest O - C semi-amplitude the star can power over {args.pmod:g} yr: {kmax:.4g} s "
            f"({kmax / 86400:.3g} d; thin-shell model of Tian et al. 2009)"
        )
        if args.oc is not None:
            dpp = mg.dp_over_p_from_oc(args.oc, args.pmod)
            e = mg.applegate_energy_ratio(
                args.mass, args.radius, args.teff, args.a, args.pbin, dpp, args.pmod
            )
            verdict = "possible" if e < 1 else "not possible (needs more energy than the star has)"
            print(f"observed {args.oc:g} s: dP/P = {dpp:.3g}, dE/E_sec = {e:.3g} -> {verdict}")
        return 0
    if args.cmd == "fit":
        from . import ocfit
        from .config import load_config as _load

        cfg = _load(args.config)
        if cfg.fit is None:
            print(f"{args.config}: no 'fit' section")
            return 2
        stage = args.stage or (2 if cfg.fit.stage2 is not None else 1)
        out = args.output or Path(cfg.base_dir or ".") / cfg.fit.output.dir
        if args.bjd_t0 is not None and not args.anchor:
            print("--bjd-t0 is the epoch for --anchor (stage 2 uses fit.stage2.bjd_t0)")
            return 2
        # --processes 0 (all cores) is the default the GUI adds to every fit job: harmless
        explicit_processes = args.processes not in (None, 0)
        if args.anchor and (args.sampler or explicit_processes or args.stage):
            print(
                "--anchor only anchors the pairs: --sampler, --processes and --stage do not apply"
            )
            return 2
        if stage == 1 and explicit_processes:
            print("--processes applies to stage 2 (the N-body Jacobian); stage 1 runs on one core")
            return 2
        if args.anchor:
            from . import dynfit

            a, t0, out = dynfit.run_anchor(args.config, args.output, args.bjd_t0)
            print(dynfit.anchor_report(a, t0))
            print(f"anchored configuration: {Path(out) / 'anchored.yaml'}")
            return 0
        if stage == 2:
            from . import dynfit

            if cfg.fit.stage2 is None:
                print(f"{args.config}: --stage 2 needs a 'fit.stage2' section")
                return 2
            import inspect

            kw = {}
            if args.sampler:
                if "sampler" not in inspect.signature(dynfit.run_stage2).parameters:
                    print("--sampler for stage 2: set fit.stage2.sampler in the configuration")
                    return 2
                kw["sampler"] = args.sampler
            res2 = dynfit.run_stage2(
                args.config, args.output, plots=not args.no_plots, processes=args.processes, **kw
            )
            print(dynfit.report(res2))
        else:
            res = ocfit.run_fit(cfg, args.output, plots=not args.no_plots, sampler=args.sampler)
            print(ocfit.report(res))
        print(f"results in {out}")
        REPORT_CONTEXT["path"] = getattr(args, "html", None) or Path(out) / "report.html"
        REPORT_CONTEXT["out_dir"] = Path(out)
        return 0
    if args.cmd == "apsidal":
        from . import apsidal
        from .config import load_config as _load

        targets = {}
        for t in args.target:
            k, _, v = t.partition("=")
            targets[k] = float(v)
        rows = apsidal.apsidal_table(_load(args.config), targets, args.nbody_years)
        print(apsidal.report(rows))
        if args.output:
            apsidal.write_csv(rows, args.output)
            print(f"table in {args.output}")
        return 0
    if args.cmd == "rv":
        from . import radvel
        from .config import load_config as _load

        cfg = _load(args.config)
        st = cfg.fit.stage2 if cfg.fit is not None else None
        t0 = args.bjd_t0 if args.bjd_t0 is not None else (st.bjd_t0 if st else None)
        if t0 is None:
            print("swingersim rv needs --bjd-t0 (or fit.stage2.bjd_t0)")
            return 2
        gamma = args.gamma
        if gamma is None:
            gamma = cfg.fit.rv.gamma if cfg.fit is not None and cfg.fit.rv is not None else 0.0
        meas = None
        if args.measured is not None:
            meas = radvel.read_rv(args.measured)
        elif cfg.fit is not None and cfg.fit.rv is not None:
            from .dynfit import load_rv

            meas = load_rv(cfg)
        t, curves = radvel.rv_curve(cfg, t0, args.start, args.end, args.step, args.stars)
        out = args.output or Path(cfg.base_dir or ".") / "output" / "rv"
        out.mkdir(parents=True, exist_ok=True)
        radvel.write_curve(out / "rv_model.csv", t, curves, gamma)
        radvel.plot_curve(out / "rv_model.png", t, curves, gamma, meas)
        print(f"radial velocities of {', '.join(curves)} in {out}")
        return 0
    if args.cmd == "gui":
        from .gui.server import serve

        serve(args.root or find_root(Path.cwd()), args.host, args.port, not args.no_browser)
        return 0
    if args.cmd != "worker":
        from . import cnoise

        cnoise.install()  # drop the REBOUNDx gr_full "Fractional Error" flood of stderr

    if args.cmd == "view":
        try:
            from .viewer import write_html
        except ImportError:
            print("the viewer needs plotly: pip install -e '.[viz]'")
            return 1
        from .runproducts import product_paths
        from .viewer import run_files

        if args.output is not None:
            if len(args.h5) != 1 or args.h5[0].is_dir():
                print("-o needs exactly one run file")
                return 2
            out = write_html(args.h5[0], args.output, args.frames, not args.cdn)
            print(f"saved {out}  (open it in a web browser)")
            return 0
        n_ok = n_skip = n_fail = 0
        for path in args.h5:
            if not path.exists():
                print(f"{path}: not found")
                n_fail += 1
                continue
            folder = path.is_dir()
            for h5 in run_files(path):
                # next to the run as <run>_viewer.html (<run>.html is the run's report)
                html = product_paths(h5)["viewer"]
                fresh = html.is_file() and html.stat().st_mtime >= h5.stat().st_mtime
                if folder and fresh and not args.force:
                    n_skip += 1
                    continue
                try:
                    write_html(h5, html, args.frames, not args.cdn)
                except Exception as exc:  # noqa: BLE001 - one broken file does not stop the rest
                    print(f"{h5}: no viewer ({type(exc).__name__}: {exc})")
                    n_fail += 1
                    continue
                print(f"saved {html}")
                n_ok += 1
        if n_ok + n_skip + n_fail > 1:
            print(f"{n_ok} viewers written, {n_skip} up to date, {n_fail} failed")
        return 1 if n_fail and not n_ok else 0

    if args.cmd == "stability":
        return _stability(args)

    if args.cmd == "emulator-train":
        from . import emulator

        bs = args.batch_size if args.batch_size == "auto" else int(args.batch_size)
        t0 = time.perf_counter()
        em = emulator.train_from_ensemble(
            args.ensemble,
            args.inputs,
            args.targets,
            args.classes,
            log_inputs=args.log_inputs,
            log_targets=args.log_targets,
            angle_inputs=args.angle_inputs,
            target_scales=_target_scales(args.scale_targets),
            hidden=tuple(args.hidden),
            n_models=args.models,
            test_fraction=args.test_fraction,
            seed=args.seed,
            max_iter=args.max_iter,
            activation=args.activation,
            alpha=args.alpha,
            learning_rate=args.learning_rate,
            batch_size=bs,
        )
        path = em.save(args.output)
        m = em.metrics
        print(f"emulator: {m['n_train']} training, {m['n_test']} test members")
        for k, v in m.items():
            if ":" in k:
                print(f"  {k:40s} {v:.4g}" if isinstance(v, float) else f"  {k:40s} {v}")
        for name, nets in [("regression", em.regressors)] + list(em.classifiers.items()):
            if nets:
                its = [m.n_iter_ for m in nets]
                sc = [m.best_validation_score_ for m in nets]
                print(
                    f"  {name}: epochs {min(its)}-{max(its)} (early stopping), best validation "
                    f"score {min(sc):.3f}-{max(sc):.3f}"
                )
        print(f"saved {path} ({time.perf_counter() - t0:.1f} s)")
        if args.plot:
            print(f"saved {emulator.diagnostic_plot(em, args.plot)}")
        return 0

    if args.cmd == "emulator-report":
        from . import emulator

        em = emulator.load(args.model)
        cfg = load_config(args.config)
        print(emulator.emulated_report(em, cfg, args.n, seed=args.seed))
        return 0

    if args.cmd == "pulsation":
        return _pulsation(args)

    if args.cmd in ("mesa-grid", "mesa-reduce"):
        from .mesa import MesaGridSpec, reduce_grid, run_grid

        spec = MesaGridSpec.load(args.config)
        if args.cmd == "mesa-grid":
            res = run_grid(spec, args.root, args.processes, args.test)
            failed = [k for k, v in res.items() if v]
            print(f"{len(res) - len(failed)} tracks ok, {len(failed)} failed")
            return 1 if failed else 0
        out = reduce_grid(spec, args.root, args.processes)
        print(f"saved {out}")
        return 0

    if args.cmd == "worker":
        from .worker import Worker, WorkerLocked

        root = args.root.resolve()
        jobdir = args.jobdir or root / "jobs"
        try:
            worker = Worker(jobdir, root, args.poll, args.cores)
        except WorkerLocked as exc:
            print(f"swingersim worker: {exc}", flush=True)
            return 1
        print(
            f"swingersim worker: jobs in {worker.jobdir}, root {root} (Ctrl+C or an empty file "
            f"{worker.jobdir / 'stop'} stops it; running jobs go back to the inbox); "
            f"{worker.cores} cores",
            flush=True,
        )
        try:
            worker.loop(once=args.once)
        except KeyboardInterrupt as exc:
            worker.heartbeat({"stopped": True})
            print(
                f"swingersim worker stopped ({str(exc) or 'Ctrl+C'}); start it again to continue",
                flush=True,
            )
        return 0

    if args.cmd in ("ensemble", "ensemble-report"):
        return _ensemble(args)

    if args.cmd == "sbi":
        return _sbi(args)

    if args.cmd == "escape-theory":
        return _escape_theory(args)

    if args.cmd == "lightcurve":
        return _lightcurve(args)

    if args.cmd == "minima":
        from . import phenom

        cfg = load_config(args.config)
        t0 = time.perf_counter()
        res = phenom.synthetic_minima(
            cfg,
            args.years,
            args.binaries,
            band=args.band,
            cadence_min=args.cadence_min,
            noise_ppm=args.noise_ppm,
            seed=args.seed,
            window=args.window,
            season_days=args.season_days,
            samples_per_period=args.samples,
            proximity=not args.no_proximity,
            light_time=not args.no_light_time,
            isolate=not args.blended,
            engine=args.engine,
            phoebe_options=_phoebe_options(args) if args.engine == "phoebe" else None,
        )
        phenom.write_minima_csv(res, args.output, config=str(args.config), bjd_t0=args.bjd_t0)
        print(f"phenomenological minima ({args.band}, cadence {args.cadence_min:g} min):")
        print(phenom.summarize(res))
        _print_activity(cfg)
        print(
            f"  {len(res.minima)} minima saved in {args.output} ({time.perf_counter() - t0:.1f} s)"
        )
        if args.seasons:
            phenom.write_seasons_csv(res, args.seasons)
            print(f"saved {args.seasons}")
        if args.plot:
            phenom.plot_minima(res, args.plot)
            print(f"saved {args.plot}")
        return 0

    if args.cmd == "timings":
        from .timing import eclipse_times, o_minus_c, plot_oc, summarize_timings, write_csv

        cfg = load_config(args.config)
        t0 = time.perf_counter()
        if args.from_years > 0:
            print("--from-years must be <= 0 (the start before t = 0)")
            return 2
        ecl, p_nom = eclipse_times(cfg, args.years, args.binaries, args.samples)
        if args.from_years < 0:
            past, _ = eclipse_times(cfg, args.from_years, args.binaries, args.samples)
            ecl = past + ecl
        rows = o_minus_c(ecl, p_nom)
        write_csv(rows, args.output, bjd_t0=args.bjd_t0)
        print(summarize_timings(rows, args.years))
        _print_activity(cfg)
        print(f"  {len(rows)} eclipses saved in {args.output} ({time.perf_counter() - t0:.1f} s)")
        if args.plot:
            plot_oc(rows, args.plot)
            print(f"saved {args.plot}")
        return 0

    if args.cmd == "chaosmap":
        import json as _json

        from .chaosmap import chaos_map, grid_values, plot_map, summarize_map

        cfg = load_config(args.config)
        if args.y and not args.y_range:
            print("--y needs --y-range")
            return 1
        xs = grid_values(*args.x_range, args.nx, args.log_x)
        ys = grid_values(*args.y_range, args.ny, args.log_y) if args.y else None
        rows = chaos_map(
            cfg,
            args.output,
            args.x,
            xs,
            args.y,
            ys,
            args.orbits,
            args.t_max,
            args.processes,
            args.max_wall,
            resume=args.resume,
            progress=lambda m: print(m, flush=True),
        )
        meta = _json.loads((args.output / "chaosmap.json").read_text())
        text = summarize_map(meta, rows)
        (args.output / "summary.txt").write_text(text + "\n")
        print(text)
        try:
            print(f"saved {plot_map(meta, rows, args.output / 'map.png')}")
        except ImportError:
            print("(matplotlib not installed: no map.png)")
        return 0

    if args.cmd == "scatter":
        from .scatter import load_scatter_config, run_scatter, summarize_scatter

        cfg = _with(load_scatter_config(args.config), **{"scatter.n": args.n})
        t0 = time.perf_counter()
        rows = run_scatter(
            cfg, args.output, args.processes, args.resume, lambda m: print(m, flush=True)
        )
        text = summarize_scatter(cfg, rows)
        (args.output / "summary.txt").write_text(text + "\n")
        print(text)
        print(f"  total wall time {time.perf_counter() - t0:.1f} s; results in {args.output}")
        return 0

    if args.cmd == "validate":
        from .validation import BENCHMARKS, format_report, run_benchmarks

        known = {b.id for b in BENCHMARKS}
        unknown = [i for i in args.ids if i not in known]
        if unknown:
            print(
                f"unknown benchmark ids: {', '.join(unknown)} (swingersim validate lists them all)"
            )
            return 2
        if args.ids and not getattr(args, "html", None):  # a partial run keeps the full report
            tag = "-".join(args.ids) if len("-".join(args.ids)) <= 80 else f"{len(args.ids)}-ids"
            REPORT_CONTEXT["path"] = Path("output") / "validation" / f"report-{tag}.html"
        outcomes: list = []
        REPORT_CONTEXT["outcomes"] = outcomes  # filled as the benchmarks finish
        run_benchmarks(args.ids or None, include_slow=args.all, out=outcomes)
        print(format_report(outcomes, command="swingersim validate " + " ".join(args.argv_rest)))
        for o in outcomes:
            print(f"  {o.benchmark.id}: {o.benchmark.source} ({o.seconds:.1f} s)")
        failed = [
            o.benchmark.id
            for o in outcomes
            if o.error or (not o.passed and not o.benchmark.explanation and not o.skipped)
        ]
        for o in outcomes:
            if o.error:
                print(f"  ERROR {o.benchmark.id}: {o.error}")
        open_ = [
            o
            for o in outcomes
            if not o.passed and o.benchmark.explanation and not o.skipped and not o.error
        ]
        for o in open_:
            print(f"  EXPLAINED {o.benchmark.id}: {o.benchmark.explanation}")
        for o in outcomes:
            if o.skipped:
                print(f"  SKIPPED {o.benchmark.id}: {o.skipped}")
        print("ALL PASSED" if not failed else f"FAILED: {', '.join(failed)}")
        return 1 if failed else 0

    cfg = _with(
        load_config(args.config), **{"run.backend": args.backend, "run.t_max_yr": args.t_max_yr}
    )
    t0 = time.perf_counter()
    if cfg.run.backend == "secular":
        from .backends.secular import run_secular

        res = run_secular(cfg)
    else:
        from .backends import nbody as nb

        partial = Path(str(args.output) + ".partial.json") if args.output else None
        nb.PARTIAL["path"] = partial
        try:
            res = run_nbody(cfg)
        finally:
            nb.PARTIAL["path"] = None
            if partial is not None:
                partial.unlink(missing_ok=True)
    text = _summary(res) + f"\n  wall time {time.perf_counter() - t0:.1f} s"
    print(text)
    if args.output:
        from .runproducts import write_products

        save_run(args.output, cfg, res)
        print(f"saved {args.output}")
        files = write_products(
            args.output, text, viewer=not args.no_viewer, plots=not args.no_plots
        )
        print("products: " + ", ".join(f.name for f in files))
    return 0
