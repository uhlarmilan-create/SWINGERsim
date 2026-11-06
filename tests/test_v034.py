"""One-sided data in fit.inclination (upper / lower limits, eclipse / no eclipse),
the primary-secondary diagnostic of stage 2 and the red-noise factor of light curves."""

from pathlib import Path

import numpy as np
import yaml

from swingersim.config import Config, load_config

ROOT = Path(__file__).resolve().parents[1]
QX = ROOT / "configs/squadra_test3/qx_cas_fit_v2.yaml"


def _spec(file, min_sigma=0.0):
    class Spec:
        pass

    Spec.file, Spec.pair, Spec.label, Spec.min_sigma_deg = str(file), "inner", None, min_sigma
    return Spec


def test_hinge_residuals():
    from swingersim.incfit import ObservedInclination, residuals

    q = ObservedInclination(
        np.zeros(6),
        np.array([80.0, 75.0, 75.0, 70.0, np.nan, np.nan]),
        np.ones(6),
        np.array([""] * 6, dtype=object),
        "inner",
        "t",
        np.array(["", "upper", "upper", "lower", "no_eclipse", "eclipse"], dtype=object),
    )
    m = np.array([[78.0, 74.0, 77.0, 72.0, 76.0, 70.0], [73.0] * 6])
    r = residuals(q, m)
    # measurement: (80 - 78) / 1; upper satisfied: 0; upper violated: 75 - 77 = -2;
    # lower satisfied: 0; no eclipse violated (76 > limit 73): -3; eclipse violated (70 < 73): +3
    assert np.allclose(r, [2.0, 0.0, -2.0, 0.0, -3.0, 3.0])


def test_load_limits_and_eclipse_limit_of_qx_cas():
    from swingersim.incfit import load, model_inclination, residuals

    q = load(_spec(ROOT / "data/observed/qx_cas/inclination_limits.csv"), None)
    assert list(q.limit).count("no_eclipse") == 4 and list(q.limit).count("eclipse") == 1
    cfg = load_config(ROOT / "configs/squadra_test3/qx_cas_fit.yaml")  # the start hypothesis
    m = model_inclination(cfg, "inner", q.time, 2444500.0, with_limit=True)
    assert m.shape == (2, len(q))
    # eclipse limit of the 6.0-d pair: arccos((R1 + R2) / r_conj), 69-74 deg as the eccentric
    # orbit turns (apsidal motion changes r_conj at the nearer conjunction)
    assert np.all((m[1] > 66) & (m[1] < 78)) and np.ptp(m[1]) > 1.0
    # the start hypothesis (i falls through the limit around 2004-2008) satisfies the eclipses of
    # the 1930s and the TESS rows; in 2006 it is still ~1.3 deg above the limit (r ~ -4.4 sigma)
    r = residuals(q, m)
    tess = np.array(["TESS" in s for s in q.source])
    assert np.all(r[tess] == 0.0) and r[0] == 0.0
    assert -6.0 < r[4] < 0.0


def test_sq5_6_solution_violates_no_eclipse():
    """The first QX Cas fit (sq5_6) kept the pair eclipsing in the TESS years: the new data
    must penalise it (tested with the configuration of the fit, outputs not needed)."""
    from swingersim.incfit import load, model_inclination, residuals

    # qx_cas_fit_v2.yaml starts from the sq5_6 solution (e cos w -0.0170, e sin w 0.2187, M0 328.51,
    # i 81.91; outer P 299.18, i 80.23; node difference +6.06 deg); sq5_6 had m3 at the MIST limit
    d = yaml.safe_load(QX.read_text())
    d["masses"]["m3"] = 6.0
    cfg = Config.model_validate(d)
    q = load(_spec(ROOT / "data/observed/qx_cas/inclination_limits.csv"), None)
    m = model_inclination(cfg, "inner", q.time, 2444500.0, with_limit=True)
    r = residuals(q, m)
    tess = np.array(["TESS" in s for s in q.source])
    assert np.sum(r[tess] ** 2) > 1.0  # at least one TESS epoch with the model eclipsing


def test_primary_secondary_split_of_v907_sco():
    """V907 Sco (SQUADRA part III): e ~ 0 gives an 11 / -13 min split of the TESS minima, the
    second-pass start (e = 0.010, omega = 301.7, anchored) none."""
    from swingersim import dynfit, ocfit

    mins = ocfit.read_minima(ROOT / "data/observed/v907_sco/minima_tess.csv")
    obs = dynfit.observed_from_minima(mins, None, None)
    d = yaml.safe_load((ROOT / "configs/squadra_test3/v907_sco_fit_v2.yaml").read_text())
    d.pop("fit")
    bad = yaml.safe_load(yaml.safe_dump(d))
    bad["orbits"]["inner"].update(e=0.0003, omega=121.7)  # sq4_5: e cos w ~ -0.00016
    good = yaml.safe_load(yaml.safe_dump(d))
    good["orbits"]["inner"].update(P=3.77475844, M0=224.7635)  # anchored (fit --anchor)
    s_good = dynfit.primary_secondary_split(
        obs, dynfit.model_times(Config.model_validate(good), obs, 2459383.9)
    )
    s_bad = dynfit.primary_secondary_split(
        obs, dynfit.model_times(Config.model_validate(bad), obs, 2459383.9)
    )
    assert not s_good["inner"]["warn"] and abs(s_good["inner"]["diff"]) < 3.0
    assert s_bad["inner"]["warn"] and abs(s_bad["inner"]["diff"]) > 10.0


def test_red_noise_beta():
    from swingersim.lcfit import red_noise_beta

    rng = np.random.default_rng(3)
    t = np.arange(0.0, 6.0, 2.0 / 1440.0)
    white = rng.normal(0.0, 1e-3, len(t))
    assert red_noise_beta(t, white) < 1.15
    red = white + 2e-3 * np.sin(2 * np.pi * t / 0.4)
    assert red_noise_beta(t, red) > 2.0
    assert red_noise_beta(t[:10], white[:10]) == 1.0


def test_new_configs_validate():
    for name in ("qx_cas_fit_v2", "v907_sco_fit_v2", "v907_sco_fit_lc_v2"):
        load_config(ROOT / f"configs/squadra_test3/{name}.yaml")


def test_gui_fit_digest(tmp_path):
    """The Outputs page of a fit: parameters, minima with the primary-secondary check,
    inclinations with limits, the other fits of the system."""
    from swingersim.gui.digest import fit_digest

    for name, red in (("fit", 0.82), ("fit_v2", 1.01)):
        f = tmp_path / name
        f.mkdir()
        (f / "stage2_report.txt").write_text(
            "swingersim fit, stage 2 - joint dynamical model (N-body) of all eclipsing pairs\n"
            "  t = 0 of the configuration: BJD 2459383.9\n"
            f"  8 minima, 1 priors, 3 free parameters, 50 model evaluations, 12 s\n"
            f"  chi^2 460.3 -> 189.3 (reduced {red}); errors of the minima scaled by 1\n"
            "  rms of the residuals [min], start -> fitted:\n"
            "    inner       8 minima      23.753 ->   12.796\n"
        )
        (f / "stage2_params.csv").write_text(
            "parameter,start,value,error,unit,prior_mean,prior_sigma\n"
            "orbits.inner.P,3.7747,3.77472,1e-07,d,,\n"
            "orbits.out.e,0.467,0.47,0.017,,0.47,0.05\n"
            "masses.m3,3.0,6.629,0.0,Msun,3.0,2.0\n"
        )
        rows = [
            "pair,type,bjd_tdb,sigma_d,model_start_bjd,model_bjd,o_minus_model_min,source,sigma_quoted_d,status"
        ]
        for k in range(8):
            r = 11.0 if k % 2 == 0 else -13.0
            rows.append(
                f"inner,{'ps'[k % 2]},{2459362.3 + k * 1.887:.6f},0.0007,{2459362.3:.6f},{2459362.3:.6f},{r + 0.1 * k:.3f},TESS S39,0.0013,fitted"
            )
        (f / "stage2_residuals.csv").write_text("\n".join(rows) + "\n")
        (f / "stage2_inclination.csv").write_text(
            "label,pair,time,i_obs,sigma,limit,i_start,i_model,i_eclipse_limit,z,source\n"
            "L,inner,2436800.0,85.4,1.0,,84.0,85.0,73.0,0.4,UBV\n"
            "L,inner,2458777.5,,0.3,no_eclipse,75.0,76.0,72.0,-13.3,TESS S17\n"
        )
    d = fit_digest(tmp_path / "fit", full=True, root=tmp_path)
    data = d["data"]
    assert d["setup"].startswith("joint N-body fit")
    assert [p["name"] for p in data["params"]] == ["orbits.inner.P", "orbits.out.e", "masses.m3"]
    assert data["params"][2]["tiny_error"]  # m3 on a bound: error 0 (bounds not known here)
    assert not data["params"][0]["tiny_error"]  # a precise period is not on a bound (0.36)
    assert any("primary and secondary minima are offset" in w for w in data["warnings"])
    assert any(q["k"] == "primary vs secondary" for q in d["quality"])
    pair = data["minima"]["pairs"][0]
    assert pair["n"] == 8 and pair["mean_p"] > 10 and pair["mean_s"] < -12
    assert [r["limit"] for r in data["inclination"]] == ["", "no_eclipse"]
    assert {s["name"] for s in data["siblings"]} == {"fit", "fit_v2"}
    card = fit_digest(tmp_path / "fit", full=False, root=tmp_path)
    assert "data" not in card and card["facts"]


def test_gui_measured_minima_against_timings(tmp_path):
    from swingersim.gui.digest import compare_minima, measured_minima_for

    P, bjd0 = 2.0, 2459000.0
    rows = ["binary,kind,epoch,t_obs_yr,o_c_d,P_d,T0_yr,impact,bjd_obs"]
    for e in range(20):
        for kind, off in (("primary", 0.0), ("secondary", 1.0)):
            t = (e * P + off) / 365.25
            rows.append(f"A,{kind},{e},{t:.9f},0.0,{P},0.0,0.1,{bjd0 + e * P + off:.6f}")
    tim = tmp_path / "timings.csv"
    tim.write_text("\n".join(rows) + "\n")
    mins = tmp_path / "minima.csv"
    mins.write_text(
        "binary,kind,bjd_tdb,sigma_d,ref\n"
        + "".join(
            f"A,p,{bjd0 + e * P + 3 / 1440:.6f},0.0005,x\nA,s,{bjd0 + e * P + 1.0 - 2 / 1440:.6f},0.0005,x\n"
            for e in range(2, 8)
        )
    )
    out = compare_minima(tim, mins, tmp_path)
    a = out["binaries"]["A"]
    assert abs(a["primary"]["mean"] - 3.0) < 0.01 and abs(a["secondary"]["mean"] + 2.0) < 0.01
    assert abs(out["bjd_t0"] - bjd0) < 1e-6
    # the measured minima of a configuration: its own fit.minima or a sibling's
    cdir = tmp_path / "configs"
    cdir.mkdir()
    (cdir / "sys_published.yaml").write_text("topology: custom\n")
    (cdir / "sys_fit.yaml").write_text("fit:\n  minima: ../minima.csv\n")
    assert measured_minima_for(tmp_path, "configs/sys_published.yaml") == mins.resolve()


def test_system_model_of_v907_sco():
    """The picture of a configuration (GUI panel, reports): hierarchy, orbits, eclipse geometry, stability."""
    from swingersim.systemview import system_model_from_file

    m = system_model_from_file(ROOT / "configs/squadra_test3/v907_sco_published.yaml")
    assert m["hierarchy"] == "((m1 + m2) + m3)"
    names = [o["name"] for o in m["orbits"]]
    assert names == ["out", "inner"]
    inner = m["orbits"][1]
    assert abs(inner["P_d"] - 3.7748) < 0.01 and abs(inner["i_deg"] - 76.61) < 0.01
    assert len(inner["path_left"]) == 240
    p = m["pairs"][0]
    assert p["eclipsing"] and 72.0 < p["i_limit_deg"] < 74.0
    assert m["triples"][0]["stable"] and abs(m["mutual"][0]["deg"] - 26.2) < 0.2


def test_gui_model_curves(tmp_path):
    from swingersim.gui import modelviews as MV

    cfg = ROOT / "configs/squadra_test3/v907_sco_published.yaml"
    (tmp_path / "output").mkdir()
    c = MV.inc_curve(tmp_path, cfg, "inner", 2459383.9, 2019.0, 2023.0, n=5)
    assert len(c["i"]) == 5 and 72 < c["i"][0] < 75 and c["limit"][0] > 70
    oc = MV.oc_curve(tmp_path, cfg, 2459383.9, 2021.3, 2021.4)
    assert oc["eclipses"]["inner"]["primary"] and oc["eclipses"]["inner"]["secondary"]
    assert MV.inc_curve(tmp_path, cfg, "inner", 2459383.9, 2019.0, 2023.0, n=5) == c  # cached
