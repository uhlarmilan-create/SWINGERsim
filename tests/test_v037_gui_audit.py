"""Outputs page: facts of stability checks, chaos maps, scattering, pulsations and viewers read
from what the jobs print; number formatting; the rms of fitted light curves."""

from swingersim.gui.digest import fmt, other_digest

STAB = """stability of configs/x.yaml (3+1)
nested triples (ratio > 1 = stable):
  mid in out (+m4)       a_in/a_out = 0.1255, e_in = 0.224, e_out = 0.274, i_mut = 1.3 deg, q_out = 0.220: Mardling-Aarseth 1.63, Vynatheya+22 1.85
  inner in mid (+m3)     a_in/a_out = 0.1405, e_in = 0.007, e_out = 0.224, i_mut = 0.4 deg, q_out = 0.477: Mardling-Aarseth 0.98, Vynatheya+22 1.86
machine-learning classifier (Vynatheya+23): stable
N-body ghost test (100 outer orbits, Newtonian): stable after 100.0 outer orbits, max relative divergence 9.89e-03  [15.8 s]
"""
CHAOS = """chaos map: orbits.out.P, 6 points, 68.4463 yr each
  regular (MEGNO < 4): 5
  chaotic (MEGNO >= 4): 1
  escape or collision: 0
  smallest orbits.out.P with regular motion from there outwards: 240
"""
SCAT = """example: 24 encounters
  outcome                                           n  fraction    1-sigma interval  sigma [AU^2]
  flyby                                            10     0.417  [0.321, 0.519]  2.788e+04
  collision                                         7     0.292  [0.209, 0.391]  1.952e+04

  energy error: median 1.8e-13, max 3.8e-11
"""


def test_stability_facts():
    f = {x["k"]: x for x in other_digest("stability", STAB)["facts"]}
    assert f["mid in out: stability ratio"]["v"] == "MA01 1.63 · V22 1.85"
    assert f["inner in mid: stability ratio"]["tone"] == "bad"  # MA01 below 1
    assert f["machine-learning classifier"]["v"] == "stable"
    assert f["N-body ghost test"]["v"] == "stable" and "9.89e-03" in f["N-body ghost test"]["note"]


def test_chaosmap_scatter_pulsation_view_facts():
    f = {x["k"]: x["v"] for x in other_digest("chaosmap", CHAOS)["facts"]}
    assert f["regular"] == "5" and f["chaotic"] == "1"
    f = {x["k"]: x["v"] for x in other_digest("scatter", SCAT)["facts"]}
    assert f["encounters"] == "24" and f["flyby"] == "41.7 %" and f["collision"] == "29.2 %"
    p = other_digest(
        "pulsation", "m1: delta_scuti (F)  P = 0.0711228 d  dP/dt = 5.35e-06 s/yr  R = 2 Rsun"
    )
    assert p["facts"][0]["v"] == "P 0.0711228 d"
    v = other_digest("view", "saved a\nsaved b\n2 viewers written, 5 up to date, 0 failed")
    assert v["facts"][0]["v"] == "2"
    assert other_digest("emulator-train", "anything") == {}


def test_fmt_hundreds_and_thousands():
    assert fmt(-1300.4, 3) == "-1300" and fmt(100, 1) == "100" and fmt(0.5, 3) == "0.5"
    assert fmt(2459383.9, 10).startswith("2.4")  # >= 1e5: scientific


def test_light_curve_rms_leaves_out_the_outliers(tmp_path):
    from swingersim.gui.digest import fit_digest

    (tmp_path / "stage2_report.txt").write_text("swingersim fit, stage 2\n")
    rows = ["source,segment,bjd_tdb,flux,sigma,sigma_quoted,model_start,model,o_minus_model,status"]
    rows += [
        f"T,S1,{2459000 + k / 100},1,0.001,0.001,1,1,{0.001 if k % 2 else -0.001},fitted"
        for k in range(100)
    ]
    rows += ["T,S1,2459002,0.6,0.001,0.001,1,1,-0.4,outlier"]
    (tmp_path / "stage2_lc.csv").write_text("\n".join(rows) + "\n")
    d = fit_digest(tmp_path, full=False)
    f = next(x for x in d["facts"] if x["k"] == "light curves: rms")
    assert f["v"] == "1 ppt" and "1 outliers left out" in f["note"]
