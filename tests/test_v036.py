"""Regression tests of the fixes of the v0.35 code review."""

import math
import warnings

import numpy as np
import pytest
import yaml

from swingersim.config import Config

warnings.simplefilter("ignore")


def _cfg(text):
    return Config.model_validate(yaml.safe_load(text))


# ---------------------------------------------------------------- backward evolving secular runs
def test_backward_evolving_secular_run_integrates():
    """segment_edges ordered the time grid ascending: a backward run with MIST evolution never
    integrated and returned the initial state at every output."""
    from swingersim.backends.evolving import segment_edges
    from swingersim.hierarchy import build_tree

    cfg = _cfg(
        """
topology: "2+2"
masses: {m1: 1.2, m2: 1.0, m3: 1.1, m4: 0.9}
orbits:
  A: {P: 3.0, e: 0.2, i: 90.0}
  B: {P: 4.0, e: 0.1, i: 90.0}
  out: {P: 2000.0, e: 0.3, i: 30.0}
stellar_model: {kind: mist, age_myr: 2000.0}
physics: {gr: true, evolution: true}
run: {t_max_yr: -1.5e9, n_output: 11, backend: secular}
"""
    )
    edges = segment_edges(cfg, build_tree(cfg)[1], -1.5e9, 0.01)
    assert edges[0] == 0.0 and edges[-1] == -1.5e9
    assert np.all(np.diff(edges) < 0)


# ---------------------------------------------------------------- periods with track mass loss
def test_period_at_t0_with_track_mass_loss_current_masses():
    from swingersim import units
    from swingersim.backends.nbody import build_simulation
    from swingersim.hierarchy import node_elements
    from swingersim.massloss import resolve_masses

    cfg = _cfg(
        """
topology: "2+2"
masses: {m1: 30.0, m2: 25.0, m3: 20.0, m4: 18.0}
orbits:
  A: {P: 5.0, e: 0.0, i: 90.0}
  B: {P: 7.0, e: 0.0, i: 90.0}
  out: {P: 3000.0, e: 0.2, i: 90.0}
stellar_model: {kind: mist, age_myr: 5.0}
physics: {gr: true, evolution: true, mass_loss: {source: track, masses: current}}
run: {t_max_yr: 1.0, n_output: 3}
"""
    )
    cfg = resolve_masses(cfg)
    sim, _, _, root, masses = build_simulation(cfg)
    ps = sim.particles
    pos = np.array([[p.x, p.y, p.z] for p in ps])
    vel = np.array([[p.vx, p.vy, p.vz] for p in ps])
    node = next(n for n in root.nodes() if n.name == "A")
    el = node_elements(node, pos, vel, masses)
    p_d = units.sma_to_period(el.a, masses[0] + masses[1]) * units.DAYS_PER_YEAR
    assert p_d == pytest.approx(5.0, rel=1e-9)


# ---------------------------------------------------------------- RV light time
RV_CFG = """
topology: '2+2'
masses: {m1: 1.2, m2: 1.0, m3: 1.1, m4: 0.9}
orbits: {A: {P: 1.5, i: 89.0, M0: 40}, B: {P: 3.0, i: 88.5, M0: 200}, out: {P: 7300.0, e: 0.3, i: 89.0, M0: 90}}
stars: {m1: {R: 1.2}, m2: {R: 1.0}, m3: {R: 1.1}, m4: {R: 0.9}}
run: {t_max_yr: 0.3, n_output: 61}
"""


def test_rv_light_time_converged_for_the_first_point():
    from swingersim.radvel import model_rv

    cfg = _cfg(RV_CFG)
    t0 = 2459000.0
    t = t0 + 10.3
    one = model_rv(cfg, ["m1"], [t], t0)[0]
    ref = model_rv(cfg, ["m1"] * 6, [t - 5, t - 4, t - 3, t - 2, t - 1, t], t0)[-1]
    mix = model_rv(cfg, ["m3", "m1"], [t, t + 1e-6], t0)[1]
    assert abs(one - ref) < 1e-3  # was 0.88 km/s
    assert abs(mix - ref) < 1e-3  # was 0.85 km/s


# ---------------------------------------------------------------- nested secular Roche
def test_nested_secular_roche_uses_the_right_semi_major_axis():
    from swingersim.backends.secular_nested import run_secular_nested

    raw = yaml.safe_load(
        """
topology: "3+1"
masses: {m1: 1.0, m2: 0.8, m3: 0.7, m4: 0.6}
orbits:
  inner: {P: 3.0,  e: 0.01, i: 80.0, Omega: 0.0,  omega: 30.0, M0: 0.0}
  mid:   {P: 60.0, e: 0.2,  i: 40.0, Omega: 30.0, omega: 60.0, M0: 90.0}
  out:   {P: 3000, e: 0.4,  i: 20.0, Omega: 70.0, omega: 0.0,  M0: 200.0}
physics: {gr: true, tides: true}
run: {t_max_yr: 3000, n_output: 100, backend: secular}
"""
    )
    raw["stars"] = {k: {"R": 1.0, "lag_s": 0.0} for k in ("m1", "m2", "m3", "m4")}
    r0 = run_secular_nested(Config.model_validate(raw)).roche
    for k in raw["stars"]:
        raw["stars"][k]["lag_s"] = 3000.0
    r1 = run_secular_nested(Config.model_validate(raw)).roche
    for k in ("m3", "m4"):  # stars of non-dissipative orbits: was 634x too large
        assert r1[k]["fill_max"] == pytest.approx(r0[k]["fill_max"], rel=0.05)


# ---------------------------------------------------------------- encounters at t = 0
def test_encounters_due_at_t0_enter_in_the_first_output_step():
    """The review suspected that encounters already approaching at t = 0 were dropped; they are
    applied by env.at(sim, 0.0) in the first output step (j = 0) - the flyby benchmark covers the
    case (its passage starts before t = 0). This test guards the add event at exactly t = 0."""
    import numpy as np

    from swingersim import environment as envmod

    enc = envmod.Encounter(t_yr=400.0, q_au=20.0, ecc=1.0, rotation=np.eye(3), m_per=1.0)
    cfg = Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"out": ["m1", "m2"]},
            "masses": {"m1": 1.0, "m2": 1.0},
            "orbits": {"out": {"a": 1.0, "e": 0.5}},
            "physics": {"gr": False, "encounters": {"n_star_pc3": 0.1}},
            "run": {"backend": "nbody", "t_max_yr": 800.0, "n_output": 3},
        }
    )
    orig = envmod.encounter_sequence
    envmod.encounter_sequence = lambda *a, **k: [enc]
    try:
        env = envmod.NBodyEnvironment(cfg, 2, 2.0, 1.0, 1.0, 800.0)
    finally:
        envmod.encounter_sequence = orig
    assert env.events[0][:2] == (0.0, "add")
    assert env.boundaries(0.0, 0.0) == [0.0]  # the j = 0 step ends at t = 0 and calls env.at


# ---------------------------------------------------------------- scatter
def test_scatter_rejects_overfilled_binary_and_parabolic_impact_parameter():
    from pydantic import ValidationError

    from swingersim.scatter import ScatterBody, ScatterSpec

    with pytest.raises(ValidationError, match="overfills its Roche lobe"):
        ScatterBody(masses=[1, 1], a_au=0.01, e_thermal=True, radii_rsun=[2, 2])
    target = ScatterBody(masses=[1, 1], a_au=1.0)
    proj = ScatterBody(masses=[1])
    with pytest.raises(ValidationError, match="parabolic"):
        ScatterSpec(target=target, projectile=proj, v_inf_kms=0.0, b_max_au=2.0)
    ScatterSpec(target=target, projectile=proj, v_inf_kms=0.0, r_peri_max_au=2.0)


def test_scatter_thermal_eccentricity_respects_the_roche_limit():
    from swingersim import units
    from swingersim.scatter import ScatterBody, _random_binary
    from swingersim.stellar import roche_lobe_eggleton

    body = ScatterBody(masses=[1, 1], a_au=0.05, e_thermal=True, radii_rsun=[1.5, 1.5])
    rng = np.random.default_rng(1)
    es = np.array([_random_binary(body, rng)["e"] for _ in range(2000)])
    e_max = 1 - 1.5 * units.RSUN_AU / (0.05 * roche_lobe_eggleton(1.0))
    assert es.max() <= e_max
    # f(e) = 2e truncated at e_max: mean 2/3 e_max
    assert es.mean() == pytest.approx(2.0 / 3.0 * e_max, rel=0.05)


# ---------------------------------------------------------------- emulator, sbi, weights
def test_emulator_classifier_observation_direction():
    from swingersim.emulator import _wanted_class

    assert _wanted_class("ghost_stable", {"value": None, "min": 0.5, "max": None}) == 1
    assert _wanted_class("ghost_stable", {"value": None, "min": None, "max": 0.5}) == 0
    assert _wanted_class("ghost_stable", {"value": 1.0}) == 1
    with pytest.raises(ValueError):
        _wanted_class("ghost_stable", {"value": None, "min": None, "max": None})


def test_limits_with_equals_and_missing_ok():
    from swingersim.ensemble import observation_weights
    from swingersim.sbi_fit import training_set

    rows = [
        {"member": 0, "direction": "forward", "fate": "intact", "x": 1.0, "y": 2.0},
        {"member": 1, "direction": "forward", "fate": "escape", "x": 1.1, "y": 2.0},
        {"member": 2, "direction": "forward", "fate": "intact", "x": math.nan, "y": 2.1},
    ]
    obs = {"fate": {"equals": "intact"}, "y": {"value": 2.0, "sigma": 0.1}}
    th, _x, cols, *_ = training_set(rows, ["x"], obs)  # crashed on equals
    assert len(th) == 1 and cols == ["y"]
    w = observation_weights(
        rows, {"fate": {"equals": "intact"}, "x": {"value": 1.0, "sigma": 0.1, "missing_ok": True}}
    )
    assert w[0] == pytest.approx(1.0) and w[1] == 0.0 and w[2] == pytest.approx(1.0)


# ---------------------------------------------------------------- mass transfer
def test_qcrit_main_sequence_and_core_helium_burning():
    from swingersim import masstransfer as mt

    assert mt.q_crit(1, 2.0, 0.0) == 3.0
    assert mt.q_crit(4, 2.0, 0.3) == 3.0
    assert mt.q_crit(1, 0.5, 0.0) == 0.695


# ---------------------------------------------------------------- jobs and GUI safety
def test_job_values_cannot_become_options(tmp_path):
    from swingersim.worker import JobError, build_command

    with pytest.raises(JobError, match="looks like an option"):
        build_command({"command": "validate", "ids": ["--html=/tmp/evil.html"]}, tmp_path)
    argv, _ = build_command({"command": "validate", "ids": ["gr-mercury"]}, tmp_path)
    assert argv[-1] == "gr-mercury"


def test_gui_job_ids_cannot_leave_the_folder(tmp_path):
    from swingersim.gui.server import GuiError, _find_job, safe_name

    with pytest.raises(GuiError):
        safe_name("../../outside/victim")
    with pytest.raises(GuiError):
        _find_job(tmp_path / "jobs", "../x")
    assert safe_name("001_run.v2") == "001_run.v2"


def test_gui_refuses_foreign_pages(tmp_path):
    import json
    import threading
    import urllib.error
    import urllib.request

    from swingersim.gui.server import make_server

    srv = make_server(tmp_path, "127.0.0.1", 0)
    port = srv.server_address[1]
    th = threading.Thread(target=srv.serve_forever, daemon=True)
    th.start()
    base = f"http://127.0.0.1:{port}"
    try:

        def post(path, ctype="application/json", origin=None, host=None):
            h = {"Content-Type": ctype}
            if origin:
                h["Origin"] = origin
            if host:
                h["Host"] = host
            req = urllib.request.Request(
                base + path, data=json.dumps({"yaml": ""}).encode(), headers=h
            )
            try:
                return urllib.request.urlopen(req).status
            except urllib.error.HTTPError as exc:
                return exc.code

        assert post("/api/validate") == 200
        assert post("/api/validate", ctype="text/plain") in (400, 500)
        assert post("/api/validate", origin="http://evil.example") == 403
        assert post("/api/validate", host="evil.example") == 403
    finally:
        srv.shutdown()
