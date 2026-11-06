"""Mergers followed by the run of the merged system (physics.mergers) and q_crit_ms."""

import math
import warnings

import numpy as np
import pytest
import yaml

from swingersim import units
from swingersim.config import Config

warnings.simplefilter("ignore")

COLLIDING = """
topology: '2+2'
masses: {m1: 1.0, m2: 0.8, m3: 0.9, m4: 0.7}
orbits:
  A: {P: 10.0, e: 0.95, i: 80, M0: 120}
  B: {P: 3.0, e: 0.0, i: 85}
  out: {P: 600.0, e: 0.2, i: 70, M0: 30}
stellar_model: {kind: mist, age_myr: 1000.0}
physics: {gr: true}
run: {t_max_yr: 10.0, n_output: 101}
"""


def _cfg(text, mergers=True):
    d = yaml.safe_load(text)
    if mergers:
        d["physics"]["mergers"] = {}
    return Config.model_validate(d)


def test_q_crit_ms_option():
    from swingersim import masstransfer as mt

    assert mt.q_crit(1, 2.0, 0.0, 1.6) == 1.6
    assert mt.q_crit(1, 0.5, 0.0, 1.6) == 0.695  # deeply convective: unchanged
    from swingersim.config import MassTransferSpec

    spec = MassTransferSpec(q_crit_ms=1.6)
    assert spec.q_crit_ms == 1.6


def test_merger_prescriptions():
    from swingersim.config import MergerSpec
    from swingersim.merger import mass_lost_fraction, rejuvenation_factor

    s = MergerSpec()
    assert mass_lost_fraction(s, 1.0, 1.0) == pytest.approx(0.075)  # Glebbeek et al. 2013
    assert mass_lost_fraction(s, 1.0, 0.25) == pytest.approx(0.3 * 0.25 / 1.25**2)
    assert rejuvenation_factor(s, 3.0) == pytest.approx(0.81)
    assert rejuvenation_factor(s, 100.0) == pytest.approx(0.64)
    assert 0.64 < rejuvenation_factor(s, 20.0) < 0.81
    assert rejuvenation_factor(MergerSpec(rejuvenation="hurley"), 2.0) == 0.1


def test_merged_main_sequence_age_is_rejuvenated():
    from swingersim.config import MergerSpec
    from swingersim.merger import _ms_window, merged_age

    sm = _cfg(COLLIDING).stellar_model
    age = merged_age(sm, (1000.0, 1000.0), (1.0, 0.8), (1.0, 0.8), 1.667, MergerSpec())
    t0n, t1n = _ms_window(sm, 1.667)
    f_new = (age - t0n) / (t1n - t0n)
    f = []
    for m in (1.0, 0.8):
        t0, t1 = _ms_window(sm, m)
        f.append((1000.0 - t0) / (t1 - t0))
    assert f_new == pytest.approx(0.81 * (1.0 * f[0] + 0.8 * f[1]) / 1.8, rel=1e-9)


def test_collision_continues_as_a_triple():
    from swingersim.backends.nbody import run_nbody

    stop = run_nbody(_cfg(COLLIDING, mergers=False))
    assert stop.event == "collision" and stop.event_detail["stars"] == ["m1", "m2"]
    res = run_nbody(_cfg(COLLIDING))
    assert res.event is None and res.times[-1] == pytest.approx(10.0)
    assert res.system["mergers_n"] == 1 and res.system["merger_1_stars"] == "m1+m2"
    assert res.hierarchy[-1] == "[1,[3,4]]"
    assert res.transitions[0]["kind"] == "merger"
    n_a = len(stop.times)
    assert np.all(np.isnan(res.orbits["A"]["a"][n_a:]))  # the merged binary is gone
    assert np.all(np.isfinite(res.orbits["B"]["a"]))
    # impulsive mass loss: the outer orbit right after the merger from the state at the
    # collision with the reduced mass (Hadjidemetriou 1963)
    st = stop.final_state
    m, p, v = st["masses"], st["pos"], st["vel"]
    m_new = res.system["merger_1_mass_msun"]
    assert m_new == pytest.approx(1.8 * (1 - 0.3 * 0.8 / 1.8**2))
    rc = (m[0] * p[0] + m[1] * p[1]) / (m[0] + m[1])
    vc = (m[0] * v[0] + m[1] * v[1]) / (m[0] + m[1])
    rb = (m[2] * p[2] + m[3] * p[3]) / (m[2] + m[3])
    vb = (m[2] * v[2] + m[3] * v[3]) / (m[2] + m[3])
    r, w = np.linalg.norm(rb - rc), np.linalg.norm(vb - vc)
    mu = units.G * (m_new + m[2] + m[3])
    a_pred = 1.0 / (2.0 / r - w * w / mu)
    assert res.orbits["out"]["a"][n_a] == pytest.approx(a_pred, rel=1e-3)
    assert res.orbits["out"]["a"][n_a] > stop.orbits["out"]["a"][0]  # the orbit widens


def test_merged_system_unbound_ends_the_run():
    from swingersim.backends.nbody import run_nbody

    d = yaml.safe_load(COLLIDING)
    d["topology"] = "custom"
    d["hierarchy"] = {"A": ["m1", "m2"]}
    d["masses"] = {"m1": 1.0, "m2": 0.8}
    d["orbits"] = {"A": d["orbits"]["A"]}
    d["physics"]["mergers"] = {}
    res = run_nbody(Config.model_validate(d))
    assert res.event == "merger"
    assert res.event_detail["merger"]["stop"] == "a single star remains"


def test_common_envelope_merger_continues_secular():
    from swingersim.backends.secular import run_secular

    cfg = Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"out": ["A", "m3"], "A": ["m1", "m2"]},
            "masses": {"m1": 1.8, "m2": 0.5, "m3": 1.0},
            "orbits": {"A": {"P": 200.0}, "out": {"P": 3650.0, "i": 10.0}},
            "stellar_model": {"kind": "sse", "age_myr": 1652.75, "z": 0.02},
            "physics": {
                "gr": False,
                "evolution": True,
                "mass_transfer": {"common_envelope": {"alpha": 0.2}},
                "mergers": {},
            },
            "run": {"backend": "secular", "t_max_yr": 2.0e5, "n_output": 21},
        }
    )
    res = run_secular(cfg)
    mg = res.system["_mergers"][0]
    assert mg["kind"] == "common_envelope_merger"
    assert mg["mass_msun"] < 2.3 and mg["hierarchy"] == "[1,3]"
    assert res.times[-1] > mg["time_yr"]
    # the envelope left at once: the outer orbit widened, a' / a ~ M / M' (circular)
    i = int(np.searchsorted(res.times, mg["time_yr"])) + 1
    ratio = res.orbits["out"]["a"][i] / res.orbits["out"]["a"][0]
    assert 1.0 < ratio < 1.3
    assert math.isnan(res.orbits["A"]["a"][-1])


def test_gui_shows_the_merger(tmp_path):
    from swingersim.backends.nbody import run_nbody
    from swingersim.gui.digest import run_digest
    from swingersim.io import save_run

    cfg = _cfg(COLLIDING)
    path = tmp_path / "run.h5"
    save_run(path, cfg, run_nbody(cfg))
    d = run_digest(path)
    fact = next(f for f in d["facts"] if f["k"].startswith("Merger"))
    assert fact["v"].startswith("m1 + m2")
    stars = {s["name"]: s for s in d["data"]["stars"]}
    assert stars["m2"]["mass"][1] is None and "merged into m1" in stars["m2"]["phase"]
    assert stars["m1"]["mass"][1] == pytest.approx(1.8 * (1 - 0.3 * 0.8 / 1.8**2))


def test_fit_digest_reads_the_mcmc_convergence():
    from swingersim.gui.digest import _TAU

    line = "sampler emcee: 72 walkers x 400 steps, first 200 discarded, seed 1, autocorrelation time tau max 37.2 steps (steps / tau 10.8)"
    assert _TAU.search(line).groups() == ("37.2", "10.8")
