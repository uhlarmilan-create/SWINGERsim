"""Hierarchy detector and run-time tracking of hierarchy changes."""

import numpy as np

from swingersim import units
from swingersim.backends.nbody import run_nbody
from swingersim.config import Config
from swingersim.diagnostics import hierarchy as hd
from swingersim.hierarchy import build_tree, initialize


def _detect(cfg):
    root, m = build_tree(cfg)
    pos, vel = initialize(root, m)
    return hd.detect(pos, vel, m, units.G), (hd.from_node(root),)


def test_examples_detected_as_configured(example_2p2, example_3p1):
    for cfg, lab in ((example_2p2, "2+2"), (example_3p1, "3+1")):
        found, configured = _detect(cfg)
        assert found == configured
        assert hd.label(found) == lab


def test_strings_and_labels():
    assert hd.to_string((((0, 1), (2, 3)),)) == "[[1,2],[3,4]]"
    assert hd.label(((((0, 1), 2), 3),)) == "3+1"
    assert hd.label((((0, 1), 2), 3)) == "triple + single"
    assert hd.label(((0, 1), (2, 3))) == "binary + binary"


def test_canonical_order():
    assert hd.canonical(((3, 2), (1, 0))) == ((0, 1), (2, 3))


def test_classify_change():
    q22 = (((0, 1), (2, 3)),)
    assert hd.classify_change(q22, (((0, 2), (1, 3)),)) == "exchange"
    assert hd.classify_change(q22, ((((0, 1), 2), 3),)) == "reconfiguration"
    assert hd.classify_change(q22, (((0, 1), 2), 3)) == "ejection"
    assert hd.classify_change(q22, ((0, 1), (2, 3))) == "breakup"
    assert hd.classify_change((((0, 1), 2), 3), ((((0, 1), 2), 3),)) == "capture"


def test_escaping_star_is_detected():
    m = np.ones(4)
    pos = np.array([[0, 0, 0], [0.1, 0, 0], [5, 0, 0], [100, 0, 0]], float)
    v_orb = np.sqrt(units.G * 2 / 0.1)
    vel = np.array([[0, -v_orb / 2, 0], [0, v_orb / 2, 0], [0, 3.0, 0], [50.0, 0, 0]])
    found = hd.detect(pos, vel, m, units.G)
    assert hd.label(found) == "triple + single"
    assert found[-1] == 3


def _unstable(stop):
    return Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 1, "m2": 1, "m3": 1, "m4": 1},
            "orbits": {
                "A": {"P": 50, "e": 0.1, "i": 30},
                "B": {"P": 50, "e": 0.1, "i": 60, "Omega": 90, "M0": 50},
                "out": {"P": 120, "e": 0.4},
            },
            "physics": {"gr": False},
            "run": {"t_max_yr": 5, "n_output": 500, "stop_on_hierarchy_change": stop},
        }
    )


def test_unstable_system_reports_transition():
    res = run_nbody(_unstable(stop=False))
    assert res.transitions, "an unstable 2+2 must change its hierarchy"
    assert res.transitions[0]["from"] == "[[1,2],[3,4]]"


def test_stop_on_hierarchy_change():
    res = run_nbody(_unstable(stop=True))
    assert res.event == "hierarchy_change"
    assert res.event_time == res.transitions[0]["time_yr"]
    assert res.times[-1] < 5.0


def test_stable_example_has_no_transition(example_2p2):
    cfg = example_2p2.model_copy(deep=True)
    cfg.run.t_max_yr, cfg.run.n_output = 5.0, 200
    res = run_nbody(cfg)
    assert res.transitions == []
    assert set(res.hierarchy) == {"[[1,2],[3,4]]"}


def test_zlk_diagnostics_example(example_2p2):
    cfg = example_2p2.model_copy(deep=True)
    cfg.run.t_max_yr, cfg.run.n_output = 2.0, 50
    diag = {d.inner: d for d in run_nbody(cfg).diagnostics}
    assert set(diag) == {"A", "B"}
    assert not diag["A"].in_zlk_window and diag["A"].stability_ratio > 1
    assert abs(diag["A"].i_mut0_deg - 18.1) < 0.5


def test_extremes_between_outputs():
    """Coarse outputs miss the ZLK peak; the internal sampling must still find it."""
    from swingersim.config import load_config

    from .conftest import EXAMPLES

    cfg = load_config(EXAMPLES / "example_2p2_zlk.yaml")
    cfg.run.t_max_yr, cfg.run.n_output = 20.0, 6
    res = run_nbody(cfg)
    b = next(d for d in res.diagnostics if d.inner == "B")
    sampled = float(np.nanmax(res.orbits["B"]["e"]))
    assert b.e_max_true > 0.79  # the first ZLK peak (e ~ 0.81) at t ~ 8 yr
    assert b.e_max_true >= sampled
    assert b.q_min_over_radii < 3.5  # ZAMS radii of the default stellar model
