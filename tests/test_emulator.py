import math

import numpy as np
import pytest

pytest.importorskip("sklearn")

from swingersim import emulator as emu


def _rows(n=600, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for m in range(n):
        x1, x2 = rng.uniform(-2, 2), rng.uniform(1, 100)
        rows.append(
            {
                "member": m,
                "direction": "forward",
                "orbits.A.e": x1,
                "masses.m1": x2,
                "rate": math.sin(x1) + 0.5 * math.log10(x2),
                "stable": float(x1 + math.log10(x2) > 1.0),
                "fate": "intact" if x1 > 0 else "collision",
            }
        )
    rows.append({"member": n, "direction": "forward", "fate": "error"})
    return rows


def test_dataset_and_training(tmp_path):
    rows = _rows()
    X, Y, C, keep = emu.dataset(
        rows, ["orbits.A.e", "masses.m1"], ["rate"], ["stable", "fate==intact"]
    )
    assert X.shape == (600, 2) and len(keep) == 600  # the error member is left out
    em = emu.train(
        X,
        Y,
        C,
        ["orbits.A.e", "masses.m1"],
        ["rate"],
        ["stable", "fate==intact"],
        log_inputs=["masses.m1"],
        hidden=(32, 32),
        n_models=3,
        seed=1,
    )
    assert em.metrics["rate:r2"] > 0.98
    assert em.metrics["stable:score"] > 0.93
    assert em.metrics["fate==intact:score"] > 0.95
    p = em.predict([[0.5, 10.0]])
    assert p["mean"][0, 0] == pytest.approx(math.sin(0.5) + 0.5, abs=0.1)
    path = em.save(tmp_path / "em.pkl")
    em2 = emu.load(path)
    assert np.allclose(em2.predict([[0.5, 10.0]])["mean"], p["mean"])
    assert (tmp_path / "em.json").exists()
    d = emu.diagnostics(em2)
    assert abs(np.nanstd(d["targets"]["rate"]["z_total"]) - 1.0) < 0.5
    assert d["classes"]["stable"]["auc"] > 0.95
    assert emu.diagnostic_plot(em2, tmp_path / "diag.png").exists()


def test_classification_scores():
    s = emu.classification_scores([1, 1, 0, 0], [1, 0, 0, 0])
    assert s["score"] == 0.75 and s["precision_1"] == 1.0 and s["recall_1"] == 0.5


def test_emulated_weights_select_the_truth():
    rows = _rows(800)
    X, Y, C, _ = emu.dataset(rows, ["orbits.A.e", "masses.m1"], ["rate"], [])
    em = emu.train(
        X,
        Y,
        C,
        ["orbits.A.e", "masses.m1"],
        ["rate"],
        [],
        log_inputs=["masses.m1"],
        hidden=(32, 32),
        n_models=2,
    )
    rng = np.random.default_rng(3)
    Xs = np.column_stack([rng.uniform(-1.4, 1.4, 5000), np.full(5000, 10.0)])
    w = emu.emulated_weights(em, Xs, {"rate": {"value": math.sin(1.0) + 0.5, "sigma": 0.02}})
    x_best = np.sum(w * Xs[:, 0]) / np.sum(w)
    assert x_best == pytest.approx(1.0, abs=0.15)


def test_ensemble_training_set_options():
    from swingersim.config import load_config
    from swingersim.ensemble import _hierarchical, _ma01_min, draw_members
    from swingersim.hierarchy import build_tree

    cfg = load_config("configs/validation/vynatheya23_2p2_ghost.yaml")
    cfg.ensemble.n = 40
    ms = draw_members(cfg)
    assert len(ms) == 40
    from swingersim.config import Config

    n_unstable = 0
    for m in ms:
        c = Config.model_validate(m.config)
        assert _hierarchical(build_tree(c)[0])
        n_unstable += _ma01_min(c) <= 1.0
    assert n_unstable < 40  # thinned
    assert [m.params for m in draw_members(cfg)] == [m.params for m in ms]  # deterministic


def test_ghost_timeout():
    from swingersim.config import load_config
    from swingersim.stability import ghost_test

    cfg = load_config("configs/validation/vynatheya23_2p2_ghost.yaml")
    cfg.ensemble = None
    g = ghost_test(cfg, n_outer=100, max_wall_s=1e-9)
    assert g.reason == "timeout" and not g.stable


def test_angle_inputs_wrap_around():
    rng = np.random.default_rng(3)
    ang = rng.uniform(0, 360, 800)
    X = ang[:, None]
    Y = np.cos(np.radians(ang))[:, None]
    em = emu.train(
        X, Y, np.zeros((800, 0)), ["w"], ["y"], angle_inputs=["w"], hidden=(16,), n_models=2, seed=0
    )
    assert em.metrics["y:r2"] > 0.99
    # 359 deg and 1 deg are neighbours
    p = em.predict([[359.0], [1.0]])["mean"][:, 0]
    assert abs(p[0] - p[1]) < 0.05
    with pytest.raises(ValueError):
        emu.train(X, Y, np.zeros((800, 0)), ["w"], ["y"], angle_inputs=["v"])


def test_target_scale_takes_out_a_power_law():
    rng = np.random.default_rng(4)
    P = 10 ** rng.uniform(1, 3.5, 1500)
    w = rng.uniform(0, 360, 1500)
    X = np.column_stack([P, w])
    Y = (np.sin(np.radians(w)) / P**2)[:, None]
    kw = {"angle_inputs": ["w"], "log_inputs": ["P"], "hidden": (32, 32), "n_models": 2}
    em = emu.train(
        X, Y, np.zeros((1500, 0)), ["P", "w"], ["y"], target_scales={"y": ("P", -2)}, **kw
    )
    assert em.metrics["y:r2"] > 0.98
    p = em.predict([[100.0, 90.0]])["mean"][0, 0]
    assert p == pytest.approx(1e-4, rel=0.1)
