"""General hierarchies: sextuples, 3+3, ten stars."""

import numpy as np
import pytest

from swingersim.config import Config
from swingersim.hierarchy import build_tree
from swingersim.stability import mlp_directory, mlp_verdicts, nested_triples, sub_quadruples


def six(topology="2+2+2", t=0.5, backend="nbody"):
    orbits = {
        "2+2+2": {
            "A": {"P": 3.0, "i": 85.0},
            "B": {"P": 2.0, "i": 80.0, "Omega": 30.0},
            "AB": {"P": 200.0, "e": 0.1, "i": 70.0},
            "C": {"P": 4.0, "i": 60.0},
            "out": {"P": 9000.0, "e": 0.2, "i": 40.0},
        },
        "3+3": {
            "A": {"P": 3.0, "i": 85.0},
            "T1": {"P": 150.0, "e": 0.1, "i": 70.0},
            "B": {"P": 2.0, "i": 80.0},
            "T2": {"P": 120.0, "e": 0.1, "i": 60.0},
            "out": {"P": 9000.0, "e": 0.2, "i": 40.0},
        },
    }[topology]
    return Config.model_validate(
        {
            "topology": topology,
            "masses": {f"m{k}": 1.0 for k in range(1, 7)},
            "orbits": orbits,
            "run": {"t_max_yr": t, "n_output": 20, "backend": backend},
        }
    )


def ten(root="top"):
    tree = {
        root: ["ABCD", "E"],
        "ABCD": ["AB", "CD"],
        "AB": ["A", "B"],
        "CD": ["C", "D"],
        "A": ["m1", "m2"],
        "B": ["m3", "m4"],
        "C": ["m5", "m6"],
        "D": ["m7", "m8"],
        "E": ["m9", "m10"],
    }
    P = {"A": 2.0, "B": 3.0, "C": 2.5, "D": 4.0, "E": 3.3, "AB": 300, "CD": 400, "ABCD": 20000}
    P[root] = 800000
    return Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": tree,
            "masses": {f"m{k}": 1.0 for k in range(1, 11)},
            "orbits": {
                k: {"P": v, "i": 60.0 + 3 * j, "e": 0.1} for j, (k, v) in enumerate(P.items())
            },
            "run": {"t_max_yr": 0.3, "n_output": 10},
        }
    )


@pytest.mark.parametrize("topology", ["2+2+2", "3+3"])
def test_sextuple_presets_run(topology):
    from swingersim.backends.nbody import run_nbody
    from swingersim.backends.secular import run_secular

    cfg = six(topology)
    assert len(cfg.masses) == 6
    res = run_nbody(cfg)
    assert np.nanmax(res.energy_error) < 1e-6
    res = run_secular(six(topology, t=100.0, backend="secular"))
    assert set(res.orbits) == set(cfg.tree)


def test_sub_quadruples():
    root, _ = build_tree(six())
    kinds = sorted((q[0], q[1].name) for q in sub_quadruples(root))
    assert kinds == [("2+2", "AB"), ("2+2", "out"), ("3+1", "A"), ("3+1", "B")]
    root, _ = build_tree(six("3+3"))
    assert sorted(q[0] for q in sub_quadruples(root)) == ["2+2", "3+1", "3+1"]


def test_ten_stars():
    from swingersim.backends.nbody import run_nbody

    cfg = ten()
    root, masses = build_tree(cfg)
    assert len(masses) == 10
    from swingersim.hierarchy import initialize

    pos, vel = initialize(root, masses)
    assert len(nested_triples(root, masses, pos, vel)) >= 8
    res = run_nbody(cfg)
    quads = res.system["quadruples"]
    assert [q["core"] for q in quads] == ["AB", "CD"]
    assert "eclipse_fraction_all" in res.system
    assert res.hierarchy[-1].count("[") == 9  # nine orbits, one bound hierarchy


def test_ghost_test_custom_root():
    from swingersim.stability import ghost_test

    g = ghost_test(six(), n_outer=0.01)
    assert g.stable


@pytest.mark.skipif(mlp_directory() is None, reason="Vynatheya+23 classifiers not installed")
def test_mlp_sub_quadruples():
    v = mlp_verdicts(six())
    assert len(v) == 4 and all(isinstance(ok, bool) for _, ok in v)


def test_ensemble_sextuple(tmp_path):
    from swingersim.ensemble import report, run_ensemble

    cfg = six(t=50.0, backend="secular")
    data = cfg.model_dump()
    data["ensemble"] = {
        "n": 3,
        "seed": 2,
        "vary": {"orbits.C.i": "isotropic", "orbits.out.Omega": "angle"},
    }
    rows = run_ensemble(Config.model_validate(data), tmp_path / "e")
    assert len(rows) == 3 and all(r["fate"] != "error" for r in rows)
    assert "C" in report(rows)


def test_viewer_sextuple(tmp_path):
    pytest.importorskip("plotly")
    from swingersim.backends.nbody import run_nbody
    from swingersim.io import save_run
    from swingersim.viewer import write_html

    cfg = six()
    save_run(tmp_path / "r.h5", cfg, run_nbody(cfg))
    assert write_html(tmp_path / "r.h5", tmp_path / "v.html", max_frames=5).exists()
