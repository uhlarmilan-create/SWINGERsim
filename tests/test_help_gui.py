"""Help system of swingersim gui (help.yaml, assist.py) and the own-MESA-grid structure option."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest
import yaml

from swingersim.config import Config
from swingersim.gui import assist

APP = Path(assist.__file__).with_name("app.html")


def test_help_yaml_complete():
    h = assist.help_texts()
    for k in (
        "tabs",
        "config",
        "compute",
        "questions",
        "templates",
        "tools",
        "results",
        "charts",
        "glossary",
        "warnings",
        "job_errors",
    ):
        assert k in h, k
    # every item of "What to compute" has a help entry
    app = APP.read_text()
    block = app[app.index("const COMPUTE = [") : app.index("const Compute = {")]
    ids = set(re.findall(r'\{ id: "(\w+)", label:', block)) - {"exact", "observer"}
    for i in ids:
        key = i if i != "oc" else "oc_exact"
        assert key in h["compute"], i
    assert "oc_observer" in h["compute"]
    # every top-level section of the configuration has an entry
    for k in Config.model_json_schema()["properties"]:
        assert k in h["config"], k
    # every tool of the Tools page has an entry
    for t in ("applegate", "escape-theory", "validate", "mesa-grid", "mesa-reduce"):
        assert t in h["tools"], t
    # questions tick existing computations
    for q in h["questions"]:
        for t in q["ticks"]:
            assert t in h["compute"], t
    # warnings and job errors are valid regular expressions
    for w in h["warnings"] + h["job_errors"]:
        re.compile(w["match"])


def test_templates_and_examples_are_valid():
    h = assist.help_texts()
    for t in h["templates"]:
        Config.model_validate(yaml.safe_load(t["yaml"]))
    base = yaml.safe_load(h["templates"][0]["yaml"])
    for key, e in h["config"].items():
        if "example" in e and "." not in key:
            d = dict(base)
            ex = yaml.safe_load(e["example"])
            if key in ("hierarchy", "masses", "orbits"):
                continue  # these replace the whole system
            d.update(ex)
            if key == "fit":
                continue
            Config.model_validate(d)


def test_quickcheck_and_methods(tmp_path):
    h = assist.help_texts()
    d = yaml.safe_load(h["templates"][0]["yaml"])
    q = assist.quickcheck(d)
    assert q["ok"]
    kinds = " ".join(i["k"] for i in q["items"])
    assert "stability" in kinds and "pair inner" in kinds
    bad = dict(d, orbits={"inner": {"P": 3.5, "i": 88}, "out": {"P": 5.0, "e": 0.3}})
    qb = assist.quickcheck(bad)
    assert any(i["tone"] == "bad" for i in qb["items"])
    assert not assist.quickcheck({"masses": "x"})["ok"]
    p = tmp_path / "c.yaml"
    p.write_text(h["templates"][0]["yaml"])
    m = assist.methods_text(p)
    assert "REBOUND" in m["text"] and any("Rein" in r for r in m["references"])


def test_applegate_from_configuration():
    from swingersim import magnetic as mg

    d = yaml.safe_load(assist.help_texts()["templates"][0]["yaml"])
    r = assist.applegate(d, "m2", 11.0, 30.0)
    v = r["values"]
    assert v["pair"] == "inner" and abs(v["pbin_d"] - 3.5) < 1e-6
    assert r["max_oc_s"] == pytest.approx(
        mg.applegate_max_oc(v["mass"], v["radius"], v["teff"], v["a_rsun"], v["pbin_d"], 11.0)
    )
    assert r["possible"] == (r["energy_ratio"] < 1)
    with pytest.raises(ValueError):
        assist.applegate(d, "m9", 11.0)


def _fake_grid(root: Path, name: str = "t1") -> Path:
    gdir = root / "mesa" / name
    gdir.mkdir(parents=True)
    rows = {k: [] for k in ("initial_mass", "z", "age_myr", "apsidal_k2", "rg2", "center_h1")}
    for z in (0.0142, 0.02):
        for m in (0.8, 1.0, 1.2):
            for a in np.linspace(0, 10000, 51):
                rows["initial_mass"].append(m)
                rows["z"].append(z)
                rows["age_myr"].append(a)
                rows["apsidal_k2"].append(0.01 * m + 1e-6 * a + (0.001 if z == 0.02 else 0.0))
                rows["rg2"].append(0.07 + 0.01 * m)
                rows["center_h1"].append(0.7 if a < 40 else 0.7 - 6e-5 * a)
    np.savez(gdir / f"{name}.npz", **{k: np.array(v) for k, v in rows.items()})
    return gdir / f"{name}.npz"


def test_mesa_structure(tmp_path, monkeypatch):
    from swingersim.stellar import find_mesa_grid, mesa_structure, star_properties

    p = _fake_grid(tmp_path)
    k2, rg2 = mesa_structure(p, 1.0, 5000.0, 0.0142)
    assert k2 == pytest.approx(0.01 + 0.005, rel=1e-3) and rg2 == pytest.approx(0.08)
    k2b, _ = mesa_structure(p, 1.1, 0.0, 0.0142)  # between two tracks (log M)
    assert 0.010 < k2b < 0.012
    k2z, _ = mesa_structure(p, 1.0, 5000.0, 0.019)  # nearest metallicity
    assert k2z == pytest.approx(0.016, rel=1e-3)
    with pytest.raises(ValueError):
        mesa_structure(p, 2.0, 100.0, 0.0142)
    monkeypatch.chdir(tmp_path)
    assert find_mesa_grid("t1") == p.resolve() or find_mesa_grid("t1").samefile(p)
    with pytest.raises(ValueError):
        find_mesa_grid("nothing")
    d = {
        "topology": "custom",
        "hierarchy": {"inner": ["m1", "m2"]},
        "masses": {"m1": 1.0, "m2": 0.9},
        "orbits": {"inner": {"P": 3.0}},
        "run": {"t_max_yr": 1},
        "stellar_model": {"kind": "mist", "age_myr": 3000, "feh": 0.0, "structure": "mesa:t1"},
    }
    cfg = Config.model_validate(d)
    sp = star_properties(cfg, "m1", 1.0)
    assert sp.apsidal_k2 == pytest.approx(0.01 + 0.003, rel=1e-3)
    from pydantic import ValidationError

    with pytest.raises(ValidationError):
        Config.model_validate({**d, "stellar_model": {"kind": "zams", "structure": "mesa:"}})
    with pytest.raises(ValidationError):
        Config.model_validate({**d, "stellar_model": {"kind": "zams", "structure": "nonsense"}})


def test_mesa_status_and_save(tmp_path):
    r = assist.mesa_save(
        tmp_path, {"name": "g1", "masses": {"from": 0.8, "to": 1.2, "step": 0.2}, "z": [0.0142]}
    )
    assert r["saved"] == "configs/mesa/g1.yaml"
    st = assist.mesa_status(tmp_path)
    g = st["grids"][0]
    assert (
        g["name"] == "g1"
        and g["tracks"] == 3
        and g["tracks_with_profiles"] == 0
        and not g["reduced"]
    )
    with pytest.raises(ValueError):
        assist.mesa_save(tmp_path, {"name": "g1", "masses": [1.0]})
    with pytest.raises(ValueError):
        assist.mesa_save(tmp_path, {"name": "bad name", "masses": [1.0]})
    _fake_grid(tmp_path, "g1")
    gd = assist.mesa_grid(tmp_path, "g1")
    assert gd["records"] == 2 * 3 * 51 and "apsidal_k2" in gd["columns"]


def test_escape_compare_needs_three_stars(tmp_path):
    (tmp_path / "scatter.yaml").write_text(
        yaml.safe_dump(
            {
                "config": {
                    "scatter": {
                        "target": {"masses": [1, 1], "a_au": 1},
                        "projectile": {"masses": [1, 1], "a_au": 0.5},
                    }
                }
            }
        )
    )
    (tmp_path / "rows.jsonl").write_text("")
    r = assist.escape_compare(tmp_path, tmp_path)
    assert not r["available"] and "3 stars" in r["why"]


def test_sibling_param_table(tmp_path):
    from swingersim.gui.digest import _param_table

    (tmp_path / "stage2_params.csv").write_text(
        "parameter,start,value,error\nmasses.m3,1,0.98,0.03\n"
    )
    assert _param_table(tmp_path / "stage2_params.csv") == {"masses.m3": [0.98, 0.03]}
    assert _param_table(tmp_path / "none.csv") == {}
