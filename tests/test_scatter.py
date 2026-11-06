"""Few-body scattering (swingersim.scatter)."""

import json
import math

import numpy as np
import pytest

from swingersim import units
from swingersim.scatter import (
    KMS_TO_AUYR,
    ScatterConfig,
    ScatterSpec,
    b_max,
    classify,
    critical_velocity_kms,
    draw_encounters,
    initial_state,
    run_encounter,
    run_scatter,
    summarize_scatter,
)
from swingersim.worker import build_command

BB = {
    "target": {"masses": [1.0, 1.0], "a_au": 1.0},
    "projectile": {"masses": [1.0, 1.0], "a_au": 1.0},
}


def test_spec_validation():
    with pytest.raises(ValueError):
        ScatterSpec(**BB, b_au=1.0)  # no velocity
    with pytest.raises(ValueError):
        ScatterSpec(**BB, v_inf_kms=1.0, b_au=1.0, b_max_au=2.0)  # two impact options
    with pytest.raises(ValueError):
        ScatterSpec(
            target={"masses": [1.0, 1.0]}, projectile={"masses": [1.0]}, v_inf_kms=1, b_au=0
        )


def test_critical_velocity_and_focusing():
    spec = ScatterSpec(**BB, v_inf_vc=1.0, r_peri_max_au=3.0)
    # mu v_c^2 / 2 = 2 G m^2 / (2 a), mu = 1  ->  v_c = sqrt(2 G) AU/yr
    assert critical_velocity_kms(spec) == pytest.approx(math.sqrt(2 * units.G) / KMS_TO_AUYR)
    v = math.sqrt(2 * units.G)
    assert b_max(spec) == pytest.approx(math.sqrt(9 + 2 * units.G * 4 * 3 / v**2))


def test_draws_are_reproducible_and_isotropic():
    spec = ScatterSpec(**BB, v_inf_vc=0.5, b_max_au=5.0, n=2000, seed=3)
    a, b = draw_encounters(spec), draw_encounters(spec)
    assert a == b
    cosi = np.array([math.cos(e["target"]["inc"]) for e in a])
    assert abs(cosi.mean()) < 0.05
    b2 = np.array([e["b_au"] for e in a]) ** 2
    assert np.mean(b2) / 25.0 == pytest.approx(0.5, abs=0.03)  # uniform in area


def test_initial_state_is_the_requested_hyperbola():
    spec = ScatterSpec(**BB, v_inf_kms=10.0, b_au=4.0)
    enc = draw_encounters(spec)[0]
    m, _, pos, vel, r0, _h = initial_state(spec, enc)
    assert np.allclose((m[:, None] * pos).sum(0), 0, atol=1e-12)
    assert np.allclose((m[:, None] * vel).sum(0), 0, atol=1e-12)
    rt, rp = pos[:2].mean(0), pos[2:].mean(0)
    vt, vp = vel[:2].mean(0), vel[2:].mean(0)
    r, v = rp - rt, vp - vt
    assert np.linalg.norm(r) == pytest.approx(r0)
    mu = units.G * 4
    v_inf = 10.0 * KMS_TO_AUYR
    assert 0.5 * v @ v - mu / np.linalg.norm(r) == pytest.approx(0.5 * v_inf**2)
    assert np.linalg.norm(np.cross(r, v)) == pytest.approx(4.0 * v_inf)
    assert r @ v < 0  # approaching


def test_classify():
    initial = ((0, 1), (2, 3))
    assert classify(((0, 1), (2, 3)), initial, False) == "flyby"
    assert classify(((0, 2), (1, 3)), initial, False) == "exchange"
    assert classify(((0, 1), 2, 3), initial, False) == "single ionization"
    assert classify(((0, 2), 1, 3), initial, False) == "exchange + ionization"
    assert classify((((0, 1), 2), 3), initial, False) == "triple + single"
    assert classify((((0, 1), (2, 3)),), initial, False) == "bound 2+2"
    assert classify((0, 1, 2, 3), initial, False) == "full ionization"
    assert classify(initial, initial, True) == "collision"


def test_wide_flyby_conserves_energy_and_keeps_binaries():
    spec = ScatterSpec(**BB, v_inf_vc=3.0, b_au=40.0)
    row = run_encounter(spec, draw_encounters(spec)[0])
    assert row["outcome"] == "flyby"
    assert row["status"] == "resolved"
    assert row["energy_error"] < 1e-9


def test_collisions_and_resume(tmp_path):
    tight = {"masses": [1.0, 1.0], "a_au": 0.05, "radii_rsun": [1.0, 1.0]}
    cfg = ScatterConfig(
        scatter={"target": tight, "projectile": tight, "v_inf_vc": 0.3, "b_max_au": 0.2, "n": 6}
    )
    rows = run_scatter(cfg, tmp_path)
    assert len(rows) == 6
    assert any(r["outcome"] == "collision" for r in rows)
    c = next(r for r in rows if r["outcome"] == "collision")
    assert len(c["collision_stars"]) == 2
    with pytest.raises(FileExistsError):
        run_scatter(cfg, tmp_path)
    cfg.scatter.n = 8
    rows = run_scatter(cfg, tmp_path, resume=True)
    assert sorted(r["id"] for r in rows) == list(range(8))
    lines = (tmp_path / "rows.jsonl").read_text().splitlines()
    assert len(lines) == 8 and all(json.loads(x)["id"] < 8 for x in lines)
    cfg.scatter.v_inf_vc = 0.5
    with pytest.raises(ValueError):
        run_scatter(cfg, tmp_path, resume=True)
    assert "collision" in summarize_scatter(cfg, rows)


def test_worker_accepts_scatter(tmp_path):
    argv, out = build_command(
        {
            "command": "scatter",
            "config": "c.yaml",
            "output": "o",
            "args": {"n": 10, "processes": 0},
        },
        tmp_path,
    )
    assert argv[3] == "scatter" and "-j" in argv and "--n" in argv
    assert out == tmp_path / "o"


def test_thermal_eccentricities_respect_roche_lobes():
    body = {"masses": [1.0, 1.0], "a_au": 0.05, "e_thermal": True, "radii_rsun": [1.0, 1.0]}
    spec = ScatterSpec(target=body, projectile=body, v_inf_kms=1.0, b_max_au=1.0, n=300)
    r_au = 1.0 * units.RSUN_AU
    for enc in draw_encounters(spec):
        for key in ("target", "projectile"):
            q = 0.05 * (1 - enc[key]["e"])
            assert q * 0.3789 > r_au  # Eggleton r_L(q = 1) = 0.3789
