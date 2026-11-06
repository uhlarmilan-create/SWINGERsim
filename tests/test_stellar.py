"""Stellar model, tides and Roche lobes (internal consistency; literature checks are benchmarks)."""

from __future__ import annotations

import numpy as np
import pytest

from swingersim.backends.nbody import run_nbody
from swingersim.config import Config
from swingersim.stellar import (
    polytrope_constants,
    star_properties,
    zams_properties,
    zams_radius_tout,
)


def _cfg(**over) -> Config:
    base = {
        "topology": "2+2",
        "masses": {"m1": 1.0, "m2": 0.8, "m3": 0.7, "m4": 0.6},
        "orbits": {
            "A": {"P": 3.0, "e": 0.2},
            "B": {"P": 2.0, "e": 0.1, "i": 10.0},
            "out": {"P": 300.0, "e": 0.3, "i": 20.0},
        },
        "physics": {"gr": True, "tides": True, "rotation": True},
        "run": {"t_max_yr": 2.0, "n_output": 20},
    }
    base.update(over)
    return Config.model_validate(base)


def test_polytrope_n1_analytic():
    k2, rg2 = polytrope_constants(1.0)
    assert k2 == pytest.approx((15.0 - np.pi**2) / (2.0 * np.pi**2), rel=1e-8)
    assert rg2 == pytest.approx(2.0 / 3.0 * (1.0 - 6.0 / np.pi**2), rel=1e-8)


def test_tout_radius_sensible():
    assert zams_radius_tout(1.0) == pytest.approx(0.89, rel=0.02)  # ZAMS Sun
    masses = np.linspace(0.1, 20.0, 200)
    radii = [zams_radius_tout(m) for m in masses]
    assert np.all(np.diff(radii) > 0)
    assert zams_properties(0.2).apsidal_k2 > zams_properties(2.0).apsidal_k2  # n=1.5 vs n=3


def test_overrides_table_and_python(tmp_path):
    table = tmp_path / "grid.csv"
    table.write_text("mass_msun,radius_rsun,apsidal_k2,rg2\n0.5,0.5,0.02,0.1\n1.5,1.5,0.01,0.08\n")
    cfg = _cfg(stellar_model={"kind": "table", "file": str(table)}, stars={"m2": {"R": 3.0}})
    p1 = star_properties(cfg, "m1", 1.0)
    assert p1.radius_rsun == pytest.approx(1.0) and p1.apsidal_k2 == pytest.approx(0.015)
    assert star_properties(cfg, "m2", 0.8).radius_rsun == 3.0  # per-star value wins
    module = tmp_path / "mymodel.py"
    module.write_text(
        "def star(m):\n    return {'radius_rsun': 2*m, 'apsidal_k2': 0.005, 'rg2': 0.07}\n"
    )
    cfg = _cfg(stellar_model={"kind": "python", "entry": f"{module}:star"})
    assert star_properties(cfg, "m3", 0.7).radius_rsun == pytest.approx(1.4)


def test_tides_conserve_energy_and_angular_momentum():
    res = run_nbody(_cfg(physics={"gr": False, "tides": True, "rotation": True}))
    assert not res.dissipative
    assert np.nanmax(res.energy_error) < 1e-9
    assert np.nanmax(res.angmom_error) < 1e-9
    assert res.spins["m1"]["P_rot_d"][0] == pytest.approx(3.0, rel=1e-6)  # synchronous


def test_dissipation_conserves_angular_momentum():
    stars = {k: {"lag_s": 100.0} for k in ("m1", "m2", "m3", "m4")}
    stars["m1"]["P_rot"] = 1.0
    res = run_nbody(_cfg(stars=stars, physics={"gr": False, "tides": True, "rotation": True}))
    assert res.dissipative
    assert np.nanmax(res.angmom_error) < 1e-7  # strong dissipation: stiffer spin equations
    p = res.spins["m1"]["P_rot_d"]
    assert p[-1] > p[0]  # the fast rotator is braked towards synchronism


def test_roche_overflow_stops_run():
    cfg = _cfg(
        stars={"m1": {"R": 9.0}},
        physics={"gr": False},
        run={"t_max_yr": 1.0, "n_output": 10, "stop_on_roche_overflow": True},
    )
    res = run_nbody(cfg)
    assert res.event == "roche_overflow"
    assert "m1" in res.event_detail["stars"]
    assert res.roche["m1"]["fill_max"] >= 1.0
