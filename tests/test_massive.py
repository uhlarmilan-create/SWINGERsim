"""Massive stars: MIST v1.2 up to 300 Msun and the structure constants of Claret (2004)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from swingersim.config import Config
from swingersim.mist import available, mist_track
from swingersim.stellar import claret2004_structure, claret_structure, star_properties


def test_grids_reach_300_msun():
    from swingersim.mist import _grid

    for feh, v in available():
        g = _grid(("p" if feh >= 0 else "m") + f"{abs(feh):.2f}", v)
        assert g["mass"][0] == pytest.approx(0.1) and g["mass"][-1] == pytest.approx(300.0)
        assert {"star_mdot", "surface_h1"} <= set(g["columns"])


@pytest.mark.parametrize("vvcrit", [0.0, 0.4])
def test_very_massive_track(vvcrit):
    tr = mist_track(150.0, 0.0, vvcrit)
    assert 2.4 < tr.t_ms < 3.2  # Myr
    zams = tr._age[201]
    assert tr.current_mass(zams) == pytest.approx(150.0, rel=0.01)
    assert tr.mdot(zams + 0.5) < -1e-6  # strong wind (Vink et al. 2001)
    end = tr.t_end * (1 - 1e-9)
    assert tr.current_mass(end) < 40.0
    assert tr.is_wolf_rayet(end)
    assert not tr.is_wolf_rayet(zams + 0.1)


def test_claret_2004_continues_2023_at_8_msun():
    """At 8 Msun the two model families give the same k2 and I/(MR^2) on the ZAMS within 2 %."""
    k23, b23 = claret_structure(8.0, 4.25, 0.0, "main")
    k04, b04 = claret2004_structure(8.0, 4.25, "main")
    assert k04 == pytest.approx(k23, rel=0.02)
    assert b04 == pytest.approx(b23, rel=0.03)


def test_claret_2004_clamps_above_125():
    assert claret2004_structure(300.0, 3.9) == claret2004_structure(125.8925, 3.9)


def _cfg(structure="auto", mass=60.0):
    return Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": mass, "m2": 30.0, "m3": 20.0, "m4": 10.0},
            "orbits": {
                "A": {"P": 10.0, "e": 0.1, "i": 85},
                "B": {"P": 12.0, "e": 0.1, "i": 80},
                "out": {"P": 5000.0, "e": 0.2, "i": 70},
            },
            "stellar_model": {
                "kind": "mist",
                "age_myr": 1.0,
                "vvcrit": 0.0,
                "structure": structure,
            },
            "run": {"t_max_yr": 1.0},
        }
    )


def test_mist_massive_star_properties():
    p = star_properties(_cfg(), "m1", 60.0)
    assert 10.0 < p.radius_rsun < 15.0
    assert 0.004 < p.apsidal_k2 < 0.02  # Claret (2004): log k2 ~ -1.8 to -2.3 on the MS
    assert p.teff_k > 40000.0


def test_mdot_follows_vink():
    from swingersim.validation import mist_mdot_over_vink

    for m in (20.0, 60.0):
        assert mist_mdot_over_vink(m, 0.3) == pytest.approx(1.0, abs=0.1)


def test_track_mass_loss_of_massive_star_in_secular_run():
    from swingersim.backends.secular import run_secular

    d = _cfg(mass=60.0).model_dump(exclude_none=True)
    d["orbits"]["A"]["P"] = 50.0  # keep the primary inside its Roche lobe
    d["physics"] = {"gr": True, "evolution": True, "mass_loss": {"source": "track"}}
    d["run"] = {"backend": "secular", "t_max_yr": 1.0e6, "n_output": 21}
    res = run_secular(Config.model_validate(d))
    m1 = res.stellar["m1"]["mass_msun"]
    ma = m1 + res.stellar["m2"]["mass_msun"]
    assert m1[-1] < m1[0] - 1.0
    am = res.orbits["A"]["a"] * ma
    assert np.max(np.abs(am / am[0] - 1.0)) < 1e-12
    assert math.isfinite(res.system["mass_loss_psi_max_out"])


def _qz(age, vvcrit=0.0, branch="lowest"):
    return Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 26.1, "m2": 32.3, "m3": 70.3, "m4": 8.8},
            "orbits": {
                "A": {"P": 5.999, "e": 0.1, "i": 85},
                "B": {"P": 20.7, "e": 0.34, "i": 80},
                "out": {"P": 14722.0, "e": 0.2, "i": 80},
            },
            "stellar_model": {"kind": "mist", "age_myr": age, "vvcrit": vvcrit},
            "physics": {
                "evolution": True,
                "mass_loss": {"source": "track", "masses": "current", "branch": branch},
            },
            "run": {"t_max_yr": 10.0, "n_output": 3},
        }
    )


def test_present_masses_are_inverted_to_initial_masses():
    from swingersim import massloss

    with pytest.warns(UserWarning, match="least evolved"):
        new = massloss.resolve_masses(_qz(2.5))
    for name, m_now in _qz(2.5).masses.items():
        tr = mist_track(new.masses[name], 0.0, 0.0)
        assert tr.current_mass(2.5) == pytest.approx(m_now, rel=1e-6)
    assert new.masses["m3"] > 80.0  # the 70 Msun supergiant has lost ~17 Msun by 2.5 Myr
    with pytest.raises(ValueError, match="reached by tracks"):
        massloss.resolve_masses(_qz(2.5, branch="unique"))
    with pytest.raises(ValueError, match="no MIST track"):
        massloss.resolve_masses(_qz(4.0))  # no 70 Msun star exists at 4 Myr


def test_mass_interpolation_is_unbiased_between_grid_masses():
    # masses are interpolated as fractions of the initial mass: no excess mass on the MS
    for m in (1.14, 8.79, 97.0):
        tr = mist_track(m, 0.0, 0.0)
        assert tr.current_mass(tr._age[250]) <= m * (1 + 1e-6)
