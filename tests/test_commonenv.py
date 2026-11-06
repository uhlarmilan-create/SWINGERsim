import math
import warnings

import pytest

from swingersim import commonenv as ce
from swingersim import masstransfer as mt
from swingersim.config import Config


def test_white_dwarf_and_helium_star_radii():
    assert ce.r_white_dwarf(0.6) == pytest.approx(
        0.0115 * math.sqrt((1.44 / 0.6) ** (2 / 3) - (0.6 / 1.44) ** (2 / 3))
    )
    assert ce.r_white_dwarf(1.0) < ce.r_white_dwarf(0.5)
    assert 0.1 < ce.r_helium_zams(1.0) < 0.3
    assert ce.core_radius(0.4, 1.5, 3) == ce.r_white_dwarf(0.4)  # degenerate core
    assert ce.core_radius(0.6, 4.0, 3) == ce.r_helium_zams(0.6)


def test_energy_formalism_survive_and_merge():
    # wide orbit: the cores survive and a_f follows eq. 71 exactly
    o = ce.outcome(1.5, 0.3, 100.0, 0.02, 0.8, 0.8, 0.75, 0.75, 1.0, 1.0, 0.5)
    e_bind = -(1.5 * 1.2 / (100.0 * 0.00465047)) / 0.5
    e_i = -0.5 * 0.3 * 0.8 / 1.0
    assert o.outcome == "survive"
    assert o.a_f == pytest.approx(-0.5 * 0.3 * 0.8 / (e_i + e_bind), rel=1e-3)
    # tight orbit, weak alpha: the companion fills its lobe -> merger with mass between core
    # and total
    m = ce.outcome(1.5, 0.3, 10.0, 0.02, 0.8, 0.8, 0.75, 0.75, 0.05, 0.2, 0.5)
    assert m.outcome == "merge"
    assert 0.3 < m.m1_f <= 2.3


def test_impulsive_mass_loss_circular():
    a, e = ce.impulsive_orbit(2.0, 1.5, 1.0, 0.0, 0.3)
    assert a == pytest.approx(1.5) and e == pytest.approx(1 / 3)
    a, e = ce.impulsive_orbit(2.0, 1.0, 1.0, 0.0, 0.0)
    assert math.isinf(a)  # half the mass lost: unbound


def _cfg(backend):
    return Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"out": ["A", "m3"], "A": ["m1", "m2"]},
            "masses": {"m1": 1.8, "m2": 0.5, "m3": 1.0},
            "orbits": {"A": {"P": 200.0}, "out": {"P": 3650.0, "i": 10.0}},
            "stellar_model": {"kind": "sse", "age_myr": 1652.75, "z": 0.02},
            "physics": {"gr": False, "evolution": True, "mass_transfer": {"common_envelope": {}}},
            "run": {"backend": backend, "t_max_yr": 4000, "n_output": 41},
        }
    )


def test_history_continues_after_common_envelope():
    cfg = _cfg("nbody")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        h = mt.histories(cfg)["A"]
    assert h.event == "post_common_envelope" and h.ce["outcome"] == "survive"
    assert h.m[0, -1] == pytest.approx(h.ce["donor_core"])
    assert h.a[-1] == pytest.approx(h.ce["a_f_au"])
    assert mt.stop(cfg) is None  # the N-body backend follows it
    s = mt.star_state(cfg, 0, 3900.0)
    assert s.radius_rsun < 0.05  # the core is a compact point mass
    assert mt.stop(_cfg("secular"))[1].startswith("post_common_envelope")
