"""Dust discs and rings as occulters (swingersim.discs)."""

import math

import numpy as np
import pytest
from pydantic import ValidationError

from swingersim.config import Config, load_config
from swingersim.discs import build_discs, path_length, transmission
from swingersim.hierarchy import build_tree, initialize
from swingersim.photometry import StarLight, light_curve
from swingersim.validation_discs import (
    disc_face_on_transmission,
    disc_inclined_transmission,
    disc_strip_transmission,
    strip_exact,
)

EDGE_ON = np.array([0.0, 1.0, 0.0])


def _binary(**disc):
    return Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"A": ["m1", "m2"]},
            "masses": {"m1": 1.0, "m2": 1.0},
            "orbits": {"A": {"P": 100.0, "i": 90.0}},
            "stars": {"m1": {"R": 1.0, "Teff": 5800}, "m2": {"R": 0.5, "Teff": 4000}},
            "discs": {"d": {"host": "m2", "r_out_au": 0.1, "thickness_au": 0.01, **disc}},
            "physics": {"gr": False},
            "run": {"t_max_yr": 1.0},
        }
    )


def test_path_length_geometry():
    n_face = np.array([0.0, 0.0, 1.0])
    assert path_length([0.0], [0.0], 10.0, n_face, 1.0, 0.0, 0.1)[0] == pytest.approx(0.2)
    i = math.radians(60.0)
    n60 = np.array([0.0, math.sin(i), math.cos(i)])
    assert path_length([0.0], [0.0], 10.0, n60, 5.0, 0.0, 0.1)[0] == pytest.approx(0.4)
    edge = path_length([0.0, 0.5, 1.5], [0.0, 0.05, 0.0], 10.0, EDGE_ON, 1.0, 0.0, 0.1)
    assert edge == pytest.approx([2.0, 2.0 * math.sqrt(0.75), 0.0])
    # inner hole and a star inside the disc (only the part in front counts)
    assert path_length([0.0], [0.0], 10.0, EDGE_ON, 1.0, 0.5, 0.1)[0] == pytest.approx(1.0)
    assert path_length([0.0], [0.0], 0.0, EDGE_ON, 1.0, 0.0, 0.1)[0] == pytest.approx(1.0)
    # behind the star: nothing
    assert path_length([0.0], [0.0], -5.0, EDGE_ON, 1.0, 0.0, 0.1)[0] == 0.0


def test_transmission_against_exact_results():
    assert disc_strip_transmission() == pytest.approx(strip_exact(), abs=5e-3)
    assert disc_face_on_transmission() == pytest.approx(math.exp(-0.7), abs=1e-12)
    assert disc_inclined_transmission() == pytest.approx(math.exp(-1.4), abs=1e-12)


def test_config_checks():
    with pytest.raises(ValidationError, match="r_in_au"):
        _binary(r_in_au=0.2)
    with pytest.raises(ValidationError, match="needs relative_to"):
        _binary(precession_period_yr=10.0)
    with pytest.raises(ValidationError, match="no star or orbit"):
        _binary(host="m7")
    with pytest.raises(ValidationError, match="no orbit"):
        _binary(relative_to="X")


def test_precession_keeps_the_tilt():
    cfg = _binary(relative_to="A", tilt=20.0, precession_period_yr=-50.0)
    d = build_discs(cfg, build_tree(cfg)[0])[0]
    for t in (0.0, 7.0, 25.0, 50.0):
        n = d.normal(t)
        assert np.linalg.norm(n) == pytest.approx(1.0)
        assert math.degrees(math.acos(n @ d.axis)) == pytest.approx(20.0)
    assert d.normal(50.0) == pytest.approx(d.normal0)  # one full period
    assert not np.allclose(d.normal(25.0), d.normal0)


def test_disc_eclipse_in_light_curve():
    """A disc around the faint star m2 dims m1 at conjunction; without the disc the curve only
    shows the stellar eclipse."""
    cfg = _binary(r_out_au=0.1, thickness_au=0.01, i=90.0, Omega=0.0)
    no = _binary(r_out_au=0.1, thickness_au=0.01)
    no.discs = {}
    t = np.linspace(0.0, 100.0, 801)
    with_disc = light_curve(cfg, t, band="V", proximity=False, light_time=False)
    without = light_curve(no, t, band="V", proximity=False, light_time=False)
    extra = without.flux - with_disc.flux
    assert extra.max() > 0.05  # the disc blocks part of m1 in front of it
    assert np.all(extra > -1e-12)  # a disc only removes light
    # the dimming lasts longer than the stellar eclipse (the disc is wider than m2)
    assert (extra > 1e-3).sum() > ((1.0 - without.flux) > 1e-3).sum()


def test_gg_tau_ring_does_not_occult():
    cfg = load_config("configs/catalog/gg_tau.yaml")
    root, masses = build_tree(cfg)
    pos, _ = initialize(root, masses)
    lights = [
        StarLight(f"m{k + 1}", 0.005, 1.0, np.zeros(4), 0.3, 1.0, 0.0) for k in range(len(masses))
    ]
    tr = transmission(build_discs(cfg, root), 0.0, pos, masses, lights)
    assert np.all(tr == 1.0)
