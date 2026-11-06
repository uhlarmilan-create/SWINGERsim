"""Validation test 1 (energy conservation) and GR periastron advance."""

import time

import numpy as np
import pytest
import rebound
import reboundx

from swingersim import units
from swingersim.backends.nbody import run_nbody
from swingersim.diagnostics import gr_periastron_advance_rate
from swingersim.kepler import state_to_elements


def _binary(m1, m2, a, e, gr):
    sim = rebound.Simulation()
    sim.G = units.G
    sim.integrator = "ias15"
    sim.add(m=m1)
    sim.add(m=m2, a=a, e=e)
    sim.move_to_com()
    rebx = None
    if gr:
        rebx = reboundx.Extras(sim)
        f = rebx.load_force("gr_full")
        rebx.add_force(f)
        f.params["c"] = units.C
    return sim, rebx  # rebx must stay referenced while integrating


def _omega(sim, mu):
    p0, p1 = sim.particles[0], sim.particles[1]
    r = np.array([p1.x - p0.x, p1.y - p0.y, p1.z - p0.z])
    v = np.array([p1.vx - p0.vx, p1.vy - p0.vy, p1.vz - p0.vz])
    return state_to_elements(mu, r, v).omega


def test_binary_energy_newton():
    sim, _ = _binary(1.0, 0.5, 0.05, 0.3, gr=False)
    e0 = sim.energy()
    sim.integrate(10.0)  # ~1100 orbits
    assert abs(sim.energy() / e0 - 1) < 1e-10


@pytest.mark.slow
def test_gr_periastron_advance():
    m1, m2, a, e = 1.0, 0.5, 0.05, 0.3
    sim, _rebx = _binary(m1, m2, a, e, gr=True)  # keep rebx alive
    mu = units.G * (m1 + m2)
    T = 100.0
    w0 = _omega(sim, mu)
    sim.integrate(T)
    measured = np.unwrap([w0, _omega(sim, mu)])[1] - w0
    expected = gr_periastron_advance_rate(m1 + m2, a, e) * T
    assert measured == pytest.approx(expected, rel=0.01)


def test_quadruple_short_run(example_2p2):
    cfg = example_2p2.model_copy(deep=True)
    cfg.run.t_max_yr = 2.0
    cfg.run.n_output = 50
    res = run_nbody(cfg)
    assert res.event is None
    assert np.nanmax(res.energy_error) < 1e-5
    assert set(res.mutual) == {"out-A", "out-B", "A-B"}


def test_quadruple_gr_hamiltonian_conserved(example_2p2):
    """With GR the 1PN Hamiltonian is checked; it already contains the Newtonian energy.
    Adding sim.energy() on top (versions <= 0.13.1) left a 1PN-sized fluctuation (~1e-6 for
    compact eccentric binaries) and false energy_error stops."""
    cfg = example_2p2.model_copy(deep=True)
    cfg.orbits["B"].e = 0.6
    cfg.run.t_max_yr = 2.0
    cfg.run.n_output = 50
    res = run_nbody(cfg)
    assert res.event is None
    assert np.nanmax(res.energy_error) < 1e-9


def test_quadruple_newton_energy(example_3p1):
    cfg = example_3p1.model_copy(deep=True)
    cfg.physics.gr = False
    cfg.run.t_max_yr = 2.0
    cfg.run.n_output = 20
    res = run_nbody(cfg)
    assert np.nanmax(res.energy_error) < 1e-10


def test_collision_detected(example_2p2):
    cfg = example_2p2.model_copy(deep=True)
    cfg.orbits["B"].P = None
    cfg.orbits["B"].a = 0.005  # ~1 Rsun: stars overlap almost immediately
    cfg.run.t_max_yr = 0.1
    res = run_nbody(cfg)
    assert res.event == "collision"
    assert set(res.event_detail["stars"]) <= {"m3", "m4"}


@pytest.mark.parametrize("sign", [1.0, -1.0])
def test_collision_detected_in_both_directions(example_2p2, sign):
    """A plunging orbit collides forward and backward in time. REBOUND's "direct" search
    skips pairs that recede in forward time, so a backward run used to miss the collision and
    hang at the point-mass singularity (job 010); max_wall_s caps the test if it regresses."""
    cfg = example_2p2.model_copy(deep=True)
    cfg.orbits["B"].P = None
    cfg.orbits["B"].a = 0.05
    cfg.orbits["B"].e = 0.9999  # periastron 5e-6 AU, far inside the stars
    cfg.physics.tides = False
    cfg.run.t_max_yr = sign * 0.1
    cfg.run.max_wall_s = 30.0
    res = run_nbody(cfg)
    assert res.event == "collision"
    assert set(res.event_detail["stars"]) <= {"m3", "m4"}


def test_wall_time_stops_inside_a_step(example_2p2, monkeypatch):
    """max_wall_s is enforced after every integrator step, not only between output steps.
    The sim.heartbeat setter is not used: in some REBOUND releases it fails ('_hb' missing
    from Simulation.__slots__), which made every member of job 010 end with an error."""

    def broken_setter(self, func):
        raise AttributeError("'Simulation' object has no attribute '_hb'")

    monkeypatch.setattr(rebound.Simulation, "heartbeat", property(None, broken_setter))
    cfg = example_2p2.model_copy(deep=True)
    cfg.run.t_max_yr = 1000.0
    cfg.run.n_output = 2  # one output step: the whole run is a single integrate() call
    cfg.run.extremes_every_orbits = 1.0e9
    cfg.run.max_wall_s = 0.2
    t0 = time.perf_counter()
    res = run_nbody(cfg)
    assert res.event == "wall_time"
    assert time.perf_counter() - t0 < 30.0


def test_line_collision_search_works_backward():
    """Why backward runs use the "line" search: a head-on plunge integrated with dt < 0 is
    found by it (the "direct" search skips the pair and the step size collapses)."""
    sim = rebound.Simulation()
    sim.integrator = "ias15"
    sim.collision = "line"
    sim.add(m=1.0, r=0.01)
    sim.add(m=1e-3, r=0.01, x=1.0, vx=1.0)  # falls onto the star going backward in time
    sim.dt = -1e-3
    with pytest.raises(rebound.Collision):
        sim.integrate(-3.0)
    assert -0.6 < sim.t < -0.5
