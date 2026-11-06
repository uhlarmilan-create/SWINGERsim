import numpy as np
import pytest

from swingersim import units
from swingersim.kepler import Elements, elements_to_state, solve_kepler, state_to_elements


@pytest.mark.parametrize("e", [0.0, 0.1, 0.5, 0.9, 0.99])
def test_solve_kepler(e):
    for M in np.linspace(0, 2 * np.pi, 37):
        E = solve_kepler(M, e)
        assert abs(E - e * np.sin(E) - np.mod(M, 2 * np.pi)) < 1e-13


def _angle_diff(x, y):
    return abs((x - y + np.pi) % (2 * np.pi) - np.pi)


def test_roundtrip_random():
    rng = np.random.default_rng(1)
    mu = units.G * 2.0
    for _ in range(500):
        el = Elements(
            a=10 ** rng.uniform(-2, 3),
            e=rng.uniform(0.01, 0.95),
            inc=rng.uniform(0.05, np.pi - 0.05),
            Omega=rng.uniform(0, 2 * np.pi),
            omega=rng.uniform(0, 2 * np.pi),
            M=rng.uniform(0, 2 * np.pi),
        )
        back = state_to_elements(mu, *elements_to_state(mu, el))
        assert back.a == pytest.approx(el.a, rel=1e-11)
        assert back.e == pytest.approx(el.e, abs=1e-11)
        for k in ("inc", "Omega", "omega", "M"):
            assert _angle_diff(getattr(back, k), getattr(el, k)) < 1e-9, k


def test_period_sma_consistency():
    # Earth around the Sun: 1 yr <-> 1 AU (to the Gaussian-year accuracy)
    assert units.period_to_sma(1.0, 1.0) == pytest.approx(1.0, rel=1e-4)
    assert units.sma_to_period(units.period_to_sma(0.3, 2.5), 2.5) == pytest.approx(0.3, rel=1e-14)
