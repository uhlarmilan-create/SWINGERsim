import pytest

from swingersim import gw
from swingersim.validation_v023 import gw_hulse_taylor_pdot, gw_secular_circular_decay


def test_peters_rates_and_hulse_taylor():
    da, de = gw.peters_rates(1.4, 1.4, 0.01, 0.0)
    assert da < 0 and de == 0
    assert gw_hulse_taylor_pdot() == pytest.approx(-2.40263e-12, rel=1e-3)


def test_secular_backend_follows_peters_decay():
    assert gw_secular_circular_decay() == pytest.approx(1.0, rel=1e-4)
