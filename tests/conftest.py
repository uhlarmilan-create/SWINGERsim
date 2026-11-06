import os

os.environ.setdefault("SWINGERSIM_NO_STDERR_FILTER", "1")  # pytest captures fd 2 itself

from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parents[1] / "configs" / "examples"


@pytest.fixture
def example_2p2():
    from swingersim.config import load_config

    return load_config(EXAMPLES / "example_2p2.yaml")


@pytest.fixture
def example_3p1():
    from swingersim.config import load_config

    return load_config(EXAMPLES / "example_3p1.yaml")
