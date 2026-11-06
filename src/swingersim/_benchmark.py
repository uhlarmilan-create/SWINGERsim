"""The benchmark record of swingersim.validation."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]  # the SWINGERsim folder (configs/, data/)


@dataclass(frozen=True)
class Benchmark:
    id: str
    tier: str  # A analytic result, B published numerical result, C other code, D observed system
    quantity: str
    unit: str
    expected: float
    tolerance: float
    relative: bool  # tolerance relative to expected (True) or absolute (False)
    source: str
    link: str
    rationale: str  # why this tolerance
    compute: Callable[[], float]
    slow: bool = False  # minutes of computing; run with `swingersim validate --all`
    manual: bool = False  # hours of computing; only run when its id is given explicitly
    explanation: str = ""  # an understood difference (input, definition, published inconsistency
    # or a limitation shared by the published models): a failure is reported as EXPLAINED

    def check(self, obtained: float) -> tuple[float, bool]:
        dev = obtained - self.expected
        if self.relative:
            dev /= abs(self.expected)
        return dev, abs(dev) <= self.tolerance
