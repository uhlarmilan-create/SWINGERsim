"""Benchmarks against published results (see swingersim.validation and docs/*/validation.md)."""

import pytest

from swingersim.validation import BENCHMARKS

pytestmark = pytest.mark.validation


def _param(b):
    marks = []
    if b.slow or b.manual:
        marks.append(pytest.mark.slow)
    if b.manual:
        marks.append(pytest.mark.skip(reason="hours of computing: swingersim validate " + b.id))
    if b.explanation:
        marks.append(pytest.mark.xfail(reason=b.explanation, strict=False))
    return pytest.param(b, marks=marks)


@pytest.mark.parametrize("bench", [_param(b) for b in BENCHMARKS], ids=[b.id for b in BENCHMARKS])
def test_benchmark(bench):
    obtained = bench.compute()
    deviation, passed = bench.check(obtained)
    assert passed, (
        f"{bench.id}: obtained {obtained:.6g} {bench.unit}, expected {bench.expected:.6g} "
        f"(deviation {deviation:+.3g}, tolerance {bench.tolerance:g}); source: {bench.source}"
    )


def test_every_benchmark_is_documented():
    for b in BENCHMARKS:
        assert b.source and b.link.startswith("http") and b.rationale, b.id
        assert b.tier in {"A", "B", "C", "D"}, b.id


def test_benchmark_detects_missing_physics():
    """Negative control: without relativity the Mercury benchmark must fail."""
    from swingersim.validation import mercury_perihelion

    bench = next(b for b in BENCHMARKS if b.id == "gr-mercury")
    obtained = mercury_perihelion(gr=False)
    assert abs(obtained) < 0.01
    assert not bench.check(obtained)[1]


def test_missing_optional_package_is_skipped(monkeypatch):
    """A benchmark whose optional package is missing is reported as skipped, not a crash
    (job 061 on VYVAR: torch not installed)."""
    from swingersim import validation as v

    b = next(x for x in v.BENCHMARKS if x.id == "disc-face-on-tau")

    def boom():
        raise ModuleNotFoundError("No module named 'torch'", name="torch")

    monkeypatch.setattr(v, "BENCHMARKS", [v.Benchmark(**{**b.__dict__, "compute": boom})])
    (o,) = v.run_benchmarks([b.id])
    assert o.skipped == "needs the optional package torch"
    assert v._verdict(o) == "SKIP"
