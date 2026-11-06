"""HTML reports of the commands (swingersim.report)."""

from pathlib import Path

import pytest

from swingersim.cli import main

EX = Path(__file__).resolve().parents[1] / "configs" / "examples"
CAT = Path(__file__).resolve().parents[1] / "configs" / "catalog"

pytest.importorskip("plotly")


def _check(path: Path, *texts: str) -> str:
    assert path.exists()
    page = path.read_text(encoding="utf-8")
    assert page.startswith("<!doctype html>")
    assert "Printed summary" in page
    for t in texts:
        assert t in page
    return page


def test_run_report_is_self_contained(tmp_path):
    cfg = tmp_path / "c.yaml"
    text = (EX / "example_2p2.yaml").read_text().replace("t_max_yr: 200", "t_max_yr: 2")
    cfg.write_text(text)
    out = tmp_path / "run.h5"
    assert main(["run", str(cfg), "-o", str(out)]) == 0
    page = _check(tmp_path / "run.html", "Orbital elements", "Configuration")
    assert "<script src=" not in page and "Plotly.newPlot" in page  # offline: plotly embedded


def test_lightcurve_report_and_no_html(tmp_path):
    out = tmp_path / "lc.csv"
    argv = [
        "lightcurve",
        str(CAT / "eps_aur.yaml"),
        "--start-day",
        "600",
        "--days",
        "200",
        "--cadence-min",
        "14400",
        "--band",
        "V",
        "-o",
        str(out),
    ]
    assert main(argv) == 0
    _check(tmp_path / "lc.html", "Light curve")
    (tmp_path / "lc.html").unlink()
    assert main([*argv, "--no-html"]) == 0
    assert not (tmp_path / "lc.html").exists()
    custom = tmp_path / "sub" / "mine.html"
    assert main([*argv, "--html", str(custom)]) == 0
    _check(custom)


def test_ensemble_report(tmp_path):
    out = tmp_path / "ens"
    assert (
        main(["ensemble", str(EX / "example_v994_ensemble.yaml"), "-o", str(out), "--n", "2"]) == 0
    )
    _check(out / "report.html", "Fates")


def test_validate_report(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert main(["validate", "disc-face-on-tau"]) == 0
    # a partial run must not overwrite the report of the full validation
    assert not (tmp_path / "output" / "validation" / "report.html").exists()
    _check(
        tmp_path / "output" / "validation" / "report-disc-face-on-tau.html",
        "disc-face-on-tau",
        "PASS",
    )
