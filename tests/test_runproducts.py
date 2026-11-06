"""Every run leaves its standard products next to the HDF5 file."""

import csv

import pytest

from swingersim.cli import main
from swingersim.runproducts import product_paths

CFG = (
    "topology: '2+2'\nmasses: {m1: 1.2, m2: 1.0, m3: 1.1, m4: 0.9}\n"
    "orbits: {A: {P: 2.0, i: 89.0}, B: {P: 3.0, i: 88.5}, out: {P: 200.0, e: 0.2, i: 89.0}}\n"
    "stars: {m1: {R: 1.2}, m2: {R: 1.0}, m3: {R: 1.1}, m4: {R: 0.9}}\n"
    "run: {t_max_yr: 0.3, n_output: 61}\n"
)


def test_run_writes_all_products(tmp_path):
    pytest.importorskip("matplotlib")
    cfg = tmp_path / "sys.yaml"
    cfg.write_text(CFG)
    h5 = tmp_path / "out" / "sys.h5"
    assert main(["run", str(cfg), "-o", str(h5)]) == 0
    p = product_paths(h5)
    for key in ("summary", "timeseries", "orbits", "inclinations", "conservation"):
        assert p[key].exists(), key
    assert not p["stars"].exists()  # no evolving stars
    assert h5.with_suffix(".html").exists()
    # the CSV starts with '#' lines (provenance, units, columns) since 0.36
    rows = list(csv.DictReader(ln for ln in p["timeseries"].open() if not ln.startswith("#")))
    assert len(rows) == 61 and {"t_yr", "e_A", "i_out_deg", "imut_A-B_deg"} <= set(rows[0])
    text = p["summary"].read_text()
    assert "sys_timeseries.csv" in text and "sys.h5" in text
    try:
        import plotly  # noqa: F401

        assert p["viewer"].exists() and p["viewer"].name in text
    except ImportError:
        assert "no 3D viewer" in text


def test_products_can_be_switched_off(tmp_path):
    cfg = tmp_path / "sys.yaml"
    cfg.write_text(CFG)
    h5 = tmp_path / "sys.h5"
    assert main(["run", str(cfg), "-o", str(h5), "--no-viewer", "--no-plots", "--no-html"]) == 0
    p = product_paths(h5)
    assert p["summary"].exists() and p["timeseries"].exists()
    assert not p["viewer"].exists() and not p["orbits"].exists()
