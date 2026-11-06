"""Eclipse timings into the past (swingersim timings --from-years) join the forward series."""

import csv

import numpy as np

from swingersim.cli import main
from swingersim.config import load_config
from swingersim.timing import eclipse_times

CFG = (
    "topology: '2+2'\nmasses: {m1: 1.2, m2: 1.0, m3: 1.1, m4: 0.9}\n"
    "orbits: {A: {P: 2.0, i: 89.0, M0: 40}, B: {P: 3.0, i: 88.5, M0: 200},"
    " out: {P: 200.0, e: 0.2, i: 89.0}}\n"
    "stars: {m1: {R: 1.2}, m2: {R: 1.0}, m3: {R: 1.1}, m4: {R: 0.9}}\n"
    "run: {t_max_yr: 0.3, n_output: 61}\n"
)


def test_backward_eclipses_continue_the_forward_ones(tmp_path):
    p = tmp_path / "s.yaml"
    p.write_text(CFG)
    cfg = load_config(p)
    fwd, pn = eclipse_times(cfg, 0.1, ["A"])
    back, _ = eclipse_times(cfg, -0.1, ["A"])
    assert all(e.t_yr < 0 for e in back) and all(e.t_yr > 0 for e in fwd)
    prim = np.array([e.t_yr for e in back + fwd if e.kind == "primary"])
    assert np.all(np.diff(prim) > 0)  # in time order
    gaps = np.diff(prim) / pn["A"]
    assert np.allclose(gaps, np.round(gaps), atol=0.01)  # one period across t = 0 as well


def test_cli_writes_observed_bjd(tmp_path):
    p = tmp_path / "s.yaml"
    p.write_text(CFG)
    out = tmp_path / "t.csv"
    rc = main(
        [
            "timings",
            str(p),
            "--from-years",
            "-0.05",
            "--years",
            "0.05",
            "--bjd-t0",
            "2459000.5",
            "-o",
            str(out),
            "--no-html",
        ]
    )
    assert rc == 0
    rows = list(csv.DictReader(out.open()))
    t = np.array([float(r["t_obs_yr"]) for r in rows])
    bjd = np.array([float(r["bjd_obs"]) for r in rows])
    assert t.min() < 0 < t.max()
    assert np.allclose(bjd, 2459000.5 + t * 365.25, atol=1e-5)
