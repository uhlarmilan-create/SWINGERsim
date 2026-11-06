"""The interactive viewer for every kind of run: secular runs, mergers, dynamically started mass
transfer, and `swingersim view` over whole folders (viewers of earlier runs)."""

import subprocess
import sys
import warnings

import numpy as np
import pytest
import yaml

warnings.simplefilter("ignore")
pytest.importorskip("plotly")

COLLIDING = """
topology: '2+2'
masses: {m1: 1.0, m2: 0.8, m3: 0.9, m4: 0.7}
orbits:
  A: {P: 10.0, e: 0.95, i: 80, M0: 120}
  B: {P: 3.0, e: 0.0, i: 85}
  out: {P: 600.0, e: 0.2, i: 70, M0: 30}
stellar_model: {kind: mist, age_myr: 1000.0}
physics: {gr: true, mergers: {}}
run: {t_max_yr: 2.0, n_output: 41}
"""
SECULAR = """
topology: '2+2'
masses: {m1: 1.2, m2: 1.0, m3: 1.1, m4: 0.9}
orbits: {A: {P: 2.0, i: 89.0}, B: {P: 3.0, i: 88.5}, out: {P: 2000.0, e: 0.3, i: 70.0}}
run: {t_max_yr: 2000.0, n_output: 61, backend: secular}
"""


def _q(cwd, *args):
    r = subprocess.run(
        [sys.executable, "-m", "swingersim", *args],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )
    assert r.returncode == 0, r.stderr[-2000:]
    return r.stdout


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    d = tmp_path_factory.mktemp("v")
    (d / "sub").mkdir()
    (d / "merge.yaml").write_text(COLLIDING)
    (d / "sec.yaml").write_text(SECULAR)
    _q(d, "run", "merge.yaml", "-o", "merge.h5", "--no-viewer", "--no-plots")
    _q(d, "run", "sec.yaml", "-o", "sub/sec.h5", "--no-viewer", "--no-plots")
    return d


def _frames(fig, trace):
    return [f.data[trace] for f in fig.frames]


def test_merged_star_disappears_and_the_product_is_labelled(runs):
    from swingersim.viewer import build_figure

    fig = build_figure(runs / "merge.h5")
    sky = [t for t in fig.frames[-1].data if getattr(t, "mode", "") == "markers+text"]
    labels = list(sky[0].text)
    assert "1+2" in labels  # the product keeps the place of star 1
    k = labels.index("2")
    xs = np.asarray(sky[0].x, float)
    assert np.isnan(xs[k])  # star 2 is gone after the merger
    first = next(t for t in fig.frames[0].data if getattr(t, "mode", "") == "markers+text")
    assert list(first.text) == ["1", "2", "3", "4"] and np.all(np.isfinite(first.x))
    assert "merged at t" in fig.layout.annotations[2].text  # zoom panel of orbit A


def test_folder_view_writes_missing_viewers_only(runs):
    out = _q(runs, "view", ".")
    assert "2 viewers written" in out
    assert (runs / "merge_viewer.html").stat().st_size > 1e6
    assert (runs / "sub" / "sec_viewer.html").exists()
    out = _q(runs, "view", ".")
    assert "0 viewers written, 2 up to date" in out
    out = _q(runs, "view", ".", "--force")
    assert "2 viewers written" in out
    html = (runs / "sub" / "sec_viewer.html").read_text(encoding="utf-8")
    assert "drag the 3D view" in html and "sec" in html


def test_zoom_follows_an_orbit_that_shrinks(tmp_path):
    """A binary circularised by the dynamic onset of mass transfer shrinks fivefold: its zoom
    panel follows it."""
    from pathlib import Path

    text = (Path(__file__).parent / "test_v037_onset.py").read_text(encoding="utf-8")
    ZLK = text.split('ZLK = """')[1].split('"""')[0]
    from swingersim.viewer import build_figure

    cfg = yaml.safe_load(ZLK)
    (tmp_path / "zlk.yaml").write_text(yaml.safe_dump(cfg))
    _q(tmp_path, "run", "zlk.yaml", "-o", "zlk.h5", "--no-viewer", "--no-plots")
    fig = build_figure(tmp_path / "zlk.h5")
    r0 = fig.frames[0].layout.xaxis2.range[1]
    r1 = fig.frames[-1].layout.xaxis2.range[1]
    assert r0 / r1 > 3
    assert "circularised" in fig.layout.annotations[-1].text
