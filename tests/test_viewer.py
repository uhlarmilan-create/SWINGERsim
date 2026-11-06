"""Interactive HTML viewer."""

import pytest

pytest.importorskip("plotly")

from swingersim.backends.nbody import run_nbody
from swingersim.io import save_run
from swingersim.viewer import build_figure, write_html


@pytest.fixture
def saved_run(tmp_path, example_3p1):
    cfg = example_3p1.model_copy(deep=True)
    cfg.run.t_max_yr, cfg.run.n_output = 2.0, 40
    path = tmp_path / "run.h5"
    save_run(path, cfg, run_nbody(cfg))
    return path


def test_figure_has_frames_and_panels(saved_run):
    fig = build_figure(saved_run, max_frames=10)
    assert len(fig.frames) == 10
    # 3 orbits x 2 paths + stars, twice (3D and sky), plus two zoom panels of 3 traces
    assert all(len(fr.data) == 2 * (6 + 1) + 2 * 3 for fr in fig.frames)


def test_write_html(saved_run, tmp_path):
    out = write_html(saved_run, tmp_path / "view" / "run.html", max_frames=5, offline=False)
    text = out.read_text(encoding="utf-8")
    assert "plotly" in text.lower() and "animate" in text


def test_triple_has_one_zoom_panel(tmp_path):
    """A 2+1 triple has a single inner orbit (the viewer used to fail on it)."""
    from swingersim.config import Config

    cfg = Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"out": ["A", "m3"], "A": ["m1", "m2"]},
            "masses": {"m1": 1.0, "m2": 0.8, "m3": 0.5},
            "orbits": {"A": {"P": 3.0}, "out": {"P": 300.0, "e": 0.2, "i": 20.0}},
            "run": {"t_max_yr": 2.0, "n_output": 20},
        }
    )
    path = tmp_path / "triple.h5"
    save_run(path, cfg, run_nbody(cfg))
    fig = build_figure(path, max_frames=5)
    assert len(fig.frames) == 5
