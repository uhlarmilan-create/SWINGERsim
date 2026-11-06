"""Regressions found on the lab computer (Python 3.14, 2+2+1 ensemble)."""

import argparse

import pytest

from swingersim.cli import main
from swingersim.diagnostics.hierarchy import classify_change, label


def _subcommands():
    from swingersim import cli

    captured = {}
    orig = argparse.ArgumentParser.parse_args

    def grab(self, args=None, namespace=None):
        captured["p"] = self
        raise SystemExit(0)

    argparse.ArgumentParser.parse_args = grab
    try:
        with pytest.raises(SystemExit):
            cli.main([])
    finally:
        argparse.ArgumentParser.parse_args = orig
    sub = next(a for a in captured["p"]._actions if isinstance(a, argparse._SubParsersAction))
    return list(sub.choices)


@pytest.mark.parametrize("cmd", _subcommands())
def test_help_strings_format(cmd, capsys):
    # argparse expands '%' in help strings; a bare '%' crashes (Python 3.14 already at start-up)
    with pytest.raises(SystemExit) as e:
        main([cmd, "--help"])
    assert e.value.code == 0


def test_label_five_stars():
    t = ((0, 1), ((2, 3), 4))
    assert label((t,)) == "3+2"
    assert label((((0, 1), (2, 3)), 4)) == "quadruple + single"
    assert label((((0, 1), ((2, 3), 4)),)) == "3+2"
    assert classify_change(((((0, 1), (2, 3)), 4),), (((0, 1), ((2, 3), 4)),)) == "reconfiguration"


def test_retry_errors(tmp_path):
    import json

    from swingersim.config import Config
    from swingersim.ensemble import run_ensemble

    cfg = Config.model_validate(
        {
            "topology": "2+2",
            "masses": {"m1": 1.0, "m2": 1.0, "m3": 1.0, "m4": 1.0},
            "orbits": {"A": {"P": 10.0}, "B": {"P": 12.0}, "out": {"P": 2000.0, "e": 0.2}},
            "run": {"backend": "secular", "t_max_yr": 1000.0, "n_output": 50},
            "ensemble": {"n": 3, "seed": 1, "vary": {"orbits.out.e": {"uniform": [0.1, 0.3]}}},
        }
    )
    out = tmp_path / "ens"
    run_ensemble(cfg, out)
    rows_file = out / "rows.jsonl"
    rows = [json.loads(x) for x in rows_file.read_text().splitlines()]
    rows[1]["fate"], rows[1]["error"] = "error", "KeyError: 5"
    rows_file.write_text("".join(json.dumps(r) + "\n" for r in rows))
    rows = run_ensemble(cfg, out, resume=True)  # plain resume keeps the error
    assert sum(r["fate"] == "error" for r in rows) == 1
    rows = run_ensemble(cfg, out, resume=True, retry_errors=True)
    assert sum(r["fate"] == "error" for r in rows) == 0
    assert len(rows) == 3


def test_worker_retry_errors_flag(tmp_path):
    from swingersim.worker import build_command

    argv, _ = build_command(
        {
            "command": "ensemble",
            "config": "c.yaml",
            "output": "o",
            "args": {"resume": True, "retry_errors": True},
        },
        tmp_path,
    )
    assert "--retry-errors" in argv and "--resume" in argv


def test_central_rates():
    from swingersim.ensemble import member_value, observation_weights

    rows = [
        {"member": 0, "direction": "forward", "B_dtiltdt0_deg_yr": -0.1, "B_e0": 0.1},
        {"member": 0, "direction": "backward", "B_dtiltdt0_deg_yr": -2.7, "B_e0": 0.1},
        {"member": 1, "direction": "forward", "B_dtiltdt0_deg_yr": -1.4, "B_e0": 0.1},
    ]
    runs = {"forward": {0: rows[0], 1: rows[2]}, "backward": {0: rows[1]}}
    assert member_value(runs, 0, "B_dtiltdt0_deg_yr") == pytest.approx(-1.4)  # central
    assert member_value(runs, 0, "forward:B_dtiltdt0_deg_yr") == -0.1
    assert member_value(runs, 1, "B_dtiltdt0_deg_yr") == -1.4  # one direction only
    assert member_value(runs, 0, "B_e0") == 0.1
    w = observation_weights(rows, {"B_dtiltdt0_deg_yr": {"value": -1.4, "sigma": 0.2}})
    assert w[0] == pytest.approx(1.0) and w[1] == pytest.approx(1.0)


def test_rate_window():
    from swingersim.backends.secular import run_secular
    from swingersim.config import Config
    from swingersim.ensemble import summarize

    base = {
        "topology": "2+2",
        "masses": {"m1": 1.0, "m2": 1.0, "m3": 1.0, "m4": 1.0},
        "orbits": {
            "A": {"P": 5.0, "i": 85.0},
            "B": {"P": 6.0, "i": 40.0, "Omega": 30.0},
            "out": {"P": 200.0, "e": 0.2, "i": 70.0},
        },
        "run": {"backend": "secular", "t_max_yr": 40.0, "n_output": 401},
    }
    rates = []
    for win in (None, 0.5, 10.0):
        d = dict(base, run=dict(base["run"], rate_window_yr=win))
        cfg = Config.model_validate(d)
        rates.append(summarize(cfg, run_secular(cfg))["B_didt0_deg_yr"])
    assert rates[0] != rates[2]  # the window changes the fitted rate
    assert abs(rates[1] - rates[0]) < abs(rates[2] - rates[0]) + 1e-12
