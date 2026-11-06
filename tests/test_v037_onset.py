"""physics.mass_transfer.dynamic_onset: a Roche-lobe overflow driven by the dynamics (here
Kozai-Lidov cycles of a triple) circularises the binary and continues the run (swingersim.merger)."""

import math
import warnings

import numpy as np
import pytest
import yaml

from swingersim import units
from swingersim.config import Config

warnings.simplefilter("ignore")

ZLK = """
topology: custom
hierarchy: {out: [A, m3], A: [m1, m2]}
masses: {m1: 1.0, m2: 0.8, m3: 1.0}
orbits:
  A: {a: 0.1, e: 0.01}
  out: {a: 1.5, e: 0.2, i: 87.0, omega: 30}
stellar_model: {kind: sse, age_myr: 1000, z: 0.02}
physics: {gr: true, evolution: true, mass_transfer: {dynamic_onset: angular_momentum}}
run: {t_max_yr: 300, n_output: 61, backend: secular}
"""


def _cfg(**over):
    d = yaml.safe_load(ZLK)
    for path, v in over.items():
        node = d
        *head, last = path.split(".")
        for k in head:
            node = node.setdefault(k, {})
        node[last] = v
    return Config.model_validate(d)


def _state(cfg):
    from swingersim.hierarchy import build_tree, initialize

    root, m = build_tree(cfg)
    pos, vel = initialize(root, m)
    return {"pos": pos, "vel": vel, "masses": m, "spins": None}


def _pair_L(state, i=0, j=1):
    m = state["masses"]
    mu = m[i] * m[j] / (m[i] + m[j])
    return mu * np.cross(state["pos"][j] - state["pos"][i], state["vel"][j] - state["vel"][i])


def _radius_au(cfg, k, t):
    from swingersim.backends.evolving import props_at
    from swingersim.hierarchy import build_tree

    return props_at(cfg, build_tree(cfg)[1], t)[k].radius_rsun * units.RSUN_AU


@pytest.mark.parametrize("mode", ["angular_momentum", "periastron"])
def test_circularisation(mode):
    """The binary at periastron of an orbit with e = 0.9: circular at a (1 - e_c^2), angular
    momentum vector of the pair kept (angular_momentum), or at a (1 - e_c) with the donor at its
    lobe and the transfer starting at once (periastron); e_c: the star touches its lobe."""
    from swingersim import masstransfer
    from swingersim.hierarchy import build_tree, initialize, node_elements
    from swingersim.merger import post_onset_config
    from swingersim.stellar import roche_lobe_eggleton

    cfg = _cfg(**{"orbits.A.e": 0.9, "physics.mass_transfer.dynamic_onset": mode})
    st = _state(cfg)
    new, info = post_onset_config(cfg, st, "A", 0, 1, 40.0, 300.0)
    assert new is not None, info
    a = 0.1
    e_c = 1.0 - _radius_au(cfg, 0, 40.0) / (roche_lobe_eggleton(1.0 / 0.8) * a)
    assert info["e_onset"] == pytest.approx(e_c, abs=1e-12) and e_c < 0.9
    o = new.orbits["A"]
    assert o.e == 0.0 and o.P is None
    if mode == "angular_momentum":
        assert o.a == pytest.approx(a * (1 - e_c**2), rel=1e-6)
        # the fill of the circular binary: R / (r_L a (1 - e_c^2)) = 1 / (1 + e_c)
        assert info["fill_after"] == pytest.approx(1.0 / (1.0 + e_c), rel=2e-3)
    else:
        assert o.a == pytest.approx(a * (1 - e_c), rel=2e-3)
        assert info["fill_after"] == pytest.approx(1.0, abs=1e-5)
        assert masstransfer.histories(new)["A"].t_onset == 0.0
    # the new configuration: the pair's angular momentum direction kept, its size that of the
    # circular orbit; the outer orbit unchanged
    root, m = build_tree(new)
    pos, vel = initialize(root, m)
    s_new = {"pos": pos, "vel": vel, "masses": m}
    L_old, L_new = _pair_L(st), _pair_L(s_new)
    assert np.allclose(L_new / np.linalg.norm(L_new), L_old / np.linalg.norm(L_old), atol=1e-12)
    mu = 0.8 / 1.8
    assert np.linalg.norm(L_new) == pytest.approx(mu * math.sqrt(units.G * 1.8 * o.a), rel=1e-10)
    if mode == "angular_momentum":  # = the angular momentum at the contact eccentricity
        assert np.linalg.norm(L_new) == pytest.approx(
            mu * math.sqrt(units.G * 1.8 * a * (1 - e_c**2)), rel=1e-6
        )
    out_old = node_elements(build_tree(cfg)[0], st["pos"], st["vel"], st["masses"])
    out_new = node_elements(root, pos, vel, m)
    for k in ("a", "e", "inc"):
        assert getattr(out_new, k) == pytest.approx(getattr(out_old, k), rel=1e-9, abs=1e-12)
    # ages advance with the run, the pair rotates synchronously
    assert new.stellar_model.age_myr == pytest.approx(1000.0 + 40e-6)
    assert new.stars["m1"].P_rot == "sync" and new.stars["m2"].obliquity == 0.0
    assert new.run.t_max_yr == pytest.approx(260.0)


def test_secular_zlk_onset_continues():
    """Kozai-Lidov cycles drive e of the inner binary to the contact of m1 with its lobe: the run
    continues to t_max with the circularised binary (P 8.6 d -> 1.3 d)."""
    from swingersim.backends.secular import run_secular

    cfg = _cfg()
    res = run_secular(cfg)
    assert res.event is None and res.times[-1] == pytest.approx(300.0)
    assert np.all(np.diff(res.times) > 0)
    (og,) = res.system["_onsets"]
    assert 55.0 < og["time_yr"] < 70.0  # first maximum of the cycle
    assert og["e_before"] == pytest.approx(og["e_onset"], abs=1e-6)  # stopped at the contact
    assert og["a_after_au"] == pytest.approx(0.1 * (1 - og["e_onset"] ** 2), rel=1e-9)
    assert og["fill_after"] == pytest.approx(1.0 / (1.0 + og["e_onset"]), rel=2e-3)
    assert res.system["onsets_n"] == 1 and res.system["onset_1_orbit"] == "A"
    k = int(np.searchsorted(res.times, og["time_yr"]))
    a, e = res.orbits["A"]["a"], res.orbits["A"]["e"]
    assert a[k - 1] == pytest.approx(0.1) and a[k + 1] == pytest.approx(og["a_after_au"], rel=1e-6)
    assert np.all(e[k + 1 :] < 0.01)  # general relativity quenches the cycles of the close pair
    a_out = res.orbits["out"]["a"]
    assert np.allclose(a_out, a_out[0], rtol=1e-6)  # the outer orbit is not touched
    tr = [t for t in res.transitions if t["kind"] == "mass_transfer_onset"]
    assert len(tr) == 1 and tr[0]["time_yr"] == og["time_yr"]


def test_overflow_without_dynamic_onset():
    """Without the option the overflow ends the secular run as before, now at the first contact
    of the star with its lobe, whose state is the last row of the run."""
    from swingersim.backends.secular import run_secular

    for stop in (False, True):
        res = run_secular(
            _cfg(
                **{"physics.mass_transfer.dynamic_onset": None, "run.stop_on_roche_overflow": stop}
            )
        )
        assert res.event == "roche_overflow" and "_onsets" not in res.system
        assert res.event_detail["stars"] == ["m1"]  # the star at its lobe (was [] in 0.37)
        assert res.event_time == pytest.approx(res.times[-1]) and 55.0 < res.event_time < 70.0
        a, e = res.orbits["A"]["a"][-1], res.orbits["A"]["e"][-1]
        r_l = 0.49 * (1.25 ** (2 / 3)) / (0.6 * 1.25 ** (2 / 3) + math.log1p(1.25 ** (1 / 3)))
        fill = _radius_au(_cfg(), 0, res.event_time) / (r_l * a * (1 - e))
        assert fill == pytest.approx(1.0, abs=1e-6)


def test_onset_not_continued():
    from swingersim.merger import post_onset_config

    cfg = _cfg(**{"orbits.A.e": 0.9, "run.t_max_yr": -300.0})
    new, info = post_onset_config(cfg, _state(cfg), "A", 0, 1, -40.0, -300.0)
    assert new is None and info["stop"] == "a backward run"
    cfg = _cfg(**{"orbits.A.e": 0.9})
    new, info = post_onset_config(cfg, _state(cfg), "A", 0, 1, 300.0, 300.0)
    assert new is None and "end of the run" in info["stop"]


@pytest.mark.slow  # ~70 s
def test_nbody_agrees_with_secular():
    """The same Kozai-Lidov onset in the N-body backend: onset time within 3 % of the secular
    one (double averaging), the circular orbit and the continuation the same."""
    from swingersim.backends.nbody import run_nbody
    from swingersim.backends.secular import run_secular

    sec = run_secular(_cfg()).system["_onsets"][0]
    res = run_nbody(_cfg(**{"run.backend": "nbody"}))
    assert res.event is None and res.times[-1] == pytest.approx(300.0)
    nb = res.system["_onsets"][0]
    assert nb["time_yr"] == pytest.approx(sec["time_yr"], rel=0.03)
    assert nb["e_onset"] == pytest.approx(sec["e_onset"], abs=2e-3)
    assert nb["a_after_au"] == pytest.approx(sec["a_after_au"], rel=5e-3)
    k = int(np.searchsorted(res.times, nb["time_yr"]))
    assert np.all(res.orbits["A"]["e"][k + 1 :] < 0.01)


def test_outputs_show_the_onset(tmp_path):
    """The run file of a continued run: the GUI digest shows a fact for the onset, does not count
    it as a change of the hierarchy, and the mass-transfer entries are those of the continued
    binary on the time axis of the whole run."""
    import subprocess
    import sys

    from swingersim.gui import digest as D

    (tmp_path / "zlk.yaml").write_text(ZLK)
    r = subprocess.run(
        [sys.executable, "-m", "swingersim", "run", "zlk.yaml", "-o", "zlk.h5", "--no-viewer"],
        cwd=tmp_path, capture_output=True, text=True, check=False,
    )  # fmt: skip
    assert r.returncode == 0, r.stderr[-2000:]
    assert "MASS TRANSFER STARTED BY THE DYNAMICS" in r.stdout
    d = D.run_digest(tmp_path / "zlk.h5")
    facts = {f["k"]: f for f in d["facts"]}
    f = facts["Mass transfer started by the dynamics (orbit A)"]
    assert f["v"].startswith("P 8.6") and "e = 0.89" in f["note"]
    assert facts["Hierarchy changes"]["v"] == "0"
    import h5py

    with h5py.File(tmp_path / "zlk.h5") as h:
        import json

        system = json.loads(h.attrs["system_json"])
    og = system["_onsets"][0]
    assert system["mass_transfer_A_t_onset_yr"] is None or (
        system["mass_transfer_A_t_onset_yr"] >= og["time_yr"]
    )
