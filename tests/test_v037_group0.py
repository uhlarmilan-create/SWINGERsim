"""Group 0 of the test campaign: the parts no reviewer had verified (Vokrouhlicky 2016 J_AB
limit, Hamers & Lai 2017 beta) checked against the dynamics and their derivation, the
radial velocities next to t = 0 (emission time before the start), and smoke tests of the
command-line tools without tests so far."""

import math
import subprocess
import sys
import warnings

import numpy as np
import pytest
import yaml

from swingersim import units
from swingersim.config import Config

warnings.simplefilter("ignore")


def _quad(iA, iB, OB, backend="secular", t_max=400000.0):
    d = yaml.safe_load(
        f"""
topology: 2+2
masses: {{m1: 1.0, m2: 0.8, m3: 1.2, m4: 0.6}}
orbits:
  A: {{a: 0.5, e: 0.0, i: {iA}}}
  B: {{a: 0.7, e: 0.0, i: {iB}, Omega: {OB}}}
  out: {{a: 20.0, e: 0.2, i: 0}}
physics: {{gr: false}}
run: {{t_max_yr: {t_max}, n_output: 2001, backend: {backend}, secular_octupole: false}}
"""
    )
    return Config.model_validate(d)


def _integrals(cfg, res):
    """C1, C2 of Vokrouhlicky (2016), eqs. (26)-(36), along a run."""
    m = [cfg.masses[f"m{k}"] for k in (1, 2, 3, 4)]
    x, y = np.cos(res.mutual["out-A"]), np.cos(res.mutual["out-B"])
    z = np.cos(res.mutual["A-B"])
    MA, MB = m[0] + m[1], m[2] + m[3]
    MAB = MA + MB
    aA, aB, a, e = 0.5, 0.7, 20.0, 0.2
    eta = math.sqrt(1 - e * e)
    n, nA, nB = (math.sqrt(units.G * M / s**3) for M, s in ((MAB, a), (MA, aA), (MB, aB)))
    s1 = -0.75 / eta**3 * n * (n / nA) * MB / MAB
    s2 = -0.75 / eta**3 * n * (n / nB) * MA / MAB
    sp = -0.75 / eta**4 * n * (m[0] / MA) * (m[1] / MA) * (aA / a) ** 2
    spp = -0.75 / eta**4 * n * (m[2] / MB) * (m[3] / MB) * (aB / a) ** 2
    return sp * x**2 + spp * y**2, (s2 / spp) * x + (s1 / sp) * y + z


@pytest.mark.parametrize("iA,iB,OB", [(10, 30, 90), (5, 35, 60), (20, 30, 180)])
def test_vokrouhlicky_integrals_and_jab_limit(iA, iB, OB):
    """Circular A and B below the Kozai-Lidov angle: the first integrals C1, C2 behind the J_AB
    limit are conserved by the secular run (1e-8), the mutual inclination of A and B varies,
    and it never exceeds the analytic limit (an upper bound)."""
    from swingersim.backends.secular import run_secular
    from swingersim.diagnostics.zlk import jab_max_vokrouhlicky
    from swingersim.hierarchy import build_tree, initialize

    cfg = _quad(iA, iB, OB)
    root, m = build_tree(cfg)
    bound = jab_max_vokrouhlicky(root, m, *initialize(root, m))
    res = run_secular(cfg)
    assert res.orbits["A"]["e"].max() < 1e-4 and res.orbits["B"]["e"].max() < 1e-4
    c1, c2 = _integrals(cfg, res)
    assert np.ptp(c1) < 1e-8 * abs(c1.mean())
    assert np.ptp(c2) < 1e-7
    J = np.degrees(res.mutual["A-B"])
    assert np.ptp(J) > 0.5 or iA + iB == 50  # (20, 30, 180): coplanar nodes, J constant
    assert J.max() <= bound + 1e-6


@pytest.mark.slow  # ~40 s
def test_vokrouhlicky_integrals_in_nbody():
    """The same integrals in the N-body model: conserved to the size of the short-period
    terms, much less than the secular variation of the mutual inclination of A and B."""
    from swingersim.backends.nbody import run_nbody

    cfg = _quad(10, 30, 90, "nbody", 30000.0)
    res = run_nbody(cfg)
    c1, c2 = _integrals(cfg, res)
    z = np.cos(res.mutual["A-B"])
    assert np.ptp(c2) < 0.2 * np.ptp(z)
    assert np.ptp(c1) < 0.01 * abs(c1.mean())


def test_hamers_lai_beta_from_the_timescales():
    """beta (Hamers & Lai 2017) = (3/4) cos i x ratio of the ZLK time of the inner orbit to that
    of the perturbing orbit, t_ZLK = (P_out^2 / P_in) (M_out / m_pert) (1 - e_out^2)^(3/2),
    written out with the periods (3+1) and with both binaries in the outer orbit (2+2)."""
    from swingersim.diagnostics.zlk import system_diagnostics
    from swingersim.hierarchy import build_tree, initialize

    def t_zlk(p_in, p_out, m_tot_out, m_pert, e_out):
        return p_out**2 / p_in * m_tot_out / m_pert * (1 - e_out**2) ** 1.5

    def per(a, m):
        return units.sma_to_period(a, m)

    # 3+1: ((m1 m2) m3) m4
    d = yaml.safe_load(
        """
topology: 3+1
masses: {m1: 1.0, m2: 0.7, m3: 0.9, m4: 1.3}
orbits:
  inner: {a: 0.3, e: 0.1}
  mid: {a: 4.0, e: 0.3, i: 50}
  out: {a: 60.0, e: 0.4, i: 0}
run: {t_max_yr: 10, n_output: 11}
"""
    )
    from swingersim.backends.nbody import run_nbody

    cfg = Config.model_validate(d)
    res = run_nbody(cfg)
    names = list(cfg.orbits)
    inner, mid, out = names[0], names[1], names[2]
    m1, m2, m3, m4 = 1.0, 0.7, 0.9, 1.3
    t12 = t_zlk(per(0.3, m1 + m2), per(4.0, m1 + m2 + m3), m1 + m2 + m3, m3, 0.3)
    t23 = t_zlk(per(4.0, m1 + m2 + m3), per(60.0, m1 + m2 + m3 + m4), m1 + m2 + m3 + m4, m4, 0.4)
    i23 = next(x.i_mut0_deg for x in res.diagnostics if x.inner == mid)
    want = 0.75 * math.cos(math.radians(i23)) * t12 / t23
    assert res.system[f"beta_hamers_lai_{inner}"] == pytest.approx(want, rel=1e-9)
    assert out  # (the outermost orbit)

    # 2+2: both binaries in the same outer orbit
    cfg = _quad(10, 30, 90, "nbody", 10.0)
    root, m = build_tree(cfg)
    pos, vel = initialize(root, m)
    res = run_nbody(cfg)
    s = system_diagnostics(root, m, res.diagnostics, (pos, vel))
    tA = t_zlk(per(0.5, 1.8), per(20.0, 3.6), 3.6, 1.8, 0.2)  # A perturbed by B (1.8 Msun)
    tB = t_zlk(per(0.7, 1.8), per(20.0, 3.6), 3.6, 1.8, 0.2)
    iB = next(x.i_mut0_deg for x in res.diagnostics if x.inner == "B")
    assert s["beta_hamers_lai_A"] == pytest.approx(
        0.75 * math.cos(math.radians(iB)) * tA / tB, rel=1e-9
    )


def test_radial_velocities_next_to_t0():
    """Velocities observed within z/c after t = 0 were emitted before it: they come from the
    backward integration and the curve is smooth across t = 0 (before: the state at t = 0
    without the acceleration, a flat piece of the curve)."""
    from swingersim.backends.nbody import build_simulation
    from swingersim.config import load_config
    from swingersim.radvel import model_rv

    cfg = load_config("configs/squadra_test/v994_her_fit_rv.yaml")
    sim = build_simulation(cfg)[0]
    lag_d = sim.particles[0].z / units.C * units.DAYS_PER_YEAR
    assert lag_d > 0.005  # the case: m1 is behind the barycentre
    t = 2459010.0 + np.linspace(-0.05, 0.05, 41)
    v = model_rv(cfg, np.repeat("m1", len(t)), t, 2459010.0)
    d2 = np.diff(v, 2)
    assert np.all(np.abs(d2 - np.median(d2)) < 2e-3)  # smooth: no kink at t = 0
    assert np.all(np.abs(np.diff(v)) > 0.1)  # no flat piece


def _cli(tmp_path, *args):
    r = subprocess.run(
        [sys.executable, "-m", "swingersim", *args],
        capture_output=True, text=True, check=False, cwd=tmp_path,
    )  # fmt: skip
    assert r.returncode == 0, r.stderr[-2000:]
    return r.stdout


ROOT = __import__("pathlib").Path(__file__).resolve().parents[1]


def test_cli_rv(tmp_path):
    """swingersim rv: the velocity curve of m1 of V994 Her A over two days has the semi-amplitude
    K1 = (2 pi G / P)^(1/3) m2 sin i / M^(2/3) / sqrt(1 - e^2) of the A pair (the motion in
    the outer orbits changes it by < 1 %)."""
    import pandas as pd

    cfg = ROOT / "configs/squadra_test/v994_her_fit_rv.yaml"
    _cli(tmp_path, "rv", str(cfg), "--start", "2459010", "--end", "2459012.1", "--step", "0.005",
         "--stars", "m1", "-o", "rv", "--no-html")  # fmt: skip
    d = pd.read_csv(tmp_path / "rv" / "rv_model.csv")
    v = d["rv_m1"].to_numpy()
    k1_run = 0.5 * (v.max() - v.min())
    from swingersim.config import load_config

    c = load_config(cfg)
    o = c.orbits["A"]
    m1, m2 = c.masses["m1"], c.masses["m2"]
    P = o.P / units.DAYS_PER_YEAR
    a = units.period_to_sma(P, m1 + m2)
    k1 = 2 * math.pi * a * m2 / (m1 + m2) * math.sin(math.radians(o.i)) / P
    k1 /= math.sqrt(1 - o.e**2)
    k1 *= units.AU_KM / units.SECONDS_PER_YEAR
    assert k1_run == pytest.approx(k1, rel=0.02)


def test_cli_apsidal(tmp_path):
    """swingersim apsidal: the relativistic part equals 6 pi G M / (c^2 a (1 - e^2)) per orbit and
    the total is the sum of its parts."""
    import pandas as pd

    from swingersim.config import load_config

    cfg = ROOT / "configs/examples/example_2p2_tides.yaml"
    _cli(tmp_path, "apsidal", str(cfg), "-o", "a.csv", "--no-html")
    d = pd.read_csv(tmp_path / "a.csv").set_index("pair")
    c = load_config(cfg)
    for pair, (i, j) in (("A", ("m1", "m2")), ("B", ("m3", "m4"))):
        M = c.masses[i] + c.masses[j]
        P = d.loc[pair, "P_d"] / units.DAYS_PER_YEAR
        a = units.period_to_sma(P, M)
        e = d.loc[pair, "e"]
        gr = 6 * math.pi * units.G * M / (units.C**2 * a * (1 - e * e)) / P  # rad/yr
        assert d.loc[pair, "gr_deg_yr"] == pytest.approx(math.degrees(gr), rel=1e-6)
        parts = d.loc[pair, ["tidal_deg_yr", "rot_configured_deg_yr", "gr_deg_yr"]].sum()
        assert d.loc[pair, "classical_gr_deg_yr"] == pytest.approx(parts, rel=1e-9)
        assert d.loc[pair, "U_classical_gr_yr"] == pytest.approx(360 / parts, rel=1e-9)


def test_cli_applegate(tmp_path):
    """swingersim applegate: the energy ratio of Tian et al. (2009) from the printed dP/P."""
    from swingersim.magnetic import applegate_energy_ratio, dp_over_p_from_oc

    out = _cli(tmp_path, "applegate", "--mass", "0.8", "--radius", "0.8", "--teff", "5000",
               "--a", "5", "--pbin", "1.5", "--pmod", "20", "--oc", "60")  # fmt: skip
    dpp = dp_over_p_from_oc(60.0, 20.0)
    ratio = applegate_energy_ratio(0.8, 0.8, 5000.0, 5.0, 1.5, dpp, 20.0)
    # by hand: 0.233 M^3 R^-10 (T/6000)^-4 a^4 (dP [s])^2 / P_mod
    by_hand = 0.233 * 0.8**3 * 0.8**-10 * (5000 / 6000) ** -4 * 5**4 * (dpp * 1.5 * 86400) ** 2 / 20
    assert ratio == pytest.approx(by_hand, rel=1e-12)
    assert f"{ratio:.3g}" in out and "not possible" in out


def test_cli_stability_escape_pulsation(tmp_path):
    out = _cli(tmp_path, "stability", str(ROOT / "configs/examples/example_2p2.yaml"),
               "--no-nbody", "--no-html")  # fmt: skip
    assert "Mardling-Aarseth" in out and "Vynatheya+22" in out
    out = _cli(tmp_path, "escape-theory", "--masses", "1", "1", "1", "--a", "1", "--v-inf", "1",
               "--no-html")  # fmt: skip
    probs = [float(x) for x in __import__("re").findall(r"probability ([\d.]+)", out)]
    assert len(probs) == 3 and sum(probs) == pytest.approx(1.0, abs=2e-3)  # equal masses: 1/3
    assert all(p == pytest.approx(1 / 3, abs=2e-3) for p in probs)
    out = _cli(tmp_path, "pulsation", str(ROOT / "configs/examples/example_pulsating_triple.yaml"))
    assert "delta_scuti" in out


def test_fit_holds_one_node_when_all_are_free():
    """A rotation of the whole system on the sky changes none of the fitted data: with the node
    of every orbit free the node of the outermost orbit is held (QX Cas fit: errors of 6.7e5
    deg and a spurious 1e-7 Msun error of m3 before)."""
    from swingersim.config import load_config
    from swingersim.dynfit import hold_one_node

    cfg = load_config("configs/squadra_test3/qx_cas_fit.yaml")
    spec = cfg.fit.stage2
    assert "orbits.inner.Omega" in spec.free and "orbits.out.Omega" in spec.free
    new, note = hold_one_node(cfg, spec, log=lambda *a: None)
    assert "orbits.out.Omega" not in new.free and "orbits.inner.Omega" in new.free
    assert "held at its configured value" in note
    one = spec.model_copy(update={"free": [n for n in spec.free if n != "orbits.inner.Omega"]})
    same, note = hold_one_node(cfg, one, log=lambda *a: None)
    assert same.free == one.free and note is None  # one node free: measured against the other


def test_model_failure_is_detected():
    from swingersim.dynfit import _model_failed

    ctx = {"obs": [1, 2], "rv": None, "inc": [], "lc": []}
    assert _model_failed({"t": None, "rv": None, "inc": None, "lc": None}, ctx)
    assert not _model_failed({"t": np.array([1.0, np.nan]), "rv": None, "inc": None,
                              "lc": None}, ctx)  # fmt: skip


def test_utc_conversion_far_ahead_and_offline():
    """Predicted minima decades ahead in geocentric UTC: the conversion needs no IERS table of
    UT1 - UTC (V994 Her fit_auto failed on the PC: 'IERS_Auto predictive values more than 30
    days old'); the round trip BJD_TDB -> UTC -> BJD_TDB is exact."""
    from astropy.utils import iers

    from swingersim.ocfit import from_bjd_tdb, to_bjd_tdb

    bjd = np.array([2440000.3, 2459010.7, 2470000.0, 2480000.0])
    with iers.conf.set_temp("auto_download", False), iers.conf.set_temp("auto_max_age", 1.0):
        utc = from_bjd_tdb(bjd, 256.6, 24.2)
        back = to_bjd_tdb(utc, "JD_UTC", 256.6, 24.2)
    assert np.all(np.abs(back - bjd) * 86400 < 1e-4)
    assert np.all(np.abs(utc - bjd) * 86400 < 600)  # Romer delay < 8.3 min, TDB - UTC ~ 1 min


@pytest.mark.slow  # ~2 min
def test_streamed_overflow_names_the_star():
    """An evolving secular run ended by a Roche-lobe overflow names the star (the extremes are
    recomputed up to the first contact, where the filling is 1 only to the bisection tolerance;
    the event listed no star in the test campaign)."""
    from swingersim.backends.secular import run_secular
    from swingersim.config import load_config

    res = run_secular(load_config("configs/catalog/hd150136.yaml"))
    assert res.event == "roche_overflow" and res.event_detail["stars"] == ["m1"]
    assert res.roche["m1"]["fill_max"] >= 1.0


def test_accreting_star_spins_at_most_at_break_up():
    """An accretor gains the angular momentum of a Keplerian orbit at its surface and spins at
    most at the break-up rate (it was spun up by the wind factor (m_new/m_old)^(2/(3 rg2)),
    beyond break-up, which made the secular run of an Algol with tides crawl)."""
    from swingersim import units as u
    from swingersim.backends import evolving as ev
    from swingersim.stellar import StarProperties

    old = new = StarProperties(radius_rsun=2.0, apsidal_k2=0.01, rg2=0.1)
    w0 = 100.0
    w = ev.scale_spin(w0, old, new, 1.0, 1.01)
    r = 2.0 * u.RSUN_AU
    j = 0.1 * 1.0 * r**2 * w0 + 0.01 * math.sqrt(u.G * 1.01 * r)
    assert w == pytest.approx(j / (0.1 * 1.01 * r**2), rel=1e-9)
    w_crit = math.sqrt(u.G * 2.0 / r**3)
    assert ev.scale_spin(w0, old, new, 1.0, 2.0) == pytest.approx(w_crit)
    assert ev.scale_spin(w0, old, new, 1.0, 1.0) == w0
    assert ev.scale_spin(w0, old, new, 1.0, 0.9) < w0


@pytest.mark.parametrize("topology", ["triple", "2+2"])
def test_locked_spins_turn_with_their_orbit(topology):
    """Aligned spins locked to their binary: the torque of the outer orbit turns the binary
    and its spins together (Correia et al. 2011, eqs. 8-15, adiabatic limit), so the total
    angular momentum (orbits + spins) is conserved even with a spin angular momentum of
    ~10 % of the binary's (it was not: the spins followed the orbit without a reaction)."""
    from swingersim.backends.secular import run_secular

    if topology == "triple":
        d = yaml.safe_load(
            """
topology: custom
hierarchy: {out: [A, m3], A: [m1, m2]}
masses: {m1: 1.0, m2: 1.0, m3: 1.0}
orbits:
  A: {a: 0.05, e: 0.05}
  out: {a: 1.0, e: 0.1, i: 50.0}
stars: {m1: {R: 3.0, P_rot: 0.4, lag_s: 1.0e-3}, m2: {R: 1.0, lag_s: 1.0e-3}}
physics: {gr: false, tides: true, rotation: true}
run: {t_max_yr: 200, n_output: 101, backend: secular}
"""
        )
    else:
        d = yaml.safe_load(
            """
topology: 2+2
masses: {m1: 1.0, m2: 1.0, m3: 1.0, m4: 1.0}
orbits:
  A: {a: 0.05, e: 0.05}
  B: {a: 0.06, e: 0.05, i: 30}
  out: {a: 1.0, e: 0.1, i: 50.0}
stars: {m1: {R: 3.0, P_rot: 0.4, lag_s: 1.0e-3}, m2: {R: 1.0, lag_s: 1.0e-3},
        m3: {R: 1.0, lag_s: 1.0e-3}, m4: {R: 1.0, lag_s: 1.0e-3}}
physics: {gr: false, tides: true, rotation: true}
run: {t_max_yr: 200, n_output: 101, backend: secular}
"""
        )
    res = run_secular(Config.model_validate(d))
    assert res.dissipative
    assert np.max(res.angmom_error) < 1e-8
    assert np.degrees(np.ptp(res.orbits["A"]["inc"])) > 5.0  # the binary does turn


@pytest.mark.slow  # ~15 s
def test_locked_spins_precess_as_the_nbody_spin_vectors():
    """A binary with a spin angular momentum of ~10 % of its orbit, tilted 50 deg to the
    tertiary: the secular node precession of the binary with the spins locked to it agrees
    with the N-body run of REBOUNDx tides_spin (full spin vectors with their reaction on the
    orbit) to 1 deg over 200 yr (without the locking: 13 deg)."""
    from swingersim.backends.nbody import run_nbody
    from swingersim.backends.secular import run_secular

    d = yaml.safe_load(
        """
topology: custom
hierarchy: {out: [A, m3], A: [m1, m2]}
masses: {m1: 1.0, m2: 1.0, m3: 1.0}
orbits:
  A: {a: 0.05, e: 0.05}
  out: {a: 1.0, e: 0.1, i: 50.0}
stars: {m1: {R: 3.0, P_rot: 0.4, lag_s: 1.0e-3}, m2: {R: 1.0, lag_s: 1.0e-3}}
physics: {gr: false, tides: true, rotation: true}
run: {t_max_yr: 200, n_output: 101, backend: secular}
"""
    )
    se = run_secular(Config.model_validate(d))
    d["run"]["backend"] = "nbody"
    nb = run_nbody(Config.model_validate(d))
    di = np.degrees(np.abs(se.orbits["A"]["inc"] - nb.orbits["A"]["inc"]))
    assert np.max(di) < 1.0


@pytest.mark.parametrize("topology", ["triple", "2+2"])
def test_locked_spins_without_dissipation(topology):
    """Without tidal dissipation the spins are constant but still locked to their binary:
    total angular momentum conserved (example_2p2_tides, secular: was 1.7e-4)."""
    from swingersim.backends.secular import run_secular
    from swingersim.config import load_config

    if topology == "2+2":
        cfg = load_config(ROOT / "configs/examples/example_2p2_tides.yaml")
        cfg.run.backend = "secular"
    else:
        cfg = Config.model_validate(
            yaml.safe_load(
                """
topology: custom
hierarchy: {out: [A, m3], A: [m1, m2]}
masses: {m1: 1.0, m2: 1.0, m3: 1.0}
orbits:
  A: {a: 0.05, e: 0.05}
  out: {a: 1.0, e: 0.1, i: 50.0}
stars: {m1: {R: 3.0, P_rot: 0.4}, m2: {R: 1.0}}
physics: {gr: false, tides: true, rotation: true}
run: {t_max_yr: 200, n_output: 101, backend: secular}
"""
            )
        )
    res = run_secular(cfg)
    assert not res.dissipative
    assert np.max(res.angmom_error) < 1e-10


def test_gui_fit_job_of_stage_one(tmp_path):
    """The GUI adds --processes 0 (all cores) to every fit job: a stage-1 fit must accept it
    (it failed with 'applies to stage 2'); an explicit number of processes is still refused."""
    cfg = str(ROOT / "configs/examples/example_fit.yaml")
    ok = subprocess.run([sys.executable, "-m", "swingersim", "fit", cfg, "--processes", "0", "-o",
                         str(tmp_path / "a"), "--no-html"], capture_output=True, text=True, check=False)  # fmt: skip
    assert ok.returncode == 0, ok.stdout + ok.stderr
    bad = subprocess.run([sys.executable, "-m", "swingersim", "fit", cfg, "--processes", "2", "-o",
                          str(tmp_path / "b")], capture_output=True, text=True, check=False)  # fmt: skip
    assert bad.returncode == 2 and "stage 2" in bad.stdout


def test_gui_serves_files_by_path(tmp_path):
    """/files/<path>: an HTML report in an iframe finds the PNG next to it (relative link; the
    GUI showed a broken image), and the path stays inside the project."""
    import threading
    import urllib.error
    import urllib.request

    from swingersim.gui.server import make_server

    d = tmp_path / "output" / "g6"
    d.mkdir(parents=True)
    (d / "lc.html").write_text('<img src="lc.png">')
    (d / "lc.png").write_bytes(b"\x89PNG fake")
    srv = make_server(tmp_path, "127.0.0.1", 0)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{srv.server_address[1]}"
    try:
        assert urllib.request.urlopen(base + "/files/output/g6/lc.png").read() == b"\x89PNG fake"
        assert b"lc.png" in urllib.request.urlopen(base + "/files/output/g6/lc.html").read()
        for bad in ("/files/..%2F..%2Fetc%2Fpasswd", "/files/../../etc/passwd"):
            with pytest.raises(urllib.error.HTTPError):
                urllib.request.urlopen(base + bad)
    finally:
        srv.shutdown()
        srv.server_close()


def test_linear_ephemeris_across_a_gap():
    """Eclipses stop for 800 cycles (the orbit precesses out of eclipse) and come back: the
    epochs after the gap are counted with the mean period of the eclipses before it, not
    with the nominal one (TIC 470710327 over 10 yr: O-C 595 min instead of 32 min)."""
    from swingersim.timing import linear_ephemeris

    P_true, P_nom = 1.10539, 1.104686
    E = np.r_[np.arange(0, 2000), np.arange(2800, 3300)]
    t = 0.3 + P_true * E
    _T0, P, Ef = linear_ephemeris(t, P_nom)
    assert P == pytest.approx(P_true, rel=1e-12) and np.array_equal(Ef, E)


@pytest.mark.parametrize("topology", ["triple", "2+2"])
def test_backward_tides_hold_fast_spins_pseudo_synchronous(topology):
    """Integrated into the past, tidal dissipation drives a spin away from synchronism
    (deviation ~ exp(|t| / t_sync)): over 100 t_sync the run failed or the spin ran away. Spins
    with t_sync < |t_max| / 3 are held pseudo-synchronous (swingersim.tides.slaved_spin)."""
    from swingersim.backends.secular import run_secular

    if topology == "triple":
        d = yaml.safe_load(
            """
topology: custom
hierarchy: {out: [A, m3], A: [m1, m2]}
masses: {m1: 1.0, m2: 0.9, m3: 0.8}
orbits:
  A: {P: 2.0, e: 0.0}
  out: {P: 3650.0, e: 0.2, i: 30.0}
stars: {m1: {lag_s: 50.0}, m2: {lag_s: 50.0}}
stellar_model: {kind: sse, age_myr: 100.0, z: 0.02}
physics: {gr: true, tides: true, rotation: true}
run: {backend: secular, t_max_yr: -3.0e5, n_output: 61}
"""
        )
    else:
        d = yaml.safe_load(
            """
topology: 2+2
masses: {m1: 1.0, m2: 0.9, m3: 1.0, m4: 0.8}
orbits:
  A: {P: 2.0, e: 0.0}
  B: {P: 2.6, e: 0.0, i: 10}
  out: {P: 3650.0, e: 0.2, i: 30.0}
stars: {m1: {lag_s: 50.0}, m2: {lag_s: 50.0}, m3: {lag_s: 50.0}, m4: {lag_s: 50.0}}
stellar_model: {kind: sse, age_myr: 100.0, z: 0.02}
physics: {gr: true, tides: true, rotation: true}
run: {backend: secular, t_max_yr: -3.0e5, n_output: 61}
"""
        )
    res = run_secular(Config.model_validate(d))
    p = res.spins["m1"]["P_rot_d"]
    assert np.all(np.isfinite(p)) and np.max(np.abs(p / p[0] - 1.0)) < 1e-6
    assert np.nanmax(res.angmom_error) < 1e-8


@pytest.mark.slow  # ~2 min
def test_backward_braking_is_the_mirror_of_forward():
    """Magnetic braking of a tidally locked binary: the period grows into the past at the rate
    it shrinks into the future (g5_mb_triple: +0.0132 d / -0.0131 d over 50 Myr)."""
    from swingersim.backends.secular import run_secular

    out = {}
    for sgn in (1, -1):
        cfg = Config.model_validate(
            yaml.safe_load(
                f"""
topology: custom
hierarchy: {{out: [A, m3], A: [m1, m2]}}
masses: {{m1: 1.0, m2: 0.9, m3: 0.8}}
orbits:
  A: {{P: 2.0, e: 0.0}}
  out: {{P: 3650.0, e: 0.2, i: 30.0}}
stellar_model: {{kind: sse, age_myr: 100.0, z: 0.02}}
physics: {{gr: true, tides: true, rotation: true, tidal_lag: hurley2002, evolution: true,
          magnetic_braking: {{}}}}
run: {{backend: secular, t_max_yr: {sgn * 2.0e7}, n_output: 21}}
"""
            )
        )
        r = run_secular(cfg)
        a = r.orbits["A"]["a"]
        out[sgn] = (a[-1] - a[0]) / a[0]
        assert np.max(np.abs(r.spins["m1"]["P_rot_d"] / r.spins["m1"]["P_rot_d"][0] - 1)) < 0.02
    assert out[1] < 0 < out[-1]
    assert out[-1] == pytest.approx(-out[1], rel=0.05)
