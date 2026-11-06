"""Roche geometry and meshes, contact binaries, tertiary mass transfer, donor and
accretor response, rotation-corrected k2, pre-main-sequence tides, rigid and shining discs."""

import math
from pathlib import Path

import numpy as np
import pytest
import yaml

from swingersim.config import Config

ROOT = Path(__file__).resolve().parents[1]


def _binary(**kw):
    d = {
        "topology": "custom",
        "hierarchy": {"A": ["m1", "m2"]},
        "masses": {"m1": 1.0, "m2": 0.8},
        "orbits": {"A": {"P": 5.0, "i": 89.0}},
        "stars": {"m1": {"R": 1.0, "Teff": 5800}, "m2": {"R": 0.8, "Teff": 5000}},
        "physics": {"gr": False},
        "run": {"t_max_yr": 0.1},
    }
    d.update(kw)
    return Config.model_validate(d)


def test_roche_lobes_and_contact_radii():
    from swingersim import roche
    from swingersim.stellar import roche_lobe_eggleton

    for q in (0.2, 1.0, 4.0):
        assert roche.roche_lobe_volume_radius(q, 0) == pytest.approx(
            roche_lobe_eggleton(1.0 / q), rel=0.01
        )
    r0 = roche.contact_radii(0.3, 0.0)
    r1 = roche.contact_radii(0.3, 1.0)
    assert r1[0] > r0[0] and r1[1] > r0[1]
    assert roche.fillout_of(0.3, roche.contact_radii(0.3, 0.4)[0]) == pytest.approx(0.4, abs=0.02)


def test_rotating_star_shape():
    from swingersim import roche

    assert float(roche.rotating_radius(np.array([math.pi / 2]), 1.0)[0]) == pytest.approx(1.5)
    s = roche.rotating_surface(0.8, 32, 64)
    assert np.linalg.norm(s.points, axis=1).min() == pytest.approx(1.0, abs=0.01)  # near the pole
    g_pol = s.gravity[np.argmax(s.points[:, 2])]
    g_eq = s.gravity[np.argmin(np.abs(s.points[:, 2]))]
    assert g_eq < g_pol  # the equator is darker (von Zeipel)


def test_mesh_matches_spheres_and_rotation_changes_flux():
    from swingersim import mesh, photometry

    cfg = _binary()
    t = np.linspace(3.65, 3.85, 21)
    a = photometry.light_curve(cfg, t, "V", proximity=False, light_time=False)
    orig = mesh.needed
    mesh.needed = lambda c: True
    try:
        b = photometry.light_curve(cfg, t, "V", proximity=False, light_time=False)
    finally:
        mesh.needed = orig
    assert np.max(np.abs(b.flux - a.flux)) < 1e-3
    fast = _binary(
        stars={
            "m1": {"R": 1.0, "Teff": 5800},
            "m2": {"R": 0.8, "Teff": 5000, "shape": "rotating", "P_rot": 0.5},
        }
    )
    c = photometry.light_curve(fast, t, "V", proximity=False, light_time=False)
    assert np.all(np.isfinite(c.flux))
    assert np.max(np.abs(c.flux - a.flux)) > 1e-4


def test_contact_binary_light_curve():
    from swingersim import photometry

    cfg = _binary(
        masses={"m1": 1.2, "m2": 0.5},
        orbits={"A": {"P": 0.4775, "i": 80.0, "contact_fillout": 0.3}},
        stars={"m1": {"R": 1.0, "Teff": 6200}, "m2": {"R": 0.7, "Teff": 6000}},
        run={"t_max_yr": 0.01},
    )
    t = np.linspace(0.0, 0.4775, 49)
    lc = photometry.light_curve(cfg, t, "V", proximity=False, light_time=False)
    mag = -2.5 * np.log10(lc.flux / lc.flux.max())
    assert 0.3 < mag.max() < 1.0  # a W UMa light curve
    # one envelope, one temperature: both lobes have the same surface brightness
    assert lc.stars[0].weight == pytest.approx(lc.stars[1].weight)


def test_k2_rotation_correction():
    from swingersim.tides import k2_rotating

    cfg = _binary(physics={"gr": False, "tides": True, "rotation": True, "k2_rotation": True})
    assert k2_rotating(cfg, 0, 0.01, 0.0, 1.0, 1.0) == 0.01
    assert k2_rotating(cfg, 0, 0.01, 50.0, 1.0, 1.0) < 0.01
    off = _binary(physics={"gr": False, "tides": True, "rotation": True})
    assert k2_rotating(off, 0, 0.01, 50.0, 1.0, 1.0) == 0.01


def test_linear_convective_reduction_is_stronger():
    from swingersim.tides import hurley_k_over_t

    args = (1.0, 1.0, 0.03, 0.05, 1.0, 0.05, 0.01, 0.0)
    q = hurley_k_over_t(*args, reduction="quadratic")
    lin = hurley_k_over_t(*args, reduction="linear")
    assert lin > q > 0


def test_rigid_disc_polar_libration_and_rates():
    from scipy.integrate import solve_ivp

    from swingersim.discs import circumstellar_rate, rigid_rhs, ring_average_rate

    e = 0.6
    rg = {"lb": np.array([0.0, 0.0, 1.0]), "eb": np.array([e, 0.0, 0.0]), "w": 1.0, "e2": e * e}
    i = math.radians(80.0)
    sol = solve_ivp(
        lambda t, l: rigid_rhs(l, rg),
        (0, 100),
        [math.sin(i), 0, math.cos(i)],
        rtol=1e-9,
        atol=1e-12,
        max_step=0.05,
    )
    assert np.all(sol.y[0] > 0)  # librates about the eccentricity vector
    assert ring_average_rate(1, 1, 1, 3, 6, 1.0) > 0
    assert circumstellar_rate(1.0, 5.0, 6.0, 0.0, 0.1, 0.8, 1.0, 0.9) < 0  # retrograde


def test_disc_soft_edge_and_light():
    from swingersim.photometry import light_curve

    d = yaml.safe_load((ROOT / "configs/catalog/eps_aur.yaml").read_text())
    d["discs"]["disc"].update({"tau": 3.0, "edge_au": 0.2, "emission": True, "albedo": 0.5})
    cfg = Config.model_validate(d)
    lc = light_curve(cfg, np.array([100.0, 900.0]), "V")
    assert lc.components["disc_disc"][0] > 0
    assert lc.flux[1] < lc.flux[0]
    with pytest.raises(ValueError):
        d["discs"]["disc"]["tau"] = None
        Config.model_validate(d)


def test_rigid_disc_config():
    from swingersim.discs import build_discs
    from swingersim.hierarchy import build_tree

    d = yaml.safe_load((ROOT / "configs/catalog/kh15d.yaml").read_text())
    d["discs"]["ring"].pop("precession_period_yr")
    d["discs"]["ring"]["dynamics"] = "rigid"
    cfg = Config.model_validate(d)
    root, _ = build_tree(cfg)
    disc = build_discs(cfg, root)[0]
    n0, n1 = disc.normal(0.0), disc.normal(1000.0)
    assert np.linalg.norm(n1) == pytest.approx(1.0)
    assert not np.allclose(n0, n1)


def test_contact_angular_momentum_loss():
    from swingersim import units
    from swingersim.masstransfer import contact_jdot, orbital_j, outer_lobe_radius

    a = units.period_to_sma(0.35 / units.DAYS_PER_YEAR, 1.5)
    t = orbital_j(1.0, 0.5, a) / abs(contact_jdot(1.0, 0.5, 1.0, 0.7, a))
    assert 1e9 < t < 5e9
    assert outer_lobe_radius(0.5, 0) > 0.4


def test_zeta_ad_giant():
    from swingersim.masstransfer import zeta_ad_giant

    assert zeta_ad_giant(1.0, 0.0) == pytest.approx(-1.0 / 3.0, abs=1e-6)
    assert zeta_ad_giant(1.0, 0.5) > zeta_ad_giant(1.0, 0.2)


def test_tertiary_mass_transfer_giant_triple():
    from swingersim import masstransfer

    d = yaml.safe_load((ROOT / "configs/examples/example_tertiary_transfer.yaml").read_text())
    cfg = Config.model_validate(d)
    hs = masstransfer.histories(cfg)
    h = hs["out"]
    assert isinstance(h, masstransfer.TertiaryHistory)
    assert h.t_onset is not None
    assert h.m[1:, -1].sum() > h.m[1:, 0].sum() - 0.01  # the inner binary accretes
    assert "out" in masstransfer.orbit_a_at(cfg, 0.0) and "A" in masstransfer.orbit_a_at(cfg, 0.0)


def test_wolf_rayet_flag_beyond_initial_track():
    """0.27.1: a donor that outlives its initial MIST track is not flagged (no exception)."""
    from swingersim.backends.evolving import _wolf_rayet
    from swingersim.config import StellarModelSpec
    from swingersim.stellar import stellar_track

    tr = stellar_track(StellarModelSpec(kind="mist", age_myr=1000.0), 2.0)
    assert _wolf_rayet(tr, 1.0e6) is False


# ------------------------------------------------------------------ v0.31: smooth and fast meshes
def _vw_lmi(**edit):
    import copy

    raw = yaml.safe_load((ROOT / "configs" / "catalog" / "vw_lmi_contact.yaml").read_text())
    raw = copy.deepcopy(raw)
    raw["orbits"]["A"]["i"] += edit.get("di", 0.0)
    cfg = Config.model_validate(raw)
    cfg._base_dir = str(ROOT / "configs" / "catalog")
    return cfg


def test_mesh_flux_is_smooth_in_the_parameters():
    """The occultation of mesh elements is anti-aliased: the flux is a smooth function of the
    inclination (finite differences with steps 0.001 and 0.01 deg agree), also in eclipse."""
    from swingersim.photometry import light_curve

    t = np.array([0.30, 0.36, 0.42])
    f0 = light_curve(_vw_lmi(), t, "TESS").flux
    d1 = (light_curve(_vw_lmi(di=0.001), t, "TESS").flux - f0) / 0.001
    d2 = (light_curve(_vw_lmi(di=0.01), t, "TESS").flux - f0) / 0.01
    assert np.allclose(d1, d2, rtol=0.05, atol=2e-5)


def test_fast_mesh_mode_matches_the_full_mesh():
    from swingersim import mesh
    from swingersim.photometry import light_curve

    t = np.linspace(0.3, 1.3, 60)
    old = mesh.FAST
    try:
        mesh.FAST = False
        full = light_curve(_vw_lmi(), t, "TESS").flux
        mesh.FAST = True
        fast = light_curve(_vw_lmi(), t, "TESS").flux
    finally:
        mesh.FAST = old
    assert np.max(np.abs(fast - full)) < 3e-4


def test_shape_policy():
    from swingersim.shapes import apply_policy, fill_factors

    raw = yaml.safe_load((ROOT / "configs" / "catalog" / "vw_lmi_contact.yaml").read_text())
    raw["orbits"]["A"].pop("contact_fillout")
    raw["orbits"]["A"]["P"] = 1.5
    raw["stars"]["m1"]["R"] = 2.6  # near its Roche lobe
    ff = fill_factors(Config.model_validate(raw))
    assert ff["A"][0] > 0.7 and ff["B"][0] < 0.3
    d, notes = apply_policy(raw, "auto", 0.5)
    assert d["stars"]["m1"]["shape"] == "roche" and d["stars"]["m2"]["shape"] == "roche"
    assert "shape" not in (d["stars"].get("m3") or {})
    assert any("Roche surfaces" in n for n in notes)
    raw["orbits"]["A"]["e"] = 0.1  # eccentric: spheres kept
    d, notes = apply_policy(raw, "auto", 0.5)
    assert "shape" not in d["stars"]["m1"] and any("circular" in n for n in notes)
    d, _ = apply_policy(raw, "spheres")
    assert all("shape" not in (s or {}) for s in d["stars"].values())


def test_detached_roche_stars_keep_their_size_and_reflect():
    """Detached Roche lobes in a slightly eccentric orbit do not breathe with the separation,
    and they get the reflection of the spheres: Roche and sphere light curves of a pair at fill
    ~0.55 agree to a few ppt after the linear calibration of the fit."""
    from swingersim import mesh
    from swingersim.photometry import light_curve

    d = {
        "topology": "custom",
        "hierarchy": {"inner": ["m1", "m2"]},
        "masses": {"m1": 1.748, "m2": 1.361},
        "orbits": {
            "inner": {"P": 3.284187, "e": 0.01, "i": 88.7, "Omega": 0, "omega": 30, "M0": 0}
        },
        "stars": {"m1": {"R": 2.996, "Teff": 6612}, "m2": {"R": 1.507, "Teff": 6724}},
        "physics": {"gr": False, "tides": False, "rotation": False},
        "run": {"t_max_yr": 1, "n_output": 11},
    }
    t = np.linspace(0.0, 3.284187, 120)
    sph = light_curve(Config.model_validate(d), t, "TESS").flux
    for k in ("m1", "m2"):
        d["stars"][k]["shape"] = "roche"
    old = (mesh.N_THETA, mesh.N_PHI, mesh.TABLE_PHASES)
    try:
        mesh.N_THETA, mesh.N_PHI, mesh.TABLE_PHASES = 32, 64, 64
        roc = light_curve(Config.model_validate(d), t, "TESS").flux
    finally:
        mesh.N_THETA, mesh.N_PHI, mesh.TABLE_PHASES = old
    A = np.column_stack([sph, np.ones_like(sph)])
    c, *_ = np.linalg.lstsq(A, roc, rcond=None)
    res = roc - A @ c
    out = sph > 0.995  # outside the eclipses (the eclipse shapes differ by up to ~1%)
    assert np.std(res) < 2.5e-3 and np.median(np.abs(res[out])) < 1e-3
