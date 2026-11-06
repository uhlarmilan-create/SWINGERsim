"""Photometry area: Roche root finding and scene cache, inclination tables,
proximity effects with meshes, Roche lobe check, grid-edge policy, CSV metadata, MIST phases,
apsidal table consistent with the N-body run, viewer size."""

import math
import warnings

import numpy as np
import pytest

from swingersim.config import Config


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


# --------------------------------------------------------------------------- Roche geometry


@pytest.mark.parametrize(
    ("q", "r", "star", "omega_bisection"),
    [
        # values of the 40-step bisection of an older version
        (0.78, 0.25, 0, 4.818285364185616),
        (0.78, 0.2, 1, 5.03441524851471),
        (0.3, 0.5, 0, 2.424627957344468),
        (2.0, 0.3, 0, 5.441619528813584),
    ],
)
def test_potential_for_radius_secant_matches_bisection(q, r, star, omega_bisection):
    from swingersim import roche

    om = roche.potential_for_radius(q, r, star)
    assert om == pytest.approx(omega_bisection, abs=1e-9)
    assert roche.volume_radius(q, om, star) == pytest.approx(r, rel=1e-12)


def test_potential_for_radius_beyond_the_outer_critical_surface_raises():
    from swingersim import roche

    with pytest.raises(ValueError, match="outer critical surface"):
        roche.potential_for_radius(0.5, 0.7, 0)


def test_cached_ray_potentials_give_identical_surfaces():
    from swingersim import roche

    a = roche.binary_surface(0.6, 4.0, 0, 20, 40)
    roche._ray_potentials_cached.cache_clear()
    b = roche.binary_surface(0.6, 4.0, 0, 20, 40)
    assert a.volume == b.volume and np.array_equal(a.points, b.points)


def _roche_pair(R1=2.0, R2=1.5, i=85.0, e=0.0):
    return Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"inner": ["m1", "m2"]},
            "masses": {"m1": 1.7, "m2": 1.4},
            "orbits": {"inner": {"P": 3.28, "e": e, "i": i, "Omega": 0, "omega": 30, "M0": 0}},
            "stars": {
                "m1": {"R": R1, "Teff": 6600, "shape": "roche"},
                "m2": {"R": R2, "Teff": 6700, "shape": "roche"},
            },
            "physics": {"gr": False, "tides": False, "rotation": False},
            "run": {"t_max_yr": 1, "n_output": 11},
        }
    )


@pytest.fixture
def coarse_mesh():
    from swingersim import mesh

    old = (mesh.N_THETA, mesh.N_PHI, mesh.TABLE_PHASES)
    mesh.N_THETA, mesh.N_PHI, mesh.TABLE_PHASES = 24, 48, 32
    yield mesh
    mesh.N_THETA, mesh.N_PHI, mesh.TABLE_PHASES = old


def test_detached_roche_star_larger_than_its_lobe_is_an_error(coarse_mesh):
    from swingersim.photometry import light_curve

    with pytest.raises(ValueError, match="larger than its Roche lobe"):
        light_curve(_roche_pair(R1=6.0), np.array([0.0, 0.1]), "TESS")


def test_scene_cache_reuses_the_roche_geometry(coarse_mesh, monkeypatch):
    """A second light curve with the same Roche parameters (e.g. another inclination step of
    a fit) builds neither the surfaces nor the phase tables again."""
    from swingersim import roche
    from swingersim.photometry import light_curve

    mesh = coarse_mesh
    mesh._SCENE_CACHE.clear()
    t = np.linspace(0.5, 1.0, 4)  # out of eclipse at i = 70
    calls = {"surface": 0, "table": 0}
    surf, table = roche.binary_surface, mesh._phase_table

    def count_surface(*a, **k):
        calls["surface"] += a[3:4] == (mesh.N_THETA,)  # mesh surfaces (not volume radii)
        return surf(*a, **k)

    def count_table(*a, **k):
        calls["table"] += 1
        return table(*a, **k)

    monkeypatch.setattr(roche, "binary_surface", count_surface)
    monkeypatch.setattr(mesh, "_phase_table", count_table)
    f1 = light_curve(_roche_pair(i=70.0), t, "TESS", light_time=False).flux
    first = dict(calls)
    assert first["surface"] == 2 and first["table"] >= 1
    f2 = light_curve(_roche_pair(i=70.0), t, "TESS", light_time=False).flux
    assert calls == first  # nothing rebuilt
    assert np.array_equal(f1, f2)


def test_phase_tables_are_interpolated_linearly_in_the_inclination(coarse_mesh):
    """Tables on a fixed inclination grid: the flux of an unocculted Roche pair is continuous
    and piecewise linear in the inclination (no steps between table rebuilds)."""
    from swingersim.photometry import light_curve

    mesh = coarse_mesh
    t = np.array([0.6])
    i0 = 70.0 + 0.3 * mesh.TABLE_DI_DEG  # inside one grid cell
    f = [
        light_curve(_roche_pair(i=i0 + d), t, "TESS", light_time=False).flux[0]
        for d in (0.0, 0.1 * mesh.TABLE_DI_DEG, 0.2 * mesh.TABLE_DI_DEG)
    ]
    assert f[1] != f[0]  # not piecewise constant
    assert f[2] - f[1] == pytest.approx(f[1] - f[0], rel=1e-3, abs=1e-12)  # linear in a cell


def test_rebuilt_table_keeps_the_lobe_factor_of_the_instant(coarse_mesh):
    """Building a phase table (a fake frame at separation 1) must not change the lobe size
    factor of the current instant, which the frames computed once per instant reuse."""
    from swingersim import roche
    from swingersim.photometry import StarLight

    mesh = coarse_mesh
    q = 0.8
    om = (roche.potential_for_radius(q, 0.15, 0), roche.potential_for_radius(q, 0.12, 1))
    geo = mesh._roche_geometry(q, om)
    b = mesh.Body("roche", [0, 1], geo["surfs"], q=q, omegas=om, extent=geo["extent"], a_ref=1.0)
    b.k = 1.03
    ld = np.array([0.5, 0.0, 0.0, 0.0])
    lights = [StarLight(f"m{k + 1}", 0.01, 1.0, ld, 0.3, 0.5, 2.0) for k in (0, 1)]
    stars = {
        k: mesh.MeshStar(k, b, sf, mesh._rel_intensity(sf, 0.3), sf.area)
        for k, sf in enumerate(geo["surfs"])
    }
    f = mesh._free_fluxes(b, stars, lights, 2, np.array([0.0, 0.5, math.sqrt(0.75)]))
    assert np.all(f > 0)
    assert b.k == 1.03 and b.a_ref == 1.0


# --------------------------------------------------------------------------- proximity


def test_mesh_keeps_the_ellipsoidal_variation_of_non_roche_pairs():
    """A close non-eclipsing pair of spheres: with a mesh star in the system (here a slowly
    rotating, nearly spherical star in the same pair) the ellipsoidal variation of
    photometry.flux_at is kept (swingersim <= 0.35 dropped it)."""
    from swingersim.photometry import light_curve

    base = {
        "masses": {"m1": 1.5, "m2": 1.2},
        "orbits": {"A": {"P": 1.2, "i": 60.0}},
        "run": {"t_max_yr": 0.01},
    }
    t = np.linspace(0.0, 1.2, 25)
    sph = light_curve(
        _binary(**base, stars={"m1": {"R": 1.6, "Teff": 6500}, "m2": {"R": 1.2, "Teff": 6000}}),
        t,
        "V",
        light_time=False,
    ).flux
    rot = light_curve(
        _binary(
            **base,
            stars={
                "m1": {"R": 1.6, "Teff": 6500},
                "m2": {"R": 1.2, "Teff": 6000, "shape": "rotating", "P_rot": 1000.0},
            },
        ),
        t,
        "V",
        light_time=False,
    ).flux
    amp_s = np.ptp(sph)
    assert amp_s > 5e-3  # a clear ellipsoidal variation
    assert np.ptp(rot) == pytest.approx(amp_s, rel=0.05)
    assert np.max(np.abs(rot / rot.mean() - sph / sph.mean())) < 1e-3


def test_proximity_terms_reproduce_flux_at():
    from swingersim.photometry import StarLight, flux_at, proximity_terms

    lights = [
        StarLight("m1", 0.01, 1.0, np.array([0.5, 0.1, 0.0, 0.0]), 0.3, 0.5, 2.0),
        StarLight("m2", 0.008, 0.7, np.array([0.6, 0.0, 0.0, 0.0]), 0.25, 1.0, 1.5),
    ]
    pos = np.array([[0.0, 0.0, 0.0], [0.03, 0.02, 0.05]])  # no eclipse
    vel = np.array([[0.0, -1.0, 2.0], [0.0, 1.5, -2.5]])
    m = np.array([1.0, 0.8])
    base = np.array([math.pi * s.radius_au**2 * s.weight for s in lights])
    ell, refl, beam = proximity_terms(pos, vel, m, lights, [(0, 1)], base)
    _, f = flux_at(pos, vel, m, lights, [(0, 1)])
    assert np.allclose(f * base.sum(), (base * (1 + ell) + refl) * beam, rtol=1e-14)
    e0, r0, _ = proximity_terms(pos, vel, m, lights, [(0, 1)], base, set(), set())
    assert not e0.any() and not r0.any()


# --------------------------------------------------------------------------- grid policy


def test_claret_grid_edges_clip_with_a_warning_naming_the_cause():
    from swingersim import gridclip, photometry
    from swingersim.photometry import claret_coefficients

    gridclip._WARNED.clear()
    photometry._LD_WARNED.clear()
    a_in, y_in = claret_coefficients(5800.0, 4.4, 0.0, "V")  # inside: plain interpolation
    direct = photometry._claret_interp(5800.0, 4.4, 0.0, "V")
    assert y_in == direct[0] and np.array_equal(a_in, direct[1:])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert np.allclose(claret_coefficients(5800.0, 4.4, 0.0, "V")[0], a_in)
    with pytest.warns(UserWarning, match=r"log g = 5\.3 is outside"):
        a, _ = claret_coefficients(3600.0, 5.3, 0.0, "V")
    assert np.allclose(a, claret_coefficients(3600.0, 5.0, 0.0, "V")[0])
    with pytest.warns(UserWarning, match=r"\[M/H\] = 1\.6 is outside"):
        a, _ = claret_coefficients(5800.0, 4.4, 1.6, "V")
    assert np.allclose(a, claret_coefficients(5800.0, 4.4, 1.0, "V")[0])
    with pytest.warns(UserWarning, match="below the coolest model"):
        claret_coefficients(3000.0, 0.5, 0.0, "V")


def test_o_stars_beyond_claret_get_von_zeipel_gravity_darkening():
    from swingersim.photometry import limb_darkening, tlusty_coefficients, von_zeipel_y

    with warnings.catch_warnings():
        warnings.simplefilter("error")  # an explicit choice, no fallback warning
        a, y = limb_darkening(45000.0, 4.0, 0.0, "V")  # Claret ends at ~39 000 K at log g 4
    assert y == pytest.approx(von_zeipel_y(45000.0, "V"))
    assert np.allclose(a, tlusty_coefficients(45000.0, 4.0, 0.0, "V"))
    _, y2 = limb_darkening(20000.0, 4.0, 0.0, "V")  # inside Claret: its y
    from swingersim.photometry import claret_coefficients

    assert y2 == pytest.approx(claret_coefficients(20000.0, 4.0, 0.0, "V")[1])


def test_mist_metallicity_tolerance():
    from swingersim.mist import MistTrack

    a = MistTrack(1.0, 0.5, 0.0)
    with pytest.warns(UserWarning, match="installed MIST grids"):
        b = MistTrack(1.0, 0.55, 0.0)
    assert b.t_ms == a.t_ms
    with pytest.raises(ValueError, match="by more than 0.1"):
        MistTrack(1.0, 0.7, 0.0)


def test_mist_phases_begin_at_the_primary_eeps():
    """Dotter (2016): ZAMS = EEP 202, TAMS = 454, RGB tip = 605, TACHeB = 707, TP-AGB = 808
    (1-based); the phase changes at these points, consistent with t_ms = age[EEP_TAMS - 1]."""
    from swingersim.mist import _phase

    assert _phase(201.5) == 0 and _phase(202.0) == 1
    assert _phase(453.9) == 1 and _phase(454.0) == 2
    assert _phase(604.9) == 3 and _phase(605.0) == 4
    assert _phase(706.9) == 4 and _phase(707.0) == 5
    assert _phase(807.9) == 5 and _phase(808.0) == 6


# --------------------------------------------------------------------------- outputs


def test_light_curve_csv_has_a_metadata_header(tmp_path):
    from swingersim.photometry import add_noise, light_curve

    lc = add_noise(light_curve(_binary(), np.linspace(0, 1, 5), "V", light_time=False), 200, 7)
    p = tmp_path / "lc.csv"
    lc.write_csv(p)
    lines = p.read_text().splitlines()
    head = [ln for ln in lines if ln.startswith("#")]
    text = "\n".join(head)
    for s in ("band: V", "light-travel time: off", "proximity", "exposure", "normalisation"):
        assert s in text
    assert "noise: 200 ppm (Gaussian), seed 7" in text and "SWINGERsim" in head[0]
    assert "a_k =" in text and "not given" in text  # LD coefficients; no BJD of t = 0
    cols = lines[len(head)]
    assert cols.startswith("t_day,flux,flux_m1")
    d = np.loadtxt(p, delimiter=",", comments="#", skiprows=len(head) + 1)
    assert np.allclose(d[:, 1], lc.flux, rtol=1e-9)


def test_light_curve_header_gives_bjd_of_t0():
    from swingersim.photometry import LightCurve

    lc = LightCurve(np.zeros(1), np.ones(1), "TESS", meta={"bjd_t0": 2459000.5})
    assert any("BJD_TDB = 2459000.500000 + t_day" in ln for ln in lc.header_lines())


def test_phoebe_engine_rejects_proximity_off_with_roche():
    from swingersim.phoebe_lc import phoebe_light_curve

    with pytest.raises(ValueError, match="no proximity switch"):
        phoebe_light_curve(_binary(), np.zeros(2), "V", proximity=False)


def test_apsidal_table_uses_the_k2_of_the_nbody_run():
    from swingersim import apsidal
    from swingersim.tides import k2_rotating

    phys = {"gr": False, "tides": True, "rotation": True, "k2_rotation": True}
    cfg = _binary(physics=phys, orbits={"A": {"P": 2.0, "e": 0.1, "i": 89.0}})
    row = apsidal.apsidal_table(cfg)[0]
    static = apsidal.apsidal_table(
        _binary(
            physics={**phys, "k2_rotation": False}, orbits={"A": {"P": 2.0, "e": 0.1, "i": 89.0}}
        )
    )[0]
    assert row["k2_1"] < static["k2_1"]  # rotation makes the star more centrally condensed
    w = 2 * math.pi / (row["P_rot_1_d"] / 365.25)
    assert row["k2_1"] == pytest.approx(k2_rotating(cfg, 0, static["k2_1"], w, 1.0, 1.0))


# --------------------------------------------------------------------------- viewer


def test_viewer_decimation_keeps_the_envelope_and_rounds():
    from swingersim.viewer import STATIC_MAX_POINTS, _compact, _decimate

    t = np.linspace(0, 100, 200_001)
    y = np.sin(t * 7.3) * (1 + 0.01 * t)
    k = _decimate(len(t), y)
    assert len(k) <= STATIC_MAX_POINTS + 2
    assert y[k].max() == y.max() and y[k].min() == y.min()
    assert k[0] == 0 and k[-1] == len(t) - 1
    c = np.asarray(_compact([1.234567891, -98765.4321, 0.0, np.nan]), float)
    assert c[0] == pytest.approx(1.23457, rel=1e-6) and c[1] == pytest.approx(-98765.4, rel=1e-6)
    assert c[2] == 0.0 and np.isnan(c[3])


# --------------------------------------------------------------------------- tides


def test_leconte_rates_zero_lag_and_zero_inertia():
    """Same rule as secular.kernel.leconte_pair: lag <= 0 contributes nothing; rg2 = 0 is a
    pseudo-synchronous star with its spin held (no NaN); normal stars are unchanged."""
    from swingersim.tides import leconte_rates

    m, r, k2 = (1.0, 0.8), (0.005, 0.004), (0.01, 0.012)
    a, e, spin = 0.05, 0.2, (40.0, 60.0)
    lag = (1e-8, 2e-8)
    full = leconte_rates(m, r, k2, (0.08, 0.1), lag, a, e, spin)
    one = leconte_rates(m, r, k2, (0.08, 0.1), (lag[0], 0.0), a, e, spin)
    alone = leconte_rates(m, r, (k2[0], 0.0), (0.08, 0.1), (lag[0], 1.0), a, e, spin)
    assert one.da_dt == pytest.approx(alone.da_dt, rel=1e-12)
    assert one.de_dt == pytest.approx(alone.de_dt, rel=1e-12)
    assert one.dspin_dt[1] == 0.0 and one.t_sync[1] == np.inf
    assert one.dspin_dt[0] == pytest.approx(full.dspin_dt[0], rel=1e-12)
    z = leconte_rates(m, r, k2, (0.0, 0.1), lag, a, e, spin)
    assert np.isfinite([z.da_dt, z.de_dt, *z.dspin_dt]).all() and z.dspin_dt[0] == 0.0
    # rg2 -> 0 limit: a tiny moment of inertia of a pseudo-synchronous star gives the same rates
    n = math.sqrt(39.47841760435743 * 1.8 / a**3)
    ps = full.spin_pseudo_sync / n
    near = leconte_rates(m, r, k2, (1e-9, 0.1), lag, a, e, (ps * n, spin[1]))
    assert z.da_dt == pytest.approx(near.da_dt, rel=1e-9)
    assert z.de_dt == pytest.approx(near.de_dt, rel=1e-9)
