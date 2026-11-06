"""Own MESA grids: inlists, driver (with a stand-in for MESA), reading and reducing profiles."""

import os
import stat
import sys
from pathlib import Path

import numpy as np
import pytest

from swingersim.inertial_waves import FULLY_CONVECTIVE_C, polytrope_profile
from swingersim.mesa import (
    G_CGS,
    MSUN_G,
    RSUN_CM,
    MesaGridSpec,
    inlist_text,
    read_mesa_table,
    reduce_grid,
    reduce_profile,
    run_grid,
)
from swingersim.stellar import polytrope_constants


def write_profile(path, n_index, conv_from=0.0, mass=0.3, radius=0.3, age=1e7, n=3000):
    """A polytrope written in the MESA profile format (surface first)."""
    r, rho, _ = polytrope_profile(n_index, n)
    m = (
        np.concatenate(
            [[0.0], np.cumsum(0.5 * (rho[1:] * r[1:] ** 2 + rho[:-1] * r[:-1] ** 2) * np.diff(r))]
        )
        * 4
        * np.pi
    )
    dpdr = -rho * m / r**2
    p = np.concatenate([[0.0], np.cumsum(-0.5 * (dpdr[1:] + dpdr[:-1]) * np.diff(r))])
    p = p[-1] - p
    keep = p > 0  # the outermost point has p = 0 (a MESA model ends at the photosphere)
    r, rho, m, p = r[keep], rho[keep], m[keep], p[keep]
    rr, mm = radius * RSUN_CM, mass * MSUN_G
    conv = (r >= conv_from).astype(int)
    cols = {
        "zone": np.arange(len(r), 0, -1),
        "radius": r * radius,
        "mass": m * mass,
        "logRho": np.log10(rho * mm / rr**3 + 1e-30),
        "logP": np.log10(p * G_CGS * mm**2 / rr**4),
        "mixing_type": conv,
        "conv_vel": 1e4 * conv,
    }
    names = list(cols)
    with open(path, "w") as f:
        f.write("1 2 3 4 5\nstar_age star_mass photosphere_L Teff center_h1\n")
        f.write(f"{age} {mass} 0.01 3300 0.7\n\n")
        f.write(" ".join(str(i + 1) for i in range(len(names))) + "\n" + " ".join(names) + "\n")
        f.writelines(
            " ".join(f"{cols[c][k]:.12e}" for c in names) + "\n" for k in range(len(r) - 1, -1, -1)
        )


def test_reduce_polytropes(tmp_path):
    for n_index in (1.5, 3.0):
        write_profile(tmp_path / "p.data", n_index)
        d = reduce_profile(*read_mesa_table(tmp_path / "p.data"))
        k2, rg2 = polytrope_constants(n_index)
        assert d["apsidal_k2"] == pytest.approx(k2, rel=1e-2)
        assert d["rg2"] == pytest.approx(rg2, rel=3e-3)
        assert d["fully_convective"]
    write_profile(tmp_path / "p.data", 1.5)
    d = reduce_profile(*read_mesa_table(tmp_path / "p.data"))
    assert d["iw_c"] == pytest.approx(FULLY_CONVECTIVE_C, rel=1e-2)
    write_profile(tmp_path / "p.data", 3.0, conv_from=0.7)
    d = reduce_profile(*read_mesa_table(tmp_path / "p.data"))
    assert not d["fully_convective"] and 0 < d["m_env"] < 0.3 and d["iw_c"] > 0
    # a stratified envelope dissipates like the two-layer model within a factor of a few
    # (Barker 2020: about 2 for solar-type stars); a sign error in the equilibrium tide gave 50x
    from swingersim.inertial_waves import two_layer

    alpha, beta = 1 - d["r_env"] / d["radius_rsun"], 1 - d["m_env"] / d["mass"]
    assert 0.5 < d["iw_c"] / two_layer(alpha, beta) < 3.0


def test_inlist():
    spec = MesaGridSpec(name="t", masses=[1.0, 1.4], z=[0.02])
    txt = inlist_text(spec, 1.4, 0.02, Path("/g"))
    for key in (
        "&star_job",
        "&controls",
        "create_pre_main_sequence_model = .true.",
        "initial_mass = 1.4",
        "xa_central_lower_limit(1)",
        "overshoot_f(1)",
    ):
        assert key in txt
    assert "overshoot" not in inlist_text(spec, 1.0, 0.02, Path("/g"))


FAKE_STAR = """#!{python}
import re, sys
sys.path.insert(0, {tests!r})
from pathlib import Path
from test_mesa import write_profile
mass = float(re.search(r"initial_mass = ([0-9.d+-]+)", open("inlist").read()).group(1).replace("d", "e"))
logs = Path("LOGS"); logs.mkdir()
(logs / "history.data").write_text("1\\nversion\\n1\\n\\n1 2\\nmodel_number star_age\\n1 1e5\\n")
for k, age in enumerate((1e5, 1e6, 1e7, 1e8, 1e9)):
    write_profile(logs / f"profile{{k + 1}}.data", 1.5, mass=mass, radius=mass, age=age, n=2000)
"""


def test_driver_with_stand_in_mesa(tmp_path, monkeypatch):
    mesa, sdk, root = tmp_path / "mesa", tmp_path / "sdk", tmp_path / "root"
    work = mesa / "star" / "work"
    work.mkdir(parents=True)
    (sdk / "bin").mkdir(parents=True)
    (sdk / "bin" / "mesasdk_init.sh").write_text("# stand-in\n")
    star = FAKE_STAR.format(python=sys.executable, tests=str(Path(__file__).parent))
    (work / "mk").write_text("#!/bin/bash\ncp star.py star\nchmod +x star\n")
    (work / "star.py").write_text(star)
    os.chmod(work / "mk", os.stat(work / "mk").st_mode | stat.S_IEXEC)
    monkeypatch.setenv("MESA_DIR", str(mesa))
    monkeypatch.setenv("MESASDK_ROOT", str(sdk))
    spec = MesaGridSpec(name="g", masses=[0.3, 0.4], z=[0.02], threads_per_run=1)
    monkeypatch.chdir(tmp_path)
    res = run_grid(spec, Path("root"), processes=2)  # relative root, as from the command line
    assert res == {"M00.3000_Z0.02000": 0, "M00.4000_Z0.02000": 0}
    assert run_grid(spec, root, processes=2) == res  # resumable: nothing to do
    out = reduce_grid(spec, root, processes=2)
    with np.load(out) as d:
        assert len(d["age_myr"]) == 10 and set(d["initial_mass"]) == {0.3, 0.4}
        assert np.allclose(d["iw_c"], FULLY_CONVECTIVE_C, rtol=0.05)
