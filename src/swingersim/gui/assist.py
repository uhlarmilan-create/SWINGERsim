"""Help and small computations that make the pages of swingersim gui self-explaining:

- help texts (help.yaml, served as JSON);
- quick check of a configuration while it is edited: stability ratios of every nested triple,
  eclipse geometry and Roche-lobe filling of every pair, plain-language warnings;
- the Applegate energy check for a star of a configuration;
- the three-body break-up theory of Stone & Leigh (2019) against a scattering experiment;
- own MESA grids: list, status, saving a grid definition, the reduced grid for plots;
- a "Methods" paragraph (with references) for a configuration, for papers.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import yaml

from .digest import _f

HELP = Path(__file__).with_name("help.yaml")
_HELP_CACHE: dict = {}


def help_texts() -> dict:
    key = HELP.stat().st_mtime_ns
    if _HELP_CACHE.get("key") != key:
        _HELP_CACHE.update(key=key, data=yaml.safe_load(HELP.read_text(encoding="utf-8")))
    return _HELP_CACHE["data"]


def _ratio(x) -> str:
    return "unknown" if x is None else f"{x:.2f}"


# --------------------------------------------------------------------------- quick check
def quickcheck(data: dict, base_dir: Path | None = None) -> dict:
    """What a configuration looks like before anything is run (a few ms; MIST once ~1 s)."""
    from ..config import Config
    from ..systemview import system_model

    if not isinstance(data, dict) or "scatter" in data:
        return {"ok": False, "items": []}
    try:
        cfg = Config.model_validate(data)
    except Exception as exc:  # noqa: BLE001 - the page shows the validation errors elsewhere
        return {"ok": False, "items": [], "error": str(exc).splitlines()[0]}
    if base_dir is not None:  # relative paths (a stellar table) are resolved against it
        cfg._base_dir = str(base_dir)
    try:
        m = system_model(cfg, n_points=8)
    except Exception as exc:  # noqa: BLE001
        return {"ok": False, "items": [], "error": f"the system could not be set up: {exc}"}
    items = []
    for t in m["triples"]:
        ma, v = _f(t["ma01"]), _f(t["v22"])
        known = [x for x in (ma, v) if x is not None]
        worst = min(known) if known else None
        tone = "" if worst is None else "good" if worst >= 1.5 else "warn" if worst >= 1 else "bad"
        items.append(
            {
                "k": f"stability of {t['label']}",
                "v": f"MA01 {_ratio(ma)} · V22 {_ratio(v)}",
                "tone": tone,
                "note": (
                    "the criteria could not be evaluated for this triple"
                    if worst is None
                    else "stable with a margin"
                    if tone == "good"
                    else "close to the stability limit – check with Stability check and a long Run"
                    if tone == "warn"
                    else "beyond the stability limit: the triple can break up"
                ),
                "help": "stability",
            }
        )
    for p in m["pairs"]:
        lim, i = p["i_limit_deg"], p["i_deg"]
        if lim is not None and i is not None:
            items.append(
                {
                    "k": f"pair {p['orbit']} ({' + '.join(p['stars'])})",
                    "v": f"i {i:.1f}° vs eclipse limit {lim:.1f}°",
                    "tone": "good" if p["eclipsing"] else "",
                    "note": "eclipses at t = 0"
                    if p["eclipsing"]
                    else "no eclipses at t = 0 (the inclination is below the limit)",
                    "help": "eclipse limit",
                }
            )
        fills = [x for x in p["fill"] if x is not None]
        if fills and max(fills) > 0.7:
            mx = max(fills)
            items.append(
                {
                    "k": f"Roche fill in {p['orbit']}",
                    "v": " / ".join(f"{x:.2f}" for x in fills),
                    "tone": "bad" if mx >= 1 else "warn",
                    "note": "a star overfills its Roche lobe – mass transfer or contact"
                    if mx >= 1
                    else "strongly deformed stars: use stars.*.shape roche for light curves "
                    "(a fit of light curves chooses Roche surfaces itself with fit.lc[].shapes: "
                    "auto)",
                    "help": "Roche fill",
                }
            )
    for o in m["orbits"]:
        e = o["e"]
        if e is not None and e > 0.9:
            items.append(
                {
                    "k": f"orbit {o['name']}",
                    "v": f"e = {e:.3f}",
                    "tone": "warn",
                    "note": "very eccentric: close periastron passages, short time steps",
                }
            )
    for mu in m["mutual"]:
        d = mu["deg"]
        if d is not None and 39.2 < d < 140.8:
            items.append(
                {
                    "k": "mutual inclination " + " / ".join(mu["orbits"]),
                    "v": f"{d:.1f}°",
                    "tone": "warn",
                    "note": "between 39° and 141°: Kozai–Lidov–Zeipel cycles can pump the inner "
                    "eccentricity",
                    "help": "mutual inclination",
                }
            )
    return {"ok": True, "items": items, "hierarchy": m["hierarchy"], "stable": m["stable"]}


# --------------------------------------------------------------------------- Applegate
def applegate(
    data: dict | None = None,
    star: str | None = None,
    p_mod_yr: float | None = None,
    oc_s: float | None = None,
    values: dict | None = None,
) -> dict:
    """Energy check of the Applegate mechanism, for a star of a configuration (data, star) or
    for given values {mass, radius, teff, a_rsun, pbin_d}."""
    from .. import magnetic as mg
    from .. import units

    v = dict(values or {})
    if data is not None and star:
        from ..backends.nbody import model_props
        from ..config import Config
        from ..hierarchy import Node, build_tree, initialize, mass_of, node_elements

        cfg = Config.model_validate(data)
        root, m_lab = build_tree(cfg)
        masses = np.asarray(m_lab, float)
        pos, vel = initialize(root, masses)
        names = [f"m{k + 1}" for k in range(len(masses))]
        if star not in names:
            raise ValueError(f"no star {star}")
        k = names.index(star)
        props = model_props(cfg, root, masses, pos, vel)
        pair = next(
            n
            for n in root.nodes()
            if not isinstance(n.left, Node)
            and not isinstance(n.right, Node)
            and k in (n.left, n.right)
        )
        el = node_elements(pair, pos, vel, masses)
        M = mass_of(pair.left, masses) + mass_of(pair.right, masses)
        v.setdefault("mass", float(masses[k]))
        v.setdefault("radius", float(props[k].radius_rsun))
        v.setdefault("teff", float(props[k].teff_k))
        v.setdefault("a_rsun", float(el.a / units.RSUN_AU))
        v.setdefault(
            "pbin_d", float(2 * math.pi * math.sqrt(el.a**3 / (units.G * M)) * units.DAYS_PER_YEAR)
        )
        v["pair"] = pair.name
    need = ("mass", "radius", "teff", "a_rsun", "pbin_d")
    if any(v.get(x) is None for x in need) or not p_mod_yr:
        raise ValueError("needs mass, radius, teff, a_rsun, pbin_d and p_mod_yr")
    kmax = mg.applegate_max_oc(
        v["mass"], v["radius"], v["teff"], v["a_rsun"], v["pbin_d"], p_mod_yr
    )
    out = {"values": v, "p_mod_yr": p_mod_yr, "max_oc_s": float(kmax)}
    if v["teff"] > 6800.0:
        out["caveat"] = (
            f"the Applegate mechanism needs a magnetically active star with a deep convective "
            f"envelope; a star of {v['teff']:.0f} K has hardly any, so the energy budget is an "
            "upper bound that does not make the mechanism plausible"
        )
    if oc_s:
        dpp = mg.dp_over_p_from_oc(oc_s, p_mod_yr)
        e = mg.applegate_energy_ratio(
            v["mass"], v["radius"], v["teff"], v["a_rsun"], v["pbin_d"], dpp, p_mod_yr
        )
        out.update(oc_s=oc_s, dp_over_p=float(dpp), energy_ratio=float(e), possible=bool(e < 1))
    return out


# --------------------------------------------------------------------------- escape theory
def escape_compare(root: Path, out_dir: Path) -> dict:
    """Stone & Leigh (2019) against a binary + single scattering experiment in out_dir."""
    from ..escape_theory import compare_with_scatter

    meta = yaml.safe_load((out_dir / "scatter.yaml").read_text())
    sc = (meta.get("config") or meta).get("scatter") or {}
    masses = list((sc.get("target") or {}).get("masses") or []) + list(
        (sc.get("projectile") or {}).get("masses") or []
    )
    if len(masses) != 3:
        return {
            "available": False,
            "why": f"the theory describes a binary meeting a single star (3 stars); this experiment has {len(masses)}",
        }
    if len((sc.get("target") or {}).get("masses") or []) == 1:  # target single, projectile binary
        masses = masses[1:] + masses[:1]
    rows = [json.loads(x) for x in (out_dir / "rows.jsonl").read_text().splitlines() if x.strip()]
    try:
        c = compare_with_scatter(rows, masses)
    except ValueError as exc:
        return {"available": False, "why": str(exc)}
    return {
        "available": True,
        "masses": masses,
        "n": c["n"],
        "observed": [c["observed"][s] for s in range(3)],
        "predicted": [c["predicted"][s] for s in range(3)],
        "observed_mean_e": c["observed_mean_e"],
        "predicted_mean_e": c["predicted_mean_e"],
    }


# --------------------------------------------------------------------------- MESA
def mesa_status(root: Path) -> dict:
    """Grid definitions in configs/mesa/ and how far each one got."""
    from ..mesa import MesaGridSpec, grid_dir, run_name

    grids = []
    for f in sorted((root / "configs" / "mesa").glob("*.yaml")):
        try:
            spec = MesaGridSpec.load(f)
        except Exception as exc:  # noqa: BLE001
            grids.append({"file": str(f.relative_to(root)), "error": str(exc)})
            continue
        gdir = grid_dir(spec, root)
        runs = spec.runs()
        done = sum(
            1
            for m, z in runs
            if any((gdir / "runs" / run_name(m, z) / "LOGS").glob("profile*.data"))
        )
        npz = gdir / f"{spec.name}.npz"
        grids.append(
            {
                "file": str(f.relative_to(root)),
                "name": spec.name,
                "masses": spec.masses,
                "z": spec.z,
                "stop": spec.stop,
                "tracks": len(runs),
                "tracks_with_profiles": done,
                "reduced": npz.is_file(),
                "spec": yaml.safe_load(f.read_text()) or {},
                "npz": str(npz.relative_to(root.resolve())) if npz.is_file() else None,
            }
        )
    mesa_dir = os.environ.get("MESA_DIR")
    return {
        "grids": grids,
        "mesa_dir": mesa_dir,
        "installed": bool(mesa_dir and Path(mesa_dir).is_dir())
        or any(g.get("error") is None and _grid_mesa_dir(root, g) for g in grids),
    }


def _grid_mesa_dir(root: Path, g: dict) -> bool:
    try:
        d = yaml.safe_load((root / g["file"]).read_text()).get("mesa_dir")
    except Exception:  # noqa: BLE001
        return False
    return bool(d and Path(d).is_dir())


def mesa_save(root: Path, spec: dict, overwrite: bool = False) -> dict:
    """Write a grid definition to configs/mesa/<name>.yaml (checked by MesaGridSpec)."""
    from ..mesa import MesaGridSpec

    name = str(spec.get("name") or "").strip()
    if not name or not all(c.isalnum() or c in "_-." for c in name):
        raise ValueError("a grid needs a name of letters, digits, _ - .")
    path = root / "configs" / "mesa" / f"{name}.yaml"
    if path.exists() and not overwrite:
        raise ValueError(f"{path.relative_to(root)} exists")
    path.parent.mkdir(parents=True, exist_ok=True)
    text = yaml.safe_dump({k: v for k, v in spec.items() if v not in (None, "")}, sort_keys=False)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(text)
    MesaGridSpec.load(tmp)  # raises on a bad definition
    tmp.replace(path)
    return {"saved": str(path.relative_to(root))}


def mesa_grid(root: Path, name: str) -> dict:
    """The reduced grid for plots: per (initial mass, Z) track the age and the quantities."""
    from ..mesa import MesaGridSpec, grid_dir
    from .server import safe_name

    spec = MesaGridSpec.load(root / "configs" / "mesa" / f"{safe_name(name, 'grid name')}.yaml")
    npz = grid_dir(spec, root) / f"{spec.name}.npz"
    if not npz.is_file():
        raise ValueError(f"{name}: not reduced yet (run mesa-reduce)")
    with np.load(npz) as d:
        g = {k: d[k] for k in d.files}
    cols = [
        c
        for c in (
            "radius_rsun",
            "teff_k",
            "luminosity_lsun",
            "apsidal_k2",
            "rg2",
            "turnover_yr",
            "iw_c",
            "center_h1",
        )
        if c in g
    ]
    tracks = []
    for z in np.unique(g["z"]):
        for m in np.unique(g["initial_mass"][g["z"] == z]):
            s = (g["z"] == z) & (g["initial_mass"] == m)
            o = np.argsort(g["age_myr"][s])
            tracks.append(
                {
                    "mass": float(m),
                    "z": float(z),
                    "age_myr": [_f(x) for x in g["age_myr"][s][o]],
                    **{c: [_f(x) for x in g[c][s][o]] for c in cols},
                }
            )
    return {"name": spec.name, "columns": cols, "tracks": tracks, "records": len(g["age_myr"])}


# --------------------------------------------------------------------------- methods text
REFS = {
    "rebound": "Rein H., Liu S.-F. (2012), A&A 537, A128",
    "ias15": "Rein H., Spiegel D. S. (2015), MNRAS 446, 1424",
    "reboundx": "Tamayo D., Rein H., Shi P., Hernandez D. M. (2020), MNRAS 491, 2885",
    "tides_spin": "Lu T., Rein H., Tamayo D., et al. (2023), ApJ 948, 41",
    "mist": "Choi J., Dotter A., Conroy C., et al. (2016), ApJ 823, 102",
    "sse": "Hurley J. R., Pols O. R., Tout C. A. (2000), MNRAS 315, 543",
    "tout": "Tout C. A., Pols O. R., Eggleton P. P., Han Z. (1996), MNRAS 281, 257",
    "claret2023": "Claret A. (2023), A&A 674, A67",
    "claret2004": "Claret A. (2004), A&A 424, 919",
    "claret_ld": "Claret A., Bloemen S. (2011), A&A 529, A75",
    "ma01": "Mardling R. A., Aarseth S. J. (2001), MNRAS 321, 398",
    "v22": "Vynatheya P., Hamers A. S., Mardling R. A., Bellinger E. P. (2022), MNRAS 516, 4146",
    "emcee": "Foreman-Mackey D., Hogg D. W., Lang D., Goodman J. (2013), PASP 125, 306",
    "hurley2002": "Hurley J. R., Tout C. A., Pols O. R. (2002), MNRAS 329, 897",
    "applegate": "Applegate J. H. (1992), ApJ 385, 621",
    "tian2009": "Tian Y.-P., Xiang F.-Y., Tao X. (2009), Ap&SS 319, 119",
    "matt2015": "Matt S. P., Brun A. S., Baraffe I., et al. (2015), ApJ 799, L23",
    "barker2020": "Barker A. J. (2020), MNRAS 498, 2270",
}


def methods_text(cfg_path: Path) -> dict:
    """A 'Methods' paragraph for a configuration, with the references it needs."""
    from ..config import load_config

    cfg = load_config(cfg_path)
    used: list[str] = []

    def ref(k):
        if k not in used:
            used.append(k)
        return f"[{used.index(k) + 1}]"

    ph = cfg.physics
    parts = []
    nstars = len(cfg.star_names)
    if (cfg.run.backend or "nbody") == "secular":
        parts.append(
            f"The {nstars}-star system was integrated with the orbit-averaged (secular) equations of swingersim."
        )
    else:
        parts.append(
            f"The {nstars}-star system was integrated as a full N-body problem with REBOUND {ref('rebound')} "
            f"and its adaptive integrator IAS15 {ref('ias15')}"
            + (
                f", with additional forces from REBOUNDx {ref('reboundx')}"
                if (ph.gr or ph.tides)
                else ""
            )
            + "."
        )
    eff = []
    if ph.gr:
        eff.append("first post-Newtonian corrections for all bodies (gr_full)")
    if ph.tides:
        eff.append(
            f"tidal and rotational deformation of the stars with spin evolution (tides_spin {ref('tides_spin')})"
        )
    if ph.evolution:
        eff.append("stellar evolution of radii and masses along the tracks")
    if ph.mass_loss is not None:
        eff.append("stellar mass loss")
    if ph.mass_transfer is not None:
        eff.append(f"mass transfer by Roche-lobe overflow {ref('hurley2002')}")
    if ph.activity:
        eff.append(
            f"magnetic activity cycles in the eclipse times {ref('applegate')}, {ref('tian2009')}"
        )
    if ph.magnetic_braking is not None:
        eff.append(f"magnetic braking {ref('matt2015')}")
    if ph.inertial_waves:
        eff.append(f"dissipation of inertial waves {ref('barker2020')}")
    if eff:
        parts.append("The model includes " + "; ".join(eff) + ".")
    sm = cfg.stellar_model
    if sm.kind == "mist":
        parts.append(
            f"Stellar radii and temperatures not given explicitly were taken from the MIST tracks {ref('mist')} "
            f"at an age of {sm.age_myr:g} Myr and [Fe/H] = {sm.feh:g}; apsidal-motion constants and moments of "
            f"inertia from the MESA models of Claret {ref('claret2023')} (Claret {ref('claret2004')} above 8 M☉)."
        )
    elif sm.kind == "sse":
        parts.append(
            f"Stellar properties follow the SSE tracks {ref('sse')} at {sm.age_myr:g} Myr (Z = {sm.z:g})."
        )
    elif sm.kind == "zams":
        parts.append(
            f"Stars not specified explicitly are zero-age main-sequence stars {ref('tout')}."
        )
    if cfg.fit is not None:
        fit = cfg.fit
        s2 = fit.stage2
        if s2 is None:
            parts.append(
                "Times of minima were fitted with analytic O − C models of every eclipsing pair "
                "(light-time effect and apsidal motion; least squares)."
            )
        else:
            data = ["times of minima"]
            if fit.rv is not None:
                data.append("radial velocities")
            if fit.lc:
                data.append(
                    f"light curves (limb and gravity darkening of Claret & Bloemen {ref('claret_ld')})"
                )
            if fit.inclination:
                data.append("orbital inclinations of the eclipsing pairs at several epochs")
            what = data[0] if len(data) == 1 else ", ".join(data[:-1]) + " and " + data[-1]
            txt = f"The {what} were fitted {'jointly ' if len(data) > 1 else ''}with the N-body model (least squares"
            if s2.sampler == "emcee":
                txt += f", errors from emcee {ref('emcee')}"
            parts.append(txt[0].upper() + txt[1:] + ").")
    parts.append(
        f"Dynamical stability was assessed with the criteria of Mardling & Aarseth {ref('ma01')} and "
        f"Vynatheya et al. {ref('v22')}."
    )
    return {
        "text": " ".join(parts),
        "references": [f"[{k + 1}] {REFS[r]}" for k, r in enumerate(used)],
    }
