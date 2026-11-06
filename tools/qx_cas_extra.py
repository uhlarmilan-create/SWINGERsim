"""QX Cas: supplementary computations for the Part III report (J. Kolar's review, 2026-10-29).

1. Eclipse depth against the inclination of the 6.0-d pair (B ~ photographic, V), for the stars
   of the fit (no third light) and the inclination at which an eclipse becomes detectable
   (0.10 mag: photographic minima of the 1930s; 0.02 mag: photoelectric).
2. Depth of the eclipses over 1900-2100 along the fitted inclination i(t), with third light
   l3 = 0, 10, 20 and 37 % of the total: when were the eclipses observable, not only possible?
3. The third light needed to bring the model depths at the published inclinations (1960:
   85.4 deg, 1980: 82 deg) to the observed ones (Sandage & Tammann; Moffett & Barnes 1983).
4. The ensemble (all third bodies that tilt the pair at the observed rate), its members weighted
   by the observations and integrated over 1900-2100: the distribution of the years when the
   eclipses began, ended and come back (the dates of one fit are one case of many).

    python tools/qx_cas_extra.py --fit <stage2_fitted.yaml> --ensemble <ensemble folder>
        [--top 150] [-o output/special/qx_cas_extra]
"""

import argparse
import copy
import json
import os
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tools"))

from _common import (
    bjd,
    eclipse_switches,
    ensemble_observed,
    ensemble_rows,
    ensemble_weights,
    fit_epoch,
)

from swingersim.config import Config

YEARS = np.arange(1900.0, 2100.5, 0.5)
L3 = (0.0, 0.10, 0.20, 0.37)
DETECT = {"photographic": 0.10, "photoelectric": 0.02}


def _cfg(d: dict) -> Config:
    d = copy.deepcopy(d)
    d.pop("fit", None)
    d.pop("ensemble", None)
    d.pop("fit_epoch_bjd_tdb", None)
    return Config.model_validate(d)


def binary_only(fit: dict, inc: float, omega: float) -> Config:
    """The 6.0-d pair alone (its stars as fitted) at the inclination inc and periastron omega."""
    o = dict(fit["orbits"]["inner"], i=float(inc), omega=float(omega), Omega=0.0, M0=0.0)
    return Config.model_validate(
        {
            "topology": "custom",
            "hierarchy": {"inner": ["m1", "m2"]},
            "masses": {k: fit["masses"][k] for k in ("m1", "m2")},
            "orbits": {"inner": o},
            "stars": {k: fit["stars"][k] for k in ("m1", "m2")},
            "physics": {"gr": False},
            "run": {"t_max_yr": 0.05, "n_output": 2},
        }
    )


def depths(fit: dict, inc: float, omega: float, band: str) -> tuple[float, float]:
    """Fractional flux drops (primary, secondary) of the eclipses of the pair over one orbit (spherical stars, no
    ellipsoidal or reflection light: the eclipse alone)."""
    from swingersim.photometry import light_curve

    P = float(fit["orbits"]["inner"]["P"])
    t = np.linspace(0.0, P, 3001)
    lc = light_curve(binary_only(fit, inc, omega), t, band=band, light_time=False, proximity=False)
    f = np.asarray(lc.flux, float)
    f = f / f.max()
    # the two minima: the deepest point and the deepest point half an orbit (in phase) away
    k1 = int(np.argmin(f))
    ph = (t - t[k1]) / P % 1.0
    far = (ph > 0.15) & (ph < 0.85)
    d1, d2 = 1.0 - f[k1], (1.0 - f[far].min()) if far.any() else 0.0
    return float(max(d1, d2)), float(min(d1, d2))


def _member(job) -> dict | None:
    """Inclination history of one ensemble member over 1900-2100 and its chi^2 per table."""
    import warnings

    from swingersim.incfit import model_inclination, residuals

    warnings.simplefilter("ignore")
    mem, t_all, t0, tabs, meta = job
    n_y = len(YEARS)
    try:
        mall = np.asarray(model_inclination(_cfg(mem), "inner", t_all, t0, with_limit=True), float)
    except Exception as exc:  # noqa: BLE001
        print(f"member {meta['member']}: {exc}", flush=True)
        return None
    ii, ll = mall[0, :n_y], mall[1, :n_y]
    chi2, j = {}, n_y
    for name, q in tabs.items():
        chi2[name] = float(np.sum(residuals(q, mall[:, j : j + len(q)]) ** 2))
        j += len(q)
    sw = eclipse_switches(YEARS, ii, ll)
    ecl = (ii > ll).astype(float)
    on = [s for s in sw if np.interp(s + 0.5, YEARS, ecl) > 0.5]
    off = [s for s in sw if s not in on]
    return {
        **meta,
        "eclipsing_1933": bool(np.interp(1933.5, YEARS, ecl) > 0.5),
        "last_end_before_2025": max([s for s in off if s < 2025], default=np.nan),
        "first_start_after_2025": min([s for s in on if s > 2025], default=np.nan),
        "first_start": min(on, default=np.nan),
        "switches": sw,
        **{f"chi2_{n}": c for n, c in chi2.items()},
    }


def mag(drop, l3=0.0):
    return -2.5 * np.log10(1.0 - (1.0 - l3) * np.asarray(drop, float))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fit", required=True, help="stage2_fitted.yaml of the QX Cas fit")
    ap.add_argument("--ensemble", required=True, help="folder of the QX Cas ensemble")
    ap.add_argument("--top", type=int, default=150, help="members with the largest weights")
    ap.add_argument(
        "--inclinations",
        nargs="*",
        default=[
            str(ROOT / f"data/observed/qx_cas/{n}.csv")
            for n in ("inclination_limits", "inclination_limits_v3")
        ],
        help="inclination tables (fit.inclination format) that also weight the ensemble members",
    )
    ap.add_argument("-j", "--processes", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--skip-ensemble", action="store_true", help="only the depths (parts 1-3)")
    ap.add_argument("-o", "--output", default=str(ROOT / "output/special/qx_cas_extra"))
    a = ap.parse_args(argv)
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    fit = yaml.safe_load(Path(a.fit).read_text())
    t0 = fit_epoch(Path(a.fit))
    res = {}

    # 1. depth against inclination (omega of the fit and its range over the apsidal cycle)
    incs = np.arange(66.0, 90.01, 0.5)
    omegas = (fit["orbits"]["inner"].get("omega", 0.0), 0.0, 90.0, 180.0, 270.0)
    table = {}
    for band in ("B", "V"):
        for w in omegas:
            table[(band, w)] = np.array([depths(fit, i, w, band) for i in incs])
    rows = []
    for band in ("B", "V"):
        for w in omegas:
            for i, (d1, d2) in zip(incs, table[(band, w)], strict=True):
                rows.append(
                    {"band": band, "omega": w, "i": i, "drop_primary": d1, "drop_secondary": d2}
                )
    import pandas as pd

    pd.DataFrame(rows).to_csv(out / "depth_vs_inclination.csv", index=False)
    w_fit = omegas[0]
    det = {}
    for name, thr in DETECT.items():
        for l3 in (0.0, 0.2):
            best = []
            for w in omegas:
                dmag = mag(table[("B", w)][:, 0], l3)
                ok = incs[dmag >= thr]
                best.append(float(ok.min()) if len(ok) else np.nan)
            det[f"{name}_l3_{l3:.1f}"] = {
                "omega_fit": best[0],
                "range": [np.nanmin(best), np.nanmax(best)],
            }
    res["detectable_inclination_deg"] = det

    # 3. third light needed at the published inclinations (V band, omega of the fit)
    need = {}
    for yr, inc_pub, obs in ((1960, 85.4, (0.32, 0.28)), (1980, 82.0, (0.12, None))):
        d1, d2 = depths(fit, inc_pub, w_fit, "V")
        r = {"model_mag": [float(mag(d1)), float(mag(d2))], "observed_mag": list(obs)}
        r["l3_needed"] = [
            None if o is None else float(1.0 - (1.0 - 10 ** (-0.4 * o)) / d)
            for o, d in zip(obs, (d1, d2), strict=True)
        ]
        need[str(yr)] = r
    res["third_light_needed"] = need

    # 2. depth over 1900-2100 along the fitted i(t)
    from swingersim.incfit import model_inclination

    m = model_inclination(_cfg(fit), "inner", bjd(YEARS), t0, with_limit=True)
    i_t, lim_t = np.asarray(m[0], float), np.asarray(m[1], float)
    d_tab = table[("B", w_fit)][:, 0]
    drop_t = np.interp(i_t, incs, d_tab, left=0.0, right=d_tab[-1])
    drop_t[i_t < lim_t] = 0.0
    hist = pd.DataFrame({"year": YEARS, "i": i_t, "limit": lim_t,
                         **{f"depth_B_l3_{l3:.2f}": mag(drop_t, l3) for l3 in L3}})  # fmt: skip
    hist.to_csv(out / "depth_history_fit.csv", index=False)
    vis = {}
    for l3 in L3:
        dm = mag(drop_t, l3)
        for name, thr in DETECT.items():
            vis[f"{name}_l3_{l3:.2f}"] = eclipse_switches(YEARS, dm, np.full_like(dm, thr))
    res["fit"] = {"geometric_switches": eclipse_switches(YEARS, i_t, lim_t), "detectable_switches": vis,
                  "depth_B_1933": {f"{l3:.2f}": float(np.interp(1933.5, YEARS, mag(drop_t, l3))) for l3 in L3}}  # fmt: skip

    if a.skip_ensemble:
        (out / "summary.json").write_text(json.dumps(res, indent=1, default=float))
        print(json.dumps(res, indent=1, default=float))
        return 0
    # 4. ensemble: dates over the weighted members
    ens = Path(a.ensemble)
    summ = ens / "summary.csv"
    d = pd.read_csv(summ)
    f_ = d[d.direction == "forward"].set_index("member")
    b_ = d[d.direction == "backward"].set_index("member")
    mm = f_.join(b_, rsuffix="_b", how="inner")
    cfg_path = (
        ens / "ensemble.yaml"
        if (ens / "ensemble.yaml").exists()
        else ROOT / "configs/squadra_test3/qx_cas_ensemble.yaml"
    )
    w = ensemble_weights(ensemble_rows(summ), ensemble_observed(cfg_path), mm.index.values)
    base = yaml.safe_load(cfg_path.read_text())
    base = base.get("config", base)  # ensemble.yaml of a run holds the configuration
    t0_ens = 2444500.0
    order = np.argsort(-w)[: a.top]
    vary = [c for c in mm.columns if c.startswith(("masses.", "orbits.")) and not c.endswith("_b")]
    from types import SimpleNamespace

    from swingersim.incfit import load as load_inc

    tabs = {Path(f).stem: load_inc(SimpleNamespace(file=f, pair="inner", label="", min_sigma_deg=0.0), None)
            for f in a.inclinations}  # fmt: skip
    t_all = np.concatenate([bjd(YEARS)] + [q.time for q in tabs.values()])
    jobs = []
    for k in order:
        if w[k] <= 0:
            continue
        mem = copy.deepcopy(base)
        for c in vary:
            v = mm.iloc[k][c]
            if not np.isfinite(v):
                continue
            node = mem
            parts = c.split(".")
            for p in parts[:-1]:
                node = node.setdefault(p, {})
            node[parts[-1]] = float(v)
        meta = {"member": int(mm.index[k]), "weight": float(w[k]),
                "P3_d": float(mm.iloc[k].get("orbits.out.P", np.nan)),
                "m3": float(mm.iloc[k].get("masses.m3", np.nan))}  # fmt: skip
        jobs.append((mem, t_all, t0_ens, tabs, meta))
    from concurrent.futures import ProcessPoolExecutor

    with ProcessPoolExecutor(max_workers=max(1, a.processes)) as ex:
        recs = [r for r in ex.map(_member, jobs) if r is not None]
    er = pd.DataFrame(recs)
    er.to_csv(out / "ensemble_dates.csv", index=False)

    def wq(col, q, ww):
        x = er[col].values
        ok = np.isfinite(x) & (ww > 0)
        if not ok.any():
            return None
        o_ = np.argsort(x[ok])
        c = np.cumsum(ww[ok][o_])
        return float(np.interp(q, c / c[-1], x[ok][o_]))

    def stats(ww):
        ww = np.asarray(ww, float)
        return {
            "effective_members": float(ww.sum() ** 2 / (ww**2).sum()) if ww.sum() > 0 else 0.0,
            "eclipsing_1933": float((ww * er.eclipsing_1933).sum() / ww.sum()),
            "end_before_2025_16_50_84": [
                wq("last_end_before_2025", q, ww) for q in (0.16, 0.5, 0.84)
            ],
            "share_with_end_1900_2025": float(
                ww[np.isfinite(er.last_end_before_2025)].sum() / ww.sum()
            ),
            "return_after_2025_16_50_84": [
                wq("first_start_after_2025", q, ww) for q in (0.16, 0.5, 0.84)
            ],
            "share_returning_by_2100": float(
                ww[np.isfinite(er.first_start_after_2025)].sum() / ww.sum()
            ),
        }

    res["ensemble"] = {
        "members_integrated": len(er),
        "weight_share_integrated": float(er.weight.sum() / w.sum()) if len(er) else 0.0,
        "weights_i0_and_rate_1980": stats(er.weight.values),
        **{
            f"weights_also_{n}": stats(er.weight.values * np.exp(-0.5 * er[f"chi2_{n}"].values))
            for n in tabs
        },
    }
    (out / "summary.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(res, indent=1, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
