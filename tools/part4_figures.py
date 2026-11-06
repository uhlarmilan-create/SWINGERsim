"""Figures and numbers of Testing SWINGERsim, Part IV (CzeV343 and CzeV4315).

    python tools/part4_figures.py [--root .] [--out docs/figures/part4]

Reads output/special/czev343/{stage1.csv, stage2/, s2tides/*/summary.parquet} and
output/special/czev4315_2p2{,p1_wide}/{rows.jsonl, ensemble.yaml}. The example histories of
one onset (draw fiducial_7) are integrated here for 1000 yr if they are missing (seconds).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6"  # categorical slots
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e0"
K2 = (0.001, 0.01, 0.05)
K2_COLOR = dict(zip(K2, (C2, C1, C3), strict=True))


def _style(plt):
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.edgecolor": MUTED,
            "axes.labelcolor": INK,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "savefig.dpi": 160,
        }
    )


def wquant(x, w, ps=(0.16, 0.5, 0.84)):
    x, w = np.asarray(x, float), np.asarray(w, float)
    o = np.argsort(x)
    cw = np.cumsum(w[o]) / w.sum()
    return [float(x[o][min(np.searchsorted(cw, p), len(x) - 1)]) for p in ps]


# ------------------------------------------------------------------ CzeV343
def czev343(root: Path, out: Path, plt, numbers: dict) -> None:
    base = root / "output/special/czev343"
    s1 = pd.read_csv(base / "stage1.csv")
    ok = s1[s1["error"].isna()]
    labels = {
        "common_envelope_merger": "merger\n(common env.)",
        "contact": "contact",
        "donor_stripped": "donor\nstripped",
        "none": "no inter-\naction",
    }
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.4))
    n1 = {}
    for k, (sol, col) in enumerate((("fiducial", C1), ("alternative", C2))):
        g = ok[ok.solution == sol]
        f = g.A_event.fillna("none").value_counts(normalize=True)
        x = np.arange(len(labels)) + (k - 0.5) * 0.38
        ax[0].bar(x, [f.get(e, 0) * 100 for e in labels], 0.36, color=col,
                  label=f"{sol} solution", edgecolor="white", linewidth=1)  # fmt: skip
        amax = g.A_a_max_rsun.dropna()
        ax[1].hist(np.log10(amax), bins=np.linspace(0.5, 2.0, 31), color=col, alpha=0.75,
                   label=f"{sol} solution", edgecolor="white", linewidth=0.5)  # fmt: skip
        uns = g[g.t_unstable_myr.notna()]
        n1[sol] = {
            "draws": len(g),
            "not_finished": int((s1.solution == sol).sum() - len(g)),
            "A_events": {e: round(float(f.get(e, 0)), 3) for e in labels},
            "a_max_gt50": round(float((amax > 50).mean()), 3),
            "unstable": round(len(uns) / len(g), 3),
            "t_unstable_myr": [
                round(float(uns.t_unstable_myr.quantile(p))) for p in (0.16, 0.5, 0.84)
            ]
            if len(uns)
            else None,
            "a_unstable_rsun": [
                round(float(uns.a_in_unstable_rsun.min()), 1),
                round(float(uns.a_in_unstable_rsun.max()), 1),
            ]
            if len(uns)
            else None,
            "ma01_unstable": round(float(uns.ma01_unstable.median()), 2) if len(uns) else None,
            "v22_unstable": [
                round(float(uns.v22_unstable.min()), 2),
                round(float(uns.v22_unstable.max()), 2),
            ]
            if len(uns)
            else None,
            "merger_t_myr": [
                round(float(g[g.A_event == "common_envelope_merger"].A_t_event_myr.quantile(p)))
                for p in (0.16, 0.5, 0.84)
            ]
            if (g.A_event == "common_envelope_merger").any()
            else None,
        }
        if len(uns):
            ax[2].hist(uns.t_unstable_myr / 1e3, bins=np.linspace(1.0, 2.2, 25), color=col,
                       alpha=0.85, label=f"{sol}: {len(uns)} of {len(g)} draws", edgecolor="white")  # fmt: skip
    ax[0].set_xticks(range(len(labels)), labels.values(), fontsize=7.5)
    ax[0].set_ylabel("share of draws [%]")
    ax[0].set_title("(a) what happens to binary A", fontsize=9, loc="left")
    ax[1].axvline(np.log10(50), color=MUTED, ls="--", lw=1)
    ax[1].text(
        np.log10(50) + 0.02, ax[1].get_ylim()[1] * 0.9, "50 R$_\\odot$", color=MUTED, fontsize=8
    )
    ax[1].set_xlabel("log largest separation of A  [R$_\\odot$]")
    ax[1].set_ylabel("draws")
    ax[1].set_title("(b) how far binary A widens", fontsize=9, loc="left")
    ax[2].set_xlabel("time from now [Gyr]")
    ax[2].set_ylabel("draws")
    ax[2].set_title("(c) onset of dynamical instability", fontsize=9, loc="left")
    for a in ax:
        a.legend(fontsize=7.5)
    fig.tight_layout()
    fig.savefig(out / "czev343_stage1.png")
    plt.close(fig)
    numbers["czev343_stage1"] = n1

    # stage 2: ensembles at the onset, three apsidal constants of the donor
    rows = []
    for f in sorted((base / "s2tides").glob("*_k2_*/summary.parquet")):
        draw, k2 = f.parent.name.rsplit("_k2_", 1)
        d = pd.read_parquet(f)
        col = d[d.fate != "intact"]
        rows.append({"draw": draw, "k2": float(k2), "n": len(d),
                     "collided": float((d.fate == "collision").mean()),
                     "t50": float(col.event_time_yr.median()) if len(col) else np.nan,
                     "t16": float(col.event_time_yr.quantile(0.16)) if len(col) else np.nan,
                     "t84": float(col.event_time_yr.quantile(0.84)) if len(col) else np.nan,
                     "eA": d.A_e_max.to_numpy(), "tcol": col.event_time_yr.to_numpy(), "imut0": float(d.imut_out_A_0.median()),
                     "err": float(d.max_energy_error.max())})  # fmt: skip
    s2 = pd.DataFrame(rows)
    ex = {}
    for k2 in K2:
        ts = base.parent / "czev343_examples" / f"fiducial_7_k2_{k2:g}_timeseries.csv"
        if not ts.exists():
            ts.parent.mkdir(parents=True, exist_ok=True)
            cfg = base / "stage2" / f"fiducial_7_k2_{k2:g}.yaml"
            subprocess.run([sys.executable, "-m", "swingersim", "run", str(cfg), "-o",
                            str(ts.with_name(ts.name.replace("_timeseries.csv", ".h5"))),
                            "--t-max-yr", "1000", "--no-viewer", "--no-plots"],
                           check=True, capture_output=True)  # fmt: skip
        ex[k2] = pd.read_csv(ts, comment="#")
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.4), gridspec_kw={"width_ratios": [1.5, 1, 1]})
    for k2 in K2:
        d = ex[k2]
        ax[0].plot(d.t_yr, d.e_A, color=K2_COLOR[k2], lw=1.2, label=f"k$_2$ = {k2:g}")
        if d.t_yr.iloc[-1] < 999:
            ax[0].plot(d.t_yr.iloc[-1], d.e_A.iloc[-1], "x", color=K2_COLOR[k2], ms=8, mew=2)
    ax[0].set_xlim(0, 400)
    ax[0].set_xlabel("time after the onset [yr]")
    ax[0].set_ylabel("eccentricity of binary A")
    ax[0].set_title("(a) one onset (fiducial draw 7): × = collision", fontsize=9, loc="left")
    ax[0].legend(fontsize=7.5)
    bins = np.geomspace(30, 1e4, 31)
    for k2 in K2:
        t = np.concatenate(s2[s2.k2 == k2].tcol.to_list())
        if len(t):
            ax[1].hist(t, bins=bins, color=K2_COLOR[k2], alpha=0.85, edgecolor="white",
                       label=f"k$_2$ = {k2:g}: {len(t)} of 600")  # fmt: skip
        parts = ax[2].violinplot([np.concatenate(s2[s2.k2 == k2].eA.to_list())], [K2.index(k2)],
                                 showmedians=True, widths=0.7)  # fmt: skip
        for b in parts["bodies"]:
            b.set_facecolor(K2_COLOR[k2])
            b.set_alpha(0.7)
        for key in ("cmedians", "cmins", "cmaxes", "cbars"):
            parts[key].set_color(MUTED)
    ax[1].set_xscale("log")
    ax[1].set_xlabel("time of the collision after the onset [yr]")
    ax[1].set_ylabel("members")
    ax[1].set_title("(b) when the collisions happen (6 onsets)", fontsize=9, loc="left")
    ax[1].legend(fontsize=7.5)
    ax[2].set_xticks(range(len(K2)), [f"{k:g}" for k in K2])
    ax[2].set_xlabel("apsidal constant k$_2$ of the donor")
    ax[2].set_ylabel("largest eccentricity of A")
    ax[2].set_title("(c) all 600 members per k$_2$", fontsize=9, loc="left")
    fig.tight_layout()
    fig.savefig(out / "czev343_stage2.png")
    plt.close(fig)
    numbers["czev343_stage2"] = {
        f"{k2:g}": {
            "collided": round(float(s2[s2.k2 == k2].collided.mean()), 3),
            "t50_range_yr": [
                round(float(s2[s2.k2 == k2].t50.min())),
                round(float(s2[s2.k2 == k2].t50.max())),
            ]
            if s2[s2.k2 == k2].t50.notna().any()
            else None,
            "t16_84_yr": [
                round(float(s2[s2.k2 == k2].t16.min())),
                round(float(s2[s2.k2 == k2].t84.max())),
            ]
            if s2[s2.k2 == k2].t50.notna().any()
            else None,
            "eA_max_median": round(
                float(np.median(np.concatenate(s2[s2.k2 == k2].eA.to_list()))), 3
            ),
            "eA_max_q84": round(
                float(np.quantile(np.concatenate(s2[s2.k2 == k2].eA.to_list()), 0.84)), 3
            ),
            "max_energy_error": float(s2[s2.k2 == k2].err.max()),
        }
        for k2 in K2
    }
    numbers["czev343_stage2"]["imut0_median_deg"] = [
        round(float(s2.imut0.min()), 1),
        round(float(s2.imut0.max()), 1),
    ]


# ------------------------------------------------------------------ CzeV4315
def _ensemble(path: Path):
    from swingersim.ensemble import observation_weights

    with open(path / "rows.jsonl") as fh:
        rows = [json.loads(line) for line in fh]
    with open(path / "ensemble.yaml") as fh:
        obs = yaml.safe_load(fh)["config"]["ensemble"]["observed"]
    obs = {
        k: {kk: vv for kk, vv in v.items() if vv is not None and vv is not False}
        for k, v in obs.items()
    }
    w = observation_weights(rows, obs)
    df = pd.DataFrame([r for r in rows if r["direction"] == "forward"]).set_index("member")
    df["w"] = pd.Series(w).reindex(df.index).fillna(0.0)
    ws = df.w.to_numpy()
    z = float(ws.mean())
    neff = float(ws.sum() ** 2 / (ws**2).sum()) if ws.sum() > 0 else 0.0
    return df, z, neff


def czev4315(root: Path, out: Path, plt, numbers: dict) -> None:
    base = root / "output/special"
    e22, z22, n22 = _ensemble(base / "czev4315_2p2")
    e5, z5, n5 = _ensemble(base / "czev4315_2p2p1_wide")
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.4), gridspec_kw={"width_ratios": [1.4, 1]})
    bins = np.linspace(-6, 3, 46)
    for d, col, lab in ((e22, C2, "2+2"), (e5, C1, "2+(2+1)")):
        x = d.B_dtiltdt0_deg_yr.clip(-6, 3)
        ax[0].hist(x, bins=bins, color=col, alpha=0.65, label=lab, edgecolor="white", linewidth=0.4)
    ax[0].axvspan(-1.4 - 0.4, -1.4 + 0.4, color=MUTED, alpha=0.18, lw=0)
    ax[0].axvline(-1.4, color=INK, lw=1)
    ax[0].text(
        -1.35, ax[0].get_ylim()[1] * 0.92, "observed\n−1.4 ± 0.2 °/yr (±2σ band)", fontsize=7.5
    )
    ax[0].set_xlabel("rate of change of |i$_B$ − 90°| today [deg/yr] (clipped)")
    ax[0].set_ylabel("members")
    ax[0].set_title("(a) tilt rate of binary B in the two hierarchies", fontsize=9, loc="left")
    ax[0].legend(fontsize=7.5, loc="upper left")
    fates = ["intact", "collision", "ejection", "exchange", "reconfiguration", "breakup"]
    for j, (d, col, lab) in enumerate(((e22, C2, "2+2"), (e5, C1, "2+(2+1)"))):
        f = d.fate.value_counts(normalize=True)
        y = np.arange(len(fates)) + (j - 0.5) * 0.38
        ax[1].barh(y, [f.get(k, 0) * 100 for k in fates], 0.36, color=col, label=lab)
    ax[1].set_yticks(range(len(fates)), fates, fontsize=8)
    ax[1].invert_yaxis()
    ax[1].set_xlabel("members [%]")
    ax[1].set_title("(b) fate in the forward run", fontsize=9, loc="left")
    ax[1].legend(fontsize=7.5, loc="lower right")
    fig.tight_layout()
    fig.savefig(out / "czev4315_hierarchies.png")
    plt.close(fig)

    g = e5[e5.w > 0.01 * e5.w.max()]
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.5))
    sc = ax[0].scatter(g["masses.m5"], g["orbits.Bt.e"], c=g.w, s=10 + 40 * g.w / g.w.max(),
                       cmap="Blues", vmin=0, vmax=g.w.max(), edgecolors=MUTED, linewidths=0.3)  # fmt: skip
    ax[0].set_xscale("log")
    ax[0].set_xlabel("mass of the third star of B, m$_5$ [M$_\\odot$]")
    ax[0].set_ylabel("eccentricity of its 119-d orbit")
    ax[0].set_title("(a) members that match the data", fontsize=9, loc="left")
    fig.colorbar(sc, ax=ax[0], label="observation weight")
    ax[1].hist(
        e5["masses.m5"], bins=np.geomspace(0.3, 6, 25), color=GRID, edgecolor="white", label="prior"
    )
    ax[1].hist(e5["masses.m5"], bins=np.geomspace(0.3, 6, 25), weights=e5.w / e5.w.sum() * len(e5),
               color=C1, alpha=0.8, edgecolor="white", label="weighted by the data")  # fmt: skip
    ax[1].set_xscale("log")
    for a_ in ax:
        a_.set_xticks([0.3, 0.5, 1, 2, 3, 5], ["0.3", "0.5", "1", "2", "3", "5"])
        a_.minorticks_off()
    ax[1].set_xlabel("m$_5$ [M$_\\odot$]")
    ax[1].set_ylabel("members (weighted, same total)")
    ax[1].set_title("(b) m$_5$: prior and posterior", fontsize=9, loc="left")
    ax[1].legend(fontsize=7.5)
    fig.tight_layout()
    fig.savefig(out / "czev4315_third_star.png")
    plt.close(fig)

    def q(c, d=e5):
        return [round(v, 3) for v in wquant(d[c].fillna(1e9), d.w)]

    numbers["czev4315"] = {
        "2+2": {
            "members": len(e22),
            "Z": z22,
            "n_eff": round(n22, 2),
            "fates": {
                k: round(float(v), 3) for k, v in e22.fate.value_counts(normalize=True).items()
            },
            "tilt_match_2sigma": int(((e22.B_dtiltdt0_deg_yr + 1.4).abs() < 0.4).sum()),
        },
        "2+(2+1)": {
            "members": len(e5),
            "Z": z5,
            "n_eff": round(n5, 2),
            "fates": {
                k: round(float(v), 3) for k, v in e5.fate.value_counts(normalize=True).items()
            },
            "tilt_match_2sigma": int(((e5.B_dtiltdt0_deg_yr + 1.4).abs() < 0.4).sum()),
            "m5": q("masses.m5"),
            "m5_95": [round(v, 2) for v in wquant(e5["masses.m5"], e5.w, (0.025, 0.975))],
            "P_m5_gt2": round(float(e5.w[e5["masses.m5"] > 2].sum() / e5.w.sum()), 3),
            "Bt_e": q("orbits.Bt.e"),
            "Bt_P": q("orbits.Bt.P"),
            "imut_Bt_B": q("imut_Bt_B_0"),
            "out_P_d": q("orbits.out.P"),
            "B_eclipse_change_yr": q("B_eclipse_change_yr"),
            "A_eclipse_change_yr": q("A_eclipse_change_yr"),
            "weighted_fates": {
                k: round(float(v), 3) for k, v in (e5.groupby("fate").w.sum() / e5.w.sum()).items()
            },
        },
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--root", default=str(ROOT))
    ap.add_argument("--out", default=str(ROOT / "docs/figures/part4"))
    a = ap.parse_args(argv)
    import warnings

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    warnings.simplefilter("ignore")
    _style(plt)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    numbers: dict = {}
    czev343(Path(a.root), out, plt, numbers)
    czev4315(Path(a.root), out, plt, numbers)
    (out / "numbers.json").write_text(json.dumps(numbers, indent=1, default=float))
    print(json.dumps(numbers, indent=1, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
