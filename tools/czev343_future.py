"""CzeV343: when and how does the 2+2 quadruple fall apart? (special test after the campaign)

Pejcha et al. (2022, A&A 667, A53, Sect. 5): with COMPAS, most realisations of either binary expand
by mass transfer to a = 50-220 Rsun, while the quadruple stays dynamically stable only for
a < ~50 Rsun (Mardling & Aarseth 2001) - "dynamical instability caused by orbital expansion ...
has not been explored so far in 2+2 quadruples". The paper gives no times.

Stage 1 (this script, minutes): for many draws of the uncertain stellar parameters, the
mass-transfer histories of both binaries (swingersim.masstransfer: MIST tracks, BSE-like stability
and rates, common envelope) and, along them, the stability ratios of both nested triples (MA01
and Vynatheya et al. 2022; the other binary as a point mass, the outer orbit adiabatic,
a_out M_tot = const, e_out and the mutual inclinations kept). Result per draw: what happens to
each binary (merger, common envelope, stable transfer, nothing), the largest separation, and the
first time the quadruple becomes dynamically unstable - before or after the binary merges.

Stage 2 (written here, run on the PC with swingersim ensemble): for every scenario that becomes
unstable, an N-body configuration at that moment (masses, separations, radii and ages of the
histories; circular binaries after the transfer; the other, tight binary as one body) with the
phases drawn, integrated for 10^4 yr (~2500 outer orbits): ejection, collision, exchange or a new
stable hierarchy.

    python tools/czev343_future.py [--n 300] [-o output/special/czev343]
"""

import argparse
import contextlib
import json
import math
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from swingersim import masstransfer as mt
from swingersim import units
from swingersim.config import Config
from swingersim.hierarchy import build_tree, initialize
from swingersim.stability import ma01_ratio, nested_triples, v22_ratio

BASE = ROOT / "configs/special/czev343.yaml"
# Pejcha et al. (2022), Table 2: isochrone solutions (A1, A2, B1, B2) [Msun] and the age [Myr]
SOLUTIONS = {
    "fiducial": {"masses": (1.27, 1.76, 1.36, 1.20), "age": 630.0},
    "alternative": {"masses": (1.42, 2.48, 1.75, 1.38), "age": 140.0},
}
MASS_SIGMA = 0.10  # "roughly 10 % uncertainty" of the binary masses (Sect. 4.1)
T_END_YR = 4.0e9


def draw(rng, sol: str) -> dict:
    s = SOLUTIONS[sol]
    fA, fB = rng.normal(1.0, MASS_SIGMA, 2).clip(0.75, 1.25)  # the binary masses
    qA, qB = (rng.normal(1.0, 0.05, 2)).clip(0.85, 1.15)  # the mass ratios within the binary
    a1, a2, b1, b2 = s["masses"]
    mA, mB = (a1 + a2) * fA, (b1 + b2) * fB
    xA = a2 / (a1 + a2) * qA
    xB = b1 / (b1 + b2) * qB
    xA, xB = min(xA, 0.8), min(xB, 0.8)
    age = s["age"] * float(np.exp(rng.normal(0.0, 0.3)))  # isochrone ages are uncertain
    return {
        "solution": sol,
        "m1": mA * xA,  # A2, the heavier star of A (Table 2), = m1 here
        "m2": mA * (1 - xA),
        "m3": mB * xB,
        "m4": mB * (1 - xB),
        "age_myr": age,
        "Omega_A": float(rng.uniform(0, 360)),
        "Omega_B": float(rng.uniform(0, 360)),
    }


def config_of(d: dict) -> Config:
    raw = yaml.safe_load(BASE.read_text())
    raw["masses"] = {k: round(d[k], 5) for k in ("m1", "m2", "m3", "m4")}
    raw["stellar_model"]["age_myr"] = round(d["age_myr"], 3)
    model = d.get("model", "mist")  # --model (passed with the draw: the workers are spawned)
    if model != "mist":  # SSE counts the age from the ZAMS (MIST from the pre-main sequence)
        raw["stellar_model"] = {"kind": model, "age_myr": round(d["age_myr"], 3)}
    raw["orbits"]["A"]["Omega"] = round(d["Omega_A"], 3)
    raw["orbits"]["B"]["Omega"] = round(d["Omega_B"], 3)
    raw["run"]["t_max_yr"] = T_END_YR
    return Config.model_validate(raw)


def stage1(cfg: Config) -> dict:
    root, m = build_tree(cfg)
    pos, vel = initialize(root, m)
    tri = {t.inner: t for t in nested_triples(root, m, pos, vel)}
    nodes = {n.name: n for n in root.nodes()}
    a_out0, e_out = float(nodes["out"].elements.a), float(nodes["out"].elements.e)
    M0 = float(m.sum())
    H = mt.histories(cfg)
    t = np.unique(np.concatenate([H[n].t for n in ("A", "B")]))
    t = t[t <= T_END_YR]
    hA, hB = H["A"], H["B"]

    def end_of(h):  # a merger ends the binary as a binary
        return h.t_event if (h.event or "").endswith("merger") else math.inf

    out = {
        "A_event": hA.event,
        "A_t_event_myr": None if hA.t_event is None else hA.t_event / 1e6,
        "A_t_onset_myr": None if hA.t_onset is None else hA.t_onset / 1e6,
        "B_event": hB.event,
        "B_t_event_myr": None if hB.t_event is None else hB.t_event / 1e6,
        "B_t_onset_myr": None if hB.t_onset is None else hB.t_onset / 1e6,
        "A_a_max_rsun": float(hA.a.max() / units.RSUN_AU),
        "B_a_max_rsun": float(hB.a.max() / units.RSUN_AU),
        "A_regimes": hA.regimes_text(),
        "B_regimes": hB.regimes_text(),
        "ma01_A0": tri["A"].ma01,
        "ma01_B0": tri["B"].ma01,
        "t_unstable_myr": None,
        "unstable_pair": None,
    }
    for tt in t:
        mA, mB = sum(hA.masses_at(tt)), sum(hB.masses_at(tt))
        a_out = a_out0 * M0 / (mA + mB)  # a M = const (isotropic loss)
        for name, h, mi, mo, end in (("A", hA, mA, mB, end_of(hA)), ("B", hB, mB, mA, end_of(hB))):
            if tt >= end:
                continue
            tr = tri[name]
            i_mut = math.radians(tr.i_mut_deg)
            a_in = h.a_at(tt)
            r1 = ma01_ratio(a_in, a_out, e_out, mo / mi, i_mut)
            r2 = v22_ratio(a_in, a_out, 0.0, e_out, mo / mi, i_mut)
            if min(r1, r2) < 1.0 and out["t_unstable_myr"] is None:
                out.update(
                    t_unstable_myr=tt / 1e6, unstable_pair=name, a_in_unstable_rsun=a_in /
                    units.RSUN_AU, a_out_unstable_au=a_out, ma01_unstable=r1, v22_unstable=r2,
                )  # fmt: skip
        if out["t_unstable_myr"] is not None:
            break
    return out


def stage2_config(cfg: Config, res: dict, path: Path, k2_giant: float = 0.01) -> None:
    """N-body configuration at the onset of the instability (members vary phases and nodes).

    Tides and rotation are on: the expanded (stripped) donor raises the tidal and rotational
    bulges whose apsidal precession can quench the von Zeipel-Lidov-Kozai cycles. Its
    apsidal-motion constant is uncertain (a giant with a thin envelope over a dense core), so it
    is a parameter (k2_giant, bracketed by the runs); rg2 comes from the stellar state."""
    tt = res["t_unstable_myr"] * 1e6
    H = mt.histories(cfg)
    base = yaml.safe_load(BASE.read_text())
    keep = res["unstable_pair"]  # the expanded binary stays a binary;
    other = "B" if keep == "A" else "A"  # the other (a_B / a_out ~ 0.005) is one body
    names = {keep: ("m1", "m2")}  # stars renumbered: the kept pair m1 + m2, the other m3
    m_k = H[keep].masses_at(tt)
    m_o = sum(H[other].masses_at(tt))
    point = "m3"
    o_k = dict(base["orbits"][keep])
    o_k.pop("P", None)
    o_k.update(a=round(H[keep].a_at(tt), 7), e=0.0 if H[keep].t_onset is not None else o_k["e"],
               Omega=round(res[f"Omega_{keep}"], 3))  # fmt: skip
    o_out = dict(base["orbits"]["out"])
    o_out.pop("P", None)
    o_out["a"] = round(res["a_out_unstable_au"], 6)
    stars = {}
    for k, nm in enumerate(names[keep]):
        st = H[keep].state_at(cfg.stellar_model, k, tt)
        stars[nm] = {"R": round(float(st.radius_rsun), 4), "Teff": round(float(st.teff_k), 1)}
        if st.teff_k < 5500.0 and st.radius_rsun > 5.0:  # the expanded donor
            stars[nm]["apsidal_k2"] = k2_giant
            # a stripped giant: dense He core + thin convective envelope; I/(M R^2) ~ 0.2 M_env/M
            # (the core adds nothing at R_core << R), at least 1e-3
            f_env = max(0.0, (m_k[k] - float(st.m_core)) / m_k[k])
            stars[nm]["rg2"] = round(min(0.4, max(1e-3, float(st.rg2), 0.2 * f_env)), 5)
    # the other binary as one body: its size is its separation (collisions with it)
    st = [H[other].state_at(cfg.stellar_model, k, tt) for k in (0, 1)]
    r_o = max(H[other].a_at(tt) / units.RSUN_AU, *(float(x.radius_rsun) for x in st))
    # a binary has no tidal bulge of its own: negligible k2 (its quadrupole is not modelled)
    stars[point] = {"R": round(r_o, 3), "Teff": round(max(float(x.teff_k) for x in st), 1),
                    "apsidal_k2": 1e-6, "lag_s": 0.0}  # fmt: skip
    raw = {
        "topology": "custom",
        "hierarchy": {
            "out": [keep, point] if keep == "A" else [point, keep],
            keep: list(names[keep]),
        },
        "masses": {
            names[keep][0]: round(m_k[0], 5),
            names[keep][1]: round(m_k[1], 5),
            point: round(m_o, 5),
        },
        "orbits": {keep: o_k, "out": o_out},
        "stars": stars,
    }
    raw["physics"] = {"gr": True, "tides": True, "rotation": True}
    raw["run"] = {"backend": "nbody", "t_max_yr": 1.0e4, "n_output": 2001, "collisions": True,
                  "max_wall_s": 1800}  # fmt: skip
    raw["ensemble"] = {
        "n": 100,
        "seed": 343,
        "directions": ["forward"],
        "vary": {
            f"orbits.{keep}.M0": "angle",
            "orbits.out.M0": "angle",
            f"orbits.{keep}.omega": "angle",
        },
    }
    head = (
        f"# CzeV343 at the onset of dynamical instability (tools/czev343_future.py, stage 2):\n"
        f"# draw {res['draw']} ({res['solution']} solution), {tt / 1e6:.1f} Myr from now; binary "
        f"{res['unstable_pair']} has a = {res['a_in_unstable_rsun']:.1f} Rsun\n"
        f"# (MA01 {res['ma01_unstable']:.2f}, V22 {res['v22_unstable']:.2f}). N-body, 1e4 yr "
        f"(~2500 outer orbits); members vary the phases and arguments of periastron.\n"
        f"# GR, tides and rotation; apsidal k2 of the expanded donor = {k2_giant:g}.\n"
        f"# Stars renumbered: m1 + m2 = binary {keep}, m3 = binary {other} as one body (R = its separation).\n"
    )
    path.write_text(head + yaml.safe_dump(raw, sort_keys=False))


def _one(d: dict) -> dict:
    warnings.simplefilter("ignore")
    try:
        r = stage1(config_of(d))
    except Exception as exc:  # noqa: BLE001 - a draw beyond the tracks
        r = {"error": f"{type(exc).__name__}: {exc}"[:200]}
    return {**d, **r}


def _run_all(draws: list, processes: int, timeout: float) -> list:
    """Each draw in its own process with a time limit (a few histories stall in the engine)."""
    import multiprocessing as mp
    import time

    ctx = mp.get_context("spawn")
    rows, pending, done = {}, list(enumerate(draws)), 0
    running = {}  # index -> (process, queue, start)
    while pending or running:
        while pending and len(running) < processes:
            i, d = pending.pop(0)
            q = ctx.Queue()
            pr = ctx.Process(target=_worker, args=(d, q), daemon=True)
            pr.start()
            running[i] = (pr, q, time.time())
        time.sleep(0.2)
        for i, (pr, q, t0) in list(running.items()):
            res = None
            with contextlib.suppress(Exception):  # an empty queue
                res = q.get_nowait()
            if res is not None:
                rows[i] = res
            elif time.time() - t0 > timeout:
                pr.kill()
                rows[i] = {**draws[i], "error": f"timeout {timeout:.0f} s"}
            elif not pr.is_alive() and pr.exitcode not in (0, None):
                rows[i] = {**draws[i], "error": f"exit {pr.exitcode}"}
            else:
                continue
            pr.join(1)
            del running[i]
            done += 1
            if done % 25 == 0:
                print(f"{done}/{len(draws)} draws", flush=True)
    return [rows[i] for i in range(len(draws))]


def _worker(d: dict, q) -> None:
    q.put(_one(d))


def _write_stage2(df, out: Path, n: int, k2s) -> None:
    """N-body configurations for a few unstable draws (spread over the onset times), one set
    per apsidal constant of the expanded donor."""
    uns = df[df.t_unstable_myr.notna()].sort_values("t_unstable_myr")
    if not len(uns):
        return
    pick = uns.iloc[np.unique(np.linspace(0, len(uns) - 1, min(n, len(uns))).astype(int))]
    s2 = out / "stage2"
    s2.mkdir(exist_ok=True)
    for k2 in k2s:
        for _, r in pick.iterrows():
            res = r.to_dict()
            stage2_config(config_of(res), res, s2 / f"{res['draw']}_k2_{k2:g}.yaml", k2)
    print(f"stage 2: {len(pick)} x {len(k2s)} N-body ensembles in {s2}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n", type=int, default=300, help="draws per isochrone solution")
    ap.add_argument("-o", "--output", default=str(ROOT / "output/special/czev343"))
    ap.add_argument("--seed", type=int, default=343)
    ap.add_argument("--stage2", type=int, default=6, help="N-body configurations to write")
    ap.add_argument("--processes", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--timeout", type=float, default=300.0, help="limit per draw [s]")
    ap.add_argument(
        "--stage2-only",
        action="store_true",
        help="reuse stage1.csv of the output folder; only write the N-body configs",
    )
    ap.add_argument(
        "--k2-giant",
        type=float,
        nargs="+",
        default=[0.001, 0.01, 0.05],
        help="apsidal-motion constants of the expanded donor (one config set each): 0.001 for "
        "a core-dominated stripped giant, 0.05 an upper bound",
    )
    ap.add_argument(
        "--model",
        choices=("mist", "sse"),
        default="mist",
        help="stellar tracks (mist: MIST v1.2 as in Part IV; sse: Hurley et al. 2000); the same "
        "seed gives the same draws, so two runs compare the track families member by member",
    )
    ap.add_argument(
        "--solutions",
        nargs="+",
        choices=tuple(SOLUTIONS),
        default=list(SOLUTIONS),
        help="isochrone solutions of Pejcha et al. (2022) to draw from",
    )
    a = ap.parse_args(argv)
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    if a.stage2_only:
        import pandas as pd

        _write_stage2(pd.read_csv(out / "stage1.csv"), out, a.stage2, a.k2_giant)
        return 0
    rng = np.random.default_rng(a.seed)
    draws = []
    for sol in SOLUTIONS:
        for k in range(a.n):
            if sol not in a.solutions:
                draw(rng, sol)  # keep the random stream of the other solutions unchanged
                continue
            d = draw(rng, sol)
            d["draw"] = f"{sol}_{k}"
            d["model"] = a.model
            draws.append(d)
    warnings.simplefilter("ignore")
    rows = _run_all(draws, a.processes, a.timeout)
    import pandas as pd

    df = pd.DataFrame(rows)
    df.to_csv(out / "stage1.csv", index=False)
    lines = []
    for sol, g in df.groupby("solution"):
        ok = g[g.get("error").isna()] if "error" in g else g
        n = len(ok)
        uns = ok[ok.t_unstable_myr.notna()]
        lines.append(
            f"{sol}: {n} draws ({len(g) - n} not finished: beyond the tracks or over the time limit)"
        )
        for p in ("A", "B"):
            ev = ok[f"{p}_event"].fillna("none").value_counts()
            lines.append(f"  binary {p}: " + ", ".join(f"{k} {v / n:.0%}" for k, v in ev.items()))
            lines.append(f"    largest separation [Rsun]: median {ok[f'{p}_a_max_rsun'].median():.1f}, "
                         f"> 50 Rsun in {(ok[f'{p}_a_max_rsun'] > 50).mean():.0%}")  # fmt: skip
        lines.append(f"  dynamically unstable: {len(uns) / n:.0%}" + (
            f", after {uns.t_unstable_myr.quantile(0.16):.0f}-{uns.t_unstable_myr.median():.0f}-"
            f"{uns.t_unstable_myr.quantile(0.84):.0f} Myr (16-50-84 %), pair "
            + ", ".join(f"{k} {v}" for k, v in uns.unstable_pair.value_counts().items())
            if len(uns) else ""))  # fmt: skip
    text = "\n".join(lines)
    print(text)
    (out / "stage1_summary.txt").write_text(text + "\n")
    _write_stage2(df, out, a.stage2, a.k2_giant)
    (out / "stage1_meta.json").write_text(
        json.dumps({"n": a.n, "seed": a.seed, "model": a.model, "solutions": a.solutions})
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
