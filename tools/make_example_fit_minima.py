"""Make configs/examples/example_fit_minima.csv: synthetic times of minima of a 2+2 for the stage-1
fit example (swingersim fit configs/examples/example_fit.yaml).

The minima are computed with swingersim.ocfit.synthetic_minima from the parameters below and then
thinned to randomly chosen cycles (seeded). Random cycles matter: minima at a fixed spacing of N
cycles cannot tell an apsidal rate dw from dw + 2 pi k / N (an alias), so a regular sampling would
not recover the apsidal motion of pair A.

Run:  python tools/make_example_fit_minima.py
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np

from swingersim import ocfit

OUT = Path(__file__).resolve().parents[1] / "configs" / "examples" / "example_fit_minima.csv"
OUTER = {"P3": 2200.0, "T3": 2455300.0, "e3": 0.35, "w3": 1.2}
PAIRS = {
    "A": {"terms": ["lite", "apsidal"], "i": 88.0, "T0": 2455000.12, "P": 6.7391, "e": 0.25,
          "w0": 0.74, "dw": 0.0015, "asini": 1.1},  # fmt: skip
    "B": {"terms": ["lite"], "i": 63.4, "T0": 2455000.40, "P": 0.9193, "asini": 0.95},
}
T_START, T_END = 2455000.0, 2461466.0  # 17.7 yr
N_CYCLES = {"A": 25, "B": 177}  # observed cycles per pair (each with a primary and a secondary)
SIGMA_D = 2e-4
SEED = 1


def main() -> None:
    rng = np.random.default_rng(SEED)
    with tempfile.TemporaryDirectory() as tmp:
        full = Path(tmp) / "all.csv"
        ocfit.synthetic_minima(full, PAIRS, OUTER, T_START, T_END, sigma=SIGMA_D, seed=SEED, every=1)
        lines = full.read_text(encoding="utf-8").splitlines()
    header, rows = lines[:2], lines[2:]
    keep = []
    for name, n_keep in N_CYCLES.items():
        own = [r for r in rows if r.startswith(name + ",")]
        cycles = len(own) // 2  # rows come as (primary, secondary) per cycle
        chosen = np.sort(rng.choice(cycles, size=n_keep, replace=False))
        for c in chosen:
            keep += own[2 * c : 2 * c + 2]
    pa, pb = PAIRS["A"], PAIRS["B"]
    with open(OUT, "w", encoding="utf-8") as f:
        f.write(
            "# synthetic minima for configs/examples/example_fit.yaml (tools/make_example_fit_minima.py,\n"
            "# swingersim.ocfit.synthetic_minima, randomly chosen cycles, seed 1, sigma 2e-4 d):\n"
            f"# outer orbit P3 = {OUTER['P3']:g} d, T3 = {OUTER['T3']:.1f}, e3 = {OUTER['e3']:g}, "
            f"w3 = {OUTER['w3']:g} rad;\n"
            f"# A: T0 = {pa['T0']}, P = {pa['P']} d, e = {pa['e']}, w0 = {pa['w0']} rad, "
            f"dw = {pa['dw']} rad/cycle, a sin i = {pa['asini']} AU;\n"
            f"# B: T0 = {pb['T0']}, P = {pb['P']} d, a sin i = {pb['asini']} AU\n"
        )
        f.write("\n".join(header + keep) + "\n")
    print(f"{len(keep)} minima written to {OUT}")


if __name__ == "__main__":
    main()
