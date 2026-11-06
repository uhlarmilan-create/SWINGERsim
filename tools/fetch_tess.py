"""Download the TESS light curves of the SQUADRA test systems for comparison with swingersim.

Needs lightkurve (``pip install lightkurve``) and internet access. For every system and
sector it writes ``data/observed/<system>/tess_s<NN>_<author>.csv`` with the columns

    bjd_tdb, flux, flux_err      (flux normalised to its median; BJD = BTJD + 2457000)

The 2-min SPOC light curves are taken when they exist, otherwise the FFI light curves of
TESS-SPOC or QLP. The published analyses (Pribulla et al. 2023, Borkovits et al. 2026,
Zasche et al. 2023) used their own FFI photometry, so depths can differ slightly from theirs.

Usage:  python tools/fetch_tess.py                 (all systems)
        python tools/fetch_tess.py bu_cmi          (one system)
        python tools/fetch_tess.py v907_sco qx_cas ah_cep   (part III)
"""

from __future__ import annotations

import sys
from pathlib import Path

SYSTEMS = {  # name: target for the MAST search; every sector found is downloaded
    "bu_cmi": "BU CMi",  # 7, 34, 61 used so far; newer ones test the predictions
    "tic120362137": "TIC 120362137",
    "v994_her": "V994 Her",  # 26, 40, 53 used so far
    # part III (2+1 with a changing inclination)
    "v907_sco": "V907 Sco",
    "qx_cas": "QX Cas",
    "ah_cep": "AH Cep",
}
AUTHORS = ["SPOC", "TESS-SPOC", "QLP"]


def fetch(name: str, root: Path) -> None:
    import lightkurve as lk
    import numpy as np

    target = SYSTEMS[name]
    out = root / "data" / "observed" / name
    out.mkdir(parents=True, exist_ok=True)
    res = lk.search_lightcurve(target, mission="TESS")
    if len(res) == 0:
        print(f"{name}: nothing found for '{target}'")
        return
    found = sorted({int(s) for s in res.table["sequence_number"]})
    print(f"{name}: sectors available and downloaded {found}")
    for sec in found:
        for author in AUTHORS:
            sel = res[(res.table["sequence_number"] == sec) & (res.table["author"] == author)]
            if len(sel) == 0:
                continue
            lc = sel[0].download(quality_bitmask="default")
            if lc is None:
                continue
            lc = lc.remove_nans().normalize()
            t = np.asarray(lc.time.value) + 2457000.0
            path = out / f"tess_s{sec:02d}_{author.lower()}.csv"
            np.savetxt(
                path,
                np.column_stack([t, np.asarray(lc.flux.value), np.asarray(lc.flux_err.value)]),
                delimiter=",",
                header="bjd_tdb,flux,flux_err",
                comments="",
                fmt="%.7f",
            )
            print(f"  sector {sec}: {author}, {len(t)} points -> {path}")
            break
        else:
            print(f"  sector {sec}: no light curve of {', '.join(AUTHORS)}")


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    names = sys.argv[1:] or list(SYSTEMS)
    for n in names:
        if n not in SYSTEMS:
            raise SystemExit(f"unknown system {n}; choose from {', '.join(SYSTEMS)}")
        fetch(n, root)


if __name__ == "__main__":
    main()
