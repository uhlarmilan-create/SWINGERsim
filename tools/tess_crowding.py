"""Crowding correction of a TESS SPOC light curve: is it why the eclipses are deeper than in QLP?

    python tools/tess_crowding.py "V994 Her" 26 data/observed/v994_her

Downloads the 2-min SPOC light curve of the sector (lightkurve, internet), prints the keywords
CROWDSAP (fraction of the flux in the aperture that belongs to the target, as SPOC estimates it)
and FLFRCSAP (fraction of the target flux inside the aperture), and writes the simple-aperture
flux (SAP, without the crowding correction) next to the PDCSAP light curve:

    <dir>/tess_sNN_spoc_sap.csv     bjd_tdb, flux, flux_err   (SAP flux normalised to its median)

PDCSAP = (SAP - (1 - CROWDSAP) * median(SAP)) / FLFRCSAP (after the systematics correction): a
CROWDSAP of 0.55 would make every eclipse 1/0.55 = 1.8 times deeper - what swingersim sees in V994 Her
sector 26 (depth scale 2.0 against 1.2 in the QLP sectors).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 3:
        print(__doc__)
        return 1
    import lightkurve as lk

    target, sector, out = argv[0], int(argv[1]), Path(argv[2])
    res = lk.search_lightcurve(target, mission="TESS", author="SPOC", sector=sector)
    if len(res) == 0:
        print(f"no SPOC light curve of {target} in sector {sector}")
        return 1
    lc = res[0].download(quality_bitmask="default", flux_column="sap_flux")
    h = lc.meta
    print(
        f"{target} sector {sector}: CROWDSAP = {h.get('CROWDSAP')}, FLFRCSAP = {h.get('FLFRCSAP')}, "
        f"TESSMAG = {h.get('TESSMAG')}, aperture pixels = {h.get('NUMPIX', '?')}"
    )
    lc = lc.remove_nans().normalize()

    def plain(x):  # astropy Masked arrays -> float array with NaN where masked
        v = getattr(x, "value", x)
        data = np.array(getattr(v, "unmasked", v), dtype=float)
        mask = np.asarray(getattr(v, "mask", False), dtype=bool)
        if mask.any():
            data[mask] = np.nan
        return data

    t = plain(lc.time) + 2457000.0
    f, e = plain(lc.flux), plain(lc.flux_err)
    ok = np.isfinite(t) & np.isfinite(f) & np.isfinite(e)
    t, f, e = t[ok], f[ok], e[ok]
    path = out / f"tess_s{sector:02d}_spoc_sap.csv"
    np.savetxt(
        path,
        np.column_stack([t, f, e]),
        delimiter=",",
        header=f"# SAP flux of SPOC, sector {sector}; CROWDSAP {h.get('CROWDSAP')}, "
        f"FLFRCSAP {h.get('FLFRCSAP')} (tools/tess_crowding.py)\nbjd_tdb,flux,flux_err",
        comments="",
        fmt="%.7f",
    )
    print(f"written {path} ({len(t)} points)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
