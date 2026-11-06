"""Gaia astrometry of systems with a suspected unseen companion (Part III: QX Cas, V907 Sco,
AH Cep; and CzeV343, CzeV4315): RUWE, astrometric excess noise and its significance, DR2 vs
DR3 parallax and proper motion, and whether DR3 has a non-single-star (NSS) solution.

Run on a computer with internet access (standard library only):
    python tools/gaia_check.py [-o output/special/gaia_check.csv] [names ...]

Positions: known ones are built in, others are resolved with CDS Sesame (three mirrors), or given
as "name=ra,dec" in degrees. The Gaia archive (ESA TAP, else the ARI Heidelberg mirror) is queried
within 3 arcsec.
A large RUWE (> 1.4) or significant excess noise (sig > 2) hints at orbital motion of the
photocentre; a DR2-DR3 parallax difference of several sigma can have the same cause."""

import argparse
import csv
import io
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

GAIA = [  # ESA archive, then the mirror of ARI Heidelberg
    "https://gea.esac.esa.int/tap-server/tap/sync",
    "https://gaia.ari.uni-heidelberg.de/tap/sync",
]
SESAME = [
    "https://cds.unistra.fr/cgi-bin/nph-sesame/-oI/SV?",
    "https://vizier.cds.unistra.fr/viz-bin/nph-sesame/-oI/SV?",
    "http://vizier.cfa.harvard.edu/viz-bin/nph-sesame/-oI/SV?",
]
DEFAULT = ["QX Cas", "V907 Sco", "AH Cep", "CzeV343", "HD 228777"]  # HD 228777 = CzeV4315
# J2000 positions already in the repository (data headers, configs): no name resolver needed
KNOWN = {
    "QX Cas": (359.67979, 61.16097),
    "AH Cep": (341.97058, 65.06217),
    "BU CMi": (119.524529, 7.213475),
    "V994 Her": (276.941222, 24.697407),
    "CzeV343": (87.10005, 30.950997),  # V849 Aur, VarAstro (var.astro.cz/en/Stars/345)
    "HD 228777": (304.417917, 39.243111),  # CzeV4315, Zasche et al. (2025), Table 1
    # V907 Sco: no precise position at hand; 2MASS J17565485-3444584 lies 9" from it (Zasche et
    # al. 2023): search 15" around that star and keep the brightest source
    "V907 Sco": (269.228542, -34.749556, 15.0, True),
}
TIMEOUT = 30
AGENT = "Mozilla/5.0 (compatible; swingersim gaia_check)"


def _value(x: str):
    """A CSV field as int, float or None (empty)."""
    if x == "":
        return None
    for f in (int, float):
        try:
            return f(x)
        except ValueError:
            pass
    return x


def _reason(body: str) -> str:
    """The error message of a TAP answer: the QUERY_STATUS INFO of a VOTable, else the HTML title
    or the start of the text."""
    m = re.search(r'<INFO[^>]*name="QUERY_STATUS"[^>]*>(.*?)</INFO>', body, re.DOTALL)
    if m and m.group(1).strip():
        return " ".join(m.group(1).split())[:400]
    m = re.search(r"<title>(.*?)</title>", body, re.DOTALL | re.IGNORECASE)
    if m:
        return " ".join(m.group(1).split())[:200]
    return " ".join(re.sub(r"<[^>]+>", " ", body).split())[:200]


def tap(adql: str) -> list[dict]:
    """An ADQL query (CSV answer) on the first Gaia TAP service that answers."""
    data = urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": adql}
    ).encode()
    err = None
    for url in GAIA:
        host = url.split("/")[2]
        for attempt in range(4):  # the ESA archive answers 503 when busy: wait and retry
            req = urllib.request.Request(url, data=data, headers={"User-Agent": AGENT})
            try:
                with urllib.request.urlopen(req, timeout=120) as r:
                    text = r.read().decode("utf-8", "replace")
            except urllib.error.HTTPError as exc:  # the body says why (bad ADQL, overload ...)
                err = f"HTTP {exc.code}: {_reason(exc.read().decode('utf-8', 'replace'))}"
                print(f"  ({host}: {err})", file=sys.stderr)
                if exc.code in (502, 503, 504) and attempt < 3:
                    time.sleep(5 * 2**attempt)
                    continue
                break
            except OSError as exc:  # timeout, refused: try the next service
                err = str(exc)
                print(f"  ({host}: {exc})", file=sys.stderr)
                break
            rows = list(csv.DictReader(io.StringIO(text)))
            return [{k: _value(v) for k, v in r.items()} for r in rows]
    raise OSError(f"no Gaia TAP service answered ({err})")


def resolve(name: str):
    """(ra, dec) in degrees: "name=ra,dec" given on the command line, a position known to the
    repository, or the CDS Sesame name resolver (SIMBAD, then VizieR; three mirrors)."""
    if "=" in name:
        ra, dec = name.split("=", 1)[1].split(",")
        return float(ra), float(dec), 3.0, False
    if name in KNOWN:
        return KNOWN[name]
    for base in SESAME:
        try:
            with urllib.request.urlopen(base + urllib.parse.quote(name), timeout=TIMEOUT) as r:
                text = r.read().decode("utf-8", "replace")
        except OSError as exc:
            print(f"  ({base.split('/')[2]}: {exc}; trying the next resolver)", file=sys.stderr)
            continue
        for line in text.splitlines():
            if line.startswith("%J "):  # %J ra dec (J2000, degrees)
                parts = line.split()
                return float(parts[1]), float(parts[2]), 3.0, False
        return None
    raise OSError(f"no name resolver answered for {name!r}; give it as 'name=ra,dec' (degrees)")


def _try(adql: str) -> list[dict]:
    """An optional query (DR2 cross-match, NSS): an empty answer if no service has the table."""
    try:
        return tap(adql)
    except OSError as exc:
        print(f"  (optional query failed: {exc})", file=sys.stderr)
        return []


def gaia(ra: float, dec: float, radius_arcsec: float = 3.0, brightest: bool = False) -> list[dict]:
    """Gaia DR3 sources within the radius (all, or only the brightest), each with its DR2
    counterpart and its NSS two-body orbit if any."""
    r = radius_arcsec / 3600.0
    q3 = (
        "SELECT g.source_id, g.phot_g_mean_mag, g.parallax, g.parallax_error, g.pmra, g.pmdec, "
        "g.ruwe, g.astrometric_excess_noise, g.astrometric_excess_noise_sig, "
        "g.ipd_frac_multi_peak, g.radial_velocity, g.radial_velocity_error, g.non_single_star, "
        f"DISTANCE(POINT('ICRS', g.ra, g.dec), POINT('ICRS', {ra:.7f}, {dec:.7f})) * 3600 AS sep_arcsec "
        "FROM gaiadr3.gaia_source AS g "
        f"WHERE CONTAINS(POINT('ICRS', g.ra, g.dec), CIRCLE('ICRS', {ra:.7f}, {dec:.7f}, {r:.7f})) = 1 "
        "ORDER BY g.phot_g_mean_mag"
    )
    out = tap(q3)
    out = out[:1] if brightest else out
    for row in out:
        sid = row["source_id"]
        nb = _try(
            "SELECT d2.source_id AS dr2_id, d2.parallax AS dr2_parallax, d2.parallax_error AS "
            "dr2_parallax_error, d2.pmra AS dr2_pmra, d2.pmdec AS dr2_pmdec, "
            "d2.astrometric_excess_noise AS dr2_aen, d2.astrometric_excess_noise_sig AS "
            "dr2_aen_sig FROM gaiadr3.dr2_neighbourhood AS n JOIN gaiadr2.gaia_source AS d2 "
            f"ON d2.source_id = n.dr2_source_id WHERE n.dr3_source_id = {sid} "
            "ORDER BY n.angular_distance",
        )
        row.update(nb[0] if nb else {})
        nss = _try(f"SELECT nss_solution_type, period, eccentricity FROM "
                        f"gaiadr3.nss_two_body_orbit WHERE source_id = {sid}")  # fmt: skip
        row["nss_two_body"] = (
            "; ".join(f"{n['nss_solution_type']} P={n['period']}" for n in nss) or ""
        )
        if row.get("dr2_parallax") is not None and row.get("parallax") is not None:
            s = (row["parallax_error"] ** 2 + row["dr2_parallax_error"] ** 2) ** 0.5
            row["dplx_sigma"] = (row["parallax"] - row["dr2_parallax"]) / s
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("names", nargs="*", default=DEFAULT)
    ap.add_argument("-o", "--output", default="output/special/gaia_check.csv")
    a = ap.parse_args(argv)
    rows = []
    for name in a.names:
        try:
            pos = resolve(name)
        except OSError as exc:
            print(f"{name}: {exc}")
            continue
        name = name.split("=")[0]
        if pos is None:
            print(f"{name}: not found in SIMBAD")
            continue
        ra, dec, radius, brightest = (*pos, 3.0, False)[:4] if len(pos) == 2 else pos
        try:
            found = gaia(ra, dec, radius, brightest)
        except OSError as exc:
            print(f"{name}: {exc}")
            continue
        if not found:
            print(f"{name}: no Gaia DR3 source within {radius:g} arcsec of {ra:.5f} {dec:+.5f}")
        for g in found:
            g["name"] = name
            rows.append(g)
            print(f"{name}: G {g['phot_g_mean_mag']:.2f}, sep {g['sep_arcsec']:.2f}\", "
                  f"RUWE {g['ruwe']}, AEN {g['astrometric_excess_noise']} mas "
                  f"(sig {g['astrometric_excess_noise_sig']}), plx {g['parallax']} +- "
                  f"{g['parallax_error']} (DR2 {g.get('dr2_parallax')}, diff "
                  f"{g.get('dplx_sigma')} sigma), NSS {g['non_single_star']} "
                  f"{g['nss_two_body']}")  # fmt: skip
    if rows:
        import os

        os.makedirs(os.path.dirname(a.output) or ".", exist_ok=True)
        keys = sorted({k for r in rows for k in r})
        with open(a.output, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["name", *[k for k in keys if k != "name"]])
            w.writeheader()
            w.writerows(rows)
        print(f"saved {a.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
