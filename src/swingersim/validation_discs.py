"""Benchmarks of the dust discs as occulters (appended to swingersim.validation.BENCHMARKS)."""

from __future__ import annotations

import math
from functools import cache

import numpy as np

from ._benchmark import ROOT, Benchmark

JD_EPSAUR_T0 = 2454466.5  # t = 0 of configs/catalog/eps_aur.yaml (2008 Jan 1)


def _star_through(normal, r_out, r_in, half, tau, offset=(0.0, 0.0)):
    """Transmission of a uniform star of radius 1 at the origin behind a disc centred at
    `offset` (sky plane), the star far behind the disc."""
    from .discs import Disc, transmission
    from .photometry import StarLight

    n = np.asarray(normal, float)
    d = Disc(
        "d",
        [1],
        r_out,
        r_in,
        half,
        math.inf if tau is None else tau / (2 * half),
        n / np.linalg.norm(n),
        None,
        0.0,
    )
    star = StarLight("m1", 1.0, 1.0, np.zeros(4), 0.3, 1.0, 0.0)
    pos = np.array([[0.0, 0.0, 100.0], [offset[0], offset[1], 0.0]])
    return float(transmission([d], 0.0, pos, np.array([1.0, 1.0]), [star, star])[0])


def disc_strip_transmission() -> float:
    """An opaque edge-on slab, much wider than the star, of half thickness 0.3 star radii
    across the centre of a uniform star: transmitted fraction."""
    return _star_through((0.0, 1.0, 0.0), 100.0, 0.0, 0.3, None)


def strip_exact(h: float = 0.3) -> float:
    return 1.0 - (2.0 / math.pi) * (math.asin(h) + h * math.sqrt(1.0 - h * h))


def disc_face_on_transmission(tau: float = 0.7) -> float:
    return _star_through((0.0, 0.0, 1.0), 100.0, 0.0, 0.05, tau)


def disc_inclined_transmission(tau: float = 0.7, inc_deg: float = 60.0) -> float:
    i = math.radians(inc_deg)
    return _star_through((0.0, math.sin(i), math.cos(i)), 100.0, 0.0, 0.05, tau)


@cache
def _epsaur_eclipse() -> tuple[float, float, float]:
    """Mid-eclipse [JD], duration I-IV [d] and depth [mag in V] of the 2009-2011 eclipse."""
    from .config import load_config
    from .photometry import light_curve

    cfg = load_config(ROOT / "configs" / "catalog" / "eps_aur.yaml")
    t = np.arange(400.0, 1500.0, 1.0)
    lc = light_curve(cfg, t, band="V")
    mag = -2.5 * np.log10(lc.flux / lc.flux[0])
    ecl = t[mag > 0.01]
    return (
        JD_EPSAUR_T0 + 0.5 * (ecl.min() + ecl.max()),
        float(ecl.max() - ecl.min()),
        float(mag.max()),
    )


def hd98800_polar_disc_inclination() -> float:
    """Inclination to the sky [deg] of a disc perpendicular to the orbit of HD 98800 B along
    its eccentricity vector (polar alignment), from the published orbit."""
    from .config import load_config
    from .discs import build_discs
    from .hierarchy import build_tree

    cfg = load_config(ROOT / "configs" / "catalog" / "hd98800.yaml")
    root, _ = build_tree(cfg)
    n = build_discs(cfg, root)[0].normal0
    return math.degrees(math.acos(abs(float(n[2]))))


EPS = "Stefanik et al. (2010) AJ 139, 1254; Hoard et al. (2010) ApJ 714, 549; Stencel (2013)"
BENCHMARKS_DISCS = [
    Benchmark(
        id="disc-strip-uniform",
        tier="A",
        quantity="opaque edge-on slab (half thickness 0.3 R) across a uniform star: transmitted "
        "fraction",
        unit="",
        expected=strip_exact(),
        tolerance=5e-3,
        relative=False,
        source="geometry: 1 - (2/pi)(asin h + h sqrt(1 - h^2))",
        link="https://en.wikipedia.org/wiki/Circular_segment",
        rationale="the stellar disk is sampled by 4000 points (sunflower lattice); a sharp edge "
        "is resolved to ~1e-3",
        compute=disc_strip_transmission,
    ),
    Benchmark(
        id="disc-face-on-tau",
        tier="A",
        quantity="face-on disc of optical depth 0.7 in front of a star: transmitted fraction",
        unit="",
        expected=math.exp(-0.7),
        tolerance=1e-9,
        relative=False,
        source="exp(-tau)",
        link="https://en.wikipedia.org/wiki/Beer%E2%80%93Lambert_law",
        rationale="every line of sight crosses the full thickness: exact",
        compute=disc_face_on_transmission,
    ),
    Benchmark(
        id="disc-inclined-tau",
        tier="A",
        quantity="disc of optical depth 0.7 inclined 60 deg in front of a star: transmitted "
        "fraction",
        unit="",
        expected=math.exp(-1.4),
        tolerance=1e-9,
        relative=False,
        source="exp(-tau / cos i)",
        link="https://en.wikipedia.org/wiki/Beer%E2%80%93Lambert_law",
        rationale="the path through the slab is thickness / cos i: exact",
        compute=disc_inclined_transmission,
    ),
    Benchmark(
        id="epsaur-mid-eclipse",
        tier="D",
        quantity="epsilon Aur: mid-eclipse of 2009-2011 from the spectroscopic orbit",
        unit="JD",
        expected=2455429.5,
        tolerance=30.0,
        relative=False,
        source=EPS + "; contacts I 2009-08-16, IV 2011-08-26",
        link="https://doi.org/10.1088/0004-6256/139/4/1254",
        rationale="the periastron time is uncertain by 80 d (Stefanik et al.); the mid-eclipse "
        "from the contacts by ~10 d",
        compute=lambda: _epsaur_eclipse()[0],
    ),
    Benchmark(
        id="epsaur-duration",
        tier="D",
        quantity="epsilon Aur: duration of the eclipse, contacts I-IV",
        unit="d",
        expected=740.0,
        tolerance=0.15,
        relative=True,
        source=EPS,
        link="https://doi.org/10.1088/0004-637X/714/1/549",
        rationale="set by the disc length (7.3 AU +-10 %) and the orbital speed (masses of the "
        "low-mass model); a fuzzy disc edge lengthens the observed ingress",
        compute=lambda: _epsaur_eclipse()[1],
    ),
    Benchmark(
        id="epsaur-depth",
        tier="D",
        quantity="epsilon Aur: depth of the eclipse in V (opaque disc, 0.6 AU thick)",
        unit="mag",
        expected=0.75,
        tolerance=0.15,
        relative=False,
        source=EPS,
        link="https://doi.org/10.1088/0004-637X/714/1/549",
        rationale="depends on the disc thickness (+-10 %), its tilt and the inclination",
        compute=lambda: _epsaur_eclipse()[2],
    ),
    Benchmark(
        id="hd98800-polar-disc",
        tier="D",
        quantity="HD 98800 B: inclination to the sky of a polar disc (along the eccentricity "
        "vector of the B orbit)",
        unit="deg",
        expected=26.0,
        tolerance=3.0,
        relative=False,
        source="Kennedy et al. (2019) Nat. Astron. 3, 230; orbit Zuniga-Fernandez et al. (2026)",
        link="https://doi.org/10.1038/s41550-018-0667-x",
        rationale="the observed disc inclination is 26 deg; a polar disc follows from the "
        "published orbit of B without any fitting",
        compute=hd98800_polar_disc_inclination,
    ),
]
