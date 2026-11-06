# Data and methods

## How SWINGERsim models a system

SWINGERsim describes a system as a hierarchy of orbits: here a close pair of stars (the eclipsing
binary) and a third star orbiting the pair. From the masses, the orbits and the stellar radii and
temperatures it integrates the motion of all stars with full gravity (N-body; REBOUND with the IAS15
integrator, Rein & Liu 2012, Rein & Spiegel 2015). The motion includes general relativity and the
tidal and rotational deformation of the stars (REBOUNDx, Tamayo et al. 2020), which slowly turn
eccentric orbits (apsidal motion). From the integrated motion SWINGERsim computes the times of the eclipses,
the radial velocities of the stars, the inclination of the pair at any date, and synthetic light curves
in a given passband (here TESS, Johnson B and V). The eclipse depths follow from the geometry and the
stellar surfaces (spheres, or Roche surfaces for close stars).

## Terms used in this report

Several words are used in a precise sense. They are defined here once.

- **Time of minimum (eclipse time).** The moment of the deepest light during an eclipse, measured by
  an observer from a light curve. These are **data**. They come from the literature, from the O − C
  gateway of the Czech Astronomical Society (var.astro.cz), from the O − C files of R. H. Nelson
  (AAVSO), and, for TESS, from our own measurement: for every TESS sector a mean eclipse profile is fitted
  to each eclipse with a free time.
- **Model eclipse time.** For every observed minimum the model computes the moment when the two stars of
  the pair are aligned with the observer (the conjunction), found in the N-body integration, plus the
  time the light needs to reach the observer from the pair's position on its orbit around the third
  star. It is computed at the epoch of every observed minimum and compared with it; it is **not a
  prediction** but the model's version of the same event. The geometric conjunction and the photometric
  minimum differ by seconds for these systems (about 1 s for AH Cep, from synthetic light curves), far
  below the errors of the data.
- **O − C** ("observed minus calculated"). The difference between an eclipse time and a strictly
  periodic clock (a linear ephemeris). It shows by how many minutes the eclipses come early or late.
  The model's O − C uses the same clock, so data and model can be drawn in one diagram.
- **Prediction.** A model eclipse time, inclination or light curve at a date **without data**. In this
  report predictions are only made outside the observed intervals, and they are marked as such.
- **Eclipse limit.** A pair with radii R₁, R₂ eclipses only while the projected distance of the two
  stars at conjunction is smaller than R₁ + R₂, i.e. while cos i < (R₁ + R₂)/r, where r is their
  separation at conjunction. For V907 Sco this is i > 72.8°, for QX Cas i > 70.7° (it changes slowly
  because the eccentric orbit of the pair turns), for AH Cep i > 50.4°. Above the limit the stars
  **can** eclipse; this does **not** mean that the eclipse was deep enough to be observed. Near the
  limit the eclipses are grazing and very shallow. Throughout the report we therefore separate a
  **geometric eclipse** (the inclination is above the limit) from a **detectable eclipse** (its depth is
  larger than what the photometry of the time could see: roughly 0.1 mag for photographic plates, a few
  hundredths of a magnitude for photoelectric photometry, a few tenths of a percent for TESS).
- **Third light.** Any light that does not take part in the eclipse: a third star of the system, a
  star in the same photometric aperture, or scattered light. If a fraction l₃ of the total light is
  third light, an eclipse that would remove a fraction d of the binary's light removes only (1 − l₃) d
  of the total, and is shallower in magnitudes. The depth of an eclipse therefore always depends on
  three things together: the inclination, the sizes and surface brightnesses of the two stars, and the
  third light. Whenever this report gives an eclipse depth, it states which third light was assumed.
- **Nodal period.** The time in which the orbit of the pair turns once around the total angular
  momentum of the system. During one turn the inclination of the pair swings between two extremes.

## Inclinations above 90°

A light curve measures only sin i: an orbit inclined by i and one inclined by 180° − i give the same
eclipses. Published inclinations are therefore always ≤ 90°, and the main panel of every inclination
figure shows the model in this **observable** form (0–90°). When the orbit passes through the edge-on
position (90°) during its precession, the observable inclination rises to 90° and falls back, and then
rises again once the orbit is past its turning point and turns back through edge-on. Over one passage the
observable curve therefore has two peaks at 90° with a dip between them, an "M" shape.

The N-body model also knows on which side of the sky the orbit is. The second panel shows the same
curve **unfolded** to 0–180°: there it is a single smooth wave, and the dip of the "M" turns out to be
the true maximum of the inclination above 90°. Both panels contain the same information. The
observable one is what a light curve can be compared with; the unfolded one shows the precession.

## How SWINGERsim fits data

The fit varies the masses and orbits of the N-body model until the model matches all data at once:

- times of minima. Every source (visual, photographic, photoelectric, CCD, each TESS sector) gets its
  own extra error, found by the fit, and outliers beyond 4–5σ are removed;
- radial velocities (each data set with its own zero point and extra error);
- inclinations from published light-curve solutions at their epochs. One-sided information is used
  as well: "the inclination was below 75°" or "there were no eclipses in this year" (the model must
  then stay below its own eclipse limit);
- light curves, where stated.

Without astrometry only the difference of the nodes of the two orbits on the sky is measurable, so
the node of the pair is fixed.

**Errors.** Unless stated otherwise the errors are from the curvature of χ² around the best fit
(linear errors). For V907 Sco the posterior was also sampled with emcee (Foreman-Mackey et al. 2013;
48 walkers × 1000 steps, the first 500 steps discarded as burn-in). The autocorrelation time of the
chain is τ ≈ 54 steps, so the whole chain covers about 19 τ and the 500 kept steps about 9 τ, i.e.
about nine independent samples per walker and about 430 in total. Foreman-Mackey et al. (2013) suggest
running for about ten τ, which the kept part nearly reaches. That is enough for the 16–84 % ranges
quoted here, but not for the 50 τ that the emcee online documentation (tutorial "Autocorrelation
analysis & convergence") recommends for a reliable estimate of τ itself. A chain of 400 steps gave errors that differ
from those of 1000 steps by a median of 10 % (at most a factor of 1.5, for m₃ and the phase of the
outer orbit), and medians that differ by less than 0.8σ. The quoted errors are therefore reliable to
a few tens of percent. 100 samples of this posterior were integrated over 1890–2070 to give the
ranges of the predicted inclination and of the eclipse switches.

Where the data are too poor for a unique fit (QX Cas), and for the question of a fourth body
(AH Cep), SWINGERsim builds an **ensemble**: thousands of possible systems drawn from broad ranges are
integrated, systems unstable by the criterion of Mardling & Aarseth (2001) are discarded, and the rest
are weighted by how well they reproduce the measurements. A single fit is one possible system; the
ensemble shows the range of systems the data allow.

## The systems and their data

**Table 1.** The three systems and the data used.

| | V907 Sco | QX Cas | AH Cep |
|---|---|---|---|
| eclipsing pair | 3.78 d, e ≈ 0.01, 2.74 + 2.56 M☉, late B (Zasche et al. 2023) | 6.00 d, e ≈ 0.21, ≈ 5.5 + 6.5 M☉, B1.5 V + B3 IV (Guinan et al. 2012a) | 1.775 d, circular, 16.1 + 13.7 M☉, B0.5 V (Pavlovski et al. 2018) |
| third body | 142 d, e = 0.47, 1.06 M☉ (Zasche et al. 2023) | not detected | light-time orbit of ~68 yr (Kim et al. 2005) |
| inclinations | 14 epochs 1900–2022 (Zasche et al. 2023, Table 2) | ~1960, ~1980, no eclipses from 2003 (Guinan et al. 2012a) | ~1983 (Pavlovski et al. 2018) |
| eclipses seen | 1899–1918, 1963–1986 (Lacy et al. 1999), since 2019 (TESS) | photographic minima 1933–1974; photoelectric ~1960 and ~1980; none in the mid-1990s, 2003–2008 (Guinan et al. 2012a,b) and TESS 2019–2024 | always |
| times of minima | 58 from TESS 2021–2025 | 96 from 1933–1974 (74 photographic) | 325 from 1901–2026 (216 TESS) |
| radial velocities | 172: CHIRON 2022 (Zasche et al. 2023), 1983–1994 (Lacy et al. 1999) | – | – |
| light curves | TESS 2021, 2023 | TESS 2019–2024 (no eclipses) | TESS 2019, 2024 |

TESS data (Ricker et al. 2015) are SPOC light curves unless stated otherwise.

## How to read the pictures of the systems

Each system is shown in three panels (Figs. 1, 8, 12), the published system dashed and orange, the
SWINGERsim result solid and blue:

- **From above.** The orbits seen from a point above the system, with Earth towards the bottom of
  the panel. The pair (label "pair m1+m2") and the third star move around the common centre of mass.
- **The eclipsing pair on the sky.** What a telescope with unlimited resolution would see: the
  orbits of the two stars of the pair projected on the sky, with the stars drawn to scale. An eclipse
  happens when one disc passes in front of the other. The flatter the projected orbit, the closer we
  look along the orbital plane and the deeper the eclipses.
- **Where the pair's orbit points.** The axis of the pair's orbit (a unit vector) projected on the
  sky. A point at distance sin i from the centre means inclination i; the rim (distance 1) means
  that we look exactly along the orbital plane. Eclipses are possible only while the axis lies in the
  yellow band between the rim and the dotted circle of radius sin(eclipse limit). As the orbit
  precesses, the axis moves around the total angular momentum of the system (×); labels are years,
  and the large dot marks 2025. The panel is zoomed on the track.
