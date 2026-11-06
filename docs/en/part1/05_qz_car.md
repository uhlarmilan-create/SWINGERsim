# QZ Car against Brož et al. (2022) {#qz-car}

QZ Car (HD 93206) in the Carina nebula is a 2+2 quadruple of O stars, among the most massive stars
known: the eclipsing binary Ac (Ac1 25.5 + Ac2 33.2 M☉, $P$ = 6.0 d), the eccentric binary Aa
(Aa1 69.8 + Aa2 8.85 M☉, $P$ = 20.7 d) and an outer orbit of 14 728 d (about 40 yr). Brož et al.
(2022) fitted 40 years of eclipses, spectra and interferometry with their N-body code Xitau and
published a complete model. It is the most massive system in our tests and a strict test of
SWINGERsim's geometry and photometry, because the orbit sizes and the light shares of the four stars
are known precisely. Their Table 2 is the input of SWINGERsim's QZ Car configurations; radii and
effective temperatures are taken from their Tables 2 and 3. In the configuration the four stars are
m1 = Ac1, m2 = Ac2, m3 = Aa1 and m4 = Aa2.

## Geometry, light and the evolutionary masses {#qzcar-benchmarks}

**Table 28.** QZ Car benchmarks. The first three rows test SWINGERsim's hierarchical initialisation
(semi-major axes from periods and masses), the next four its photometry (V-band light shares from radii
and temperatures), the last three the MIST tracks (the mass of a single-star model with the observed
radius and temperature). Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `qzcar-broz-a1` | own | 54.032 R☉ | 53.9946 R☉ | −0.07 % | 0.2 % | PASS |
| `qzcar-broz-a2` | own | 136.132 R☉ | 136.081 R☉ | −0.04 % | 0.2 % | PASS |
| `qzcar-broz-a3` | own | 13 053 R☉ | 13 046 R☉ | −0.05 % | 0.2 % | PASS |
| `qzcar-broz-lv-ac1` | own | 0.22 | 0.220775 | +0.000775 | 0.005 | PASS |
| `qzcar-broz-lv-ac2` | own | 0.133 | 0.133856 | +0.000856 | 0.005 | PASS |
| `qzcar-broz-lv-aa1` | own | 0.635 | 0.635024 | +2.4·10⁻⁵ | 0.005 | PASS |
| `qzcar-broz-lv-aa2` | own | 0.01 | 0.0103458 | +0.000346 | 0.005 | PASS |
| `qzcar-evolutionary-mass-ac1` | own/tab | 25.5 M☉ | 27.0258 M☉ | +5.98 % | 15 % | PASS |
| `qzcar-evolutionary-mass-ac2` | own/tab | 33.2 M☉ | 26.0243 M☉ | −21.61 % | 15 % | **EXPLAINED** |
| `qzcar-evolutionary-mass-aa1` | own/tab | 69.8 M☉ | 42.5172 M☉ | −39.09 % | 15 % | **EXPLAINED** |

Table 28 lists the ten benchmarks; the paragraphs below read it group by group.

**Geometry (level B).** The three semi-major axes follow from the published periods and masses to
0.04–0.07 %; the tolerance of 0.2 % covers the rounding of the masses in Table 2 to three digits.

**Light (level C).** Black bodies with the Claret/TLUSTY limb darkening and the published radii and
temperatures give the shares of the four stars in the V-band light – 0.221, 0.134, 0.635 and 0.010 – to
0.001 of the values that Brož et al. derived from synthetic spectra. The tolerance of 0.005 is the
difference between a black body and a model atmosphere at these temperatures.

**Evolutionary masses (level D).** The **evolutionary mass** of a star is the mass of the MIST track
that reaches its observed radius and temperature; the **dynamical mass** comes from the orbits. For a
normal single star they agree. Ac1 is a 27.0 M☉ star of 5.6 Myr on the MIST tracks, within 6 % of its
dynamical mass. Ac2 (26.0 M☉, 5.1 Myr) and above all Aa1 (42.5 M☉, 3.9 Myr) do not fit their dynamical
masses of 33.2 and 69.8 M☉. These two benchmarks are EXPLAINED: both stars are likely mass gainers. The
nominal model of Brož et al. (2022, Table 2) fits the temperature of Ac2 (32 979 K) and leaves its
rotation free; the fit gives 410 km/s, the signature of a star that has accreted angular momentum in the Ac1 → Ac2
transfer that Brož et al. describe. For Aa1, Brož et al. note that it is extremely massive compared
with Aa2 and suggest past transfer Aa2 → Aa1. A gainer is rejuvenated (fresh hydrogen mixed into its
core) and looks younger and less massive than a single star of the same mass; single-star tracks
cannot describe it. Placed at their dynamical masses, the MIST tracks reach the observed radii at
6.0 Myr (Ac1), 3.5 Myr (Ac2) and 2.5 Myr (Aa1), and the rotating tracks never give Aa1 its observed
28 R☉ at 70 M☉: single-star tracks cannot make the system coeval. SWINGERsim therefore lets a star
start with an effective age (`stars.<name>.age_offset_myr`) so that every component begins with its
observed radius, and it rejuvenates gainers during a run ([Chapter 3](#mass-transfer)).

## Synthetic minima {#qzcar-minima}

**Input and method.** The minima of the eclipsing pair Ac were computed over one outer period
(40.5 yr) with the spherical model from a V-band light curve with 10-min cadence and yearly observing
seasons (4933 minima), and over one year with the PHOEBE hybrid, which deforms the stars into their
Roche shapes (4000 triangles per star).

![**Figure 7.** O−C diagram of the synthetic minima of the eclipsing pair Ac of QZ Car over one outer period of about 40 years. Horizontal axis: time; vertical axis: O−C in days, the time of each minimum minus the time predicted with a constant period. Blue: primary minima; orange: secondary minima. The wave is the light-travel-time effect of the outer orbit.](../../figures/showcase/qzcar_oc.png)

**Result.** Figure 7 shows the O−C of the 4933 minima. If Ac were alone, all points would lie on a
horizontal line. Instead they follow a wave of about ±0.2 d (an O−C range of 0.38 d): while Ac moves
away from us on the outer orbit its minima come late, while it approaches they come early. The wave is
not a sinusoid because the outer orbit is eccentric. Its amplitude equals $a_{12}\sin i/c$ of the
published outer orbit. Primary and secondary minima lie on the same curve, so both are timed equally
well. Measured against the barycentric conjunctions, the primary minima come 16 s late and the
secondary minima 10 s early: the Rømer delay inside Ac is ±16 s (Kaplan 2010), and the slightly
eccentric orbit of Ac adds the rest. The PHOEBE minima over one year agree with the spherical-model
minima to 1 s for the secondary and lie within +6 s and −5 s of the conjunctions over three weeks. The
timing of a massive, tidally distorted eclipsing binary is thus reproduced to a few seconds by two
independent light-curve models.

## Present-day dynamics {#qzcar-dynamics}

**Input.** An N-body ensemble of 2000 members with masses and orbital elements drawn within the
uncertainties of Brož et al. (2022) and the unknown nodes of the binaries drawn isotropically,
integrated 100 yr forward and 100 yr backward, with the condition that Ac has eclipsed since 1955.

**Result.** All 2000 members stay intact in both directions, and all are stable by the Mardling–Aarseth
and Vynatheya et al. (2022) criteria. Over a century the orbits hardly change: the inclination of Ac
drifts by a median of 9·10⁻⁴ °/yr (16–84 % range −0.0011 to +0.0010 °/yr), about 0.04° in 40 yr, which
Brož et al. could not have measured. Ac therefore eclipses for centuries; Aa eclipses in 24 % of the
isotropic orientations of its node. QZ Car is a stable laboratory for timing studies.

## Evolution to the first Roche-lobe overflow {#qzcar-evolution}

**Question.** The stars of QZ Car are young (2–6 Myr) but massive stars evolve fast. Which star fills
its Roche lobe first, and when?

**Input.** Three ensembles of 1000 members each with the evolving backend and MIST tracks, masses drawn
within the published errors: (a) all stars coeval at 2.5 Myr on non-rotating tracks; (b) the same with
rotating tracks ($v/v_\mathrm{crit}$ = 0.4), which have stronger winds; (c) effective ages chosen so that
every star starts with its observed radius (the stars are then not coeval).

![**Figure 8.** QZ Car: time from today to the first Roche-lobe overflow (horizontal axis, Myr) in three ensembles of 1000 members (vertical axis: number of members per 0.05-Myr bin). Colours show which star overflows first: blue Aa1 (the 69.8 M☉ supergiant), green Aa2, orange Ac1. The legend gives the number of members per star and, at its top, the members that end otherwise. (a) Coeval stars at 2.5 Myr, non-rotating MIST tracks. (b) The same with rotating tracks. (c) Effective ages that give every star its observed radius.](../../figures/part1/qzcar_overflow.png)

**Result.** Figure 8 shows the three ensembles. (a) With coeval non-rotating stars, all 1000 members
end with a Roche-lobe overflow of Aa1, after a median of 0.44 Myr (16–84 %: 0.30–0.53 Myr); the outer
orbit widens slightly by the winds (median from 60.5 to 62.1 au). (b) With rotating tracks the stronger
winds widen the orbits: only 338 members overflow (Aa1 first in 244, median 0.92 Myr; Aa2 first in 94,
median 1.10 Myr), while 616 reach the end of their tracks without overflow, with Aa widened to about
2.3 au and the outer orbit to about 105 au; 46 members are refused, because an Aa1 heavier than about
73.3 M☉ does not exist at 2.5 Myr on the rotating grid. (c) With the observed radii, 803 members end in
an overflow – Aa1 first in 626 (median 0.38 Myr) and Ac1 first in 177 (median 0.12 Myr) – and 197 are
refused at the start, because an Ac1 heavier than 26.4 M☉ at its effective age already overfills its
Roche lobe in the observed 6-d orbit. This last number is a result in itself: the observed radius of
Ac1 and the fact that Ac is detached bound the mass of Ac1 from above, a constraint the dynamical fit
does not contain. In panel (c) the two binaries race: Ac1 wins only when the overflow comes very soon,
in members in which Ac1 is drawn close to that upper bound.

There is no published simulation of the future of QZ Car to compare with; the evolution is a
prediction. Its ingredients are validated separately: the MIST tracks of massive stars against Geneva
models and leave-one-out interpolation, the winds against Vink et al. (2001) ([Chapter 3](#massive-giants)),
the widening of orbits by winds against Jeans' law (Chapters 3 and [7](#catalogue)). The spread
between the three ensembles – a factor of two in the time to the overflow – is the uncertainty of the
stellar models of very massive stars, not of the dynamics.
