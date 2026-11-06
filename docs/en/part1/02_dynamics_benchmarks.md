# Dynamics: orbits, secular evolution and stability {#dynamics}

This chapter tests how SWINGERsim moves the stars: the N-body backend with relativity, the secular
backends, the detection of instability, and the reproduction of four classic studies of
hierarchical systems. It also tests hierarchies of more than four stars, scattering experiments and
the influence of the environment (gas, the Galactic tide, passing stars, gravitational waves).

Two effects recur in the whole report. **Apsidal precession** is the slow turning of an elliptical
orbit in its own plane (the argument of periastron $\omega$ grows); its causes are the other stars,
general relativity and the tidal and rotational bulges of the stars. **Nodal precession** is the
turning of the orbital plane around the total angular momentum of the system; it changes the
inclination to the sky and so makes eclipses appear and disappear. When the mutual inclination of an
inner and an outer orbit is large (roughly 40°–140°), the outer body exchanges angular momentum with
the inner orbit and pumps its eccentricity up and down periodically: the **von Zeipel–Lidov–Kozai
(ZLK) cycles** (Kozai 1962; Lidov 1962; Naoz 2016). The characteristic time of these cycles is the
**Kozai time** $t_K$, of order $P_\mathrm{out}^2/P_\mathrm{in}$ times a mass ratio.

## Relativity, ZLK cycles and nodal precession in the N-body backend {#nbody-basics}

**Reference and input.** Three analytic limits test the N-body backend. (i) The first post-Newtonian
(1PN) perihelion advance of Mercury, 42.98″ per century for the real orbit ($a$ = 0.387 au,
$e$ = 0.206), one of the classic tests of general relativity (Will 2014). (ii) The ZLK cycles of a
test particle at 1 au around a 1 M☉ star, perturbed by a 1 M☉ star on a circular orbit at 20 au. In
this quadrupole, test-particle limit the theory gives $e_\mathrm{max} = \sqrt{1 - \tfrac53\cos^2 i_0}$
for an initial mutual inclination $i_0$, a mutual inclination of 39.23° at the eccentricity maximum,
no eccentricity growth below 39.2°, and an exactly conserved quantity $\sqrt{1-e^2}\cos i$ (Kozai
1962; Lidov 1962; Naoz 2016). (iii) The nodal precession period of the same set-up at $i$ = 20°
from the quadrupole formula of Söderhjelm (1975). The secular backends are tested in the same limits,
and in addition on the relativistic apsidal period of an isolated binary ($a$ = 0.05 au, $e$ = 0.3).

**Table 3.** Relativity, ZLK cycles and precession: N-body backend (class R) and secular backends
(class own). Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `gr-mercury` | R | 42.98″/century | 42.9803″/century | +0.00 % | 0.5 % | PASS |
| `zlk-emax-65` | R | 0.838047 | 0.84023 | +0.26 % | 1 % | PASS |
| `zlk-emax-80` | R | 0.974548 | 0.976155 | +0.16 % | 1 % | PASS |
| `zlk-imin-65` | R | 39.2315° | 39.3441° | +0.113° | 0.5° | PASS |
| `zlk-below-critical` | R | 0 | 0.00173037 | +0.00173 | 0.01 | PASS |
| `zlk-kozai-constant` | R | 0 | 5.16811·10⁻⁵ | +5.17·10⁻⁵ | 0.001 | PASS |
| `nodal-precession` | R | 11351.4 yr | 11339.8 yr | −0.10 % | 1 % | PASS |
| `zlk-period-nbody` | own | 0 | 0.00660493 | +0.0066 | 0.02 | PASS |
| `quad-tp-emax` | own | 0.838047 | 0.838082 | +0.00 % | 0.1 % | PASS |
| `secular-kozai-emax` | own | 0.838047 | 0.838082 | +0.00 % | 0.5 % | PASS |
| `secular-3p1-kozai-emax` | own | 0.838047 | 0.838082 | +0.00 % | 0.5 % | PASS |
| `secular-nodal-precession` | own | 11351.4 yr | 11351.4 yr | −0.00 % | 0.5 % | PASS |
| `secular-gr-apsidal` | own | 7113.77 yr | 7113.77 yr | −0.00 % | 0.5 % | PASS |
| `brouwer-law` | R | 0.5 | 0.59301 | +0.093 | 0.15 | PASS |

**Result.** Table 3 lists the comparisons. The N-body backend gives Mercury's advance to better than
10⁻⁴; with relativity switched off it gives 0″, so the benchmark detects missing physics. The ZLK
maxima of the N-body runs lie 0.16–0.26 % above the formula. This small excess is expected and has
the right sign: the formula neglects terms of order $a_\mathrm{in}/a_\mathrm{out}$ (here 1/20), and
the N-body orbit starts at $e$ = 0.01 rather than exactly zero. The quantity $\sqrt{1-e^2}\cos i$,
averaged over the first and the last quarter of a 20 000-yr run, drifts by 5·10⁻⁵ while $e$ swings
between 0 and 0.84; its instantaneous value fluctuates by about 0.5 % because the elements are
osculating (computed from the instantaneous positions and velocities). Below the critical
inclination the eccentricity stays below 0.002. The ZLK period of the N-body run is 0.66 % longer
than the quadrupole secular period. The nodal precession period agrees with Söderhjelm's formula to
0.10 %.

The three secular backends agree with the same analytic limits to better than 10⁻⁴ (columns
"Difference" of 0.00 %). This is expected and is a strong test of their implementation: they solve
exactly the averaged equations of these limits, without the short-period terms that make the N-body
values deviate slightly. The test-particle $e_\mathrm{max}$ is reproduced by the triple code, by the
3+1 branch of the nested code and by the 2+2 code with the perturber replaced by a binary of
negligible size.

**Numerical quality of the integrator.** The energy error of a high-order integrator with
round-off errors grows like a random walk, $\propto t^{1/2}$ (Brouwer's law; Rein & Spiegel 2015).
Six realisations over 3·10⁴ orbits give a growth exponent of 0.59 (`brouwer-law`), within the
tolerance of 0.15 set by the scatter of six realisations. This is a test of IAS15 as SWINGERsim uses
it, not of SWINGERsim's physics.

## Stability: hierarchy changes and chaos indicators {#stability-basics}

A hierarchical configuration is **dynamically stable** if it keeps its hierarchy (who orbits whom)
for a very long time. The **Mardling–Aarseth criterion** (MA01; Mardling & Aarseth 2001) gives the
smallest ratio of outer periastron to inner semi-major axis for which a triple is stable. **MEGNO**
(Mean Exponential Growth factor of Nearby Orbits; Cincotta & Simó 2000) measures how fast
neighbouring trajectories diverge: it tends to 2 for regular (quasi-periodic) motion and grows
without bound for chaotic motion. SWINGERsim computes MEGNO with REBOUND's variational equations.

**Table 4.** Detection of instability and chaos (class R).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `stability-unstable` | R | 1 | 1 | +0 | 0 | PASS |
| `stability-stable` | R | 0 | 0 | +0 | 0 | PASS |
| `megno-regular` | R | 2 | 1.99959 | −0.00041 | 0.05 | PASS |
| `megno-chaotic` | R | 1 | 1 | +0 | 0 | PASS |

An equal-mass coplanar 2+2 is placed at half and at 2.5 times the MA01 limit and integrated for 100
outer periods. The first must change its hierarchy (1 = yes) and does, ejecting a star after about two
outer orbits; the second must not (0 = no change) and does not. The margin of a factor of two on
either side is needed because MA01 was derived for triples, not for quadruples. The MEGNO of a
regular hierarchical triple converges to 1.998 (expected 2), and a 2+2 inside the unstable region
reaches MEGNO > 4 (1 = yes). The stability tools as a whole are tested against a published neural
classifier in [Chapter 9](#stability).

## Reproduction of Naoz et al. (2013) {#naoz2013}

Naoz et al. (2013) showed that octupole-order terms – the next order beyond the quadrupole in the
expansion in $a_\mathrm{in}/a_\mathrm{out}$, important when the outer orbit is eccentric – can flip an
inner orbit from prograde to retrograde. Two of their figures serve as level-B benchmarks, with the
input values taken from the captions and the text.

**Input.** *Fig. 3* (eccentric ZLK flip): a 1 M☉ star with a Jupiter-mass planet at 6 au and a brown
dwarf of 40 Jupiter masses at 100 au with $e$ = 0.6, initial mutual inclination 65°. *Fig. 4* (the
PSR B1620−26 set-up): a massive inner binary flipping at quadrupole order. The text gives the range of
the mutual inclination as 57.5°–106.7°. The caption states initial inclinations of the two orbits to
the total angular momentum of 6.75° and 63.25°; this corresponds to an angular-momentum ratio
$G_1/G_2 = \sin 63.25°/\sin 6.75° = 7.60$, whereas the quoted masses and $a_2$ = 50 au give 6.79. An
independent quadrupole calculation reproduces the published range with 7.60 (57.53°–106.70°) and not
with 6.79 (57.3°–105.0°). The benchmark therefore uses $a_2$ = 39.9 au, which gives the published
ratio.

**Table 5.** Naoz et al. (2013), Figs. 3 and 4: N-body backend (class R) and secular backends (own).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `naoz13-fig4-imin` | R | 57.5° | 57.2391° | −0.261° | 0.5° | PASS |
| `naoz13-fig4-imax` | R | 106.7° | 106.93° | +0.23° | 0.5° | PASS |
| `quad-naoz13-fig4-imin` | own | 57.5° | 57.5376° | +0.0376° | 0.1° | PASS |
| `quad-naoz13-fig4-imax` | own | 106.7° | 106.721° | +0.0207° | 0.1° | PASS |
| `naoz13-fig3-flip-time` | R | 3.8 Myr | 3.9868 Myr | +4.92 % | 10 % | PASS |
| `naoz13-fig3-imax` | R | 142.6° | 143.01° | +0.408° | 2° | PASS |
| `secular-naoz13-fig3-flip-time` | own | 3.8 Myr | 4.04981 Myr | +6.57 % | 10 % | PASS |
| `secular-naoz13-fig3-imax` | own | 142.6° | 142.556° | −0.0435° | 2° | PASS |

![**Figure 2.** SWINGERsim reproducing Naoz et al. (2013). Both panels show the mutual inclination of the inner and outer orbit (vertical axis, degrees) against time; the horizontal black line marks 90°, the boundary between prograde and retrograde. *Left:* their Fig. 3, the eccentric ZLK flip, N-body run over 5 Myr (blue line). The dashed red vertical line is the first flip read from the published figure (3.8 Myr), the dotted red horizontal line the published maximum after the flip (142.6°). *Right:* their Fig. 4, the PSR B1620−26 set-up over 500 kyr; blue with the angular-momentum ratio of the paper ($G_1/G_2$ = 7.60), grey with the quoted $a_2$ = 50 au ($G_1/G_2$ = 6.79); the dotted red lines mark the published range 57.5°–106.7°.](../../figures/part1/naoz2013.png)

**Result.** Figure 2 shows both set-ups. In the left panel, look at the ZLK minima before the flip:
they lie at 38°, and their envelope rises slowly towards 90°, as in the published figure. The orbit
first crosses 90° at 3.99 Myr in the N-body run and 4.05 Myr in the secular backend, against 3.8 Myr
read from the published figure (+4.9 % and +6.6 %; the tolerance of 10 % reflects the reading of a
time from a figure and the sensitivity of the flip time to the phase of the cycles). After the flip
the inclination reaches 143.0° (N-body) and 142.6° (secular), against 142.6°. In the right panel, the
blue curve oscillates exactly between the two dotted lines: the secular quadrupole code reproduces
the published range to 0.04° (57.54°–106.72°), the N-body backend to 0.26° (57.24°–106.93°), the
difference being the short-period terms that the averaged equations omit. The grey curve shows that
the quoted $a_2$ = 50 au would give a narrower range; the published curve belongs to the ratio in
the caption.

## Reproduction of Vokrouhlický (2016) {#vokrouhlicky2016}

Vokrouhlický (2016) developed the secular theory of a **2+2 quadruple** – two binaries A and B on a
common outer orbit – and applied it to a V994 Her-like system. His
quantity of interest is $J_{AB}$, the angle between the orbital planes of A and B: it controls whether
both binaries can eclipse at the same time. He derived an analytic upper limit of $J_{AB}$ and
illustrated four set-ups (his Figs. 2, 3, 4 and 7).

**Input.** The V994 Her-like system of the paper, set up from sky-frame elements: Fig. 2 with circular
A and B and nodes $\Omega_A$ = 10°, $\Omega_B$ = 15°; Fig. 3 with the observed $e_A$ = 0.031,
$e_B$ = 0.125; Fig. 4 with $\Omega_A$ = 35°, $\Omega_B$ = 40°; Fig. 7 as Fig. 4 plus the tidal
quadrupole distortion of the stars. For Fig. 7 the paper tunes one deformability per binary to the
observed apsidal periods of about 116 yr (A) and 111 yr (B), with radii $R\propto M^{0.8}$ (2.4, 2.1,
1.6 and 1.7 R☉). SWINGERsim uses the same radii, tides without rotation and dissipation, and apsidal
constants $k_2$ = 0.00754 (A) and 0.00333 (B), for which Sterne's formula gives exactly these apsidal
periods for the isolated binaries. The four set-ups were integrated with the N-body backend (25–60
minutes each) and with the 2+2 secular backend (seconds each).

**Table 6.** Vokrouhlický (2016): angles and eclipse probabilities (all class own). Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `vok16-v994-JA` | own | 11.6° | 11.6462° | +0.0462° | 0.1° | PASS |
| `vok16-v994-JB` | own | 15.8° | 15.7932° | −0.00678° | 0.1° | PASS |
| `vok16-jab-max-analytic` | own | 5.2° | 5.20074° | +0.000736° | 0.05° | PASS |
| `vok16-fig2-jab-min` | own | 3.5° | 3.47152° | −0.0285° | 0.15° | PASS |
| `vok16-fig2-jab-max` | own | 5.2° | 5.20069° | +0.000687° | 0.15° | PASS |
| `vok16-fig2-jab-1myr` | own | 5.2° | 5.2006° | +0.000633° | 0.15° | PASS |
| `vok16-fig3-jab-max` | own | 6° | 6.20874° | +0.209° | 0.5° | PASS |
| `vok16-eclipse-prob-A` | own | 0.366 | 0.369265 | +0.89 % | 2 % | PASS |
| `vok16-eclipse-prob-B` | own | 0.397 | 0.399075 | +0.52 % | 2 % | PASS |

The first two rows of Table 6 check SWINGERsim's conversion of sky-frame elements into the tilts $J_A$, $J_B$
of the binaries to the outer orbit – an independent check of the angle conventions. The analytic
limit of $J_{AB}$ is 5.201° (paper 5.2°). The secular backend reproduces the ranges of $J_{AB}$
quoted in the text for Figs. 2 and 3 to 0.03° and 0.21°, and over 1 Myr the maximum stays at the
analytic limit (the paper finds the same over 100 Myr). SWINGERsim's eclipse criterion applied to
isotropically oriented orbits gives the chances that A and B eclipse, 0.369 and 0.399, against the
published 0.366 and 0.397.

**Table 7.** Vokrouhlický (2016), Figs. 2, 3, 4 and 7: published behaviour against the N-body
backend (class R) and the 2+2 secular backend (class own). Each row is one set-up; the second
column gives what the paper shows, the third and fourth the ranges SWINGERsim computes.

| Set-up | Paper | SWINGERsim N-body | SWINGERsim secular |
|---|---|---|---|
| Fig. 2: circular A and B | $J_{AB}$ 3.5°–5.2°, $e_A$, $e_B$ near zero | $J_{AB}$ 3.41°–5.19° over 50 kyr, $e$ < 0.005 | $J_{AB}$ 3.47°–5.20°, $e_A$ < 0.0051 |
| Fig. 3: observed $e_A$, $e_B$ | $J_{AB}$ ≈ 1°–6°, small oscillations of $e$ | $J_{AB}$ 0.94°–6.30°, $e_A$ 0.029–0.036, $e_B$ 0.123–0.140 | $J_{AB}$ 0.78°–6.21°, $e_A$ 0.029–0.036, $e_B$ 0.123–0.138 |
| Fig. 4: $\Omega_A$ = 35°, $\Omega_B$ = 40° | $J_A$ ≈ 35°, $J_B$ ≈ 40°, $e_B$ above 0.6, $J_{AB}$ from 0° to more than 70° | $J_A$ = 35.4°, $J_B$ = 40.3°, $e_B$ up to 0.631, $J_{AB}$ 3.5°–73.4° | $e_B$ up to 0.610, $J_{AB}$ 3.2°–71.1° |
| Fig. 7: Fig. 4 + tidal distortion | order restored: small oscillations of $e$ and $i$ | $e_B$ 0.124–0.136, $e_A$ 0.029–0.032, $J_{AB}$ 0.4°–5.8°, apsidal periods 113 and 108 yr | $e_B$ 0.123–0.134, $J_{AB}$ 0.29°–5.91°, apsidal periods 113.1 and 108.4 yr |

![**Figure 3.** SWINGERsim reproducing Vokrouhlický (2016), Figs. 2, 3, 4 and 7 (one column per set-up). *Top row:* the angle $J_{AB}$ between the orbital planes of the two binaries (degrees) against time (kyr); black: N-body backend; red dashed: 2+2 secular backend; grey dashed: the analytic limit 5.2° of the paper. *Bottom row:* the eccentricities of A (blue) and B (orange; red dashed for the secular run) against time.](../../figures/vokrouhlicky2016.png)

**Result.** Figure 3 and Table 7 compare the four set-ups. In the first column, $J_{AB}$ oscillates
between 3.4° and 5.2° and never crosses the grey analytic limit; the black and red curves overlap.
In the second column the observed eccentricities let $J_{AB}$ swing between about 1° and 6°, as in
the paper. The third column is the interesting one: with larger tilts, B enters ZLK cycles and
reaches $e_B$ = 0.63, $J_{AB}$ wanders up to 73°, and the N-body and secular curves agree for about
4 kyr and then separate. The paper reports exactly this for its own two integrators ("they start to
significantly deviate after ≃5 kyr"): the motion is chaotic, and two correct integrations of it
diverge. In this set-up the closest periastron of B reaches 1.06 times the sum of the stellar radii,
which supports the paper's conclusion that the point-mass model is incomplete here. The fourth
column adds the tidal bulges: the fast apsidal precession they cause (113 and 108 yr) detunes the
ZLK cycles of B, its eccentricity stays at 0.124–0.136 instead of reaching 0.63, and the planes of A
and B stay within 6° of each other – the order that the paper restores in its Fig. 7. Energy is
conserved to 2.5·10⁻¹² in the N-body runs. Quantitatively, the published ranges are reproduced to
0.1–0.3° wherever the motion is regular.

## Reproduction of Pejcha et al. (2013): an ensemble {#pejcha2013}

Pejcha et al. (2013) integrated thousands of equal-mass 2+2 quadruples with an N-body code and
counted how often a binary reaches an extreme eccentricity, $1-e<10^{-3}$, within 300 Kozai times of
binary A. The inclinations of A and B to the outer orbit were not drawn isotropically but set on a
grid in $(\cos i_A, \cos i_B)$, with 100 random set-ups of the other angles at every grid point; SWINGERsim draws them uniformly in
$\cos i$ (isotropically), the continuous counterpart of a uniform grid. For their set F6 ($a_A$ = 0.573 au, $e_A$ = 0.3, $a_B$ = 0.351 au, $e_B$ = 0.1,
$a_{AB}$ = 6 au, $e_{AB}$ = 0.3, $t_K$ = 64.9 yr) they found 35.8 % for A and 53.0 % for B – 4.8 and 10
times more than for the equivalent triples. This "enhancement" in 2+2 quadruples is the main result
of the paper. Vokrouhlický (2016, Sect. 4) repeated the experiment with his secular model
($a_{AB}$ = 5 au, 11 kyr) and found 31.0 % and 48.1 %.

An **ensemble** in SWINGERsim is a set of runs (**members**) whose unknown parameters are drawn from a
**prior** (a distribution of plausible values); the result is the fraction of members with a given
outcome, with a 68 % Wilson interval as its statistical error.

**Table 8.** Pejcha et al. (2013) and Vokrouhlický (2016): fractions of binaries reaching
$1-e<10^{-3}$. The N-body benchmarks are manual (class R), the secular ones slow (class own). Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `pejcha13-F6-fA` | R | 0.358 | 0.24 | −0.118 | 0.14 | PASS |
| `pejcha13-F6-fB` | R | 0.530 | 0.41 | −0.12 | 0.15 | PASS |
| `vok16-pejcha-fA` | own | 0.31 | 0.21 | −0.10 | 0.07 | **EXPLAINED** |
| `vok16-pejcha-fB` | own | 0.481 | 0.4025 | −0.0785 | 0.07 | **EXPLAINED** |

**Table 9.** The F6 ensembles in detail. Column 2 is the number of members, columns 3–4 the fractions
of A and B reaching $1-e<10^{-3}$ with their 1σ statistical errors, column 5 the published values.

| Set-up | Members | Fraction A | Fraction B | Published |
|---|---|---|---|---|
| F6, N-body, 300 $t_K$ | 100 | 0.24 ± 0.04 | 0.41 ± 0.05 | 0.358, 0.530 (N-body) |
| F6, secular, 300 $t_K$ | 400 | 0.215 ± 0.021 | 0.358 ± 0.024 | – |
| F6, secular, 900 $t_K$ | 200 | 0.325 ± 0.033 | 0.505 ± 0.035 | – |
| Vokrouhlický set-up, secular, 11 kyr (≈ 300 $t_K$) | 400 | 0.21 | 0.40 | 0.310, 0.481 (secular) |
| the same, secular, 50 kyr | 400 | 0.36 | 0.57 | 0.310, 0.481 (secular) |
| F6 with B replaced by a point mass (triple), secular | 400 | 0.020 ± 0.007 | – | about 0.07 (N-body), 0.03 analytic |
| the same triple, N-body | 100 | 0.06 ± 0.02 | – | about 0.07 |

**What agrees** (Tables 8 and 9; in Table 9 each row is one ensemble, and the fractions are given with their 1σ errors where the ensemble was large enough to state them). (i) *The triple.* With B replaced by a point mass, the N-body backend gives
6 ± 2 % of members reaching $1-e<10^{-3}$, as Pejcha et al. (about 7 %). The secular backend gives
2 % (the analytic quadrupole estimate is about 3 %): in a triple of this mild hierarchy the
short-period terms raise the eccentricity peaks (Luo, Katz & Dong 2016). (ii) *The enhancement.* In
N-body, binary A of the quadruple reaches the limit 4.0 times more often than in the triple (paper:
4.8), B more often than A, and no member gets there from outside the ZLK window – the main result of
both papers. (iii) *The two SWINGERsim backends agree with each other* within the statistical errors:
for the same 100 members 0.24 vs 0.24 (A) and 0.41 vs 0.39 (B), although individual chaotic
members end differently.

**What differs, and why.** At the nominal 300 Kozai times SWINGERsim finds about two thirds of the
published quadruple fractions, with both backends and in both set-ups, while the triple agrees. Four
tests narrow the cause down. An independent secular code, SecularMultiple (Hamers & Portegies Zwart
2016), on the same 200 F6 members gives 0.19 (A) and 0.37 (B) at the same multipole orders, and
0.225 and 0.385 with all orders, against SWINGERsim's 0.23 and 0.37 for these members. SWINGERsim's
secular equations integrated with a fixed 0.1-yr fourth-order Runge–Kutta step, as in Vokrouhlický
(2016), give the same fractions as the adaptive integrator (0.275 and 0.30 for 40 members). And the
fractions keep growing with the length of the run (Figure 4).

![**Figure 4.** Set F6 of Pejcha et al. (2013): the fraction of members that have reached $1-e<10^{-3}$ at least once (vertical axis) against the time in Kozai times of binary A (horizontal axis, logarithmic); SWINGERsim secular backend, 200 members. Blue: binary A; orange: binary B. The dotted horizontal lines are the published fractions; the dashed vertical line marks 300 $t_K$, the stated length of the published runs.](../../figures/part1/pejcha_time_curve.png)

Figure 4 shows the cumulative fraction against time. The curves only rise, because a member that has
reached the limit stays counted. They do not saturate: they grow roughly logarithmically, by 0.2 (A)
and 0.3 (B) per factor of ten in time, and reach the published values at about 1500 $t_K$ (A) and
1200 $t_K$ (B). Fang, Thompson & Hirata (2018) found the same non-saturation in wider 2+2 systems
(about 10 % per factor of ten). With such growth, the published numbers depend strongly on the exact
length of the runs. Running the set of Vokrouhlický (2016) to 50 kyr gives 0.36 and 0.57: the
published 0.31 and 0.48 lie between SWINGERsim at 11 kyr and at 50 kyr (Figure 6, bottom right, in
[Chapter 4](#discs)). Three independent integrations (SWINGERsim secular, SWINGERsim N-body with
REBOUND, SecularMultiple) agree with each other, and the integration method of the secular paper is
ruled out as the cause. The difference therefore lies in the effective length of the published runs
(the published fractions correspond to 4–5 times the stated length) and possibly in the sampling of
their inclination grid, which Pejcha et al. note "might not be distributed uniformly". The N-body
benchmark passes at three standard deviations; the secular benchmark is EXPLAINED by the length of
the published runs. The morphology of the high-eccentricity region in the plane of
$(\cos J_A, \cos J_B)$ agrees with Fig. 11 of Vokrouhlický (2016): a diagonal band for A and a broad
band $|\cos J_B|\lesssim 0.6$ for B.

## Comparison with MSE {#mse}

**Reference.** MSE (Multiple Stellar Evolution; Hamers et al. 2021) is an independent code for
hierarchical multiples. It was built from its public repository and run on the same 3+1 and 2+2
quadruples as SWINGERsim's secular backends, with the same Hamiltonian (MSE's binary–triplet,
hexadecupole and dotriacontupole terms switched off). A third comparison is a *fate*: a 3+1 whose
inner binary (0.3 au) is driven by ZLK cycles into Roche-lobe overflow – MSE's event "RLOF start"
against SWINGERsim's first Roche-lobe filling at periastron, with the same radii.

**Table 10.** SWINGERsim's secular backends against MSE (level C, class own). Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `mse-3p1-inner-emax` | own | 0.9707 | 0.970714 | +1.37·10⁻⁵ | 0.002 | PASS |
| `mse-3p1-inner-t-emax` | own | 49 325 yr | 49 322.8 yr | −0.00 % | 0.2 % | PASS |
| `mse-3p1-mid-emax` | own | 0.2245 | 0.224461 | −3.86·10⁻⁵ | 0.002 | PASS |
| `mse-2p2-A-emax` | own | 0.9494 | 0.949492 | +9.2·10⁻⁵ | 0.002 | PASS |
| `mse-2p2-B-emax` | own | 0.7758 | 0.776227 | +0.000427 | 0.002 | PASS |
| `mse-3p1-rlof-time` | own | 380.03 yr | 380.1 yr | +0.02 % | 0.2 % | PASS |

**Result (Table 10).** The eccentricity maxima agree to 10⁻⁴ or better and the time of the maximum of the 3+1
inner binary to 2 yr in 49 kyr. For the 2+2, the first three ZLK cycles of A and all of B agree, after
which the two runs drift apart in phase, as expected in a chaotic regime. MSE's extra multipole terms
change $e_\mathrm{max}$ by less than 0.001 in these systems. The onset of Roche-lobe overflow agrees to
0.07 yr in 380 yr. Two practical differences are worth knowing: MSE switches to direct N-body
integration when a system becomes dynamically unstable, whereas SWINGERsim's secular backends only
report the stability diagnostics (the user then switches to the N-body backend); and MSE loses the
argument of periastron of an orbit with an inclination of exactly 0°, so the comparison uses 0.5–2°.

## Scattering and the statistical escape theory {#scattering}

SWINGERsim can set up encounters of binaries and single stars (scattering experiments) and integrate
them with REBOUND until the outcome is clear. Antognini & Thompson (2016, Table 3) give the
cross-sections of binary–binary encounters (equal masses, $a_1 = a_2$, circular orbits, velocity at
infinity equal to the critical velocity $v_c$ at which the total energy is zero, impact parameters up
to 7$a$). Stone & Leigh (2019) derived the statistics of the break-up of a resonant
("democratic") three-body system from ergodic arguments; SWINGERsim implements their distributions.

**Table 11.** Scattering experiments (class R) and the statistical escape theory (own). Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `scatter-hyperbola-deflection` | own | 173.548° | 173.548° | +1.4·10⁻¹²° | 10⁻⁶° | PASS |
| `at16-binbin-single-ionization` | R | 2.645 | 2.6873 | +1.60 % | 10 % | PASS |
| `at16-binbin-exchange-ionization` | R | 3.259 | 3.2305 | −0.88 % | 10 % | PASS |
| `at16-binbin-exchange` | R | 0.302 | 0.37165 | +23.06 % | 30 % | PASS |
| `escape-theory-energy-slope-L0` | own | −3.5 | −3.5026 | −0.0026 | 0.05 | PASS |
| `escape-theory-energy-slope-highL` | own | −4 | −4.02479 | −0.0248 | 0.1 | PASS |
| `escape-theory-mean-e-L0` | own | 0.8133 | 0.812817 | −0.000483 | 0.01 | PASS |
| `escape-theory-mean-e-highL` | own | 0.7 | 0.710776 | +0.0108 | 0.02 | PASS |
| `escape-theory-vs-nbody-light-escaper` | own | 0.769 | 0.75369 | −0.0153 | 0.07 | PASS |

**Result (Table 11).** The set-up of a single hyperbolic passage gives the exact two-body deflection angle. The
cross-sections (in units of $\pi a^2$) come from 12 000 encounters, against 10⁶ in the paper. Single
ionization and exchange + ionization agree to 1.6 % and 0.9 %. The paper states the normalisation
$\pi(a_1^2+a_2^2)$, but its numbers agree with $\pi a^2$ of one binary; the benchmarks use the
latter (a factor-two ambiguity in the paper). Exchanges come out 23 % more frequent than published;
with 91 exchange events this is 1.8 standard deviations, and the tolerance of 30 % is the
statistical error of so few events. Stable triples are rare in both (0.008 against 0.020 ± 0.002).
Two of the 12 000 encounters had no final outcome within the integration limit (0.02 %) and are left
out of the cross-sections. The slopes of the energy distribution of the escaping binary and its mean
eccentricity reproduce the closed forms of Stone & Leigh (2019) to 0.0026–0.025 (slopes) and 0.0005–
0.011 (mean $e$). In resonant binary–single encounters (1 + 0.6 M☉ binary, 0.8 M☉ star) integrated
with N-body, the lightest star escapes in 75.4 % of the cases against 76.9 % from the theory. The 76.9 % is not a
published number: it is SWINGERsim's own evaluation of the escape distribution of Stone & Leigh
(2019) for these masses, so this benchmark compares SWINGERsim's N-body integrations with its own
implementation of the theory.

## Hierarchies of more than four stars {#larger-hierarchies}

SWINGERsim describes any hierarchy of nested binaries (2+2+1, 3+1+1, 2+(2+1), 2+2+2, 3+3 or a
custom nesting) and is tested up to ten stars. Two kinds of evidence cover the larger hierarchies.

**Secular against N-body for a sextuple (level C).** A 2+2+2 sextuple whose binary A is inclined by
60° to the orbit of the quadruple AB, and whose AB orbit is inclined by 50° to the outer orbit, is
integrated with the N-body backend (IAS15 + 1PN, 400 yr, 98 s) and with the nested secular backend
(binary pairs to octupole order; Hamers et al. 2015; Hamers & Portegies Zwart 2016; 1 s). A 2+2+1
quintuple tests the ZLK cycles of the orbit of the two binaries driven by a fifth star.

**V994 Her as a sextuple (level D).** V994 Her is a triply eclipsing sextuple (Zasche et al. 2023): two
binaries A and B on a mutual orbit and a third binary C on a wide orbit around them. That orbit is
unknown: the median of the period distribution that Zasche et al. estimate is close to 3000 yr,
while their photodynamical model uses a fixed circular orbit seen face-on with a period of about
30 kyr. SWINGERsim takes 3000 yr ($e$ = 0.3) for it. With the published masses, radii and the
other orbits (all coplanar), SWINGERsim's N-body run without relativity and tides gives
the apsidal advance that the other stars cause in each binary; with relativity it gives the
relativistic share. Both are expressed as fractions of the observed rates (1.86, 3.66 and 1.68 °/yr).

**Table 12.** Larger hierarchies. Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `sextuple-secular-emax` | own | 0.75207 | 0.751821 | −0.03 % | 2 % | PASS |
| `sextuple-secular-zlk` | own | 70.2 yr | 71.04 yr | +1.20 % | 5 % | PASS |
| `sextuple-nbody-emax` | R | 0.75207 | 0.75207 | +0.00 % | 0.2 % | PASS |
| `quint-out-zlk-emax` | R | 0.81272 | 0.8152 | +0.31 % | 1 % | PASS |
| `v994s-dyn-A` | R | 0.12 | 0.0999365 | −16.72 % | 25 % | PASS |
| `v994s-dyn-B` | R | 0.05 | 0.0519361 | +3.87 % | 25 % | PASS |
| `v994s-dyn-C` | R | 0 | 1.23·10⁻⁶ | +1.23·10⁻⁶ | 0.01 | PASS |
| `v994s-gr-A` | R | 0.10 | 0.0938886 | −6.11 % | 10 % | PASS |
| `v994s-gr-B` | R | 0.08 | 0.0748582 | −6.43 % | 10 % | PASS |
| `v994s-gr-C` | R | 0.06 | 0.0811201 | +35.20 % | 10 % | **EXPLAINED** |

**Result (Table 12).** The nested secular backend gives the largest eccentricity of A to 0.03 % and its ZLK
period to 1.2 % of the N-body run, at a hundredth of the cost. The ranges of the other orbits agree as
well (AB: $e$ = 0.190–0.327 secular, 0.179–0.324 N-body; mutual inclination of B 6.6°–23.5° and
6.6°–23.2°). The N-body result itself is reproducible to the last digit (`sextuple-nbody-emax`). For
V994 Her the dynamical shares are 10.0 % (A), 5.2 % (B) and 0 (C), against the published 12 %, 5 %
and "negligible". The difference for A lies within the effect of the uncertainty of the outer
eccentricity ($e_{AB}$ = 0.687 +0.050/−0.037), which changes the rate by up to 25 % – hence the
tolerance. The relativistic shares are 9.4 % (A) and 7.5 % (B), against the quoted 10 % and 8 %. For C
SWINGERsim gives 8.1 %, the paper 6 %. This difference is EXPLAINED: the 1PN rate computed by hand
from the paper's own Table 3 (1.81 + 1.077 M☉, $P$ = 1.960 d, $e$ = 0.189) is 0.136 °/yr, which is
8.1 % of the observed 1.68 °/yr – SWINGERsim's value. The 6 % is probably a slip in the paper (the
primary's mass alone gives 6 %); A and B, computed the same way, match.

The stability of larger hierarchies is judged by the algebraic criteria applied to every nested
triple and by the N-body ghost test ([Chapter 9](#stability)); the neural classifiers were trained
on quadruples only, and SWINGERsim applies them to every 2+2 and 3+1 sub-quadruple, which is an
approximation outside their training set.

## The environment: gas, the Galactic tide, passing stars and gravitational waves {#environment}

Four effects act on a multiple star from outside or over very long times.

- A **parametric gas phase** (for young systems) damps the semi-major axis, the eccentricity and the
  inclination to the gas with prescribed time scales, $\dot a = -a/\tau_a$, $\dot e = -e/\tau_e$, and
  an alignment of the orbit normal with $\tan(i/2) = \tan(i_0/2)\,e^{-t/\tau_i}$ (exponential migration,
  e.g. Papaloizou & Larwood 2000). The benchmarks check that the integrated orbits follow these
  definitions.
- The **Galactic tide** of the disc (Heisler & Tremaine 1986) drives a ZLK-like oscillation of very
  wide orbits; for an orbit inclined by $i$ to the disc and without the radial terms the maximum
  eccentricity is $e_\mathrm{max} = \sqrt{1 - \tfrac54\cos^2 i}$ (Grishin & Perets 2022, eq. 6).
- **Passing stars** kick a binary; for a distant parabolic passage the first-order change of the
  eccentricity vector is given by Hamers & Samsing (2019, eq. 9).
- **Gravitational waves** shrink and circularise close orbits at the orbit-averaged rates of Peters
  (1964). The observed orbital decay of the binary pulsar PSR B1913+16 (Weisberg & Huang 2016) is the
  level-D reference.

**Table 13.** The environment and gravitational waves. Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `gas-migration-a` | own | 5.48812 au | 5.48812 au | −0.00 % | 10⁻⁴ % | PASS |
| `gas-damping-e` | own | 0.111565 | 0.111565 | +0.00 % | 10⁻⁴ % | PASS |
| `gas-alignment-i` | own | 5.63999° | 5.63943° | −0.000561° | 0.01° | PASS |
| `galactic-tide-emax-60` | own | 0.829156 | 0.829181 | +2.5·10⁻⁵ | 0.01 | PASS |
| `galactic-tide-emax-60-nbody` | R | 0.82916 | 0.82908 | −0.01 % | 1 % | PASS |
| `flyby-first-order-kick` | own | 0 | 0.00594917 | +0.00595 | 0.03 | PASS |
| `flyby-nbody-backend` | R | 0 | 0.00738403 | +0.00738 | 0.05 | PASS |
| `gw-hulse-taylor-pdot` | own | −2.40263·10⁻¹² | −2.4021·10⁻¹² | +0.02 % | 0.1 % | PASS |
| `gw-secular-circular-decay` | own | 1 | 1 | −0.00 % | 0.01 % | PASS |
| `gw-nbody-decay` | R | 1 | 1.0003 | +0.03 % | 1 % | PASS |

**Result (Table 13).** An orbit of 2·10⁴ au inclined by 60° to the disc reaches $e$ = 0.82918 in the secular
backend and 0.82908 in the N-body backend, against the analytic 0.82916; the energy of the averaged
system is conserved to 10⁻¹⁰. The first-order kick formula agrees with a direct IAS15 integration of
a parabolic passage at a periastron distance of 20 $a$ to 0.6 % (0.2 % at 40 $a$); the same passage
integrated as a particle in the N-body backend agrees with it to 0.7 %. The Peters rates give the
orbital decay of PSR B1913+16 (masses 1.438 and 1.390 M☉, $P$ = 0.323 d, $e$ = 0.617) to 0.02 % of the
general-relativity value; the secular backend follows the analytic decay of a circular binary to
better than 10⁻⁴, the N-body backend (where SWINGERsim applies the orbit-averaged loss between
integration steps) the rate to 3·10⁻⁴.
