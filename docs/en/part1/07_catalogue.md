# A catalogue of published multiple stars {#catalogue}

The benchmarks of Chapters 2–4 test one effect at a time. A real multiple star tests them together,
with stars from every part of the Hertzsprung–Russell diagram. The SWINGERsim distribution contains
29 configurations of 23 published systems (Table 34): Wolf–Rayet and O stars, β Cep, δ Sct and Cepheid
pulsators, giants, young T Tauri stars, white dwarfs, a millisecond pulsar, and stars hidden by dust
discs. Each configuration names its sources (paper and table), every conversion (for example, the
relative orbit uses $\omega$ = published $\omega$ of the primary + 180°; the mean anomaly follows from
the time of periastron) and every assumption, marked ASSUMED or UNKNOWN; where a quantity is unknown, an
ensemble draws it from a prior. Every quantity that SWINGERsim computes and that has also been published
is compared below. A difference is called *agrees* if it lies within the published uncertainty or
about 10 %, *explained* if it has been traced to the input, a definition or the published numbers, and
*prediction* if nothing has been published to compare with.

**Table 34.** The catalogue. Column 2 gives the hierarchy in SWINGERsim's notation (2+1 = a binary with
a third star; 1+2 = a single star orbited by a binary), column 3 the stars, column 4 the main source of
the configuration.

| System | Hierarchy | Stars | Source |
|---|---|---|---|
| Algol | 2+1 | B8 V + K0 IV (fills its lobe) + Am | Baron et al. 2012 |
| BU CMi | 2+2 | four A stars | Pribulla et al. 2023 |
| VW LMi | 2+2 | contact W UMa pair + G pair | Pribulla et al. 2020 |
| ξ Tau | 3+1 | B9 V + B9 V, B5 V, F5 V | Nemravová et al. 2016 |
| TIC 120362137 | 3+1 | F + F + F + G | Borkovits et al. 2026 |
| TIC 454140642 | 2+2 | four F stars | Kostov et al. 2021 |
| TIC 470710327 | 2+1 | B + B + O9.5 V | Eisner et al. 2022 |
| δ Ori A | 2+1 | O9.5 II + B2 V + B0 IV | Oplištilová et al. 2023 |
| σ Ori | 2+1 | O9.5 V + B0.5 V + B | Schaefer et al. 2016 |
| λ Sco | 2+1 | β Cep + pre-main-sequence star + B2 V | Tango et al. 2006 |
| Polaris | 2+1 | first-overtone Cepheid + F6 V + F3 V | Evans et al. 2024 |
| HD 181068 | 1+2 | red giant + G + K | Borkovits et al. 2013 |
| PSR J0337+1715 | 2+1 | millisecond pulsar + two white dwarfs | Ransom et al. 2014 |
| KIC 4150611 | 2+(2+1) | δ Sct/γ Dor + M + M, G + G | Hełminiak et al. 2017 |
| HD 98800 | 2+2 + disc | four T Tauri stars, 10 Myr, polar disc | Zúñiga-Fernández et al. 2026; Kennedy et al. 2019 |
| HD 5980 | 2+2 | WR/LBV + WN4, O supergiant + unseen | Koenigsberger et al. 2014 |
| HD 150136 | 2+1 | O3 V + O5.5 V + O6.5 V | Mahy et al. 2018 |
| α Cen | 2+1 | G2 V + K1 V + M5.5 V | Akeson et al. 2021; Kervella et al. 2017 |
| GG Tau | 2+(2+1) + ring | five T Tauri stars / brown dwarf | Köhler 2011; Keppler et al. 2020 |
| TIC 168789840 | 2+2+2 | three eclipsing F + K/M pairs | Powell et al. 2021 |
| ε Aur | 1+1 + disc | F supergiant + B star in a dust disc | Stefanik et al. 2010; Hoard et al. 2010 |
| EE Cep | 1+1 + disc | B5 III + dark disc (illustrative) | Pieńkowski et al. 2020 |
| KH 15D | 2 + ring | T Tauri binary behind a ring (illustrative) | Winn et al. 2006 |

The dust-disc systems are discussed in [Chapter 4](#discs); the others follow, grouped by the physics
they test. Classes in the sense of Table 1: precession and timing runs use the N-body backend (R) or the
secular backends (own), the evolving runs and pulsation timing use SWINGERsim's own code.

## Apsidal and nodal precession {#catalogue-precession}

**Method.** The turning of the ellipses (apsidal) and of the orbital planes (nodal) is measured directly
from the integrated orbital elements – with all stars, relativity, tides and rotation – and compared
with the published rates, which the authors derived from eclipse timing, photodynamical models or
analytic formulae. Runs cover 30–1000 yr.

**Table 35.** Apsidal and nodal precession in the catalogue. Column 4 is SWINGERsim's value, column 5
the published one, column 6 their difference.

| System | Quantity | SWINGERsim | Published | Difference | Status |
|---|---|---|---|---|---|
| BU CMi | apsidal period of A | 27.8 yr | 28.7 yr (Pribulla et al. 2023) | −3 % | agrees |
| BU CMi | apsidal period of B | 24.8 yr | 25.1 yr | −1 % | agrees |
| ξ Tau | apsidal rate of the 145-d orbit | 2.14 °/yr | 2.02 ± 0.31 °/yr | +6 % (0.4σ) | agrees |
| TIC 120362137 | apsidal period of the 51-d orbit | 32.4 yr | 34.79 yr (Borkovits et al. 2026) | −7 % | agrees |
| TIC 120362137 | nodal period of the outer planes | 282 yr | 265.8 yr | +6 % | agrees |
| TIC 454140642 | apsidal periods of A, B | 69.7, 90.9 yr | 70, 97 yr (Kostov et al. 2021) | −0.4 %, −6 % | agrees |
| HD 181068 | nodal period | 12.8 yr | 13–14 yr (Borkovits et al. 2013) | −5 % | agrees |
| VW LMi | apsidal period of the detached pair | 95.1 yr | 95.5 yr (the paper's own N-body model) | −0.4 % | agrees |
| δ Ori A | apsidal rate of the 5.73-d orbit | 2.42 °/yr | 1.45 °/yr | +67 % | explained: input |
| TIC 470710327 | nodal period of the eclipsing pair | 14.1 yr | – | – | prediction |

![**Figure 13.** Precession measured in four catalogue runs; in every panel an angle grows steadily with time, the slope is the precession rate and 360° divided by the rate is the period. *Top left:* longitude of periastron of the two pairs of BU CMi (A and B). *Top right:* the same for the detached pair of VW LMi (blue) against the dashed line of the 78.6-yr period fitted by Pribulla et al. (2020) to three years of data. *Bottom left:* nodes of the three orbits of TIC 120362137 in the invariable plane (the plane perpendicular to the total angular momentum) over 1000 yr. *Bottom right:* inclinations of the two inner orbits of ξ Tau to the sky over 300 yr.](../../figures/studies/catalogue_precession.png)

**Result (Table 35).** Figure 13 shows four of the runs. In BU CMi, 360° are reached after about 28 yr (A) and
25 yr (B), 3 % and 1 % faster than published. In VW LMi the detached pair turns in 95 yr, slower than
the dashed 78.6-yr line: the 78.6 yr of Pribulla et al. (2020) is a fit to only three years of radial
velocities, whereas the paper's own analytic and N-body models give 96.2 and 95.5 yr (and DDO data of
1998–2008 give 91 ± 11 yr); SWINGERsim agrees with the paper's dynamical model to 0.4 %. The nodes of
TIC 120362137 all turn together, backwards, by 360° in about 280 yr: the planes precess as one block.
For this system, Table 1 of Borkovits et al. (2026) gives 34.79 yr for the eccentric 51-d orbit,
against SWINGERsim's 32.4 yr (−7 %). The nearly circular 3.3-d pair ($e$ = 0.007) is pushed by the
tertiary to $e$ between 0.0001 and 0.016; its osculating $\omega$ jumps, and a measured apsidal period is
not defined (the paper gives only a nominal theoretical value for it). The inclination of the eclipsing
inner pair of ξ Tau wobbles by ±0.6° with a period of about 20 yr, so its eclipse depth should vary with
that period – a prediction.

The one large difference, δ Ori A, comes from the input. Its apsidal motion is dominated by tides and
rotation (the binary alone gives an apsidal period of 155 yr), and the tidal part scales as $k_2R^5$ of
the primary. With the radius of Oplištilová et al. (2023) and $k_2$ from Claret (2004), SWINGERsim is
1.7 times faster than the rate of Pablo et al. (2015); the published 1.45 °/yr requires $R_1\approx$
11.5 R☉ (11 % smaller) or $\log k_2\approx$ −2.85. The published rate is therefore a measurement of
$k_2R^5$ of Aa1 rather than a test of the dynamics. For TIC 470710327, with a mutual inclination of 17°
the 1.1-d pair precesses in 14 yr, and its eclipses should change within years (a prediction).

## Eclipse-timing variations {#catalogue-timings}

**Method.** `swingersim timings` computes the exact times of conjunction from the N-body integration,
including the light travel time, for TIC 470710327 (5 yr), TIC 454140642 (5 yr), BU CMi (30 yr) and
ξ Tau (10 yr). The resulting O−C diagrams are what an observer would obtain.

![**Figure 14.** Synthetic O−C diagrams of four catalogue systems; horizontal axes: time; vertical axes: how much earlier or later each eclipse comes than a strictly periodic clock predicts. *Top left* (TIC 470710327, minutes): the 1.1-d pair with the 52-d orbit of the third star. *Top right* (ξ Tau, minutes): the 7.1-d pair with the 145-d orbit. *Bottom left* (TIC 454140642, minutes): the primary minima of both pairs, A (blue) and B (orange), with the 432-d outer orbit. *Bottom right* (BU CMi, hours): the primary (blue) and secondary (orange) minima of pair A over 30 yr.](../../figures/part1/catalogue_timings.png)

**Table 36.** Eclipse-timing variations. "Semi-amplitude" is half the peak-to-peak O−C of a fitted
periodic curve.

| System | Quantity | SWINGERsim | Published | Difference | Status |
|---|---|---|---|---|---|
| TIC 470710327 | ETV semi-amplitude of the 1.1-d pair (52-d orbit) | 4.97 min | ≈ 6 min (read from a figure) | −17 % | agrees |
| ξ Tau | ETV semi-amplitude of the 7.1-d pair (145-d orbit) | 19.7 min | – | – | prediction |
| TIC 454140642 | ETV semi-amplitudes of A and B (432-d orbit) | 51.4 and 40.2 min, in phase | dynamical ETV dominates (Kostov et al. 2021) | – | agrees |
| BU CMi | apsidal period of A from the O−C crossings | 27.78 yr | 28.7 yr | −3 % | agrees |
| BU CMi | apsidal period of B from the O−C crossings | 24.79 yr | 25.1 yr | −1 % | agrees |

**Result (Table 36).** In Figure 14 the eclipses of the 1.1-d pair of TIC 470710327 run about ±5 min early and late
with the 52-d period of the third star; the published value, read from a figure, is about 6 min. ξ Tau
shows ±20 min with the 145-d orbit; the pure light-time effect of that orbit would be at most 8.6 min,
so the dynamical ETV dominates – the tertiary is close enough to change the period of the pair within
each outer orbit. The bottom-left panel is the instructive one. If only the light travel time mattered,
the two pairs of a 2+2 would move in antiphase (A late when B is early), here with about ±8 min (7.6 and
7.9 min). The computed ETVs are five to seven times larger and *in phase*: in this compact quadruple
($P_\mathrm{out}/P_\mathrm{in}\approx$ 32–42) the gravitational pull of each pair on the other dominates,
as Kostov et al. (2021) describe. Antiphase is therefore not a general signature of a 2+2, and a
comparison with the minima of a real compact 2+2 must use the full dynamical model, as
`swingersim timings` does. In BU CMi the primary and secondary minima of A drift apart and cross again
as the ellipse turns; the time between two crossings is half the apsidal period, which gives 27.78 and
24.79 yr – the same as measured from $\omega$ (27.76 and 24.82 yr) by an independent route.

## ZLK cycles and long-term survival {#catalogue-zlk}

**Method.** Three systems with strongly inclined orbits: σ Ori (secular, 3 Myr), HD 98800 (secular,
10 Myr) and α Cen with Proxima (secular with the Galactic tide, 1 Gyr).

**Table 37.** ZLK cycles and survival. Column 3 gives the result of the run, column 4 the published
observation.

| System | Run | SWINGERsim | Published | Status |
|---|---|---|---|---|
| σ Ori | secular, 3 Myr | $e$ of Aa–Ab 0.69–0.90, $i_\mathrm{mut}$ 116°–134°, ZLK period 39.5 kyr; closest periastron 3.4 $(R_1+R_2)$ | observed $e$ = 0.78 (Schaefer et al. 2016) | agrees; survives |
| HD 98800 | secular, 10 Myr | $e_A$ 0.27–0.63, $e_B$ 0.42–0.80, ZLK periods 46 and 29 kyr | observed 0.48 and 0.78 (Zúñiga-Fernández et al. 2026) | agrees; survives |
| α Cen | secular + Galactic tide, 1 Gyr | Proxima stays bound; $e_\mathrm{out}$ 0.50 → 0.92, periastron down to 724 au; $i_\mathrm{mut}$ 64°–83° | bound (Kervella et al. 2017) | agrees; the growth of $e$ is a prediction |

![**Figure 15.** ZLK cycles and the Galactic tide: eccentricity (vertical axes) against time. *Left:* σ Ori Aa–Ab, the first 0.15 of the 3-Myr run. *Middle:* the two pairs of HD 98800, A (blue) and B (orange), the first 0.5 of the 10-Myr run. *Right:* the orbit of Proxima around α Cen AB over 1 Gyr.](../../figures/studies/zlk_systems.png)

**Result (Table 37).** In Figure 15 (left) the eccentricity of σ Ori Aa–Ab swings regularly between 0.69 and 0.90
every 40 kyr, and the observed 0.78 lies in the middle of the swing; the closest approach stays 3.4 times
the sum of the radii, so the system survives its age. In HD 98800 both pairs oscillate (0.27–0.63 and
0.42–0.80), and the observed 0.48 and 0.78 lie within these ranges: the large observed eccentricities are
naturally explained by ZLK cycles, as Zúñiga-Fernández et al. (2026) suggest. In α Cen the Galactic tide
raises the eccentricity of Proxima's orbit slowly from 0.50 to 0.92 over a billion years without unbinding
it; the ensemble of [Ensembles of published systems](#catalogue-ensembles) confirms this with passing
stars. The tide itself is validated against the analytic maximum eccentricity of a wide binary
([Chapter 2](#environment)).

## Growing stars and winds {#catalogue-growing}

**Method.** Systems in which the stars change during the run: HD 181068 (N-body over 30 yr; evolving
backend with MIST from 362.2 Myr), HD 150136 (evolving, MIST), HD 5980 (secular with winds, 100 kyr) and
PSR J0337+1715 (N-body over 30 yr) as a check of exact coplanarity.

**Table 38.** Growing stars and winds.

| System | Run | SWINGERsim | Published | Status |
|---|---|---|---|---|
| HD 181068 | N-body, 30 yr | coplanar ($i_\mathrm{mut}$ 0.8°), nodal period 12.8 yr | 0.8 ± 1.4°; 13–14 yr (Borkovits et al. 2013) | agrees |
| HD 181068 | evolving, MIST | the giant grows from 12.4 to 38.3 R☉ and fills its Roche lobe after 2.1 Myr; tertiary mass transfer ends in a triple common envelope at 2.1055 Myr (age 364.3 Myr) | – | prediction |
| HD 150136 | evolving, MIST | the O3 V primary fills its lobe after 2.45 Myr (age 3.45 Myr) | – | prediction |
| HD 5980 | secular with winds, 100 kyr | $a_\mathrm{out}$ 41.3 → 42.9 au = $a_0M_0/M$ exactly; masses 61 → 56.5 and 66 → 64 M☉ | Jeans (1924) | agrees |
| PSR J0337+1715 | N-body, 30 yr | $i_\mathrm{mut}$ stays 0.00° | coplanar to 0.01° (Ransom et al. 2014) | agrees |

![**Figure 16.** Growing stars. *Left and middle:* the radius of a growing star (blue) and the size of its Roche lobe (dashed) against time, for the giant of HD 181068 (whose lobe is around the inner pair) and for the O3 V primary of HD 150136; where the blue line reaches the dashed one the star starts to transfer mass. *Right:* the semi-major axis of the outer orbit of HD 5980 divided by its initial value (blue) and the ratio of the initial to the present total mass (dashed) against time.](../../figures/studies/growing_stars.png)

**Result (Table 38).** In Figure 16 (left) the giant of HD 181068 grows from 12 to 38 R☉ in 2.1 Myr until it meets
its Roche lobe of 38.3 R☉. Because the giant (3.0 M☉) is heavier than the inner pair (0.92 + 0.87 M☉),
the transfer onto the pair is unstable (Hamers et al. 2022), and SWINGERsim ends the run in a triple
common envelope at the moment of contact. HD 181068 is the textbook case for mass transfer from a
tertiary onto a binary (de Vries et al. 2014; Tokovinin & Moe 2020). The O3 V primary of HD 150136 grows
from 9.5 to 14.6 R☉ in about 2.5 Myr while its lobe grows a little, because its wind widens the orbit; it
will become a very young, massive Algol-like system. In the right panel the two lines for HD 5980
coincide: the outer orbit widens exactly as fast as the system loses mass, $a\propto 1/M$, which is
Jeans' law for slow mass loss. PSR J0337+1715, a millisecond pulsar with two white dwarfs, stays exactly
coplanar, as observed.

## Pulsation timing {#catalogue-pulsations}

**Method.** SWINGERsim's light curves contain the pulsations, and the delay is measured as an observer
does it: a sinusoid at the pulsation frequency is fitted in short windows and its phase is followed
([Chapter 4](#pulsations)). Three systems: λ Sco (β Cep pulsator, 1100 d, V band), Polaris
(first-overtone Cepheid with the observed period of 3.97 d, 30 yr, V) and KIC 4150611 (δ Sct, 200 d, TESS
band).

**Table 39.** Pulsation timing. Column 3 is the semi-amplitude of the delays measured from the
synthetic light curve, column 4 the light-time semi-amplitude of the orbit used in the run (Irwin 1952),
column 5 the published value.

| System | Pulsation frequency | Measured | Orbit of the run | Published | Difference | Status |
|---|---|---|---|---|---|---|
| λ Sco | 4.679 c/d | 1067 s | 1069 s | 1193 s | −11 % | explained: input |
| Polaris | 0.252 c/d | 1309 s | 1336 s | 1335 s | −2 % | agrees |
| KIC 4150611 | 20.243 c/d | 97.1 s | 97.9 s | 98 ± 9.5 s | −1 % | agrees |

![**Figure 17.** Pulsation timing in three catalogue systems. *Top:* a few days of each synthetic light curve (relative flux against time); the regular waves are the pulsations, for λ Sco also two eclipses. *Bottom:* the delay of the pulsation phase in seconds (points), measured in short windows along the whole light curve, with the light-time curve of the orbit that carries the pulsating star (line); the dashed lines mark ± the published semi-amplitude. The text above each panel gives the measured, model and published values.](../../figures/part1/pulsator_delays.png)

**Result (Table 39).** In Figure 17 the measured delays (points) follow the light-time curve of the run's own orbit
(line) to within 2 % in all three systems: the measurement is accurate. The remaining differences to the
published numbers come from the input. For λ Sco the configuration uses the interferometric orbit and
masses of Tango et al. (2006), which give $a_1\sin i$ = 2.15 au; the published 1193 s belong to the older
orbit of Uytterhoeven et al. (2004). For Polaris the 1430 s often quoted is $a_1\sin i/c$; the
semi-amplitude of the light-time curve is smaller by $\sqrt{1-e^2\cos^2\omega}$, 1335 s (Evans et al.
2024), which SWINGERsim reproduces to 2 % (only one 30-yr orbit is computed, so the curve is not centred
on zero; its height from lowest to highest point is twice the semi-amplitude). The first-overtone period
that the stellar model of Polaris predicts from the period–luminosity–mass relation is 3.62 d against the
observed 3.97 d (−9 %), within the scatter of such relations; the run uses the observed period. For
KIC 4150611 the configuration takes the companion mass $M_\mathrm{Ab}$ = 0.86 M☉, consistent with the
size and period ratios of Kemp et al. (2024); their dynamical mass ratio of 3.61 contradicts these ratios
(which imply 1.70) and would give only 54 s. With consistent masses the orbit gives the published 98 s
(Hełminiak et al. 2017 from radial velocities; Shibahashi & Kurtz 2012 from frequency modulation).
Pulsation timing is thus a sensitive test of the masses in a configuration.

## Ensembles of published systems {#catalogue-ensembles}

**Method.** Where a parameter of a real system is unknown (an outer orbit, a node, a mutual inclination),
an ensemble draws it from a prior and asks what the observed system implies.

**Table 40.** Ensembles of published systems. Column 3 is the number of members, column 4 the question,
column 5 the answer.

| System | Backend | Members | Question | Answer |
|---|---|---|---|---|
| HD 150136 | N-body | 500 | does the unknown mutual inclination matter within a century? | $i_\mathrm{mut}$ changes by < 0.04°; the 2.7-d pair stays eclipsing; its inclination drifts by 0–0.57° per century depending on $i_\mathrm{mut}$ |
| HD 5980 | N-body | 1000 | is the 2+2 bound and safe? | 29.7 % collide (median largest $e_B$ 0.91), mostly for outer periods below 20 yr |
| GG Tau | N-body | 500 | fate of the young quintuple in 20 kyr | intact 80.8 %, reconfiguration 9.2 %, collision 7.6 %, ejection 2.4 % |
| TIC 168789840 | secular | 500 | how often does each pair eclipse? | A 31 %, B 91 %, C 35 % of the time; all three together 11 % |
| KIC 4150611 | secular | 500 | which outer orbits keep the Ab pair circular, and is the system stable? | with the consistent masses all members stay intact and stable by MA01 and Vynatheya et al. (2022); for unknown outer orbits 42.6 % keep $e$ < 0.05 over 1 Myr, the others are pumped up to 0.7 |
| α Cen | secular + Galactic tide + passing stars | 300 | does Proxima stay bound for 1 Gyr? | 100 % bound; median largest $e$ 0.90; median closest periastron 849 au, minimum 70 au |

![**Figure 18.** Ensembles of published systems, one panel each. *HD 5980:* fraction of members that end in a collision against the outer period. *GG Tau:* fates (blue intact, orange collision, green reconfiguration, violet ejection) per range of the Aa–Ab semi-major axis. *TIC 168789840:* mean fraction of time each of the three pairs eclipses. *KIC 4150611:* each point is a member, the largest eccentricity of the Ab pair against the size of the outer orbit; dashed: $e$ = 0.05. *α Cen:* histogram of the largest eccentricity of Proxima's orbit over 1 Gyr. *HD 150136:* each point is a member, the rate of change of the inclination of the eclipsing pair against the unknown mutual inclination.](../../figures/part1/catalogue_ensembles.png)

**Result (Table 40).** Figure 18 shows the six ensembles. For HD 5980 most members with outer periods below 20 yr
collide, long ones are safe: if the two binaries form a bound 2+2, the outer period is likely above 20 yr
– a dynamical constraint for the spectroscopists. GG Tau stays intact in about 80 % of the members
everywhere; 74 % of the members that MA01 calls unstable stay intact for 20 kyr, which only shows that
20 kyr is short compared with the time scale of MA01 (MA01 calls 34.8 % stable). In TIC 168789840 pair B
eclipses almost always, A and C a third of the time: the triple eclipses seen by TESS are a 1-in-9
coincidence. The KIC 4150611 members either keep the Ab pair circular or pump it to 0.6–0.7 at any outer
size. Proxima stays bound in all 300 members, mostly with largest eccentricities of 0.85–0.97 (Kervella
et al. 2017 find it bound). For HD 150136 the drift of the inclination vanishes at mutual inclinations
of 0°, 90° and 180°, so a measured change of the eclipse depth would determine the mutual inclination;
on observable time scales the effects are negligible, as Mahy et al. (2018) state.

## Summary of the catalogue {#catalogue-summary}

All 23 systems ran without failure. Every quantity that can be compared with a publication – apsidal
and nodal periods, eclipse-timing amplitudes, ZLK ranges, coplanarities, light-time amplitudes,
Roche-lobe sizes and disc geometries – agrees within the published uncertainty or about 10 %, or the
difference has been traced to the input of the configuration (δ Ori A, λ Sco) or to a different
definition of the published quantity (VW LMi, TIC 120362137, Polaris). The figure of all comparisons on
one scale is in [Chapter 11](#conclusions).
