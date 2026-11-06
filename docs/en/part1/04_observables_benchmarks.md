# Observables: light curves, timings and pulsations {#observables}

Astronomers do not observe orbits; they observe light. This chapter tests how SWINGERsim turns its
orbits and stars into what a telescope records: light curves of eclipsing stars (spherical,
distorted and in contact), the dimming by dust discs, radial velocities, the times of eclipse minima
and their deviations from a strict period, and the timing of pulsations. All modules in this chapter
are SWINGERsim's own code, except the PHOEBE hybrid model, which calls the external code PHOEBE 2.5
(Prša et al. 2016) for the light of Roche-distorted stars.

## Light curves of spherical stars {#light-curves}

**Definitions.** SWINGERsim's basic photometric model (called **tier A** in the software) treats stars
as spheres with a radially symmetric limb darkening (the four-coefficient law of Claret 2000) and
computes the flux blocked in an eclipse by Gauss–Legendre integration over the visible disk. The limb
darkening comes from the tables of Claret & Bloemen (2011, ATLAS atmospheres) and, for stars hotter
than 15 000 K, of Reeve & Howarth (2016, non-LTE TLUSTY atmospheres). Ellipsoidal variation – the
light change of a tidally elongated star – is added as the leading term of Morris (1985).

**References.** The analytic occultation of a limb-darkened star by a dark disk (Mandel & Agol 2002)
as implemented in batman (Kreidberg 2015); the eclipse depths of a detached binary of spheres in ellc
(Maxted 2016); the ellipsoidal variation of ellc's Roche model; PHOEBE 2.5 with the same physics
(spheres, black bodies, Claret limb darkening); the published TLUSTY coefficients and the overlap of
TLUSTY and ATLAS at 15 000–30 000 K.

**Table 23.** Light curves of spherical stars and limb darkening. Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `lc-occultation-small` | own | 0.987891 | 0.987891 | −1.9·10⁻¹² | 10⁻⁷ | PASS |
| `lc-occultation-limb` | own | 0.945236 | 0.945236 | +3.7·10⁻⁹ | 10⁻⁷ | PASS |
| `lc-occultation-large` | own | 0.128979 | 0.128979 | −5.5·10⁻⁹ | 10⁻⁷ | PASS |
| `lc-detached-primary` | own | 0.380821 | 0.380868 | +4.7·10⁻⁵ | 10⁻⁴ | PASS |
| `lc-detached-secondary` | own | 0.144658 | 0.144654 | −3.9·10⁻⁶ | 10⁻⁴ | PASS |
| `lc-ellipsoidal` | own | 0.0027949 | 0.00260619 | −6.75 % | 10 % | PASS |
| `phoebe-sphere-vs-tier-a` | own | 0 | 0.00037 | +0.00037 | 0.001 | PASS |
| `ld-tlusty-row` | own/tab | 1.18262 | 1.18262 | +0 | 10⁻⁴ | PASS |
| `ld-tlusty-vs-atlas` | own/tab | 1 | 0.997587 | −0.24 % | 2 % | PASS |

**Result (Table 23).** For small (radius ratio 0.1), grazing and large (1.2) occulters the flux agrees with the
analytic model to better than 10⁻⁸: the integration is exact for spheres in practice. The eclipse
depths of a detached binary agree with ellc to 5·10⁻⁵, which is the accuracy of ellc's integration
grid. Against PHOEBE 2.5 with the same physics, the largest difference of the normalised light curves
over a whole orbit ($e$ = 0.1) is 3.7·10⁻⁴ of the flux. The small-distortion ellipsoidal term is
6.8 % below ellc's Roche model at $R/a$ = 0.15; at $R/a$ = 0.10, 0.20 and 0.25 the deficit is 3, 11
and 17 %. This benchmark marks the limit of the spherical model: closer pairs need the surface meshes
of the next section. The TLUSTY coefficients are reproduced at the grid points, and the disk-averaged
intensities of TLUSTY and ATLAS agree to 0.2–0.7 % where both apply, so the switch between them at
15 000 K is smooth.

## Distorted and contact stars on surface meshes {#meshes}

**Definitions.** For stars that are far from spherical SWINGERsim builds **surface meshes**: the
surface of a rigidly rotating oblate star or of a star filling (or overfilling) its Roche lobe is
divided into elements of 64 × 128 directions, each with its own temperature (gravity darkening), limb
darkening and projected area; the eclipses are computed element by element. The same meshes give the
**centre of light** of a star: in a contact binary the radial velocity measured from spectral lines
is the flux-weighted velocity of the visible surface, which differs from the velocity of the centre of
mass (and during an eclipse shows the Rossiter–McLaughlin effect; Rossiter 1924; McLaughlin 1924).

**Table 24.** Surface meshes, contact binaries and centre-of-light velocities. Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `mesh-sphere-eclipse` | own | 0 | 0.00026 | +0.00026 | 0.001 | PASS |
| `contact-vs-phoebe` | own | 0 mag | 0.0030 mag | +0.0030 mag | 0.03 mag | PASS |
| `rv-light-vs-phoebe` | own | 0 | 0.022 | +0.022 | 0.05 | PASS |

**Result (Table 24).** Meshes of spheres reproduce the exact sphere eclipses of the previous section to
2.6·10⁻⁴ of the flux; the edges of the stars in front are anti-aliased by the projected size of each
element. A contact binary ($q$ = 0.42, fill-out factor 0.3, $i$ = 80°) has the same V amplitude as in
PHOEBE 2.5 with black bodies and the same limb darkening to 0.003 mag. The centre-of-light velocity of
its secondary differs from PHOEBE's flux-weighted velocity by at most 2.2 % of the 52-km/s amplitude;
the largest differences occur inside the eclipses, where the two codes' meshes resolve the occulted
edge differently. Automatic tests add that out of eclipse the offset of a spherical star is exactly
zero, that the Rossiter–McLaughlin signal is antisymmetric about mid-eclipse with the expected size
$v_\mathrm{rot}(R_\mathrm{occ}/R)^2$, and that velocities of a contact binary computed with the centre of
light but fitted with the centre of mass give a $\chi^2$ more than twice as large and masses several per
cent too large, whereas fitting with `fit.rv.centre: light` returns the masses within 3σ.

## Dust discs as occulters {#discs}

**Definitions.** Some eclipses are caused not by a star but by a dark dust disc around a companion, or
by a ring around a whole binary. SWINGERsim represents a disc as a **slab** – a cylinder of outer radius
$r_\mathrm{out}$, an optional inner hole $r_\mathrm{in}$, thickness $H$ and optical depth $\tau$. A star
behind it is dimmed by $e^{-\tau L/H}$, where $L$ is the length of the line of sight inside the slab.
The edge can be sharp or soft, the disc can emit and scatter light, and it can be fixed, precess at a
prescribed rate or precess as a rigid body under the torque of its binary or of a companion. A disc
around an eccentric binary can settle **polar**, perpendicular to the binary orbit; above a critical
inclination it librates around the polar orientation (Farago & Laskar 2010; Martin & Lubow 2017).

**References.** Exact geometry for the transmission; the test-particle precession rate for a narrow
rigid ring (Zanazzi & Lai 2018, eq. 2); $\cos I_\mathrm{crit} = \sqrt{5e^2/(1+4e^2)}$ for the polar
libration; and three observed systems: ε Aur (the F supergiant eclipsed every 27.1 yr for about two
years by the disc around its B-star companion; orbit of Stefanik et al. 2010, low-mass model of Hoard et
al. 2010, disc of Kloppenborg et al. 2015), HD 98800 B (a polar circumbinary disc, inclined 26° to the
sky; Kennedy et al. 2019; orbit of Zúñiga-Fernández et al. 2026) and EE Cep (a B star eclipsed every
5.6 yr with depths that change with a period of 11–12 orbits, 62–67 yr; Gałan et al. 2012; Pieńkowski
et al. 2020).

**Table 25.** Dust discs. Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `disc-strip-uniform` | own | 0.623838 | 0.62475 | +0.000912 | 0.005 | PASS |
| `disc-face-on-tau` | own | 0.496585 | 0.496585 | +0 | 10⁻⁹ | PASS |
| `disc-inclined-tau` | own | 0.246597 | 0.246597 | −4·10⁻¹⁶ | 10⁻⁹ | PASS |
| `disc-ring-precession` | own | 1 | 1.00018 | +0.02 % | 0.1 % | PASS |
| `disc-polar-critical` | own | 37.7612° | 37.7612° | +6·10⁻⁹° | 0.2° | PASS |
| `epsaur-mid-eclipse` | own | JD 2455429.5 | JD 2455414.5 | −15 d | 30 d | PASS |
| `epsaur-duration` | own | 740 d | 638 d | −13.78 % | 15 % | PASS |
| `epsaur-depth` | own | 0.75 mag | 0.83134 mag | +0.0813 mag | 0.15 mag | PASS |
| `hd98800-polar-disc` | own | 26° | 28.4735° | +2.47° | 3° | PASS |
| `eecep-disc-precession` | own | 65 yr | 183.163 yr | +181.79 % | 15 % | **EXPLAINED** |

**Result: geometry and dynamics (Table 25).** An opaque edge-on slab of half-thickness 0.3 $R$ across a uniform
star transmits 0.6248 of the light against the exact 0.6238: the stellar disk is sampled on a
sunflower lattice of 4000 points, which resolves a sharp edge to about 10⁻³ of the flux. Face-on and
inclined discs of optical depth 0.7 transmit exactly $e^{-0.7}$ and $e^{-1.4}$. A narrow rigid ring at
four binary separations precesses at the test-particle rate to 0.02 %, and the rigid-disc equation
reproduces the critical inclination of polar libration for $e$ = 0.5 (37.76°) exactly.

![**Figure 5.** Dust discs as occulters, with a sharp-edged opaque disc. In all panels the vertical axis is the dimming in magnitudes, growing downwards as astronomers plot it. *Left* (ε Aur): the computed light curve of the 2009–2011 eclipse against Julian date; dashed vertical lines: the observed first and last contacts; dotted horizontal line: the observed depth. *Middle* (EE Cep, illustrative disc): the depth of each eclipse (points) over 70 yr as the disc precesses in a prescribed 65 yr; the grey band is the observed range of depths. *Right* (KH 15D, illustrative ring): 25 years of the light curve of a young binary disappearing behind the edge of a precessing circumbinary ring and re-emerging.](../../figures/discs.png)

**Result: observed systems.** Figure 5 shows the three discs with a sharp edge. For ε Aur (left), the
time of mid-eclipse follows from the spectroscopic orbit alone and falls 15 d before the observed
mid-point of the contacts (2009 August 16 to 2011 August 26). The eclipse is 638 d long against the
observed 740 d (−14 %) and 0.83 mag deep against 0.75 mag. The shorter duration is the sharp edge:
the observed ingress takes about 145 d (Stencel 2013), so the real edge is soft. The depth is
0.08 mag too large because the model omits the light of the B star inside the disc and the light
scattered by the disc. For HD 98800 B, a disc put perpendicular to the orbit of B and along its
eccentricity vector – the polar alignment predicted for discs around eccentric binaries – is inclined
28.5° to the sky, computed from the published orbit without any fit; the observed disc is inclined 26°.
In the middle panel, an illustrative EE Cep disc (1.6 au across, $\tau$ = 1, tilted 15°) whose
precession is *prescribed* at 65 yr produces eclipse depths of 0.19–1.73 mag, close to the observed
0.5–2 mag. In the right panel, an illustrative KH 15D ring (opaque, 3–5 au, tilted 20°, precessing in
600 yr) reproduces the observed sequence: first one star of the binary peeks out near apastron
(spikes), then both are hidden for years, then the binary re-emerges.

![**Figure 6.** Discs with soft edges and rigid dynamics, contact evolution and the Pejcha fractions. *Top left* (ε Aur): dimming against Julian date for the opaque sharp-edged disc (orange) and for a disc with a soft edge and its own light (blue); grey area: the observed first to last contact. *Top right* (EE Cep): the depth of each eclipse (magnitudes, points connected by lines) over 65 years, for the disc of outer radius 0.8 au whose 65-yr precession is prescribed (orange) and for the rigid disc of outer radius 1.6 au whose torque-driven precession has 65 yr (blue); grey band: observed depths. *Bottom left:* masses of the two stars (solid, left axis) and their separation (dashed, right axis) of a 1.0 + 0.5 M☉ binary evolving into contact and merging (see [Chapter 8](#life-cycles)). *Bottom right:* the fraction of the Pejcha et al. (2013) set (Vokrouhlický set-up, 400 members) whose binary A (blue) or B (orange) has reached $1-e<10^{-3}$ after secular runs of 11 kyr and 50 kyr; dotted: the published fractions.](../../figures/part1/v027_jobs.png)

**Result: soft edges and torques.** Figure 6 shows what happens when the discs get a soft edge, their
own light, and a precession computed from the torque instead of prescribed. For ε Aur (top left) the
soft edge makes the ingress take 141 d, as observed (about 145 d; the sharp edge gives 89 d), but the
eclipse becomes 829 d long, 12 % longer than observed; the width of the soft edge is assumed, not
fitted, and the disc adds only 0.1 % of light. For EE Cep the illustrative disc of outer radius 0.8 au, made rigid,
precesses in 183 yr under the torque of the B star – not in 65 yr. The benchmark is EXPLAINED as a
constraint on the disc rather than an error: the rate scales as $r_\mathrm{out}^{3/2}$, and 65 yr
requires $r_\mathrm{out}\approx$ 1.6 au (or a heavier B star). With that size and the illustrative
thickness and optical depth, however, the eclipses become 3.7–14 mag deep and 36–158 d long (top right,
blue), far outside the observed 0.5–2 mag. A disc that is rigid and precesses with the observed period
must therefore also be thinner or less opaque than the illustrative one. The KH 15D ring of 3–5 au,
made rigid, turns by only 2.6° in 60 yr (period about 6500 yr) and hides the binary 98 % of the time;
the change observed within decades needs 1000–3000 yr, so the ring must be closer to the binary. These
are constraints on the discs, which SWINGERsim derives because its disc dynamics passes the analytic
benchmarks above. The bottom panels of Figure 6 belong to [Chapter 8](#life-cycles) and to the Pejcha
ensembles of [Chapter 2](#pejcha2013): the published fractions (dotted) lie between the values after
11 kyr and after 50 kyr.

## Eclipse timings and O−C diagrams {#timings}

**Definitions.** An observer predicts the times of minima of an eclipsing binary with a constant
period, $T_n = T_0 + nP$, and plots the **O−C diagram** (observed minus computed time against cycle
number $n$). Any structure in it is a signal. The **light-travel-time effect** (LTTE) of a third body
moves the minima by up to $a_{12}\sin i/c$, where $a_{12}$ is the semi-major axis of the binary's orbit
around the common centre of mass (Irwin 1952). The **dynamical** eclipse-timing variation (ETV) is the
change of the binary's period by the pull of the third body. **Apsidal motion** moves the primary and
secondary minima in opposite directions. Inside the binary itself, the secondary eclipse is delayed by
the **Rømer delay** $(P/\pi c)(K_2-K_1)$, because the two stars are at different distances from us at
the two conjunctions (Kaplan 2010). A steady period change $\dot P$ bends the O−C into a parabola,
$\tfrac12\dot P P n^2$ (Kwee 1958; Hilditch 2001).

SWINGERsim offers three ways to compute minima: exact conjunctions of the N-body integration,
including the light travel time (`swingersim timings`); minima measured from a synthetic light curve
as an observer would; and the PHOEBE hybrid, which computes the light of Roche-distorted stars with
PHOEBE 2.5 in chunks of the N-body orbit. For observed minima it provides two fits: **stage 1**, the
analytic O−C model used by the SQUADRA group in OCFit (Gajdoš & Parimucha 2019) – LTTE, apsidal motion
in the series of Giménez & Bastero (1995) and a quadratic term for every pair – and **stage 2**, the
joint N-body model of all eclipsing pairs (tested in Parts II–IV).

**Table 26.** Eclipse timings and the stage-1 fit. Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `timing-ltte-amplitude` | own | 41.5837 min | 41.5863 min | +0.01 % | 0.1 % | PASS |
| `phenom-romer-delay` | own | 27.4517 s | 27.3578 s | −0.34 % | 1.5 % | PASS |
| `phoebe-romer-delay` | ext | 27.452 s | 27.287 s | −0.60 % | 2 % | PASS |
| `pdot-parabola` | R | 1 | 1.00005 | +0.00 % | 0.1 % | PASS |
| `ocfit-gimenez-bastero` | own | 0 s | 0.00063 s | +0.00063 s | 0.01 s | PASS |
| `ocfit-lite-amplitude` | own | 1 | 1 | −0.00 % | 10⁻⁴ % | PASS |
| `ocfit-nbody-asini-A` | own | 1 | 0.99864 | −0.14 % | 1 % | PASS |
| `ocfit-nbody-asini-B` | own | 1 | 1.0024 | +0.24 % | 1 % | PASS |
| `ocfit-nbody-P3` | own | 1 | 0.99996 | −0.004 % | 0.5 % | PASS |

**Result (Table 26).** The synthetic O−C of a binary on a circular edge-on outer orbit of 10 au (equal masses)
has the LTTE semi-amplitude of Irwin (1952) to 10⁻⁴; the dynamical delays of this wide system are below
0.01 min. In a 2 + 1 M☉ ZAMS binary with $P$ = 5 d and $i$ = 89°, the secondary minima measured from the
light curve come 27.36 s after half the period, against the Rømer delay of 27.45 s (−0.34 %); the
PHOEBE hybrid gives 27.29 s (−0.60 %). A binary whose period changes by 4·10⁻⁶ d/yr over 60 yr,
integrated with N-body and period kicks, gives the quadratic coefficient of its minima to 5·10⁻⁵; a
stage-2 fit of 160 such minima recovers $\dot P$ within its error (below 5 % of the value). In the stage-1
fit, SWINGERsim's apsidal term computed from exact conjunctions agrees with the nine-term series of
Giménez & Bastero (1995) used by OCFit to 0.6 ms at $e$ = 0.3 (at $e$ = 0.6 the truncated series is off
by 1.6 s), and its LTTE term equals OCFit's. Finally, the eclipse timings of an N-body 2+2 with a
1500-d outer orbit ($e$ = 0.3, $i$ = 60°), fitted in stage 1 with one shared outer orbit, return the
Keplerian $a\sin i$ of both pairs to 0.14 % and 0.24 % and the outer period to 4·10⁻⁵. The remaining
difference is the dynamical delay, of order $(3/8\pi)P_\mathrm{in}^2/P_\mathrm{out}$, which the analytic
stage-1 model does not contain.

## Pulsation timing {#pulsations}

**Definitions.** A pulsating star is a clock. In a binary its pulsations arrive early and late by the
light travel time across its orbit, exactly as eclipses do. SWINGERsim includes the pulsations in its
light curves, with periods from published relations (period–luminosity of δ Scuti stars, pulsation
constant $Q$, Cepheid models), and measures the delays as an observer does: it fits the pulsation
phase in short segments of the light curve (Murphy et al. 2014).

**References.** For a circular orbit the amplitude of the time delays is $a_1\sin i/c$ (Murphy et al.
2014, Sect. 2.7, eq. 10); for an eccentric orbit the ratio of the second to the first harmonic of the delays
depends only on the eccentricity; KIC 9651065, a δ Scuti star in a 272.7-d orbit (Murphy et al. 2014,
Table A1); the period–luminosity relation of 2375 Galactic δ Scuti stars (Poro et al. 2024, eq. 10)
with the fundamental pulsation constant $Q$ = 0.0329 d (Breger 1979); and the period–radius relation of
rotating Geneva Cepheid models (Anderson et al. 2016, Table 5) with Cepheid periods from Bono et al.
(2000).

**Table 27.** Pulsation timing and pulsation periods. Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `pulsation-time-delay-circular` | own | 53.859 s | 53.859 s | +0.00 % | 0.01 % | PASS |
| `murphy2014-kic9651065-a1` | own | 165.57 s | 167.767 s | +1.33 % | 6 % | PASS |
| `murphy2014-kic9651065-a2-a1` | own | 0.231986 | 0.219651 | −0.0123 | 0.025 | PASS |
| `delta-scuti-pl-1p8` | own/tab | 0 mag | −0.0252 mag | −0.0252 mag | 0.25 mag | PASS |
| `delta-scuti-pl-2p0` | own/tab | 0 mag | −0.115 mag | −0.115 mag | 0.25 mag | PASS |
| `cepheid-period-radius-6` | own/tab | 0 dex | −0.0077 dex | −0.0077 dex | 0.05 dex | PASS |
| `cepheid-period-radius-7` | own/tab | 0 dex | +0.0039 dex | +0.0039 dex | 0.05 dex | PASS |

**Result (Table 27).** A δ Scuti-like pulsator (1.7 M☉, pulsation period 0.1 d) with a 0.5 M☉ companion
($P$ = 100 d, $i$ = 60°), sampled like Kepler long-cadence data and analysed in 2-d segments, gives the
delay amplitude $a_1\sin i/c$ = 53.859 s to 3·10⁻⁷ after the segment averaging is corrected by its
sinc factor. For KIC 9651065, with the published orbit ($a_1\sin i$ = 0.37 au, $e$ = 0.468,
$\varpi$ = 2.01 rad), the first harmonic of the delays is 167.8 s against the measured 165.57 ± 2.03 s
(+1.3 %; the tolerance of 6 % is the uncertainty of $a_1\sin i$), and the ratio of the harmonics 0.220
against 0.232 ± 0.013, which checks the eccentric-orbit structure of the delays. MIST stars of 1.8
and 2.0 M☉, where they cross $T_\mathrm{eff}$ = 7500 K on the main sequence, lie 0.03 and 0.12 mag from the
δ Scuti period–luminosity relation (the tolerance of 0.25 mag reflects the uncertainty of $Q$, up to
18 %). MIST Cepheids of 6 and 7 M☉ on their blue loops have radii within 0.008 and 0.004 dex of the
Geneva period–radius relation; the tolerance of 0.05 dex is the difference that Anderson et al. (2016)
find between predicted and empirical radii, and the two chains of models share no code. The light
curves and delays of real pulsators in multiple systems are shown in [Chapter 7](#catalogue-pulsations)
and [Chapter 8](#life-pulsators).
