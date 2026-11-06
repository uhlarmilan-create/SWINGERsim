# Stars: structure, tides, evolution and mass exchange {#stellar}

The stars of a multiple system are not points. They have a size and an internal structure, they are
distorted by tides and rotation, they evolve, lose mass in winds and exchange mass with their
companions. This chapter tests every module of SWINGERsim that describes this, from the structure
constants of a single star to a binary that evolves into contact and merges.

## Structure constants, Roche lobes and rotating stars {#structure}

**Definitions.** The **apsidal-motion constant** $k_2$ measures how strongly a star is distorted by a
tidal or centrifugal force; with the radius $R$ it sets the extra apsidal precession of an eccentric
binary, which scales as $k_2R^5$. The **radius of gyration** $\beta$ (or $r_g^2 = I/MR^2$) sets the
moment of inertia $I$. The **Roche lobe** is the region around a star in a binary within which its own
gravity dominates; its volume-equivalent radius $R_L$ is given to better than 1 % by the formula of
Eggleton (1983). A star rotating rigidly with angular velocity $\Omega$ has, in the Roche model, an
equatorial radius that grows with $\Omega$ and reaches 1.5 times the polar radius at the critical
rotation.

**References.** The classical apsidal constants of polytropes of index $n$ = 1.5 and 3 (Brooker & Olle
1955); the published tables of $\log k_2$ and $\beta$ along stellar tracks of Claret (2004, up to
125 M☉) and Claret (2023); the exact Roche potential, integrated directly by the benchmark, for the
Eggleton formula and for SWINGERsim's surface mesh; the closed forms of the Roche model of a rotating
star; and the rotational correction of $k_2$, $\delta\log k_2 = -0.87\,\lambda_s$ with
$\lambda_s = 2V^2/(3gR)$, that is a factor $10^{-0.87\cdot 2V^2/3gR}$ (Claret 1999).

**Table 14.** Structure constants, Roche lobes and rotating stars. Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `polytrope-k2-n1.5` | own | 0.1433 | 0.143279 | −0.01 % | 0.2 % | PASS |
| `polytrope-k2-n3` | own | 0.0144 | 0.014443 | +0.30 % | 0.5 % | PASS |
| `claret2023-k2-2msun` | own/tab | −2.478 | −2.47794 | +6·10⁻⁵ | 0.002 | PASS |
| `claret2023-beta-1msun` | own/tab | 0.2265 | 0.226499 | −6.3·10⁻⁷ | 0.0005 | PASS |
| `claret2004-log-k2` | own/tab | −2.81523 | −2.81569 | −0.000463 | 0.002 | PASS |
| `claret2004-beta` | own/tab | 0.161706 | 0.161699 | −6.6·10⁻⁶ | 0.001 | PASS |
| `k2-rotation` | own | 1 | 1 | +0.00 % | 10⁻⁷ % | PASS |
| `roche-eggleton-q1` | own | 0.379861 | 0.378921 | −0.25 % | 1 % | PASS |
| `roche-eggleton-q0.1` | own | 0.2054 | 0.206773 | +0.67 % | 1 % | PASS |
| `roche-lobe-volume` | own | 1 | 1.00509 | +0.51 % | 1 % | PASS |
| `rotating-star-critical` | own | 1.5 | 1.5 | −0.00 % | 10⁻⁴ % | PASS |
| `rotating-star-equator` | own | 0 | −2.2·10⁻¹⁶ | −2.2·10⁻¹⁶ | 10⁻⁹ | PASS |

**Result (Table 14).** SWINGERsim integrates the Lane–Emden and Radau equations for its default stellar model;
the apsidal constants of the polytropes agree with Brooker & Olle (1955) to 0.01 % and 0.3 % (the
automatic tests also check the analytic $n$ = 1 values to 10⁻⁸). The interpolation of the Claret
tables in $\log g$ along the tracks reproduces the published rows to their last digit (6·10⁻⁵ in
$\log k_2$, 6·10⁻⁷ in $\beta$). The Eggleton formula, which SWINGERsim uses for all Roche-lobe
filling factors, differs from the exact volume-equivalent radius by −0.25 % at $q$ = 1 and +0.67 % at
$q$ = 0.1, within the 1 % that Eggleton states for all mass ratios. The surface mesh of a
Roche-filling star (used for light curves of distorted stars, [Chapter 4](#meshes)) encloses a volume
whose equivalent radius is 0.51 % above Eggleton's; the residual is the resolution of the mesh.

## Apsidal motion and tidal dissipation {#tides}

**Definitions.** Tides have a conservative part – the static bulge, which adds apsidal precession –
and a dissipative part – the lag of the bulge behind the line joining the stars, which transfers
angular momentum between the orbit and the spins and dissipates energy. SWINGERsim's N-body backend
uses the REBOUNDx effect `tides_spin` (the constant-time-lag equilibrium tide of Eggleton et al. 1998;
Lu et al. 2023) with stellar parameters set by SWINGERsim; the secular backends use SWINGERsim's own
orbit-averaged equations (Hut 1981).

**References.** The apsidal period from the tidal and rotational distortion of both stars (Sterne
1939); the orbit-averaged rates $\dot a$, $\dot e$ and $\dot\Omega$ of the constant-time-lag model
(Hut 1981; Leconte et al. 2010); the circularisation rate of Hut (1981, eq. 10) for the nested secular
backend; and the three contributions to the apsidal motion of the eclipsing pair A of the observed
quadruple CzeV343 (tidal, rotational and relativistic) computed by Pejcha et al. (2022, Table 3).

**Table 15.** Apsidal motion and tidal dissipation. Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `apsidal-sterne-static` | R | 318.227 yr | 318.201 yr | −0.01 % | 0.5 % | PASS |
| `apsidal-sterne-rotating` | R | 291.213 yr | 291.185 yr | −0.01 % | 0.5 % | PASS |
| `tides-leconte-a` | R | 1 | 0.999285 | −0.07 % | 1 % | PASS |
| `tides-leconte-e` | R | 1 | 0.999356 | −0.06 % | 1 % | PASS |
| `tides-leconte-spin` | R | 1 | 1.00113 | +0.11 % | 1 % | PASS |
| `tides-leconte-e-eccentric` | R | 1 | 0.996923 | −0.31 % | 1 % | PASS |
| `secular-nested-tides-hut81` | own | 1 | 1 | −0.00 % | 10⁻⁴ % | PASS |
| `czev343-apsidal-tidal` | own | 4.14·10⁻⁴ rad/d | 4.22272·10⁻⁴ rad/d | +2.00 % | 3 % | PASS |
| `czev343-apsidal-rotational` | own | 6.4·10⁻⁵ rad/d | 6.32533·10⁻⁵ rad/d | −1.17 % | 5 % | PASS |
| `czev343-apsidal-gr` | own | 1.47·10⁻⁵ rad/d | 1.47874·10⁻⁵ rad/d | +0.59 % | 2 % | PASS |
| `czev343-apsidal-period` | R | 33.4 yr | 34.491 yr | +3.27 % | 6 % | PASS |

**Result (Table 15).** An isolated eccentric binary of extended stars ($R/a\approx0.14$, $k_2$ = 0.015,
$e$ = 0.3), first without rotation and then rotating synchronously, precesses in a 60-yr N-body run
with the apsidal period of Sterne's formula to 10⁻⁴. With a time lag, the N-body rates of change of
$a$, $e$ and the spin agree with the orbit-averaged theory to 0.06–0.31 %, both at $e$ = 0.3 and at
$e$ = 0.6 with a star spinning at 1.2 times the mean motion. The rates are taken from the difference
between runs with and without dissipation, which removes the short-period terms of the conservative
tide. These benchmarks also fix a convention: REBOUNDx's `tau` is twice SWINGERsim's time lag;
without SWINGERsim's factor ½ the rates would be twice as large. The nested secular backend reproduces
Hut's circularisation rate exactly (the same equation). For CzeV343 the three contributions agree with
Pejcha et al. (2022) to 0.6–2 %, and the full N-body quadruple with tides and rotation gives an
apsidal period of pair A of 34.5 yr against the observed 33.4 ± 0.2 yr (Pejcha et al. 2022, Sect. 4.3) (+3.3 %; tolerance 6 %, the
uncertainty of the stellar parameters). CzeV343 is studied in full in Part IV.

## Stellar tracks of low- and intermediate-mass stars {#tracks}

**Definitions.** SWINGERsim follows stellar evolution with tabulated tracks. The **SSE tracks** are
the analytic formulae of Hurley, Pols & Tout (2000), evaluated once with the published SSE Fortran
routines on a grid of 320 masses × 5 metallicities and interpolated by SWINGERsim along
"equivalent evolutionary points" (fractions of each evolutionary phase); they end at helium ignition.
The **MIST tracks** (Choi et al. 2016; Dotter 2016) are MESA models converted from the published EEP
archives, from 0.1 to 300 M☉, with and without rotation, up to the end of the asymptotic giant branch
(AGB). The envelope quantities (mass and radius of the convective envelope) follow BSE (Hurley et al.
2002).

**References.** The ZAMS radii of Tout et al. (1996); the present Sun; the helioseismic base of the
solar convection zone, 0.713 $R_\odot$ (Christensen-Dalsgaard et al. 1991); the published MIST EEP ages;
MIST against SSE (two families of models from different codes); the Reimers wind on the red giant
branch; and the radius of the fully convective M dwarf CM Dra A, which standard models underpredict
(Feiden & Chaboyer 2014).

**Table 16.** Stellar tracks of low- and intermediate-mass stars (class own/tab: SWINGERsim's
interpolation of published tracks). Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `sse-zams-radius-1` | own/tab | 0.88825 R☉ | 0.888312 R☉ | +0.01 % | 0.1 % | PASS |
| `sse-zams-radius-2` | own/tab | 1.79434 R☉ | 1.79435 R☉ | +0.00 % | 0.1 % | PASS |
| `sse-sun-radius` | own/tab | 1 R☉ | 0.984843 R☉ | −1.52 % | 5 % | PASS |
| `sse-sun-luminosity` | own/tab | 1 L☉ | 0.95786 L☉ | −4.21 % | 5 % | PASS |
| `sse-sun-convection-zone` | own/tab | 0.713 | 0.733903 | +0.0209 | 0.03 | PASS |
| `mist-grid-tams-age-1` | own/tab | 9919.42 Myr | 9919.42 Myr | +0.00 % | 0.001 % | PASS |
| `mist-interpolation` | own/tab | 1 | 1.01068 | +1.07 % | 2 % | PASS |
| `mist-ms-lifetime-1` | own/tab | 9878 Myr | 11 003.5 Myr | +11.39 % | 15 % | PASS |
| `mist-ms-lifetime-2` | own/tab | 1054 Myr | 1164.17 Myr | +10.45 % | 15 % | PASS |
| `mist-tams-radius-2` | own/tab | 3.913 R☉ | 3.56583 R☉ | −8.87 % | 15 % | PASS |
| `sse-wind-rgb-1msun` | own/tab | 0.04158 M☉ | 0.0421789 M☉ | +1.44 % | 10 % | PASS |
| `cmdra-a-radius-excess-solar` | own/tab | 0.035 | 0.0335989 | −0.0014 | 0.015 | PASS |
| `cmdra-a-radius-excess-lowz` | own/tab | 0.06 | 0.0641756 | +0.00418 | 0.015 | PASS |

**Result (Table 16).** The SSE grid reproduces the ZAMS radii of Tout et al. (1996) to 10⁻⁴, also for 2.5 M☉,
which lies between grid masses. A 1 M☉, $Z$ = 0.02 star at 4.57 Gyr has $R$ = 0.985 R☉ and
$L$ = 0.958 L☉, within the 5 % accuracy that Hurley et al. (2000) quote for their formulae; the base of
its convective envelope lies at 0.734 $R$, 0.021 above the helioseismic value. The MIST grid reproduces
the published EEP ages exactly, and a track interpolated from neighbours at twice the grid spacing
(1.14 M☉ from 1.12 and 1.16 M☉) has a radius at the end of the main sequence within 1.1 % of the
tabulated one. The main-sequence lifetimes of 1 and 2 M☉ are 11 % and 10 % longer in MIST than in SSE,
and the radius of 2 M☉ at the end of the main sequence is 9 % smaller. These are the known differences
between the two families (convective overshooting, solar abundance scale $Z$ = 0.02 against 0.0142);
they set the systematic uncertainty of ages and of the onset of Roche-lobe overflow in SWINGERsim, and
the tolerance of 15 %. The SSE Reimers wind with $\eta$ = 0.1 removes 0.0422 M☉ from a 1 M☉ star up to
helium ignition, against 0.0416 M☉ in the MIST model with the same law. CM Dra A is 3.4 % larger than
the MIST radius of its mass at solar metallicity and 6.4 % at [Fe/H] = −0.3, as Feiden & Chaboyer
(2014) find for standard models (3–4 % and 5–7 %): SWINGERsim reproduces the known radius inflation
problem of M dwarfs rather than hiding it.

## Massive stars and giants to the end of the AGB {#massive-giants}

**References.** For massive stars: leave-one-out interpolation of the MIST grid (a track computed from
its neighbours and compared with the tabulated one); the wind recipe of Vink et al. (2001, eq. 24); the
hydrogen-burning lifetimes of the Geneva models (Ekström et al. 2012), an independent code with and
without rotation. For giants: the white-dwarf initial–final mass relation of Cummings et al. (2018),
measured for 80 white dwarfs in 13 clusters with MIST progenitor ages (scatter 0.06 M☉). The mass of
the helium–carbon core at the end of the MIST track is almost the mass of the future white dwarf.

**Table 17.** Massive stars and the end of the AGB. Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `mist-interpolation-massive-r` | own/tab | 1 | 0.999431 | −0.06 % | 1 % | PASS |
| `mist-interpolation-massive-m` | own/tab | 1 | 1.00107 | +0.11 % | 3 % | PASS |
| `mist-vink2001` | own/tab | 1 | 1.05511 | +5.51 % | 10 % | PASS |
| `mist-geneva-tauh-15-0` | own/tab | 11.015 Myr | 12.3538 Myr | +12.15 % | 15 % | PASS |
| `mist-geneva-tauh-15-4` | own/tab | 13.447 Myr | 12.6825 Myr | −5.69 % | 15 % | PASS |
| `mist-geneva-tauh-60-0` | own/tab | 3.53 Myr | 3.71791 Myr | +5.32 % | 15 % | PASS |
| `mist-geneva-tauh-60-4` | own/tab | 4.465 Myr | 3.94657 Myr | −11.61 % | 15 % | PASS |
| `mist-geneva-tauh-120-0` | own/tab | 2.672 Myr | 2.78283 Myr | +4.15 % | 15 % | PASS |
| `mist-geneva-tauh-120-4` | own/tab | 3.178 Myr | 3.33128 Myr | +4.82 % | 15 % | PASS |
| `agb-final-mass-1p5` | own/tab | 0.609 M☉ | 0.565426 M☉ | −0.0436 M☉ | 0.1 M☉ | PASS |
| `agb-final-mass-2p5` | own/tab | 0.689 M☉ | 0.610396 M☉ | −0.0786 M☉ | 0.1 M☉ | PASS |
| `agb-final-mass-3p0` | own/tab | 0.745 M☉ | 0.653775 M☉ | −0.0912 M☉ | 0.1 M☉ | PASS |
| `agb-final-mass-4p0` | own/tab | 0.899 M☉ | 0.835342 M☉ | −0.0637 M☉ | 0.1 M☉ | PASS |

**Result (Table 17).** The radius of a 60 M☉ star in the middle of its main sequence, interpolated from 55 and
65 M☉, is within 0.06 % of the tabulated one, and the present mass of 100 M☉ at the end of the main
sequence within 0.11 %. The limits of the interpolation lie later in the evolution: from the late main
sequence on, the inflated envelopes of 40–70 M☉ stars are not resolved by the grid, with radius
errors up to 20–40 % and mass errors up to 7 % (for example −24 % in radius for 60 M☉ at the end of the
main sequence without rotation). The MIST mass-loss rate of a 40 M☉ star follows Vink et al. (2001)
with an offset of +5.5 %, which is the metallicity passed to the recipe. The hydrogen-burning lifetimes
of 15, 60 and 120 M☉ agree with the Geneva models to 4–12 %, with and without rotation. The masses at
the end of the main sequence differ more (60 M☉ without rotation: 46.6 against 36.3 M☉; 120 M☉ with
rotation: 25.3 against 34.6 M☉): this spread between model families, not the interpolation, is the
uncertainty of the mass history of very massive stars.

The white-dwarf masses of MIST stars of 1.5–4 M☉ lie 0.04–0.09 M☉ below the relation of Cummings et
al. (2018), the offset the authors themselves note for these models; the tolerance of 0.1 M☉ is that
offset plus the scatter of the relation. The benchmark compares the core, because the MIST tracks stop
with an envelope of up to 0.15–0.2 M☉ that the post-AGB phase would remove. Derived quantities are
checked by the automatic tests: the Sun's convective envelope has 0.034 M☉ and 0.29 $R$, with a
turnover time of 21 d.

## Inertial waves and tidal circularisation {#circularisation}

**Definitions.** In convective envelopes the tidal flow is damped by turbulent viscosity (the
**equilibrium tide**). In rotating stars the tide also excites **inertial waves**, waves restored by
the Coriolis force, which can dissipate much more energy (Ogilvie 2013; Barker 2020, 2022). The
efficiency is expressed as the modified tidal quality factor $Q'$ (larger $Q'$ = weaker dissipation).
The **circularisation period** of a population of binaries of the same age is the period below which
the orbits are circular; observers define it from the final orbits (Meibom & Mathieu 2005).

**References.** The closed formula for the inertial-wave dissipation of a two-layer star (Ogilvie
2013, in the form of Barker 2020, eq. 39); $Q' \approx 10^7$ for 0.3 and 1 M☉ main-sequence stars
rotating in 10 d (Barker 2022, from MESA profiles); circularisation on the pre-main sequence, 7.2–8.5 d
(Zahn & Bouchet 1989); the circularisation period of pre-main-sequence binaries, 7.1 ± 1.2 d, and
that of the 4-Gyr open cluster M67, 12.1 +1.0/−1.5 d (both Meibom & Mathieu 2005, Table 3; quoted by
Mazeh 2008).

**Table 18.** Inertial waves and tidal circularisation. Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `iw-two-layer` | own | 1 | 1.00031 | +0.03 % | 0.5 % | PASS |
| `iw-q-solar` | own | 7 | 7.0122 | +0.0122 | 0.3 | PASS |
| `iw-q-m-dwarf` | own | 7 | 6.94697 | −0.053 | 0.3 | PASS |
| `tides-circularization-pms` | own | 7.85 d | 7.9307 d | +1.03 % | 10 % | PASS |
| `tides-circularization-pms-population` | own | 7.1 d | 6.5831 d | −0.517 d | 1.2 d | PASS |
| `tides-circularization-m67` | own | 12.1 d | 3.7847 d | −68.72 % | 30 % | **EXPLAINED** |
| `tides-circularization-m67-iw` | own | 12.1 d | 7.4572 d | −38.37 % | 30 % | **EXPLAINED** |
| `tides-circularization-m67-final` | own | 12.1 d | 7.5337 d | −4.57 d | 1.2 d | **EXPLAINED** |

**Result (Table 18).** SWINGERsim's full-profile solver of the inertial-wave dissipation reproduces the closed
two-layer formula to 0.03 % for a smoothed two-layer star; this checks the solver and its
normalisation. For 0.3 and 1 M☉ main-sequence stars rotating in 10 d it gives $\log Q'$ = 6.95 and 7.01,
as Barker (2022). Twin 1 M☉ binaries followed on MIST pre-main-sequence tracks from 1 Myr with the
equilibrium tide and inertial waves are circular up to 7.9 d when they reach the main sequence (Zahn
& Bouchet: 7.2–8.5 d), and a population of such binaries at 3 Myr has a circularisation period of
6.6 d in the observers' definition (final period, $e$ < 0.01), against the observed 7.1 ± 1.2 d of
pre-main-sequence binaries.

The three M67 benchmarks are EXPLAINED as an open point of SWINGERsim's tidal model. With the
equilibrium tide of Hurley et al. (2002) alone, twin solar-type binaries started on the main
sequence with $e$ = 0.3 are circular after 4 Gyr only below 3.8 d; with inertial waves below 7.5 d;
and with the whole history from the pre-main sequence 7.5 d in the final period. The
pre-main-sequence phase is reproduced, but the main sequence adds only about 1 d, whereas M67 requires
growth to 12 d. Barker (2022) finds that inertial waves can explain the circularisation periods of
solar-type binaries, with a circularisation period that grows with main-sequence age as observed.
SWINGERsim's estimate of the inertial-wave dissipation is also frequency-averaged, but it uses a
two-layer approximation (Ogilvie 2013) of the SSE and MIST structures and follows the spins of the
stars, whereas Barker computes the dissipation from full MESA profiles and assumes synchronised
spins. The difference
therefore lies in the stellar structure and the spin history; which of the two causes it is not
resolved here. (Resonance locking, which SWINGERsim does not model, acts mainly on the pre-main
sequence, with a negligible main-sequence contribution; Zanazzi & Wu 2021.) The tolerance is not widened; the difference is stated as a limit of SWINGERsim's
tidal model.

## Mass loss by winds {#mass-loss}

**References.** A binary whose star loses mass slowly (compared with the orbital period) keeps the
adiabatic invariant $aM$ (Jeans 1924; Hadjidemetriou 1963). Within one orbit the eccentricity follows
$e = e_0 + \Psi(1-e_0^2)^{3/2}\sin f/(1+e_0\cos f)$, obtained by integrating eq. 4 of Veras et al.
(2011) over the true anomaly $f$, where $\Psi$ is the mass lost per orbit in units of the total mass;
an initially circular orbit reaches $e = 2\Psi$. An instantaneous loss of the fraction $\delta$ of the
mass gives $a'/a = (1-\delta)/(1-2\delta)$ and $e' = \delta/(1-\delta)$ (Blaauw 1961; Hills 1983). In a
triple whose inner primary loses mass, the expanding orbits can cross the stability limit (the
"triple evolution dynamical instability" of Perets & Kratter 2012).

In the N-body backend the masses change through the REBOUNDx operator `modify_mass` on SWINGERsim's
mass-loss schedule (class R); in the secular and evolving backends the orbits are rescaled by
SWINGERsim's own equations (class own).

**Table 19.** Mass loss. Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `ml-adiabatic-invariant` | R | 0 | 3.57·10⁻⁹ | +3.57·10⁻⁹ | 10⁻⁶ | PASS |
| `ml-osculating-e-03` | R | 0 | 0.0126471 | +0.0126 | 0.05 | PASS |
| `ml-osculating-e-07` | R | 0 | 0.0123538 | +0.0124 | 0.05 | PASS |
| `ml-circular-emax` | R | 2 | 2.00265 | +0.13 % | 2 % | PASS |
| `ml-impulsive-a` | R | 1.21429 | 1.21491 | +0.05 % | 1 % | PASS |
| `ml-impulsive-e` | R | 0.176471 | 0.176892 | +0.24 % | 1 % | PASS |
| `ml-zlk-nbody-secular` | own | 1 | 1.0021 | +0.21 % | 1 % | PASS |
| `ml-tedi-onset` | R | 1 | 1.06475 | +6.47 % | 15 % | PASS |
| `wind-transfer-jeans` | own | 0 | 1.69·10⁻¹¹ | +1.69·10⁻¹¹ | 10⁻⁸ | PASS |

**Result (Table 19).** A binary whose star loses 0.3 M☉ slowly keeps $aM$ to 4·10⁻⁹ in the N-body backend; the
secular backends keep it exactly by construction. Within two orbits the eccentricity follows the
formula above to 1.3 % of $\Psi$ for $e_0$ = 0.3 and 0.7, and a circular orbit reaches $e/\Psi$ = 2.003.
Veras et al. print the denominator of their eq. 17 as $1-e_0\cos f$; that form does not follow from
their eq. 4 and disagrees with the N-body runs by about $\Psi$, so we take it for a misprint. Losing
15 % of the mass within 10⁻³ of an orbit gives the instantaneous results to 0.05 % and 0.24 %; the
residual is the finite duration of the loss. While a ZLK perturber loses half its mass, the N-body and
the secular backend agree on the time of the fourth ZLK maximum to 0.2 %. For the masses of Perets &
Kratter (2012) (5 + 2 M☉ inner binary, 2.25 M☉ tertiary) and the inner primary losing 3 M☉ in
8000 yr, the inner eccentricity starts to grow 6.5 % after the adiabatically expanding orbits cross
the MA01 limit; the system had not broken up after a further 2000 outer orbits, in line with the long
survival on the edge found by Toonen et al. (2022). With winds and mass transfer together, a
25 + 18 M☉ binary keeps $a(M_1+M_2)$ to 2·10⁻¹¹ before its transfer begins.

## Mass transfer and rejuvenation {#mass-transfer}

**Definitions.** A star that fills its Roche lobe transfers mass to its companion (**Roche-lobe
overflow**). How the orbit responds depends on how much of the transferred mass the accretor keeps and
how much angular momentum the lost mass carries away. SWINGERsim's mass-transfer module is its own
code: it follows the donor on its track, computes the rate, the response of the orbit and of the
Roche lobe, and the stability of the transfer. A main-sequence accretor that gains mass is
**rejuvenated**: its fresh hydrogen makes it look younger (Hurley et al. 2000).

**References.** The closed forms of the orbit response with constant fractions of the mass lost from
the donor and from the accretor (Soberman, Phinney & van den Heuvel 1997, eqs. 20–28), and for
conservative transfer and isotropic re-emission (Postnov & Yungelson 2014, eqs. 37 and 43); the
response $\zeta_L = d\ln R_L/d\ln M_d$ of the Roche lobe, approximated by $2.13q - 1.67$ (Tout et al.
1997, quoted by Hurley et al. 2002); the **BSE Algol example** of Hurley et al. (2002, Sect. 3.1) –
2.9 + 0.9 M☉, $Z$ = 0.02, started at $P$ = 8 d and $e$ = 0.7, which tides have circularised to
$P$ = 2.7 d when the primary fills its Roche lobe; the benchmarks start from this state – with BSE run independently in COSMIC 4.2.1 (Breivik
et al. 2020); and the age correction of a rejuvenated accretor (Hurley et al. 2000).

**Table 20.** Mass transfer and rejuvenation (all class own). Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `mt-soberman-general` | own | 1.29573 $a/a_0$ | 1.29573 $a/a_0$ | −0.00 % | 10⁻⁵ % | PASS |
| `mt-conservative` | own | 2.25 $a/a_0$ | 2.25 $a/a_0$ | −0.00 % | 10⁻⁵ % | PASS |
| `mt-isotropic-reemission` | own | 2.03003 $a/a_0$ | 2.03003 $a/a_0$ | −0.00 % | 10⁻⁵ % | PASS |
| `mt-zeta-lobe-q2` | own | 2.59 | 2.60707 | +0.0171 | 0.05 | PASS |
| `mt-bse-algol-donor` | own | 0.42 M☉ | 0.404336 M☉ | −0.0157 M☉ | 0.05 M☉ | PASS |
| `mt-cosmic-algol-donor` | own | 0.401 M☉ | 0.404336 M☉ | +0.00334 M☉ | 0.03 M☉ | PASS |
| `mt-bse-algol-accretor` | own | 2.85 M☉ | 2.51105 M☉ | −0.339 M☉ | 0.15 M☉ | **EXPLAINED** |
| `mt-bse-algol-period` | own | 30 d | 21.9301 d | −26.90 % | 15 % | **EXPLAINED** |
| `mt-bse-algol-giant-phase` | own | 30 d | 30.3561 d | +1.19 % | 3 % | PASS |
| `rejuvenation-hurley` | own | 1 | 0.999174 | −0.08 % | 0.2 % | PASS |

**Result (Table 20).** The numerically integrated orbit response agrees with the closed forms to 3·10⁻¹³
(general fractions: 30 % of the transferred mass lost from the donor, 40 % re-emitted from the
accretor, 30 % accreted) and 2·10⁻¹² (conservative, isotropic re-emission). $\zeta_L$ at $q$ = 2 is
2.607 against the fit's 2.59; the fit is accurate to 0.03 for $q$ = 0.5–3. In the BSE Algol example
SWINGERsim, re-modelling the stars as BSE does, ends the first transfer with a donor of 0.404 M☉ (paper
0.42 M☉, BSE in COSMIC 0.401 M☉; both PASS), an accretor of 2.51 M☉ and $P$ = 21.9 d. The paper gives
2.85 M☉ and 30 d. These two benchmarks are EXPLAINED because the published history is not
self-consistent: from 2.7 d to 17 d at the end of the Hertzsprung gap the orbital angular momentum
would have to *grow* by at least 8 % (an upper bound without any loss gives 13.5 d), whereas the
following giant phase, from 17 d to 30 d, is exactly conservative. Started from the published state at
the end of the Hertzsprung gap (0.53 + 2.74 M☉, 17 d), SWINGERsim reproduces this giant phase (30.4 d
against 30 d). The BSE code itself, run in COSMIC 4.2.1 with the same input and the BSE accretion
limit, ends at 3.2 d with an accretor of 1.07 M☉: how much mass a thermal-time-scale accretor keeps
depends on its thermal time scale (COSMIC uses 3.13·10⁷ instead of 10⁷ yr in the corresponding
equation) and on the implementation of BSE. It is the least certain quantity of Algol evolution, and
SWINGERsim lets the user set it (`accretion: fractions`). The age of a 2 M☉ accretor after gaining 1 %
of its mass follows the rejuvenation formula of Hurley et al. (2000) to 0.08 %.

## Common envelope, magnetic braking, the Applegate mechanism and tidal capture {#ce-magnetism}

**Definitions.** If mass transfer becomes dynamically unstable, the envelope of the donor engulfs both
stars (**common envelope**); in the energy formalism the orbital energy released by the spiral-in
unbinds the envelope with an efficiency $\alpha_\mathrm{CE}$. **Magnetic braking** is the loss of spin
angular momentum of a cool star through its magnetised wind. The **Applegate mechanism** explains
cyclic period changes of close binaries by a magnetic cycle that changes the quadrupole moment of one
star; the question is whether the star has enough energy for it. **Tidal capture**: two stars that
pass close enough lose enough orbital energy to tides to become bound.

**References.** The 35-h orbit after the common envelope of a cataclysmic-variable progenitor of
Tout et al. (1997) – 7 + 1 M☉, $P$ = 1200 d, whose primary fills its Roche lobe on the AGB at 5.36 M☉
with a core of 1.014 M☉ (this is not the cataclysmic-variable example of Hurley et al. 2002, Sect. 3.2,
which starts from 6.0 + 1.3 M☉ in a circular 630-d orbit) – and the 0.04-d orbit after the second
common envelope of the BSE Algol example; the impulsive response of an outer orbit to the sudden loss
of the envelope (Hadjidemetriou 1963; Veras et al. 2011); the wind torque of Matt et al. (2015) and the
Skumanich law $\Omega\propto t^{-1/2}$ (Skumanich 1972); the thin-shell energy ratios of HU Aqr, NN Ser and V471 Tau
(Völschow et al. 2016, Table 4, column of the thin-shell model of Tian et al. 2009); and the capture radii of Ray et al. (1987), quoted by Portegies Zwart
& Meinen (1993), and of McMillan et al. (1987).

**Table 21.** Common envelope, magnetism and tidal capture. Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `ce-tout97-period` | own | 35 h | 34.84 h | −0.46 % | 10 % | PASS |
| `ce-bse-algol-second-period` | own | 0.04 d | 0.0365638 d | −8.59 % | 15 % | PASS |
| `ce-impulsive-outer-a` | R | 1 | 1.01274 | +1.27 % | 3 % | PASS |
| `ce-impulsive-outer-e` | R | 0 | 0.00442155 | +0.00442 | 0.02 | PASS |
| `matt15-sun-period` | own | 27.97 d | 27.7635 d | −0.74 % | 2 % | PASS |
| `matt15-unsaturated-powerlaw` | own | 0.5 | 0.491949 | −0.00805 | 0.02 | PASS |
| `applegate-huaqr` | own | 0.1 | 0.105659 | +5.66 % | 10 % | PASS |
| `applegate-nnser` | own | 3.3 | 3.30338 | +0.10 % | 10 % | PASS |
| `applegate-v471tau` | own | 0.014 | 0.0143688 | +2.63 % | 10 % | PASS |
| `tidal-capture-radius-n15-mcmillan` | own | 3.4 | 3.33454 | −1.93 % | 5 % | PASS |
| `tidal-capture-radius-n3` | own | 2.72 | 2.36599 | −13.01 % | 5 % | **EXPLAINED** |
| `tidal-capture-radius-n15` | own | 3.84 | 3.33454 | −13.16 % | 5 % | **EXPLAINED** |

**Result (Table 21).** The energy formalism reproduces the 35-h orbit (34.8 h; the donor fills its lobe in the
orbit widened by its wind) and the 0.04-d orbit (0.037 d with $\alpha_\mathrm{CE}$ = 3). A 1.8 + 0.5
M☉ binary whose giant goes through a common envelope inside a circular 10-yr outer orbit (N-body)
widens the outer orbit 1.3 % more than the impulsive prediction and makes it eccentric to within 0.004
of the prediction; the residual is the finite duration of the ejection. A solar model started at
$P$ = 1 d spins at 27.76 d at the solar age (the solar value $2\pi/\Omega_\odot$ is 27.97 d), and the spin
follows the Skumanich law (ratio 0.492 against 0.5 between 1.14 and 4.55 Gyr). The Applegate energy
ratios agree with Völschow et al. to 0.1–5.7 %.

The two tidal-capture benchmarks are EXPLAINED as a difference of definition. SWINGERsim uses the fits
of Portegies Zwart & Meinen (1993) (as implemented in NBODY6); for two 0.6 M☉ stars meeting at
10 km/s they give capture radii of 2.37 (polytrope $n$ = 3) and 3.33 (n = 1.5) stellar radii, 13 %
less than the values of Ray et al. (1987) quoted in that paper. At the quoted radii, however, the
tidal energy is 5 and 13 times smaller than the orbital energy that must be removed, so the quoted
values cannot belong to the same velocity and stellar structure; they belong to other conventions of
Ray et al. The same fits agree with McMillan et al. (1987), about 3.4 $R$ for $n$ = 1.5, to 1.9 %
(PASS).

## Contact binaries and the dynamical onset of mass transfer {#contact}

**Definitions.** A **contact binary** is a pair whose stars both overfill their Roche lobes and share
a common envelope (W UMa stars). It loses angular momentum (Gazeas & Stępień 2008) and merges when it
reaches the Darwin instability or loses mass through the outer Lagrange point L2 (Rasio 1995). A
**dynamical onset** of mass transfer happens when ZLK cycles drive the eccentricity of an inner binary
so high that a star touches its Roche lobe at periastron; if tides then circularise the orbit while
conserving angular momentum, the circular separation is $a(1-e_c^2)$, and the donor fills
$(1+e_c)$ times less of its lobe than at the contact (Hut 1981).

**Table 22.** Contact binaries and the dynamical onset. Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `contact-jdot` | own | 2.18893 Gyr | 2.18898 Gyr | +0.00 % | 0.1 % | PASS |
| `dynamic-onset-circularisation` | own | 1 | 1 | −0.00 % | 0.01 % | PASS |
| `dynamic-onset-nbody-vs-secular` | own | 1 | 1.0173 | +1.73 % | 3 % | PASS |

**Result (Table 22).** The angular-momentum loss time of a 1.0 + 0.5 M☉ contact binary with $P$ = 0.35 d,
computed in SWINGERsim's units, equals the value computed in cgs units from Gazeas & Stępień (2008)
to 2·10⁻⁵. The filling of the donor after circularisation at $e_c$ = 0.89 times $(1+e_c)$ is 1 to
2·10⁻⁹. The time of the contact driven by ZLK cycles (about 61 yr) differs by 1.7 % between the N-body
and the secular backend; the difference is the short-period terms of the outer orbit
($a_\mathrm{out}/a_\mathrm{in}$ = 15), which the double-averaged equations omit. Both backends place the
onset at the eccentricity at which $a(1-e_c) = R/r_L$, with $r_L$ the Eggleton lobe radius in units of
$a$. The full evolution of a binary into contact and its merger is shown in
[Chapter 8](#life-cycles).
