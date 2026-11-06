\appendix

# All 219 benchmarks {#appendix-benchmarks}

This appendix lists every benchmark registered in SWINGERsim with its result from the final code, by level
of evidence. Each level has two tables with the same rows in the same order: the first says what each
benchmark compares and which code it tests (Tables A1, A3, A5, A7), the second gives its result (Tables A2,
A4, A6, A8); a benchmark is found in both by its Id. Within a level the benchmarks appear in the order in
which the software registers them, so that a row described as "the same for …" follows the row it refers
to. The columns are:

- **Id** – the identifier; `swingersim validate <Id>` runs the benchmark and prints its source (paper and
  equation, figure or table), its inputs and the reason for its tolerance;
- **Code** – the class of code it tests (Table 1): *own* – equations and numerics written for SWINGERsim;
  *own/tab* – SWINGERsim's interpolation of published tables; *R* – REBOUND IAS15 with REBOUNDx effects,
  driven by SWINGERsim's set-up; *ext* – an external library through a SWINGERsim wrapper. A benchmark of
  own code is the stronger evidence, because nothing outside SWINGERsim could make it pass;
- **Quantity** – what is compared, as the benchmark describes it (symbols as in the source code: `e` the
  eccentricity, `i0` an initial inclination, `MA01` the Mardling–Aarseth criterion, `tier A` the spherical
  light-curve model);
- **Expected** – the reference value; **Computed** – the value SWINGERsim computes, with all printed digits;
- **Difference** – SWINGERsim minus expected, relative (%) or absolute, as defined by the benchmark;
- **Tolerance** – the largest difference that counts as agreement, fixed before the comparison;
- **Verdict** – PASS, or **EXPLAINED** (outside the tolerance for a documented reason; Table 48).

To read a row, compare the column *Difference* with the column *Tolerance*: a PASS has an absolute
difference not larger than the tolerance. Where *Expected* is 0 or 1, the benchmark compares a ratio or a
difference with its ideal value. The fast benchmarks come from the last run of `swingersim validate`, the
slow ones from a full run of `swingersim validate --all` on the production machine, and the nine manual
ones (`phoebe-*`, `sbi-*`, `pejcha13-*`, `vyn23-*`) from separate runs; all with the final code. Totals:
93 benchmarks of level A (93 PASS), 60 of level B (55 PASS, 5 EXPLAINED), 31 of level C (29 PASS, 2
EXPLAINED) and 35 of level D (28 PASS, 7 EXPLAINED); 165 of the 219 test code of the classes own and
own/tab.

\footnotesize

## Level A: analytic results {#appendix-level-a}

93 benchmarks: 93 PASS, 0 EXPLAINED.

**Table A1.** Level A: what each benchmark compares and which code it tests (columns *Id*, *Code*, *Quantity*). The results follow in Table A2.

| Id | Code | Quantity |
|---|---|---|
| `gr-mercury` | R | 1PN perihelion advance of Mercury |
| `zlk-emax-65` | R | ZLK maximum eccentricity, i0 = 65 deg, test particle |
| `zlk-emax-80` | R | ZLK maximum eccentricity, i0 = 80 deg, test particle |
| `zlk-imin-65` | R | mutual inclination at maximum eccentricity, i0 = 65 deg |
| `zlk-below-critical` | R | maximum eccentricity below the critical angle, i0 = 30 deg |
| `zlk-kozai-constant` | R | secular drift of sqrt(1-e\^2) cos i during ZLK cycles, i0 = 65 deg |
| `stability-unstable` | R | 2+2 at 0.5 x Mardling-Aarseth limit changes hierarchy (1 = yes) |
| `stability-stable` | R | 2+2 at 2.5 x Mardling-Aarseth limit keeps hierarchy (0 = no change) |
| `quad-tp-emax` | own | secular quadrupole code: ZLK e\_max, test particle, i0 = 65 deg |
| `zlk-period-nbody` | own | ZLK period, N-body relative to the quadrupole secular period (i0 = 65 deg) |
| `megno-regular` | R | MEGNO of a regular (quasi-periodic) hierarchical triple |
| `megno-chaotic` | R | MEGNO > 4 for a 2+2 inside the unstable region (1 = yes) |
| `brouwer-law` | R | growth exponent of the IAS15 energy error (Brouwer's law) |
| `nodal-precession` | R | nodal precession period, circular inner orbit, i = 20 deg |
| `polytrope-k2-n1.5` | own | apsidal-motion constant of a polytrope n = 1.5 |
| `polytrope-k2-n3` | own | apsidal-motion constant of a polytrope n = 3 |
| `roche-eggleton-q1` | own | Roche-lobe radius / separation, q = 1 (Eggleton formula vs exact) |
| `roche-eggleton-q0.1` | own | Roche-lobe radius / separation, q = 0.1 (Eggleton formula vs exact) |
| `apsidal-sterne-static` | R | apsidal period from the tidal distortion, e = 0.3, no rotation |
| `apsidal-sterne-rotating` | R | apsidal period from tidal + rotational distortion, e = 0.3, Omega = n |
| `tides-leconte-a` | R | tidal da/dt, N-body / theory, e = 0.3, non-rotating stars |
| `tides-leconte-e` | R | tidal de/dt, N-body / theory, e = 0.3, non-rotating stars |
| `tides-leconte-spin` | R | tidal spin-up dOmega/dt, N-body / theory, e = 0.6, Omega = 1.2 n |
| `tides-leconte-e-eccentric` | R | tidal de/dt, N-body / theory, e = 0.6, Omega = 1.2 n |
| `secular-kozai-emax` | own | secular backend: ZLK maximum eccentricity, i0 = 65 deg, test particle |
| `secular-nodal-precession` | own | secular backend: nodal precession period, i = 20 deg \[yr\] |
| `secular-gr-apsidal` | own | secular backend: relativistic apsidal period, a = 0.05 au, e = 0.3 |
| `secular-3p1-kozai-emax` | own | 3+1 secular backend: ZLK e\_max of the inner orbit, test particle, i0 = 65 deg |
| `gas-migration-a` | own | gas phase: semi-major axis after 3000 yr with tau\_a = 5000 yr \[au\] |
| `gas-damping-e` | own | gas phase: eccentricity after 3000 yr with tau\_e = 2000 yr |
| `gas-alignment-i` | own | gas phase: inclination to the disc after 3000 yr, i0 = 40 deg, tau\_i = 1500 yr |
| `scatter-hyperbola-deflection` | own | scattering set-up: deflection angle of a 1 + 0.5 M☉ hyperbolic passage, v\_inf = 5 km/s, b = 3 au |
| `escape-theory-energy-slope-L0` | own | statistical escape theory, equal masses, L0 = 0: slope of dP/dE\_B for /E\_B/ >> /E0/ |
| `escape-theory-energy-slope-highL` | own | statistical escape theory, L0 = L\_c: slope of dP/dE\_B |
| `escape-theory-mean-e-L0` | own | statistical escape theory, L0 = 0: mean eccentricity of the binary |
| `escape-theory-mean-e-highL` | own | statistical escape theory, L0 = L\_c: mean eccentricity of the binary |
| `timing-ltte-amplitude` | own | synthetic O-C: light-time semi-amplitude of binary A, circular edge-on outer orbit, a = 10 au, equal masses |
| `phenom-romer-delay` | own | phenomenological minima: delay of the secondary minimum after half the period, 2 + 1 M☉ ZAMS, P = 5 d, circular, i = 89 deg |
| `phoebe-romer-delay` | ext | PHOEBE hybrid: delay of the secondary minimum after half the period, 2 + 1 M☉, P = 5 d, circular (with minus without light time) |
| `pulsation-time-delay-circular` | own | amplitude of the light-travel-time delays of a pulsating star in a circular binary (1.7 + 0.5 M☉, P = 100 d, i = 60 deg), measured from the synthetic light curve like Kepler data (30-min cadence, pulsation phases in 2-d segments) |
| `mt-soberman-general` | own | separation after transfer with constant fractions (donor 3 -> 1.5 M☉, accretor 2 M☉; 30 % donor wind, 40 % isotropic re-emission, 30 % accreted): numerical orbit response against the closed form |
| `mt-conservative` | own | separation after conservative transfer (donor 3 -> 1 M☉, accretor 2 -> 4) |
| `mt-isotropic-reemission` | own | separation after transfer with complete isotropic re-emission (donor 3 -> 1 M☉, accretor 2 M☉) |
| `ce-impulsive-outer-a` | R | outer semi-major axis after a common envelope inside a triple (N-body) relative to before, divided by the impulsive-mass-loss prediction for a circular orbit |
| `ce-impulsive-outer-e` | R | outer eccentricity after the common envelope of the same triple minus the impulsive-mass-loss prediction |
| `matt15-unsaturated-powerlaw` | own | Omega(4.55 Gyr) / Omega(1.14 Gyr) of the braked solar model |
| `galactic-tide-emax-60` | own | largest eccentricity of a 2e4 au orbit inclined by 60 deg to the Galactic disc (secular backend, physics.galactic\_tide) |
| `flyby-first-order-kick` | own | relative difference of the eccentricity-vector change of a binary in a parabolic passage at Q = 20 a: first-order formula against direct N-body |
| `sbi-linear-gaussian-std` | ext | SBI (neural posterior estimation): posterior width of a linear Gaussian model |
| `quint-out-zlk-emax` | R | 2+2+1: ZLK e\_max of the orbit of the two binaries driven by a fifth star, i0 = 65 deg |
| `mist-grid-tams-age-1` | own/tab | MIST grid: age of 1 M☉ at the end of the main sequence (EEP 454), \[Fe/H\] = 0, v/vcrit = 0 |
| `mist-interpolation` | own/tab | MIST: radius of 1.14 M☉ at the TAMS interpolated from 1.12 and 1.16 M☉ / tabulated |
| `claret2023-k2-2msun` | own/tab | Claret (2023) grid: log k2 of 2.0 M☉, \[Fe/H\] = 0, log g = 4.069 |
| `claret2023-beta-1msun` | own/tab | Claret (2023) grid: radius of gyration beta of 1.0 M☉, \[Fe/H\] = 0, log g = 4.216 |
| `sse-zams-radius-1` | own/tab | SSE track grid: ZAMS radius of 1 M☉, Z = 0.02 |
| `sse-zams-radius-2` | own/tab | SSE track grid: ZAMS radius of 2.5 M☉, Z = 0.02 |
| `ml-adiabatic-invariant` | R | max abs(a M / (a M)\_0 - 1) while m1 loses 0.3 M☉ slowly (Psi = 5.6e-6), N-body |
| `ml-osculating-e-03` | R | max abs(e - e\_theory) / Psi over two orbits, e0 = 0.3, Psi = 5.6e-5, N-body |
| `ml-osculating-e-07` | R | max abs(e - e\_theory) / Psi over two orbits, e0 = 0.7, Psi = 5.6e-5, N-body |
| `ml-circular-emax` | R | largest e / Psi of an initially circular binary under slow mass loss |
| `ml-impulsive-a` | R | a'/a after a circular binary loses 15 % of its mass within 1e-3 of a period |
| `ml-impulsive-e` | R | eccentricity after a circular binary loses 15 % of its mass within 1e-3 P |
| `ml-zlk-nbody-secular` | own | time of the 4th ZLK rise (e > 0.6) while the perturber loses half its mass, N-body / secular quadrupole |
| `mist-interpolation-massive-r` | own/tab | MIST: radius of 60 M☉ in the middle of the main sequence (EEP 330) interpolated from 55 and 65 M☉ / tabulated |
| `mist-interpolation-massive-m` | own/tab | MIST: current mass of 100 M☉ at the TAMS interpolated from 95 and 105 M☉ / tabulated |
| `mist-vink2001` | own/tab | MIST wind of 40 M☉ at 30 % of its main sequence / Vink et al. (2001) eq. 24 |
| `claret2004-log-k2` | own/tab | Claret (2004): log k2 of 39.8 M☉ at log g = 3.416 (3.65 Myr) |
| `claret2004-beta` | own/tab | Claret (2004): gyration radius beta of 39.8 M☉ at log g = 3.416 |
| `ld-tlusty-row` | own/tab | TLUSTY limb darkening: a1 of OG40000g400v10 in Bessell V |
| `ocfit-gimenez-bastero` | own | apsidal term of the times of minima, e = 0.3: largest difference between the exact conjunction times and the series of Gimenez & Bastero (1995) used by OCFit |
| `ocfit-lite-amplitude` | own | half peak-to-peak light-time effect over one orbit / K of Irwin (1952), e3 = 0.55 |
| `ocfit-nbody-asini-A` | own | stage-1 fit of N-body eclipse timings of a 2+2 (outer P = 1500 d, e = 0.3, i = 60 deg): fitted / Keplerian a sin i of pair A |
| `ocfit-nbody-asini-B` | own | the same for pair B (anti-phase, shared outer orbit) |
| `ocfit-nbody-P3` | own | the same fit: fitted / configured outer period |
| `wind-transfer-jeans` | own | 25 + 18 M☉ MIST binary with winds and mass transfer: largest relative change of a (M1 + M2) before the transfer starts |
| `secular-nested-tides-hut81` | own | 3+1 in the nested secular backend: tidal de/dt of the inner binary (P = 4 d, e = 0.05, synchronous) / Hut (1981) eq. 10 |
| `galactic-tide-emax-60-nbody` | R | largest eccentricity of a 2e4 au binary inclined by 60 deg to the disc, N-body backend with the Galactic tide (from e = 0.01) |
| `gw-secular-circular-decay` | own | circular 1.4 + 1.4 M☉ binary in a triple, secular backend with physics.gw: a / a0 (1 - t/tau)\^(1/4) at 90 % of the merger time |
| `gw-nbody-decay` | R | the same binary with e = 0.3, N-body backend: change of a over 10 yr / Peters rate x 10 yr |
| `disc-strip-uniform` | own | opaque edge-on slab (half thickness 0.3 R) across a uniform star: transmitted fraction |
| `disc-face-on-tau` | own | face-on disc of optical depth 0.7 in front of a star: transmitted fraction |
| `disc-inclined-tau` | own | disc of optical depth 0.7 inclined 60 deg in front of a star: transmitted fraction |
| `roche-lobe-volume` | own | volume-equivalent radius of the Roche lobe (q = M2/M1 = 0.3) from the surface mesh of the Roche module / Eggleton's formula |
| `rotating-star-critical` | own | equatorial / polar radius of a Roche-model star at critical rotation |
| `rotating-star-equator` | own | R\_eq/R\_pol - (1 + w\^2/2) with w = Omega / Omega\_K(R\_eq), Roche model, Omega / Omega\_crit = 0.6 |
| `mesh-sphere-eclipse` | own | largest difference of the light curves of two spheres through both eclipses: surface meshes (the mesh module) against the exact tier-A eclipses |
| `disc-ring-precession` | own | precession rate of a narrow rigid ring at 4 a\_b around a circular binary / test-particle rate |
| `disc-polar-critical` | own | critical inclination \[deg\] for polar libration of a disc around a binary with e = 0.5, from integrations of SWINGERsim's rigid-disc equation |
| `rejuvenation-hurley` | own | main-sequence age of a 2 M☉ accretor after gaining 1 % of its mass / (t'\_MS / t\_MS)(M / M') t |
| `k2-rotation` | own | rotational correction of k2 (2 M☉, 2 R☉, V = 200 km/s) / 10\^(-0.87 × 2 V\^2 / 3 g R) (Claret 1999) |
| `contact-jdot` | own | angular-momentum loss time \[Gyr\] of a 1.0 + 0.5 M☉ contact binary, P = 0.35 d (SWINGERsim units vs cgs) |
| `pdot-parabola` | R | quadratic coefficient of the eclipse times of a binary with dP/dt = 4e-6 d/yr over 60 yr (N-body with period kicks) / analytic (1/2) P dP/dt |
| `dynamic-onset-circularisation` | own | Roche-lobe filling of the donor after the circularisation of a Kozai-Lidov triple at the contact (e\_c = 0.89), times (1 + e\_c) |

**Table A2.** Level A: results (columns *Expected*, *Computed*, *Difference*, *Tolerance*, *Verdict*), in the same order as Table A1.

| Id | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|
| `gr-mercury` | 42.98 ″/century | 42.9803 ″/century | +0.00 % | 0.5 % | PASS |
| `zlk-emax-65` | 0.838047 | 0.84023 | +0.26 % | 1 % | PASS |
| `zlk-emax-80` | 0.974548 | 0.976155 | +0.16 % | 1 % | PASS |
| `zlk-imin-65` | 39.2315° | 39.3441° | +0.113° | 0.5° | PASS |
| `zlk-below-critical` | 0 | 0.00173037 | +0.00173 | 0.01 | PASS |
| `zlk-kozai-constant` | 0 | 5.16811·10⁻⁵ | +5.17·10⁻⁵ | 0.001 | PASS |
| `stability-unstable` | 1 | 1 | +0 | 0 | PASS |
| `stability-stable` | 0 | 0 | +0 | 0 | PASS |
| `quad-tp-emax` | 0.838047 | 0.838082 | +0.00 % | 0.1 % | PASS |
| `zlk-period-nbody` | 0 | 0.00660493 | +0.0066 | 0.02 | PASS |
| `megno-regular` | 2 | 1.99959 | −0.00041 | 0.05 | PASS |
| `megno-chaotic` | 1 | 1 | +0 | 0 | PASS |
| `brouwer-law` | 0.5 | 0.59301 | +0.093 | 0.15 | PASS |
| `nodal-precession` | 11351.4 yr | 11339.8 yr | −0.10 % | 1 % | PASS |
| `polytrope-k2-n1.5` | 0.1433 | 0.143279 | −0.01 % | 0.2 % | PASS |
| `polytrope-k2-n3` | 0.0144 | 0.014443 | +0.30 % | 0.5 % | PASS |
| `roche-eggleton-q1` | 0.379861 | 0.378921 | −0.25 % | 1 % | PASS |
| `roche-eggleton-q0.1` | 0.2054 | 0.206773 | +0.67 % | 1 % | PASS |
| `apsidal-sterne-static` | 318.227 yr | 318.201 yr | −0.01 % | 0.5 % | PASS |
| `apsidal-sterne-rotating` | 291.213 yr | 291.185 yr | −0.01 % | 0.5 % | PASS |
| `tides-leconte-a` | 1 | 0.999285 | −0.07 % | 1 % | PASS |
| `tides-leconte-e` | 1 | 0.999356 | −0.06 % | 1 % | PASS |
| `tides-leconte-spin` | 1 | 1.00113 | +0.11 % | 1 % | PASS |
| `tides-leconte-e-eccentric` | 1 | 0.996923 | −0.31 % | 1 % | PASS |
| `secular-kozai-emax` | 0.838047 | 0.838082 | +0.00 % | 0.5 % | PASS |
| `secular-nodal-precession` | 11351.4 yr | 11351.4 yr | −0.00 % | 0.5 % | PASS |
| `secular-gr-apsidal` | 7113.77 yr | 7113.77 yr | −0.00 % | 0.5 % | PASS |
| `secular-3p1-kozai-emax` | 0.838047 | 0.838082 | +0.00 % | 0.5 % | PASS |
| `gas-migration-a` | 5.48812 au | 5.48812 au | −0.00 % | 0.0001 % | PASS |
| `gas-damping-e` | 0.111565 | 0.111565 | +0.00 % | 0.0001 % | PASS |
| `gas-alignment-i` | 5.63999° | 5.63943° | −0.000561° | 0.01° | PASS |
| `scatter-hyperbola-deflection` | 173.548° | 173.548° | +1.36·10⁻¹²° | 1·10⁻⁶° | PASS |
| `escape-theory-energy-slope-L0` | −3.5 | −3.5026 | −0.0026 | 0.05 | PASS |
| `escape-theory-energy-slope-highL` | −4 | −4.02479 | −0.0248 | 0.1 | PASS |
| `escape-theory-mean-e-L0` | 0.8133 | 0.812817 | −0.000483 | 0.01 | PASS |
| `escape-theory-mean-e-highL` | 0.7 | 0.710776 | +0.0108 | 0.02 | PASS |
| `timing-ltte-amplitude` | 41.5837 min | 41.5863 min | +0.01 % | 0.1 % | PASS |
| `phenom-romer-delay` | 27.4517 s | 27.3578 s | −0.34 % | 1.5 % | PASS |
| `phoebe-romer-delay` | 27.452 s | 27.287 s | −0.60 % | 2 % | PASS |
| `pulsation-time-delay-circular` | 53.859 s | 53.859 s | +0.00 % | 0.01 % | PASS |
| `mt-soberman-general` | 1.29573 a/a0 | 1.29573 a/a0 | −0.00 % | 1·10⁻⁵ % | PASS |
| `mt-conservative` | 2.25 a/a0 | 2.25 a/a0 | −0.00 % | 1·10⁻⁵ % | PASS |
| `mt-isotropic-reemission` | 2.03003 a/a0 | 2.03003 a/a0 | −0.00 % | 1·10⁻⁵ % | PASS |
| `ce-impulsive-outer-a` | 1 | 1.01274 | +1.27 % | 3 % | PASS |
| `ce-impulsive-outer-e` | 0 | 0.00442155 | +0.00442 | 0.02 | PASS |
| `matt15-unsaturated-powerlaw` | 0.5 | 0.491949 | −0.00805 | 0.02 | PASS |
| `galactic-tide-emax-60` | 0.829156 | 0.829181 | +2.5·10⁻⁵ | 0.01 | PASS |
| `flyby-first-order-kick` | 0 | 0.00594917 | +0.00595 | 0.03 | PASS |
| `sbi-linear-gaussian-std` | 0.3092 | 0.291 | −5.8 % | 10 % | PASS |
| `quint-out-zlk-emax` | 0.81272 | 0.8152 | +0.31 % | 1.0 % | PASS |
| `mist-grid-tams-age-1` | 9919.42 Myr | 9919.42 Myr | +0.00 % | 0.001 % | PASS |
| `mist-interpolation` | 1 | 1.01068 | +1.07 % | 2 % | PASS |
| `claret2023-k2-2msun` | −2.478 | −2.47794 | +6·10⁻⁵ | 0.002 | PASS |
| `claret2023-beta-1msun` | 0.2265 | 0.226499 | −6.34·10⁻⁷ | 0.0005 | PASS |
| `sse-zams-radius-1` | 0.88825 R☉ | 0.888312 R☉ | +0.01 % | 0.1 % | PASS |
| `sse-zams-radius-2` | 1.79434 R☉ | 1.79435 R☉ | +0.00 % | 0.1 % | PASS |
| `ml-adiabatic-invariant` | 0 | 3.56804·10⁻⁹ | +3.57·10⁻⁹ | 1·10⁻⁶ | PASS |
| `ml-osculating-e-03` | 0 | 0.0126471 | +0.0126 | 0.05 | PASS |
| `ml-osculating-e-07` | 0 | 0.0123538 | +0.0124 | 0.05 | PASS |
| `ml-circular-emax` | 2 | 2.00265 | +0.13 % | 2 % | PASS |
| `ml-impulsive-a` | 1.21429 | 1.21491 | +0.05 % | 1 % | PASS |
| `ml-impulsive-e` | 0.176471 | 0.176892 | +0.24 % | 1 % | PASS |
| `ml-zlk-nbody-secular` | 1 | 1.0021 | +0.21 % | 1.0 % | PASS |
| `mist-interpolation-massive-r` | 1 | 0.999431 | −0.06 % | 1 % | PASS |
| `mist-interpolation-massive-m` | 1 | 1.00107 | +0.11 % | 3 % | PASS |
| `mist-vink2001` | 1 | 1.05511 | +5.51 % | 10 % | PASS |
| `claret2004-log-k2` | −2.81523 | −2.81569 | −0.000463 | 0.002 | PASS |
| `claret2004-beta` | 0.161706 | 0.161699 | −6.61·10⁻⁶ | 0.001 | PASS |
| `ld-tlusty-row` | 1.18262 | 1.18262 | +0 | 0.0001 | PASS |
| `ocfit-gimenez-bastero` | 0 s | 0.000629837 s | +0.00063 s | 0.01 s | PASS |
| `ocfit-lite-amplitude` | 1 | 1 | −0.00 % | 0.0001 % | PASS |
| `ocfit-nbody-asini-A` | 1 | 0.99864 | −0.14 % | 1.0 % | PASS |
| `ocfit-nbody-asini-B` | 1 | 1.0024 | +0.24 % | 1.0 % | PASS |
| `ocfit-nbody-P3` | 1 | 0.99996 | −0.00 % | 0.5 % | PASS |
| `wind-transfer-jeans` | 0 | 1.68718·10⁻¹¹ | +1.69·10⁻¹¹ | 1·10⁻⁸ | PASS |
| `secular-nested-tides-hut81` | 1 | 1 | −0.00 % | 0.0001 % | PASS |
| `galactic-tide-emax-60-nbody` | 0.82916 | 0.82908 | −0.01 % | 1.0 % | PASS |
| `gw-secular-circular-decay` | 1 | 1 | −0.00 % | 0.01 % | PASS |
| `gw-nbody-decay` | 1 | 1.0003 | +0.03 % | 1.0 % | PASS |
| `disc-strip-uniform` | 0.623838 | 0.62475 | +0.000912 | 0.005 | PASS |
| `disc-face-on-tau` | 0.496585 | 0.496585 | +0 | 1·10⁻⁹ | PASS |
| `disc-inclined-tau` | 0.246597 | 0.246597 | −4.16·10⁻¹⁶ | 1·10⁻⁹ | PASS |
| `roche-lobe-volume` | 1 | 1.00509 | +0.51 % | 1 % | PASS |
| `rotating-star-critical` | 1.5 | 1.5 | −0.00 % | 0.0001 % | PASS |
| `rotating-star-equator` | 0 | −2.22045·10⁻¹⁶ | −2.22·10⁻¹⁶ | 1·10⁻⁹ | PASS |
| `mesh-sphere-eclipse` | 0 | 0.000259643 | +0.00026 | 0.001 | PASS |
| `disc-ring-precession` | 1 | 1.00018 | +0.02 % | 0.1 % | PASS |
| `disc-polar-critical` | 37.7612° | 37.7612° | +6.29·10⁻⁹° | 0.2° | PASS |
| `rejuvenation-hurley` | 1 | 0.999174 | −0.08 % | 0.2 % | PASS |
| `k2-rotation` | 1 | 1 | +0.00 % | 1·10⁻⁷ % | PASS |
| `contact-jdot` | 2.18893 Gyr | 2.18898 Gyr | +0.00 % | 0.1 % | PASS |
| `pdot-parabola` | 1 | 1.00005 | +0.00 % | 0.1 % | PASS |
| `dynamic-onset-circularisation` | 1 | 1 | −0.00 % | 0.01 % | PASS |

## Level B: published numerical results {#appendix-level-b}

60 benchmarks: 55 PASS, 5 EXPLAINED.

**Table A3.** Level B: what each benchmark compares and which code it tests (columns *Id*, *Code*, *Quantity*). The results follow in Table A4.

| Id | Code | Quantity |
|---|---|---|
| `naoz13-fig4-imin` | R | Naoz et al. 2013 Fig. 4 (PSR B1620-26 setup): minimum mutual inclination |
| `naoz13-fig4-imax` | R | Naoz et al. 2013 Fig. 4: maximum mutual inclination (orbit flip to retrograde) |
| `naoz13-fig3-flip-time` | R | Naoz et al. 2013 Fig. 3 (eccentric ZLK): time of the first flip through 90 deg |
| `naoz13-fig3-imax` | R | Naoz et al. 2013 Fig. 3: maximum inclination after the flip |
| `quad-naoz13-fig4-imin` | own | secular quadrupole code: Naoz et al. 2013 Fig. 4, minimum inclination |
| `quad-naoz13-fig4-imax` | own | secular quadrupole code: Naoz et al. 2013 Fig. 4, maximum inclination |
| `vok16-v994-JA` | own | Vokrouhlicky 2016, V994 Her setup: initial tilt J\_A of binary A to the outer orbit |
| `vok16-v994-JB` | own | Vokrouhlicky 2016, V994 Her setup: initial tilt J\_B of binary B to the outer orbit |
| `vok16-jab-max-analytic` | own | Vokrouhlický 2016, analytic limit J\_AB\^max for V994 Her (Fig. 2 setup) |
| `vok16-fig2-jab-min` | own | secular backend: Vokrouhlický 2016 Fig. 2, minimum J\_AB over 100 kyr |
| `vok16-fig2-jab-max` | own | secular backend: Vokrouhlický 2016 Fig. 2, maximum J\_AB over 100 kyr |
| `vok16-fig2-jab-1myr` | own | secular backend: Fig. 2 setup, maximum J\_AB over 1 Myr |
| `vok16-fig3-jab-max` | own | secular backend: Vokrouhlický 2016 Fig. 3, maximum J\_AB over 20 kyr |
| `vok16-eclipse-prob-A` | own | V994 Her: chance that binary A eclipses, isotropic orientation |
| `vok16-eclipse-prob-B` | own | V994 Her: chance that binary B eclipses, isotropic orientation |
| `secular-naoz13-fig3-flip-time` | own | 3+1 secular backend: Naoz et al. 2013 Fig. 3, time of the first flip \[yr\] |
| `secular-naoz13-fig3-imax` | own | 3+1 secular backend: Naoz et al. 2013 Fig. 3, maximum inclination after the flip |
| `at16-binbin-single-ionization` | R | binary-binary scattering, equal masses, a1 = a2, e = 0, v\_inf = v\_crit: cross-section of 'single ionization' / (pi a\^2) |
| `at16-binbin-exchange-ionization` | R | binary-binary scattering, equal masses, a1 = a2, e = 0, v\_inf = v\_crit: cross-section of 'exchange + ionization' / (pi a\^2) |
| `at16-binbin-exchange` | R | binary-binary scattering, equal masses, a1 = a2, e = 0, v\_inf = v\_crit: cross-section of 'exchange' / (pi a\^2) |
| `chaosmap-circumbinary-limit` | R | MEGNO map: inner stability limit of a circumbinary orbit, equal-mass circular binary, all of 8 phases regular \[a\_bin\] |
| `qzcar-broz-a1` | own | QZ Car: semi-major axis of Ac1+Ac2 (P = 6.0 d) from the Table 2 elements and masses |
| `qzcar-broz-a2` | own | QZ Car: semi-major axis of Aa1+Aa2 (P = 20.7 d) from the Table 2 elements and masses |
| `qzcar-broz-a3` | own | QZ Car: semi-major axis of Ac+Aa (P = 14 728 d) from the Table 2 elements and masses |
| `mt-zeta-lobe-q2` | own | response zeta\_L = dln R\_L / dln M\_d of the donor's Roche lobe (Eggleton 1983) to conservative transfer at q = 2 |
| `ce-tout97-period` | own | orbital period after the common envelope of the CV progenitor of Tout et al. (1997): 5.36 M☉ AGB star (1.014 M☉ core) with a 1 M☉ companion, alpha\_CE = 1 |
| `ce-bse-algol-second-period` | own | orbital period after the second common envelope of the BSE Algol example (2.85 M☉ giant with a 0.416 M☉ core, 0.417 M☉ white dwarf, P = 18 d, alpha\_CE = 3) |
| `matt15-sun-period` | own | rotation period of a solar model at 4.55 Gyr with the wind torque of Matt et al. (2015), started at P = 1 d |
| `applegate-huaqr` | own | energy needed by the Applegate mechanism over the energy radiated in one modulation period, HU Aqr (thin-shell model) |
| `applegate-nnser` | own | energy needed by the Applegate mechanism over the energy radiated in one modulation period, NN Ser (thin-shell model) |
| `applegate-v471tau` | own | energy needed by the Applegate mechanism over the energy radiated in one modulation period, V471 Tau (thin-shell model) |
| `vyn23-2p2-ma01-score` | own | overall score of the MA01 criterion against the N-body ghost test on SWINGERsim's 2+2 training set (sampled as in Vynatheya et al. 2023, Sect. 3.1) |
| `vyn23-2p2-mlp-score` | ext | held-out score of SWINGERsim's MLP stability emulator (2+2; 11 parameters, 4 x 50 neurons) against the N-body ghost test |
| `vyn23-3p1-ma01-score` | own | overall score of the MA01 criterion against the N-body ghost test on SWINGERsim's 3+1 training set (sampled as in Vynatheya et al. 2023, Sect. 3.1) |
| `vyn23-3p1-mlp-score` | ext | held-out score of SWINGERsim's MLP stability emulator (3+1; 11 parameters, 4 x 50 neurons) against the N-body ghost test |
| `escape-theory-vs-nbody-light-escaper` | own | resonant binary-single encounters (1 + 0.6 M☉ binary, 0.8 M☉ star): N-body fraction in which the 0.6 M☉ star escapes vs the statistical theory (expected value: SWINGERsim's own evaluation of the theory of Stone & Leigh 2019, not a published number) |
| `czev343-apsidal-tidal` | own | CzeV343 binary A: tidal apsidal rate \[rad/day\] |
| `czev343-apsidal-rotational` | own | CzeV343 binary A: rotational apsidal rate, pseudo-synchronous spins \[rad/day\] |
| `czev343-apsidal-gr` | own | CzeV343 binary A: relativistic apsidal rate \[rad/day\] |
| `vok16-pejcha-fA` | own | Pejcha-like 2+2 population, secular: binary A reaches 1 - e < 1e-3 |
| `vok16-pejcha-fB` | own | Pejcha-like 2+2 population, secular: binary B reaches 1 - e < 1e-3 |
| `pejcha13-F6-fA` | R | Pejcha et al. F6 population, N-body: binary A reaches 1 - e < 1e-3 |
| `pejcha13-F6-fB` | R | Pejcha et al. F6 population, N-body: binary B reaches 1 - e < 1e-3 |
| `mist-ms-lifetime-1` | own/tab | main-sequence lifetime of 1 M☉ (ZAMS to TAMS), SSE vs MIST |
| `mist-ms-lifetime-2` | own/tab | main-sequence lifetime of 2 M☉ (ZAMS to TAMS), SSE vs MIST |
| `mist-tams-radius-2` | own/tab | radius of 2 M☉ at the end of the main sequence, SSE vs MIST |
| `iw-two-layer` | own | inertial-wave dissipation of a two-layer star (alpha 0.6, beta 0.95): impulsive formalism / closed formula |
| `tides-circularization-pms` | own | period below which twin 1 M☉ binaries (e0 = 0.3) are circular when they reach the main sequence (MIST pre-main sequence from 1 Myr to 50 Myr; equilibrium tide and inertial waves) |
| `ml-tedi-onset` | R | onset of the inner instability (e\_in > 0.1) / crossing of the MA01 limit, inner primary losing 3 of 5 M☉ in 8000 yr |
| `mist-geneva-tauh-15-0` | own/tab | hydrogen-burning lifetime of 15 M☉, non-rotating, MIST vs Geneva |
| `mist-geneva-tauh-15-4` | own/tab | hydrogen-burning lifetime of 15 M☉, rotating, MIST vs Geneva |
| `mist-geneva-tauh-60-0` | own/tab | hydrogen-burning lifetime of 60 M☉, non-rotating, MIST vs Geneva |
| `mist-geneva-tauh-60-4` | own/tab | hydrogen-burning lifetime of 60 M☉, rotating, MIST vs Geneva |
| `mist-geneva-tauh-120-0` | own/tab | hydrogen-burning lifetime of 120 M☉, non-rotating, MIST vs Geneva |
| `mist-geneva-tauh-120-4` | own/tab | hydrogen-burning lifetime of 120 M☉, rotating, MIST vs Geneva |
| `flyby-nbody-backend` | R | parabolic passage at Q = 20 a integrated by the N-body backend (physics.encounters): relative difference of the eccentricity kick from the first-order result of Hamers & Samsing (2019) |
| `tidal-capture-radius-n3` | own | largest periastron distance for tidal capture / R, two n = 3 stars (0.6 M☉, 0.65 R☉) at 10 km/s |
| `tidal-capture-radius-n15` | own | the same for two n = 1.5 stars |
| `tidal-capture-radius-n15-mcmillan` | own | largest periastron distance for tidal capture / R, two equal n = 1.5 stars at 10 km/s |
| `mt-bse-algol-giant-phase` | own | period at the end of the BSE Algol example's transfer, computed from its state at the end of the Hertzsprung gap (0.53 + 2.74 M☉, 17 d; conservative) |

**Table A4.** Level B: results (columns *Expected*, *Computed*, *Difference*, *Tolerance*, *Verdict*), in the same order as Table A3.

| Id | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|
| `naoz13-fig4-imin` | 57.5° | 57.2391° | −0.261° | 0.5° | PASS |
| `naoz13-fig4-imax` | 106.7° | 106.93° | +0.23° | 0.5° | PASS |
| `naoz13-fig3-flip-time` | 3.8·10⁶ yr | 3.9868·10⁶ yr | +4.92 % | 10.0 % | PASS |
| `naoz13-fig3-imax` | 142.6° | 143.01° | +0.408° | 2° | PASS |
| `quad-naoz13-fig4-imin` | 57.5° | 57.5376° | +0.0376° | 0.1° | PASS |
| `quad-naoz13-fig4-imax` | 106.7° | 106.721° | +0.0207° | 0.1° | PASS |
| `vok16-v994-JA` | 11.6° | 11.6462° | +0.0462° | 0.1° | PASS |
| `vok16-v994-JB` | 15.8° | 15.7932° | −0.00678° | 0.1° | PASS |
| `vok16-jab-max-analytic` | 5.2° | 5.20074° | +0.000736° | 0.05° | PASS |
| `vok16-fig2-jab-min` | 3.5° | 3.47152° | −0.0285° | 0.15° | PASS |
| `vok16-fig2-jab-max` | 5.2° | 5.20069° | +0.000687° | 0.15° | PASS |
| `vok16-fig2-jab-1myr` | 5.2° | 5.2006° | +0.000633° | 0.15° | PASS |
| `vok16-fig3-jab-max` | 6° | 6.20874° | +0.209° | 0.5° | PASS |
| `vok16-eclipse-prob-A` | 0.366 | 0.369265 | +0.89 % | 2 % | PASS |
| `vok16-eclipse-prob-B` | 0.397 | 0.399075 | +0.52 % | 2 % | PASS |
| `secular-naoz13-fig3-flip-time` | 3.8·10⁶ yr | 4.04981·10⁶ yr | +6.57 % | 10 % | PASS |
| `secular-naoz13-fig3-imax` | 142.6° | 142.556° | −0.0435° | 2° | PASS |
| `at16-binbin-single-ionization` | 2.645 | 2.6873 | +1.60 % | 10.0 % | PASS |
| `at16-binbin-exchange-ionization` | 3.259 | 3.2305 | −0.88 % | 10.0 % | PASS |
| `at16-binbin-exchange` | 0.302 | 0.37165 | +23.06 % | 30.0 % | PASS |
| `chaosmap-circumbinary-limit` | 2.27 | 2.3 | +0.03 | 0.06 | PASS |
| `qzcar-broz-a1` | 54.032 R☉ | 53.9946 R☉ | −0.07 % | 0.2 % | PASS |
| `qzcar-broz-a2` | 136.132 R☉ | 136.081 R☉ | −0.04 % | 0.2 % | PASS |
| `qzcar-broz-a3` | 13053 R☉ | 13046 R☉ | −0.05 % | 0.2 % | PASS |
| `mt-zeta-lobe-q2` | 2.59 | 2.60707 | +0.0171 | 0.05 | PASS |
| `ce-tout97-period` | 35 h | 34.84 h | −0.46 % | 10 % | PASS |
| `ce-bse-algol-second-period` | 0.04 d | 0.0365638 d | −8.59 % | 15 % | PASS |
| `matt15-sun-period` | 27.97 d | 27.7635 d | −0.74 % | 2 % | PASS |
| `applegate-huaqr` | 0.1 | 0.105659 | +5.66 % | 10 % | PASS |
| `applegate-nnser` | 3.3 | 3.30338 | +0.10 % | 10 % | PASS |
| `applegate-v471tau` | 0.014 | 0.0143688 | +2.63 % | 10 % | PASS |
| `vyn23-2p2-ma01-score` | 0.83 | 0.869 | +0.039 | 0.05 | PASS |
| `vyn23-2p2-mlp-score` | 0.94 | 0.9375 | −0.0025 | 0.03 | PASS |
| `vyn23-3p1-ma01-score` | 0.56 | 0.7896 | +0.23 | 0.05 | **EXPLAINED** |
| `vyn23-3p1-mlp-score` | 0.93 | 0.9195 | −0.0105 | 0.03 | PASS |
| `escape-theory-vs-nbody-light-escaper` | 0.769 | 0.75369 | −0.0153 | 0.07 | PASS |
| `czev343-apsidal-tidal` | 0.000414 rad/d | 0.000422272 rad/d | +2.00 % | 3 % | PASS |
| `czev343-apsidal-rotational` | 6.4·10⁻⁵ rad/d | 6.32533·10⁻⁵ rad/d | −1.17 % | 5 % | PASS |
| `czev343-apsidal-gr` | 1.47·10⁻⁵ rad/d | 1.47874·10⁻⁵ rad/d | +0.59 % | 2 % | PASS |
| `vok16-pejcha-fA` | 0.31 | 0.21 | −0.1 | 0.07 | **EXPLAINED** |
| `vok16-pejcha-fB` | 0.481 | 0.4025 | −0.0785 | 0.07 | **EXPLAINED** |
| `pejcha13-F6-fA` | 0.358 | 0.24 | −0.118 | 0.14 | PASS |
| `pejcha13-F6-fB` | 0.530 | 0.41 | −0.12 | 0.15 | PASS |
| `mist-ms-lifetime-1` | 9878 Myr | 11003.5 Myr | +11.39 % | 15 % | PASS |
| `mist-ms-lifetime-2` | 1054 Myr | 1164.17 Myr | +10.45 % | 15 % | PASS |
| `mist-tams-radius-2` | 3.913 R☉ | 3.56583 R☉ | −8.87 % | 15 % | PASS |
| `iw-two-layer` | 1 | 1.00031 | +0.03 % | 0.5 % | PASS |
| `tides-circularization-pms` | 7.85 d | 7.9307 d | +1.03 % | 10 % | PASS |
| `ml-tedi-onset` | 1 | 1.06475 | +6.47 % | 15 % | PASS |
| `mist-geneva-tauh-15-0` | 11.015 Myr | 12.3538 Myr | +12.15 % | 15 % | PASS |
| `mist-geneva-tauh-15-4` | 13.447 Myr | 12.6825 Myr | −5.69 % | 15 % | PASS |
| `mist-geneva-tauh-60-0` | 3.53 Myr | 3.71791 Myr | +5.32 % | 15 % | PASS |
| `mist-geneva-tauh-60-4` | 4.465 Myr | 3.94657 Myr | −11.61 % | 15 % | PASS |
| `mist-geneva-tauh-120-0` | 2.672 Myr | 2.78283 Myr | +4.15 % | 15 % | PASS |
| `mist-geneva-tauh-120-4` | 3.178 Myr | 3.33128 Myr | +4.82 % | 15 % | PASS |
| `flyby-nbody-backend` | 0 | 0.00738403 | +0.00738 | 0.05 | PASS |
| `tidal-capture-radius-n3` | 2.72 | 2.36599 | −13.01 % | 5 % | **EXPLAINED** |
| `tidal-capture-radius-n15` | 3.84 | 3.33454 | −13.16 % | 5 % | **EXPLAINED** |
| `tidal-capture-radius-n15-mcmillan` | 3.4 | 3.33454 | −1.93 % | 5 % | PASS |
| `mt-bse-algol-giant-phase` | 30 d | 30.3561 d | +1.19 % | 3 % | PASS |

## Level C: other codes {#appendix-level-c}

31 benchmarks: 29 PASS, 2 EXPLAINED.

**Table A5.** Level C: what each benchmark compares and which code it tests (columns *Id*, *Code*, *Quantity*). The results follow in Table A6.

| Id | Code | Quantity |
|---|---|---|
| `sextuple-secular-emax` | own | 2+2+2 sextuple: maximum eccentricity of A, secular backend vs N-body |
| `sextuple-secular-zlk` | own | 2+2+2 sextuple: ZLK period of A, secular backend vs N-body |
| `sextuple-nbody-emax` | R | 2+2+2 sextuple: maximum eccentricity of A, N-body (reproducibility) |
| `mse-3p1-inner-emax` | own | 3+1 (a = 1, 15, 300 au, i\_mid = 85 deg): maximum e of the inner binary in 1e5 yr |
| `mse-3p1-inner-t-emax` | own | 3+1: time of the maximum e of the inner binary |
| `mse-3p1-mid-emax` | own | 3+1: maximum e of the middle orbit |
| `mse-2p2-A-emax` | own | 2+2 (a = 1, 1.5, 40 au): maximum e of binary A in 8e4 yr |
| `mse-2p2-B-emax` | own | 2+2: maximum e of binary B in 8e4 yr |
| `mse-3p1-rlof-time` | own | 3+1 (a = 0.3, 5, 200 au, i = 89 deg): onset of Roche-lobe overflow of the inner binary driven by ZLK cycles |
| `phoebe-sphere-vs-tier-a` | own | largest difference of the normalised light curves, tier A vs PHOEBE 2.5 with the same physics (spheres, black bodies, Claret limb darkening), e = 0.1 |
| `qzcar-broz-lv-ac1` | own | QZ Car: share of Ac1 in the V-band light from its radius and Teff |
| `qzcar-broz-lv-ac2` | own | QZ Car: share of Ac2 in the V-band light from its radius and Teff |
| `qzcar-broz-lv-aa1` | own | QZ Car: share of Aa1 in the V-band light from its radius and Teff |
| `qzcar-broz-lv-aa2` | own | QZ Car: share of Aa2 in the V-band light from its radius and Teff |
| `cepheid-period-radius-6` | own/tab | log R of a 6 M☉ MIST Cepheid (blue loop, middle of the strip, period from Bono et al. 2000) minus the period-radius relation of rotating Geneva Cepheid models |
| `cepheid-period-radius-7` | own/tab | as cepheid-period-radius-6 for a 7 M☉ Cepheid |
| `mt-bse-algol-donor` | own | donor mass at the end of the first mass transfer of the BSE Algol example (2.9 + 0.9 M☉, Z = 0.02, started at P = 8 d, e = 0.7; from the onset of overflow, circular, P = 2.7 d) |
| `mt-bse-algol-accretor` | own | accretor mass at the end of the first mass transfer of the BSE Algol example |
| `mt-bse-algol-period` | own | orbital period at the end of the first mass transfer of the BSE Algol example |
| `mt-cosmic-algol-donor` | own | donor mass when the first mass transfer of the BSE Algol example ends, against the BSE code in COSMIC 4.2.1 (2.9 + 0.9 M☉, circular 2.7 d, Z = 0.02, acc\_lim = -1) |
| `lc-occultation-small` | own | flux of a limb-darkened star behind a dark disk, rp = 0.1, central |
| `lc-occultation-limb` | own | flux of a limb-darkened star behind a dark disk, rp = 0.3 at the limb (z = 0.9) |
| `lc-occultation-large` | own | flux of a limb-darkened star behind a larger dark disk, rp = 1.2, z = 0.5 |
| `lc-detached-primary` | own | primary eclipse depth of a detached binary of spheres |
| `lc-detached-secondary` | own | secondary eclipse depth of a detached binary of spheres |
| `lc-ellipsoidal` | own | semi-amplitude of the ellipsoidal variation, R/a = 0.15, q = 0.8, i = 60 deg |
| `ld-tlusty-vs-atlas` | own/tab | disk-averaged intensity <I>/I(1) in V at 25 000 K, log g 4: TLUSTY non-LTE / ATLAS LTE |
| `sse-wind-rgb-1msun` | own/tab | mass lost by a 1 M☉ star (Z = 0.02) on the red giant branch with the Reimers wind, eta = 0.1: SSE track with the SSE winds against the MIST v1.2 (MESA) model at the tip of the giant branch |
| `contact-vs-phoebe` | own | V amplitude of a contact binary (q = 0.42, fill-out 0.3, i = 80 deg): SWINGERsim mesh minus PHOEBE 2.5 with the same physics |
| `rv-light-vs-phoebe` | own | centre-of-light minus centre-of-mass velocity of the secondary of a contact binary (q = 0.42, fill-out 0.3, i = 80 deg): largest difference SWINGERsim - PHOEBE / its amplitude (52 km/s) |
| `dynamic-onset-nbody-vs-secular` | own | time of the Roche-lobe contact driven by Kozai-Lidov cycles (first maximum, \~61 yr): N-body / secular (double-averaged, octupole) |

**Table A6.** Level C: results (columns *Expected*, *Computed*, *Difference*, *Tolerance*, *Verdict*), in the same order as Table A5.

| Id | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|
| `sextuple-secular-emax` | 0.75207 | 0.751821 | −0.03 % | 2 % | PASS |
| `sextuple-secular-zlk` | 70.2 yr | 71.04 yr | +1.20 % | 5 % | PASS |
| `sextuple-nbody-emax` | 0.75207 | 0.75207 | +0.00 % | 0.2 % | PASS |
| `mse-3p1-inner-emax` | 0.9707 | 0.970714 | +1.37·10⁻⁵ | 0.002 | PASS |
| `mse-3p1-inner-t-emax` | 49325 yr | 49322.8 yr | −0.00 % | 0.2 % | PASS |
| `mse-3p1-mid-emax` | 0.2245 | 0.224461 | −3.86·10⁻⁵ | 0.002 | PASS |
| `mse-2p2-A-emax` | 0.9494 | 0.949492 | +9.2·10⁻⁵ | 0.002 | PASS |
| `mse-2p2-B-emax` | 0.7758 | 0.776227 | +0.000427 | 0.002 | PASS |
| `mse-3p1-rlof-time` | 380.03 yr | 380.1 yr | +0.02 % | 0.2 % | PASS |
| `phoebe-sphere-vs-tier-a` | 0 | 0.00037 | +0.00037 | 0.001 | PASS |
| `qzcar-broz-lv-ac1` | 0.22 | 0.220775 | +0.000775 | 0.005 | PASS |
| `qzcar-broz-lv-ac2` | 0.133 | 0.133856 | +0.000856 | 0.005 | PASS |
| `qzcar-broz-lv-aa1` | 0.635 | 0.635024 | +2.36·10⁻⁵ | 0.005 | PASS |
| `qzcar-broz-lv-aa2` | 0.01 | 0.0103458 | +0.000346 | 0.005 | PASS |
| `cepheid-period-radius-6` | 0 dex | −0.00765757 dex | −0.00766 dex | 0.05 dex | PASS |
| `cepheid-period-radius-7` | 0 dex | 0.00392438 dex | +0.00392 dex | 0.05 dex | PASS |
| `mt-bse-algol-donor` | 0.42 M☉ | 0.404336 M☉ | −0.0157 M☉ | 0.05 M☉ | PASS |
| `mt-bse-algol-accretor` | 2.85 M☉ | 2.51105 M☉ | −0.339 M☉ | 0.15 M☉ | **EXPLAINED** |
| `mt-bse-algol-period` | 30 d | 21.9301 d | −26.90 % | 15 % | **EXPLAINED** |
| `mt-cosmic-algol-donor` | 0.401 M☉ | 0.404336 M☉ | +0.00334 M☉ | 0.03 M☉ | PASS |
| `lc-occultation-small` | 0.987891 | 0.987891 | −1.92·10⁻¹² | 1·10⁻⁷ | PASS |
| `lc-occultation-limb` | 0.945236 | 0.945236 | +3.66·10⁻⁹ | 1·10⁻⁷ | PASS |
| `lc-occultation-large` | 0.128979 | 0.128979 | −5.52·10⁻⁹ | 1·10⁻⁷ | PASS |
| `lc-detached-primary` | 0.380821 | 0.380868 | +4.72·10⁻⁵ | 0.0001 | PASS |
| `lc-detached-secondary` | 0.144658 | 0.144654 | −3.85·10⁻⁶ | 0.0001 | PASS |
| `lc-ellipsoidal` | 0.0027949 | 0.00260619 | −6.75 % | 10 % | PASS |
| `ld-tlusty-vs-atlas` | 1 | 0.997587 | −0.24 % | 2 % | PASS |
| `sse-wind-rgb-1msun` | 0.04158 M☉ | 0.0421789 M☉ | +1.44 % | 10 % | PASS |
| `contact-vs-phoebe` | 0 mag | 0.0030377 mag | +0.00304 mag | 0.03 mag | PASS |
| `rv-light-vs-phoebe` | 0 | 0.022028 | +0.022 | 0.05 | PASS |
| `dynamic-onset-nbody-vs-secular` | 1 | 1.0173 | +1.73 % | 3.0 % | PASS |

## Level D: observed systems {#appendix-level-d}

35 benchmarks: 28 PASS, 7 EXPLAINED.

**Table A7.** Level D: what each benchmark compares and which code it tests (columns *Id*, *Code*, *Quantity*). The results follow in Table A8.

| Id | Code | Quantity |
|---|---|---|
| `v994s-dyn-A` | R | V994 Her sextuple: share of the apsidal motion of A caused by B (and C) |
| `v994s-dyn-B` | R | V994 Her sextuple: share of the apsidal motion of B caused by A (and C) |
| `v994s-dyn-C` | R | V994 Her sextuple: share of the apsidal motion of C caused by A and B |
| `v994s-gr-A` | R | V994 Her: relativistic share of the apsidal motion of A |
| `v994s-gr-B` | R | V994 Her: relativistic share of the apsidal motion of B |
| `v994s-gr-C` | R | V994 Her: relativistic share of the apsidal motion of C |
| `qzcar-evolutionary-mass-ac1` | own/tab | QZ Car: present mass of the single-star MIST model with the radius and Teff of Ac1, against its dynamical mass |
| `qzcar-evolutionary-mass-ac2` | own/tab | the same for Ac2 |
| `qzcar-evolutionary-mass-aa1` | own/tab | the same for Aa1 |
| `agb-final-mass-1p5` | own/tab | white-dwarf mass (core at the end of the AGB) of a 1.5 M☉ MIST star (solar \[Fe/H\]) against the white-dwarf initial-final mass relation |
| `agb-final-mass-2p5` | own/tab | white-dwarf mass (core at the end of the AGB) of a 2.5 M☉ MIST star (solar \[Fe/H\]) against the white-dwarf initial-final mass relation |
| `agb-final-mass-3p0` | own/tab | white-dwarf mass (core at the end of the AGB) of a 3.0 M☉ MIST star (solar \[Fe/H\]) against the white-dwarf initial-final mass relation |
| `agb-final-mass-4p0` | own/tab | white-dwarf mass (core at the end of the AGB) of a 4.0 M☉ MIST star (solar \[Fe/H\]) against the white-dwarf initial-final mass relation |
| `murphy2014-kic9651065-a1` | own | first orbital harmonic A1 of the time delays of KIC 9651065 (delta Sct in a 272.7-d eccentric orbit) from its orbit (a1 sin i = 0.37 au, e = 0.468, varpi = 2.01) |
| `murphy2014-kic9651065-a2-a1` | own | ratio of the second to the first orbital harmonic of the time delays of KIC 9651065 (sets the eccentricity; e = 0.468, varpi = 2.01 rad) |
| `delta-scuti-pl-1p8` | own/tab | M\_V of a 1.8 M☉ MIST star at Teff = 7500 K on the main sequence minus M\_V of the observed fundamental-mode delta Scuti P-L relation at its period from Q = 0.0329 d |
| `delta-scuti-pl-2p0` | own/tab | as delta-scuti-pl-1p8 for a 2.0 M☉ star |
| `cmdra-a-radius-excess-solar` | own/tab | observed radius of CM Dra A over the MIST radius of its mass (4 Gyr, \[Fe/H\] = 0), minus 1 |
| `cmdra-a-radius-excess-lowz` | own/tab | the same for \[Fe/H\] = -0.3 (the metallicity of Terrien et al. 2012) |
| `czev343-apsidal-period` | R | CzeV343: observed apsidal period of binary A, full N-body quadruple \[yr\] |
| `sse-sun-radius` | own/tab | Sun: radius of a 1 M☉, Z = 0.02 SSE star at 4.57 Gyr |
| `sse-sun-luminosity` | own/tab | Sun: luminosity of a 1 M☉, Z = 0.02 SSE star at 4.57 Gyr |
| `sse-sun-convection-zone` | own/tab | Sun: base of the convective envelope r\_bcz/R (BSE envelope model) |
| `iw-q-solar` | own | log10 Q' of inertial waves, 1 M☉ main-sequence star, P\_rot = 10 d |
| `iw-q-m-dwarf` | own | log10 Q' of inertial waves, fully convective 0.3 M☉ star, P\_rot = 10 d |
| `tides-circularization-m67` | own | circularization period of solar-type binaries at 4 Gyr (M67), Hurley et al. (2002) tides (observed: 12.1 +1.0/−1.5 d, Meibom & Mathieu 2005, Table 3) |
| `tides-circularization-m67-iw` | own | circularization period at 4 Gyr (M67) with inertial waves, main sequence only (equilibrium tide only: 3.8 d; observed as above) |
| `gw-hulse-taylor-pdot` | own | orbital decay dP/dt of PSR B1913+16 by gravitational waves (masses 1.438 and 1.390 M☉, P = 0.323 d, e = 0.617) |
| `epsaur-mid-eclipse` | own | epsilon Aur: mid-eclipse of 2009-2011 from the spectroscopic orbit |
| `epsaur-duration` | own | epsilon Aur: duration of the eclipse, contacts I-IV |
| `epsaur-depth` | own | epsilon Aur: depth of the eclipse in V (opaque disc, 0.6 au thick) |
| `hd98800-polar-disc` | own | HD 98800 B: inclination to the sky of a polar disc (along the eccentricity vector of the B orbit) |
| `eecep-disc-precession` | own | precession period \[yr\] of the EE Cep disc from the torque of the B star on a rigid disc (illustrative disc of configs/catalog/ee\_cep.yaml, r\_out 0.8 au) |
| `tides-circularization-pms-population` | own | circularization period (final period, e < 0.01) of twin 1 M☉ binaries at 3 Myr (MIST pre-main sequence from 1 Myr, equilibrium tide with linear damping and inertial waves; observed pre-main-sequence binaries: 7.1 ± 1.2 d, Meibom & Mathieu 2005, Table 3) |
| `tides-circularization-m67-final` | own | circularization period (final period) of twin 1 M☉ binaries at 4 Gyr (MIST from 1 Myr, equilibrium tide with linear damping and inertial waves; observed M67: 12.1 +1.0/−1.5 d) |

**Table A8.** Level D: results (columns *Expected*, *Computed*, *Difference*, *Tolerance*, *Verdict*), in the same order as Table A7.

| Id | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|
| `v994s-dyn-A` | 0.12 | 0.0999365 | −16.72 % | 25 % | PASS |
| `v994s-dyn-B` | 0.05 | 0.0519361 | +3.87 % | 25 % | PASS |
| `v994s-dyn-C` | 0 | 1.23366·10⁻⁶ | +1.23·10⁻⁶ | 0.01 | PASS |
| `v994s-gr-A` | 0.1 | 0.0938886 | −6.11 % | 10 % | PASS |
| `v994s-gr-B` | 0.08 | 0.0748582 | −6.43 % | 10 % | PASS |
| `v994s-gr-C` | 0.06 | 0.0811201 | +35.20 % | 10 % | **EXPLAINED** |
| `qzcar-evolutionary-mass-ac1` | 25.5 M☉ | 27.0258 M☉ | +5.98 % | 15 % | PASS |
| `qzcar-evolutionary-mass-ac2` | 33.2 M☉ | 26.0243 M☉ | −21.61 % | 15 % | **EXPLAINED** |
| `qzcar-evolutionary-mass-aa1` | 69.8 M☉ | 42.5172 M☉ | −39.09 % | 15 % | **EXPLAINED** |
| `agb-final-mass-1p5` | 0.609 M☉ | 0.565426 M☉ | −0.0436 M☉ | 0.1 M☉ | PASS |
| `agb-final-mass-2p5` | 0.689 M☉ | 0.610396 M☉ | −0.0786 M☉ | 0.1 M☉ | PASS |
| `agb-final-mass-3p0` | 0.745 M☉ | 0.653775 M☉ | −0.0912 M☉ | 0.1 M☉ | PASS |
| `agb-final-mass-4p0` | 0.899 M☉ | 0.835342 M☉ | −0.0637 M☉ | 0.1 M☉ | PASS |
| `murphy2014-kic9651065-a1` | 165.57 s | 167.767 s | +1.33 % | 6 % | PASS |
| `murphy2014-kic9651065-a2-a1` | 0.231986 | 0.219651 | −0.0123 | 0.025 | PASS |
| `delta-scuti-pl-1p8` | 0 mag | −0.0251584 mag | −0.0252 mag | 0.25 mag | PASS |
| `delta-scuti-pl-2p0` | 0 mag | −0.115268 mag | −0.115 mag | 0.25 mag | PASS |
| `cmdra-a-radius-excess-solar` | 0.035 | 0.0335989 | −0.0014 | 0.015 | PASS |
| `cmdra-a-radius-excess-lowz` | 0.06 | 0.0641756 | +0.00418 | 0.015 | PASS |
| `czev343-apsidal-period` | 33.4 yr | 34.491 yr | +3.27 % | 6.0 % | PASS |
| `sse-sun-radius` | 1 R☉ | 0.984843 R☉ | −1.52 % | 5 % | PASS |
| `sse-sun-luminosity` | 1 L☉ | 0.95786 L☉ | −4.21 % | 5 % | PASS |
| `sse-sun-convection-zone` | 0.713 | 0.733903 | +0.0209 | 0.03 | PASS |
| `iw-q-solar` | 7 | 7.0122 | +0.0122 | 0.3 | PASS |
| `iw-q-m-dwarf` | 7 | 6.94697 | −0.053 | 0.3 | PASS |
| `tides-circularization-m67` | 12.1 d | 3.7847 d | −68.72 % | 30.0 % | **EXPLAINED** |
| `tides-circularization-m67-iw` | 12.1 d | 7.4572 d | −38.37 % | 30 % | **EXPLAINED** |
| `gw-hulse-taylor-pdot` | −2.40263·10⁻¹² | −2.4021·10⁻¹² | +0.02 % | 0.1 % | PASS |
| `epsaur-mid-eclipse` | 2.45543·10⁶ JD | 2.45541·10⁶ JD | −15 JD | 30 JD | PASS |
| `epsaur-duration` | 740 d | 638 d | −13.78 % | 15 % | PASS |
| `epsaur-depth` | 0.75 mag | 0.83134 mag | +0.0813 mag | 0.15 mag | PASS |
| `hd98800-polar-disc` | 26° | 28.4735° | +2.47° | 3° | PASS |
| `eecep-disc-precession` | 65 yr | 183.163 yr | +181.79 % | 15 % | **EXPLAINED** |
| `tides-circularization-pms-population` | 7.1 d | 6.5831 d | −0.517 d | 1.2 d | PASS |
| `tides-circularization-m67-final` | 12.1 d | 7.5337 d | −4.57 d | 1.2 d | **EXPLAINED** |

\normalsize
