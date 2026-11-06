| Id | Tier | Quantity | Expected | quadsim | Difference | Result |
| --- | --- | --- | --- | --- | --- | --- |
| apsidal-sterne-rotating | A | apsidal period from tidal + rotational distortion, e = 0.3, Omega = n | 291.21 yr | 291.18 yr | -0.01 % | PASS |
| apsidal-sterne-static | A | apsidal period from the tidal distortion, e = 0.3, no rotation | 318.23 yr | 318.2 yr | -0.01 % | PASS |
| brouwer-law | A | growth exponent of the IAS15 energy error (Brouwer's law) | 0.5 | 0.59301 | +0.093 | PASS |
| ce-impulsive-outer-a | A | outer semi-major axis after a common envelope inside a triple (N-body) relative to before, divided by the impulsive-mass-loss prediction for a circular orbit | 1 | 1.0127 | +1.27 % | PASS |
| ce-impulsive-outer-e | A | outer eccentricity after the common envelope of the same triple minus the impulsive-mass-loss prediction | 0 | 0.0044215 | +0.00442 | PASS |
| claret2004-beta | A | Claret (2004): gyration radius beta of 39.8 Msun at log g = 3.416 | 0.16171 | 0.1617 | -6.61e-06 | PASS |
| claret2004-log-k2 | A | Claret (2004): log k2 of 39.8 Msun at log g = 3.416 (3.65 Myr) | -2.8152 | -2.8157 | -0.000463 | PASS |
| claret2023-beta-1msun | A | Claret (2023) grid: radius of gyration beta of 1.0 Msun, [Fe/H] = 0, log g = 4.216 | 0.2265 | 0.2265 | -6.34e-07 | PASS |
| claret2023-k2-2msun | A | Claret (2023) grid: log k2 of 2.0 Msun, [Fe/H] = 0, log g = 4.069 | -2.478 | -2.4779 | +6e-05 | PASS |
| contact-jdot | A | angular-momentum loss time [Gyr] of a 1.0 + 0.5 Msun contact binary, P = 0.35 d (quadsim units vs cgs) | 2.1889 Gyr | 2.189 Gyr | +0.00 % | PASS |
| disc-face-on-tau | A | face-on disc of optical depth 0.7 in front of a star: transmitted fraction | 0.49659 | 0.49659 | +0 | PASS |
| disc-inclined-tau | A | disc of optical depth 0.7 inclined 60 deg in front of a star: transmitted fraction | 0.2466 | 0.2466 | -4.16e-16 | PASS |
| disc-polar-critical | A | critical inclination [deg] for polar libration of a disc around a binary with e = 0.5, from integrations of quadsim's rigid-disc equation | 37.761 deg | 37.761 deg | +6.29e-09 deg | PASS |
| disc-ring-precession | A | precession rate of a narrow rigid ring at 4 a_b around a circular binary / test-particle rate | 1 | 1.0002 | +0.02 % | PASS |
| disc-strip-uniform | A | opaque edge-on slab (half thickness 0.3 R) across a uniform star: transmitted fraction | 0.62384 | 0.62475 | +0.000912 | PASS |
| escape-theory-energy-slope-L0 | A | statistical escape theory, equal masses, L0 = 0: slope of dP/dE_B for/E_B/>>/E0/ | -3.5 | -3.5026 | -0.0026 | PASS |
| escape-theory-energy-slope-highL | A | statistical escape theory, L0 = L_c: slope of dP/dE_B | -4 | -4.0248 | -0.0248 | PASS |
| escape-theory-mean-e-L0 | A | statistical escape theory, L0 = 0: mean eccentricity of the binary | 0.8133 | 0.81282 | -0.000483 | PASS |
| escape-theory-mean-e-highL | A | statistical escape theory, L0 = L_c: mean eccentricity of the binary | 0.7 | 0.71078 | +0.0108 | PASS |
| flyby-first-order-kick | A | relative difference of the eccentricity-vector change of a binary in a parabolic passage at Q = 20 a: first-order formula against direct N-body | 0 | 0.0059492 | +0.00595 | PASS |
| galactic-tide-emax-60 | A | largest eccentricity of a 2e4 AU orbit inclined by 60 deg to the Galactic disc (secular backend, physics.galactic_tide) | 0.82916 | 0.82918 | +2.5e-05 | PASS |
| galactic-tide-emax-60-nbody | A | largest eccentricity of a 2e4 AU binary inclined by 60 deg to the disc, N-body backend with the Galactic tide (from e = 0.01) | 0.82916 | 0.82908 | -0.01 % | PASS |
| gas-alignment-i | A | gas phase: inclination to the disc after 3000 yr, i0 = 40 deg, tau_i = 1500 yr | 5.64 deg | 5.6394 deg | -0.000561 deg | PASS |
| gas-damping-e | A | gas phase: eccentricity after 3000 yr with tau_e = 2000 yr | 0.11157 | 0.11157 | +0.00 % | PASS |
| gas-migration-a | A | gas phase: semi-major axis after 3000 yr with tau_a = 5000 yr [AU] | 5.4881 AU | 5.4881 AU | -0.00 % | PASS |
| gr-mercury | A | 1PN perihelion advance of Mercury | 42.98 arcsec/century | 42.98 arcsec/century | +0.00 % | PASS |
| gw-nbody-decay | A | the same binary with e = 0.3, N-body backend: change of a over 10 yr / Peters rate x 10 yr | 1 | 0.99992 | -0.01 % | PASS |
| gw-secular-circular-decay | A | circular 1.4 + 1.4 Msun binary in a triple, secular backend with physics.gw: a / a0 (1 - t/tau)^(1/4) at 90 % of the merger time | 1 | 1 | -0.00 % | PASS |
| k2-rotation | A | rotational correction of k2 (2 Msun, 2 Rsun, V = 200 km/s) / 10^(-2 V^2 / 3 g R) | 1 | 1 | +0.00 % | PASS |
| ld-tlusty-row | A | TLUSTY limb darkening: a1 of OG40000g400v10 in Bessell V | 1.1826 | 1.1826 | +0 | PASS |
| matt15-unsaturated-powerlaw | A | Omega(4.55 Gyr) / Omega(1.14 Gyr) of the braked solar model | 0.5 | 0.49195 | -0.00805 | PASS |
| megno-chaotic | A | MEGNO > 4 for a 2+2 inside the unstable region (1 = yes) | 1 | 1 | +0 | PASS |
| megno-regular | A | MEGNO of a regular (quasi-periodic) hierarchical triple | 2 | 2.0016 | +0.00162 | PASS |
| mesh-sphere-eclipse | A | largest difference of the light curves of two spheres through both eclipses: surface meshes (quadsim.mesh) against the exact tier-A eclipses | 0 | 0.00025964 | +0.00026 | PASS |
| mist-grid-tams-age-1 | A | MIST grid in quadsim: age of 1 Msun at the end of the main sequence (EEP 454), [Fe/H] = 0, v/vcrit = 0 | 9919.4 Myr | 9919.4 Myr | +0.00 % | PASS |
| mist-interpolation | A | MIST: radius of 1.14 Msun at the TAMS interpolated from 1.12 and 1.16 Msun / tabulated | 1 | 1.0107 | +1.07 % | PASS |
| mist-interpolation-massive-m | A | MIST: current mass of 100 Msun at the TAMS interpolated from 95 and 105 Msun / tabulated | 1 | 1.0011 | +0.11 % | PASS |
| mist-interpolation-massive-r | A | MIST: radius of 60 Msun in the middle of the main sequence (EEP 330) interpolated from 55 and 65 Msun / tabulated | 1 | 0.99943 | -0.06 % | PASS |
| mist-vink2001 | A | MIST wind of 40 Msun at 30 % of its main sequence / Vink et al. (2001) eq. 24 | 1 | 1.0551 | +5.51 % | PASS |
| ml-adiabatic-invariant | A | max abs(a M / (a M)_0 - 1) while m1 loses 0.3 Msun slowly (Psi = 5.6e-6), N-body | 0 | 3.5681e-09 | +3.57e-09 | PASS |
| ml-circular-emax | A | largest e / Psi of an initially circular binary under slow mass loss | 2 | 2.0026 | +0.13 % | PASS |
| ml-impulsive-a | A | a'/a after a circular binary loses 15 % of its mass within 1e-3 of a period | 1.2143 | 1.2149 | +0.05 % | PASS |
| ml-impulsive-e | A | eccentricity after a circular binary loses 15 % of its mass within 1e-3 P | 0.17647 | 0.17689 | +0.24 % | PASS |
| ml-osculating-e-03 | A | max abs(e - e_theory) / Psi over two orbits, e0 = 0.3, Psi = 5.6e-5, N-body | 0 | 0.012647 | +0.0126 | PASS |
| ml-osculating-e-07 | A | max abs(e - e_theory) / Psi over two orbits, e0 = 0.7, Psi = 5.6e-5, N-body | 0 | 0.012354 | +0.0124 | PASS |
| ml-zlk-nbody-secular | A | time of the 4th ZLK rise (e > 0.6) while the perturber loses half its mass, N-body / secular quadrupole | 1 | 1.0021 | +0.21 % | PASS |
| mt-conservative | A | separation after conservative transfer (donor 3 -> 1 Msun, accretor 2 -> 4) | 2.25 a/a0 | 2.25 a/a0 | -0.00 % | PASS |
| mt-isotropic-reemission | A | separation after transfer with complete isotropic re-emission (donor 3 -> 1 Msun, accretor 2 Msun) | 2.03 a/a0 | 2.03 a/a0 | -0.00 % | PASS |
| mt-soberman-general | A | separation after transfer with constant fractions (donor 3 -> 1.5 Msun, accretor 2 Msun; 30 % donor wind, 40 % isotropic re-emission, 30 % accreted): numerical orbit response against the closed form | 1.2957 a/a0 | 1.2957 a/a0 | -0.00 % | PASS |
| nodal-precession | A | nodal precession period, circular inner orbit, i = 20 deg | 11351 yr | 11340 yr | -0.10 % | PASS |
| ocfit-gimenez-bastero | A | apsidal term of the times of minima, e = 0.3: largest difference between the exact conjunction times and the series of Gimenez & Bastero (1995) used by OCFit | 0 s | 0.00062984 s | +0.00063 s | PASS |
| ocfit-lite-amplitude | A | half peak-to-peak light-time effect over one orbit / K of Irwin (1952), e3 = 0.55 | 1 | 1 | -0.00 % | PASS |
| ocfit-nbody-P3 | A | the same fit: fitted / configured outer period | 1 | 0.99996 | -0.00 % | PASS |
| ocfit-nbody-asini-A | A | stage-1 fit of N-body eclipse timings of a 2+2 (outer P = 1500 d, e = 0.3, i = 60 deg): fitted / Keplerian a sin i of pair A | 1 | 0.99864 | -0.14 % | PASS |
| ocfit-nbody-asini-B | A | the same for pair B (anti-phase, shared outer orbit) | 1 | 1.0024 | +0.24 % | PASS |
| phenom-romer-delay | A | phenomenological minima: delay of the secondary minimum after half the period, 2 + 1 Msun ZAMS, P = 5 d, circular, i = 89 deg | 27.452 s | 27.347 s | -0.38 % | PASS |
| phoebe-romer-delay | A | PHOEBE hybrid: delay of the secondary minimum after half the period, 2 + 1 Msun, P = 5 d, circular (with minus without light time) | 27.452 s | 27.287 s | -0.60 % | PASS |
| polytrope-k2-n1.5 | A | apsidal-motion constant of a polytrope n = 1.5 | 0.1433 | 0.14328 | -0.01 % | PASS |
| polytrope-k2-n3 | A | apsidal-motion constant of a polytrope n = 3 | 0.0144 | 0.014443 | +0.30 % | PASS |
| pulsation-time-delay-circular | A | amplitude of the light-travel-time delays of a pulsating star in a circular binary (1.7 + 0.5 Msun, P = 100 d, i = 60 deg), measured from the synthetic light curve like Kepler data (30-min cadence, pulsation phases in 2-d segments) | 53.859 s | 53.859 s | +0.00 % | PASS |
| quad-tp-emax | A | secular quadrupole code: ZLK e_max, test particle, i0 = 65 deg | 0.83805 | 0.83808 | +0.00 % | PASS |
| quint-out-zlk-emax | A | 2+2+1: ZLK e_max of the orbit of the two binaries driven by a fifth star, i0 = 65 deg | 0.81272 | 0.8152 | +0.31 % | PASS |
| rejuvenation-hurley | A | main-sequence age of a 2 Msun accretor after gaining 1 % of its mass / (t'_MS / t_MS)(M / M') t | 1 | 0.99917 | -0.08 % | PASS |
| roche-eggleton-q0.1 | A | Roche-lobe radius / separation, q = 0.1 (Eggleton formula vs exact) | 0.2054 | 0.20677 | +0.67 % | PASS |
| roche-eggleton-q1 | A | Roche-lobe radius / separation, q = 1 (Eggleton formula vs exact) | 0.37986 | 0.37892 | -0.25 % | PASS |
| roche-lobe-volume | A | volume-equivalent radius of the Roche lobe (q = M2/M1 = 0.3) from the surface mesh of quadsim.roche / Eggleton's formula | 1 | 1.0051 | +0.51 % | PASS |
| rotating-star-critical | A | equatorial / polar radius of a Roche-model star at critical rotation | 1.5 | 1.5 | -0.00 % | PASS |
| rotating-star-equator | A | R_eq/R_pol - (1 + w^2/2) with w = Omega / Omega_K(R_eq), Roche model, Omega / Omega_crit = 0.6 | 0 | -2.2204e-16 | -2.22e-16 | PASS |
| sbi-linear-gaussian-std | A | SBI (neural posterior estimation): posterior width of a linear Gaussian model | 0.30923 | 0.29122 | -5.82 % | PASS |
| scatter-hyperbola-deflection | A | scattering set-up: deflection angle of a 1 + 0.5 Msun hyperbolic passage, v_inf = 5 km/s, b = 3 AU | 173.55 deg | 173.55 deg | +1.36e-12 deg | PASS |
| secular-3p1-kozai-emax | A | 3+1 secular backend: ZLK e_max of the inner orbit, test particle, i0 = 65 deg | 0.83805 | 0.83808 | +0.00 % | PASS |
| secular-gr-apsidal | A | secular backend: relativistic apsidal period, a = 0.05 AU, e = 0.3 | 7113.8 yr | 7113.8 yr | -0.00 % | PASS |
| secular-kozai-emax | A | secular backend: ZLK maximum eccentricity, i0 = 65 deg, test particle | 0.83805 | 0.83808 | +0.00 % | PASS |
| secular-nested-tides-hut81 | A | 3+1 in the nested secular backend: tidal de/dt of the inner binary (P = 4 d, e = 0.05, synchronous) / Hut (1981) eq. 10 | 1 | 1 | -0.00 % | PASS |
| secular-nodal-precession | A | secular backend: nodal precession period, i = 20 deg [yr] | 11351 yr | 11351 yr | -0.00 % | PASS |
| sse-zams-radius-1 | A | SSE track grid: ZAMS radius of 1 Msun, Z = 0.02 | 0.88825 Rsun | 0.88831 Rsun | +0.01 % | PASS |
| sse-zams-radius-2 | A | SSE track grid: ZAMS radius of 2.5 Msun, Z = 0.02 | 1.7943 Rsun | 1.7944 Rsun | +0.00 % | PASS |
| stability-stable | A | 2+2 at 2.5 x Mardling-Aarseth limit keeps hierarchy (0 = no change) | 0 | 0 | +0 | PASS |
| stability-unstable | A | 2+2 at 0.5 x Mardling-Aarseth limit changes hierarchy (1 = yes) | 1 | 1 | +0 | PASS |
| tides-leconte-a | A | tidal da/dt, N-body / theory, e = 0.3, non-rotating stars | 1 | 0.99929 | -0.07 % | PASS |
| tides-leconte-e | A | tidal de/dt, N-body / theory, e = 0.3, non-rotating stars | 1 | 0.99936 | -0.06 % | PASS |
| tides-leconte-e-eccentric | A | tidal de/dt, N-body / theory, e = 0.6, Omega = 1.2 n | 1 | 0.99692 | -0.31 % | PASS |
| tides-leconte-spin | A | tidal spin-up dOmega/dt, N-body / theory, e = 0.6, Omega = 1.2 n | 1 | 1.0011 | +0.11 % | PASS |
| timing-ltte-amplitude | A | synthetic O-C: light-time semi-amplitude of binary A, circular edge-on outer orbit, a = 10 AU, equal masses | 41.584 min | 41.586 min | +0.01 % | PASS |
| wind-transfer-jeans | A | 25 + 18 Msun MIST binary with winds and mass transfer: largest relative change of a (M1 + M2) before the transfer starts | 0 | 1.6527e-11 | +1.65e-11 | PASS |
| zlk-below-critical | A | maximum eccentricity below the critical angle, i0 = 30 deg | 0 | 0.0017304 | +0.00173 | PASS |
| zlk-emax-65 | A | ZLK maximum eccentricity, i0 = 65 deg, test particle | 0.83805 | 0.84023 | +0.26 % | PASS |
| zlk-emax-80 | A | ZLK maximum eccentricity, i0 = 80 deg, test particle | 0.97455 | 0.97615 | +0.16 % | PASS |
| zlk-imin-65 | A | mutual inclination at maximum eccentricity, i0 = 65 deg | 39.232 deg | 39.344 deg | +0.113 deg | PASS |
| zlk-kozai-constant | A | secular drift of sqrt(1-e^2) cos i during ZLK cycles, i0 = 65 deg | 0 | 5.1681e-05 | +5.17e-05 | PASS |
| zlk-period-nbody | A | ZLK period, N-body relative to the quadrupole secular period (i0 = 65 deg) | 0 | 0.0066049 | +0.0066 | PASS |
| applegate-huaqr | B | energy needed by the Applegate mechanism over the energy radiated in one modulation period, HU Aqr (thin-shell model) | 0.1 | 0.10566 | +5.66 % | PASS |
| applegate-nnser | B | energy needed by the Applegate mechanism over the energy radiated in one modulation period, NN Ser (thin-shell model) | 3.2 | 3.3034 | +3.23 % | PASS |
| applegate-v471tau | B | energy needed by the Applegate mechanism over the energy radiated in one modulation period, V471 Tau (thin-shell model) | 0.014 | 0.014369 | +2.63 % | PASS |
| at16-binbin-exchange | B | binary-binary scattering, equal masses, a1 = a2, e = 0, v_inf = v_crit: cross-section of 'exchange' / (pi a^2) | 0.302 | 0.37158 | +23.04 % | PASS |
| at16-binbin-exchange-ionization | B | binary-binary scattering, equal masses, a1 = a2, e = 0, v_inf = v_crit: cross-section of 'exchange + ionization' / (pi a^2) | 3.259 | 3.2299 | -0.89 % | PASS |
| at16-binbin-single-ionization | B | binary-binary scattering, equal masses, a1 = a2, e = 0, v_inf = v_crit: cross-section of 'single ionization' / (pi a^2) | 2.645 | 2.6868 | +1.58 % | PASS |
| ce-bse-algol-second-period | B | orbital period after the second common envelope of the BSE Algol example (2.85 Msun giant with a 0.416 Msun core, 0.417 Msun white dwarf, P = 18 d, alpha_CE = 3) | 0.04 d | 0.036564 d | -8.59 % | PASS |
| ce-tout97-period | B | orbital period after the common envelope of the CV progenitor of Tout et al. (1997): 5.36 Msun AGB star (1.014 Msun core) with a 1 Msun companion, alpha_CE = 1 | 35 h | 34.84 h | -0.46 % | PASS |
| chaosmap-circumbinary-limit | B | MEGNO map: inner stability limit of a circumbinary orbit, equal-mass circular binary, all of 8 phases regular [a_bin] | 2.27 | 2.3 | +0.03 | PASS |
| czev343-apsidal-gr | B | CzeV343 binary A: relativistic apsidal rate [rad/day] | 1.47e-05 rad/d | 1.4787e-05 rad/d | +0.59 % | PASS |
| czev343-apsidal-rotational | B | CzeV343 binary A: rotational apsidal rate, pseudo-synchronous spins [rad/day] | 6.4e-05 rad/d | 6.3253e-05 rad/d | -1.17 % | PASS |
| czev343-apsidal-tidal | B | CzeV343 binary A: tidal apsidal rate [rad/day] | 0.000414 rad/d | 0.00042227 rad/d | +2.00 % | PASS |
| escape-theory-vs-nbody-light-escaper | B | resonant binary-single encounters (1 + 0.6 Msun binary, 0.8 Msun star): N-body fraction in which the 0.6 Msun star escapes vs the statistical theory | 0.769 | 0.75369 | -0.0153 | PASS |
| flyby-nbody-backend | B | parabolic passage at Q = 20 a integrated by the N-body backend (physics.encounters): relative difference of the eccentricity kick from the first-order result of Hamers & Samsing (2019) | 0 | 0.007384 | +0.00738 | PASS |
| iw-two-layer | B | inertial-wave dissipation of a two-layer star (alpha 0.6, beta 0.95): impulsive formalism / closed formula | 1 | 1.0003 | +0.03 % | PASS |
| matt15-sun-period | B | rotation period of a solar model at 4.55 Gyr with the wind torque of Matt et al. (2015), started at P = 1 d | 27.97 d | 27.759 d | -0.75 % | PASS |
| mist-geneva-tauh-120-0 | B | hydrogen-burning lifetime of 120 Msun, non-rotating, MIST vs Geneva | 2.672 Myr | 2.7828 Myr | +4.15 % | PASS |
| mist-geneva-tauh-120-4 | B | hydrogen-burning lifetime of 120 Msun, rotating, MIST vs Geneva | 3.178 Myr | 3.3313 Myr | +4.82 % | PASS |
| mist-geneva-tauh-15-0 | B | hydrogen-burning lifetime of 15 Msun, non-rotating, MIST vs Geneva | 11.015 Myr | 12.354 Myr | +12.15 % | PASS |
| mist-geneva-tauh-15-4 | B | hydrogen-burning lifetime of 15 Msun, rotating, MIST vs Geneva | 13.447 Myr | 12.682 Myr | -5.69 % | PASS |
| mist-geneva-tauh-60-0 | B | hydrogen-burning lifetime of 60 Msun, non-rotating, MIST vs Geneva | 3.53 Myr | 3.7179 Myr | +5.32 % | PASS |
| mist-geneva-tauh-60-4 | B | hydrogen-burning lifetime of 60 Msun, rotating, MIST vs Geneva | 4.465 Myr | 3.9466 Myr | -11.61 % | PASS |
| mist-ms-lifetime-1 | B | main-sequence lifetime of 1 Msun (ZAMS to TAMS), SSE vs MIST | 9878 Myr | 11003 Myr | +11.39 % | PASS |
| mist-ms-lifetime-2 | B | main-sequence lifetime of 2 Msun (ZAMS to TAMS), SSE vs MIST | 1054 Myr | 1164.2 Myr | +10.45 % | PASS |
| mist-tams-radius-2 | B | radius of 2 Msun at the end of the main sequence, SSE vs MIST | 3.913 Rsun | 3.5658 Rsun | -8.87 % | PASS |
| ml-tedi-onset | B | onset of the inner instability (e_in > 0.1) / crossing of the MA01 limit, inner primary losing 3 of 5 Msun in 8000 yr | 1 | 1.0647 | +6.47 % | PASS |
| mt-bse-algol-giant-phase | B | period at the end of the BSE Algol example's transfer, computed from its state at the end of the Hertzsprung gap (0.53 + 2.74 Msun, 17 d; conservative) | 30 d | 30.356 d | +1.19 % | PASS |
| mt-zeta-lobe-q2 | B | response zeta_L = dln R_L / dln M_d of the donor's Roche lobe (Eggleton 1983) to conservative transfer at q = 2 | 2.59 | 2.6071 | +0.0171 | PASS |
| naoz13-fig3-flip-time | B | Naoz et al. 2013 Fig. 3 (eccentric ZLK): time of the first flip through 90 deg | 3.8e+06 yr | 3.9868e+06 yr | +4.92 % | PASS |
| naoz13-fig3-imax | B | Naoz et al. 2013 Fig. 3: maximum inclination after the flip | 142.6 deg | 143.01 deg | +0.408 deg | PASS |
| naoz13-fig4-imax | B | Naoz et al. 2013 Fig. 4: maximum mutual inclination (orbit flip to retrograde) | 106.7 deg | 106.93 deg | +0.23 deg | PASS |
| naoz13-fig4-imin | B | Naoz et al. 2013 Fig. 4 (PSR B1620-26 setup): minimum mutual inclination | 57.5 deg | 57.239 deg | -0.261 deg | PASS |
| pejcha13-F6-fA | B | Pejcha et al. F6 population, N-body: binary A reaches 1 - e < 1e-3 | 0.358 | 0.24 | -0.118 | PASS |
| pejcha13-F6-fB | B | Pejcha et al. F6 population, N-body: binary B reaches 1 - e < 1e-3 | 0.53 | 0.41 | -0.12 | PASS |
| quad-naoz13-fig4-imax | B | secular quadrupole code: Naoz et al. 2013 Fig. 4, maximum inclination | 106.7 deg | 106.72 deg | +0.0207 deg | PASS |
| quad-naoz13-fig4-imin | B | secular quadrupole code: Naoz et al. 2013 Fig. 4, minimum inclination | 57.5 deg | 57.538 deg | +0.0376 deg | PASS |
| qzcar-broz-a1 | B | QZ Car: semi-major axis of Ac1+Ac2 (P = 6.0 d) from the Table 2 elements and masses | 54.032 Rsun | 53.995 Rsun | -0.07 % | PASS |
| qzcar-broz-a2 | B | QZ Car: semi-major axis of Aa1+Aa2 (P = 20.7 d) from the Table 2 elements and masses | 136.13 Rsun | 136.08 Rsun | -0.04 % | PASS |
| qzcar-broz-a3 | B | QZ Car: semi-major axis of Ac+Aa (P = 14 728 d) from the Table 2 elements and masses | 13053 Rsun | 13046 Rsun | -0.05 % | PASS |
| secular-naoz13-fig3-flip-time | B | 3+1 secular backend: Naoz et al. 2013 Fig. 3, time of the first flip [yr] | 3.8e+06 yr | 4.0498e+06 yr | +6.57 % | PASS |
| secular-naoz13-fig3-imax | B | 3+1 secular backend: Naoz et al. 2013 Fig. 3, maximum inclination after the flip | 142.6 deg | 142.56 deg | -0.0435 deg | PASS |
| tidal-capture-radius-n15 | B | the same for two n = 1.5 stars | 3.84 | 3.3345 | -13.16 % | EXPLAINED |
| tidal-capture-radius-n15-mcmillan | B | largest periastron distance for tidal capture / R, two equal n = 1.5 stars at 10 km/s | 3.4 | 3.3345 | -1.93 % | PASS |
| tidal-capture-radius-n3 | B | largest periastron distance for tidal capture / R, two n = 3 stars (0.6 Msun, 0.65 Rsun) at 10 km/s | 2.72 | 2.366 | -13.01 % | EXPLAINED |
| tides-circularization-pms | B | period below which twin 1 Msun binaries (e0 = 0.3) are circular when they reach the main sequence (MIST pre-main sequence from 1 Myr to 50 Myr; equilibrium tide and inertial waves) | 7.85 d | 7.9209 d | +0.90 % | PASS |
| vok16-eclipse-prob-A | B | V994 Her: chance that binary A eclipses, isotropic orientation | 0.366 | 0.36927 | +0.89 % | PASS |
| vok16-eclipse-prob-B | B | V994 Her: chance that binary B eclipses, isotropic orientation | 0.397 | 0.39908 | +0.52 % | PASS |
| vok16-fig2-jab-1myr | B | secular backend: Fig. 2 setup, maximum J_AB over 1 Myr | 5.2 deg | 5.2006 deg | +0.000633 deg | PASS |
| vok16-fig2-jab-max | B | secular backend: Vokrouhlický 2016 Fig. 2, maximum J_AB over 100 kyr | 5.2 deg | 5.2007 deg | +0.000687 deg | PASS |
| vok16-fig2-jab-min | B | secular backend: Vokrouhlický 2016 Fig. 2, minimum J_AB over 100 kyr | 3.5 deg | 3.4715 deg | -0.0285 deg | PASS |
| vok16-fig3-jab-max | B | secular backend: Vokrouhlický 2016 Fig. 3, maximum J_AB over 20 kyr | 6 deg | 6.2087 deg | +0.209 deg | PASS |
| vok16-jab-max-analytic | B | Vokrouhlický 2016, analytic limit J_AB^max for V994 Her (Fig. 2 setup) | 5.2 deg | 5.2007 deg | +0.000736 deg | PASS |
| vok16-pejcha-fA | B | Pejcha-like 2+2 population, secular: binary A reaches 1 - e < 1e-3 | 0.31 | 0.21 | -0.1 | EXPLAINED |
| vok16-pejcha-fB | B | Pejcha-like 2+2 population, secular: binary B reaches 1 - e < 1e-3 | 0.481 | 0.4025 | -0.0785 | EXPLAINED |
| vok16-v994-JA | B | Vokrouhlicky 2016, V994 Her setup: initial tilt J_A of binary A to the outer orbit | 11.6 deg | 11.646 deg | +0.0462 deg | PASS |
| vok16-v994-JB | B | Vokrouhlicky 2016, V994 Her setup: initial tilt J_B of binary B to the outer orbit | 15.8 deg | 15.793 deg | -0.00678 deg | PASS |
| vyn23-2p2-ma01-score | B | overall score of the MA01 criterion against the N-body ghost test on quadsim's 2+2 training set (sampled as in Vynatheya et al. 2023, Sect. 3.1) | 0.83 | 0.869 | +0.039 | PASS |
| vyn23-2p2-mlp-score | B | held-out score of quadsim's MLP stability emulator (2+2; 11 parameters, 4 x 50 neurons) against the N-body ghost test | 0.94 | 0.9375 | -0.0025 | PASS |
| vyn23-3p1-ma01-score | B | overall score of the MA01 criterion against the N-body ghost test on quadsim's 3+1 training set (sampled as in Vynatheya et al. 2023, Sect. 3.1) | 0.56 | 0.7896 | +0.23 | EXPLAINED |
| vyn23-3p1-mlp-score | B | held-out score of quadsim's MLP stability emulator (3+1; 11 parameters, 4 x 50 neurons) against the N-body ghost test | 0.93 | 0.9195 | -0.0105 | PASS |
| cepheid-period-radius-6 | C | log R of a 6 Msun MIST Cepheid (blue loop, middle of the strip, period from Bono et al. 2000) minus the period-radius relation of rotating Geneva Cepheid models | 0 dex | -0.0076576 dex | -0.00766 dex | PASS |
| cepheid-period-radius-7 | C | as cepheid-period-radius-6 for a 7 Msun Cepheid | 0 dex | 0.0039244 dex | +0.00392 dex | PASS |
| contact-vs-phoebe | C | V amplitude of a contact binary (q = 0.42, fill-out 0.3, i = 80 deg): quadsim mesh minus PHOEBE 2.5 with the same physics | 0 mag | 0.0037365 mag | +0.00374 mag | PASS |
| lc-detached-primary | C | primary eclipse depth of a detached binary of spheres | 0.38082 | 0.38087 | +4.72e-05 | PASS |
| lc-detached-secondary | C | secondary eclipse depth of a detached binary of spheres | 0.14466 | 0.14465 | -3.85e-06 | PASS |
| lc-ellipsoidal | C | semi-amplitude of the ellipsoidal variation, R/a = 0.15, q = 0.8, i = 60 deg | 0.0027949 | 0.0026062 | -6.75 % | PASS |
| lc-occultation-large | C | flux of a limb-darkened star behind a larger dark disk, rp = 1.2, z = 0.5 | 0.12898 | 0.12898 | -5.52e-09 | PASS |
| lc-occultation-limb | C | flux of a limb-darkened star behind a dark disk, rp = 0.3 at the limb (z = 0.9) | 0.94524 | 0.94524 | +3.66e-09 | PASS |
| lc-occultation-small | C | flux of a limb-darkened star behind a dark disk, rp = 0.1, central | 0.98789 | 0.98789 | -1.92e-12 | PASS |
| ld-tlusty-vs-atlas | C | disk-averaged intensity <I>/I(1) in V at 25 000 K, log g 4: TLUSTY non-LTE / ATLAS LTE | 1 | 0.99759 | -0.24 % | PASS |
| mse-2p2-A-emax | C | 2+2 (a = 1, 1.5, 40 AU): maximum e of binary A in 8e4 yr | 0.9494 | 0.94949 | +9.2e-05 | PASS |
| mse-2p2-B-emax | C | 2+2: maximum e of binary B in 8e4 yr | 0.7758 | 0.77623 | +0.000427 | PASS |
| mse-3p1-inner-emax | C | 3+1 (a = 1, 15, 300 AU, i_mid = 85 deg): maximum e of the inner binary in 1e5 yr | 0.9707 | 0.97071 | +1.37e-05 | PASS |
| mse-3p1-inner-t-emax | C | 3+1: time of the maximum e of the inner binary | 49325 yr | 49323 yr | -0.00 % | PASS |
| mse-3p1-mid-emax | C | 3+1: maximum e of the middle orbit | 0.2245 | 0.22446 | -3.86e-05 | PASS |
| mse-3p1-rlof-time | C | 3+1 (a = 0.3, 5, 200 AU, i = 89 deg): onset of Roche-lobe overflow of the inner binary driven by ZLK cycles | 380.03 yr | 380.04 yr | +0.00 % | PASS |
| mt-bse-algol-accretor | C | accretor mass at the end of the first mass transfer of the BSE Algol example | 2.85 Msun | 2.5111 Msun | -0.339 Msun | EXPLAINED |
| mt-bse-algol-donor | C | donor mass at the end of the first mass transfer of the BSE Algol example (2.9 + 0.9 Msun, P = 2.7 d, Z = 0.02) | 0.42 Msun | 0.40434 Msun | -0.0157 Msun | PASS |
| mt-bse-algol-period | C | orbital period at the end of the first mass transfer of the BSE Algol example | 30 d | 21.93 d | -26.90 % | EXPLAINED |
| mt-cosmic-algol-donor | C | donor mass when the first mass transfer of the BSE Algol example ends, against the BSE code in COSMIC 4.2.1 (2.9 + 0.9 Msun, circular 2.7 d, Z = 0.02, acc_lim = -1) | 0.401 Msun | 0.40434 Msun | +0.00334 Msun | PASS |
| phoebe-sphere-vs-tier-a | C | largest difference of the normalised light curves, tier A vs PHOEBE 2.5 with the same physics (spheres, black bodies, Claret limb darkening), e = 0.1 | 0 | 0.00037195 | +0.000372 | PASS |
| qzcar-broz-lv-aa1 | C | QZ Car: share of Aa1 in the V-band light from its radius and Teff | 0.635 | 0.63502 | +2.36e-05 | PASS |
| qzcar-broz-lv-aa2 | C | QZ Car: share of Aa2 in the V-band light from its radius and Teff | 0.01 | 0.010346 | +0.000346 | PASS |
| qzcar-broz-lv-ac1 | C | QZ Car: share of Ac1 in the V-band light from its radius and Teff | 0.22 | 0.22078 | +0.000775 | PASS |
| qzcar-broz-lv-ac2 | C | QZ Car: share of Ac2 in the V-band light from its radius and Teff | 0.133 | 0.13386 | +0.000856 | PASS |
| sextuple-nbody-emax | C | 2+2+2 sextuple: maximum eccentricity of A, N-body (reproducibility) | 0.75207 | 0.75207 | +0.00 % | PASS |
| sextuple-secular-emax | C | 2+2+2 sextuple: maximum eccentricity of A, secular backend vs N-body | 0.75207 | 0.75182 | -0.03 % | PASS |
| sextuple-secular-zlk | C | 2+2+2 sextuple: ZLK period of A, secular backend vs N-body | 70.2 yr | 71.04 yr | +1.20 % | PASS |
| sse-wind-rgb-1msun | C | mass lost by a 1 Msun star (Z = 0.02) on the red giant branch with the Reimers wind, eta = 0.1: SSE track with the SSE winds against the MIST v1.2 (MESA) model at the tip of the giant branch | 0.04158 Msun | 0.042072 Msun | +1.18 % | PASS |
| agb-final-mass-1p5 | D | white-dwarf mass (core at the end of the AGB) of a 1.5 Msun MIST star (solar [Fe/H]) against the white-dwarf initial-final mass relation | 0.609 Msun | 0.56543 Msun | -0.0436 Msun | PASS |
| agb-final-mass-2p5 | D | white-dwarf mass (core at the end of the AGB) of a 2.5 Msun MIST star (solar [Fe/H]) against the white-dwarf initial-final mass relation | 0.689 Msun | 0.6104 Msun | -0.0786 Msun | PASS |
| agb-final-mass-3p0 | D | white-dwarf mass (core at the end of the AGB) of a 3.0 Msun MIST star (solar [Fe/H]) against the white-dwarf initial-final mass relation | 0.745 Msun | 0.65377 Msun | -0.0912 Msun | PASS |
| agb-final-mass-4p0 | D | white-dwarf mass (core at the end of the AGB) of a 4.0 Msun MIST star (solar [Fe/H]) against the white-dwarf initial-final mass relation | 0.899 Msun | 0.83534 Msun | -0.0637 Msun | PASS |
| cmdra-a-radius-excess-lowz | D | the same for [Fe/H] = -0.3 (the metallicity of Terrien et al. 2012) | 0.06 | 0.064176 | +0.00418 | PASS |
| cmdra-a-radius-excess-solar | D | observed radius of CM Dra A over the MIST radius of its mass (4 Gyr, [Fe/H] = 0), minus 1 | 0.035 | 0.033599 | -0.0014 | PASS |
| czev343-apsidal-period | D | CzeV343: observed apsidal period of binary A, full N-body quadruple [yr] | 33.4 yr | 34.491 yr | +3.27 % | PASS |
| delta-scuti-pl-1p8 | D | M_V of a 1.8 Msun MIST star at Teff = 7500 K on the main sequence minus M_V of the observed fundamental-mode delta Scuti P-L relation at its period from Q = 0.0329 d | 0 mag | -0.025158 mag | -0.0252 mag | PASS |
| delta-scuti-pl-2p0 | D | as delta-scuti-pl-1p8 for a 2.0 Msun star | 0 mag | -0.11527 mag | -0.115 mag | PASS |
| eecep-disc-precession | D | precession period [yr] of the EE Cep disc from the torque of the B star on a rigid disc (illustrative disc of configs/catalog/ee_cep.yaml, r_out 0.8 AU) | 65 yr | 183.16 yr | +181.79 % | EXPLAINED |
| epsaur-depth | D | epsilon Aur: depth of the eclipse in V (opaque disc, 0.6 AU thick) | 0.75 mag | 0.83134 mag | +0.0813 mag | PASS |
| epsaur-duration | D | epsilon Aur: duration of the eclipse, contacts I-IV | 740 d | 638 d | -13.78 % | PASS |
| epsaur-mid-eclipse | D | epsilon Aur: mid-eclipse of 2009-2011 from the spectroscopic orbit | 2.4554e+06 JD | 2.4554e+06 JD | -15 JD | PASS |
| gw-hulse-taylor-pdot | D | orbital decay dP/dt of PSR B1913+16 by gravitational waves (masses 1.438 and 1.390 Msun, P = 0.323 d, e = 0.617) | -2.4026e-12 | -2.4021e-12 | +0.02 % | PASS |
| hd98800-polar-disc | D | HD 98800 B: inclination to the sky of a polar disc (along the eccentricity vector of the B orbit) | 26 deg | 28.474 deg | +2.47 deg | PASS |
| iw-q-m-dwarf | D | log10 Q' of inertial waves, fully convective 0.3 Msun star, P_rot = 10 d | 7 | 6.947 | -0.053 | PASS |
| iw-q-solar | D | log10 Q' of inertial waves, 1 Msun main-sequence star, P_rot = 10 d | 7 | 7.0122 | +0.0122 | PASS |
| murphy2014-kic9651065-a1 | D | first orbital harmonic A1 of the time delays of KIC 9651065 (delta Sct in a 272.7-d eccentric orbit) from its orbit (a1 sin i = 0.37 au, e = 0.468, varpi = 2.01) | 165.57 s | 167.77 s | +1.33 % | PASS |
| murphy2014-kic9651065-a2-a1 | D | ratio of the second to the first orbital harmonic of the time delays of KIC 9651065 (sets the eccentricity; e = 0.468, varpi = 2.01 rad) | 0.23199 | 0.21965 | -0.0123 | PASS |
| qzcar-evolutionary-mass-aa1 | D | the same for Aa1 | 69.8 Msun | 42.517 Msun | -39.09 % | EXPLAINED |
| qzcar-evolutionary-mass-ac1 | D | QZ Car: present mass of the single-star MIST model with the radius and Teff of Ac1, against its dynamical mass | 25.5 Msun | 27.026 Msun | +5.98 % | PASS |
| qzcar-evolutionary-mass-ac2 | D | the same for Ac2 | 33.2 Msun | 26.024 Msun | -21.61 % | EXPLAINED |
| sse-sun-convection-zone | D | Sun: base of the convective envelope r_bcz/R (BSE envelope model) | 0.713 | 0.7339 | +0.0209 | PASS |
| sse-sun-luminosity | D | Sun: luminosity of a 1 Msun, Z = 0.02 SSE star at 4.57 Gyr | 1 Lsun | 0.95786 Lsun | -4.21 % | PASS |
| sse-sun-radius | D | Sun: radius of a 1 Msun, Z = 0.02 SSE star at 4.57 Gyr | 1 Rsun | 0.98484 Rsun | -1.52 % | PASS |
| tides-circularization-m67 | D | circularization period of solar-type binaries at 4 Gyr (M67), Hurley et al. (2002) tides | 12.1 d | 2.6431 d | -78.16 % | EXPLAINED |
| tides-circularization-m67-final | D | circularization period (final period) of twin 1 Msun binaries at 4 Gyr (MIST from 1 Myr, equilibrium tide with linear damping and inertial waves) | 12.1 d | 7.5256 d | -4.57 d | EXPLAINED |
| tides-circularization-m67-iw | D | circularization period at 4 Gyr (M67) with inertial waves, main sequence only | 12.1 d | 7.4471 d | -38.45 % | EXPLAINED |
| tides-circularization-pms-population | D | circularization period (final period, e < 0.01) of twin 1 Msun binaries at 3 Myr (MIST pre-main sequence from 1 Myr, equilibrium tide with linear damping and inertial waves) | 7.1 d | 6.5754 d | -0.525 d | PASS |
| v994s-dyn-A | D | V994 Her sextuple: share of the apsidal motion of A caused by B (and C) | 0.12 | 0.099937 | -16.72 % | PASS |
| v994s-dyn-B | D | V994 Her sextuple: share of the apsidal motion of B caused by A (and C) | 0.05 | 0.051936 | +3.87 % | PASS |
| v994s-dyn-C | D | V994 Her sextuple: share of the apsidal motion of C caused by A and B | 0 | 1.2337e-06 | +1.23e-06 | PASS |
| v994s-gr-A | D | V994 Her: relativistic share of the apsidal motion of A | 0.1 | 0.093889 | -6.11 % | PASS |
| v994s-gr-B | D | V994 Her: relativistic share of the apsidal motion of B | 0.08 | 0.074858 | -6.43 % | PASS |
| v994s-gr-C | D | V994 Her: relativistic share of the apsidal motion of C | 0.06 | 0.08112 | +35.20 % | EXPLAINED |
