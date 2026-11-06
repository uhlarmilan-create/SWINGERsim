# Research and student projects {#projects}

SWINGERsim is meant for modelling specific systems from their data. The projects below use what it
already does. Each has a question, the tools it would use and the data it needs.

## Projects with SQUADRA data {#projects-squadra}

**Table 3.** Projects on systems observed by the SQUADRA group. Columns: the project and its
question, what SWINGERsim provides, and what is still needed (usually data).

| project | question | SWINGERsim provides | still needed |
|:--|:--|:--|:--|
| A. Systems with changing inclination | which outer orbits explain the measured inclination changes (CzeV4315, ASASSN-V J075203 and similar), and when do their eclipses appear or disappear? | ensembles with observation weights and evidences for alternative hierarchies; the dynamical fit of minima, radial velocities and light curves | times of minima of both pairs; light-curve inclinations at several epochs |
| B. Eclipse forecasts for the SQUADRA list | for each target, the probability that the eclipses of either pair disappear within 10–50 yr: which systems to observe first | configurations from published solutions, ensembles, the worker for batches | a complete table of periods, masses and outer orbits |
| C. CzeV343 and the 3:2 resonance | is the system in the resonance or passing through it, and how will its O−C change? | apsidal motion, chaos maps (MEGNO), long integrations, the fit chain | resonance-angle diagnostics (Table 2) |
| D. Compact 2+2 systems with outer periods below a year | secular theory fails; what does a full N-body model of their O−C give? | the dynamical fit (stage 2) with all pairs | minima of the group's targets |

## Questions of populations and theory {#projects-theory}

**Table 4.** Wider questions. Columns: project, question, what SWINGERsim provides.

| project | question | SWINGERsim provides |
|:--|:--|:--|
| E. Alignment of 2+2 orbits | are 2+2 orbits mostly aligned (formation in a disc) or isotropic (dynamical capture), judged from the statistics of doubly eclipsing systems? | eclipse probabilities of populations (validated against Vokrouhlický 2016); missing: a survey selection model |
| F. The excess of the 3:2 period ratio | why do 2+2 systems show an excess near 3:2? Among 781 candidates, 16 lie within 1 % of 3:2 against 11 expected for uniform periods (27 against 21 within 2 %; Kolář et al. 2025), and 8 against 12 lie at 2–3 % from it (Kolář 2025, PhD thesis): an excess at the resonance and a deficit next to it, as capture would produce, but in small numbers. Capture into resonance during orbital migration, as known for planets? The present CzeV343 is not in resonance (Part IV). | chaos maps, a parametric gas phase, tides |
| G. Mergers and contact binaries from ZLK | how many 2+2 systems end as mergers, contact binaries or mass-transferring binaries, and when? | ZLK with tides, stellar evolution, mass transfer, common envelope, contact binaries; Part IV is one example |
| H. Instability after mass transfer | how often does mass transfer destabilise a quadruple, and what follows: merger by ZLK, as in CzeV343, or ejection? | the two-stage method of Part IV |
| I. Formation by binary–binary scattering | cross-sections for forming stable 2+2 and 3+1 systems in young clusters | scattering experiments, escape theory, tidal capture |
| J. Gas phase and alignment | which distribution of mutual inclinations does a gas phase leave behind? | parametric gas phase |
| K. Pulsating stars as clocks | do the light-time delays of δ Scuti, β Cephei or γ Dor pulsations together with the eclipse times give the masses of the outer orbit without spectroscopy? | pulsations evaluated at the emission time; minima from synthetic light curves |
| L. Mass transfer in triples | does a tertiary change the outcome of mass transfer (Algols, stripped stars), and does the transfer change the outer orbit and the eclipses? | mass transfer in inner binaries and from a tertiary, with its effect on the other orbits |
| M. Common envelope in multiples | how often does a common envelope unbind or destabilise the outer orbit? | energy prescription with the response of the outer orbits |
| N. Young massive multiples | which star of a young O-star quadruple interacts first (QZ Car-like systems)? | MIST to high masses, winds, effective ages (Part I, QZ Car) |
| O. Hidden multiples | how many eclipsing binaries in TESS, ZTF and ASAS-SN start or stop eclipsing within the survey baseline? | eclipse appearance statistics of a population |

## Student projects {#projects-students}

**Table 5.** Student projects. Columns: level, topic (with the project letter of Tables 3–4), task,
difficulty.

| level | topic | task | difficulty |
|:--|:--|:--|:--|
| Bachelor | eclipse forecasts for ten SQUADRA systems (B) | configurations from published solutions, ensembles, a summary of when eclipses vanish | low |
| Bachelor | binary–binary scattering (I) | scan velocities and orbit ratios, compare with Antognini & Thompson (2016) | low to medium |
| Bachelor | time delays of a δ Scuti star in a TESS eclipsing binary (K) | measure the delays in TESS, compare with a SWINGERsim model | low to medium |
| Master | dynamical model of CzeV4315 (A) | fit of the minima of both pairs with the 2+(2+1) model; a paper with the observers | medium |
| Master | QX Cas: radii, third light and the old eclipses (open question of Section 3) | joint fit of depths, TESS heartbeat and times of minima; DASCH plates | medium |
| Master | chaos maps and the 3:2 resonance in CzeV343 (C, F) | MEGNO maps, comparison with the resonant secular theory of Breiter & Vokrouhlický (2018) | medium to high |
| Master | Algols with tertiaries (L) | ensembles of mass-transferring binaries with tertiaries, comparison with known systems | medium |
| PhD | alignment of 2+2 orbits from the eclipsing population (E) | survey selection model, Bayesian inversion | high |
| PhD | origin of the 3:2 excess (F) | N-body with migration and tides, comparison with the observed systems | high |
| PhD | fates of 2+2 systems (G, H) | ZLK with stellar evolution and mass transfer for a population, comparison with population-synthesis codes | high |

Ensembles of $10^4$ members run in hours on a 12-core computer, so a student can concentrate on
the physics.
