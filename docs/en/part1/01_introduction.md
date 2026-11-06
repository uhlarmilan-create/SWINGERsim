---
title: "Testing SWINGERsim"
subtitle: "Part I: Validation against the literature"
---

```{=latex}
\makeatletter
\addtocontents{toc}{\protect\renewcommand*{\protect\numberline}{\protect\hb@xt@ 3.4em}}
\makeatother
```

**Abstract.** SWINGERsim simulates hierarchical multiple stars – triples, quadruples and larger
nested systems – together with the evolution of their stars, tides, mass transfer, eclipses, light
curves, eclipse timings, fits to observations and ensembles of possible systems. This first of four
test reports asks whether every physical module of SWINGERsim reproduces independent published
results. The evidence has four levels: analytic results (A), published numerical results (B), other
codes run on the same input (C) and observed systems (D). SWINGERsim carries 219 registered
benchmarks, each with its source, its expected value and a tolerance fixed before the comparison.
With the final code all 219 were run: 205 lie within their tolerance and 14 lie outside it for a
documented reason that is not an error of SWINGERsim (a slip or an inconsistency in the paper, a
different definition, the length of the published runs, or a limitation shared by all published
models). No benchmark fails. Three quarters of the benchmarks (165 of 219) test code written for
SWINGERsim itself – the secular backends, the coupling of stellar evolution, tides and mass
transfer to the orbits, the light curves on spheres and surface meshes, the eclipse timings, the
fits and the ensembles – and not only the external N-body integrator it uses. Beyond the benchmarks,
SWINGERsim reproduces the published model of the massive quadruple QZ Car, the dynamics of the SQUADRA
candidates CzeV4315 and ASASSN-V J075203, 23 published multiple stars across the
Hertzsprung–Russell diagram and the stability scores of a published neural classifier. The report
ends with the problems met during testing, written as properties of the software that a user has
to know, and with a list of what is validated, to which accuracy, and what is not.

# Introduction {#introduction}

## Aim of this report {#aim}

A simulation code is useful to astronomers only if its physics is right. That a code runs without
errors, or that it agrees with itself, proves little: a wrong force can conserve energy perfectly.
The aim of this report is therefore narrow and strict. For every physical effect in SWINGERsim we
ask: **does SWINGERsim reproduce an independent, published result for it, and to what accuracy?**
Where it does not, we ask why, and we trace the difference to its source.

Part I covers the physics modules one by one ([Chapters 2–4](#dynamics)), then whole published
systems: the massive quadruple QZ Car ([Chapter 5](#qz-car)), the SQUADRA candidates
([Chapter 6](#squadra-systems)), a catalogue of 23 published multiple stars ([Chapter 7](#catalogue)),
examples that follow a multiple star through a phase of its life ([Chapter 8](#life-cycles)) and the
stability tools ([Chapter 9](#stability)). [Chapter 10](#problems-found) lists the problems that the
tests revealed and how SWINGERsim handles them, [Chapter 11](#conclusions) summarises what is
validated and what is not, and [Appendix A](#appendix-benchmarks) lists all 219 benchmarks with
their results.

## What SWINGERsim is, and which parts are its own code {#own-code}

SWINGERsim follows hierarchical multiple stars: systems of three to about ten stars in which close
pairs orbit each other and the pairs orbit common centres of mass on wider orbits. It integrates the
orbits with one of three **backends** (the engines that advance the system in time):

- the **N-body backend** integrates the equations of motion of all stars directly. It uses the
  integrator IAS15 of the REBOUND library (Rein & Liu 2012; Rein & Spiegel 2015) and three effects of
  the REBOUNDx library (Tamayo et al. 2020): `gr_full` (first post-Newtonian relativity for all
  bodies), `tides_spin` (the equilibrium tide with spins, Eggleton et al. 1998; Lu et al. 2023) and
  `modify_mass` (masses that change with time);
- the **secular backends** average the motion over the orbital periods and follow only the slow
  changes of the orbits (precession, ZLK cycles, defined in [Chapter 2](#dynamics)). They are written
  for SWINGERsim: a quadrupole–octupole code for triples, a 2+2 code after Vokrouhlický (2016) and a
  nested code for any hierarchy after Hamers et al. (2015) and Hamers & Portegies Zwart (2016);
- the **evolving backend** couples a secular backend to evolving stars: it cuts the run into
  segments within which no stellar radius changes by more than 1 % and applies winds, mass transfer,
  common envelopes, tides and magnetic braking between them. It is written for SWINGERsim.

The question "which code does a benchmark test?" matters for the weight of the evidence. A
benchmark of the N-body backend that reproduces, say, the relativistic perihelion advance of
Mercury mainly confirms that REBOUNDx is used correctly – the force law itself is REBOUNDx's. A
benchmark of the secular backend, of the mass-transfer coupling or of the light-curve model tests
equations and numerics that exist only in SWINGERsim; it is the stronger evidence, because nothing
outside SWINGERsim could have made it pass. We read the source code (`src/swingersim/`) and
classify every benchmark into one of four classes (Table 1). The class is given for every benchmark
in the tables of Chapters 2–9 and in [Appendix A](#appendix-benchmarks).

**Table 1.** Classes of benchmarks by the code they test. Column 1 is the label used in all later
tables; column 2 says which code computes the tested number; column 3 lists the main modules; the
last column gives how many of the 219 benchmarks belong to the class.

| Class | What computes the tested number | Modules (examples) | Benchmarks |
|---|---|---|---|
| own | equations and numerics written for SWINGERsim | secular backends, evolving backend, mass transfer and common envelope, tides and apsidal formulae, inertial waves, magnetic braking, Roche geometry and surface meshes, light curves, dust discs, eclipse timings, pulsation timing, fits of minima, ensembles, stability criteria, escape theory, gravitational waves, Galactic tide and passing stars | 126 |
| own/tab | SWINGERsim's interpolation and use of published tables | MIST and SSE tracks, Claret (2004, 2023) structure constants, Claret & Bloemen (2011) and TLUSTY limb darkening | 39 |
| R | REBOUND IAS15 with REBOUNDx effects, driven by SWINGERsim's set-up | hierarchical initialisation, element conversion, collision and hierarchy detection, mass-loss schedule, time-lag convention, operator splitting of SWINGERsim's own forces | 50 |
| ext | an external library through a SWINGERsim wrapper | PHOEBE 2.5 hybrid light curves, the `sbi` package, the scikit-learn neural network trained on SWINGERsim's own training sets | 4 |

Even an "R" benchmark is not a test of REBOUND alone: SWINGERsim builds the initial positions and
velocities from the hierarchical orbital elements of the configuration, converts the integrated
state back to elements, detects collisions and changes of the hierarchy, and sets the parameters of
the REBOUNDx effects (for example the factor between the time lag of SWINGERsim and the `tau` of
`tides_spin`). A benchmark of class R therefore tests this set-up and the external integrator
together. Where a benchmark compares a SWINGERsim backend with an N-body run of the same system (for
example the secular backend against N-body), the class names the tested side, which is SWINGERsim's
own code.

Of the 219 benchmarks, 165 (75 %) test SWINGERsim's own code (classes own and own/tab), 50 the
N-body set-up around REBOUND, and 4 an external library through a wrapper.

## The validation rule {#rule}

Every physical effect in SWINGERsim satisfies five rules.

1. **An independent reference.** The effect is checked against a result that SWINGERsim did not
   produce: an equation from the literature, a figure or table of a paper, another code, or an
   observed system. Internal checks (round trips of coordinate conversions, energy conservation)
   are part of the automatic tests too, but they are not counted as validation.
2. **A recorded source.** Each benchmark stores the paper and the equation, figure or table it
   reproduces, the input values taken from it, the expected value and the tolerance. The command
   `swingersim validate` prints the source next to every result.
3. **A physical tolerance.** The tolerance follows from the approximation behind the reference (for
   example the neglected order of a series, the accuracy of a fitting formula, the error of a
   measurement or the statistical error of an ensemble). It is fixed before the comparison and
   never widened to make a benchmark pass.
4. **Every difference is explained or it is an error.** A result outside its tolerance is either an
   error of SWINGERsim or it has a documented reason. The reason is stored with the benchmark and
   printed with the result.
5. **Results are reproducible.** Every benchmark runs with one command, and its number in this
   report is the number printed by that command with the final code.

## Levels of evidence and verdicts {#levels}

The references fall into four **levels of evidence** (Table 2). The levels answer different
questions. Level A tests the physics in limits in which the exact answer is known. Levels B and C
test regimes in which no formula exists – chaotic dynamics, coupled evolution of stars and orbits.
Level D tests whether the model as a whole describes real stars, including everything that the
model leaves out.

**Table 2.** Levels of evidence. Column 2 says what the reference is, column 3 gives typical
examples, column 4 the number of benchmarks of that level.

| Level | Reference | Examples | Benchmarks |
|---|---|---|---|
| A | analytic result from the literature | ZLK maximum eccentricity, relativistic perihelion advance, Sterne's apsidal formula, Roche-lobe radius, closed forms of mass transfer | 93 |
| B | published numerical result (a figure or table of a paper) | orbital flips of Naoz et al. (2013), mutual inclinations of Vokrouhlický (2016), scattering cross-sections of Antognini & Thompson (2016) | 60 |
| C | another code on the same input | MSE, BSE in COSMIC, PHOEBE 2.5, batman, ellc, the Xitau model of QZ Car, Geneva Cepheid models | 31 |
| D | observed system | V994 Her, CzeV343, QZ Car, the Sun, M67, white dwarfs, PSR B1913+16, ε Aur, HD 98800 | 35 |

Every benchmark ends with one of four **verdicts**:

- **PASS** – the result lies within the tolerance;
- **EXPLAINED** – the result lies outside the tolerance for a documented reason that is not an error
  of SWINGERsim: the input or the definition used by the paper, published numbers that contradict
  each other, the length or sampling of the published runs, or a stated limitation of
  SWINGERsim's physical model that is an open point rather than a coding error. The tolerance is not changed to hide the difference;
- **FAIL** – outside the tolerance without such a reason; this would be an error of SWINGERsim;
- **SKIP** – an optional package that the benchmark needs is not installed.

![**Figure 1.** The 219 benchmarks of SWINGERsim by level of evidence (rows A–D, see Table 2) with the final code. The length of a bar is the number of benchmarks of that level; the blue part counts PASS, the orange part EXPLAINED. The numbers at the end of each bar repeat the counts. There is no FAIL and no SKIP.](../../figures/part1/benchmark_tiers.png)

Figure 1 summarises the result. All 93 analytic benchmarks (level A) pass. The 14 explained
differences lie in levels B (5), C (2) and D (7) – that is, where the reference is a published
simulation, another code or a real star. This is what one expects if SWINGERsim's physics is right:
an analytic limit has no input to disagree about, whereas a paper can quote a number from a
different definition, and a real star contains physics that no model has. Each of the 14 is
discussed where its physics is treated, and they are collected in [Chapter 11](#conclusions).

## Where the numbers come from {#numbers}

The benchmarks are part of the software. `swingersim validate` runs the fast ones (seconds each),
`swingersim validate --all` adds the slow ones (minutes each, for example ensembles of hundreds of
scattering experiments) and `swingersim validate <id>` runs one benchmark by name. Nine *manual*
benchmarks run only when named: they need hours of computing (the N-body ensemble of Pejcha et al.
2013), training sets of 10⁴ N-body runs (the stability classifier of Vynatheya et al. 2023), or an
optional package (PHOEBE, `sbi`).

```bash
swingersim validate                  # the 184 fast benchmarks
swingersim validate --all            # all 210 automatic benchmarks
swingersim validate gr-mercury       # one benchmark by name
swingersim validate pejcha13-F6-fA   # a manual benchmark (hours), runs only when named
```

Every number in this report was computed with the final code of SWINGERsim: the fast benchmarks
in the last run of the suite, the slow benchmarks in a full run on the production machine (a
12-core workstation), and the manual benchmarks in separate runs of the same code. The command also
writes an HTML report with the table of all results and the reasons of the EXPLAINED ones.

The studies of whole systems ([Chapters 5–9](#qz-car)) are simulations with the configuration files
of the SWINGERsim distribution; their numbers come from the outputs of these runs. Where a study is
an ensemble, its size is given with the result.

## How to read this report {#how-to-read}

Each section of Chapters 2–4 follows the same order: the physical effect and its reference, the
input, the result (a table, often a figure), and the agreement stated as a number. The tables of
benchmarks have the same columns throughout:

- **Benchmark** – the identifier by which `swingersim validate` runs it;
- **Class** – which code it tests (Table 1);
- **Expected** – the reference value (analytic, published or from another code);
- **Computed** – the value SWINGERsim computes;
- **Difference** – SWINGERsim minus expected, either relative (in %) or absolute (in the unit of
  the quantity), as defined by the benchmark;
- **Tolerance** – the largest difference that counts as agreement, in the same form;
- **Verdict** – PASS or EXPLAINED.

Terms are defined at first use. The notation follows the other SWINGERsim documents: $P$ is a
period, $a$ a semi-major axis, $e$ an eccentricity, $i$ an inclination to the sky, $i_\mathrm{mut}$
the mutual inclination of two orbits, $\omega$ the argument of periastron and $\Omega$ the longitude
of the node. The orbits of a 2+2 quadruple are called A and B (the two pairs) and "out" (their mutual
orbit); in a 3+1 they are "inner", "mid" and "out". The equations behind each module are derived in
*SWINGERsim – Physics and Methods*; this report quotes only what is needed to judge a comparison.

## Parts II–IV {#parts}

Part I asks whether the physics of SWINGERsim reproduces the literature. The other three parts ask
whether SWINGERsim, used as an astronomer would use it, reproduces real observations of systems
studied by the SQUADRA group:

- **Part II: Three published multiple stars** – the 2+2 quadruple BU CMi, the 3+1 quadruple
  TIC 120362137 and the sextuple V994 Her: their published solutions integrated forward and
  backward, synthetic light curves and eclipse times against TESS and ground-based data, and joint
  fits of all data;
- **Part III: Eclipsing binaries with a changing inclination** – V907 Sco, QX Cas and AH Cep, whose
  eclipses appear, deepen or vanish because a third body tilts their orbits;
- **Part IV: Special cases** – the 2+2 quadruples CzeV343 and CzeV4315, which test the limits of the
  models (resonances, unstable configurations, chaotic members of ensembles).
