---
title: "SWINGERsim"
subtitle: "User Guide"
---

# Introduction {#introduction}

SWINGERsim is a program for modelling **specific, observed multiple-star systems**: a close
binary with a third star, two binaries that orbit each other (a 2+2 quadruple), a triple with a
distant fourth star (3+1), and larger hierarchies of up to about ten stars, such as the triply
eclipsing sextuple V994 Her. You describe one system – the masses of the stars, how they are
paired, the orbits on which they move and the physics that matters – and SWINGERsim computes
how the system evolves and what an observer would measure: eclipses that appear and vanish,
times of minima and their O−C diagram, light curves, radial velocities, and whether the system
survives or two of its stars merge.

This User Guide explains what the software is for, how its parts fit together, how to install
and control it, and what goes in and comes out of every workflow. It does not derive the
physics; that is the subject of the companion document *Physics and Methods*.

## What SWINGERsim is for {#what-for}

A **hierarchical multiple star** is a system whose orbits have very different sizes, so that
each orbit is close to a Kepler ellipse that is slowly perturbed by the others. Every orbit of
such a system joins two components, each of which is a star or a smaller subsystem; the whole
system is a tree of pairs (Figure 1). The slow perturbations are what makes these systems
interesting for an observer: the planes of the orbits precess, so the inclination of an
eclipsing binary to the line of sight changes and its eclipses become deeper, shallower or
disappear; eccentricities oscillate; the times of minima shift because the binary moves on its
outer orbit (the light-travel-time effect, LTTE) and because the other stars pull on it.

![**Figure 1.** The hierarchies SWINGERsim describes, drawn as trees. Orange circles are stars
(named `m1`, `m2`, …), green squares are orbits (named `A`, `B`, `out`, …); every orbit joins
the two components drawn below it. The panels show a binary, a triple (2+1), the two kinds of
quadruples (2+2 and 3+1), three quintuples and sextuples that have ready-made names in the
configuration (2+2+1, 2+(2+1), 2+2+2, 3+3) and a ten-star tree written by hand (`custom`). Read
a tree from the bottom: in 3+1, `m1` and `m2` form the orbit `inner`, `inner` and `m3` form
`mid`, and `mid` and `m4` form the outermost orbit `out`.](../../figures/concepts/hierarchies.png)

Figure 1 is the key to every configuration file and every output table: the names of the
orbits in the tree (`A`, `B`, `out`, `inner`, `mid`, …) are the names under which the elements
of the orbits are given in the input and reported in the output. The presets in the top two
rows and the two sextuples have fixed names; a custom tree uses the names you write.

SWINGERsim works with one system at a time and offers three ways of asking questions about it.

| Way of working | Question it answers | Typical result |
| --- | --- | --- |
| forward model | "What does a system with these parameters do, and what would we observe?" | the evolution of the orbits over years to millions of years, eclipse windows, an O−C diagram, a light curve |
| fit of measured data | "Which parameters of this system reproduce the measured minima, radial velocities, light curves and inclinations?" | fitted elements with errors, residuals, predictions of future minima, a fitted configuration that runs like any other |
| ensemble | "Some parameters are unknown – what is possible, and what do the measurements allow?" | distributions of outcomes over thousands of drawn systems, weighted by the measurements |

**Table 1.** The three ways of working with SWINGERsim. Each row names the question, and the
right-hand column what you get. The three are connected: a fit produces a configuration that is
then run forward, and an ensemble is a set of forward models whose members are weighted by the
same kinds of measurements that a fit uses.

## Questions it answers {#questions}

The following questions are typical for the systems SWINGERsim is made for. Each is answered
by one of the workflows of chapter [Workflows](#workflows).

- When did the eclipses of this binary begin, and when will they end?
- Is the observed change of the eclipse depth explained by the precession driven by the
  known third body, or does it need another, unseen companion?
- Which outer orbits and masses of an unseen companion reproduce a measured inclination and
  its rate of change?
- Does an N-body model of the whole system reproduce the measured times of minima of all
  eclipsing pairs at once, and with which elements?
- Is the system dynamically stable, and for how long? Is it near a resonance or in a chaotic
  region?
- What does the system look like in a TESS light curve, and how precisely would the minima be
  measured?
- How do the stars and orbits change over millions of years: when does a star fill its Roche
  lobe, and what follows (mass transfer, a contact phase, a common envelope, a merger)?
- Can the wave in an O−C diagram be magnetic activity instead of a third body?

## What it does not do {#not-for}

SWINGERsim is deliberately narrow. It is not a population-synthesis code: it does not draw
whole galaxies of multiple stars and does not model the selection effects of a survey (an
ensemble can draw broad distributions, but the interpretation as a population is the user's).
It does not reduce photometry or spectra and does not measure times of minima on real data;
it starts from tables of measured quantities. Its light curves are computed from the positions
of the stars in the dynamical model (a photodynamical model) with spherical, rotating or
Roche-shaped stars; it is not a replacement for a dedicated light-curve code when the detailed
surface physics of one binary is the aim (spots, detailed atmospheres). Supernovae and the
evolution of neutron stars and black holes are outside its scope, and non-hierarchical
(democratic) systems appear only as the outcome of scattering experiments. The limits of each
piece of physics are listed in *Physics and Methods*, chapter Limitations.

## The documentation set {#documentation-set}

The documentation of SWINGERsim consists of seven documents.

| Document | Content | Read it when |
| --- | --- | --- |
| *User Guide* (this document) | what the software does, installation, the web application, the command line, workflows, configuration, outputs | you want to use the software |
| *Physics and Methods* | the equations and algorithms: hierarchies and orbits, the stars, what the observer sees, the N-body and secular dynamics, stellar evolution, tides, mass transfer, stability, observables, fitting, ensembles, limitations | you need to know what is computed and how accurately |
| *Testing, Part I: Validation against the literature* | benchmarks: published results reproduced by SWINGERsim | you want to know how far the physics is verified |
| *Testing, Part II: Three published multiple stars* | BU CMi, TIC 120362137 and V994 Her modelled from their data | you want worked applications to real systems |
| *Testing, Part III: Eclipsing binaries with a changing inclination* | V907 Sco, QX Cas and AH Cep | you study binaries whose eclipses come and go |
| *Testing, Part IV: Special cases* | questions without a published answer: the future of CzeV343 (mass transfer, instability, ZLK merger), the hierarchy of CzeV4315 (2+2 against 2+(2+1) by Bayesian evidence), and the convergence of the MCMC fits | you combine evolution, dynamics and ensembles, or compare hypotheses |
| *Open questions and future work* | what is not solved and what could come next | you plan further work |

**Table 2.** The seven documents of the SWINGERsim documentation. The middle column says what
each contains, the right column when to read it. References between documents are given in the
text by the title of the document and of its chapter.

## How to read this guide {#how-to-read}

The guide is written for astronomers and students of astrophysics; no programming beyond
editing a text file and typing a command is needed. Every technical term is explained at its
first use or in the [Glossary](#glossary).

| If you want to … | read |
| --- | --- |
| know how the program is built and which part does what | [How SWINGERsim works](#how-it-works) |
| install it | [Installation](#installation) |
| see results in a few minutes | [Quick start](#quick-start) |
| work with forms in a web browser | [The web application](#gui) |
| work from a terminal or a script | [The command line](#command-line) |
| do a concrete task (fit minima, run an ensemble, follow evolution) | [Workflows](#workflows) and [Worked examples](#examples) |
| look up a key of the configuration file | [The configuration file](#configuration) |
| understand an output file or column | [Outputs](#outputs) |
| solve a problem | [Troubleshooting and FAQ](#troubleshooting) |
| look up a term | [Glossary](#glossary) |

**Table 3.** Where to find what in this guide. The left column is a task, the right column the
chapter that covers it; the chapters can be read in any order after the first two.

Commands are written as they are typed, e.g. `swingersim run config.yaml`; configuration keys
are written with dots for nesting, e.g. `run.t_max_yr` is the key `t_max_yr` inside the section
`run`, and `orbits.<name>.P` is the period of any orbit. Units follow the conventions of the
field: masses in M☉, radii in R☉, periods in days (d), semi-major axes in au, angles in degrees
and times in years (yr) unless stated otherwise.
