# How SWINGERsim works {#how-it-works}

This chapter describes what happens between a configuration file and the finished results: the
processing chain, the parts of the program that carry out each step, how they are connected,
and the two ways of computing the motion (the N-body and the secular backend). It ends with a
map of all modules and a table that separates what SWINGERsim computes with its own code from
what it takes from external libraries and published data.

## The processing chain {#chain}

Every computation follows the same chain (Figure 2). A **configuration** – a text file in the
YAML format – describes one system. SWINGERsim checks it, builds the initial state of the
system, integrates the motion with one of two **backends**, analyses the result with a set of
**diagnostics**, and writes **products** (tables, figures, an HTML report) and **observables**
(eclipse times, light curves, radial velocities). Ensembles, fits and maps are layers on top of
this chain: they call it many times with different parameters. The command line and the web
application start the computations and show their results.

![**Figure 2.** The processing chain of SWINGERsim. Boxes 1–5 (top) are the path of one
configuration: the input file (grey) is checked and turned into an initial state (blue,
box 2), integrated by the N-body backend (3a) or the secular backend (3b), analysed (violet,
box 4) and turned into run products and observables (orange, 5a and 5b). The green boxes at the
bottom (6a–6c) are the methods that run this chain many times – ensembles, fits of measured
data, maps and experiments; the arrows from them point into the chain they call. The yellow
box (7) is the control layer: the command line, the web application and the worker that runs
job files. Each box lists the main tasks of that step; the modules that carry them out are in
Table 5.](../../figures/guide/pipeline.png)

Figure 2 is read from left to right for one run. The arrows from box 2 split because each run
uses exactly one backend, chosen in the configuration; both backends hand the same kind of
result to the diagnostics, which is why everything to the right of them works with either.
The bottom row shows that ensembles, fits and maps are not separate programs: each of them
produces many configurations (or many variants of one) and sends them through boxes 2–5.

The steps in detail:

1. **Input.** The configuration names the hierarchy (`topology`), the masses, the Keplerian
   elements of every orbit at the start ($t = 0$), optional properties of the stars, the
   physics to include and the length of the run. Optional sections describe an ensemble
   (`ensemble`), a fit of measured data (`fit`) and occulting dust discs (`discs`). Measured
   data (times of minima, radial velocities, light curves, inclinations) are CSV tables that
   the configuration refers to. The complete reference is the chapter
   [The configuration file](#configuration).
2. **Check and build.** The file is validated against a schema before anything is computed:
   every key must be known, every value must have the right type and range, the orbits must
   match the hierarchy, and combinations that do not make sense are refused with a message
   naming the key. The hierarchy becomes a binary tree; the elements of every orbit are
   converted to positions and velocities in Jacobi coordinates (each orbit is the relative
   motion of the centres of mass of its two components). Every star receives a radius, an
   effective temperature, an apsidal-motion constant $k_2$, a moment of inertia and, if
   requested, a tidal time lag, from a stellar model at the age of the system or from values
   given in the file. With mass transfer switched on, the history of every interacting binary is
   computed from the stellar tracks at this point (see *Physics and Methods*, chapter Mass
   transfer, common envelope and mergers).
3. **Backend.** The N-body backend integrates the positions and velocities of all stars; the
   secular backend integrates orbit-averaged elements. Both are described
   [below](#backends). With evolving stars the run is cut into segments within which the stars
   are held fixed, and the stars are updated between segments (operator splitting).
4. **Diagnostics.** The time series of the elements are turned into numbers that answer
   questions: which stars are bound to which at every output (and whether the hierarchy
   changed), the ranges of eccentricity and mutual inclination compared with the analytic
   theory of von Zeipel–Lidov–Kozai (ZLK) cycles, the precession and apsidal periods, the
   fraction of time each binary eclipses, the largest Roche-lobe filling, and the stability
   ratios of every nested triple.
5. **Products and observables.** A run writes its summary, an HDF5 file with all time series,
   a CSV table, PNG figures, an interactive 3D viewer and an HTML report. Observables are
   computed from the same dynamics on request: eclipse times and the O−C diagram, synthetic
   light curves, times of minima measured on them, radial velocities and pulsation periods.
6. **Layers.** Ensembles draw unknown parameters and run every member; fits adjust free
   parameters until the model reproduces measured data; chaos maps and scattering experiments
   run grids or random sets of configurations.
7. **Control.** The command `swingersim` with one sub-command per task, the web application
   (`swingersim gui`) and the worker, which runs job files on a computer that stays on for
   long computations.

## Two backends {#backends}

The motion of a hierarchical system has two very different time scales: the revolutions of
the orbits (days to years) and their slow secular changes (precession, ZLK cycles; decades to
millions of years). SWINGERsim offers two backends for this reason.

| | N-body (`run.backend: nbody`) | secular (`run.backend: secular`) |
| --- | --- | --- |
| what is integrated | positions and velocities of all stars | orbit-averaged elements (angular-momentum and eccentricity vectors of every orbit) |
| integrator | IAS15 of REBOUND, adaptive, error at machine precision | DOP853 of SciPy; the 2+2 equations compiled with Numba when it is installed |
| relativity | 1PN for all pairs of bodies (REBOUNDx `gr_full`) | periastron advance of every orbit |
| stars | tides and rotation as forces on extended bodies (REBOUNDx `tides_spin`), tidal dissipation, mass change | distortion as apsidal and nodal precession, tidal dissipation, evolving stars, mass loss and mass transfer |
| events | collisions, changes of the hierarchy, escapes, Roche-lobe overflow, energy-error limit, wall-clock limit | Roche-lobe overflow and mass-transfer events; the hierarchy cannot change |
| outputs | osculating elements | mean elements |
| observables (eclipse times, light curves, radial velocities) and fits | yes | no |
| cost | seconds to hours, set by the number of revolutions of the shortest orbit | about a second per run, independent of the inner periods |
| use it for | decades to thousands of years; mild hierarchies; close encounters; instability; everything an observer measures | millions of years; ensembles of thousands of members; strongly hierarchical systems |

**Table 4.** The two backends. Each row compares one property; the last row says when to use
which. "Osculating" elements are the instantaneous Kepler ellipse computed from positions and
velocities; "mean" elements are averaged over the orbital motion and lack the short-period
wiggles.

Both backends return the same structure – time series of the elements of every orbit, mutual
inclinations, energy and angular-momentum errors, events – so every diagnostic, output file
and result page works with either, and a result can always be checked by running it with the
other backend. Figure 3 shows the difference on one system.

![**Figure 3.** The same hierarchical system computed with both backends. Horizontal axis: time
in years; vertical axis: eccentricity of the inner binary. The thin grey line is the N-body
backend, the thick blue line the secular backend. Both show two ZLK cycles in which the
eccentricity grows from 0.05 to 0.77 and returns. The N-body curve carries small fast wiggles –
the osculating elements respond to every revolution of the outer orbit – which the averaged
equations leave out; the slow cycle, which matters for the long-term evolution, is the same in
both.](../../figures/concepts/nbody_vs_secular.png)

Figure 3 illustrates the rule of thumb: when the period ratios are large and the system stays
hierarchical, the averaged equations follow the slow evolution closely at a tiny fraction of
the cost. They lose accuracy when the hierarchy is mild (the outer period only a few tens of
inner periods), near mean-motion resonances and in chaotic systems; the summary of a run then
warns that the "secular estimates" are approximate (the parameter $\varepsilon_\mathrm{SA}$,
see [Outputs](#outputs)). The secular equations of a 2+2 quadruple follow Vokrouhlický (2016);
every other hierarchy uses the general formulation for nested binaries of Hamers et al. (2015),
of which the 2+2 equations are a special case. The equations, their validity and the tests
against the N-body backend are in *Physics and Methods*, chapter Dynamics: N-body, relativity
and the secular backend.

**Which backend to choose.** Start with N-body for anything an observer measures or for runs of
up to a few thousand years. Use the secular backend for evolution over millions of years, for
long ensembles, and when the inner binaries are so short that N-body would need billions of
steps. When in doubt, run both over a short time and compare the summaries.

## How the parts are connected {#connections}

The parts communicate only through files and one internal result structure, which keeps them
independent:

- **Configuration → everything.** Every command reads the same YAML file and validates it with
  the same schema. A configuration written in the web application, by hand or by a fit
  (`stage2_fitted.yaml`) is the same kind of file.
- **Backend → diagnostics → files.** Both backends return the same result structure; the
  diagnostics read only this structure; the HDF5 file stores the result together with the
  complete configuration, so a run file always says how it was made.
- **Observables from the dynamics.** Eclipse times, light curves, radial velocities and the
  fitted model observables all integrate the same N-body equations as a run; a light curve is
  therefore consistent with the O−C diagram of the same configuration.
- **Layers call the chain.** An ensemble writes one summary row per member and direction; a
  fit calls the N-body model once per evaluation of its free parameters; a chaos map runs one
  integration per grid point.
- **Jobs.** The web application does not compute anything itself: it writes job files (small
  YAML files that name a command, a configuration and options) into a folder, and the worker
  executes them with the same command line that a user would type. Results are read back from
  the output files.

## Module map {#module-map}

The program is a Python package (`swingersim`, source in `src/swingersim/`). Table 5 lists
every module, grouped by its role, with one line on what it does and whether it is SWINGERsim's
own code or a wrapper around an external library. "Own" means that the method is implemented
in SWINGERsim (it may use NumPy and SciPy for arrays, root finding and integrators);
"wrapper" means that the essential computation is done by the named external package and the
module prepares its input and reads its output; "own + external" means that SWINGERsim's own
algorithm uses an external package for one well-defined part.

| Module | What it does | Own code or wrapper |
| --- | --- | --- |
| **Configuration and set-up** | | |
| `config.py` | schema of the configuration, units, presets of hierarchies, all consistency checks | own (schema written with pydantic) |
| `units.py` | internal units (M☉, au, yr) and conversions | own |
| `kepler.py` | Kepler's equation, elements ↔ position and velocity | own |
| `hierarchy.py` | hierarchy tree, Jacobi initialisation, elements of every node | own |
| `tables.py`, `rootpath.py` | one rule for reading the CSV tables and for paths of jobs | own |
| **Stars** | | |
| `stellar.py` | radius, $T_\mathrm{eff}$, $k_2$, moment of inertia, tidal lag of every star; ZAMS relations, polytropes, Claret tables, user tables, Roche-lobe radius | own; data of Tout et al. (1996), Claret (2004, 2023) |
| `evolution.py` | SSE tracks (radius, luminosity, envelope) at any mass, age and metallicity by interpolation | own interpolation; grid computed from the published SSE routines of Hurley et al. (2000) |
| `mist.py` | MIST tracks from the pre-main sequence to the end of the AGB, interpolated at equivalent evolutionary points | own interpolation; data of Choi et al. (2016) |
| `gridclip.py` | one policy for values outside the grids of tabulated models | own |
| `massloss.py` | stellar winds in the dynamics: masses from tracks or rates, the mass-loss index Ψ | own |
| `masstransfer.py` | Roche-lobe overflow, stability, accretion, orbit response, rejuvenation, contact binaries, transfer from a tertiary | own (prescriptions of Hurley et al. 2002 and others) |
| `commonenv.py` | outcome of a common envelope (energy formalism) and the response of the outer orbits | own |
| `merger.py` | continuation of a run after two stars merge | own |
| `tides.py`, `apsidal.py` | analytic apsidal motion (tidal, rotational, relativistic) and tidal time scales | own |
| `inertial_waves.py` | dissipation of inertial waves in convective envelopes | own |
| `magnetic.py` | magnetic braking and the Applegate energy test | own |
| `pulsation.py` | pulsation periods from the stellar models and the flux modulation | own |
| `mesa.py` | own stellar-model grids: writes MESA input, runs MESA, reduces profiles to $k_2$, envelopes and dissipation | wrapper of MESA (external, installed separately) + own reduction |
| `gw.py`, `environment.py`, `tidalcapture.py`, `pdot.py` | gravitational-wave decay, Galactic tide and passing stars, tidal capture in encounters, an imposed period change | own |
| **Dynamics** | | |
| `backends/nbody.py` | the N-body backend: sets up the simulation, events, extremes between outputs, spins | own + external (REBOUND, REBOUNDx) |
| `backends/secular.py`, `backends/secular_nested.py`, `backends/secular_common.py` | the secular backends for 2+2 and for any nested hierarchy | own (SciPy integrator) |
| `backends/evolving.py`, `backends/chunked.py` | evolving stars and mass loss by operator splitting; integration in chunks that fit into memory | own |
| `secular/twotwo.py`, `secular/kernel.py`, `secular/nested.py`, `secular/quadrupole.py` | the secular equations (2+2, any nested tree), their compiled kernel, quadrupole predictions | own (kernel compiled by Numba if installed) |
| **Diagnostics and stability** | | |
| `diagnostics/hierarchy.py`, `zlk.py`, `binary.py`, `extremes.py`, `orbits.py` | hierarchy detection, ZLK and precession, apsidal motion and eclipse fractions, extremes, mutual inclinations | own |
| `stability.py` | Mardling–Aarseth and Vynatheya et al. (2022) criteria, the N-body ghost test, the optional neural classifiers of Vynatheya et al. (2023) | own; classifiers are the authors' published networks (optional) |
| `chaosmap.py` | MEGNO maps over one or two parameters | own driver; MEGNO computed by REBOUND |
| `scatter.py`, `escape_theory.py` | few-body scattering experiments; statistical theory of three-body break-up (Stone & Leigh 2019) | own (scattering integrated with REBOUND) |
| **Observables** | | |
| `timing.py` | eclipse times from the dynamics and their O−C | own |
| `photometry.py` | passbands, limb darkening, occultation of limb-darkened discs, synthetic light curves | own; limb-darkening tables of Claret & Bloemen (2011) and Reeve & Howarth (2016) |
| `roche.py`, `mesh.py`, `shapes.py` | Roche geometry, light of distorted stars on surface meshes, choice of the star shapes | own |
| `discs.py` | dust discs and rings as occulters, their precession, thermal and scattered light | own |
| `phenom.py` | phenomenological eclipse profiles (Mikulášek 2015) and times of minima measured on synthetic light curves | own |
| `radvel.py`, `rvlight.py` | radial velocities of the centres of mass and of the centres of light | own |
| `phoebe_lc.py` | the light of the inner binaries computed by PHOEBE 2 in a hybrid with the N-body positions | wrapper of PHOEBE 2 (optional) |
| `systemview.py` | the model system as data for pictures | own |
| **Fitting of measured data** | | |
| `ocfit.py` | stage 1: analytic O−C model of every pair, least squares, MCMC, predictions | own + external (SciPy least squares, emcee, astropy time conversion) |
| `dynfit.py` | stage 2: joint N-body model of all pairs, error scaling, outliers, priors, anchoring, MCMC, predictions | own + external (SciPy least squares, emcee) |
| `obsmodel.py`, `lcfit.py`, `incfit.py` | model observables of one integration; light curves and inclinations in the joint fit | own |
| **Ensembles and machine learning** | | |
| `ensemble.py`, `parallel.py` | drawing members, running them in parallel, summary rows, weights by measurements, reports | own |
| `emulator.py` | neural emulators and classifiers trained on an ensemble | own features and evaluation; networks of scikit-learn |
| `sbi_fit.py` | neural posterior estimation on an ensemble | wrapper of `sbi` and PyTorch |
| **Products, control and display** | | |
| `io.py`, `runproducts.py` | HDF5 file of a run; summary, CSV, PNG figures next to it | own (h5py, matplotlib) |
| `viewer.py`, `report.py` | interactive 3D viewer of a run; self-contained HTML report of every command | own (plotly for the plots) |
| `cli.py`, `__main__.py` | the command line | own |
| `worker.py` | the job worker | own |
| `gui/server.py`, `gui/digest.py`, `gui/assist.py`, `gui/modelviews.py`, `gui/app.html`, `gui/help.yaml` | the web application: server, result pages, help texts and quick checks, model curves, the page itself | own (Python standard-library web server; plotly.js in the browser) |
| `cnoise.py` | filter of messages that compiled libraries print | own |
| **Validation** | | |
| `validation.py` and its companions, `_benchmark.py` | the benchmarks run by `swingersim validate` | own |

**Table 5.** Module map of SWINGERsim. Bold rows are the groups. The middle column states the
task of each module, the right column whether the computation is SWINGERsim's own code, a
wrapper, or own code with one external part (named in brackets). Paths are relative to
`src/swingersim/`.

The package contains about 45 000 lines of Python plus the web page. Its own code covers the
physics that couples dynamics, stars and observables; the external libraries are used for
well-defined numerical tasks.

## What SWINGERsim implements and what it uses {#own-vs-external}

Table 6 separates the two explicitly. It answers the question a referee asks first: which
results depend on code written for SWINGERsim, and which on established external software and
published data.

| Task | Implemented in SWINGERsim | External software or data used |
| --- | --- | --- |
| hierarchies | builder of any tree of pairs, Jacobi initialisation, detection of the hierarchy at every output and of its changes | – |
| N-body integration | set-up, events (collisions, Roche overflow, escapes, wall clock), sampling of extremes between outputs, coupling to evolving stars and mass change, spins | REBOUND with IAS15 (Rein & Liu 2012; Rein & Spiegel 2015); REBOUNDx forces `gr_full`, `tides_spin` and operator `modify_mass` (Tamayo et al. 2020; Lu et al. 2023) |
| secular integration | the equations for 2+2 (Vokrouhlický 2016) and for any nested tree (Hamers et al. 2015), tides, gas damping, Galactic tide and passing stars, chunked integration | SciPy integrator DOP853; Numba compiler (optional) |
| stars | interpolation of tracks, structure constants, envelopes, tidal lags, rejuvenation, pulsation periods | tracks: MIST (Choi et al. 2016), SSE grid from the routines of Hurley et al. (2000); structure: Claret (2004, 2023); ZAMS: Tout et al. (1996); optional own grids computed with MESA (Paxton et al. 2011; Jermyn et al. 2023) |
| interacting binaries | coupling of mass transfer, contact phases, common envelopes, tertiary transfer and mergers to the dynamics | prescriptions from the literature (Hurley et al. 2002 and others; *Physics and Methods*) |
| eclipses and light curves | eclipse geometry, occultation of limb-darkened discs, ellipsoidal, beaming and reflection terms, Roche and rotating surfaces on meshes, dust discs, light-travel time, pulsations | limb-darkening tables of Claret & Bloemen (2011) and the TLUSTY coefficients of Reeve & Howarth (2016); passband curves; PHOEBE 2 (Prša et al. 2016) as an optional engine for the inner binaries |
| eclipse timing and O−C | conjunctions found in the integration, light-time, O−C, phenomenological minima (Mikulášek 2015) | – |
| radial velocities | velocities at the emission time, centre of light with the Rossiter–McLaughlin effect | – |
| fit chain | stage 1 and stage 2 models, error scaling and jitter, outlier rejection, priors, red-noise factors, anchoring, predictions with uncertainties | SciPy least squares; emcee (Foreman-Mackey et al. 2013) for MCMC; astropy for time conversions |
| ensembles | drawing, parallel runs, resumable rows, weights, effective number of members, Wilson intervals, re-weighting | pandas and pyarrow only to write Parquet (optional) |
| machine learning | features, training sets, evaluation, calibration plots | scikit-learn networks; `sbi` and PyTorch (optional) |
| stability | Mardling–Aarseth and Vynatheya et al. (2022) criteria for every nested triple, N-body ghost test, chaos maps, scattering, escape theory | MEGNO of REBOUND; the neural classifiers published by Vynatheya et al. (2023) (optional, downloaded separately) |
| configuration | schema, defaults, checks, error messages | pydantic (validation framework), PyYAML |
| outputs | HDF5 layout, products, 3D viewer, HTML reports | h5py, matplotlib, plotly |
| web application, worker | server, result pages, help, job queue, pause and stop | Python standard library; plotly.js in the browser |

**Table 6.** What SWINGERsim implements itself (middle column) and which external libraries and
published data it relies on (right column), task by task. A dash means that no external
component is involved beyond NumPy and SciPy.

## Design principles {#design}

- **Forward model first.** The core only predicts. Fits, ensemble weights and machine learning
  are layers that call the core and never change it.
- **Two backends, one result.** The N-body backend is the reference; the secular backend is
  the fast approximation; both produce the same structure, so one can check the other.
- **Any hierarchy as a tree.** Nothing in the chain assumes four stars: presets are trees,
  diagnostics walk the tree, stability is tested on every nested triple.
- **Validated physics.** Every physical ingredient is compared with a published result
  (*Testing, Part I*).
- **Reproducibility.** Every run file stores its complete configuration and the versions of the
  software that computed it; an ensemble is fixed by its seed; resuming an ensemble checks a
  fingerprint of the configuration, so rows of different set-ups are never mixed.
- **Safe remote computing.** The worker runs only SWINGERsim commands with the arguments their
  parser allows, and only with paths inside the SWINGERsim folder; a job file cannot run shell
  commands.
