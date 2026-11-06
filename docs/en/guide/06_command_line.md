# The command line {#command-line}

Everything SWINGERsim computes is started by one command, `swingersim`, followed by a
sub-command and its arguments. `swingersim --help` lists the sub-commands and
`swingersim <sub-command> --help` lists the options of one. The web application and the worker
build their forms and job command lines from the same definitions, so an option described here
is available everywhere.

## Common conventions {#cli-conventions}

- The first argument of most sub-commands is a configuration file (YAML). Relative paths inside
  a configuration (tables of minima, light curves, output folders of a fit) are relative to the
  configuration file.
- `-o` / `--output` names the output file or folder; missing folders are created.
- Every sub-command that produces results also writes an **HTML report**: one self-contained file
  with the command line, the configuration, the printed summary, interactive plots and the
  figures, which opens offline in any browser. It is written next to the output; `--html PATH`
  chooses another name and `--no-html` switches it off. `gui`, `worker`, `view`, `pulsation`,
  `applegate`, `mesa-grid` and `mesa-reduce` write no report.
- `-j` / `--processes` sets the number of parallel processes of the commands that can use
  several cores; `0` means all cores (of the job, when run by the worker).

| Sub-command | Purpose | Main input | Main output |
| --- | --- | --- | --- |
| `run` | one simulation | configuration | summary, HDF5, CSV, PNG, 3D viewer, report |
| `view` | 3D viewer of saved runs | HDF5 files or folders | HTML viewer |
| `stability` | stability criteria and the N-body ghost test | configuration | report |
| `timings` | eclipse times from the dynamics and their O−C | configuration | CSV, PNG |
| `minima` | minima measured on a synthetic light curve | configuration | CSV, PNG |
| `lightcurve` | synthetic light curve | configuration | CSV, PNG |
| `rv` | radial-velocity curves | configuration | CSV, PNG |
| `apsidal` | apsidal motion split into its parts | configuration | table |
| `pulsation` | periods of pulsating stars | configuration | printed table |
| `fit` | fit of measured minima, radial velocities, light curves, inclinations | configuration with `fit` | parameters, residuals, fitted configuration, predictions |
| `ensemble` | Monte Carlo ensemble | configuration with `ensemble` | folder of rows and summary |
| `ensemble-report` | statistics of a finished ensemble, weighted by measurements | ensemble folder | report, figure |
| `emulator-train`, `emulator-report` | neural emulators of an ensemble | ensemble folder | network, posterior summary |
| `sbi` | simulation-based inference on an ensemble | ensemble folder | posterior samples, corner plot |
| `chaosmap` | MEGNO map over one or two parameters | configuration | map, table |
| `scatter` | few-body scattering experiments | scattering configuration | outcome table |
| `escape-theory` | statistical break-up of a three-body system | numbers | outcome distributions |
| `applegate` | can magnetic activity explain an O−C wave? | numbers | energy ratio |
| `validate` | benchmarks against published results | benchmark ids | verdicts, report |
| `mesa-grid`, `mesa-reduce` | own stellar-model grids with MESA | grid file | MESA runs, reduced grid |
| `worker` | runs job files from a folder | job folder | results of the jobs |
| `gui` | the web application | SWINGERsim folder | – |

**Table 17.** All sub-commands. Columns: the name, what it is for, what it reads and what it
writes. The sections below give the options and an example of each.

## Simulations {#cli-sims}

### `swingersim run` {#cli-run}

Integrates one configuration over `run.t_max_yr`, prints a summary and, with `-o`, writes the
run file and its standard products (chapter [Outputs](#out-run)).

| Option | Meaning |
| --- | --- |
| `config` | the configuration file |
| `-o`, `--output` | HDF5 output file; the products are written next to it |
| `--backend nbody` or `secular` | replace `run.backend` of the configuration |
| `--t-max-yr T` | replace `run.t_max_yr`; a negative value integrates into the past |
| `--no-viewer` | do not write the 3D viewer `<output>_viewer.html` |
| `--no-plots` | do not write the PNG figures `<output>_*.png` |

**Table 18.** Options of `swingersim run`.

```bash
swingersim run configs/examples/example_2p2.yaml -o output/example.h5
swingersim run configs/examples/example_2p2.yaml --backend secular --t-max-yr -1000 -o output/past.h5
```

### `swingersim view` {#cli-view}

Writes the interactive 3D viewer of saved runs: one run, several, or every run in folders.
`run -o` already writes it; `view` rebuilds it, e.g. with more animation frames.

| Option | Meaning |
| --- | --- |
| `h5 …` | HDF5 files of `run`, or folders searched for them |
| `-o`, `--output` | HTML file of a single run (default `<run>_viewer.html`) |
| `--frames N` | number of animation frames (default 150) |
| `--force` | folders: rebuild every viewer (default: only missing ones and those older than their run) |
| `--cdn` | load plotly.js from the internet: a smaller file that needs a connection |

**Table 19.** Options of `swingersim view`.

```bash
swingersim view output/example.h5 --frames 300
swingersim view output                # every run under output/ that has no up-to-date viewer
```

### `swingersim stability` {#cli-stability}

Evaluates the stability of every nested triple of the hierarchy with the criteria of Mardling &
Aarseth (2001) and Vynatheya et al. (2022), the neural classifiers of Vynatheya et al. (2023)
for every sub-quadruple when they are installed, and runs the N-body ghost test of Vynatheya et
al. (2023) (*Physics and Methods*, chapter Stability and fates).

| Option | Meaning |
| --- | --- |
| `--orbits N` | outer orbits of the N-body test (default 100) |
| `--no-nbody` | skip the N-body test |

**Table 20.** Options of `swingersim stability`.

```bash
swingersim stability configs/examples/example_2p2_zlk.yaml
```

## Observables {#cli-observables}

### `swingersim timings` {#cli-timings}

Finds every eclipse of the inner binaries in an N-body integration (the conjunction of the two
stars of a pair, with the light-travel time of the pair's centre of mass) and writes the times
and their O−C against a linear ephemeris fitted to the primary eclipses.

| Option | Meaning |
| --- | --- |
| `--years Y` | end of the series [yr after $t = 0$] (required) |
| `--from-years Y0` | start of the series [yr, ≤ 0]; negative values add the eclipses before $t = 0$, integrated into the past |
| `--bjd-t0 BJD` | BJD of model time $t = 0$: adds the column `bjd_obs` |
| `-o`, `--output` | CSV file (required) |
| `--binaries …` | orbit names (default: all binaries of two stars) |
| `--samples N` | samples per binary period (default 40) |
| `--plot` | PNG of the O−C diagram |

**Table 21.** Options of `swingersim timings`.

```bash
swingersim timings configs/examples/example_timings_czev4315.yaml --years 5 \
  -o output/timings_czev4315.csv --plot output/timings_czev4315.png
swingersim timings configs/squadra_test3/v907_sco_published.yaml --from-years -10 --years 10 \
  --bjd-t0 2459360.0 -o output/v907_timings.csv
```

### `swingersim lightcurve` {#cli-lightcurve}

Computes the flux of all stars in a passband from the N-body positions: limb darkening, mutual
eclipses of every pair of stars, ellipsoidal variation, Doppler beaming, reflection, the
light-travel time, pulsations, rotating and Roche-shaped stars and dust discs as configured
(*Physics and Methods*, chapter Observables).

| Option | Meaning |
| --- | --- |
| `--days D` | length of the series [d] (required) |
| `--start-day S` | first time [d after $t = 0$] (default 0) |
| `--cadence-min C` | sampling [min] (default 2) |
| `--band B` | passband: `U`, `B`, `V`, `R`, `I`, `TESS`, `Gaia_G`, `Gaia_BP`, `Gaia_RP` (default `TESS`) |
| `--exposure-min E`, `--supersample N` | exposure time [min] and sub-exposures per point, for the smearing of the eclipses |
| `--noise-ppm N`, `--seed S` | white noise [ppm] and its random seed |
| `--ld FILE` | YAML file: star → four limb-darkening coefficients, instead of the tables |
| `--no-proximity` | no ellipsoidal variation, beaming and reflection |
| `--no-light-time` | ignore the light-travel time |
| `--engine swingersim` or `phoebe` | the light of the inner binaries from SWINGERsim's own model (default) or from PHOEBE 2 |
| `--phoebe-atm`, `--phoebe-distortion`, `--phoebe-ntriangles`, `--no-irradiation` | PHOEBE settings: atmosphere (`auto`, `ck2004`, `blackbody`), shape (`roche`, `sphere`), mesh size (default 1500), no reflection |
| `-o`, `--output` | CSV file (required) |
| `--plot` | PNG of the light curve |

**Table 22.** Options of `swingersim lightcurve`.

```bash
swingersim lightcurve configs/squadra_test3/v907_sco_published.yaml --days 30 --cadence-min 10 \
  -o output/v907_lc.csv --plot output/v907_lc.png
```

### `swingersim minima` {#cli-minima}

Computes a synthetic light curve and measures every minimum on it with the phenomenological
eclipse profile of Mikulášek (2015), with an error, as an observer would measure real data. It
takes the light-curve options of Table 22 (`--band`, `--cadence-min`, `--noise-ppm`, `--seed`,
`--no-proximity`, `--no-light-time`, the PHOEBE options) and:

| Option | Meaning |
| --- | --- |
| `--years Y` | length of the series [yr] (required) |
| `-o`, `--output` | CSV file of the minima (required) |
| `--binaries …`, `--samples N` | orbit names; samples per binary period |
| `--window W` | width of the fitted window in eclipse durations (default 2.5) |
| `--season-days S`, `--seasons FILE` | fit one eclipse shape per season of this length and write the shapes to a CSV |
| `--blended` | measure on the total light (default: each binary on its own light plus third light) |
| `--bjd-t0 BJD` | BJD_TDB of $t = 0$ for the time column (default `fit.stage2.bjd_t0`, else J2000) |
| `--plot` | PNG of the O−C of the minima and of the conjunctions |

**Table 23.** Options of `swingersim minima` beyond those of `lightcurve`.

```bash
swingersim minima configs/benchmarks/qzcar_minima.yaml --years 2 --binaries A \
  -o output/qzcar_minima.csv --plot output/qzcar_minima.png
```

### `swingersim rv` {#cli-rv}

Radial-velocity curves of the stars from the N-body model, at the emission times, with measured
velocities for comparison.

| Option | Meaning |
| --- | --- |
| `--start`, `--end` | first and last time [BJD_TDB] (required) |
| `--bjd-t0` | BJD_TDB of $t = 0$ (default `fit.stage2.bjd_t0`) |
| `--step` | step of the curves [d] (default 0.01) |
| `--stars …` | stars (default all) |
| `--gamma` | systemic velocity [km/s] (default `fit.rv.gamma`) |
| `--measured FILE` | CSV of measured velocities (default `fit.rv.file`) |
| `-o`, `--output` | output folder (default `output/rv`) |

**Table 24.** Options of `swingersim rv`.

```bash
swingersim rv configs/squadra_test3/v907_sco_fit_lc_v2.yaml --start 2459700 --end 2459710 -o output/v907_rv
```

### `swingersim apsidal` {#cli-apsidal}

The apsidal motion of every binary of two stars split into its contributions: the tidal and the
rotational distortion of both stars after Sterne (1939) for synchronous, pseudo-synchronous and
the configured rotation, and the relativistic advance; optionally the total advance measured in
an N-body run, which adds the perturbations by the other stars; and the $k_2$ that a target
apsidal period would need.

| Option | Meaning |
| --- | --- |
| `--target A=194 B=98` | target apsidal periods [yr] per pair |
| `--nbody-years Y` | also run the N-body model for Y years and measure the total advance |
| `-o`, `--output` | CSV file of the table |

**Table 25.** Options of `swingersim apsidal`.

```bash
swingersim apsidal configs/squadra_test3/v907_sco_published.yaml -o output/v907_apsidal.csv
```

The table lists per pair the period, eccentricity, the two $k_2$, the tidal advance, the
rotational advance for the three rotation assumptions, the relativistic advance, the total in
deg/yr and the apsidal period $U$ in years.

### `swingersim pulsation` {#cli-pulsation}

Prints for every star with `stars.<m>.pulsation` its type and mode, the period of the main mode,
$dP/dt$ from the evolving stellar model, its radius and $T_\mathrm{eff}$, and notes on validity.

```bash
swingersim pulsation configs/examples/example_pulsating_triple.yaml
# m1: delta_scuti (F)  P = 0.0711228 d  dP/dt = 5.35e-06 s/yr  R = 2.03 Rsun  Teff = 7931 K
# m3: lpv (F)  P = 334.992 d  dP/dt = 135 s/yr  R = 348 Rsun  Teff = 2867 K
```

## Fits of measured data {#cli-fit}

### `swingersim fit` {#cli-fit-cmd}

Fits the model to the measurements listed in the `fit` section: stage 1 is an analytic O−C model
of every eclipsing pair, stage 2 the joint N-body model of the whole system fitted to all minima
at once, optionally with radial velocities, light curves and inclinations (workflow in
[Workflows](#wf-fit), keys in [The configuration file](#cfg-fit)).

| Option | Meaning |
| --- | --- |
| `config` | configuration with a `fit` section |
| `-o`, `--output` | output folder (default `fit.output.dir`) |
| `--stage 1` or `2` | 1: analytic O−C (default without `fit.stage2`); 2: joint N-body model (default with it) |
| `--processes N` | stage 2: cores for the derivatives (0 = all; default `fit.stage2.processes`) |
| `--sampler none` or `emcee` | replace the sampler of the stage that runs |
| `--anchor` | only anchor the period and phase of every eclipsing pair to the minima and write `anchored.yaml` |
| `--bjd-t0 BJD` | with `--anchor`: BJD_TDB of $t = 0$ (default `fit.stage2.bjd_t0`) |
| `--no-plots` | do not draw the figures |

**Table 26.** Options of `swingersim fit`. `--sampler`, `--processes` and `--stage` do not
apply to `--anchor`; `--processes` applies to stage 2 only.

```bash
swingersim fit configs/examples/example_fit.yaml                       # stage 1, seconds
swingersim fit configs/examples/example_fit.yaml --sampler emcee        # stage 1 with MCMC
swingersim fit configs/squadra_test/tic120362137_fit.yaml --processes 0 # stage 2
swingersim fit configs/squadra_test/bu_cmi_fit.yaml --anchor            # anchoring only
```

### `swingersim applegate` {#cli-applegate}

Can the magnetic activity cycle of a star explain a periodic O−C variation of its binary
(Applegate 1992)? The command gives the largest O−C semi-amplitude the star can power over the
modulation period in the thin-shell model of Tian et al. (2009) and, for an observed amplitude,
the ratio of the energy needed to the energy available (below 1: possible).

| Option | Meaning |
| --- | --- |
| `--mass`, `--radius`, `--teff` | the active star [M☉, R☉, K] |
| `--a`, `--pbin` | separation [R☉] and period [d] of the binary |
| `--pmod` | modulation period [yr] |
| `--oc` | observed O−C semi-amplitude [s] (optional) |

**Table 27.** Options of `swingersim applegate`; all except `--oc` are required.

```bash
swingersim applegate --mass 0.93 --radius 0.96 --teff 5040 --a 3.3 --pbin 0.522 --pmod 30.5 --oc 137.2
# largest O - C semi-amplitude the star can power over 30.5 yr: 1145 s (0.0132 d; ...)
# observed 137.2 s: dP/P = 1.79e-06, dE/E_sec = 0.0144 -> possible
```

The example is V471 Tau with the values of Völschow et al. (2016).

## Ensembles and inference {#cli-ensembles}

### `swingersim ensemble` {#cli-ensemble}

Draws the members described in the `ensemble` section, runs each forward and/or backward in
parallel and writes one summary row per member and direction (workflow in
[Workflows](#wf-ensemble), output in [Outputs](#out-ensemble)).

| Option | Meaning |
| --- | --- |
| `-o`, `--output` | output folder (required) |
| `-j`, `--processes` | parallel processes (0 = all cores) |
| `--n N` | replace `ensemble.n` |
| `--backend` | replace `run.backend` |
| `--resume` | continue an interrupted ensemble (finished members are kept) |
| `--retry-errors` | with `--resume`: run again the members that ended with an error |

**Table 28.** Options of `swingersim ensemble`.

```bash
swingersim ensemble configs/examples/example_v994_ensemble.yaml -o output/v994 --n 40 -j 0
swingersim ensemble configs/examples/example_v994_ensemble.yaml -o output/v994 --n 400 -j 0 --resume
```

### `swingersim ensemble-report` {#cli-ensemble-report}

Statistics of a finished ensemble: the fractions of the fates, the distribution of the largest
eccentricities, mutual inclinations, eclipses, rates of change of the inclinations and the
stability verdicts, separately forward and backward, with 68 % Wilson intervals; with
observations, every statistic is weighted and the drawn parameters are summarised before and
after the conditioning.

| Option | Meaning |
| --- | --- |
| `directory` | the ensemble folder |
| `--observed FILE` | YAML file: summary column → `{value, sigma}` or limits; replaces `ensemble.observed` |
| `--prior` | ignore the observations (unconditioned report) |
| `--thresholds …` | report the fraction of members with $1-e_\mathrm{max}$ below these values (default 0.1, 0.01, 0.001) |
| `--plot FILE` | an overview figure |
| `--x`, `--y` | drawn parameter and summary column of the last panel of the figure |

**Table 29.** Options of `swingersim ensemble-report`.

```bash
swingersim ensemble-report output/v994 --plot output/v994/overview.png --x orbits.out.i --y B_didt0_deg_yr
```

### `swingersim emulator-train` and `emulator-report` {#cli-emulator}

`emulator-train` trains a neural emulator (a deep ensemble of multi-layer perceptrons of
scikit-learn) that predicts summary columns of an ensemble from its drawn parameters, or a
classifier of a yes/no column, and prints its scores on held-out members. `emulator-report` draws
the emulator inputs from the distributions of `ensemble.vary`, weights them by
`ensemble.observed` and prints the posterior of every input (median, 68 % and 95 % intervals) –
with millions of draws at negligible cost.

| Option | Meaning |
| --- | --- |
| `ensemble`, `-o` | the ensemble folder and the model file (`.pkl`) |
| `--inputs …` | input columns (drawn parameters) |
| `--targets …`, `--classes …` | regression columns; binary columns or `column==value` (e.g. `fate==intact`) |
| `--log-inputs`, `--log-targets` | use $\log_{10}$ |
| `--angle-inputs` | angles in degrees, used as cos and sin |
| `--scale-targets target=input^power` | fit the target divided by a known power law |
| `--hidden`, `--models`, `--activation`, `--alpha`, `--learning-rate`, `--batch-size`, `--max-iter`, `--test-fraction`, `--seed` | network settings (defaults: four layers of 50 units, 5 networks, ReLU, L2 $10^{-4}$, rate $10^{-3}$, 2000 iterations, 20 % held out) |
| `--plot` | PNG of the training and held-out diagnostics |
| `emulator-report model config`, `--n`, `--seed` | the model, the ensemble configuration, prior draws (default $10^6$) |

**Table 30.** Options of the emulator commands.

```bash
swingersim emulator-train output/v994 -o output/v994/emulator.pkl \
  --inputs orbits.out.i orbits.out.Omega orbits.B.Omega \
  --angle-inputs orbits.out.i orbits.out.Omega orbits.B.Omega \
  --targets B_didt0_deg_yr --plot output/v994/emulator.png
swingersim emulator-report output/v994/emulator.pkl my_v994_with_observed.yaml
```

The inputs must be columns that the ensemble drew (`ensemble.vary`); `emulator-report` needs a
configuration whose `ensemble` section also has `observed` entries. An emulator is only as good
as its training set: trained on the 40 members of the first `swingersim ensemble` example above,
the held-out score of this example is $R^2 = -0.25$, i.e. worse than predicting the mean; it
needs thousands of members.

### `swingersim sbi` {#cli-sbi}

Simulation-based inference: a neural posterior estimator (the package `sbi`) is trained on a
finished ensemble and returns the posterior of chosen parameters for measured values; a trained
model is reused for new measurements with the same errors.

| Option | Meaning |
| --- | --- |
| `directory` | the ensemble folder (training set) |
| `--observed FILE` | YAML: column → `{value, sigma}` and/or `{min, max}` |
| `--params …`, `--log …` | parameters to infer; those inferred in $\log_{10}$ |
| `--model FILE`, `--retrain` | trained model file (reused if it matches); train anew |
| `--epochs`, `-n` | maximum training epochs (2000); posterior samples (20 000) |
| `--save-samples`, `--plot` | CSV of the samples; PNG corner plot |
| `--coverage`, `--diagnostics` | calibration test on 20 % held-out members; PNG of loss, coverage and calibration |

**Table 31.** Options of `swingersim sbi` (needs the `sbi` group, Table 8).

```bash
cat > output/v994/observed.yaml <<'EOF'
B_didt0_deg_yr: {value: 0.05, sigma: 0.02}
EOF
swingersim sbi output/v994 --observed output/v994/observed.yaml \
  --params orbits.out.i orbits.out.Omega --model output/v994/sbi_model.pt \
  --plot output/v994/corner.png --coverage
```

The observed value in this example is illustrative; a posterior is only as good as the number of
members the network was trained on (*Physics and Methods*, chapter Ensembles and inference).

## Stability, chaos and scattering {#cli-dynamics}

### `swingersim chaosmap` {#cli-chaosmap}

Runs one N-body integration per point of a grid of one or two configuration parameters and
computes the chaos indicator MEGNO (Mean Exponential Growth factor of Nearby Orbits; 2 for
regular motion, larger for chaos).

| Option | Meaning |
| --- | --- |
| `-o`, `--output` | output folder (required) |
| `--x PATH`, `--x-range MIN MAX`, `--nx N`, `--log-x` | first parameter (e.g. `orbits.out.P`), its range, points (30), logarithmic spacing |
| `--y`, `--y-range`, `--ny`, `--log-y` | optional second parameter |
| `--orbits N` | length in periods of the widest orbit (default 1000) |
| `--t-max T` | length in years (replaces `--orbits`) |
| `--max-wall S` | seconds per grid point (default 300) |
| `-j`, `--resume` | parallel processes; continue an interrupted map |

**Table 32.** Options of `swingersim chaosmap`.

```bash
swingersim chaosmap configs/validation/quarles18_circumbinary.yaml -o output/chaosmap_cb \
  --x orbits.out.a --x-range 1.8 4.0 --nx 23 --y orbits.out.e --y-range 0 0.5 --ny 11 --t-max 3000 -j 0
```

### `swingersim scatter` {#cli-scatter}

Few-body scattering experiments: a binary or a single star meets another binary or single star;
the encounters are integrated until the outcome is clear (flyby, exchange, ionisation, a bound
triple or quadruple, a collision). The configuration has its own section `scatter`
([The configuration file](#cfg-scatter)).

| Option | Meaning |
| --- | --- |
| `-o`, `--output` | output folder (required) |
| `-j`, `--processes`, `--n`, `--resume` | parallel processes, number of encounters (replaces `scatter.n`), continue |

**Table 33.** Options of `swingersim scatter`.

```bash
swingersim scatter configs/examples/example_scatter_2p2.yaml -o output/scatter_2p2 --n 200 -j 0
```

### `swingersim escape-theory` {#cli-escape}

The statistical theory of the break-up of a non-hierarchical three-body system (Stone & Leigh
2019): which star escapes and the eccentricity of the binary left behind, for given masses,
binary, velocity at infinity and impact parameter (or angular momentum); optionally compared with
a `scatter` folder of binary–single runs.

```bash
swingersim escape-theory --masses 1 1 1 --a 1 --v-inf 1 --b 2
```

Options: `--masses m1 m2 m3` (binary star 1, binary star 2, single star, M☉), `--a` [au], `--e`,
`--v-inf` [km/s], `--b` [au] or `--L0-ratio`, `--alpha`, `--criterion AE|SE`, `--compare FOLDER`.

## Validation and own stellar models {#cli-other}

### `swingersim validate` {#cli-validate}

Reproduces published results and prints PASS, EXPLAINED or FAIL for each benchmark; `--all`
includes the slow ones (minutes each); benchmark ids select some. The report is written to
`output/validation/`. The benchmarks are the subject of *Testing, Part I*.

```bash
swingersim validate            # the fast benchmarks
swingersim validate --all      # all of them
```

### `swingersim mesa-grid` and `mesa-reduce` {#cli-mesa}

Compute an own grid of stellar models with MESA from a grid file (`configs/mesa/*.yaml`) and
reduce its profiles to the quantities SWINGERsim needs ([Installation](#install-mesa)). Options:
`--root` (the SWINGERsim folder), `-j` (parallel MESA runs), `--test` (one 1 M☉ track only).

## Running and controlling {#cli-control}

### `swingersim worker` {#cli-worker}

Runs job files dropped into a folder ([Installation](#install-worker), [Workflows](#wf-worker)).

| Option | Meaning |
| --- | --- |
| `jobdir` | job folder (default `<root>/jobs`) |
| `--root` | the SWINGERsim folder (default `.`) |
| `--poll S` | seconds between checks (default 10) |
| `--once` | run the waiting jobs and stop |
| `--cores N` | cores the jobs may use together (default all) |

**Table 34.** Options of `swingersim worker`.

### `swingersim gui` {#cli-gui}

Starts the web application; its options are in Table 11.
