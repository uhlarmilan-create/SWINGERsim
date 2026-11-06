# Workflows {#workflows}

This chapter describes the typical tasks from start to finish. Each workflow begins with a table
of its inputs, steps and outputs and ends with a short worked example. The commands are those of
chapter [The command line](#command-line); every step can also be done in the web application
(chapter [The web application](#gui)).

## Forward model of a published system {#wf-forward}

**Question.** A paper gives the periods, eclipses, masses and perhaps an outer orbit of a
multiple system. What does this system do, and what will observers see next?

| Inputs | Steps | Outputs |
| --- | --- | --- |
| published masses, radii, temperatures; elements of every orbit; the epoch they refer to; measured apsidal motion or inclination changes for comparison | 1. write the configuration; 2. check stability; 3. run with N-body and, for long times, secular; 4. compute eclipse times and a light curve; 5. compare with the measurements | run files, eclipse windows, O−C diagram, light curve, stability verdict; the differences from the measurements |

**Table 35.** The forward-model workflow. The left column lists what you need, the middle the
steps, the right what you get.

1. **Collect the parameters** with their sources: inner periods, eccentricities, inclinations,
   radii (from light-curve solutions), masses (spectroscopic if available), the outer orbit
   (from the LTTE, radial velocities or astrometry), apsidal periods and measured rates of
   change. Note the epoch of the elements: it becomes $t = 0$.
2. **Write the configuration.** Start from an example in `configs/` (a catalog of published
   systems is in `configs/catalog/`). Mark every value that is not measured with a comment
   (`# ASSUMED`, `# ESTIMATED`). Give the radii (`stars.<m>.R`) when they are known; otherwise
   zero-age main-sequence radii are used, which are too small for evolved stars. Unknown nodes:
   keep `Omega = 0` for the eclipsing pair and set the mutual inclination through `i` and
   `Omega` of the other orbits. For close binaries switch on `physics.gr`, `tides` and
   `rotation`; they set the apsidal motion and can suppress ZLK cycles.
3. **Check stability**: `swingersim stability config.yaml`. An unstable nominal system means
   that a value is wrong or the system is not what was assumed.
4. **Run** with the N-body backend over a few precession periods and look at the result page or
   the viewer. For evolution over millions of years see [below](#wf-evolution).
5. **What would observers see?** `swingersim timings` gives the O−C of every binary with the
   light-travel time and the dynamical effects; `--bjd-t0` gives each minimum its BJD, and a
   negative `--from-years` adds the minima before the epoch of the elements. `swingersim
   lightcurve` gives the light curve in a passband. If measured minima disagree with the model,
   fit them ([next section](#wf-fit)).

**Worked example.** V907 Sco is a 3.78-d eclipsing binary with a third star on a 142-d orbit;
its eclipses stopped and restarted during the 20th century. The configuration
`configs/squadra_test3/v907_sco_published.yaml` takes the solution of Zasche et al. (2023):

```bash
swingersim stability configs/squadra_test3/v907_sco_published.yaml
swingersim run configs/squadra_test3/v907_sco_published.yaml -o output/v907/run.h5
swingersim timings configs/squadra_test3/v907_sco_published.yaml --from-years -100 --years 100 \
  -o output/v907/timings.csv
```

The stability check gives ratios of 1.67 (Mardling–Aarseth) and 2.2 (Vynatheya et al. 2022) and
a stable ghost test over 100 outer orbits. The 200-yr run (27 s) shows a mutual inclination of
25.6–26.3°, a nodal precession period of 63 yr and eclipses of the inner pair during 50.5 % of
the time, in windows of about 30 years (Figures 9–11); the eclipse times show the same windows
(Figure 13). The result pages of these three jobs are the ones shown in chapter
[The web application](#gui).

## Fitting measured data {#wf-fit}

**Question.** Which parameters of the system reproduce the measured times of minima of all
eclipsing pairs – and, if available, radial velocities, light curves and inclinations measured
at different epochs?

| Inputs | Steps | Outputs |
| --- | --- | --- |
| a configuration with start values; a table of minima; optionally radial velocities, light-curve segments, inclinations; the BJD of $t = 0$; priors from published solutions | 1. prepare the tables; 2. stage 1 (analytic O−C of each pair); 3. anchoring; 4. stage 2 (joint N-body model) with error scaling, outlier rejection, priors, red noise; 5. check residuals and warnings; 6. MCMC for reliable errors; 7. predictions and the fitted configuration | fitted parameters with errors and correlations, residuals of every measurement, predicted minima with their uncertainty, `stage2_fitted.yaml` |

**Table 36.** The fitting workflow ("the fitting chain"). The steps are described below; every
option is a key of the `fit` section ([The configuration file](#cfg-fit)).

### The data {#wf-fit-data}

The **times of minima** are a CSV table with one row per minimum and the columns `pair` (the
name of the eclipsing orbit), `type` (`p` primary, `s` secondary), `time`, `sigma` (both in
days), `source` and `method`. The column names of the measured-minima files of the web
application (`binary`, `kind`, `bjd_tdb`, `sigma_d`, `ref`) are accepted as well. Times are
BJD_TDB; `fit.time_system: HJD_UTC` or `JD_UTC` converts them (with `fit.coordinates`). A
*primary* minimum is an eclipse of the first star of the pair (`m1` of `inner`), whatever its
depth; if the observer's primary is the other star, swap the types or the stars.

**Radial velocities** (`fit.rv`): columns `star`, `time` (BJD_TDB), `rv` and `sigma` (km/s) and
`source`. Each source can have its own zero point (`fit.rv.offsets`, free as
`rv.offset.<source>`); the systemic velocity is `rv.gamma`. For semi-detached and contact
binaries `centre: light` uses the flux-weighted velocity of the visible surface, including the
Rossiter–McLaughlin effect in eclipse, instead of the velocity of the centre of mass.

**Light curves** (`fit.lc`): columns `time`, `flux`, `sigma` and optionally `segment`. Each
segment gets its own scale, zero point and polynomial trend, solved linearly at every evaluation;
with `scale: false` the depths follow from the stars as given, so the light curve constrains
the inclination and radii.

**Inclinations at epochs** (`fit.inclination`): columns `mid_hjd` (or `time`), `i_deg`,
`sigma_deg` and optionally `limit` and `source`, e.g. inclinations from light-curve solutions of
archival data. Only $\lvert\cos i\rvert$ is compared, because eclipses do not tell the sense of
revolution. The column `limit` makes a row one-sided: `upper`, `lower`, `no_eclipse` (the model
must be below its eclipse limit then) or `eclipse` (it must be above).

### Stage 1: the analytic O−C of each pair {#wf-fit-stage1}

Stage 1 describes the O−C of every pair with analytic terms: the light-time effect of an outer
orbit (Irwin 1952), apsidal motion of an eccentric pair, and a quadratic term (`fit.stage1.terms`
or per pair `fit.stage1.pairs`). With `shared_outer: true` all pairs of a 2+2 share one outer
orbit (common period, epoch and eccentricity, opposite pericentres). It is fast (seconds) and
gives start values and a first look: the outer period, the amplitude and mass function, the
apsidal period. `sampler: emcee` adds an MCMC run (Foreman-Mackey et al. 2013).

### Anchoring {#wf-fit-anchor}

Stage 2 converges only from a start that numbers the minima correctly: an error of the inner
period accumulated over the span of the data must stay well below half a period. Published
elements from a photodynamical solution often drift in an N-body model with slightly different
physics. *Anchoring* adjusts only the period $P$ and the mean anomaly $M_0$ of every eclipsing
pair so that the model minima fall on the measured ones; `swingersim fit --anchor` writes the
result as `anchored.yaml`, and `fit.stage2.anchor: true` does it automatically before the fit.

### Stage 2: the joint N-body model {#wf-fit-stage2}

Stage 2 computes the times of all minima of all pairs from one N-body integration of the whole
system – the same physics as `run` and `timings` – forward and backward from $t = 0$ (the BJD
`fit.stage2.bjd_t0`), and fits the parameters listed in `fit.stage2.free` by weighted least
squares. Light-time effect, apsidal motion and dynamical perturbations are in the model at once;
no term has to be chosen. Free parameters are elements of orbits (`orbits.<o>.P`, `M0`, `e`,
`omega`, `Omega`, `i`, or `ecosw` and `esinw` for nearly circular orbits, and `Pdot`), masses
(`masses.<m>`), properties of stars (`stars.<m>.R`, `Teff`, `apsidal_k2`, `rg2`) and the
radial-velocity zero points (`rv.gamma`, `rv.offset.<source>`). The derivatives are computed by
finite differences, two N-body integrations per free parameter and iteration, in parallel on
`processes` cores. The method is described in *Physics and Methods*, chapter Fitting observed
data.

What the minima constrain: the periods, phases, eccentricities and arguments of pericentre of
all orbits and, through the light-time amplitude and the dynamical perturbations, the masses of
the components relative to each other. The inclinations of the eclipsing pairs and the absolute
masses need light curves and radial velocities; keep them fixed at the values of such a
solution, fit them together with those data, or give them priors. Without astrometry only the
difference of the nodes of two orbits is measurable: keep the node of one orbit fixed.

**Errors and jitter.** Minima from many sources rarely have realistic errors.
`error_scale: per_source` gives every source (column `source`) its own correction so that its
reduced $\chi^2$ is 1, never below its quoted errors; `global` uses one factor for all. The
correction is a factor on the errors (`error_model: scale`) or an error floor added in quadrature
(`error_model: jitter`), suited to precise timings with optimistic errors. Sources with fewer
than `min_source_minima` minima share one correction. The fit is then repeated from its solution
until the corrections settle, at most `refits` times.

**Outliers.** `clip_sigma: 5` removes, after each fit, the minima that deviate by more than five
robust standard deviations of their source (at most `max_clip_fraction` of them) and fits again;
they stay in the residual table with the status `outlier`. Alternatively a robust loss
(`loss: soft_l1`, `huber` or `cauchy`, from `f_scale` σ on) reduces their weight without removing
them.

**Priors.** `priors: {orbits.out.e: [0.47, 0.05]}` adds a Gaussian term
$((x-\mu)/\sigma)^2$ to $\chi^2$ for a free parameter, e.g. from a published light-curve or
radial-velocity solution. Priors keep weakly constrained parameters (an outer orbit seen through
one eclipsing pair, a light-time orbit covered once) from wandering to solutions that fit the
minima equally well but contradict other data. The report gives the deviation of every prior in
σ; several σ mean a real tension between the data and the prior.

**Red noise.** Residuals of a light curve are often correlated (instrumental trends, spots,
oscillations), so the formal errors are too small. `red_noise: true` on a light-curve entry
multiplies the errors of each segment by the factor β of the time-averaging method (Pont et al.
2006; Winn et al. 2008): how much more the residuals averaged over 20–120 min scatter than white
noise would. It is applied between the fits in a row, so it needs `error_scale` or
`clip_sigma`.

**MCMC.** The errors from the covariance matrix assume a straight valley of $\chi^2$. When the
report lists strongly correlated parameters ($\lvert r\rvert > 0.95$), `sampler: emcee` runs an
affine-invariant MCMC (Goodman & Weare 2010; Foreman-Mackey et al. 2013) around the solution
with `walkers` × `mcmc_steps` model evaluations; the first half is discarded. With
`mcmc_tau_target: 50` the chain is extended until it is at least 50 integrated autocorrelation
times long, as the emcee authors recommend, or reaches `mcmc_max_steps`. The fit page shows the
ratio of steps to autocorrelation time.

**Predictions.** `fit.output.predict: {years: 5}` writes the next minima of every pair with
their uncertainty: parameter sets are drawn from the MCMC samples or the covariance matrix and
each is integrated, so the errors include the correlations of the parameters.

**Checking the result.** Read, in this order: the warnings; $\chi^2$ from the start to the fitted
model and the reduced $\chi^2$; the rms of every data type; the mean residual of primary and
secondary minima (a difference of several σ means a wrong $e\cos\omega$); systematic patterns in
the residuals (a parameter that should be free is fixed, or the fit sits in a local minimum);
parameters on their bounds; deviations from priors.

**Worked example.** The fit of V907 Sco in *Testing, Part III* combines 58 TESS minima, 172
radial velocities from two sources, 2562 points of two TESS light-curve segments and 14
inclinations measured between 1900 and 2022, with 16 free parameters and 3 priors:

```bash
swingersim fit configs/squadra_test3/v907_sco_fit_lc_v2.yaml
```

Its configuration anchors the inner pair (`anchor: true`), corrects the errors per source with a
jitter, removes outliers beyond 5σ and applies red-noise factors to the light curves. The fit
lowers $\chi^2$ from 11 920 to 1052 in 4384 model evaluations; the rms of the minima is 2.6 min,
of the radial velocities 3.85 km/s and of the light curves 0.666 ppt (Figures 19–23).

## Ensembles: what we do not know {#wf-ensemble}

**Question.** Some parameters of the system are unknown (an unseen companion, the orientation of
an outer orbit, the sense of revolution). Which values are possible, what do the measurements
allow, and what follows for the future and the past?

| Inputs | Steps | Outputs |
| --- | --- | --- |
| the nominal configuration; distributions of the unknown parameters (`ensemble.vary`); measured values or limits (`ensemble.observed`); stability requirements | 1. write the `ensemble` section; 2. run a few hundred members locally, more on the worker; 3. condition on the measurements; 4. read $n_\mathrm{eff}$, the conditioned distributions and the outcomes; 5. open the best members as configurations | `summary.csv` (one row per member and direction), the report, the weighted distributions and fractions, best-fitting members |

**Table 37.** The ensemble workflow.

![**Figure 25.** The idea of an ensemble, schematic. Left: the unknown parameter (here an outer
period in days) is drawn many times from its distribution (histogram of the draws). Middle: every
member is integrated; the curves are the eccentricity of the inner binary against time in years
for a few members. Right: the distribution of the answer over all members (here the largest
eccentricity reached). The weighting by measurements described below turns the left histogram
into a posterior and the right one into a conditioned prediction.](../../figures/concepts/ensemble_idea.png)

1. **What we do not know** (Figure 25). Draw only what is unknown; keep measured values at their
   values or draw them within their errors (`normal`). Typical choices: an unknown orientation
   is `isotropic` (uniform in $\cos i$); an unknown node or phase is `angle` (uniform 0–360°); a
   scale parameter such as an outer period or a companion mass is `loguniform`; an eclipsing
   pair whose sense of revolution is unknown gets `mirror: true` ($i$ or $180° - i$).
   `require: [stable_ma01]` redraws members that a stability criterion classifies as unstable.
   Write down every choice: the result depends on it.
2. **Run** forward and backward (`directions: [forward, backward]`): the rates of change "now"
   are best measured from both sides of $t = 0$, and the past is half of the answer.
3. **Observation weights.** Every entry of `ensemble.observed` is a summary column of the
   ensemble with either a measurement `{value, sigma}`, limits `{min, max}`, or a required text
   value `{equals: intact}`. A member gets the weight
   $w = \exp(-\chi^2/2)$ with $\chi^2 = \sum ((x - \mathrm{value})/\sigma)^2$, and $w = 0$ when it
   violates a limit; the weighted members approximate the posterior distribution of the drawn
   parameters (importance weighting with the prior as the proposal). Present rates of change are
   by default *central* – the mean of the forward and backward fits – as an observer measures
   them; prefix a column with `forward:`, `backward:` or `central:` to choose.
4. **How much is left.** The effective number of members $n_\mathrm{eff} = (\sum w)^2/\sum w^2$
   (Kish 1965) says how many equally weighted members the weighted ensemble is worth. Below about
   20 the conditioned distributions are noise: run more members, or narrow the drawn ranges to
   the region the measurements allow.
5. **Evidence.** The mean weight of the members, $Z = \langle w \rangle$, is the probability of
   the measurements under the hypothesis of the ensemble (the evidence). The ratio of $Z$ of two
   ensembles compares two hypotheses (e.g. 2+2 against 2+(2+1)) – but only when the shared
   parameters are drawn from the same priors; otherwise the ratio measures the priors.
   `ensemble-report` prints $Z$ with its Monte Carlo error $\mathrm{std}(w)/\sqrt{n}$ next to
   $n_\mathrm{eff}$; when fewer than two members match, it asks you to quote the upper bound
   $Z < 1/n$. The same numbers in a few lines of Python:

   ```python
   import numpy as np
   from swingersim.config import load_config
   from swingersim.ensemble import load_rows, observation_weights

   cfg = load_config("configs/squadra_test3/qx_cas_ensemble.yaml")
   rows = load_rows("output/squadra_test3/qx_cas/ensemble")
   w = np.array(list(observation_weights(rows, cfg.ensemble.observed).values()))
   print("evidence Z =", w.mean(), " n_eff =", w.sum() ** 2 / (w ** 2).sum())
   ```

6. **Use the members.** On the result page, the members with the largest weight can be opened as
   configurations and run with a light curve or eclipse times. For a fast posterior of new
   measurements, train an emulator or a neural posterior on the ensemble (`emulator-train`,
   `sbi`).

**Worked example.** The eclipsing binary QX Cas (6.0 d) changes its inclination by
−0.23 ± 0.08 °/yr; the cause is an unseen third body. `configs/squadra_test3/qx_cas_ensemble.yaml`
draws its mass (0.3–6 M☉, log-uniform), its period (30–11 000 d, log-uniform), eccentricity
(0–0.6) and the orientation of its orbit, runs 2000 members 30 years forward and backward, and
weights them by the present inclination (82 ± 1.5°), its central rate and the survival of the
system both ways:

```bash
swingersim ensemble configs/squadra_test3/qx_cas_ensemble.yaml -o output/qx_cas_ensemble -j 0
swingersim ensemble-report output/qx_cas_ensemble
```

The conditioned ensemble has $n_\mathrm{eff} = 186.7$; it narrows the outer period from the drawn
118–4429 d (16–84 %) to 193–1028 d, leaves the mass and the angles almost unconstrained, and
predicts that in 19.5 % of the systems that eclipse now the eclipses end within 30 years
(Figures 15–18).

## Evolution over millions of years {#wf-evolution}

**Question.** How do the stars and orbits change over the life of the system: when does a star
fill its Roche lobe, and what follows?

| Inputs | Steps | Outputs |
| --- | --- | --- |
| masses, orbits, the age and metallicity of the system; the physics of the stars and their interaction | 1. choose the stellar model and age; 2. check the stars with a short run; 3. switch on evolution, tides, winds, mass transfer; 4. run with the secular backend forward or backward; 5. read the events and compare light curves at different ages | runs with evolving stars (`_stars.png`), events of the mass-transfer history, mergers, the system at later ages |

**Table 38.** The evolution workflow.

1. **Stellar model.** `stellar_model: {kind: mist, age_myr: …, feh: …}` uses the MIST tracks
   from the pre-main sequence to the end of the asymptotic giant branch (AGB); `kind: sse` the
   faster SSE formulae up to helium ignition. All stars are coeval; `stars.<m>.age_offset_myr`
   gives a star an apparent age (a stripped donor looks older, a rejuvenated gainer younger).
   A short run (small `t_max_yr`) lists $T_\mathrm{eff}$, radius, phase and age of every star in
   the summary.
2. **Physics.** `physics.evolution: true` lets the stars follow their tracks during the run;
   `tides: true` with `tidal_lag: hurley2002` takes the tidal dissipation from the envelope of
   every star; `mass_loss` adds winds; `mass_transfer` follows Roche-lobe overflow (stable
   transfer, contact binaries, common envelopes, transfer from a tertiary); `mergers` continues
   a run after two stars merge; `magnetic_braking` spins down stars with convective envelopes.
   The models are described in *Physics and Methods*, chapters Stellar models and evolution,
   Tides and rotation, and Mass transfer, common envelope and mergers.
3. **Run** with `run.backend: secular` (the N-body backend cannot follow millions of years of a
   close binary); a negative `t_max_yr` integrates back towards the birth of the stars. Without
   `mass_transfer` the run stops when a star fills its Roche lobe; it also stops at the end of
   the tracks, and the summary says why (`EVENT`).
4. **Look at the system at a later age**: a configuration with the later age and the orbits at
   that time gives its light curve and eclipse times.

Keep in mind the systematic uncertainty of the stellar tracks; when a conclusion depends on the
age, draw `stellar_model.age_myr` in an ensemble. The worked example is in
[Worked examples](#ex-evolution).

## Stability, chaos maps and scattering {#wf-stability}

**Question.** Is the system stable? Where in the space of its unknown parameters is it stable or
chaotic? Can it have formed by a capture?

| Inputs | Steps | Outputs |
| --- | --- | --- |
| a configuration; one or two parameters and their ranges; for scattering, the two colliding systems | 1. `stability` (seconds); 2. a long `run`; 3. `chaosmap` over the uncertain parameters; 4. `scatter` and `escape-theory` for formation | stability ratios and verdicts, MEGNO maps with the stability limit, outcome fractions and cross-sections |

**Table 39.** The stability workflow.

1. **Criteria.** `swingersim stability` applies the criterion of Mardling & Aarseth (2001) and
   its update by Vynatheya et al. (2022) to every nested triple (an orbit and the body that
   perturbs it); a ratio above 1 means stable. The N-body ghost test integrates the system for a
   number of outer orbits and compares it with a slightly displaced copy (Vynatheya et al. 2023).
   The criteria are approximate near 1; a long run decides.
2. **Chaos maps.** `swingersim chaosmap` integrates a grid of configurations and computes MEGNO
   (2 for regular motion, growing for chaos); escapes and collisions are marked. The map shows
   where the system sits relative to unstable regions and resonances.
3. **Scattering.** `swingersim scatter` sends a binary towards another binary or a single star
   from a range of impact parameters and records the outcomes; `escape-theory` gives the
   statistical prediction of Stone & Leigh (2019) for comparison.

**Worked example.** The chaos map of a circumbinary orbit is in [Worked examples](#ex-chaos), a
binary–binary scattering experiment in [Worked examples](#ex-scatter).

## Long computations on the worker {#wf-worker}

Ensembles of $10^4$ members, N-body ensembles, stage-2 fits with MCMC and MESA grids take hours
to days. They run as **jobs** on the worker ([Installation](#install-worker)). A job is a YAML
file in `jobs/inbox/`:

```yaml
command: ensemble                      # any sub-command except gui and worker
config: configs/squadra_test3/qx_cas_ensemble.yaml
output: output/qx_cas_ensemble
args: {n: 2000, processes: 0, resume: true}
note: QX Cas, unseen third body, 2000 members forward and backward
```

| Key | Meaning |
| --- | --- |
| `command` | the sub-command |
| `config` (and other positional arguments by their name: `ensemble`, `directory`, `h5`, `ids`, `model`) | the positional arguments of the command |
| `output` | the output argument (`-o`, or the folder of `ensemble-report` and `view`) |
| `args` | the options by their long name, with `-` written as `_` (`retry_errors: true`, `t_max_yr: 3000`) |
| `cores` | cores the job occupies (a number or `all`); default from `processes` |
| `note` | free text, shown in the web application |

**Table 40.** Keys of a job file. Every path must lie inside the SWINGERsim folder; the worker
runs only SWINGERsim commands with the arguments their parser allows.

Rules that matter in practice:

- One-core jobs (`run`, `lightcurve`, `timings`, `minima`, a stage-1 `fit`) run side by side up
  to the number of cores; a job that needs all cores runs alone. Every job receives the number of
  cores it was booked for (environment variable `SWINGERSIM_CORES`), so `processes: 0` means the
  job's cores, not the whole computer.
- Use `resume: true` for ensembles, scattering and chaos maps: an interrupted job continues where
  it stopped. A job interrupted by stopping the worker or by a restart of the computer goes back
  to the inbox automatically. `retry_errors: true` runs again only the members that ended with an
  error.
- A file with another suffix (e.g. `.yaml.hold`) stays in the inbox without being run.
- When the free memory falls below 8 % of the total (at least 1.5 GB), the worker stops the job
  that uses the most memory and marks it failed with an explanation; its results so far stay.
- Restart the worker after every update of SWINGERsim.

## Reproducibility {#wf-reproducibility}

Every result can be traced back to its input and recomputed:

- the HDF5 file of a run stores the complete configuration (with every default filled in), the
  versions of SWINGERsim, Python, NumPy, SciPy, h5py, REBOUND, REBOUNDx and emcee, and the time
  of its creation; the CSV and fit tables carry the same provenance in their header lines;
- an ensemble folder stores the configuration and a fingerprint of it; `--resume` refuses to mix
  rows of a changed configuration, and the drawn members are fixed by `ensemble.seed`;
- a fit writes the fitted configuration with its epoch (`fit_epoch_bjd_tdb`), the history of
  the optimisation and, with MCMC, the samples (`seed` fixes the walkers);
- the worker keeps every job file, its log and its status in `jobs/done/` or `jobs/failed/`, and
  the web application keeps a frozen copy of every submitted configuration in `jobs/submitted/`.

To share a result, share the configuration, the data tables it refers to and the version of
SWINGERsim; with the seeds in the configuration anyone can recompute it. The *Methods paragraph*
of a result page writes the description of the computation for a paper.
