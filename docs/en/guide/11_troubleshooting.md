# Troubleshooting and FAQ {#troubleshooting}

The problems met most often, ordered by where they appear. Each table gives the symptom, its
usual cause and the fix. General questions are answered at the end.

## The configuration is refused {#ts-config}

SWINGERsim checks the whole configuration before it starts and stops with a message that names
the key. In the web application the same message appears next to the field.

| Message (shortened) | Cause | Fix |
| --- | --- | --- |
| `give exactly one of P (days) or a (AU)` | an orbit has both or neither | keep one |
| `hierarchy: unknown child 'x'` | a name in `hierarchy` is neither a star (`m1`, …) nor an orbit | check the spelling; stars are numbered from 1 without gaps |
| `hierarchy: needs exactly one outermost orbit` | two orbits are not joined | add the orbit that joins them |
| `topology … needs orbits [...]`, `masses must be exactly [...]` | the orbits or masses do not match the hierarchy | use *Add stars and orbits* in the web application, or compare with Table 46 |
| `stellar_model.kind = mist needs 'age_myr'` | evolving models need an age | give `stellar_model.age_myr` |
| `physics.rotation needs physics.tides` | rotation and tides share one model of the stellar distortion | add `tides: true` |
| `physics.evolution: stars [...] have a fixed radius R or Teff` | with evolution, radii and temperatures follow from mass and age | remove `R` and `Teff`, or switch evolution off |
| `physics.mass_transfer needs physics.evolution: true and stellar_model.kind sse or mist` | transfer histories come from the tracks | switch on evolution with `sse` or `mist` |
| `physics.gas needs run.backend: secular` | the gas phase is implemented in the secular backend | set `run.backend: secular` |
| `run.megno assumes constant masses` | MEGNO with mass loss or transfer | switch one of them off |
| `give value and sigma, or min and/or max (or equals)` | an `observed` entry is incomplete | an observation is a value with an error, a limit, or a text value |
| `fit.stage2.priors.<name>: not a free parameter` | a prior on a fixed parameter has no effect | free it or remove the prior |
| `fit.lc red_noise: true … needs fit.stage2.error_scale … or clip_sigma` | red-noise factors are applied between fits in a row | add `error_scale` or `clip_sigma` |
| `Extra inputs are not permitted` | a typo in a key; keys are case-sensitive | compare with [The configuration file](#configuration) |

**Table 76.** Common messages when a configuration is refused.

## The run is slow or stops early {#ts-run}

| Symptom | Cause | Fix |
| --- | --- | --- |
| an N-body run takes hours | every revolution of the shortest orbit needs tens of steps | use the secular backend if the hierarchy stays intact; shorten `t_max_yr`; install Numba (`fast`) for the secular backend |
| `EVENT: collision` early | two stars touch at periastron | expected for high eccentricities; check the radii and whether tides and relativity suppress ZLK cycles |
| `EVENT: roche_overflow` | a star fills its Roche lobe | switch on `physics.mass_transfer` to follow it, or take the time of the event as the answer |
| `EVENT: stellar_evolution_end` | a star reached the end of its track | the run cannot continue past the models |
| `EVENT: wall_time` | `run.max_wall_s` reached | raise it or use the secular backend |
| `EVENT: energy_error` | $\lvert\Delta E/E\rvert$ exceeded `run.max_energy_error`, usually at a very close approach | look for near-collisions; check with the stability tools |
| the secular backend refuses a run with mass loss | the mass-loss index Ψ exceeds `max_psi`: the orbit-averaged adiabatic response is not valid | use the N-body backend |

**Table 77.** Problems during a run.

## The results look wrong {#ts-results}

| Symptom | Cause | Fix |
| --- | --- | --- |
| no eclipses although the binary is observed to eclipse | the default zero-age main-sequence radii are too small for the inclination | give the measured `stars.<m>.R` or an evolved stellar model; the Quick look shows the eclipse limit |
| apsidal motion much too fast | the polytropic $k_2$ of `zams` is too large for stars of 2–3 M☉ | use `mist` or `sse`, or give measured `apsidal_k2` |
| ZLK cycles reach $e \approx 1$ and the stars collide | point masses: no tidal or relativistic precession | switch on `gr`, `tides` and `rotation` |
| osculating elements wiggle | they are instantaneous ellipses, perturbed on every orbit | normal; the secular backend reports the averages |
| a rate of change differs from a published one | different time window or definition | set `run.rate_window_yr` to about half of the observing baseline |
| primary and secondary eclipses swapped against a paper | SWINGERsim's primary is the eclipse of the first star of the pair | swap the stars or the types |
| the O−C of a model has the wrong phase against measured minima | the model has no calendar | give the BJD of $t = 0$ (`--bjd-t0`) |
| different ensemble result each time | the seed changed | fix `ensemble.seed` |

**Table 78.** Results that look wrong.

## Ensembles {#ts-ensembles}

| Symptom | Cause | Fix |
| --- | --- | --- |
| an interrupted ensemble starts from zero | `--resume` missing | run the same command with `--resume` |
| members ended with `error` | a configuration that fails for some draws | read the message in `rows.jsonl`; after a fix use `--resume --retry-errors` |
| the conditioned result rests on a few members | precise observations, broad prior | run more members, narrow the prior, or use an emulator or SBI |
| fractions with wide intervals | too few members | more members narrow the 68 % interval as $1/\sqrt{n}$ |
| `require` rejects most draws | the prior produces mostly unstable systems | narrow the drawn ranges |
| `observed columns not in the summary` | a typo in a column name, or a column the run does not produce | compare with Table 71 |

**Table 79.** Problems with ensembles.

## Fits of measured data {#ts-fits}

| Symptom | Cause | Fix |
| --- | --- | --- |
| stage 2 jumps to the wrong cycle: residuals of whole periods, $\chi^2$ above the start | the start is too far from the solution, so the minima are numbered differently | anchor first (`anchor: true` or `--anchor`); free only `P` and `M0` in a first fit and the rest in a second one started from its `stage2_fitted.yaml`; narrow `bounds` |
| warning "strongly correlated parameters" | the data do not separate them (e.g. $M_0$, $\omega$ and $P$ of an outer orbit covered less than twice) | the large, correlated errors are correct; check them with `sampler: emcee`, or fix one parameter |
| stage 2 is slow | each evaluation integrates the span of the data; each iteration needs two more per free parameter | raise `processes`; leave out imprecise minima (`max_sigma_d`); fit fewer parameters |
| a few minima have residuals of many σ | wrong type or pair in the table, or outliers | check those rows in `stage2_residuals.csv`; `clip_sigma: 5`, a robust `loss`, or `exclude` |
| reduced $\chi^2$ of tens to thousands | quoted errors far below the real scatter | `error_scale: per_source` with `error_model: jitter` |
| warning "the errors / outliers still changed after N fits" | the corrections did not settle within `refits` | raise `refits`, or start the next fit from `stage2_fitted.yaml` |
| primary and secondary residuals differ by several σ | a wrong eccentricity vector of the pair (a local minimum) | restart with $e\cos\omega$ of the opposite sign, or free `ecosw` and `esinw` |
| a parameter on its bound | the solution is probably not the true minimum, its error is meaningless | widen `bounds` or add a prior; check whether the data constrain it |
| the fitted outer orbit disagrees with a published one | the minima constrain it weakly | add priors from the published solution; the report gives each prior's deviation in σ |
| a stage-1 apsidal rate far from the expected value, high reduced $\chi^2$ | minima sampled at regular cycle intervals allow aliases of the rate | add minima at other cycles, or start or constrain the rate near the expected value ([Worked examples](#ex-fit1)) |
| `sampler: emcee needs the package emcee` | the `fit` group is not installed | `pip install -e ".[fit]"` |

**Table 80.** Problems with fits.

## The worker and the web application {#ts-worker}

| Symptom | Cause | Fix |
| --- | --- | --- |
| jobs stay in the inbox | the worker is not running, or the jobs are on hold (`.yaml.hold`) | start the worker (Jobs page or `swingersim worker`); queue the jobs |
| a job runs with old code | the worker was started before an update | stop it (Jobs page or an empty file `jobs/stop`) and start it again |
| a job failed with an out-of-memory message | the worker stopped the job that used the most memory | fewer parallel `processes`, fewer outputs, a shorter run |
| a file is ignored in the inbox | its name does not end in `.yaml` | rename it |
| the plots of the web application do not appear | plotly is not installed and there is no internet connection | `pip install -e ".[viz]"` |
| the page cannot be reached from another computer | the server listens on 127.0.0.1 | `--host 0.0.0.0` (only on a trusted network) |

**Table 81.** Problems with the worker and the web application.

## Frequently asked questions {#faq}

**How many stars can a system have?** Any hierarchy of nested pairs; systems of up to ten stars
are tested. The cost grows with the number of short inner orbits, not with the number of stars.

**Can SWINGERsim fit my light curve?** It fits light-curve segments together with the times of
minima in stage 2 of `swingersim fit`, with spherical, rotating or Roche-shaped stars and a
calibration per segment; it is a photodynamical fit of the whole system, not a replacement for a
detailed model of one binary's surface (spots, detailed atmospheres).

**Can it run backwards in time?** Yes, with a negative `run.t_max_yr` or
`ensemble.directions: [backward]`. Dissipative processes (tides, winds, mass transfer) are not
time-reversible, and chaotic systems make sense backwards only as ensembles (*Physics and
Methods*, chapter Limitations).

**Which backend should I use?** N-body for anything an observer measures and for up to thousands
of years; secular for millions of years and large ensembles ([Two backends](#backends)).

**Which stellar model should I use?** `zams` for quick dynamical studies of unevolved stars;
`mist` for real systems and anything that ages; `sse` for fast evolution with mass transfer up to
helium ignition. Measured radii and temperatures in `stars` always take precedence.

**Why do the eclipse times "exact" and "as an observer" differ?** The first are the conjunctions
of the dynamics; the second are minima measured on a synthetic light curve. Asymmetric eclipses,
the light of the other stars, pulsations, noise and sampling move a measured minimum away from
the conjunction. Compare measured minima with the second.

**How accurate is SWINGERsim?** Every physical ingredient is compared with published results in
*Testing, Part I*; applications to real systems are in Parts II–IV.

**How do I cite the results?** Cite the papers of the methods used (the *Methods paragraph* of a
result page lists them) and give the configuration and the version of SWINGERsim, which every
output file records.
