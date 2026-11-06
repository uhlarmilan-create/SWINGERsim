# Outputs {#outputs}

This chapter describes every file SWINGERsim writes, the meaning of its columns, how to
interpret the values and how to read the files in Python. All tables are CSV files whose first
lines may start with `#`: they record the software versions, the units, the time system and,
for some, the configuration; the first line that does not start with `#` is the header.

## The summary of a run {#out-summary}

`swingersim run` prints a summary, which is also saved as `<name>_summary.txt`. Example 3 of
[Worked examples](#ex-tides) (a 2+2 with tides and rotation) gives, shortened:

```text
integrated to t = 300 yr (1500 outputs, osculating elements)
  out    a = 2.442..2.4431 AU  e = 0.3998..0.4002
  A      a = 0.041764..0.041767 AU  e = 0.0494..0.0502
  i_mut A-B           8.44..  25.08 deg
  max |dE/E| = 2.22e-11
  max |dL/L| = 1.83e-07
hierarchy: configured [[1,2],[3,4]]
  no change detected
orbit pairs:
  A in out: i_mut0 = 11.6 deg, g1_0 = 98 deg
    quadrupole theory : e = 0.048..0.050, i_mut = 11.6..11.6 deg, e-cycle 323 yr (outside ZLK window)
    simulation        : e = 0.049..0.050, i_mut = 11.6..13.6 deg
    precession period : 1.29e+03 yr estimated, not completed in the run, closest periastron = 2.78 x (R1+R2)
    stability ratio Mardling-Aarseth = 7.97 (stable), Vynatheya+22 = 10.17 (stable), eps_oct = 0.00054, eps_SA = 0.0017
    GR apsidal period ~ 2.31e+03 yr
    tidal+rotational apsidal period 64 yr, ZLK time 949 yr -> faster than ZLK, ZLK suppressed
binaries A and B: P_A/P_B = 0.6207 (nearest 3:5, +3.45 %), precession ratio A/B = 1.97 estimated
  mutual inclination of A and B: 8.44..25.08 deg, analytic limit 104.21 deg (Vokrouhlicky 2016, circular A and B)
  eclipses (observer along z): A 100 %, B 91 %, both at once 91 % of the run
secular coupling beta (Hamers & Lai 2017): A 1.45, B 0.365
stars:
  m1  R = 1.7 Rsun, k2 = 0.0144, rg2 = 0.0754, lag = 0 s, P_rot = 1.8 -> 1.8 d, obliquity 0.0 -> 0.1 deg, Roche lobe filled up to 0.51
  binary A (m1+m2): apsidal period 62.29 yr from the binary itself, tides+rotation 64.01 yr, GR 2311 yr; measured 59.43 yr (with the other stars)
  wall time 152.9 s
```

| Line | Content |
| --- | --- |
| `integrated to` | time reached, number of outputs and the kind of elements (osculating for N-body, mean for secular); shorter than `t_max_yr` if an event stopped the run |
| `EVENT` | only when the run stopped early: `collision` (with the stars), `roche_overflow`, `hierarchy_change`, `energy_error`, `wall_time`, `e_reached`, `stellar_evolution_end`, `stellar_birth`, `gw_unresolved`, or `mass_transfer: <kind>` (below) |
| `MERGER`, `MASS TRANSFER STARTED BY THE DYNAMICS` | with `physics.mergers` or `mass_transfer.dynamic_onset`: what happened and how the run continued |
| one line per orbit | minimum..maximum of the semi-major axis [au] and the eccentricity |
| `i_mut X-Y` | minimum..maximum of the mutual inclination of two orbits [deg] |
| max \|dE/E\|, max \|dL/L\| | largest relative error of the energy and of the total angular momentum; a measure of numerical accuracy (with dissipation the energy changes physically and a note says so) |
| `hierarchy` | the configured hierarchy and every confirmed change: time, kind, old → new structure |
| `X in Y` | each inner orbit in its parent orbit: initial mutual inclination and argument of periastron measured from the line of nodes |
| `quadrupole theory` | prediction of the orbit-averaged quadrupole equations: ranges of $e$ and $i_\mathrm{mut}$ and the period of the eccentricity cycle; whether the orbit is inside the ZLK window (39.2–140.8°) |
| `simulation` | the same ranges measured in the run (the maximum of $e$ includes peaks between outputs) and the number of completed ZLK cycles |
| `precession period` | estimated and measured nodal precession period; closest periastron in units of the sum of the stellar radii |
| `stability ratio`, `eps_oct`, `eps_SA` | stability ratios of Mardling & Aarseth (2001) and Vynatheya et al. (2022) for this nested triple (above 1 stable); octupole strength; single-averaging parameter – above about 0.01 the averaged estimates lose accuracy (Luo et al. 2016) and a note says so |
| `GR apsidal period`, `tidal+rotational apsidal period` | relativistic and (with `physics.tides`) tidal and rotational apsidal periods, compared with the ZLK time scale; a note says when ZLK cycles are suppressed |
| `binaries A and B` | 2+2 cores: period ratio and the nearest small-integer commensurability, ratio of the precession periods; the range of the mutual inclination of A and B and its analytic limit (Vokrouhlický 2016) |
| `eclipses` | fraction of the run during which each binary, and both at once, eclipse for an observer along $z$ |
| `secular coupling beta` | parameter $\beta$ of Hamers & Lai (2017) for each inner orbit of a 2+2 |
| `stars` | radius, $k_2$, $I/(MR^2)$, lag, rotation period and obliquity at start and end, largest Roche-lobe filling; with evolving stars the age, radius, $T_\mathrm{eff}$ and phase at start and end |
| `binary X` | each binary of two stars: apsidal period expected from the binary itself and measured in the run (which includes the perturbation by the other stars) |
| `tidal time scales` | with a tidal lag: circularisation, orbital decay and synchronisation times, pseudo-synchronous rotation period |
| `MEGNO` | with `run.megno`: the chaos indicator and the Lyapunov time |
| `wall time` | computing time |

**Table 65.** Lines of the summary of a run. Left: the start of the line; right: its content.

The kinds of `mass_transfer` events are, for a binary, `common_envelope`, `contact`,
`contact_merger_darwin`, `contact_merger_l2`, `common_envelope_merger`, `post_common_envelope`
and `donor_stripped`; for transfer from a tertiary, `triple_common_envelope`,
`inner_roche_overflow`, `dynamical_instability` and `donor_stripped`.

## Files of a run {#out-run}

`swingersim run config.yaml -o out/sys.h5` writes, next to each other:

| File | Content |
| --- | --- |
| `sys.h5` | all data of the run (below) |
| `sys.html` | the HTML report: command line, configuration, summary, interactive plots, the figures |
| `sys_summary.txt` | the summary and the list of the files of the run |
| `sys_timeseries.csv` | one row per output time (at most 20 000 rows; every k-th output for longer runs, all of them in the HDF5 file) |
| `sys_viewer.html` | the interactive 3D viewer (needs plotly) |
| `sys_orbits.png` | $a/a_0$, $e$ and the inclination to the sky of every orbit over time |
| `sys_inclinations.png` | the mutual inclinations over time |
| `sys_conservation.png` | relative errors of the energy and of the angular momentum |
| `sys_stars.png` | mass, radius and $T_\mathrm{eff}$ of the stars over time (only with evolving stars) |

**Table 66.** The standard products of a run. Every run started by the worker or the web
application writes the same files; `--no-viewer`, `--no-plots` and `--no-html` leave out the
viewer, the figures and the report.

The columns of `sys_timeseries.csv` (orbit names as in the configuration; angles in degrees, unlike
the radians of the HDF5 file):

| Column | Content |
| --- | --- |
| `t_yr` | time [yr] from $t = 0$ (negative for backward runs) |
| `a_<orbit>_au`, `e_<orbit>` | semi-major axis and eccentricity of every orbit |
| `i_<orbit>_deg`, `i_inv_<orbit>_deg` | inclination to the sky and to the invariable plane (the plane perpendicular to the total angular momentum) |
| `Omega_<orbit>_deg`, `Omega_inv_<orbit>_deg`, `omega_<orbit>_deg` | node on the sky and on the invariable plane, argument of periastron |
| `imut_<o1>-<o2>_deg` | mutual inclination of every pair of orbits |
| `dE_over_E`, `dL_over_L` | relative errors of the energy and of the angular momentum |
| `M_<star>_msun`, `R_<star>_rsun`, `Teff_<star>_k` | mass, radius and $T_\mathrm{eff}$ (only with evolving stars) |

**Table 67.** Columns of the time-series CSV.

### The HDF5 file {#out-hdf5}

HDF5 is a standard binary format for numerical data, readable from Python (`h5py`) and tools such
as HDFView, MATLAB or Julia.

| Path | Content | Unit |
| --- | --- | --- |
| `time_yr` | output times | yr |
| `energy_error`, `angmom_error` | relative errors of the energy and of the total angular momentum (orbits + spins) | – |
| `orbits/<orbit>/a`, `e` | semi-major axis, eccentricity | au, – |
| `orbits/<orbit>/inc`, `Omega`, `omega`, `M` | inclination to the sky, node, argument of periastron, mean anomaly | **rad** |
| `orbits/<orbit>/inc_inv`, `Omega_inv` | inclination and node with respect to the invariable plane | **rad** |
| `mutual_inclination/<o1>-<o2>` | mutual inclination | **rad** |
| `stars/<m>` (attributes) | values used: `radius_rsun`, `apsidal_k2`, `rg2`, `lag_s`, `mass_msun`, `mass_msun_config`, Roche-lobe filling (`roche_fill_max`, its time), and from the stellar model `luminosity_lsun`, `teff_k`, `age_myr`, `phase`, envelope mass and radius, convective turnover time | R☉, –, –, s, M☉ … |
| `stars/<m>/P_rot_d`, `obliquity_deg` | rotation period and spin–orbit angle at each output (with `physics.tides`) | d, deg |
| `stars/<m>/radius_rsun`, `luminosity_lsun`, `teff_k`, `phase`, `age_myr`, `mass_msun` | the evolving star at each output (with `physics.evolution` / `mass_loss`); phase 0 pre-main sequence, 1 main sequence, 2 Hertzsprung gap, 3 giant branch, 4 core helium burning, 5 early AGB, 6 thermally pulsing AGB | R☉, L☉, K, –, Myr, M☉ |
| `megno` | MEGNO at each output (with `run.megno`) | – |
| `extremes/<orbit>` (attributes) | `e_max`, `t_e_max_yr`, `q_min_au`, `t_q_min_yr`: extremes found between outputs | –, yr, au |
| `hierarchy/detected` | the detected hierarchy at each output, e.g. `[[1,2],[3,4]]` | text |
| `hierarchy` (attributes) | `configured`, `final`, `transitions_json` (confirmed changes: time, kind, from, to), `detected_changes_json` | – |
| `diagnostics/<inner orbit>` (attributes) | the diagnostics of each orbit pair (Table 69) | – |

**Table 68.** Content of the HDF5 file of a run. Angles are stored in radians.

File attributes: `config_yaml` (the complete configuration with all defaults, so that the run can
be repeated), `swingersim_version` and the versions of Python, NumPy, SciPy, h5py, REBOUND,
REBOUNDx and emcee, `created_utc`, `wall_s`, `units_json` (the conventions: Julian years from
$t = 0$, au, M☉, radians, the sky frame), `backend`, `event`, `event_time_yr` and the stars and
age of the event, `dissipative`, `gas_phase`, `mass_loss` (true when the energy is not conserved
physically), `binaries_json` (per binary of two stars: the expected and measured apsidal periods,
the tidal time scales, the eclipse fraction), `system_json` (period ratios and commensurabilities
of 2+2 cores, the analytic limit of the mutual inclination, eclipse fractions, the coupling
parameters, the mass-loss indices, the mass-transfer and common-envelope history of every orbit,
the numbers of passing-star encounters) and `lyapunov_time_yr` (with `run.megno`).

**Hierarchy notation.** Brackets group bound pairs: `[[1,2],[3,4]]` is a 2+2, `[[[1,2],3],4]` a
3+1; unbound parts are separated by `|`, e.g. `[[1,2],3] | 4` is a triple with an escaping star.
The kinds of change are `exchange` (same type, different partners), `reconfiguration` (the type
changes, e.g. 2+2 → 3+1), `ejection` (one star becomes unbound), `breakup` (several unbound
parts) and `capture` (unbound parts become bound).

| Field | Meaning |
| --- | --- |
| `i_mut0_deg`, `g1_0_deg` | initial mutual inclination to the parent orbit; initial argument of periastron from the node on the local invariable plane |
| `in_zlk_window`, `e_max_quadrupole` | whether $i_\mathrm{mut,0}$ lies between 39.2° and 140.8°; the test-particle estimate $\sqrt{1-\frac53\cos^2 i_0}$ |
| `t_zlk_yr`, `t_prec_yr`, `t_gr_yr`, `gr_suppresses_zlk` | ZLK time scale, nodal precession period, relativistic apsidal period, and whether relativity is faster than ZLK |
| `epsilon_oct`, `epsilon_sa` | octupole strength (Naoz 2016) and single-averaging parameter (Luo et al. 2016) |
| `stability_ratio`, `stability_ratio_v22` | Mardling–Aarseth and Vynatheya et al. (2022) ratios |
| `e_min`, `e_max`, `i_mut_min_deg`, `i_mut_max_deg` | ranges in the run |
| `e_max_true`, `t_e_max_yr`, `q_min_au`, `q_min_over_radii` | largest eccentricity including internal samples, its time; smallest periastron, also in units of $R_1+R_2$ |
| `quad_e_min`, `quad_e_max`, `quad_i_min_deg`, `quad_i_max_deg`, `quad_period_yr` | prediction of the averaged quadrupole equations (Naoz et al. 2013) |
| `n_zlk_cycles`, `zlk_period_yr`, `prec_period_measured_yr` | completed eccentricity cycles and their period; precession period measured about the invariable plane |

**Table 69.** Diagnostics of an orbit pair (attributes of `diagnostics/<inner orbit>`).

## Ensembles {#out-ensemble}

`swingersim ensemble config.yaml -o folder` writes into the folder:

| File | Content |
| --- | --- |
| `ensemble.yaml` | the full configuration and its fingerprint (a resumed ensemble must use the same one) |
| `rows.jsonl` | one line per finished member and direction, written as members finish; `--resume` continues from it |
| `summary.csv` | the same rows as a table, sorted by member and direction |
| `summary.parquet` | the same table for pandas (with the `ensemble` group installed) |
| `runs/member_<n>_<direction>.h5` | the full run of every member (only with `save_runs: true`) |
| `report.html` | the report |

**Table 70.** Files of an ensemble.

| Column | Meaning |
| --- | --- |
| `member`, `direction` | number of the member; `forward` or `backward` |
| `masses.m3`, `orbits.out.P`, … | the drawn values of the `vary` parameters |
| `fate` | `intact` (no event, hierarchy unchanged), a hierarchy change (`hierarchy_change`, `exchange`, `ejection`, `breakup`, `reconfiguration`, `capture`), an event (`collision`, `merger`, `roche_overflow`, `energy_error`, `wall_time`, `e_reached`, `gw_unresolved`) or `error` |
| `event`, `event_time_yr`, `t_end_yr` | the event, its time, the end of the run |
| `n_transitions`, `first_transition_yr`, `final_hierarchy` | changes of the hierarchy |
| `<orbit>_e0`, `_e_min`, `_e_max`, `_log1me_min` | initial, smallest, largest eccentricity; $\log_{10}(1-e_\mathrm{max})$ |
| `<orbit>_t_e_max_yr`, `_q_min_au`, `_a_min`, `_a_max` | time of the largest eccentricity, closest periastron, range of $a$ |
| `<orbit>_t_1me_lt_1e-1_yr` … `_1e-4_yr` | first time at which $1-e$ fell below $10^{-1}$ … $10^{-4}$ (empty if never) |
| `<orbit>_U_measured_yr` | apsidal period measured over the run (empty when the periastron advanced too little) |
| `imut_<o1>_<o2>_0`, `_min`, `_max` | initial, smallest, largest mutual inclination [deg] |
| `<binary>_eclipse_fraction`, `eclipse_fraction_both`, `eclipse_fraction_all` | fraction of the run during which the binary (both binaries of the first 2+2; all binaries) eclipse |
| `<binary>_i0`, `_didt0_deg_yr`, `_dtiltdt0_deg_yr` | inclination to the sky now; its rate of change; the rate of change of $\lvert i - 90°\rvert$ (negative: turning towards edge-on, deeper eclipses); fitted one-sided from $t = 0$ over `run.rate_window_yr` |
| `<binary>_eclipsing0`, `_eclipse_change_yr` | 1 if the binary eclipses now; the first time the eclipse state changes |
| `<star>_roche_fill_max`, `roche_first`, `roche_first_kind` | largest Roche-lobe filling; with `roche_overflow`, the star that filled its lobe first and whether its companion is a star or a subsystem |
| `mass_transfer_*` | the scalar mass-transfer and common-envelope entries of the run |
| `ma01_min`, `v22_min` | smallest stability ratio of all nested triples at $t = 0$ |
| `mlp_stable`, `ghost_stable`, `ghost_delta_max` | with `mlp` / `ghost_test`: 1 stable, 0 unstable |
| `max_energy_error`, `wall_s` | quality and cost of the run |

**Table 71.** Columns of the ensemble summary. These names are the ones used in
`ensemble.observed`.

`swingersim ensemble-report` prints, separately forward and backward, the fraction of every fate,
the distribution of the largest eccentricities and the fractions with $1-e$ below 0.1, 0.01 and
0.001, the mutual inclinations, the present rates of change and when eclipses end or began, the
eclipse fractions and the stability verdicts; fractions come with their 68 % Wilson interval
(Wilson 1927). With observations every number is weighted, the header gives the effective number
of members, and the drawn parameters are listed as prior → conditioned median and 16–84 % range.

## Eclipse times, minima, light curves and radial velocities {#out-observables}

**`swingersim timings`** writes one row per eclipse:

| Column | Meaning |
| --- | --- |
| `binary`, `kind` | orbit name; `primary` (the first star of the binary eclipsed) or `secondary` |
| `epoch` | cycle number of the linear ephemeris (half-integer for secondary eclipses) |
| `t_yr`, `t_obs_yr` | time of mid-eclipse [yr from $t = 0$], without and with the light-travel time |
| `o_c_d` | observed minus calculated [d] against the linear ephemeris of the primary eclipses |
| `impact` | projected separation at mid-eclipse / $(R_1+R_2)$: 0 central, near 1 grazing |
| `T0_yr`, `P_d` | the fitted linear ephemeris |
| `bjd_obs` | with `--bjd-t0`: the BJD of the minimum |

**Table 72.** Columns of the eclipse-times CSV. The summary on screen gives the number of
eclipses, the period, the range of the O−C and the intervals without eclipses; with
`physics.activity` it also reports the energy check of every active star.

**`swingersim minima`** writes one row per measured minimum:

| Column | Meaning |
| --- | --- |
| `binary`, `kind`, `epoch` | as for the eclipse times |
| `time`, `t_min_d`, `sigma_d` | time of minimum (BJD_TDB, and days after $t = 0$) from the Mikulášek (2015) profile and its standard error |
| `t_geom_d` | the geometric conjunction with light time [d] |
| `o_c_d`, `o_c_geom_d` | O−C of the minima and of the conjunctions against their own linear ephemerides [d] |
| `depth_mag`, `width_d`, `rms_mag` | depth and width of this eclipse, residuals of its fit |

**Table 73.** Columns of the minima CSV. `--seasons` writes the shape per season (`t_mid_d`, `n`,
`depth_mag`, `width_d`, `gamma`, `corr`, `rms_mag`).

**`swingersim lightcurve`** writes one row per time: `t_day` (days after $t = 0$, observer's
time), `flux` (total flux normalised to the sum of the stars outside eclipse), `flux_m1`, … (the
contribution of every star; pulsating stars include their pulsation, stars behind a disc its
absorption) and `flux_disc_<name>` (the light of discs with emission or scattering). With
`--engine phoebe` the columns are `flux_A`, `flux_B`, … (each inner binary), the single stars and
`flux_mutual` (light lost in eclipses between different binaries). The header lines give the band,
the time convention and the BJD of $t = 0$ if known.

**`swingersim rv`** writes `rv_model.csv` (columns `bjd_tdb`, `rv_m1`, `rv_m2`, …, km/s) and a
figure `rv_model.png` with the measured velocities if given.

**`swingersim apsidal`** writes per pair the period, eccentricity, $k_2$, radii and rotation
periods of both stars, the tidal advance, the rotational advance for synchronous,
pseudo-synchronous and configured rotation, the relativistic advance, the total (deg/yr) and the
apsidal period $U$ (yr).

## Fits {#out-fit}

`swingersim fit` writes into `fit.output.dir`. Stage 1:

| File | Content |
| --- | --- |
| `stage1_parameters.csv`, `.tex` | fitted parameters with 1σ errors (`<pair>.T0`, `.P`, `.Q`, `.e`, `.w0`, `.dw` [rad/cycle], `.asini` [au], `.P3`, `.T3`, `.e3`, `.w3`; `out.*` for a shared outer orbit) and derived values: `.K_s` (semi-amplitude of the light-time effect), `.f_m3` (mass function), `.P3_yr`, `.U_yr` (apsidal period), `.Pa_d` (anomalistic period), `mass_ratio_A_over_B` |
| `stage1_oc.csv` | every minimum: pair, type, epoch, time, sigma, source, O−C against the linear ephemeris, the model, the residual [d] |
| `stage1_summary.json` | $\chi^2$, degrees of freedom, parameters, derived values, error scaling, outliers, ephemerides, software versions |
| `stage1_predictions.csv` | with `output.predict`: future minima and their spread |
| `stage1_oc.png`, `stage1_corner.png` | O−C with model and residuals and the terms separately; corner plot of the MCMC samples |

**Table 74.** Files of a stage-1 fit.

Stage 2:

| File | Content |
| --- | --- |
| `stage2_params.csv` | every free parameter: `parameter`, `start`, `value`, `error` (1σ; with MCMC half the 16–84 % range), `unit`, `prior_mean`, `prior_sigma`, `error_linear` (from the covariance matrix), `p16`, `p50`, `p84` (MCMC percentiles) |
| `stage2_covariance.csv` | the covariance matrix of the linear model |
| `stage2_residuals.csv` | every minimum: `pair`, `type`, `bjd_tdb`, `sigma_d` (the error used), `model_start_bjd`, `model_bjd`, `o_minus_model_min`, `source`, `sigma_quoted_d` (the error of the table), `status` (`fitted` or `outlier`) |
| `stage2_rv.csv` | every velocity: `star`, `bjd_tdb`, `rv`, `sigma`, `sigma_quoted`, `model_start`, `model`, `o_minus_model`, `source`, `status` |
| `stage2_lc.csv` | every light-curve point: `source`, `segment`, `bjd_tdb`, `flux`, `sigma`, `sigma_quoted`, `model_start`, `model`, `o_minus_model` |
| `stage2_inclination.csv` | every inclination epoch: `label`, `pair`, `time`, `i_obs`, `sigma`, `limit`, `i_start`, `i_model`, `i_eclipse_limit`, `z` (normalised residual), `source` |
| `stage2_fitted.yaml` | the configuration with the fitted values and `fit_epoch_bjd_tdb`, without the `fit` section; its header gives $\chi^2$, the rms per pair and the epoch; it runs as it is |
| `stage2_report.txt` | $\chi^2$ from start to fit, rms per pair and data type, parameters with errors, strongly correlated pairs ($\lvert r\rvert > 0.95$), error factors or jitters per source, priors and their deviations, the fits in a row and the outliers; the same text is printed |
| `stage2_history.json` | $\chi^2$ and the parameters at every model evaluation |
| `stage2_samples.csv` | with MCMC: the samples after burn-in |
| `stage2_predictions.csv`, `.png` | predicted minima: `pair`, `type`, `bjd_tdb`, `sigma_min`, `bjd_16`, `bjd_84`, `duration_h`, `impact`, the observing window and, with coordinates, UTC times |
| `stage2_oc.png`, `stage2_rv.png`, `stage2_lc.png`, `stage2_inclination.png` | figures of the O−C, the velocities, the light curves and the inclinations with start and fitted model and residuals |
| `anchored.yaml` | the anchored start configuration (with `anchor` or `--anchor`) |
| `report.html` | the report |

**Table 75.** Files of a stage-2 fit.

## Scattering experiments and chaos maps {#out-scatter}

`swingersim scatter` writes `scatter.yaml` (the configuration and the critical velocity),
`rows.jsonl` (one line per encounter: impact parameter, outcome, status, final hierarchy, end time,
energy error, the bound orbits at the end, stability of a bound subsystem, colliding stars, total
energy and angular momentum, the drawn initial conditions) and `summary.txt` (number, fraction
with its 68 % Wilson interval and cross-section of every outcome).

`swingersim chaosmap` writes `chaosmap.json` (the grid and the configuration), `grid.jsonl` (one
line per grid point: indices, values, `megno`, `lyapunov_yr`, `t_end_yr`, `status` – `ok`,
`escape`, `collision` or `wall_time`), `summary.txt` (numbers of regular, chaotic and escaping
points and the stability limit along x) and `map.png` (MEGNO over the grid, clipped at 10;
escapes and collisions drawn as 10). `lyapunov_yr` (and `lyapunov_time_yr` of a run) is the Lyapunov
time $1/\lambda$: REBOUND's `sim.lyapunov()` is the slope of the mean MEGNO, $\langle Y\rangle \approx \lambda t/2$,
so SWINGERsim reports $1/(2\,\mathrm{slope})$. The initial tangent vector of MEGNO has a fixed seed, so a map
is the same on every run.

## Machine learning {#out-ml}

`emulator-train` writes the model (`.pkl`) and `model.json` next to it (columns, scores,
settings) and, with `--plot`, its diagnostics; `emulator-report` prints the posterior of every
input with the effective sample size. `sbi` prints the posterior (median, 68 % and 95 % intervals)
and optionally writes the network (`--model`), the samples (`--save-samples`), a corner plot and
the calibration diagnostics.

## How to interpret the values {#out-interpret}

- **Osculating elements** of the N-body backend wiggle on orbital time scales (typically below
  1 % in $a$ and 0.01 in $e$); look at trends and at extremes over many outer orbits. The secular
  backend gives mean elements without these wiggles.
- **Semi-major axes** stay constant on average in point-mass dynamics (with mass loss $aM$ stays
  constant). A drift together with a growing energy error points to a numerical problem; a sudden
  jump to a close encounter.
- **Inclination to the sky** decides whether a binary eclipses; **mutual inclination** decides
  the dynamics (precession, ZLK).
- **Ω and ω** can wrap around 360°; unwrap them before plotting.
- **NaN** means that an orbit was not bound at that output, e.g. after a break-up.
- **Apsidal periods** measured in a run include the perturbation by the other stars and can
  differ strongly from the value of the isolated binary.
- **Roche-lobe filling** near 1 means that a star nearly fills its lobe at periastron.
- **Energy error.** Without relativity $10^{-14}$–$10^{-12}$ is normal; with relativity the 1PN
  energy is checked and $10^{-13}$–$10^{-11}$ is typical; values approaching `max_energy_error`
  ($10^{-5}$) mean that the result should not be trusted. With dissipation, gas or mass loss the
  energy changes physically.
- **Primary eclipse** in SWINGERsim is the eclipse of the first star of the binary, whatever the
  depth; an observer usually calls the deeper eclipse primary.
- **Model time** starts at $t = 0$, the epoch of the elements; measured BJDs are compared after
  giving the BJD of $t = 0$ (`--bjd-t0`, `fit.stage2.bjd_t0`).

## Reading the files in Python {#out-python}

```python
import h5py
import numpy as np
import pandas as pd

with h5py.File("output/example.h5") as f:
    t = f["time_yr"][:]
    e_A = f["orbits/A/e"][:]
    i_AB = np.degrees(f["mutual_inclination/A-B"][:])      # radians in the file
    print(f.attrs["event"] or "no event")
    print(f.attrs["config_yaml"])

ts = pd.read_csv("output/example_timeseries.csv", comment="#")
ens = pd.read_csv("output/v994/summary.csv")
fwd = ens[ens.direction == "forward"]
print(fwd.fate.value_counts(normalize=True))
res = pd.read_csv("output/example_fit/stage1_oc.csv", comment="#")
```
