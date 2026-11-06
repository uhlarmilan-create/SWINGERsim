# The configuration file {#configuration}

A configuration is a YAML text file that describes one system and what to do with it. This
chapter is its complete reference: every section, every key with its type, default, allowed
values and meaning, as checked by the program. A key that is not listed here is refused, and so
is a value outside its range; the error message names the key.

## Units and conventions {#cfg-units}

| Quantity | Unit in the file |
| --- | --- |
| masses | M☉ |
| periods of orbits and stars | days |
| semi-major axes, disc radii | au |
| stellar radii | R☉ |
| angles (inclination, node, argument of pericentre, mean anomaly, obliquity) | degrees |
| times of runs, time scales | years; ages in Myr |
| times of measurements | BJD_TDB (unless `fit.time_system` says otherwise) |
| velocities | km/s |

**Table 44.** Units of the configuration file. Internally SWINGERsim works in M☉, au and years.

**Reference frame.** Inclinations and nodes are measured from the plane of the sky; the $z$ axis
is along the line of sight, pointing away from the observer, so $i = 90°$ is edge-on (eclipses
possible) and $i = 0°$ face-on. Values of $i$ above 90° describe the same orbit revolving in the
opposite sense; this matters for mutual inclinations, not for eclipses. Each orbit is the
relative orbit of its second component around its first, both taken as their centres of mass
(Jacobi coordinates). The elements given are osculating elements at $t = 0$, the epoch of the
configuration.

**Names.** Stars are `m1`, `m2`, … without gaps. Orbits have the names of the hierarchy: `A`,
`B`, `out` for 2+2; `inner`, `mid`, `out` for 3+1; your names for `custom`.

## Top level {#cfg-top}

| Key | Type, default | Meaning |
| --- | --- | --- |
| `topology` | text, required | `2+2`, `3+1`, `2+2+1`, `3+1+1`, `2+(2+1)`, `2+2+2`, `3+3`, or `custom` |
| `hierarchy` | map orbit → [child 1, child 2]; only with `custom` | the tree: each orbit has two children, a star or another orbit; exactly one orbit is nobody's child |
| `masses` | map star → M☉, required | one positive mass per star |
| `orbits` | map orbit → orbit, required | one entry per orbit of the hierarchy (below) |
| `stars` | map star → star, default empty | properties beyond the mass (below) |
| `discs` | map name → disc, default empty | dust discs and rings that occult stars in light curves (below) |
| `stellar_model` | section | where radii, temperatures and structure come from (below) |
| `physics` | section | the physics beyond Newtonian point masses (below) |
| `run` | section, required (except in a configuration with only a `fit`) | how the motion is computed (below) |
| `ensemble` | section, optional | a Monte Carlo ensemble (`swingersim ensemble`); ignored by `run` |
| `fit` | section, optional | a fit of measured data (`swingersim fit`); ignored by `run` |
| `fit_epoch_bjd_tdb` | number, optional | written by a stage-2 fit into `stage2_fitted.yaml`: the BJD_TDB of $t = 0$ of the fitted elements; ignored by runs |

**Table 45.** Keys at the top level of a configuration.

The presets of `topology` are shortcuts for these trees (Figure 1):

| Preset | Tree |
| --- | --- |
| `2+2` | `out: [A, B]`, `A: [m1, m2]`, `B: [m3, m4]` |
| `3+1` | `out: [mid, m4]`, `mid: [inner, m3]`, `inner: [m1, m2]` |
| `2+2+1` | `out5: [out, m5]`, `out: [A, B]`, `A: [m1, m2]`, `B: [m3, m4]` |
| `3+1+1` | `out5: [out, m5]`, `out: [mid, m4]`, `mid: [inner, m3]`, `inner: [m1, m2]` |
| `2+(2+1)` | `out: [A, Bt]`, `A: [m1, m2]`, `Bt: [B, m5]`, `B: [m3, m4]` |
| `2+2+2` | `out: [AB, C]`, `AB: [A, B]`, `A: [m1, m2]`, `B: [m3, m4]`, `C: [m5, m6]` |
| `3+3` | `out: [T1, T2]`, `T1: [A, m3]`, `A: [m1, m2]`, `T2: [B, m6]`, `B: [m4, m5]` |

**Table 46.** The preset hierarchies. Each orbit is written as `name: [first child, second
child]`. A triple is written as `custom`, e.g. `hierarchy: {out: [inner, m3], inner: [m1, m2]}`.
Trees of up to ten stars are tested.

## `orbits.<name>` {#cfg-orbits}

| Key | Type, default | Meaning |
| --- | --- | --- |
| `P` | number > 0 | period [d]; give exactly one of `P` and `a` |
| `a` | number > 0 | semi-major axis [au] |
| `e` | 0 ≤ e < 1, default 0 | eccentricity |
| `i` | 0–180°, default 0 | inclination to the plane of the sky |
| `Omega` | degrees, default 0 | longitude of the ascending node |
| `omega` | degrees, default 0 | argument of pericentre |
| `M0` | degrees, default 0 | mean anomaly at $t = 0$ |
| `Pdot` | d/yr, default 0 | imposed secular period change of an orbit of two stars, $P(t) = P + \dot P\,t$ (mass transfer or angular-momentum loss of semi-detached and contact binaries); N-body backend; free in a fit as `orbits.<name>.Pdot` |
| `Pdot_mode` | `angular_momentum` (default) or `conservative` | what changes the period: angular-momentum loss at constant masses, or conservative transfer between the two stars ($\dot M_1 = (\dot P/P)\,M_1 M_2 / (3(M_1 - M_2))$) |
| `contact_fillout` | 0–1, optional | light curves: the two stars of this binary are a contact binary in a common envelope with fill-out factor $f$ (Lucy 1968); their radii follow from $f$ and they share one temperature |

**Table 47.** Keys of an orbit.

## `stars.<name>` {#cfg-stars}

Everything left out comes from the stellar model.

| Key | Type, default | Meaning |
| --- | --- | --- |
| `R` | R☉ > 0 | radius |
| `Teff` | K > 0 | effective temperature (light curves; with `R` it fixes the luminosity) |
| `apsidal_k2` | 0 < k2 < 0.75 | apsidal-motion constant (Love number = 2 $k_2$) |
| `rg2` | 0 < rg2 ≤ 0.4 | moment of inertia $I/(MR^2)$ |
| `lag_s` | s ≥ 0 | constant tidal time lag (Mignard 1979); 0 = no dissipation |
| `P_rot` | days, `sync` (default) or `pseudosync` | rotation period; `sync`: the mean motion of its binary; `pseudosync`: pseudo-synchronous rotation of an eccentric binary (Hut 1981) |
| `obliquity` | 0–180°, default 0 | angle between the spin and the orbit normal |
| `spin_azimuth` | degrees, default 0 | azimuth of the spin around the orbit normal, from the node |
| `age_offset_myr` | Myr, default 0 | `sse`, `mist`: apparent age of this star minus `stellar_model.age_myr` (products of binary interaction) |
| `pulsation` | section, optional | pulsations in light curves (below) |
| `shape` | `sphere` (default), `rotating`, `roche` | light curves: a sphere; an oblate, gravity-darkened rigid rotator; the Roche equipotential in its circular synchronised binary |
| `radius_inflation` | −0.5 to 1, default 0 | fractional radius change against the stellar model at fixed luminosity (e.g. 0.05 for magnetically active M dwarfs) |

**Table 48.** Keys of a star.

### `stars.<name>.pulsation` {#cfg-pulsation}

| Key | Type, default | Meaning |
| --- | --- | --- |
| `type` | `manual` (default), `lpv`, `delta_scuti`, `gamma_dor`, `beta_cep`, `spb`, `rr_lyrae`, `cepheid` | kind of pulsator; the period of `lpv`, `delta_scuti`, `beta_cep`, `cepheid` and `rr_lyrae` follows from published relations, the others need `period_d` or `Q` |
| `mode` | `F` (default), `1O`, `2O`, `3O` | radial mode |
| `period_d` | days > 0 | period; overrides the relation |
| `Q` | days > 0 | pulsation constant, $P = Q/\sqrt{\rho/\rho_\odot}$ |
| `amplitude_mag` | mag ≥ 0, default 0.01 | semi-amplitude in the band |
| `phase` | cycles, default 0 | phase at emission time $t = 0$ |
| `period_spacing_s`, `n_modes` | s; 1–100, default 1 | gravity modes (`spb`, `gamma_dor`): a series of `n_modes` modes $P_0 + k\,\Delta P$ |
| `extra_modes` | list | further modes, each with `period_d` or `frequency_ratio`, `amplitude_mag` (required) and `phase` |

**Table 49.** Keys of a pulsation.

## `discs.<name>` {#cfg-discs}

A dust disc or ring around a star or a binary that occults stars behind it: a slab of outer
radius `r_out_au`, optional inner hole and full thickness, dimming starlight by
$\exp(-\tau L/\mathrm{thickness})$ along each line of sight.

| Key | Type, default | Meaning |
| --- | --- | --- |
| `host` | star or orbit, required | the disc is centred on the star or on the centre of mass of the stars of that orbit |
| `r_out_au`, `r_in_au` | au; `r_in_au` default 0 | outer radius (required) and inner hole |
| `thickness_au` | au, required | full thickness |
| `tau` | > 0, default none (opaque) | optical depth through the full thickness, perpendicular to the disc |
| `i`, `Omega` | degrees, defaults 90 and 0 | orientation of the disc plane, as for orbits |
| `relative_to`, `tilt`, `tilt_node` | orbit; degrees | tilt the disc from the plane of an orbit (then `i`, `Omega` are ignored) |
| `precession_period_yr` | years | with `relative_to`: precession of the disc normal about the orbit normal (negative = retrograde) |
| `dynamics` | `prescribed` (default) or `rigid` | the precession as given, or computed for a rigid disc torqued by the binary or by the companion of its host (needs `relative_to`) |
| `sigma_index`, `damping_yr` | 0–2.5, default 1; years | `rigid`: surface-density exponent; alignment time scale by viscous dissipation |
| `edge_au` | au, optional | soft outer edge (needs a finite `tau`) |
| `emission` | default false | add the thermal light of a flat reprocessing disc |
| `albedo`, `g` | 0–1, default 0; −1 to 1, default 0.5 | scattered starlight: single-scattering albedo and asymmetry parameter |

**Table 50.** Keys of a disc. Discs affect light curves and minima only, not the dynamics.

## `stellar_model` {#cfg-stellar-model}

| Key | Type, default | Meaning |
| --- | --- | --- |
| `kind` | `zams` (default), `sse`, `mist`, `table`, `python` | `zams`: zero-age main sequence (Tout et al. 1996) and polytropes; `sse`: stars of age `age_myr` after Hurley et al. (2000), up to helium ignition, 0.1–8 M☉; `mist`: MIST tracks (Choi et al. 2016) from the pre-main sequence to the end of the AGB; `table`: own CSV table; `python`: own function |
| `age_myr` | Myr ≥ 0 | `sse`, `mist` (required there): age of the stars at $t = 0$; `sse` counts from the ZAMS, `mist` from the start of the pre-main sequence |
| `feh` | default 0 | `mist`: metallicity [Fe/H] of an installed grid (−1.50 to +0.50 in steps of 0.25) |
| `vvcrit` | 0 or 0.4 (default) | `mist`: initial rotation of the models |
| `z` | 0.004–0.03, default 0.02 | `sse`: metallicity $Z$ |
| `structure` | `auto` (default), `polytrope`, `claret2023`, `mesa:<grid>` | source of $k_2$ and $I/(MR^2)$: the Claret (2023) models at the star's mass and $\log g$, polytropes of index 1.5 or 3, or an own MESA grid; `auto` = `claret2023` for `sse`, polytropes for `zams` |
| `file` | path | `table`: the CSV file (relative to the configuration) |
| `entry` | `file.py:function` or `module:function` | `python`: the function |

**Table 51.** Keys of the stellar model. `mist` always uses the Claret structure (Claret 2004
above 8 M☉).

## `physics` {#cfg-physics}

| Key | Type, default | Meaning |
| --- | --- | --- |
| `gr` | default true | 1PN relativity (N-body), periastron advance (secular) |
| `gw` | default false | gravitational-wave emission of orbits of two stars (Peters 1964) |
| `tides` | default false | tidal distortion of the stars (apsidal motion; dissipation with a lag) |
| `rotation` | default false | rotational flattening and spin–orbit coupling (needs `tides`) |
| `evolution` | default false | the stars follow their tracks during the run (needs `stellar_model.kind` `sse` or `mist`; no fixed `R` or `Teff`) |
| `tidal_lag` | `fixed` (default) or `hurley2002` | the lag from `stars.<m>.lag_s`, or from the envelope of every star (Hurley et al. 2002; needs `tides`) |
| `tidal_efficiency` | > 0, default 1 | `hurley2002`: factor on the dissipation rate |
| `tidal_convective_reduction` | `quadratic` (default) or `linear` | `hurley2002`: reduction of the turbulent viscosity for tidal periods shorter than the convective turnover time |
| `k2_rotation` | default false | correct $k_2$ for the rotation of the stars (Claret 1999) |
| `inertial_waves` | default false | add the dissipation of inertial waves in convective envelopes (needs `hurley2002`) |
| `gas` | section | parametric gas phase (secular backend): `t_end_yr`, `t_decay_yr`, `plane` (`invariable` or `{i, Omega}`), `orbits: {<orbit>: {tau_a_yr, tau_e_yr, tau_i_yr}}` |
| `mass_loss` | section | stellar winds in the dynamics (below) |
| `mass_transfer` | section | Roche-lobe overflow in two-star binaries (below) |
| `mergers` | section | continue after two stars merge (below) |
| `galactic_tide` | section | tide of the Galactic disc: `rho0_msun_pc3` (0.185), `pole_i_deg` (0), `pole_Omega_deg` (0) |
| `encounters` | section | passing stars: `n_star_pc3` (0.1), `sigma_kms` (40), `m_perturber` (0.5), `b_max_au` (default 30 × the outermost $a$), `seed` (1) |
| `activity` | map star → cycle | magnetic activity cycle (Applegate 1992) in the eclipse times: `p_mod_yr`, `amplitude_s` (required), `phase_deg` (0) |
| `magnetic_braking` | section | wind braking of spins (Matt et al. 2015; needs `evolution` and `rotation`): `chi` (10), `p` (2), `t0_erg` ($9.5\times10^{30}$), `stars` (default all), `max_step_yr` ($10^6$) |

**Table 52.** Keys of `physics`. Default values of sub-keys in brackets.

### `physics.mass_loss` {#cfg-mass-loss}

| Key | Type, default | Meaning |
| --- | --- | --- |
| `source` | `track` (default) or `rates` | masses from the tracks (needs `evolution`), or constant rates |
| `rates` | map star → M☉/yr | `rates`: $dM/dt$ (negative = loss) |
| `reimers_eta` | default 0.5 | SSE stars: coefficient of the Reimers wind |
| `masses` | `initial` (default) or `current` | `track`: the configured masses are initial masses, or present masses at `age_myr` |
| `branch` | `lowest` (default) or `unique` | `current`: which initial mass to take when several tracks reach the present mass |
| `max_psi` | default 0.05 | secular backend: refuse runs whose mass-loss index $\Psi = \lvert\dot M\rvert/(nM)$ exceeds this |

**Table 53.** Keys of `physics.mass_loss`.

### `physics.mass_transfer` {#cfg-mass-transfer}

Needs `physics.evolution: true` with `stellar_model.kind` `sse` or `mist`.

| Key | Type, default | Meaning |
| --- | --- | --- |
| `rejuvenation` | `hurley` (default) or `fractional_age` | rejuvenation of accreting main-sequence stars |
| `stability` | `bse` (default) or `zeta` | critical mass ratios of Hurley et al. (2002), or the adiabatic exponents of giant donors |
| `q_crit_ms` | default 3.0 | critical mass ratio of main-sequence and core-helium-burning donors |
| `accretion` | `bse` (default) or `fractions` | accretion limited by the thermal time scale of the accretor, or fixed fractions |
| `alpha`, `beta` | 0–1, defaults 0 | `fractions`: mass lost from the donor's and the accretor's vicinity ($\alpha + \beta \leq 1$) |
| `dynamic_onset` | `angular_momentum`, `periastron` or none (default) | continue when the dynamics drives a star over its Roche lobe: circularise the binary and start the transfer |
| `tertiary` | section (default present) or null | transfer from a star onto the binary it orbits (secular backend): `accretion` (`auto`), `split` (`auto`), `inner_orbit` (`energy` or `angular_momentum`), `alpha_lambda` (5), `h` (0) |
| `contact` | section or none (default) | follow contact binaries: `braking` (true), `q_min` (0.09) |
| `common_envelope` | section or none (default) | follow a common envelope with the energy formalism: `alpha` (1), `lam` (0.5), `duration_yr` (1) |

**Table 54.** Keys of `physics.mass_transfer`. The prescriptions are in *Physics and Methods*,
chapter Mass transfer, common envelope and mergers.

### `physics.mergers` {#cfg-mergers}

| Key | Type, default | Meaning |
| --- | --- | --- |
| `mass_loss` | `glebbeek2013` (default) or a fraction | mass ejected by the merger |
| `rejuvenation` | `glebbeek` (default), `hurley` or a number | fractional main-sequence age of the product relative to the mass-weighted ages of the two stars |
| `spin_fraction` | 0–1, default 0.3 | rotation of the product as a fraction of break-up |

**Table 55.** Keys of `physics.mergers`; `mergers: {}` uses the defaults.

## `run` {#cfg-run}

| Key | Type, default | Meaning |
| --- | --- | --- |
| `t_max_yr` | years, required, ≠ 0 | end of the run; negative = into the past |
| `backend` | `nbody` (default) or `secular` | the backend ([How SWINGERsim works](#backends)) |
| `integrator` | `ias15` | the N-body integrator |
| `n_output` | ≥ 2, default 1000 | number of stored outputs |
| `collisions` | default true | stop when two stars touch (N-body) |
| `stop_on_hierarchy_change` | default false | stop at the first confirmed change of the hierarchy |
| `hierarchy_persistence` | ≥ 1, default 3 | outputs a new hierarchy must persist before it is confirmed |
| `stop_on_roche_overflow` | default false | stop when a star fills its Roche lobe |
| `stop_when_e_reached` | map orbit → e | N-body: stop once every listed orbit has reached its eccentricity |
| `max_wall_s` | seconds | stop after this computing time (event `wall_time`) |
| `max_energy_error` | default $10^{-5}$ | abort if $\lvert\Delta E/E\rvert$ exceeds this |
| `extremes_every_orbits` | ≥ 0, default 10 | internal sampling of $e_\mathrm{max}$ and the closest periastron, in periods of the shortest orbit; 0 = only at outputs |
| `rate_window_yr` | years | length of the linear fits for the present rates of change of the inclinations, on each side of $t = 0$ (default: two periods of the perturbing orbit, at least 6 outputs) |
| `megno` | default false | compute MEGNO (Newtonian; not with mass loss or transfer) |
| `secular_octupole` | default true | secular backend: include octupole terms |
| `secular_rtol` | default $10^{-10}$ | secular backend: relative tolerance |
| `evolution_max_dr` | 0–0.2, default 0.01 | `evolution`: largest relative change of a radius within one segment |

**Table 56.** Keys of `run`.

## `ensemble` {#cfg-ensemble}

| Key | Type, default | Meaning |
| --- | --- | --- |
| `n` | ≥ 1, required | number of members |
| `seed` | default 1 | random seed of the draws |
| `vary` | map parameter → value or distribution | the drawn parameters: `masses.<m>`, `orbits.<o>.<key>`, `stars.<m>.<key>` (numeric keys only) |
| `directions` | list of `forward`, `backward`; default `[forward]` | integrate each member to $+\lvert t_\mathrm{max}\rvert$, $-\lvert t_\mathrm{max}\rvert$ or both |
| `require` | list of `stable_ma01`, `stable_v22`, `hierarchical` | redraw members that a criterion classifies as unstable (`hierarchical`: every inner apocentre inside the pericentre of the orbit around it) |
| `max_tries` | default 1000 | draws allowed per member |
| `ma01_unstable_keep`, `ma01_window` | 0–1, default 1; [low, high] | keep only a fraction of the draws classified unstable, or only draws whose smallest MA01 ratio lies in a window (training sets for classifiers) |
| `mlp` | default false | add the verdict of the classifiers of Vynatheya et al. (2023) to every member |
| `ghost_test`, `ghost_orbits`, `ghost_max_wall_s` | default false; 100; none | add the N-body ghost test to every member |
| `save_runs` | default false | also store the HDF5 file of every member |
| `observed` | map column → observation | measurements that weight the members (below) |

**Table 57.** Keys of `ensemble`.

A value of `vary` is a number (fixed), `isotropic` (an inclination with $\cos i$ uniform),
`angle` (uniform 0–360°) or a distribution with exactly one of the kinds of Table 58 and the
optional keys `min`, `max` (truncation: the value is drawn again) and `mirror` (an inclination is
replaced by $180° - i$ with probability ½, for eclipsing binaries of unknown sense of
revolution).

| Kind | Form | Meaning |
| --- | --- | --- |
| `uniform` | `[low, high]` | uniform |
| `loguniform` | `[low, high]`, positive | uniform in the logarithm |
| `normal` | `[mean, sigma]` | Gaussian |
| `isotropic` | `[i1, i2]` | inclination with $\cos i$ uniform between two limits |
| `choice` | list | one of the listed values |

**Table 58.** Distributions of `ensemble.vary`.

An entry of `ensemble.observed` names a column of the ensemble summary ([Outputs](#out-ensemble))
and is one of:

- a measurement `{value, sigma}` – Gaussian weight;
- limits `{min, max}` (either or both) – weight 0 outside; `missing_ok: true` lets a missing
  value pass;
- a text value `{equals: intact}` or a list of values – for text columns such as `fate`.

The column name may start with `forward:`, `backward:` or `central:` (the mean of both) to take
the value from that run; the default is `central` for the present rates of change
(`*_didt0_deg_yr`, `*_dtiltdt0_deg_yr`) and `forward` for everything else.

## `fit` {#cfg-fit}

| Key | Type, default | Meaning |
| --- | --- | --- |
| `minima` | path, required | CSV table of minima: `pair`, `type` (`p`/`s`), `time`, `sigma`, `source`, `method` (also `binary`, `kind`, `bjd_tdb`, `sigma_d`, `ref`) |
| `time_system` | `BJD_TDB` (default), `HJD_UTC`, `JD_UTC` | time system of the table; the others are converted (needs `coordinates`; times before 1960 are taken as Universal Time, TT = UT + ΔT of Espenak & Meeus 2006) |
| `coordinates` | `{ra_deg, dec_deg}` | ICRS coordinates of the system |
| `ephemerides` | map pair → `{T0, P}` | linear ephemeris for the epochs and the O−C (default: fitted to the primary minima) |
| `exclude` | `{before, after, sources, pairs}` | drop minima by time, source or pair |
| `error_scale` | `none` (default), `global`, `per_source` | rescale the errors so that the reduced $\chi^2$ is 1 |
| `rv` | section | radial velocities (Table 60) |
| `lc` | list | light-curve segments (Table 61) |
| `inclination` | list | inclinations at epochs: `file` (columns `mid_hjd` or `time`, `i_deg`, `sigma_deg`, optional `limit` and `source`), `pair`, `label`, `min_sigma_deg` (0) |
| `stage1` | section | the analytic O−C model (Table 62) |
| `stage2` | section | the joint N-body model (Table 63) |
| `output` | section | `dir` (`output/fit`), `tables` (`[csv]`, also `latex`), `plots` (`[oc, oc_components]`, also `corner`), `predict` (`years` 5, `pairs`, `start`, `samples` 100) |

**Table 59.** Keys of `fit`. Paths are relative to the configuration.

| Key of `fit.rv` | Type, default | Meaning |
| --- | --- | --- |
| `file` | path, required | CSV: `star`, `time` (BJD_TDB), `rv`, `sigma` [km/s], `source` |
| `gamma` | km/s, default 0 | systemic velocity (start value; free as `rv.gamma`) |
| `offsets` | map source → km/s | zero-point offsets (start values; free as `rv.offset.<source>`) |
| `stars` | list | stars used (default all) |
| `max_sigma` | km/s | leave out velocities with a larger error |
| `centre` | `mass` (default) or `light` | velocity of the centre of mass, or flux-weighted velocity of the visible surface (distortion, rotation, Rossiter–McLaughlin effect) |
| `band` | default `V` | passband of the flux weights for `centre: light` |

**Table 60.** Keys of `fit.rv`.

| Key of `fit.lc[]` | Type, default | Meaning |
| --- | --- | --- |
| `file` | path, required | CSV: `time` (BJD_TDB), `flux` (normalised), `sigma`, optional `segment` |
| `label` | text | name in the report |
| `band` | default `TESS` | passband |
| `exposure_min`, `n_sub` | min ≥ 0, default 0; 1–15, default 1 | exposure or binning time and model points per exposure |
| `baseline_degree` | 0–2, default 1 | polynomial trend per segment |
| `scale` | default true | scale the model per segment (dilution); false: depths from the stars as given |
| `gap_d` | default 0.5 | without a `segment` column: split segments at gaps longer than this [d] |
| `segments` | list | use only these segments |
| `shapes` | `auto` (default), `as_configured`, `spheres` | star shapes: Roche surfaces where a star fills its lobe beyond `roche_threshold`, as in the configuration, or spheres |
| `roche_threshold` | 0–1, default 0.5 | volume radius / Roche-lobe radius above which `auto` uses a Roche surface |
| `mesh_resolution`, `coarse_resolution` | default (64, 128) and (24, 48) | mesh directions of the accurate and the coarse model |
| `fidelity` | `mapped` (default) or `direct` | fit on the coarse mesh corrected by the accurate one, or the accurate mesh in every evaluation |
| `red_noise` | default false | inflate the errors per segment by the red-noise factor β (needs `error_scale` or `clip_sigma`) |

**Table 61.** Keys of an entry of `fit.lc`.

| Key of `fit.stage1` | Type, default | Meaning |
| --- | --- | --- |
| `terms` | list of `lite`, `apsidal`, `quadratic`; default `[lite]` | terms of every pair |
| `pairs` | map pair → terms | terms per pair |
| `shared_outer` | default false | one outer orbit for all pairs |
| `P3_range_d` | [min, max] | search range of the outer period |
| `clip_sigma`, `max_clip_fraction` | none; 0.1 | outlier rejection |
| `sampler` | `none` (default) or `emcee` | MCMC after the least squares |
| `walkers`, `steps`, `seed` | 48, 4000, 1 | MCMC settings |

**Table 62.** Keys of `fit.stage1`.

| Key of `fit.stage2` | Type, default | Meaning |
| --- | --- | --- |
| `bjd_t0` | BJD_TDB, required | the epoch of $t = 0$ of the configuration |
| `free` | list, required | free parameters: `orbits.<o>.<key>` with key `P`, `M0`, `e`, `omega`, `Omega`, `i`, `ecosw`, `esinw` or `Pdot`; `masses.<m>`; `stars.<m>.<key>` with key `R`, `Teff`, `apsidal_k2` or `rg2`, `rv.gamma`, `rv.offset.<source>` |
| `pairs` | list | pairs whose minima are fitted (default all) |
| `max_sigma_d` | days | leave out minima with a larger error |
| `anchor` | default false | first anchor $P$ and $M_0$ of every eclipsing pair to the minima |
| `priors` | map parameter → `[mean, sigma]` or `{normal: [mean, sigma]}` | Gaussian priors on free parameters (on `e`, `omega` also when `ecosw`, `esinw` are free, and the reverse) |
| `bounds` | map parameter → [low, high] | default: $P$ within 0.5–2 × the start, $e$ 0–0.95, $i$ 0–180°, masses 0.05–20 × |
| `fd_steps` | map parameter → step | finite-difference steps (default: $P$ $10^{-6}$ relative, angles 0.01°, $e$ $10^{-4}$, masses $10^{-4}$ relative) |
| `loss`, `f_scale` | `linear` (default), `soft_l1`, `huber`, `cauchy`; 3 | robust loss from `f_scale` σ |
| `error_scale` | `none`, `global`, `per_source`; default `fit.error_scale` | error correction after a first fit |
| `error_model` | `scale` (default) or `jitter` | factor on the errors or an error floor in quadrature |
| `min_source_minima` | default 5 | sources with fewer minima share one correction |
| `clip_sigma`, `max_clip_fraction` | none; 0.1 | outlier rejection |
| `refits` | 1–10, default 4 | fits in a row while corrections or outliers change |
| `max_iterations` | default 30 | iterations of the optimiser |
| `processes` | default 0 (all) | cores for the derivatives |
| `sampler` | `none` (default) or `emcee` | MCMC after the least squares |
| `walkers`, `mcmc_steps` | 0 (= 2 × free + 2), 300 | MCMC size; the first half of the chain is discarded |
| `mcmc_tau_target`, `mcmc_max_steps` | 0, 5000 | extend the chain until it is this many autocorrelation times long |
| `seed` | default 1 | seed of the MCMC |

**Table 63.** Keys of `fit.stage2`. `bounds`, `fd_steps` and `priors` must name free
parameters; MCMC keys have no effect with `sampler: none` (a warning says so).

## The `scatter` configuration {#cfg-scatter}

A scattering experiment (`swingersim scatter`) uses a configuration of its own with the keys
`name`, `description` and `scatter`:

| Key of `scatter` | Type, default | Meaning |
| --- | --- | --- |
| `target`, `projectile` | body, required | a single star or a binary: `masses` (one or two, M☉), `a_au` (binaries), `e` (0), `e_thermal` (draw $e$ from $f(e) = 2e$), `radii_rsun` (enables collisions), `polytrope_n` (for tidal capture) |
| `v_inf_kms` or `v_inf_vc` | one required | velocity at infinity [km/s] or in units of the critical velocity |
| `b_au`, `b_max_au` or `r_peri_max_au` | one required | fixed impact parameter, largest one, or largest pericentre (with gravitational focusing) |
| `b_distribution` | `area` (default) or `linear` | uniform in $b^2$ or in $b$ |
| `r_start_factor` | default 50 | start distance in units of the largest semi-major axis |
| `t_max_yr`, `max_wall_s`, `check_periods` | none, 600, 5 | limits and the interval of hierarchy checks |
| `gr`, `tidal_capture` | default false | 1PN relativity; tidal energy loss at close passages |
| `n`, `seed` | 100, 1 | number of encounters and random seed |

**Table 64.** Keys of a scattering configuration.

## A complete example {#cfg-example}

The following configuration is a triple with an eclipsing inner binary and an unseen third star.
It runs as it is with `run`, `timings`, `lightcurve` and `stability`, and its `ensemble` and
`fit` sections are used by `ensemble` and `fit` (the table `minima.csv` must exist next to it).

```yaml
topology: custom
hierarchy:
  out: [inner, m3]          # the inner binary and a third star
  inner: [m1, m2]
masses: {m1: 2.0, m2: 1.8, m3: 1.0}
orbits:
  inner: {P: 3.5, e: 0.0, i: 88.0, Omega: 0.0}
  out:   {P: 150.0, e: 0.3, i: 70.0, Omega: 20.0, omega: 45.0, M0: 10.0}
stars:
  m1: {R: 2.1, Teff: 9000}
  m2: {R: 1.8, Teff: 8200}
stellar_model: {kind: mist, age_myr: 100, feh: 0.0}
physics: {gr: true, tides: true, rotation: true}
run: {backend: nbody, t_max_yr: 100, n_output: 2001}
ensemble:
  n: 500
  seed: 7
  directions: [forward, backward]
  require: [stable_ma01]
  vary:
    masses.m3: {loguniform: [0.3, 3.0]}
    orbits.out.P: {loguniform: [50.0, 2000.0]}
    orbits.out.e: {uniform: [0.0, 0.6]}
    orbits.out.i: isotropic
    orbits.out.Omega: angle
  observed:
    inner_i0: {value: 88.0, sigma: 1.0}
    inner_didt0_deg_yr: {value: -0.05, sigma: 0.02}
    "backward:fate": {equals: intact}
fit:
  minima: minima.csv
  time_system: BJD_TDB
  stage2:
    bjd_t0: 2459000.0
    anchor: true
    free: [orbits.inner.P, orbits.inner.M0, orbits.out.P, orbits.out.M0,
           orbits.out.e, orbits.out.omega, masses.m3]
    error_scale: per_source
    error_model: jitter
    clip_sigma: 5.0
  output:
    dir: output/example_triple_fit
    predict: {years: 5}
```
