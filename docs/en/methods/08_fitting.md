# Fitting observed data {#fitting}

The previous chapters run a model forward: from orbits and stars to light curves, velocities and
times of minima. This chapter describes the inverse step – adjusting the model to measured data –
which SWINGERsim performs in a sequence of steps called the *fit chain*. The chain starts from an
analytic description of every eclipsing pair (stage 1), optionally anchors published elements to the
measured minima, and ends with one N-body model of the whole system fitted to all data at once
(stage 2), with robust error treatment, optional Markov-chain Monte Carlo (MCMC) uncertainties and
predictions of future minima.

![**Figure 19.** The fit chain (new figure for this volume). Boxes are steps, arrows the flow of
information. Top row: the data (times of minima by source, radial velocities, light-curve segments,
inclinations or eclipse limits at given epochs) enter stage 1, an analytic model per pair, whose
outer orbit can provide the starting configuration together with published elements and stellar
parameters. Right: optional anchoring of the period and phase of each eclipsing pair. Middle: stage
2, the least-squares fit of one N-body model to all data, alternates with the "fits in a row" (left)
that rescale the errors per source, reject outliers and estimate the red-noise factor, until nothing
changes. Bottom: uncertainties from the covariance matrix and optionally from emcee, and the products
– a fitted configuration usable by every other command and predicted minima with their bands. Blue
boxes are analytic or linear steps; orange boxes call the N-body model many times.](../../figures/methods/fit_chain.png)

Figure 19 is the map of this chapter. Each box corresponds to one of the following sections, and the
loop between the two orange boxes is the robust part of the fit: the error model is not fixed in
advance but estimated from the residuals of the previous fit.

## Data {#fit-data}

The fit uses up to four kinds of data, all on one time axis (BJD in the TDB system; tables in HJD or
JD in UTC are converted with the sky coordinates of the target):

- **Times of minima** (table `fit.minima`): pair, type (primary or secondary), time, error, source
  (instrument, observer or paper) and method. The *source* matters: errors and outliers are treated
  per source.
- **Radial velocities** (`fit.rv`): star, time, velocity, error, source; the model velocity is
  computed at the emission time and, optionally, as the centre-of-light velocity
  ([Observables](#radial-velocities)).
- **Light curves** (`fit.lc`): segments of normalised flux with errors, for example the windows
  around eclipses or binned sectors, in a given passband and exposure time.
- **Inclinations against time** (`fit.inclination`): inclinations of an eclipsing pair measured at
  several epochs, or one-sided constraints ([Inclination data](#inclination-data)).

A *pair* is normally a binary of two stars, but it can be any orbit: the conjunctions of a wider orbit
– the centres of mass of its two components aligned along the line of sight – are the centres of
third-body eclipse events, so a tertiary eclipse can be fitted like a binary minimum.

## Stage 1: analytic O−C model per pair {#stage-1}

Stage 1 is the analytic O−C model that the SQUADRA group fits with OCFit (Gajdoš & Parimucha 2019),
one pair at a time or several pairs with one outer orbit:

$$T(n, j) = T_0 + P\,(n + j/2) + Q\,n^2 + \Delta_\mathrm{LTTE}(t) + D_j(\omega_n),$$

with the cycle number $n$, the minimum type $j$ (0 primary, 1 secondary), the reference epoch $T_0$,
the period $P$, the quadratic coefficient $Q$ (a constant period change, $Q = \tfrac12 P\dot P$),
the light-travel-time term of an outer orbit ([Observables](#ltte)) and the apsidal term $D_j$. The
terms are chosen per pair (`fit.stage1.terms`, `pairs`).

**Apsidal term.** The apsidal motion of an eccentric pair, $\omega_n = \omega_0 + \dot\omega\,n$,
enters through the time of minimum projected separation. At $i = 90°$ this is the conjunction
$\nu + \omega = \pi/2$ (primary) or $3\pi/2$ (secondary), and

$$ D_j = \frac{P}{2\pi}\left[M(\nu_j) + \omega - \left(\frac{\pi}{2} + j\pi\right)\right], $$

the time from the mean to the true conjunction ($M(\nu)$ the mean anomaly belonging to the true
anomaly $\nu$). This is the exact form of the series of Giménez & Bastero (1995) used by OCFit; for
$i < 90°$ the minimum of the projected separation is found numerically (the $e\cot^2 i$ terms of
Giménez & Garcia-Pelayo 1983). The apsidal period is $U = 2\pi P/\dot\omega$.

**Shared outer orbit.** With `shared_outer` the pairs share one outer orbit (common $P_3$,
periastron time $T_3$ and $e_3$, the periastron of the second pair opposite to the first) with
separate $a\sin i$; the ratio of the $a\sin i$ values gives the mass ratio of the pairs,
$m_A/m_B = (a\sin i)_B/(a\sin i)_A$.

**Method.** A weighted least-squares fit started from a grid of circular LTTE fits over the outer
period (and of the anti-phase of primary and secondary minima for the apsidal term); times are fitted
relative to a reference epoch near the middle of the data so that the derivatives keep their
precision. Errors come from the covariance matrix, scaled by $\sqrt{\chi^2_\mathrm{red}}$ when the
reduced $\chi^2$ exceeds 1, or from an MCMC run (below). Outliers can be rejected and errors rescaled
per source as in stage 2. The model agrees with OCFit's series to $10^{-3}$ s at $e = 0.3$, and a fit
of SWINGERsim's own N-body timings of a 2+2 recovers the Keplerian $a\sin i$ of both pairs to 0.3 % and
the outer period to $4\times10^{-5}$ (see *Testing, Part I*).

**Limits.** Stage 1 treats each pair as a Keplerian binary on a Keplerian outer orbit. It cannot
describe dynamical perturbations, precession or the coupling of the pairs in a compact system; its
outer orbit and $a\sin i$ are a starting point for stage 2.

## Anchoring to the present {#anchoring}

Published photodynamical solutions give osculating elements at an epoch. Started from them, the
N-body model of SWINGERsim has a mean inner period that differs slightly from the one measured in the
minima (the elements are osculating, the period seen in the minima is a mean; the other code may
differ in physics), and the computed minima drift away from the measured ones over the years.
*Anchoring* (`fit.stage2.anchor`) removes this drift before any comparison: the period $P$ and the
mean anomaly $M_0$ of every eclipsing pair are adjusted so that the mean linear ephemeris of the
model's eclipses equals that of the measured minima. A coarse step scales $P$ by
$P_\mathrm{data}/P_\mathrm{model}$ and shifts $M_0$ by the delay of the model eclipse nearest
$t = 0$; a fine step fits a straight line through measured − model at the epochs of the minima of
each pair (unbiased by uneven sampling of a curved O−C), repeated until no minimum moves by more
than 0.1 min. Everything else stays as configured. Stage 2 removes the same difference by
construction, because it fits the osculating elements that reproduce the minima; anchoring provides
a good starting point and is useful alone when only a prediction is wanted.

## Stage 2: the joint N-body model {#stage-2}

Stage 2 (`fit.stage2`) fits one N-body model of the whole system to all data at once. The model is
the same as that of a normal run (Newtonian gravity, 1PN relativity, tides, imposed period changes as
configured), so the light-travel-time effect, apsidal motion, dynamical perturbations and precession
come out of one model, and the fitted configuration can be run directly (stability, evolution
forward and backward, light curves, predictions).

**Model times.** For every measured minimum the system is integrated to shortly before the expected
conjunction; the two-body motion of the pair from that state gives the nearest conjunction of the
right type (minimum projected separation; primary when the first component is the farther one); the
N-body integration continues to it and refines it. The observed time adds the light time of the pair's
centre of mass, $t_\mathrm{obs} = t + z_\mathrm{cm}/c$. Minima before $t = 0$ come from a second,
backward integration. All requests of all data types (minima, velocities at their emission times,
inclinations at their epochs) go into one queue ordered by time, so one model evaluation costs one
integration over the span of the data in each direction. The times agree with those of the synthetic
timings at 400 samples per period to 0.06 s.

### Free parameters {#free-parameters}

Table 13 lists the parameters that can be fitted.

**Table 13.** Parameters that stage 2 can fit (`fit.stage2.free`). *Name* is the form in the
configuration; *Meaning* the quantity; *Default bounds* the box within which the optimiser searches
unless `bounds` gives another. Everything not listed as free stays as configured.

| Name | Meaning | Default bounds |
|---|---|---|
| `orbits.<o>.P` | period [d] | 0.5–2 × start |
| `orbits.<o>.M0`, `omega`, `Omega` | mean anomaly, argument of periastron, node [deg] | none |
| `orbits.<o>.e` | eccentricity | 0–0.95 |
| `orbits.<o>.ecosw`, `esinw` | $e\cos\omega$, $e\sin\omega$ (replace $e$, $\omega$ for small $e$) | ±0.95, $e \le 0.95$ |
| `orbits.<o>.i` | inclination [deg] | 0–180 |
| `orbits.<o>.Pdot` | imposed period change [d yr⁻¹] | none |
| `masses.<m>` | mass [M☉] | 0.05–20 × start |
| `stars.<m>.R`, `Teff` | radius [R☉], effective temperature [K] | 0.2–5 ×, 0.5–2 × start |
| `stars.<m>.apsidal_k2`, `rg2` | structure constants | $10^{-5}$–0.74, 0.005–0.4 |
| `rv.gamma`, `rv.offset.<source>` | systemic velocity, zero point of a source [km s⁻¹] | none |

Fitting $k_2$ (or the radii) turns the apsidal motion into a constraint on stellar structure;
fitting `Pdot` lets a parabola in the O−C be fitted together with the LTTE
([Mass transfer and evolution](#semidetached-contact)). A rotation of the whole system about the line
of sight changes none of the data, so with the node of every orbit free only their differences are
measured; SWINGERsim then holds the node of the outermost orbit at its configured value and says so.

**What the data constrain.** Minima alone determine the periods, phases, eccentricities and arguments
of periastron of all orbits and, through the light-time amplitude and the dynamical perturbations, the
mass ratios. The inclinations of the eclipsing pairs and the absolute masses need light curves and
radial velocities; without them they should be fixed at the values of a light-curve or spectroscopic
solution or constrained by priors. Mutual inclinations and nodes are constrained only in compact
systems with a measurable dynamical (not light-time) signal.

### Likelihood {#likelihood}

The fit minimises

$$ \chi^2 = \sum_k \left(\frac{d_k - m_k(\mathbf{p})}{\sigma_k}\right)^2 + \sum_j \left(\frac{p_j - \mu_j}{s_j}\right)^2 , $$

where $d_k$ are the data (times, velocities, fluxes, inclinations), $m_k(\mathbf{p})$ the model at the
parameters $\mathbf{p}$, $\sigma_k$ the errors after the corrections below, and the second sum the
Gaussian priors on parameters $p_j$ with means $\mu_j$ and widths $s_j$. This is $-2\ln$ of a Gaussian
likelihood times Gaussian priors, up to a constant. Angles in priors are wrapped to ±180°. A data
point for which the model does not exist (no conjunction of that type, a star beyond its track)
counts $10^3\sigma$, so the optimiser moves away from such regions.

**Light curves** enter with a calibration per segment, solved linearly at every evaluation and not
part of the free parameters:

$$ F_\mathrm{obs}(t) = a\,F_\mathrm{model}(t) + b + c\,(t - t_\mathrm{mid}) + d\,(t - t_\mathrm{mid})^2 , $$

where $a$ scales the eclipse depths (dilution by third light or a different normalisation), $b$ is the
zero point and $c$, $d$ describe slow trends (polynomial degree 0–2); with `scale: false` the depths
must follow from the stars as given. The degrees of freedom of the fit are reduced by these
coefficients. The model light curve uses the exposure time of the data. Star shapes are chosen per
light curve: with `shapes: auto` a star of a circular two-star binary gets its Roche surface when its
volume radius exceeds a fraction `roche_threshold` (default 0.5) of its Roche lobe – spheres err by
about 2 ppt there ([Observables](#distorted-stars)) – and stays a sphere otherwise.

**Deformed stars in the fit.** Meshes are expensive. With `fidelity: mapped` (default) the fit runs on
a coarse mesh (24 × 48 directions) plus the difference *accurate − coarse* evaluated at the current
solution with the accurate mesh (64 × 128); the difference is recomputed at each new solution until
it changes by less than 0.3σ (at most five passes). This *output space mapping* converges to the
solution of the accurate model with one accurate evaluation per pass; $\chi^2$ and the light curves
reported at the solution are those of the accurate mesh.

### Inclination data {#inclination-data}

A table of inclinations of an eclipsing pair measured at several epochs (light-curve solutions of
separate data sets, as published for V907 Sco by Zasche et al. 2023a) is compared with the inclination of the
pair's orbit in the N-body model at the same times. Eclipses measure only $|\cos i|$, so the model
inclination is folded to 0°–90°, and so are measured values above 90°
([Hierarchies and orbits](#orbital-elements)). A row can also be one-sided:

- `upper` / `lower`: $i_\mathrm{model} \le i$ or $\ge i$ (for example "about 75° or less" from a light
  curve without eclipses);
- `no_eclipse`: the pair did not eclipse at that epoch, i.e. $i_\mathrm{model}$ lies below the
  eclipse limit $\arccos[(R_1+R_2)/r_\mathrm{conj}]$ of the model, with the model radii and the
  separation $r_\mathrm{conj}$ at the nearer conjunction;
- `eclipse`: the pair did eclipse, $i_\mathrm{model}$ above that limit.

A one-sided row contributes $(\mathrm{bound} - i_\mathrm{model})/\sigma$ only when the model violates
it (a *hinge*) and nothing otherwise; $\sigma$ sets how sharply. Satisfied one-sided rows are not
counted as data in the degrees of freedom. A lower limit `min_sigma_deg` on the errors represents the
systematics of the light-curve solutions.

### Optimiser and covariance {#optimiser}

$\chi^2$ is minimised by a trust-region reflective least-squares method with bounds (SciPy
`least_squares`), optionally with a robust loss (`soft_l1`, `huber`, `cauchy`, which reduce the
influence of residuals beyond `f_scale` σ). The Jacobian comes from **central differences**, two
N-body integrations per free parameter and iteration, computed in parallel; forward differences are
biased by the curvature of the residuals (for example the cumulative effect of a period change).
At a bound, or where the model fails on one side (a star beyond the end of its track, an unstable
orbit), a one-sided difference is used and the parameter is reported as lying at the edge of the
model's validity, where its error is not valid. The optimiser stops when $\chi^2$ changes by less
than $10^{-5}$ relative or 0.05 absolute per iteration (far below the change of 1 that moves a
parameter by one σ), or after `max_iterations` (default 30).

The covariance matrix is computed from the singular-value decomposition of the column-scaled Jacobian
rather than by inverting $J^\mathsf{T}J$, which squares the condition number: with parameters of very
different scale (a period known to $10^{-5}$ d, angles to 0.1°) and correlations near 0.999, a
pseudo-inverse of $J^\mathsf{T}J$ silently drops the weakest direction and gives errors far too small.
Directions with a singular value below $10^{-12}$ of the largest are reported as unconstrained. The
errors of the data rows are scaled by $\sqrt{\max(\chi^2_\mathrm{red}, 1)}$; the priors keep their
widths. Strongly correlated parameters ($|r| > 0.95$) are listed in the report – typically $M_0$,
$\omega$ and $P$ of an outer orbit covered by less than about two periods – because their linear
errors hold only while the $\chi^2$ valley is straight.

A closure test on TIC 120362137 (synthetic minima at the 147 epochs of the real TESS minima, 12 free
parameters perturbed) converges from $\chi^2 = 1.9\times10^8$ to 107 ($\chi^2_\mathrm{red} = 0.79$)
in 331 model evaluations and recovers all 12 parameters within 1.5σ of the linear errors (see
*Testing, Part I*).

## Errors per source, jitter and outliers {#errors-outliers}

Times of minima from different sources have errors of different quality. The quoted errors of old
photographic, visual or survey minima are often far below their real scatter, and even formal errors
of TESS minima are below what any model reproduces (reduced $\chi^2$ of 86–1432 in a test on SQUADRA
data). Left as they are, a few hundred old minima with optimistic errors outweigh the precise ones.
SWINGERsim therefore estimates the error model from the residuals in **fits in a row**: fit, correct
the errors and reject outliers, fit again from the solution, until neither the rejected points nor
the corrections change by more than 2 % (at most `refits` fits, default 4). Minima, velocities
(groups per source) and light curves (groups per file and segment) are treated alike.

**Error correction per source** (`error_scale: per_source`; `global` uses one correction for all
minima). Sources with fewer than `min_source_minima` (default 5) points share one group. Each group
$s$ gets either

- a **factor** (`error_model: scale`), $\sigma_k \to f_s\,\sigma_k$ with
  $f_s = \sqrt{\chi^2_s/\nu_s} \ge 1$, or
- a **jitter** added in quadrature (`error_model: jitter`), $\sigma_k \to \sqrt{\sigma_k^2 + j_s^2}$,
  with $j_s \ge 0$ the root of $\sum_{k \in s} r_k^2/(\sigma_k^2 + j_s^2) = \nu_s$,

where $\chi^2_s$ is the contribution of the group, $r_k$ its residuals and $\nu_s$ its degrees of
freedom: the free parameters are shared among the groups in proportion to their points (a light-curve
segment also subtracts its calibration coefficients). A factor keeps the relative weights within a
source; a jitter acts as an error floor – model systematics and timing systematics add in quadrature
to the measurement error – which suits sources of very precise minima. Per source the errors are
never reduced below the quoted ones.

**Outlier rejection** (`clip_sigma`, e.g. 5). A point is left out when its normalised residual exceeds
`clip_sigma` robust standard deviations of its group, $1.4826\,\mathrm{median}|r_k/\sigma_k|$ (at
least 1, so that a source that fits better than its errors loses no good points). The worst points go
first, at most `max_clip_fraction` (default 10 %) of all; a point left out earlier comes back when the
present model fits it within 0.8 × the threshold. Points left out remain in the outputs.

**Red noise** (`fit.lc.red_noise`). Residuals of light curves are often correlated in time
(instrumental trends, oscillations, spots), and white-noise errors then underestimate the parameter
errors. The *red-noise factor* $\beta$ of the time-averaging method (Pont, Zucker & Queloz 2006; Winn
et al. 2008) compares the scatter of residuals averaged in bins of $N$ points with the white-noise
expectation,

$$ \beta_N = \frac{\sigma_N}{\sigma_1 N^{-1/2}\sqrt{M/(M-1)}}, $$

where $\sigma_1$ is the scatter of the unbinned residuals, $\sigma_N$ that of the bin means and $M$ the
number of bins; $\beta$ is the median over bin widths of 20, 40, 60 and 120 min, at least 1. The
errors of each segment are multiplied by $\beta$ in the fits in a row.

## MCMC uncertainties {#mcmc}

The linear errors rely on a locally linear model. In a long, curved valley of $\chi^2$ – typical for
an outer orbit covered by less than about two periods – they come out too small. With
`sampler: emcee` SWINGERsim samples the posterior around the least-squares solution with the
affine-invariant ensemble sampler of Goodman & Weare (2010) as implemented in emcee
(Foreman-Mackey et al. 2013).

- **Posterior.** $\ln p = -\tfrac12[\chi^2_\mathrm{data}/\max(\chi^2_\mathrm{red}, 1) + \chi^2_\mathrm{prior}]$
  inside the bounds, so the data errors are inflated as for the linear errors while the priors keep
  theirs. For orbits fitted with $e\cos\omega$ and $e\sin\omega$ a Jacobian $1/e$ is added: a flat
  prior in $(e\cos\omega, e\sin\omega)$ is a prior $\propto e$ that pushes nearly circular orbits away
  from $e = 0$ (Lucy & Sweeney 1971; Ford 2006; Eastman, Gaudi & Agol 2013), and the correction makes
  the prior flat in $e$ and $\omega$.
- **Walkers** – the parallel chains of the ensemble sampler – number at least $2n + 2$ for $n$ free
  parameters, and even (emcee moves half of them at a time). The number does not depend on the
  computer, so the same seed gives the same chain on any number of cores. They start in a ball of half the covariance matrix around the solution, i.e. spread along
  the valley.
- **Steps.** Each walker makes `mcmc_steps` steps (default 300); the first half is discarded as
  burn-in. Errors are half of the 16–84 % range of the samples.
- **Autocorrelation time.** The integrated autocorrelation time $\tau$ is the number of steps after
  which a chain forgets its state; a chain of $N$ steps contains about $N/\tau$ independent samples per
  walker. SWINGERsim reports $\tau$ of every parameter, its maximum and the ratio steps/$\tau$, and
  warns when the chain after burn-in is shorter than 10 $\tau$.

**Convergence rule of thumb.** A chain at least 50 $\tau$ long is the guidance of the emcee authors
(Foreman-Mackey et al. 2013) for a reliable estimate of $\tau$ and of the posterior quantiles. With
`mcmc_tau_target: 50` the chain is extended by a quarter of its length at a time until it reaches
50 $\tau$ or `mcmc_max_steps`. Every model evaluation is an N-body integration, so 50 $\tau$ can be
expensive. **When the chain is short**, report it as such: give the number of walkers, steps, the
discarded burn-in, $\tau_\mathrm{max}$, steps/$\tau$ and the acceptance fraction; quote the MCMC
intervals as indicative (they are usually too narrow, because the walkers have not yet explored the
tails), and quote the linear errors and the list of strongly correlated parameters alongside.

## Predictions {#predictions}

With `fit.output.predict` the fitted model predicts the eclipses of the chosen pairs for a number of
years ahead (from now, or from the last measured minimum if it is later). Each predicted minimum
carries an uncertainty band: parameter sets are drawn from the emcee samples (when the fit ran emcee)
or from the covariance matrix, each set is integrated once (in parallel), and the band is the 16–84 %
range of the times of the same eclipse across the sets. Only eclipses that occur in the fitted model
(impact parameter below 1) are listed. Bands from the covariance matrix inherit its linear
approximation; in curved valleys the MCMC samples give the more honest band.

## Limits {#fitting-limits}

- The fit finds a local minimum near the starting point; a poor start (a wrong outer period) can end
  in a wrong minimum. Stage 1, anchoring and published solutions provide the start.
- Linear errors are reliable only in straight $\chi^2$ valleys; MCMC errors only for converged
  chains.
- Systematics of the model (missing physics, fixed stellar parameters) enter the error corrections
  per source, not the parameters; a jitter that is large compared with the quoted errors is a sign of
  missing physics.
- The light-curve model inherits the limits of chapter [Observables](#observables): no spots, Roche
  shapes only for circular synchronised binaries.

## Implementation {#implementation-fitting}

**External libraries.** SciPy's `least_squares` (trust-region reflective) as optimiser; emcee
(Foreman-Mackey et al. 2013) as MCMC sampler and for the autocorrelation times; REBOUND/REBOUNDx in
every model evaluation; Astropy for time-system conversions. **Own code.** The stage-1 model
(LTTE, exact apsidal term, shared outer orbit, start grids); anchoring; the stage-2 model with the
event queue that serves minima, velocities and inclinations from one integration per direction; the
parameter mapping, bounds, $e\cos\omega$/$e\sin\omega$ handling and the fixed outermost node; the
central-difference Jacobian with one-sided steps at bounds and at the edge of model validity; the
SVD covariance; the fits in a row with per-source factors and jitters, robust clipping with
re-inclusion and the red-noise factor; the light-curve calibration and the output-space-mapped mesh
fidelity; the hinge rows of one-sided inclination data; the extension of emcee chains to a target
number of $\tau$ and the $1/e$ prior correction; and the predictions with their bands.
