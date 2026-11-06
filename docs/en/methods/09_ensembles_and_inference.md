# Ensembles and inference {#ensembles-and-inference}

The parameters of a real multiple star are never all known. Photometry fixes the inclination of an
eclipsing binary but not the orientation of its node on the sky; the outer orbit is often known only
from a timing signal; masses carry uncertainties. A single simulation answers "what if"; to answer
"how likely", SWINGERsim runs **Monte Carlo ensembles** – many members, each a full run with
parameters drawn from distributions – and conditions them on measurements. This chapter describes
the ensembles, the observation weights and their interpretation as Bayesian quantities, and the
auxiliary machine-learning tools built on ensembles.

SWINGERsim is aimed at **specific systems**: a configuration describes one observed multiple star,
and an ensemble explores the uncertainty of that one system. Emulators and simulation-based inference
(end of the chapter) are auxiliary tools that make such ensembles cheaper to use; they are not
population-synthesis machinery.

![**Figure 20.** A Monte Carlo ensemble in three steps (illustrative). Left: histogram of an unknown
parameter drawn for the members (here an outer period in days; vertical axis: number of members).
Middle: the eccentricity of the inner binary against time in years for a sample of members, one
coloured line each – every member is a full run. Right: histogram of the answer to the question asked
(here the largest eccentricity reached; vertical axis: number of members). The distribution on the
right, not any single run, is the result.](../../figures/concepts/ensemble_idea.png)

Figure 20 shows the logic: uncertain inputs (left) become a bundle of possible histories (middle) and
then a distribution of the quantity of interest (right). Every statement drawn from an ensemble is a
statement about a fraction or a quantile of such a distribution, with a statistical error.

## What is varied {#what-is-varied}

The configuration describes the nominal system; `ensemble.vary` maps parameters (`masses.m1`,
`orbits.out.P`, `orbits.B.Omega`, `stars.m2.R`, …) to distributions, which form the **prior** of the
ensemble. Table 14 lists the distributions.

**Table 14.** Distributions of `ensemble.vary`. *Form* is how it is written; *Distribution* what is
drawn; *Typical use* the parameters it suits. Every distribution can be truncated with `min`, `max`
(a value outside is drawn again).

| Form | Distribution | Typical use |
|---|---|---|
| a number | fixed value | – |
| `uniform: [a, b]` | uniform between $a$ and $b$ | eccentricities, phases |
| `loguniform: [a, b]` | uniform in $\log x$ | periods, separations spanning decades |
| `normal: [μ, σ]` | Gaussian | measured masses, periods, radii |
| `isotropic` or `isotropic: [i1, i2]` | $\cos i$ uniform (between the limits) | inclination of an orbit of unknown orientation |
| `angle` | uniform in 0°–360° | nodes, arguments of periastron, mean anomalies |
| `choice: [x1, x2, …]` | one of the values with equal probability | discrete alternatives |
| `mirror: true` (with an inclination) | $i \to 180° - i$ with probability 1/2 | sense of revolution of an eclipsing binary |

**Angles in the sky frame.** For an orbit of unknown orientation, an isotropic inclination together
with a uniform node gives an isotropic distribution of its normal, and the mutual inclinations follow
from the geometry ([Hierarchies and orbits](#mutual-inclination)). For a measured eclipsing binary,
$i$ is drawn around its measured value, and `mirror: true` adds the unknown sense of revolution;
only the relative sense of two orbits matters, so the first binary can be kept at $i < 90°$. In a
population study (for example Pejcha et al. 2013) the reference plane is the outer orbit instead.

**Requirements.** `require` redraws members that fail a criterion: `stable_ma01`, `stable_v22`
([Stability and fates](#criteria)) or `hierarchical` (every inner apocentre inside the pericentre of
the orbit around it). A requirement truncates the prior; results are then conditional on it.
`ma01_unstable_keep` and `ma01_window` shape training sets for stability classifiers (as Vynatheya et
al. 2023 kept 0.3 of the MA01-unstable draws for 2+2 and 0.2 for 3+1).

**Directions.** With `directions: [forward, backward]` every member is integrated both ways from the
same present state: the backward run tells for how long the system has looked as it does now (for
example since when both binaries eclipse), the forward run what comes next, and survival in both
directions is a strong constraint ([Stability and fates](#survival)).

**Summary and reproducibility.** Every member is reduced to one row: fate, eccentricity extremes,
mutual inclinations and their ranges, eclipse fractions and eclipse changes, apsidal periods,
Roche-lobe filling, stability verdicts (MA01, V22, optionally the MLP classifier and the ghost test)
and the present rates of change. The members depend only on the configuration and `seed`, so an
ensemble can be extended or rerun member by member with another backend.

### Rates of change today {#rates-today}

The most informative measurements of a precessing system are its present rates of change, for
example of the inclination of an eclipsing binary measured from changing eclipse depths (Zasche et al.
2025). SWINGERsim computes `<orbit>_didt0_deg_yr`, the rate of the sky inclination, and
`<orbit>_dtiltdt0_deg_yr`, the rate of the *tilt* $|i - 90°|$, by a linear fit from $t = 0$ over two
periods of the perturbing orbit (at least six outputs) or over `run.rate_window_yr`, so that the
oscillations of the osculating elements with the outer period average out. For an eclipsing binary
the tilt rate is the useful one: eclipses that deepen mean a negative tilt rate whatever the unknown
sense of revolution, whereas the sign of $di/dt$ depends on it. The forward run gives a one-sided rate
forward, the backward run backward; an observer measures a rate over a baseline centred on now, so
the conditioning uses by default the *central* rate, the mean of the two. When an inclination
oscillates quickly the one-sided rates can differ by several degrees per year; `rate_window_yr`
should be about half the observing baseline.

## Statistical error of a fraction {#wilson}

A fraction $p$ estimated from $n$ members has a binomial uncertainty of about $\sqrt{p(1-p)/n}$:
±2.5 % for $p$ = 0.3 and $n$ = 400. SWINGERsim reports the Wilson score interval (Wilson 1927) at
68 %, which stays inside 0–1 also for fractions near 0 or 1. Quantiles (median, 16 % and 84 %) are
reported for continuous columns. Members with a mild hierarchy are often chaotic: two members with
almost the same parameters end differently, so only distributions over many members are meaningful.

## Conditioning on observations {#weights}

Without measurements of how the system changes, forward and backward ensembles of an observed system
are statistically alike. Measurements break this symmetry. Every member $i$ gets the weight of a
Gaussian likelihood

$$ w_i = \exp\left[-\frac{1}{2}\sum_k \left(\frac{x_{ik} - x_k^\mathrm{obs}}{\sigma_k}\right)^2\right] , $$

where $x_{ik}$ is the value of summary column $k$ for member $i$ (for example the tilt rate of B or an
apsidal period), and $x_k^\mathrm{obs} \pm \sigma_k$ the measurement (`ensemble.observed`, entries
`{value, sigma}`). Hard limits (`{min, max}`, or `{equals: …}` for text columns such as the fate) set
$w_i = 0$ for members that contradict what is known for certain – for example
`backward:m1_roche_fill_max: {max: 1}` (the stars of an observed detached binary have not filled their
Roche lobes in the past) or `backward:fate: {equals: intact}`. A member without a finite value of an
observed column gets weight 0. Column names can carry `forward:`, `backward:` or `central:` to say
from which run the value comes.

All statistics of the report are then weighted, both directions use the weight of the same member,
and the drawn parameters are listed before (prior) and after (conditioned) the weighting. Weighted
quantiles are computed from the cumulative weights.

**Interpretation.** This is *importance sampling* of the posterior distribution with the drawn
distributions as the prior: by Bayes' theorem $p(\theta \mid x^\mathrm{obs}) \propto
L(x^\mathrm{obs} \mid \theta)\,\pi(\theta)$, and members drawn from the prior $\pi$ and weighted by the
likelihood $L$ represent the posterior. Two assumptions are built in: the measurement errors are
Gaussian and independent, and the simulated summary values are taken without model error (the scatter
between members represents the unknown parameters, not a noise of the model).

**Effective number of members.** Weighting concentrates the information in fewer members. The
effective number (Kish 1965)

$$ n_\mathrm{eff} = \frac{\left(\sum_i w_i\right)^2}{\sum_i w_i^2} $$

replaces $n$ in the statistical errors (the Wilson intervals are computed with it). If it is small
(below about 20), the measurement selects only a small part of the ensemble and the results are noisy:
more members should be drawn, or the prior narrowed around the region that fits.

![**Figure 21.** Conditioning by importance weights in a toy model (new figure for this volume;
synthetic numbers, not a real system). A parameter θ – here playing the role of an outer period in
days – is drawn log-uniformly for 4000 members; each member's observable is a power law of θ with
30 % random scatter (standing for the unknown angles). (a) Simulated observable against θ, each member
a dot coloured by its weight $w$ (white 0, dark blue 1); the orange line and band are the measurement
−1.4 ± 0.2. (b) Probability density of θ for the prior (grey histogram, all members with equal weight)
and after conditioning (blue outline, weighted). The text gives $n$, $n_\mathrm{eff}$ and the mean
weight $Z$.](../../figures/methods/weights_toy.png)

Figure 21 shows what the weights do. In panel (a) only members whose simulated value falls inside the
orange band receive a high weight; because of the scatter they come from a range of θ, not one value.
Panel (b) shows the result: the broad prior (grey) becomes a narrow conditioned distribution (blue)
around θ ≈ 100–150. Of the 4000 members about 245 carry the weight ($n_\mathrm{eff}$), so the
conditioned statistics have the precision of an ensemble of about 245 members, not 4000.

### The mean weight as Bayesian evidence {#evidence}

The mean weight of an ensemble,

$$ Z = \frac{1}{n}\sum_{i=1}^{n} w_i , $$

is a Monte Carlo estimate of the integral of the likelihood over the prior,
$\int L(x^\mathrm{obs}\mid\theta)\,\pi(\theta)\,d\theta$, i.e. of the *Bayesian evidence* (marginal
likelihood) of the hypothesis that the ensemble represents, up to the constant normalisation
$\prod_k (2\pi\sigma_k^2)^{-1/2}$ of the Gaussian likelihood. Members removed by hard limits count with
$w = 0$; members redrawn by `require` do not count at all, so a requirement is part of the hypothesis.
Its statistical error is about $\mathrm{std}(w)/\sqrt{n}$, and it is noisy when $n_\mathrm{eff}$ is
small.

**Comparing hypotheses.** For two hypotheses $H_1$, $H_2$ – for example a 2+2 quadruple and a
2+(2+1) quintuple, or two ranges of an unknown outer period – conditioned on **the same observations
with the same errors and the same hard limits**, the normalisation cancels and the ratio

$$ B_{12} = \frac{Z_1}{Z_2} $$

is the *Bayes factor*: the factor by which the data change the odds of $H_1$ against $H_2$. With equal
prior odds it is the posterior odds. In the application to CzeV4315 in *Testing SWINGERsim, Part IV*
(one binary with a measured inclination rate, the other with persistent eclipses), a 2+(2+1) model
reached $Z = 0.022$ with $n_\mathrm{eff} = 46$ (1500 N-body members). A 2+2 model with the 119-d
period as its mutual orbit gave no member above a weight of $10^{-12}$ among 1000, i.e. $Z \lesssim
10^{-3}$. The Bayes factor in favour of the quintuple is therefore at least about 20, formally
$10^{13}$. When one hypothesis has (almost) no matching member, quote the upper bound $Z <
1/n$, not the formal value. Three cautions
apply: (1) the evidence depends on the prior – a broad prior that wastes members on regions excluded
by the data lowers $Z$, which is the Bayesian penalty for vagueness, but the comparison is only as
meaningful as the priors are physically motivated; (2) a small $n_\mathrm{eff}$ makes $Z$ itself
uncertain by a large factor; (3) $Z$ from a secular ensemble can be overestimated when members that
fit are chaotic or unstable on longer time scales, which the secular backend does not see – an N-body
check of the best members is advisable. The ensemble report prints $n_\mathrm{eff}$ and $Z$ with
its error $\mathrm{std}(w)/\sqrt{n}$, and it flags the case of fewer than two matching members, for
which the upper bound applies.

**Reporting.** A conditioned result should state the prior (the `vary` distributions and
requirements), the observations and their errors, $n$, $n_\mathrm{eff}$, the weighted median and the
16–84 % range of each quantity, fractions with their Wilson intervals, and, when hypotheses are
compared, both values of $Z$ with their $n_\mathrm{eff}$.

## Emulators {#emulators}

An ensemble maps the parameters of its members to summary columns, at a cost of seconds to hours per
member. A neural *emulator* learns this map from a finished ensemble and evaluates it in microseconds.
The networks are multi-layer perceptrons (scikit-learn, as in Vynatheya et al. 2022, 2023); inputs and
regression targets are standardised (optionally after $\log_{10}$), angles can be mapped to
$(\cos\varphi, \sin\varphi)$, and a known power law can be divided out of a target. A *deep ensemble*
of networks with different seeds (default 5; Lakshminarayanan et al. 2017) gives the spread of the
prediction; 20 % of the members are held out and give the scores. Two uses:

- **Inference.** Draw $10^6$–$10^7$ parameter sets from the ensemble's distributions, predict the
  observed columns and weight every draw with the likelihood of the measurements, in which the variance
  of the emulator (spread of the deep ensemble plus its held-out residual) is added to that of the
  measurement. This is the conditioning above with far more members; precise measurements no longer
  collapse $n_\mathrm{eff}$.
- **Classification.** A binary column (the ghost-test verdict, `fate==intact`) is learned by a
  classifier that returns a probability – for example the stability of a quadruple from its
  scale-free parameters, like the classifiers of Vynatheya et al. (2022, 2023).

An emulator is only as good as its ensemble: outside the sampled ranges and where members are sparse
it extrapolates, which the deep-ensemble spread shows only partly.

## Simulation-based inference {#sbi}

Neural posterior estimation (NPE; Papamakarios & Murray 2016; Greenberg et al. 2019), with the `sbi`
package (Tejero-Cantero et al. 2020), trains a conditional density $q(\theta \mid x)$ on the pairs of
drawn parameters $\theta$ and simulated observables $x$ of an ensemble. Measurement errors enter as
Gaussian noise added to the simulated observables, and hard limits select the training members.
Afterwards the posterior for any measured $x$ is sampled in under a second, without new simulations.
It answers the same question as the importance weights, with the ensemble's distributions as the prior;
weighting needs no training but its $n_\mathrm{eff}$ collapses for precise measurements, while NPE
interpolates between members. A trained network can be wrong where the ensemble has few members, so it
is checked against the weighting and by a coverage test: the fraction of held-out members whose true
values lie inside their 68 % intervals (ideally 0.68), and simulation-based calibration (Talts et al.
2018), in which the rank of the true value among posterior samples is uniform for a calibrated
posterior.

## Implementation {#implementation-ensembles}

**External libraries.** scikit-learn multi-layer perceptrons for the emulators and classifiers; the
`sbi` package (Tejero-Cantero et al. 2020) with PyTorch for neural posterior estimation; the trained
classifiers of Vynatheya et al. (2023) when installed. **Own code.** The drawing of members (all
distributions, truncation, mirroring, requirements, training-set shaping) with deterministic seeds;
the parallel execution with resumption; the reduction of every run to its summary row (including the
central rates of change); the Wilson intervals; the observation weights with Gaussian terms, hard
limits and direction prefixes; $n_\mathrm{eff}$ and weighted quantiles; the emulator layer (input and
target transforms, deep ensembles, retraining of networks stuck on a plateau, emulator variance in the
likelihood); and the SBI interface with noise augmentation, coverage and calibration tests.
