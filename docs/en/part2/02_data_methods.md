# Data and methods

## The three systems

**Table 1.** The test systems and the data used. *Minima*: published times of minima (number,
years, sources). *TESS*: sectors of the light-curve comparison (SPOC: pipeline of the TESS Science
Processing Operations Center; QLP: Quick-Look Pipeline).

| | BU CMi | TIC 120362137 | V994 Her |
|---|---|---|---|
| hierarchy | 2+2: (Aa+Ab)+(Ba+Bb) | 3+1: ((Aa+Ab)+B)+C | (2+2)+2: ((A)+(B))+(C) |
| eclipsing pairs | A 2.94 d, B 3.26 d | Aa–Ab 3.28 d | A 2.08 d, B 1.42 d, C 1.96 d |
| outer orbits | A–B 121.8 d | A–B 51.3 d, AB–C 1046 d | A–B 1062 d; C around AB (not measured) |
| source of the solution | Pribulla et al. (2023), Table 7 | Borkovits et al. (2026), Tables 1–2 | Zasche et al. (2023), Table 3 |
| epoch of the elements (BJD) | 2458491.5 | 2458683.3 | 2459010.0 |
| minima | 183, 2005–2023; HAT, KELT, MASCARA, TESS, ground follow-up, Ogloza et al. (2008), Volkov et al. (2021) | 147, 2019–2024; Supplementary Table S3 / Supplementary Data 6 of the paper | 343, 2002–2022; TESS, observers RU, SD & UM, FRAM, MM, GP |
| radial velocities | 305 of all four stars (Pribulla et al. 2023) | 380 of all four stars (Borkovits et al. 2026) | 128 of the stars of A and B (Lee et al. 2008) |
| TESS sectors | 7 (SPOC 30 min), 34, 61 (SPOC 2 min) | 14, 26, 40–41, 53–54, 74 (SPOC) | 26 (SPOC 2 min), 40, 53 (QLP 10 min) |

The published solutions of all three systems were obtained with the code Lightcurvefactory
(Borkovits et al. 2019, 2020, 2021), which fits the light curves, the radial velocities, the eclipse
timing variations and the spectral energy distribution together, with the radii and temperatures of
the stars taken from stellar isochrones (PARSEC). The orbits are computed in two different ways. For
BU CMi and TIC 120362137 the code integrates the motion of the stars numerically (a *photodynamical*
model); its integrator can include tidal and rotational effects (Borkovits et al. 2019, Appendix A),
but neither paper states whether these terms were used. For V994 Her the perturbations between the
pairs are small, and Zasche et al. (2023, Sect. 3) used analytic Keplerian orbits with a constant,
freely fitted apsidal advance rate for each pair instead of an N-body integration. The orbital elements
in the papers are *instantaneous osculating elements* at the epoch in Table 1, i.e. the Keplerian orbit
that the stars would follow from that moment if the other stars disappeared. SWINGERsim uses the same
kind of elements as initial conditions, with the argument of periastron ω of the secondary
(relative orbit), the convention of all three papers. For V994 Her the papers give no nodes; we
assumed Ω = 0 for all orbits (the authors consider the system close to coplanar) and an orbit of C
around A+B with P = 3000 yr and e = 0.3 (the projected separation is ~300 au); the ensemble varies
these assumptions.

## Terms used in this report

- **Time of minimum (eclipse time).** The moment of the deepest light loss in an eclipse, measured
  by an observer from a light curve: these are **data**. The papers list them with errors, the
  *primary* minimum being the eclipse of the first star of a pair.
- **Model eclipse time.** For every measured minimum SWINGERsim computes the moment of minimum
  projected separation of the two stars of the pair in the N-body integration, plus the time the
  light needs to cross the system to the observer (Irwin 1952). It is the model's version of the
  same event, not a prediction.
- **O − C** ("observed minus calculated"). The measured time minus the time of a strictly periodic
  (linear) ephemeris T = T0 + P E, where E is the cycle number. A strictly periodic binary gives a
  horizontal line. In a multiple star the line bends: the *apsidal motion* moves the primary and
  secondary minima in opposite directions (two mirror-image waves); the *light-travel-time effect*
  (LTTE) of an outer orbit shifts all minima of a pair together; and in compact systems the third
  star changes the period directly (*dynamical* eclipse timing variations, ETV). In this report the
  linear ephemeris of each pair is the mean ephemeris of the model primaries, so that measured points
  and model curves share one reference.
- **Prediction.** A model eclipse time or light curve at a date without data.
- **Osculating and mean period.** The osculating period is that of the Keplerian orbit at one instant;
  the mean period is the average interval between eclipses. In a compact multiple star they differ by
  10⁻⁴–10⁻³ of the period, because the orbit of a pair is perturbed continuously.
- **Anchoring.** Adjusting only the period and the phase (mean anomaly at the epoch) of each eclipsing
  pair so that its mean eclipse period and phase match the minima.
- **Apsidal motion.** Rotation of the line of apsides (the long axis) of an eccentric orbit; its
  period U is the time for one full turn. It is driven by the other stars (*dynamical*), by the tidal
  and rotational bulges of the stars themselves (*classical*), and by general relativity.
- **Stability ratio.** For a pair inside an outer orbit, the ratio of the outer periastron distance to
  the critical distance of Mardling & Aarseth (2001, "MA01") or Vynatheya et al. (2022, "V22"); above
  1 the triple is stable. For quadruples SWINGERsim also uses the machine-learning classifier of
  Vynatheya et al. (2023) and a *ghost test*: two N-body runs differing by 10⁻⁶ au in one semi-major
  axis must stay close.
- **Extra light.** Light in the photometric aperture that does not take part in the eclipses: other
  stars of the system, neighbouring stars, or light removed or added by the pipeline's correction
  for crowding. It changes all eclipse depths by the same factor.
- **ppt**: parts per thousand of the total flux. **rms**: root mean square of the residuals.
  **Median |O − C|**: the typical residual, less sensitive to a few bad points.

## Three model configurations

Every system was computed with three configurations:

- **published**: exactly the elements, masses, radii and temperatures of the paper;
- **anchored**: the published configuration with the period and phase of each eclipsing pair adjusted
  to the minima; all other numbers as published;
- **fitted**: the result of SWINGERsim's fit (below).

Why anchoring is needed is one of the results of the test: osculating periods are not portable. Two
codes that describe the same motion with slightly different coordinates or forces map the same
osculating numbers to slightly different mean periods, and over a few thousand cycles a difference of
10⁻⁴ in the period becomes hours. Anchoring removes this difference with two numbers per pair and
leaves the dynamics untouched.

## The model

SWINGERsim integrates the motion of all stars with full gravity (REBOUND with the IAS15 integrator, Rein &
Liu 2012, Rein & Spiegel 2015), including general relativity and the tidal and rotational deformation
of the stars (REBOUNDx, Tamayo et al. 2020), with the apsidal constants from stellar models (Claret
2023) at the age given by MIST tracks (Choi et al. 2016; Dotter 2016). The light curve sums the light of
all stars: limb-darkened spheres (or Roche surfaces) with the four-coefficient law of Claret (2000) and
coefficients of Claret & Bloemen (2011), with the eclipses of every star by every other star in front of
it, the ellipsoidal variation, reflection and Doppler beaming of the close pairs, and the light-travel
time across the system.

## Comparing light curves

The light curves of the TESS satellite (Ricker et al. 2015) were downloaded with lightkurve (Lightkurve Collaboration 2018) from the SPOC
(Jenkins et al. 2016) and QLP (Huang et al. 2020) pipelines. The model was computed at the observed
times and averaged over the exposure (2, 10 or 30 min). Two things in the TESS photometry are not
properties of the stars and have to be allowed for: extra light in the aperture changes all depths by
the same factor, and slow instrumental trends bend the baseline. We therefore compared the observations
with `a × model + baseline`, where `a` is one scale factor per sector and the baseline is a quadratic
polynomial in time for each continuous segment. These are the only adjusted numbers; the eclipses come
from the model unchanged. The figure of merit is the rms of observed minus model in ppt, after removing
points more than five robust standard deviations away.

## The fit

The fit varies the elements of the N-body model until it matches the data. For every measured minimum
the system is integrated to just before the expected eclipse, and the model eclipse is found from the
motion of the pair. The fit minimises χ² = Σ ((data − model)/σ)² over all data at once:

- **times of minima**, every source (each survey, each observer, each TESS sector) with its own extra
  error ("jitter"), added in quadrature to the quoted errors and determined by the fit so that every
  source has a reduced χ² of 1; minima more than 5 robust standard deviations from the model are left
  out (and taken back if a later iteration brings them closer);
- **radial velocities**, with the masses of the stars, the systemic velocity and a zero point per
  instrument free;
- **light curves**: windows of ±0.3 d around the eclipses of one or two TESS sectors, with the
  inclinations and radii free and a linear calibration per window (scale for extra light, zero point,
  trend). Neighbouring points of a light curve are not independent (red noise); their errors are
  multiplied by the factor β of Pont et al. (2006), measured from the residuals;
- **priors**: Gaussian priors from the papers on the elements that the data constrain only weakly
  (eccentricities and ω of the pairs, the outer orbits), never on the periods of the pairs (osculating
  periods are not portable).

The minimisation uses a trust-region least-squares method (SciPy, Virtanen et al. 2020). The posterior
was then sampled with emcee (Foreman-Mackey et al. 2013; 72 walkers × 1000 steps, the first half
discarded as burn-in). The autocorrelation time of the chains is τ ≈ 52–53 steps, so the whole chain
covers about 19 τ and the kept second half about ten τ per walker. That meets the rule of thumb of
Foreman-Mackey et al. (2013), who suggest running for about ten τ, and is enough for the 16–84 % ranges
quoted here (about 700 independent samples), but it is short of the 50 τ that the emcee online
documentation (tutorial "Autocorrelation analysis & convergence") recommends for a reliable estimate of
τ itself. A run of 400 steps gave errors
that differ from those of 1000 steps by a median of 3–5 % (at most a factor of 2, for the eccentricity
vector of the almost circular pair of TIC 120362137) and medians that differ by less than 0.4σ. The
quoted errors are therefore reliable to a few tens of percent.

To see what each kind of data contributes, the fit was also run step by step: times of minima alone,
then with the extra error per source, with priors, starting automatically from the published elements,
with radial velocities, and with light curves. The steps are summarised for each system.

## Ensembles

An ensemble is a set of systems drawn at random within the published uncertainties (masses,
eccentricities, inclinations, nodes, the outer orbits; for V994 Her also the unknown orbit of C).
Each system is integrated forward and backward, and SWINGERsim records its fate (intact, unstable,
merged), its precession periods and whether its pairs eclipse. Systems are then weighted by how
well they reproduce the published apsidal periods and the condition that the pairs eclipse now. The
effective number of members n_eff measures how many systems carry the weight.
