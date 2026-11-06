# Discussion

## Published elements are not portable without anchoring

The most visible result is that the published elements, entered literally into SWINGERsim, produce eclipses
that are off by hours to days after a few years. The orbits are not wrong: the same osculating numbers mean
slightly different motions in two codes that use different coordinates for the hierarchy and different
models for the tidal and rotational forces of the stars. A relative difference of a few 10⁻⁴ in the mean
period is invisible in the light curve of one TESS sector but adds up to hours over thousands of cycles. Two
numbers per pair (period and phase) remove the effect completely. Anyone who transfers a photodynamical
solution into an N-body code should expect this and anchor the solution to the minima first; SWINGERsim does
it automatically when a fit starts from published elements.

## The dynamics is right; the classical apsidal motion tests stellar structure

Where the apsidal motion is driven by the other stars, SWINGERsim reproduces the published integration: BU CMi
to 1 % with point masses, and the full ETV of TIC 120362137 as well as the authors' own model. Where the
stars themselves contribute, the models differ in what they assume. In BU CMi SWINGERsim adds the classical
advance that stellar models predict (2–4 % of the total), and the fit compensates it through the outer
eccentricity (2.6σ). The published integration behaves like SWINGERsim's point-mass model, although the
paper's own analytic estimate of the tidal advance (0.65 and 0.39°/yr) is as large as SWINGERsim's
(0.56 and 0.33°/yr) and the published code can include tides and rotation; whether its tidal term was
smaller or switched off is not stated, and we leave it open. With the apsidal constants free, the minima
want an even weaker classical advance ($k_2 \le 0.0023$ for pair A at Δχ² = 4 against 0.0042 of the
models) but leave the outer eccentricity at 0.266, so the constants do not explain the difference. In V994 Her the paper fitted the apsidal rates
freely, while SWINGERsim predicts them; B and A would need apsidal constants changed by about 20 % and
10 %, plausible for stellar models, C by a factor of 1.76, which no model gives. A companion of C of
0.13–0.31 M☉ on an 80–182-d orbit would do it instead; the present minima of C neither show nor exclude it. For B, Zasche et al.
(2023) found the same themselves: their measured k₂ is 0.83 times the theoretical value, which they
relate to a possibly higher metallicity; SWINGERsim's factor is 0.81. In both systems, therefore, the
differences are not errors of either code but a measurement of the classical apsidal motion: the timings
of the next decade will tell how large it is.

## What the fit can and cannot do

- **Extra errors per source** are the single most effective element of the fit. They bring the reduced χ²
  to 1 in all three systems, measure how much each observer's minima can be trusted (from 1 min for TESS to
  10–18 min for surveys and early ground minima) and stop a few old minima from pulling the solution.
- **Times of minima alone** determine the eclipsing pairs well (V994 Her: within 2σ of the published
  values), but not orbits seen only indirectly: the outer orbit of TIC 120362137, seen through one eclipsing
  pair, and the A–B orbit of V994 Her, seen through light-time waves. There the fit wanders to solutions that
  fit the timings equally well but contradict other data.
- **Priors** do what they should for an orbit that the minima do not determine (TIC 120362137), and they
  reveal the orbits that the minima do determine differently from the papers (BU CMi, V994 Her).
- **Radial velocities** give masses, but they test the velocities as much as the model: in TIC 120362137 they
  agree to their errors, in BU CMi they scatter far beyond their quoted errors, in V994 Her they agree with the
  published masses on their own but are partly absorbed by the A–B orbit in the joint fit.
- **Light curves** fix inclinations and radii where the eclipses are clean, leave the ratio of the radii of
  partially eclipsing stars undetermined without a spectroscopic light ratio, and need any additional signal
  (the oscillation of V994 Her) removed first. Their correlated noise must be taken into account; without it
  the errors were several times too small.
- **MCMC with 1000 steps** gives errors that agree with the linear ones to about 10 % (median), and changed
  the errors of a 400-step run by a few percent. The chains cover about 19 autocorrelation times in total
  and about ten after the burn-in half is dropped, enough for the quoted ranges but not a formal proof of
  convergence.

## The role of the observers

The TESS minima cover at most five years. The apsidal motion of BU CMi (25–29 yr) and of V994 Her (a century
and more), the 2.9-year A–B orbit of V994 Her and the outer orbit of TIC 120362137 are visible only thanks to
ground-based minima that extend the baseline to 15–20 years. In V994 Her the 64 minima of R. Uhlař and those of
FRAM agree with the model to 2–5 minutes, close to TESS. The oldest minima (2002–2004 in V994 Her, the surveys
of 2008–2016 in BU CMi) are the least precise and the most decisive for the long-term terms. Light curves of
the observers were not available to us, only the times of minima derived from them, so the light-curve
comparison uses TESS only; ground-based light curves, e.g. from the SQUADRA observers processed with the VYVAR
pipeline, would allow the same comparison for the eclipses between the TESS sectors.

## Limitations of this test

- Temperatures were taken from the papers; for V994 Her the nodes and the orbit of C are assumptions.
- SWINGERsim's light curves use limb-darkened spheres or Roche surfaces; the pairs in these systems are detached
  (Roche lobes filled to ≤ 57 %), and Roche surfaces changed the results by less than the errors.
- The oscillation of V994 Her was removed as a sum of sine waves, not modelled.
