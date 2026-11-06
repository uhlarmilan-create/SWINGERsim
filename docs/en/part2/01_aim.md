---
title: "Testing SWINGERsim"
subtitle: "Part II: Three published multiple stars – BU CMi, TIC 120362137 and V994 Her"
---

**Abstract.** We test whether the N-body simulator SWINGERsim reproduces three well-studied compact
multiple stars whose solutions were published by members and collaborators of the SQUADRA group:
the 2+2 quadruple BU CMi (Pribulla et al. 2023), the 3+1 quadruple TIC 120362137 (Borkovits et al.
2026) and the triply eclipsing sextuple V994 Her (Zasche et al. 2023). For every system SWINGERsim
integrated the published solution forward and backward in time, tested its stability, computed
synthetic TESS light curves and the times of all eclipses, fitted its own N-body model to the
measured times of minima, radial velocities and TESS light curves, and integrated an ensemble of
400–600 systems drawn within the published uncertainties.

All three systems are stable in every test. The published elements, entered literally, give
eclipses that drift by hours to days within a few years, because osculating periods are not the
periods with which eclipses repeat; after the period and phase of each eclipsing pair are adjusted
to the minima ("anchored"), the synthetic TESS light curves agree with the observed ones to
2–4 parts per thousand (ppt) of the flux in the best sectors, about one per cent of the eclipse
depth, and the eclipse times to the scatter of the data. A joint fit of all data, with a realistic
error for every source of minima and the posterior sampled by MCMC, reproduces the published masses
to 1–5 % (BU CMi, TIC 120362137) and the inclinations to 0.15°. It also shows three differences, each
of which we traced as far as the published information allows. (1) In BU CMi the published apsidal
periods agree to 1 % with SWINGERsim's point-mass integration. When SWINGERsim adds the apsidal advance
caused by the tidal and rotational bulges of the stars, the fit reproduces the observed apsidal motion by
lowering the outer eccentricity by 2.6σ; with point masses the tension disappears. The published code can
include tidal and rotational effects, and the paper's own analytic estimate of the tidal advance is even
somewhat larger than SWINGERsim's; whether the published integration used a smaller tidal term or none is
not stated, and remains open. (2) In TIC 120362137 the published apsidal period of the 51-d orbit is a
theoretical value, not a measurement; SWINGERsim's 32.4 yr is therefore not in conflict with the data.
(3) In V994 Her the apsidal motion is mostly classical: SWINGERsim predicts it from stellar structure,
while the paper fitted the rates freely. For pair B the paper itself finds an apsidal constant about
0.83 times the theoretical one, close to the factor 0.81 that SWINGERsim needs; pair C would need apsidal
constants 1.8 times larger than stellar models give.

# Aim of the test

SWINGERsim is a forward model: from masses, radii and orbits it computes how a hierarchical
multiple star moves, what its light curve looks like and when its eclipses happen. SQUADRA works
the other way round: it starts from observations, light curves and times of minima measured by
satellites and by observers, and asks which system produced them. Before SWINGERsim can help with
that, it has to pass a simple test: **given a system that has already been solved carefully by
others, does SWINGERsim reproduce what was observed?**

We chose three systems with complete published solutions, long series of times of minima and
TESS light curves, one of each hierarchy that matters for SQUADRA:

- **BU CMi**: two eclipsing binaries orbiting each other (2+2), the most compact known flat doubly
  eclipsing quadruple (outer period 121.8 d);
- **TIC 120362137**: an eclipsing binary with a third star on a 51-d orbit and a fourth on a
  1046-d orbit (3+1), the most compact known quadruple of this type;
- **V994 Her**: three eclipsing binaries, A and B in a 2+2 and the visual companion C
  ((2+2)+2), the first known triply eclipsing sextuple.

For each system we ask five questions:

1. **Stability.** Is the published system dynamically stable, by analytic criteria and by a direct
   N-body test?
2. **Dynamics in time.** What happens to the orbits over a century forward and backward, and do
   the precession periods that SWINGERsim computes agree with the published ones?
3. **Light curves.** Do the synthetic light curves of SWINGERsim, of the whole system and of the
   individual pairs, match the TESS light curves?
4. **Times of minima.** Do the eclipses of SWINGERsim happen when the observers saw them, over the whole
   span of the published minima?
5. **Fit.** Can SWINGERsim's own fit to all these data find a configuration that reproduces them, and how
   does it compare with the published solution? An ensemble of systems within the published errors
   shows how robust the conclusions are.

Success means: the published system is stable in SWINGERsim; the synthetic light curves agree with
TESS to the level of the residuals of the published models (a few ppt); the eclipse times agree
with the measured minima to minutes once the periods are anchored; and the fit reproduces the data
without leaving the range of physically sensible parameters. Where SWINGERsim disagrees with a paper,
we look for the reason rather than for agreement.
