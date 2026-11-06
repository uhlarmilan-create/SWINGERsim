# Accuracy and limits {#accuracy-and-limits}

This chapter collects what the previous chapters say about accuracy and scope: how precisely the
numerical methods solve their equations, which physics is not modelled, in which regimes the models
stop being valid, and how every piece of physics was validated.

## Numerical accuracy {#numerical-accuracy}

Table 15 lists the numerical accuracy of the main methods – the error with which SWINGERsim solves
its own equations, not the error of the physics.

**Table 15.** Numerical accuracy. *Method* is the step; *Measure* the quantity compared; *Accuracy*
the achieved level; *Chapter* where the method is described.

| Method | Measure | Accuracy | Chapter |
|---|---|---|---|
| elements → state → elements | elements after the round trip | $10^{-12}$ | Hierarchies and orbits |
| N-body integration (IAS15) | relative energy error, Newtonian | ~$10^{-14}$ | Dynamics |
| N-body with 1PN | 1PN Hamiltonian | $10^{-13}$–$10^{-11}$ | Dynamics |
| forward and backward run | initial elements recovered | $10^{-12}$ | Dynamics |
| secular backends | relative tolerance of DOP853; conserved $H$ checked | $10^{-10}$ (default) | Dynamics |
| SSE interpolation | radius against direct SSE (95th percentile) | 0.8 % MS, 1.7 % HG, 4 % GB | Stars |
| MIST interpolation | radius at the end of the MS, ages (twice the grid spacing) | 1 %, 1–2 % | Stars |
| sphere occultations | flux against Mandel & Agol (2002) | $10^{-8}$ | Observables |
| surface meshes (64 × 128) | flux against exact sphere eclipses | $3\times10^{-4}$ | Observables |
| disc occultation (4000 points) | flux at a sharp edge | ~$10^{-3}$ | Observables |
| synthetic minima in the fit | times against dense timings | 0.06 s | Fitting |
| imposed period change | deviation from the parabola | < 0.05 s | Mass transfer and evolution |
| radial velocities | emission-time mismatch | < 1 ms | Observables |

Numerical accuracy is rarely the limiting factor. The limits come from the physics and from the data.

## Regime limits {#regime-limits}

- **Secular backends** require strong hierarchies ($\varepsilon_\mathrm{SA}$ well below 0.01), no
  mean-motion resonances and no close encounters; in chaotic regions they agree with the N-body
  solution only for a limited time ([Dynamics](#secular-validity)).
- **Chaos** limits every single run to a few Lyapunov times, forward and backward; beyond, only
  ensembles are meaningful ([Dynamics](#megno)).
- **Stellar models**: SSE stops at helium ignition and covers 0.1–8 M☉; MIST covers 0.1–300 M☉ up to the
  departure from the AGB (massive stars: carbon burning); remnants are not evolved; the masses and
  radii of evolved massive stars are model-dependent at the level of tens of per cent
  ([Stars](#mist-model)).
- **Structure constants** from polytropes are rough (a factor of about three for 2–3 M☉ stars); Claret's
  grids cover 0.2–8 M☉ (Claret 2023) and up to 125 M☉ at one metallicity (Claret 2004).
- **Tides**: the equilibrium-tide dissipation is uncertain by large factors; main-sequence
  circularisation of old clusters is not reproduced ([Tides and rotation](#dissipation-from-structure)).
- **Mass transfer** follows BSE recipes for circular binaries; helium stars and wind accretion are not
  modelled ([Mass transfer and evolution](#mass-transfer-limits)).
- **Light curves**: Roche shapes only for circular synchronised binaries; the small-distortion terms of
  the sphere model err by 3–17 % in the ellipsoidal amplitude at $R/a$ = 0.10–0.25 and should be
  replaced by meshes for close pairs ([Observables](#distorted-stars)).
- **Stability tools** refer to about 100 outer orbits; the MLP classifiers are approximations outside
  their training set ([Stability and fates](#ghost-test)).
- **Fits** find local minima; their linear errors hold in straight $\chi^2$ valleys and MCMC errors only
  for converged chains ([Fitting](#fitting-limits)).

## What is not modelled {#not-modelled}

Table 16 lists the effects that lie outside the scope of SWINGERsim, with what replaces them where
anything does.

**Table 16.** Effects outside the scope of SWINGERsim. *Effect* is the physics; *Consequence* what its
absence means for the results and what replaces it where something does.

| Effect | Consequence |
|---|---|
| hydrodynamics (gas flows, common-envelope hydrodynamics, collisions as fluid events) | gas enters only parametrically (damping of orbits during formation); common envelopes and mergers through energy and mass recipes |
| magnetohydrodynamics (dynamos, resolved magnetic fields) | magnetism enters through parametrised braking and activity laws calibrated on MHD and observations |
| gas drag in scattering experiments | tidal capture is included, drag is not |
| eccentric or asynchronous distorted binaries in the mesh light curves | Roche geometry assumes circular synchronised orbits; the PHOEBE hybrid handles eccentric distorted binaries |
| post-AGB stars, white dwarfs, neutron stars, black holes as evolving stars | they can be point masses; a common-envelope core is a point mass of constant mass |
| helium stars, wind accretion | a donor stripped to its core ends the history; winds are isotropic and leave the system |
| main-sequence tidal dissipation beyond the equilibrium tide and the two-layer estimate of inertial waves | main-sequence circularisation periods of old clusters too short (M67: 7.5 d instead of 12.1 d; Barker 2022 explains M67 with inertial waves computed from full stellar profiles) |
| stellar spots in the light curves | activity cycles enter the eclipse timings, spots not the light |
| pulsation shapes, amplitudes from the structure, gravity-mode spacings from the structure | modes are sinusoids of fixed amplitude; SPB and γ Dor stars need $P_0$ and $\Delta P$ from the user |
| disc mass and gravity, warps, forward scattering in thin parts | discs are massless slabs or rigidly precessing planes |
| astrometry (photocentre motion, Gaia astrometric orbits) | positions are in the output, but no astrometric observable is computed or fitted |
| supernovae and natal kicks | not modelled |

## The validation principle {#validation-principle}

Every physics module of SWINGERsim is compared with published results before it is used: analytic
formulae, published numerical results, other codes (REBOUND/REBOUNDx examples, MSE, PHOEBE, BSE/COSMIC,
MESA, batman, ellc) and observed systems. Each comparison is one number with a tolerance fixed before
the comparison, and its verdict is recorded: passed, explained (outside the tolerance for a documented
reason that is not an error of SWINGERsim – for example published numbers that contradict each other,
or a stated limitation of the physical model), or failed. The tolerance is not
widened to hide a disagreement. The one-sentence summaries in the previous chapters ("reproduces … to
…") are the outcomes of these comparisons; the comparisons themselves, with their sources, inputs and
results, are documented in *Testing SWINGERsim, Part I: Validation against the literature*. Parts II–IV
of the testing documents apply the full chain to published multiple stars, to eclipsing binaries with
changing inclinations and to special cases.

## Implementation {#implementation-accuracy}

The accuracy figures of Table 15 are measured by SWINGERsim's own test suite and benchmark collection
(reproduction of published results), which are part of the code. The integrators whose round-off
behaviour sets the floor are external: IAS15 of REBOUND (Rein & Spiegel 2015) and DOP853 of SciPy.
