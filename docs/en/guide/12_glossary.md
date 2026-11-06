# Glossary {#glossary}

Terms used in SWINGERsim and in this guide, explained for a reader from another field of
astrophysics or a student starting with multiple stars. Within each group the terms are sorted
alphabetically; names in code font are configuration keys or output columns.

## Multiple stars and their hierarchy {#gl-hierarchy}

| Term | Meaning |
| --- | --- |
| 2+2 quadruple | Two binaries (A and B) orbiting their common centre of mass on a wide outer orbit (`out`). |
| 2+(2+1), 2+2+1, 3+1+1 | Quintuples: in 2+(2+1) the fifth star orbits one of the binaries; in 2+2+1 and 3+1+1 it orbits the whole quadruple. |
| 2+2+2, 3+3 | Sextuples: three binaries in a hierarchy (e.g. V994 Her), or two triples orbiting each other. |
| 3+1 quadruple | A close binary (`inner`) orbited by a third star (orbit `mid`), the triple orbited by a fourth star (`out`). |
| centre of mass | The mass-weighted mean position of a group of stars; in a hierarchy each subsystem moves as if its mass were concentrated there. |
| democratic system | Three or more stars without a clear hierarchy, all interacting strongly; short-lived. |
| doubly eclipsing | A system in which two binaries eclipse. |
| hierarchical system | A multiple star whose orbits have very different sizes, so that each is close to a Kepler ellipse slowly perturbed by the others. |
| hierarchy tree | The description of a hierarchy as nested pairs: each orbit joins two children, a star or another orbit. |
| nested triple | An orbit together with the body that perturbs it (the other component of its parent orbit); stability criteria of triples are applied to every nested triple. |
| preset | A named hierarchy (`2+2`, `3+1`, …) that stands for a full tree (Table 46). |
| sub-quadruple | Four stars inside a larger hierarchy that form a 2+2 or 3+1 on their own. |
| triple (2+1) | A binary with a third star on a wider orbit. |

## Orbits and geometry {#gl-orbits}

| Term | Meaning |
| --- | --- |
| argument of periastron (ω, `omega`) | Angle in the orbital plane from the ascending node to the periastron. |
| eccentricity (e) | Shape of the orbit: 0 a circle, close to 1 a very elongated ellipse. |
| edge-on, face-on | An orbit seen from its plane ($i = 90°$, eclipses possible) or from its pole ($i = 0°$). |
| inclination (i) | Angle between the orbital plane and the plane of the sky. |
| invariable plane | The plane perpendicular to the total angular momentum of the system; fixed in time. |
| Jacobi coordinates | Coordinates in which each orbit is the relative motion of two subsystems represented by their centres of mass. |
| longitude of the ascending node (Ω, `Omega`) | Angle in the plane of the sky to the point where the orbit crosses it moving away from the observer; only differences of Ω between orbits are measurable without astrometry. |
| mean anomaly (M, `M0`) | Phase of the motion along the orbit, growing uniformly in time; `M0` is its value at $t = 0$. |
| mean (averaged) elements | Elements averaged over the orbital motion (secular backend). |
| mutual inclination ($i_\mathrm{mut}$) | Angle between two orbital planes; the key quantity of the dynamics. |
| osculating elements | The Kepler ellipse computed from the instantaneous positions and velocities; they oscillate slightly on orbital time scales. |
| periastron, apoastron | Closest and farthest point of an orbit; periastron distance $q = a(1-e)$. |
| prograde, retrograde | Two orbits revolving in the same sense ($i_\mathrm{mut} < 90°$) or in opposite senses. |
| semi-major axis (a) | Half of the longest diameter of the orbital ellipse. |

## Dynamics {#gl-dynamics}

| Term | Meaning |
| --- | --- |
| 1PN, 2.5PN | Orders of the post-Newtonian expansion of general relativity: 1PN gives the periastron advance, 2.5PN the loss of energy by gravitational waves. |
| apsidal motion | Rotation of the line of apsides of an eccentric orbit, caused by stellar distortion, relativity and other bodies; seen as diverging primary and secondary minima. |
| backend | The part of SWINGERsim that integrates the motion: `nbody` or `secular`. |
| chaos, MEGNO | Sensitive dependence on initial conditions; MEGNO (Mean Exponential Growth factor of Nearby Orbits; Cincotta & Simó 2000) tends to 2 for regular motion and grows for chaos. |
| Galactic tide | The tidal field of the Milky Way's disc acting on very wide orbits. |
| hierarchy change | A change of which stars are bound to which (exchange, ejection, break-up, reconfiguration, capture). |
| mean-motion resonance | A near-integer ratio of two orbital periods. |
| N-body integration | Direct numerical solution of the equations of motion for all stars, without averaging. |
| nodal precession | Slow rotation of an orbital plane around the total angular momentum; it changes the inclination to the sky and therefore the eclipses. |
| octupole, quadrupole | Terms of the expansion of the interaction in the ratio of inner to outer orbit size; quadrupole is the leading one. |
| secular dynamics | The slow evolution that remains after averaging over the orbital motions. |
| ZLK cycles (von Zeipel–Lidov–Kozai) | Periodic exchange between eccentricity and inclination of an inner orbit perturbed by a distant inclined body, for mutual inclinations between about 39.2° and 140.8°. |

## Stars, tides and interaction {#gl-stars}

| Term | Meaning |
| --- | --- |
| AGB | Asymptotic giant branch: a late stage of low- and intermediate-mass stars with a huge, cool envelope and a strong wind. |
| apsidal-motion constant ($k_2$) | Measure of the central concentration of a star; sets how strongly its tidal and rotational bulges turn the periastron. |
| Applegate mechanism | A magnetic activity cycle that changes the shape of a star and thereby the period of its binary (Applegate 1992). |
| common envelope | Unstable mass transfer in which gas engulfs both stars; the orbit shrinks until the envelope is ejected or the stars merge. |
| contact binary | A binary whose two stars both overfill their Roche lobes and share one envelope (W UMa stars); fill-out factor $f$ from 0 (inner critical surface) to 1 (outer). |
| effective temperature ($T_\mathrm{eff}$) | Temperature of a black body with the same surface flux; $L = 4\pi R^2\sigma T_\mathrm{eff}^4$. |
| gravity darkening | A rotating or distorted star is brighter where its surface gravity is higher. |
| limb darkening | Decrease of the surface brightness of a star from the centre of its disc to the edge. |
| mass transfer | Flow of gas from a star that fills its Roche lobe (donor) to its companion (accretor). |
| MIST, SSE, BSE | Stellar models: MIST tracks (Choi et al. 2016); single- and binary-star evolution formulae (Hurley et al. 2000, 2002). |
| moment of inertia ($I/(MR^2)$, `rg2`) | Squared gyration radius of a star. |
| pre-main sequence | The youngest stage of a star, still contracting towards the main sequence. |
| pulsating star | A star whose radius and brightness change periodically (δ Scuti, γ Doradus, β Cephei, SPB, RR Lyrae, Cepheids, long-period variables). |
| rejuvenation | An accreting or merged main-sequence star mixes fresh hydrogen into its core and behaves as a younger star. |
| Roche lobe | The region around a star in a binary within which matter is bound to it; its radius after Eggleton (1983). |
| tidal time lag (`lag_s`) | Delay of the tidal bulge behind the companion; sets the tidal dissipation (Mignard 1979; Hut 1981). |
| tertiary mass transfer | Transfer from a third star that fills its Roche lobe onto the inner binary it orbits. |
| Wolf–Rayet star | A hot, massive star that has lost most of its hydrogen envelope and drives a dense wind. |
| ZAMS | Zero-age main sequence: the start of hydrogen burning; the default radii of SWINGERsim. |

## Stability and fates {#gl-stability}

| Term | Meaning |
| --- | --- |
| fate | The end state of a run: intact, a hierarchy change, collision, merger, Roche overflow, or stopped by a limit. |
| ghost test | N-body stability test of Vynatheya et al. (2023): a system is compared with a slightly displaced copy; fast divergence signals instability. |
| Mardling–Aarseth criterion (MA01) | Analytic stability limit of a hierarchical triple (Mardling & Aarseth 2001); ratio above 1 = stable. |
| stability ratio | The ratio of the outer periastron distance to its critical value of a criterion; above 1 stable. |
| Vynatheya+22 criterion (V22) | Updated triple stability criterion (Vynatheya et al. 2022). |

## Observations {#gl-observations}

| Term | Meaning |
| --- | --- |
| BJD_TDB, HJD, JD | Julian dates: barycentric in Barycentric Dynamical Time (the time system of SWINGERsim's measured data); heliocentric and geocentric in UTC (converted with `fit.time_system`). |
| centre of light | The flux-weighted velocity of the visible surface of a star, which a spectrograph measures (`fit.rv.centre: light`). |
| eclipse limit | The smallest inclination at which a pair still eclipses: $\lvert\cos i\rvert < (R_1+R_2)/a$ for a circular orbit. |
| ellipsoidal variation | Brightness change of a binary with tidally distorted stars, two maxima per orbit. |
| ephemeris | Formula predicting times of minima, $T = T_0 + nP$. |
| light curve | Brightness against time in a passband. |
| LTTE (LiTE) | Light-travel-time effect: eclipses come earlier or later as the binary moves on its outer orbit (Irwin 1952). |
| mass function | Combination of masses and inclination measured from the LTTE or radial velocities. |
| minimum (primary, secondary) | Time of an eclipse. In SWINGERsim the primary is the eclipse of the first star of the pair; observers usually call the deeper one primary. |
| O−C, ETV | Observed minus calculated time of minimum (eclipse timing variation) against a linear ephemeris. |
| passband | The wavelength sensitivity of a filter and detector: UBVRI, TESS, Gaia G, BP, RP. |
| photodynamical model | A light curve computed directly from the positions of the stars in an N-body integration. |
| Rossiter–McLaughlin effect | The distortion of the measured radial velocity during an eclipse, when part of the rotating star is hidden. |
| third light | Light of stars that do not take part in an eclipse; it makes the eclipse shallower. |

## Fitting and statistics {#gl-statistics}

| Term | Meaning |
| --- | --- |
| anchoring | Adjusting the period and phase of every eclipsing pair so that the model minima fall on the measured ones, before a stage-2 fit. |
| autocorrelation time (τ) | Number of MCMC steps after which a walker has "forgotten" its position; a chain should be many τ long. |
| $\chi^2$, reduced $\chi^2$ | Sum of squared normalised residuals; divided by the degrees of freedom it is about 1 when the model fits within the errors. |
| conditioning, weights | Weighting ensemble members by how well they reproduce measurements, $w = \exp(-\chi^2/2)$, 0 outside limits. |
| central rate | A present rate of change taken as the mean of the forward and backward fits around $t = 0$. |
| effective number of members ($n_\mathrm{eff}$) | $(\sum w)^2/\sum w^2$ (Kish 1965): how many equally weighted members a weighted ensemble is worth. |
| error scaling, jitter | Correction of quoted errors per source: a factor, or an error floor added in quadrature. |
| evidence (Z) | The mean weight of an ensemble: the probability of the measurements under its hypothesis; ratios compare hypotheses with equal priors. |
| isotropic | Random orientation with $\cos i$ uniform. |
| MCMC, emcee | Markov-chain Monte Carlo sampling of the posterior; emcee is the affine-invariant ensemble sampler (Goodman & Weare 2010; Foreman-Mackey et al. 2013). |
| outlier rejection | Removing measurements that deviate by more than a number of robust standard deviations (`clip_sigma`). |
| prior, posterior | Distribution of parameters before and after taking the measurements into account. |
| red noise (β) | Correlated residuals; β is the factor by which binned residuals scatter more than white noise would (Pont et al. 2006; Winn et al. 2008). |
| seed | Number that fixes the random draws. |
| stage 1, stage 2 | The analytic O−C model of each pair, and the joint N-body model of the whole system, of `swingersim fit`. |
| Wilson interval | Confidence interval of a fraction that behaves well also near 0 and 1 (Wilson 1927). |

## Machine learning {#gl-ml}

| Term | Meaning |
| --- | --- |
| coverage, calibration | Whether stated intervals contain the truth as often as they claim; checked on held-out members. |
| deep ensemble | Several networks trained from different random starts; their spread estimates the uncertainty. |
| emulator | A neural network trained on an ensemble that predicts its summary columns from the drawn parameters. |
| neural network, MLP | A flexible function with many adjustable weights; a multi-layer perceptron is the simplest kind. |
| R² | Score of a regression on held-out members: 1 perfect, 0 no better than the mean. |
| SBI, NPE | Simulation-based inference with neural posterior estimation: a network learns the posterior of parameters from an ensemble. |

## Software {#gl-software}

| Term | Meaning |
| --- | --- |
| configuration | The YAML file that describes one system and what to compute. |
| draft | A configuration written in the web application, kept in `jobs/drafts/`. |
| HDF5, CSV, Parquet | File formats of the run data and of the tables. |
| HTML report | The self-contained web page that every command writes with its results. |
| IAS15 | The adaptive high-order integrator of REBOUND (Rein & Spiegel 2015). |
| job, worker | A job file describes one computation; the worker runs job files from a folder. |
| MESA | The stellar-evolution code used for own grids (Paxton et al. 2011). |
| Numba | A compiler of Python functions; it speeds up the secular equations. |
| operator splitting | Alternating two processes on different time scales: the stars evolve between segments, the orbits are integrated within them. |
| PHOEBE | A code for light curves of eclipsing binaries (Prša et al. 2016), an optional engine of SWINGERsim. |
| REBOUND, REBOUNDx | The N-body code and its library of additional forces (Rein & Liu 2012; Tamayo et al. 2020). |
| YAML | The human-readable text format of configurations and job files. |

## References {#references}

- Applegate J. H. (1992). A mechanism for orbital period modulation in close binaries. *ApJ* 385, 621.
- Choi J., Dotter A., Conroy C., Cantiello M., Paxton B., Johnson B. D. (2016). Mesa Isochrones and Stellar Tracks (MIST). I. Solar-scaled models. *ApJ* 823, 102.
- Cincotta P. M., Simó C. (2000). Simple tools to study global dynamics in non-axisymmetric galactic potentials – I. *A&AS* 147, 205.
- Claret A. (1999). *A&A* 350, 56.
- Claret A. (2004). New grids of stellar models including tidal-evolution constants up to carbon burning. I. From 0.8 to 125 M☉ at Z = 0.02. *A&A* 424, 919.
- Claret A. (2023). Theoretical tidal evolution constants for stellar models from the pre-main sequence to the white dwarf stage. *A&A* 674, A67.
- Claret A., Bloemen S. (2011). Gravity and limb-darkening coefficients for the Kepler, CoRoT, Spitzer, uvby, UBVRIJHK, and Sloan photometric systems. *A&A* 529, A75.
- Eggleton P. P. (1983). Approximations to the radii of Roche lobes. *ApJ* 268, 368.
- Foreman-Mackey D., Hogg D. W., Lang D., Goodman J. (2013). emcee: the MCMC hammer. *PASP* 125, 306.
- Goodman J., Weare J. (2010). Ensemble samplers with affine invariance. *Commun. Appl. Math. Comput. Sci.* 5, 65.
- Hamers A. S. et al. (2015). Secular dynamics of hierarchical quadruple systems: the case of a triple system orbited by a fourth body. *MNRAS* 449, 4221.
- Hamers A. S., Lai D. (2017). Secular chaotic dynamics in hierarchical quadruple systems. *MNRAS* 470, 1657.
- Hamers A. S. et al. (2021). Multiple Stellar Evolution (MSE). *MNRAS* 502, 4479.
- Hurley J. R., Pols O. R., Tout C. A. (2000). Comprehensive analytic formulae for stellar evolution as a function of mass and metallicity. *MNRAS* 315, 543.
- Hurley J. R., Tout C. A., Pols O. R. (2002). Evolution of binary stars and the effect of tides on binary populations. *MNRAS* 329, 897.
- Hut P. (1981). Tidal evolution in close binary systems. *A&A* 99, 126.
- Irwin J. B. (1952). The determination of a light-time orbit. *ApJ* 116, 211.
- Jermyn A. S. et al. (2023). Modules for Experiments in Stellar Astrophysics (MESA): time-dependent convection, energy conservation, automatic differentiation, and infrastructure. *ApJS* 265, 15.
- Kish L. (1965). *Survey Sampling*. Wiley, New York.
- Lacy C. H. S., Helt B. E., Vaz L. P. R. (1999). *AJ* 117, 541 (V907 Sco).
- Lu T., Rein H., Tamayo D., Hadden S., Mardling R., Millholland S. C., Laughlin G. (2023). Self-consistent spin, tidal, and dynamical equations of motion in the REBOUNDx framework. *ApJ* 948, 41.
- Lucy L. B. (1968). The structure of contact binaries. *ApJ* 151, 1123.
- Luo L., Katz B., Dong S. (2016). Double-averaging can fail to characterize the long-term evolution of Lidov–Kozai cycles and derivation of an analytical correction. *MNRAS* 458, 3060.
- Mardling R. A., Aarseth S. J. (2001). Tidal interactions in star cluster simulations. *MNRAS* 321, 398.
- Matt S. P. et al. (2015). The mass-dependence of angular momentum evolution in Sun-like stars. *ApJL* 799, L23.
- Mignard F. (1979). The evolution of the lunar orbit revisited, I. *Moon Planets* 20, 301.
- Mikulášek Z. (2015). Phenomenological modelling of eclipsing system light curves. *A&A* 584, A8.
- Naoz S. (2016). The eccentric Kozai–Lidov effect and its applications. *ARA&A* 54, 441.
- Naoz S., Farr W. M., Lithwick Y., Rasio F. A., Teyssandier J. (2013). Secular dynamics in hierarchical three-body systems. *MNRAS* 431, 2155.
- Paxton B. et al. (2011). Modules for Experiments in Stellar Astrophysics (MESA). *ApJS* 192, 3.
- Peters P. C. (1964). Gravitational radiation and the motion of two point masses. *Phys. Rev.* 136, B1224.
- Pont F., Zucker S., Queloz D. (2006). The effect of red noise on planetary transit detection. *MNRAS* 373, 231.
- Prša A. et al. (2016). Physics of eclipsing binaries. II. Toward the increased model fidelity. *ApJS* 227, 29.
- Reeve D. C., Howarth I. D. (2016). Non-LTE limb-darkening coefficients for hot stars. *MNRAS* 456, 1294.
- Rein H., Liu S.-F. (2012). REBOUND: an open-source multi-purpose N-body code for collisional dynamics. *A&A* 537, A128.
- Rein H., Spiegel D. S. (2015). IAS15: a fast, adaptive, high-order integrator for gravitational dynamics. *MNRAS* 446, 1424.
- Sterne T. E. (1939). Apsidal motion in binary stars. *MNRAS* 99, 451.
- Stone N. C., Leigh N. W. C. (2019). A statistical solution to the chaotic, non-hierarchical three-body problem. *Nature* 576, 406.
- Tamayo D., Rein H., Shi P., Hernandez D. M. (2020). REBOUNDx: a library for adding conservative and dissipative forces to otherwise symplectic N-body integrations. *MNRAS* 491, 2885.
- Tian Y. P., Xiang F. Y., Tao X. (2009). *Ap&SS* 319, 119.
- Tout C. A., Pols O. R., Eggleton P. P., Han Z. (1996). Zero-age main-sequence radii and luminosities as analytic functions of mass and metallicity. *MNRAS* 281, 257.
- Vokrouhlický D. (2016). Secular dynamics of gravitationally bound pair of binaries. *MNRAS* 461, 3964.
- Völschow M. et al. (2016). Eclipsing time variations in close binary systems: planetary hypothesis vs. Applegate mechanism. *A&A* 587, A34.
- Vynatheya P., Hamers A. S., Mardling R. A., Bellinger E. P. (2022). Algebraic and machine learning approach to hierarchical triple-star stability. *MNRAS* 516, 4146.
- Vynatheya P., Mardling R. A., Hamers A. S. (2023). Quadruple-star systems are not always nested triples: a machine learning approach to dynamical stability. *MNRAS* 525, 2388.
- Wilson E. B. (1927). Probable inference, the law of succession, and statistical inference. *J. Am. Stat. Assoc.* 22, 209.
- Winn J. N. et al. (2008). The transit light curve project. IX. Evidence for a smaller radius of the exoplanet XO-3b. *ApJ* 683, 1076.
- Zasche P., Vokrouhlický D., Barlow B. N., Mašek M. (2023). *AJ* 165, 81 (V907 Sco).
- Zasche P. et al. (2023b). V994 Herculis: a unique triply eclipsing sextuple star system. *MNRAS* 520, 3127.
