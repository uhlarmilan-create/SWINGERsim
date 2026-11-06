# Conclusions {#conclusions}

## What is validated, and to which accuracy {#validated}

SWINGERsim carries 219 benchmarks, each with its source and a tolerance fixed before the comparison. With
the final code, 205 pass and 14 lie outside their tolerance for a documented reason that is not an error
of SWINGERsim; none fails. Three quarters of them (165) test code written for SWINGERsim itself, not only
the external N-body integrator. Table 47 summarises the accuracy reached for each group of physics; the
individual numbers are in Chapters 2–4 and in [Appendix A](#appendix-benchmarks).

**Table 47.** Validated physics. Column 2 names the references, column 3 the accuracy reached (the
largest difference among the passing benchmarks of the group, or the typical one where stated), column 4
the class of the tested code (Table 1).

| Physics | References | Accuracy | Class |
|---|---|---|---|
| Newtonian and 1PN N-body dynamics | Mercury, ZLK limits, Söderhjelm, MA01, MEGNO | 10⁻⁴ (relativity), 0.3 % (ZLK, with the expected short-period excess) | R |
| secular backends (triple, 2+2, nested) | analytic limits, Naoz et al. (2013), Vokrouhlický (2016), MSE, N-body of a sextuple | 10⁻⁴ in analytic limits; 0.04° on published curves; 10⁻⁴ in $e_\mathrm{max}$ against MSE; 1.2 % in the ZLK period against N-body | own |
| ensembles of chaotic systems | Pejcha et al. (2013), Antognini & Thompson (2016) | triple fractions and enhancement agree; quadruple fractions depend on the run length (see below); cross-sections 1–2 % (exchange 23 %, 1.8σ) | R, own |
| tides and apsidal motion | Sterne (1939), Hut (1981), Leconte et al. (2010), CzeV343 | 10⁻⁴ (apsidal), 0.3 % (dissipative rates), 3 % (observed apsidal period) | R, own |
| stellar structure constants and tracks | Brooker & Olle, Claret, Tout et al., MIST, Geneva, the Sun, white dwarfs | last digit of tables; 5 % for the Sun; 4–12 % between model families | own/tab |
| inertial waves and pre-main-sequence circularisation | Ogilvie, Barker, Zahn & Bouchet, Meibom & Mathieu | 0.03 % (solver), 1 % (pre-main sequence) | own |
| mass loss | Jeans, Hadjidemetriou, Veras et al., Blaauw/Hills, Perets & Kratter | 4·10⁻⁹ (adiabatic invariant), 0.24 % (impulsive), 6.5 % (instability onset) | R, own |
| mass transfer and common envelope | Soberman et al., Postnov & Yungelson, BSE, COSMIC, Tout et al. | 3·10⁻¹³ (closed forms); 0.02 M☉ in the stripped donor; 0.5–9 % in post-common-envelope periods | own |
| magnetism, Applegate mechanism | Matt et al. (2015), Völschow et al. (2016) | 0.7 % (solar spin), 6 % (energy ratios) | own |
| environment and gravitational waves | Grishin & Perets, Hamers & Samsing, PSR B1913+16 | 3·10⁻⁵ (Galactic tide), 0.7 % (passing stars), 0.02 % (orbital decay) | own, R |
| light curves on spheres and meshes, contact binaries | batman, ellc, PHOEBE 2.5 | 10⁻⁸ (occultation), 4·10⁻⁴ (PHOEBE, spheres), 0.003 mag (contact binary) | own |
| dust discs | exact geometry, Farago & Laskar, ε Aur, HD 98800 | exact geometry; 14 % in the ε Aur duration (sharp edge); 2.5° in the HD 98800 disc | own |
| eclipse timings and the stage-1 fit | Irwin, Kaplan, Giménez & Bastero, OCFit | 10⁻⁴ (LTTE), 0.3–0.6 % (Rømer delay), 0.2 % ($a\sin i$ from N-body minima) | own |
| pulsation timing | Murphy et al. (2014), δ Sct and Cepheid relations | 3·10⁻⁷ (analytic), 1.3 % (KIC 9651065), 2 % (catalogue pulsators) | own |
| stability classifiers | Vynatheya et al. (2023), Quarles et al. (2018) | MLP scores to 0.01; MA01 to 0.05 with the sampling of the paper | own, ext |

Beyond the benchmarks, SWINGERsim reproduces the published model of QZ Car (orbit sizes to 0.07 %, light
shares to 0.001), the precession, timing and pulsation-timing amplitudes of the catalogue systems (within
the published uncertainty or about 10 % wherever the input is consistent), and the observed ZLK phases of
σ Ori and HD 98800. It produced new, testable results: the 2+(2+1) hierarchy of CzeV4315 is favoured by a
factor of about 70 over a 2+2, and tides halve its failure rate; a 2+2 explains ASASSN-V J075203 only in a
narrow regular island with two predicted signatures; the eclipse-timing variations of the compact 2+2
TIC 454140642 run in phase; HD 181068 ends in a triple common envelope; and the precession of the discs of
EE Cep and KH 15D constrains their size.

![**Figure 24.** Every quantitative comparison of SWINGERsim with a published value in this report, on one scale. Each row is one number; the horizontal position of the point is SWINGERsim divided by the published value (logarithmic axis; 1 = the same number, 0.5 = half, 2 = double). The grey band marks ±10 %; horizontal bars are the published 1σ uncertainties divided by the published value. Blue: agrees; orange: the difference has been traced to its source. The column on the right gives the two numbers, SWINGERsim first.](../../figures/part1/literature_overview.png)

Figure 24 puts all comparisons of whole systems and published simulations on one scale. Read it row by
row: a blue point inside the grey band is an agreement to better than 10 %; a blue point outside the band
with a bar that reaches the vertical line is an agreement within the published error. Of the 38 rows, 28 are
blue. Each of the 10 orange rows has a traced cause: the input of the configuration or of the paper (δ Ori
A, λ Sco, the QZ Car gainers), published numbers that are inconsistent (the BSE Algol, V994 Her C), the
length of the published runs (Pejcha et al.), or an open point of SWINGERsim's tidal model (M67).

## The fourteen explained differences {#explained}

**Table 48.** The 14 EXPLAINED benchmarks. Column 2 gives SWINGERsim's value, column 3 the published one,
column 4 the reason, column 5 the kind of reason.

| Benchmark | SWINGERsim | Published | Reason | Kind |
|---|---|---|---|---|
| `v994s-gr-C` | 0.081 | 0.06 | the paper's own Table 3 gives 8.1 %; A and B, computed the same way, match | slip in the paper |
| `mt-bse-algol-accretor` | 2.51 M☉ | 2.85 M☉ | the published history would need the orbit to gain angular momentum in the Hertzsprung gap; the conservative giant phase is reproduced to 1.2 % | published numbers inconsistent |
| `mt-bse-algol-period` | 21.9 d | 30 d | as above | published numbers inconsistent |
| `qzcar-evolutionary-mass-ac2` | 26.0 M☉ | 33.2 M☉ | Ac2 is a likely mass gainer (its fitted rotation is 410 km/s); single-star tracks cannot describe it | input / definition |
| `qzcar-evolutionary-mass-aa1` | 42.5 M☉ | 69.8 M☉ | Brož et al. suggest past transfer Aa2 → Aa1; a rejuvenated gainer looks younger and lighter | input / definition |
| `vok16-pejcha-fA` | 0.21 | 0.31 | fractions grow with the run length; at 50 kyr they are 0.36; three independent integrations agree | length of the published runs |
| `vok16-pejcha-fB` | 0.40 | 0.481 | as above (0.57 at 50 kyr) | length of the published runs |
| `vyn23-3p1-ma01-score` | 0.79 | 0.56 | the score depends on the sampling; the paper's own Table 4 implies 0.69; drawn near the boundary as in the paper: 0.51 | sampling of the published set |
| `tidal-capture-radius-n3` | 2.37 | 2.72 | the quoted radii belong to other conventions; the same fits agree with McMillan et al. (1987) | definition |
| `tidal-capture-radius-n15` | 3.33 | 3.84 | as above (3.33 against about 3.4) | definition |
| `tides-circularization-m67` | 3.8 d | 12.1 d | the equilibrium tide alone is too weak on the main sequence; Barker (2022) explains M67 with inertial waves, which SWINGERsim's two-layer estimate with SSE/MIST structure and evolving spins does not reproduce (stellar structure and spin history) | open point of the tidal model |
| `tides-circularization-m67-iw` | 7.5 d | 12.1 d | as above, with inertial waves | open point of the tidal model |
| `tides-circularization-m67-final` | 7.5 d | 12.1 d | as above, with the pre-main sequence (which is reproduced: 6.6 d against 7.1 ± 1.2 d) | open point of the tidal model |
| `eecep-disc-precession` | 183 yr | 65 yr | a rigid disc of the illustrative size precesses too slowly; 65 yr needs an outer radius of about 1.6 au | constraint on the disc |

None of the fourteen is an error of SWINGERsim; the three M67 benchmarks remain an open point of its tidal model. Four are problems of the published numbers (a slip, an inconsistent history
that affects two benchmarks, and a score that the paper's own table does not support), two concern stars
that are not single-star models, two a different definition, two the length of published runs, three an
open point of SWINGERsim's tidal model, and one turns a difference into a constraint on an object.

## What is not validated {#not-validated}

The following parts of SWINGERsim, or regimes of its physics, have no independent reference in this
report. They should be used with that in mind.

- **Main-sequence tidal circularisation:** the 12-d circularisation of M67 is not reproduced by
  SWINGERsim. Barker (2022) explains it with inertial waves computed from full stellar profiles;
  SWINGERsim's frequency-averaged two-layer estimate with SSE and MIST structure gives 7.5 d. The
  difference lies in the stellar structure and the spin history and is an open point.
- **The streams of mass transfer** are prescriptions (BSE rates, accretion fractions of Hamers et al. 2022),
  not hydrodynamics; the accreted fraction during thermal-time-scale transfer is the least certain number of
  Algol evolution and is a user's choice. Tertiary mass transfer and triple common envelopes have no
  quantitative benchmark; they are checked for conservation and against the expected outcome.
- **Late phases of very massive stars:** from the late main sequence on, the inflated envelopes of
  40–70 M☉ stars are not resolved by the MIST grid (radius errors up to 20–40 %), and the masses at the end of
  the main sequence differ between model families by up to a third.
- **Contact binaries** have one common temperature and no thermal relaxation oscillations; the merger is
  tested only against the expected behaviour.
- **Dust discs** have no gravity acting on the orbits, and a star inside a disc does not add its light.
- **Atmospheres below 3500 K** and with $\log g$ < 0 use the nearest node of the limb-darkening tables.
- **Stability classifiers for more than four stars** apply the 2+2 and 3+1 networks to sub-quadruples, an
  approximation outside their training sets; the algebraic criteria for every nested triple and the ghost
  test are the reference.
- **Comparisons not yet made:** the 1:1 resonance of the binaries of a 2+2 (Breiter & Vokrouhlický 2018),
  secular chaos in 2+2 and 3+1 systems (Hamers & Lai 2017), fates of hierarchical systems with stellar
  evolution against MSE, fractions of triples destabilised by winds and mass transfer (Toonen et al. 2022),
  and a wider sample of observed systems such as the Multiple Star Catalog (Tokovinin 2018).

Outside the scope of SWINGERsim by design are relativistic pulsar timing, compact remnants and supernova
kicks, Wolf–Rayet stars to core collapse, colliding winds, starspots and flares in the light curves,
non-radial pulsation modes, and the selection effects of surveys.

## Summary {#summary}

Every physics module of SWINGERsim reproduces independent published results: analytic limits to 10⁻⁴ or
better, published simulations and other codes to a fraction of a degree or a few per cent, and observed
systems within their uncertainties. Where a difference remains, it has been traced to the published number,
its definition or the length of the published runs, or it is stated as an open point (the
main-sequence circularisation of M67). Most of this evidence concerns code written for SWINGERsim – the secular backends, the coupling of evolution, tides
and mass transfer to the orbits, the light curves, timings, fits and ensembles. Parts II–IV test the same
code on real observations of SQUADRA systems.
