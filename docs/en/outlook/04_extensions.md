# Possible extensions of the software {#extensions}

SWINGERsim is complete within its stated scope: hierarchical systems of three to ten stars, from the
pre-main sequence to the end of the AGB, with tides, mass transfer, common envelopes, discs and the
observables of eclipsing systems (*User Guide*, chapter Introduction). The tests suggest the
following extensions. None of them is needed for the results of the test reports.

**Table 2.** Extensions suggested by the tests. Columns:

- *extension*: what is missing;
- *motivated by*: the test in which it would have helped;
- *effort*: a rough estimate of the work — small (days), medium (weeks) or large (months).

| extension | motivated by | effort |
|:--|:--|:--|
| fit of stellar radii, light ratio and third light together with historical eclipse depths (not only times) | QX Cas | medium |
| history-conditioned search: MCMC directly on an inclination history | QX Cas dates | small (existing fit chain, configuration and computing time) |
| synthetic Gaia astrometry: light-centre track sampled with the Gaia scanning law, RUWE and an epoch-astrometry likelihood for DR4 | V907 Sco, QX Cas | medium |
| Hipparcos–Gaia acceleration as an observable of an ensemble or fit | AH Cep | small |
| apsidal constants from detailed structure of stripped and evolved stars along a mass-transfer history | CzeV343 | medium |
| a stage-2 integration that continues mass transfer and a common envelope inside the N-body run | CzeV343 | large |
| resonance-angle diagnostics for near-commensurate 2+2 systems | CzeV343 (3:2) | small |
| selection model of surveys (TESS, ZTF) for population statistics of eclipsing multiples | projects E and O below | large |
