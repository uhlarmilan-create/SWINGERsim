# Independent checks with Gaia

Gaia measures the position of the light centre of a star many times over three years. If the light
centre moves on an orbit, the single-star solution fits worse: the renormalised unit weight error
(RUWE) grows above about 1.4 (the threshold proposed by Lindegren 2018), the astrometric excess noise becomes significant, the parallaxes of DR2
and DR3 (different time spans) may disagree, and for clear cases the Non-Single Star catalogue gives an
astrometric orbit. These are independent of everything used above. A third star of mass m₃ and light
fraction l₃ on an orbit of size a around a pair moves the light centre by about
a (m₃/M − l₃), where M is the total mass.

**Table 12.** The expected motion of the light centre from the models of this report. The distances are
rough values from the brightness of the stars; Gaia's parallax will replace them.

| | V907 Sco | QX Cas | AH Cep |
|---|---|---|---|
| outer orbit | 142 d, a ≈ 0.98 au | 190–1030 d, a ≈ 1.6–4.8 au (fit: 288 d, 2.2 au) | 70 yr, a ≈ 56 au |
| m₃/M − l₃ | 0.16 − 0.02 = 0.14 | 0.04–0.17 (fit: 0.28 − 0.13 = 0.15) | ≥ 0.19 − 0.07 ≈ 0.12 |
| light centre moves by | ≈ 0.14 au | ≈ 0.06–0.8 au (fit: 0.3 au) | ≈ 7 au over 70 yr |
| at a distance of | 0.3–0.5 kpc: ≈ 0.3–0.5 mas | ~3.3 kpc (NGC 7790): 0.02–0.25 mas (fit: 0.1 mas) | ~1 kpc: ≈ 7 mas; the proper motion changes by ~1 mas/yr in 25 yr |
| expected in Gaia | RUWE > 1.4 likely; NSS orbit possible | at the limit | a proper-motion difference Hipparcos–Gaia |

V907 Sco should show the third star clearly; QX Cas is at the limit of what Gaia can see, so a normal
RUWE would not exclude the ensemble; AH Cep, bright enough for Hipparcos, should show the 70-yr orbit as
a change of its proper motion between Hipparcos (1991) and Gaia (2016).

## What Gaia DR3 shows {#gaia-results}

We queried the Gaia archive (DR3, Gaia Collaboration et al. 2023, and for comparison DR2) at the
positions of the three systems with `tools/gaia_check.py`. The astrometry of DR3 is the solution
described by Lindegren et al. (2021).

**Table 13.** Gaia DR3 astrometry of the three systems. Columns:

- *G*: Gaia magnitude;
- *parallax*: DR3 parallax with its error, and the corresponding distance (1/parallax);
- *DR2 parallax*: the value in DR2 and the difference DR3 − DR2 in units of the combined error;
- *RUWE*: values above about 1.4 indicate a poor single-star fit;
- *excess noise*: the astrometric excess noise (AEN, the extra scatter needed to fit the positions)
  and its significance;
- *NSS*: the solution type in the Non-Single Star catalogue, if any.

| system | G | parallax [mas] (distance) | DR2 parallax (difference) | RUWE | excess noise (signif.) | NSS |
|---|---|---|---|---|---|---|
| V907 Sco | 8.55 | 2.279 ± 0.033 (439 pc) | 2.337 ± 0.064 (−0.8σ) | 0.81 | 0.16 mas (11) | SB2, P = 3.978 d |
| QX Cas | 10.12 | 0.268 ± 0.019 (3.7 kpc) | 0.115 ± 0.031 (+4.2σ) | 1.18 | 0.17 mas (39) | – |
| AH Cep | 6.84 | 1.139 ± 0.024 (878 pc) | 1.102 ± 0.034 (+0.9σ) | 1.18 | 0.23 mas (52) | eclipsing binary, P = 1.774737 d |

**How to read Table 13.** Each row is one system; the columns are explained above. None of the three
has a RUWE above 1.4. The excess noise is formally significant in all three, but for bright stars this
is the rule rather than the exception: among the five-parameter solutions of Gaia EDR3 (the same
astrometry as DR3) with G = 9–12, about 97 % have a significant excess noise, with a median of about
0.11 mas (Gaia EDR3 online documentation, Sect. 4.5.2). Much of it is likely to come from imperfect
calibration of bright stars rather than from orbital motion. A significant excess noise of 0.16–0.23 mas,
as here, is therefore weak evidence on its own; the RUWE is the more robust indicator.

- **V907 Sco: a tension.** At the Gaia distance of 439 pc our model moves the light centre by about
  0.14 au, i.e. about 0.32 mas, on the 142-d orbit. Gaia covered about seven such orbits within the
  DR3 time span. A motion of 0.3 mas at G = 8.5 usually raises the RUWE above 1.4, yet it is 0.81. The
  excess noise (0.16 mas) is of the right order, but it does not by itself show an orbit. The NSS
  catalogue gives a double-lined spectroscopic (SB2) solution with P = 3.978 d. This is 5 % longer
  than the eclipse period of 3.776 d, so it probably reflects the inner pair seen with a poorly
  determined period, not the 142-d orbit. Possible reasons for the low RUWE are a smaller light-centre
  motion than modelled (more light from the third star, or a smaller m₃/M), or a cancellation by the
  scanning geometry. The third star itself is established from the eclipse times and radial velocities
  (Zasche et al. 2023), so this does not question it. But an astrometric orbit in a future Gaia
  release (DR4, with epoch astrometry) is a direct test of our m₃.
- **QX Cas: consistent, not decisive.** The RUWE is normal. The excess noise of 0.17 mas is about the
  size of the light-centre motion of our fit (0.1 mas at 3.3 kpc), but for a G = 10 star in a cluster
  that is not evidence. The parallax changed by 4.2σ between DR2 and DR3. The DR2 value (0.115 mas,
  i.e. 8.7 kpc) is incompatible with the cluster NGC 7790 at about 3.3 kpc, while DR3 agrees with it.
  Orbital motion of the light centre can bias a parallax measured over a shorter time span, so this is
  a hint, but DR2 systematics are an equally likely cause.
- **AH Cep: as expected.** The RUWE is normal, as expected for an outer orbit of 70 years that Gaia
  sees as nearly linear motion. The NSS catalogue recognises the eclipsing binary with exactly the
  period of the inner pair (1.774737 d). The test proposed above, a change of the proper motion
  between Hipparcos and Gaia, uses the Hipparcos–Gaia catalogue of accelerations (Brandt 2021). We have
  not yet made it; it is listed in *Open questions and future work*.
