# Worked examples {#examples}

The examples of this chapter go from one simple system to ensembles, fits, scattering and
evolving stars. The configurations are in `configs/examples/` (the chaos map in
`configs/validation/`); the parameters of examples 1–3 are illustrative, not those of a real
system. The numbers quoted are those printed by the commands; computing times are for one core
of a desktop computer.

## Example 1: a precessing 2+2 quadruple {#ex-2p2}

The configuration `configs/examples/example_2p2.yaml` is the system of the
[Quick start](#quick-start): binary A (2.2 M☉, 6.74 d) is nearly edge-on to us ($i = 85°$),
binary B (1.7 M☉, 2.68 d) is inclined 60° to the sky, and the two orbit each other every 500 d
on an eccentric orbit ($e = 0.3$). From the angles, the mutual inclinations are about 18°
between `out` and A and 29° between `out` and B.

```bash
swingersim run configs/examples/example_2p2.yaml -o output/example_2p2.h5
```

| Line of the summary | Value | Meaning |
| --- | --- | --- |
| `out a`, `out e` | 1.940–1.942 au; 0.2995–0.3007 | the 500-d period with 3.9 M☉ gives $a_\mathrm{out} = 1.94$ au (Kepler's third law); the 0.1 % spread is the osculating wiggle, not a real change |
| `A a`, `A e` | 0.0908 au; 0.044–0.054 | 19.5 R☉, constant to 0.07 %; a small oscillation of $e$ forced by the outer orbit; no ZLK because 18° is below 39.2° |
| `i_mut out-A`, `out-B` | 16.3–18.1°; 28.5–32.9° | each binary precesses around the outer orbit keeping its tilt |
| `i_mut A-B` | 25.6–46.7° | A and B precess at different rates, so the angle between their planes changes |
| max \|dE/E\| | $3.2\times10^{-13}$ | the integration is accurate |
| `precession period` | 283 yr (A), 598 yr (B) estimated | the 200-yr run is shorter than one cycle; a 1000-yr run measures 255 yr and 590 yr |
| `stability ratio` | 3.67 (Mardling–Aarseth), 4.52 (Vynatheya+22) for A; 6.80 and 8.46 for B | comfortably stable |
| `binaries A and B` | $P_A/P_B = 2.515$, nearest 5:2 (+0.6 %) | the periods are not in resonance |
| `eclipses` | A 2 %, B 43 %, both 0 % | in these 200 years A eclipses only briefly at the start, B for almost half of the time |
| `binary A` | apsidal period 25 660 yr from the binary itself, 270 yr measured | without tides only relativity turns the periastron of an isolated binary; in the quadruple the outer orbit makes it circulate a hundred times faster |
| wall time | 80 s | it scales with the shortest period and the length of the run |

**Table 41.** The summary of example 1, line by line. Left: the line of the summary; middle: the
printed value; right: its meaning.

![**Figure 26.** Mutual inclinations of the three pairs of orbits of example 1 over 1000 years
(`run.t_max_yr: 1000`). Horizontal axis: time in years; vertical axis: mutual inclination in
degrees. Blue: A–B; orange: out–B; violet: out–A. The angle between A and B oscillates between
17.8° and 46.7° with a period of about 450 years; the angle of B to the outer orbit oscillates
between about 29° and 34° and that of A between about 16° and 18°, in the same rhythm.](../../figures/example_2p2_mutual.png)

The 450-year cycle of the angle between A and B (Figure 26) is the beat of the two precession
periods, $1/(1/255 - 1/590) \approx 450$ yr: the two planes are alternately on the same side and
on opposite sides of the outer orbit, so the angle between them swings between the difference
and the sum of their tilts to the outer orbit.

![**Figure 27.** Inclination of A (blue) and B (orange) to the plane of the sky over 1000 years.
Horizontal axis: time in years; vertical axis: inclination in degrees. The shaded bands are the
ranges in which each binary eclipses ($\lvert i - 90°\rvert < 6.0°$ for A, $< 9.1°$ for B, from
the default radii). A swings between 54° and 88° and eclipses only near its maxima; B swings
between 42° and 102°.](../../figures/example_2p2_sky_inclination.png)

Figure 27 is what an observer sees of the precession: a binary eclipses only while its
inclination is inside its band, so A eclipses in windows a few decades long every 255 years and
is inside such a window today ($t = 0$), and B eclipses about a third of the time. Inside each
window the eclipse depth changes continuously. The same run with `--backend secular` takes
2.5 s instead of 80 s and gives the mutual inclinations to a few hundredths of a degree.

## Example 2: ZLK cycles {#ex-zlk}

`configs/examples/example_2p2_zlk.yaml` has all nodes aligned, so the mutual inclinations are
differences of the inclinations: B is tilted by 65° to the outer orbit (200 d), inside the ZLK
window of 39.2–140.8°; A by only 5°.

```bash
swingersim run configs/examples/example_2p2_zlk.yaml -o output/zlk.h5
```

```text
  B in out: i_mut0 = 65.0 deg, g1_0 = 90 deg
    quadrupole theory : e = 0.010..0.807, i_mut = 39.7..65.0 deg, e-cycle 18.9 yr (inside ZLK window)
    simulation        : e = 0.007..0.811, i_mut = 40.6..65.1 deg, 6 ZLK cycles of ~14 yr
    stability ratio Mardling-Aarseth = 1.64 (stable), Vynatheya+22 = 1.73 (stable), eps_oct = 0.002, eps_SA = 0.032  (mild hierarchy: secular estimates approximate)
```

![**Figure 28.** ZLK cycles of binary B in example 2. Top: eccentricity of B (orange) and of A
(blue) against time in years; the dashed line is the largest eccentricity of the test-particle
quadrupole approximation, $\sqrt{1-\frac53\cos^2 65°} = 0.84$. Bottom: mutual inclination of B
(orange) and of A (violet) to the outer orbit in degrees; the dotted line marks 39.2°, the lower
edge of the ZLK window. In every cycle the eccentricity of B rises to 0.81 while its inclination
drops to about 41°, then both return.](../../figures/example_zlk.png)

The eccentricity range agrees with the orbit-averaged quadrupole theory, but the cycle is about
25 % shorter than predicted (14 yr against 18.9 yr; Figure 28). The summary warns about this:
$\varepsilon_\mathrm{SA} = 0.032$ means a mild hierarchy (the outer period is only 20 times the
period of B), for which the averaged equations lose accuracy (Luo et al. 2016); with a
two-times wider outer orbit the averaged period agrees to 1.5 %. At the eccentricity maximum the
stars of B approach to $a(1-e) = 4.4$ R☉, three times the sum of their radii.

## Example 3: stars with structure {#ex-tides}

`configs/examples/example_2p2_tides.yaml` is a compact 2+2 of two eclipsing binaries with periods
of 1.8 and 2.9 days, with `physics: {gr: true, tides: true, rotation: true}`: the stars are
extended bodies distorted by tides and rotation, rotating synchronously with their binary.

```bash
swingersim run configs/examples/example_2p2_tides.yaml -o output/example_tides.h5
```

| Line | Value | Meaning |
| --- | --- | --- |
| `tidal+rotational apsidal period` (A) | 64 yr vs ZLK time 949 yr | the distortion turns the periastron of A fifteen times faster than the outer orbit could build up eccentricity, so ZLK cycles of A are suppressed |
| `stars` | $k_2 = 0.0144$, $P_\mathrm{rot}$ 1.8 → 1.8 d | structure from the default polytropes; without dissipation the spin rates do not change |
| `Roche lobe filled up to` | 0.29–0.51 | detached binaries |
| `binary A` | 62.3 yr expected, 59.4 yr measured | dominated by its own stars; the perturbation by B and the outer orbit adds a few per cent |
| `binary B` | 1116 yr expected, 505 yr measured | B is wider, its own distortion is weak; the other binary doubles the rate |
| max \|dE/E\| | $2.2\times10^{-11}$ | conserved because the tidal lag is zero |

**Table 42.** Key lines of the summary of example 3.

In a quadruple, the observed apsidal period of an eclipsing binary need not agree with the value
expected from the stellar structure alone (Table 42). Giving a tidal time lag
(`stars: {m3: {lag_s: 10}}`) switches on tidal dissipation; the summary then lists the
circularisation and synchronisation time scales. For apsidal-motion studies of stars of
2–3 M☉ use `stellar_model.kind: mist` or `sse` (the Claret 2023 values of $k_2$) instead of the
default polytropes, whose $k_2$ is about three times too large for such stars.

## Example 4: into the past, and stability {#ex-past}

The same system as example 1 integrated 200 years into the past needs only a negative end time,
`run: {t_max_yr: -200}` or `--t-max-yr -200`:

```text
integrated to t = -200 yr (2000 outputs, osculating elements)
  i_mut A-B          17.82..  45.32 deg
  eclipses (observer along z): A 27 %, B 0 %, both at once 0 % of the run
```

In the 200 years before today A eclipsed about a quarter of the time and B not at all; in the
next 200 years it is the other way round. The stability check of example 2:

```bash
swingersim stability configs/examples/example_2p2_zlk.yaml
```

```text
nested triples (ratio > 1 = stable):
  A in out (+B)   a_in/a_out = 0.0706, e_in = 0.010, e_out = 0.300, i_mut = 5.0 deg, q_out = 0.773: Mardling-Aarseth 2.38, Vynatheya+22 3.04
  B in out (+A)   a_in/a_out = 0.1029, e_in = 0.010, e_out = 0.300, i_mut = 65.0 deg, q_out = 1.294: Mardling-Aarseth 1.64, Vynatheya+22 1.73
N-body ghost test (100 outer orbits, Newtonian): stable after 100.0 outer orbits, max relative divergence 1.70e-03
```

Binary B, tilted by 65° and in ZLK cycles, is the more critical one; both criteria and the direct
test agree that the system is stable.

## Example 5: an ensemble – V994 Her {#ex-ensemble}

For V994 Her, photometry fixes the periods, eccentricities and inclinations of both binaries
(84° and 85°), but not how the orbits are turned around the line of sight.
`configs/examples/example_v994_ensemble.yaml` draws 200 systems: masses within 2 %,
inclinations within 1°, binary B at 85° or 95° (sense of revolution unknown), its node relative
to A anywhere, the outer orbit oriented at random; every member is integrated 10 000 years
forward and backward with the secular backend.

```bash
swingersim ensemble configs/examples/example_v994_ensemble.yaml -o output/v994/ensemble -j 4
swingersim ensemble-report output/v994/ensemble --plot output/v994/overview.png --x orbits.out.i --y B_didt0_deg_yr
```

```text
forward: 200 members (fractions with 68 % Wilson intervals)
  fate intact           100.0 % (99.5-100.0) [200/200]
  i_mut A_B     now 90.3 deg (16-84 %: 27.7-150.9), range over the run: median 40.0...141.9 deg, members that cross 90 deg: 69.5 %
  A      di/dt now +0.00462 deg/yr (16-84 %: -0.0765..+0.0723); of the members eclipsing now, eclipses end within the run in 92 % (median after 440.4 yr)
  eclipses (time averaged over members): A 48.9 %, B 51.6 %, both 27.9 %
  stable by Vynatheya+22       100.0 % (99.5-100.0) [200/200]
```

The 400 runs take about 10 minutes on one core with Numba. The system is safe in every member,
and its eclipses are temporary: in 92 % of the members the eclipses of A stop within 10 000
years, typically after about 400 years, and both binaries eclipse at once only 28 % of the time –
the doubly eclipsing state observed today is a phase. The present rate of change of the
inclination depends on the unknown orientation of the outer orbit, so a measured rate constrains
it; that is what the conditioning of [Workflows](#wf-ensemble) does. The past looks statistically
like the future, as it must without a measured direction of change.

## Example 6: a fit of times of minima, stage 1 {#ex-fit1}

`configs/examples/example_fit.yaml` is a 2+2 whose minima are synthetic: the table
`example_fit_minima.csv` (50 minima of A, 354 of B over 17.7 years; made by
`tools/make_example_fit_minima.py`) was made with an outer orbit of 2200 d ($e_3 = 0.35$,
$\omega_3 = 1.2$ rad), an apsidal motion of A of 0.0015 rad per cycle and $a\sin i$ of 1.1 au (A)
and 0.95 au (B), at randomly chosen cycles, with errors of 2·10⁻⁴ d. The fit uses the light-time
effect for both pairs, apsidal motion for A and one shared outer orbit:

```bash
swingersim fit configs/examples/example_fit.yaml
```

```text
stage 1: chi^2 = 385.8 for 391 degrees of freedom (reduced 0.99)
  out.P3        2199.8 ± 1.2       d
  out.e3        0.3508 ± 0.0058
  out.w3         1.191 ± 0.015     rad
  A.e          0.24986 ± 0.00040
  A.dw       0.0015010 ± 0.0000028 rad/cycle
  A.asini       1.0965 ± 0.0073    AU
  B.asini       0.9461 ± 0.0032    AU
```

All values are recovered within 1.3 standard deviations, and the reduced $\chi^2$ is 0.99, as it
must be for errors that are correct by construction. The cycles are random on purpose: minima
spaced by a fixed number $N$ of cycles cannot tell an apsidal rate $\dot\omega$ from
$\dot\omega + 2\pi k/N$ (an alias). An earlier version of this table had the minima of A every
40 cycles, and the fit returned 0.1586 rad per cycle ($0.0015 + 2\pi/40$) with a reduced
$\chi^2$ of 14. Real data with a regular cadence (one minimum per season, say) can fall into the
same trap; a reduced $\chi^2$ far above 1 for trusted errors is the warning sign.

## Example 7: two binaries meet {#ex-scatter}

Can a 2+2 quadruple be assembled from two binaries that meet in a young cluster?
`configs/examples/example_scatter_2p2.yaml` sends a binary with $a = 0.5$ au at 1 km/s towards
one with $a = 1$ au (all stars 1 M☉ and 1 R☉, eccentricities drawn from the thermal
distribution, pericentres up to 3 au):

```bash
swingersim scatter configs/examples/example_scatter_2p2.yaml -o output/scatter_2p2 --n 200 -j 0
```

```text
  outcome                               n  fraction    1-sigma interval  sigma [AU^2]
  collision                            97     0.485  [0.450, 0.520]  3.245e+04
  flyby                                56     0.280  [0.249, 0.313]  1.874e+04
  exchange + ionization                19     0.095  [0.076, 0.118]  6357
  triple + single (stable)             12     0.060  [0.045, 0.079]  4015
  single ionization                    10     0.050  [0.037, 0.068]  3346
  exchange                              2     0.010  [0.005, 0.020]  669.2
  triple + single (unstable)            2     0.010  [0.005, 0.020]  669.2
  bound 3+1 (unstable)                  1     0.005  [0.002, 0.013]  334.6
  bound 2+2 (unstable)                  1     0.005  [0.002, 0.013]  334.6
```

The columns are the outcome, the number of encounters, the fraction with its 68 % interval and
the cross-section. At 1 km/s the two binaries are strongly bound once they meet: almost every
close encounter becomes a long resonant interaction, and in half of them two stars touch. A
stable triple with an escaping star is the most common dynamical product (6 %); a bound 2+2 or
3+1 appears only as a temporary, unstable state. Captures into stable quadruples need a third
party or dissipation (tides, gas). The 200 encounters took 22 minutes on two cores.

## Example 8: where is a system stable? A chaos map {#ex-chaos}

`configs/validation/quarles18_circumbinary.yaml` puts a test particle around an equal-mass binary
($a = 1$ au); the map scans its semi-major axis and eccentricity:

```bash
swingersim chaosmap configs/validation/quarles18_circumbinary.yaml -o output/chaosmap_cb \
  --x orbits.out.a --x-range 1.8 4.0 --nx 23 --y orbits.out.e --y-range 0 0.5 --ny 11 --t-max 3000 -j 0
```

![**Figure 29.** MEGNO map of a circumbinary orbit after 3000 years. Horizontal axis: semi-major
axis of the outer orbit in au; vertical axis: its eccentricity. Colour: MEGNO (scale on the
right), 2 for regular motion (blue), above 4 chaotic, 10 marks an escape or a collision (dark
red). The boundary between the regular region (lower right) and the escapes (upper left) follows
a nearly constant pericentre distance $a(1-e)$.](../../figures/chaosmap_circumbinary.png)

The 253 runs take 7 minutes on two cores. At $e = 0$ the orbit is regular from 2.0 au outwards
for this single initial phase; the published limit, which considers all phases, is reproduced in
*Testing, Part I*. The same command maps any parameter pair of a quadruple, e.g. `orbits.out.P`
against `orbits.B.P` near a 3:2 ratio of the periods of the binaries.

## Example 9: what would observers measure? Eclipse timings {#ex-timings}

`configs/examples/example_timings_czev4315.yaml` is the doubly eclipsing quadruple CzeV4315
with one possible mutual orbit (120 d, close to the 119-d signal in the timings of B; the
orientation is illustrative):

```bash
swingersim timings configs/examples/example_timings_czev4315.yaml --years 5 \
  -o output/timings_czev4315.csv --plot output/timings_czev4315.png
```

```text
eclipse timings over 5 yr
  A: 45 primary and 73 secondary eclipses, P = 6.74571156 d, O-C of primaries -840.62 .. +787.62 min (semi-amplitude 814.12 min)
    eclipses missing between 0.330 and 3.158 yr
    no primary eclipses after t = 3.635 yr
  B: 1987 primary and 1988 secondary eclipses, P = 0.91885817 d, O-C of primaries -4.30 .. +3.59 min (semi-amplitude 3.94 min)
```

![**Figure 30.** Synthetic O−C of CzeV4315 for a 120-d mutual orbit over 5 years. Top: binary A,
bottom: binary B; horizontal axis: time in years from $t = 0$; vertical axis: O−C in minutes;
blue dots are primary, orange dots secondary eclipses. B shows a wave of ±4 min with the 120-d
period. A eclipses only in two intervals; its primaries and secondaries are separated by thousands
of minutes and drift.](../../figures/timings_czev4315.png)

B shows the expected 120-day oscillation (light time plus the dynamical effect of A; Figure 30).
A, eccentric and only 18 times faster than the mutual orbit, reacts strongly: its apsides turn
within years, its secondary eclipses move by about ±3 days relative to the primaries, and for
almost three years it does not eclipse at all. CzeV4315 A does eclipse in the data, so this
particular orbit is excluded – the kind of test that an ensemble applies to thousands of orbits.

## Example 10: six stars – V994 Her {#ex-sextuple}

V994 Her is a triply eclipsing sextuple: the doubly eclipsing quadruple A+B (mutual orbit
1062 d) and a third eclipsing binary C on a wide orbit (Zasche et al. 2023b). The preset `2+2+2`
describes it; `configs/examples/example_v994_sextuple.yaml` takes the masses, radii and orbits
of the published solution and marks what is assumed.

```bash
swingersim stability configs/examples/example_v994_sextuple.yaml --no-nbody
swingersim run configs/examples/example_v994_sextuple.yaml -o output/v994_sextuple.h5
```

The stability check applies the criteria to the three nested triples (AB in `out` with C:
Mardling–Aarseth 21.00; A in AB with B: 4.43; B in AB with A: 5.65) and, with the classifiers
installed, to every sub-quadruple. The run (40 s for 30 years) measures apsidal periods of A, B
and C that correspond to 1.42, 3.51 and 0.71 °/yr against the observed 1.86, 3.66 and 1.68 °/yr:
B, whose apsidal constant was derived from these data, agrees to 4 %; A and C turn too slowly
with the assumed theoretical $k_2$, so their structure constants or radii are what the data would
constrain.

## Example 11: an evolving quadruple and its light curve {#ex-evolution}

`configs/examples/example_evolution.yaml` is an illustrative 2+2 of age 1.1 Gyr: binary A with a
2 M☉ primary (4 d), binary B of two solar-type stars (5 d). With `stellar_model.kind: sse` and
`physics.evolution: true` the radii, luminosities and structure constants follow the age of the
stars, and `tidal_lag: hurley2002` takes the tidal dissipation of every star from its envelope.

```bash
swingersim run configs/examples/example_evolution.yaml -o output/evolution.h5
swingersim lightcurve configs/examples/example_evolution.yaml --days 30 -o output/evolution_lc.csv
```

```text
EVENT: roche_overflow at t = 7.445e+07 yr {'age_myr': 1174.45, 'stars': ['m1']}
  max |dE/E| = 2.22e-01  (tidal dissipation: energy decreases)
  max |dL/L| = 1.51e-05
  m1  ... P_rot = 4 -> 3.509 d, Roche lobe filled up to 1.01
      evolving: age 1100 -> 1174.3 Myr, R 3.72 -> 5.96 Rsun, Teff 6761 -> 5057 K, main sequence -> giant branch
```

![**Figure 31.** Example 11. Top: radii of the four stars in R☉ against the time of the run in
Myr (age 1100 Myr at $t = 0$); the 2 M☉ primary m1 (red) leaves the main sequence after 64 Myr
and expands to 6 R☉. Middle: orbital period of A (blue) and rotation period of m1 (dashed red)
in days; the steps are the operator-splitting segments, the drop at the end is the tidal
synchronisation that takes angular momentum from the orbit. Bottom: TESS light curves of 12 days
at the start (age 1100 Myr, grey) and at the end (1174 Myr, red, shifted by −0.25), relative flux
against time in days.](../../figures/example_evolution.png)

The primary of A crosses the Hertzsprung gap and fills its Roche lobe at 1174 Myr; there the run
stops (`roche_overflow`), because this example does not switch on `physics.mass_transfer`
(Figure 31). While the star expands, its spin slows down; once its envelope becomes convective,
the tides synchronise it within a few hundred thousand years and take the angular momentum from
the orbit: the period of A shrinks from 4.0 to 3.5 d. The energy is therefore not conserved (the
tides dissipate it), while the total angular momentum is kept to $1.5\times10^{-5}$. The run takes
3 minutes for 74 Myr. The large ellipsoidal variation of the late light curve is at the edge of
the small-distortion approximation of spherical stars; `stars.<m>.shape: roche` gives the star its
Roche surface.

## Example 12: more kinds of systems {#ex-zoo}

`configs/examples/` holds a ready configuration for each kind of system; run any of them with
`swingersim run <file> -o <name>.h5`.

| Kind of system | Configuration | What happens |
| --- | --- | --- |
| young stars in their gas | `example_gas_2p2.yaml` | the gas shrinks, circularises and aligns the outer orbit, then disperses (secular backend) |
| a triple with a fourth star | `example_3p1.yaml` | precession of the nested orbits |
| a giant beside a binary | `example_giant_triple.yaml` | an AGB tertiary loses mass in its wind; the outer orbit widens |
| a giant feeding a binary | `example_tertiary_transfer.yaml` | the tertiary transfers mass onto the inner binary |
| mass transfer (Algol) | `example_algol_triple.yaml` | the donor overflows onto its companion; the mass ratio reverses |
| common envelope | `example_ce_triple.yaml` | unstable transfer engulfs the binary; the orbit shrinks |
| contact binary | `example_contact_evolution.yaml` | a close binary comes into contact, shrinks and merges |
| pulsating stars | `example_pulsating_triple.yaml`, `example_cepheid_triple.yaml` | δ Scuti, long-period and Cepheid pulsations in the light curve |
| massive stars, Wolf–Rayet | `example_wr_triple.yaml` | a massive tertiary becomes a Wolf–Rayet star; its wind widens the outer orbit |
| winds that destabilise | `example_mass_loss.yaml` | mass loss pushes a triple over the stability limit |

**Table 43.** Further example configurations. Left: the kind of system; middle: the file in
`configs/examples/`; right: what the run shows. The physics each needs is described in *Physics
and Methods*, chapter Stars in multiple systems.
