# Quick start {#quick-start}

This chapter takes you from an installed SWINGERsim to the first results in a few minutes, once
with the web application and once from the command line. Both compute the same thing: the
example 2+2 quadruple `configs/examples/example_2p2.yaml` over 200 years.

## The example system {#qs-system}

```yaml
topology: "2+2"
masses: {m1: 1.2, m2: 1.0, m3: 0.9, m4: 0.8}        # solar masses
orbits:
  A:   {P: 6.74, e: 0.05, i: 85.0, Omega: 0.0,  omega: 90.0,  M0: 0.0}
  B:   {P: 2.68, e: 0.02, i: 60.0, Omega: 40.0, omega: 10.0,  M0: 120.0}
  out: {P: 500,  e: 0.3,  i: 70.0, Omega: 10.0, omega: 120.0, M0: 45.0}
physics: {gr: true}                                  # general relativity on
run: {t_max_yr: 200, n_output: 2000}                 # 200 years, 2000 snapshots
```

Four Sun-like stars form two binaries, A (periods in days: 6.74) and B (2.68), which orbit each
other on the outer orbit `out` (500 d, eccentricity 0.3). Each orbit has a period `P` in days,
an eccentricity `e` and three angles in degrees: the inclination `i` to the plane of the sky
($i = 90°$ is edge-on), the longitude of the ascending node `Omega` and the argument of
periastron `omega`; `M0` is the mean anomaly at the start, i.e. where on the orbit the stars
are at $t = 0$. The stars have no radii given, so they get zero-age main-sequence radii.

## The web application way {#qs-gui}

```bash
swingersim gui
```

The command starts a small web server and opens `http://127.0.0.1:8765` in the browser.

1. On the **Configuration** page choose `configs/examples/example_2p2.yaml` in the list
   *Open…*. The form and the YAML text show the system; the badge at the top right says
   *valid*.
2. Type a name for the draft (e.g. `my_2p2`) and press **Save**. The configuration is now a
   draft in `jobs/drafts/`.
3. In the panel **What to compute** tick **Run** (and, if you like, **Eclipse times (O − C)**),
   then press **Add and run**. The jobs are queued and the worker starts in the background.
4. The **Jobs** page shows the job running with a progress bar. After about a minute it moves
   to *Finished*.
5. On the **Outputs** page open the card of the run. The result page shows the model system,
   when each binary eclipses, the orbital elements over time and the stability ratios.

The pages are described panel by panel in [The web application](#gui); Figures 4, 7, 8 and 9
there show what you will see.

## The command-line way {#qs-cli}

```bash
swingersim run configs/examples/example_2p2.yaml -o output/example.h5
```

After about a minute and a half the program prints a summary (the tutorial run of this
configuration took 80 s on one core). The most important lines are:

| Line of the summary | What it tells you |
| --- | --- |
| `out a = 1.94..1.9418 AU e = 0.2995..0.3007` | the outer orbit keeps its size (about twice the Earth–Sun distance) and shape |
| `i_mut A-B 25.59..46.67 deg` | the angle between the planes of A and B swings between 26° and 47°: the binaries tilt against each other |
| max \|dE/E\| = 3.23e-13 | the energy is conserved to 13 digits: the integration is accurate |
| `hierarchy: configured [[1,2],[3,4]]`, `no change detected` | no star was ejected or changed partner |
| `precession period : 283 yr estimated` | the plane of binary A turns around once in about 283 years |
| `stability ratio Mardling-Aarseth = 3.67 (stable)` | the orbits are far enough apart for a long life (above 1 is stable) |
| `eclipses (observer along z): A 2 %, B 43 %` | B eclipses during 43 % of the run, A only briefly at the start |
| `m1 R = 1.14 Rsun, …, Teff = 6144 K, main sequence` | the stars: radius, temperature, evolutionary phase |

**Table 10.** The key lines of the summary of the example run. The left column quotes the
line as printed, the right column says what it means. Every line of the summary is explained
in [Outputs](#out-summary).

The run has also written its standard products next to the HDF5 file:

```bash
ls output/example*
# example.h5  example.html  example_summary.txt  example_timeseries.csv
# example_orbits.png  example_inclinations.png  example_conservation.png  example_viewer.html
```

Open `output/example.html` (the report) and `output/example_viewer.html` (an animated 3D view of
the run) in a web browser. Both are single files that work without SWINGERsim and without an
internet connection.

## Change one thing {#qs-change}

The typical way of working is to change one parameter or switch one effect on or off and
compare. A prepared variant tilts binary B by 65° against the outer orbit:

```bash
swingersim run configs/examples/example_2p2_zlk.yaml -o output/zlk.h5
```

The eccentricity of B now swings between 0.01 and 0.81 every 14 years: these are ZLK cycles
(von Zeipel–Lidov–Kozai), the exchange of eccentricity and inclination of an orbit perturbed by
a distant inclined body ([Worked examples](#ex-zlk)).

## Eclipse times of the same system {#qs-timings}

```bash
swingersim timings configs/examples/example_2p2.yaml --years 20 -o output/example_timings.csv \
  --plot output/example_oc.png
```

The command finds every eclipse of A and B in the next 20 years, including the light-travel
time across the outer orbit, and writes their O−C against a linear ephemeris (a table and a
figure). Continue with [Workflows](#workflows) for complete tasks, or with
[Worked examples](#examples) for more systems.
