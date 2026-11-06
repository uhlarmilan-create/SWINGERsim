# The web application {#gui}

`swingersim gui` starts a small web server on your computer and opens a page in the browser
from which everything SWINGERsim does can be controlled: writing configurations, choosing what
to compute, running and controlling jobs, and reading the results on pages made for each kind
of result. The page is in English. It does not compute anything itself: it writes job files
that the worker executes (chapter [Installation](#install-worker)) and reads the files the jobs
write. A result computed from the command line therefore looks the same as one started from the
page, and the browser can be closed while jobs run.

The screenshots of this chapter show a real session. The system is V907 Sco, a 3.78-d
eclipsing binary with a third star on a 142-d orbit whose eclipses come and go (Zasche et al.
2023); its runs, light curve, eclipse times and fit are those of *Testing, Part III*. The
ensemble is that of QX Cas from the same document.

## Starting the application {#gui-start}

```bash
swingersim gui                                  # http://127.0.0.1:8765, opens the browser
swingersim gui --port 8800 --no-browser         # another port, no browser
swingersim gui --root ~/swingersim              # another SWINGERsim folder
swingersim gui --host 0.0.0.0                   # reachable from other computers
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--root` | the current folder, or the nearest folder above it that contains `configs/` | the SWINGERsim folder whose `configs/`, `jobs/` and `output/` are used |
| `--host` | `127.0.0.1` | address the server listens on; `0.0.0.0` makes it reachable from the network |
| `--port` | `8765` | port of the server |
| `--no-browser` | off | do not open a browser |

**Table 11.** Options of `swingersim gui`. With `--host 0.0.0.0` anyone on the network can
read the outputs, write configurations and start jobs; use it only on a trusted network. The
server serves only files below `--root`.

## The frame of every page {#gui-frame}

The bar at the top carries the name of the program and four tabs – **Configuration**, **Jobs**,
**Outputs** and **Tools** – which are the four pages of the application. The work goes in one
direction: on *Configuration* you describe a system and tick what to compute, *Jobs* runs the
computations, *Outputs* shows what they gave; *Tools* holds the calculators that need no
configuration. On the right of the bar, **?** opens the help drawer for the current page and
**◐** switches between a light and a dark colour scheme. Below the bar each page starts with a
blue box that says in two sentences what the page is for; **×** closes it. Every section,
field, computation, chart and warning has its own **?**, which opens the help drawer with an
explanation (Figure 5).

## Configuration {#gui-config}

The Configuration page writes and edits configuration files and says what to compute with
them (Figure 4).

![**Figure 4.** The Configuration page with the configuration of V907 Sco opened. Left: the
list of drafts (empty here). Top: the toolbar – *New…*, *Open…*, the topology selector
(`custom`), *Add stars and orbits* and the status badge (*valid*, green). Below it the name bar:
the draft will be saved as `jobs/drafts/v907_sco_published.yaml`; *Save as …*, *Save as new*
and *Delete*; the line under it says where the configuration came from. Middle: the form,
here with the sections *hierarchy* (collapsed), *masses* (three stars) and *orbits* (the orbit
`inner` with its period, semi-major axis, eccentricity and inclination fields; each field has a
**?**, a unit and its default below it, and **×** removes it). Right: the YAML text of the same
configuration (*preview*), and the top of the *Quick look*
panel.](../../figures/guide/gui_config.png)

**Drafts.** The list on the left shows the configurations in `jobs/drafts/`. A click opens a
draft. A configuration from elsewhere (an example in `configs/`, a fitted configuration, a
member of an ensemble) becomes a draft when it is saved.

**Toolbar.** *New…* starts a configuration, empty or from one of five templates (an eclipsing
binary with a third body, a 2+2 quadruple, a 3+1 quadruple, an ensemble for an unseen companion
and a fit of times of minima). *Open…* loads any configuration under `configs/`, a draft or a
submitted job. The **topology** selector sets the hierarchy (the presets of Figure 1 or
`custom`); *Add stars and orbits* fills in the masses and orbits that the topology needs. The
badge on the right shows the state of the configuration: *valid*, *not checked* or the number of
errors.

**Name bar.** The name of the draft (`jobs/drafts/<name>.yaml`). *Save* writes the draft under
this name (the button says *Save as …* while the name differs from the opened file), *Save as
new* asks for a new name and never overwrites an existing draft, *Delete* moves the draft to
`jobs/drafts/.deleted/`. A dot next to the name marks unsaved changes; jobs always use the
saved version.

**Form.** The form is generated from the schema of the configuration, so it always offers every
key with its unit, default and a one-line description; it is the same schema that the command
line uses (chapter [The configuration file](#configuration)). Sections that are not used show
*not set* and say which computations use them (e.g. *stars – used by Light curve, Eclipse times,
Tides*). A section opens with a click; *new key* and *Add* add an entry to a dictionary (a star
to `masses`, an orbit to `orbits`); **×** removes a field. Every change is validated at once
with the same checks as the command line, and errors are shown next to the field and listed
under the YAML text.

**YAML.** The same configuration as text. *preview* shows it; *edit as text* turns it into an
editor, and *Apply text to the form* reads the edited text back into the form. Comments in the
text are not kept by the form; keep annotated files under `configs/`.

![**Figure 5.** The help drawer, opened with the **?** of the section *orbits*. It covers the
right third of the page and has the parts *What it is*, *When to use it*, *What it affects*,
*Typical values* and *Example*; *Insert the example into the configuration* copies the example
into the form. The form on the left stays usable; **×** closes the
drawer.](../../figures/guide/gui_help.png)

The help drawer (Figure 5) is the built-in manual. For a configuration key it explains what
the key means, when it is needed, what it affects, typical values, an example and, where
relevant, the literature; for a computation it says what is computed, what the configuration
needs, how long it takes and what the result page shows; for a chart it explains how to read
it; for a warning it says what it means and what to do.

![**Figure 6.** The right column of the Configuration page for V907 Sco, in two parts placed
side by side. Left part, from the top: the *Quick look* panel with the hierarchy
`((m1 + m2) + m3)`, the stability ratios of the only nested triple (MA01 1.67, V22 2.20,
"stable with a margin") and the eclipse geometry of the pair `inner` (inclination 76.6° against
the eclipse limit 72.8°: "eclipses at t = 0"); then *What to compute* with the six questions
of *What do I want to know?* – the question *What would a photometer see?* is selected and its
answer is shown in the blue box – and the items *Run* (not ticked) and *Light curve* (ticked,
with its options band, start, length, cadence and noise). Right part: *Eclipse times (O − C)*,
ticked, switched to *as an observer* with its options (span, band, cadence, noise; *shapes per
season* and *also a PNG figure* as extra files); *Stability check* and *Chaos map* (not ticked);
*Pulsation periods*, *Ensemble*, *Fit the measurements* and *Scattering experiment* greyed out
with the reason (the configuration has no pulsating star, no `ensemble`, `fit` or `scatter`
section); at the bottom *Add to the inbox* and *Add and run* (disabled until the draft is saved)
and *From the command line*.](../../figures/guide/gui_config_compute.png)

**Quick look** (Figure 6, top left) is computed from the configuration as you type, before
anything is run: the hierarchy in bracket notation, the stability ratios of every nested triple
after Mardling & Aarseth (2001, MA01) and Vynatheya et al. (2022, V22) – above 1 is stable,
green from 1.5 – and, for every pair of two stars, its inclination against the eclipse limit
$\lvert\cos i\rvert < (R_1+R_2)/a$ and the Roche-lobe filling of its stars. It warns early, e.g.
that an observed eclipsing binary has radii too small to eclipse at the given inclination.

**What to compute** lists every computation a configuration can give. Each ticked item becomes
one job of the saved draft, named `<number>_<draft>_<item>`; its results appear on the Outputs
page under the name of the configuration. Items the configuration cannot give are greyed out and
say what is missing. The questions under *What do I want to know?* tick the items that answer
them and say in one sentence how (Table 12).

| Question | Items it ticks | Note |
| --- | --- | --- |
| Is it stable, and for how long? | Stability check, Run | the check gives the verdict of three criteria in a minute; a long run shows what actually happens |
| When will it eclipse? | Run, Eclipse times (exact) | the run shows the eclipsing geometry over time, the eclipse times every eclipse and its O−C |
| What would a photometer see? | Light curve, Eclipse times (as an observer) | a light curve and minima measured on it like real data |
| I have measurements – what system made them? | Fit the measurements | opens the `fit` section to fill in |
| Something unseen perturbs it – what can it be? | Ensemble | opens the `ensemble` section to fill in |
| How did it look in the past? | Run | with a negative `run.t_max_yr` |

**Table 12.** The questions of *What do I want to know?*. Each row is a question (left), the
computations it selects (middle) and the hint shown with it (right).

| Item | Command | Options in the panel | Needs |
| --- | --- | --- | --- |
| Run | `run` | backend (empty: as configured) | masses and orbits |
| Light curve | `lightcurve` | band, start [d after $t = 0$], length [d] (default 27), cadence [min] (default 2), noise [ppm]; a PNG figure | radii and temperatures (from `stars` or the stellar model) |
| Eclipse times (O − C), exact | `timings` | span [yr] (default 10); a PNG figure | radii for the eclipse limit |
| Eclipse times (O − C), as an observer | `minima` | span [yr] (default 3), band, cadence, noise; shapes per season, a PNG figure | as the light curve |
| Stability check | `stability` | outer orbits of the N-body test (default 100) | masses and orbits |
| Chaos map | `chaosmap` | x parameter, x range, x points (30); optional y parameter, range, points | masses and orbits |
| Pulsation periods | `pulsation` | – | a star with `pulsation` |
| Ensemble | `ensemble` | members (empty: `ensemble.n`) | an `ensemble` section |
| Fit the measurements | `fit` | stage (empty: 2 if `fit.stage2` exists) | a `fit` section |
| Scattering experiment | `scatter` | – | a scattering configuration |

**Table 13.** The items of *What to compute*. Columns: the label in the panel, the command the
job runs (chapter [The command line](#command-line)), the options shown when the item is ticked,
and what the configuration must contain.

The two kinds of eclipse times answer different questions. *Exact* takes the moments of
conjunction from the integration, with the light-travel time: the true O−C of the model, without
a light curve or noise, in seconds of computing. *As an observer* computes a synthetic light
curve and measures every minimum on it with the phenomenological model of Mikulášek (2015), with
an error: an O−C of the same kind as one from real minima, in minutes. They differ by what moves
a measured minimum away from the conjunction (asymmetric eclipses of eccentric orbits, the light
of the other stars, pulsations, noise and sampling).

*Add to the inbox* writes the jobs as waiting jobs (`.yaml.hold`), *Add and run* queues them and
starts the worker. *From the command line* shows the commands the ticked items correspond to,
for use in a terminal or a script.

## Jobs {#gui-jobs}

The Jobs page is the control panel of the worker and of the job queue (Figure 7).

![**Figure 7.** The Jobs page while a job runs. Top: the worker bar – *worker running* (green),
the computer (`vm`), its 2 CPUs, the time of the last check, *cores in use: 1 of 2* – and
*Stop worker*. Left: *Drafts* with the draft `v907_sco` (tick box; *run* says which item would
be run), *To the inbox* and *Run*. Right, from the top: the *Inbox* with two waiting jobs
(columns job, command, cores, configuration, output, state) and the buttons *Select one-core
jobs*, *Run selected*, *Hold*, *Remove* and ↻; *Running* with the job
`010_v907_sco_long_run` – a 3000-yr N-body run – its progress bar (*t = 784.5 of 3000 yr ·
elapsed 1 min 1 s · about 2 min 52 s left*), the command it runs and the buttons *Pause*, *Stop*
and *Intermediate results*; *Finished* with the tabs *Done (7)*, *Failed (0)*, *Cancelled (0)*,
*Removed (0)* and one line per finished job with its command, finishing time, duration and the
buttons *HTML report*, *Results* and *Run again*.](../../figures/guide/gui_jobs.png)

**Worker bar.** Whether the worker runs, on which computer, with how many cores, when it last
reported and how many cores the running jobs use. *Start worker* starts it in the background,
so that it keeps running when the browser or the web server is closed; *Stop worker* stops it,
and a running job goes back to the inbox (an ensemble with `resume` continues where it stopped).
A worker started with an older version of SWINGERsim is marked: restart it after an update.

**Drafts.** Tick drafts and press *To the inbox* (waiting jobs) or *Run* (queued, worker
started). This is the same as *Add to the inbox* and *Add and run* on the Configuration page.

**Inbox.** Waiting (`.yaml.hold`) and queued (`.yaml`) jobs. The *cores* column says whether a
job needs one core or all of them. *Run selected* runs only the ticked jobs: other queued jobs
are put on hold, and the page says which. *Select one-core jobs* ticks every job that needs one
core, so that they run side by side. *Hold* puts ticked jobs on hold, *Remove* moves them to
`jobs/cancel/`, ↻ refreshes the list.

**Running.** One block per running job: its command and cores, the progress bar with the
estimated time left (paused time is not counted), the command line, and three buttons. *Pause*
stops the processes of the job (SIGSTOP) and *Continue* lets them go on exactly where they were.
*Stop* ends the job, deletes its results so far (only inside `output/`) and moves the job to
`jobs/cancel/`. *Intermediate results* opens what the job has written so far: a long N-body run
writes its state every 20 s, an ensemble writes one row per finished member.

**Finished.** Done, failed, cancelled and removed jobs on four tabs. *HTML report* opens the
report of the job, *Results* its result page on the Outputs page, *Run again* puts a copy of the
job back into the queue. A failed job shows the end of its log and, for common errors (a missing
file, an invalid configuration, too little memory, the energy-error limit), what the error means
and what to do.

Jobs sent from the page that can use several cores (`ensemble`, `scatter`, `chaosmap`, a stage-2
`fit`, the MESA commands) get `processes: 0` (all cores of the job) unless another value is
set. A single run, light curve or set of eclipse times is one calculation on one core.

## Outputs {#gui-outputs}

Every finished, failed or running job is a **result**. The Outputs page lists them as cards and
opens each on a page made for its kind (Figure 8).

![**Figure 8.** The list of results. Left: *Configurations* – the results grouped by the
configuration they were computed from, with their number (*all configurations 7*,
`example_fit` 1, `qx_cas_ensemble` 1, `v907_sco_fit_lc_v2` 1, `v907_sco_published` 4), and
*everything under output/*. Top: the title *Results*, the search field and ↻. Below it the
filter chips *All*, *Runs*, *Ensembles*, *Light curves*, *Timings & minima*, *Fits*,
*Benchmarks* and *Other*. Each card names its kind (coloured label), job number and date, its
state (*done*), the configuration, one line on what was computed, a small plot (O−C of the
eclipses, the inclination of the inner binary over the run, the light curve, the measured minus
model residuals of a fit, the conditioned distribution of an ensemble parameter) and two to
four facts read from the output, e.g. *Eclipsing (fraction of run) inner 50.5 %*. *Open* opens
the result page.](../../figures/guide/gui_outputs.png)

The facts on the cards are read from the output files; they state what the output contains
(ranges, fractions, events, stability ratios, effective numbers of members) and never interpret
the system. Colours mark the quality of a number consistently on all pages (Table 14).

| Quantity | green | amber | red |
| --- | --- | --- | --- |
| largest relative energy error $\lvert\Delta E/E\rvert$ of a run (not marked for dissipative runs) | $< 10^{-6}$ | $< 10^{-3}$ | $\geq 10^{-3}$ |
| stability ratio (Mardling–Aarseth, Vynatheya et al. 2022) | $\geq 1.5$ | 1–1.5 | $< 1$ |
| effective number of members $n_\mathrm{eff}$ of an ensemble | $\geq 100$ | 20–100 | $< 20$ |
| reduced $\chi^2$, MCMC length, priors of a fit | as stated on the fit page | | |

**Table 14.** The colour code of the result pages. Each row is a quantity, the columns the
ranges shown in green (good), amber (look at it) and red (do not trust the number without
checking). The ranges are printed under every *Data quality* box.

Every result page has the same frame: a header with the name of the configuration, the kind of
result, its state, a line on what was computed and the job's note; the main facts in boxes; a
*Data quality* box on the right; a blue box *What am I looking at?*; the parts made for this
kind of result; *Other results of this configuration*; *Compute next*; *For a paper*; and
*Advanced*. *Compute next* offers what usually follows (a light curve or the O−C of a run, more
members of an ensemble, the fitted configuration of a fit) and opens the Configuration page with
the item ticked; greyed-out actions say why they are not available. *For a paper* explains that
every figure can be saved as a vector graphic (the camera icon above it) and its data as CSV (the
disk icon), and *Methods paragraph* writes a paragraph that describes how the result was computed,
with the references it needs. *Advanced* holds the files of the job, its log, its configuration
and job file, free plots of any columns, all time series of a run and the 3D viewer.

### The result page of a run {#gui-run}

![**Figure 9.** Top of the result page of a run (V907 Sco, 200 yr N-body with relativity, tides
and rotation). Header: `v907_sco_published`, *run*, *finished*, "3 stars · orbits out, inner ·
nbody · GR, tides, rotation · 0 → 200 yr · 4001 outputs · 39 s". Fact boxes: *Stopping event*
none, *Hierarchy changes* 0, *Eclipsing (fraction of run)* inner 50.5 %, *Mutual inclination
out–inner* 25.6–26.3°, *Lowest stability ratio* 1.67 (MA01), *Max |ΔE/E|*
$1.3\times10^{-11}$. *Data quality*: energy conservation $1.3\times10^{-11}$, angular momentum
$4.1\times10^{-7}$, stars evolve no, backend nbody. *The model system*: the summary boxes
(hierarchy, dynamical stability, eclipses of inner at $t = 0$); *Hierarchy* – the tree with the
orbits as boxes and the stars as circles coloured by temperature; *On the sky* – the orbits of
both components of every orbit projected on the sky (x, y in au, the node $\Omega$ measured from x towards y; blue: out, orange:
inner); *From above* – the same seen from the pole of the sky plane, the observer at the top; *Pair
inner on the sky, stars to scale* – the two stars of the eclipsing pair drawn with their radii
along their apparent orbit, with the inclination, the eclipse limit, $a$ and the impact
parameter. Below: tables of the orbits ($P$, $e$, $i$) and of the nested triple ($a_\mathrm{in}/a_\mathrm{out}$,
mutual inclination, MA01 and V22 ratios).](../../figures/guide/gui_run_top.png)

The top of a run page (Figure 9) answers "what is this system and did the computation work".
The fact boxes say whether the run ended early, whether the hierarchy changed, which binaries
eclipse and for how much of the run, the range of the mutual inclinations, the lowest stability
ratio and the energy error. *The model system* draws the configuration at $t = 0$ four ways. In
the sky view, a pair eclipses when its two stars overlap on the sky at conjunction; in the
zoomed panel the stars are drawn to scale, so the eclipse geometry can be read directly (here
the inclination 76.6° is just above the limit 72.85°, and the impact parameter is 0.785 of
$R_1+R_2$: grazing, partial eclipses).

![**Figure 10.** The middle of the run page. *System view* (dark panel): the stars at the time
of the slider, as seen from Earth (*as seen from Earth*) or face-on to the outer orbit; on the
right a close-up of the pair `inner` with the stars to scale and the label *out of eclipse* or
*eclipse now*; *Play*, the time slider and the playback speed (*inner-orbit speed* or
*outer-orbit speed*). *Hierarchy*: the configured hierarchy and whether it changed during the
run. *Stars*: mass, radius, $T_\mathrm{eff}$ and the largest Roche-lobe filling of every star at
the start and the end. *Interactive viewer*: opens or saves the 3D viewer of the run as one HTML
file. *Eclipses over the run*: the inclination of the binary `inner` to the sky (blue curve,
degrees) against time (years), its eclipse limit (dashed line, 72.85°) and the time windows in
which it can eclipse (blue bars at the top); here 50.5 % of the run in four
windows.](../../figures/guide/gui_run_system.png)

The *System view* (Figure 10, top) animates the run: positions come from the stored osculating
elements, with Kepler motion between outputs. *Eclipses over the run* is the central
observational result: the inclination of the inner binary oscillates between about 53° and 89°
with the nodal precession period (63 yr measured in this run), and the binary eclipses only while
the inclination is above its eclipse limit $\arccos[(R_1+R_2)/a]$. The bars at the top are these
windows; a gap is a time when the orbit has turned away and the eclipses are absent.

![**Figure 11.** The bottom of the run page. *Orbits over the run*: tabs *Eccentricity*,
*Inclination to the sky*, *Mutual inclinations*, *Semi-major axis* and *Energy error*; one small
panel per orbit, each with its own vertical range (min–max printed at the top right); here the
eccentricity of `inner` (0.0041–0.0131) and of `out` (0.4663–0.4712) against time in years.
*Secular dynamics: theory vs. the run*: for every inner orbit the mutual inclination at $t = 0$,
the eccentricity range of the quadrupole theory and of the run, the precession period estimated
and measured, the apsidal period measured and expected for the binary alone, and the MA01 and V22
ratios. *Other results of this configuration* (eclipse times, stability check, light curve) and
*Compute next* (Light curve, Eclipse times exact or as an observer, Stability check; *Ensemble
around it* and *Pulsation periods* greyed out with the reason; *Change and run
again*).](../../figures/guide/gui_run_orbits.png)

*Orbits over the run* (Figure 11) shows the osculating elements, sampled at the outputs; the fast
oscillations of the outer eccentricity (orange) are the response to the inner binary on every
outer revolution, not a secular change. The table *Secular dynamics* compares the run with the
orbit-averaged quadrupole theory: agreement says that the system behaves as the analytic theory
predicts; large differences point to octupole, resonant or non-hierarchical effects (*Physics
and Methods*, chapter Dynamics).

### Light curves {#gui-lightcurve}

![**Figure 12.** The light-curve page (V907 Sco, TESS band, 30 d at 10-min cadence). *Light
curve*: the relative flux (1 = all stars outside eclipse) against time in days after $t = 0$,
with *show:* buttons for the total light, the binary *inner (m1 + m2)* and each star, and the
choice *relative flux* or *magnitude*; the bar under the plot selects a time range. *Folded on
the period of a binary*: the light of `inner` folded on its period (3.77487 d, the period in the configuration; phase 0 = its
deepest point), with the choice *light of this binary* or *total light*. *Light of each star and
binary*: the share of the light outside eclipse and the deepest loss of its own light (inner
98.3 % and 6.18 %; m1 46.5 % and 12.1 %; m2 52.1 % and 10.0 %; m3 1.69 % and 0.0413 %). *Compare
with a measured light curve*: a CSV file (time in BJD or days, flux or magnitude) drawn over the
model, with the BJD of model $t = 0$.](../../figures/guide/gui_lightcurve.png)

The light-curve page (Figure 12) separates the light into its parts, so that each feature can be
traced to the star that causes it. The folded curve shows the primary and secondary eclipses of
the binary on its own period; scatter around the mean shape comes from the other stars and from
changes of the orbit during the window. Measured data are drawn over the model after giving the
BJD of model $t = 0$, because the model has no calendar: its time starts at $t = 0$, the epoch
of the elements.

### Eclipse times and minima {#gui-timings}

![**Figure 13.** The eclipse-times page (exact times from the dynamics, V907 Sco, 100 yr before
and after $t = 0$). *O − C of every eclipsing binary*: the O−C in minutes against a linear
ephemeris fitted to the primary eclipses, against time in years from $t = 0$; filled dots:
primary eclipses (eclipses of m1), open circles: secondary eclipses, shifted by their mean offset
(−8.44 min, given in the legend). The eclipses exist only in four windows; within each window
the O−C runs downward over about 200 min. *Per binary*: the numbers of primary and secondary
eclipses (8953, 8924), the period of the ephemeris (3.776279318 d), the O−C range (−104.2 to 103.8
min), its half range and the impact parameter. *Compare with measured minima*: a CSV file of
measured minima and the BJD of model $t = 0$. *Other results* and *Compute
next*.](../../figures/guide/gui_timings.png)

Figure 13 shows how the dynamics appears in an O−C diagram. A linear ephemeris fitted to all
eclipses cannot follow the slow drift of the eclipse period inside each window, which is caused
by the precession of the orbit with respect to the line of sight; within a window, primary and
secondary eclipses run almost parallel, so the apsidal motion of the nearly circular binary is
small. For a binary with a light-time orbit, the O−C instead shows a wave with the outer period.
The page of minima measured *as an observer* has the same layout and adds the error of every
minimum, the shape of the eclipses (depth, width and kurtosis of the Mikulášek profile, per
season) and the shift between the measured minimum and the geometric conjunction.

### Stability checks and other results {#gui-other}

![**Figure 14.** The page of a stability check, an example of the pages of the commands without
a dedicated layout. Fact boxes: *inner in out: stability ratio* MA01 1.67 · V22 2.2 ("with m3,
i_mut 26.2° · above 1: stable") and *N-body ghost test* stable ("100 outer orbits, largest
divergence $9.60\times10^{-4}$"). *Report*: the text the job printed – the nested triple with
$a_\mathrm{in}/a_\mathrm{out}$, the eccentricities, the mutual inclination, the outer mass ratio
and both ratios, and the result of the N-body test. *Other results of this configuration*,
*Compute next* (*Open the configuration*, *Tools*), *For a paper* and
*Advanced*.](../../figures/guide/gui_stability.png)

Commands without a page of their own (stability checks, chaos maps, scattering experiments,
pulsation periods, emulators) show the text they printed and the facts recognised in it
(Figure 14), with their files under *Advanced*.

### The result page of an ensemble {#gui-ensemble}

![**Figure 15.** Top of an ensemble page (QX Cas: 2000 members, forward and backward, 30 yr
N-body each; the mass and orbit of the third star drawn). Header: 3 stars, the backend, the run
length, "7 varied parameters · 4 observations · required: stable_ma01". Fact boxes: *Members*
2 000, *Effective members (n_eff)* 186.7 ("conditioned on 4 observations"), *Most common fate*
intact 100 %, *Stable (Mardling–Aarseth)* 100 %, *Eclipsing (time average)* inner 94 %. *Data
quality*: $n_\mathrm{eff}$ with the observations (green), backend, members with an error
(0).](../../figures/guide/gui_ensemble_top.png)

An ensemble page (Figures 15–18) answers "which systems are possible, and which do the
measurements allow". The model system at the top (not shown) is the nominal configuration.

![**Figure 16.** The conditioning and the parameters of the ensemble. *Observations used for
conditioning*: one row per entry of `ensemble.observed` with a tick box *use*, the quantity
(`inner_i0`, `inner_didt0_deg_yr`, `fate` backward, `fate` forward), the direction (forward,
central, backward) and the constraint (82 ± 1.5, −0.23 ± 0.08, = intact); the effective number
of members with the ticked observations (186.7) and *Re-weight*. *Parameters: drawn vs.
conditioned*: one histogram per drawn parameter and for the initial mutual inclination, grey =
drawn (prior), blue = weighted by the observations; above each, the median and 16–84 % range
before → after; the panels are sorted by how much the observations narrowed the range (the outer
period first, its range ×0.46; logarithmic axis for a log-uniform draw).](../../figures/guide/gui_ensemble_params.png)

The observations table (Figure 16) is interactive: untick a measurement and press *Re-weight* to
see what it contributes; nothing is integrated again. The weight of a member is
$w = \exp(-\chi^2/2)$ for Gaussian measurements (1 for an exact match, 0.61 at 1σ, 0.14 at 2σ)
and 0 when it violates a limit; the effective number of members
$n_\mathrm{eff} = (\sum w)^2/\sum w^2$ (Kish 1965) says how many members carry the weight. Below
about 20 the conditioned histograms are noise and are hidden. Where the blue histogram is
narrower than the grey one, the measurements constrain the parameter: here they prefer outer
periods of 190–1030 d over the drawn 120–4400 d, and leave the angles almost unconstrained.

![**Figure 17.** Outcomes of the ensemble, weighted by the observations. Tabs *Forward
(future)* and *Backward (past)*. *Fates*: intact 100 % (99.5–100 %, 68 % Wilson interval).
*Stability*: the fractions stable by Mardling–Aarseth and Vynatheya et al. (2022). *Orbits*: per
orbit the largest eccentricity (median and 16–84 %), the fraction with $1-e < 0.1$, the present
rate of change of the inclination $di/dt$ at $t = 0$, the time-averaged eclipsing fraction and the
fraction of members whose eclipses end within the run (19.5 %, median 23 yr). *Mutual
inclinations*: at $t = 0$ (84.7°, 44–135°) and the fraction of members that cross 90°.
Histograms (grey drawn, blue weighted): largest $e$ of inner and of out, eclipsing fraction of
inner, lowest MA01 stability ratio.](../../figures/guide/gui_ensemble_outcomes.png)

The outcomes (Figure 17) are the answer of the ensemble: fractions with their 68 % Wilson
intervals (Wilson 1927), which behave well also for fractions close to 0 or 1, and medians with
16–84 % ranges, all weighted. The *Backward* tab gives the same for the runs into the past.

![**Figure 18.** *Members with the largest weight*. The grey box explains the table: each row is
one drawn system; *weight* says how well it reproduces the ticked observations; *share* is its
part of the total weight; the stability columns are the MA01 and V22 ratios (green ≥ 1.5, amber
1–1.5); "59 members have a weight above 0.5: they reproduce the observations equally well". The
columns: member number, weight, $\chi^2$, share, the drawn values ($m_3$, outer $P$, $e$, $i$,
$\Omega$, $\omega$, $M_0$), the fate forward/backward, the two stability ratios and *open as
configuration*.](../../figures/guide/gui_ensemble_members.png)

*Open as configuration* (Figure 18) turns a member into an ordinary draft with its drawn values,
ready to be run, with a light curve or an O−C diagram. The answer of an ensemble is the range of
the parameters, not one member: members with similar weights are equally good examples.

### The result page of a fit {#gui-fit}

![**Figure 19.** Top of the page of a stage-2 fit (V907 Sco: 58 minima, 172 radial velocities,
2562 light-curve points, 14 inclinations, 3 priors, 16 free parameters, 4384 model evaluations).
Fact boxes: *χ² start → fitted* 11 920 → 1052 ("reduced 0.3765 – below 0.5: the errors used are
larger than the scatter, the parameter errors are conservative"), *priors > 2σ off* 1
(`stars.m2.R` −4.3σ), *inner: rms of the minima* 2.6 min (start 3.5 min), *radial velocities:
rms* 3.85 km/s, *light curves: rms* 0.666 ppt (2 segments), *inclinations: χ²* 222 (14 epochs).
*Data quality*: reduced χ² and two warnings (the error factors still changed after 4 fits; a
strongly correlated pair of parameters). *The model system*: as for a run, with the start
configuration dashed and the fitted one solid.](../../figures/guide/gui_fit_top.png)

The fit page (Figures 19–23) shows the model fitted to the measurements and everything needed to
judge it. Read the warnings first: they say what needs attention.

![**Figure 20.** Warnings, parameters and the O−C of the fit. *Check these first*: the
warnings of the fit with their **?**. *Parameters*: one row per free parameter with its start
value, fitted value, linear error from the covariance matrix, unit, prior (mean ± σ) and the
deviation of the fitted value from the prior, and the shift (fitted − start)/error; deviations
from a prior above 2σ are red. *Eclipse timing (O − C) with the model*: the measured minima
(dots: primary filled, secondary open, with error bars) and every eclipse of the fitted model
(solid: primaries, dashed: secondaries), as O−C in minutes against one linear ephemeris over the
years 2021–2025.](../../figures/guide/gui_fit_params.png)

In the O−C panel of Figure 20, apsidal motion moves the primary and secondary curves apart, while
the light-time effect and the perturbation by the third star move both together; here the
secondary curve lies about 15–25 min below the primary one and both carry the 142-d wave of the
outer orbit.

![**Figure 21.** Residuals of the minima and the radial velocities. *Times of minima*: measured
minus fitted model in minutes against time, filled primary, open secondary, crosses for minima
left out as outliers; the buttons switch between the fitted and the start model; the error bars
are the errors the fit used (with the error floor of the source). Tables: the rms per source
(TESS sectors) and the mean residual of primary and secondary minima with their difference
(0.435 min). *Radial velocities*: measured points and the model at the same times for m1 and m2,
against the phase of the binary (tabs *phase of its own binary*, *against time*, *phase of
inner*, *phase of out*); blue CHIRON 2022, orange Lacy (1999); the grey ticks are the model,
which includes the motion around the third star and is therefore not one sine curve; below, the
residuals in km/s.](../../figures/guide/gui_fit_minima_rv.png)

The difference between the mean residuals of primary and secondary minima (Figure 21) is a
sensitive check: a difference of several σ means that the model puts one kind of eclipse at the
wrong phase, usually a wrong $e\cos\omega$ of the pair.

![**Figure 22.** Light curves and inclinations of the fit. *Light curves*: one tab per segment
(TESS S39, S66), *fitted* or *+ start model*; the measured flux (dots) and the calibrated model
(line, scale, zero point and trend solved per segment) in windows around the eclipses, against
BJD − 2459361; below, the residuals in parts per thousand (ppt); the table gives the points and
the rms per segment. *Inclinations against time*: the inclination of the pair in the N-body
model (blue: fitted, orange: start) and the measured epochs (black dots with error bars) from
1900 to 2022; the table lists every epoch with the measured value, its error, the start and
fitted model values, the eclipse limit (for one-sided data), the normalised residual and the
source (DASCH, APPLAUSE, TESS, FRAM).](../../figures/guide/gui_fit_lc_inc.png)

The inclination panel (Figure 22) shows a property that only a dynamical model can fit: the
eclipses of V907 Sco disappeared and reappeared during the 20th century because the inclination
of its orbit changes, and the N-body model of the whole system follows the measured epochs.

![**Figure 23.** The end of the fit page. *Predicted minima* (above, cut): the next minima in
BJD_TDB and as a calendar date, with their uncertainty σ in minutes, the eclipse duration and
the impact parameter. *Other fits of this system*: every fit folder next to this one with its
reduced χ², number of free parameters, rms of the minima and data used; *Compare the parameters
with* adds a second fit to the parameter table. *Report*: the text report of the fit.
*Compute next*: *Open the fitted configuration* (run, light curve, eclipse times of the
solution), *Change the fit and run again*, *HTML report* (greyed out here) and *Could a residual
wave be magnetic activity?* (the Applegate test). *For a paper* and
*Advanced*.](../../figures/guide/gui_fit_bottom.png)

The predictions (Figure 23) are computed from many parameter sets drawn from the MCMC samples or
the covariance matrix, one N-body integration each, so their errors include the correlations of
the parameters. *Other fits of this system* makes variants (with and without light curves,
priors or radial velocities) comparable at a glance.

## Tools {#gui-tools}

![**Figure 24.** The Tools page with the Applegate calculator selected. Left: the tools grouped
as *Calculators* (Applegate mechanism, three-body break-up theory), *Checks* (validation against
published results), *Results* (interactive viewers of the runs) and *Stellar models* (own models
with MESA); each with the command it runs. Right: the help of the tool (what it does, when to use
it, what it needs, what you get, references), its form (mass, radius, $T_\mathrm{eff}$, $a$,
$P_\mathrm{bin}$, $P_\mathrm{mod}$ of the active star and binary; the observed O−C semi-amplitude),
the job name and *Run now*, *Add to the inbox*, *Add and run*. Bottom: *Where the former commands
are* – for every command the page and the item where it is found.](../../figures/guide/gui_tools.png)

The Tools page (Figure 24) holds what belongs neither to a configuration nor to a result. Each
tool has a form generated from the options of its command (chapter
[The command line](#command-line)). *Run now* runs the command at once and shows what it prints
(at most 5 minutes); longer computations go to the inbox. The table at the bottom says where each
command is found in the application (Table 15).

| Command | Where in the application |
| --- | --- |
| `run`, `lightcurve`, `timings`, `minima`, `stability`, `chaosmap`, `ensemble`, `pulsation`, `scatter`, `fit` | Configuration › What to compute |
| `view` | Outputs › a run › Interactive viewer; all runs of a folder: Tools › Interactive viewers |
| `ensemble-report` | Outputs › an ensemble (text report under Advanced) |
| `emulator-train`, `sbi` | Outputs › an ensemble › Compute next |
| `emulator-report` | Outputs › an emulator › Compute next |
| `validate`, `applegate`, `escape-theory`, `mesa-grid`, `mesa-reduce` | Tools |
| `rv`, `apsidal`, `worker`, `gui` | command line only (the worker is started on the Jobs page) |

**Table 15.** Where every command of SWINGERsim is found in the web application. Left: the
sub-command; right: the page and the place on it.

## Folders used by the application {#gui-folders}

| Folder | Content | Written by |
| --- | --- | --- |
| `configs/` | configurations that come with SWINGERsim and your own annotated files | you |
| `jobs/drafts/` | configurations written in the application, members of ensembles opened as configurations | Configuration, Outputs |
| `jobs/submitted/` | frozen copies of the drafts as they were sent | Configuration, Jobs |
| `jobs/inbox/`, `running/`, `done/`, `failed/`, `cancel/` | the job queue (Table 9) | the application and the worker |
| `output/<job>/` | the results of a job sent from the application | the command of the job |
| `output/.gui_cache/` | cached cards, pages and viewers; rebuilt when a job or its output changes | the application |

**Table 16.** The folders the web application reads and writes, all inside the SWINGERsim folder
given by `--root`.
