# The SQUADRA systems {#squadra-systems}

SQUADRA is the survey of 2+2 quadruple candidates of J. Kolář, P. Zasche and collaborators: systems in
which TESS sees two eclipsing pairs in one star image (Kolář et al. 2025: 781 candidates). The question
for SWINGERsim is whether the observed behaviour – eclipse depths that change, eclipse times that
deviate from a constant period – can be explained by a dynamical model, and which hierarchy is the
right one. This chapter does not reproduce a published model, because none exists for these systems;
it tests the ensemble machinery of SWINGERsim on real targets: conditioning on observations, the
evidence of competing models, the agreement of the secular and N-body backends, and the symmetry of the
dynamics in time. The detailed analyses of the SQUADRA systems with real minima and light curves are
the subject of Parts II–IV.

**Definitions.** In an ensemble, each member gets a **weight** – the likelihood of the observations
given its parameters, for example $\exp(-\chi^2/2)$ for a measured rate of change with its error. The
mean weight over all members drawn from the prior is the **evidence** $Z$ of the model; the ratio of
the evidences of two models (the Bayes factor) says which explains the data better. The **effective
number of members** $n_\mathrm{eff} = (\sum w)^2/\sum w^2$ (Kish 1965) counts how many members actually
carry the weight; an $n_\mathrm{eff}$ near 1 means that one member carries everything and the evidence is
not reliable. The **fate** of a member says how its run ended: intact, collision, ejection,
reconfiguration (a change of hierarchy), exchange or breakup.

## CzeV4315: 2+2 or 2+(2+1)? {#czev4315}

**The system.** CzeV4315 shows two eclipsing pairs: A (6.74 d, eccentric, $e$ = 0.64), whose eclipses
have been stable for years, and B (0.92 d), whose inclination changes by −1.4 ± 0.2 °/yr – its eclipses
are getting shallower quickly. An eclipse-timing variation with a period of about 119 d is seen. Is it a
plain 2+2 quadruple, or a 2+(2+1) quintuple in which B has its own close companion on a 119-d orbit?

**Input.** A 2+2 ensemble of 20 000 members with the secular backend (±300 yr; outer period
log-uniform between 40 and 3000 d, isotropic orientations, estimated masses), and a 2+(2+1) ensemble of
10 000 members with the N-body backend and relativity (±60 yr; companion of B on an orbit of
119 ± 3 d). The observations used as conditions: the tilt rate of B, no Roche-lobe overflow, and
optionally A eclipsing for at least 5 yr in the past and at least 0.25 yr ahead.

**Table 29.** Evidence of the two models of CzeV4315. Column 3 lists the observations that condition
the ensemble; columns 4–5 the evidence and the effective number of members.

| Model | Members | Conditions | Evidence $Z$ | $n_\mathrm{eff}$ |
|---|---|---|---|---|
| 2+2, broad prior of the outer period | 20 000 | tilt rate of B | 0.0063 | 177 |
| 2+2, broad prior of the outer period | 20 000 | + A eclipsing for years | 0.00013 | 4.8 |
| 2+(2+1) | 10 000 | tilt of B + A eclipsing + detached | 0.0104 | 165 |
| 2+(2+1) | 10 000 | the same + intact in both directions of time | 0.0095 | 153 |

**Result.** Table 29 shows that the 2+(2+1) quintuple is favoured over the 2+2 by a factor of about 70
under the same conditions. The observation that discriminates most is that A has eclipsed for years:
it removes almost all 2+2 members. The 2+2 members that fit lie at the edge of stability (weighted
median MA01 ratio 1.43, outer period about 116 d). Here the secular backend, which assumes a stable
hierarchy, cannot be trusted alone. An N-body check of the 40 best 2+2 members over 100 yr with MEGNO
finds 4 collisions, 14 chaotic members (MEGNO 4–35, Lyapunov times about 1.5–25 yr) and 22 regular ones: about
half of the 2+2 solutions cannot describe a system that exists today, so the 2+2 evidence is an
overestimate and the preference for 2+(2+1) is conservative.

![**Figure 9.** The best 2+2 configuration of CzeV4315 in the secular ensemble: the impact parameter of both pairs over 600 years. Horizontal axis: time from 300 yr in the past to 300 yr ahead (dashed line: today). Vertical axis: the smallest projected separation of the two stars of a pair at conjunction, in units of the sum of their radii; below 1 (shaded band) the pair eclipses. Blue: pair A; orange: pair B.](../../figures/showcase/czev4315.png)

Figure 9 shows why a 2+2 has difficulty. Pair A (blue) stays far below 1 – it eclipses all the time.
Pair B (orange) swings up and down as its orbital plane precesses; it crosses 1 every few decades, so
its eclipses vanish for years and return (B eclipses 80 % of the time, in windows of about 20 yr).
Today B is on its way down – the observed fast change of its eclipse depth. A 2+2 reproduces
this tilt rate only with outer periods of 50–170 d (median 116 d), close to the stability limit; in the 2+(2+1) model the 119-d companion of B drives the precession of B directly.

**Time symmetry.** Conservative dynamics is symmetric in time: a configuration that is unstable will
show it in the past as often as in the future. In 20 000 runs of the 2+(2+1) ensemble (10 000 members
in both directions), 5.99 % of the members fail forward and 5.97 % backward. 3.93 % fail in both
directions within the window, where independent failures would give 0.36 %: instability is a property
of the configuration, and the starting phase only decides when it shows. The cumulative failure
fraction grows without saturating (0.4 % at 20 yr, 2.4 % at 40 yr, 6.0 % at 60 yr), as expected near the
stability boundary, where the disruption time grows roughly exponentially with the hierarchy (Hayashi
et al. 2022; Mushkin & Katz 2020). The 4.2 % of runs that end in a collision do so by three routes: ZLK
cycles of B driven by its companion (mutual inclinations 60°–120°; B reaches $e\approx0.45$, where its
stars touch), an eccentricity of the companion pumped to about 0.9 by pair A on the outer orbit, so that
its periastron comes within 1–3 semi-major axes of B and the inner triple becomes unstable, and – rarely –
a growth of the eccentricity of A from 0.64 to about 0.86.

## CzeV4315 with tides and rotation {#czev4315-tides}

**Question.** The 2+(2+1) ensemble used point masses with relativity. Real stars are distorted by tides
and rotation; the distortion adds a fast apsidal precession and dissipates energy. With the apsidal
constants of these stars (Claret 2023) the tidal and rotational precession of B has a period of about
9 yr, shorter than the ZLK periods of 14–118 yr. Does this change which members survive?

**Input.** The same prior and random seed as before, so that members 0–4999 are exactly the same draws,
plus tides and rotation; N-body, 100 yr forward and backward (10 000 runs).

**Table 30.** CzeV4315 with tides and rotation: fates of the 5000 members after 100 yr in each direction.

| Fate (100 yr) | Forward | Backward |
|---|---|---|
| intact | 94.68 % | 94.34 % |
| collision | 2.50 % | 2.52 % |
| ejection | 1.26 % | 1.38 % |
| reconfiguration | 0.92 % | 1.12 % |
| exchange | 0.36 % | 0.38 % |
| breakup | 0.28 % | 0.26 % |
| not intact in both directions | 4.20 % (0.30 % if independent) | |

![**Figure 10.** CzeV4315 with tides and rotation. *Left:* for every fate other than intact, the percentage of the 5000 members that ended so within 100 yr; blue bars: runs into the future, orange bars: runs into the past. *Middle:* the percentage of the same members that did not stay intact within the first 60 yr, with point masses and relativity (dark grey) and with tides and rotation added (green), forward and backward. *Right:* histogram (both axes logarithmic) of the largest eccentricity that the 0.92-d pair B reached in each run, with point masses (grey filled) and with tides and rotation (green outline).](../../figures/part1/job010_fates.png)

**Result.** Figure 10 and Table 30 show three things. In the left panel the blue and orange bars have
almost the same length for every fate (collisions 2.50 % forward, 2.52 % backward): the time symmetry
holds with tides as well. In the middle panel the failure rate within the first 60 yr drops from 6.2 %
to 3.3–3.6 % when tides and rotation are added; 223 runs that collided with point masses stay intact.
The right panel shows why: with point masses the eccentricity of B is pumped to anything up to 1
(median 0.45 in the runs that collided), whereas with tides almost all runs stay near $10^{-3}$ (median
0.0015) – the tidal and rotational bulges make the orbit of B precess so fast that its companion can no
longer pump its eccentricity. This suppression of ZLK cycles by short-range forces is known (Fabrycky &
Tremaine 2007; Liu, Muñoz & Lai 2015); here it appears in a real candidate. Rerun with tides over 100 yr,
examples of all three collision routes survive, while a deep plunge of the companion (periastron at
1.6 $a_B$) still collides: short-range forces remove the secular collisions and part of the dynamical
ones. The practical conclusion: point masses overestimate the failures of compact 2+(2+1) systems by a
factor of two, and the ensembles of such systems need tides and rotation. CzeV4315 is analysed with its
real minima in Part IV.

## ASASSN-V J075203 {#j075203}

**The system.** In this 2+2 candidate, B tilts at 1.03 ± 0.2 °/yr, A keeps a constant inclination
(0 ± 0.3 °/yr) and A shows a fast apsidal motion (apsidal period $U\approx3.5\pm1.5$ yr). Can a 2+2
reproduce all three at once?

**Table 31.** ASASSN-V J075203: the four steps of the analysis. Column 3 gives what each step found.

| Step | Members | Result |
|---|---|---|
| isotropic prior, secular | 20 000 | no member fits all three observations (best $\chi^2$ 6.2; $n_\mathrm{eff}$ = 1.0) |
| outer orbit near the plane of A ($i\sim N(80.7°, 3°)$, $P_\mathrm{out}$ 40–600 d), secular ±100 yr | 40 000 | $n_\mathrm{eff}$ = 53; 32 members with $\chi^2$ < 3; $P_\mathrm{out}$ 50–80 d, A almost coplanar with the outer orbit |
| N-body ±30 yr of the 32 | 32 | 2 collide; 7 fit all three ($\chi^2$ < 3.5, $P_\mathrm{out}$ 55–83 d) |
| N-body with MEGNO ±1000 yr of the 7 | 7 | 4 regular, 1 chaotic (Lyapunov time about 13.5 yr), 2 collide within 50–150 yr |

**Table 32.** The seven members that fit all three observations in N-body. Columns: member number, outer
period, MA01 stability ratio (above 1 = stable by the criterion), MEGNO (2 = regular) and verdict.

| Member | $P_\mathrm{out}$ [d] | MA01 ratio | MEGNO | Verdict |
|---|---|---|---|---|
| 642 | 55.3 | 1.06 | 2 | regular |
| 12744 | 62.3 | 1.16 | 2 | regular |
| 39269 | 66.3 | 1.22 | 2 | regular |
| 38234 | 77.5 | 1.21 | 2 | regular |
| 1043 | – | 1.08 | 40–55 | chaotic, intact for 1000 yr |
| 2654, 30610 | – | 0.87–0.89 | – | collide, excluded |

**Result.** Tables 31 and 32 trace the search. With an isotropic prior one member carries the whole
weight ($n_\mathrm{eff}$ = 1.0): the evidence of such an ensemble is meaningless and the prior must be
focused, which the second step does. The members that fit lie at the stability limit (MA01 0.9–1.2,
period ratio $P_\mathrm{out}/P_A\approx$ 6–9). There, the orbit averaging of the secular backend fails:
the N-body apsidal period of A is about half the secular one, as expected when the period ratio is small
(Luo, Katz & Dong 2016). Only 4 regular members survive all tests, with outer periods of 55–78 d. A 2+2 is
therefore possible only in this narrow regular island at the stability limit, and it makes two
predictions that existing minima can test: the eccentricity of A oscillates by 0.13–0.17 within years,
and the minima of A show an eclipse-timing variation with the 55–78-d outer period. If these are absent,
the 2+2 is excluded. The lesson for the software is general: near the stability limit the secular result
is a pre-selection, and the N-body backend gives the answer.

## Eclipse statistics of 2+2 populations {#populations}

**Question.** How often do both pairs of a 2+2 eclipse, and how fast do eclipses appear and disappear?
This calibrates what SQUADRA can expect to see.

**Input.** 50 000 synthetic 2+2 systems with the secular backend over 3000 yr. Assumed prior: masses
uniform in 0.8–3 M☉, inner periods log-uniform in 0.3–30 d, inner $e$ uniform in 0–0.2, outer periods
log-uniform from 100 d to 3000 yr, outer $e$ uniform in 0–0.8, isotropic orientations, stable by the
criterion of Vynatheya et al. (2022).

**Table 33.** Eclipse statistics of the synthetic 2+2 population. Rows: fraction of pairs eclipsing
today, and the fractions whose eclipses stop (among those eclipsing today) or start (among those not
eclipsing today) within 10, 100 and 1000 yr. Columns: pair A and pair B.

| Quantity | A | B |
|---|---|---|
| eclipsing now | 31.7 % | 32.3 % |
| eclipses stop within 10 / 100 / 1000 yr | 3.4 / 10.2 / 19.0 % | 3.0 / 10.2 / 19.4 % |
| eclipses start within 10 / 100 / 1000 yr | 2.6 / 12.9 / 24.7 % | 2.6 / 13.0 / 24.3 % |

![**Figure 11.** Synthetic 2+2 population (50 000 systems): the fraction of eclipsing pairs whose eclipses stop (left) and of non-eclipsing pairs whose eclipses begin (right) within a given time from today (horizontal axis, logarithmic, in years). Blue: pair A; orange: pair B.](../../figures/showcase/population.png)

![**Figure 12.** The same population sorted by the outer period (horizontal axis, logarithmic). Each point is the fraction of systems in that period bin; the shaded band is its statistical uncertainty. Red: pair A tilts faster than 0.1°/yr. Blue: pair A starts or stops eclipsing within 100 yr. Green: both pairs eclipse at some time within 3000 yr; the dotted line is the fraction expected at a random moment for independent pairs.](../../figures/population_eclipses.png)

**Result (Table 33).** In Figure 11, read for example at 100 years: about 10 % of the eclipsing pairs stop
eclipsing and about 13 % of the hidden ones start. The curves of A and B lie on top of each other, as
they must for pairs drawn from the same prior, and P(B eclipses | A eclipses) equals P(B eclipses) – a
check that the eclipse statistics treat isotropic orbits correctly (cf. Vokrouhlický 2016, Sect. 3.3).
Both pairs eclipse together 10 % of the time on average. Figure 12 shows that fast changes need short
outer periods: among eclipsing pairs a tilt faster than 0.1 °/yr occurs in 92 % of the systems with
outer periods below 0.6 yr, in 27 % at about 2 yr, 5 % at 4 yr and never beyond 10 yr. Rates above
1 °/yr – the rates measured for CzeV4315 and ASASSN-V J075203 – need outer periods of 105–460 d, which
agrees with the observed fast-tilting systems (outer periods of 0.3–3 yr; Kolář et al. 2025; V994 Her;
CzeV4315). The probability that a pair eclipses, 32 % for this prior, agrees with the geometric
estimate $(R_1+R_2)/a$ averaged over the prior. The population is assumed, not fitted; a comparison
with the 781 candidates needs the selection function of TESS, which is outside the scope of SWINGERsim.

## CzeV343 {#czev343}

The fourth SQUADRA system of Part I is CzeV343, a 2+2 quadruple in a period resonance studied
dynamically by Pejcha et al. (2022). Its apsidal-motion contributions are four of SWINGERsim's
benchmarks ([Chapter 3](#tides)): the tidal, rotational and relativistic rates of pair A agree with the
paper to 0.6–2 %, and the full N-body quadruple gives the observed apsidal period to 3.3 %. CzeV343 is
analysed in full in Part IV.
