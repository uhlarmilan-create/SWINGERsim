---
title: "SWINGERsim"
subtitle: "Physics and Methods"
---

# Hierarchies and orbits {#hierarchies-and-orbits}

This volume describes the physics and the numerical methods of SWINGERsim, a simulator of
hierarchical multiple stars. For each piece of physics it states what is modelled, the
equations, the assumptions, the numerical method, the configuration keys that control it, how
it was validated and where it stops being valid. The validation itself – the comparison of
every module with published results – is documented in *Testing SWINGERsim, Part I:
Validation against the literature*; here each method carries a one-sentence summary of its
test. How to write a configuration, run the commands and read the outputs is the subject of
the *User Guide* (chapters Configuration, Command line and GUI).

This first chapter defines the geometry: how a multiple star is described as a tree of
orbits, which coordinates and orbital elements SWINGERsim uses, how the observer enters, and
which slow (secular) phenomena change the orbits over time.

## Hierarchical multiple stars {#hierarchical-multiple-stars}

A *multiple star* is a gravitationally bound system of three or more stars. About half of the
Sun-like stars and nearly all massive stars have at least one companion, and many systems
contain three, four or more stars (Tokovinin 2018). Such a system survives for a long time
only if it is *hierarchical*: it can be split into pairs whose orbits have very different
sizes. Two stars close together form a *binary*; seen from a distance much larger than its
separation, the binary acts on other bodies almost like a single star of the combined mass.
It can therefore orbit another star or another binary on a much wider orbit. Each orbit is
then nearly a fixed ellipse (a *Keplerian* orbit), slowly perturbed by the others. Without this
separation of scales the stars exchange partners and eject each other within a few orbits
(chapter [Stability and fates](#stability-and-fates)).

A hierarchy is a **binary tree**: every node is an orbit that joins two components, and each
component is either a star (a leaf) or another orbit. SWINGERsim accepts any such tree up to
about ten stars. The common ones are presets of the key `topology`; any other tree is given
explicitly with `topology: custom` and `hierarchy`, a map from the name of each orbit to its
two components; Table 1 lists the presets.

**Table 1.** Hierarchies available as presets. *Preset* is the value of `topology`; *Stars*
the number of stars; *Orbits* the names under which the orbits appear in the configuration
and in all outputs; *Description* the arrangement. The notation "2+2" means two binaries
orbiting each other, "3+1" a triple (a binary with a close third star) orbited by a distant
fourth star; parentheses mark a sub-system inside a larger one.

| Preset | Stars | Orbits | Description |
|---|---|---|---|
| `custom` | 2–3 | free names, e.g. `A`, `out` | a binary, or a triple ("2+1") |
| `2+2` | 4 | `A`, `B`, `out` | two binaries A = (m1, m2) and B = (m3, m4) on a mutual orbit |
| `3+1` | 4 | `inner`, `mid`, `out` | binary (m1, m2), star m3 around it, star m4 around the triple |
| `2+2+1`, `3+1+1` | 5 | as above plus `out5` | a quadruple with a distant fifth star |
| `2+(2+1)` | 5 | `A`, `B`, `Bt`, `out` | binary A and a triple made of binary B and a close star m5 |
| `2+2+2` | 6 | `A`, `B`, `C`, `AB`, `out` | a 2+2 quadruple AB and a third binary C |
| `3+3` | 6 | `A`, `B`, `T1`, `T2`, `out` | two triples on a mutual orbit |
| `custom` | up to ~10 | free names | any tree of nested binaries |

![**Figure 1.** Hierarchies drawn as trees. Each panel is one hierarchy; its name (binary,
triple 2+1, 2+2, 3+1, 2+2+1, 2+(2+1), 2+2+2, 3+3 and a custom tree of ten stars) is the title.
Orange circles at the bottom are the stars m1, m2, …; green squares are the orbits, labelled
with their names (A, B, out, inner, mid, …); the grey lines connect each orbit with its two
components. The height of an orbit in the tree indicates its level: an orbit higher up is wider
and slower. The figure shows that every hierarchy, however large, is built from the same
element – an orbit of two components – which is why one set of equations serves all of
them.](../../figures/concepts/hierarchies.png)

Figure 1 shows the trees of the presets. Read each panel from the bottom: two stars join into
an orbit (for example A), the orbit joins with another component into a wider orbit (out), and
so on to the top. In the 2+2 panel both components of `out` are binaries; in the 3+1 panel the
tree is a chain, each orbit containing the previous one. The last panel shows that a custom
tree can mix both patterns. Because the model works with the tree and never with a fixed
number of stars, the dynamics, the diagnostics, the observables and the fits are the same for
every hierarchy; the presets only provide names.

Doubly eclipsing 2+2 systems, in which both binaries eclipse, are a particularly valuable class
because the eclipses constrain the orientation of both inner orbits (Zasche et al. 2019, 2023b,
2024; Kolář et al. 2026). Many examples in this volume refer to them, but nothing in the methods
is restricted to them.

## Jacobi (hierarchical) coordinates {#jacobi-coordinates}

Every orbit of the tree describes the **relative motion of its second component around its
first**, each component taken as its centre of mass. These are *Jacobi* (or *hierarchical*)
coordinates. In a 2+2 system the orbit `out` is the motion of the centre of mass of B around
the centre of mass of A; in a 3+1 system the orbit `mid` is the motion of m3 around the centre
of mass of (m1, m2). The gravitational parameter of an orbit is $\mu = G(M_1 + M_2)$, where
$M_1$ and $M_2$ are the total masses of its two components and $G$ the gravitational constant.

**Initial conditions.** The configuration gives Keplerian elements for every orbit. SWINGERsim
converts them to positions and velocities recursively from the outermost orbit inwards: the
outermost orbit places the centres of mass of its two components, then each inner orbit places
its two components around its own centre of mass, down to the stars. Finally the whole system
is shifted so that its barycentre is at rest at the origin. Converting the resulting state back
to elements returns the input to $10^{-12}$ (see *Testing, Part I*).

**Osculating and mean elements.** The elements reported by the N-body backend are *osculating*:
the Keplerian ellipse the two components would follow if all other stars vanished at that
instant. Because the other stars do not vanish, osculating elements oscillate on the time scale
of the orbits even when nothing changes on average. The secular backends
([Dynamics](#dynamics)) work with *mean* (orbit-averaged) elements, which contain no such
oscillations. The difference matters when a published photodynamical solution is used as a
starting point: its osculating elements at an epoch give, in another code, a slightly different
mean period, which accumulates into a drift of the eclipse times
([Fitting](#anchoring)).

## Orbital elements and the reference frame {#orbital-elements}

Each orbit is described by six elements, listed in Table 2.

**Table 2.** The orbital elements of one orbit. *Symbol* is the notation of this volume;
*Key* the name under `orbits.<name>` in the configuration; *Meaning* the definition. The size
is given either as a period or as a semi-major axis, never both.

| Symbol | Key | Meaning |
|---|---|---|
| $P$ or $a$ | `P` [d] or `a` [au] | period, or semi-major axis; linked by Kepler's third law |
| $e$ | `e` | eccentricity, $0 \le e < 1$ |
| $i$ | `i` [deg] | inclination of the orbital plane to the reference plane, 0°–180° |
| $\Omega$ | `Omega` [deg] | longitude of the ascending node in the reference plane |
| $\omega$ | `omega` [deg] | argument of periastron, measured in the orbital plane from the node |
| $M_0$ | `M0` [deg] | mean anomaly at $t = 0$ |

Kepler's third law links the period and the semi-major axis,

$$ a^3 = \frac{G\,(M_1 + M_2)}{4\pi^2}\,P^2 , $$

where $M_1 + M_2$ is the total mass of the two components of the orbit. Two Suns with a period
of one year are $2^{1/3} \approx 1.26$ au apart; with a period of 10 d about 0.11 au, or
24 R☉.

![**Figure 2.** Orbital elements of one orbit. The grey plane is the reference plane (the plane
of the sky when the elements are given relative to the sky). The black vertical axis is its
normal $z$, which points **away from the observer**; the black line in the plane is the
reference direction ($x$ axis) from which the node is measured. The blue ellipse is the orbit,
solid above the reference plane and dashed below it. The dotted line is the line of nodes; the
orange arc $\Omega$ is the longitude of the ascending node, the green line the orbit normal
tilted by the inclination $i$ from $z$, and the violet arc $\omega$ the argument of periastron
from the ascending node to the periastron (violet dot). Read together, $i$ and $\Omega$ fix the
orientation of the orbital plane, $\omega$ the orientation of the ellipse within that plane, and
$a$, $e$ its size and shape.](../../figures/concepts/orbit_elements.png)

Figure 2 shows the geometry. SWINGERsim uses one **reference frame** for all orbits: the
$x$–$y$ plane is the reference plane, $z$ its normal, and the rotation from the orbital frame to
the reference frame is $R = R_z(\Omega)\,R_x(i)\,R_z(\omega)$. When the elements are given
relative to the sky – the normal case for an observed system – the reference plane is the
plane of the sky and **the line of sight is the $z$ axis, pointing away from the observer**.
Consequences of this convention, used throughout the volume:

- the radial velocity of a star is $v_z$, positive when the star recedes;
- light emitted at depth $z$ arrives later by $z/c$, so the observed time of an event at the
  dynamical time $t$ is $t_\mathrm{obs} = t + z/c$, with the barycentre of the system at
  $z = 0$ ([Observables](#light-travel-time));
- $i = 90°$ means the orbit is seen edge-on and eclipses are possible; $i = 0°$ or $180°$ means
  face-on.

**Inclination range 0°–180° and folding to 0°–90°.** In the configuration $i$ lies between 0°
and 180°: an orbit with $i < 90°$ revolves counter-clockwise as seen from the observer's side
of the sky plane, one with $i > 90°$ clockwise. The two senses matter dynamically, because the
mutual inclination of two orbits depends on them (next section). An observer, however, cannot
tell them apart from photometry: an eclipsing binary at $i$ and at $180° - i$ has the same
light curve, since eclipses depend only on $|\cos i|$. SWINGERsim therefore

- keeps the full range 0°–180° in the dynamics;
- folds inclinations to 0°–90° whenever they are compared with photometric measurements
  ($i \to 180° - i$ for $i > 90°$; [Fitting](#inclination-data));
- in ensembles, draws both senses with equal probability when `mirror: true` is set for an
  inclination ([Ensembles](#what-is-varied));
- reports the *tilt* $|i - 90°|$ of an eclipsing binary, whose rate of change has the same sign
  for both senses of revolution ([Ensembles](#rates-today)).

Only the *relative* sense of revolution of two orbits is physical; one orbit (by convention the
first binary) can always be kept at $i < 90°$.

**The invariable plane.** Besides the elements relative to the sky, SWINGERsim stores the
inclination and node of every orbit relative to the *invariable plane*, the plane
perpendicular to the conserved total angular momentum of the system (output columns
`inc_inv`, `Omega_inv`). Secular theories are written in this plane, and the nodal precession
periods in the run summary are measured in it.

## Mutual inclination {#mutual-inclination}

The **mutual inclination** $i_\mathrm{mut}$ of two orbits is the angle between their orbital
planes, computed from their orbital angular-momentum vectors $\mathbf{h}_1$ and
$\mathbf{h}_2$:

$$ \cos i_\mathrm{mut} = \frac{\mathbf{h}_1 \cdot \mathbf{h}_2}{|\mathbf{h}_1|\,|\mathbf{h}_2|} . $$

Unlike $i$, it does not depend on the reference plane, and it is the angle that governs the
dynamics. Values below 90° mean that the orbits revolve in the same sense (prograde), above 90°
in opposite senses (retrograde). Every pair of orbits has its own mutual inclination;
SWINGERsim reports all of them (columns `imut_<pair>`).

For two orbits with sky inclinations $i_1$, $i_2$ and nodes $\Omega_1$, $\Omega_2$,

$$ \cos i_\mathrm{mut} = \cos i_1 \cos i_2 + \sin i_1 \sin i_2 \cos(\Omega_1 - \Omega_2) . $$

Eclipses measure $i_1$ and $i_2$ (up to the folding above) but not the nodes on the sky, so
the mutual inclination of two eclipsing binaries is known only within the range allowed by the
unknown $\Omega_1 - \Omega_2$. This is the main reason why observed systems are studied with
ensembles ([Ensembles and inference](#ensembles-and-inference)).

## Nodal precession {#nodal-precession}

When two orbits are not coplanar, the distant companion exerts a torque on an inner binary.
The torque makes the inner orbital plane **precess**: its angular-momentum vector moves on a
cone around the total angular momentum, keeping $i_\mathrm{mut}$ nearly constant. To lowest
(quadrupole) order, for a circular inner orbit, the precession period is

$$ P_\mathrm{prec} \approx \frac{4}{3}\,\frac{P_\mathrm{out}^2}{P_\mathrm{in}}\,
\frac{m_\mathrm{tot}}{m_\mathrm{pert}}\,\left(1-e_\mathrm{out}^2\right)^{3/2}
\frac{1}{\cos i_\mathrm{mut}} , $$

where $P_\mathrm{in}$ and $P_\mathrm{out}$ are the inner and outer periods, $e_\mathrm{out}$ the
outer eccentricity, $m_\mathrm{tot}$ the total mass of the inner orbit and its perturber, and
$m_\mathrm{pert}$ the mass of the perturber (for a 2+2 system the other binary). The derivation
for triples goes back to Söderhjelm (1975); the 2+2 case is treated by Vokrouhlický (2016).
SWINGERsim prints this estimate for every orbit pair (`t_prec`) next to the period measured in
the run from the node in the invariable plane. The estimate is an orientation; the simulation
is the result.

![**Figure 3.** Nodal precession and its observable consequence (illustrative 2+2 system).
Left: geometry. The green arrow is the total angular momentum; the blue arrows are three
positions of the normal of binary B, which turns on a cone (dashed ellipse) around it.
Right: the inclination of B to the sky (vertical axis, degrees) against time (horizontal axis,
years). The blue curve oscillates with the precession period; the orange band, about 83°–97°,
marks the inclinations at which this binary eclipses. Eclipses are seen only while the curve
crosses the band – here for a few decades twice per precession period – and are absent in
between.](../../figures/concepts/precession_eclipses.png)

Figure 3 shows why precession matters for observers. The cone on the left keeps its opening
angle (the mutual inclination), but the orientation of B relative to the line of sight
changes, so the sky inclination on the right oscillates between $90° - \theta$ and
$90° + \theta$ (for a cone axis in the sky plane; $\theta$ the opening angle). The binary
eclipses only inside the narrow orange band, whose width is set by the eclipse condition of
chapter [Observables](#eclipse-geometry). A slowly precessing binary therefore shows eclipses
that deepen, become shallower and disappear over years to centuries; such changes are observed
in several 2+2 systems (Zasche et al. 2025). In a 2+2 system both binaries precess, generally
at different rates, so their mutual inclination oscillates as well. In larger hierarchies the
same happens on every level, more slowly for wider orbits.

## ZLK cycles {#zlk-cycles}

When the mutual inclination of an inner and an outer orbit lies between about 39.2° and
140.8°, the torque does more than turn the plane: it exchanges inclination for eccentricity.
The inner eccentricity grows to a maximum while $i_\mathrm{mut}$ falls towards 39.2°, then
both return. These are the **von Zeipel–Lidov–Kozai (ZLK) cycles** (Kozai 1962; Lidov 1962;
review by Naoz 2016). The window follows from $\cos^2 i_\mathrm{mut} < 3/5$, i.e.
$39.23° < i_\mathrm{mut} < 140.77°$.

For an initially circular inner orbit, a distant perturber and quadrupole order, the maximum
eccentricity depends only on the initial mutual inclination $i_0$,

$$ e_\mathrm{max} = \sqrt{1 - \tfrac{5}{3}\cos^2 i_0} , $$

and the cycles repeat on a time scale of order

$$ t_\mathrm{ZLK} \sim \frac{P_\mathrm{out}^2}{P_\mathrm{in}}\,
\frac{m_\mathrm{tot}}{m_\mathrm{pert}}\,\left(1-e_\mathrm{out}^2\right)^{3/2} , $$

with the symbols of the previous section. The exact period depends on the initial conditions
(Antognini 2015). SWINGERsim reports for every orbit pair whether it starts in the ZLK window,
the quadrupole $e_\mathrm{max}$, $t_\mathrm{ZLK}$, and, from the double-averaged equations of
Naoz et al. (2013) valid for any mass ratio, the predicted range of $e$ and $i_\mathrm{mut}$
and the period of the cycle ([Dynamics](#secular-prediction)).

![**Figure 4.** ZLK cycles of an inner binary computed with the secular backend. Top:
eccentricity of the inner binary against time in years. Bottom: mutual inclination in degrees
against the same time axis. Blue: initial mutual inclination 70°; grey: 20°. The red dotted
line marks 39.2°, the lower edge of the ZLK window. At 70° the eccentricity climbs from 0.05 to
about 0.86 while the inclination drops to the 39.2° line, and both return – inclination is
traded for eccentricity. At 20°, outside the window, both stay almost constant.](../../figures/concepts/zlk.png)

Figure 4 shows the exchange. The two panels share the time axis, and the eccentricity maxima
of the blue curve coincide with the inclination minima; the inclination minimum touches the
39.2° line, as the quadrupole theory predicts. The maximum eccentricity agrees with
$\sqrt{1 - (5/3)\cos^2 70°} = 0.81$ to within the effect of the small initial eccentricity and
of the finite outer angular momentum. ZLK cycles matter because at maximum eccentricity the
stars approach each other at periastron: tides can then shrink the orbit (Fabrycky & Tremaine
2007), the stars can fill their Roche lobes or collide.

**Octupole order.** When the masses of the inner binary differ and the outer orbit is
eccentric, the next (octupole) term of the expansion breaks the regularity of the cycles and
can drive the eccentricity to values close to 1 and flip the orbit (Naoz 2016). Its strength
is measured by

$$ \varepsilon_\mathrm{oct} = \frac{|m_a - m_b|}{m_a + m_b}\,\frac{a_\mathrm{in}}{a_\mathrm{out}}\,
\frac{e_\mathrm{out}}{1 - e_\mathrm{out}^2} , $$

where $m_a$, $m_b$ are the masses of the two components of the inner orbit and $a_\mathrm{in}$,
$a_\mathrm{out}$ the inner and outer semi-major axes. SWINGERsim reports $\varepsilon_\mathrm{oct}$
for every pair; the N-body backend contains all orders automatically.

**Suppression by short-range forces.** Any additional apsidal precession of the inner orbit
that is faster than the ZLK cycle – relativity, tidal and rotational distortion of the stars, a
close third star – suppresses the cycles, because the periastron turns before the torque can
build up eccentricity (Liu, Muñoz & Lai 2015; Naoz 2016). The run summary flags a pair when the
relativistic apsidal period is shorter than $t_\mathrm{ZLK}$. Tidal and rotational distortion
are often stronger: in the V994 Her setup of Vokrouhlický (2016), binary B reaches $e = 0.63$
with point masses but stays near its initial $e = 0.13$ with the stellar distortion included
(see *Testing, Part I*). A pure point-mass model therefore overestimates eccentricity growth in
close binaries; `physics.tides` should be switched on for them
([Tides and rotation](#tides-and-rotation)).

## Apsidal motion {#apsidal-motion}

The periastron of an eccentric orbit turns slowly in its plane: **apsidal motion**, measured
as the rate $\dot\omega$ or as the apsidal period $U = 2\pi/\dot\omega$. In an eclipsing binary
it is directly observable, because primary and secondary minima move in opposite directions in
the O−C diagram ([Observables](#apsidal-motion-in-o-c)). Three causes contribute:

- **General relativity.** To first post-Newtonian (1PN) order,
  $$ \dot\omega_\mathrm{GR} = \frac{3\,G\,(m_1+m_2)\,n}{c^2\,a\,(1-e^2)} , \qquad n = \frac{2\pi}{P}, $$
  with $m_1$, $m_2$ the masses, $a$ the semi-major axis, $e$ the eccentricity, $n$ the mean
  motion and $c$ the speed of light.
- **Tidal and rotational distortion of the stars** (Sterne 1939), which depends on the
  apsidal-motion constants $k_2$ and the radii of the stars ([Tides and rotation](#classical-apsidal-motion)).
- **The other stars** of the system, whose secular torque turns the periastron as well.

For every binary of two stars SWINGERsim gives the classical expectation (Sterne plus
relativity) and the rate measured in the run from the longitude of periastron
$\varpi = \Omega + \omega$; the difference is the contribution of the other stars. In a compact
quadruple this contribution can be as large as the classical one. The command
`swingersim apsidal` splits the rate into its tidal, rotational, relativistic and dynamical parts
and gives the factor by which $k_2$ would have to change to match a measured rate.

## Resonances {#resonances}

Two kinds of near-commensurability are typical of 2+2 systems:

- **Mean-motion (period) resonance.** When the periods of A and B are close to a ratio of small
  integers, the binaries can be locked in a resonance. The 1:1 case is treated by Breiter &
  Vokrouhlický (2018); the 3:2 resonance is observed in CzeV343 (Pejcha et al. 2022) and its
  statistical excess among 2+2 systems is discussed by Kolář et al. (2025). SWINGERsim reports
  the nearest ratio of small integers of the two periods. Mean-motion resonances are absent from
  the secular backends by construction and are followed only by the N-body backend.
- **Precession (secular) resonance.** When A and B precess at similar rates, their mutual
  torques add up coherently and can drive much larger eccentricity oscillations than in either
  binary alone (Pejcha et al. 2013; Hamers & Lai 2017).

The parameter $\beta$ of Hamers & Lai (2017, their eqs. 32–33) compares the precession rate of
the perturber's orbit with the ZLK rate of an inner orbit. For $\beta \ll 1$ the classical ZLK
picture holds; for $\beta \gg 1$ the inner orbit feels a perturber averaged over its
precession; for $\beta \sim 0.3$–3 the inner eccentricity can be excited strongly and
chaotically, even at initial inclinations well below the ZLK window. SWINGERsim prints $\beta$
for A and B (2+2) or for the inner orbit (3+1).

## Coupling of the orbital planes of A and B {#plane-coupling}

If the outer orbit held all the angular momentum, the planes of A and B would precess
independently around it, and the angle $J_{AB}$ between them would range from $|J_A - J_B|$ to
$J_A + J_B$, where $J_A$ and $J_B$ are their tilts to the outer orbit. In reality the outer
orbit precesses too, and the binaries are coupled through it. For circular A and B at
quadrupole order Vokrouhlický (2016, Sect. 2.3) found two integrals of motion that limit
$J_{AB}$ to a maximum $J_{AB}^\mathrm{max}$, often much smaller than $J_A + J_B$: in his V994 Her
setup $J_A = 11.6°$ and $J_B = 15.8°$ but $J_{AB}^\mathrm{max} = 5.2°$. The planes of two
eclipsing binaries can therefore stay nearly parallel for a long time, which raises the chance
that both eclipse. SWINGERsim prints $J_{AB}^\mathrm{max}$ for every 2+2 core where the bound is
informative (`J_AB_max_analytic_deg`). It assumes circular binaries; eccentric binaries can
exceed it slightly: for the observed eccentricities of V994 Her (Fig. 3 of the paper) Vokrouhlický
finds $J_{AB}$ between 1° and 6°, and SWINGERsim's N-body integration of the same set-up reaches
6.3° (see *Testing, Part I*).

## Larger hierarchies and a fifth star {#larger-hierarchies}

Real quadruples are often parts of quintuples, and an unseen fifth star can change the
dynamics in both directions:

- **Destabilising.** A fifth star on a wide inclined orbit drives ZLK cycles of the orbit of the
  two binaries; its eccentricity grows and its periastron approaches them, which can push a
  quadruple that is stable on its own across the stability limit. The stability report lists
  the nested triple of `out` inside `out5` for this purpose.
- **Stabilising.** Any precession faster than the ZLK cycle suppresses it (Liu, Muñoz & Lai
  2015). A close third star around one binary makes that binary precess quickly and shields it;
  likewise the quadrupole of the two binaries makes their mutual orbit precess and limits a
  ZLK cycle driven by a fifth star. In a test setup with real binaries (1 + 1 M☉, $a$ = 1 au)
  instead of point masses, the eccentricity of `out` reaches 0.73 instead of the point-mass
  value 0.80.
- **Observable consequences.** A companion of only one binary changes only that binary's
  inclination and eclipse timings, as proposed for some doubly eclipsing systems whose binaries
  behave differently (Zasche et al. 2025).

Every method of this volume works with any tree. SWINGERsim remains a code for *hierarchical*
systems: a non-hierarchical (democratic) configuration is followed by the N-body backend until
it breaks up, as in the scattering experiments ([Stability and fates](#scattering)), but it is
not the target of the diagnostics.

## Implementation {#implementation-orbits}

Everything in this chapter is SWINGERsim's own code; no external package is involved apart
from NumPy and SciPy for linear algebra and ODE integration. Own implementations: the
hierarchy tree and the recursive conversion of the elements of all orbits into positions and
velocities and back (Kepler's equation solved by Newton iteration to machine precision); the
mutual inclinations and the invariable plane; the analytic diagnostics of every orbit pair
(quadrupole $e_\mathrm{max}$, $t_\mathrm{ZLK}$, $t_\mathrm{prec}$, $\varepsilon_\mathrm{oct}$,
$\varepsilon_\mathrm{SA}$, the $\beta$ parameter of Hamers & Lai 2017, period ratios) coded
from the cited formulae; the integration of the double-averaged quadrupole equations of Naoz et
al. (2013, eqs. A27, A29, A36), with the SciPy integrator DOP853; and the bound
$J_{AB}^\mathrm{max}$ coded from Vokrouhlický (2016, eqs. 26–36).
