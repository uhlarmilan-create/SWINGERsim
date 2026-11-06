# Dynamics {#dynamics}

SWINGERsim computes the motion of the stars in two complementary ways, called *backends*. The
**N-body backend** integrates the equations of motion of all stars directly; it is exact up to
the included physics, sees every orbit and every close encounter, and is slow when the inner
orbits are short and the run is long. The **secular backends** average every orbit over its
revolution and follow only the slow changes of the averaged ellipses; they are thousands of
times faster and reach billions of years, but they assume that the hierarchy stays intact. The
backend is chosen with `run.backend` (`nbody` or `secular`). This chapter describes both, the
diagnostics computed during a run, and the weak external perturbations (gas, Galactic tide,
passing stars, gravitational waves).

## Direct N-body integration {#nbody}

All stars (2 to about 10) are integrated as point masses – extended bodies when tides are on
([Tides and rotation](#tides-and-rotation)) – with the IAS15 integrator of REBOUND (Rein & Liu
2012; Rein & Spiegel 2015). IAS15 is a 15th-order Gauss–Radau integrator with an adaptive step.
Its error stays at the level of floating-point round-off, about $10^{-15}$ per step, and it
handles close encounters and high eccentricities without special settings. The price is speed:
the step follows the shortest orbit, so the computing time grows with
$t_\mathrm{max}/P_\mathrm{in}$, the run length divided by the shortest period.

Because the full N-body problem is solved, all orders of the interaction (quadrupole,
octupole and higher), the short-period terms that orbit averaging removes, mean-motion
resonances, the coupling between the two binaries of a 2+2 system and any break-down of the
hierarchy are included automatically.

## Relativity {#relativity}

With `physics.gr: true` (the default) SWINGERsim adds the first post-Newtonian (1PN) equations of
motion for all bodies with the `gr_full` effect of REBOUNDx (Tamayo et al. 2020). The simpler
REBOUNDx effects `gr` and `gr_potential` assume one dominant central mass, which is not true in
a multiple star, and are not used. The main consequence is the periastron advance
$\dot\omega_\mathrm{GR}$ of each orbit ([Hierarchies and orbits](#apsidal-motion)); the
simulated advance of an isolated binary (1.5 M☉, $a$ = 0.05 au, $e$ = 0.3) agrees with the
analytic formula to 1 % (see *Testing, Part I*).

Relativity is a small correction for most stellar binaries, with apsidal periods of thousands
to tens of thousands of years at periods of days. It becomes decisive for binaries with periods
of about a day or shorter and a slowly acting distant companion, where it can suppress ZLK
cycles completely ([Hierarchies and orbits](#zlk-cycles)).

## Conservation checks {#energy-error}

The relative energy error $|E/E_0 - 1|$ is recorded at every output. Without relativity it stays
near $10^{-14}$. With relativity the checked quantity is the conserved 1PN Hamiltonian of all
bodies, which already contains the Newtonian kinetic and potential energy; it stays at
$10^{-13}$–$10^{-11}$, also at high eccentricities and with tides. With tides the spin and
tidal energies are included. The total angular momentum is checked in its Newtonian form, which
the 1PN equations do not conserve exactly: with relativity it fluctuates at about $10^{-7}$
without growing. A run stops with the event `energy_error` when the error exceeds
`run.max_energy_error` (default $10^{-5}$), far above both levels. When energy is lost
physically – tidal dissipation, mass loss, gravitational waves – the energy is reported but no
longer used as a stop criterion; with tidal dissipation the conservation of the total angular
momentum of orbits and spins is checked instead (`angmom_error`).

## Collisions {#collisions}

Every star has a radius, from the stellar model or the configuration. With `run.collisions`
(default on) REBOUND tests at every step whether the surfaces of two stars touch; the run then
stops with the event `collision`, the time and the stars involved – or, with `physics.mergers`,
continues with the merged star ([Mass transfer and evolution](#mergers)). The default radii of
the `zams` model are those of zero-age main-sequence stars; for any system in which collisions
matter, evolved radii should be used.

## Hierarchy detection {#hierarchy-detection}

At every output SWINGERsim determines which stars are bound to which. Starting with every star
as a separate object, it repeatedly joins the pair whose two-body relative orbit is bound and
has the smallest semi-major axis, replaces the pair by its centre of mass, and stops when no
bound pair remains. A stable 2+2 system gives one tree, `((m1, m2), (m3, m4))`. A new structure
must persist for `run.hierarchy_persistence` consecutive outputs (default 3) before it is
confirmed, so a single unusual snapshot near a close passage does not count. Confirmed changes
are classified as exchange, reconfiguration, ejection, break-up or capture
([Stability and fates](#fates)). The detection is purely energetic (bound or not); it does not
judge long-term stability, which the subsequent dynamics shows.

## Diagnostics of orbit pairs {#pair-diagnostics}

For every orbit that has a parent in the tree (A and B in `out` for a 2+2; `inner` in `mid` and
`mid` in `out` for a 3+1) SWINGERsim compares analytic expectations with the run. Table 5 lists
them.

**Table 5.** Diagnostics of an orbit pair (inner orbit and the orbit around it). *Quantity* is
what is reported; *Basis* the formula or theory; *Meaning* how to read it. All analytic values
are computed from the initial elements and are estimates for orientation; the simulated values
are the result.

| Quantity | Basis | Meaning |
|---|---|---|
| ZLK window, $e_\mathrm{max}$, $t_\mathrm{ZLK}$ | quadrupole test-particle theory ([Hierarchies and orbits](#zlk-cycles)) | whether and how strongly ZLK cycles are expected |
| range of $e$, $i_\mathrm{mut}$, cycle period | double-averaged quadrupole equations, any mass ratio (Naoz et al. 2013) | prediction of the secular cycle |
| $t_\mathrm{prec}$ | quadrupole nodal precession ([Hierarchies and orbits](#nodal-precession)) | expected precession period |
| GR apsidal period, ZLK suppression flag | 1PN periastron advance | relativity suppresses ZLK if faster |
| $\varepsilon_\mathrm{oct}$ | Naoz (2016) | strength of the octupole term |
| $\varepsilon_\mathrm{SA}$ | Luo, Katz & Dong (2016) | validity of orbit averaging; warning above 0.01 |
| Mardling–Aarseth and Vynatheya et al. ratios | stability criteria ([Stability and fates](#criteria)) | > 1 stable |
| $\beta$ | Hamers & Lai (2017) | coupling of precession and ZLK |
| $J_{AB}^\mathrm{max}$, period ratio | Vokrouhlický (2016); nearest small-integer ratio | 2+2 plane coupling, resonances |
| measured $e$ range, $i_\mathrm{mut}$ range, number and period of ZLK cycles, measured precession period | the run | the simulated result |

### Secular prediction at quadrupole order {#secular-prediction}

For every orbit pair SWINGERsim integrates the orbit-averaged quadrupole equations of Naoz et
al. (2013, their eqs. A27, A29 and A36), with the perturber treated as a point mass. Unlike the
test-particle formula, these equations hold for any mass ratio, conserve the total angular
momentum exactly and let the outer orbit change its inclination. They predict the range of $e$
and $i_\mathrm{mut}$ and the period of the eccentricity cycle. The implementation reproduces the
published curve of Naoz et al. (2013, Fig. 4) to 0.04° and the Kozai limit to $5\times10^{-5}$
(see *Testing, Part I*).

### Extremes between outputs {#extremes}

Eccentricity peaks of ZLK cycles can be much shorter than the spacing of the outputs; in an
eccentric ZLK test $1-e$ drops to about $10^{-5}$ for a few centuries. The N-body backend
therefore samples the state internally every `run.extremes_every_orbits` periods of the shortest
inner orbit (default 10) and keeps for every orbit the largest eccentricity and the smallest
periastron distance $q = a(1-e)$ with their times, and for every binary of two stars the largest
Roche-lobe filling factor at periastron. A spacing of a few orbits suffices because these peaks
are secular and last many inner orbits; collisions are caught separately at every integrator
step. The extra cost is a few per cent.

### Chaos indicator MEGNO {#megno}

With `run.megno: true` SWINGERsim computes the Mean Exponential Growth factor of Nearby Orbits
(MEGNO; Cincotta & Simó 2000) from REBOUND's variational equations. MEGNO tends to 2 for
regular (quasi-periodic) motion and grows linearly with time for chaotic motion; its slope gives
the Lyapunov time, the e-folding time of the separation of neighbouring trajectories. Chaos means
that a single run cannot be predicted beyond a few Lyapunov times and only statistics over many
runs are meaningful ([Ensembles and inference](#ensembles-and-inference)). REBOUND's variational
equations cover only Newtonian gravity, so MEGNO is computed in a parallel Newtonian point-mass
copy of the run, without relativity, tides or mass loss (MEGNO is refused together with mass
loss); the cost roughly doubles.

## Backward integration {#backward}

A negative `run.t_max_yr` integrates into the past, with both backends. Gravity and the
relativistic terms are time-reversible: integrating a run forward and its end state back
returns the initial elements to $10^{-12}$. Three limits apply to the past:

- **Chaos.** Beyond a few Lyapunov times the uncertainties of the observed elements grow so
  much that a single past has no meaning; then only an ensemble of runs within the
  uncertainties gives the distribution of possible pasts.
- **Tidal dissipation** run backward makes orbits wider and more eccentric; the result depends
  directly on the poorly known dissipation rate.
- **Stellar evolution** runs backward too: the stars become smaller towards the ZAMS (SSE) or
  the start of the pre-main sequence (MIST), where the run stops (`stellar_birth`).

Survival in both directions is a constraint used in ensembles ([Stability and fates](#survival)).

## Secular backends {#secular}

Orbit averaging replaces each Keplerian orbit by a ring of mass spread along it in proportion
to the time spent there. The interaction energy averaged over the inner and the outer orbit
(*double averaging*) depends only on the shapes and orientations of the orbits, not on the
positions of the stars, so the semi-major axes are constant (in the point-mass problem) and the
equations contain only the slow, secular time scales. The integrator step then follows these
time scales instead of the orbital periods.

### Vectorial elements and equations of motion {#vector-elements}

Each orbit X is described by two vectors, the normalised angular momentum
$\mathbf{K}_X = \sqrt{1-e_X^2}\,\mathbf{k}_X$ and the Laplace vector
$\boldsymbol{\pi}_X = e_X\,\mathbf{p}_X$, where $\mathbf{k}_X$ is the unit orbit normal and
$\mathbf{p}_X$ the unit vector to the periastron (Vokrouhlický 2016). These elements have no
singularities at zero eccentricity or inclination. Their equations of motion (his eqs. 12–17)
are

$$\frac{d\mathbf{K}_X}{dt} = -\frac{1}{\Lambda_X}\left(\boldsymbol{\pi}_X\times\frac{\partial H}{\partial\boldsymbol{\pi}_X} + \mathbf{K}_X\times\frac{\partial H}{\partial\mathbf{K}_X}\right), \qquad
\frac{d\boldsymbol{\pi}_X}{dt} = -\frac{1}{\Lambda_X}\left(\mathbf{K}_X\times\frac{\partial H}{\partial\boldsymbol{\pi}_X} + \boldsymbol{\pi}_X\times\frac{\partial H}{\partial\mathbf{K}_X}\right),$$

with $\Lambda_X = \mu_X n_X a_X^2$ the circular angular momentum of the orbit ($\mu_X$ its
reduced mass, $n_X$ its mean motion, $a_X$ its semi-major axis) and $H$ the averaged interaction
energy of the whole system. The equations are integrated with the eighth-order Runge–Kutta
method DOP853 of SciPy with relative tolerance `run.secular_rtol` (default $10^{-10}$).

### 2+2 quadruples {#secular-2p2}

For a 2+2 system $H$ is the sum of the quadrupole and octupole terms of two triples: binary A
perturbed by the centre of mass of B, and B perturbed by the centre of mass of A
(Vokrouhlický 2016, eqs. 18–25). The octupole terms can be switched off with
`run.secular_octupole: false`. Added to $H$ are terms that depend only on the eccentricity of
each orbit and reproduce the relativistic apsidal advance (with `physics.gr`) and the apsidal
advance from the tidal and rotational distortion of the stars after Sterne (1939) (with
`physics.tides`). Because all these terms are conservative, $H$ stays exactly conserved and
serves as a check of the integration. With tidal dissipation the semi-major axes,
eccentricities and spin rates of A and B evolve in addition by the constant-time-lag theory
([Tides and rotation](#tidal-evolution)). The Fig. 2 run of Vokrouhlický (2016) over 50 kyr
takes about 5 s instead of about an hour with the N-body backend.

### 3+1 and other hierarchies {#secular-nested}

For every other hierarchy (3+1, and trees with five or more stars) SWINGERsim uses the
general nested-binary formulation of Hamers et al. (2015) and Hamers & Portegies Zwart (2016),
written in the same vectors and equations of motion. $H$ is a sum of *binary-pair* terms, one
for each orbit $k$ and each orbit $l$ that contains it:

$$H_{kl} = -\frac{G\,\mu_k M_o a_k^2}{8\,a_l^3}\,f_2(\mathbf{K}_k,\boldsymbol{\pi}_k,\mathbf{K}_l)
- \frac{15}{64}\,\frac{G\,\mu_k M_o a_k^3}{a_l^4}\,\frac{m_{k1}-m_{k2}}{M_k}\,s\,f_3(\mathbf{K}_k,\boldsymbol{\pi}_k,\mathbf{K}_l,\boldsymbol{\pi}_l),$$

where $\mu_k$ is the reduced mass of orbit $k$, $a_k$ and $a_l$ the semi-major axes, $M_o$ the
mass of the component of $l$ that does not contain $k$, $m_{k1}$, $m_{k2}$ the masses of the two
components of $k$ and $M_k$ their sum, $s = \pm 1$ says in which component of $l$ the orbit $k$
lies, and $f_2$, $f_3$ are the quadrupole and octupole functions of the vectorial elements
(Hamers et al. 2015). For a 2+2 quadruple these are exactly the terms of the 2+2 backend, and
both backends give the same result. Not included are the *binary-triplet* terms (octupole order,
$(a_k/a_m)^2(a_m/a_l)^3$ for $k$ inside $m$ inside $l$) and terms of hexadecupole and higher
order; the comparison with the code MSE (Hamers et al. 2021), which has them, shows differences
below 0.001 in $e_\mathrm{max}$ for the tested systems. Relativity, the stellar distortion and
tidal dissipation in every orbit of two single stars are added as for the 2+2 backend; for a 2+2
system both backends agree to the integration tolerance and the circularisation rate equals that
of Hut (1981).

### Output of the secular backends {#secular-output}

The output has the same structure as that of the N-body backend. The elements are mean
elements; the mean anomaly is not followed and is stored as 0. The hierarchy cannot change and
collisions cannot occur. Stellar spins stay aligned with their orbits; with rotating stars this
causes a small error in the angular-momentum check, of the order of $10^{-5}$ in a typical
example.

### When orbit averaging fails {#secular-validity}

Orbit averaging assumes that the outer orbit is much slower than the inner one. The size of the
neglected terms is measured by the single-averaging parameter (Luo, Katz & Dong 2016)

$$ \varepsilon_\mathrm{SA} = \left(\frac{a_\mathrm{in}}{a_\mathrm{out}(1-e_\mathrm{out}^2)}\right)^{3/2}
\frac{m_\mathrm{pert}}{\sqrt{m_\mathrm{in}\,m_\mathrm{tot}}} , $$

where $a_\mathrm{in}$, $a_\mathrm{out}$ are the semi-major axes, $e_\mathrm{out}$ the outer
eccentricity, $m_\mathrm{in}$ the mass of the inner orbit, $m_\mathrm{pert}$ that of the
perturber and $m_\mathrm{tot}$ their sum. Table 6 shows the effect for one binary and perturber
as the outer orbit widens.

**Table 6.** ZLK cycle period of the same binary and perturber computed with the N-body backend
and with the double-averaged theory, for three outer semi-major axes. *$a_\mathrm{out}$* is
given relative to the closest case, which has $\varepsilon_\mathrm{SA}$ = 0.03; *difference* is
(N-body − averaged)/averaged. The difference falls quickly as the hierarchy becomes stronger.

| $a_\mathrm{out}$ | N-body period | averaged period | difference |
|---|---|---|---|
| 1 × | 12.8 yr | 18.9 yr | −32 % |
| 2 × | 153.8 yr | 151.4 yr | +1.5 % |
| 4 × | 1211 yr | 1213 yr | −0.2 % |

At $\varepsilon_\mathrm{SA}$ = 0.03 the eccentricity range still agrees but the cycle period is
a third shorter; at twice the separation ($\varepsilon_\mathrm{SA}$ eight times smaller) the
difference is at the per-cent level. The run summary therefore warns when
$\varepsilon_\mathrm{SA} > 0.01$. Near very high eccentricities the short-period terms neglected
by the averaging can push the osculating eccentricity higher than the mean one (Luo, Katz & Dong
2016); in a reproduction of the population study of Pejcha et al. (2013) the N-body backend found
three times more triples reaching $1-e < 10^{-3}$ than the secular one, while for the quadruples
both agreed.

**Validity in summary.** The secular backends are first-order secular theories. They are
reliable for strong hierarchies ($\varepsilon_\mathrm{SA}$ well below 0.01) far from mean-motion
resonances. They neglect the short-period terms and, in a 2+2, the direct coupling of A and B of
order $\varepsilon^4$; they cannot describe close encounters, exchanges, ejections or
collisions. In chaotic regions they follow the N-body solution only for a while: for the Fig. 4
setup of Vokrouhlický (2016) both agree for about 4 kyr and then diverge, as in the paper.

![**Figure 7.** The same triple computed by both backends. Horizontal axis: time in years;
vertical axis: eccentricity of the inner binary. The thin grey curve is the N-body solution,
which shows every short-period oscillation of the osculating eccentricity; the thick blue curve
is the secular (orbit-averaged) solution. Both show the same ZLK cycles of about 70 yr with a
maximum near 0.77; the secular curve follows the average of the N-body curve.](../../figures/concepts/nbody_vs_secular.png)

Figure 7 shows the relation of the two backends for a well-separated triple. The grey curve
oscillates on the time scale of the outer orbit around the blue one but neither drifts away:
in this regime the secular backend reproduces the slow evolution at a small fraction of the
cost. Table 7 summarises when to use which.

**Table 7.** Choice of backend. *Situation* is the question or the system; *Backend* the
recommended one; *Reason* the physics behind the choice.

| Situation | Backend | Reason |
|---|---|---|
| stability, break-up, exchanges, collisions | N-body | only it follows non-hierarchical motion |
| orbits not well separated ($\varepsilon_\mathrm{SA} > 0.01$, stability ratio near 1) | N-body | averaging fails for mild hierarchies |
| eclipse timings, light curves, radial velocities | N-body | needs the positions on each orbit |
| inner periods of days, runs of Myr or longer | secular | an N-body run would take years of computing |
| evolution over the stellar lifetime (tides, winds, mass transfer) | secular | billions of years |
| large ensembles of stable systems | secular | seconds per member |
| unsure | both | run one system with both and compare |

## External perturbations {#external}

### Gas phase of formation {#gas}

Multiple stars form in gas, and their orbits keep changing as long as the gas is present.
SWINGERsim describes this phase with the parametric damping of Papaloizou & Larwood (2000): for
every orbit listed in `physics.gas.orbits`

$$\frac{da}{dt} = -\frac{f(t)\,a}{\tau_a}, \qquad \frac{de}{dt} = -\frac{f(t)\,e}{\tau_e}, \qquad
\frac{d\mathbf{k}}{dt} = \frac{f(t)}{\tau_i}\left[\mathbf{n} - (\mathbf{n}\cdot\mathbf{k})\,\mathbf{k}\right],$$

where $\tau_a$, $\tau_e$, $\tau_i$ are the migration, circularisation and alignment time scales
(`tau_a_yr`, `tau_e_yr`, `tau_i_yr`), $\mathbf{k}$ the orbit normal, $\mathbf{n}$ the normal of
the gas disc (the initial invariable plane or a given plane), and $f = 1$ until `t_end_yr`,
$f = \exp[-(t - t_\mathrm{end})/t_\mathrm{decay}]$ afterwards. The last equation turns the orbit
plane towards the gas plane without changing $e$, so that $\tan(i/2)$ to the gas plane decays as
$e^{-t/\tau_i}$. The semi-major axes become state variables. The phase runs in the secular
backend. The time scales are inputs: they depend on the mass and structure of the gas, which
SWINGERsim does not model.

### Galactic tide {#galactic-tide}

For orbits wider than about 1000 au the tide of the Galactic disc matters. In the solar
neighbourhood its vertical part dominates, with the potential $2\pi G\rho_0 z_\mathrm{G}^2$
($z_\mathrm{G}$ the height along the Galactic pole, $\rho_0$ the local disc density, default
0.185 M☉ pc⁻³; Heisler & Tremaine 1986; Grishin & Perets 2022, eqs. 3–4). Averaged over an orbit
of reduced mass $\mu$,

$$H = \pi G\rho_0\,\mu\,a^2\left[(1-e^2) - (\mathbf{j}\cdot\hat{\mathbf{z}})^2 + 5(\mathbf{e}\cdot\hat{\mathbf{z}})^2\right],$$

with $\mathbf{j} = \sqrt{1-e^2}\,\mathbf{k}$, the eccentricity vector $\mathbf{e}$ and the
Galactic pole $\hat{\mathbf{z}}$. It drives ZLK-like cycles above the critical inclination
$\cos^2 i < 4/5$ to $e_\mathrm{max} = \sqrt{1 - (5/4)\cos^2 i}$, on the time scale
$\tau_\mathrm{GT} = 0.89\ \mathrm{Gyr}\,(M/\mathrm{M}_\odot)^{1/2}(a/10^4\,\mathrm{au})^{-3/2}$
(Grishin & Perets 2022, eqs. 6–7). In the secular backends these terms are added to $H$ of every
orbit; in the N-body backend the tide is applied as velocity kicks
$\Delta\mathbf{v}_i = -4\pi G\rho_0\,[(\mathbf{r}_i - \mathbf{r}_\mathrm{cm})\cdot\hat{\mathbf{z}}]\,\hat{\mathbf{z}}\,\Delta t$
every 1/400 of the outer period (operator splitting), and its potential energy enters the energy
check. A 60°-inclined wide binary reaches $e_\mathrm{max}$ within $10^{-4}$ of the analytic value.

### Passing stars {#encounters}

A star passing at a periapsis distance $Q$ much larger than an orbit changes the orbit's
eccentricity and orientation but hardly its energy. To first order in the perturbation the
changes of $\mathbf{e}$ and $\mathbf{j}$ follow from Hamers & Samsing (2019, eq. 9), with the
small parameter $\epsilon = [M^2/(m(m+M))\,(a/Q)^3]^{1/2}(1+E)^{-3/2}$ ($M$ the perturber mass,
$m$ the mass of the orbit, $E$ the eccentricity of the perturber's hyperbola). With
`physics.encounters` SWINGERsim draws a Poisson sequence of encounters (rate
$n\pi b_\mathrm{max}^2\langle v\rangle$, relative speeds from a Maxwellian weighted by $v$,
impact parameters uniform in area, isotropic orientations). The secular backends apply the
first-order kicks to every orbit for which the passage is distant ($Q > 3a(1+e)$) and count
closer passages without applying them. The N-body backend integrates every passing star as a
particle from 20 periapsis distances before to 20 after its periapsis, so close passages are
followed too. At $Q = 20a$ the first-order kick agrees with a direct N-body passage to 0.6 %.

### Gravitational waves {#gravitational-waves}

With `physics.gw` every orbit of two single stars loses energy and angular momentum by
gravitational radiation, orbit-averaged after Peters (1964):

$$\dot a = -\frac{64}{5}\frac{G^3 m_1 m_2 M}{c^5 a^3 (1-e^2)^{7/2}}\left(1+\tfrac{73}{24}e^2+\tfrac{37}{96}e^4\right),\qquad
\dot e = -\frac{304}{15}\frac{e\,G^3 m_1 m_2 M}{c^5 a^4 (1-e^2)^{5/2}}\left(1+\tfrac{121}{304}e^2\right),$$

with $m_1$, $m_2$ the masses and $M = m_1 + m_2$. The secular backends integrate these rates
with the orbits; the N-body backend applies them once per output step at fixed orientation and
mean anomaly. The orbital decay of the Hulse–Taylor pulsar is reproduced to 0.02 % (Weisberg &
Huang 2016; see *Testing, Part I*). The effect matters only for compact remnants and very close
binaries.

## Implementation {#implementation-dynamics}

**External libraries.** REBOUND (Rein & Liu 2012) with its IAS15 integrator (Rein & Spiegel
2015), collision search and variational equations for MEGNO; REBOUNDx (Tamayo et al. 2020) with
the effects `gr_full` (1PN for all bodies), `tides_spin` (Lu et al. 2023; next chapter) and the
operator `modify_mass` (mass loss); SciPy's DOP853 integrator for the secular equations; Numba,
when installed, compiles SWINGERsim's 2+2 secular kernel. **Own code.** The secular equations of
motion and the averaged energies for 2+2 (Vokrouhlický 2016, eqs. 12–25) and for any nested tree
(binary-pair terms after Hamers et al. 2015 and Hamers & Portegies Zwart 2016) were coded for
SWINGERsim, together with the eccentricity-only Hamiltonian terms of relativity and stellar
distortion that keep $H$ conserved, the constant-time-lag dissipation in the secular equations,
the gas damping, the Galactic tide (secular terms and N-body kicks), the encounter sequence and
first-order kicks, the orbit-averaged gravitational-wave rates, the hierarchy detection, the
tracking of extremes between outputs, the pair diagnostics, the energy book-keeping (1PN
Hamiltonian, spin and tidal energy, angular momentum) and the backward runs.
