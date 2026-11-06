# Stability and fates {#stability-and-fates}

A hierarchical multiple star is *dynamically stable* when its orbits stay well separated for a
very long time; an unstable one breaks up, exchanges partners or ejects a star, usually within
a number of outer orbits that is small compared with its age. Whether an observed configuration
is stable is therefore a constraint on its unknown parameters: a solution that disrupts in a few
thousand years cannot describe a system that exists today. SWINGERsim offers tools from cheap to
definitive – algebraic criteria, machine-learning classifiers, a direct N-body test, chaos maps –
and classifies the fates of the runs.

## Algebraic criteria for nested triples {#criteria}

Each criterion applies to a *nested triple*: an inner orbit and the orbit around it, the rest of
the system reduced to point masses at the centres of mass. A 2+2 quadruple has two (A with B as a
point mass, B with A); a 3+1 also two (`inner` with m3, `mid` with m4); larger trees one per orbit
pair. Both criteria are written as a ratio that exceeds 1 for stable systems.

**Mardling & Aarseth (2001)** (MA01), in the form of Vynatheya et al. (2022, eq. 2):

$$ \frac{a_\mathrm{out}(1-e_\mathrm{out})}{a_\mathrm{in}} >
2.8\left[\frac{(1+q_\mathrm{out})(1+e_\mathrm{out})}{(1-e_\mathrm{out})^{1/2}}\right]^{2/5}
\left(1 - 0.3\,\frac{i_\mathrm{mut}}{180°}\right), \qquad
q_\mathrm{out} = \frac{m_\mathrm{out}}{m_\mathrm{in}} , $$

where $a_\mathrm{in}$, $a_\mathrm{out}$ are the semi-major axes of the inner and outer orbit,
$e_\mathrm{out}$ the outer eccentricity, $m_\mathrm{in}$ the mass of the inner orbit,
$m_\mathrm{out}$ that of the outer component and $i_\mathrm{mut}$ the mutual inclination. The left
side is the outer periastron distance in units of the inner semi-major axis.

**Vynatheya et al. (2022)** (V22, their eqs. 3–4) add the eccentricity the inner orbit reaches in
ZLK cycles and a better inclination dependence:

$$\widetilde Y_\mathrm{crit} = 2.4\left[\frac{1+q_\mathrm{out}}{(1+\tilde e_\mathrm{in})(1-e_\mathrm{out})^{1/2}}\right]^{2/5}\left[\frac{1-0.2\tilde e_\mathrm{in}+e_\mathrm{out}}{8}(\cos i_\mathrm{mut}-1)+1\right],$$

with $\widetilde Y = a_\mathrm{out}(1-e_\mathrm{out})/[a_\mathrm{in}(1+\tilde e_\mathrm{in})]$ (outer
periastron over inner apastron), $\tilde e_\mathrm{in} = \max(e_\mathrm{in}, \tfrac12 e_\mathrm{max}^2)$
and $e_\mathrm{max}$ the quadrupole ZLK maximum ([Hierarchies and orbits](#zlk-cycles)). The ratio
$\widetilde Y/\widetilde Y_\mathrm{crit}$ is reported.

![**Figure 10.** The Mardling & Aarseth (2001) stability limit of a coplanar triple. Horizontal
axis: eccentricity of the outer orbit; vertical axis: ratio of the outer to the inner semi-major
axis $a_\mathrm{out}/a_\mathrm{in}$ (logarithmic). The black curve is the limit for three equal masses ($q_\mathrm{out} = 0.5$);
above it (green) the triple is stable for very long times, below it (red) it breaks up or
exchanges stars. At $e_\mathrm{out} = 0$ the limit is $a_\mathrm{out}/a_\mathrm{in} \approx 3.3$;
it rises steeply with the outer eccentricity because the outer periastron, not the semi-major
axis, decides.](../../figures/concepts/stability.png)

Figure 10 shows the content of the criterion: what matters is how close the outer component comes
at periastron, so an eccentric outer orbit must be several times wider than a circular one. For
inclined orbits the factor $(1 - 0.3\,i_\mathrm{mut}/180°)$ lowers the limit slightly (retrograde
orbits are more stable), while the V22 form raises it near the ZLK window, where the inner
eccentricity grows.

**Limits.** Both criteria are fits to triple simulations. For quadruples, reducing a binary to a
point mass is not always reliable, especially in 3+1 systems, where the fourth star excites the
eccentricity of the middle orbit (Vynatheya et al. 2023). The run summary prints both ratios for
every orbit pair as a first check.

## Machine-learning classifiers {#mlp}

Vynatheya et al. (2023) trained multi-layer-perceptron classifiers on $5\times10^5$ N-body
integrations each of 2+2 and 3+1 quadruples; they classify 94 % and 93 % of the systems correctly,
against 88 % and 66 % for their triple MLP (Vynatheya et al. 2022) applied to the nested triples
and 83 % and 56 % for the criterion of Mardling & Aarseth (2001) (Vynatheya et al. 2023,
Tables 3–4). SWINGERsim uses the networks published by the
authors when their repository is installed. Hierarchies of more than four stars are classified
sub-quadruple by sub-quadruple: every orbit with two orbits as components gives a 2+2, every orbit
inside an orbit inside an orbit a 3+1, deeper sub-systems being point masses at their centres of
mass (a 2+2+2 gives four sub-quadruples). The system counts as stable if all are; outside the
training set this is an approximation. SWINGERsim can also train its own classifiers on ensembles
([Ensembles and inference](#emulators)).

## The N-body ghost test {#ghost-test}

The definitive answer for one configuration is a direct integration with the definition of
stability of Vynatheya et al. (2023, Sect. 2). The system and a *ghost* copy, in which the
semi-major axis of the lighter inner binary (2+2) or of the innermost binary (3+1) is larger by
$10^{-6}$ au, are integrated for 100 outer orbits as Newtonian point masses. The system is
unstable if an orbit becomes unbound (for at least three consecutive checks) or if the relative
difference of the two semi-major axes exceeds $10^{-2}$ – the divergence of neighbouring
trajectories that signals chaos. The test takes seconds to a minute.

**Limits.** All three tools refer to dynamical stability over about 100 outer orbits. A system that
passes can still evolve secularly (ZLK cycles, chaotic diffusion) on longer time scales; MEGNO and
long runs check that. The ghost criterion is sensitive to the short-period oscillations of the
osculating semi-major axis: for moderately hierarchical systems the two copies, which drift in
phase, differ by about 1 % even when both stay bound for 1000 outer orbits. "Diverged" verdicts with
a difference only slightly above $10^{-2}$ are therefore marginal.

## Chaos maps {#chaos-maps}

`swingersim chaosmap` integrates a grid of configurations in which one or two parameters vary
(for example `orbits.out.a` and `orbits.out.e`) and records MEGNO ([Dynamics](#megno)) at the end
of each run. A point counts as chaotic when MEGNO ≥ 4; a point where a star escapes (beyond 20
times the largest initial semi-major axis) or two stars collide is marked separately. The runs are
Newtonian point masses because REBOUND's variational equations cover only Newtonian gravity.

![**Figure 11.** MEGNO chaos map of a test particle around an equal-mass circular binary of
separation 1 au, after 3000 yr. Horizontal axis: semi-major axis of the particle's orbit
(`orbits.out.a`, au, i.e. in units of the binary separation); vertical axis: its eccentricity
(`orbits.out.e`). Colour: MEGNO, from blue (2, regular) through white (4) to dark red (10, used
for escapes). The blue region at large distances and small eccentricities is regular; the red
region close to the binary and at high eccentricity is chaotic or unstable; the boundary rises
with eccentricity, as the stability criteria predict.](../../figures/chaosmap_circumbinary.png)

Figure 11 reads like Figure 10 with the axes exchanged: stable orbits lie far from the binary and
on nearly circular orbits, and the boundary moves outwards as the eccentricity grows. The isolated
red cells inside the blue region and blue cells inside the red one are chaotic islands and
resonances, which a single criterion cannot represent. With the initial phase (`orbits.<orbit>.M0`)
as the second parameter, the smallest distance beyond which every phase is regular is a stability
limit: for a test particle around an equal-mass circular binary SWINGERsim finds 2.30 binary
separations, against 2.27 in the published grid of Quarles et al. (2018) (see *Testing, Part I*).

## Fates and their classification {#fates}

Every run, and every member of an ensemble, ends with a *fate*. Table 9 lists them.

**Table 9.** Fates of a run. *Fate* is the value stored in the outputs (column `fate` of
ensembles); *Meaning* what happened; *Detected by* the mechanism that recognises it.

| Fate | Meaning | Detected by |
|---|---|---|
| `intact` | the hierarchy survived the whole run | no event |
| `exchange` | a star swapped partner | hierarchy detection |
| `reconfiguration` | the tree changed shape without losing a star | hierarchy detection |
| `ejection` | a star left the system | hierarchy detection |
| `breakup` | the system split into unbound parts | hierarchy detection |
| `capture` | a new bound structure formed | hierarchy detection |
| `hierarchy_change` | another confirmed change of the tree | hierarchy detection |
| `collision` | two stellar surfaces touched | REBOUND collision search |
| `merger` | two stars merged (common envelope, contact binary) | mass-transfer history |
| `roche_overflow` | a star filled its Roche lobe and transfer was not followed | Roche filling factor |
| `e_reached` | a requested eccentricity was reached | `run.stop_when_e_reached` |
| `energy_error`, `wall_time`, `gw_unresolved`, `error` | the run stopped for a numerical or computing reason | run control |

A change of the hierarchy is confirmed only after it has persisted for
`run.hierarchy_persistence` outputs ([Dynamics](#hierarchy-detection)). The last row is not
physics: such members should be counted separately and, if frequent, the settings changed.

## Survival in both directions of time {#survival}

A system observed today has existed since its formation, so a configuration that disrupts in the
recent past cannot describe it. Because the dynamics without dissipation is time-reversible, the
same holds for the near future: an unstable configuration is unstable in both directions. In a
10 000-member N-body ensemble of CzeV4315 integrated ±60 yr, 6.0 % of the members disrupt forward
and as many backward, and two thirds of those that fail in one direction fail within the same 60
years in the other; the rest only because the window is short. A one-sided survival over a short
run is therefore no evidence of stability. Survival is imposed in ensembles as a condition on both
runs, `backward:fate: {equals: intact}` and `forward:fate: {equals: intact}`
([Ensembles and inference](#weights)).

A run of decades only catches the fastest instabilities. Near the stability boundary the time to
disruption grows roughly exponentially with $(r_\mathrm{p,out}/a_\mathrm{in})^{3/2}$
($r_\mathrm{p,out}$ the outer periastron distance) and scatters by one to two orders of magnitude
between neighbouring initial phases, a signature of chaos (Hayashi et al. 2022, 2023; Mushkin &
Katz 2020); for strongly inclined orbits (60°–150°) ZLK oscillations shorten it further, and
retrograde coplanar systems are the most stable. A configuration consistent with an observed system
must survive its age – usually millions of outer orbits – so the criteria and MEGNO are the
complementary tests for members close to the boundary.

**What changes the stability.** Without dissipation and with constant masses the stability of a
hierarchy does not change by itself; slow processes move it across the boundary (Table 10).

**Table 10.** Processes that move a hierarchy across the stability boundary. *Destabilising* and
*Stabilising* list the processes with the literature; the rows are not paired.

| Destabilising | Stabilising |
|---|---|
| mass loss of the inner binary (winds, supernova): the hierarchy weakens (Perets & Kratter 2012; Toonen et al. 2022) | tidal shrinking of an inner binary, also through ZLK cycles with tidal friction (Fabrycky & Tremaine 2007) |
| expansion of a star: Roche-lobe overflow, mass transfer, common envelope | precession by tides, rotation and relativity that quenches ZLK cycles |
| ZLK cycles driving a high eccentricity | dissipation in gas during formation (damped eccentricities, aligned orbits) |
| encounters and the Galactic tide on wide orbits (Kaib et al. 2013; Grishin & Perets 2022) | accretion and migration in the formation disc (Tokovinin & Moe 2020) |
| sweeping through resonances when periods change | retrograde coplanar configurations (Hayashi et al. 2022) |

All processes of the table except supernovae are modelled in SWINGERsim; gas damping is
parametric.

## Scattering experiments {#scattering}

`swingersim scatter` studies dynamical formation: a binary (or a single star) meets another binary
(or a single star). Each encounter is set up as in Hut & Bahcall (1983) and FEWBODY (Fregeau et al.
2004): isotropic orientations and random phases of the binaries, fixed or thermal eccentricities
($f(e) = 2e$), the centres of mass on a hyperbola with velocity at infinity $v_\infty$ and impact
parameter $b$ ($b^2$ uniform up to $b_\mathrm{max}$; with a maximum periastron $q$,
$b_\mathrm{max}^2 = q^2 + 2GMq/v_\infty^2$ from gravitational focusing). The critical velocity $v_c$
is the one at which the total energy is zero, $\tfrac12\mu_\infty v_c^2 = \sum G m_1 m_2/(2a)$
($\mu_\infty$ the reduced mass of the two incoming objects, the sum over the binaries). The
integration (IAS15, point masses, optional 1PN, collisions when radii are given) starts at a
multiple of the largest semi-major axis and ends when the top-level objects found by the hierarchy
detection are pairwise unbound, receding and farther apart than the start distance, with every
bound sub-system of three or more stars stable by the MA01 criterion. A system that stays bound
counts as captured once it has been stable for 10 of its outer orbits; otherwise the run stops at
its time or computing limit as *unresolved*.

The outcomes follow Antognini & Thompson (2016): *flyby*, *exchange*, *single ionization*,
*exchange + ionization*, *triple + single*, *full ionization*, *bound 2+2* or *bound 3+1* (a
captured quadruple) and *collision*. Cross-sections are $\sigma = \pi b_\mathrm{max}^2 \times$ the
fraction of the outcome. With tidal capture, every close periastron passage of two stars removes
the energy of the tides it raises (Press & Teukolsky 1977),
$E = \sum G M_\mathrm{o}^2/R\,[(R/R_p)^6 T_2(\eta) + (R/R_p)^8 T_3(\eta)]$ with
$\eta = (M/(M+M_\mathrm{o}))^{1/2}(R_p/R)^{3/2}$, where $M$, $R$ are the mass and radius of the
perturbed star, $M_\mathrm{o}$ the mass of the other, $R_p$ the periastron distance and $T_2$,
$T_3$ the polytrope fits of Portegies Zwart & Meinen (1993). The relative velocity is reduced at
periastron, so unbound stars can be captured into binaries: two 0.6 M☉ stars ($n = 3$, 0.65 R☉)
meeting at 10 km s⁻¹ are captured within 2.37 stellar radii, $n = 1.5$ stars within 3.33 radii,
close to the about 3.4 of McMillan et al. (1987). Gas drag during the encounter is not modelled.

## Statistical escape theory {#escape-theory}

When three stars interact resonantly (no hierarchy, a "democratic" phase) the outcome is
practically random, and Stone & Leigh (2019) showed that its distribution is the phase-space
volume allowed by the conserved energy and angular momentum (their eq. 2, with the strong-
interaction radius $R = 2a_B(1+e_B)$, $a_B$ and $e_B$ the elements of the final binary).
`swingersim escape-theory` evaluates this integral numerically for each of the three possible
escapers and gives the escape probabilities and the distributions of the binary's energy,
eccentricity and orientation, in a second instead of thousands of N-body runs. It reproduces the
limits of the paper (energy distribution $\propto |E_B|^{-7/2}$ for zero angular momentum and
$|E_B|^{-4}$ near the maximum; their eccentricity distributions) and agrees with SWINGERsim's own
scatterings: for resonant encounters of a 1 + 0.6 M☉ binary with a 0.8 M☉ star the lightest star
escapes in 75 % of the N-body runs and in 77 % in the theory. The theory covers only resonant
three-body encounters; prompt exchanges and fly-bys lie outside it.

## Implementation {#implementation-stability}

**External code and data.** The MLP classifiers are the trained networks of Vynatheya et al.
(2023), used as published; the N-body integrations of the ghost test, chaos maps and scattering
use REBOUND/IAS15 and its MEGNO (Rein & Liu 2012; Rein & Spiegel 2015). **Own code.** The MA01
and V22 criteria for every nested triple of any tree; the decomposition of large hierarchies into
sub-quadruples and the computation of the classifiers' input parameters; the ghost test; the chaos
maps and their classification; the fate classification; the scattering set-up, termination and
outcome classification; the tidal-capture energy loss with the Portegies Zwart & Meinen (1993)
fits; and the numerical evaluation of the Stone & Leigh (2019) phase-space integral.
