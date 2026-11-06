# Problems found during testing and how they are handled {#problems-found}

Running the benchmarks, the catalogue and ensembles of up to 5·10⁴ members tested the software as much as
the physics. This chapter collects what the tests revealed about the numerics and the practical use of
SWINGERsim, written as properties of the software that a user needs to know: what can go wrong, how
SWINGERsim detects or prevents it, and what the user should do. A second group of findings concerns the
published references themselves; it is listed at the end, because it affects how any code should be
compared with the literature.

Table 45 summarises the findings; each row is detailed in the section named in its last column.

**Table 45.** Summary of the findings. Column 2 says what can happen, column 3 how SWINGERsim handles it,
column 4 where in this report the finding appears.

| Area | What can happen | How SWINGERsim handles it | Seen in |
|---|---|---|---|
| unstable members | a close encounter makes the time step collapse and a run never ends | collision search valid in both directions of time; wall-clock limit checked after every step; fate `wall_time` | [Ch. 6](#czev4315) |
| short evolutionary phases | segments of the evolving backend step over thermal pulses | segment boundaries at the nodes of the tracks, refined to 1 % in radius | [Ch. 8](#giant-tertiary) |
| osculating elements | osculating periods are not the periods of eclipses | mean motion over an orbit for light-curve chunks; anchoring of eclipse periods in fits | [Ch. 5](#qzcar-minima) |
| mass transfer in eccentric orbits | a tertiary pumps the eccentricity of a transferring pair | warning; tides circularise; onset check allows for forced eccentricity | [Ch. 8](#algol-triple) |
| chaotic and ill-conditioned runs | two correct integrations diverge | compare statistics of ensembles, not single members | [Ch. 2](#vokrouhlicky2016) |
| secular backend near the stability limit | orbit averaging fails | stability diagnostics in every run; N-body check of the best members | [Ch. 6](#j075203) |
| ensemble weights | one member carries the whole weight | effective number of members reported with every evidence | [Ch. 6](#czev4315) |
| undefined quantities | apsidal period of a nearly circular orbit | reported as not measurable | [Ch. 7](#catalogue-precession) |
| stellar tracks | a present mass is reached by two tracks; a mass outside the grid; a star overfilling its lobe at the start | warning and the least evolved track; member refused with the reason | [Ch. 5](#qzcar-evolution) |
| tables of atmospheres | temperatures or gravities outside the limb-darkening grids | nearest node with a warning | [Ch. 8](#life-pulsators) |
| conserved quantities | the Newtonian energy is not conserved with relativity | errors measured against the 1PN Hamiltonian | [Ch. 2](#nbody-basics) |

## Unstable members and wall-clock limits {#wall-clock}

An ensemble drawn from a broad prior always contains members close to or beyond the stability limit. In
such a member two stars can fall towards each other. If the collision is detected, the run ends with the
fate *collision*. The tests showed that REBOUND's default ("direct") collision search misses pairs that
approach each other when the integration runs backward in time; two stars then pass through each other,
the adaptive time step of IAS15 shrinks without bound (to 10⁻⁴⁸ yr) and the run does not finish; with the direct
search, eleven of the 10⁴ runs of the CzeV4315 ensemble with tides behave like this. SWINGERsim therefore uses a collision
search that works in both directions of time ("line" search) for backward runs, and checks a wall-clock
limit (`run.max_wall_s`) after *every* integration step, not only between output times. A member that
reaches the limit ends with the fate `wall_time` and is reported as such; it is not silently dropped. Integrated
with the line search, all eleven of these members end in a collision, which is their correct fate. The ghost
tests of [Chapter 9](#stability) have their own limit of 30 minutes per member; none of the 3·10⁴ ghost
tests reached it. What a user should do: keep the limit (it protects the ensemble, not the physics), read
the fates, and treat `wall_time` members as unresolved, not as stable.

## Step size of the evolving backend {#segments}

The evolving backend advances the orbits in segments within which the stellar radii change by less than
1 %. On the thermally pulsing AGB a star's radius jumps by a factor of 2–4 at each pulse for only a few
hundred years, every 10⁴–10⁵ yr. A uniform time grid steps over the pulses: on such a grid the giant
tertiary of [Chapter 8](#giant-tertiary) reaches at most 313 R☉ and 73 % of its Roche lobe and ends at the
end of its track, while the MIST track itself reaches about 500 R☉. SWINGERsim therefore puts the segment
boundaries at the nodes of the MIST track and refines them until no radius changes by more than 1 %
between two boundaries; with this the same giant overflows its lobe at 452 R☉ on a pulse peak. The same
rule protects every short evolutionary phase (the Hertzsprung gap, helium flashes). The price is run time
on the AGB; there is no switch to turn the refinement off.

## Osculating elements and eclipse periods {#osculating}

SWINGERsim's orbital elements are **osculating**: computed from the instantaneous positions and
velocities. In a multiple system, and in a binary of tidally distorted stars, the osculating period is not
the period with which eclipses repeat. Two consequences appeared in the tests. (i) The PHOEBE hybrid
computes the light of distorted stars in chunks of the N-body orbit; if each chunk used the osculating
period of a tidally distorted binary (137 s per orbit off for QZ Car Ac), its minima would drift by tens
of seconds. SWINGERsim uses the mean motion over one orbit, and the PHOEBE minima agree with the spherical
model to a few seconds ([Chapter 5](#qzcar-minima)). (ii) Published elements entered literally into an
N-body model give eclipses that drift by hours within years; Part II shows how SWINGERsim *anchors* the
period and phase of each eclipsing pair to the observed minima. A related point: the light travel time must
be applied to the positions consistently. With exact retarded positions the minima of QZ Car follow the
barycentric conjunctions displaced only by the Rømer delay inside the binary (±16 s); a first-order
correction of the positions would shift them by about 200 s. The Rømer benchmarks of
[Chapter 4](#timings) guard this.

## Mass transfer in an eccentric orbit {#eccentric-transfer}

SWINGERsim's mass transfer is formulated for circular orbits, as in BSE. A tertiary can, however, force a
small eccentricity on a transferring pair, or pump it by ZLK cycles. Two effects followed. First, a star in
a slightly eccentric orbit fills its lobe at periastron earlier than in a circular one (in the Algol triple
by 11 kyr); the onset check allows for the forced eccentricity, so the run does not stop before the
transfer. Second, without tides the tertiary of the Algol triple pumps the eccentricity of the transferring
pair to 0.44 (0.33 at the end of the run), outside the validity of the circular formulation; SWINGERsim
prints a warning that names the binary and the eccentricity. With tides switched on, as in the run of
[Chapter 8](#algol-triple), the pair stays circular ($e$ < 4·10⁻⁴) and the masses at the end are the same
(0.454 and 3.172 M☉). What a user should do: switch on tides for any binary that transfers mass in a
multiple system.

## Chaotic runs and ill-conditioned starts {#chaos}

In chaotic configurations two correct integrations diverge after a finite time: the N-body and secular
runs of the Fig. 4 set-up of Vokrouhlický (2016) after about 4 kyr (the paper reports the same for its two
integrators), SWINGERsim and MSE after three ZLK cycles of the 2+2, and individual members of the Pejcha
ensembles. A similar sensitivity appears for regular orbits started next to a separatrix: a ZLK system
started at $e_0$ = 0.01 differs by 10 % in the timing of its cycles after two cycles between the N-body and
the secular backend (0.5 % for constant masses away from the separatrix). Neither is an error. What a user
should do: compare chaotic systems statistically (fractions of ensembles, as in [Chapter 2](#pejcha2013),
where the two backends agree within the errors), and avoid comparing single members or times beyond the
Lyapunov time.

## The secular backend near the stability limit {#secular-limit}

The secular backends average over the orbits and therefore assume a stable hierarchy with well-separated
periods. Near the stability limit this fails: in ASASSN-V J075203 (period ratio 6–9) the N-body apsidal
period of A is about half the secular one, and in CzeV4315 half of the secular 2+2 solutions turn out to be
chaotic or to collide in N-body ([Chapter 6](#squadra-systems)). A secular run cannot see instability; it
reports the MA01 and V22 ratios and the octupole and "single-averaging" parameters of every nested triple
in its summary and marks mild hierarchies ("secular estimates approximate"). What a user should do: use
the secular backend to pre-select, then check the members that matter with the N-body backend and MEGNO.

## What an ensemble weight of one member means {#weights}

The evidence of a model is the mean weight of its members, and its reliability is measured by the
effective number of members $n_\mathrm{eff}$ ([Chapter 6](#squadra-systems)). When the observations are
precise and the prior broad, almost all weight can fall on one member. ASASSN-V J075203 with an isotropic
prior gives $n_\mathrm{eff}$ = 1.0: one member out of 20 000 carries everything, so the "evidence" is the
weight of a single lucky draw and says nothing about the model. The 2+2 of CzeV4315 conditioned on the
eclipses of A has $n_\mathrm{eff}$ = 4.8, and a comparison of models with the same 119-d prior had
$n_\mathrm{eff}$ = 1.6 and is therefore only a noisy estimate. SWINGERsim prints $n_\mathrm{eff}$ next to every
evidence and every weighted fraction. What a user should do: trust an evidence only with $n_\mathrm{eff}$ of
tens or more; otherwise run more members, focus the prior on the region that fits (as the targeted ensemble
of J075203 did, $n_\mathrm{eff}$ = 53), or use an emulator or simulation-based inference.

## Quantities that are not defined {#undefined}

Some quantities a user may ask for do not exist for a given orbit. The apsidal period of a nearly circular
orbit is the clearest case: in TIC 120362137 the tertiary forces the eccentricity of the 3.3-d pair to
oscillate between 0.0001 and 0.016, the osculating $\omega$ jumps, and a "measured" apsidal period is
meaningless (the paper gives only a theoretical value). The same holds for Algol ($e\approx0.001$) and for
orbits with $e$ < 10⁻⁶. SWINGERsim reports such a period as "not measurable in the run" rather than
printing a number. Long-term extremes of strongly interacting members, in which small numerical
differences grow, are similarly ill-defined and are not used as results.

## Stellar tracks at their limits {#track-limits}

Three situations occur when a configuration asks the stellar tracks for something they cannot give.
(i) For very massive stars with strong winds the present mass of a star at a given age can be reached by
two tracks of different initial mass (for example 79 M☉ at 2.5 Myr by tracks of about 108 and 139 M☉ in
QZ Car ensembles); SWINGERsim takes the least evolved star and warns. (ii) A mass that does not exist at the
requested age on the chosen grid (Aa1 of QZ Car above 73.3 M☉ on the rotating tracks at 2.5 Myr) cannot be
started; the member ends with fate *error* and the reason. (iii) A star that already overfills its Roche
lobe at the start is refused with the reason. The refused members are information, not failures: in QZ Car
they bound the mass of Ac1 from above ([Chapter 5](#qzcar-evolution)). Errors are counted separately in
every ensemble summary, and a user can rerun them after a change with `--resume --retry-errors`.

## Tables of model atmospheres {#atmospheres}

The limb-darkening tables of Claret & Bloemen (2011) stop at 3500 K and cover $\log g$ from 0 to 5. A
2867 K long-period variable or a supergiant with $\log g$ < 0 uses the nearest node, with a warning that names
the star and the values. For the relative light of such a star in a multiple system this is adequate; for
its own eclipse profile it is not, and tables for cooler stars are outside the scope of SWINGERsim.

## Conserved quantities as diagnostics {#conservation}

Every run reports its relative energy and angular-momentum errors. With relativity the conserved quantity
is the 1PN Hamiltonian, not the Newtonian energy; SWINGERsim measures the error against the Hamiltonian of
REBOUNDx's `gr_full`, otherwise a correct relativistic run would show a spurious drift. With dissipation or
mass loss the energy is not conserved by the physics, and the diagnostic is then the adiabatic invariant
($aM$ for winds, tested to 4·10⁻⁹ in [Chapter 3](#mass-loss)). Time symmetry is a further diagnostic of
conservative dynamics: failures forward and backward in time agree to 0.02 % in the CzeV4315 ensembles
([Chapter 6](#czev4315)).

## Problems in the published references {#published-problems}

Tracing every difference back to its paper found problems in the references themselves. They are listed
here because any code compared with these papers will meet them.

Table 46 lists them; each row names the reference, what is wrong or ambiguous in it, and how the benchmark deals with it.

**Table 46.** Problems in the published references found while reproducing them. Column 2 describes the
problem, column 3 what SWINGERsim's benchmark uses.

| Reference | Problem | What the benchmark uses |
|---|---|---|
| Naoz et al. (2013), Fig. 4 | the caption (angles 6.75° and 63.25°) and the quoted $a_2$ = 50 au imply different angular-momentum ratios (7.60 and 6.79); the published curve belongs to 7.60 | $a_2$ = 39.9 au, which gives 7.60 |
| Veras et al. (2011), eq. 17 | denominator $1-e_0\cos f$ does not follow from their eq. 4 | the form derived from eq. 4, $1+e_0\cos f$ |
| Antognini & Thompson (2016), Table 3 | stated normalisation $\pi(a_1^2+a_2^2)$; the numbers correspond to $\pi a^2$ | $\pi a^2$ |
| Zasche et al. (2023), V994 Her C | quoted relativistic share 6 %; their Table 3 gives 8.1 % | EXPLAINED with 8.1 % |
| Hurley et al. (2002), Algol example | the published history would need the orbit to gain angular momentum in the Hertzsprung gap | EXPLAINED; the conservative giant phase is a separate benchmark |
| Portegies Zwart & Meinen (1993) | capture radii quoted from Ray et al. (1987) belong to other conventions | EXPLAINED; McMillan et al. (1987) as the reference |
| Vynatheya et al. (2023), MA01 for 3+1 | the quoted score 0.56 does not follow from their own Table 4 (0.69) | EXPLAINED; set drawn near the boundary |
| Pejcha et al. (2013), Vokrouhlický (2016) | the fractions grow with the run length, which is not stated precisely enough | N-body PASS at 3σ; secular EXPLAINED |
| Kemp et al. (2024), KIC 4150611 | the dynamical mass ratio contradicts the size and period ratios of the same paper | masses consistent with the ratios |
| Pribulla et al. (2020), VW LMi | the quoted 78.6-yr apsidal period is a fit to three years of data; the paper's own models give 95.5–96.2 yr | comparison with the paper's N-body model |
| MSE (Hamers et al. 2021) | the argument of periastron is lost at an inclination of exactly 0° | inclinations of 0.5–2° |

The lesson is simple but important: a published number is a measurement of something, and before a
difference is called an error of either side, one has to establish what exactly was computed, from which
input, with which definition, and over which time.
