# Mass transfer and evolution {#mass-transfer-and-evolution}

Stars change their masses: they lose mass in winds, give mass to a companion when they fill
their Roche lobes, share an envelope in contact binaries, engulf a companion in a common
envelope, or merge with it. Every such change acts on all orbits of the system that contain the
star. This chapter describes how SWINGERsim models these processes and how it couples the
evolution of the stars (chapter [Stars](#stars)) to the dynamics (chapter [Dynamics](#dynamics)).

## Coupling stellar evolution to the dynamics {#coupling}

Stars change on nuclear time scales (Myr–Gyr), in most systems much more slowly than the
secular dynamics. With `physics.evolution: true` the coupling is an **operator splitting** as in
the code MSE (Hamers et al. 2021):

- **Secular backends.** The run is cut into segments within which no stellar radius, apsidal
  constant, moment of inertia (and, with mass loss, mass) changes by more than
  `run.evolution_max_dr` (default 1 %). Within a segment the stellar structure (radius, $k_2$,
  $r_g^2$, tidal lag) is held at its value in the middle of the segment; the orbital state is
  carried from segment to segment, and spins conserve their angular momentum when the moment of
  inertia changes. On the thermally pulsing AGB, where the radius peaks for only a few centuries
  at every pulse, the segment grid includes the nodes of the MIST tracks and is refined by
  bisection until no radius changes by more than `evolution_max_dr` between two grid points.
  Runs of millions to billions of years are streamed in chunks and reduced to their extremes,
  so they fit in memory.
- **N-body backend.** The stellar structure is updated at every output step (radius for
  collisions and Roche lobes, Love number, moment of inertia and lag for the tidal forces, spins
  at constant spin angular momentum), held at its value in the middle of the step. The output
  steps (`run.n_output`) must therefore be short compared with the evolution.

A run ends at the end of the track of any star (`stellar_evolution_end`: helium ignition with
SSE, the departure from the AGB or the end of carbon burning with MIST), backwards in time at the
birth of the stars (`stellar_birth`), and, unless mass transfer is followed, when a star fills
its Roche lobe (`roche_overflow`).

## Winds {#winds}

With `physics.mass_loss` the masses change during the run, from one of two sources:

- `source: track` – the mass of every star follows its track. MIST masses contain the wind
  prescriptions of the MESA models behind MIST (Vink et al. 2001 for hot stars with full
  strength, Nugis & Lamers 2000 for Wolf–Rayet stars, Reimers 1975 on the giant branch and
  Blöcker 1995 on the AGB, boosted for rotating stars by $(1-\Omega/\Omega_\mathrm{crit})^{-0.43}$;
  Choi et al. 2016). The SSE tracks are computed without winds; for them SWINGERsim adds the wind
  of the SSE code (Hurley, Pols & Tout 2000, Sect. 7.1): the larger of
  $9.6\times10^{-15}(Z/0.02)^{1/2}R^{0.81}L^{1.24}M^{0.16}$ M☉ yr⁻¹ (Nieuwenhuijzen & de Jager
  1990) above 4000 L☉ and $4\times10^{-13}\,\eta\,LR/M$ M☉ yr⁻¹ (Reimers 1975) from the
  Hertzsprung gap on, with $\eta$ = `reimers_eta` (default 0.5); $R$, $L$, $M$ in solar units.
- `source: rates` – constant rates $\dot M$ per star, for parametric studies and benchmarks.

The wind is isotropic in the frame of the star and carries the star's momentum (the *Jeans mode*;
Jeans 1924; Hadjidemetriou 1963). Its effect on an orbit X of total mass $M_X$ is measured by the
mass-loss index (Veras et al. 2011, eq. 15)

$$ \Psi = \frac{|\dot M_X|}{n_X M_X}, $$

the ratio of the orbital period to the mass-loss time scale divided by $2\pi$ ($n_X$ the mean
motion).

- **N-body backend.** The masses change during the integration with the REBOUNDx operator
  `modify_mass` ($\dot m = m/\tau$ with velocities unchanged), with $\tau$ of every star set at
  every output step so that the mass reaches its value at the end of the step. This is exact for
  any $\Psi$, from the adiabatic regime ($\Psi \ll 1$: $aM$ constant, $e$ constant on average) to
  the impulsive one ($\Psi \gg 1$: position and velocity are kept while the mass drops; a
  circular binary that loses half its mass at once is unbound; Blaauw 1961; Hills 1983).
- **Secular backends.** Orbit-averaged equations need the adiabatic response: every orbit keeps
  its eccentricity vector and orientation and $a_X M_X$ is constant (Hadjidemetriou 1963). A run
  is refused when $\Psi$ of any orbit exceeds `physics.mass_loss.max_psi` (default 0.05; the
  transition to the non-adiabatic regime lies at $\Psi \sim 0.01$–0.1, Veras et al. 2011); the
  largest $\Psi$ of every orbit is stored. The test-particle ZLK cycles with a perturber of
  changing mass then agree with the N-body backend to 0.2 % in timing.
- **Spins.** The wind carries away spin angular momentum at the rate
  $\dot J = \tfrac{2}{3}\dot M R^2\Omega$ (a thin spherical shell; Hurley et al. 2002, eq. 11), so
  $J \propto M^{2/(3 r_g^2)}$.

**Triple evolution dynamical instability.** Mass loss from an inner binary widens the inner orbit
by $M_{\mathrm{in},0}/M_\mathrm{in}$ and the outer one only by the smaller factor
$M_{\mathrm{tot},0}/M_\mathrm{tot}$, so the hierarchy weakens and a stable triple can cross the
stability limit (Perets & Kratter 2012; Toonen et al. 2022). In a benchmark the inner eccentricity
starts to grow within 7 % of the time at which the adiabatically expanding orbits cross the
Mardling & Aarseth (2001) limit; the break-up itself takes many outer orbits after that, as
Toonen et al. found.

## Roche-lobe overflow and stable mass transfer {#roche-lobe-overflow}

A star of a two-star binary that fills its Roche lobe ([Stars](#roche-lobes)) gives mass to its
companion. With `physics.mass_transfer` (which needs `physics.evolution` with SSE or MIST stars)
SWINGERsim follows stable transfer with the binary-star evolution algorithm (BSE) of Hurley, Tout
& Pols (2002, Sect. 2.6). Every two-star binary is treated as circular – tides circularise close
binaries before the overflow – and its **history** (masses, separation, stellar states) is
computed before the dynamics from the stellar tracks and the initial separation. The dynamics
then follows that history (section [Histories in the dynamics](#histories)).

![**Figure 9.** Mass transfer, schematically. The donor (orange, left) fills its Roche lobe; a
stream of gas (dark orange line) leaves it through the inner Lagrange point L1 and reaches the
accretor (blue, right), either directly or through an accretion disc (light orange ellipse).
The grey note at the top recalls that part of the transferred mass may leave the binary. The
masses of both stars, the orbit and both stars' structure change together.](../../figures/concepts/mass_transfer.png)

Figure 9 names the ingredients that the equations below quantify: the rate at which the donor
loses mass, the fraction the accretor keeps, the angular momentum carried by the mass that
leaves, and the response of both stars.

**Rate** (Hurley et al. 2002, eqs. 58–63). A donor of radius $R$ overfilling its lobe $R_L$
transfers

$$ \dot M = 3\times10^{-6}\,[\min(M, 5)]^2\,[\ln(R/R_L)]^3\ \mathrm{M_\odot\,yr^{-1}} $$

($M$ in M☉), so it stays just above its lobe while its expansion or the shrinking lobe drives
the transfer. The rate is limited by the thermal rate $M/\tau_\mathrm{KH}$ of giant-like donors
(Hertzsprung gap and later; $\tau_\mathrm{KH} = 10^7\,M M'/(RL)$ yr with the envelope mass
$M' = M - M_c$, $M_c$ the core mass, solar units) or by the dynamical rate.

**Stability.** Transfer is dynamically unstable and leads to a common envelope when the mass
ratio $q = M_\mathrm{donor}/M_\mathrm{accretor}$ exceeds a critical value $q_\mathrm{crit}$
(Table 8).

**Table 8.** Critical mass ratios $q_\mathrm{crit} = M_\mathrm{donor}/M_\mathrm{accretor}$ above
which transfer is dynamical (Hurley et al. 2002, Sect. 2.6). *Donor* is the evolutionary phase
of the donor; $x = 0.3$ is the exponent of the giant mass–radius relation and $M_c/M$ the core
mass fraction.

| Donor | $q_\mathrm{crit}$ | Note |
|---|---|---|
| main sequence below 0.7 M☉ (deeply convective) | 0.695 | |
| main sequence above 0.7 M☉, core helium burning | 3.0 | `q_crit_ms`; e.g. 1.6 of binary_c (de Mink et al. 2007) |
| Hertzsprung gap | 4 | |
| giant branch, AGB | $(1.67 - x + 2(M_c/M)^5)/2.13$ | eq. 57 |

With `stability: zeta` giant donors are judged instead by the adiabatic response of a condensed
polytrope (Hjellming & Webbink 1987; Soberman et al. 1997, eqs. 60–61),

$$\zeta_\mathrm{ad} = \frac{2}{3}\frac{m_c}{1-m_c} - \frac{1}{3}\frac{1-m_c}{1+2m_c} - 0.03\,m_c + \frac{0.2\,m_c}{1+(1-m_c)^{-6}},\qquad m_c = \frac{M_c}{M},$$

where $\zeta_\mathrm{ad} = d\ln R/d\ln M$ of the donor on the dynamical time scale; transfer is
unstable when $\zeta_\mathrm{ad}$ falls below the exponent $\zeta_L = d\ln R_L/d\ln M$ of the Roche
lobe for the actual accretion fractions.

**Accretion** (eqs. 64–65). A main-sequence, Hertzsprung-gap or core-helium-burning accretor takes
at most $10\,\tau_{\dot M}/\tau_\mathrm{KH,2}$ of the transferred mass ($\tau_{\dot M} = M_2/\dot M$,
$\tau_\mathrm{KH,2}$ its thermal time scale); the rest leaves with the donor's specific orbital
angular momentum. Giant accretors take everything. With `accretion: fractions` the fractions
`alpha` and `beta` below are constant.

**Orbit.** Of the mass $|dM_d|$ lost by the donor a fraction $\alpha$ leaves from the donor's
vicinity, a fraction $\beta$ from the accretor's (isotropic re-emission), and $1-\alpha-\beta$ is
accreted. With the orbital angular momentum $J = M_d M_a\sqrt{Ga/M}$ of a circular binary
($M_d$, $M_a$ the donor and accretor masses, $M = M_d + M_a$),

$$\frac{d\ln a}{dM_d} = 2\,\frac{\alpha M_a^2 + \beta M_d^2}{M M_d M_a} + \frac{\alpha+\beta}{M}
- \frac{2}{M_d} + \frac{2(1-\alpha-\beta)}{M_a}.$$

This reproduces the closed forms of Soberman et al. (1997) and Postnov & Yungelson (2014) to
$10^{-12}$: conservative transfer keeps $a M_d^2 M_a^2$ constant, the Jeans mode $aM$.

**Response of donors and accretors.** Every star carries its mass, its age on the track it
follows and the mass of that track (Hurley et al. 2000, Sect. 7.1; 2002, Sect. 2.6.6). A
main-sequence star, or a Hertzsprung-gap star whose new mass would not shrink its helium core,
follows the track of its present mass at the same evolutionary state (EEP for MIST, phase
fraction for SSE): a donor shrinks with its mass and crosses the gap on the longer time scale of
its new mass, an accretor brightens. With `rejuvenation: hurley` (default) an accretor with a
convective core (above 1.25 M☉) or a fully convective one (below 0.35 M☉) mixes the fresh
hydrogen into its burning region and is rejuvenated further,
$t' = (M/M')\,(t'_\mathrm{MS}/t_\mathrm{MS})\,t$, where $t$ and $t'$ are the ages before and after,
$M$, $M'$ the masses and $t_\mathrm{MS}$, $t'_\mathrm{MS}$ the main-sequence lifetimes (Hurley et
al. 2000). A gainer therefore looks younger and less evolved than a single star of its mass.
Giants keep the track of the mass they had when they left the gap (their structure follows the
core), and their radius follows the envelope mass with $\zeta_\mathrm{eq} = -x + 2(M_c/M)^5$,
$x = 0.3$ (eq. 56). The helium core never shrinks.

**Events** that end a history: `common_envelope` (unless the common envelope is followed),
`contact` (the accretor fills its own lobe, unless contact binaries are followed),
`donor_stripped` (the donor is down to its helium core; helium stars are not followed).

**Check.** For the BSE Algol example (2.9 + 0.9 M☉, 2.7 d) SWINGERsim ends with the donor mass of
BSE (0.404 M☉; BSE in COSMIC, Breivik et al. 2020: 0.401; the paper 0.42); the accreted mass is
the least certain number of that example, and the published numbers are not self-consistent (see
*Testing, Part I*).

## Histories in the dynamics {#histories}

In the dynamics the masses of the two stars of a transferring binary follow its history. The
separation is scaled at every operator-splitting step (secular backends) or output step (N-body:
the relative orbit $\mathbf{r} \to g\mathbf{r}$, $\mathbf{v} \to \mathbf{v}/\sqrt g$ at fixed
eccentricity vector and phase) by the change beyond the adiabatic response $aM$ = const, which
the backends apply themselves. Outer orbits respond adiabatically to the mass that leaves the
binary. During the transfer the donor's radius is its Roche-lobe radius (for eclipses, tides and
diagnostics), and the overflow of the transferring star does not end the run. The re-modelled
stars carry their radius, luminosity and temperature into the dynamics and the light curves.
Winds and transfer combine: the stars lose mass by their winds before, during and after the
transfer, and the binary widens by $aM$ = const for the wind part.

### Transfer started by the dynamics {#dynamic-onset}

The histories are computed for circular binaries of the configured separations. The dynamics can
drive a star over its lobe that no history expects: ZLK cycles raise the eccentricity of an inner
binary until a star touches its lobe at periastron (Eggleton & Kisseleva-Eggleton 2006; Fabrycky &
Tremaine 2007). Tides then circularise the binary quickly; triple codes assume that the orbit is
circular at the onset (Toonen et al. 2016). With `physics.mass_transfer.dynamic_onset` SWINGERsim
ends the run at the contact – the eccentricity $e_c$ at which $a(1-e_c)\,r_L = R$, with $r_L$ the
lobe in units of the separation – and continues with the circularised binary: at constant orbital
angular momentum, $a_c = a(1-e_c^2)$ (option `angular_momentum`; Hut 1981), or at the periastron,
$a_c = a(1-e_c)$ (option `periastron`, the donor stays at its lobe). The pair keeps its centre of
mass, plane and phase, its stars rotate synchronously, and a new mass-transfer history starts
there. In a triple driven to $e_c = 0.89$ the N-body and the secular backend find the contact within
1.7 % in time and the same circular orbit.

## Semi-detached and contact binaries {#semidetached-contact}

**Imposed period change.** Algols and W UMa binaries change their periods by mass transfer and by
losing angular momentum; their minima then lie on a parabola,
$O - C = \tfrac12 \dot P P E^2$ ($E$ the cycle number; Kwee 1958; Hilditch 2001). When the cause
need not be modelled, `orbits.<pair>.Pdot` imposes a constant $\dot P$ (d yr⁻¹) on a two-star
binary inside the N-body run, $P(t) = P + \dot P\,t$. The time axis is cut into intervals of length
$\Delta t$, and where the integration crosses a boundary (in either direction) the relative orbit
is scaled to the period of the new interval, $a \to a(P_\mathrm{new}/P_\mathrm{old})^{2/3}$, at
fixed eccentricity vector and phase and with the centre of mass of the pair unchanged. The eclipse
times then deviate from the parabola by at most $(\dot P/P)\Delta t^2/8$, kept below 0.05 s.
`Pdot_mode: conservative` changes the masses with the period as in conservative transfer,
$\dot P/P = 3\dot M_1 (M_1 - M_2)/(M_1 M_2)$. All observables and the third body see the same
drifting orbit, so the parabola and the light-time effect can be fitted together
([Fitting](#free-parameters)). The quadratic coefficient of the computed minima equals
$\tfrac12 P\dot P$ to $5\times10^{-5}$. This option exists only in the N-body backend.

**Contact binaries.** When both stars fill their lobes, `physics.mass_transfer.contact` lets them
share a common envelope (Lucy 1968) instead of ending the history: the less massive star stays at
its lobe and feeds the primary through the envelope (conservatively), while the magnetised wind
of the envelope removes angular momentum,

$$\dot J = -4.9\times10^{41}\,\frac{R_1^2 M_1 + R_2^2 M_2}{P}\ \mathrm{g\,cm^2\,s^{-1}\,yr^{-1}}$$

($R$, $M$ in solar units, $P$ in days; Gazeas & Stępień 2008, eq. 19). The orbit shrinks and the
mass ratio falls; the binary merges when it reaches the Darwin instability ($q < q_\mathrm{min}$,
default 0.09; Rasio 1995) or overfills its outer critical surface through L2 (events
`contact_merger_darwin`, `contact_merger_l2`). Thermal relaxation oscillations and energy transfer
between the components beyond one common temperature are not modelled.

## Tertiary mass transfer {#tertiary-transfer}

A star that orbits a two-star binary (for example the giant tertiary of HD 181068) and fills its
Roche lobe transfers mass onto that binary (`physics.mass_transfer.tertiary`, on by default,
secular backend). SWINGERsim computes the history of the three stars before the dynamics, like that
of a two-star binary: the donor sees the inner binary as a point mass, so its lobe is that of a
binary with $q = M_3/(M_1+M_2)$ at the periastron of the outer orbit (Eggleton 1983), and the rate
and its limits follow BSE (Hurley et al. 2002; Hamers et al. 2022). The stream forms a circumbinary
disc when the inner separation is below the circularisation distance
$a_c = a_\mathrm{peri}(1 - R_L)^2$, where $R_L$ is the Roche-lobe radius of the donor in units of
the outer separation (Eggleton 1983). Leigh et al. (2020, eq. 5) write $a_c = a_\mathrm{out}(1 - R_L)^2$
with the outer semi-major axis; SWINGERsim evaluates it at the periastron $a_\mathrm{peri}$ of the
outer orbit, where the donor overflows. With a disc, 90 % of the mass is accreted (two thirds by the more massive inner star), otherwise 20 %
(Hamers et al. 2022). The rest leaves the inner binary isotropically; for the outer orbit it is
re-emitted by the accretor (Soberman et al. 1997), and its ejection costs inner orbital energy
(de Vries et al. 2014):

$$\mathrm{d}\ln a_\mathrm{in} = \frac{\mathrm{d}M_1}{M_1} + \frac{\mathrm{d}M_2}{M_2} - \frac{2 M_\mathrm{in}\,\mathrm{d}M_\mathrm{ej}}{\alpha\lambda\, M_1 M_2},$$

with the inner masses $M_1$, $M_2$, their sum $M_\mathrm{in}$, the ejected mass
$\mathrm{d}M_\mathrm{ej}$ and the efficiency $\alpha\lambda$ = 5 (de Vries et al. found 3–7 from
hydrodynamics); `inner_orbit: angular_momentum` instead lets the accreted gas bring the specific
angular momentum $h\sqrt{G M_\mathrm{in} a_\mathrm{in}}$. The history stops at a triple common
envelope, when an inner star fills its own lobe, when the triple becomes dynamically unstable
(Mardling & Aarseth 2001), when the donor is stripped to its core or at the end of its track. Not
modelled: the drag of the stream on the inner orbit beyond the energy formalism, the change of
the outer eccentricity during the transfer, and episodes at every periastron of an eccentric outer
orbit (the rate is orbit-averaged). Outside the secular backend, or without this option, the
overflow of a tertiary ends the run; ensembles record which star overflowed first and whether its
companion is a star or a binary.

## Common envelope {#common-envelope}

When the transfer becomes dynamically unstable, the donor's envelope engulfs both stars. With
`physics.mass_transfer.common_envelope` SWINGERsim computes the outcome with the energy formalism of
Hurley, Tout & Pols (2002, Sect. 2.7.1; after Webbink 1984 and Tout et al. 1997) instead of
stopping:

- binding energy of the envelopes (eq. 69),
  $E_\mathrm{bind} = -\frac{G}{\lambda}\left(\frac{M_1 M_\mathrm{env,1}}{R_1} + \frac{M_2 M_\mathrm{env,2}'}{R_2}\right)$,
  with the structure parameter $\lambda$ (`lam`, default 0.5) and a main-sequence companion
  counted as all core ($M_{c,2}' = M_2$, $M'_\mathrm{env,2} = 0$);
- the cores spiral in and give orbital energy to the envelope with the efficiency
  $\alpha_\mathrm{CE}$ (`alpha`, default 1; eqs. 70–72),
  $E_\mathrm{bind} = \alpha_\mathrm{CE}\left(-\frac{GM_{c,1}M_{c,2}'}{2a_f} + \frac{GM_{c,1}M_{c,2}'}{2a_i}\right)$,
  which gives the final separation $a_f$ from the initial $a_i$;
- if a core would fill its Roche lobe at $a_f$, the cores merge where this happens, and the merger
  keeps the part of the envelope the released energy could not unbind (eqs. 73–77). Core radii:
  the white-dwarf radius for degenerate helium cores and AGB cores, the zero-age helium-star radius
  otherwise (Hurley et al. 2000).

The envelope leaves isotropically within `duration_yr`. In the N-body backend the binary's centre
of mass keeps its position and velocity while the stars take their new masses and a circular orbit
of separation $a_f$; the outer orbits feel an impulsive mass loss and widen, become eccentric
($e' = f/(1-f)$ for a circular orbit losing the fraction $f$ of its mass) or unbound ($f \ge 1/2$;
Hadjidemetriou 1963; Veras et al. 2011). The secular backends stop at the end of the envelope,
because the impulsive response depends on the orbital phase. After the envelope the core is a
point mass of constant mass (remnants are not evolved); the companion continues on the track of
its present mass. The outcome is validated on published BSE examples and on the analytic impulsive
response (see *Testing, Part I*). A tertiary that fills its lobe in an eccentric outer orbit can
also lead to a triple common envelope (Glanz & Perets 2021; Hamers et al. 2022), which ends the run.

## Mergers {#mergers}

With `physics.mergers` a run continues after two stars merge – a collision in the N-body backend,
or a common envelope or contact binary that ends in a merger – with one merged star and the
remaining orbits (a quadruple becomes a triple). Collisions in bound few-body systems are slow
(relative velocity about the escape velocity), so the stars merge rather than destroy each other.
The merger ejects the fraction $\phi = 0.3\,q/(1+q)^2$ of the total mass, $q = M_2/M_1 \le 1$
(head-on collisions of main-sequence stars, Glebbeek et al. 2013: 2–7.5 %; Lombardi et al. 2002:
at most about 8 %); a common-envelope merger keeps the mass of the energy formalism. The product
sits at the centre of mass of the two stars with their total momentum; the mass leaves at once,
so the outer orbits feel an impulsive mass loss. It mixes fresh hydrogen into its core: its
fractional main-sequence age is $f = r\,(M_1 f_1 + M_2 f_2)/(M_1 + M_2)$, with $r$ from 0.81
(M ≤ 6 M☉) to 0.64 (M ≥ 60 M☉) (Glebbeek & Pols 2008; Schneider et al. 2016) or 0.1 of BSE
(Hurley et al. 2002). Born near break-up, the product starts at a fraction `spin_fraction`
(default 0.3) of the break-up rate and spins down quickly (Schneider et al. 2019).

## Magnetic braking and activity {#magnetic}

**Magnetic braking.** With `physics.magnetic_braking` (requires `evolution` and `rotation`) the
spins of stars with convective envelopes slow down by the torque of Matt et al. (2015),

$$T = -T_0\left(\frac{\tau_\mathrm{cz}}{\tau_\mathrm{cz,\odot}}\right)^p\left(\frac{\Omega}{\Omega_\odot}\right)^{p+1}
\ \text{(unsaturated)},\qquad T = -T_0\,\chi^p\,\frac{\Omega}{\Omega_\odot}\ \text{(saturated)},$$

with $T_0 = 9.5\times10^{30}\,\mathrm{erg}\,(R/\mathrm{R}_\odot)^{3.1}(M/\mathrm{M}_\odot)^{0.5}$,
$\chi = 10$, $p = 2$, $\Omega_\odot = 2.6\times10^{-6}$ s⁻¹, saturation below the Rossby number
$\mathrm{Ro}_\odot/\chi$ ($\mathrm{Ro} = 1/(\Omega\tau_\mathrm{cz})$), and the convective turnover
time of Cranmer & Saar (2011),
$\tau_\mathrm{cz} = 314.24\,\exp[-T_\mathrm{eff}/1952.5\,\mathrm{K} - (T_\mathrm{eff}/6250\,\mathrm{K})^{18}] + 0.002$ d
(12.9 d for the Sun, nearly zero above 6500 K). Within each segment the moment of inertia is
constant and $d\Omega/dt = T/I$ is integrated exactly. With tides, a braked star in a close binary
is spun up again by its companion and the binary loses orbital angular momentum: magnetic braking
shrinks close binaries, a path to contact binaries in triples. The solar spin at the solar age and
the Skumanich law $\Omega \propto t^{-1/2}$ (Skumanich 1972) are reproduced; for stars below 0.6 M☉ the torque is
known to be too weak (Douglas et al. 2017).

**Applegate mechanism.** A magnetically active star whose quadrupole moment changes through its
activity cycle modulates the period of its binary (Applegate 1992). `swingersim applegate` compares
the energy the modulation needs with the energy the star radiates in one modulation period
(thin-shell model of Tian et al. 2009 as used by Völschow et al. 2016,
$\Delta P/P = 2\pi\,(O-C)/P_\mathrm{mod}$ with the full O−C amplitude) and gives the largest O−C
semi-amplitude that activity can produce. A timing variation larger than that cannot be activity;
the two-zone models of Völschow et al. need one to two orders of magnitude more energy, so the
thin-shell limit is generous. With `physics.activity` a cycle becomes part of the forward model:
the eclipses of the binary are shifted by $K\sin(2\pi t/P_\mathrm{mod} + \varphi)$ in the synthetic
timings and minima.

## Limits {#mass-transfer-limits}

- Mass-transfer histories assume circular binaries; eccentric transfer is entered only through the
  dynamically started onset with instantaneous circularisation.
- Helium stars, white dwarfs, neutron stars and black holes are not evolved; a common-envelope core
  is a point mass of constant mass.
- Winds are isotropic and leave the system; wind accretion by a companion is not modelled.
- Rates and stability follow BSE recipes and carry their uncertainties (critical mass ratios,
  $\alpha_\mathrm{CE}\lambda$, accretion efficiency).

## Implementation {#implementation-mass-transfer}

**External library and data.** The N-body mass change uses the REBOUNDx operator `modify_mass`
(Tamayo et al. 2020); MIST winds come with the MIST tracks (Choi et al. 2016). **Own code.**
Everything else was written for SWINGERsim from the cited papers: the operator splitting of
stellar evolution and dynamics with adaptive segments and streaming; the SSE wind prescription
integrated along the SSE tracks; the adiabatic response and the $\Psi$ check of the secular
backends; the BSE mass-transfer algorithm (rates, limits, critical mass ratios, accretion
efficiency, the $d\ln a/dM_d$ equation, re-modelling of donors and accretors, rejuvenation); the
$\zeta$ stability criterion; the dynamically started onset with circularisation; the imposed
period change; contact binaries with angular-momentum loss and their merger criteria; tertiary
mass transfer; the common-envelope energy formalism; mergers; magnetic braking; and the Applegate
energy check and activity modulation.
