# Stars {#stars}

The dynamics of a multiple star needs more from each star than its mass. Its radius decides
collisions, Roche-lobe overflow and eclipses; its internal structure, through the
apsidal-motion constant $k_2$ and the moment of inertia, decides how strongly it is distorted
by tides and rotation; its luminosity and effective temperature decide its light. All of these
change as the star evolves. This chapter describes where SWINGERsim takes these quantities from,
how they evolve, and how a star's size compares with its Roche lobe. Mass loss and mass transfer,
which change the masses themselves, are the subject of chapter
[Mass transfer and evolution](#mass-transfer-and-evolution).

## The Hertzsprung–Russell diagram {#hr-diagram}

The *Hertzsprung–Russell (HR) diagram* plots the luminosity $L$ of a star (the energy it
radiates per second, in units of the solar luminosity L☉) against its effective temperature
$T_\mathrm{eff}$ (the temperature of a black body with the same surface flux; hot stars on the
left). The two are linked to the radius $R$ by $L = 4\pi R^2 \sigma T_\mathrm{eff}^4$, with
$\sigma$ the Stefan–Boltzmann constant. Most stars lie on the **main sequence**, where they burn
hydrogen in their cores: the heavier, the hotter, brighter and shorter-lived. Before it, a
young star contracts on the **pre-main sequence**; after it, the core contracts and the
envelope swells into a **giant**, later an **asymptotic-giant-branch (AGB)** star. Massive stars
lose their hydrogen envelopes in strong winds and become hot **Wolf–Rayet** stars.

![**Figure 5.** Where the members of a modelled system can lie in the HR diagram. Horizontal
axis: effective temperature in kelvin, decreasing to the right (logarithmic); vertical axis:
luminosity in L☉ (logarithmic). Black lines: MIST evolutionary tracks of stars of 0.5, 1, 2, 5,
15 and 60 M☉ (labelled at their zero-age main-sequence points), from the main sequence to the
end of the track; grey dashed lines: the same stars on the pre-main sequence; red: the
Wolf–Rayet phase of the 60 M☉ star. The orange band is the classical instability strip in which
stars pulsate. Boxes name the regions. Any star of a modelled system can sit anywhere on these
tracks, so one system can combine, for example, a young binary with a giant.](../../figures/concepts/hr_zoo.png)

Figure 5 summarises the stellar side of the model. Each black track shows how one star moves
in the diagram with age: nearly fixed on the main sequence, then to the right and upwards as it
becomes a giant. The radius, which matters most for the orbits, grows by a factor of 100 and
more along the giant branch and the AGB. The instability strip marks the region where
pulsations make a star a clock ([Observables](#pulsating-stars)).

## Stellar models {#stellar-models}

A *stellar model* maps the mass of a star (and its age and metallicity, for evolving models) to
its radius, luminosity, effective temperature, apsidal-motion constant $k_2$ and gyration
constant $r_g^2 = I/(MR^2)$ ($I$ the moment of inertia). It is chosen with
`stellar_model.kind`; Table 3 lists the options.

**Table 3.** Stellar models. *Kind* is the value of `stellar_model.kind`; *Source* the model or
data behind it; *Range* the masses, metallicities and phases covered; *Structure* where $k_2$
and $r_g^2$ come from. Values given for a single star under `stars.<name>` (`R`, `Teff`,
`apsidal_k2`, `rg2`, `lag_s`) always override the model.

| Kind | Source | Range | Structure ($k_2$, $r_g^2$) |
|---|---|---|---|
| `zams` (default) | zero-age main-sequence radius of Tout et al. (1996) | 0.1–100 M☉, solar metallicity, no ageing | polytrope, $n = 1.5$ below 0.35 M☉, $n = 3$ above |
| `sse` | single-star evolution formulae of Hurley, Pols & Tout (2000) | 0.1–8 M☉, $Z$ = 0.004–0.03; main sequence, Hertzsprung gap, giant branch | Claret (2023) for 0.2–8 M☉, polytrope outside |
| `mist` | MIST v1.2 tracks (Choi et al. 2016; Dotter 2016) | 0.1–300 M☉, [Fe/H] −1.5 to +0.5, $v/v_\mathrm{crit}$ = 0 or 0.4; pre-main sequence to the end of the AGB | Claret (2023) for 0.2–8 M☉, Claret (2004) above 8 M☉ |
| `table` | a CSV table of the user (mass, radius, $k_2$, $r_g^2$, optional lag) | as tabulated, linear in mass | from the table |
| `python` | a function of the user returning the same quantities | any | from the function |

With `structure: mesa:<grid>` the structure constants come from a grid of the user's own MESA
models (Jermyn et al. 2023) instead (section [Own MESA grids](#own-mesa-grids)).

**Radius, temperature and luminosity of a single star.** When only the radius of a star is
given, the luminosity of the model is kept and the effective temperature follows from
$L = 4\pi R^2 \sigma T_\mathrm{eff}^4$; when $R$ and $T_\mathrm{eff}$ are both given, they fix
the luminosity. `radius_inflation` enlarges a star relative to the model at fixed luminosity.
Magnetically active low-mass stars are larger than standard models predict – CM Dra by 3–4 %
at solar metallicity and 5–7 % at a lower metallicity (Feiden & Chaboyer 2014) – and MIST
shows the same deficit.

### The default: zero-age main-sequence stars {#zams-model}

The `zams` model gives every star the zero-age main-sequence (ZAMS) radius of the analytic fit
of Tout et al. (1996) for solar metallicity (1 M☉ → 0.89 R☉, 3 M☉ → 1.98 R☉). The stars do not
age. The structure constants are those of a polytrope – a star in which pressure and density are
related by $P \propto \rho^{1 + 1/n}$ – with index $n = 1.5$ for fully convective stars below
0.35 M☉ and $n = 3$ above. SWINGERsim integrates the Lane–Emden equation for the polytrope
together with Radau's equation for its tidal distortion and obtains $k_2$ and $r_g^2$; it
reproduces the analytic values for $n = 1$ and the tabulated $k_2 = 0.1433$ ($n = 1.5$) and
0.0144 ($n = 3$) of Brooker & Olle (1955) (see *Testing, Part I*).

The polytropes are a rough approximation. Real main-sequence stars above about 1.5 M☉ are more
centrally condensed, with $k_2$ of a few $10^{-3}$, and $k_2$ falls further as the star evolves;
the $n = 3$ polytrope overestimates $k_2$ of 2–3 M☉ stars by about a factor of three. For any
quantitative apsidal motion the structure should come from stellar models (`sse`, `mist`, a table
or per-star values).

### Evolving stars: SSE {#sse-model}

With `kind: sse` the stars have an age `age_myr` (counted from the ZAMS, all stars coeval unless
`stars.<name>.age_offset_myr` shifts one) and a metallicity `z` ($Z = 0.02$ is solar in this
calibration). Radius, luminosity, effective temperature and the convective envelope follow from
the formulae of Hurley, Pols & Tout (2000), fitted to the detailed models of Pols et al. (1998),
for the main sequence, the Hertzsprung gap and the first giant branch. The mass and depth of the
convective envelope come from the binary-code extension of Hurley, Tout & Pols (2002), which is
calibrated for $Z = 0.02$.

The formulae are not evaluated at run time. They were evaluated once with the published SSE
routines on a grid of 320 masses and 5 metallicities; every phase is sampled at the same
fractions of the phase (*equivalent evolutionary points*, as in MIST), and a star of any mass
and metallicity is interpolated at fixed phase fraction, linearly in $\log M$ and $\log Z$. Against
direct SSE evaluations the interpolation is accurate to 0.8 % in radius on the main sequence,
1.7 % in the Hertzsprung gap and 4 % on the giant branch (95th percentiles for 0.8–3 M☉,
$Z$ = 0.01–0.03) and to 0.2 % in the phase boundaries; for massive metal-poor giants the errors
reach 30 %. Beyond helium ignition (core helium burning, AGB, remnants) SSE stars are not
followed; a run stops there. For the Sun (1 M☉, $Z$ = 0.02, 4.57 Gyr) the model gives
$R$ = 0.985 R☉, $L$ = 0.958 L☉ and the base of the convective envelope at 0.734 $R$, against
the helioseismic 0.713 $R$ (Christensen-Dalsgaard et al. 1991).

### Evolving stars: MIST {#mist-model}

With `kind: mist` SWINGERsim uses the MESA tracks of MIST v1.2 (Choi et al. 2016; Dotter 2016) as
an independent model family. They cover 0.1–300 M☉, [Fe/H] from −1.5 to +0.5 in steps of 0.25
and the two initial rotations of MIST, $v/v_\mathrm{crit}$ = 0 and 0.4 (the default, the MIST
standard). Every track is stored at the *equivalent evolutionary points* (EEPs) 1–1409:
pre-main sequence (1–201), main sequence (202–453), subgiant and red giant branch (454–604),
helium flash and core helium burning (605–707), early AGB (708–807) and thermally pulsing AGB to
the departure from the AGB (808–1409); massive stars end with core carbon burning. Post-AGB
stars and white dwarfs are not followed. A star is interpolated at fixed EEP between the
neighbouring masses (linearly in $\log M$) and metallicities, as in the MIST isochrones, and in
age between the EEPs. With twice the grid spacing the interpolation reproduces the radius at the
end of the main sequence to 1 % and the ages to 1–2 %. Ages count from the start of the MIST
pre-main-sequence models, not from the ZAMS as in SSE.

MIST tracks contain no convective-envelope structure, which the tidal dissipation needs
([Tides and rotation](#dissipation-from-structure)). SWINGERsim derives it: giants on the giant
branch and the AGB, and core-helium-burning stars cooler than 5500 K, have a convective
envelope containing all the mass above the helium core, as deep as the star (Hurley et al.
2002), shrinking linearly to zero between 5500 and 7000 K on the blue loop; main-sequence and
Hertzsprung-gap stars take the envelope mass and depth fractions of the SSE star of the same
mass at the same fraction of its phase; pre-main-sequence stars of 0.1–2.5 M☉ are fully
convective on the Hayashi track, and from the minimum of their effective temperature (where the
radiative core appears) the envelope shrinks linearly in EEP to its ZAMS size. For the Sun the
[Fe/H] = 0 track gives $R$ = 1.025 R☉ and $L$ = 1.11 L☉ at 4.57 Gyr, the SSE track 0.985 and
0.958: the two families bracket the Sun.

**Systematic differences between the families.** Against MIST the SSE main-sequence lifetimes
are about 10 % longer and the radii at the end of the main sequence up to 10 % smaller. This is
the systematic uncertainty of any age derived with one family.

**Massive stars.** Above 8 M☉ three properties of MIST matter for the dynamics. (1) The winds
are strong and depend on rotation: a 100 M☉ star ([Fe/H] = 0) leaves the main sequence with
68 M☉ without rotation and with 25 M☉ at $v/v_\mathrm{crit}$ = 0.4, so `vvcrit` is a choice of
the mass history ([Mass transfer and evolution](#winds)). (2) The hydrogen-burning lifetimes
agree with the Geneva models of Ekström et al. (2012) to 4–12 %, but the masses at the end of
the main sequence differ by up to 30 %; Choi et al. (2016) call 1D models of the
radiation-dominated envelopes of very massive stars uncertain. The masses of evolved massive
stars are therefore a model choice, not a prediction. (3) From the late main sequence on, the
envelopes of 40–70 M☉ stars inflate and contract within a few M☉ of initial mass; there the
grid does not resolve the radius (errors of 20–40 %), and grid masses should be used where the
radius matters.

**Present masses.** Masses from an orbit solution are present masses, whereas tracks are
labelled by initial masses. With `physics.mass_loss: {source: track, masses: current}`
SWINGERsim finds the initial mass of every star whose track has the given mass at `age_myr`.
Massive stars can be ambiguous (at 2.5 Myr a 26 M☉ star is either a 26.4 M☉ main-sequence star
or a 257 M☉ star stripped to a Wolf–Rayet star); by default the least evolved solution is taken
with a warning, and a mass that no track reaches at that age is refused as a contradiction.

## Structure constants: $k_2$ and $r_g^2$ {#structure-constants}

The **apsidal-motion constant** $k_2$ measures how strongly a star responds to a tidal or
rotational perturbation; it is half the potential *Love number* $k_{2,\mathrm{Love}} = 2k_2$ used
in tidal theory. A centrally condensed star (mass concentrated in the core) has a small $k_2$:
the $n = 3$ polytrope has 0.0144, a homogeneous body 0.75. The **gyration constant**
$r_g^2 = I/(MR^2)$ sets how much angular momentum a spinning star stores. Both enter the tidal
forces, the spin evolution and the apsidal motion
([Tides and rotation](#tides-and-rotation)).

With `sse` and `mist` both come from the MESA models of Claret (2023): 69 tracks, 0.2–8 M☉,
[Fe/H] = −0.5, 0, +0.5, from the pre-main sequence to the white dwarf. As Claret recommends, the
tables are interpolated in the surface gravity $\log g$ of the star – which carries the
evolutionary state independently of the track code – along a monotonic segment of each track
(pre-main sequence, main sequence and giant branch, core helium burning), then linearly in
$\log M$ and [Fe/H]; the [Fe/H] of an SSE star is $\log(Z/0.02)$. The values reproduce the table
rows to the published precision. On the main sequence $k_2$ decreases as the star becomes more
centrally condensed (for 2 M☉ from 0.005 to 0.003); on the giant branch the deep convective
envelope raises it by more than an order of magnitude. Above 8 M☉ (MIST) the constants come
from Claret (2004) (0.8–125 M☉, $Z$ = 0.02, ZAMS to carbon exhaustion), interpolated in the same
way and agreeing with Claret (2023) to 1 % at 8 M☉; above 125 M☉ the 125 M☉ track is used. AGB
stars take the values at the tip of the giant branch, because Claret's grids have no AGB models.

**Rotation.** A rotating star is more centrally condensed than a static one. With
`physics.k2_rotation` SWINGERsim corrects the model constants,
$\delta\log k_2 = -0.87\,\lambda_s$ with $\lambda_s = 2\Omega^2R^3/(3GM) = 2V^2/(3gR)$, where
$\Omega$ is the angular velocity of rotation, $V$ the equatorial rotation velocity, $g$ the
surface gravity, $R$ the radius and $M$ the mass (Claret 1999; Claret & Giménez 2010 apply this
correction to observed binaries), with the
spin of each star at that time.

### Own MESA grids {#own-mesa-grids}

For structure quantities that published grids do not provide, SWINGERsim can drive a local
installation of MESA (Paxton et al. 2011; Jermyn et al. 2023). A grid file lists masses,
metallicities and the physics (mixing length, overshooting, end point); every track starts on
the pre-main sequence and writes structure profiles. The reduction computes from each profile
the base of the convective envelope, the convective turnover time, $k_2$ from Radau's equation
(the method of Claret 2023), $I/(MR^2)$ and the inertial-wave dissipation coefficient
([Tides and rotation](#inertial-waves)). On polytropes written in the MESA format the reduction
reproduces $k_2$ to 1 %, $I/(MR^2)$ to 0.3 % and the inertial-wave coefficient of a fully
convective star to 1 %.

## Roche lobes {#roche-lobes}

In a binary that rotates with its orbit, the effective potential (gravity of both stars plus
the centrifugal term) has a saddle point between the stars, the inner Lagrange point L1. The
equipotential surface through L1 encloses two teardrop-shaped regions, the **Roche lobes**:
material inside a lobe belongs to its star. A star that swells beyond its lobe pours gas
through L1 onto its companion.

![**Figure 6.** Roche lobes of a circular synchronised binary with mass ratio $q = 0.5$, cut in
the orbital plane. The black dots are the centres of the heavier (left) and the lighter (right)
star. Grey dashed lines are equipotentials of the effective potential; the thick red dashed
figure-of-eight is the critical equipotential through the inner Lagrange point L1 (red cross),
which bounds the two Roche lobes. The large curves at the edges are equipotentials far from the
binary, where the centrifugal term dominates. A star whose surface reaches its lobe passes mass
to the companion through L1.](../../figures/concepts/roche.png)

Figure 6 shows the geometry. Close to each star the equipotentials are nearly spheres, so a
star much smaller than its lobe is almost spherical. Near the lobe the surfaces become pointed
towards L1: a star that nearly fills its lobe is visibly distorted, which matters for its light
curve ([Observables](#distorted-stars)). The lobe of the heavier star is larger.

SWINGERsim uses the *volume-equivalent* lobe radius $R_L$ (the radius of a sphere with the
lobe's volume) of Eggleton (1983),

$$\frac{R_L}{d} = \frac{0.49\,q^{2/3}}{0.6\,q^{2/3} + \ln(1 + q^{1/3})}, \qquad
q = \frac{M}{M_\mathrm{companion}},$$

where $d$ is the separation, $M$ the mass of the star and $M_\mathrm{companion}$ that of its
companion; for equal masses $R_L/d \approx 0.38$. For every binary of two stars the *filling
factor* $R/R_L$ is evaluated at periastron, $d = a(1-e)$ of the osculating orbit, as is usual
for eccentric orbits (Sepinsky et al. 2007). A star that orbits a sub-system – a giant tertiary
around an inner binary – has a lobe against the total mass of that sub-system (Eggleton 1983).
The largest filling factor of every star is stored; with `run.stop_on_roche_overflow` a run
stops at $R/R_L = 1$, and with `physics.mass_transfer` the transfer is followed instead
([Mass transfer and evolution](#roche-lobe-overflow)).

**Limits.** The formula assumes a circular synchronised binary; in an eccentric or asynchronous
orbit the lobe changes along the orbit and its value at periastron is approximate (Sepinsky et
al. 2007). For the shapes of distorted stars in light curves SWINGERsim solves the full Roche
potential instead ([Observables](#distorted-stars)).

## Kinds of systems {#kinds-of-systems}

The same machinery describes systems at very different stages of life. Table 4 is a guide to
the physics each kind needs; the chapters named in the last column describe it.

**Table 4.** Kinds of multiple stars and the physics they need. *Kind* describes the system;
*Model* the recommended `stellar_model.kind`; *Physics* the switches under `physics` (and per-star
keys) that matter; *Backend* the integration method of chapter [Dynamics](#dynamics); *Chapter*
where the physics is described.

| Kind | Model | Physics | Backend | Chapter |
|---|---|---|---|---|
| young, still in gas | `mist`, age < 1 Myr | `gas`, `tides`, `tidal_lag: hurley2002` | secular | Dynamics, Tides |
| Sun-like, main sequence | `zams` or `mist` | `gr`, `tides`, `rotation`; `magnetic_braking` | either | Tides, Mass transfer |
| close eclipsing binaries in a multiple | `mist` or `sse` | `gr`, `tides`, `rotation` | N-body for timings | Tides, Observables |
| giants, AGB | `mist` | `evolution`, `mass_loss`, `tides`, `tidal_lag: hurley2002` | secular | Mass transfer |
| Algols, giant tertiaries | `sse` or `mist` | `evolution`, `mass_transfer` | secular | Mass transfer |
| contact binaries, rotating stars | `mist` or `sse` | `mass_transfer.contact`; `shape`, `contact_fillout` | secular for evolution, N-body for light | Mass transfer, Observables |
| pulsating components | `mist` | `stars.<m>.pulsation`, `evolution` | N-body for light | Observables |
| massive stars, Wolf–Rayet | `mist`, `vvcrit: 0.4` | `evolution`, `mass_loss`, `tides` | secular | Mass transfer |
| very wide (> 1000 au) | any | `galactic_tide`, `encounters` | either | Dynamics |
| compact remnants | point masses | `gr`, `gw` | either | Dynamics |

A useful habit, borrowed from the PHOEBE documentation: switch the effects on one at a time and
compare the runs. Each switch shows what that piece of physics does to the system, and a result
that changes much with an effect thought negligible deserves a closer look.

## Implementation {#implementation-stars}

**External data.** The SSE formulae were evaluated with the published SSE Fortran routines of
Hurley, Pols & Tout (2000) once, offline, into a stored grid; the MIST v1.2 tracks (Choi et al.
2016; Dotter 2016) are converted from the published EEP archives; the structure constants are
the tables of Claret (2023) and Claret (2004); MESA (Jermyn et al. 2023) is an external program
run by the user for own grids. **Own code.** The ZAMS radius fit of Tout et al. (1996) and the
polytropes (Lane–Emden and Radau equations integrated with SciPy); the EEP and phase-fraction
interpolation of both track families; the $\log g$ interpolation of the Claret tables along
monotonic segments; the convective envelopes of MIST stars (giants, blue loop, pre-main
sequence) derived as described above; the search for initial masses from present masses; the
rotation correction of $k_2$; the reduction of MESA profiles ($k_2$ by Radau's equation,
$I/(MR^2)$, turnover time, inertial-wave coefficient); and the Eggleton (1983) lobes evaluated
at periastron for binaries and for stars orbiting sub-systems.
