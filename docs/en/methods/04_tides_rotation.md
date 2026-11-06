# Tides and rotation {#tides-and-rotation}

A star in a close binary is not a point. The companion pulls harder on the near side than on
the far side, so the star is stretched along the line of centres (a *tidal bulge*); if it
rotates, it is also flattened at the poles (a *rotational bulge*). Both distortions add a
quadrupole to the star's gravity, which turns the periastron (apsidal motion) and couples the
spin of the star to the orbit. Internal friction makes the tidal bulge lag behind the
companion; the lag transfers angular momentum between spin and orbit and dissipates energy, so
that spins synchronise, orbits circularise and shrink. In a multiple star tides have a second
role: the fast apsidal motion they cause in close binaries suppresses ZLK cycles
([Hierarchies and orbits](#zlk-cycles)).

![**Figure 8.** Tides, schematically. Left: a star (orange) distorted by its companion (blue);
the black line is the axis of the tidal bulge, which lags behind the line of centres (grey
dotted) because of internal friction. Right: evolution under tidal friction in units of the
tidal time scale (horizontal axis). Orange: spin frequency of the star divided by the orbital
frequency, rising from 0.3 to 1 (synchronisation); blue: orbital eccentricity, falling from 0.4
towards 0 (circularisation). Synchronisation is reached several times faster than
circularisation because a star stores much less angular momentum in its spin than the orbit
does.](../../figures/concepts/tides.png)

Figure 8 summarises the order of events. The right panel is illustrative (arbitrary time unit),
but its message holds for real binaries: the spin adjusts first, the orbit later. The
equations behind it are given below.

## Equilibrium tides in the N-body backend {#equilibrium-tides}

With `physics.tides: true` the stars become extended bodies. SWINGERsim uses the REBOUNDx effect
`tides_spin` (Lu et al. 2023), which implements the equilibrium-tide model of Eggleton,
Kiseleva & Hut (1998) for any number of bodies:

- every star is distorted by its own rotation and by the tidal field of **all** other stars,
  and the distorted stars attract each other with the corresponding quadrupole forces; this
  conservative part makes the periastra advance and the orbital planes precess about the
  stellar spins;
- the spin vector of each star evolves under the torques of all other stars (for example the
  precession of a spin tilted to its orbit);
- with a time lag $\Delta t > 0$ the tidal bulge lags behind, energy is dissipated, orbits
  circularise and shrink and the spins approach pseudo-synchronism.

The strength of the distortion is set by the potential Love number $2k_2$ and by $R^5$, the
spin evolution by the moment of inertia $I = r_g^2 M R^2$ ([Stars](#structure-constants)).
`physics.rotation` (which requires `physics.tides`) gives the stars a spin: synchronous with
their binary by default (`stars.<name>.P_rot: sync`), pseudo-synchronous (`pseudosync`), or a
given period in days, with an `obliquity` (spin–orbit angle) and a `spin_azimuth`. Without
`physics.rotation` the stars start without spin, which leaves only the tidal distortion (the
model used by Vokrouhlický 2016, his Fig. 7).

**Time-lag convention.** SWINGERsim expects the time lag `lag_s` (seconds) in the convention of
Mignard (1979), Hut (1981) and Leconte et al. (2010). Reading the source of `tides_spin` and
comparing its rates with these formulae showed that the same number given to REBOUNDx as its
parameter $\tau$ produces rates exactly twice as large; SWINGERsim therefore passes
$\tau = \Delta t/2$, and benchmarks against Leconte et al. (2010) check the result (see
*Testing, Part I*). For orientation, a constant time lag corresponds to a tidal quality factor
$Q \approx 1/(2\,|\Omega - n|\,\Delta t)$ at the tidal frequency, where $\Omega$ is the spin and
$n$ the orbital angular frequency.

**Conservation.** Without dissipation the energy including spin and tidal energy is conserved
to $|dE/E| \sim 10^{-13}$ (without relativity). With dissipation the energy decreases
physically, and the run checks the total angular momentum of orbits and spins instead.

## Classical apsidal motion of a binary {#classical-apsidal-motion}

For every binary of two stars SWINGERsim gives the apsidal advance expected from the binary
itself, the classical tidal and rotational term of Sterne (1939):

$$\frac{\dot\omega}{n} = \sum_{i} k_{2,i}\left(\frac{R_i}{a}\right)^5
\left[15\,\frac{m_j}{m_i}\,f(e) + \left(1 + \frac{m_j}{m_i}\right)
\left(\frac{\Omega_i}{n}\right)^2 g(e)\right],$$

$$ f(e) = \frac{1 + \tfrac32 e^2 + \tfrac18 e^4}{(1-e^2)^5}, \qquad g(e) = \frac{1}{(1-e^2)^2}, $$

where the sum runs over the two stars $i$ (with companion $j$), $k_{2,i}$ is the apsidal-motion
constant, $R_i$ the radius, $m_i$ the mass, $\Omega_i$ the spin angular frequency (aligned with
the orbit), $a$ the semi-major axis, $e$ the eccentricity and $n$ the mean motion. The first
term in the brackets is the tidal, the second the rotational contribution; the relativistic
term ([Hierarchies and orbits](#apsidal-motion)) is added separately. The steep dependence on
$(R/a)^5$ makes apsidal motion a probe of the radii and of the internal structure through $k_2$
(Claret & Giménez 2010).

In a multiple system the measured advance contains in addition the perturbation by the other
stars; the run summary gives both, and `swingersim apsidal` splits the total into tidal,
rotational, relativistic and dynamical parts and gives the factor by which $k_2$ of both stars
would have to change to reach a target rate. Published photodynamical solutions often fit the
apsidal rate as a free parameter instead of computing it from stellar structure; a difference
against SWINGERsim then measures the tension between the measured rate and stellar models, not a
difference between two codes.

**Pseudo-synchronisation.** In an eccentric orbit tides cannot synchronise a spin with the
orbit at every point; they drive it to the *pseudo-synchronous* rate at which the
orbit-averaged tidal torque vanishes (Hut 1981, eq. 42), faster than the mean motion
(1.24 times at $e = 0.2$). SWINGERsim uses it for `P_rot: pseudosync` and reports the
pseudo-synchronous period in the tidal time scales below.

**Rotation and $k_2$.** A rotating star is more centrally condensed than a static one. With
`physics.k2_rotation` the model constants are corrected by
$\delta\log k_2 = -0.87\,\lambda_s$, $\lambda_s = 2\Omega^2R^3/(3GM) = 2V^2/(3gR)$ (Claret 1999;
Claret & Giménez 2010 apply this rotational correction to observed binaries), with the angular velocity $\Omega$, equatorial velocity $V$, surface
gravity $g$, radius $R$ and mass $M$ of each star at that time, in all backends.

## Tidal evolution: the constant-time-lag model {#tidal-evolution}

The secular backends and the analytic time scales use the constant-time-lag theory (Mignard
1979; Hut 1981) in the form of Leconte et al. (2010, eqs. 1–3) for aligned spins. With the
potential Love number $2k_{2,i}$ and the time lag $\Delta t_i$ of star $i$,

$$ K_i = \tfrac32\,(2k_{2,i})\,\Delta t_i\,\frac{G m_j^2}{R_i}\left(\frac{R_i}{a}\right)^6 n^2 , $$

$$ \frac{da}{dt} = \frac{4a^2}{G m_1 m_2}\sum_i K_i\left[N(e)\frac{\Omega_i}{n} - N_{a1}(e)\right], \qquad
\frac{de}{dt} = \frac{11\,a\,e}{G m_1 m_2}\sum_i K_i\left[\Omega_e(e)\frac{\Omega_i}{n} - \frac{18}{11}N_{e1}(e)\right], $$

$$ \frac{d\Omega_i}{dt} = -\frac{2K_i}{C_i\,n}\left[\Omega(e)\frac{\Omega_i}{n} - N(e)\right], \qquad C_i = r_{g,i}^2 m_i R_i^2 , $$

where $N$, $N_{a1}$, $N_{e1}$, $\Omega_e$ and $\Omega(e)$ are the polynomial functions of the
eccentricity of Leconte et al. (2010) (the same functions appear in Hut 1981, eqs. 9–11), $C_i$
the moment of inertia and the other symbols as above. The spin equation follows from the balance
of the tidal torque with the change of the orbital angular momentum; for a circular orbit it
reduces to Murray & Dermott (1999, eq. 4.159). The secular backends integrate these equations
for every orbit of two single stars together with the gravitational dynamics.

For every binary the run summary lists the time scales at the start: circularisation
$e/|\dot e|$, orbital decay $a/|\dot a|$, the e-folding time of each spin towards
pseudo-synchronism and the pseudo-synchronous period. Realistic stellar lags give time scales of
millions of years and longer, beyond direct N-body integrations; in N-body runs dissipation is
mainly useful for strongly eccentric phases (ZLK peaks, where it shrinks orbits; Fabrycky &
Tremaine 2007) and for tests.

## Dissipation from the stellar structure {#dissipation-from-structure}

The time lag is the least known input of tidal evolution. With `physics.tidal_lag: hurley2002`
SWINGERsim derives it for every star of a two-star binary from its structure after Hurley, Tout
& Pols (2002, Sect. 2.3), instead of a fixed `lag_s`:

- a star with a **convective envelope** dissipates through turbulent viscosity (equilibrium
  tide; Zahn 1977; Rasio et al. 1996),
  $$ \left(\frac{k}{T}\right)_\mathrm{c} = \frac{2}{21}\,\frac{f_\mathrm{conv}}{\tau_\mathrm{conv}}\,\frac{M_\mathrm{env}}{M}, \qquad f_\mathrm{conv} = \min\left[1, \left(\frac{P_\mathrm{tid}}{2\tau_\mathrm{conv}}\right)^2\right], $$
  with the envelope mass $M_\mathrm{env}$, the eddy turnover time
  $\tau_\mathrm{conv} = 0.4311\,[M_\mathrm{env} R_\mathrm{env}(R - R_\mathrm{env}/2)/(3L)]^{1/3}$ yr
  ($R_\mathrm{env}$ the depth of the envelope, $L$ the luminosity, solar units) and the tidal
  pumping period $P_\mathrm{tid}$, $1/P_\mathrm{tid} = |1/P_\mathrm{orb} - 1/P_\mathrm{spin}|$
  (Hurley et al. 2002, eq. 33; the factor 2 in $f_\mathrm{conv}$, their eq. 32, is there because
  the tide comes around twice in this period). The factor $f_\mathrm{conv}$ reduces the viscosity when the tidal period is shorter than the turnover time
  – quadratically (Goldreich & Nicholson 1977; Hurley et al. 2002) by default, or linearly
  (Zahn 1966, 1989; Zahn & Bouchet 1989) with `physics.tidal_convective_reduction: linear`;
- a star with a **radiative envelope** dissipates through the dynamical tide (Zahn 1975, 1977),
  $(k/T)_\mathrm{r} = 1.9782\times10^4\,(M R^2/a^5)(1+q)^{5/6} E_2$ yr⁻¹ with
  $E_2 = 1.592\times10^{-9} M^{2.84}$ ($q$ the mass ratio, solar units).

Hut's (1981) $(k/T)$ corresponds to the constant time lag through
$T = R^3/(G M \tau_\mathrm{Hut})$ and $\tau_\mathrm{Hut} = 2\,\Delta t$ (SWINGERsim uses the Love
number $2k_2$ of Leconte et al. 2010; the factor was checked numerically). The lag is computed
from the initial orbit and spin and, for evolving stars, again in every segment; a `lag_s` given
for a star is kept, and `physics.tidal_efficiency` multiplies the dissipation.

**Calibration.** This recipe is widely used (BSE, MSE, COSMIC), but its calibration is the weakest
link of long-term evolution. On the main sequence alone, solar-type binaries circularise in
4 Gyr only up to a period of about 3.8 d, whereas the cluster M67 shows 12.1 +1.0/−1.5 d (Meibom &
Mathieu 2005, Table 3); a `tidal_efficiency` of about 150 would be needed on the main sequence
alone (100 gives 11.2 d, 200 gives 12.9 d).

### Inertial waves {#inertial-waves}

In a rotating convective envelope a tidal forcing with frequency below $2\Omega$ excites
*inertial waves*, restored by the Coriolis force. Their dissipation – the dynamical tide of
convective zones – is much stronger than the damping of the equilibrium tide by turbulent
viscosity, which is weak because the tidal period is short compared with the convective turnover
time (Ogilvie 2013; Mathis 2015; Barker 2020, 2022). Averaged over the tidal frequency it
depends only on the structure and the spin,

$$ \frac{1}{\langle Q'_\mathrm{IW}\rangle} = \epsilon^2\,C, \qquad \epsilon^2 = \frac{\Omega^2 R^3}{GM}, $$

where $Q'$ is the modified tidal quality factor, $\Omega$ the spin rate and $C$ a structure
coefficient. With `physics.inertial_waves: true` (requires `tidal_lag: hurley2002`) SWINGERsim
computes $C$ for a star with a radiative core of radius $\alpha R$ and mass $\beta M$ (from the
envelope model) with the two-layer formula of Ogilvie (2013), and for a fully convective star
(below 0.35 M☉) for an $n = 1.5$ polytrope, $C = 1/130.9$, computed once with the method of Barker
(2020, eqs. 12–38; 2022, App. A) for continuous profiles. The eccentricity damping
$d\ln e/dt = -(225/16)(M_2/M)(R/a)^5\,n/Q'$ (Barker 2022, eq. 1) is converted into an equivalent
constant time lag by matching the damping rate of the lag model for a synchronous star,
$-(21/2)(k/T)\,q(1+q)(R/a)^8$ (Hut 1981), and added to the equilibrium tide. Inertial waves are
switched off when the star rotates slower than half the orbital frequency (no forcing frequency
below $2\Omega$).

The numerical method reproduces the two-layer formula in the form of Barker (2020, eq. 39) for a
smoothed two-layer profile to 0.2 %, a check of the solver and its normalisation. For
main-sequence stars of 0.3 and 1 M☉ rotating in 10 d it gives $Q' = 0.9$ and
$1.0\times10^7$, as Barker (2022) finds with MESA profiles. For realistic profiles of solar-type
stars the two-layer model under-predicts the dissipation by about a factor of two, and for F stars
above 1.2 M☉ by an order of magnitude or more (Barker 2020).

### Tides before the main sequence {#pms-tides}

Young stars are large and fully convective, and most of the circularisation of solar-type
binaries happens before they reach the main sequence (Zahn & Bouchet 1989). With MIST stars on
the pre-main sequence ([Stars](#mist-model)) the equilibrium tide and the inertial waves act
there too. Twin 1 M☉ binaries started at 1 Myr are circular up to an initial period of 7.9 d at
50 Myr (Zahn & Bouchet: 7.2–8.5 d). In the period observers use – the final period of the
binaries just circularised – they are circular up to 6.6 d at 3 Myr (pre-main-sequence binaries:
7.1 ± 1.2 d; Meibom & Mathieu 2005, Table 3) and up to 7.5 d at 4 Gyr (M67: 12.1 +1.0/−1.5 d).
With inertial waves the main sequence alone circularises solar-type twins up to 7.5 d in 4 Gyr
(3.8 d with the equilibrium tide alone). The pre-main sequence is thus reproduced; the further
growth on the main sequence is not. Barker (2022), with frequency-averaged inertial-wave
dissipation computed from MESA profiles and synchronised spins, finds circularisation periods
that grow with main-sequence age as observed. SWINGERsim's two-layer estimate with SSE and MIST
structures and evolving spins does not reproduce this; the difference lies in the stellar
structure (the two-layer model under-predicts the dissipation of realistic solar-type profiles,
see above) and in the spin history, and it remains an open point. Resonance locking of gravity
modes, which SWINGERsim does not model, acts mainly on the pre-main sequence: Zanazzi & Wu (2021)
find circular periods of about 4–6 d from it and a negligible main-sequence contribution. The result also
depends on the start of the calculation: a younger, larger star circularises more (started at
0.15 Myr, 11.6 d already at 3 Myr). On the pre-main sequence the inertial waves dominate, so the
linear and quadratic reductions of the turbulent viscosity agree there within 1 %.

## Limits {#tides-limits}

- The equilibrium-tide model assumes small, quasi-static distortions. Dynamical tides of
  radiative stars enter only through the Zahn recipe in the lag; resonance locking and
  non-linear tides are not modelled.
- In the secular backends spins stay aligned with their orbits; spin–orbit misalignment and spin
  precession are followed only by the N-body backend.
- The dissipation rates carry the calibration uncertainty discussed above; the circularisation
  periods of old clusters are not reproduced by any equilibrium-tide model, including this one.
- Without `physics.tides` the eccentricity growth in ZLK cycles of close binaries is
  overestimated (Vokrouhlický 2016).

## Implementation {#implementation-tides}

**External library.** The tidal and rotational forces and the spin evolution of the N-body
backend are the REBOUNDx effect `tides_spin` (Lu et al. 2023; Tamayo et al. 2020).
**Own code.** The conversion of the time-lag convention ($\tau = \Delta t/2$, found by reading
the REBOUNDx source and verified against Leconte et al. 2010); the initial spins (synchronous,
pseudo-synchronous after Hut 1981, given periods and obliquities); the Sterne (1939) apsidal
rates and the split into contributions; the constant-time-lag equations of Leconte et al. (2010)
in the secular backends and the tidal time scales; the structure-based lag of Hurley, Tout &
Pols (2002) with both viscosity reductions; the rotation correction of $k_2$; and the
inertial-wave dissipation – the two-layer formula of Ogilvie (2013) and an own solver of the
continuous-profile method of Barker (2020, 2022) for polytropes and MESA profiles.
