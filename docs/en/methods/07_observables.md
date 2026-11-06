# Observables {#observables}

Nobody watches a multiple star for a million years. Astronomers have snapshots: the brightness
of the system over nights to years, the times at which its eclipses happen, velocities from
spectra. This chapter describes how SWINGERsim turns the orbits and stars of the previous
chapters into these measurements: eclipses and their detectability, light curves, radial
velocities, times of minima and O−C diagrams, pulsating stars as clocks, and dust discs as
occulters. All observables are computed from the N-body backend, because they need the position
of every star on its orbit.

## Eclipse geometry {#eclipse-geometry}

A binary eclipses when, seen along the line of sight, the disk of one star passes over the disk
of the other. For spherical stars of radii $R_1$, $R_2$ on a circular orbit of separation $a$
the condition is

$$ |\cos i| < \frac{R_1 + R_2}{a} , $$

where $i$ is the inclination to the sky ([Hierarchies and orbits](#orbital-elements)). For an
eccentric orbit $a$ is replaced by the separation at conjunction. For an isotropic orientation
the probability that a binary eclipses is $(R_1+R_2)/a$; for the two binaries of V994 Her this
gives 37 % and 40 %, as found by Vokrouhlický (2016) (see *Testing, Part I*).

![**Figure 12.** An eclipsing binary. Left and middle: the binary seen on the sky (the observer
looks along $z$); the grey ellipse is the projected orbit of the blue star around the orange one.
At $i = 88°$ (left) the projected orbit is nearly a line and the disks overlap at conjunction; at
$i = 60°$ (middle) they pass each other without overlapping. Right: the light curve of the
eclipsing case, brightness (normalised to 1 outside eclipse) against orbital phase. The deeper
minimum at phase 0 is the primary eclipse, when the hotter star is hidden; the shallower one at
phase 0.5 the secondary eclipse.](../../figures/concepts/eclipse_geometry.png)

Figure 12 shows the two cases. The depth of each minimum is the fraction of the total light that
is hidden, which depends on the overlap of the disks and on the surface brightness of the hidden
star: hiding the hotter star costs more light, so the primary minimum is deeper.

During a run SWINGERsim reports for every binary of two stars the fraction of the time in which
the condition above holds (spherical stars, circular-orbit form, as in Vokrouhlický 2016, Sect.
3.3), with the radii changing as the stars evolve; for 2+2 systems it also gives the fraction in
which both binaries eclipse at once – the defining property of doubly eclipsing systems.
Ensembles report whether a binary eclipses at $t = 0$ and when its eclipses start or end.

### Geometric eclipses and detectable eclipses {#detectable-eclipses}

The geometric condition says only that the disks overlap. Whether an observer *detects* the
eclipse depends on its depth in the total light. In a multiple star the light of all other stars
dilutes the eclipses of a binary. If a binary alone would lose the fraction $d$ of its own light
at minimum, and the other stars contribute the fraction $l_3$ of the total light (the *third
light*), the observed depth in magnitudes is

$$ \Delta m = -2.5\,\log_{10}\left[1 - (1 - l_3)\,d\right] . $$

A grazing eclipse ($d$ small) in a bright multiple ($l_3$ close to 1) can therefore be far below
the photometric precision even though the geometric condition holds.

![**Figure 13.** Eclipse limit and detectability (new figure for this volume, computed from the
two equations above). (a) The largest tilt from edge-on, $|90° - i|_\mathrm{max} = \arcsin[(R_1 + R_2)/a]$,
for which a circular binary still eclipses, against its period $P$ in days (both axes logarithmic).
Blue: two 1 M☉ ZAMS stars of 0.89 R☉; orange: two 2 M☉ stars of 1.70 R☉; green: two 1 M☉
subgiants of 5 R☉ (shown only where the stars are well detached). (b) Observed depth $\Delta m$
in magnitudes against the intrinsic fractional depth $d$ of the binary (both logarithmic), for
third light $l_3$ = 0 (blue), 0.5 (orange), 0.9 (green) and 0.99 (yellow). The dashed line marks
an example detection limit of 1 mmag.](../../figures/methods/eclipse_detectability.png)

Figure 13 separates the two questions. Panel (a): a pair of Sun-like stars with a period of 1 d
eclipses within about 20° of edge-on, at 10 d within about 4°, at 100 d within less than 1°;
larger stars widen the window in proportion to $R_1 + R_2$. This is why long-period binaries are
rarely eclipsing and why a precessing binary gains and loses its eclipses
([Hierarchies and orbits](#nodal-precession)). Panel (b): for small $d$ the observed depth is
$\Delta m \approx 1.086\,(1 - l_3)\,d$, so each line is the blue one shifted down by the factor
$1 - l_3$. An intrinsic depth of 1 % appears as 11 mmag alone, 5 mmag with $l_3$ = 0.5, 1.1 mmag with
$l_3$ = 0.9 and 0.11 mmag with $l_3$ = 0.99 – the last one below the example detection limit. When
SWINGERsim is used to predict when eclipses become visible or disappear, the synthetic light curve
(below), not the geometric condition, decides.

## Light curves {#light-curves}

`swingersim lightcurve` computes the light an observer records, from the N-body positions of all
stars. Three levels of fidelity are available (Table 11).

**Table 11.** Light-curve models. *Model* names the level; *Stars* the shapes; *Physics* what is
included; *Use* where it is appropriate. All levels share the positions, the light-travel time and
the normalisation (1 = all stars outside eclipse).

| Model | Stars | Physics | Use |
|---|---|---|---|
| spheres (default) | limb-darkened spheres | exact occultations of any number of stars, analytic ellipsoidal, reflection and beaming terms | detached systems, $R/R_L$ well below 1 |
| meshes | Roche surfaces, rotating stars, contact binaries; all stars meshes | gravity and limb darkening per surface element, occultation along every ray | close, semi-detached and contact binaries, fast rotators |
| PHOEBE hybrid | PHOEBE 2 for each inner binary | Roche geometry, model atmospheres, irradiation of PHOEBE; mutual eclipses of different binaries as spheres | selected solutions, cross-checks, training sets |

### Spheres: intensity, limb darkening and passbands {#limb-darkening}

The surface brightness of a star in a passband is that of a black body integrated over the
passband with photon-counting weights, $w = \int B_\lambda(T_\mathrm{eff})\,S(\lambda)\,\lambda\,d\lambda$,
where $B_\lambda$ is the Planck function and $S(\lambda)$ the passband response. Available
passbands are Johnson–Cousins U, B, V, R, I (Bessell 1990), TESS (Ricker et al. 2015) and Gaia G,
BP, RP (Riello et al. 2021). Model atmospheres change the surface-brightness ratio of stars of
different temperature by a few per cent; the PHOEBE hybrid uses model atmospheres where this
matters.

The intensity across the disk follows the four-coefficient limb-darkening law of Claret (2000),

$$ \frac{I(\mu)}{I(1)} = 1 - \sum_{k=1}^{4} a_k\left(1-\mu^{k/2}\right), $$

where $\mu$ is the cosine of the angle between the line of sight and the surface normal and $a_k$
are the coefficients; linear and quadratic laws are special cases, and coefficients can be given
per star. Below 15 000 K the coefficients and the gravity-darkening exponents come from Claret &
Bloemen (2011), interpolated in $T_\mathrm{eff}$, $\log g$ and [M/H]; for TESS and Gaia the nearest
Johnson–Cousins band of those tables is used (TESS and RP: I; G: R; BP: V). From 15 000 to 55 000 K
the coefficients come from the non-LTE TLUSTY atmospheres (OSTAR2002, BSTAR2006) of Reeve & Howarth
(2016); where both grids apply, the disk-averaged intensities agree to 0.2–0.7 %. Beyond the
gravity-darkening grid of Claret & Bloemen (O stars above about 39 000 K) the exponent follows von
Zeipel's law in the band, $y = \tfrac14\,x e^x/(e^x - 1)$, $x = hc/(\lambda k T_\mathrm{eff})$. Stars
outside the tables (for example cool AGB stars) take the nearest tabulated node, with a warning.

**Occultations.** Every star is occulted by all stars nearer to the observer whose disks overlap
it. The hidden flux is $\int_0^R I(r)\,\alpha(r)\,2\pi r\,dr$, where $\alpha(r)$ is the fraction of
the circle of radius $r$ around the centre of the star covered by the union of the occulting disks;
the radial integral is computed with Gauss–Legendre quadrature between the radii at which a disk
edge crosses (as in ellc, Maxted 2016). This is exact for any number of spheres, so mutual eclipses
of the two binaries of a 2+2 system and triple overlaps need no special treatment. It agrees with
the analytic model of Mandel & Agol (2002), as implemented in batman (Kreidberg 2015), to
$10^{-8}$.

**Proximity effects** of the stars of two-star binaries, in the small-distortion approximation:
the ellipsoidal variation $-\alpha_\mathrm{ell}\,q\,(R/d)^3(2\cos^2\psi - \sin^2 i)$ with
$\alpha_\mathrm{ell} = 0.15(15+u)(1+\tau)/(3-u)$ (Morris 1985; Morris & Naftilan 1993), where $q$ is
the mass ratio, $d$ the separation, $\psi$ the angle between the line of centres and the line of
sight, $u$ the linear limb-darkening and $\tau$ the gravity-darkening coefficient; Doppler beaming
$-(3-\alpha_\nu)\,v_z/c$ with the spectral index $\alpha_\nu$ (Loeb & Gaudi 2003); and reflection of
the companion's light by a Lambert sphere with Bond albedo 0.5 for convective and 1 for radiative
envelopes (Rucinski 1969). The ellipsoidal term underestimates the variation by 3, 7, 11 and 17 % at
$R/a$ = 0.10, 0.15, 0.20 and 0.25 compared with the Roche model of ellc; closer pairs need meshes.

**Exposure and noise.** Finite exposures are modelled by averaging the model over sub-exposures
(exposure time and number of sub-points per exposure); Gaussian white noise can be added.

![**Figure 14.** Synthetic light curve of a 2+2 system computed with spheres. Top: total flux
(normalised to 1 outside eclipses) against time in days. Bottom: the flux of each star divided by
its own unocculted flux, for m1 (dark blue), m2 (light blue), m3 (orange) and m4 (light orange).
The eclipses of binary A (m1 and m2 hidden in turn) appear in the total light diluted by the light
of B; the small slow wiggles are the ellipsoidal and reflection variations.](../../figures/concepts/light_curve.png)

Figure 14 shows how the signal of a multiple star is built. The bottom panel separates the stars:
m2 loses about 12 % of its own light at its eclipse, but in the total light (top) this becomes a
3 % dip, because the other three stars keep shining. This dilution is the third light of the
previous section, and it is why a model of the whole system, not of one binary, is needed to
interpret eclipse depths.

### Light-travel time {#light-travel-time}

Every star is seen at its own *emission time*: the light observed at time $t$ left a star at depth
$z$ (positive away from the observer, barycentre at $z = 0$) at $t_e = t - z/c$. Every two-star
binary is integrated separately to the emission time of its centre of mass,
$t_e = t - z_\mathrm{com}(t_e)/c$ (solved iteratively), and only the offsets of its two stars from
the centre are corrected to first order, $-\mathbf{v}\,\Delta z/c$ (the Rømer delay inside the
binary); other stars use $\mathbf{x}(t) - \mathbf{v}\,z/c$. A first-order correction of the whole
delay is not enough when $z/c$ is a sizeable fraction of the binary period: for QZ Car
($z/c \approx 0.2$ d, $P$ = 6 d) it would shift the eclipses by about 200 s.

### Distorted stars {#distorted-stars}

Stars can be given a Roche shape (`stars.<m>.shape: roche`, both stars of a circular synchronised
binary), the oblate shape of a rigidly rotating star (`shape: rotating`), or an orbit can be a
contact binary with fill-out factor $f$ (`orbits.<X>.contact_fillout`). The binary model uses the
Roche potential (Kopal 1959; Wilson & Devinney 1971)

$$ \Omega(x,y,z) = \frac{1}{r_1} + q\left(\frac{1}{r_2} - x\right) + \frac{1+q}{2}\,(x^2+y^2) , $$

in units of the separation, with star 1 at the origin, star 2 at $(1, 0, 0)$, $r_1$, $r_2$ the
distances to the two stars, $q = M_2/M_1$ and $z$ along the orbit normal. The inner critical
surface $\Omega_\mathrm{in} = \Omega(\mathrm{L1})$ bounds the Roche lobes, the outer one
$\Omega_\mathrm{out}$ (through L2 or L3) the common envelope; the fill-out factor is
$f = (\Omega_\mathrm{in}-\Omega)/(\Omega_\mathrm{in}-\Omega_\mathrm{out})$, negative for detached
and between 0 and 1 for contact binaries (Lucy & Wilson 1979; Rucinski 1973). A rotating star
follows the Roche model $1/r + (4/27)\,w^2 r^2\sin^2\theta = 1$ in units of its polar radius, with
$w = \Omega/\Omega_\mathrm{crit}$ and $\theta$ the colatitude (Maeder 2009); its equator reaches 1.5
times the polar radius at $w = 1$. The configured radius is the volume-equivalent radius.

Each surface is sampled on 64 × 128 directions. An element emits
$I = I_\mathrm{ref}\,(g/\bar g)^y\,L(\mu)$, with the local effective gravity $g$, its surface mean
$\bar g$, the band gravity-darkening exponent $y$ (von Zeipel 1924; Claret & Bloemen 2011) and the
limb darkening $L(\mu)$; a contact binary's two lobes share one effective temperature. An element is
hidden when its normal points away from the observer, when a sphere lies in front of it
(anti-aliased by its projected size), or when its line of sight passes through another body or the
other lobe of its own Roche body. When any star is a mesh, all stars are meshes. Accuracy: $3\times10^{-4}$
of the flux against the exact sphere eclipses; the amplitude of a contact binary agrees with PHOEBE
2.5 to 0.004 mag (see *Testing, Part I*). The error of spheres against Roche surfaces grows steeply
with the filling of the lobe: about 0.5 ppt rms at $R/R_L$ = 0.3, 2 ppt at 0.5, 5 ppt at 0.7, 9 ppt
at 0.85 and 17 ppt at 0.95, and 15–45 ppt for contact binaries. The Roche model assumes a circular
synchronised orbit; eccentric distorted binaries need the PHOEBE hybrid.

### The PHOEBE hybrid {#phoebe}

With the PHOEBE engine the light of every inner binary of two stars comes from PHOEBE 2 (Prša et
al. 2016; Horvat et al. 2018; Conroy et al. 2020): Roche geometry, gravity darkening, irradiation,
Castelli–Kurucz model atmospheres or black bodies, and the light-travel time inside the binary.
The multiple system comes from SWINGERsim:

1. in chunks of one binary period the N-body integration gives the osculating elements of the
   binary at the chunk centre and the motion of its centre of mass along the line of sight; period
   and time of periastron come from the mean longitude over one orbit, because the osculating
   Keplerian period of a tidally distorted binary is wrong by minutes per orbit;
2. PHOEBE computes the binary through the chunk at the emission times;
3. stars outside two-star binaries shine with their passband luminosity;
4. eclipses between stars of different binaries are taken from the sphere model and subtracted.

The conventions are converted ($z$ away from the observer; PHOEBE's argument of periastron is that
of the primary, `per0` = $\omega$ + 180°). With the same physics the hybrid reproduces the sphere
model to a few $10^{-4}$. PHOEBE needs 0.1–0.2 s per point and binary, about $10^4$ times more
than the sphere model, and is meant for selected solutions. Rotation is synchronous, Doppler
beaming is absent in PHOEBE 2.5, and discs are ignored. The native treatment of multiple systems by
Brož et al. (2025) is an alternative when all light must come from one code.

## Radial velocities {#radial-velocities}

The radial velocity of star $k$ is its velocity along the line of sight relative to the barycentre
plus the systemic velocity $\gamma$: $\mathrm{RV}_k = \gamma + v_{z,k}$, positive away from the
observer. A velocity measured at time $t_\mathrm{obs}$ belongs to the emission time $t_e$ with
$t_e + z_k(t_e)/c = t_\mathrm{obs}$; SWINGERsim finds it by Newton steps until the mismatch is below
1 ms. The light time of the star (up to tens of minutes in a wide outer orbit) is therefore included,
which matters for the steep velocity curves of close binaries (a semi-amplitude of 100 km s⁻¹ in a
3-day orbit changes by 0.15 km s⁻¹ per minute).

**Centre of light and the Rossiter–McLaughlin effect.** A spectrograph measures the Doppler shift
of the light of a star – the flux-weighted mean line-of-sight velocity of its visible surface – not
the velocity of its centre of mass. For spherical, slowly rotating stars the two agree. A tidally
distorted star is brighter and larger on some sides and its synchronous rotation adds surface
velocity, so its centre of light moves on a different orbit (tens of km s⁻¹ in contact binaries;
Wilson & Sofia 1976); during an eclipse the covered part of a rotating star is missing from the
mean, which shifts the measured velocity (the Rossiter–McLaughlin effect; Rossiter 1924; McLaughlin
1924). With `fit.rv.centre: light` the model velocity is the flux-weighted line-of-sight velocity
of the visible surface elements of the meshes above – Roche stars corotate with their orbit, other
stars rotate with their spin – in the passband `fit.rv.band`. The offset from the centre-of-mass
velocity is added to the exact centre-of-mass velocity, so the light-travel time and the dynamics
stay those of the N-body model. Irradiation and beaming are not included in the weights. For a
contact binary ($q$ = 0.42, fill-out 0.3) the offsets reach ±34 and ±52 km s⁻¹ and agree with
PHOEBE 2 to 2 % of their amplitude; a fit of such velocities with the centre of mass gives masses
too large by several per cent.

## Times of minima and O−C diagrams {#times-of-minima}

If a binary were alone, its eclipses would follow a linear *ephemeris*, $T_n = T_0 + nP$, with the
cycle number $n$. Observers compare each measured time of minimum (O) with the time computed from
such an ephemeris (C) and plot the difference O−C against time. A straight line means a wrong
period; anything else means that something moves the binary or changes its orbit.

**Synthetic timings.** `swingersim timings` integrates the system (N-body with the configured
physics) and finds every eclipse of each inner binary: the minima of the projected separation of
its two stars, sampled 40 times per binary period and refined with the two-body motion between
samples. A conjunction is an eclipse when the projected separation is smaller than $R_1 + R_2$ (the
ratio is reported as the impact parameter). A *primary* eclipse is an eclipse of the binary's first
star (m1 of A, m3 of B), i.e. the first star is the farther one. Observed times include the
light-travel time of the binary's centre of mass, $t_\mathrm{obs} = t + z_\mathrm{com}/c$; the
dynamical delays come from the integration itself. A linear ephemeris fitted to the primary
eclipses gives the O−C; secondary eclipses use half-integer epochs.

### Light-travel-time effect {#ltte}

The binary orbits the barycentre of the whole system, so it is sometimes farther from the observer
and sometimes closer. Its eclipses arrive late and early in turn: the *light-travel-time effect*
(LTTE). For an outer orbit of period $P_3$, eccentricity $e_3$, argument of periastron $\omega_3$
and true anomaly $\nu_3$, the delay is (Irwin 1952)

$$\Delta_\mathrm{LTTE} = \frac{a_{12}\sin i_3}{c}\left[\frac{1-e_3^2}{1+e_3\cos\nu_3}\sin(\nu_3+\omega_3) + e_3\sin\omega_3\right],$$

where $a_{12}\sin i_3$ is the projected semi-major axis of the binary's orbit around the barycentre.
Its semi-amplitude is $K = (a_{12}\sin i_3/c)\sqrt{1-e_3^2\cos^2\omega_3}$ and the mass function
$f(m_3) = (a_{12}\sin i_3)^3/P_3^2$ (au, yr, M☉). The N-body timings contain this term exactly,
together with every other effect.

![**Figure 15.** The light-travel-time effect. Left: a binary (blue dot) on its outer orbit (green
ellipse) around the barycentre (cross); the arrow points to the observer. On the far side of the
orbit its eclipses are seen later, on the near side earlier. Right: the resulting O−C in minutes
against time in days, a periodic curve with the outer period and an amplitude equal to the light
time across the binary's projected orbit (here about 10 min; illustrative).](../../figures/concepts/lite_oc.png)

Figure 15 links the geometry to the diagram. The O−C curve repeats with the outer period; its
semi-amplitude is $a_{12}\sin i_3/c$ (about 8.3 min per au), and its shape departs from a sine for
an eccentric outer orbit. In a 2+2 quadruple the two binaries orbit each other, so their LTTE curves
are in anti-phase with amplitudes in the inverse ratio of their masses; in compact systems, with
outer periods of a year or less, the dynamical perturbations of the period add a part that is not
in anti-phase and can be larger than the LTTE.

### Apsidal motion in the O−C diagram {#apsidal-motion-in-o-c}

In an eccentric binary the time between primary and secondary minimum depends on the orientation
of the orbit. When the periastron turns (apsidal motion), the primary and secondary minima move in
opposite directions in the O−C diagram, with the apsidal period $U$.

![**Figure 16.** Apsidal motion (illustrative). Left: an eccentric orbit (blue ellipses, shades for
successive times) around the orange star; the dots mark the periastron, which slowly turns. Right:
O−C in minutes against time in years for the primary (blue) and secondary (orange) minima: two
curves in anti-phase with the apsidal period, here about 25 yr.](../../figures/concepts/apsidal.png)

Figure 16 shows the signature used to measure $\dot\omega$. The amplitude of the two curves grows
with the eccentricity, and their anti-phase distinguishes apsidal motion from the LTTE, which shifts
both minima together. The analytic form used in the stage-1 fit is given in chapter
[Fitting](#stage-1).

### Phenomenological minima {#phenomenological-minima}

Observers do not fit physics to an eclipse; they measure its time of minimum, depth and width with a
phenomenological profile. `swingersim minima` reduces the synthetic light curve in the same way,
with the profile of Mikulášek (2015, eq. 13),

$$\Delta m(t) = A\,(1 + C x^2)\left\{1 - \left[1 - e^{1 - \cosh x}\right]^{\Gamma}\right\},
\qquad x = \frac{t - t_\mathrm{min}}{D},$$

with central depth $A$ (mag), width $D$, kurtosis $\Gamma$ and correction $C$, on a quadratic
baseline (the light of the other stars, the ellipsoidal and reflection variation, which is
symmetric about the minimum to second order, and a slow drift). The eclipses are located with the
synthetic timings; the window around each is 2.5 times its duration (at most half the period). A
*shape* ($A$, $D$, $\Gamma$, $C$) is fitted per binary and kind (primary, secondary) to all its
eclipses together, optionally per season (which follows the slow change of depth and width as the
orbit precesses); then the time of every eclipse is fitted with that shape and free depth, width
and baseline, its error from the covariance matrix scaled by the residuals. By default every
binary is measured on its own light plus the mean light of the others (ideal disentangling); with
the blended option the total light is used as observed, eclipses of other binaries in the window
are masked, and their proximity variation shifts the minima (in a CzeV4315-like example the scatter
of B's minima grows from 0.2 s to 8 s).

The difference between a minimum and the barycentric conjunction is physical. Light-travel time
inside the binary delays the eclipse of the farther star: the secondary minimum of a circular orbit
comes $(P/\pi c)(K_2 - K_1) = 2a\sin i\,(m_1 - m_2)/(Mc)$ after half the period (Kaplan 2010; $K_1$,
$K_2$ the velocity semi-amplitudes), reproduced to 0.4 %. For a noise-free curve the Mikulášek
profile fits detached eclipses to about 1 % of their depth.

## Pulsating stars as clocks {#pulsating-stars}

Many components of multiple systems pulsate: δ Scuti and γ Doradus stars near 1.5–2.5 M☉, β Cephei
and slowly pulsating B (SPB) stars, RR Lyrae stars, Cepheids and long-period variables on the AGB.
With `stars.<m>.pulsation` a star's own flux varies as

$$\Delta m(t) = \sum_k A_k \sin\left[2\pi\left(\frac{t_e}{P_k} + \phi_k\right)\right],$$

with amplitudes $A_k$, periods $P_k$ and phases $\phi_k$ of the modes, evaluated at the star's
*emission time* $t_e$. An orbiting pulsator is therefore a clock: its phases lag by the light time
across its orbits, the phase modulation of Murphy et al. (2014). For a circular orbit the delays are
$\tau(t) = (a_1\sin i/c)\sin(2\pi t/P_\mathrm{orb} + \varphi)$ ($a_1$ the semi-major axis of the
pulsator's orbit around the barycentre); in a multiple system the delays of the inner orbit add to
those of the outer one. SWINGERsim reproduces delays measured from a light curve sampled like Kepler
data to $3\times10^{-7}$.

![**Figure 17.** A pulsating component as a clock (illustrative). Top: brightness of the pulsator
(arbitrary units) against time in hours, on the near side of its orbit (grey) and half an orbit
later (orange): the same pulses arrive later. Bottom: the delay of the pulsation phase in minutes
against time in years (orange line) and the O−C of the eclipses of the same binary (blue dots):
both trace the same orbital light-time curve.](../../figures/concepts/pulsator_clock.png)

Figure 17 shows why a pulsator in an eclipsing binary carries two clocks that must agree: the
eclipses sample the light-time curve once or twice per binary period, the pulsations continuously.

Periods follow from the stellar model (current mass, radius, luminosity, $T_\mathrm{eff}$) unless
`period_d` is given (Table 12).

**Table 12.** Pulsation types and the source of their periods. *Type* is the value of
`pulsation.type`; *Period relation* the formula used; *Source* the reference and its range.

| Type | Period relation | Source |
|---|---|---|
| `lpv` | non-linear fundamental-mode period–mass–radius relation, broken power law with saturation | Trabucchi et al. (2021), eqs. 1–2; O-rich AGB, 0.6–7 M☉ |
| `delta_scuti` | $Q = P\sqrt{\bar\rho/\bar\rho_\odot}$ = 0.0329 d (F), 0.0251, 0.0203, 0.0172 d (1st–3rd overtone) | Breger (1979); consistent with Poro et al. (2024) |
| `rr_lyrae` | $\log P_0 = 11.627 + 0.823\log L - 3.506\log T_\mathrm{eff} - 0.582\log M$; $P_1/P_0 = 0.744$ | Bono et al. (1996) as quoted by Alcock et al. (1997) |
| `cepheid` | $\log P = a + b\log T_\mathrm{eff} + c\log M + d\log L$ (F), first overtone separately | Bono, Castellani & Marconi (2000), Z = 0.004–0.02 |
| `beta_cep` | radial fundamental mode, $Q$ = 0.033 d | Handler (2013) |
| `spb`, `gamma_dor` | series $P_0 + k\,\Delta P$ of gravity modes; $P_0$, $\Delta P$ given by the user | Wu, Li & Deng (2018) |
| `manual` | `period_d` or $Q$ given | – |

Every type warns when the star's $T_\mathrm{eff}$ lies outside its class (a plausibility check).
`swingersim pulsation` also reports $dP/dt$ from the evolving stellar model, the evolutionary period
change that O−C diagrams of pulsators measure. MIST stars in the middle of the δ Scuti strip lie on
the observed period–luminosity relation of Poro et al. (2024) to 0.1 mag, and MIST Cepheids on their
blue loops match the period–radius relation of Anderson et al. (2016) to 0.01 dex. Amplitudes are
fixed per band, mode shapes are sinusoids (no combination frequencies or non-sinusoidal RR Lyrae and
Mira shapes), and pulsation does not change the radius used for eclipses.

## Discs as occulters {#discs}

Some eclipses are caused by a dark disc of dust: the two-year eclipses of ε Aur by a disc around its
companion, the eclipses of EE Cep whose depth changes from cycle to cycle as its disc precesses, or
KH 15D, a young binary disappearing behind the edge of a precessing circumbinary ring. A disc in
SWINGERsim is a slab: a cylinder of outer radius $r_\mathrm{out}$, optional hole $r_\mathrm{in}$ and
full thickness $H$, centred on a star or on the centre of mass of an orbit. For every line of sight
from a point of a star towards the observer the length $L$ inside the slab is computed exactly, and
the transmission is $e^{-\tau L/H}$, with $\tau$ the optical depth perpendicular to the disc (opaque
without $\tau$). An edge-on disc is a dark band of height $H$, an inclined one an ellipse with a rim.
The stellar disk is sampled on 4000 equal-area points weighted by the limb darkening; a sharp edge
crossing a star is resolved to about $10^{-3}$ of its flux. A disc can be fixed on the sky, tilted
from an orbit, precess at a given rate, or precess as a rigid body under the torque of its binary
(`dynamics: rigid`; Zanazzi & Lai 2018, eq. 3; Larwood et al. 1996 for circumstellar discs), with
polar alignment around eccentric binaries (Farago & Laskar 2010; Martin & Lubow 2017). Discs can
shine as flat reprocessing discs, $T = (2/3\pi)^{1/4}(R_*/r)^{3/4}T_*$ (Chiang & Goldreich 1997), and
scatter starlight with an albedo and a Henyey–Greenstein phase function (Henyey & Greenstein 1941).
Not modelled: the disc's mass and gravity on the orbits, forward scattering through optically thin
parts, warps.

![**Figure 18.** Light curves of disc-eclipsing systems computed with SWINGERsim. Left: ε Aur
2009–2011, dimming in magnitudes against year; dashed vertical lines: the observed contacts I and
IV; dotted horizontal line: the observed depth. Middle: EE Cep, depth of every eclipse in magnitudes
against time in years over one precession of its disc (illustrative model; shaded: the observed
range 0.5–2 mag). Right: KH 15D, dimming in magnitudes against time in years as the edge of its
opaque ring sweeps over the binary orbit (illustrative model; phases in which the binary is fully
hidden lie off the scale).](../../figures/discs.png)

Figure 18 shows three regimes. ε Aur: an opaque slab passing in front of a star gives a flat-bottomed
eclipse whose contacts and depth the model matches. EE Cep: the precessing disc changes the depth
from eclipse to eclipse, so the depths trace the precession. KH 15D: the ring edge progressively
covers the orbit of the binary, so the fraction of each orbit in which a star is hidden changes over
years. For the illustrative EE Cep disc (0.8 au) the torque of the B star gives a precession period
of 183 yr against the observed 62–67 yr; a disc of about 1.6 au would match – the precession
constrains the size of the disc.

## What is not modelled {#observables-not-modelled}

Stellar spots, combination frequencies of pulsations, the light and gravity of disc material beyond
the slab model, and astrometry: SWINGERsim computes no photocentre motion and no Gaia astrometric
signal. Relative and absolute positions of the stars on the sky are available in the run output and
could serve for such a model, but no astrometric observable is computed or fitted.

## Implementation {#implementation-observables}

**External libraries and data.** Limb-darkening and gravity-darkening tables of Claret & Bloemen
(2011) and the TLUSTY coefficients of Reeve & Howarth (2016); passband curves of Bessell (1990), TESS
(Ricker et al. 2015) and Gaia DR3 (Riello et al. 2021); PHOEBE 2 (Prša et al. 2016; Conroy et al.
2020) in the optional hybrid. batman, ellc and PHOEBE serve only as references in the tests.
**Own code.** The black-body passband integration and the table interpolation; the exact occultation
of any number of limb-darkened spheres; the light-travel-time solution per binary; the analytic
proximity terms; the Roche and rotating-star surfaces, their meshes, gravity darkening and ray-traced
occultations including self-occultation of contact binaries; the coupling to PHOEBE (chunks, element
conversion, mutual eclipses); radial velocities at emission time and the centre-of-light velocities
on the meshes; the synthetic eclipse timings and O−C; the phenomenological minima with the Mikulášek
profile; the pulsation periods and the phase modulation; and the disc geometry, transmission, light
and rigid precession.
