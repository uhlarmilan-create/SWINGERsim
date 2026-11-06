# Open questions from the tests {#open}

Each item states the finding, why it matters, and what would settle it. The first group concerns
individual systems, the second the physics and methods.

## Individual systems {#open-systems}

**QX Cas: the old eclipses and the stellar parameters** (Part III, QX Cas)

- *Finding:* With the adopted radii (3.5 and 5.1 R☉) and temperatures, the photographic minima of
  1933–1942 need an inclination of at least 76–78°. A fit that enforces this keeps the eclipses
  detectable until 1999.0, while none were seen in the mid-1990s. The depths of 1960 and 1980 need
  13–47 % of third light, more than a single value can supply.
- *Why it matters:* The dates of the eclipses, past and future, depend on it.
- *What would settle it:*
  - a joint fit of the radii, the light ratio and the third light together with the old depths and
    the TESS heartbeat light curve;
  - the DASCH plates of 1900–1950;
  - the unpublished light-curve solutions of Bonaro et al. (2009).

**QX Cas: a range of dates** (Part III, QX Cas)

- *Finding:* A random ensemble weighted afterwards keeps only one or two systems consistent with the
  whole inclination history. An MCMC of the dynamical fit with the inclination history as data has
  now been run (24 walkers × 6000 steps, 30 h). The eclipses ended in 2000–2006.5 and return in
  2072–2130 (16–84 %, 300 integrated posterior samples).
- *Why it is still open:* The times of minima do not fix the period of the third star. χ² is flat
  from 270 to 400 d (the best sample in every 20-d bin lies within 3 of the optimum), and the
  autocorrelation time of the chain is 311 steps, so even 6000 steps are only 9.6 τ after the burn-in.
  The median return date moved from 2096 (1500 steps) to 2110 (6000 steps). More computing does not
  fix this; the data do not contain the information.
- *What would settle it:* An independent measurement of the outer orbit: radial velocities of the
  third star or of the binary over a few hundred days, or Gaia epoch astrometry.

**V907 Sco: Gaia sees less motion than modelled** (Part III, Gaia)

- *Finding:* The RUWE is 0.81, whereas the model moves the light centre by about 0.3 mas.
- *Why it matters:* It is the only independent check of $m_3$ and of the light of the third star.
- *What would settle it:* Gaia DR4 epoch astrometry. Until then, a simulation of the Gaia scanning of
  this source would show whether 0.3 mas must raise the RUWE.

**AH Cep: inclination and third light** (Part III, AH Cep)

- *Finding:* TESS prefers $i \approx 72.5°$ with about 5 % third light. The published value is 69.8°,
  and the spectra allow less than 2 % third light.
- *What would settle it:* An independent light-curve analysis, for example with PHOEBE, of the same
  TESS sectors. Also the Hipparcos–Gaia acceleration for the 70-yr orbit.

**BU CMi: the classical apsidal advance** (Part II, BU CMi)

- *Finding:* With the apsidal advance that stellar models predict (0.56 and 0.33°/yr), the fit moves
  the outer eccentricity by 2.6σ. Without it the fit agrees with the paper. A fit with the apsidal
  constants free and a χ² profile in $k_2$ of pair A (errors held fixed) give $k_{2,A} \le 0.0012$
  (Δχ² = 1) and $\le 0.0023$ (Δχ² = 4), against 0.0042 of the models (Δχ² = +13.5 there, with
  $k_{2,B}$ pushed to its lower bound). The two pairs trade off ($k_{2,A} + k_{2,B} \approx 0.0045$).
  Freeing the constants moves $e_\mathrm{out}$ only from 0.2631 to 0.2663 (published 0.2728 ± 0.0037).
- *Why it matters:* The data want a weaker tidal and rotational apsidal advance than models of
  synchronously rotating stars give. Lower $k_2$, slower (sub-synchronous) rotation or missing physics
  would each do it; this is a tension, not a measurement of $k_2$.
- *What would settle it:* Projected rotational velocities ($v \sin i$) of the four stars, which fix
  the rotational part of the advance, and apsidal constants from models of the actual stars.

**V994 Her: pair C** (Part II, V994 Her)

- *Finding:* The published apsidal rate of pair C needs apsidal constants 1.8 times larger than
  stellar models give, and the A–B eccentricity is 0.51–0.59 from the minima against 0.69 published.
- *What would settle it:* A closer, unknown companion of C would explain the first point. An
  ensemble of 3000 companions (0.1–2 M☉, 8–3000 d) shows which ones give the measured apsidal
  period. The existing 56 primary minima of C limit a timing signal of 80–400 d to below 2.4 min
  (2σ, median over the period), and that excludes 55 % of those companions. The rest have
  0.13–0.31 M☉ on 80–182-d orbits. The highest peak in the times of minima (80.9 d, 2.0 ± 0.5 min)
  is not a detection: aliases at 49–84 d are almost as high. Dense timing of C over at least
  two such periods is the test. The crowding correction of TESS sector 26 must be checked for the
  depths.

**CzeV343: how often does binary A widen?** (Part IV, CzeV343)

- *Finding:* COMPAS (Pejcha et al. 2022) widens binary A in ≳ 85 % of realisations, SWINGERsim with
  MIST tracks in 10–29 %. The same 300 draws run with SSE tracks (the Hurley et al. formulae on which
  COMPAS is built) leave a stripped donor in 82 % of them and a separation above 50 R☉ in 83 %; the
  outcome of the two track families agrees in only 27 % of the draws. With SSE 75 % of the quadruples
  become unstable, with MIST 14 %.
- *What remains:* The difference lies in the stellar tracks (the response of the donor and so the
  stability of the mass transfer), not in the dynamics. Which family is closer to the truth needs a
  detailed binary-evolution code such as MESA on identical inputs.

**V907 Sco: the spin of the B stars** (Part III, V907 Sco)

- *Finding:* The nodal precession of the binary (63 yr) and the spin precession of its stars
  (~36 yr at synchronous rotation) are nearly equal: $\eta = g/\alpha \approx 0.58$ against the
  bifurcation of the Cassini states at 0.587. 600 N-body realisations with free rotation periods
  and initial obliquities show obliquity changes of 81–82° (median) over 600 yr for rotation periods of
  3.5–4.5 d, 13° at 0.5–1 d and 42–44° at 10–30 d. If tides have damped the spin into its Cassini
  state (6.7°; synchronisation time 9–940 Myr against an age of 400 Myr), the obliquity is small and
  constant.
- *What would settle it:* The Rossiter–McLaughlin effect during the eclipses, which last until
  about 2052, and $v \sin i$ of the two stars.

**CzeV343: is it in the 3:2 resonance?** (Part IV, CzeV343)

- *Finding:* No. A scan of 141 periods of binary A around the exact commensurability (1000 yr each)
  finds the 3:2 Lindblad resonance of Tremaine (2020), with the angle $3\lambda_A - 2\lambda_B -
  \varpi_A$ and a width of about $10^{-5}$ d in $P_A$. The present system lies $7.9 \times 10^{-4}$ d,
  about 80 widths, away. None of its satellites split off by the outer orbit
  ($k = -6 \ldots +6$) is active either. Why 2+2 quadruples cluster near such ratios (Kolář et al. 2025) remains open.
- *What would settle it:* A capture experiment with slow migration through the resonance.

**CzeV4315: the mass of the third star of B** (Part IV, CzeV4315)

- *Finding:* The mass is degenerate with the eccentricity and inclination of its orbit (1.4–4.9 M☉,
  16–84 %).
- *What would settle it:*
  - the amplitude of the 119-day timing signal, or radial velocities of B;
  - a fit of the eclipse times of B with the 2+(2+1) model, which SWINGERsim can do once the minima
    are available.

## Physics and methods {#open-physics}

- **Apsidal constants of stripped and evolved stars.** The future of CzeV343 hinges on $k_2$ of a
  stripped giant: about $10^{-3}$ by our two-zone estimate. SWINGERsim takes $k_2$ from tables or
  polytropes for normal stars. For stripped donors, giants with thin envelopes and mass gainers, a
  structure model, for example from MESA profiles (the `mesa-reduce` tool), is the reliable source.
- **MA01 vs V22 near the boundary (answered for one configuration).** At the onset of instability
  in CzeV343 the criterion of Mardling & Aarseth gives 1.00, that of Vynatheya et al. 1.19–1.24.
  1200 N-body realisations of CzeV343 after the mass transfer, followed for 200 outer orbits and the
  243 doubtful ones for $10^4$ outer orbits, give the answer (Part IV, CzeV343). Below a Mardling–Aarseth
  ratio of 1, 81, 59, 33 and 11 % of the members break up in the bins 0.6–0.7 … 0.9–1.0; above 1.1 none
  does. The Mardling–Aarseth criterion finds 98 % of the break-ups. That of Vynatheya et al. (2022)
  finds 51 % at its threshold 1 and works best at about 1.17. The neural network of Vynatheya et al.
  (2023) finds 25 %, and only 32 % even of the systems its own training test (a ghost orbit) calls
  unstable. These systems lie just outside its training range (mass ratio 10.7 > 10).
- **Convergence of MCMC fits.** The 1000-step chains are about 19 autocorrelation times long; the
  half kept after the burn-in is about 9.4 τ. For QX Cas, whose outer period the data do not fix,
  τ is 311 steps, so even 6000 steps are only 9.6 τ after the burn-in (see above). An MCMC also checks
  the least-squares fit: for V907 Sco its best sample has a χ² 2.3 below the optimum (with the final
  jitters of the fit). A least-squares fit started at that sample, with the jitters re-estimated as in
  Part III, returned to the same optimum (χ² 222.73, every parameter within 0.03σ), so the optimum does
  not depend on the start. Δχ² = 2.3 lies far inside the 1σ region of 14 parameters (Δχ² ≈ 16): the
  valley is flat, and fit and MCMC agree. Comparing the lowest sample with the optimum requires removing
  the term −2 ln e that the ecosw/esinw parametrisation adds to −2 ln p. The medians are stable, and the widths of the 95 %
  ranges are less certain than those of the 1σ intervals. About 1100 steps would reach 10 τ after the
  burn-in (Foreman-Mackey et al. 2013), and about 5300 steps the 50 τ of the emcee documentation.
- **Main-sequence tides.** The equilibrium tide with inertial waves reproduces the circularisation of
  young binaries. It gives 7.5 d instead of the 12.1 ± 1.2 d of M67 at 4 Gyr (Part I). Stronger
  dissipation (resonance locking, non-linear dynamical tides) is missing from all published
  equilibrium-tide models and lies outside the scope of SWINGERsim.
- **Stellar tracks.** SSE and MIST differ by about 10 % in main-sequence lifetimes and radii at its
  end. When the timing of an event matters, both families should be compared (Part I).
