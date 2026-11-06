---
title: "Testing SWINGERsim"
subtitle: "Part III: Eclipsing binaries with a changing inclination – V907 Sco, QX Cas and AH Cep"
---

**Abstract.** When an eclipsing binary has a third star on an inclined orbit, the orbit of the binary
slowly turns (precesses) like a spinning top. Its inclination to our line of sight changes, so over
decades the eclipses deepen, become shallower, or switch off and on. We test whether SWINGERsim, an
N-body simulator for bound systems of three to ten stars, can take such a system from the published
data, find a model that reproduces a century of observations, and say what comes next. Three binaries
stand for three situations.

**V907 Sco** has a complete published solution (Zasche et al. 2023). From 58 TESS times of minima,
172 radial velocities and 14 inclinations measured since 1900, SWINGERsim finds a third star of
0.98 (+0.05/−0.04) M☉ on a 142.09 ± 0.04-d orbit, inclined by 25.8 ± 0.8° to the binary, with a nodal period of 65.5 (+0.9/−1.0) yr.
All of these agree with the published values within about one standard deviation (the published
ranges of m₃, the mutual inclination and the nodal period are 90 % limits). The model switches the eclipses off
and on in the same years as the plates, the old photometry and TESS show. It predicts the depths of the
2021 and 2023 TESS eclipses to within 5 %, and it puts the end of the present eclipses in 2052–2057
(published: about 2049). The radial velocities of Lacy et al. (1999) fit the same 142-d orbit; a fit
started at their 99-d period leaves it and ends in a clearly worse solution.

**QX Cas** stopped eclipsing in the 1990s, and its third star has never been seen (Guinan et al.
2012a,b). SWINGERsim shows which third stars can turn the orbit fast enough (most likely 0.5–3 years and
0.5–4 M☉) and that the inclination must have fallen at least twice as fast as the published linear
rate. It also finds the binary in the TESS light curves of 2019–2024 through its light variation at the
orbital period, at the amplitude the model predicts. Computing the eclipse **depth**, not only whether
an eclipse is geometrically possible, shows a weak point: in our best model the eclipses of the 1930s
are far too shallow for the photographic plates that recorded them. A model with deep eclipses in the
1930s exists, but it keeps the eclipses until the late 1990s, when none were seen. The old eclipse
depths also require 13–47 % of third light, more than any third star of the model provides. Both point
to the adopted stellar parameters, not to the dynamics.

**AH Cep** is a massive binary for which a change of the eclipse depth was once claimed. SWINGERsim
reproduces 125 years of eclipse times with a third body on a 70-yr orbit of at least 6.6 M☉. Such a
distant companion changes the inclination by only 0.02° in 400 years, so it cannot cause the claimed
change. A fourth body that tilted the pair noticeably would show in the eclipse times or in the TESS
eclipse depths, and the 9-yr fourth body once proposed would be dynamically unstable. The TESS light
curves are best described by an inclination of about 72.50 ± 0.21° (published 69.8°) with about 5 % of extra light;
this is in tension with the < 2 % third light allowed by the spectra and depends on the light-curve
model.

**Gaia DR3** shows no RUWE above 1.4 for the three systems. For V907 Sco, where our model moves the
light centre by about 0.3 mas, this is a tension worth following up with Gaia DR4.

# Aim

Earlier parts of the SQUADRA test checked SWINGERsim on compact multiple stars with complete published
solutions. Here the question is harder: **given what has been measured over a century — a few
inclinations, a list of eclipse times of very different quality, some radial velocities — can
SWINGERsim find the system that produced it, check the published conclusions, and say something new?**

For each system we ask four things:

1. What do the published papers conclude, and from which data?
2. What does SWINGERsim obtain from the same data?
3. Does SWINGERsim **confirm** the papers, **add** to them, or **disagree** with them, and why?
4. What does SWINGERsim predict that can be checked with data that already exist (here mostly TESS and
   Gaia)?
