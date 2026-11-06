# Discussion

## What makes these results reliable, and where they are weaker

- **V907 Sco** rests on the best data: inclinations from four independent sources, velocities from
  two epochs 40 years apart, precise TESS minima, and a posterior sampled with MCMC. The SWINGERsim
  result agrees with the published one within the errors, and SWINGERsim reproduces the published
  predictions when given the published parameters. The orientation of the orbits on the sky (the
  nodes) is not given in the paper; only the difference of the nodes enters any observable.
- **QX Cas** rests on two inclinations, one upper limit and the absence of eclipses. The ensemble is
  therefore the main result, and the fit is an example. Computing the eclipse depths showed that the fit
  satisfies the 1930s only geometrically, not observably, and that the old depths need third light and
  different stellar parameters. The TESS detection of the binary is an independent check of the model's
  amplitude, but not yet of its orbit.
- **AH Cep** rests on the light-time effect alone for the third body. The third-light result depends
  on the light-curve model.

## A lesson: "the model eclipses" is not "the eclipses were seen"

For all three systems the switching dates were first compared as geometric crossings of the eclipse
limit. For V907 Sco this was enough, because its eclipses deepen quickly after the crossing and the
model depths at the observed switch dates were checked (Table 2). For QX Cas the same comparison hid a
real problem: a geometric eclipse in 1933 that is 0.0006 mag deep is no eclipse for a photographic
plate. Whenever an observation says "eclipses were seen", the model has to give a depth above the
detection limit of that observation, with the third light stated. We now enter such observations in that
form.

## The light-curve model

In all three systems SWINGERsim's synthetic light curves were compared with TESS. They reproduce the
eclipses of V907 Sco to within 5 % in depth without being fitted, and the heartbeat amplitude of QX Cas.
Where depths are fitted, SWINGERsim's eclipses for given parameters are systematically about 5 %
shallower than PHOEBE's (V907 Sco, AH Cep). This shifts fitted inclinations upwards by a fraction of a
degree for V907 Sco and probably by more for the shallow, partial eclipses of AH Cep. Results that rest
on eclipse depths (the light-curve fit of V907 Sco, the inclination and third light of AH Cep, the
depths and detectability of the QX Cas eclipses, the shape of the QX Cas signal) should be read with
this in mind. The dynamical results (eclipse times, velocities, geometric switching dates) do not depend
on it.

## Errors

The errors are linear errors from the curvature of χ², except for V907 Sco, where the MCMC posterior
errors are given; for V907 Sco both agree within 50 %. The MCMC chain covers about 19 autocorrelation
times in total and about 9 after the burn-in half is discarded (see Data and methods); a chain of
400 steps gave errors within a median of 10 % of these, so we consider them correct to a few tens of
percent. For strongly correlated parameters (for example mass and orbital inclination of
the third body of AH Cep) linear errors understate the true range.
