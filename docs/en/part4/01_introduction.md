---
title: "Testing SWINGERsim"
subtitle: "Part IV: Special cases – the future of CzeV343, the hierarchy of CzeV4315 and the convergence of the posteriors"
---

# Aim {#aim}

Parts I–III of this test series asked whether SWINGERsim reproduces what is already known: published
benchmarks of the physics (Part I), three well-observed multiple stars (Part II) and three eclipsing
binaries whose inclination changes because of a third body (Part III). Part IV asks questions to
which no published answer exists. These are the questions for which the software was built.

1. **How and when does CzeV343 fall apart?** CzeV343 is a doubly eclipsing 2+2 quadruple.
   Pejcha et al. (2022) found that the quadruple becomes dynamically unstable if one of its binaries
   widens beyond about 50 R☉ by mass transfer. They did not say when that happens or what follows
   ("dynamical instability caused by orbital expansion … has not been explored so far in 2+2
   quadruples"). We follow the stellar evolution and mass transfer of both binaries to the onset of
   the instability. We then integrate the four-body motion from that moment on
   ([Section 2](#czev343)).
2. **Is CzeV4315 a 2+2 quadruple or a 2+(2+1) quintuple?** In CzeV4315 the inclination of
   binary B changes by 1.4° per year (Zasche et al. 2025). B also shows a 119-day cycle in its
   eclipse times. We compare two hierarchies with the same data, weigh them by the Bayesian
   evidence, and ask what the preferred one predicts ([Section 3](#czev4315)).
3. **Are the posteriors of the fits converged?** The Markov-chain Monte Carlo (MCMC) fits of
   Parts II and III were run for 1000 steps. We check how far they are from convergence and what
   this means for the quoted uncertainties ([Section 4](#mcmc)).

A word on the two kinds of result in this part. A **forward model** starts from the present system
and integrates it into the future; its answer is a distribution over the unknown parameters, not a
single date. An **ensemble** is a set of such models whose unknown parameters are drawn from stated
prior ranges. Each member gets an **observation weight**: the likelihood of the observed data given
that member. The mean weight over the ensemble estimates the **Bayesian evidence** $Z$ of the
hypothesis behind the ensemble; the ratio of two evidences, the **Bayes factor**, says how much
better one hypothesis explains the data than the other. The **effective sample size**
$n_\mathrm{eff} = (\sum w)^2/\sum w^2$ counts how many members carry the result. A result resting
on a handful of members is a hint, not a measurement. The methods are described in *Physics and
Methods*, chapters Ensembles and inference and Fitting observed data.

All numbers in this part come from the output files listed in
[Data and reproducibility](#data). The figures are made by `tools/part4_figures.py`.
