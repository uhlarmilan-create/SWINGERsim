# Are the posteriors converged? A 1000-step MCMC {#mcmc}

## The question {#mcmc-question}

The uncertainties quoted in Parts II and III come from Markov-chain Monte Carlo (MCMC) sampling of
the posterior probability of the fitted parameters. We use the affine-invariant ensemble sampler
emcee (Foreman-Mackey et al. 2013); the fit chain is described in *Physics and Methods*, chapter
Fitting observed data. An MCMC run is a set of **walkers** that move through parameter space. Each
**step** costs one full N-body integration per walker. After a **burn-in** the positions of the
walkers sample the posterior.

The positions of one walker are correlated from step to step. The **integrated autocorrelation
time** $\tau$ is the number of steps after which a walker has "forgotten" where it was. A chain of
$N$ steps contains roughly $N/\tau$ independent samples per walker. Foreman-Mackey et al. (2013)
suggest running "a few (say 10) autocorrelation times"; the emcee documentation recommends at least
50 τ before the estimate of $\tau$ itself is reliable. SWINGERsim discards the first half of the
chain as burn-in, measures $\tau$ on the second half, and warns when the kept half is shorter than
10 τ.

The first fits used 400 steps, for which emcee reported $\tau \approx 21$ (19 $\tau$). To test
whether the quoted uncertainties are reliable, we re-ran three fits with 1000 steps: one 2+2 system
(BU CMi), one 3+1 system (TIC 120362137) and one 2+1 system with a changing inclination (V907 Sco).

## Results {#mcmc-results}

**Table 7.** The 1000-step runs. Columns:

- *run time*: wall-clock time on 12 cores;
- *$\tau$*: the mean integrated autocorrelation time over the parameters, in steps;
- *steps/$\tau$*: the chain length in units of $\tau$;
- *error ratio*: 1σ uncertainties from 1000 steps divided by those from 400 steps, as the median
  over all parameters and the range;
- *largest shift*: the largest change of a median value between the two runs, in units of its
  uncertainty.

| system | run time | $\tau$ [steps] | steps/$\tau$ | error ratio (median; range) | largest shift |
|:--|--:|--:|--:|:--|:--|
| TIC 120362137 | 10.5 h | 53.3 | 18.7 | 1.05; 0.86–2.02 | < 0.4σ |
| V907 Sco | 19.3 h | 53.9 | 18.6 | 1.10; 0.86–1.53 | < 0.8σ |
| BU CMi | 19.4 h | 52.3 | 19.1 | 1.03; 0.87–1.36 | < 0.25σ |

Three things stand out.

1. **The autocorrelation time was underestimated.** It grew from about 21 to about 53 steps when
   the chain was made 2.5 times longer. This is the known behaviour of short chains: $\tau$ cannot
   be measured from a chain shorter than a few tens of $\tau$, and the estimate grows with the chain
   (Foreman-Mackey et al. 2013; Sokal 1997). The chains are therefore about 19 $\tau$ long in total, and the
   half kept after the burn-in is about 9.4 $\tau$. That is just below the 10 τ of Foreman-Mackey et
   al. and well short of the 50 τ of the emcee documentation.
2. **The results did not move.** No median moved by more than 0.8σ, and most by much less. The
   typical uncertainty grew by 3–10 %. In a few weakly constrained parameters with broad,
   non-Gaussian posteriors it grew by up to a factor of two (TIC 120362137) or 1.5 (V907 Sco).
3. **MCMC and linearised uncertainties agree.** The covariance matrix of the least-squares fit (the
   "linear" errors) agrees with the MCMC widths to about 10 % in the median. Both describe the same
   well-defined minimum.

## What this means for the quoted numbers {#mcmc-meaning}

- The **medians** quoted in Parts II and III are robust. Doubling the chain length would not move
  them by more than a fraction of their uncertainty.
- The **1σ uncertainties** are reliable to a few tens of per cent. They are slightly too small
  rather than too large, by about 5–10 % typically and more for the weakly constrained elements.
- The **tails** (95 % ranges) are less certain than the 1σ widths. With about 9–10 independent
  samples per walker they rest on fewer samples. This matters most for the dates of future events
  quoted as 95 % ranges, such as the end of the eclipses of V907 Sco in 2052.5–2056.8. The width of
  such ranges is itself uncertain, more than that of the 1σ intervals.
- Ten τ in the kept half would need about 1100 steps, just above the 1000 used. The 50 τ of the emcee
  documentation would need about 5300 steps, i.e. 55–100 h per system on 12 cores. The comparison of
  400 and 1000 steps shows that longer chains refine the widths but do not change the conclusions.
  We report the 1000-step results.

The values in Parts II and III are those of the 1000-step runs.
