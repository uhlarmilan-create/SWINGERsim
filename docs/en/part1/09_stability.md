# Stability: the ghost test, criteria and classifiers {#stability}

An ensemble of multiple stars contains many members that would not survive: their orbits are too close,
and they break up within a few outer orbits. Before the long-term evolution of a member means anything,
SWINGERsim has to decide whether the member is stable. It offers three tools of increasing cost:

1. **Algebraic criteria** applied to the nested triples of a hierarchy: MA01 (Mardling & Aarseth 2001)
   and the criterion of Vynatheya et al. (2022, V22), which includes the dependence on the mutual
   inclination. They cost microseconds and are SWINGERsim's own code.
2. **Neural classifiers** (multi-layer perceptrons, MLP) trained on 2+2 and 3+1 quadruples. SWINGERsim
   builds the training sets itself (N-body ghost tests of 10⁴ quadruples per hierarchy) and trains the
   network with the scikit-learn library; the network library is external, the training data and the
   use are SWINGERsim's.
3. The **ghost test**, the definition of stability of Vynatheya et al. (2023, Sect. 2): the system and
   a "ghost" copy whose inner semi-major axis is larger by 10⁻⁶ au are integrated with N-body for 100
   outer orbits (Newtonian point masses); the system is unstable if it becomes unbound or if the inner
   semi-major axes of the two copies differ by more than 1 %. This is the reference.

## Reproduction of Vynatheya et al. (2023) {#vynatheya2023}

**Reference.** Vynatheya et al. (2023) trained MLP classifiers on 4·10⁵ quadruples labelled by the ghost
test and published the **scores** – the fraction of test members classified correctly – of MA01 and of
their MLP for 2+2 and 3+1 quadruples (their Tables 3–4). Two further measures appear below: the
**recall** of the stable class (the fraction of truly stable members the classifier calls stable) and
the **precision** (the fraction of members called stable that are truly stable).

**Input.** Two training sets of 10⁴ quadruples each, drawn from the prior of the paper (their Sect.
3.1: masses, period ratios, eccentricities, inclinations), classified by the ghost test with a 30-min
limit per member; an MLP with the 11 input parameters and the architecture of the paper (4 hidden layers
of 50 logistic units, Adam optimiser, L2 regularisation 10⁻⁴; their Sect. 3.2), five networks averaged,
trained on 80 % of the members and scored on the remaining 20 %. A third set of 10⁴ 3+1 quadruples was
drawn, as in the paper, only near the MA01 stability boundary.

**Table 43.** Vynatheya et al. (2023): scores against the ghost test. These four benchmarks are manual;
the MLP is class ext (the scikit-learn network on SWINGERsim's training sets), MA01 class own. Columns as defined in [How to read this report](#how-to-read).

| Benchmark | Class | Expected | Computed | Difference | Tolerance | Verdict |
|---|---|---|---|---|---|---|
| `vyn23-2p2-ma01-score` | own | 0.83 | 0.869 | +0.039 | 0.05 | PASS |
| `vyn23-2p2-mlp-score` | ext | 0.94 | 0.9375 | −0.0025 | 0.03 | PASS |
| `vyn23-3p1-ma01-score` | own | 0.56 | 0.7896 | +0.23 | 0.05 | **EXPLAINED** |
| `vyn23-3p1-mlp-score` | ext | 0.93 | 0.9195 | −0.0105 | 0.03 | PASS |

**Table 44.** The three training sets. Columns 2–4: the fraction of members that are stable by the ghost
test, by MA01 and by V22; columns 5–7: the scores of MA01, V22 and the MLP against the ghost test.

| Set | Stable (ghost) | Stable (MA01) | Stable (V22) | Score MA01 | Score V22 | Score MLP |
|---|---|---|---|---|---|---|
| 2+2, prior of the paper | 28.2 % | 39.3 % | 34.9 % | 0.869 | 0.895 | 0.938 |
| 3+1, prior of the paper | 9.5 % | 28.5 % | 26.3 % | 0.790 | 0.812 | 0.920 |
| 3+1, near the MA01 boundary | 21.0 % | – | – | 0.509 | 0.568 | – |

![**Figure 22.** Scores of the three classifiers – MA01 (orange), V22 (green) and the neural network MLP (blue) – against the N-body ghost test, for the 2+2 set (left) and the first 3+1 set (right). The vertical axis is the score, the fraction of test quadruples whose stability is predicted correctly (1 = perfect); the number above each bar repeats it. The black dashed lines are the scores published by Vynatheya et al. (2023), who give none for V22.](../../figures/part1/stability.png)

**Result (Tables 43 and 44).** In Figure 22 the blue bars reach the dashed published lines in both sets: the MLP scores
0.938 (2+2) and 0.920 (3+1) against the published 0.94 and 0.93 – within 0.01, although it was trained
on 40 times fewer systems. MA01 scores 0.869 on the 2+2 set against the published 0.83 (+0.04, within the
tolerance of 0.05 set by the difference between training sets). On the first 3+1 set, however, MA01
scores 0.79 against 0.56. This benchmark is EXPLAINED by the sampling of the test set, for three reasons.
First, with only 9.5 % stable systems a classifier that calls *everything* unstable already scores
0.905; MA01 calls 28.5 % stable – too many – but recognises most of the many unstable members, so its
score is high. Second, the paper drew its members close to the MA01 boundary, and its own Table 4
(precision and recall) implies a score of 0.69, not 0.56. Third, the set drawn near the boundary as in the
paper gives 0.51, close to the published 0.56 (and V22 0.57). This window set has fewer stable members
(21 %) than the about 36 % that the published precisions and recalls imply, which lowers the precision of
the stable class (0.29 against the published 0.54; recall 0.94 against 0.95). A score depends on how the
test set is drawn and is meaningful only together with the sampling. No member of the 3·10⁴ ghost tests
reached the time limit.

![**Figure 23.** Diagnostics of the 3+1 neural network. *Top left:* training loss of each of the five networks against the training epoch. *Top middle:* score on the validation members kept aside during training against the epoch. *Bottom left:* ROC curve – for every threshold of the predicted probability, the fraction of truly stable members called stable (vertical) against the fraction of unstable members called stable (horizontal); dotted: a random classifier; the area under the curve is 0.914. *Bottom middle:* reliability – of the members with a predicted probability $p$, the fraction that is really stable (dotted: ideal). *Bottom right:* histogram of the predicted probabilities of the held-out members by their true class (blue stable, orange unstable).](../../figures/studies/vyn23_3p1_mlp.png)

**Calibration.** Figure 23 shows that the 3+1 network has learnt (the loss falls and levels off; the
validation score rises from 0.90, the score of calling everything unstable, to 0.92–0.93) and that its
probabilities can be taken at face value: the reliability points lie on the diagonal. The histogram shows
the limitation: the few stable 3+1 systems rarely get a probability above 0.5, so the recall of the stable
class at the threshold 0.5 is only 0.20 (precision 0.65; for 2+2: recall 0.90, precision 0.88). For 3+1
systems the probability itself should be used, not the label at 0.5.

**Findings about the criteria.** On the broad sets MA01 and V22 call more systems stable than the ghost
test (2+2: 39.3 % and 34.9 % against 28.2 %; 3+1: 28.5 % and 26.3 % against 9.5 %). They are optimistic:
a system they call stable is not necessarily stable for 100 outer orbits. In ensembles SWINGERsim
therefore reports the algebraic verdicts as diagnostics and uses the N-body backend (or the ghost test)
where stability matters, as in [Chapter 6](#czev4315).

## Chaos maps and simulation-based inference {#chaos-sbi}

**Chaos maps.** SWINGERsim computes MEGNO on a grid of two orbital parameters to map regular and chaotic
regions. The benchmark is the inner stability limit of a circumbinary orbit around an equal-mass circular
binary, published as a grid by Quarles et al. (2018): 2.27 binary separations (the fit of Holman &
Wiegert 1999 gives 2.39). SWINGERsim, requiring all of 8 orbital phases to be regular, finds 2.30
(`chaosmap-circumbinary-limit`, class R; tolerance 0.06, the grid step of 0.05 plus the coarser sampling of
8 phases and 10⁴ binary orbits instead of 10⁵).

**Simulation-based inference.** For posteriors that an ensemble cannot reach by weighting, SWINGERsim can
train a neural posterior estimator with the external `sbi` package (Tejero-Cantero et al. 2020) on its own
simulations. The benchmark is a linear Gaussian model with an analytic posterior: with 2500 training
simulations the posterior width is 0.291 against the exact 0.309 (−5.8 %, tolerance 10 %;
`sbi-linear-gaussian-std`, class ext). This benchmark checks that SWINGERsim passes its simulations to the
package correctly; the inference algorithm itself is external.
