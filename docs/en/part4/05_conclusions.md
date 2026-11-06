# Conclusions {#conclusions}

Part IV used SWINGERsim on three questions that the literature leaves open.

1. **CzeV343.**
   - The quadruple becomes dynamically unstable only if the donor of binary A is stripped and the
     binary widens to about 57–62 R☉. With the mass-transfer prescriptions of SWINGERsim this
     happens in 14 % of the futures of the fiducial solution (1.3–1.6 Gyr from now) and in 1 % of
     the alternative solution; in most futures binary A merges or becomes a contact binary first.
     How often A widens depends on the stellar tracks: the same draws with SSE instead of MIST
     widen it in 82 % (COMPAS, built on the same formulae: ≳ 85 %; Pejcha et al. 2022), and 75 %
     of them become unstable. The onset separation and what follows do not depend on this.
   - Near the onset the Mardling–Aarseth criterion is the safe one for such systems: in 1200
     N-body realisations it finds 98 % of the break-ups within $4 \times 10^4$ yr, the criterion of
     Vynatheya et al. (2022) 51 % and their neural network 25 %.
   - CzeV343 is not in the 3:2 resonance between its binaries; it lies about 80 resonance widths
     away.
   - After the onset the system is not torn apart. With the apsidal constant computed for the
     stripped giant ($k_2 \approx 10^{-3}$), ZLK cycles driven by binary B raise the eccentricity
     of A to 0.7, and the giant collides with its companion within 50–550 years in all
     600 integrations.
   - Only an implausibly large envelope ($k_2 \gtrsim 0.01$) would quench ZLK and keep the system
     intact for $10^4$ yr.
   - The answer to "when and how" is therefore: if binary A widens (rarely in our models, often in
     COMPAS), the instability starts at 57–62 R☉, in our models about 1.3–1.6 Gyr from now, and ends
     within centuries in a merger inside the quadruple driven by ZLK, not in an ejection.
2. **CzeV4315.**
   - The fast inclination change of binary B cannot come from binary A. In a 2+2 the perturbation
     that turns B also turns A by 2–50°/yr, or the system collides.
   - A third star on the 119-day orbit around B explains all the observations with an evidence
     larger by a factor of at least 20 (formally $10^{13}$). Zasche et al. (2025) left the
     arrangement open, with a weak hint in the same direction.
   - Its mass is degenerate with the eccentricity and inclination of its orbit. 74 % of the weight
     lies above 2 M☉, which spectroscopy or the light-travel-time amplitude of the 119-day signal
     can test.
   - A further physical companion at 1.84″ makes the system at least sextuple.
3. **Convergence of the fits.**
   - The 1000-step MCMC runs are about 19 autocorrelation times long; the half kept after the
     burn-in is about 9.4 τ.
   - Their medians agree with the 400-step runs within 0.8σ, and their uncertainties are larger by
     3–10 % typically.
   - The values quoted in Parts II and III are robust. Their 1σ uncertainties are reliable to a few
     tens of per cent, and their 95 % ranges somewhat less.

**What these tests show about the software.** Each question needed several parts of SWINGERsim
working together:

- stellar evolution and mass transfer with stability tracking (CzeV343, stage 1);
- N-body integration with relativity, tides and rotation, set up automatically from an evolved
  state (stage 2);
- ensembles with observation weights and evidences, across two different hierarchies (CzeV4315);
- the fit chain with MCMC (Section 4).

The answers depend on physical inputs that the observations do not fix: the apsidal constant of a
stripped giant, the mass of an unseen star. SWINGERsim makes these dependences explicit and
quantifies them instead of hiding them in a single run.

**Open points**, collected in *Open questions and future work*:

- a detailed binary-evolution code on the CzeV343 inputs, to decide between the MIST and SSE
  behaviour of the donor of binary A;
- spectroscopy of CzeV4315, to detect or exclude a massive third star around B, and a measurement
  of the amplitude of its 119-day timing signal;
- longer MCMC chains (≈ 1100 steps for 10 τ after the burn-in, ≈ 5300 for 50 τ) where the tails of a
  posterior matter.

# Data and reproducibility {#data}

**Table 8.** Inputs and outputs of Part IV. Each row names a configuration or tool (input) and the
output folder it produces; paths are relative to the SWINGERsim directory.

| item | input | output |
|:--|:--|:--|
| CzeV343 stage 1 and stage-2 set-ups | `tools/czev343_future.py`, `configs/special/czev343.yaml` | `output/special/czev343/stage1.csv`, `stage2/*.yaml` |
| CzeV343 stage 2 | `swingersim ensemble output/special/czev343/stage2/<onset>_k2_<k2>.yaml` | `output/special/czev343/s2tides/` |
| apsidal constant of the stripped giant | `tools/k2_condensed.py` | Table 2 |
| CzeV4315, 2+2 | `configs/special/czev4315_2p2_nbody_ensemble.yaml` | `output/special/czev4315_2p2/` |
| CzeV4315, 2+(2+1) | `configs/squadra/czev4315_2p2p1_wide_ensemble.yaml` | `output/special/czev4315_2p2p1_wide/` |
| 1000-step MCMC | fit configurations of Parts II and III with 1000 steps | the fit folders of Parts II and III |
| figures and numbers | `tools/part4_figures.py` | `docs/figures/part4/` |

The commands are:

```bash
python tools/czev343_future.py -o output/special/czev343          # stage 1 + stage-2 set-ups
for f in output/special/czev343/stage2/*_k2_*.yaml; do
  swingersim ensemble "$f" -o "output/special/czev343/s2tides/$(basename "$f" .yaml)" -j 12
done
swingersim ensemble configs/special/czev4315_2p2_nbody_ensemble.yaml -o output/special/czev4315_2p2 -j 12
swingersim ensemble configs/squadra/czev4315_2p2p1_wide_ensemble.yaml -o output/special/czev4315_2p2p1_wide -j 12
python tools/k2_condensed.py
python tools/part4_figures.py
```
