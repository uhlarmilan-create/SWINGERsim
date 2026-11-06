#!/usr/bin/env bash
# Part II of the SQUADRA test (docs/<lang>/squadra_test.md): figures and PDFs from the outputs of the
# test campaign. Links output/squadra_test/<system>/<variant> to output/campaign/..., integrates the
# published systems forward and backward if missing, then makes figures, numbers and PDFs.
#     bash tools/partII_rebuild.sh
set -euo pipefail
cd "$(dirname "$0")/.."
C=output/campaign
link() { # link <target under output/campaign> <link name under output/squadra_test>
  if [ -e "$C/$1" ]; then mkdir -p "$(dirname "output/squadra_test/$2")"; ln -sfn "$(realpath "$C/$1")" "output/squadra_test/$2"
  else echo "missing: $C/$1"; fi
}
for s in bu_cmi tic120362137 v994_her; do
  for v in fit fit_robust fit_priors fit_auto fit_rv; do link "g2/${s}_$v" "$s/$v"; done
done
link mcmc1000/bu_cmi_fit_final_mcmc_1000 bu_cmi/fit_mcmc
link mcmc1000/bu_cmi_fit_final_mcmc_1000 bu_cmi/fit_lc
link mcmc1000/tic120362137_fit_final_mcmc_1000 tic120362137/fit_mcmc
link mcmc1000/tic120362137_fit_final_mcmc_1000 tic120362137/fit_lc
link g2/tic120362137_fit_lc_roche tic120362137/fit_lc_roche
link g2/v994_her_fit_mcmc v994_her/fit_mcmc
link g2/v994_her_fit_lc_noosc v994_her/fit_lc
link g7b/bu_cmi bu_cmi/ensemble
link g7b/tic120362137 tic120362137/ensemble2
link g7b/v994_her v994_her/ensemble
for s in bu_cmi tic120362137 v994_her; do
  D=output/squadra_test/$s
  [ -f $D/published_forward.h5 ] || swingersim run configs/squadra_test/${s}_published.yaml -o $D/published_forward.h5
  t=$(grep -o "t_max_yr: [0-9.]*" configs/squadra_test/${s}_published.yaml | awk '{print $2}')
  [ -f $D/published_backward.h5 ] || swingersim run configs/squadra_test/${s}_published.yaml -o $D/published_backward.h5 --t-max-yr=-$t
done
python tools/squadra_report.py --root .
python tools/squadra_fits.py --root .
python tools/build_docs_pdf.py part2
echo "PART II DONE"
