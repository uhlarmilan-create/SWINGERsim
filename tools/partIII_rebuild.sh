#!/usr/bin/env bash
# Part III of the tests (docs/<lang>/part3/): figures, numbers and PDFs from the outputs of the test
# campaign. Links output/squadra_test3/<system>/<variant> to output/campaign/... and
# output/campaign_final/... (the QX Cas and AH Cep fits redone with the pre-1960 minima converted with Delta T,
# tools/campaign_final/run_all.sh; the V907 Sco fit v2 repeated from the best sample of its MCMC, run_finish.sh), then runs tools/squadra3_report.py and builds the Part III PDFs.
#     bash tools/partIII_rebuild.sh [--no-png]   (--no-png: JSON figures only, no kaleido needed)
set -euo pipefail
cd "$(dirname "$0")/.."
O=output/squadra_test3
link() { # link <target under output/> <link name under output/squadra_test3>
  if [ -e "output/$1" ]; then mkdir -p "$(dirname "$O/$2")"; ln -sfn "$(realpath "output/$1")" "$O/$2"
  else echo "missing: output/$1"; fi
}
link campaign_final/v907_sco_fit_v2_polished v907_sco/fit_v2
link campaign/g2/v907_sco_fit_lc_v2 v907_sco/fit_lc_v2
link campaign/mcmc1000/v907_sco_fit_v2_mcmc_1000 v907_sco/fit_v2_mcmc
link campaign/g2/v907_sco_fit_p99 v907_sco/fit_p99
link campaign_final/qx_cas_fit_v2_best qx_cas/fit_v2
link campaign/g7b/qx_cas qx_cas/ensemble
link campaign_final/qx_cas_extra qx_cas/extra
link campaign_final/ah_cep_fit ah_cep/fit
link campaign/g2/ah_cep_fit_lc ah_cep/fit_lc
link campaign/g7b/ah_cep_4body ah_cep/ensemble_4body
P=$O/v907_sco/fit_v2_mcmc/posterior_inc.json
[ -f $P ] || python tools/squadra3_posterior.py 100 $P --fit-dir $O/v907_sco/fit_v2_mcmc
python tools/squadra3_report.py --root . "$@"
if ls docs/figures/squadra_test3/*.png >/dev/null 2>&1; then python tools/build_docs_pdf.py part3
else echo "no PNG figures (run without --no-png, needs kaleido >= 1 with Chrome): PDFs not built"; fi
echo "PART III DONE"
