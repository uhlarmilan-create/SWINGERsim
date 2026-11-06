#!/usr/bin/env bash
# Final campaign, last part (after tools/campaign_final/run_all.sh): what its analysis on 2026-11-06 left to do.
#
#   1. QX Cas depths and dates (tools/qx_cas_extra.py) again with ALL weighted members of the ensemble (--top 5000):
#      run_all.sh used the default --top 150, the published version integrated all 1654 weighted members.
#   2. V907 Sco: the fit v2 of Part III repeated from the best point of its MCMC (2.3 lower in chi^2 than the optimum;
#      tools/campaign_final/make_v907_polish.py). Result: it returned to the same optimum (all parameters within 0.03 sigma).
#   3. tools/partIII_rebuild.sh (figures, numbers, PDFs of Part III) with the new outputs. It needs Kaleido with
#      Chrome for the PNG figures: once, before this script, run `source .venv/bin/activate && plotly_get_chrome -y` (Chrome for Kaleido).
#
#   cd ~/Public/SWINGERsim && git pull && nohup tools/campaign_final/run_finish.sh > output/campaign_finish.log 2>&1 &
#
# About 1.5 h. Steps write output/campaign_final/.done_<step> and a line to output/campaign_final/status.txt, like run_all.sh.
set -u
cd "$(dirname "$0")/../.." || exit 1
[ -f .venv/bin/activate ] && source .venv/bin/activate
OUT=output/campaign_final
CFG=configs/campaign_final
mkdir -p "$OUT"
echo $$ > "$OUT/run_finish.pid"
NPROC=$(python -c "import os; print(os.cpu_count() or 2)")  # not nproc: it follows OMP_NUM_THREADS
J=${SWINGERSIM_CORES:-0}
JP=$(( NPROC > 1 ? NPROC - 1 : 1 ))
export PYTHONWARNINGS=ignore

stamp() { date '+%Y-%m-%d %H:%M:%S'; }
step() {
  local name=$1; shift
  if [ -f "$OUT/.done_$name" ]; then echo "[$(stamp)] $name: already done, skipped"; return 0; fi
  echo; echo "[$(stamp)] ===== $name: $*"
  local t0=$SECONDS
  if "$@"; then
    touch "$OUT/.done_$name"; echo "$name  $(stamp)  $(( (SECONDS - t0) / 60 )) min  ok" >> "$OUT/status.txt"
    echo "[$(stamp)] $name: done in $(( (SECONDS - t0) / 60 )) min"
  else
    echo "$name  $(stamp)  $(( (SECONDS - t0) / 60 )) min  FAILED" >> "$OUT/status.txt"
    echo "[$(stamp)] $name: FAILED"
  fi
}

echo "[$(stamp)] final campaign, last part, $NPROC cores"; git log -1 --oneline 2>/dev/null
# Kaleido must be able to write a PNG, otherwise the Part III rebuild fails at its end (as on 2026-11-04)
if ! python -c "import plotly.graph_objects as go, plotly.io as pio; pio.write_image(go.Figure(), '/tmp/swingersim_kaleido_check.png')" 2>/dev/null; then
  echo "[$(stamp)] Kaleido cannot write PNG figures (Chrome missing). Run: source .venv/bin/activate && plotly_get_chrome -y   - then start this script again."
  rm -f "$OUT/run_finish.pid"; exit 1
fi
rm -f /tmp/swingersim_kaleido_check.png

# the run of 2026-11-04 with --top 150 is kept for comparison
if [ -f "$OUT/qx_cas_extra/summary.json" ] && [ ! -d "$OUT/qx_cas_extra_top150" ]; then
  mv "$OUT/qx_cas_extra" "$OUT/qx_cas_extra_top150"
fi
step f2_qx_extra python tools/qx_cas_extra.py --fit "$OUT/qx_cas_fit_v2_best/stage2_fitted.yaml" \
  --ensemble output/campaign/g7b/qx_cas --top 5000 -o "$OUT/qx_cas_extra" -j "$JP"
step f2_v907_config python tools/campaign_final/make_v907_polish.py
step f2_v907_fit swingersim fit "$CFG/v907_sco_fit_v2_polished.yaml" --processes "$J"
rm -f "$OUT/.done_p3_rebuild"
step p3_rebuild bash tools/partIII_rebuild.sh

echo; echo "[$(stamp)] last part finished"; tail -5 "$OUT/status.txt"
rm -f "$OUT/run_finish.pid"
