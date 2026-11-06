#!/usr/bin/env bash
# Final campaign of the SWINGERsim validation: everything whose published numbers depend on the corrections of
# 2026-11-03 (pre-1960 times converted with Delta T, MEGNO seed, corrected published configurations), and the last
# open questions of the campaign "open questions". All runs in one go, on all cores, about three days.
#
#   cd ~/Public/SWINGERsim && git pull && nohup tools/campaign_final/run_all.sh > output/campaign_final.log 2>&1 &
#
# Every step writes output/campaign_final/.done_<step> when it has finished; started again, the script skips the
# finished steps (the long steps also continue where they stopped). Progress: tail -f output/campaign_final.log, or
# output/campaign_final/status.txt (one line per step: end, wall time, result). Stop: kill the PID that
# `cat output/campaign_final/run_all.pid` prints (not pkill), then start the script again later.
#
#   V    swingersim validate --all (Part I: the benchmark tables; MEGNO with its fixed seed)
#   P3   Part III: the QX Cas and AH Cep fits again (their pre-1960 minima now converted with Delta T), the QX Cas
#        depths (tools/qx_cas_extra.py), then tools/partIII_rebuild.sh (figures, numbers, PDFs; corrected published
#        configurations of V907 Sco and QX Cas)
#   S5   BU CMi: profile of chi^2 in the apsidal constant of pair A (7 fits with the errors of S5 frozen)
#   S3   CzeV343: 1000-yr runs on a fine grid of P_A across the 3:2 resonance (141 runs), resonant angles
#   S2   the ghost test of Vynatheya et al. (2023) for all 1200 members of the S2 ensemble
#   S4   QX Cas: MCMC of the refitted fit v3 until 50 autocorrelation times, then the eclipse dates of 300 samples
#   S2L  the 243 members of S2 that were chaotic or below the Mardling-Aarseth limit, 10^4 outer orbits each
set -u
cd "$(dirname "$0")/../.." || exit 1
ROOT=$(pwd)
[ -f .venv/bin/activate ] && source .venv/bin/activate
OUT=output/campaign_final
mkdir -p "$OUT"
echo $$ > "$OUT/run_all.pid"
CFG=configs/campaign_final
NPROC=$(python -c "import os; print(os.cpu_count() or 2)")  # not nproc: it follows OMP_NUM_THREADS
J=${SWINGERSIM_CORES:-0}            # 0 = all cores
JP=$(( NPROC > 1 ? NPROC - 1 : 1 ))  # tools with their own pools: leave one core for the desktop
export PYTHONWARNINGS=ignore

stamp() { date '+%Y-%m-%d %H:%M:%S'; }
step() {
  local name=$1; shift
  if [ -f "$OUT/.done_$name" ]; then
    echo "[$(stamp)] $name: already done, skipped"
    return 0
  fi
  echo; echo "[$(stamp)] ===== $name: $*"
  local t0=$SECONDS
  if "$@"; then
    touch "$OUT/.done_$name"
    echo "$name  $(stamp)  $(( (SECONDS - t0) / 60 )) min  ok" >> "$OUT/status.txt"
    echo "[$(stamp)] $name: done in $(( (SECONDS - t0) / 60 )) min"
  else
    echo "$name  $(stamp)  $(( (SECONDS - t0) / 60 )) min  FAILED" >> "$OUT/status.txt"
    echo "[$(stamp)] $name: FAILED (the campaign goes on with the next step)"
  fi
}
need() {  # need <file>: stop the campaign if an input of an earlier campaign is missing
  if [ ! -e "$1" ]; then echo "[$(stamp)] missing input: $1 - campaign stopped"; rm -f "$OUT/run_all.pid"; exit 1; fi
}

echo "[$(stamp)] final campaign in $ROOT, $NPROC cores"
python -c "import swingersim; print('SWINGERsim', swingersim.__version__)"
git log -1 --oneline 2>/dev/null
need output/campaign_oq/s2_stability/summary.csv
need output/campaign_oq/s5_bu_cmi_k2/stage2_residuals.csv
need output/campaign_oq/s5_bu_cmi_k2/stage2_rv.csv
need output/campaign_oq/s5_bu_cmi_k2/stage2_lc.csv
need output/campaign/g7b/qx_cas/summary.csv

# --- V: benchmarks ----------------------------------------------------------------------------------------------
mkdir -p "$OUT/validation"
validate_all() { swingersim validate --all --html "$OUT/validation/report.html" 2>&1 | tee "$OUT/validation/validate.txt"; return "${PIPESTATUS[0]}"; }
step v_all validate_all

# --- P3: Part III fits with the Delta-T minima -------------------------------------------------------------------
step p3_qx_v2best swingersim fit configs/squadra_test3/qx_cas_fit_v2_best.yaml -o "$OUT/qx_cas_fit_v2_best" --processes "$J"
step p3_qx_v3 swingersim fit configs/squadra_test3/qx_cas_fit_v3.yaml -o "$OUT/qx_cas_fit_v3" --processes "$J"
step p3_ahcep swingersim fit configs/squadra_test3/ah_cep_fit.yaml -o "$OUT/ah_cep_fit" --processes "$J"
step p3_qx_extra python tools/qx_cas_extra.py --fit "$OUT/qx_cas_fit_v2_best/stage2_fitted.yaml" \
  --ensemble output/campaign/g7b/qx_cas -o "$OUT/qx_cas_extra" -j "$JP"
step p3_rebuild bash tools/partIII_rebuild.sh

# --- S5: BU CMi, profile in k2 of pair A -------------------------------------------------------------------------
step s5_data python tools/campaign_final/make_s5_profile.py
for k in 0.0003 0.0010 0.0020 0.0030 0.0042 0.0060 0.0090; do
  step "s5_k2A_$k" swingersim fit "$CFG/s5_k2A_$k.yaml" --processes "$J"
done

# --- S3: CzeV343, fine scan across 3:2 ---------------------------------------------------------------------------
step s3_scan python tools/campaign_final/s3_resonance_scan.py "$OUT/s3_resonance_scan" -j "$J"

# --- S2: ghost test of all members --------------------------------------------------------------------------------
step s2_ghost python tools/campaign_final/s2_extend.py ghost output/campaign_oq/s2_stability "$OUT/s2_ghost" -j "$J"

# --- S4: QX Cas, MCMC until 50 tau -----------------------------------------------------------------------------------
step s4_config python tools/campaign_oq/make_s4_config.py --start "$OUT/qx_cas_fit_v3/stage2_fitted.yaml" \
  --tau-target 50 --max-steps 6000 --iterations 40 --output-dir ../../output/campaign_final/s4_qx_cas_mcmc \
  -o "$CFG/s4_qx_cas_mcmc.yaml"
step s4_mcmc swingersim fit "$CFG/s4_qx_cas_mcmc.yaml" --processes "$J"
step s4_dates python tools/campaign_oq/posterior_dates.py "$OUT/s4_qx_cas_mcmc" 300 --years 1900 2150 \
  --processes "$JP"

# --- S2L: the long runs --------------------------------------------------------------------------------------------
step s2_long python tools/campaign_final/s2_extend.py long output/campaign_oq/s2_stability "$OUT/s2_long" -j "$J"

echo; echo "[$(stamp)] final campaign finished"; cat "$OUT/status.txt"
rm -f "$OUT/run_all.pid"
