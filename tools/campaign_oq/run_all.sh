#!/usr/bin/env bash
# Campaign "open questions" (S1-S7): all runs in one go, on all cores, about a day.
#
#   cd ~/Public/SWINGERsim && nohup tools/campaign_oq/run_all.sh > output/campaign_oq.log 2>&1 &
#
# Every step writes output/campaign_oq/.done_<step> when it has finished; started again, the script skips the finished
# steps (ensembles and chaos maps also continue where they stopped). Progress: tail -f output/campaign_oq.log, or
# output/campaign_oq/status.txt (one line per step: start, end, wall time, result).
# Stop: kill the PID that `cat output/campaign_oq/run_all.pid` prints (not pkill), then start it again later.
#
#   S1  V907 Sco: Cassini states of the spin of the eclipsing pair (reference run 2000 yr + ensemble of 600)
#   S5  BU CMi: joint fit with the apsidal constants free
#   S7  CzeV343: the future with SSE against MIST tracks (stage 1, 300 + 300 draws)
#   S6  V994 Her: a hidden companion of pair C (ensemble of 3000)
#   S3  CzeV343: the 3:2 commensurability (present system 400 yr + MEGNO scan of P_A across 3:2)
#   S2  stability criteria (MA01, V22, MLP) against N-body for a 2+2 with one widened binary (ensemble of 1200)
#   S4  QX Cas: MCMC of the dynamical fit + eclipse dates of 300 posterior samples
set -u
cd "$(dirname "$0")/../.." || exit 1
ROOT=$(pwd)
[ -f .venv/bin/activate ] && source .venv/bin/activate
OUT=output/campaign_oq
mkdir -p "$OUT"
echo $$ > "$OUT/run_all.pid"
CFG=configs/campaign_oq
NPROC=$(nproc 2>/dev/null || echo 2)
J=${SWINGERSIM_CORES:-0}          # 0 = all cores
JP=$(( NPROC > 1 ? NPROC - 1 : 1 )) # tools with their own pools: leave one core for the desktop
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

echo "[$(stamp)] campaign 'open questions' in $ROOT, $NPROC cores"
python -c "import swingersim; print('SWINGERsim', swingersim.__version__)"
git log -1 --oneline 2>/dev/null

# --- S1 V907 Sco: spin ------------------------------------------------------------------------------------------
mkdir -p "$OUT/s1_v907_spin"
step s1a swingersim run "$CFG/s1a_v907_spin_reference.yaml" -o "$OUT/s1_v907_spin/reference.h5"
step s1b swingersim ensemble "$CFG/s1b_v907_spin_ensemble.yaml" -o "$OUT/s1_v907_spin/ensemble" -j "$J" --resume

# --- S5 BU CMi: k2 free --------------------------------------------------------------------------------------------
step s5 swingersim fit "$CFG/s5_bu_cmi_k2_free.yaml" --processes "$J"

# --- S7 CzeV343: SSE against MIST -----------------------------------------------------------------------------------
step s7_sse python tools/czev343_future.py --model sse --n 300 --stage2 0 --processes "$JP" -o "$OUT/s7_tracks/sse"
step s7_mist python tools/czev343_future.py --model mist --solutions fiducial --n 300 --stage2 0 \
  --processes "$JP" -o "$OUT/s7_tracks/mist"

# --- S6 V994 Her: companion of C ------------------------------------------------------------------------------------
step s6 swingersim ensemble "$CFG/s6_v994_her_c_companion.yaml" -o "$OUT/s6_v994_her_c" -j "$J" --resume

# --- S3 CzeV343: 3:2 ------------------------------------------------------------------------------------------------
mkdir -p "$OUT/s3_czev343_resonance"
step s3a swingersim run "$CFG/s3_czev343_resonance.yaml" -o "$OUT/s3_czev343_resonance/present.h5"
step s3b swingersim chaosmap "$CFG/s3_czev343_resonance.yaml" -o "$OUT/s3_czev343_resonance/megno_PA" \
  --x orbits.A.P --x-range 1.2070 1.2136 --nx 121 --t-max 400 --max-wall 3600 -j "$J" --resume

# --- S2 stability criteria ------------------------------------------------------------------------------------------
step s2 swingersim ensemble "$CFG/s2_stability_criteria_2p2.yaml" -o "$OUT/s2_stability" -j "$J" --resume

# --- S4 QX Cas: MCMC + dates ----------------------------------------------------------------------------------------
step s4_config python tools/campaign_oq/make_s4_config.py
step s4_mcmc swingersim fit "$CFG/s4_qx_cas_mcmc.yaml" --processes "$J"
step s4_dates python tools/campaign_oq/posterior_dates.py "$OUT/s4_qx_cas_mcmc" 300 --years 1900 2150 \
  --processes "$JP"

echo; echo "[$(stamp)] campaign finished"; cat "$OUT/status.txt"
rm -f "$OUT/run_all.pid"
