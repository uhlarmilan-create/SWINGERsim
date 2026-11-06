#!/usr/bin/env bash
# Capture into the 3:2 resonance: the grid of tools/capture/capture_grid.py (3 epsilon x 5 migration rates x 5
# eccentricities x 6 phase sets + 6 divergent control runs = 456 runs, about 11 million years of N-body in total,
# about 30 h on 12 cores). Restartable: started again, it skips the finished runs.
#
#   cd ~/Public/SWINGERsim && git pull && nohup tools/capture/run_grid.sh > output/capture_grid.log 2>&1 &
#
# Progress: tail -f output/capture_grid.log. Stop: kill the PID in output/capture/grid/run.pid (not pkill).
set -u
cd "$(dirname "$0")/../.." || exit 1
[ -f .venv/bin/activate ] && source .venv/bin/activate
OUT=output/capture/grid
mkdir -p "$OUT"
echo $$ > "$OUT/run.pid"
NPROC=$(python -c "import os; print(os.cpu_count() or 2)")  # not nproc: it follows OMP_NUM_THREADS
JP=$(( NPROC > 1 ? NPROC - 1 : 1 ))
export PYTHONWARNINGS=ignore OMP_NUM_THREADS=1
echo "[$(date '+%Y-%m-%d %H:%M:%S')] capture grid on $JP processes"; git log -1 --oneline 2>/dev/null
python tools/capture/capture_grid.py "$OUT" -j "$JP"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] capture grid finished"
rm -f "$OUT/run.pid"
