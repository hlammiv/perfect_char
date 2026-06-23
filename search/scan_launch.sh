#!/bin/bash
# Generic detached scan launcher for lenore. Waits for any in-flight flows, runs a
# points scan + phase-aware analysis, signals done. Run detached:
#   setsid bash search/scan_launch.sh OUTDIR POINTSFILE [NT N K NTHERM NCONC THREADS] </dev/null >/dev/null 2>&1 &
cd ~/perfect_char || exit 1
OUT=$1; PF=$2; NT=${3:-12}; N=${4:-60}; K=${5:-40}; NTH=${6:-1500}; NCONC=${7:-8}; THR=${8:-4}
while pgrep -x wilson_flow >/dev/null 2>&1; do sleep 15; done
rm -f "$OUT.DONE"
bash search/run_points.sh "$OUT" "$PF" 4 "$NT" "$NT" "$N" "$K" "$NTH" "$NCONC" "$THR" > "$OUT.log" 2>&1
python3 search/analyze_scan.py "$OUT" >> "$OUT.log" 2>&1
touch "$OUT.DONE"
