#!/bin/bash
# Phase 1: parallel bisection of the freezing boundary beta_adj*(beta_f) for a set of beta_f.
# Each beta_f runs its own sequential bisection chain; chains run concurrently. Detached:
#   setsid bash search/bisect_all.sh OUTDIR NT N K NTHERM THREADS </dev/null >OUT.log 2>&1 &
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
OUT=$1; NT=$2; N=$3; K=$4; NTH=$5; THR=$6
mkdir -p "$OUT"
# per-beta_f starting adjoint (near the rough boundary from the 8^4 phase map; 8,10 start frozen and expand)
declare -A ADJ0=( [3]=-1.25 [4]=-1.75 [6]=-2.25 [8]=-3.0 [10]=-3.0 )
for bf in 3 4 6 8 10; do
  python3 search/bisect_boundary.py "$bf" "$OUT/bf$bf" 4 "$NT" "$NT" "$N" "$K" "$NTH" "$THR" 1111 "${ADJ0[$bf]}" 0.5 0.06 \
      > "$OUT/bf$bf.log" 2>&1 &
done
wait
echo "=== freezing boundary beta_adj*(beta_f), pure fund+adj ==="
grep -h "^BOUNDARY" "$OUT"/bf*.log | sort -t= -k1
touch "$OUT.DONE"
