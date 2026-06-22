#!/bin/bash
# Generate perfect-action ladder ensembles on S1080ctm and save NERSC configs.
# Polite by default (nice + bounded concurrency/threads) so it can share a busy box.
# Does NOT flow -- run scripts/flow_all.sh afterwards.  Run from the repo root.
#
# usage:
#   scripts/run_ladder.sh OUTDIR "TRUNCS" "BVLIST" D NT NX "SEEDS" K N NTHERM [NCONC] [NTHREADS] [NICE]
#
#   TRUNCS   perfect-action truncations, e.g. "1t 2t 3t 4t 5t"
#   BVLIST   bare couplings,            e.g. "7.6 8.0 8.4"
#   D NT NX  lattice (D=4 required for the Wilson flow); 12^3x4 => D=4 NT=4 NX=12
#   SEEDS    e.g. "1 2 3"
#   K N NTHERM   decorrelation sweeps / #configs / thermalization sweeps
#   NCONC    concurrent points        (default 4)
#   NTHREADS OMP threads per point    (default 4)   -- keep NCONC*NTHREADS <= free cores
#   NICE     niceness                 (default 19)
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
[ $# -ge 10 ] || { sed -n '2,18p' "$0"; exit 1; }
OUT=$1; TRUNCS=$2; BVS=$3; D=$4; NT=$5; NX=$6; SEEDS=$7; K=$8; N=$9; NTHERM=${10}
NCONC=${11:-4}; NTHREADS=${12:-4}; NICE=${13:-19}

[ -x ./dym-mod-metro-savecfg ] || { echo "build first:  make"; exit 1; }
[ -f ./groups/S1080ctm ]      || { echo "missing groups/S1080ctm (run generators/build_s1080_matrix_file.py)"; exit 1; }
[ "$D" = 4 ] || echo "WARNING: D=$D; the Wilson flow needs 4D (D=4) configs."
mkdir -p "$OUT"
echo "# ladder: TRUNCS=[$TRUNCS] BV=[$BVS] D=$D Nt=$NT Nx=$NX SEEDS=[$SEEDS] K=$K N=$N Ntherm=$NTHERM"
echo "# nice=$NICE  concurrency=$NCONC x ${NTHREADS}thr -> out=$OUT"

running=0
for tr in $TRUNCS; do for bv in $BVS; do for sd in $SEEDS; do
  betas=$(BETAS_ONLY=1 python3 generators/couplings_for_ctm.py "$tr" "$bv" "$D" "$NT" "$NX") \
      || { echo "coupling generation failed for $tr bv=$bv"; exit 1; }
  pdir="$OUT/$tr/b$bv/seed$sd"; mkdir -p "$pdir"
  echo "launch  $tr  bv=$bv  seed=$sd  -> $pdir"
  nice -n "$NICE" env OMP_NUM_THREADS="$NTHREADS" \
    ./dym-mod-metro-savecfg ./groups/S1080ctm "$D" "$NT" "$NX" $betas "$sd" "$pdir/" "$K" "$N" "$NTHERM" \
    > "$pdir/run.log" 2>&1 &
  running=$((running+1))
  [ "$running" -ge "$NCONC" ] && { wait -n; running=$((running-1)); }
done; done; done
wait
echo "# all ensemble points complete -> $OUT"
echo "# next:  scripts/flow_all.sh $OUT   then  scripts/w0_table.py $OUT"
