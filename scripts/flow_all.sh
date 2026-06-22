#!/bin/bash
# Wilson-flow every saved NERSC config under OUTDIR and collect a per-point flow.dat.
# Run from the repo root, after scripts/run_ladder.sh.
#
# usage: scripts/flow_all.sh OUTDIR [EPS] [TMAX] [NCONC]
#   EPS    flow step    (default 0.01)
#   TMAX   max flow time(default 10.0; lower it for small lattices to save time)
#   NCONC  concurrent flows (default 6)   -- each flow is single-threaded
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
[ $# -ge 1 ] || { sed -n '2,9p' "$0"; exit 1; }
OUT=$1; EPS=${2:-0.01}; TMAX=${3:-10.0}; NCONC=${4:-6}
WF="$ROOT/tools/wilsonflow/wilson_flow"
[ -x "$WF" ] || { echo "build the flow tools:  ( cd tools/wilsonflow && ./build.sh )"; exit 1; }

find "$OUT" -name flow.dat -delete 2>/dev/null
flow_one() {
  local cfg=$1 d b
  d=$(dirname "$cfg"); b=$(basename "$cfg")
  ( cd "$d" && OMP_NUM_THREADS=1 "$WF" -f nersc -e "$EPS" -t "$TMAX" "$b" >/dev/null 2>&1 )
  cat "$d/flow.$b" >> "$d/flow.dat" 2>/dev/null && rm -f "$d/flow.$b"
}
export -f flow_one; export WF EPS TMAX
n=$(find "$OUT" -name 'nersc-*' | wc -l)
echo "# flowing $n configs (eps=$EPS tmax=$TMAX) ${NCONC}-way ..."
find "$OUT" -name 'nersc-*' | xargs -P "$NCONC" -I{} bash -c 'flow_one "$@"' _ {}
echo "# done.  summarize:  scripts/w0_table.py $OUT"
