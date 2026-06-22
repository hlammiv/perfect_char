#!/bin/bash
# Wilson-flow each saved config in a base-ensemble dir to a PER-CONFIG file
# flow_<num>.dat (keyed by the config index), so search/search_scale.py can pair
# each flow curve with its pchar.dat row for reweighting.
#
# usage: search/flow_configs.sh ENS_DIR [eps] [tmax] [nconc]
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
[ $# -ge 1 ] || { sed -n '2,8p' "$0"; exit 1; }
DIR=$1; EPS=${2:-0.02}; TMAX=${3:-6.0}; NCONC=${4:-6}
WF="$ROOT/tools/wilsonflow/wilson_flow"
[ -x "$WF" ] || { echo "build the flow tools:  ( cd tools/wilsonflow && ./build.sh )"; exit 1; }

flow_one() {
  local cfg=$1 d b num
  d=$(dirname "$cfg"); b=$(basename "$cfg")
  num=$(echo "$b" | sed -E 's/.*num([0-9]+)$/\1/')
  ( cd "$d" && OMP_NUM_THREADS=1 "$WF" -f nersc -e "$EPS" -t "$TMAX" "$b" >/dev/null 2>&1 )
  mv "$d/flow.$b" "$d/flow_${num}.dat" 2>/dev/null
}
export -f flow_one; export WF EPS TMAX
find "$DIR" -maxdepth 1 -name 'nersc-*' | xargs -P "$NCONC" -I{} bash -c 'flow_one "$@"' _ {}
echo "flowed $(ls "$DIR"/flow_*.dat 2>/dev/null | wc -l) configs in $DIR (eps=$EPS tmax=$TMAX)"
