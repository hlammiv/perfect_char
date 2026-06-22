#!/bin/bash
# Flow every config in each trajectory-point dir and extract w0 (BMW w0_scale).
# Usage: flow_w0.sh SCANDIR [EPS] [TMAX] [NC]
#   SCANDIR is the OUTDIR from deploy_scaling.sh
set -e
SCAN=$1; EPS=${2:-0.01}; TMAX=${3:-10.0}; NC=${4:-4}
HERE="$(cd "$(dirname "$0")/.." && pwd)"
WF="$HERE/tools/wilsonflow/wilson_flow"; W0="$HERE/tools/wilsonflow/w0_scale"
flow_one() {  # flow a single config, append its WFLOW block to the point's allflow.dat
  local cfg=$1 d
  d="$(dirname "$(dirname "$cfg")")"            # .../b1.._b2.. (cfg is in cfg/)
  ( cd "$(dirname "$cfg")" && OMP_NUM_THREADS=1 "$WF" -f nersc -e "$EPS" -t "$TMAX" "$(basename "$cfg")" >/dev/null 2>&1 )
  cat "$(dirname "$cfg")/flow.$(basename "$cfg")" >> "$d/allflow.dat" 2>/dev/null || true
  rm -f "$(dirname "$cfg")/flow.$(basename "$cfg")"
}
export -f flow_one; export WF EPS TMAX
for d in "$SCAN"/b1*; do [ -d "$d/cfg" ] && : > "$d/allflow.dat"; done
# flow all configs across all points, NC-parallel
find "$SCAN"/b1*/cfg -name 'nersc-*' 2>/dev/null | xargs -P "$NC" -I{} bash -c 'flow_one "$@"' _ {}
echo "# beta1 beta2 nconf  w0[a]  (BMW w0_scale)"
for d in "$SCAN"/b1*; do
  [ -s "$d/allflow.dat" ] || continue
  b1=$(basename "$d" | sed -E 's/b1([0-9.]+)_b2.*/\1/')
  b2=$(basename "$d" | sed -E 's/.*_b2(-?[0-9.]+)/\1/')
  nc=$(grep -c '^WFLOW 0\.0[0-9]* ' "$d/allflow.dat" 2>/dev/null || echo '?')
  w0=$("$W0" "$d/allflow.dat" 2>&1 | tail -1)
  echo "$b1 $b2 $nc  $w0"
done
