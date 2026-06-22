#!/bin/bash
# Scan the fundamental coupling to locate the scaling window (w0/a ~ O(1), single
# phase below freezing). Per beta: generate a small ensemble, flow it, report
# <plaq>, |Polyakov|, and w0/a. Look for w0/a ~ 1-2 with intermediate <plaq>
# (not ~3 = frozen, not small = strong coupling).
#
# usage: search/scan_window.sh OUTDIR "BV_LIST" D NT NX N K NTHERM [eps tmax nthreads]
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
[ $# -ge 8 ] || { sed -n '2,10p' "$0"; exit 1; }
OUT=$1; BVS=$2; D=$3; NT=$4; NX=$5; N=$6; K=$7; NTH=$8
EPS=${9:-0.04}; TMAX=${10:-8.0}; THR=${11:-6}
mkdir -p "$OUT"
echo "# beta_fund   <plaq>    |poly|     w0/a"
for bv in $BVS; do
  betas=$(python3 -c "print('0 $bv $bv'+' 0'*14)")
  d="$OUT/b$bv"; mkdir -p "$d"
  nice -n 19 env OMP_NUM_THREADS="$THR" ./dym-mod-metro-savecfg ./groups/S1080ctm "$D" "$NT" "$NX" \
       $betas 7 "$d/" "$K" "$N" "$NTH" >"$d/run.log" 2>&1
  nice -n 19 bash search/flow_configs.sh "$d" "$EPS" "$TMAX" "$THR" >/dev/null 2>&1
  python3 - "$d" "$bv" <<'PY'
import sys; sys.path.insert(0, 'search')
import numpy as np
from reweight import load_pchar, load_flows, align_flows
from cost_flow import reweighted_t2E, w0_from
d, bv = sys.argv[1], sys.argv[2]
e = load_pchar(d + '/pchar.dat'); fl = load_flows(d)
plaq = e['obs']['simpleplaq'].mean()
poly = abs(e['obs']['repoly'].mean() + 1j*e['obs']['impoly'].mean())
if fl:
    t, E1, E2, idx = align_flows(e['nums'], fl)
    t2E1, t2E2 = reweighted_t2E(t, E1, E2, np.ones(len(idx)))
    w0 = w0_from(t, 0.5*(t2E1 + t2E2))
else:
    w0 = float('nan')
print('  %-9s %8.4f %8.4f %8.4f' % (bv, plaq, poly, w0))
PY
done
echo "# done -> $OUT"
