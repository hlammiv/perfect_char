#!/bin/bash
# 2D scan over (beta_fund, beta_adjoint) to soften the first-order freezing
# transition and locate a single-phase scaling window (w0/a ~ 1-2, smooth plaquette).
# fund -> chars 1,2 (3 + 3bar);  adjoint -> char 3 (the 8).
# A SMOOTH variation of <plaq> across beta_fund at fixed beta_adj (no jump) marks the
# softened/crossover region near the transition endpoint.
#
# usage: search/scan_plane.sh OUTDIR "BF_LIST" "BA_LIST" D NT NX N K NTHERM [eps tmax nthreads]
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
[ $# -ge 9 ] || { sed -n '2,12p' "$0"; exit 1; }
OUT=$1; BFS=$2; BAS=$3; D=$4; NT=$5; NX=$6; N=$7; K=$8; NTH=$9
EPS=${10:-0.04}; TMAX=${11:-8.0}; THR=${12:-6}
mkdir -p "$OUT"
echo "# beta_fund  beta_adj   <plaq>    |poly|     w0/a"
for bf in $BFS; do for ba in $BAS; do
  betas=$(python3 -c "v=['0']*17; v[1]='$bf'; v[2]='$bf'; v[3]='$ba'; print(' '.join(v))")
  d="$OUT/bf${bf}_ba${ba}"; mkdir -p "$d"
  nice -n 19 env OMP_NUM_THREADS="$THR" ./dym-mod-metro-savecfg ./groups/S1080ctm "$D" "$NT" "$NX" \
       $betas 7 "$d/" "$K" "$N" "$NTH" >"$d/run.log" 2>&1
  nice -n 19 bash search/flow_configs.sh "$d" "$EPS" "$TMAX" "$THR" >/dev/null 2>&1
  python3 - "$d" "$bf" "$ba" <<'PY'
import sys; sys.path.insert(0, 'search')
import numpy as np
from reweight import load_pchar, load_flows, align_flows
from cost_flow import reweighted_t2E, w0_from
d, bf, ba = sys.argv[1], sys.argv[2], sys.argv[3]
e = load_pchar(d + '/pchar.dat'); fl = load_flows(d)
plaq = e['obs']['simpleplaq'].mean()
poly = abs(e['obs']['repoly'].mean() + 1j*e['obs']['impoly'].mean())
if fl:
    t, E1, E2, idx = align_flows(e['nums'], fl)
    t2E1, t2E2 = reweighted_t2E(t, E1, E2, np.ones(len(idx)))
    w0 = w0_from(t, 0.5*(t2E1 + t2E2))
else:
    w0 = float('nan')
print('  %-9s %-9s %8.4f %8.4f %8.4f' % (bf, ba, plaq, poly, w0))
PY
done; done
echo "# done -> $OUT"
