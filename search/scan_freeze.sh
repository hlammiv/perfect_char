#!/bin/bash
# Map the FREEZING line in (beta_fund, beta_X) at fixed beta_adj, to see whether
# character X pushes freezing OUT (=> finer accessible lattice spacing before freezing).
# Generation-only (no flow): the confined<->frozen transition is a sharp plaquette jump,
# so <plaq> and |Polyakov| identify the phase cheaply. Higher beta_fund still-confined
# = finer a reachable.  fund -> chars 1,2; adjoint -> char 3.
#
# usage: search/scan_freeze.sh OUTDIR CHAR "BF_LIST" "BX_LIST" BADJ D NT NX N K NTHERM [threads]
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
[ $# -ge 11 ] || { sed -n '2,11p' "$0"; exit 1; }
OUT=$1; CH=$2; BFS=$3; BXS=$4; BADJ=$5; D=$6; NT=$7; NX=$8; N=$9; K=${10}; NTH=${11}; THR=${12:-16}
mkdir -p "$OUT"
echo "# char X = $CH (extra), beta_adj = $BADJ"
echo "# beta_fund  beta_X     <plaq>    |poly|    phase"
for bf in $BFS; do for bx in $BXS; do
  betas=$(python3 -c "v=['0']*17; v[1]='$bf'; v[2]='$bf'; v[3]='$BADJ'; v[$CH]='$bx'; print(' '.join(v))")
  d="$OUT/bf${bf}_x${bx}"; mkdir -p "$d"
  nice -n 5 env OMP_NUM_THREADS="$THR" ./dym-mod-metro-savecfg ./groups/S1080ctm "$D" "$NT" "$NX" \
       $betas 7 "$d/" "$K" "$N" "$NTH" >"$d/run.log" 2>&1
  python3 - "$d" "$bf" "$bx" <<'PY'
import sys; sys.path.insert(0,'search'); import numpy as np
from reweight import load_pchar
d,bf,bx=sys.argv[1],sys.argv[2],sys.argv[3]
e=load_pchar(d+'/pchar.dat'); sp=e['obs']['simpleplaq']
poly=abs(e['obs']['repoly'].mean()+1j*e['obs']['impoly'].mean())
ph='FROZEN' if sp.mean()>2.3 else ('confined' if poly<0.1 else 'ordered')
print('  %-9s %-9s %8.4f %8.4f   %s'%(bf,bx,sp.mean(),poly,ph))
PY
done; done
echo "# done -> $OUT"