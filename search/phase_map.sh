#!/bin/bash
# Gen-only phase map (confined vs frozen) — no Wilson flow, fast. For each point:
# generate a short ensemble, classify by <plaq> and |Polyakov|. Detached:
#   setsid bash search/phase_map.sh OUTDIR POINTSFILE D NT NX N K NTHERM NCONC THREADS </dev/null >OUT.log 2>&1 &
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
OUT=$1; PF=$2; D=$3; NT=$4; NX=$5; N=$6; K=$7; NTH=$8; NCONC=$9; THR=${10}
mkdir -p "$OUT"
one() {
  local label=$1; local seed=$2; shift 2; local betas="$*"
  local d="$OUT/$label"; mkdir -p "$d"
  env OMP_NUM_THREADS="$THR" ./dym-mod-metro-savecfg ./groups/S1080ctm "$D" "$NT" "$NX" \
      $betas "$seed" "$d/" "$K" "$N" "$NTH" > "$d/gen.log" 2>&1
  rm -f "$d"/nersc-*          # phase map: configs not needed
  python3 - "$d" "$label" <<'PY' > "$d/phase.txt" 2>/dev/null
import sys; sys.path.insert(0,'search'); import numpy as np
from reweight import load_pchar
d, lab = sys.argv[1], sys.argv[2]
e = load_pchar(d+'/pchar.dat'); sp = e['obs']['simpleplaq'].mean()
po = abs(e['obs']['repoly'].mean()+1j*e['obs']['impoly'].mean())
ph = 'confined' if (sp < 2.3 and po < 0.1) else 'FROZEN'
print('%-22s plaq=%.3f poly=%.3f %s' % (lab, sp, po, ph))
PY
}
export -f one; export OUT D NT NX N K NTH THR
running=0
while read -r line; do
  [ -z "$line" ] && continue; case "$line" in \#*) continue;; esac
  bash -c "one $line" &
  running=$((running+1)); [ "$running" -ge "$NCONC" ] && { wait -n; running=$((running-1)); }
done < "$PF"
wait
echo "=== phase map ($OUT, ${D}D ${NT}x${NX}) ==="; cat "$OUT"/*/phase.txt 2>/dev/null | sort
touch "$OUT.DONE"
