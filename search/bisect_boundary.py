#!/usr/bin/env python3
"""Bisection search for the freezing boundary beta_adj*(beta_f) in pure fund+adj.
For a fixed beta_f: locate the beta_adj where the system transitions confined<->frozen.
Confined = more-negative adjoint; frozen = less-negative. Expands beta_adj more negative
WITHOUT CAP (safety stop -20) until a confined point is found, then bisects to tol.
Gen-only (plaq/poly phase), hot start, Ntherm for metastability.

usage: bisect_boundary.py BF OUTDIR D NT NX N K NTHERM THREADS SEED ADJ0 STEP TOL
"""
import sys, os, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from reweight import load_pchar

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
a = sys.argv
bf = float(a[1]); OUT = a[2]; D, NT, NX = a[3], a[4], a[5]; N, K, NTH = a[6], a[7], a[8]
THR = a[9]; SEED = a[10]; adj0 = float(a[11]); step = float(a[12]); tol = float(a[13])
os.makedirs(OUT, exist_ok=True)

def phase(adj):
    d = os.path.join(OUT, f"a{adj:.3f}"); os.makedirs(d, exist_ok=True)
    betas = ['0'] * 17; betas[1] = betas[2] = f"{bf}"; betas[3] = f"{adj:.4f}"
    env = dict(os.environ, OMP_NUM_THREADS=THR)
    with open(d + "/gen.log", "w") as log:
        subprocess.run([os.path.join(ROOT, "dym-mod-metro-savecfg"), os.path.join(ROOT, "groups/S1080ctm"),
                        D, NT, NX, *betas, SEED, d + "/", K, N, NTH], cwd=ROOT, env=env, stdout=log, stderr=log)
    for f in os.listdir(d):
        if f.startswith("nersc"):
            try: os.remove(os.path.join(d, f))
            except OSError: pass
    try:
        e = load_pchar(d + "/pchar.dat"); sp = float(e['obs']['simpleplaq'].mean())
        po = abs(e['obs']['repoly'].mean() + 1j * e['obs']['impoly'].mean())
        return ('confined' if (sp < 2.3 and po < 0.1) else 'FROZEN'), sp, po
    except Exception:
        return 'ERR', float('nan'), float('nan')

def ev(adj):
    ph, sp, po = phase(adj)
    print(f"  bf={bf:>4} adj={adj:+.3f} -> {ph:8s} (plaq={sp:.3f} poly={po:.3f})", flush=True)
    return ph

# establish bracket: conf (more negative, confined) and froz (less negative, frozen)
conf = froz = None
if ev(adj0) == 'confined':
    conf = adj0; x = adj0
    while x + step < -1e-6:
        x += step
        if ev(x) == 'FROZEN': froz = x; break
else:
    froz = adj0; x = adj0
    while x - step >= -20.0:           # no cap (safety -20)
        x -= step
        if ev(x) == 'confined': conf = x; break

if conf is None:
    print(f"BOUNDARY bf={bf}: NO confined point down to adj=-20 (freezes for all adjoint)"); sys.exit()
if froz is None:
    print(f"BOUNDARY bf={bf}: confined up to adj~0 (no freezing) — boundary above {conf:+.3f}"); sys.exit()

lo, hi = conf, froz   # lo more-negative (confined), hi less-negative (frozen)
while abs(hi - lo) > tol:
    mid = 0.5 * (lo + hi)
    if ev(mid) == 'confined': lo = mid
    else: hi = mid
print(f"BOUNDARY bf={bf}: beta_adj* = {0.5*(lo+hi):+.3f}  (confined adj<={lo:+.3f}, frozen adj>={hi:+.3f})")
