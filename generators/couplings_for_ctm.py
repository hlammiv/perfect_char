#!/usr/bin/env python3
"""Emit a perfect-action multi-character beta-vector for groups/S1080ctm.

couplings_Nt.py produce the character-expansion (perfect-action) couplings in
S1080ct's irrep-row order.  S1080ctm (which also carries the SU(3) matrices needed
for the Wilson flow) uses a different row order, so the couplings must be permuted
by pi (generators/s1080ct_to_ctm.json, from map_s1080ct_to_ctm.py) before they are
applied -- otherwise a coupling lands on the wrong same-dimension irrep.

Usage:
    python3 generators/couplings_for_ctm.py TRUNC bv D Nt Nx [seed] [K N Ntherm]
      TRUNC = 1t | 2t | 3t | 4t | 5t   (which perfect-action truncation)
      bv    = bare coupling passed to couplings_<TRUNC>.py
prints the beta-vector and a ready dym-mod-metro-savecfg command for S1080ctm.
"""
import sys, os, io, json, contextlib

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

if len(sys.argv) < 6:
    print(__doc__); sys.exit(1)
trunc = sys.argv[1].lower().lstrip("couplings_").rstrip(".py")
bv, D, Nt, Nx = sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
seed = sys.argv[6] if len(sys.argv) > 6 else "12345"
K   = sys.argv[7] if len(sys.argv) > 7 else "200"
N   = sys.argv[8] if len(sys.argv) > 8 else "100"
NTH = sys.argv[9] if len(sys.argv) > 9 else "2000"

cfile = os.path.join(HERE, f"couplings_{trunc}.py")
if not os.path.exists(cfile):
    sys.exit(f"no such truncation: {cfile}")

# run couplings_<trunc>.py to obtain g[] and mask[] (S1080ct row order); suppress its print
ns = {}
saved = sys.argv
sys.argv = [cfile, bv, D, Nt, Nx, "S1080", "data/"]
try:
    with contextlib.redirect_stdout(io.StringIO()):
        exec(open(cfile).read(), ns)
finally:
    sys.argv = saved
g, mask = ns["g"], ns["mask"]
beta_ct = [float(mask[i]) * float(g[i]) for i in range(len(g))]

pi = json.load(open(os.path.join(HERE, "s1080ct_to_ctm.json")))["pi_ct_to_ctm"]
C = len(pi)
assert len(beta_ct) == C, f"coupling vector length {len(beta_ct)} != {C}"
beta_ctm = [0.0] * C
for s in range(C):
    beta_ctm[pi[s]] = beta_ct[s]

betastr = " ".join("%.10g" % b for b in beta_ctm)
if os.environ.get("BETAS_ONLY"):      # machine-readable: just the S1080ctm beta-vector
    print(betastr); sys.exit(0)
print(f"# perfect-action truncation couplings_{trunc}, bv={bv}")
print(f"# beta (S1080ct row order):  {['%.5g'%b for b in beta_ct]}")
print(f"# beta (S1080ctm row order): {['%.5g'%b for b in beta_ctm]}")
print()
print("# multi-character run (no config save):")
print(f"./dym-mod-metro ./groups/S1080ctm {D} {Nt} {Nx} {betastr} {seed}")
print()
print("# config-saving run for Wilson flow (outprefix data/, K N Ntherm):")
print(f"./dym-mod-metro-savecfg ./groups/S1080ctm {D} {Nt} {Nx} {betastr} {seed} data/ {K} {N} {NTH}")
