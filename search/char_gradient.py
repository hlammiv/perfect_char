#!/usr/bin/env python3
"""Rank characters by their (within-ensemble, reweighted) flow-cost gradient.

The action is linear in beta, so from ONE ensemble we can reweight to beta +/- delta*e_c
for every character c and measure d(cost)/d(beta_c) cheaply and with LOW noise (same
configs, same t0 baseline -> the cost DIFFERENCE is correlated, unlike re-anchoring).
The characters with the steepest cost-reducing gradient are the ones worth turning on.

usage: search/char_gradient.py ENS_DIR [--ref FILE] [--delta 0.02] [--window 0.2 1.5]
"""
import sys, os, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from reweight import load_pchar, load_flows, align_flows, weights, ess
from cost_flow import reweighted_t2E, t0_from, RATIO_T0_W0SQ

ap = argparse.ArgumentParser()
ap.add_argument("ensdir")
ap.add_argument("--ref", default="search/continuum_t2E_su3.dat")
ap.add_argument("--delta", type=float, default=0.02)
ap.add_argument("--window", nargs=2, type=float, default=[0.2, 1.5])
a = ap.parse_args()
win = tuple(a.window)
ref = None
if a.ref and os.path.exists(a.ref):
    r = np.loadtxt(a.ref); ref = (r[:, 0], r[:, 1])

e = load_pchar(os.path.join(a.ensdir, "pchar.dat")); fl = load_flows(a.ensdir)
t, E1, E2, idx = align_flows(e["nums"], fl)
P, base = e["P"][idx], e["base_beta"]
C = e["C"]

def cost(b):
    w = weights(P, base, b); esf = ess(w) / len(w)
    t2E1, t2E2 = reweighted_t2E(t, E1, E2, w); te = 0.5 * (t2E1 + t2E2); t0 = t0_from(t, te)
    if not np.isfinite(t0) or t0 <= 0: return np.nan, esf
    x = t / t0; sel = (x >= win[0]) & (x <= win[1])
    if sel.sum() < 3: return np.nan, esf
    d = te[sel] - np.interp(RATIO_T0_W0SQ * x[sel], ref[0], ref[1]) if ref is not None else t2E1[sel] - t2E2[sel]
    return float(np.trapezoid(d * d, x[sel]) / (win[1] - win[0])), esf

# representative dims per row (S1080ctm): row->dim
dims = {0:1,1:3,2:3,3:8,4:3,5:3,6:5,7:5,8:6,9:6,10:8,11:9,12:9,13:9,14:10,15:15,16:15}
c0, es0 = cost(base)
# per-character step sized to a FIXED overlap budget: delta_c = target / std(P_c)
# (so each character is moved by the same ESS cost; then rank by actual cost reduction)
sig = P.std(axis=0)
print(f"ensemble: {a.ensdir}   base cost = {c0:.4e}   (overlap budget per step ~{a.delta})")
print(f"{'char':>4} {'dim':>4} {'delta_c':>9} {'best dcost':>12} {'dir':>5} {'ess':>5}   note")
rows = []
for c in range(1, C):
    if sig[c] < 1e-9:
        continue
    dc = a.delta / sig[c]
    bp = base.copy(); bp[c] += dc; cp, esp = cost(bp)
    bm = base.copy(); bm[c] -= dc; cm, esm = cost(bm)
    cand = [(cp, "+", esp), (cm, "-", esm)]
    cand = [x for x in cand if np.isfinite(x[0])]
    if not cand:
        continue
    cbest, dirn, esb = min(cand, key=lambda x: x[0])
    rows.append((c, dims[c], dc, cbest - c0, dirn, esb))
rows.sort(key=lambda r: r[3])           # most negative dcost (biggest reduction) first
for c, dim, dc, dcost, dirn, esb in rows:
    note = "lowers cost" if dcost < 0 else "(at optimum / raises)"
    print(f"{c:>4} {dim:>4} {dc:>9.4f} {dcost:>+12.3e} {dirn:>5} {100*esb:>4.0f}%   {note}")
print("\n[most-negative dcost at fixed overlap = the character to turn on next]")
