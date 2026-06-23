#!/usr/bin/env python3
"""FREE reweighting screen: d(w0/a)/d(beta_c) at an anchor, per REDUCED character,
with jackknife error. CAVEAT: this is the response at FIXED beta_f -- it does NOT
capture the freezing-edge-pushing benefit (the main way characters help). Known
edge-pushers (adjoint, 6) are expected to show NEGATIVE direct gradient here.
Use this only to find DIRECT improvers; a null is not 'no character helps'.

Reduced characters = rows of ReChar that are identical (conjugate pairs); only the
sum enters the action, so we vary one representative per group.

usage: w0_gradient.py ENS_DIR [--group GROUPFILE] [--budget 0.3]
"""
import sys, os, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from reweight import load_pchar, load_flows, align_flows, weights, ess
from cost_flow import reweighted_t2E, t0_from, RATIO_T0_W0SQ

ap = argparse.ArgumentParser()
ap.add_argument("ens")
ap.add_argument("--group", default=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "groups", "S1080ctm"))
ap.add_argument("--budget", type=float, default=0.3)
a = ap.parse_args()

def load_rechar(fn):
    it = (x for line in open(fn) for x in line.split())
    int(next(it)); C = int(next(it))
    return np.array([[float(next(it)) for _ in range(C)] for _ in range(C)])

Re = load_rechar(a.group)
e = load_pchar(a.ens + "/pchar.dat"); fl = load_flows(a.ens)
t, E1, E2, idx = align_flows(e["nums"], fl); P = e["P"][idx]; base = e["base_beta"]; C = e["C"]
dims = {0:1,1:3,2:3,3:8,4:3,5:3,6:5,7:5,8:6,9:6,10:8,11:9,12:9,13:9,14:10,15:15,16:15}
ON = {1,2,3,8,9}   # already on: fund(1,2), adjoint(3), 6(8,9)

# reduced groups: rows with identical Re chi
groups = []; used = set()
for r in range(1, C):
    if r in used: continue
    g = [r] + [r2 for r2 in range(r+1, C) if r2 not in used and np.max(np.abs(Re[r]-Re[r2])) < 1e-6]
    used.update(g); groups.append(g)

def w0a_jack(b):
    w = weights(P, base, b)
    te = 0.5*np.array(reweighted_t2E(t, E1, E2, w)).sum(0); t0 = t0_from(t, te)
    if not np.isfinite(t0): return np.nan, np.nan, ess(w)/len(w)
    n = len(idx); jk = []
    for i in range(n):
        k = np.arange(n) != i
        ww = weights(P[k], base, b); tt = 0.5*np.array(reweighted_t2E(t, E1[k], E2[k], ww)).sum(0); t0k = t0_from(t, tt)
        if np.isfinite(t0k): jk.append((t0k/RATIO_T0_W0SQ)**0.5)
    jk = np.array(jk)
    err = np.sqrt((len(jk)-1)/len(jk)*np.sum((jk-jk.mean())**2)) if len(jk) > 1 else np.nan
    return (t0/RATIO_T0_W0SQ)**0.5, err, ess(w)/len(w)

sig = P.std(0)
w0b, eb, _ = w0a_jack(base)
print(f"anchor {a.ens}:  w0/a = {w0b:.3f} +- {eb:.3f}   eff_fund={base[1]+base[2]:.2f} adj={base[3]:.2f} 6={base[8]:.2f}")
print(f"{'group':14} {'dim':>3} {'d(w0/a)/dbeta':>16} {'ESS':>5}  tag  signif")
rows = []
for g in groups:
    rep = g[0]
    if sig[rep] < 1e-9: continue
    d = a.budget / sig[rep]
    bp = base.copy(); bp[rep] += d; wp, ep, esp = w0a_jack(bp)
    bm = base.copy(); bm[rep] -= d; wm, em, esm = w0a_jack(bm)
    if not (np.isfinite(wp) and np.isfinite(wm)): continue
    grad = (wp - wm) / (2*d); gerr = np.hypot(ep, em) / (2*d)
    rows.append((g, dims[rep], grad, gerr, min(esp, esm), rep not in ON))
rows.sort(key=lambda r: -r[2])
for g, dim, grad, gerr, es, isnew in rows:
    tag = "NEW" if isnew else "on "
    s = "++(>2sig)" if grad-2*gerr > 0 else ("--(<-2sig)" if grad+2*gerr < 0 else "~0")
    print(f"{str(g):14} {dim:>3} {grad:>+9.4f} +-{gerr:6.4f} {100*es:>4.0f}%  {tag}  {s}")
print("\n[direct improver = NEW with ++(>2sig). Reminder: edge-pushers (which is how adj & 6 help)")
print(" show -- here; their benefit needs the freezing-edge scan, not this fixed-beta_f gradient.]")
