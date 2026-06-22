#!/usr/bin/env python3
"""Summarize a perfect_char ladder run: Wilson-flow scales (w0/a, sqrt(t0)/a) and
the mean fundamental plaquette, per (truncation, bare coupling), averaged over seeds.

Walks  OUTDIR/<trunc>/b<bv>/seed<seed>/{flow.dat, run.log}.

usage: scripts/w0_table.py OUTDIR [ref]      (ref = flow reference value, default 0.3)
"""
import sys, os, glob, re
import numpy as np

if len(sys.argv) < 2:
    print(__doc__); sys.exit(1)
OUT = sys.argv[1]
REF = float(sys.argv[2]) if len(sys.argv) > 2 else 0.3

def split_configs(fn):
    rows = []
    for l in open(fn):
        if l.startswith("WFLOW"):
            p = l.split(); rows.append((float(p[1]), float(p[2]), float(p[3])))
    if not rows:
        return []
    rows = np.array(rows); cfgs = []; start = 0
    for i in range(1, len(rows)):
        if rows[i, 0] < rows[i-1, 0]:        # flow time reset -> new config
            cfgs.append(rows[start:i]); start = i
    cfgs.append(rows[start:]); return cfgs

def cross(x, y, ref):
    for i in range(len(y)-1):
        if (y[i]-ref)*(y[i+1]-ref) < 0:
            return x[i] + (ref-y[i])/(y[i+1]-y[i])*(x[i+1]-x[i])
    return np.nan

def scales(cfg):                              # t0 from t^2E=ref, w0 from t d/dt(t^2E)=ref
    t = cfg[:, 0]; E = (cfg[:, 1]+cfg[:, 2])/2.0; t2E = t*t*E
    W = t[1:-1]*(t2E[2:]-t2E[:-2])/(t[2:]-t[:-2])
    t0 = cross(t, t2E, REF); w0sq = cross(t[1:-1], W, REF)
    return (np.sqrt(t0) if np.isfinite(t0) and t0 > 0 else np.nan,
            np.sqrt(w0sq) if np.isfinite(w0sq) and w0sq > 0 else np.nan)

def mean_se(vals):
    v = np.array([x for x in vals if np.isfinite(x)])
    if len(v) == 0: return (np.nan, np.nan)
    if len(v) == 1: return (v[0], np.nan)
    return v.mean(), v.std(ddof=1)/np.sqrt(len(v))

def plaq_from_log(fn):
    vals = []
    try:
        for l in open(fn):
            if l.startswith("GMES:"):
                vals.append(float(l.split()[4]))     # simpleplaq column
    except FileNotFoundError:
        pass
    return np.mean(vals) if vals else np.nan

pts = {}
for fd in sorted(glob.glob(os.path.join(OUT, "*", "b*", "seed*", "flow.dat"))):
    m = re.search(r"/([^/]+)/b([^/]+)/seed([^/]+)/flow\.dat$", fd)
    if not m: continue
    tr, bv = m.group(1), m.group(2)
    cfgs = split_configs(fd)
    w0, _ = mean_se([scales(c)[1] for c in cfgs])     # per-point w0 (mean over its configs)
    t0, _ = mean_se([scales(c)[0] for c in cfgs])
    plq = plaq_from_log(os.path.join(os.path.dirname(fd), "run.log"))
    pts.setdefault((tr, bv), []).append((w0, t0, plq))

print(f"# ref={REF}  OUTDIR={OUT}")
print(f"# {'trunc':6} {'bv':>7} {'nseed':>5} {'w0/a':>9} {'+-':>8} {'sqrt(t0)/a':>11} {'<plaq>':>9}")
for (tr, bv), rows in sorted(pts.items()):
    mw, ew = mean_se([r[0] for r in rows])
    mt, _  = mean_se([r[1] for r in rows])
    plq    = np.nanmean([r[2] for r in rows])
    print(f"  {tr:6} {bv:>7} {len(rows):5d} {mw:9.4f} {ew:8.4f} {mt:11.4f} {plq:9.4f}")
if not pts:
    print("# no flow.dat found -- run scripts/flow_all.sh first")
