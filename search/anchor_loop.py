#!/usr/bin/env python3
"""Re-anchor loop: walk the multi-character couplings down the gradient-flow cost.

Each iteration: generate a base ensemble at the current beta (single phase), flow it,
reweight-search the active characters within the overlap (ESS) ball, step beta to the
best, re-anchor, repeat. Stops if a step pushes into the frozen phase, or when the cost
stops improving. Logs cost & beta vs iteration.

Run from the repo root.
usage:
  search/anchor_loop.py OUTDIR --beta "0 2.7 2.7 -1.0 0 ..." --active "1 3" [--iters 6]
     [--D 4 --nt 8 --nx 8 --N 60 --K 40 --ntherm 1500 --eps 0.04 --tmax 8.0 --threads 16]
     [--ref search/continuum_t2E_su3.dat --window 0.2 1.0 --essmin 0.25 --maxstep 0.15 --seed 4242]

Note: chars 1,2 are the conjugate 3+3bar (identical real characters), so only beta1+beta2
enters the action -- put char 1 (not 2) in --active and leave beta2 fixed to avoid the flat
direction; the effective fundamental coupling is then beta1+beta2.
"""
import sys, os, argparse, subprocess
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from reweight import load_pchar, load_flows, align_flows, weights, ess
from cost_flow import reweighted_t2E, t0_from, RATIO_T0_W0SQ
from scipy.optimize import minimize

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GROUP = "./groups/S1080ctm"

ap = argparse.ArgumentParser()
ap.add_argument("outdir")
ap.add_argument("--beta", required=True)
ap.add_argument("--active", required=True)
ap.add_argument("--iters", type=int, default=6)
ap.add_argument("--D", type=int, default=4); ap.add_argument("--nt", type=int, default=8); ap.add_argument("--nx", type=int, default=8)
ap.add_argument("--N", type=int, default=60); ap.add_argument("--K", type=int, default=40); ap.add_argument("--ntherm", type=int, default=1500)
ap.add_argument("--eps", type=float, default=0.04); ap.add_argument("--tmax", type=float, default=8.0); ap.add_argument("--threads", type=int, default=16)
ap.add_argument("--ref", default="search/continuum_t2E_su3.dat")
ap.add_argument("--window", nargs=2, type=float, default=[0.2, 1.0])
ap.add_argument("--essmin", type=float, default=0.25); ap.add_argument("--maxstep", type=float, default=0.15)
ap.add_argument("--seed", type=int, default=4242)
a = ap.parse_args()
os.chdir(ROOT)
beta = np.array([float(x) for x in a.beta.split()]); C = len(beta)
active = [int(x) for x in a.active.split()]
win = tuple(a.window)
ref = None
if a.ref and os.path.exists(a.ref):
    r = np.loadtxt(a.ref); ref = (r[:, 0], r[:, 1])
os.makedirs(a.outdir, exist_ok=True)
logf = open(os.path.join(a.outdir, "loop.log"), "w")
def log(*m):
    s = " ".join(str(x) for x in m); print(s, flush=True); logf.write(s + "\n"); logf.flush()

def gen_and_flow(d, b, seed):
    os.makedirs(d, exist_ok=True)
    cmd = ["./dym-mod-metro-savecfg", GROUP, str(a.D), str(a.nt), str(a.nx)] + \
          ["%.10g" % x for x in b] + [str(seed), d + "/", str(a.K), str(a.N), str(a.ntherm)]
    env = dict(os.environ, OMP_NUM_THREADS=str(a.threads))
    subprocess.run(cmd, env=env, stdout=open(d + "/gen.log", "w"), stderr=subprocess.STDOUT, check=True)
    subprocess.run(["bash", "search/flow_configs.sh", d, str(a.eps), str(a.tmax), str(a.threads)],
                   stdout=open(d + "/flow.log", "w"), stderr=subprocess.STDOUT, check=True)

def load(d):
    e = load_pchar(d + "/pchar.dat"); fl = load_flows(d)
    t, E1, E2, idx = align_flows(e["nums"], fl)
    return e, t, E1, E2, e["P"][idx], e["base_beta"], e["obs"]

def cost_at(P, base, t, E1, E2, b):
    w = weights(P, base, b); esf = ess(w) / len(w)
    if esf < a.essmin: return np.inf, esf, np.nan
    t2E1, t2E2 = reweighted_t2E(t, E1, E2, w); te = 0.5 * (t2E1 + t2E2); t0 = t0_from(t, te)
    if not np.isfinite(t0) or t0 <= 0: return np.inf, esf, np.nan
    x = t / t0; sel = (x >= win[0]) & (x <= win[1])           # scale by t0 (stable value crossing)
    if sel.sum() < 3: return np.inf, esf, t0
    if ref is None:
        d = t2E1[sel] - t2E2[sel]
    else:
        d = te[sel] - np.interp(RATIO_T0_W0SQ * x[sel], ref[0], ref[1])  # ref is in t/w0^2 = R*(t/t0)
    return float(np.trapezoid(d * d, x[sel]) / (win[1] - win[0])), esf, t0

log("# iter  eff_fund(b1+b2)  beta_adj(b3)  w0/a   base_cost   best_cost   ess%   status")
prev_cost = None
for it in range(a.iters):
    d = os.path.join(a.outdir, f"iter{it:02d}")
    gen_and_flow(d, beta, a.seed + it)
    e, t, E1, E2, P, base, obs = load(d)
    sp = obs["simpleplaq"]
    if sp.mean() > 2.3 or sp.std() > 0.06:
        log(f"  {it:3d}  STOP: ensemble not single confined phase (plaq mean={sp.mean():.3f} std={sp.std():.3f}) -> beta hit freezing")
        break
    c0, es0, t0s = cost_at(P, base, t, E1, E2, base)
    def obj(z):
        b = base.copy()
        for k, c in enumerate(active): b[c] = base[c] + np.clip(z[k], -a.maxstep, a.maxstep)
        return cost_at(P, base, t, E1, E2, b)[0]
    res = minimize(obj, np.zeros(len(active)), method="Nelder-Mead",
                   options={"xatol": 1e-3, "fatol": 1e-10, "maxiter": 4000})
    best = base.copy()
    for k, c in enumerate(active): best[c] = base[c] + np.clip(res.x[k], -a.maxstep, a.maxstep)
    cb, esb, t0b = cost_at(P, base, t, E1, E2, best)
    w0disp = (t0s / RATIO_T0_W0SQ) ** 0.5 if np.isfinite(t0s) else float("nan")
    log(f"  {it:3d}    {base[1]+base[2]:.4f}        {base[3]:+.4f}    {w0disp:.3f}  {c0:.4e}  {cb:.4e}  {100*esb:.0f}   ok")
    log(f"        best beta active {active}: {[round(best[c],4) for c in active]}")
    beta = best
    if prev_cost is not None and cb > 0.98 * prev_cost:
        log(f"  converged: cost not improving ({cb:.4e} vs prev {prev_cost:.4e})"); break
    prev_cost = cb
log("# done")
