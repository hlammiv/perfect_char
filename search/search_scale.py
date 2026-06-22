#!/usr/bin/env python3
"""Search the multi-character couplings at a SINGLE scale (fixed active set) to
minimize the gradient-flow cost, by REWEIGHTING one base ensemble.

The action is linear in beta, so an ensemble at base beta0 (generated on the
weak-coupling / ordered side) scores any nearby beta by reweighting its per-config
character sums P_c and flow curves -- no new Monte Carlo.  Only the user-specified
active characters are varied; the rest stay at their base value.  An ESS floor keeps
the optimizer inside the overlap region (and inside one phase -- reweighting cannot
cross the freezing transition).

usage:
  search/search_scale.py ENS_DIR --active "1 2 3" [--ref FILE] [--essmin 0.1]
                         [--window 0.2 1.5] [--maxstep 0.05]

ENS_DIR must contain pchar.dat (from dym-mod-metro-savecfg) and flow_<num>.dat
(from search/flow_configs.sh).
"""
import sys, os, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from reweight import load_pchar, load_flows, align_flows, weights, ess
from cost_flow import flow_cost
from scipy.optimize import minimize

ap = argparse.ArgumentParser()
ap.add_argument("ensdir")
ap.add_argument("--active", required=True, help='space-separated character indices, e.g. "1 2 3"')
ap.add_argument("--ref", default=None, help="continuum reference curve file (cols: t/w0^2  t^2E)")
ap.add_argument("--essmin", type=float, default=0.1, help="minimum ESS fraction (overlap floor)")
ap.add_argument("--window", nargs=2, type=float, default=[0.2, 1.5])
ap.add_argument("--maxstep", type=float, default=0.05, help="max |dbeta| per active char")
a = ap.parse_args()
active = [int(x) for x in a.active.split()]
win = tuple(a.window)

e = load_pchar(os.path.join(a.ensdir, "pchar.dat"))
flows = load_flows(a.ensdir)
if not flows:
    sys.exit("no flow_<num>.dat in ENS_DIR -- run search/flow_configs.sh first")
t, E1, E2, idx = align_flows(e["nums"], flows)
P, base = e["P"][idx], e["base_beta"]
ref = None
if a.ref:
    r = np.loadtxt(a.ref); ref = (r[:, 0], r[:, 1])
else:
    print("# no --ref: using self-contained anisotropy proxy (t^2E1 - t^2E2)^2")

def evaluate(beta):
    w = weights(P, base, beta)
    esf = ess(w) / len(w)
    if esf < a.essmin:
        return np.inf, esf, np.nan
    c, w0 = flow_cost(t, E1, E2, w, ref, win)
    return c, esf, w0

def objective(x):
    beta = base.copy()
    for k, c in enumerate(active):
        beta[c] = base[c] + np.clip(x[k], -a.maxstep, a.maxstep)
    return evaluate(beta)[0]

c0, es0, w00 = evaluate(base)
print(f"base beta      : {np.round(base,4)}")
print(f"active chars   : {active}")
print(f"base cost      : {c0:.6e}   (ESS {100*es0:.0f}%, w0/a={w00:.4f})\n")

res = minimize(objective, np.zeros(len(active)), method="Nelder-Mead",
               options={"xatol": 1e-3, "fatol": 1e-8, "maxiter": 3000})
beta = base.copy()
for k, c in enumerate(active):
    beta[c] = base[c] + np.clip(res.x[k], -a.maxstep, a.maxstep)
c, es, w0 = evaluate(beta)

print(f"best beta      : {np.round(beta,5)}")
print(f"best cost      : {c:.6e}   (was {c0:.6e},  {100*(1-c/c0):.0f}% lower)")
print(f"ESS at optimum : {100*es:.0f}%   w0/a={w0:.4f}")
at_edge = any(abs(abs(res.x[k]) - a.maxstep) < 1e-3 for k in range(len(active)))
if es <= a.essmin*1.2 or at_edge:
    print("NOTE: optimum hit the overlap/step boundary -> re-anchor: generate a new base"
          " ensemble at this beta and search again.")
