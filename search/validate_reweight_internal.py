#!/usr/bin/env python3
"""Transition-immune validation of the reweighting ENGINE: check it against the
analytic linear response computed from the SAME ensemble (no second run, no phase
ambiguity).  W(U) ∝ exp(-sum_c beta_c P_c)  =>  d<O>/dbeta_c = -cov(O, P_c).

usage: search/validate_reweight_internal.py ENS_DIR [active_chars...]
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from reweight import load_pchar, weights, rw_mean, ess

e = load_pchar(os.path.join(sys.argv[1], "pchar.dat"))
active = [int(a) for a in sys.argv[2:]] or [1, 2]
P, O, b0 = e["P"], e["obs"]["simpleplaq"], e["base_beta"]

# (1) dbeta = 0  -> weights uniform -> reweighted mean == plain mean (exact)
w0 = weights(P, b0, b0)
print(f"[plumbing] dbeta=0:  rw_mean - plain_mean = {rw_mean(O, w0) - O.mean():.2e}  (expect 0)")

# (2) small delta on the active chars: reweighted vs analytic linear response
Pa = P[:, active].sum(axis=1)
cov = np.mean(O * Pa) - O.mean() * Pa.mean()      # connected correlation
print(f"\nactive chars {active};  susceptibility d<O>/dbeta = -cov = {-cov:.5f}")
print(f"{'delta':>8} {'reweighted':>12} {'linear pred':>12} {'diff':>10} {'~delta^2?':>10} {'ESS%':>6}")
for delta in (0.002, 0.005, 0.01, 0.02):
    b = b0.copy()
    for c in active:
        b[c] += delta
    w = weights(P, b0, b)
    orw = rw_mean(O, w)
    lin = O.mean() - delta * cov
    print(f"{delta:8.3f} {orw:12.6f} {lin:12.6f} {orw-lin:+10.2e} {(orw-lin)/delta**2:10.2f} {100*ess(w)/len(w):6.0f}")
print("\n[PASS if diff -> 0 like delta^2 (curvature), i.e. reweighting matches the analytic")
print(" linear slope as delta->0. A wrong sign would give diff ~ 2*delta*cov (linear in delta).]")
