#!/usr/bin/env python3
"""Validate the reweighting engine: reweight a BASE ensemble to the coupling of a
DIRECT ensemble and check the predicted observable matches the directly-measured one.

usage: search/validate_reweight.py BASE_DIR DIRECT_DIR
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from reweight import load_pchar, weights, ess, rw_mean_err

base   = load_pchar(os.path.join(sys.argv[1], "pchar.dat"))
direct = load_pchar(os.path.join(sys.argv[2], "pchar.dat"))
target = direct["base_beta"]

w = weights(base["P"], base["base_beta"], target)
print(f"base   beta = {base['base_beta']}")
print(f"target beta = {target}")
print(f"ESS = {ess(w):.1f} / {len(w)} configs ({100*ess(w)/len(w):.0f}% overlap)\n")
print(f"{'obs':10} {'base(beta0)':>12} {'reweighted->tgt':>18} {'direct(tgt)':>16} {'pull(sigma)':>11}")
for k in ("simpleplaq", "repoly", "impoly"):
    mr, er = rw_mean_err(base["obs"][k], w)
    d = direct["obs"][k]; md = d.mean(); ed = d.std(ddof=1) / np.sqrt(len(d))
    mb = base["obs"][k].mean()
    pull = (mr - md) / np.hypot(er, ed) if np.hypot(er, ed) > 0 else np.nan
    print(f"{k:10} {mb:12.4f} {mr:11.4f}+-{er:5.4f} {md:10.4f}+-{ed:5.4f} {pull:11.2f}")
print("\n[pass if reweighted≈direct (|pull|<~2-3) AND base differs from direct (reweight moved it)]")
