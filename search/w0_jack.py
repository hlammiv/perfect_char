#!/usr/bin/env python3
"""w0/a with a statistical error from an ensemble's per-config Wilson-flow files.
Per config: t0 (t^2<E>=0.3 crossing) -> w0/a = sqrt(t0/0.9596); then mean +/- std/sqrt(N).
The per-config spread also flags finite-volume noise (large spread = flow hitting the box).

usage: search/w0_jack.py ENS_DIR     (ENS_DIR contains flow_<num>.dat)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from reweight import load_flows
from cost_flow import t0_from, RATIO_T0_W0SQ

d = sys.argv[1]
fl = load_flows(d)
w0 = []
for num, c in fl.items():
    t = c[:, 0]; te = 0.5 * t * t * (c[:, 1] + c[:, 2])
    t0 = t0_from(t, te)
    if np.isfinite(t0) and t0 > 0:
        w0.append((t0 / RATIO_T0_W0SQ) ** 0.5)
w0 = np.array(w0); n = len(w0)
lab = os.path.basename(d.rstrip('/')) or d
if n < 2:
    print(f"{lab:14s} w0/a: insufficient resolved configs ({n})")
else:
    m = w0.mean(); err = w0.std(ddof=1) / np.sqrt(n)
    print(f"{lab:14s} w0/a = {m:.3f} +- {err:.3f}   (N={n} cfgs, per-cfg spread {100*w0.std()/m:.0f}%)")
