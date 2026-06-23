#!/usr/bin/env python3
"""w0/a (BMW convention) with jackknife error from an ensemble's Wilson-flow files.
Matches tools/wilsonflow/w0_scale: energy density E = Et+Es (the flow code's
f(t)=t^2(Es+Et)), and w0 = sqrt(t) where t d/dt(t^2 E) = 0.3 (the DIRECT derivative
crossing, NOT the t0/0.9596 proxy). Also reports t0/w0^2 (continuum SU(3)=0.9596) as
an SU(3)-fidelity diagnostic, and the phase (a frozen config flows to an artifact w0/a).

usage: search/w0_jack.py ENS_DIR     (ENS_DIR contains flow_<num>.dat [+ pchar.dat])
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from reweight import load_flows, load_pchar
from cost_flow import w0_from, t0_from, RATIO_T0_W0SQ

d = sys.argv[1]
phase = "?"
try:
    e = load_pchar(d.rstrip('/') + '/pchar.dat')
    sp = e['obs']['simpleplaq'].mean(); po = abs(e['obs']['repoly'].mean() + 1j*e['obs']['impoly'].mean())
    phase = "confined" if (sp < 2.3 and po < 0.1) else "FROZEN"
except Exception:
    pass

fl = load_flows(d)
nums = sorted(fl)
lab = os.path.basename(d.rstrip('/')) or d
if len(nums) < 2:
    print(f"{lab:14s} w0/a: insufficient configs ({len(nums)})"); sys.exit()
t = fl[nums[0]][:, 0]
TE = np.array([t * t * (fl[n][:, 1] + fl[n][:, 2]) for n in nums])   # E = Et+Es
mean = TE.mean(0)
w0 = w0_from(t, mean); t0 = t0_from(t, mean)
ratio = t0 / (w0 * w0) if (np.isfinite(w0) and w0 > 0 and np.isfinite(t0)) else np.nan
# jackknife error on the direct w0 (delete-1 over configs, on the mean curve)
n = len(TE); jk = []
for i in range(n):
    wj = w0_from(t, TE[np.arange(n) != i].mean(0))
    if np.isfinite(wj): jk.append(wj)
jk = np.array(jk)
err = np.sqrt((len(jk) - 1) / len(jk) * np.sum((jk - jk.mean()) ** 2)) if len(jk) > 1 else np.nan
flag = "" if phase == "confined" else f"  [{phase}: ARTIFACT]"
print(f"{lab:14s} [{phase}] w0/a = {w0:.3f} +- {err:.3f}   t0/w0^2={ratio:.3f} (SU3={RATIO_T0_W0SQ})   N={n}{flag}")
