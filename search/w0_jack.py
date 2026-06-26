#!/usr/bin/env python3
"""w0/a and t0/a (BMW convention) with jackknife error from an ensemble's flow files.
Matches tools/wilsonflow/w0_scale: energy density E = Et+Es (the flow code's
f(t)=t^2(Es+Et)); w0 = sqrt(t) where t d/dt(t^2 E)=0.3; sqrt(t0)/a = sqrt(t) where
t^2 E = 0.3. Reports BOTH dimensionless scales, t0/w0^2 (continuum SU(3)=0.9596) as an
SU(3)-fidelity diagnostic, and the phase (a frozen config flows to an artifact scale).

SCALE-SETTING (read this): pure SU(3) Yang-Mills has NO intrinsic 'fm'. Convert with a
PURE-GAUGE reference: a[fm] = sqrt(t0)_phys / (sqrt(t0)/a), with sqrt(t0)_phys ~ 0.145-0.166 fm
(convention-dependent), OR a = w0_phys/(w0/a) with the PURE-GAUGE w0. Do NOT use the
2+1f-QCD value w0=0.1755 fm (that was an earlier bug). The literature (Alexandru et al)
sets the scale with t0 -- prefer t0/a for comparison.

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
st0 = np.sqrt(t0) if (np.isfinite(t0) and t0 > 0) else np.nan   # sqrt(t0)/a (the literature's scale)
ratio = t0 / (w0 * w0) if (np.isfinite(w0) and w0 > 0 and np.isfinite(t0)) else np.nan
# jackknife error on the direct w0 (delete-1 over configs, on the mean curve)
n = len(TE); jk = []
for i in range(n):
    wj = w0_from(t, TE[np.arange(n) != i].mean(0))
    if np.isfinite(wj): jk.append(wj)
jk = np.array(jk)
err = np.sqrt((len(jk) - 1) / len(jk) * np.sum((jk - jk.mean()) ** 2)) if len(jk) > 1 else np.nan
flag = "" if phase == "confined" else f"  [{phase}: ARTIFACT]"
print(f"{lab:14s} [{phase}] w0/a={w0:.3f}+-{err:.3f} sqrt(t0)/a={st0:.3f} t0/w0^2={ratio:.3f} (SU3={RATIO_T0_W0SQ})   N={n}{flag}")
