#!/usr/bin/env python3
"""Validate that a discrete-action ensemble reproduces continuum SU(3) gradient-flow.
Two non-circular checks per point:
  (1) t0/w0^2 measured INDEPENDENTLY (t0 = value crossing t^2<E>=0.3; w0 = derivative
      crossing t d/dt(t^2<E>)=0.3) vs the continuum SU(3) value 0.9596 (FlowQCD).
  (2) overlay t^2<E> vs t/w0^2 on the digitized continuum SU(3) curve; RMS deviation
      in the physical window t/w0^2 in [0.2,1.5].
A discrete action that is genuinely SU(3)-like at this a sits on the continuum curve.

usage: validate_su3.py "label:dir1,dir2" ["label2:dir3,dir4" ...]  [--out FILE.png]
"""
import sys, os, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from reweight import load_flows
from cost_flow import t0_from, w0_from, RATIO_T0_W0SQ

ap = argparse.ArgumentParser()
ap.add_argument("specs", nargs="+", help='each "label:dir1,dir2" (dirs pooled = one action point)')
ap.add_argument("--out", default="validate_su3.png")
a = ap.parse_args()

ref = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), "continuum_t2E_su3.dat"))
plt.figure(figsize=(7.5, 5.2))
plt.plot(ref[:, 0], ref[:, 1], "k-", lw=2.5, label="continuum SU(3) (Nf=0)", zorder=10)

for spec in a.specs:
    label, dirs = spec.split(":", 1); dirs = dirs.split(",")
    curves = []
    for d in dirs:
        for num, c in load_flows(d).items():
            t = c[:, 0]; te = t * t * (c[:, 1] + c[:, 2]); curves.append(te)   # E = Et+Es (code: f=t^2(Es+Et))
    t = c[:, 0]  # common flow-time grid
    TE = np.array(curves); mean = TE.mean(0); err = TE.std(0) / np.sqrt(len(TE))
    t0 = t0_from(t, mean); w0 = w0_from(t, mean); w0sq = w0 * w0
    ratio = t0 / w0sq if (np.isfinite(t0) and np.isfinite(w0sq) and w0sq > 0) else np.nan
    x = t / w0sq
    win = (x >= 0.2) & (x <= 1.5)
    rms = np.sqrt(np.mean((mean[win] - np.interp(x[win], ref[:, 0], ref[:, 1])) ** 2))
    w0a = np.sqrt(t0 / RATIO_T0_W0SQ)
    print(f"{label:18s} N={len(curves):3d}  t0={t0:.3f} w0^2={w0sq:.3f}  "
          f"t0/w0^2={ratio:.4f} (SU3={RATIO_T0_W0SQ})  w0/a={w0a:.3f}  RMS_dev={rms:.4f}")
    plt.errorbar(x, mean, yerr=err, fmt="o", ms=3.5, capsize=2,
                 label=f"{label}: w0/a={w0a:.2f}, t0/w0^2={ratio:.3f}, RMS={rms:.3f}")

plt.axhline(0.3, ls=":", c="gray", lw=1)
plt.xlabel(r"$t / w_0^2$"); plt.ylabel(r"$t^2\langle E\rangle$")
plt.xlim(0, 1.6); plt.ylim(0, 0.45)
plt.legend(fontsize=8, loc="lower right")
plt.title("Gradient-flow validation: discrete Σ(1080) action vs continuum SU(3)")
plt.tight_layout(); plt.savefig(a.out, dpi=120)
print("saved", a.out)
