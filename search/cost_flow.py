"""Gradient-flow cost functions for the perfect-action search.

Given an ensemble's per-config flow curves (E1=gact4i, E2=gactij) and reweighting
weights w (to a trial coupling), form the reweighted t^2<E>(t) curve, set the scale
w0, and score how far the action is from continuum SU(3):

  - cost_ref : squared deviation of t^2<E>(t/w0^2) from a supplied continuum SU(3)
               reference curve, over a flow-time window.  (the Holland-style target;
               requires a reference, --ref)
  - cost_aniso : self-contained proxy = integrated (t^2 E1 - t^2 E2)^2 over the window;
               the electric-magnetic (an)isotropy of the flowed action density, which
               a perfect action / isotropic continuum drives to zero.  (default)
"""
import numpy as np

REF_W0 = 0.3      # reference value defining w0:  t d/dt(t^2<E>) = 0.3


def reweighted_t2E(t, E1, E2, w):
    """Reweighted t^2<E1>(t), t^2<E2>(t) with weights w (aligned to E rows)."""
    W = w / w.sum()
    return t*t*(W @ E1), t*t*(W @ E2)


def w0_from(t, t2E, ref=REF_W0):
    """w0 from W(t)=t d/dt(t^2<E>)=ref; returns sqrt(t_ref) or nan."""
    W = t[1:-1]*(t2E[2:] - t2E[:-2])/(t[2:] - t[:-2])
    tt = t[1:-1]
    for i in range(len(W)-1):
        if (W[i]-ref)*(W[i+1]-ref) < 0:
            tc = tt[i] + (ref-W[i])/(W[i+1]-W[i])*(tt[i+1]-tt[i])
            return np.sqrt(tc) if tc > 0 else np.nan
    return np.nan


RATIO_T0_W0SQ = 0.9596     # pure SU(3) t0/w0^2 (FlowQCD 1503.06516): converts t/t0 <-> t/w0^2

def t0_from(t, t2E, ref=REF_W0):
    """t0 = flow time where t^2<E>=ref (a VALUE crossing -- far more stable than the
    w0 derivative crossing, especially near finite-volume saturation)."""
    for i in range(len(t2E) - 1):
        if (t2E[i] - ref) * (t2E[i+1] - ref) < 0:
            return t[i] + (ref - t2E[i]) / (t2E[i+1] - t2E[i]) * (t[i+1] - t[i])
    return np.nan


def cost_aniso(t, t2E1, t2E2, w0, window=(0.2, 1.5)):
    if not np.isfinite(w0) or w0 <= 0:
        return np.inf
    x = t/w0**2
    sel = (x >= window[0]) & (x <= window[1])
    if sel.sum() < 3:
        return np.inf
    d = t2E1[sel] - t2E2[sel]
    return float(np.trapz(d*d, x[sel]) / (window[1]-window[0]))


def cost_ref(t, t2E1, t2E2, w0, ref_xy, window=(0.2, 1.5)):
    if not np.isfinite(w0) or w0 <= 0:
        return np.inf
    x = t/w0**2
    sel = (x >= window[0]) & (x <= window[1])
    if sel.sum() < 3:
        return np.inf
    cont = np.interp(x[sel], ref_xy[0], ref_xy[1])
    d = 0.5*(t2E1[sel] + t2E2[sel]) - cont
    return float(np.trapz(d*d, x[sel]) / (window[1]-window[0]))


def flow_cost(t, E1, E2, w, ref_xy=None, window=(0.2, 1.5)):
    """Full cost at given weights: returns (cost, w0).  Uses cost_ref if a reference
    curve is supplied, else the self-contained cost_aniso proxy."""
    t2E1, t2E2 = reweighted_t2E(t, E1, E2, w)
    w0 = w0_from(t, 0.5*(t2E1 + t2E2))
    if ref_xy is not None:
        return cost_ref(t, t2E1, t2E2, w0, ref_xy, window), w0
    return cost_aniso(t, t2E1, t2E2, w0, window), w0
