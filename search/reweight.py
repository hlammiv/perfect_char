"""Reweighting engine for the multi-character action.

The action is LINEAR in the couplings,  S(U) = -sum_c beta_c P_c(U),  where
P_c(U) = sum_plaq Re chi_c(U_p) is logged per config by dym-mod-metro-savecfg
(pchar.dat).  Hence the Monte Carlo weight is exp(+sum_c beta_c P_c) ... up to
the sign convention fixed by the update: the sampled weight is

    W(U) ∝ exp( - sum_c beta_c P_c(U) ),

so observables at any beta' reweight EXACTLY from an ensemble at beta0:

    <O>_{beta'} = sum_n O_n w_n / sum_n w_n,   w_n = exp( -sum_c (beta'_c-beta0_c) P_c[n] ).

Validity is limited only by ensemble overlap (effective sample size).
"""
import os, glob, re
import numpy as np


def load_flows(dirpath):
    """Per-config Wilson-flow curves: {num: array[(t, gact4i, gactij), ...]}.
    Flow files are flow_<num>.dat (search/flow_configs.sh)."""
    flows = {}
    for f in glob.glob(os.path.join(dirpath, "flow_*.dat")):
        num = int(re.search(r"flow_(\d+)\.dat", f).group(1))
        rows = [tuple(map(float, l.split()[1:4])) for l in open(f) if l.startswith("WFLOW")]
        if rows:
            flows[num] = np.array(rows)
    return flows


def align_flows(nums, flows):
    """Return (t_grid, E1[nconf,nt], E2[nconf,nt], idx) aligned to the pchar rows in `nums`.
    E1=gact4i, E2=gactij; idx selects the P/obs rows that have a flow file (same order)."""
    keys = [int(n) for n in nums if int(n) in flows]
    t = flows[keys[0]][:, 0]
    E1 = np.array([flows[n][:, 1] for n in keys])
    E2 = np.array([flows[n][:, 2] for n in keys])
    idx = np.array([i for i, n in enumerate(nums) if int(n) in flows])
    return t, E1, E2, idx


def load_pchar(path):
    """Parse a pchar.dat -> dict(base_beta, C, nums, P[n,C], obs{simpleplaq,repoly,impoly})."""
    base = None; C = None; rows = []
    for l in open(path):
        if l.startswith("# BASE_BETA"):
            base = np.array([float(x) for x in l.split()[2:]])
        elif l.startswith("# NCHAR"):
            C = int(l.split()[2])
        elif l.startswith("#"):
            continue
        else:
            rows.append([float(x) for x in l.split()])
    rows = np.array(rows)
    return {"base_beta": base, "C": C, "nums": rows[:, 0].astype(int),
            "P": rows[:, 1:1+C],
            "obs": {"simpleplaq": rows[:, 1+C], "repoly": rows[:, 2+C], "impoly": rows[:, 3+C]}}


def weights(P, base_beta, beta):
    """Reweighting weights from base_beta to beta (numerically stabilized)."""
    db = np.asarray(beta, float) - np.asarray(base_beta, float)
    x = -(P @ db)                      # log weight = -sum_c dbeta_c P_c
    x -= x.max()
    return np.exp(x)


def ess(w):
    """Effective sample size (Kish): (sum w)^2 / sum w^2.  ESS/N near 1 = good overlap."""
    return (w.sum() ** 2) / np.sum(w ** 2)


def rw_mean(vals, w):
    return np.sum(w * vals) / np.sum(w)


def rw_mean_err(vals, w):
    """Reweighted mean with delete-1 jackknife error over configs."""
    vals = np.asarray(vals); n = len(vals)
    m = rw_mean(vals, w)
    jk = np.array([rw_mean(np.delete(vals, i), np.delete(w, i)) for i in range(n)])
    return m, np.sqrt((n - 1) / n * np.sum((jk - m) ** 2))
