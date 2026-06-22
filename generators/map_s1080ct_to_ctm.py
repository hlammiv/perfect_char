#!/usr/bin/env python3
"""Resolve the irrep-row correspondence between S1080ct and S1080ctm.

couplings_Nt.py emit a beta-vector in S1080ct's irrep-row order; S1080ctm uses a
different row order (and a different element/class labeling).  To apply the
perfect-action couplings to the *correct* characters on S1080ctm we need the row
permutation pi with  ReChar_ct[s][c] == ReChar_ctm[pi[s]][sigma[c]]  for all s,c,
where sigma is the class correspondence.

Found by constructing the group isomorphism phi: S1080ct -> S1080ctm from the two
Cayley tables, ANCHORED on the defining (fundamental, row 1) representation so the
isomorphism respects the SU(3) content (rules out outer/Galois automorphisms that
would send a coupling to the wrong same-dimension irrep).  The result is verified
by full character-table equality and written to generators/s1080ct_to_ctm.json.
"""
import sys, os, json
import numpy as np
from collections import Counter, deque

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CT  = os.path.join(ROOT, "groups", "S1080ct")
CTM = os.path.join(ROOT, "groups", "S1080ctm")
OUT = os.path.join(HERE, "s1080ct_to_ctm.json")

def toks(fn):
    with open(fn) as f:
        for line in f:
            for t in line.split():
                yield t

def load(fn):
    it = toks(fn); P = int(next(it)); C = int(next(it))
    Re = np.array([[float(next(it)) for _ in range(C)] for _ in range(C)])
    Im = np.array([[float(next(it)) for _ in range(C)] for _ in range(C)])
    mult = np.array([int(next(it)) for _ in range(P*P)], dtype=np.int32).reshape(P, P)
    conclass = np.array([int(next(it)) for _ in range(P)], dtype=np.int32)
    return P, C, Re, Im, mult, conclass

def derive(P, Re, Im, mult, conclass):
    idn = next(n for n in range(P) if np.array_equal(mult[n], np.arange(P)))
    order = np.ones(P, dtype=int)
    for e in range(P):
        x, o = e, 1
        while x != idn:
            x = int(mult[x][e]); o += 1
        order[e] = o
    cnt = Counter(int(c) for c in conclass)
    csize = np.array([cnt[int(conclass[e])] for e in range(P)])
    fund = np.array([complex(-Re[1][conclass[e]], -Im[1][conclass[e]]) for e in range(P)])
    return idn, order, csize, fund

print("loading files ...")
Pa, C, ReA, ImA, multA, ccA = load(CT)
Pb, Cb, ReB, ImB, multB, ccB = load(CTM)
assert Pa == Pb and C == Cb, "size mismatch"
idA, ordA, szA, fundA = derive(Pa, ReA, ImA, multA, ccA)
idB, ordB, szB, fundB = derive(Pb, ReB, ImB, multB, ccB)

def keyset(order, sz, fund):
    return [(int(order[e]), int(sz[e]), round(fund[e].real, 4), round(fund[e].imag, 4)) for e in range(Pa)]

def closure(gens, mult, idn):
    seen = {idn}; dq = deque([idn])
    while dq:
        a = dq.popleft()
        for g in gens:
            c = int(mult[a][g])
            if c not in seen:
                seen.add(c); dq.append(c)
    return seen

def find_iso(conj):
    """Try to find phi: ct->ctm preserving (order,csize,fund); conjugate fund_a if conj."""
    fA = np.conj(fundA) if conj else fundA
    keyA = keyset(ordA, szA, fA)
    keyB = keyset(ordB, szB, fundB)
    bucket = {}
    for e in range(Pb):
        bucket.setdefault(keyB[e], []).append(e)
    # candidate count per ct element
    cand = [len(bucket.get(keyA[e], [])) for e in range(Pa)]
    if any(c == 0 for e, c in enumerate(cand) if ordA[e] > 1):
        return None  # some element has no image under this anchor -> wrong conj
    # greedy generating set, preferring small candidate buckets
    order_by_cand = sorted((e for e in range(Pa) if ordA[e] > 1), key=lambda e: cand[e])
    gens, cur = [], {idA}
    for e in order_by_cand:
        if e in cur:
            continue
        gens.append(e); cur = closure(gens, multA, idA)
        if len(cur) == Pa:
            break
    assert len(cur) == Pa, "generators do not generate the group"
    buckets = [bucket[keyA[g]] for g in gens]
    print(f"  conj={conj}: generators={gens}  bucket sizes={[len(b) for b in buckets]}")

    import itertools
    phi = np.empty(Pa, dtype=np.int32)
    attempts = 0
    for images in itertools.product(*buckets):
        attempts += 1
        phi[:] = -1; phi[idA] = idB
        dq = deque([idA]); ok = True
        while dq and ok:
            a = dq.popleft(); pa = int(phi[a])
            for g, im in zip(gens, images):
                c = int(multA[a][g]); pc = int(multB[pa][im])
                if phi[c] != -1:
                    if phi[c] != pc: ok = False; break
                else:
                    if keyB[pc] != keyA[c]: ok = False; break
                    phi[c] = pc; dq.append(c)
        if ok and np.all(phi >= 0) and len(set(phi.tolist())) == Pa:
            print(f"  found isomorphism after {attempts} generator-image combos")
            return phi
    return None

phi = find_iso(conj=False)
if phi is None:
    phi = find_iso(conj=True)
assert phi is not None, "no SU(3)-respecting isomorphism found"

# class map sigma: ct class -> ctm class
sigma = np.empty(C, dtype=int)
for c in range(C):
    e = int(np.where(ccA == c)[0][0])
    sigma[c] = int(ccB[int(phi[e])])
assert len(set(sigma.tolist())) == C, "sigma not a bijection"

# row map pi: ct row s -> ctm row with matching character (columns aligned by sigma)
pi = np.full(C, -1, dtype=int)
for s in range(C):
    matches = [r for r in range(C)
               if np.max(np.abs(ReB[r][sigma] - ReA[s])) < 1e-6
               and np.max(np.abs(ImB[r][sigma] - ImA[s])) < 1e-6]
    assert len(matches) == 1, f"row {s}: {len(matches)} matches"
    pi[s] = matches[0]
assert len(set(pi.tolist())) == C, "pi not a bijection"

# full verification
maxerr = max(max(np.max(np.abs(ReA[s] - ReB[pi[s]][sigma])),
                 np.max(np.abs(ImA[s] - ImB[pi[s]][sigma]))) for s in range(C))
print(f"full char-table match error = {maxerr:.2e}")
assert maxerr < 1e-6

# dims via |chi(id)|: identity class of ct
cidA = int(ccA[idA])
dims_ct = [int(round(abs(ReA[s][cidA]))) for s in range(C)]
print("\n S1080ct row -> S1080ctm row   (dim)")
for s in range(C):
    print(f"   {s:2d} -> {pi[s]:2d}    (dim {dims_ct[s]})")

json.dump({"pi_ct_to_ctm": pi.tolist(), "sigma_class_ct_to_ctm": sigma.tolist(),
           "dims_ct_row": dims_ct},
          open(OUT, "w"), indent=1)
print(f"\nwrote {OUT}")
