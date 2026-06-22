#!/usr/bin/env python3
"""Build a unified Sigma(1080) group file carrying BOTH the multi-character table
(so dym-mod-metro can run the multi-character action) AND the appended defining-rep
SU(3) matrices (so dym-mod-metro-savecfg can dump NERSC configs for the Wilson
flow), all in ONE element ordering.

Ordering + matrices come from mys1080-v4.  The full character table is computed
directly from v4's Cayley table by the Burnside-Dixon method, so it is already in
v4's conjugacy-class order (no fragile alignment).  Rows are then permuted to match
groups/char_tabs/char_table.py so the beta-vector / irrep indices agree with S1080ct.

Output (perfect_char group-file format, same header as S1080ct):
    P C
    ReChar[C][C]   (= -Re chi, negated convention used by group.c)
    ImChar[C][C]   (= -Im chi)
    mult[P][P]
    conclass[P]
    <appended>  18 reals/element: defining-rep 3x3 matrix, data[2*c + 6*r + part]

Usage: python3 generators/build_s1080_matrix_file.py [V4] [CHAR_TABLE_PY] [OUT]
"""
import sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
V4   = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "..", "improved", "groups", "mys1080-v4")
CHAR = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "groups", "char_tabs", "char_table.py")
OUT  = sys.argv[3] if len(sys.argv) > 3 else os.path.join(ROOT, "groups", "S1080ctm")
CMP  = os.path.join(ROOT, "groups", "S1080ct")

def toks(fn):
    with open(fn) as f:
        for line in f:
            for t in line.split():
                yield t

# ---------------------------------------------------------------- parse v4
print(f"reading {V4}")
it = toks(V4)
P = int(next(it))
ReTr = np.array([float(next(it)) for _ in range(P)])   # stored = -Re Tr
ImTr = np.array([float(next(it)) for _ in range(P)])   # stored = -Im Tr
mult = np.array([int(next(it)) for _ in range(P * P)], dtype=int).reshape(P, P)
rest = [float(x) for x in it]
assert len(rest) == 18 * P, f"appended reals {len(rest)} != {18*P}"
mdata = np.array(rest).reshape(P, 18)

def mat(e):
    d = mdata[e]
    return np.array([[d[2*c + 6*r] + 1j*d[2*c + 6*r + 1] for c in range(3)] for r in range(3)])

# ---------------------------------------------------------------- group structure
idn = next(n for n in range(P) if np.array_equal(mult[n], np.arange(P)))
inv = np.array([int(np.where(mult[a] == idn)[0][0]) for a in range(P)], dtype=int)
print(f"identity index = {idn}")

seen = np.zeros(P, dtype=bool)
classes = []
for g in range(P):
    if seen[g]:
        continue
    orb = sorted(set(int(mult[mult[x][g]][inv[x]]) for x in range(P)))
    for e in orb:
        seen[e] = True
    classes.append(orb)
nC = len(classes)
classof = np.empty(P, dtype=int)
for ci, cl in enumerate(classes):
    for e in cl:
        classof[e] = ci
h = np.array([len(cl) for cl in classes])           # class sizes
reps = [cl[0] for cl in classes]
i0 = int(classof[idn])
print(f"#conjugacy classes = {nC};  sizes sum = {h.sum()}")

# ---------------------------------------------------------------- Burnside-Dixon
# class-multiplication constants a[i][j][k]: # of (x,y), x in C_i, y in C_j, x*y=rep(C_k)
a = np.zeros((nC, nC, nC), dtype=np.int64)
for k in range(nC):
    z = reps[k]
    for x in range(P):
        a[classof[x], classof[mult[inv[x]][z]], k] += 1
# M_i (j,k) = a[i,j,k]; common right-eigenvectors omega satisfy M_i omega = omega_i omega
rng = np.random.default_rng(12345)
Msum = np.tensordot(rng.random(nC), a, axes=([0], [0]))     # sum_i t_i M_i  -> (j,k)
evals, evecs = np.linalg.eig(Msum)
assert evecs.shape[1] == nC
chi = np.zeros((nC, nC), dtype=complex)                     # chi[r][class]
for r in range(nC):
    w = evecs[:, r]
    w = w / w[i0]                                           # normalize identity-class component to 1 -> omega
    d = np.sqrt(P / np.sum(np.abs(w) ** 2 / h))             # dimension
    chi[r] = d * w / h
# order rows by dimension then value, for determinism
order = sorted(range(nC), key=lambda r: (round(chi[r, i0].real), tuple(np.round(chi[r].view(float), 3))))
chi = chi[order]

# ================================================================ VERIFY (computed table)
ok = True
def check(name, cond):
    global ok
    print(("  OK  " if cond else " FAIL ") + name); ok = ok and cond

dims = chi[:, i0].real
check("dimensions integer", np.max(np.abs(dims - np.round(dims))) < 1e-6)
check("sum d^2 = |G|", abs(np.sum(np.round(dims) ** 2) - P) < 1e-4)
orth = max(abs(np.sum(h * chi[r] * np.conj(chi[s])) / P - (1.0 if r == s else 0.0))
           for r in range(nC) for s in range(nC))
check("row orthonormality", orth < 1e-6)
colorth = max(abs(np.sum(chi[:, i] * np.conj(chi[:, j])) - (P / h[i] if i == j else 0.0))
              for i in range(nC) for j in range(nC))
check("column orthogonality", colorth < 1e-4)

# fundamental row = the one reproducing v4's per-element fundamental character
fund_per_class = np.array([complex(-ReTr[reps[c]], -ImTr[reps[c]]) for c in range(nC)])  # actual chi_fund
frow = min(range(nC), key=lambda r: np.max(np.abs(chi[r] - fund_per_class)))
check("a computed row equals v4 fundamental character", np.max(np.abs(chi[frow] - fund_per_class)) < 1e-6)
print(f"  (computed fundamental row = {frow}, dim={dims[frow]:.0f})")

# ---- canonical irrep (row) order: pin the indices that matter.
#   0 = trivial (1)        1 = fundamental (3, the defining rep matching the matrices)
#   2 = conj-fundamental(3) 3 = adjoint (8 = |chi_fund|^2 - 1)   then the rest by (dim, value)
trivial = min(range(nC), key=lambda r: np.max(np.abs(chi[r] - 1.0)))
fund    = frow
cfund   = min(range(nC), key=lambda r: np.max(np.abs(chi[r] - np.conj(fund_per_class))))
adj_t   = np.abs(fund_per_class) ** 2 - 1.0               # 3 (x) 3bar = 1 + adjoint
adj     = min(range(nC), key=lambda r: np.max(np.abs(chi[r] - adj_t)))
check("adjoint row = |chi_fund|^2 - 1 (dim 8)",
      np.max(np.abs(chi[adj] - adj_t)) < 1e-6 and round(dims[adj]) == 8)
head = [trivial, fund, cfund, adj]
assert len(set(head)) == 4, f"degenerate head rows {head}"
rest = sorted((r for r in range(nC) if r not in head),
              key=lambda r: (round(dims[r]), tuple(np.round(chi[r].view(float), 3))))
chi = chi[head + rest]                                    # row1=fund, row3=adjoint
check("row1=fundamental (group.c), row3=adjoint (run_file beta convention)",
      np.max(np.abs(chi[1] - fund_per_class)) < 1e-6 and np.max(np.abs(chi[3] - adj_t)) < 1e-6)
print("  row dims:", [int(round(d)) for d in chi[:, i0].real])

# negated convention; columns are in v4 class order
ReChar = -chi.real
ImChar = -chi.imag
conclass = classof.copy()                                 # element -> v4 class index = char column
check("ReChar[1][conclass[n]] == v4 ReTr[n]", np.max(np.abs(ReChar[1][conclass] - ReTr)) < 1e-6)

# cross-check: the computed table really is Sigma(1080)'s (matches char_table.py up to perm)
ns = {}; exec(open(CHAR).read(), ns); ct = np.asarray(ns["char_tab"])
def rowkey(v): return tuple(sorted((round(z.real, 3), round(z.imag, 3)) for z in v))
check("computed table == char_table.py up to row/col permutation",
      sorted(rowkey(chi[r]) for r in range(nC)) == sorted(rowkey(ct[s]) for s in range(nC)))

# matrices: SU(3) and a faithful representation of the Cayley table
import random
pr = random.Random(0)
unit = max(np.max(np.abs(mat(e).conj().T @ mat(e) - np.eye(3))) for e in range(0, P, 37))
det = max(abs(np.linalg.det(mat(e)) - 1) for e in range(0, P, 37))
homo = max(np.max(np.abs(mat(a_) @ mat(b_) - mat(int(mult[a_][b_]))))
           for a_, b_ in ((pr.randrange(P), pr.randrange(P)) for _ in range(400)))
check("matrices unitary (sample)", unit < 1e-6)
check("matrices det=1 (sample)", det < 1e-5)
check("homomorphism M(a)M(b)=M(ab) (400 samples)", homo < 1e-6)

# the irrep set (real characters) matches S1080ct -> same physics / beta semantics
if os.path.exists(CMP):
    jt = toks(CMP); Pc = int(next(jt)); Cc = int(next(jt))
    ReC = np.array([float(next(jt)) for _ in range(Cc * Cc)]).reshape(Cc, Cc)
    msa = sorted(tuple(np.round(np.sort(ReChar[r]), 5)) for r in range(nC))
    msb = sorted(tuple(np.round(np.sort(ReC[r]), 5)) for r in range(Cc))
    check("ReChar irrep-row set matches S1080ct", msa == msb)

assert ok, "VERIFICATION FAILED -- not writing file"

# ================================================================ WRITE
print(f"writing {OUT}")
def row(vals):                       # full-precision, round-trippable plain doubles
    return " ".join("%.17g" % float(x) for x in vals)
with open(OUT, "w") as f:
    f.write(f"{P} {nC}\n")
    for r in range(nC):
        f.write(row(ReChar[r]) + "\n")
    for r in range(nC):
        f.write(row(ImChar[r]) + "\n")
    for n in range(P):
        f.write(" ".join(str(int(x)) for x in mult[n]) + "\n")
    f.write(" ".join(str(int(x)) for x in conclass) + "\n")
    for e in range(P):
        f.write(row(mdata[e]) + "\n")
print("done:", OUT)
