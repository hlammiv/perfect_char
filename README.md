# perfect_char

**Toward a perfect lattice action for SU(N) gauge theory using a discrete subgroup and multiple characters.**

This repository contains a lattice gauge-theory Monte Carlo in which the gauge
group is a *finite subgroup* of SU(N) (e.g. the 1080-element subgroup Σ(1080) ⊂
SU(3)) and the single-plaquette action is a sum over **multiple group characters
(irreducible representations)**, each with its own independent coupling:

```
S = - Σ_plaquettes  Σ_c  β_c · Re χ_c(U_plaquette)          (Boltzmann weight e^{-S})
```

where `c` runs over the conjugacy-class characters of the finite group. A single
character (the fundamental) reproduces the ordinary Wilson action restricted to a
discrete group, which suffers a weak-coupling **freezing transition**. The
research question here is whether *tuning several character couplings `β_c`* (a
truncated character expansion of a perfect / fixed-point action) lets a fixed
discrete subgroup approximate continuum SU(N) over a much wider range — i.e. a
**perfect action realized on a discrete group**.

See [`LITERATURE_REVIEW.md`](LITERATURE_REVIEW.md) for the background survey
(perfect/fixed-point actions, discrete-subgroup Monte Carlo, the freezing
problem, character expansions, mixed/multi-representation actions, and the
quantum-simulation motivation).

---

## Quick start

```bash
make                       # builds ./dym-mod-metro (needs g++ with OpenMP)

# Fundamental-only Wilson run on the 108-element group, 2D 4x4 lattice:
#   args: groupfile  D  Nt  Nx  beta_0 ... beta_{C-1}  [seed]
./dym-mod-metro ./groups/S108ct 2 4 4   0 0.8 0 0 0 0 0 0 0 0 0 0 0 0   12345

# Two-character (fundamental + adjoint) action on Sigma(1080), 4D 32^4:
./dym-mod-metro ./groups/S1080ct 4 32 32  0 8.407 0 -1.65 0 0 0 0 0 0 0 0 0 0 0 0 0  14931
```

The number of `beta` arguments **must equal the number of conjugacy classes `C`**
of the group (the first line of the group file is `P C`). Couplings left at `0`
are skipped at run time, so a multi-character action is just a vector with a few
non-zero entries.

> In the Σ(1080) example, `β[1]=8.407` couples the 3-dimensional **fundamental**
> irrep and `β[3]=-1.65` couples the 8-dimensional **adjoint** irrep — a concrete
> fundamental+adjoint mixed action (cf. the Bhanot–Creutz plane).

---

## Repository layout

```
perfect_char/
├── README.md, LITERATURE_REVIEW.md
├── Makefile
├── dym-mod-metro.cpp        # main: multi-character Metropolis Monte Carlo
├── dym-mod-metro-savecfg.cpp# same MC, also dumps SU(3) configs in NERSC format (for Wilson flow)
├── lattice.c / lattice.h    # geometry, plaquette, Wilson loop, Polyakov line
├── group.c   / group.h      # group-file loader (mult table + character table)
├── timer.cpp / timer.h      # lightweight timers
├── run_file                 # example list of run commands
├── now_that_we_dont_talk.sh # build + launch helper
├── groups/                  # gauge-group data (the "data needed to run")
│   ├── S1080ct              #   Σ(1080) ⊂ SU(3): 1080 elements, 17 irreps   ← main target
│   ├── S1080ctm             #   Σ(1080) + appended defining-rep matrices (for NERSC dump/flow)
│   ├── S108ct              #   order-108 subgroup: 108 elements, 14 irreps  ← fast for testing
│   ├── char_tabs/          #   character tables + coupling result dumps (.nb derivations: ../chars)
│   └── character-expansion/#   standalone character-expansion utility
├── generators/              # scripts that produce the data above
│   ├── genS108ct.py        #   build a group file in the loader's format (template)
│   ├── build_s1080_matrix_file.py  # build S1080ctm (char table + matrices) via Burnside–Dixon
│   ├── couplings_*.py      #   emit run commands with perfect-action couplings β_c(β_bare)
│   └── group_matrices/     #   per-subgroup matrix + multiplication-table generators
├── tools/wilsonflow/        # standalone BMW gradient-flow + w0 scale tools (build.sh)
├── scripts/                 # flow_w0.sh, w0_extract.py, string_tension.py, ...
└── data_s1080/              # example output log + analysis scripts
```

---

## The group-file format

`group.c` (`load_group`) reads a single text file. Whitespace is ignored; tokens
appear in this order:

```
P  C                       # group order P, number of conjugacy classes C
ReChar[0][0..C-1]          # real part of the character table, C×C
...                        #   (row = irrep, column = conjugacy class)
ReChar[C-1][0..C-1]
ImChar[0][0..C-1]          # imaginary part, C×C
...
ImChar[C-1][0..C-1]
mult[0][0..P-1]            # group multiplication table, P×P  (entries are element indices)
...
mult[P-1][0..P-1]
conclass[0..P-1]           # conjugacy class of each element (0..C-1)
```

Total = `2 + 2·C² + P² + P` tokens. The identity and inverses are found
automatically from the multiplication table.

**Sign convention.** The character entries are stored *negated* (`-χ_c`), matching
the generators (`genS108ct.py` prints `-chartab`); the Metropolis ratio in
`dym-mod-metro.cpp` undoes the sign, so the Boltzmann weight is
`exp(+Σ β_c Re χ_c)`. Keep this convention when adding new groups.

**⚠ Hard-coded fundamental row.** `group.c` defines the per-element trace used for
the plaquette/Polyakov observables as `ReTr[n] = ReChar[1][conclass[n]]` — i.e. it
assumes **irrep row 1 is the fundamental**. This is true for `S1080ct`. For
`S108ct` the fundamental is **row 4** (see the commented alternative in
`group.c:72`). The multi-character *action* uses every row and is unaffected, but
the reported `plaq`/Polyakov observables are only physical if row 1 is the
fundamental — adjust `group.c` (or reorder the irreps) for other groups.

### Two formats in the wild
- `*ct` files (`S1080ct`, `S108ct`) — the **character-table** format above; these
  are what the current code loads.
- The `generators/group_matrices/gen*.py` scripts emit a **legacy single-trace**
  format (`P`, then per-element Re/Im trace, then the multiplication table — no
  class count, no full character table). Those define the subgroup matrices and
  multiplication law but must be combined with a character table (see
  `groups/char_tabs/`) to produce a loadable `*ct` file.

---

## Output format

Each measurement sweep prints:

```
GMES: 999  <Re Polyakov>  <Im Polyakov>  <plaq>     # plaq = mean fundamental plaquette
CC m WL k:  v1 v2 ... vNt                            # character-m Wilson loop, spatial size k+1, temporal 1..Nt
...
ACC: <acceptance rate>
CONFIGS: <V·D element indices>                       # full configuration dump
```

`is_it_over_now.sh` (in `data_s1080/`) averages the plaquette from `GMES` lines;
`say_dont_go.sh` reshapes the `CC … WL …` character-Wilson-loop lines for fitting.

Decorrelation/measurement counts are compile-time `#define K` / `#define N` at the
top of `dym-mod-metro.cpp` (default `K=1000`, `N=10000`); lower them for quick
tests and rebuild.

---

## Data-generation pipeline

```
 (1) subgroup matrices + multiplication table     generators/group_matrices/gen*.py
            │
            ▼
 (2) conjugacy classes + character table          groups/char_tabs/char_table.py
            │                                       ../chars/groups/char_tabs/*.nb (Mathematica)
            ▼
 (3) loadable group file  (P C / ReChar / ImChar / mult / conclass)
            │                                       generators/genS108ct.py  (template)
            ▼
 (4) perfect-action character couplings β_c(β_bare)  generators/couplings_*.py
            │                                       ../chars/groups/char_tabs/final_result_gamma_in_beta*.nb
            ▼
 (5) Monte Carlo run                               ./dym-mod-metro ...
```

- **`generators/genS108ct.py`** — worked example of step (3): closes the group,
  assigns conjugacy classes, and prints a loadable file. Reproduces `groups/S108ct`
  exactly. Use it as the template for new groups (`python3 generators/genS108ct.py > groups/S108ct`).
- **`generators/couplings_*.py`** — step (4): print a `dym-mod-metro` command whose
  `β_c` are the **character expansion of an (approximately) perfect action**,
  expressed as polynomials in the bare coupling. `couplings_1t.py … couplings_5t.py`
  are truncations keeping 1…5 expansion terms; `couplings_3rd.py` and `couplings.py`
  are higher-order / full versions; `couplings_scale*.py` handle scale setting.
- **`groups/char_tabs/`** — the Σ(1080) character table (`char_table.py`, 17×17) plus
  the strong-coupling result dumps (`c5*`, `result_class_element.npy`, read by
  `read_results.py`). The heavy **Mathematica notebooks** that derive the
  perfect-action couplings to 5th order (`*-order.nb`, `final_result_gamma_in_beta*.nb`)
  are **not** kept in this repo to keep it light — they live in the sibling tree
  `../chars/groups/char_tabs/`. Nothing in the run path depends on them.

---

## Wilson-flow scale setting (w₀, t₀)

To set the lattice scale and to use **gradient-flow observables** as a clean,
renormalized probe of improvement (the methodology of the classically-perfect
gradient-flow program — Holland, Ipp, Müller & Wenger, arXiv:2504.15870), the
discrete-group configurations are embedded back into SU(3) and flowed:

```
multi-character MC ──► save SU(3) configs (NERSC) ──► Wilson flow ──► w₀ / t₀
   dym-mod-metro-savecfg          (each config)        tools/wilsonflow/    scripts/
```

1. **A group file with appended defining-rep matrices** is required so each link
   index can be written as its SU(3) matrix. `groups/S1080ctm` is exactly this:
   the Σ(1080) character table **and** the 1080 defining-rep 3×3 matrices in one
   element ordering. It is built by
   `python3 generators/build_s1080_matrix_file.py` (ordering + matrices taken from
   the `dym_symanzik_action` repo's `mys1080-v4`; the full character table is
   recomputed from the group's Cayley table by the **Burnside–Dixon** method, then
   verified: SU(3)-ness, the homomorphism `M(a)M(b)=M(ab)`, character orthogonality,
   and that row 1 = fundamental / row 3 = adjoint as in `S1080ct`).

2. **Dump configs** (note the extra `outprefix [K] [N] [Ntherm]` args; needs D=4):
   ```
   ./dym-mod-metro-savecfg groups/S1080ctm 4 Nt Nx  β_0 … β_{C-1}  seed  out/  K N Ntherm
   ```
   writes `out/nersc-bf<β1>-ba<β3>-…-num####` in NERSC `4D_SU3_GAUGE_3X3` format.
   `wilson_flow` re-validates the plaquette from the matrices on load, so a wrong
   embedding is caught immediately.

3. **Flow + extract** (`tools/wilsonflow/build.sh` builds `wilson_flow` + `w0_scale`,
   the BMW standalone tools).

The whole 1t→5t ladder is orchestrated by three scripts (see
[`scripts/README.md`](scripts/README.md) for laptop/lenore presets):
```
scripts/run_ladder.sh OUTDIR "1t 2t 3t 4t 5t" "BVLIST" D NT NX "SEEDS" K N NTHERM [NCONC NTHREADS NICE]
scripts/flow_all.sh   OUTDIR [eps] [tmax] [nconc]      # Wilson-flow every saved config
scripts/w0_table.py   OUTDIR [ref]                     # w0/a, sqrt(t0)/a, <plaq> per (trunc, bv)
```
`run_ladder.sh` pulls the correct remapped β-vector per truncation via
`couplings_for_ctm.py`, launches one `nice`'d `dym-mod-metro-savecfg` per point, and
writes `OUTDIR/<trunc>/b<bv>/seed<seed>/{nersc-*, run.log, flow.dat}`.
`scripts/string_tension.py` fits √σ from the Wilson-loop grid; along a coupling
trajectory the dimensionless ratio √σ·w₀ should stay constant where the discrete
theory mimics continuum SU(3) (watch for the freezing line pre-empting the window).

> On a tiny lattice (e.g. 4⁴) `w0_scale` will report "scale determination failed" —
> there is no scale window; use a thermalized ensemble on a larger lattice.

### Perfect-action couplings on `S1080ctm`

`couplings_*.py` emit the character-expansion (perfect-action) β-vector in
`S1080ct`'s irrep-row order, which differs from `S1080ctm`'s. To apply them to the
flow-enabled file:

```
python3 generators/couplings_for_ctm.py  TRUNC bv D Nt Nx [seed]    # TRUNC = 1t..5t
```

prints the β-vector **remapped to `S1080ctm`** and ready `dym-mod-metro` /
`dym-mod-metro-savecfg` commands. The permutation π is computed once and verified
to machine precision by `generators/map_s1080ct_to_ctm.py`
(→ `s1080ct_to_ctm.json`): it is the group isomorphism between the two element
labelings, **anchored on the defining representation** so each coupling lands on
the correct irrep — e.g. the 6-dim couplings go to `S1080ctm` rows 8,9, *not* its
rows 4,5 (which are 3-dims). Verified that the remapped couplings reproduce the
identical action (`|S_ct − S_ctm∘σ| ≈ 2·10⁻¹³`). The truncation ladder turns on
successively more characters in perfect-action order: `1t` = the defining 3⊕3̄
(≈ discrete Wilson) through `5t` = all 17 irreps.

---

## Building

```bash
make            # -> dym-mod-metro  and  dym-mod-metro-savecfg  (g++ -O3 -fopenmp)
( cd tools/wilsonflow && ./build.sh )   # -> wilson_flow, w0_scale
make clean
```

Requirements: a C/C++ toolchain with OpenMP; Python 3 + NumPy for the generators;
Mathematica only if you want to re-derive the perfect-action couplings.

---

## Notes & caveats

- The Metropolis proposal updates a link by multiplying with a randomly chosen
  *near-identity* group element ("small group" — the elements with the
  largest `Re Tr`), giving a local, ergodic, detailed-balance move.
- `groups/character-expansion/` is a vendored copy of an external utility (its
  upstream `.git` was removed); see its own `README.md`/`LICENSE`.
- Larger / additional groups (`S216`, `S648`, `S60`, binary `BI`/`BO`/`BT`,
  dihedral `Dn`, cyclic `Zn`, `A6`, Clifford groups) have matrix generators under
  `generators/group_matrices/` but need the character-table step before they load.

---

## Origin

Extracted and cleaned from the `chars_clean` working tree of the `dym` (discrete
Yang–Mills) project, reorganized into a self-contained repository focused on the
multi-character / perfect-action investigation.
