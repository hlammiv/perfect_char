# Perfect-action search (single scale, fixed active set, by reweighting)

Search the multi-character couplings `{β_c}` that minimize gradient-flow cutoff
effects — i.e. find the *right trajectory* directly, NOT from the strong-coupling
`couplings_Nt` expansion.

## Why this is cheap

The action is **linear** in the couplings, `S(U) = -Σ_c β_c P_c(U)` with
`P_c(U)=Σ_plaq Reχ_c(U_p)`. So one ensemble at base `β⁰` scores *any* nearby `β`
by **exact reweighting** of its per-config `P_c` (logged in `pchar.dat`) and its
per-config flow curves — no new Monte Carlo. The search varies only the
**user-specified active characters**; the rest stay fixed.

## Workflow

```
1. base ensemble in the SCALING WINDOW (weak side, single phase below freezing):
     ./dym-mod-metro-savecfg ./groups/S1080ctm 4 Nt Nx  <beta vector>  seed  ENS/  K N Ntherm
   -> ENS/nersc-*  +  ENS/pchar.dat   (P_c, simpleplaq, polyakov per config)

2. flow every config (per-config, keyed by num):
     search/flow_configs.sh ENS [eps] [tmax] [nconc]      -> ENS/flow_<num>.dat

3. search the active set:
     search/search_scale.py ENS --active "1 2 3" [--ref FILE] [--essmin 0.1]
                                [--window 0.2 1.5] [--maxstep 0.05]
   -> best β (only the active chars moved), cost, ESS, w0/a; flags if it hit the
      overlap/step boundary (then re-anchor: regenerate a base ensemble at the new β).
```

## Cost (what defines "right")

Gradient-flow curve, single scale:
- `--ref FILE` (cols `t/w0^2  t^2E`): squared deviation of the lattice
  `t²⟨E⟩(t/w0²)` from your continuum SU(3) reference curve (the Holland target).
  **Supply your own continuum-extrapolated curve.**
- default (no `--ref`): self-contained **anisotropy** proxy `∫(t²E1−t²E2)²`
  (electric `gact4i` vs magnetic `gactij`); zero in the isotropic continuum.

## Critical caveats (learned the hard way)

- The base ensemble **must be in the scaling window**: `w0/a ~ O(1)`. Over-ordered
  / nearly frozen couplings give a trivial smooth field with no resolved scale
  (`t²E` never reaches 0.3 → `w0=nan` → cost `inf`). Too fine for the box does the
  same (scale doesn't fit). The search reports this rather than guessing.
- Reweighting (and any local search) **cannot cross the first-order freezing
  transition** — stay inside one phase. The `--essmin` floor enforces overlap;
  re-anchor with a fresh ensemble to move along a path.

## Validation (engine correctness)

- `validate_reweight_internal.py ENS [chars...]` — transition-immune: checks the
  reweighting against the analytic linear response `d⟨O⟩/dβ=-cov(O,P)` from the
  same ensemble (confirms sign + math; residual scales like δ²). **Passes.**
- `validate_reweight.py BASE DIRECT` — checks reweighted ≈ a direct run; only valid
  inside one phase with good overlap (fails across transitions — by design).

## Files

| file | role |
|---|---|
| `reweight.py` | load `pchar.dat` + `flow_<num>.dat`; reweighting weights, ESS, means |
| `cost_flow.py` | reweighted `t²⟨E⟩`, `w0`, flow-curve cost (ref or anisotropy) |
| `search_scale.py` | fixed-active-set optimizer (Nelder–Mead, ESS floor) |
| `flow_configs.sh` | per-config Wilson flow → `flow_<num>.dat` |
| `validate_reweight*.py` | engine validation |
