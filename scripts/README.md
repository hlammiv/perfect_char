# Perfect-action ladder + Wilson-flow pipeline

Three steps, all run **from the repo root**. They generate multi-character
ensembles on Σ(1080) (`S1080ctm`), Wilson-flow the saved SU(3) configs, and
tabulate the scale `w0/a`.

```
1.  scripts/run_ladder.sh  OUTDIR "TRUNCS" "BVLIST" D NT NX "SEEDS" K N NTHERM [NCONC] [NTHREADS] [NICE]
2.  scripts/flow_all.sh    OUTDIR [EPS] [TMAX] [NCONC]
3.  scripts/w0_table.py    OUTDIR [ref]
```

Output tree: `OUTDIR/<trunc>/b<bv>/seed<seed>/{nersc-*, run.log, flow.dat}`.

## Prerequisites (once, on each machine)

```
make                                   # dym-mod-metro + dym-mod-metro-savecfg
( cd tools/wilsonflow && ./build.sh )  # wilson_flow + w0_scale
# groups/S1080ctm is committed; regenerate only if needed:
#   python3 generators/build_s1080_matrix_file.py   (needs ../improved/groups/mys1080-v4)
```

The β-vectors are the perfect-action character couplings, **remapped to
`S1080ctm`'s irrep order** automatically (via `generators/couplings_for_ctm.py`
+ `s1080ct_to_ctm.json`). The ladder `1t..5t` turns on successively more
characters (1t = defining 3⊕3̄ ≈ discrete Wilson … 5t = all 17).

## What to watch

- `w0/a` should grow smoothly along a coupling trajectory; where the discrete
  theory mimics continuum SU(3), dimensionless ratios (e.g. √σ·w0) stay constant.
- A first-order **freezing** line can pre-empt the scaling window — watch the
  plaquette/acceptance in `run.log`.
- Tiny lattices have no scale window (`w0_scale` → "failed"); use thermalized
  ensembles on a real volume.

## Presets

**Laptop pilot** (12³×4, full ladder, one coupling, polite on ~6 of the free
cores — won't disturb other jobs). ~2–3 h, then flow + table:
```
scripts/run_ladder.sh runs/pilot "1t 2t 3t 4t 5t" "8.4" 4 4 12 "1" 200 40 2000  3 2 19
scripts/flow_all.sh   runs/pilot 0.01 6.0 6
scripts/w0_table.py   runs/pilot
```

**lenore production** (12³×4 scan, several couplings × seeds; set NCONC×NTHREADS
to the free thread budget — 8×4=32 once its gh_string batch finishes, lower while
it's busy):
```
scripts/run_ladder.sh runs/scan "1t 2t 3t 4t 5t" "7.6 8.0 8.4 8.8" 4 4 12 "1 2 3" 300 60 3000  8 4 10
scripts/flow_all.sh   runs/scan 0.01 10.0 12
scripts/w0_table.py   runs/scan
```
For 16⁴ (needs lenore's 125 GB): use `D=4 NT=16 NX=16` (≈ 1.3 s/sweep, ~3 h per point).

On lenore: `git clone https://github.com/hlammiv/perfect_char.git` (S1080ctm is
committed), then the prerequisites above. No scheduler — these scripts just launch
`nice`'d background processes, like your existing runs.

## Cost (measured, this i7-12700H @ 14 threads)

| lattice | s/sweep 1t→5t | one point (10k sweeps) |
|---|---|---|
| 12³×4 | 0.10 → 0.14 | ~17–23 min |
| 16⁴   | 0.9 → 1.3   | ~2.6–3.6 h |

Wilson flow adds ~20–50% wall, embarrassingly parallel.
