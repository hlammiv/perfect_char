# Scaling pipeline: w0 (Wilson flow) + string tension

Both physical scales along a trajectory parallel to the freezing line, from one
ensemble. All dependency-free (no gwu-qcd, no MPI; the broken precompiled
`kentucky2nersc` is bypassed -- `dym-mod-metro-savecfg` writes NERSC directly).

## Tools (build: `./build.sh`)
- `wilson_flow` / `w0_scale` -- Budapest-Marseille-Wuppertal standalone code.
  Canonical source: `~/Dropbox/discrete groups/thermo-tuning/wilsonflow/`.

## Pipeline
1. **Generate** ensembles (configs + Wilson loops):
   ```
   scripts/deploy_scaling.sh OUTDIR NC NT NX DELTA "B1LIST" K N NTHERM
   # 8^4, one trajectory at Delta=1.0, beta1=2..6, ~30 configs/point:
   scripts/deploy_scaling.sh /tmp/scan8 6 8 8 1.0 "2 3 4 5 6" 300 30 3000
   ```
   Trajectory: beta2 = (0.308 - DELTA) - 0.296*beta1  (parallel to the MC freezing line).
   Writes `OUTDIR/b1*/cfg/nersc-*` and `OUTDIR/b1*/loops.dat`.

2. **w0** (flow every config, jackknife):
   ```
   scripts/flow_w0.sh /tmp/scan8 0.01 10.0 6     # eps tmax NC
   ```
3. **string tension** (static potential -> Cornell fit -> sqrt(sigma)):
   ```
   scripts/string_tension.py /tmp/scan8/b1X_b2Y/loops.dat NT NX
   ```

## Scaling test
Along the trajectory, w0/a should run smoothly and dimensionless ratios
(e.g. sqrt(sigma)*w0) stay constant where the theory mimics continuum SU(3).
A strongly first-order freezing line can pre-empt the window -- watch for it.

## Cost / placement
8^4 ~ 0.33 s/sweep single-thread (~16x a 4^4 sweep). Run single-threaded,
one process per point (plain Metropolis -- cache-friendly). Needs a genuinely
free machine; do NOT oversubscribe.
