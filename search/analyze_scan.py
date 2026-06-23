#!/usr/bin/env python3
"""Group a run_points scan by coupling combo (averaging over seeds), rank by w0/a.
Strips a trailing _s<seed> from each point label to group seeds. Frozen points count
as failures for that seed.

usage: search/analyze_scan.py SCAN_OUTDIR
"""
import sys, os, glob, re, collections
import numpy as np

out = sys.argv[1]
g = collections.defaultdict(list)
for f in glob.glob(os.path.join(out, '*', 'result.txt')):
    txt = open(f).read().strip()
    lab = os.path.basename(os.path.dirname(f))
    combo = re.sub(r'_s\d+$', '', lab)
    m = re.search(r'w0/a=([0-9.]+)', txt)
    if 'FROZEN' in txt or not m:
        g[combo].append(np.nan)
    else:
        g[combo].append(float(m.group(1)))

rows = []
for combo, vals in g.items():
    v = [x for x in vals if np.isfinite(x)]
    if v:
        rows.append((combo, np.mean(v), (np.std(v, ddof=1)/np.sqrt(len(v)) if len(v) > 1 else 0.0), len(v), len(vals)))
    else:
        rows.append((combo, np.nan, np.nan, 0, len(vals)))
rows.sort(key=lambda r: -(r[1] if np.isfinite(r[1]) else -1))
print(f"# scan: {out}")
print(f"{'combo':24s} {'w0/a':>8} {'+-':>7} {'nconf':>6}")
for combo, m, e, nc, nt in rows:
    w = f"{m:.3f}" if np.isfinite(m) else "frozen"
    print(f"{combo:24s} {w:>8} {e:7.3f} {nc}/{nt}")
