import numpy as np, glob, re, os
def split_configs(fn):
    """allflow.dat -> list of (t,E) arrays, one per config (t resets each config)."""
    rows=[]
    for l in open(fn):
        if l.startswith("WFLOW"):
            p=l.split(); rows.append((float(p[1]),float(p[2]),float(p[3])))
    rows=np.array(rows); cfgs=[]; start=0
    for i in range(1,len(rows)):
        if rows[i,0] < rows[i-1,0]:   # t reset -> new config
            cfgs.append(rows[start:i]); start=i
    cfgs.append(rows[start:])
    return cfgs
def cross(x,y,ref):
    for i in range(len(y)-1):
        if (y[i]-ref)*(y[i+1]-ref)<0:
            return x[i]+(ref-y[i])/(y[i+1]-y[i])*(x[i+1]-x[i])
    return np.nan
def scales(cfg):
    t=cfg[:,0]; E=(cfg[:,1]+cfg[:,2])/2.0; t2E=t*t*E
    W=t[1:-1]*(t2E[2:]-t2E[:-2])/(t[2:]-t[:-2])
    t0=cross(t,t2E,0.3); w0sq=cross(t[1:-1],W,0.3)
    return (np.sqrt(t0) if t0>0 else np.nan), (np.sqrt(w0sq) if w0sq>0 else np.nan)
def jack(vals):
    v=np.array([x for x in vals if np.isfinite(x)]); n=len(v)
    if n<2: return np.nan,np.nan
    m=v.mean(); jk=np.array([np.delete(v,i).mean() for i in range(n)])
    return m, np.sqrt((n-1)/n*np.sum((jk-m)**2))

print("# beta1 beta2  nconf   sqrt(t0)/a        w0/a")
res=[]
for d in sorted(glob.glob("/tmp/scan8/b1*")):
    f=f"{d}/allflow.dat"
    if not os.path.exists(f) or os.path.getsize(f)==0: continue
    b1=float(re.search(r"b1([\d.]+)_b2",d).group(1)); b2=float(re.search(r"_b2(-?[\d.]+)",d).group(1))
    cfgs=split_configs(f)
    st0=[scales(c)[0] for c in cfgs]; sw0=[scales(c)[1] for c in cfgs]
    mt,et=jack(st0); mw,ew=jack(sw0)
    res.append((b1,b2,mw,ew))
    print(f"  {b1:.0f}   {b2:+.3f}   {len(cfgs):3d}   {mt:.4f}+-{et:.4f}   {mw:.4f}+-{ew:.4f}")
print("\n(w0/a < 1 => flow scale below one lattice spacing: coarse/strong-coupling, no scale separation)")
