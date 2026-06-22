#!/usr/bin/env python3
"""
String tension from the Wilson loops saved by dym-mod-metro-savecfg (GMES lines).

  W(r,t) = -wloop/3   (normalized 1/3 ReTr; our ReTr is negated, so actual = -wloop)
  a V(r) = ln[ W(r,t) / W(r,t+1) ]   (plateau in t)
  a V(r) = V0 + (a^2 sigma) r - (pi/12)/r     (Cornell)  -> sqrt(sigma) in 1/a

Jackknife errors over configurations.
Usage: string_tension.py loops.dat NT NX [--tmin T] [--tmax T] [--rmin R] [--rmax R]
"""
import sys, argparse, numpy as np
ap=argparse.ArgumentParser(); ap.add_argument("file"); ap.add_argument("NT",type=int); ap.add_argument("NX",type=int)
ap.add_argument("--tmin",type=int,default=1); ap.add_argument("--tmax",type=int,default=0)
ap.add_argument("--rmin",type=int,default=1); ap.add_argument("--rmax",type=int,default=0)
a=ap.parse_args(); NT,NX=a.NT,a.NX

rows=[]
for l in open(a.file):
    if not l.startswith("GMES:"): continue
    t=l.split()
    w=np.array([float(x) for x in t[5:5+NX*NT]])   # after 999,rep,imp,plaq
    if len(w)==NX*NT: rows.append(w.reshape(NX,NT))  # W[r-1][t-1]
G=np.array(rows)  # (ncfg, NX, NT)
if len(G)<2: sys.exit(f"need >=2 configs, got {len(G)}")
W = -G/3.0       # normalized Wilson loops
nc=len(W)
tmax=a.tmax or NT-1; rmax=a.rmax or NX-1

def sigma_from(Wmean):
    """effective potential + Cornell fit -> (sqrt_sigma, V0, [r,aV])."""
    r_list=[]; V_list=[]
    for r in range(a.rmin,rmax+1):
        # plateau of ln(W(r,t)/W(r,t+1)) over t in [tmin,tmax]
        vs=[]
        for t in range(a.tmin,tmax):
            num,den=Wmean[r-1,t-1],Wmean[r-1,t]
            if num>0 and den>0: vs.append(np.log(num/den))
        if vs: r_list.append(r); V_list.append(np.median(vs))
    r=np.array(r_list,float); V=np.array(V_list)
    if len(r)<3: return np.nan,np.nan,(r,V)
    # Cornell V = V0 + s*r + c/r ; fit linear in [1, r, 1/r] (c ~ -pi/12 if Coulombic)
    A=np.vstack([np.ones_like(r), r, 1.0/r]).T
    coef,*_=np.linalg.lstsq(A,V,rcond=None)
    s=coef[1]
    return (np.sqrt(s) if s>0 else np.nan), coef[0], (r,V)

ss,V0,(r,V)=sigma_from(W.mean(0))
# jackknife
jk=[]
for i in range(nc):
    Wj=np.delete(W,i,0).mean(0)
    sj,_,_=sigma_from(Wj)
    if np.isfinite(sj): jk.append(sj)
jk=np.array(jk)
err=np.sqrt((nc-1)/nc*np.sum((jk-jk.mean())**2)) if len(jk)>1 else np.nan
print(f"# configs={nc}  NT={NT} NX={NX}  t in [{a.tmin},{tmax}]  r in [{a.rmin},{rmax}]")
print("# r   aV(r)")
for rr,vv in zip(r,V): print(f"  {rr:.0f}  {vv:.5f}")
print(f"\nsqrt(sigma)*a = {ss:.5f} +/- {err:.5f}   V0 = {V0:.4f}")
if np.isfinite(ss): print(f"(1/sqrt(sigma) = {1/ss:.3f} a  -- want >> 1 and a stable ratio along the trajectory)")
