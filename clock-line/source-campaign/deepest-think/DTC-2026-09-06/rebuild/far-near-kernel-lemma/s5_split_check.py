#!/usr/bin/env python3
"""S5 - KILL K4 and K5: evaluate the far / collar / inner pieces EXACTLY (term-by-term e-fold
integrals of the shell operator) and compare with the proved constants C1, (4/s+pi/8), C2in.
Also: the phi-profile of the inner-edge plateau constant c_edge(phi) and its behaviour as phi -> 0.
Bang-bang cap w = -M sgn(z), rho0 = 1, R = 4096, M = 1."""
import numpy as np, math, json, hashlib, sys, time
sys.path.insert(0, '.')
from kern import geg_table, Nl, bangbang_g, cesaro

t0=time.time(); LMAX = 1001
g = bangbang_g(LMAX, nquad=2500)
ls = np.arange(1, LMAX+1, 2)
def pieces(rho, phi, rho0=1.0, R=4096.0):
    t = math.cos(phi); C = geg_table(np.array([t]), LMAX+1)[:, 0]
    cin  = -((ls+2)*g[ls]/(2*ls+3))*C[ls-1]
    cout =  ((ls+1)*g[ls]/(2*ls+3))*C[ls+1]
    eps, U = rho/R, rho/rho0                       # u = rho/rho' ranges over (eps, U)
    # FAR  : u in (eps, 1/2)   [rho' > 2 rho]
    Ifar = np.where(ls == 1, math.log(0.5/eps), (0.5**np.maximum(ls-1,1) - eps**np.maximum(ls-1,1))/np.maximum(ls-1,1))
    far_x   = cesaro(cin*Ifar)[1]
    far_0   = float(cin[0]*Ifar[0])                # only l=1 survives at x = 0
    # COLLAR: u in (1/2, min(1,U)) interior + (1, min(2,U)) exterior
    hi_in  = min(1.0, U)
    Icol_in = np.where(ls == 1, math.log(hi_in/0.5),
                       (hi_in**np.maximum(ls-1,1) - 0.5**np.maximum(ls-1,1))/np.maximum(ls-1,1))
    Icol_in = np.where(hi_in > 0.5, Icol_in, 0.0)
    lo_out, hi_out = max(1.0, eps), min(2.0, U)
    Icol_out = np.where(hi_out > lo_out, (lo_out**(-(ls+4)) - hi_out**(-(ls+4)))/(ls+4), 0.0)
    collar = cesaro(cin*Icol_in)[1] + cesaro(cout*Icol_out)[1]
    # INNER : u in (2, U)   [rho' < rho/2]
    Iin = np.where(U > 2.0, (2.0**(-(ls+4)) - U**(-(ls+4)))/(ls+4), 0.0)
    inner = cesaro(cout*Iin)[1]
    return far_x, far_0, collar, inner

C1odd, C2inodd = 0.291999, 0.014754     # from s2_constants.json (z-odd column)
print("exact piecewise split of a(x) for the bang-bang cap (rho0=1, R=4096, M=1)")
print(f"{'rho':>7} {'phi(deg)':>8} {'far(x)':>10} {'far(0)':>10} {'|far-far0|':>11} {'C1':>7} {'collar':>9} {'bnd 4/s+pi/8':>13} {'inner':>10} {'C2in':>9} {'sum':>10} {'(M/2)lnR/2rho':>14}")
k4 = k5 = False; rows=[]
for rho in (8.0, 64.0, 512.0):
    for phid in (10.0, 45.0, 90.0):
        phi = math.radians(phid); s = math.sin(phi)
        fx, f0, col, inn = pieces(rho, phi)
        bnd = 4.0/s + math.pi/8
        tot = fx + col + inn
        rows.append(dict(rho=rho, phi_deg=phid, far=fx, far0=f0, collar=col, inner=inn, total=tot, bound_collar=bnd))
        print(f"{rho:7.1f} {phid:8.1f} {fx:10.5f} {f0:10.5f} {abs(fx-f0):11.6f} {C1odd:7.4f} {col:9.5f} {bnd:13.3f} {inn:10.2e} {C2inodd:9.5f} {tot:10.5f} {0.5*math.log(4096/(2*rho)):14.5f}")
        if abs(fx-f0) > C1odd: k4 = True
        if abs(col) > bnd or abs(inn) > C2inodd: k5 = True
print("KILL K4", "FIRED" if k4 else "did not fire", "   KILL K5", "FIRED" if k5 else "did not fire")

# K5': does the NEAR part (collar+inner) grow with log(rho/rho0)?  fit over 6 octaves at fixed phi
print("\nK5': near part (collar+inner) vs log(rho/rho0) at phi = 45 deg (must be flat)")
xs=[]; ys=[]
for j in range(1, 7):
    rho = 2.0**j; fx, f0, col, inn = pieces(rho, math.radians(45.0))
    xs.append(math.log(rho)); ys.append(col+inn); print(f"   rho={rho:7.1f}  near = {col+inn:+.6f}")
sl = float(np.polyfit(xs, ys, 1)[0]); print(f"   fitted d(near)/d log(rho/rho0) = {sl:+.5f}  (threshold 0.02)")
k5p = abs(sl) > 0.02; print("KILL K5'", "FIRED" if k5p else "did not fire")

# --- the phi-profile of the inner-edge plateau constant, and phi -> 0 -----------------------------
print("\ninner-edge plateau constant c_edge(phi) = a(rho0+,phi) - (M/2) log(R/rho0),  R -> inf")
def c_edge(phi):
    t = math.cos(phi); C = geg_table(np.array([t]), LMAX+1)[:, 0]
    l3 = np.arange(3, LMAX+1, 2)
    return cesaro(-((l3+2)*g[l3]/((2*l3+3)*(l3-1)))*C[l3-1])[1]
prof = {}
for phid in (1,2,3,5,7.5,10,15,20,30,45,60,75,90):
    prof[phid] = c_edge(math.radians(phid)); print(f"   phi={phid:5.1f} deg   c_edge = {prof[phid]:+.6f}")
print("   (c_edge grows without bound as phi -> 0: the bang-bang cap's inner-edge CORNER on the axis,")
print("    where eta = omega^theta/r is unbounded -- the same inadmissibility the sharp-clock refuter found.)")
json.dump(dict(rows=rows, near_slope=sl, c_edge_profile=prof), open('s5_results.json','w'), indent=1)
print(f"\nelapsed {time.time()-t0:.1f}s  SCRIPT-SHA256 {hashlib.sha256(open(__file__,'rb').read()).hexdigest()}")
