#!/usr/bin/env python3
"""G15 — gate G-b of GENERATION_1: exact sign structure of blob-dipole and compressor-pair strains (numerics, own 5D-lift Biot–Savart).
Elements: Gaussian ω^θ blobs W exp(-((r-r_i)^2+(z-z_i)^2)/σ^2); η = ω^θ/r; a(r,z) = (3/(2π))∫∫ η r'^3 (z-z') J dr'dz',
J = ∫_0^π sin²θ dθ /(r²+r'²-2rr'cosθ+(z-z')²)^{5/2}  (same kernel as G2, calibrated there to 1e-5).
Blob dipole (what the swirl source ∂_z(q²) creates around a swirl maximum at z_c=0): ring at (r_G, -Z) with ω^θ = +W, ring at (r_G, +Z) with ω^θ = −W.
Claims: (1) a(0,0) > 0 (expansion at the blob center); (2) at the midplane outside the dipole, a(ρ,0) > 0 for ρ ∈ {2,4,8} r_G (the blob stretches an outer compressor);
(3) a 'separating' compressor pair (ω^θ = −W at (ρ_c,−Z_c), +W at (ρ_c,+Z_c)) gives a(0,0) < 0 (inflow at the blob), the 'colliding' orientation gives > 0.
Control: reversing the dipole sign must reverse every sign; the single-ring axis value must match the analytic thin-ring limit a(0,0) = (3/4)·Γ_ring·r_G²·(−z_ring)/(r_G²+z_ring²)^{5/2}·(normalization) — checked by comparing the Gaussian ring at small σ with the δ-ring formula."""
import numpy as np, sys, hashlib
def gl(n,a,b):
    x,w=np.polynomial.legendre.leggauss(n); return 0.5*(b-a)*x+0.5*(b+a), 0.5*(b-a)*w
TH=[gl(48,0,np.pi/64),gl(48,np.pi/64,np.pi/8),gl(64,np.pi/8,np.pi)]
th=np.concatenate([p[0] for p in TH]); s2=np.sin(th)**2*np.concatenate([p[1] for p in TH])
def strain_from_blob(r,z,rc,zc,W,sig,n=96):
    rp,wr=gl(n,max(rc-4*sig,1e-9),rc+4*sig); zp,wz=gl(n,zc-4*sig,zc+4*sig)
    RP,ZP=np.meshgrid(rp,zp,indexing='ij'); W2=np.outer(wr,wz)
    om=W*np.exp(-((RP-rc)**2+(ZP-zc)**2)/sig**2); eta=om/RP; dz=z-ZP
    D0=r**2+RP**2+dz**2; c=2*r*RP
    J=np.zeros_like(RP)
    for k in range(0,th.size,32):
        t=th[k:k+32]; ws=s2[k:k+32]; D=D0[...,None]-c[...,None]*np.cos(t); J+=np.einsum('ijk,k->ij',D**(-2.5),ws)
    return (3/(2*np.pi))*np.sum(eta*RP**3*dz*J*W2)
def a_total(r,z,elems):
    return sum(strain_from_blob(r,z,*e) for e in elems)
rG,Z,W,sig=1.0,1.0,1.0,0.15
# control 1: thin ring vs delta-ring formula at the axis: Gamma = ∫ω dr dz = W π σ² ; a(0,0) = (3/4) ∫∫ η r'^3 (−z')/(r'^2+z'^2)^{5/2} dr'dz' ≈ (3/4) Γ r_G^2 (−z_ring)/(r_G^2+z_ring^2)^{5/2}
Gam=W*np.pi*sig**2
a_num=strain_from_blob(0,0,rG,-Z,W,sig); a_delta=(3/4)*Gam*rG**2*(Z)/(rG**2+Z**2)**2.5
print(f'CTRL thin-ring axis value: numeric {a_num:.6f} vs delta-ring {a_delta:.6f} (rel {abs(a_num-a_delta)/a_delta:.3f}; σ/r_G={sig})')
dipole=[(rG,-Z,+W,sig),(rG,+Z,-W,sig)]
a0=a_total(0,0,dipole); print(f'(1) blob dipole at its center a(0,0) = {a0:+.5f}  (claim > 0: expansion)')
signs=[]
for rho in [2,4,8]:
    v=a_total(rho*rG,0,dipole); signs.append(v>0); print(f'(2) blob dipole at midplane rho={rho} r_G: a = {v:+.6f}')
Zc=2.0; rc=4.0
sep=[(rc,-Zc,-W,sig),(rc,+Zc,+W,sig)]; col=[(rc,-Zc,+W,sig),(rc,+Zc,-W,sig)]
a_sep=a_total(0,0,sep); a_col=a_total(0,0,col)
print(f'(3) separating compressor pair at rho_c=4, Z_c=2: a(0,0) = {a_sep:+.6f} (claim < 0: inflow); colliding pair: {a_col:+.6f} (claim > 0)')
rev=[(rG,-Z,-W,sig),(rG,+Z,+W,sig)]; a0r=a_total(0,0,rev)
print(f'CTRL reversed dipole a(0,0) = {a0r:+.5f} (must be −{a0:.5f})')
ok = abs(a_num-a_delta)/a_delta < 0.05 and a0>0 and all(signs) and a_sep<0 and a_col>0 and abs(a0r+a0)<1e-9*max(1,abs(a0))
print('RESULT','PASS' if ok else 'FAIL','| SCRIPT-SHA256',hashlib.sha256(open(__file__,'rb').read()).hexdigest())
sys.exit(0 if ok else 1)
