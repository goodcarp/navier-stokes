#!/usr/bin/env python3
"""x4 -- (a) independent re-derivation of the seat's exact identities (2.1)-(2.5),
            Lemma V.3 propagator, Lemma V.3', Prop 3.1, by my own sympy route;
        (b) apples-to-apples inset comparison: the seat's sec-5 f-table vs the same
            search run against Theorem V.4's OWN E_tail (+E_4), with the seat's own
            d = f rho0 convention and the seat's own target eps_bulk."""
import json, math
import numpy as np, sympy as sp
from scipy.stats import chi2
from scipy.special import roots_hermitenorm
R = {}

# ---------------------------- (a) identities, my own route -------------------------------
y1,y2,y3,y4,z = sp.symbols('y1 y2 y3 y4 z', real=True)
Y = sp.Matrix([y1,y2,y3,y4]); r = sp.sqrt((Y.T*Y)[0,0])
def lap5(f):
    return sum(sp.diff(f,v,2) for v in (y1,y2,y3,y4,z))
# (2.1) r^3 Delta_5 (1/r) = -1
e21 = sp.simplify(r**3*lap5(1/r)); R['A_r3_Lap5_inv_r'] = str(e21)
# cylindrical form Delta_5 = d_rr + (3/r) d_r + d_zz on f(r)
rr = sp.symbols('rho', positive=True)
e21b = sp.simplify(rr**3*(sp.diff(1/rr,rr,2) + 3/rr*sp.diff(1/rr,rr)))
R['A_r3_Lap5_inv_r_cyl'] = str(e21b)
# (2.2) Delta_5 eta_P = -eta_P/r^2 for eta_P = -M/r (z>0)
M = sp.symbols('M', positive=True); etaP = -M/r
R['A_2_2_residual'] = str(sp.simplify(lap5(etaP) + etaP/r**2))
# (2.3) Hessian of 1/r in R^5
H = sp.hessian(1/r, (y1,y2,y3,y4,z))
Hs = sp.simplify(H.subs({y1:sp.Symbol('R0',positive=True), y2:0,y3:0,y4:0}))
R['B_hess_eigs'] = str(sp.simplify(sp.Matrix(Hs).eigenvals()))
R['B_hess_trace_r3'] = str(sp.simplify(sp.trace(H)*r**3))
# (2.4) sup_|v|=1 |d_v^k (1/r)| = k!/r^{k+1}  via the Legendre generating function
s_, mu = sp.symbols('s mu'); Rv = sp.Symbol('Rv', positive=True)
gen = 1/sp.sqrt(Rv**2 + 2*s_*mu*Rv + s_**2)
res = {}
for k in range(1,7):
    coeff = sp.simplify(sp.series(gen, s_, 0, k+1).removeO().coeff(s_, k))
    res[k] = str(sp.simplify(coeff - sp.legendre(k, -mu)/Rv**(k+1)))
R['C_legendre_residuals'] = res
# (2.5) Delta_y^k (1/r) coefficients in R^4 (y only)
def lapy(f): return sum(sp.diff(f,v,2) for v in (y1,y2,y3,y4))
g = 1/r; cs=[]
for k in range(1,5):
    g = sp.simplify(lapy(g)); cs.append(str(sp.simplify(g*r**(2*k+1))))
R['D_Delta_y_k_coeffs'] = cs
# Lemma V.3 propagator: Psi solves Psi' = -S(tau-s) Psi with S = a diag(1,1,1,1,-2)
t_, u_, p_ = sp.symbols('t u p', positive=True); a = sp.Function('a')
tau = sp.Symbol('tau', positive=True)
mu_s = sp.Integral(a(t_), (t_, tau-sp.Symbol('s',positive=True), tau))
R['E_propagator_note'] = ("Psi_s = diag(e^{-mu(s)}x4, e^{2mu(s)}), mu(s)=int_{tau-s}^{tau} a; "
                          "Psi_tau Psi_p^{-1} = diag(lam(tau-p)^{-1} x4, lam(tau-p)^2), "
                          "lam(u)=exp int_0^u a  -- verified by direct substitution below")
S_ = sp.Symbol('s', positive=True)
lamf = sp.Function('lam')
# check d/ds [e^{-mu(s)}] = -a(tau-s) e^{-mu(s)}  (the y-block of Psi' = -S(tau-s)Psi)
mu_expr = sp.Integral(a(t_), (t_, tau-S_, tau))
lhs = sp.diff(sp.exp(-mu_expr), S_)
R['E_prop_residual_y'] = str(sp.simplify(sp.expand(lhs + a(tau-S_)*sp.exp(-mu_expr))))
lhs2 = sp.diff(sp.exp(2*mu_expr), S_)
R['E_prop_residual_z'] = str(sp.simplify(sp.expand(lhs2 - 2*a(tau-S_)*sp.exp(2*mu_expr))))
# Lemma V.3': (1/2) C : Hess(eta_P) with C = diag(2 sy I4, 2 sz)
sy, sz, R0 = sp.symbols('sigma_y sigma_z R0', positive=True)
Cmat = sp.diag(2*sy,2*sy,2*sy,2*sy,2*sz)
Hp = sp.hessian(etaP, (y1,y2,y3,y4,z))
contr = sp.Rational(1,2)*sum(Cmat[i,j]*Hp[i,j] for i in range(5) for j in range(5))
contr0 = sp.simplify(contr.subs({y1:R0,y2:0,y3:0,y4:0}))
R['E_half_C_contract'] = str(contr0)
R['E_identification_residual'] = str(sp.simplify(contr0 - (M*sy/R0**3)))
R['E_dsigma_z'] = str(sp.simplify(sp.diff(contr0, sz)))
# Prop 3.1: V^k = D^{-1} d_l(D G^{lk}) vanishes for affine X  (my own n=3 check)
al = sp.symbols('alpha1:4'); A3 = sp.Matrix(3,3, lambda i,j: sp.Rational((i+1)*(j+2), (i+j+3)))
X3 = A3*sp.Matrix(al) + sp.Matrix([1,2,3])
J3 = X3.jacobian(sp.Matrix(al)); D3 = sp.simplify(J3.det()); G3 = sp.simplify((J3.T*J3).inv())
V3 = [sp.simplify(sum(sp.diff(D3*G3[l,k], al[l]) for l in range(3))/D3) for k in range(3)]
R['F_V_affine_n3'] = [str(v) for v in V3]
# curved map in n=3 -> V != 0 at first order
eps = sp.Symbol('epsilon')
Xc = sp.Matrix([al[0]+eps*al[1]**2, al[1]+eps*al[0]*al[1], al[2]])
Jc = Xc.jacobian(sp.Matrix(al)); Dc = sp.simplify(Jc.det()); Gc = sp.simplify((Jc.T*Jc).inv())
Vc = [sp.simplify(sum(sp.diff(Dc*Gc[l,k], al[l]) for l in range(3))/Dc) for k in range(3)]
R['F_V_curved_n3_linear_in_eps'] = [str(sp.simplify(sp.series(v, eps, 0, 2).removeO())) for v in Vc]
print("(a) identities, my own route:")
for k in ['A_r3_Lap5_inv_r','A_r3_Lap5_inv_r_cyl','A_2_2_residual','B_hess_eigs','B_hess_trace_r3',
          'D_Delta_y_k_coeffs','E_prop_residual_y','E_prop_residual_z','E_half_C_contract',
          'E_identification_residual','E_dsigma_z','F_V_affine_n3','F_V_curved_n3_linear_in_eps']:
    print(f"   {k:34s} {R[k]}")
print(f"   C_legendre_residuals               {R['C_legendre_residuals']}")

# ---------------------------- (b) inset comparison ---------------------------------------
TH = 4*(1-math.sqrt(2/3)); I2 = 4/5 - 16*math.sqrt(6)/135; I4 = -4/7 + 27*math.sqrt(6)/28
Cw = math.sqrt(1.5)*TH; n = 5; sd = math.sin(math.radians(7.5))
_hx,_hw = roots_hermitenorm(400); _hw=_hw/_hw.sum()
def P_gt(a, sy_, sz_):
    if a<=0: return 1.0
    g=math.sqrt(2*sz_)*_hx; rem=a*a-g*g
    return float((_hw*np.where(rem<=0,1.0,chi2.sf(np.maximum(rem,0.0)/(2*sy_),df=4))).sum())
def bound(L,f,phi0=30.0,K2hat=0.0,Rm_is_r0=False):
    nu=sd**2; tau=TH/L; r0=(1+f)*math.sin(math.radians(phi0))
    sy=nu*I2/L; sz=nu*I4/L; sC=sy/r0**2; N=r0/sd
    d=min(f, 0.999*r0)
    Rm = r0 if Rm_is_r0 else max(r0-d,1e-9)
    K2=K2hat; V1=n*nu*tau*(math.exp(2*Cw)-1)/Cw; EW=0.5*K2*math.exp(Cw)*tau*V1
    trC=8*sy+2*sz; trC2=4*(2*sy)**2+(2*sz)**2; m2=trC; m4=trC**2+2*trC2
    g=math.sqrt(2*sz)*_hx; a2=2*sy
    m6=float((_hw*(a2**3*192+3*a2**2*24*g*g+3*a2*4*g**4+g**6)).sum())
    P=2*P_gt(d/2,sy,sz)+P_gt(d,sy,sz)
    Eh=(r0/Rm**2)*EW + (4*N*EW/d if d>0 else float('inf'))
    E4=(r0/Rm**5)*m4
    Et=2*N*P+P+sum(r0**(-k)*math.sqrt([m2,m4,m6][k-1])*math.sqrt(P) for k in (1,2,3))
    return dict(s_C=sC,E_hess=Eh,E_4=E4,E_tail=Et,total=Eh+E4+Et)
def scan(L, target, K2hat=0.0):
    """scan f and phi0; return the smallest f (over the best phi0) meeting the target,
       respecting (H3): d = f rho0 must satisfy d < r0 = (1+f) sin phi0."""
    best=None
    for phi0 in [30.0,45.0,60.0,75.0]:
        sp0=math.sin(math.radians(phi0))
        fmax = 0.98*sp0/(1-sp0) if sp0<1 else 50.0
        for f in np.linspace(1e-4, min(fmax,50.0), 4000):
            b=bound(L,f,phi0=phi0,K2hat=K2hat)
            tg=(1.0/L) if target=='1/L' else b['s_C']
            if b['E_tail']+b['E_4']+b['E_hess'] <= tg:
                if best is None or f<best[0]: best=(f,phi0,b)
                break
    return best
print("\n(b) inset f, d = f rho0 (the seat's own convention), K2hat = 0 (tail+E_4 only);")
print("    (H3) forces d < r0, i.e. f < sin(phi0)/(1-sin(phi0)); phi0 scanned over 30/45/60/75 deg")
print(" L    seat f (vs eps_bulk)   theorem f (vs eps_bulk)   theorem f (vs 1/L)   gap-V f (vs 1/L)")
seat_f={10:0.25972,20:0.19228,40:0.14189,80:0.10441,160:0.07663}
gapV_f={10:0.9131,20:0.6590,40:0.4694,80:0.3382,160:0.2437}
rows=[]
for L in [10,20,40,80,160,320,640]:
    be=scan(L,'eps'); b1=scan(L,'1/L')
    rows.append(dict(L=L,f_seat=seat_f.get(L),
                     f_theorem_eps=(be[0] if be else None), phi_eps=(be[1] if be else None),
                     f_theorem_1L=(b1[0] if b1 else None), phi_1L=(b1[1] if b1 else None),
                     f_gapV=gapV_f.get(L)))
    fe = f"{be[0]:.4f}@{be[1]:.0f}deg" if be else "NOT ACHIEVABLE"
    f1 = f"{b1[0]:.4f}@{b1[1]:.0f}deg" if b1 else "NOT ACHIEVABLE"
    print(f"{L:5d}   {str(seat_f.get(L,'-')):>10s}           {fe:>16s}      {f1:>16s}      {str(gapV_f.get(L,'-')):>8s}")
R['inset_comparison']=rows
json.dump(R, open('x4_results.json','w'), indent=1, default=str)
print("\nWROTE x4_results.json")
