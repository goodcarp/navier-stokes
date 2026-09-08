"""
t3_maps_check.py -- numerical test of LEMMA T on two explicit map pairs.

a[eta0 o Phi^{-1}](0) = INT_S eta0(x) Ktil(Phi(x)) J_Phi(x) dx5 ,
Ktil(y) = (3/(8 pi^2)) (-y_z) |y|^{-5} ,  dx5 = 2 pi^2 rho^4 sin^3(phi) drho dphi.

Quadrature written from scratch (tensor Gauss-Legendre in u = log(rho/rho0) and phi,
phi split at pi/2 because eta0 has a sign jump there; also split at the taper angle).

Datum: eta0 = -M sgn(z) h_delta(phi) / r ,  h_delta = min(1, phi_ax/delta)  (delta=0 -> bang-bang)
       so |omega^theta| = M h_delta <= M, supported in S = {rho0<|x|<R}.

MAP PAIR A (uniform axial displacement; |Phi-T_lam| = mu*|x| EXACTLY):
   Phi = T_lambda + mu*rho*(0,-1)   [(r,z) components];  J = lambda^2 (1 - mu lambda^2 cos phi)
MAP PAIR B (shell-dependent strain -- the physically relevant one for GAP T):
   Phi = T_{Lambda(rho)},  Lambda(rho) = lambda (1 + eps sin(pi log(rho/rho0)/L))
MAP PAIR C (angular displacement, unit field (sin 2phi, cos 2phi)) -- kept as a third case;
   its first-order response happens to vanish, so it probes the quadratic regime.
"""
import json, numpy as np, sympy as sp, mpmath as mp
mp.mp.dps = 40

M, RHO0, LL = 1.0, 1.0, 6.0
R = RHO0*np.exp(LL)

# ---------------------------------------------------------------- symbolic maps -> callables
r_,z_,lam_,mu_,eps_ = sp.symbols('r z lam mu eps', positive=True)
rho_ = sp.sqrt(r_**2+z_**2)

def build(F,G):
    Fr,Fz,Gr,Gz = sp.diff(F,r_),sp.diff(F,z_),sp.diff(G,r_),sp.diff(G,z_)
    J = (F/r_)**3*(Fr*Gz-Fz*Gr)
    args = (r_,z_,lam_,mu_)
    return (sp.lambdify(args,F,'numpy'), sp.lambdify(args,G,'numpy'), sp.lambdify(args,J,'numpy'))

# T_lambda
FT, GT = lam_*r_, z_/lam_**2
T_F,T_G,T_J = build(FT,GT)
# pair A : uniform axial displacement
FA = lam_*r_
GA = z_/lam_**2 - mu_*rho_
A_F,A_G,A_J = build(FA,GA)
# pair C : angular unit field (sin 2phi, cos 2phi)
FC = lam_*r_ + mu_*rho_*(2*r_*z_/rho_**2)
GC = z_/lam_**2 + mu_*rho_*((z_**2-r_**2)/rho_**2)
C_F,C_G,C_J = build(FC,GC)
# pair B  (mu_ plays the role of eps here)
u_  = sp.log(rho_/RHO0)/LL
LamB = lam_*(1+mu_*sp.sin(sp.pi*u_))
FB, GB = LamB*r_, z_/LamB**2
B_F,B_G,B_J = build(FB,GB)

# ---------------------------------------------------------------- quadrature
def nodes(a,b,n):
    x,w = np.polynomial.legendre.leggauss(n)
    return 0.5*(b-a)*x+0.5*(a+b), 0.5*(b-a)*w

def a_functional(Ff,Gf,Jf,lam,mu,delta=0.0,nu=120,np_=120):
    """a[eta0 o Phi^{-1}](0) by my own tensor Gauss-Legendre."""
    U,WU = nodes(0.0,LL,nu)
    brk = [0.0]
    if delta>0: brk += [delta, np.pi-delta]
    brk += [np.pi/2]; brk = sorted(set(brk+[np.pi]))
    P=[];WP=[]
    for a_,b_ in zip(brk[:-1],brk[1:]):
        p,w = nodes(a_,b_,np_); P.append(p); WP.append(w)
    P=np.concatenate(P); WP=np.concatenate(WP)
    UU,PP = np.meshgrid(U,P,indexing='ij'); WW = np.outer(WU,WP)
    rho = RHO0*np.exp(UU); s,c = np.sin(PP), np.cos(PP)
    r,z = rho*s, rho*c
    h = np.ones_like(PP) if delta<=0 else np.minimum(1.0,np.minimum(PP,np.pi-PP)/delta)
    eta0 = -M*np.sign(c)*h/(rho*s)
    Fv,Gv,Jv = Ff(r,z,lam,mu),Gf(r,z,lam,mu),Jf(r,z,lam,mu)
    Fv = np.broadcast_to(np.asarray(Fv,dtype=float),PP.shape)
    Gv = np.broadcast_to(np.asarray(Gv,dtype=float),PP.shape)
    Jv = np.broadcast_to(np.asarray(Jv,dtype=float),PP.shape)
    absPhi = np.sqrt(Fv**2+Gv**2)
    integ = 0.75*(-Gv)*absPhi**-5*eta0*Jv*rho**5*s**3
    return float(np.sum(WW*integ))

def map_diagnostics(Ff,Gf,Jf,lam,mu,n=400):
    """measured mu (sup |Phi-T_lam|/|x|) and muJ (sup |J/lam^2 - 1|), on a fine grid."""
    U = np.linspace(0,LL,n); P = np.linspace(1e-6,np.pi-1e-6,n)
    UU,PP = np.meshgrid(U,P,indexing='ij')
    rho = RHO0*np.exp(UU); r,z = rho*np.sin(PP), rho*np.cos(PP)
    Fv,Gv,Jv = Ff(r,z,lam,mu),Gf(r,z,lam,mu),Jf(r,z,lam,mu)
    Fv=np.broadcast_to(np.asarray(Fv,float),PP.shape); Gv=np.broadcast_to(np.asarray(Gv,float),PP.shape)
    Jv=np.broadcast_to(np.asarray(Jv,float),PP.shape)
    d = np.sqrt((Fv-lam*r)**2+(Gv-z/lam**2)**2)/rho
    return float(d.max()), float(np.max(np.abs(Jv/lam**2-1))), float(Jv.min()), float((np.sqrt(Fv**2+Gv**2)/rho).min())

def I1(lam,mu):
    if mu >= min(lam,lam**-2): return float('inf')
    f = lambda p: mp.sin(p)**2*(mp.sqrt(lam**2*mp.sin(p)**2+lam**-4*mp.cos(p)**2)-mu)**-5
    return float(mp.quad(f,[0,mp.pi/2,mp.pi]))

def rhs(lam,mu,muJ):
    return M*LL*(3*(1+muJ)*lam**2*I1(lam,mu)*mu + lam/2*muJ)

out={}
# ---- calibration: Phi = T_lambda must reproduce Lemma 1 exactly
cal=[]
for lam in (1.0,1.1,1.25,1.5):
    for delta in (0.0, np.deg2rad(7.5)):
        v = a_functional(T_F,T_G,T_J,lam,0.0,delta)
        cal.append(dict(lam=lam,delta_deg=np.rad2deg(delta),a=v,lemma1=(M/2)*lam*LL if delta==0 else None,
                        rel_err=(abs(v-(M/2)*lam*LL)/((M/2)*lam*LL) if delta==0 else None)))
out['calibration_T_lambda_vs_Lemma1']=cal

# ---- convergence of the quadrature
conv=[]
for n in (60,120,240):
    conv.append(dict(n=n, a_T=a_functional(T_F,T_G,T_J,1.25,0.0,0.0,n,n),
                     a_A=a_functional(A_F,A_G,A_J,1.25,0.05,0.0,n,n)))
out['quadrature_convergence']=conv

# ---- PAIR A : three (lambda,mu)
rows=[]
for lam,mu in ((1.0,0.10),(1.25,0.05),(1.5,0.02)):
    for delta in (0.0, np.deg2rad(7.5)):
        aT = a_functional(T_F,T_G,T_J,lam,0.0,delta)
        aP = a_functional(A_F,A_G,A_J,lam,mu,delta)
        mm,mJ,Jmin,qmin = map_diagnostics(A_F,A_G,A_J,lam,mu)
        L_,R_ = abs(aP-aT), rhs(lam,mm,mJ)
        rows.append(dict(pair='A',lam=lam,mu_target=mu,mu_meas=mm,muJ_meas=mJ,Jmin=Jmin,
                         min_absPhi_over_rho=qmin,delta_deg=round(np.rad2deg(delta),2),
                         a_T=aT,a_Phi=aP,LHS=L_,RHS=R_,ratio_LHS_over_RHS=L_/R_,
                         LHS_over_mu_M_L=L_/(mm*M*LL),HOLDS=bool(L_<=R_)))
out['pairA']=rows

# ---- PAIR B : three (lambda,mu); solve eps so that measured mu hits the target
def solve_eps(lam,target):
    lo,hi=1e-8,0.9
    for _ in range(60):
        mid=0.5*(lo+hi)
        m,_,_,_ = map_diagnostics(B_F,B_G,B_J,lam,mid,200)
        if m<target: lo=mid
        else: hi=mid
    return 0.5*(lo+hi)
rowsB=[]
for lam,mu in ((1.0,0.10),(1.25,0.05),(1.5,0.02)):
    eps = solve_eps(lam,mu)
    for delta in (0.0, np.deg2rad(7.5)):
        aT = a_functional(T_F,T_G,T_J,lam,0.0,delta)
        aP = a_functional(B_F,B_G,B_J,lam,eps,delta)
        mm,mJ,Jmin,qmin = map_diagnostics(B_F,B_G,B_J,lam,eps)
        L_,R_ = abs(aP-aT), rhs(lam,mm,mJ)
        rowsB.append(dict(pair='B',lam=lam,eps=eps,mu_target=mu,mu_meas=mm,muJ_meas=mJ,Jmin=Jmin,
                          min_absPhi_over_rho=qmin,delta_deg=round(np.rad2deg(delta),2),
                          a_T=aT,a_Phi=aP,LHS=L_,RHS=R_,ratio_LHS_over_RHS=L_/R_,
                          LHS_over_mu_M_L=L_/(mm*M*LL),HOLDS=bool(L_<=R_)))
out['pairB']=rowsB

# ---- small-mu slope: is the true Lipschitz constant finite and << the proved one?
slope=[]
for lam in (1.0,1.25,1.5):
    for mu in (0.02,0.01,0.005,0.0025):
        aT=a_functional(T_F,T_G,T_J,lam,0.0,0.0); aP=a_functional(A_F,A_G,A_J,lam,mu,0.0)
        mm,mJ,_,_ = map_diagnostics(A_F,A_G,A_J,lam,mu,200)
        if True: pass
        slope.append(dict(lam=lam,mu=mu,ratio=abs(aP-aT)/(mm*M*LL),proved_C=3*(1+mJ)*lam**2*I1(lam,mm)+lam/2*(mJ/mm)))
out['pairA_small_mu_slope']=slope

# ---- PAIR C (quadratic-response case), same three (lambda,mu)
rowsC=[]
for lam,mu in ((1.0,0.10),(1.25,0.05),(1.5,0.02)):
    for delta in (0.0, np.deg2rad(7.5)):
        aT = a_functional(T_F,T_G,T_J,lam,0.0,delta)
        aP = a_functional(C_F,C_G,C_J,lam,mu,delta)
        mm,mJ,Jmin,qmin = map_diagnostics(C_F,C_G,C_J,lam,mu)
        L_,R_ = abs(aP-aT), rhs(lam,mm,mJ)
        rowsC.append(dict(pair='C',lam=lam,mu_target=mu,mu_meas=mm,muJ_meas=mJ,Jmin=Jmin,
                          min_absPhi_over_rho=qmin,delta_deg=round(np.rad2deg(delta),2),
                          a_T=aT,a_Phi=aP,LHS=L_,RHS=R_,ratio_LHS_over_RHS=L_/R_,
                          LHS_over_mu_M_L=L_/(mm*M*LL),HOLDS=bool(L_<=R_)))
out['pairC']=rowsC


def show(key):
    print('\n==',key)
    hdr=('pair','lam','mu_meas','muJ','minJ','min|Phi|/rho','delta','a_T','a_Phi','LHS','RHS','LHS/RHS','LHS/(mu M L)','HOLDS')
    print(('{:>4}{:>6}{:>10}{:>10}{:>9}{:>13}{:>7}{:>10}{:>10}{:>12}{:>10}{:>10}{:>14}{:>7}').format(*hdr))
    for x in out[key]:
        print(('{:>4}{:>6.2f}{:>10.5f}{:>10.5f}{:>9.4f}{:>13.5f}{:>7.1f}{:>10.5f}{:>10.5f}{:>12.3e}{:>10.4f}{:>10.2e}{:>14.5f}{:>7}').format(
            x['pair'],x['lam'],x['mu_meas'],x['muJ_meas'],x['Jmin'],x['min_absPhi_over_rho'],x['delta_deg'],
            x['a_T'],x['a_Phi'],x['LHS'],x['RHS'],x['ratio_LHS_over_RHS'],x['LHS_over_mu_M_L'],x['HOLDS']))
print('calibration (Phi=T_lambda) vs Lemma 1:')
for c in out['calibration_T_lambda_vs_Lemma1']: print('  ',c)
print('quadrature convergence:',out['quadrature_convergence'])
for k in ('pairA','pairB','pairC'): show(k)
print('\nsmall-mu slope (pair A):')
for r_ in out['pairA_small_mu_slope']: print('  ',r_)

json.dump(out,open('t3_results.json','w'),indent=1)
