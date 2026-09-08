"""
t6_stepwise_check.py -- verify EACH step of LEMMA T separately, not just the final inequality.

S2+S3 (pointwise):  |Ktil(Phi(x)) - Ktil(T_lam x)| <= C_K (|x|(q_lam(phi)-mu))^{-5} * mu |x|
S4:                 |Term I|  = |INT eta0 [Ktil(Phi)-Ktil(T_lam)] J_Phi|  <= 3(1+muJ) lam^2 I1 mu M L
S5:                 |Term II| = |INT eta0 Ktil(T_lam)(J_Phi - lam^2)|      <= (lam/2) muJ M L
S6:                 Term I + Term II == Delta   (identity, checks the split itself)
"""
import json, numpy as np, mpmath as mp
mp.mp.dps=30
exec(open('t3_maps_check.py').read().split("out={}")[0])   # same maps, same quadrature, verbatim
CK = 3/(2*np.pi**2)

def grid(delta,nu=160,np_=160):
    U,WU = nodes(0.0,LL,nu)
    brk=[0.0]+([delta,np.pi-delta] if delta>0 else [])+[np.pi/2,np.pi]; brk=sorted(set(brk))
    P=[];WP=[]
    for a_,b_ in zip(brk[:-1],brk[1:]):
        p,w=nodes(a_,b_,np_); P.append(p); WP.append(w)
    P=np.concatenate(P);WP=np.concatenate(WP)
    UU,PP=np.meshgrid(U,P,indexing='ij'); WW=np.outer(WU,WP)
    return UU,PP,WW

def pieces(Ff,Gf,Jf,lam,mu,delta=0.0):
    UU,PP,WW = grid(delta)
    rho=RHO0*np.exp(UU); s,c=np.sin(PP),np.cos(PP); r,z=rho*s,rho*c
    h=np.ones_like(PP) if delta<=0 else np.minimum(1.0,np.minimum(PP,np.pi-PP)/delta)
    eta0=-M*np.sign(c)*h/(rho*s)
    def bc(v): return np.broadcast_to(np.asarray(v,float),PP.shape)
    Fp,Gp,Jp = bc(Ff(r,z,lam,mu)),bc(Gf(r,z,lam,mu)),bc(Jf(r,z,lam,mu))
    Ft,Gt    = lam*r, z/lam**2
    Kp = 0.375/np.pi**2*(-Gp)*np.hypot(Fp,Gp)**-5
    Kt = 0.375/np.pi**2*(-Gt)*np.hypot(Ft,Gt)**-5
    meas = 2*np.pi**2*rho**4*s**3*rho          # dx5 with drho = rho du
    TI  = float(np.sum(WW*eta0*(Kp-Kt)*Jp*meas))
    TII = float(np.sum(WW*eta0*Kt*(Jp-lam**2)*meas))
    # pointwise step S2+S3
    q  = np.sqrt(lam**2*s**2+lam**-4*c**2)
    ptw_lhs = np.abs(Kp-Kt)
    ptw_rhs = CK*(rho*(q-mu))**-5 * mu*rho
    return TI,TII,float(np.max(ptw_lhs/ptw_rhs)),float(np.min(q-mu))

rows=[]
for name,(Ff,Gf,Jf) in (('A',(A_F,A_G,A_J)),('B',(B_F,B_G,B_J)),('C',(C_F,C_G,C_J))):
    for lam,mu in ((1.0,0.10),(1.25,0.05),(1.5,0.02)):
        p = mu
        if name=='B':   # eps chosen in t3 so that measured mu hits the target
            p = {(1.0,0.10):0.0555157470703125,(1.25,0.05):0.03691432530973254,(1.5,0.02):0.013334012241165191}[(lam,mu)]
        # measure mu, muJ on a fine grid (same as t3)
        UU,PP,_ = grid(0.0,400,400)
        rho=RHO0*np.exp(UU); s,c=np.sin(PP),np.cos(PP); r,z=rho*s,rho*c
        def bc(v): return np.broadcast_to(np.asarray(v,float),PP.shape)
        Fp,Gp,Jp = bc(Ff(r,z,lam,p)),bc(Gf(r,z,lam,p)),bc(Jf(r,z,lam,p))
        mm = float(np.max(np.hypot(Fp-lam*r,Gp-z/lam**2)/rho)); mJ=float(np.max(np.abs(Jp/lam**2-1)))
        for delta in (0.0,np.deg2rad(7.5)):
            TI,TII,ptw,qmin = pieces(Ff,Gf,Jf,lam,p,delta)
            aT = a_functional(T_F,T_G,T_J,lam,0.0,delta); aP=a_functional(Ff,Gf,Jf,lam,p,delta)
            bI  = 3*(1+mJ)*lam**2*I1(lam,mm)*mm*M*LL
            bII = lam/2*mJ*M*LL
            rows.append(dict(pair=name,lam=lam,mu=mm,muJ=mJ,delta_deg=round(np.rad2deg(delta),2),
                             pointwise_max_ratio_S2S3=ptw, min_q_minus_mu=qmin,
                             TermI=TI,bound_TermI=bI,TermI_ok=bool(abs(TI)<=bI),
                             TermII=TII,bound_TermII=bII,TermII_ok=bool(abs(TII)<=bII),
                             split_residual=abs((TI+TII)-(aP-aT))))
print(json.dumps(dict(rows=rows,
      all_pointwise_ok=bool(all(r['pointwise_max_ratio_S2S3']<=1 for r in rows)),
      all_TermI_ok=bool(all(r['TermI_ok'] for r in rows)),
      all_TermII_ok=bool(all(r['TermII_ok'] for r in rows)),
      max_split_residual=max(r['split_residual'] for r in rows)),indent=1))
json.dump(rows,open('t6_results.json','w'),indent=1)
