"""
t4_controls.py -- the two controls for GAP T.

CONTROL C1 (harness control -- MUST FIRE).  Re-run the pair-A/B/C checks against the
deliberately FALSE claim |Delta| <= (C/50) mu M L.  If nothing fails, the checker is not
diagnostic and the pass in t3 means nothing.

CONTROL C2 (sharpness -- the estimate MUST FAIL when mu is not small).
   Phi_mu := T_lambda - mu rho0 e_z      (a rigid translation composed with T_lambda)
   * smooth, injective, J_Phi = lambda^2 EXACTLY  (muJ = 0: volume-exact in the 5D lift)
   * sup_{x in S} |Phi_mu(x) - T_lambda(x)| / |x| = mu   (attained at |x| = rho0)
   * dist(0, Phi_mu(S-bar)) = (lambda^{-2} - mu) rho0  for mu < lambda^{-2} = min_phi q_lambda(phi)
   For the extremal datum eta0 = -M sgn(z)/r (|omega^theta| = M, admitted by the hypothesis)
   the perturbed functional blows up logarithmically as mu -> lambda^{-2} and is +infinity at
   mu = lambda^{-2}.  Hence NO bound |Delta| <= C mu M L with a fixed C survives, and the
   threshold mu < min_phi q_lambda(phi) of LEMMA T is SHARP.

Exact evaluator used for C2 (no singular quadrature).  With y = Phi_mu(x), y = (t sin th, t cos th):
   a = (3/4) M lambda INT_0^pi cos th sin^2 th * W(th) dth ,
   W(th) = INT sgn(z_x) 1_{x in S} dlog t ,      x(t,th) = (t sin th/lambda, lambda^2(t cos th + mu rho0))
W(th) is computed in closed form (the t-set is cut out by two quadratics in t).
"""
import json, numpy as np, mpmath as mp
mp.mp.dps = 30
M, RHO0, LL = 1.0, 1.0, 6.0
R = RHO0*np.exp(LL)

# ---------------------------------------------------------------- exact evaluator for C2
def t_roots(alpha,beta,gamma,c2):
    disc = beta*beta - 4*alpha*(gamma-c2)
    if disc <= 0: return []
    s = mp.sqrt(disc)
    return sorted([r for r in ((-beta-s)/(2*alpha), (-beta+s)/(2*alpha)) if r > 0])

def W_of_theta(th, lam, mu):
    """INT sgn(z_x) 1_{x in S} dlog t , in closed form."""
    th = mp.mpf(th); lam = mp.mpf(lam); mu = mp.mpf(mu)
    s,c = mp.sin(th), mp.cos(th)
    alpha = s**2/lam**2 + lam**4*c**2
    beta  = 2*lam**4*mu*RHO0*c
    gamma = lam**4*mu**2*RHO0**2
    inn = t_roots(alpha,beta,gamma,mp.mpf(RHO0)**2)   # |x| = rho0
    out = t_roots(alpha,beta,gamma,mp.mpf(R)**2)      # |x| = R
    # build the set {t>0 : rho0^2 < alpha t^2 + beta t + gamma < R^2}
    pts = sorted(set([mp.mpf(0)]+inn+out))
    segs=[]
    for a_,b_ in zip(pts[:-1],pts[1:]):
        m = (a_+b_)/2
        v = alpha*m*m+beta*m+gamma
        if RHO0**2 < v < R**2: segs.append((a_,b_))
    if pts:
        m = pts[-1]*2+1
        v = alpha*m*m+beta*m+gamma
        if RHO0**2 < v < R**2: segs.append((pts[-1], mp.inf))
    tot = mp.mpf(0)
    tstar = mu*RHO0/abs(c) if c < 0 else None   # sgn(z_x) flips here when cos th < 0
    for a_,b_ in segs:
        if b_ == mp.inf: raise RuntimeError('unbounded segment')
        cuts = [a_,b_] if (tstar is None or not (a_ < tstar < b_)) else [a_,tstar,b_]
        for u,v in zip(cuts[:-1],cuts[1:]):
            if u <= 0: raise RuntimeError('t=0 in segment (origin interior to image)')
            mid = (u+v)/2
            sg = mp.sign(mid*c + mu*RHO0)
            tot += sg*mp.log(v/u)
    return tot

def a_C2(lam, mu):
    f = lambda th: mp.cos(th)*mp.sin(th)**2*W_of_theta(th,lam,mu)
    br = [0, mp.pi/2, mp.pi]
    return float(mp.mpf(0.75)*M*lam*mp.quad(f, br))

out = {}
CLAM = {1.0:5.212388980384690, 1.25:6.657497095904176, 1.5:9.142247938917823}   # C(lambda,0) from t2

# 1. calibration of the exact evaluator: mu = 0 must give (M/2) lambda L
out['C2_evaluator_calibration'] = [dict(lam=l, a=a_C2(l,0.0), lemma1=0.5*l*LL) for l in (1.0,1.25,1.5)]

# 2. blow-up as mu -> lambda^{-2} = min_phi q_lambda(phi)
rows=[]
for lam in (1.0,1.5):
    thr = lam**-2
    for f in ('0.5','0.9','0.99','0.999','0.9999','0.99999'):
        mu = mp.mpf(f)*thr
        a = a_C2(lam,mu); a0 = 0.5*lam*LL
        rows.append(dict(lam=lam,threshold=thr,mu=float(mu),mu_over_thr=float(f),a_Phi=a,a_T=a0,
                         LHS=abs(a-a0), LHS_over_mu_M_L=abs(a-a0)/float(mu*M*LL),
                         C_lambda_0=CLAM[lam],
                         VIOLATES_frozen_constant=bool(abs(a-a0) > CLAM[lam]*float(mu)*M*LL)))
out['C2_blowup'] = rows

# 3. log-slope fit  a ~ (lambda/4) log(1/(thr-mu)) + const
fits={}
for lam in (1.0,1.5):
    thr=lam**-2; xs=[];ys=[]
    for k in (6,8,10,12):
        mu = thr*(1-mp.mpf(10)**(-k))
        xs.append(float(mp.log(1/(thr-mu)))); ys.append(a_C2(lam,mu))
    A,B = np.polyfit(xs,ys,1)
    fits[str(lam)] = dict(slope=float(A), predicted_slope=lam/4, intercept=float(B))
out['C2_log_slope_fit'] = fits

# 4. explicit violation  (needs high precision: mu = lambda^{-2}(1-1e-k)) of the frozen-constant claim |Delta| <= C(lambda,0) mu M L
vio=[]
for lam in (1.0,):
    thr=lam**-2
    for k in (20,60,100,140,160):
        mp.mp.dps = max(60, 2*k+40)
        mu = thr*(1-mp.mpf(10)**(-k))
        a = a_C2(lam,mu); a0=0.5*lam*LL; L_=abs(a-a0); R_=CLAM[lam]*float(mu)*M*LL
        vio.append(dict(lam=lam,k=k,mu_str=f'{float(thr)}*(1-1e-{k})',LHS=L_,frozen_RHS=R_,
                        LHS_over_mu_M_L=L_/float(mu*M*LL), FAILS=bool(L_>R_)))
out['C2_frozen_constant_violation']=vio
out['C2_frozen_constant_control_fires']=any(v['FAILS'] for v in vio)
mp.mp.dps = 30

# 5. the SAME map with an ADMISSIBLE tapered datum stays bounded (minimal-hypothesis remark)
def a_C2_taper_np(lam,mu,delta,tmin=1e-14,n=600):
    lam=float(lam); mu=float(mu)
    def W(th):
        s,c=np.sin(th),np.cos(th)
        alpha=s*s/lam**2+lam**4*c*c; beta=2*lam**4*mu*RHO0*c; gamma=lam**4*mu*mu*RHO0**2
        rts=[]
        for c2 in (RHO0**2,R**2):
            d=beta*beta-4*alpha*(gamma-c2)
            if d>0:
                sq=np.sqrt(d); rts += [r for r in ((-beta-sq)/(2*alpha),(-beta+sq)/(2*alpha)) if r>0]
        pts=sorted(set([tmin]+[r for r in rts if r>tmin])); segs=[]
        for a_,b_ in zip(pts[:-1],pts[1:]):
            m=0.5*(a_+b_); v=alpha*m*m+beta*m+gamma
            if RHO0**2<v<R**2: segs.append((a_,b_))
        tot=0.0
        x,w=np.polynomial.legendre.leggauss(n)
        for a_,b_ in segs:
            la,lb=np.log(a_),np.log(b_)
            lt=0.5*(lb-la)*x+0.5*(la+lb); ww=0.5*(lb-la)*w
            t=np.exp(lt); rx=t*s/lam; zx=lam**2*(t*c+mu*RHO0); rr=np.hypot(rx,zx)
            ph=np.arccos(np.clip(zx/rr,-1,1)); phax=np.minimum(ph,np.pi-ph)
            tot+=np.sum(ww*np.sign(zx)*np.minimum(1.0,phax/delta))
        return tot
    x,w=np.polynomial.legendre.leggauss(400); tot=0.0
    for a_,b_ in ((0,np.pi/2),(np.pi/2,np.pi)):
        th=0.5*(b_-a_)*x+0.5*(a_+b_); ww=0.5*(b_-a_)*w
        tot+=np.sum(ww*np.cos(th)*np.sin(th)**2*np.array([W(t) for t in th]))
    return 0.75*M*lam*tot
d75=np.deg2rad(7.5)
tap=[]
for f in (0.5,0.9,0.99,0.9999,1.0,1.5,3.0):
    mu=f*1.0
    tap.append(dict(lam=1.0,mu=mu,mu_over_thr=f,
                    a_taper_7p5deg=a_C2_taper_np(1.0,mu,d75),
                    a_taper_tmin1em20=a_C2_taper_np(1.0,mu,d75,tmin=1e-20)))
out['C2_tapered_datum_stays_bounded']=tap

# 6. C1 harness control (MUST FIRE)
t3 = json.load(open('t3_results.json'))
fired=[]
for key in ('pairA','pairB','pairC'):
    for x in t3[key]:
        false_rhs = x['RHS']/50.0
        fired.append(dict(pair=x['pair'],lam=x['lam'],mu=x['mu_meas'],delta_deg=x['delta_deg'],
                          LHS=x['LHS'],false_RHS=false_rhs,FIRES=bool(x['LHS']>false_rhs)))
out['C1_false_claim_rows']=fired
out['C1_control_fires']=any(r['FIRES'] for r in fired)
out['C1_n_fired']=sum(r['FIRES'] for r in fired)

print(json.dumps(out,indent=1))
json.dump(out,open('t4_results.json','w'),indent=1)
