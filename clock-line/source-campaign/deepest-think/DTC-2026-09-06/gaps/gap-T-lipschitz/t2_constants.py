"""
t2_constants.py -- the explicit constants of LEMMA T.

  q_l(p) = sqrt(lambda^2 sin^2 p + lambda^{-4} cos^2 p)   (so |T_lambda x| = |x| q_l(phi))
  m(lambda) = min_p q_l(p) = min(lambda, lambda^{-2})
  I1(lambda,mu) = INT_0^pi sin^2 p (q_l(p) - mu)^{-5} dp        (finite iff mu < m(lambda))

LEMMA T bound:
  |Delta| <= [ 3 (1+muJ) lambda^2 I1(lambda,mu) ] mu M L  +  (lambda/2) muJ M L
Corollary (muJ <= mu):
  |Delta| <= C(lambda,mu) mu M L ,  C = 3(1+mu) lambda^2 I1(lambda,mu) + lambda/2 .
Everything here is computed with mpmath (50 dps) and cross-checked with a coarse/fine
Gauss-Legendre pair written from scratch.
"""
import json, numpy as np, mpmath as mp
mp.mp.dps = 50

def qlam(p, lam): return mp.sqrt(lam**2*mp.sin(p)**2 + lam**-4*mp.cos(p)**2)
def mlam(lam):    return min(mp.mpf(lam), mp.mpf(lam)**-2)

def I1(lam, mu):
    lam = mp.mpf(lam); mu = mp.mpf(mu)
    if mu >= mlam(lam): return mp.inf
    f = lambda p: mp.sin(p)**2*(qlam(p,lam)-mu)**-5
    return mp.quad(f, [0, mp.pi/2, mp.pi])

def I1_gl(lam, mu, n):                      # independent quadrature, my own nodes
    x,w = np.polynomial.legendre.leggauss(n); tot = 0.0
    for a,b in ((0.0,np.pi/2),(np.pi/2,np.pi)):
        p = 0.5*(b-a)*x+0.5*(a+b); ww = 0.5*(b-a)*w
        q = np.sqrt(lam**2*np.sin(p)**2 + lam**-4*np.cos(p)**2)
        tot += np.sum(ww*np.sin(p)**2*(q-mu)**-5)
    return tot

def C_of(lam, mu, muJ=None):
    muJ = mu if muJ is None else muJ
    return float(3*(1+muJ)*lam**2*I1(lam,mu) + lam/2*(muJ/mu if mu>0 else 0))

out = {}
out['I1(1,0)_vs_pi/2'] = [float(I1(1,0)), float(mp.pi/2), float(abs(I1(1,0)-mp.pi/2))]
tab = []
for lam in (1.0,1.05,1.1,1.25,1.5):
    m = float(mlam(lam))
    for mu in (0.0,0.01,0.02,0.05,0.10,0.20):
        if mu >= m: continue
        i1 = float(I1(lam,mu))
        row = dict(lam=lam, mu=mu, m_lambda=m, I1=i1,
                   I1_gl400=I1_gl(lam,mu,400), I1_gl800=I1_gl(lam,mu,800),
                   C_disp_coef=3*(1+mu)*lam**2*i1,      # multiplies mu
                   C_jac_coef=lam/2,                     # multiplies muJ
                   C_total_if_muJ_eq_mu=3*(1+mu)*lam**2*i1 + lam/2)
        tab.append(row)
out['table'] = tab

# blow-up as mu -> m(lambda)^-  (the sharpness mechanism)
bl = []
for lam in (1.0,1.5):
    m = float(mlam(lam))
    for frac in (0.5,0.8,0.9,0.99,0.999):
        mu = frac*m
        bl.append(dict(lam=lam, mu=mu, mu_over_m=frac, I1=float(I1(lam,mu)),
                       C_total=3*(1+mu)*lam**2*float(I1(lam,mu))+lam/2))
out['blowup_as_mu_to_m'] = bl
out['headline'] = {
  'C(1,0)=3*pi/2+1/2': 3*float(mp.pi)/2+0.5,
  'C(1.25,0)': 3*1.25**2*float(I1(1.25,0))+1.25/2,
  'C(1.5,0)':  3*1.5**2*float(I1(1.5,0))+1.5/2,
  'C_K=3/(2pi^2)': 3/(2*float(mp.pi)**2),
}
print(json.dumps(out, indent=1))
json.dump(out, open('t2_results.json','w'), indent=1)
