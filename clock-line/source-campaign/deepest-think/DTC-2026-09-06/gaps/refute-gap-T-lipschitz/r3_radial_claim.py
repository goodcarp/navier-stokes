"""r3 - two secondary claims of the seat's NOTE, tested.

(i) NOTE sec 3(N) "Consequence for the campaign: the O(c/L) SHEAR BETWEEN SHELLS in
    prove-lagrangian sec 4(2) is a RADIAL REARRANGEMENT -- the least dangerous of the three error
    terms listed there."
    A radial rearrangement is Psi(u) = s(|u|) u : it moves points ALONG the ray, polar angle fixed.
    The shear between shells is lambda = lambda(rho), i.e. Lambda(x) = T_{lambda(rho)} x, and
    T_lambda moves the polar angle by  tan(phi') = lambda^3 tan(phi).  Measured below: the angular
    part of the displacement, and the fraction of the total displacement that is NOT radial.

(ii) NOTE sec 3(N)/(A) "if s is log-periodic ... Delta = 0 EXACTLY -- for every lambda, every datum,
     and every amplitude mu (however large)" and sec 5 "unchanged past the diffeomorphism limit
     mu = 0.9342".  The seat's own s5 run has a k=2, mu=3.0 row.  Re-run here independently.
"""
import json, numpy as np
KTH = 2*(1-np.sqrt(2.0/3.0))
def lam_of_sigma(s): return (1-(1-s)*KTH/2)**-2

out = {}
# ---- (i) how much of the shell-dependent strain is angular? ----
sig = np.linspace(0,1,2001); ph = np.linspace(1e-6, np.pi-1e-6, 2001)
S,P = np.meshgrid(sig, ph, indexing='ij')
lam = lam_of_sigma(S)
lb = 1.2247448713915889                       # sqrt(3/2), the log-mean; best-value single lambda
# image of the unit sphere direction under T_lam vs under T_lb
def img(l, p): return l*np.sin(p), l**-2*np.cos(p)
r1,z1 = img(lam,P); r0,z0 = img(lb,P)
n1 = np.hypot(r1,z1); n0 = np.hypot(r0,z0)
disp = np.hypot(r1-r0, z1-z0)
radial = np.abs(n1-n0)                         # purely-radial part of the displacement
angular = np.sqrt(np.maximum(disp**2 - radial**2, 0.0))
out['sup_rel_displacement'] = float(np.max(disp/n0))
out['sup_rel_radial_part']  = float(np.max(radial/n0))
out['sup_rel_angular_part'] = float(np.max(angular/n0))
out['max_angular_fraction'] = float(np.max(angular/np.maximum(disp,1e-300)))
# polar-angle change: tan phi' = lambda^3 tan phi
ph2 = np.arctan2(lam*np.sin(P), lam**-2*np.cos(P))
ph2b = np.arctan2(lb*np.sin(P), lb**-2*np.cos(P))
out['max_polar_angle_shift_deg'] = float(np.max(np.abs(ph2-ph2b))*180/np.pi)
out['lambda_cubed_ratio_extremes'] = [float(lam_of_sigma(0.0)**3), float(lam_of_sigma(1.0)**3)]

# ---- (ii) the log-periodic null direction past the diffeomorphism limit ----
M, rho0, R = 1.0, 1.0, 4096.0
L = np.log(R/rho0)
def gl(a,b,n):
    x,w = np.polynomial.legendre.leggauss(n); return 0.5*(b-a)*x+0.5*(a+b), 0.5*(b-a)*w
def panels(pts,n):
    xs,ws=[],[]
    for a,b in zip(pts[:-1],pts[1:]):
        x,w=gl(a,b,n); xs.append(x); ws.append(w)
    return np.concatenate(xs), np.concatenate(ws)
def a_ripple(mu,k,n_phi=300,n_tau=600):
    ph,wph = panels([0.0,np.pi/2,np.pi],n_phi); ta,wta = gl(0,1,n_tau)
    PH,TA = np.meshgrid(ph,ta,indexing='ij'); WP,WT = np.meshgrid(wph,wta,indexing='ij')
    vr = rho0*np.exp(L*TA)                                   # lambda = 1 so |T x| = rho
    th = k*np.pi*np.log(vr/rho0)/L
    s = 1+mu*np.sin(th); vs = mu*(k*np.pi/L)*np.cos(th)
    om = -M*np.sign(np.cos(PH))
    # (Z/|Phi|^5) J = z_u (s + vs)/varrho^5  (the s^4 cancellation)
    integ = np.cos(PH)*np.sin(PH)**2*om*(s+vs)
    return -(3.0/4.0)*L*np.sum(integ*WP*WT)
aT = M/2*L
rows=[]
for mu in (0.05,0.5,0.9,0.9342,1.0,1.2,3.0):
    for k in (2,4):
        a = a_ripple(mu,k)
        smin = 1-mu   # min of s over the shell
        vsmin = -mu*(k*np.pi/L)
        rows.append(dict(mu=mu,k=k,a=a,aT=aT,rel=abs(a-aT)/aT,
                         is_diffeo=bool(mu*np.sqrt(1+(k*np.pi/L)**2)<1),
                         s_positive=bool(smin>0)))
out['null_rows']=rows
out['max_rel_diffeo']    = max(r['rel'] for r in rows if r['is_diffeo'])
out['max_rel_nondiffeo'] = max(r['rel'] for r in rows if not r['is_diffeo'])
json.dump(out, open('r3_results.json','w'), indent=1)
for k,v in out.items():
    if k!='null_rows': print(f"{k:32s} {v}")
print("\nlog-periodic ripple, lambda=1, bang-bang cap:")
print(f"{'mu':>7s} {'k':>3s} {'diffeo':>7s} {'s>0':>5s} {'a':>12s} {'rel |a-aT|/aT':>15s}")
for r in rows:
    print(f"{r['mu']:7.4f} {r['k']:3d} {str(r['is_diffeo']):>7s} {str(r['s_positive']):>5s} "
          f"{r['a']:12.7f} {r['rel']:15.3e}")
