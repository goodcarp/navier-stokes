"""r1 -- independent recomputation of every [K#] number in gaps/kukavica-at-source,
by routes different from the original script, plus controls that CAN fail.

n = 3.  G(x,t) = (4 pi t)^{-3/2} exp(-|x|^2/(4t)).
"""
import math, json
import numpy as np
from scipy import integrate, optimize

out = {}
n = 3
G   = lambda r,t: (4*math.pi*t)**(-n/2)*math.exp(-r*r/(4*t))
aG  = lambda r,t: (4*math.pi*t)**(-n/2)*math.exp(-r*r/(4*t))*(r/(2*t))

# ---------- 1. K1 : masses -------------------------------------------------
# analytic: int G dx = 1 ; int |grad G| dx = 2/sqrt(pi t).
# check by a DIFFERENT quadrature: substitution r = a*tan(theta), theta in (0,pi/2)
def quad_tan(f, a=1.0):
    g = lambda th: f(a*math.tan(th))*a/math.cos(th)**2
    return integrate.quad(g, 0, math.pi/2, limit=800)[0]
for tau in (1.0, 1e-2, 1e-4):
    a = math.sqrt(tau)
    out[f"int_G_tau{tau:g}_tanquad"]      = quad_tan(lambda r: 4*math.pi*r*r*G(r,tau), a)
    out[f"int_absgradG_tau{tau:g}_tanquad"]= quad_tan(lambda r: 4*math.pi*r*r*aG(r,tau), a)
    out[f"analytic_2_over_sqrt_pi_tau{tau:g}"] = 2/math.sqrt(math.pi*tau)

# ---------- 1b. CONTROL THAT CAN FAIL: mean-zero of d_1 G ------------------
# original script checked E[-X_1/(2tau)] on SYMMETRIC Gauss-Hermite nodes -> 0 by
# construction.  Here: tensor Gauss-Legendre on a DELIBERATELY ASYMMETRIC box.
def box_int(f, lo, hi, m=140):
    x,w = np.polynomial.legendre.leggauss(m)
    def mp(l,h): return 0.5*(h-l)*x+0.5*(h+l), 0.5*(h-l)*w
    X,WX = mp(lo[0],hi[0]); Y,WY = mp(lo[1],hi[1]); Z,WZ = mp(lo[2],hi[2])
    XX,YY,ZZ = np.meshgrid(X,Y,Z, indexing='ij')
    WW = WX[:,None,None]*WY[None,:,None]*WZ[None,None,:]
    return float(np.sum(WW*f(XX,YY,ZZ)))
t = 1.0
gfun  = lambda X,Y,Z: (4*math.pi*t)**(-1.5)*np.exp(-(X**2+Y**2+Z**2)/(4*t))
d1fun = lambda X,Y,Z: -(X/(2*t))*gfun(X,Y,Z)
LO, HI = (-6.0,-6.5,-7.0), (7.0,6.5,6.0)      # asymmetric box, no cancellation by symmetry
out["asym_box_int_d1G"] = box_int(d1fun, LO, HI)
out["asym_box_int_G_CONTROL_should_be_1"] = box_int(gfun, LO, HI)   # same machinery, must give 1

# ---------- 2. K2 : sharp constants, by an optimizer not a grid ------------
f_grad = lambda r: -aG(r,1.0)*(r+1)**(n+1)
f_G    = lambda r: -G(r,1.0)*(r+1)**(n+1)
r1 = optimize.minimize_scalar(f_grad, bracket=(0.5,2.0,6.0), method='brent', options={'xtol':1e-14})
r2 = optimize.minimize_scalar(f_G,    bracket=(0.5,2.0,6.0), method='brent', options={'xtol':1e-14})
out["sharp_C_gradG_brent"] = -r1.fun; out["argmax_gradG"] = r1.x
out["sharp_c_G_brent"]     = -r2.fun; out["argmax_G"]     = r2.x
out["grid_value_reported_gradG"] = 0.9231098778782552
out["grid_value_reported_G"]     = 0.7109819926171765
out["digits_agree_gradG"] = -math.log10(abs(-r1.fun-0.9231098778782552)/abs(r1.fun))
out["digits_agree_G"]     = -math.log10(abs(-r2.fun-0.7109819926171765)/abs(r2.fun))

# ---------- 3. K3 : Duhamel scalings, CLOSED FORM derivation ---------------
# int_{R^3}(|z|+a)^{-4}dz = 4pi int r^2/(r+a)^4 dr = 4pi/(3a)  [substitute r=a u]
# => I_K(t) = int_0^t 4pi/(3 sqrt(t-s)) ds = 8 pi sqrt(t)/3
#    I_B(t) = int_0^t sqrt(t-s) 4pi/(3 sqrt(t-s)) ds = 4 pi t/3
out["inner_check_a1_numeric"] = quad_tan(lambda r: 4*math.pi*r*r/(r+1.0)**4, 1.0)
out["inner_check_a1_closed_4pi_over_3"] = 4*math.pi/3
for T in (0.25,1.0,4.0):
    IK = integrate.quad(lambda s: 4*math.pi/(3*math.sqrt(T-s)), 0, T, limit=400)[0]
    IB = integrate.quad(lambda s: math.sqrt(T-s)*4*math.pi/(3*math.sqrt(T-s)), 0, T, limit=400)[0]
    out[f"I_K_t{T:g}_closedinner"] = IK; out[f"I_K_closedform_t{T:g}"] = 8*math.pi*math.sqrt(T)/3
    out[f"I_B_t{T:g}_closedinner"] = IB; out[f"I_B_closedform_t{T:g}"] = 4*math.pi*T/3

# ---------- 4. K4 : weight equivalence, exactly ---------------------------
w = lambda r: (r+1)**4/(r**4+1)
rs = optimize.minimize_scalar(lambda r: -w(r), bracket=(0.2,1.0,3.0), method='brent')
out["weight_sup_brent"] = -rs.fun; out["weight_argmax"] = rs.x
out["weight_at_0_exact_inf"] = w(0.0)
out["weight_limit_at_infinity"] = w(1e12)
out["weight_inf_reported_in_NOTE"] = 1.0
out["weight_inf_actually_in_s1_results_json"] = 1.000400060004

# ---------- 5. K5 : the h_R table, independent quadrature + limits --------
def hR(r,R):
    v = 1 - math.log1p(r)/math.log(R)
    return v if v>0 else 0.0
kappa = 2*math.log(2) - 5/6            # = int_0^1 3 s^2 log(1+s) ds, closed form
out["kappa_closed_2log2_minus_5_6"] = kappa
out["kappa_numeric"] = integrate.quad(lambda s: 3*s*s*math.log1p(s),0,1,limit=400)[0]
rows=[]
for R in (1e2,1e4,1e8,1e16,1e32,1e64):
    lr = math.log(R); supp = R-1.0
    avg = integrate.quad(lambda s: 3*s*s*hR(s,R),0,1,limit=400)[0]
    # independent quadrature: split at the support radius, tan-substitution on the tail
    J0 = (integrate.quad(lambda r: 4*math.pi*r*r*hR(r,R)/(r**4+1),0,min(supp,1e6),limit=900)[0]
          + (integrate.quad(lambda r: 4*math.pi*r*r*hR(r,R)/(r**4+1),1e6,supp,limit=900)[0] if supp>1e6 else 0.0))
    J1 = (integrate.quad(lambda r: 4*math.pi*r*r*abs(hR(r,R)-avg)/(r**4+1),0,min(supp,1e6),limit=900)[0]
          + (integrate.quad(lambda r: 4*math.pi*r*r*abs(hR(r,R)-avg)/(r**4+1),1e6,supp,limit=900)[0] if supp>1e6 else 0.0)
          + integrate.quad(lambda r: 4*math.pi*r*r*avg/(r**4+1),supp,np.inf,limit=900)[0])
    rows.append(dict(R=R,log_R=lr,support_radius=supp,avg=avg,
                     avg_predicted_1_minus_kappa_over_logR=1-kappa/lr,
                     J0=J0,J0_logR=J0*lr,J1=J1,J1_logR=J1*lr))
out["h_table_independent"]=rows
out["J0_limit_closed_4pi*pi/(2sqrt2)"] = 4*math.pi*math.pi/(2*math.sqrt(2))
out["J0_limit_numeric"] = quad_tan(lambda r: 4*math.pi*r*r/(r**4+1),1.0)
out["J1_logR_limit_predicted"] = quad_tan(lambda r: 4*math.pi*r*r*abs(math.log1p(r)-kappa)/(r**4+1),1.0)

# ---------- 6. the one-line counterexample the NOTE never states ----------
# f == 1 : ||f||_BMO = 0 but the unsubtracted LHS is 4pi/3 > 0.
out["unsub_LHS_for_f_equal_1"] = quad_tan(lambda r: 4*math.pi*r*r/(r+1)**4,1.0)
out["unsub_LHS_for_f_equal_1_closed_4pi_over_3"] = 4*math.pi/3
# and h_R IS compactly supported, so it satisfies BFG's stated decay hypothesis.

# ---------- 7. BMO norm of h_R: the hypothesis the NOTE only asserts ------
# radial h(x)=phi(|x|).  Ball B(a e1, rho).  Radial density of Lebesgue measure
# restricted to that ball: 2 pi r^2 (1 - mu),  mu = clip((r^2+a^2-rho^2)/(2 a r),-1,1).
def mean_osc(phi, a, rho, m=400):
    lo = max(0.0, a-rho); hi = a+rho
    x,wq = np.polynomial.legendre.leggauss(m)
    r = 0.5*(hi-lo)*x+0.5*(hi+lo); wr = 0.5*(hi-lo)*wq
    if a < 1e-14:
        dens = 4*math.pi*r**2
    else:
        mu = np.clip((r**2+a**2-rho**2)/(2*a*r),-1,1)
        dens = 2*math.pi*r**2*(1-mu)
    vals = np.array([phi(float(rr)) for rr in r])
    Z = float(np.sum(wr*dens)); 
    if Z <= 0: return 0.0
    mean = float(np.sum(wr*dens*vals))/Z
    return float(np.sum(wr*dens*np.abs(vals-mean)))/Z
bmo_rows=[]
for R in (1e2,1e4,1e8,1e16):
    phi = lambda r,R=R: hR(r,R)
    best=0.0; arg=None
    A = [0.0]+[10**e for e in np.linspace(-2, math.log10(R), 60)]
    P = [10**e for e in np.linspace(-2, math.log10(R)+0.5, 60)]
    for a in A:
        for rho in P:
            v = mean_osc(phi,a,rho,m=200)
            if v>best: best, arg = v,(a,rho)
    bmo_rows.append(dict(R=R,log_R=math.log(R),bmo_lower_bound=best,
                         times_logR=best*math.log(R),argmax_a_rho=arg))
out["bmo_norm_search_h_R"]=bmo_rows

print(json.dumps(out, indent=1, default=float))
json.dump(out, open("~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/gaps/refute-kukavica-at-source/scripts/r1_results.json","w"), indent=1, default=float)
