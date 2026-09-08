"""s1 -- the kernel facts that separate Kukavica's Lemma 3.1 from BFG's (10).

Every number in NOTE.md tagged [K#] comes from here.  n = 3 throughout unless said.
G(x,t) = (4 pi t)^{-n/2} exp(-|x|^2/(4t)) is the heat kernel.
"""
import numpy as np
from scipy import integrate
import json, math

out = {}
n = 3

# ---- K1: mean of the two kernels over R^3 at fixed tau -------------------
# G is a probability density; dG/dx_j is odd in x_j hence has integral 0.
def G(r, t):  return (4*math.pi*t)**(-n/2)*math.exp(-r*r/(4*t))
def dG_r(r, t):  # radial derivative magnitude |grad G| = -dG/dr
    return (4*math.pi*t)**(-n/2)*math.exp(-r*r/(4*t))*(r/(2*t))

for tau in (1.0, 1e-2, 1e-4):
    I_G   = integrate.quad(lambda r: 4*math.pi*r*r*G(r,tau), 0, np.inf, limit=400)[0]
    I_aG  = integrate.quad(lambda r: 4*math.pi*r*r*dG_r(r,tau), 0, np.inf, limit=400)[0]
    # signed integral of d_1 G over R^3, done as a 1-d radial reduction:
    # int d_1G dx = int (-x_1/(2t)) G dx = 0 by oddness; check to machine eps by 3-d quadrature on a big box
    out[f"int_G_tau{tau:g}"]      = I_G
    out[f"int_absgradG_tau{tau:g}"] = I_aG
    out[f"two_over_sqrt_pi_tau_{tau:g}"] = 2/math.sqrt(math.pi*tau)

# signed mean-zero check for d_1 G by direct 3-d Monte-Carlo-free quadrature (Gauss-Hermite)
xg, wg = np.polynomial.hermite_e.hermegauss(120)   # weight exp(-x^2/2)
tau = 1.0
s = math.sqrt(2*tau)          # G ~ N(0, 2 tau I)
# int d_1 G dx = E[ -X_1/(2 tau) ] under the N(0,2tau) law  = 0
E = float(np.sum(wg*(-(xg*s)/(2*tau)))/np.sum(wg))
out["signed_int_d1G_gausshermite"] = E

# ---- K2: sharp constants in the pointwise bounds -------------------------
# Kukavica p.6:  |grad G(x,t)| <= C /(|x|+t^{1/2})^{n+1}
# BFG p.10:      G(x,t)        <= c sqrt(t)/(|x|+sqrt t)^{4}
# By parabolic scaling both sups may be evaluated at t = 1.
rr = np.linspace(0, 80, 800001)[1:]
C_grad = float(np.max(dG_r(1.0*rr, 1.0)*(rr+1)**(n+1))) if False else float(np.max(
    [(4*math.pi)**(-n/2)*math.exp(-r*r/4)*(r/2)*(r+1)**(n+1) for r in rr]))
C_G    = float(np.max([(4*math.pi)**(-n/2)*math.exp(-r*r/4)*(r+1)**(n+1) for r in rr]))
out["sharp_C_gradG_weight_(|x|+sqrt t)^4"] = C_grad
out["sharp_c_G_weight_sqrt_t_over_(|x|+sqrt t)^4"] = C_G

# ---- K3: the Duhamel time-scalings the two kernels produce ---------------
# Kukavica's kernel bound:  I_K(t) = int_0^t int_{R^3} (|z|+sqrt(t-s))^{-4} dz ds
# BFG's kernel bound:       I_B(t) = int_0^t int_{R^3} sqrt(t-s) (|z|+sqrt(t-s))^{-4} dz ds
def inner(a):   # int_{R^3} (|z|+a)^{-4} dz , a>0
    return integrate.quad(lambda r: 4*math.pi*r*r/(r+a)**4, 0, np.inf, limit=400)[0]
for t in (0.25, 1.0, 4.0):
    IK = integrate.quad(lambda s: inner(math.sqrt(t-s)), 0, t, limit=400)[0]
    IB = integrate.quad(lambda s: math.sqrt(t-s)*inner(math.sqrt(t-s)), 0, t, limit=400)[0]
    out[f"I_K_t{t:g}"] = IK; out[f"I_K_over_sqrt_t_t{t:g}"] = IK/math.sqrt(t)
    out[f"I_B_t{t:g}"] = IB; out[f"I_B_over_t_t{t:g}"] = IB/t
out["closed_form_I_K_coeff_8pi_over_3"] = 8*math.pi/3
out["closed_form_I_B_coeff_4pi_over_3"] = 4*math.pi/3

# ---- K4: the two weights in the BMO inequality are equivalent ------------
# Kukavica p.7 uses 1/(|x|^{n+1}+1);  BFG p.10 uses 1/(|x|+1)^{n+1}.
ratio = [( (r+1)**4 )/( r**4 + 1 ) for r in rr]
out["weight_ratio_sup_(|x|+1)^4_over_(|x|^4+1)"] = float(np.max(ratio))
out["weight_ratio_inf"] = float(np.min(ratio))

# ---- K5: subtraction is exactly what carries the inequality --------------
# h_R(x) = (1 - log(1+|x|)/log R)_+ ,  ||h_R||_BMO <= C/log R.
# unsubtracted:  J0(R) = int |h_R| /(|x|^4+1) dx     -- must blow up rel. to 1/log R
# subtracted:    J1(R) = int |h_R - avg_{B_1} h_R|/(|x|^4+1) dx  -- must stay O(1/log R)
def hR(r, R):
    v = 1 - math.log1p(r)/math.log(R)
    return v if v > 0 else 0.0
rows = []
for R in (1e2, 1e4, 1e8, 1e16, 1e32, 1e64):
    avg = integrate.quad(lambda r: 3*r*r*hR(r,R), 0, 1, limit=400)[0]   # mean over unit ball
    J0 = integrate.quad(lambda r: 4*math.pi*r*r*hR(r,R)/(r**4+1), 0, np.inf, limit=600)[0]
    J1 = integrate.quad(lambda r: 4*math.pi*r*r*abs(hR(r,R)-avg)/(r**4+1), 0, np.inf, limit=600)[0]
    lr = math.log(R)
    rows.append(dict(R=R, log_R=lr, avg_B1=avg, J0=J0, J0_times_logR=J0*lr,
                     J1=J1, J1_times_logR=J1*lr))
out["subtraction_table"] = rows

print(json.dumps(out, indent=1))
json.dump(out, open("scripts/s1_results.json","w"), indent=1)
