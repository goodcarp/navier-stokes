"""
REFUTER independent recomputation. Deliberately different machinery from both prior scripts:
mpmath adaptive quadrature at 30 dps (not numpy Gauss-Legendre / Gauss-Hermite), plus exact
closed forms derived by hand and checked symbolically-by-value.

[R1] heat-kernel means in R^3 (analytic + quadrature)
[R2] radial closed form int_{R^3}(|z|+a)^{-4} dz and the two Duhamel scalings (analytic)
[R3] the weighted-BMO family h_R: unsubtracted J0 and subtracted J1, times log R
[R4] weight equivalence sup (r+1)^4/(r^4+1)
[R5] BMO-norm control of the family: is C in ||h_R||_BMO <= C/log R really R-independent?
     Numeric lower bound on ||h_R||_BMO by maximising the mean oscillation over a grid of
     concentric balls -- a control that CAN fail (if ||h_R||_BMO*log R grew, the family dies).
"""
import json, math
from mpmath import mp, mpf, quad, exp, sqrt, pi, log, inf, log10

mp.dps = 30
out = {}

# ---- [R1] ----
r1 = []
for tau in ('1', '0.01', '0.0001'):
    tau = mpf(tau)
    G   = lambda r: (4*pi*tau)**mpf('-1.5')*exp(-r*r/(4*tau))
    r1.append(dict(
        tau=str(tau),
        int_G      = str(4*pi*quad(lambda r: G(r)*r*r, [0, 20*sqrt(tau), inf])),
        int_absgrad= str(4*pi*quad(lambda r: G(r)*(r/(2*tau))*r*r, [0, 20*sqrt(tau), inf])),
        analytic_2_over_sqrt_pi_tau = str(2/sqrt(pi*tau)),
    ))
out['R1_heat_kernel'] = r1
# analytic: int |grad G| = 2/sqrt(pi tau) exactly; int d_1 G = 0 exactly (odd integrand)
out['R1_analytic_note'] = ("int G = 1; int d_1 G = 0 exactly by oddness in z_1; "
                           "int |grad G| dz = 2/sqrt(pi tau) derived by hand and matched below")

# ---- [R2] ----
r2 = []
for a in ('0.1', '1', '3'):
    a = mpf(a)
    num = 4*pi*quad(lambda r: r*r/(r+a)**4, [0, a, 10*a, inf])
    r2.append(dict(a=str(a), quad=str(num), closed_4pi_over_3a=str(4*pi/(3*a)),
                   rel_err=str(abs(num-4*pi/(3*a))/(4*pi/(3*a)))))
out['R2_radial'] = r2
# Duhamel, done exactly: int_0^t 4pi/(3 sqrt(t-s)) ds = 8pi sqrt(t)/3 ;
#                        int_0^t sqrt(t-s) 4pi/(3 sqrt(t-s)) ds = 4pi t/3
duh = []
for t in ('0.25', '1', '4'):
    t = mpf(t)
    A = quad(lambda s: 4*pi/(3*sqrt(t-s)), [0, t])
    B = quad(lambda s: sqrt(t-s)*4*pi/(3*sqrt(t-s)), [0, t])
    duh.append(dict(t=str(t), gradG_over_sqrt_t=str(A/sqrt(t)), G_over_t=str(B/t),
                    eight_pi_over_3=str(8*pi/3), four_pi_over_3=str(4*pi/3)))
out['R2_duhamel'] = duh

# ---- [R3] ----
r3 = []
for k in (2, 4, 8, 16, 32, 64):
    R = mpf(10)**k
    lg = log(R)
    h  = lambda r: max(mpf(0), 1 - log(1+r)/lg)
    avg = 3*quad(lambda r: h(r)*r*r, [0, 1])          # (1/|B1|) int_{B1} h = 3 int_0^1 h r^2 dr
    br = [0, 1, 10, R-1, inf]
    J0 = 4*pi*quad(lambda r: h(r)/(r**4+1)*r*r, br)
    J1 = 4*pi*quad(lambda r: abs(h(r)-avg)/(r**4+1)*r*r, br)
    r3.append(dict(R=f'1e{k}', log_R=str(lg), avg_B1=str(avg),
                   J0=str(J0), J0_logR=str(J0*lg), J1=str(J1), J1_logR=str(J1*lg)))
out['R3_weighted_bmo_family'] = r3

# ---- [R4] ----
f = lambda r: (r+1)**4/(r**4+1)
from mpmath import findroot, diff
# stationary points of log f: 4/(r+1) - 4r^3/(r^4+1) = 0  <=>  (r^4+1) = r^3 (r+1) <=> 1 = r^3
# so r = 1 is the unique positive critical point; f(1) = 16/2 = 8, f(0) = 1, f(inf) = 1
crit = findroot(lambda r: 4/(r+1) - 4*r**3/(r**4+1), mpf('1.3'))
out['R4_weight_equivalence'] = dict(unique_positive_critical_point=str(crit),
                                    f_at_crit=str(f(crit)),
                                    f_at_1=str(f(mpf(1))), f_at_0=str(f(mpf(0))),
                                    f_at_1e6=str(f(mpf(10)**6)),
                                    algebra_note='log f stationary  <=>  4/(r+1)=4r^3/(r^4+1) <=> r^4+1=r^4+r^3 <=> r^3=1 <=> r=1')

# ---- [R5] a control that can fail: does ||h_R||_BMO * log R stay bounded? ----
# lower bound: mean oscillation of h_R over balls B(0,rho), rho = 10^j, j = 0..k
r5 = []
for k in (2, 4, 8, 16, 32, 64):
    R = mpf(10)**k; lg = log(R)
    h = lambda r: max(mpf(0), 1 - log(1+r)/lg)
    worst = mpf(0); arg = None
    for j in [0] + [k*m//8 for m in range(1,9)] + [k+1]:
        rho = mpf(10)**j
        m = 3*quad(lambda r: h(r)*r*r, [0, min(rho,mpf(1)), rho])/rho**3
        osc = 3*quad(lambda r: abs(h(r)-m)*r*r, [0, min(rho,mpf(1)), rho])/rho**3
        if osc > worst: worst, arg = osc, rho
    r5.append(dict(R=f'1e{k}', lower_bound_BMO=str(worst), times_logR=str(worst*lg),
                   argmax_rho=str(arg)))
out['R5_bmo_lower_bound'] = r5
print(json.dumps(out, indent=1))
