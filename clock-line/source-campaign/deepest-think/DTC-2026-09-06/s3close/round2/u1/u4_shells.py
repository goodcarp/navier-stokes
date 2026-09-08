"""
u4 -- the closed integro-differential system on a discretised shell, and the comparison.

The lower bound proved in PROOF.md sec.4 is, for every label shell sigma = log(rho/rho_0)/L,

  d/dtheta log lambda(sigma,theta)  >=  (1 - epsT) * (1/2) INT_{sigma+d}^{1} P_h(lambda(s',theta)) ds'
                                        -  C_visc/L ,        d := ell_shell/L ,

with lambda(.,0) = 1.  P_h is increasing on [1, lam_max] (u2), so the right-hand side is
monotone in the profile and the discretised system is QUASIMONOTONE:  d/dtheta x_i = f_i(x)
with  d f_i / d x_j >= 0  for every j != i.  Hence (comparison theorem for quasimonotone
systems; Walter, "Differential and Integral Inequalities", Thm 10.XII) any sub-solution stays
below any super-solution.  The two-line proof for this system is in PROOF.md sec.5.1.

Three comparison objects at the tracked shell sigma = 0:
  (A) the discretised system above, N = 200 log-shells;
  (B) the same with P_h(lambda) replaced by r_h * lambda  (a genuine lower bound: P_h >= r_h lam);
  (C) the CONSERVATIVE exponential exp( kappa_delta (1-epsT) (1-d) theta )  -- what the theorem
      of PROOF.md actually claims;
  (D) prove-lagrangian sec.4's closed form (1 - (1-sigma) kappa_c theta/2)^{-2} at the reduced,
      shifted rate kappa_c = 2 kappa_delta (1-epsT)(1-d)/2 -- reported for comparison only.
"""
import json, math
import numpy as np
from scipy.integrate import solve_ivp, quad

U2 = json.load(open("u2_results.json"))
DEG = math.pi/180.0
DELTA = 7.5*DEG
KAPPA = U2["kappa_delta_7.5"]
CSTAR = U2["c_star"]
LAM_MAX = U2["lam_max_apriori"]
R_H = U2["r_h_1_to_lam_max"]

# --- P_h on a table, interpolated (the ODE needs it thousands of times) --------------
def P_h_exact(lam, delta=DELTA):
    A = lam**2 - lam**-4; B = lam**-4
    f = lambda v: min(1.0, math.asin(min(1.0, v))/delta)*v*v*(A*v*v + B)**-2.5
    w = math.sin(delta)
    return 3.0*(quad(f, 0.0, w, epsabs=1e-13, epsrel=1e-13, limit=300)[0]
                + quad(f, w, 1.0, epsabs=1e-13, epsrel=1e-13, limit=300)[0])

LAM_TAB = np.linspace(1.0, 2.05, 1051)   # P_h increases on [1, 2.075]; the window caps lam at exp(3c/4)=1.838
P_TAB = np.array([P_h_exact(l) for l in LAM_TAB])
assert (np.diff(P_TAB) > 0).all(), "P_h must be increasing on the table range"
def P_h(lams):
    return np.interp(np.clip(lams, 1.0, 2.05), LAM_TAB, P_TAB)

OUT = {"P_h_table_increasing_on_[1,2.05]": bool((np.diff(P_TAB) > 0).all()),
       "P_h_argmax_above": 2.075,
       "P_h_interp_max_abs_err_vs_exact":
           float(max(abs(float(P_h(np.array([l]))[0]) - P_h_exact(l))
                     for l in [1.0, 1.13, 1.5, 1.8377, 2.0, 2.05]))}

N = 200
edges = np.linspace(0.0, 1.0, N+1)
sig = 0.5*(edges[:-1] + edges[1:])
dsig = 1.0/N

def run_system(L, d, epsT, C_visc, kernel, theta_max):
    """kernel(lams) -> per-shell rate density; d = shift in sigma units."""
    def rhs(th, y):
        lam = np.exp(y)
        k = kernel(lam)
        # cumulative tail integral INT_{sigma_i + d}^{1} k(sigma') dsigma'
        tail = np.concatenate([np.cumsum(k[::-1])[::-1]*dsig, [0.0]])   # tail[i] = INT_{edges[i]}^1
        val = np.interp(sig + d, edges, tail)
        return (1.0 - epsT)*val - C_visc/L
    sol = solve_ivp(rhs, (0.0, theta_max), np.zeros(N), rtol=1e-9, atol=1e-12,
                    max_step=theta_max/400.0, dense_output=True)
    return sol

def theta_to_lambda(sol, target, theta_max):
    """first theta with lambda(sigma=0) >= target, by bisection on the dense output."""
    f = lambda th: math.exp(sol.sol(th)[0]) - target
    if f(theta_max) < 0:
        return None
    lo, hi = 0.0, theta_max
    for _ in range(80):
        mid = 0.5*(lo + hi)
        if f(mid) >= 0: hi = mid
        else: lo = mid
    return hi

RES = {}
for L in [40.0, 160.0, 640.0]:
    ell = U2["ell_shell_mu0"]
    d = ell/L
    epsT = 0.0            # the comparison is run at eps_T' = 0; the budget applies (1-epsT)
    Cv = 0.0
    theta_max = 1.1
    solA = run_system(L, d, epsT, Cv, lambda lam: 0.5*P_h(lam), theta_max)
    solB = run_system(L, d, epsT, Cv, lambda lam: 0.5*R_H*lam, theta_max)
    thA = theta_to_lambda(solA, 1.5, theta_max)
    thB = theta_to_lambda(solB, 1.5, theta_max)
    kap_c = KAPPA*(1.0 - d - sig[0])   # the tracked cell sits at sigma = sig[0], not 0
    thC = math.log(1.5)/kap_c                                   # conservative exponential
    thD = 2.0*(1.0 - math.sqrt(2.0/3.0))/(kap_c)                # closed form (1-kap_c th/2)^-2 = 3/2
    ths = np.linspace(0.0, min(theta_max, 1.2*thC), 200)
    lamA = np.array([math.exp(solA.sol(t)[0]) for t in ths])
    lamB = np.array([math.exp(solB.sol(t)[0]) for t in ths])
    lamC = np.exp(kap_c*ths)
    lamD = (1.0 - kap_c*ths/2.0)**-2
    RES["L=%g" % L] = {
        "d_shift": d, "kappa_c": kap_c, "sigma_of_tracked_cell": float(sig[0]),
        "kappa_c_continuum": KAPPA*(1.0 - d),
        "theta_A_lambda=1.5": thA, "theta_B_lambda=1.5": thB,
        "theta_C_conservative": thC, "theta_D_closed_form": thD,
        "A_ge_B_everywhere": bool((lamA >= lamB - 1e-12).all()),
        "A_ge_C_everywhere": bool((lamA >= lamC - 1e-12).all()),
        "B_ge_C_everywhere": bool((lamB >= lamC - 1e-12).all()),
        "A_ge_D_everywhere": bool((lamA >= lamD - 1e-12).all()),
        "min(A-C)": float((lamA - lamC).min()), "min(A-D)": float((lamA - lamD).min()),
        "lambda_A_at_theta_C": float(math.exp(solA.sol(thC)[0])),
        "lambda_max_reached_A": float(max(math.exp(solA.sol(t)[j])
                                          for t in ths for j in [0])),
        "profile_decreasing_in_sigma": bool((np.diff(solA.sol(thC)) <= 1e-12).all()),
        "acceleration_ratio_C_over_A": thC/thA if thA else None,
    }
OUT["comparison"] = RES

# --- the a priori bound lambda <= exp(3c/4): check it is respected -------------------
solchk = run_system(640.0, U2["ell_shell_mu0"]/640.0, 0.0, 0.0,
                    lambda lam: 0.5*P_h(lam), CSTAR)
OUT["lambda_at_c_star_sigma0_L640"] = float(math.exp(solchk.sol(CSTAR)[0]))
OUT["lambda_apriori_cap"] = LAM_MAX
OUT["apriori_cap_respected"] = bool(math.exp(solchk.sol(CSTAR)[0]) <= LAM_MAX + 1e-9)

# --- conservative clock: theta_* and t_* --------------------------------------------
OUT["theta_star_conservative_noloss"] = math.log(1.5)/KAPPA
OUT["c2_conservative"] = 2.0*math.log(1.5)/KAPPA        # M T_d <= c2/(2L) -> T(Re) <= 2*c2/logRe
OUT["c2_times_2"] = 4.0*math.log(1.5)/(2.0*KAPPA)
# the brief's t_* formula, evaluated, against the two defensible ones
I_brief = quad(lambda l: 1.0/(l*P_h_exact(l)), 1.0, 1.5, epsabs=1e-12, epsrel=1e-12)[0]
OUT["brief_t_star_theta"] = 2.0*I_brief
OUT["closed_form_theta_max"] = 4.0*(1.0 - math.sqrt(2.0/3.0))
OUT["conservative_theta_max"] = math.log(1.5)/KAPPA
OUT["brief_formula_vs_closed_form_ratio"] = 2.0*I_brief/(4.0*(1.0 - math.sqrt(2.0/3.0)))

with open("u4_results.json", "w") as fh:
    json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
print(json.dumps(OUT, indent=1, sort_keys=True, default=str))
