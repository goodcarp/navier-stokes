"""
r5 -- three smaller checks.

 (1) The feedback exponent.  ASSEMBLY section 3.2 quotes
        p/(ML) = 3/2 + C'/L + (3/2)(2 kappa)(15 pi/4)(9/4)/r_h = 41.974 (L -> inf)
     with 15 pi/4 = 2 pi (3/2 + 3/8), i.e. it sets mu_J = mu in Lemma T''s bracket.  The
     bootstrap actually integrated in a2 is a 2-state linear-plus-nonlinear system in
     (m, mu_J^Lambda); its linearised growth rate is the dominant eigenvalue of the Jacobian
     of a2's rhs at the origin.  Computed here by finite differences of a2's own rhs (from
     the copy) and compared with the quoted p.

 (2) The L-scaling column of the term table (section 3.1): evaluate a2's viscous_budget at
     L and 4L (fixed c, f, C_K) and report the fitted exponents.

 (3) BFG bookkeeping.  Section 4 says 'both sides are invariant under u -> lambda u(lambda x,
     lambda^2 t) at fixed nu' as the reason nu = 1 can be assumed.  That scaling does not change
     nu; the map that sets nu = 1 is  u~(x,t) = u(x, t/nu)/nu, checked here symbolically
     (1D template of the NS operator, which is what the scaling acts on), and M0 T_d is invariant
     under it.  Also the strict/non-strict subtlety in 'T_d >= T'.
"""
import json, math, os, time
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
COPIES = os.path.join(HERE, "copies")
os.chdir(COPIES)
while not os.path.exists(os.path.join(COPIES, "a1_results.json")):
    time.sleep(5)
src = open("a2_budget.py").read().split('if __name__ == "__main__":')[0]
G = {"__name__": "r5"}
exec(src, G)
KAPPA = G["KAPPA"]; R_H = G["R_H"]; KAPPA_S = G["KAPPA_S"]; RECORD = G["RECORD"]
eps_Tprime = G["eps_Tprime"]; SENS = G["SENS_MEASURED_FACTOR"]; C_a_proved = G["C_a_proved"]; PHI0 = G["PHI0"]
OUT = {}

# ---------------------------------------------------------------- (1) Jacobian
def rhs_factory(L, Cprime, C_a, sens, c, geom_lo=4.0/9.0, geom_hi=1.5, lamfac_grad=1.5):
    GoML = lamfac_grad + Cprime/L
    c_G = GoML*c
    def rhs(th, y):
        m, muJL = y
        mu = m/geom_lo
        muJ = muJL*(1.0+2.0*KAPPA_S/L) + 2.0*KAPPA_S/L
        epsT = eps_Tprime(mu, muJ, R_H, sens)
        lam = min(math.exp(KAPPA*th), 1.5)
        drive = 2.0*(KAPPA*epsT + KAPPA*c_G/L) + C_a*lam/L
        return np.array([GoML*m + geom_hi*drive,
                         2.0*(KAPPA*epsT + KAPPA*c_G/L + C_a*lam/L + GoML*mu)])
    return rhs

jac = {}
cstar = math.log(1.5)/KAPPA
for col, Cp, Ca, sens in [("proved", RECORD["Cprime_A1_computedV"], C_a_proved(math.sin(PHI0)), 1.0),
                          ("measured", RECORD["Cprime_measured"], RECORD["material_offset_M"], SENS)]:
    for L in [160.0, 1e6, 1e12]:
        f = rhs_factory(L, Cp, Ca, sens, cstar)
        h = 1e-7
        y0 = np.array([0.0, 0.0])
        J = np.zeros((2, 2))
        for j in range(2):
            e = np.zeros(2); e[j] = h
            J[:, j] = (f(0.0, y0+e) - f(0.0, y0))/h
        ev = np.linalg.eigvals(J)
        pq = 1.5 + Cp/L + 1.5*2.0*KAPPA*sens*(15.0*math.pi/4.0)/(R_H*(4.0/9.0))
        lam_max = float(max(ev.real))
        jac["%s_L=%g" % (col, L)] = {"Jacobian": J.tolist(), "eigenvalues": [float(x.real) for x in ev],
                                     "dominant_rate": lam_max, "assembly_p_over_ML": pq,
                                     "ratio_assembly_over_true": pq/lam_max,
                                     "pc_star_true": lam_max*cstar, "pc_star_assembly": pq*cstar,
                                     "exp_pc_true": math.exp(lam_max*cstar), "exp_pc_assembly": math.exp(pq*cstar)}
OUT["feedback_exponent"] = jac

# ---------------------------------------------------------------- (2) scaling exponents
vb = G["viscous_budget"]
sc = {}
for c in [0.25, cstar]:
    for (L1, L2) in [(160.0, 640.0), (640.0, 2560.0)]:
        a = vb(L1, c, 1.0, 1.0, Cprime=119.33); b = vb(L2, c, 1.0, 1.0, Cprime=119.33)
        row = {}
        for k in ["eps_bulk", "E_hess", "E_4", "q_half"]:
            row[k + "_exponent"] = math.log(b[k]/a[k])/math.log(L2/L1)
        row["E_tail_L1"] = a["E_tail"]; row["E_tail_L2"] = b["E_tail"]
        row["c_G_L1"] = a["c_G"]; row["c_G_L2"] = b["c_G"]
        sc["c=%g_L=%g->%g" % (c, L1, L2)] = row
OUT["scaling_exponents"] = sc
# E_hess at fixed L vs C_K, and vs c (e^{3 c_G} claim)
L = 640.0
r1 = vb(L, 0.5, 1.0, 1.0, Cprime=119.33); r2 = vb(L, 1.0, 1.0, 1.0, Cprime=119.33)
OUT["E_hess_c_dependence"] = {"E_hess(c=1)/E_hess(c=0.5)": r2["E_hess"]/r1["E_hess"],
                              "e^{3(cG(1)-cG(0.5))}": math.exp(3*(r2["c_G"]-r1["c_G"])),
                              "c-prefactor (tau*V1 ~ c^2 (e^{2cG}-1)/cG)": (1.0**2/0.5**2)*((math.exp(2*r2["c_G"])-1)/r2["c_G"])/((math.exp(2*r1["c_G"])-1)/r1["c_G"])*math.exp(r2["c_G"]-r1["c_G"])}
OUT["E_hess_CK_linear"] = vb(L, 0.5, 66.6622, 1.0, Cprime=119.33)["E_hess"]/r1["E_hess"]

# ---------------------------------------------------------------- (3) BFG scaling
x, t, nu = sp.symbols('x t nu', positive=True)
u = sp.Function('u')
ut = u(x, t/nu)/nu                                   # candidate solution of the nu = 1 equation
NS_nu = lambda w, n, tv: sp.diff(w, tv) + w*sp.diff(w, x) - n*sp.diff(w, x, 2)   # 1D template of the operator
lhs = NS_nu(ut, 1, t)
s = sp.symbols('s')
rhs = (NS_nu(u(x, s), nu, s)/nu**2).subs(s, t/nu)
OUT["bfg_viscosity_rescaling_residual"] = str(sp.simplify(lhs - rhs))
OUT["bfg_M0Td_invariance"] = "omega~ = omega/nu so M~ = M/nu ; T~ = nu T ; M~ T~ = M T"
OUT["bfg_fixed_nu_scaling"] = "u -> lambda u(lambda x, lambda^2 t): omega -> lambda^2 omega, T_d -> T_d/lambda^2, nu unchanged; it does not reach nu = 1"
OUT["bfg_strictness"] = ("Theorem 10 gives sup_{t<=T} ||omega(t)||_inf <= M ||omega_0||_inf; with T_d := first time "
                         ">= (3/2)M0, take any M in (1, 3/2): then T_d > T >= 1/(c(M) M0), so M0 T_d >= 1/c(M) with M fixed, "
                         "an absolute constant; the assembly's M = 3/2 reading needs strictness it does not state")

with open(os.path.join(HERE, "r5_results.json"), "w") as fh:
    json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
print(json.dumps(OUT, indent=1, default=str))
