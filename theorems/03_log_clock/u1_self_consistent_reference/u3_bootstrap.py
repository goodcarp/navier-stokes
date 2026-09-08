"""
u3 -- the bootstrap with the SLAVED reference.  The l=1 rate error is gone from the driver,
so eps_T' no longer appears on the right-hand side of the majorant system.

Two forms of the same estimate:

(G) GROENWALL form (the form the brief asks for).  m := sup_x |Phi-Lambda|/|x| :
      d m/dtheta   = G m + lam_max * lam_om * (C_R + 2 log lam_max)/L ,   m(0)=0
      G            = (3/2)(1 + C_1/L) + C''/L                (hypothesis (Gamma-off))
      mu           = lam_max^2 m
    linearised exponent  p_G c = G c .

(C) CONJUGATED form.  Phi_s(x) = T_{lam(|x|,s)}(x + e(s)),  ehat := sup |e|/|x| :
      d ehat/dtheta = (lam_om/L) [ 2 log lam_max + log(1/(1-ehat)) + lam_max^3 C_R ] (1+ehat)
      mu            = lam_max^3 ehat
    linearised exponent  p_C c = (lam_om/L)[1 + 2 log lam_max + lam_max^3 C_R] c  =  O(c/L) .
    (Gamma-off) is NOT used here at all.

Jacobian (both forms, no Corollary-3 composition needed because lambda^2 is what Lemma T'
(H3) compares against and  d/dtheta log(J_Phi/lambda^2) = (2/(ML))[a(Phi) - frak_a(|x|)] ):
      d w/dtheta = (2 lam_om/L) [ C_a + (1/2) Lrad(m) ] ,  mu_J = exp(w) - 1 .

For reference the ASSEMBLY exponent is recomputed here from its own formula.
"""
import json, math
import numpy as np
from scipy.integrate import solve_ivp

U2 = json.load(open("u2_results.json"))
REC = U2["RECORD"]
KAPPA = U2["kappa_delta_7.5"]
CSTAR = U2["c_star"]
LAM_MAX = U2["lam_max_apriori"]
LAM_OM = 1.5                       # M_s/M <= 3/2 : the contradiction hypothesis
R_H_SLAVED = U2["r_h_1_to_lam_max"]
R_H_ASSEMBLY = U2["r_h_1_to_1.5"]
C_A = U2["C_a_phi0_30deg"]

def eps_Tprime(mu, muJ, r_h, sens=1.0):
    if mu >= 0.999 or not np.isfinite(mu) or not np.isfinite(muJ):
        return float("inf")
    return sens*(2.0*math.pi/r_h)*(3.0*mu*(1.0+muJ)/(2.0*(1.0-mu)**5) + 3.0*muJ/8.0)

def Lrad_from_m(m, lam_mx):
    """|log(|Phi|/|x|)| with |Phi| in [lam^-2|x| - m|x|, lam|x| + m|x|]."""
    lo = lam_mx**-2 - m
    if lo <= 1e-9:
        return float("inf")
    return max(math.log(lam_mx + m), math.log(1.0/lo))

def bootstrap_groenwall(L, c, Cpp, C_R, C_a=C_A, C1=1.0, lam_mx=LAM_MAX, lam_om=LAM_OM):
    G = 1.5*(1.0 + C1/L) + Cpp/L
    drive = lam_mx*lam_om*(C_R + 2.0*math.log(lam_mx))/L

    def rhs(th, y):
        m, w = y
        Lr = Lrad_from_m(m, lam_mx)
        if not np.isfinite(Lr):
            return [0.0, 0.0]
        return [G*m + drive, (2.0*lam_om/L)*(C_a + 0.5*Lr)]

    def ev(th, y):
        return y[0] - (lam_mx**-2 - 1e-6)
    ev.terminal = True; ev.direction = 1
    sol = solve_ivp(rhs, (0.0, c), [0.0, 0.0], rtol=1e-10, atol=1e-14,
                    events=ev, max_step=c/200.0)
    m = float(sol.y[0, -1]); w = float(sol.y[1, -1])
    blew = (sol.t[-1] < c - 1e-12) or (not sol.success)
    mu = lam_mx**2*m if not blew else float("inf")
    muJ = math.exp(w) - 1.0 if (not blew and w < 50) else float("inf")
    return {"m": m, "mu": mu, "muJ": muJ, "w": w, "p_c": G*c, "G": G,
            "drive": drive, "closed": (not blew) and mu < 1.0}

def bootstrap_conjugated(L, c, C_R, C_a=C_A, lam_mx=LAM_MAX, lam_om=LAM_OM):
    K = lam_mx**3*C_R

    def rhs(th, y):
        eh, w = y
        if eh >= 0.999:
            return [0.0, 0.0]
        Lr = 2.0*math.log(lam_mx) + math.log(1.0/(1.0 - eh))
        de = (lam_om/L)*(Lr + K)*(1.0 + eh)
        # for the Jacobian, |Phi|/|x| in [lam^-2(1-eh), lam(1+eh)]
        Lr2 = max(math.log(lam_mx*(1.0 + eh)), 2.0*math.log(lam_mx) + math.log(1.0/(1.0 - eh)))
        return [de, (2.0*lam_om/L)*(C_a + 0.5*Lr2)]

    def ev(th, y):
        return y[0] - 0.99
    ev.terminal = True; ev.direction = 1
    sol = solve_ivp(rhs, (0.0, c), [0.0, 0.0], rtol=1e-10, atol=1e-14,
                    events=ev, max_step=c/200.0)
    eh = float(sol.y[0, -1]); w = float(sol.y[1, -1])
    blew = (sol.t[-1] < c - 1e-12) or (not sol.success)
    mu = lam_mx**3*eh if not blew else float("inf")
    muJ = math.exp(w) - 1.0 if (not blew and w < 50) else float("inf")
    p_c = (LAM_OM/L)*(1.0 + 2.0*math.log(lam_mx) + K)*c
    return {"ehat": eh, "mu": mu, "muJ": muJ, "w": w, "p_c": p_c,
            "closed": (not blew) and mu < 1.0}

def assembly_exponent(L, c, Cprime, r_h=R_H_ASSEMBLY, sens=1.0):
    """ASSEMBLY sec.3.2, recomputed from its own formula."""
    p = 1.5 + Cprime/L + 1.5*2.0*KAPPA*sens*(15.0*math.pi/4.0)/(r_h*(4.0/9.0))
    return p*c

OUT = {"inputs": {"kappa": KAPPA, "c_star": CSTAR, "lam_max": LAM_MAX,
                  "lam_omega": LAM_OM, "r_h_slaved": R_H_SLAVED,
                  "r_h_assembly": R_H_ASSEMBLY, "C_a": C_A}}

CR = {"proved_Gc0": U2["C_R_proved_Gc0"],
      "proved_Gc_value": U2["C_R_proved_Gc_value"],
      "proved_Gc_5x": U2["C_R_proved_Gc_5x"],
      "measured": U2["C_R_measured_Gc0"],
      "assembly_convention": U2["C_R_assembly_convention"]}
OUT["C_R_variants"] = CR

# ---- the three feedback exponents at c = c_* ---------------------------------------
ex = {}
for L in [40.0, 160.0, 640.0, 2560.0, 1e18]:
    ex["L=%g" % L] = {
        "assembly_proved_C'=151.15": assembly_exponent(L, CSTAR, REC["Cpp_proved"]),
        "assembly_proved_C'=119.33": assembly_exponent(L, CSTAR, REC["Cpp_computedV"]),
        "assembly_measured_sens": assembly_exponent(L, CSTAR, REC["Cpp_measured"],
                                                    sens=REC["lemmaT_sens_proved"]/REC["lemmaT_sens_measured"]),
        "slaved_groenwall_C''=151.15": bootstrap_groenwall(L, CSTAR, REC["Cpp_proved"], CR["proved_Gc0"])["p_c"],
        "slaved_groenwall_C''=11.74": bootstrap_groenwall(L, CSTAR, REC["Cpp_measured"], CR["measured"])["p_c"],
        "slaved_conjugated_proved": bootstrap_conjugated(L, CSTAR, CR["proved_Gc0"])["p_c"],
        "slaved_conjugated_measured": bootstrap_conjugated(L, CSTAR, CR["measured"])["p_c"],
    }
OUT["feedback_exponents"] = ex

# ---- mu(c), mu_J(c), eps_T' over an L-grid -----------------------------------------
grid = {}
for name, (Cpp, cr, sens) in {
        "proved":            (REC["Cpp_proved"],   CR["proved_Gc0"], 1.0),
        "proved_Ca_measured": (REC["Cpp_proved"],  CR["measured"],   1.0),
        "measured":          (REC["Cpp_measured"], CR["measured"],
                              REC["lemmaT_sens_proved"]/REC["lemmaT_sens_measured"]),
}.items():
    for L in [40.0, 160.0, 640.0, 2560.0, 10240.0]:
        g = bootstrap_groenwall(L, CSTAR, Cpp, cr)
        cj = bootstrap_conjugated(L, CSTAR, cr)
        grid["%s_L=%g" % (name, L)] = {
            "groenwall_mu": g["mu"], "groenwall_muJ": g["muJ"],
            "groenwall_epsT": eps_Tprime(g["mu"], g["muJ"], R_H_SLAVED, sens),
            "groenwall_mu_times_L": g["mu"]*L,
            "conjugated_mu": cj["mu"], "conjugated_muJ": cj["muJ"],
            "conjugated_epsT": eps_Tprime(cj["mu"], cj["muJ"], R_H_SLAVED, sens),
            "conjugated_mu_times_L": cj["mu"]*L,
            "p_c_groenwall": g["p_c"], "p_c_conjugated": cj["p_c"]}
OUT["grid"] = grid

# ---- which form wins, and the mu(c) coefficient in units of 1/L ---------------------
OUT["mu_c_over_1_over_L"] = {}
for name, cr in [("proved", CR["proved_Gc0"]), ("proved_Ca_measured", CR["measured"]),
                 ("measured", CR["measured"])]:
    Cpp = REC["Cpp_proved"] if name != "measured" else REC["Cpp_measured"]
    big = 1e9
    g = bootstrap_groenwall(big, CSTAR, Cpp, cr)
    cj = bootstrap_conjugated(big, CSTAR, cr)
    OUT["mu_c_over_1_over_L"][name] = {"groenwall": g["mu"]*big, "conjugated": cj["mu"]*big,
                                       "groenwall_muJ_L": g["muJ"]*big,
                                       "conjugated_muJ_L": cj["muJ"]*big}

with open("u3_results.json", "w") as fh:
    json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
print("C_R variants:", json.dumps(CR, indent=1))
print("\nfeedback exponents p*c at c = c_* = %.7f:" % CSTAR)
print(json.dumps(ex["L=640"], indent=1))
print("\nmu(c) * L  (L -> infinity):")
print(json.dumps(OUT["mu_c_over_1_over_L"], indent=1))
print("\ngrid:")
for k in sorted(grid):
    r = grid[k]
    print("%-28s muG=%10.4g muJG=%10.4g epsTG=%10.4g | muC=%10.4g epsTC=%10.4g"
          % (k, r["groenwall_mu"], r["groenwall_muJ"], r["groenwall_epsT"],
             r["conjugated_mu"], r["conjugated_epsT"]))
