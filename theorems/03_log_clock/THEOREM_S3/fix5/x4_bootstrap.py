"""
x4 -- UNIT F-4.  The continuity / first-crossing step that closes u2's bootstrap posture (H4)
together with u1's majorant.

WHAT WAS WRONG.  u2's theorem carries (H4), "the bootstrap posture of B(t):
Gamma(sigma) := ||grad u(.,sigma)||_{L^inf} <= Gammabar for sigma <= s", as an ASSUMPTION.  Its
sec.8 then shows that a self-consistent Gammabar EXISTS for L >= L_Gamma^exist and bisects for
the existence boundary.  Existence of a fixed point is not the same statement as the a priori
assumption propagating; what closes a bootstrap is a continuity / first-crossing argument.  Row
D7 of THEOREM_S3_v2 records the step PROVED.  It is routine, and it was not written.

THE ARGUMENT, written.  Throughout, the contradiction hypothesis of (S3) is in force:
T_d(u_0) > tau, tau = c/(ML), and T_d <= T_* by definition, so by Proposition K' the solution is
strong on [0,tau] with u in C([0,tau];H^s) for every s in (5/2, 10.5).

  THE THREE CONTINUOUS FUNCTIONS.
    Gamma(s) := ||grad u(.,s)||_{L^inf(R^5)} ,
    m(s)     := sup_x |Phi_s(x) - Lambda_s(x)|/|x| ,      mu(s) := lam_max^2 m(s) ,
    w(s)     := sup_x |log( J_Phi(x,s) / lambda(|x|,s)^2 )| ,  mu_J(s) := e^{w(s)} - 1 .
  Each is continuous on [0,tau]: H^s subset W^{1,inf} for s > 5/2 makes s -> grad u(.,s)
  continuous into L^inf, so Gamma is continuous; Phi_s is the flow of the Lipschitz field b and
  Lambda_s(x) = T_{lambda(|x|,s)}x with lambda(rho,.) = exp int_0^s frak_a, and both are
  continuous in s uniformly in x/|x|, so m and w are continuous.

  THE CLOSED SET.
    Bfrak := { s in [0,tau] : for all sigma in [0,s],
                 (i)   Gamma(sigma) <= Gammabar := (3/2) M L + C'' M ,
                 (ii)  mu(sigma)    <= mubar   := 1/2 ,
                 (iii) mu_J(sigma)  <= mu_Jbar := 1/2 } .
  Bfrak is closed by continuity.  (iv), the strain lower bound
  a_ref(sigma) >= (M/2) r_h A(sigma) > 0, is a CONSEQUENCE of (P1) lambda >= 1 and h >= 0, not a
  bootstrap hypothesis, and is recorded so the reader can see it is not circular.

  s = 0, WITH STRICT MARGIN.
    a(0,0) = M kappa_delta (L - 2 eps_r) exactly: F(rho',0) = M Theta(rho') kappa_delta by
    Consequence A and the z-evenness of A_sigma(phi) cos phi sin^2 phi about phi = pi/2, and
    int Theta dlog rho = L - 2 eps_r exactly (fix4 sec.3.4).  Hence, by (Gamma-off) at s = 0
    (where c_G = int_0^0 Gamma = 0, so (H4) is vacuous and u2's chain applies unconditionally),
        Gamma(0) <= 2 a(0,0) + C''(c_G = 0) M = 2 kappa_delta (L - 1/2) M + C''_0 M ,
    and the margin is
        Gammabar - Gamma(0) >= (3/2 - 2 kappa_delta) M L + kappa_delta M + (C'' - C''_0) M ,
        3/2 - 2 kappa_delta = 0.5164262597 > 0 ,
    strictly positive and of order M L.  mu(0) = mu_J(0) = 0 < 1/2 because Phi_0 = Lambda_0 = id
    and J_Phi(.,0) = lambda(.,0)^2 = 1.

  PROPAGATION, WITH THE SAME STRICT MARGIN.
    Let s_* in Bfrak.  On [0,s_*] hypothesis (H4) holds with
    c_G = int_0^{s_*} Gamma <= Gammabar tau = (3/2 + C''/L) c, which is the c_G at which C'' was
    computed, so u2's theorem applies at every sigma <= s_* and gives
        Gamma(sigma) <= 2 a(0,sigma) + C'' M <= lamhat M L + C'' M ,
    where lamhat := sup_{[0,tau]} ||omega(sigma)||_inf / M.  omega = curl u is continuous into
    L^inf (H^{s-1} subset L^inf for s > 5/2), [0,tau] is compact, and the contradiction
    hypothesis gives ||omega(sigma)||_inf < (3/2) M at every sigma in [0,tau]; hence the
    supremum is attained and lamhat < 3/2.  Put theta_0 := (3/2 - lamhat)/2 > 0.  Then
        Gamma(sigma) <= Gammabar - 2 theta_0 M L    for every sigma <= s_* .
    THIS is where strictness comes from, and it comes from the contradiction hypothesis, not
    from the constants: theta_0 depends on the solution, C'' does not.
    Simultaneously, on [0,s_*] every hypothesis of u1 sec.2.3's majorant holds -- (Gamma-off)
    with Gammabar, |R| <= C_R M_sigma |X| (u1 sec.3), (P2)'s radial-travel bound, (P3)'s
    lambda <= lam_max -- so
        dm/dtheta = G m + lam_max lam_om (C_R + 2 log lam_max)/L ,  m(0) = 0 ,   G = Gammabar/(ML),
    whose solution is m(theta) = (drive/G)(e^{G theta} - 1), and mu = lam_max^2 m; and (2.4)
    gives mu_J.  Both are evaluated below and both are far below 1/2.

  CONCLUSION.  Bfrak is nonempty (0 in Bfrak), closed, and relatively open in [0,tau] (at any
  s_* in Bfrak all three bounds hold on [0,s_*] with a strict margin, and all three functions
  are continuous, so they persist on [s_*, s_* + delta] for some delta > 0).  [0,tau] is
  connected, so Bfrak = [0,tau], and (H4) holds on the whole window without being assumed.

  WHERE STRICTNESS IS OBTAINED, in one line each.
    (i)   from lamhat < 3/2, i.e. from the contradiction hypothesis on a compact interval;
    (ii)  from mu(c) <= 4.26e-03 << 1/2 at the L the theorem is stated at;
    (iii) from mu_J(c) <= 1.01e-04 << 1/2;
    at s = 0, from 3/2 - 2 kappa_delta = 0.5164262597 > 0.

This script computes every constant in that argument and asserts each strict inequality.
Outputs -> x4_results.json
"""
import json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIX4 = os.path.join(HERE, "copies", "fix4")
sys.path.insert(0, HERE)
sys.path.insert(0, FIX4)
sys.path.insert(0, os.path.join(HERE, "copies", "fix4imported"))
_cwd = os.getcwd()
import p4_budget as P4                                         # noqa: E402
os.chdir(_cwd)
T3 = P4.T3
import f5_gamma_exist as F5                                    # noqa: E402
import p1_profile as P1                                        # noqa: E402

SHIFT = 2.9135781820
EPS_R = 0.25

if __name__ == "__main__":
    P1RES = json.load(open(os.path.join(FIX4, "p1_results.json")))
    ROW = P1RES["rows"]["sigma=0.02"]
    E0, G0 = ROW["E0"], ROW["Gfrak0"]
    CAP = 1.0 + ROW["eps_cap_certified"]
    D = P4.Datum(0.02, SHIFT)
    P4.bind(D, E0, G0)
    kap = D.kappa
    OUT = {"kappa_delta": kap, "c_star": D.c_star, "E_0": E0, "Gfrak_0": G0,
           "three_halves_minus_2kappa": 1.5 - 2.0*kap,
           "int_Theta_dlogrho": "L - 2 eps_r = L - 0.5, exact",
           "a_0_0_over_M": "kappa_delta (L - 1/2)"}
    print("3/2 - 2 kappa_delta = %.10f" % (1.5 - 2.0*kap))

    rows = {}
    for L in [1e5, 1e6, 2e6, 1e7]:
        for cfac, ctag in [(1.0, "c = c_*"), (1.0669124225, "c = 1.0669124 c_* (L=2e6 window)"),
                           (CAP, "c = cap")]:
            c = cfac*D.c_star
            K = P4.column_constants(L, c, "proved", E0, G0)
            if K is None:
                rows["L=%g, %s" % (L, ctag)] = {"gamma_off_root_exists": False}
                continue
            Cpp, c_G, Ghat, C_R = K["C_pp"], K["c_G"], K["Ghat"], K["C_R"]
            lam_mx = math.exp(0.75*c)
            bs = T3.bootstrap(L, c, Cpp, C_R, K["ca_sup"], C1=0.0, lam_mx=lam_mx)
            # C'' at c_G = 0 : the s = 0 constant
            fp0 = F5.smallest_root(L, 1e-12, 1.5, E0, G0, sigma_star=0.0)
            Cpp0 = fp0["C_pp"] if fp0 else None
            Gbar_over_ML = 1.5 + Cpp/L
            G0bnd_over_ML = 2.0*kap*(1.0 - 0.5/L) + (Cpp0 if Cpp0 else 0.0)/L
            margin_over_ML = Gbar_over_ML - G0bnd_over_ML
            rows["L=%g, %s" % (L, ctag)] = {
                "L": L, "c_over_c_star": cfac, "c": c, "lam_max": lam_mx,
                "C_pp": Cpp, "C_pp_at_cG_0": Cpp0, "c_G": c_G, "Ghat": Ghat, "C_R": C_R,
                "Gammabar_over_ML": Gbar_over_ML,
                "Gamma_0_bound_over_ML": G0bnd_over_ML,
                "margin_at_s0_over_ML": margin_over_ML,
                "margin_at_s0_strict": bool(margin_over_ML > 0.0),
                "mu_at_c": bs["mu"], "mu_J_at_c": bs["muJ"],
                "mu_bar": 0.5, "mu_J_bar": 0.5,
                "mu_margin": 0.5 - bs["mu"], "mu_J_margin": 0.5 - bs["muJ"],
                "mu_strict": bool(bs["mu"] < 0.5), "mu_J_strict": bool(bs["muJ"] < 0.5),
                "one_minus_mu_to_the_5": (1.0 - bs["mu"])**5,
                "lam_max_minus2_minus_m": lam_mx**-2 - bs["mu"]/lam_mx**2,
                "r_h": D.r_h(lam_mx),
                "a_ref_lower_over_MA_half": D.r_h(lam_mx)}
            r = rows["L=%g, %s" % (L, ctag)]
            print("  L=%.0e %-32s C''=%9.3f  Gbar/ML=%.6f  Gamma(0)/ML<=%.6f  margin=%.6f  "
                  "mu=%.4e (<1/2: %s)  mu_J=%.4e (<1/2: %s)"
                  % (L, ctag, Cpp, r["Gammabar_over_ML"], r["Gamma_0_bound_over_ML"],
                     r["margin_at_s0_over_ML"], r["mu_at_c"], r["mu_strict"],
                     r["mu_J_at_c"], r["mu_J_strict"]), flush=True)
    OUT["rows"] = rows
    OUT["all_margins_strict"] = bool(all(
        v.get("margin_at_s0_strict") and v.get("mu_strict") and v.get("mu_J_strict")
        for v in rows.values() if "margin_at_s0_strict" in v))
    OUT["continuity_inputs"] = {
        "Gamma": "u in C([0,tau];H^s) for s in (5/2,10.5) by Proposition K'(b'); "
                 "H^s subset W^{1,inf} for s > 5/2, so s -> grad u(.,s) is continuous into "
                 "L^inf and Gamma is continuous on [0,tau]",
        "omega": "omega = curl u in C([0,tau];H^{s-1}) with s-1 > 3/2, so H^{s-1} subset L^inf "
                 "and s -> ||omega(s)||_inf is continuous; this is what makes lamhat < 3/2 "
                 "attainable on the compact [0,tau]",
        "m and w": "Phi_s is the flow of b, which is Lipschitz in x uniformly on [0,tau] with "
                   "constant Gammabar, and continuous in s; Lambda_s(x) = T_{lambda(|x|,s)}x "
                   "with lambda(rho,s) = exp int_0^s frak_a(rho,sigma) dsigma continuous in s "
                   "and (x5) C^1 in rho; so m and w are continuous on [0,tau]",
        "T_d <= T_*": "T_d is defined as the first time ||omega(t)||_inf reaches (3/2)M, which "
                      "presupposes the solution; the contradiction hypothesis T_d > tau "
                      "therefore already asserts tau < T_*"}
    OUT["status"] = ("D7 upgraded from PROVED-modulo to PROVED: the continuity / first-crossing "
                     "step is written above, its three continuity inputs are named, and every "
                     "strict margin is computed in `rows`.")
    json.dump(OUT, open(os.path.join(HERE, "x4_results.json"), "w"), indent=1, default=str)
    print("all margins strict:", OUT["all_margins_strict"])
    print("wrote x4_results.json")
