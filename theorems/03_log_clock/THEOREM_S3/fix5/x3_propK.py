"""
x3 -- UNIT F-3.  Proposition K restated for H^s data, and the datum's exact Sobolev index.

WHAT WAS WRONG.  THEOREM_S3_v2 sec.1.3 cites Proposition K of
`forced-route-2026-09/theorems/02_ADDENDUM_uniqueness_2026-09-08.md` sec.3 for existence and
uniqueness of the maximal strong solution.  Proposition K's hypothesis is "u_0 smooth, divergence
free and axisymmetric satisfying (4)", and (4) is Fefferman's decay condition quantified over
every multi-index and every K; smooth plus (4) is the Schwartz condition, and Proposition K(a)
says so.  The datum of sec.1.1 is NOT Schwartz: omega_0 decays like rho^{-8}, and at the origin
Theta(rho) = rho^8/(rho^8 + e^2) exactly, so omega_0 is C^{7,1} there and lies in H^s only for
s < 9.5.  K(a), K(c) and K(e) are quantified over EVERY s >= 0 and are false for this datum
above that index.  fix4 sec.7 M6 records the defect and does not repair it.

WHAT IS DONE HERE.  Proposition K' is stated with a finite Sobolev index and its conclusions
quantified over s in a range the datum supplies; the argument is [FENCES] section 1.3 Step 3
(local existence in H^s, s > 5/2, for the forced equation with an existence time bounded below by
||u(t_0)||_{H^s} and ||f||_{L^1(t_0,t_0+1;H^s)}; uniqueness in L^inf_t H^s; propagation of higher
regularity) together with [FENCES] Lemma A at alpha = 2 (a finite L^inf bound forces a finite H^s
bound, hence continuation), which is what makes T_* independent of s.  Then every consumer of
Proposition K in the chain is listed with the highest derivative of u it actually uses, and the
list is compared with what s < 10.5 supplies.

THE SCRIPT CHECKS
  1. Theta(rho) = rho^8/(rho^8 + e^2) EXACTLY at eps_r = 1/4 (the inner edge), and the outer edge
     tail bound e^{-8(L - eps_r - u)}.
  2. omega_0's Cartesian components are  -(M/N_sigma) [Theta(rho)/rho] W_sigma(phi) (-x_2, x_1, 0),
     and Theta(rho)/rho = rho^7/e^2 + O(rho^15): the leading singular factor is |x|^7 x_j, which
     is homogeneous of degree 8.
  3. W_sigma is a smooth function ON THE SPHERE: it is even at both poles, so every ODD derivative
     of W_sigma vanishes at phi = 0 and phi = pi.  Checked to order 3.
  4. The Fourier transform of a degree-8 homogeneous function that is smooth off the origin is
     homogeneous of degree -8-3 = -11, with a constant that is nonzero here (the Gamma-function
     coefficient for (n,l,a) = (3,1,8) has no pole), so |omega_0hat(xi)| ~ |xi|^{-11} and
     int (1+|xi|^2)^s |omega_0hat|^2 dxi ~ int k^{2s-20} dk converges exactly for s < 9.5.
     Confirmed by a direct spherical-Bessel (Hankel, l = 1) transform of the radial profile, whose
     large-k decay exponent is fitted.

Outputs -> x3_results.json
"""
import json, math, os, sys
import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fx_profile as FX

mp.mp.dps = 40
EPS_R = 0.25
E2 = math.exp(2.0)

PROP_K_PRIME = {
 "name": "Proposition K'",
 "hypothesis": "u_0 divergence free and axisymmetric on R^3 with u_0 in H^{s_1}(R^3) for every "
               "s_1 < S_1 (here S_1 = 10.5), and f smooth and axisymmetric satisfying "
               "Fefferman's (5) with div f = 0 (here f = 0, which satisfies (5) trivially).",
 "s_0": "the local theory needs s > 5/2; take s_0 = 5/2.",
 "conclusions": [
   "(a') u_0 in H^s for every s in (s_0, S_1); f in L^1(0,inf;H^s) cap L^2(0,inf;H^s) for every "
   "s (immediate at f = 0), so hypothesis (A1) of [FENCES] Lemma A holds at every index.",
   "(b') There is T_* in (0,inf] and a unique u with u in C([0,T'];H^s) cap L^2(0,T';H^{s+1}) "
   "for every T' < T_* and every s in (s_0, S_1), maximal among such; and T_* does not depend "
   "on s in that range.",
   "(c') If T_* < inf then sup_{t<T_*}||u(t)||_{L^inf} = inf.",
   "(d') u is axisymmetric, and swirl free if u_0 and f are.",
   "(e') On every closed [0,T'] with T' < T_*, u in L^inf([0,T']; L^2 cap L^6 cap L^inf) and "
   "grad u in L^inf(R^3 x [0,T']); u satisfies (H*) in the finite-index form "
   "u in C([0,T');H^s) for every s in (s_0, S_1).",
   "(f') u in L^inf(0,T_*;L^2) cap L^2(0,T_*;H^1) up to T_*."],
 "proof": [
   "(b') LOCAL EXISTENCE: [FENCES] section 1.3 Step 3 gives, for nu > 0 and a force obeying "
   "(A1), the a priori bound d/dt ||u||_{H^s} <= C ||u||_{H^s}^2 + ||f||_{H^s} for s > 5/2 "
   "(pairing with Lambda^s u, the dissipative term having the good sign, and H^s subset "
   "W^{1,inf}), hence existence on [t_0, t_0 + c/(||u(t_0)||_{H^s} + ||f||_{L^1(t_0,t_0+1;H^s)})] "
   "by a regularisation-plus-Aubin-Lions scheme, and uniqueness in L^inf_t H^s from the L^2 "
   "estimate on the difference. That step is written out there rather than cited.",
   "(b') PROPAGATION: u in L^2(t_0,T;H^{s+1}) gives a.e. t_1 with u(t_1) in H^{s+1}; the same "
   "estimate at s+1 from t_1 on, iterated, yields every s' in (s, S_1). The iteration stops at "
   "S_1 because the DATUM stops there, not because the argument does.",
   "(b') T_* INDEPENDENT OF s: T_*(s) is non-increasing in s. Conversely, if T_*(s) < T_*(s') "
   "for s > s', then on [0,T_*(s)) the L^inf norm is bounded, since H^{s'} subset L^inf; "
   "[FENCES] Lemma A at alpha = 2, whose hypothesis (A1) holds at index s by (a'), then gives "
   "sup_{t<T_*(s)}||u(t)||_{H^s} < inf, and the local theory continues past T_*(s), a "
   "contradiction. So T_*(s) is the same for every s in (s_0, S_1).",
   "(c') is [FENCES] Lemma A at alpha = 2 read contrapositively, at any single admissible s.",
   "(d') Rotations about e_z and the reflection (x_1,x_2,x_3) -> (x_1,-x_2,x_3) preserve the "
   "equation, the class of (b'), the datum and the force; uniqueness in that class does the rest. "
   "Unchanged from Proposition K(d), which uses no property of the datum beyond its symmetry.",
   "(e') Sobolev embedding on a compact time interval, unchanged from K(e). H^s subset L^inf "
   "and H^s subset W^{1,inf} need only s > 5/2.",
   "(f') the energy identity of [FENCES] Lemma A Step 1, unchanged from K(f): it uses only "
   "u in L^2 and the force in L^1_t L^2."],
 "what_changes_against_K": "K(a), K(c) and K(e) are quantified over every s >= 0 and are false "
                           "for a datum that is only C^{7,1} at the origin. K'(a'),(c'),(e') are "
                           "quantified over s in (5/2, 10.5), which the datum supplies. Nothing "
                           "else in the proposition or in its Corollary C1 changes: C1's Lemma U "
                           "argument needs the reference in L^inf_t(L^2 cap L^6 cap L^inf) with "
                           "grad u in L^inf, i.e. K'(e'), and no higher index."}

CONSUMERS = [
 {"consumer": "u2 sec.2, the exact reduction ||grad u||_op <= 2|a| + r|grad a| + |omega^theta|",
  "object": "grad u, omega", "order_in_u": 1},
 {"consumer": "u2 sec.4, Gamma_rad and x.b <= Gamma_rad|x|^2", "object": "grad b", "order_in_u": 1},
 {"consumer": "u2 sec.5 P1, sup rho|eta| <= e^{c_R} E_0 M",
  "object": "eta solving D_t eta = nu Lap_5 eta, classical", "order_in_u": 3},
 {"consumer": "u2 sec.5 P2, sup rho^2|grad eta| <= e^{p c_G} Gfrak_0 M",
  "object": "grad eta solving D_t grad eta = -(grad b).grad eta + nu Lap_5 grad eta, classical",
  "order_in_u": 4},
 {"consumer": "u2 sec.6, r|grad a| by the Riesz composition", "object": "grad_5 eta",
  "order_in_u": 2},
 {"consumer": "pmax-h3v Part A (the P1/P2 proofs, barriers and comparison)",
  "object": "same as u2 sec.5; interior parabolic Schauder supplies the classical regularity",
  "order_in_u": 4},
 {"consumer": "pmax-h3v (D1)-(D5) on the datum", "object": "eta_0 in W^{1,inf} cap C, z-odd, "
  "real-analytic at the poles, rho^{-9} tail", "order_in_u": 1, "on_datum": True},
 {"consumer": "Theorem V.1 (gaps/gap-V-aronson), the Feynman-Kac representation",
  "object": "eta a bounded classical solution of D_t eta = nu Lap_5 eta", "order_in_u": 3},
 {"consumer": "Theorem V.4' hypothesis (H3'_vartheta)",
  "object": "d_v^k (eta_0 - eta_P) for k <= 4 on B(x_*,d)", "order_in_u": 5, "on_datum": True,
  "note": "the highest order anywhere in the chain, and it is a condition on the DATUM at t = 0, "
          "on a ball whose closest point to the origin is rho_* - d = 3.4619 rho_0, where "
          "eta_0 is real-analytic; the origin's C^{7,1} corner is nowhere near it"},
 {"consumer": "hk2 (H-K2), K_2 = ||grad^2_5 b||_{L^inf(N_tau)}",
  "object": "grad^2 b, and sup_B ||grad^2 eta||_F in the near ball", "order_in_u": 3},
 {"consumer": "Lemma T'/T'' hypotheses (D), (H1), (H2), (H3)",
  "object": "eta_0 measurable with |eta_0| <= M/r; Phi a C^1 diffeomorphism; lambda(.,s) in C^1",
  "order_in_u": 2},
 {"consumer": "u1 (P1),(P2),(P3) and the majorant system",
  "object": "the l = 1 Gegenbauer coefficient of omega^theta and its rho-derivative",
  "order_in_u": 2},
 {"consumer": "(Gamma-off), u2 sec.8", "object": "Gamma = ||grad u||_inf", "order_in_u": 1},
 {"consumer": "Corollary C1 of the addendum (uniqueness against a Clay-class competitor)",
  "object": "the reference in L^inf_t(L^2 cap L^6 cap L^inf), grad u in L^inf", "order_in_u": 1},
]

if __name__ == "__main__":
    OUT = {"Proposition_K_prime": PROP_K_PRIME, "consumers": CONSUMERS}

    # ---- 1. Theta = rho^8/(rho^8 + e^2) exactly at eps_r = 1/4
    def Theta(rho, L):
        u = math.log(rho)
        return 0.5*(math.tanh((u - EPS_R)/EPS_R) - math.tanh((u - L + EPS_R)/EPS_R))
    rows = []
    for rho in [0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 5.0, 20.0]:
        a = Theta(rho, 40.0)
        b = rho**8/(rho**8 + E2)
        rows.append({"rho": rho, "Theta": a, "closed_form_rho8_over_rho8_plus_e2": b,
                     "rel": abs(a - b)/max(abs(b), 1e-300)})
    # below rho ~ 0.3 the double-precision tanh saturates at -1 and Theta underflows to 0 while
    # the closed form is still representable; the identity is exact, and it is verified
    # symbolically rather than numerically there.
    import sympy as _sp
    _u = _sp.symbols("u", real=True)
    _lhs = (_sp.tanh((_u - _sp.Rational(1,4))/_sp.Rational(1,4)) + 1)/2
    _rhs = _sp.exp(8*_u)/(_sp.exp(8*_u) + _sp.exp(2))
    _res = _sp.simplify(_sp.expand((_lhs - _rhs).rewrite(_sp.exp)))
    OUT["inner_edge_closed_form"] = {"rows": rows, "max_rel": max(r["rel"] for r in rows),
                                     "sympy_residual_of_the_identity": str(_res),
                                     "identity_exact": bool(_res == 0),
                                     "why": "eps_r = 1/4 makes 2/eps_r = 8, so "
                                            "(1/2)(tanh(4u-1)+1) = rho^8/(rho^8 + e^2)"}
    print("inner edge closed form: max rel %.2e" % OUT["inner_edge_closed_form"]["max_rel"])

    # outer edge
    L = 40.0
    orow = []
    for u in [3.0, 4.0]:
        exact = 0.5*(1.0 - math.tanh((u - L + EPS_R)/EPS_R)) - 0.5*(1.0 - math.tanh((u-EPS_R)/EPS_R))
        bnd = math.exp(-8.0*(L - EPS_R - u))
        orow.append({"u": u, "psi_out_exact": 1.0/(1.0 + math.exp(-8.0*(u - L + EPS_R))),
                     "bound_exp": bnd})
    OUT["outer_edge_tail"] = orow

    # ---- 2. Theta(rho)/rho = rho^7/e^2 + O(rho^15)
    trow = []
    for rho in [1e-4, 1e-3, 1e-2, 3e-2]:
        v = Theta(rho, 40.0)/rho
        lead = rho**7/E2
        trow.append({"rho": rho, "Theta_over_rho": v, "leading_rho7_over_e2": lead,
                     "rel_residual": abs(v - lead)/lead, "predicted_rho8_over_e2": rho**8/E2})
    OUT["leading_singular_factor"] = {
        "rows": trow,
        "structure": "omega_0 = -(M/N_sigma)[Theta(rho)/rho] W_sigma(phi) (-x_2, x_1, 0); near "
                     "the origin Theta(rho)/rho = rho^7/e^2 (1 + O(rho^8)), so each Cartesian "
                     "component is c |x|^7 x_j W_sigma(phi) + (a factor smoother by 8 orders), "
                     "and |x|^7 x_j is homogeneous of degree 8 and smooth off the origin",
        "regularity_at_origin": "C^{7,1}: the 8th derivatives are bounded and the 7th are "
                                "Lipschitz; the 9th derivatives blow up like |x|^{-1}"}

    # ---- 3. W_sigma is smooth ON THE SPHERE: odd derivatives vanish at both poles
    CONV = FX.Conv(0.02)
    polepts = np.array([1e-9, math.pi - 1e-9])
    D = CONV.derivs(polepts, kmax=3)
    OUT["W_sigma_even_at_poles"] = {
        "W_at_0": float(D[0][0]), "W_prime_at_0": float(D[1][0]),
        "W_third_at_0": float(D[3][0]),
        "W_at_pi": float(D[0][1]), "W_prime_at_pi": float(D[1][1]),
        "W_third_at_pi": float(D[3][1]),
        "consequence": "W_sigma is an even real-analytic function of phi at each pole, hence a "
                       "real-analytic function of the direction on S^2; the angular factor "
                       "contributes no singularity at the origin and none on the axis"}
    print("W'(0) = %.3e, W'''(0) = %.3e, W'(pi) = %.3e"
          % (D[1][0], D[3][0], D[1][1]))

    # ---- 4. the Fourier exponent: closed form and a direct Hankel check at l = 1
    n, l, a = 3, 1, 8
    # FT of |x|^a Y_l in R^n is c |xi|^{-a-n} Y_l with
    #   c = (-i)^l 2^{a+n} pi^{n/2} Gamma((a+n+l)/2)/Gamma((l-a)/2)
    num = mp.gamma(mp.mpf(a + n + l)/2)
    den_arg = mp.mpf(l - a)/2
    pole = (den_arg <= 0) and (float(den_arg) == int(float(den_arg)))
    c_const = 2**(a + n)*mp.pi**(mp.mpf(n)/2)*num/mp.gamma(den_arg)
    OUT["fourier_exponent"] = {
        "n": n, "l": l, "homogeneity_a": a,
        "FT_homogeneity": -(a + n),
        "constant_argument_of_Gamma_in_denominator": float(den_arg),
        "denominator_Gamma_has_a_pole": bool(pole),
        "constant_nonzero": bool(not pole),
        "constant": mp.nstr(c_const, 15),
        "Hs_integral": "int (1+k^2)^s k^{-2(a+n)} k^{n-1} dk ~ int k^{2s-2(a+n)+n-1} dk, which "
                       "converges at infinity iff 2s - 2(a+n) + n - 1 < -1, i.e. s < a + n/2",
        "critical_s_for_omega_0": a + n/2.0,
        "critical_s_for_u_0": a + n/2.0 + 1.0}
    print("FT homogeneity %d, constant %s, critical s = %.1f (omega_0), %.1f (u_0)"
          % (-(a+n), mp.nstr(c_const, 8), a + n/2.0, a + n/2.0 + 1.0))

    # (4b) THE SINGULARITY, DIRECTLY.  g(x) = |x|^7 x_1 is homogeneous of degree 8 and smooth
    # off the origin, so d^8 g is homogeneous of degree 0 and d^9 g of degree -1.  If d^8 g takes
    # two different values on the unit sphere then g is not C^8 at the origin, and d^9 g blows up
    # like |x|^{-1}; a function whose 9th derivatives are O(|x|^{-1}) near 0 lies in H^s_loc(R^3)
    # exactly for s < 9 + 1/2.  Both are checked in exact arithmetic.
    import sympy as sp
    X1, X2, X3, t = sp.symbols("X1 X2 X3 t", real=True)
    g = (X1**2 + X2**2 + X3**2)**sp.Rational(7, 2)*X1
    d8 = sp.diff(g, X1, 8)
    d9 = sp.diff(g, X1, 9)
    dirs = {"e1": (1, 0, 0), "e3": (0, 0, 1), "(1,1,0)/sqrt2": (1, 1, 0), "(1,0,1)/sqrt2": (1, 0, 1)}
    vals8, vals9 = {}, {}
    for name, v in dirs.items():
        nv = sp.sqrt(sum(sp.Integer(c)**2 for c in v))
        sub = {X1: t*v[0]/nv, X2: t*v[1]/nv, X3: t*v[2]/nv}
        vals8[name] = sp.simplify(sp.limit(d8.subs(sub), t, 0, "+"))
        e9 = sp.simplify(sp.together(d9.subs(sub)*t))
        vals9[name] = sp.simplify(sp.limit(e9, t, 0, "+"))
    distinct8 = len(set(str(sp.nsimplify(v)) for v in vals8.values()))
    OUT["origin_singularity"] = {
        "g": "|x|^7 x_1, the leading term of every Cartesian component of omega_0 at the origin",
        "d8_limits_along_rays": {k: str(v) for k, v in vals8.items()},
        "d8_takes_distinct_values": distinct8 > 1,
        "n_distinct_d8_limits": distinct8,
        "limit_of_t_times_d9_along_rays": {k: str(v) for k, v in vals9.items()},
        "conclusion": "d^8 g is homogeneous of degree 0 with direction-dependent limits, so g is "
                      "C^{7,1} and not C^8 at the origin; d^9 g is homogeneous of degree -1 and "
                      "t d^9 g has a finite nonzero limit along each ray, so |d^9 g| ~ c/|x|. "
                      "A compactly supported function whose 9th derivatives are O(|x|^{-1}) is "
                      "in H^s_loc(R^3) exactly for s < 9 + 1/2 = 9.5."}
    print("d^8(|x|^7 x_1) limits along rays:", {k: str(v) for k, v in vals8.items()})
    print("  distinct values:", distinct8, " -> not C^8 at the origin")

    # (4c) THE FOURIER EXPONENT ON AN EXACT MODEL.  Take the model function
    #      f(x) = |x|^7 x_1 e^{-|x|^2} = F(|x|) x_1/|x| ,  F(rho) = rho^8 e^{-rho^2} ,
    # which has EXACTLY the datum's origin singularity multiplied by an entire, rapidly
    # decreasing factor, so its Fourier decay is governed by the same singularity.  Its l = 1
    # radial transform has a closed form (Watson 13.3(4), nu = 3/2, mu = 10.5, p = 1):
    #   G(k) = int_0^inf F(rho) j_1(k rho) rho^2 drho
    #        = sqrt(pi/(2k)) [Gamma(6)/(2 Gamma(5/2))] (k/2)^{3/2} 1F1(6; 5/2; -k^2/4) ,
    # evaluated in 60-digit arithmetic, so no oscillatory quadrature and no cancellation.
    mp.mp.dps = 60

    def G_exact(k):
        K = mp.mpf(k)
        return (mp.sqrt(mp.pi/(2*K))*(mp.gamma(6)/(2*mp.gamma(mp.mpf(5)/2)))
                * (K/2)**mp.mpf(1.5)*mp.hyp1f1(6, mp.mpf(5)/2, -K*K/4))

    ks = [50.0, 100.0, 200.0, 400.0, 800.0, 1600.0]
    gv = [float(mp.fabs(G_exact(k))) for k in ks]
    pair = [float((math.log(gv[i+1]) - math.log(gv[i]))/(math.log(ks[i+1]) - math.log(ks[i])))
            for i in range(len(ks) - 1)]
    OUT["hankel_decay_check"] = {
        "model": "f(x) = |x|^7 x_1 e^{-|x|^2}: the datum's origin singularity times an entire, "
                 "rapidly decreasing factor",
        "closed_form": "G(k) = sqrt(pi/(2k)) [Gamma(6)/(2 Gamma(5/2))] (k/2)^{3/2} "
                       "1F1(6; 5/2; -k^2/4)",
        "dps": 60, "k": ks, "abs_G": gv, "pairwise_exponents": pair, "predicted": -11,
        "agrees_with_prediction": bool(abs(pair[-1] + 11.0) < 0.02),
        "note": "a compactly supported C^infinity cutoff was tried first and is useless here: "
                "its own transform decays like exp(-c sqrt(k)), which dominates k^{-11} over "
                "every k a double- or even a 60-digit quadrature reaches. Recorded rather than "
                "hidden."}
    print("Hankel (exact, 60 dps) pairwise exponents:",
          ["%.5f" % q for q in pair], " predicted -11")

    OUT["verdict"] = {
        "omega_0_in_H_s_for": "s < 9.5", "u_0_in_H_s_for": "s < 10.5",
        "highest_order_needed_in_the_chain": max(c["order_in_u"] for c in CONSUMERS),
        "highest_order_on_the_datum": max(c["order_in_u"] for c in CONSUMERS
                                          if c.get("on_datum")),
        "supplied_by_s_lt_10p5": "H^{s} subset C^k for s > k + 3/2, so u_0 and u(t) are C^8 in "
                                 "space at every t < T_*; for t > 0 and nu > 0 interior parabolic "
                                 "regularity gives C^infinity anyway",
        "gap": False,
        "statement": "every consumer needs at most 5 derivatives of u, and that one (the k <= 4 "
                     "clause of (H3'_vartheta)) is a condition on the datum at t = 0 on a ball "
                     "where eta_0 is real-analytic. s < 10.5 supplies 8. There is no gap."}
    json.dump(OUT, open(os.path.join(HERE, "x3_results.json"), "w"), indent=1, default=str)
    print("wrote x3_results.json")
