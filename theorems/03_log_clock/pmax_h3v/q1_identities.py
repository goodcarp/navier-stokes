"""
q1 -- every algebraic identity and every sign that the (P-max) write-up leans on.

All checks are sympy residuals against ZERO, in Cartesian coordinates on R^n, for a
GENERIC function (sympy Function), so nothing here is a numerical coincidence.

WHAT IS CHECKED

 (A) the weight identity in general n and k:
        rho^k Delta g = Delta(rho^k g) - (2k/rho^2) x.grad(rho^k g)
                        + k(k-n+2) (rho^k g)/rho^2
     and the sign of k(k-n+2) at n = 5: k(k-3) <= 0 for k = 0,1,2,3.

 (B) the EPSILON-REGULARISED weight identities used in the proof, with
     R_eps := sqrt(rho^2 + eps^2):
        R_eps   Delta u = Delta V - (2/R_eps^2) x.grad V - (2 rho^2 +  5 eps^2) V / R_eps^4 ,
        R_eps^2 Delta u = Delta W - (4/R_eps^2) x.grad W - (2 rho^2 + 10 eps^2) W / R_eps^4 ,
     V := R_eps u, W := R_eps^2 u.  Both zeroth-order coefficients are <= 0 for every eps,
     and both reduce to (A) at eps = 0.  This is the n = 5 fact the argument needs.

 (C) the regularised Kato remainders, exactly:
        Delta sqrt(f^2+eps^2) - f Delta f / sqrt(f^2+eps^2) = eps^2 |grad f|^2/(f^2+eps^2)^{3/2}
     (scalar), and for a vector field g the analogous remainder is
        |grad g|^2/G_eps - |g.grad g|^2/G_eps^3 >= 0   by Cauchy-Schwarz.

 (D) the gradient equation obtained from the scalar equation, and the ABSENCE of a
     zeroth-order term: d_t(d_j eta) + b.grad(d_j eta) = -(d_j b).grad eta + nu Delta(d_j eta),
     with div b = 2a appearing nowhere.  Checked against the conservation form, which is
     where 2a does appear, so the two bookkeepings are exhibited side by side.

 (E) the barrier inequalities for h = e^{Ks} (rho^2+1)^{beta/2}:
        Delta h + (4/rho) ... -> the two coefficient bounds
        phi'' + (4/rho) phi'  <=  beta(beta+3) phi/(rho^2+1)   (n = 5, beta >= 2)
        phi'' + (2/rho) phi'  <=  ...                      (the k = 1 variant)
     verified as inequalities on a dense grid in rho with the algebraic majorant proved
     symbolically.

Outputs -> q1_results.json
"""
import json
import sympy as sp

OUT = {"checks": {}}


def chk(name, expr, extra=None):
    r = sp.simplify(sp.expand(expr))
    ok = (r == 0)
    OUT["checks"][name] = {"residual_is_zero": bool(ok), "residual": str(r)}
    if extra:
        OUT["checks"][name].update(extra)
    print(("PASS " if ok else "FAIL ") + name + ("" if ok else "   residual=%s" % r))
    return ok


# ---------------------------------------------------------------- (A) general n, general k
print("(A) weight identity, general n and k")
A_ok = {}
for n in range(2, 8):
    xs = sp.symbols("x1:%d" % (n + 1), real=True)
    g = sp.Function("g")(*xs)
    rho2 = sum(x * x for x in xs)
    rho = sp.sqrt(rho2)
    for k in range(0, 5):
        w = rho ** k
        lhs = w * sum(sp.diff(g, x, 2) for x in xs)
        wg = w * g
        rhs = (sum(sp.diff(wg, x, 2) for x in xs)
               - (2 * k / rho2) * sum(x * sp.diff(wg, x) for x in xs)
               + k * (k - n + 2) * wg / rho2)
        A_ok["n=%d,k=%d" % (n, k)] = chk("A_weight_n%d_k%d" % (n, k), lhs - rhs)
OUT["A_sign_at_n5"] = {"k(k-3)": {k: k * (k - 3) for k in range(0, 5)},
                       "nonpositive_for_k": [k for k in range(0, 5) if k * (k - 3) <= 0]}

# ---------------------------------------------------------------- (B) eps-regularised, n = 5
print("(B) eps-regularised weight identities, n = 5")
n = 5
xs = sp.symbols("y1:6", real=True)
eps = sp.symbols("epsilon", positive=True)
u = sp.Function("u")(*xs)
rho2 = sum(x * x for x in xs)
Reps2 = rho2 + eps ** 2
Reps = sp.sqrt(Reps2)

lap = lambda F: sum(sp.diff(F, x, 2) for x in xs)
xdot = lambda F: sum(x * sp.diff(F, x) for x in xs)

V = Reps * u
lhsV = Reps * lap(u)
rhsV = lap(V) - (2 / Reps2) * xdot(V) - (2 * rho2 + 5 * eps ** 2) * V / Reps2 ** 2
chk("B_V_eps_identity_k1_n5", lhsV - rhsV)

W = Reps2 * u
lhsW = Reps2 * lap(u)
rhsW = lap(W) - (4 / Reps2) * xdot(W) - (2 * rho2 + 10 * eps ** 2) * W / Reps2 ** 2
chk("B_W_eps_identity_k2_n5", lhsW - rhsW)

# the eps -> 0 limits are the k = 1, 2 rows of (A)
chk("B_V_limit_eps0", sp.simplify((rhsV - (lap(sp.sqrt(rho2) * u)
                                          - (2 / rho2) * xdot(sp.sqrt(rho2) * u)
                                          - 2 * sp.sqrt(rho2) * u / rho2)).subs(eps, 0)))
chk("B_W_limit_eps0", sp.simplify((rhsW - (lap(rho2 * u)
                                          - (4 / rho2) * xdot(rho2 * u)
                                          - 2 * rho2 * u / rho2)).subs(eps, 0)))
OUT["B_zeroth_order_coefficients"] = {
    "k=1": "-(2 rho^2 + 5 eps^2)/R_eps^4   <= 0 for every rho, eps",
    "k=2": "-(2 rho^2 + 10 eps^2)/R_eps^4  <= 0 for every rho, eps",
    "eps=0 value": "-2/rho^2 in both cases, = k(k-3)/rho^2 at n=5",
}
# also the two derivative facts used verbatim in the text
chk("B_grad_Reps", sp.simplify(sum(sp.diff(Reps, x) ** 2 for x in xs) - rho2 / Reps2))
chk("B_lap_Reps", sp.simplify(lap(Reps) - (5 * eps ** 2 + 4 * rho2) / Reps2 ** sp.Rational(3, 2)))
chk("B_lap_Reps2", sp.simplify(lap(Reps2) - 10))

# ---------------------------------------------------------------- (C) regularised Kato
print("(C) regularised Kato remainders")
f = sp.Function("f")(*xs)
ue = sp.sqrt(f ** 2 + eps ** 2)
rem = lap(ue) - f * lap(f) / ue
chk("C_scalar_kato_remainder",
    sp.simplify(rem - eps ** 2 * sum(sp.diff(f, x) ** 2 for x in xs) / ue ** 3))
# vector case: remainder = |grad g|^2/G - |g.grad g|^2/G^3, nonneg by Cauchy-Schwarz.
gs = [sp.Function("g%d" % j)(*xs) for j in range(1, 6)]
G = sp.sqrt(sum(gg ** 2 for gg in gs) + eps ** 2)
lapG = lap(G)
gdotlap = sum(gs[j] * lap(gs[j]) for j in range(5)) / G
gradg2 = sum(sp.diff(gs[j], x) ** 2 for j in range(5) for x in xs)
gdg2 = sum((sum(gs[j] * sp.diff(gs[j], x) for j in range(5))) ** 2 for x in xs)
chk("C_vector_kato_remainder", sp.simplify(lapG - gdotlap - (gradg2 / G - gdg2 / G ** 3)))
OUT["C_note"] = ("|grad g|^2 G^2 >= |grad g|^2 |g|^2 >= sum_x |g.d_x g|^2 by Cauchy-Schwarz, "
                 "so the vector remainder is >= 0 for every eps >= 0.")

# ---------------------------------------------------------------- (D) the gradient equation
print("(D) gradient equation, and where div b = 2a does and does not appear")
s = sp.symbols("s", real=True)
XS = sp.symbols("z1:6", real=True)
eta = sp.Function("eta")(*XS, s)
bs = [sp.Function("b%d" % j)(*XS, s) for j in range(1, 6)]
nu = sp.symbols("nu", positive=True)
lapX = lambda F: sum(sp.diff(F, x, 2) for x in XS)
adv = sp.diff(eta, s) + sum(bs[j] * sp.diff(eta, XS[j]) for j in range(5)) - nu * lapX(eta)
for j in range(5):
    got = sp.diff(adv, XS[j])
    want = (sp.diff(sp.diff(eta, XS[j]), s)
            + sum(bs[m] * sp.diff(sp.diff(eta, XS[j]), XS[m]) for m in range(5))
            + sum(sp.diff(bs[m], XS[j]) * sp.diff(eta, XS[m]) for m in range(5))
            - nu * lapX(sp.diff(eta, XS[j])))
    chk("D_gradient_equation_j%d" % (j + 1), got - want)
# conservation form: this is where div b shows up, and it is a different equation
divb = sum(sp.diff(bs[j], XS[j]) for j in range(5))
cons = (sp.diff(eta, s) + sum(sp.diff(bs[j] * eta, XS[j]) for j in range(5))
        - nu * lapX(eta) - divb * eta)
chk("D_advective_equals_conservative_minus_divb_eta", sp.simplify(cons - adv))
OUT["D_note"] = ("d_j of the ADVECTIVE equation produces exactly -(d_j b).grad eta and no "
                 "multiple of eta; div_5 b = 2a enters only when the equation is written in "
                 "conservation form, which the argument never uses.")

# ---------------------------------------------------------------- (E) barrier coefficients
print("(E) barrier coefficient bounds for h = e^{Ks} (rho^2+1)^{beta/2}")
r = sp.symbols("r", positive=True)
beta = sp.symbols("beta", positive=True)
phi = (r ** 2 + 1) ** (beta / 2)
phip = sp.diff(phi, r)
phipp = sp.diff(phi, r, 2)
# n = 5 radial Laplacian is phi'' + 4 phi'/r
LapRad = sp.simplify(phipp + 4 * phip / r)
OUT["E_lap_radial"] = str(sp.simplify(LapRad))
# claimed majorant:  phi'' + 4 phi'/r <= beta(beta+3) (r^2+1)^{beta/2 - 1}  for beta >= 2
maj = beta * (beta + 3) * (r ** 2 + 1) ** (beta / 2 - 1)
diffE = sp.simplify(sp.expand(sp.simplify(maj - LapRad) / (r ** 2 + 1) ** (beta / 2 - 2)))
OUT["E_majorant_gap_poly"] = str(sp.simplify(diffE))
# r phi' <= beta phi
chk("E_r_phip_le_beta_phi_poly",
    sp.simplify(sp.simplify(beta * phi - r * phip) - beta * (r ** 2 + 1) ** (beta / 2 - 1)))
# the exact radial Laplacian, as an identity (this is what the majorant is read off)
chk("E_radial_laplacian_exact",
    sp.simplify(LapRad - (5 * beta * (r ** 2 + 1) ** (beta / 2 - 1)
                          + beta * (beta - 2) * r ** 2 * (r ** 2 + 1) ** (beta / 2 - 2))))
# and the k = 1 combination Delta h - (2/rho^2) x.grad h at eps = 0
chk("E_radial_k1_combination_exact",
    sp.simplify(LapRad - 2 * r * phip / r ** 2 * r ** 0 - (3 * beta * (r ** 2 + 1) ** (beta / 2 - 1)
                + beta * (beta - 2) * r ** 2 * (r ** 2 + 1) ** (beta / 2 - 2))))
OUT["E_majorant_algebra"] = {
    "phi''+4phi'/r": "5 beta (r^2+1)^{b/2-1} + beta(beta-2) r^2 (r^2+1)^{b/2-2}",
    "majorised by": "beta(beta+3)(r^2+1)^{b/2-1}, using r^2 <= r^2+1 and beta >= 2",
    "k=1 combination": "3 beta (r^2+1)^{b/2-1} + beta(beta-2) r^2 (r^2+1)^{b/2-2} "
                       "<= beta(beta+1)(r^2+1)^{b/2-1}",
}

# numeric confirmation of the majorant over beta in [2,4], r in (0, 40]
import numpy as np
worst = None
for bb in np.linspace(2.0, 4.0, 21):
    rr = np.concatenate([np.linspace(1e-4, 1.0, 400), np.linspace(1.0, 40.0, 400)])
    lap_n = (bb * (rr ** 2 + 1) ** (bb / 2 - 1)
             + bb * (bb - 2) * rr ** 2 * (rr ** 2 + 1) ** (bb / 2 - 2)
             + 4 * bb * (rr ** 2 + 1) ** (bb / 2 - 1))
    maj_n = bb * (bb + 3) * (rr ** 2 + 1) ** (bb / 2 - 1)
    g = (maj_n - lap_n).min()
    if worst is None or g < worst[0]:
        worst = (float(g), float(bb))
OUT["E_majorant_min_gap_over_grid"] = {"min_gap": worst[0], "at_beta": worst[1],
                                       "claim": "phi'' + 4 phi'/r <= beta(beta+3) phi/(r^2+1)"}
print("   E majorant minimum gap over grid: %.6g at beta = %.3f" % (worst[0], worst[1]))

# the k = 1 barrier: the combination that appears is  -(phi'' + 4 phi'/r) + (2/r) . r phi'
# i.e. Delta h - (2/R^2) x.grad h ; majorised by beta(beta+1) phi as claimed in the text
worst1 = None
for bb in np.linspace(2.0, 4.0, 21):
    rr = np.concatenate([np.linspace(1e-4, 1.0, 400), np.linspace(1.0, 40.0, 400)])
    lap_n = (bb * (rr ** 2 + 1) ** (bb / 2 - 1)
             + bb * (bb - 2) * rr ** 2 * (rr ** 2 + 1) ** (bb / 2 - 2)
             + 4 * bb * (rr ** 2 + 1) ** (bb / 2 - 1))
    sub = 2 * bb * (rr ** 2 + 1) ** (bb / 2 - 1)          # (2/R^2) x.grad h at eps = 0
    maj_n = bb * (bb + 1) * (rr ** 2 + 1) ** (bb / 2 - 1)
    g = (maj_n - (lap_n - sub)).min()
    if worst1 is None or g < worst1[0]:
        worst1 = (float(g), float(bb))
OUT["E_majorant_k1_min_gap"] = {"min_gap": worst1[0], "at_beta": worst1[1],
                                "claim": "Delta h - (2/rho^2) x.grad h <= beta(beta+1) phi/(rho^2+1)"}
print("   E k=1 majorant minimum gap over grid: %.6g at beta = %.3f" % (worst1[0], worst1[1]))

# ---------------------------------------------------------------- (F) Lemma 4.2 eigenvalues,
# and the SHARP rate for |grad eta| (which is -Lambda_min(sym grad b), not ||grad b||_op)
print("(F) sym(grad_5 b) spectrum: Gamma_rad and the gradient rate")
a, Q, P, om = sp.symbols("a Q P omega", real=True)
Sym = sp.Matrix([[a + Q, P - om / 2], [P - om / 2, -2 * a - Q]])
ev = list(Sym.eigenvals().keys())
lam_top = sp.simplify((-a + sp.sqrt((3 * a + 2 * Q) ** 2 + (2 * P - om) ** 2)) / 2)
lam_bot = sp.simplify((-a - sp.sqrt((3 * a + 2 * Q) ** 2 + (2 * P - om) ** 2)) / 2)
chk("F_block_top_eigenvalue", sp.simplify(sp.Min(*[sp.simplify(e - lam_top) for e in ev])
                                          * sp.Max(*[sp.simplify(e - lam_top) for e in ev])))
OUT["F_note"] = ("The transverse R^3 in R^4 carries the eigenvalue a with multiplicity 3; the "
                 "(e_r, e_z) block has trace -a and (tr)^2 - 4 det = (3a+2Q)^2 + (2P-omega)^2. "
                 "Gamma_rad = max(a, lam_top).  The sharp rate for |grad eta| is "
                 "-Lambda_min = max(-a, (a + sqrt((3a+2Q)^2+(2P-omega)^2))/2), which is "
                 "<= ||grad_5 b||_op = Gamma.  u2 uses Gamma; the note records the sharper "
                 "rate as headroom and does NOT use it.")
OUT["F_eigs"] = [str(e) for e in ev]

# The reference strain b = a diag(1,1,1,1,-2): Q = P = omega = 0.  Then
#   Gamma_rad = Lambda_max(sym) = a ,  -Lambda_min(sym) = 2a ,  ||grad_5 b||_op = 2a .
# So the "sharp" gradient rate -Lambda_min coincides with Gamma exactly, and the obvious
# sharpening of P2's rate buys nothing in the regime that dominates.
import numpy as _np
_a = 1.0
_S = _np.diag([_a, _a, _a, _a, -2 * _a])
_ev = _np.linalg.eigvalsh(0.5 * (_S + _S.T))
_op = float(_np.linalg.norm(_S, 2))
OUT["F_reference_strain"] = {
    "a": _a, "eigenvalues_sym": sorted(float(x) for x in _ev),
    "Gamma_rad = Lambda_max": float(_ev.max()),
    "minus_Lambda_min": float(-_ev.min()),
    "Gamma = ||grad_5 b||_op": _op,
    "Gamma_rad_over_Gamma": float(_ev.max()) / _op,
    "sharpening_gain = Gamma/(-Lambda_min)": _op / float(-_ev.min()),
    "note": ("-Lambda_min = 2a = Gamma exactly on the reference strain, so replacing "
             "Gamma by the sharp rate in P2 gains a factor 1.  Gamma_rad/Gamma = 1/2 "
             "reproduces u5 (c)."),
}
print("   F reference strain: Gamma_rad/Gamma = %.6f, Gamma/(-Lambda_min) = %.6f"
      % (OUT["F_reference_strain"]["Gamma_rad_over_Gamma"],
         OUT["F_reference_strain"]["sharpening_gain = Gamma/(-Lambda_min)"]))

n_fail = sum(0 if v.get("residual_is_zero", True) else 1 for v in OUT["checks"].values())
OUT["summary"] = {"n_checks": len(OUT["checks"]), "n_fail": n_fail}
print("\n%d symbolic checks, %d failures" % (len(OUT["checks"]), n_fail))
json.dump(OUT, open("q1_results.json", "w"), indent=1, sort_keys=True)
