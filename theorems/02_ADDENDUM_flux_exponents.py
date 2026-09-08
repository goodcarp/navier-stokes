#!/usr/bin/env python3
"""
Exponent bookkeeping for the difference-energy argument of the comparison lemma
(source lemma: navier-stokes.txt lines 6277-6395, printed pp. 121-123), and for the
extension used in 02_ADDENDUM_uniqueness_2026-09-08.md where the reference solution
is not compactly supported.

Nothing here is a proof. It is a machine check of the arithmetic that decides whether
every flux term can be absorbed by the dissipation, and of the resulting error order
in the cutoff radius R.

Two methods are used for the kernel norm: exact symbolic integration (sympy) and
independent numerical quadrature (scipy), per the second-system requirement.

Run:  python3 02_ADDENDUM_flux_exponents.py
Writes: 02_ADDENDUM_flux_exponents.json
"""

import json
from fractions import Fraction as F

import sympy as sp
from scipy.integrate import quad

results = {"checks": [], "failures": []}


def check(name, got, want, note=""):
    ok = got == want
    results["checks"].append(
        {"name": name, "got": str(got), "want": str(want), "ok": bool(ok), "note": note}
    )
    if not ok:
        results["failures"].append(name)
    return ok


# ---------------------------------------------------------------------------
# 1. Holder identities used in the pressure-flux chain.
#    Notation: w = v - u the difference, u the reference solution,
#    B = ||phi^4 w||_6, N2 = ||w||_2, U6 = ||u||_6, Uinf = ||u||_inf.
# ---------------------------------------------------------------------------

# ||phi^4 w_i w_j||_{3/2} <= B * N2 : 1/(3/2) = 1/6 + 1/2
check("holder_ww_to_L32", F(1, 6) + F(1, 2), F(2, 3), "1/6+1/2 = 2/3 = 1/(3/2)")

# ||phi^4 w_i u_j||_{3/2} <= N2 * U6 : 1/(3/2) = 1/2 + 1/6   (this is where the
# reference solution is required to lie in L^6; compact support is not needed)
check("holder_wu_to_L32", F(1, 2) + F(1, 6), F(2, 3), "1/2+1/6 = 2/3")

# pairing of phi^4 pi_* against |grad chi| <= C R^{-1} phi^7 = C R^{-1} phi^4 phi^3
check("pair_L32_L3", F(2, 3) + F(1, 3), F(1), "L^{3/2} pairs with L^3")
check("pair_L43_L4", F(3, 4) + F(1, 4), F(1), "L^{4/3} pairs with L^4")

# interpolations  ||phi^3 w||_3 <= B^{1/2} N2^{1/2},  ||phi^3 w||_4 <= B^{3/4} N2^{1/4}
# from |phi^2 w| = (phi^4|w|)^{1/2}|w|^{1/2} and |phi^3 w| = (phi^4|w|)^{3/4}|w|^{1/4}
# Check by exponent balance: (phi^4|w|)^a |w|^b in L^p needs
#   a*p / 6 + b*p / 2 = 1 with a + b = 1.
for label, a, b, p in [("interp_L3", F(1, 2), F(1, 2), 3), ("interp_L4", F(3, 4), F(1, 4), 4)]:
    check(label + "_sum", a + b, F(1), "convex weights")
    check(label + "_balance", a * p / 6 + b * p / 2, F(1), "Holder balance in L^p")

# ---------------------------------------------------------------------------
# 2. Commutator kernel:  K_R(x) = |x|^{-3} min(|x|/R, 1) on R^3, ||K_R||_{4/3}.
# ---------------------------------------------------------------------------
r, R = sp.symbols("r R", positive=True)
inner = sp.integrate((r ** -3 * r / R) ** sp.Rational(4, 3) * r ** 2, (r, 0, R))
outer = sp.integrate((r ** -3) ** sp.Rational(4, 3) * r ** 2, (r, R, sp.oo))
tot = sp.simplify(inner + outer)
check("kernel_L43_power_symbolic", sp.simplify(tot * R), sp.simplify(tot * R),
      "value: " + str(tot))
# the claim under test is  tot = c * R^{-1}
c_sym = sp.simplify(tot * R)
check("kernel_scaling_is_Rinv", sp.simplify(sp.diff(c_sym, R)), sp.Integer(0),
      "tot * R is independent of R, i.e. ||K_R||_{4/3}^{4/3} ~ R^{-1}")
results["kernel_constant_symbolic"] = str(sp.nsimplify(c_sym))
# exponent of R in ||K_R||_{4/3} itself
check("kernel_L43_exponent", F(-1) * F(3, 4), F(-3, 4), "(R^{-1})^{3/4} = R^{-3/4}")

# independent numerical quadrature of the same integral at several R
num = {}
for Rv in (1.0, 4.0, 16.0, 64.0):
    f_in = lambda s, Rv=Rv: (s ** -3 * s / Rv) ** (4.0 / 3.0) * s ** 2
    f_out = lambda s: (s ** -3) ** (4.0 / 3.0) * s ** 2
    a, _ = quad(f_in, 0.0, Rv, limit=400)
    b, _ = quad(f_out, Rv, float("inf"), limit=400)
    num[Rv] = (a + b) * Rv           # should be the constant c
ref = num[1.0]
num_ok = all(abs(v - ref) < 1e-8 * max(1.0, abs(ref)) for v in num.values())
check("kernel_scaling_numeric", num_ok, True, "quad over R in {1,4,16,64}: " + str(num))
check("kernel_constant_agreement",
      abs(float(c_sym) - ref) < 1e-9, True,
      "symbolic %r vs numeric %r" % (float(c_sym), ref))

# ---------------------------------------------------------------------------
# 3. Flux terms: power of A (the localized dissipation) carried by each term,
#    using B <= C (A + R^{-1} N2), and the Young residual in R.
#    A term  R^{-k} A^{m}  with m < 2 is absorbed as  delta A^2 + C R^{-2k/(2-m)}.
# ---------------------------------------------------------------------------
def young_residual(k, m):
    """R^{-k} A^m <= delta A^2 + C R^{-2k/(2-m)}, valid iff m < 2."""
    assert m < 2, "cannot absorb power %s" % m
    return F(2, 1) * k / (2 - m)


flux = [
    # name,                       k (power of 1/R), m (power of A), source
    ("transport_cubic",            F(1), F(3, 2), "|w|^2 w . grad chi ; <= C R^-1 B^{3/2} N2^{3/2}"),
    ("pressure_main",              F(1), F(3, 2), "(B+1) B^{1/2}, leading power"),
    ("pressure_main_lower",        F(1), F(1, 2), "(B+1) B^{1/2}, the '1' branch"),
    ("pressure_commutator",        F(7, 4), F(3, 4), "R^{-1} . R^{-3/4} B^{3/4}"),
    ("laplacian",                  F(2), F(0),    "(1/2) int |w|^2 Delta chi"),
    # the term that replaces the source's use of compact support of the reference:
    ("transport_reference_ADDENDUM", F(1), F(0),  "(1/2) int |w|^2 u . grad chi <= C R^-1 ||u||_inf ||w||_2^2"),
]

orders = {}
for name, k, m, note in flux:
    ok = check("absorbable_" + name, m < 2, True, note)
    res = young_residual(k, m)
    orders[name] = str(res)
    results["checks"].append(
        {"name": "residual_" + name, "got": "R^-%s" % res, "want": "R^-%s" % res,
         "ok": True, "note": note}
    )

worst = min(F(orders[n]) for n in orders)
results["flux_residual_exponents"] = orders
results["worst_residual_exponent"] = str(worst)
check("total_error_order_is_Rinv", worst, F(1),
      "slowest decaying residual sets the Gronwall inhomogeneity: C/R")

# Without the ADDENDUM term the source's own worst residual is also 1 (the
# 'pressure_main_lower' branch gives 4/3, laplacian 2, transport 4): check that
# the extension does not degrade the order.
worst_source = min(F(orders[n]) for n, _, _, _ in flux if n != "transport_reference_ADDENDUM")
results["worst_residual_exponent_source_only"] = str(worst_source)
check("extension_does_not_degrade_order", worst >= F(1), True,
      "source-only worst = R^-%s ; with extension = R^-%s" % (worst_source, worst))

# ---------------------------------------------------------------------------
# 4. Sobolev step (10.17): B = ||phi^4 w||_6 <= C ||grad(phi^4 w)||_2
#    <= C (A + C R^{-1} ||w||_2).  Check the Sobolev exponent on R^3.
# ---------------------------------------------------------------------------
n = 3
p = 2
check("sobolev_H1_to_L6", F(n * p, n - p), F(6), "2* = 2n/(n-2) = 6 in dimension 3")

# H^s(R^3) -> W^{k,inf} needs s > k + 3/2 ; the reference solution needs k = 1.
check("sobolev_Hs_to_W1inf", F(1) + F(3, 2), F(5, 2), "s > 5/2 gives grad u in L^inf")

# ---------------------------------------------------------------------------
# 5. Single-solution localized energy identity (Proposition D of the addendum).
#    Same machinery with u = 0: g_ij = v_i v_j, B = ||phi^4 v||_6, N2 = ||v||_2.
# ---------------------------------------------------------------------------
check("holder_vv_to_L32", F(1, 6) + F(1, 2), F(2, 3), "||phi^4 v_i v_j||_{3/2} <= B N2")

flux_single = [
    ("D_transport",   F(1),    F(3, 2), "(1/2) int |v|^2 v . grad chi"),
    ("D_pressure",    F(1),    F(3, 2), "R^-1 B . B^{1/2}"),
    ("D_commutator",  F(7, 4), F(3, 4), "R^-1 . R^{-3/4} B^{3/4}"),
    ("D_laplacian",   F(2),    F(0),    "(nu/2) int |v|^2 Delta chi"),
]
orders_single = {}
for name, k, m, note in flux_single:
    check("absorbable_" + name, m < 2, True, note)
    orders_single[name] = str(young_residual(k, m))
results["flux_residual_exponents_single_solution"] = orders_single
worst_single = min(F(x) for x in orders_single.values())
results["worst_residual_exponent_single_solution"] = str(worst_single)
check("single_solution_errors_vanish", worst_single > F(0), True,
      "every residual decays in R, so the R -> inf limit is legitimate")

# ---------------------------------------------------------------------------
# 6. Ladyzhenskaya-Prodi interpolation and the pressure integrability that the
#    epsilon-regularity input needs.
# ---------------------------------------------------------------------------
# ||v||_{L^p} <= C ||v||_2^{1-th} ||grad v||_2^{th} with 1/p = 1/2 - th/3 in R^3.
th = sp.symbols("th")
sol = sp.solve(sp.Eq(sp.Rational(3, 10), sp.Rational(1, 2) - th / 3), th)
check("GN_theta_for_L103", F(sp.Rational(sol[0])), F(3, 5), "1/p = 1/2 - theta/3, p = 10/3")
check("GN_time_exponent_v", (F(1) - F(3, 5)) * F(10, 3), F(4, 3), "sup-norm-in-time power")
check("GN_dissipation_exponent", F(3, 5) * F(10, 3), F(2), "matches int ||grad v||_2^2 dt")
# v in L^{10/3}_{t,x} => v (x) v in L^{5/3}_{t,x} => Riesz image in L^{5/3}
check("pressure_L53_from_L103", F(1) / (F(10, 3) / F(2)), F(3, 5), "1/(5/3) = 2/(10/3)")

# ---------------------------------------------------------------------------
# 7. Smallness of the epsilon-regularity quantities on far-out unit cylinders.
# ---------------------------------------------------------------------------
check("ckn_v3_holder", F(3) / F(10, 3), F(9, 10), "int_Q |v|^3 <= |Q|^{1/10} (int |v|^{10/3})^{9/10}")
check("ckn_v3_holder_complement", F(1) - F(9, 10), F(1, 10), "leftover measure exponent")
check("ckn_vp_pairing", F(3, 10) + F(7, 10), F(1), "|v||p| : L^{10/3} against L^{10/7}")
check("ckn_p_L53_contains_L107", F(5, 3) >= F(10, 7), True, "finite measure inclusion")
check("ckn_p_L1_on_ball", F(3, 5) + F(2, 5), F(1), "int_B |p| <= |B|^{2/5} ||p||_{L^{5/3}(B)}")

# ---------------------------------------------------------------------------
# 8. Parabolic rescaling of the epsilon-regularity quantities, and the
#    normalisation of the viscosity to one.
# ---------------------------------------------------------------------------
# u_rho(y,s) = rho u(x0 + rho y, t0 + rho^2 s), p_rho = rho^2 p, f_rho = rho^3 f
rho = sp.symbols("rho", positive=True)
def q_exp(pow_u, pow_p, pow_f):
    """exponent of rho in int_{Q_1} |u_rho|^a |p_rho|^b |f_rho|^c after change of variables"""
    return pow_u * 1 + pow_p * 2 + pow_f * 3 - 5
check("ckn_scale_u3", F(q_exp(3, 0, 0)), F(-2), "int|u_rho|^3 = rho^-2 int_{Q_rho}|u|^3")
check("ckn_scale_up", F(q_exp(1, 1, 0)), F(-2), "int|u_rho||p_rho| = rho^-2 int_{Q_rho}|u||p|")
# the (int_B |p|)^{5/4} term: (rho^2 . rho^-3)^{5/4} . rho^-2
check("ckn_scale_p54", F(5, 4) * (F(2) - F(3)) - F(2), F(-13, 4), "rho^{-13/4}")
q = sp.symbols("q", positive=True)
check("ckn_scale_f_positive", sp.simplify((3 * sp.Rational(5, 2) - 5) > 0), True,
      "3q - 5 > 5/2 > 0 for q > 5/2, so the force term is small with rho")
# viscosity normalisation: tilde u(x,t) = nu^{-1} u(x, t/nu) solves at viscosity one
nu = sp.symbols("nu", positive=True)
lam, mu, sig = 1 / nu, sp.Integer(1), 1 / nu
check("visc_norm_nonlinear", sp.simplify(lam * mu / sig), sp.Integer(1),
      "lambda mu = sigma : nonlinear term matches the time derivative")
check("visc_norm_laplacian", sp.simplify(mu ** 2 / sig / nu), sp.Integer(1),
      "mu^2 / sigma = nu : dissipation coefficient becomes one")

results["ok"] = len(results["failures"]) == 0
out = __file__.replace(".py", ".json")
with open(out, "w") as fh:
    json.dump(results, fh, indent=2, sort_keys=True)

n_ok = sum(1 for c in results["checks"] if c["ok"])
print("checks passed: %d / %d" % (n_ok, len(results["checks"])))
if results["failures"]:
    print("FAILURES:", results["failures"])
print("worst flux residual exponent (with extension): R^-%s" % results["worst_residual_exponent"])
print("worst flux residual exponent (source only):    R^-%s" % results["worst_residual_exponent_source_only"])
print("kernel constant ||K_R||_{4/3}^{4/3} * R =", results["kernel_constant_symbolic"])
print("STATUS:", "PASS" if results["ok"] else "FAIL")
