"""Check auxiliary viscosity-audit identities with SymPy.

This is algebra verification, not a PDE solver or a Lean-proof verification.
No source code from fluid_lean is imported or executed.
"""

import sympy as s


def check(name, expression):
    assert s.simplify(expression) == 0, (name, s.simplify(expression))
    print(f"PASS: {name}")


r = s.symbols("r", positive=True)
z, y = s.symbols("z y", real=True)
G = s.Function("G")
g = G(z, r**2 / 2)
g_zz = s.Subs(s.diff(G(z, y), z, 2), y, r**2 / 2)
g_yy = s.Subs(s.diff(G(z, y), y, 2), y, r**2 / 2)
g_y = s.Subs(s.diff(G(z, y), y), y, r**2 / 2)

check(
    "circulation diffusion in volume coordinates",
    s.diff(g, r, 2) - s.diff(g, r) / r + s.diff(g, z, 2)
    - (g_zz + r**2 * g_yy),
)
check(
    "reduced-vorticity diffusion in volume coordinates",
    s.diff(g, r, 2) + 3 * s.diff(g, r) / r + s.diff(g, z, 2)
    - (g_zz + r**2 * g_yy + 4 * g_y),
)

zeta1, zeta2 = s.symbols("zeta1 zeta2", real=True)
kappa = zeta1**2 / r**2 + zeta2**2
check("physical diffusion symbol q = r^2 kappa", r**2 * kappa - zeta1**2 - r**2 * zeta2**2)

a, b, damping, eigenvalue = s.symbols("a b damping eigenvalue", real=True)
matrix = s.Matrix([[-damping, a], [b, -damping]])
check(
    "equal damping shifts the characteristic polynomial",
    (eigenvalue * s.eye(2) - matrix).det() - ((eigenvalue + damping)**2 - a*b),
)

T, O, c0, ratio = s.symbols("T O c0 ratio", nonzero=True)
T_dot = a * O - damping * T
O_dot = b * T - damping * O
ratio_dot = (c0 * (O_dot * T - O * T_dot) / T**2).subs(O, ratio * T / c0)
check("equal damping cancels in the projective ratio", ratio_dot - (c0*b - a*ratio**2/c0))
check("absolute amplitude retains damping", (T_dot/T).subs(O, ratio*T/c0) - (a*ratio/c0-damping))

alpha, Qm, Qprev, log_Nm = s.symbols("alpha Qm Qprev log_Nm", positive=True)
beta = s.Rational(1, 8)
exponent = 2*alpha - beta/(2*Qm*Qprev)
check(
    "two-generation frequency conversion",
    2*alpha*log_Nm - (beta/2)*(log_Nm/(Qm*Qprev)) - exponent*log_Nm,
)
check("ordinary viscosity exponent", exponent.subs(alpha, 1) - (2-s.Rational(1,16)/(Qm*Qprev)))

print("\nExact scale exponent:", exponent)
print("Once Qm*Qprev >= 1/(16*alpha), this exponent is at least alpha > 0.")
print("The analytic limit also uses log(Lambda_m) = o(log(N_m)); see viscosity-audit.md.")
