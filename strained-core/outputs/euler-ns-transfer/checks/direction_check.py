"""Algebra checks for ../REPORT.md; not a PDE existence verification."""
import sympy as sp

r, z, y2 = sp.symbols("r z y2", positive=True)
G = sp.Function("G")
g = G(z, r**2 / 2)
gamma_operator = sp.diff(g, z, 2) + sp.diff(g, r, 2) - sp.diff(g, r) / r
xi_operator = sp.diff(g, z, 2) + sp.diff(g, r, 2) + 3 * sp.diff(g, r) / r
gamma_expected = (
    sp.diff(G(z, y2), z, 2) + r**2 * sp.diff(G(z, y2), y2, 2)
).subs(y2, r**2 / 2)
xi_expected = gamma_expected + 4 * sp.diff(G(z, y2), y2).subs(y2, r**2 / 2)
assert sp.simplify(gamma_operator - gamma_expected) == 0
assert sp.simplify(xi_operator - xi_expected) == 0

zz, zr, gz, gr, circulation = sp.symbols("zz zr gz gr circulation", real=True)
q = zz**2 + r**2 * zr**2
kappa = q / r**2
d = -zr * gz + zz * gr / r
c = 2 * circulation / r**4
lambda_squared = -d * c * zz / kappa
physical_lambda_squared = -2 * circulation / r**3 * (zz**2 * gr - zz * r * zr * gz) / q
assert sp.simplify(lambda_squared - physical_lambda_squared) == 0

sigma, decay, eig = sp.symbols("sigma decay eig", nonzero=True)
# Both diagonal viscous terms have the same principal decay.
B = sp.Matrix([[-decay, -d / (sigma * kappa)], [sigma * c * zz, -decay]])
assert sp.simplify((B - eig * sp.eye(2)).det() - ((eig + decay)**2 - lambda_squared)) == 0

alpha = sp.Rational(1, 4)
beta = sp.Rational(1, 16)
assert 0 < 2 * beta < alpha < sp.Rational(1, 2)
assert alpha - 2 * beta > 0
assert 1 - 2 * alpha > 0
assert 1 - alpha - 2 * beta > 0

print("PASS: both exact reduced diffusion operators")
print("PASS: physical-coordinate growth exponent identity")
print("PASS: frozen viscous characteristic polynomial")
print("PASS: rational scale example and positive summability exponents")
print("These checks establish algebra only; they do not construct an NS stage or singularity.")
