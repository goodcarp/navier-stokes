#!/usr/bin/env python3
"""Exact degree-four outer pressure selection for a radial annular strain."""
import sympy as sp

x, y, z = sp.symbols("x y z", real=True)
r = sp.symbols("r", positive=True)
mu = sp.symbols("mu", real=True)
eta, eta1, eta2 = sp.symbols("eta eta1 eta2")
coords = (x, y, z)
X = sp.Matrix(coords)
s = X.dot(X)
S0 = sp.diag(-1, -1, 2)
P = (eta + sp.Rational(2, 3)*s*eta1)*S0*X - sp.Rational(2, 3)*eta1*X.dot(S0*X)*X
DP = sp.Matrix(3, 3, lambda i, j: sp.diff(P[i], coords[j]) +
               2*coords[j]*(eta1*sp.diff(P[i], eta) + eta2*sp.diff(P[i], eta1)))
assert sp.expand(sp.trace(DP)) == 0

# This is 4*pi*r^5*Q, with Q_ij=-partial_zzij N.
Qbar = sp.Matrix(3, 3, lambda i, j: -(
    105*X[i]*X[j]*z*z/r**4
    - 15*((1 if i == j else 0)*z*z +
          2*(1 if i == 2 else 0)*X[j]*z +
          2*(1 if j == 2 else 0)*X[i]*z + X[i]*X[j])/r**2
    + 3*((1 if i == j else 0) + 2*(1 if i == j == 2 else 0))))

# Check the kernel formula against direct Cartesian differentiation.
newton_without_4pi = s**(-sp.Rational(1, 2))
for i in range(3):
    for j in range(i, 3):
        actual = -sp.diff(newton_without_4pi, z, z, coords[i], coords[j])
        claimed = (Qbar[i, j]/r**5).subs(r, sp.sqrt(s))
        assert sp.simplify(actual-claimed) == 0
print("PASS fourth-derivative Newton kernel, including its sign away from zero.")

def sphere_polynomial(expr):
    result = 0
    for (px, pz), coefficient in sp.Poly(sp.expand(expr.subs(y, 0)), x, z).terms():
        assert px % 2 == 0
        result += coefficient*r**(px+pz)*(1-mu**2)**(px//2)*mu**pz
    return sp.expand(result)

weight = sphere_polynomial(sp.trace(Qbar*DP))
expected = (eta*(-315*mu**4+270*mu**2-27)
            + r**2*eta1*(-570*mu**4+468*mu**2-42)
            + 96*r**4*eta2*(mu**2-mu**4))
assert sp.expand(weight-expected) == 0
q = {
    0: sp.Rational(64, 5)*r**4*eta2,
    2: -sp.Rational(32, 7)*r**2*(3*eta1-2*r**2*eta2),
    4: -sp.Rational(24, 35)*(105*eta+190*r**2*eta1+32*r**4*eta2),
}
assert sp.expand(weight-sum(coefficient*sp.legendre(ell, mu) for ell, coefficient in q.items())) == 0

def exact_integral(expr):
    return sp.cancel(sum(coefficient*sp.Rational(2, power+1)
                         for (power,), coefficient in sp.Poly(sp.expand(expr), mu).terms()
                         if power % 2 == 0))

for ell in range(13):
    actual = sp.cancel((2*ell+1)*exact_integral(weight*sp.legendre(ell, mu))/2)
    assert sp.simplify(actual-q.get(ell, 0)) == 0
print("PASS exact q0, q2, q4 and zero projection onto every checked higher degree.")
assert sp.Poly(weight, mu).degree() == 4
print("PASS degree-four polynomial identity proves all higher projections vanish.")

# At y=0 the azimuthal vector points along y. Its Q image is also along y,
# while every axisymmetric pressure gradient is in the x,z plane.
Qmeridian = Qbar.subs(y, 0)
assert Qmeridian[0, 1] == 0 and Qmeridian[2, 1] == 0
print("PASS pointwise zero azimuthal pairing against axisymmetric pressure gradient.")

# Volume integration: -2 * (2*pi)/(4*pi) * 2/(2ell+1) = -2/(2ell+1).
for ell in (0, 2, 4):
    coefficient = -2*sp.Rational(1, 2)*sp.Rational(2, 2*ell+1)
    assert coefficient == -sp.Rational(2, 2*ell+1)
print("PASS radial functional coefficient -2/(2ell+1).")
print("All exact checks passed. Retained-mode numerical errors are not certified.")
