#!/usr/bin/env python3
"""Exact initial core response to exterior harmonic pressure; no quadrature."""
import sympy as sp

x, y, z = sp.symbols("x y z", real=True)
r = sp.symbols("r", positive=True)
mu = sp.symbols("mu", real=True)
p, p1, p2 = sp.symbols("p p1 p2")
b, H2, H4, moment = sp.symbols("b H2 H4 moment")
coords = (x, y, z)
s = x*x + y*y + z*z
X = sp.Matrix(coords)
S0 = sp.diag(-1, -1, 2)
J = sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])
strain = (p + sp.Rational(2, 3)*s*p1)*S0*X - sp.Rational(2, 3)*p1*(X.dot(S0*X))*X
swirl = (p + sp.Rational(2, 3)*s*p1)*J*X

def derivative(expr, coordinate):
    return sp.diff(expr, coordinate) + 2*coordinate*(p1*sp.diff(expr, p) + p2*sp.diff(expr, p1))

grad_strain = sp.Matrix(3, 3, lambda i, j: derivative(strain[i], coords[j]))
grad_swirl = sp.Matrix(3, 3, lambda i, j: derivative(swirl[i], coords[j]))
assert sp.expand(sp.trace(grad_strain)) == 0
assert sp.expand(sp.trace(grad_swirl)) == 0

def solid_harmonic(degree):
    result = 0
    for (power,), coefficient in sp.Poly(sp.legendre(degree, mu), mu).terms():
        result += coefficient*z**power*s**((degree-power)//2)
    result = sp.expand(result)
    assert sp.expand(sum(sp.diff(result, coordinate, 2) for coordinate in coords)) == 0
    return result

def angular_polynomial(expr):
    """Axisymmetry permits y=0; all remaining x powers in the scalar are even."""
    result = 0
    for (x_power, z_power), coefficient in sp.Poly(sp.expand(expr.subs(y, 0)), x, z).terms():
        assert x_power % 2 == 0
        result += coefficient*r**(x_power+z_power)*(1-mu**2)**(x_power//2)*mu**z_power
    return sp.expand(result)

def exact_angular_integral(expr):
    return sp.expand(sum(coefficient*sp.Rational(2, power+1)
                         for (power,), coefficient in sp.Poly(sp.expand(expr), mu).terms()
                         if power % 2 == 0))

expected = {
    0: 0,
    2: sp.Rational(32, 35)*r**4*p2 + sp.Rational(16, 5)*r**2*p1,
    4: sp.Rational(64, 7)*r**6*p2 + sp.Rational(1536, 35)*r**4*p1 + sp.Rational(144, 5)*r**2*p,
    6: 0,
    8: 0,
}
angular = {}
contacts = {}
for degree, expected_value in expected.items():
    harmonic = solid_harmonic(degree)
    hessian = sp.hessian(harmonic, coords)
    assert sp.expand(sp.trace(grad_swirl*hessian)) == 0
    contraction = sp.expand(sp.trace(grad_strain*hessian))
    angular[degree] = exact_angular_integral((3*mu**2-1)*angular_polynomial(contraction))
    assert sp.expand(angular[degree]-expected_value) == 0
    hessian_at_zero = hessian.subs({x: 0, y: 0, z: 0})
    # h=-2b tr(grad S Hess phi); the distributional Hessian contact is -h(0)/3.
    contacts[degree] = sp.Rational(2, 3)*b*sp.trace(S0*hessian_at_zero)
    assert contacts[degree] == (4*b if degree == 2 else 0)
    print(f"PASS degree {degree}: angular contraction = {angular[degree]}, contact = {contacts[degree]}")

# The radial integrals below follow by s=r^2 and compact-cutoff integration
# by parts, using psi(0)=1. All derivatives here are with respect to s.
# Dictionary entries are integral_0^infty r^(power-1) psi^(derivative)(r^2) dr.
radial_moments = {
    (4, 2): sp.Rational(1, 2),       # 1/2 integral s psi'' = 1/2
    (2, 1): -sp.Rational(1, 2),      # 1/2 integral psi' = -1/2
    (6, 2): moment,                 # 1/2 integral s^2 psi'' = integral psi
    (4, 1): -moment/2,              # 1/2 integral s psi' = -integral psi/2
    (2, 0): moment/2,               # 1/2 integral psi
}

def integrate_radially(expr):
    result = 0
    for cutoff_symbol, derivative_order in ((p, 0), (p1, 1), (p2, 2)):
        for (power,), coefficient in sp.Poly(sp.expand(expr).coeff(cutoff_symbol), r).terms():
            if coefficient:
                result += coefficient*radial_moments[(power, derivative_order)]
    return sp.expand(result)

radial2 = integrate_radially(angular[2])
radial4 = integrate_radially(angular[4])
assert radial2 == -sp.Rational(8, 7)
assert radial4 == sp.Rational(8, 5)*moment
response = sp.expand(H2*(contacts[2]-b*radial2) + H4*(contacts[4]-b*radial4))
assert response == sp.Rational(36, 7)*b*H2 - sp.Rational(8, 5)*b*moment*H4
print("PASS radial integration and contact: response =", response)

# Independently verify the radial integrations on two compactly supported
# polynomial profiles. These are finite-regularity algebraic probes; the note
# proves the identities for every smooth compact profile constant near zero.
t = sp.symbols("t", nonnegative=True)
for exponent in (4, 6):
    cutoff = (1-t)**exponent
    exact_moment = sp.integrate(cutoff, (t, 0, 1))
    for (power, derivative_order), claimed in radial_moments.items():
        direct = sp.integrate(t**(sp.Rational(power, 2)-1)*sp.diff(cutoff, t, derivative_order)/2, (t, 0, 1))
        assert sp.simplify(direct-claimed.subs(moment, exact_moment)) == 0
print("PASS independent polynomial-profile checks of all five radial moments.")

pzz, pzzzz = sp.symbols("pzz pzzzz")
derivative_response = response.subs({H2: pzz/2, H4: pzzzz/24})
assert sp.expand(derivative_response-(sp.Rational(18, 7)*b*pzz-b*moment*pzzzz/15)) == 0
print("PASS derivative form: (18/7)b p_zz - (b moment/15) p_zzzz.")
print("All exact algebraic checks passed. No numerical pressure sign is certified.")
