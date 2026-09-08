"""Symbolic checks for the compact-core/remote-packet initial pressure jet."""
import sympy as s


def check(name, expression):
    result = s.simplify(expression)
    assert result == 0, (name, result)
    print(f"PASS: {name}")


x, y, z, Omega = s.symbols("x y z Omega", real=True)
q = x*x + y*y + z*z
psi = s.Function("psi")
u = s.Matrix([-Omega*y*psi(q), Omega*x*psi(q), 0])
J = u.jacobian(s.Matrix([x,y,z]))
check("compact swirl incompressibility", s.trace(J))
check(
    "compact swirl exact pressure source",
    s.trace(J*J) + Omega**2*(2*psi(q)**2 + 4*(x*x+y*y)*psi(q)*s.Subs(s.diff(psi(s.Symbol('q')), s.Symbol('q')), s.Symbol('q'), q)),
)

# Kernel below omits the common factor 1/(4*pi); pointwise formulas apply off zero.
N = q**(-s.Rational(1,2))
check("Newtonian Hessian regular part", s.diff(N,z,2) - (3*z*z-q)/q**s.Rational(5,2))
mu = s.symbols("mu", real=True)
angular = s.integrate((3*mu*mu-1)*(1-mu*mu), (mu,-1,1))/2
check("angular contraction", angular + s.Rational(4,15))
R = s.symbols("R", positive=True)
check(
    "radial integral is a boundary term",
    s.diff(psi(R*R)**2,R)/4
    - R*psi(R*R)*s.Subs(s.diff(psi(s.Symbol('q')), s.Symbol('q')),s.Symbol('q'),R*R),
)
check("core pressure Hessian coefficient", s.Rational(2,3) - 4*(-s.Rational(1,4))*angular - s.Rational(2,5))

horizontal = [x,y]
for i in range(2):
    for j in range(2):
        delta = 1 if i == j else 0
        expected = (
            (3-15*z*z/q)*delta
            + (105*z*z/q-15)*horizontal[i]*horizontal[j]/q
        )/q**s.Rational(5,2)
        check(f"horizontal pressure fourth derivative {i},{j}", s.diff(N,z,2,horizontal[i],horizontal[j])-expected)

m = s.symbols("m", real=True)
tangent = 3-15*m
radial = -12+105*m*(1-m)
check("cone tangent endpoint", tangent.subs(m,s.Rational(16,17))+s.Rational(189,17))
check("cone radial endpoint", radial.subs(m,s.Rational(16,17))+s.Rational(1788,289))
assert -s.Rational(189,17) < -6
assert -s.Rational(1788,289) < -6
# Both eigenvalues decrease on 16/17 <= m <= 1.
assert s.diff(tangent,m) < 0
assert s.diff(radial,m).subs(m,s.Rational(16,17)) < 0
assert s.diff(radial,m,2) < 0
print("PASS: both horizontal kernel eigenvalues are below -6 throughout the cone")
print("All algebra checks passed; this is an initial time-jet result, not a finite-gain or blowup proof.")
