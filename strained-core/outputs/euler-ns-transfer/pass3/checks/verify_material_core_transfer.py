"""Exact algebra behind material-core-transfer.md; no NS trajectory is simulated."""
import sympy as s


def check(name, expression):
    value = s.simplify(expression)
    assert value == 0, (name, value)
    print(f"PASS: {name}")


t, Omega, K, nu = s.symbols("t Omega K nu", real=True)
h11, h12, h22 = s.symbols("h11 h12 h22", real=True)
L2 = s.Matrix([[0, -Omega], [Omega, 0]])
H = s.Matrix([[h11, h12], [h12, h22]])
F = s.eye(2) + t*L2 + t**2*(H + L2**2)/2
detF = s.expand(F.det())
check("no first-order material area change", detF.coeff(t,1))
check("second-order material area coefficient", detF.coeff(t,2) - (h11+h22)/2)
check("incompressible axial/transverse area cancellation", detF.coeff(t,2).subs(h22,-K-h11) + K/2)
metric = s.expand(F.T*F)
for i in range(2):
    for j in range(2):
        check(f"transverse metric second coefficient {i},{j}", s.expand(metric[i,j]).coeff(t,2)-H[i,j])
flux = s.expand(2*Omega*(1+K*t*t/2)*(1-K*t*t/2))
check("central second-order vorticity-area flux cancellation", flux.coeff(t,2))

x,y,z = s.symbols("x y z", real=True)
variables = (x,y,z)
xx = s.Matrix(variables)
L = s.Matrix([[0,-Omega,0],[Omega,0,0],[0,0,0]])
base = L*xx
f = s.Function("f")(*variables)
Phi = s.Function("Phi")(*variables)
Pi = s.Function("Pi")(*variables)


def grad(g):
    return s.Matrix([s.diff(g,q) for q in variables])


def lap(g):
    return sum(s.diff(g,q,2) for q in variables)


def lapv(v):
    return s.Matrix([lap(g) for g in v])


def rotate(g):
    return sum(base[i]*s.diff(g,variables[i]) for i in range(3))


def rotatev(v):
    return s.Matrix([rotate(g) for g in v])


check("rigid rotation commutes with scalar Laplacian", lap(rotate(f))-rotate(lap(f)))
check("skew pressure pairing vanishes for potential acceleration", s.trace(L*s.hessian(Phi,variables)))

u1 = grad(Phi)
u2 = -L*u1 - rotatev(u1) - grad(Pi)
expected_lapu2 = -L*grad(lap(Phi)) - rotatev(grad(lap(Phi))) - grad(lap(Pi))
for i in range(3):
    check(f"second velocity jet harmonic reduction {i}", lap(u2[i])-expected_lapu2[i])

omega2 = 2*Omega*u1.diff(z)
omega3 = -rotatev(omega2) + L*omega2 + 2*Omega*u2.diff(z)
for i in range(3):
    check(f"second vorticity jet harmonic reduction {i}", lap(omega2[i])-2*Omega*s.diff(lap(Phi),variables[i],z))
expected_lapomega3 = -rotatev(lapv(omega2)) + L*lapv(omega2) + 2*Omega*lapv(u2).diff(z)
for i in range(3):
    check(f"third vorticity jet harmonic reduction {i}", lap(omega3[i])-expected_lapomega3[i])

# Explicit nonquadratic harmonic inputs provide a separate nontrivial instance.
a,b,c,d = s.symbols("a b c d", real=True)
quartic = z**4-3*z*z*(x*x+y*y)+s.Rational(3,8)*(x*x+y*y)**2
mixed = z*z*(x*x-y*y)-(x**4-y**4)/6
phi_poly = a*(2*z*z-x*x-y*y)+b*(x*x-y*y)+c*quartic+d*mixed
pi_poly = x*(x*x-3*y*y)
check("chosen nonquadratic acceleration potential is harmonic",lap(phi_poly))
check("chosen first pressure jet is harmonic",lap(pi_poly))
u1_poly = grad(phi_poly)
u2_poly = -L*u1_poly-rotatev(u1_poly)-grad(pi_poly)
omega2_poly = 2*Omega*u1_poly.diff(z)
omega3_poly = -rotatev(omega2_poly)+L*omega2_poly+2*Omega*u2_poly.diff(z)
for field_name, field in [("u1",u1_poly),("u2",u2_poly),("omega2",omega2_poly),("omega3",omega3_poly)]:
    for i in range(3):
        check(f"nonquadratic local harmonic jet {field_name}[{i}]",lap(field[i]))

delta4,j3 = s.symbols("delta4 j3",real=True)
area_series=1-K*t*t/2+j3*t**3/6
viscous_flux_rate=s.expand(nu*area_series*delta4*t**4/s.factorial(4))
for order in range(4):
    check(f"viscous flux rate coefficient below order four: {order}",viscous_flux_rate.coeff(t,order))
check("first not-excluded circulation Taylor coefficient", s.integrate(viscous_flux_rate,t).coeff(t,5)-nu*delta4/s.factorial(5))

print("All algebra checks passed. Circulation is controlled through order four only; no finite-time conservation is asserted.")
