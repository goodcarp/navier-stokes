"""Exact model identities for affine-stretch-gate.md, not a finite-energy PDE proof."""
import sympy as s


def check(name, expression):
    value = s.simplify(expression)
    assert value == 0, (name, value)
    print(f"PASS: {name}")


t, r = s.symbols("t r", positive=True)
nu, circulation = s.symbols("nu circulation", positive=True)
width = s.Function("width")(t)
strain = s.Function("strain")(t)
omega = circulation / (s.pi * width) * s.exp(-r**2 / width)
swirl = circulation / (2 * s.pi * r) * (1 - s.exp(-r**2 / width))
width_equation = {s.diff(width, t): -strain * width + 4 * nu}

check(
    "Gaussian stretched-vorticity full scalar residual",
    (
        s.diff(omega, t) - strain * r * s.diff(omega, r) / 2
        - strain * omega - nu * (s.diff(omega, r, 2) + s.diff(omega, r) / r)
    ).subs(width_equation),
)
check(
    "Gaussian stretched-swirl full azimuthal momentum residual",
    (
        s.diff(swirl, t) - strain * r * s.diff(swirl, r) / 2
        - strain * swirl / 2
        - nu * (s.diff(swirl, r, 2) + s.diff(swirl, r) / r - swirl / r**2)
    ).subs(width_equation),
)
w_positive = s.symbols("w_positive", positive=True)
check(
    "Gaussian circulation normalization",
    s.integrate(2 * s.pi * r * omega.subs(width, w_positive), (r, 0, s.oo)) - circulation,
)

# Check an anisotropic Kelvin wave directly in the full unforced NS equation.
x, y, z = s.symbols("x y z", real=True)
a1, a2 = s.symbols("a1 a2", positive=True)
k, amplitude = s.Function("k")(t), s.Function("amplitude")(t)
velocity = s.Matrix([
    -a1 * x + amplitude * s.sin(k * y),
    -a2 * y,
    (a1 + a2) * z,
])
coordinates = s.Matrix([x, y, z])
pressure = -(a1**2 * x**2 + a2**2 * y**2 + (a1 + a2)**2 * z**2) / 2
substitutions = {
    s.diff(k, t): a2 * k,
    s.diff(amplitude, t): (a1 - nu * k**2) * amplitude,
}
check("Kelvin incompressibility", s.trace(velocity.jacobian(coordinates)))
residual = (
    velocity.diff(t) + velocity.jacobian(coordinates) * velocity
    + s.Matrix([s.diff(pressure, q) for q in coordinates])
    - nu * s.Matrix([sum(s.diff(component, q, 2) for q in coordinates) for component in velocity])
).subs(substitutions)
for index, entry in enumerate(residual):
    check(f"exact anisotropic Kelvin full NS component {index + 1}", entry)

vorticity_amplitude = -k * amplitude
check(
    "Kelvin axial vorticity stretching with evolving wavevector",
    s.diff(vorticity_amplitude, t).subs(substitutions)
    - (a1 + a2 - nu * k**2) * vorticity_amplitude,
)
check("affine material volume conservation", s.exp(-a1*t) * s.exp(-a2*t) * s.exp((a1+a2)*t) - 1)

X, eta = s.symbols("X eta", positive=True)
vorticity_gain = X * s.exp(-eta * (X - 1))
velocity_gain = s.sqrt(X) * s.exp(-eta * (X - 1))
check("vorticity gain derivative", s.diff(vorticity_gain, X) - s.exp(-eta*(X-1)) * (1-eta*X))
check("velocity gain derivative", s.diff(velocity_gain, X) - s.exp(-eta*(X-1)) / s.sqrt(X) * (s.Rational(1,2)-eta*X))
check("vorticity gain maximizer", s.diff(vorticity_gain, X).subs(X, 1/eta))
check("velocity gain maximizer", s.diff(velocity_gain, X).subs(X, 1/(2*eta)))

U0, Lx0, k0, A1, A2, D = s.symbols("U0 Lx0 k0 A1 A2 D", positive=True)
U = U0 * s.exp(A1 - D)
kt = k0 * s.exp(A2)
Lx = Lx0 * s.exp(-A1)
check("anisotropic Reynolds-localization product", U/(nu*kt) * (kt*Lx) - U0*Lx0/nu*s.exp(-D))
check("isotropic carrier Reynolds has only viscous loss", (U/(nu*kt)).subs(A1, A2) - U0/(nu*k0)*s.exp(-D))

print("\nAll identities passed. Global examples have infinite energy; localization and regeneration are unproved.")
