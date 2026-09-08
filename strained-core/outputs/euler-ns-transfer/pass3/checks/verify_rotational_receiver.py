#!/usr/bin/env python3
"""Exact finite checks for rotational-receiver.md; no PDE simulation or kernel proof.

Requires SymPy. The analytic lemma still depends on local smooth NS existence,
the full initial pressure calculation, radial integration and the stated norm bound.
"""

import sympy as s

x, y, z = s.symbols("x y z", real=True)
xyz = (x, y, z)
Omega = s.symbols("Omega", positive=True)
J = s.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])
X = s.Matrix(xyz)
checks = 0


def zero(expr, name):
    global checks
    values = list(expr) if isinstance(expr, s.MatrixBase) else [expr]
    assert all(s.simplify(v) == 0 for v in values), name
    checks += 1
    print(f"PASS {checks}: {name}")


def lap(p):
    return sum(s.diff(p, v, 2) for v in xyz)


def grad(p):
    return s.Matrix([s.diff(p, v) for v in xyz])


def homogeneous_harmonic(degree):
    """A symbolic linear combination of the entire harmonic space in this degree."""
    mons = [x**i * y**j * z**(degree-i-j)
            for i in range(degree+1) for j in range(degree+1-i)]
    low = [x**i * y**j * z**(degree-2-i-j)
           for i in range(degree-1) for j in range(degree-1-i)]
    matrix = s.Matrix([[s.Poly(lap(m), *xyz).coeff_monomial(n)
                        for m in mons] for n in low])
    basis = matrix.nullspace()
    assert len(basis) == 2 * degree + 1
    coeffs = s.symbols(f"h{degree}_0:{len(basis)}")
    return s.expand(sum(a * sum(v[i] * mons[i] for i in range(len(mons)))
                        for a, v in zip(coeffs, basis)))


def sphere_average(poly):
    """Exact normalized S² polynomial integral via even monomial moments."""
    ans = 0
    for powers, coefficient in s.Poly(s.expand(poly), *xyz).terms():
        if any(p % 2 for p in powers):
            continue
        numerator = s.prod(s.factorial2(p-1) for p in powers)
        denominator = s.factorial2(sum(powers)+1)
        ans += coefficient * numerator / denominator
    return s.simplify(ans)


psi_parts = {d: homogeneous_harmonic(d) for d in (2, 4, 6)}
for degree, psi in psi_parts.items():
    zero(lap(psi), f"generic degree-{degree} harmonic potential ({2*degree+1} coefficients)")
    # Degree 4 and 6 have nonconstant Hessians. Their angular first moments
    # must vanish, so no radial weight or radius can introduce a correction.
    moment_errors = []
    for i in range(3):
        for j in range(3):
            central = s.diff(psi, xyz[i], xyz[j]).subs({x: 0, y: 0, z: 0})
            moment_errors.append(sphere_average(xyz[i] * s.diff(psi, xyz[j]))
                                 - central / 3)
    zero(s.Matrix(moment_errors), f"degree-{degree} radial first-moment angular identities")

psi = sum(psi_parts.values())
u0 = Omega * J * X
ut = grad(psi)
zero(s.Matrix([lap(v) for v in ut]), "initial time derivative has zero local viscous Laplacian")
direct_utt = -u0.jacobian(xyz) * ut - ut.jacobian(xyz) * u0
projectable_utt = -2 * Omega * J * ut - grad(u0.dot(ut))
zero(direct_utt-projectable_utt, "second-time nonlinear identity before removing pressure gradient")

N_sphere = sphere_average((J * X).dot(J * X))
zero(N_sphere-s.Rational(2, 3), "radial rotational normalization is two-thirds of radial second moment")
K = s.diff(psi, z, 2).subs({x: 0, y: 0, z: 0})
modal_second = sum(sphere_average((-2 * Omega * J * grad(p)).dot(J * X))
                   for p in psi_parts.values())
zero(modal_second/(Omega*N_sphere)-K, "all generic harmonic orders give normalized a''(0)=K")

E, P, N, ar, omega_r = s.symbols("E P N ar omega_r", real=True)
# Weighted orthogonal decomposition: integral chi|u-omega_r Jx|².
zero((2*E-2*omega_r*P+omega_r**2*N).subs(P, omega_r*N)
     -(2*E-omega_r**2*N), "weighted projection decomposition underlying energy lower bound")
E0 = Omega**2*N/2
zero((omega_r**2*N/2).subs(omega_r, Omega*ar)-ar**2*E0,
     "initial rigid rotation saturates weighted energy normalization")

q = s.symbols("q", positive=True)
a = s.symbols("a", real=True)
zero(q*q**(-1-a)-q**(-a), "RMS velocity scale match implies Reynolds scale match")
zero(q**3*q**(-2*(1+a))-q**(1-2*a), "weight mass and RMS scaling give exact energy ratio")
M3, t = s.symbols("M3 t", positive=True)
zero((s.Rational(3, 2)*M3/Omega)*t**3/6-M3*t**3/(4*Omega),
     "uniform third derivative yields the recorded cubic remainder constant")

print(f"\n{checks} exact checks passed. No PDE evolution, full-pressure kernel proof, "
      "uniform Sobolev estimate, or profile regeneration is certified by this script.")
