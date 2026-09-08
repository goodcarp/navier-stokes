#!/usr/bin/env python3
"""Exact scalar identities behind the circulation return obstruction."""
import sympy as sp

r = sp.symbols("r", positive=True)
ur, uz, w, wr, wz, wrr, wzz, wt = sp.symbols("ur uz w wr wz wrr wzz wt")
nu, force = sp.symbols("nu force", nonnegative=True)
Gamma_r, Gamma_z = w+r*wr, r*wz
Gamma_rr, Gamma_zz = 2*wr+r*wrr, r*wzz
swirl_equation = wt+ur*wr+uz*wz+ur*w/r-nu*(wrr+wr/r+wzz-w/r**2)-force
Gamma_equation = r*wt+ur*Gamma_r+uz*Gamma_z-nu*(Gamma_rr-Gamma_r/r+Gamma_zz)-r*force
assert sp.expand(Gamma_equation-r*swirl_equation) == 0
print("PASS exact Gamma equation, cylindrical diffusion, and weighted force source.")

q = sp.symbols("q", positive=True)
alpha = sp.symbols("alpha", positive=True)
u_value = sp.symbols("u_value")
old_Gamma = q*r*u_value
new_Gamma = r*q**(1+alpha)*u_value
assert sp.simplify(new_Gamma-q**alpha*old_Gamma) == 0
assert sp.simplify((q**alpha*sp.Symbol("G"))/(q**alpha*sp.Symbol("nu0"))-sp.Symbol("G")/sp.Symbol("nu0")) == 0
print("PASS Gamma and viscosity scale by the same q^alpha factor.")

# Receiver mass scales R^3; its rho^2 moment scales R^5.
R, mass1, moment1, G, Omega = sp.symbols("R mass1 moment1 G Omega", positive=True)
massR, momentR = R**3*mass1, R**5*moment1
receiver_upper = G*massR/momentR
assert sp.simplify(R**2*receiver_upper-G*mass1/moment1) == 0
lower_G = Omega*momentR/massR
assert sp.simplify(lower_G-(moment1/mass1)*R**2*Omega) == 0

# A finite-regularity weight is enough for this integral inequality;
# the actual theorem also permits any nonzero smooth nonnegative radial weight.
x = sp.symbols("x", nonnegative=True)
weight = (1-x*x)**2
dchi = sp.Rational(2,3)*sp.integrate(x**4*weight,(x,0,1))/sp.integrate(x*x*weight,(x,0,1))
assert dchi == sp.Rational(2,9)
print("PASS receiver-radius cancellation and an exact sample d_chi=2/9.")

# Exact finite-stage example: G0/g*=16, q=1/2, alpha=1/4.
contraction = sp.Rational(1,2)**sp.Rational(1,4)
assert sp.simplify(16*contraction**16) == 1
assert sp.simplify(16*contraction**17) < 1
print("PASS finite cap example: 16 returns possible under the bound; stage 17 excluded.")

# Taylor core lower bound with a<=omega/M.
a, omega, M = sp.symbols("a omega M", positive=True)
error = M*a*a/2
assert sp.simplify(error.subs(M,omega/a)-omega*a/2) == 0
assert sp.simplify(a*(omega*a-omega*a/2)-omega*a*a/2) == 0
print("PASS fixed-core Taylor lower bound G >= omega*a^2/2.")

# Forced recurrence telescopes with the physical source budget when
# normalized source on stage n is Q_n^alpha F_n.
Qn, qn, G0, budget, increment = sp.symbols("Qn qn G0 budget increment", positive=True)
forced_step = qn**alpha*(Qn**alpha*(G0+budget)+Qn**alpha*increment)
assert sp.simplify(forced_step-(Qn*qn)**alpha*(G0+budget+increment)) == 0
print("PASS finite-time forcing-budget scaling.")
print("The note proves the maximum principle and return implications; this is an algebraic checker.")
