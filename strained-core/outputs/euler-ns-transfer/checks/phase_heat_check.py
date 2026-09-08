"""Exact checks of the auxiliary phase-heat and return calculations.

These checks verify algebraic examples/assumptions, not an NS solution.
Requires SymPy; no downloaded research code is imported.
"""
import json
import sympy as s

t, theta = s.symbols('t theta', real=True)
checks = []


def zero(label, value):
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(s.simplify(v) == 0 for v in entries), (label, value)
    checks.append(label)


# A noncommuting, explicitly solvable time-dependent component system.
Phi = s.Matrix([[1+t**3, t], [t**2, 1]])
assert Phi.det() == 1
B = Phi.diff(t) * Phi.inv()
assert B.subs(t, 0)*B.subs(t, 1) != B.subs(t, 1)*B.subs(t, 0)
D = t + t**3/3
diffusion = s.diff(D, t)
heat_data = s.Matrix([
    s.exp(-D)*s.sin(theta) + s.exp(-4*D)*s.cos(2*theta),
    2*s.exp(-D)*s.cos(theta) - 3*s.exp(-9*D)*s.sin(3*theta),
])
W = Phi*heat_data
zero('three-mode PDE residual with noncommuting B(t)', W.diff(t)-B*W-diffusion*W.diff(theta,2))
zero('initial data are preserved', W.subs(t,0)-s.Matrix([s.sin(theta)+s.cos(2*theta),2*s.cos(theta)-3*s.sin(3*theta)]))

# Generic scalar quotient, including the diagnostic unequal-damping term.
T, O, a, b, d, d1, d2 = s.symbols('T O a b d d1 d2', nonzero=True)
ratio_rhs = ((b*T-d2*O)*T-O*(a*O-d1*T))/T**2
zero('unequal-damping ratio identity', ratio_rhs-(b+(d1-d2)*O/T-a*(O/T)**2))
zero('equal damping cancels from ratio', ratio_rhs.subs({d1:d,d2:d})-(b-a*(O/T)**2))

# Exact finite pulse: the reset survives and the retained amplitude decreases.
pulse = s.exp(-2*t)*s.Matrix([s.cos(t)+s.sin(t), s.cos(t)-s.sin(t)])
zero('damped pulse solves its ODE', pulse.diff(t)-(s.Matrix([[0,1],[-1,0]])-2*s.eye(2))*pulse)
zero('same exact return time', pulse[1].subs(t,s.pi/4))
zero('retained amplitude includes heat loss',pulse[0].subs(t,s.pi/4)-s.sqrt(2)*s.exp(-s.pi/2))
assert float(pulse[0].subs(t,s.pi/4)) < 1 < float(s.sqrt(2))
checks.append('return-only success would accept a contracting pulse')

# Required negative control: the phase-independence assumption is material.
# For phase-dependent B(theta), B exp(D L) and exp(D L) B need not commute.
g=s.sin(theta)
zero_mode_counterexample = s.diff(s.cos(theta)*g,theta,2)-s.cos(theta)*s.diff(g,theta,2)
assert s.simplify(zero_mode_counterexample) != 0
checks.append('phase-dependent coefficient commutation is correctly rejected')

print(json.dumps({'status':'PASS','scope':'exact auxiliary-model algebra; no PDE blowup certification',
                  'checks':checks,
                  'pulse_retained_amplitude':float(pulse[0].subs(t,s.pi/4)),
                  'undamped_retained_amplitude':float(s.sqrt(2))},indent=2))
