#!/usr/bin/env python3
"""Exact algebra checks for c4-adapted-reference.md, not a fluid simulation.

The moving-coordinate manufactured field is polynomial and need not have
finite energy. It tests the local tensor identities used by the note.
"""
import json
from pathlib import Path
import sympy as s

results = []


def check(label, value):
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    residuals = [s.simplify(s.expand(v)) for v in entries]
    ok = all(v == 0 for v in residuals)
    results.append({'name': label, 'passed': ok})
    if not ok:
        raise AssertionError((label, residuals))


t, r, I = s.symbols('t r I', positive=True)
x = s.Matrix(s.symbols('x0:3'))
y = s.Matrix(s.symbols('y0:3'))
J = s.Matrix([[0,-1,0],[1,0,0],[0,0,0]])
S = s.diag(-1,-1,2)
Q = s.Matrix([[0,-1,0],[1,0,0],[0,0,1]])
b, o, b1, o1 = s.symbols('b o b1 o1', real=True)
A = b*S+o*J
check('central C4 commutant', A*Q-Q*A)
check('central incompressibility', s.trace(A))
check('central matrices at different times commute', A*(b1*S+o1*J)-(b1*S+o1*J)*A)
G = s.Matrix(3,3,s.symbols('g0:9'))
sol = s.linsolve(list(G*Q-Q*G)+[s.trace(G)], list(G))
assert len(sol) == 1
G_sol = s.Matrix(3,3,next(iter(sol)))
check('C4 commutant has only the two central parameters',
      G_sol-((-G_sol[0,0])*S+G_sol[1,0]*J))
results.append({'name':'central commutant dimension is two', 'passed':True})

# Local quartic pressure freedoms are allowed by central C4 symmetry.
z = x[2]
rho2 = x[0]**2+x[1]**2
p4 = x[0]**4-6*x[0]**2*x[1]**2+x[1]**4
p40 = z**4-3*z*z*rho2+s.Rational(3,8)*rho2**2
def lap(f, coords):
    return sum(s.diff(f,a,2) for a in coords)
def subvec(expr, old, new):
    return expr.subs(dict(zip(old,new)), simultaneous=True)
for label,p in [('mode4',p4),('axisymmetric',p40)]:
    check(label+' quartic pressure is harmonic',lap(p,x))
    check(label+' quartic pressure is C4 invariant',subvec(p,x,Q*x)-p)
    check(label+' quartic pressure Hessian vanishes at center',
          s.hessian(p,x).subs(dict.fromkeys(x,0)))

# Fully noncommuting, volume-preserving time-dependent matrix reference.
F = s.Matrix([[1,t,0],[t,1+t*t,0],[0,0,1]])
Fi = F.inv()
At = F.diff(t)*Fi
Bt = Fi*F.diff(t)
C = Fi*Fi.T
check('manufactured determinant one',F.det()-1)
check('manufactured trace A zero',s.trace(At))
assert At.subs(t,0)*At.subs(t,1)-At.subs(t,1)*At.subs(t,0) != s.zeros(3)
results.append({'name':'manufactured reference requires time ordering','passed':True})

# A divergence-free cubic relative field with time dependence.
potential = s.Matrix([y[1]*y[2]**3, t*y[0]**2*y[2]**2, y[0]**2*y[1]**2])
v = s.Matrix([s.diff(potential[2],y[1])-s.diff(potential[1],y[2]),
              s.diff(potential[0],y[2])-s.diff(potential[2],y[0]),
              s.diff(potential[1],y[0])-s.diff(potential[0],y[1])])
Dv = v.jacobian(y)
check('manufactured relative divergence',s.trace(Dv))
u_y = F.diff(t)*y+F*v
u_x = subvec(u_y,y,Fi*x)
grad_u = u_x.jacobian(x)
check('full gradient transform',subvec(grad_u,x,F*y)-(At+F*Dv*Fi))
check('physical divergence',s.trace(grad_u))

acc_x = u_x.diff(t)+grad_u*u_x
acc_y = F.diff(t,2)*y+2*F.diff(t)*v+F*(v.diff(t)+Dv*v)
check('complete material acceleration and factor 2',subvec(acc_x,x,F*y)-acc_y)
Lv = s.Matrix([sum(C[i,j]*s.diff(v[k],y[i],y[j])
                     for i in range(3) for j in range(3)) for k in range(3)])
check('anisotropic vector diffusion',
      subvec(s.Matrix([lap(ui,x) for ui in u_x]),x,F*y)-F*Lv)
check('full pressure source after affine subtraction',
      subvec(s.trace(grad_u*grad_u),x,F*y)-s.trace(At*At)
      -2*s.trace(Bt*Dv)-s.trace(Dv*Dv))

p = x[0]**2*x[1]**2+x[1]**2*x[2]**2+t*x[0]**3*x[2]
py = subvec(p,x,F*y)
check('pressure gradient metric',
      Fi*subvec(s.Matrix([s.diff(p,a) for a in x]),x,F*y)
      -C*s.Matrix([s.diff(py,a) for a in y]))
check('pressure Laplacian metric',subvec(lap(p,x),x,F*y)
      -sum(C[i,j]*s.diff(py,y[i],y[j]) for i in range(3) for j in range(3)))

# Exact comparison majorants and integrated variational exponents.
z3 = r/s.sqrt(1-2*r*r*I)
z2 = r/(1-r*I)
check('cubic radius solves dz/dI=z cubed',s.diff(z3,I)-z3**3)
check('quadratic radius solves dz/dI=z squared',s.diff(z2,I)-z2**2)
check('cubic variational exponent',
      s.diff(-s.Rational(3,2)*s.log(1-2*r*r*I),I)-3*z3*z3)
check('quadratic variational exponent',
      s.diff(-2*s.log(1-r*I),I)-2*z2)

# H4-compatible endpoint modulus and the mixed 5/2-power comparison.
h, radial, tau = s.symbols('h radial tau', positive=True)
radial_bound = h*h*(2/h)+4/(2/h)
check('Fourier endpoint radial bound integral', radial_bound-4*h)
check('endpoint D4 Holder constant squared',
      (2*s.pi)**(-3)*24*16*s.pi-(4*s.sqrt(3)/s.pi)**2)
check('second Taylor Holder factor',
      s.integrate((1-tau)*s.sqrt(tau),(tau,0,1))-s.Rational(4,15))
check('gradient Taylor Holder factor',
      s.integrate(s.sqrt(tau),(tau,0,1))-s.Rational(2,3))
Dm=1-s.Rational(3,2)*r**s.Rational(3,2)*I
zm=r*Dm**(-s.Rational(2,3))
check('mixed radius Bernoulli equation',s.diff(zm,I)-r**s.Rational(5,2)*Dm**(-s.Rational(5,3)))
check('mixed variational exponent',
      s.diff(-2*s.log(Dm),I)-3*r**s.Rational(3,2)/Dm)
xi=s.symbols('xi0:3')
Hpoly=s.Poly((1+sum(q*q for q in xi))**4,*xi)
assert max(Hpoly.coeffs()) == 24
results.append({'name':'Fourier H4 versus D4 coefficient bound 24','passed':True})

# The outer mean maximum has shear, even though its value rotates rapidly.
c, om = s.symbols('c om', real=True)
outer = s.Matrix([[c,-om,0],[-om,c,0],[0,0,-2*c]])
check('outer maximum gradient is symmetric',outer-outer.T)
check('rotating-basis shear generator',
      outer-om*J-s.Matrix([[c,0,0],[-2*om,c,0],[0,0,-2*c]]))
check('outer eigenvector plus',outer*s.Matrix([1,1,0])-(c-om)*s.Matrix([1,1,0]))
check('outer eigenvector minus',outer*s.Matrix([1,-1,0])-(c+om)*s.Matrix([1,-1,0]))

report={'status':'PASS','scope':'exact manufactured/local algebra, no NS trajectory or error enclosure',
        'checks':len(results),'results':results}
Path(__file__).with_name('c4-adapted-reference-check.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS:',len(results),'exact C4, coordinate, pressure, viscosity and comparison identities')
