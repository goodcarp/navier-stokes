#!/usr/bin/env python3
"""Exact identities behind the independent initial-evolution audit."""
import json
import sympy as s

x=s.symbols('x y z',real=True)
u=s.Matrix([s.Function(f'u{i}')(*x) for i in range(3)])
v=s.Matrix([s.Function(f'v{i}')(*x) for i in range(3)])
Gu=u.jacobian(x);Gv=v.jacobian(x)
C=Gu*v+Gv*u
div=lambda w:sum(s.diff(w[i],x[i]) for i in range(3))
grad=lambda f:s.Matrix([s.diff(f,xx) for xx in x])
assert s.simplify(div(C)-2*s.trace(Gu*Gv)
                  -(u.dot(grad(div(v))))-v.dot(grad(div(u))))==0

b,o,b1,o1,d,h,hz=s.symbols('b Omega b1 Omega1 d h hz',real=True)
L=s.Matrix([[-b,-o,0],[o,-b,0],[0,0,2*b]])
L1=s.Matrix([[-b1+d/2,-o1,0],[o1,-b1+d/2,0],[0,0,2*b1]])
L2=-(L*L1+L1*L)-s.diag(h,h,hz)
b2=L2[2,2]/2
o2=(L2[1,0]-L2[0,1])/2
assert s.simplify(b2+4*b*b1+hz/2)==0
assert s.simplify(o2-(2*b*o1+2*b1*o-d*o))==0
beta2=b2/o-b*o2/o**2-2*b1*o1/o**2+2*b*o1**2/o**3
beta1=s.symbols('beta1',real=True)
affine=s.simplify(beta2.subs({b:1,o:1,o1:2,b1:beta1+2}))
assert s.simplify(affine-(-(hz+32)/2-10*beta1+d))==0

t=s.symbols('t',real=True)
B0,B1,B2,O0,O1,O2=s.symbols('B0 B1 B2 O0 O1 O2',real=True,nonzero=True)
q=(B0+B1*t+B2*t*t/2)/(O0+O1*t+O2*t*t/2)
assert s.simplify(s.diff(q,t).subs(t,0)-(B1/O0-B0*O1/O0**2))==0
assert s.simplify(s.diff(q,t,2).subs(t,0)
       -(B2/O0-B0*O2/O0**2-2*B1*O1/O0**2+2*B0*O1**2/O0**3))==0

m,rate=s.symbols('m rate',real=True)
U=s.Function('U')(t)
rotated=s.exp(s.I*m*rate*t)*U
assert s.simplify(s.diff(rotated,t).subs(t,0)
                  -s.diff(U,t).subs(t,0)-s.I*m*rate*U.subs(t,0))==0
assert s.simplify(s.diff(rotated,t,2).subs(t,0)
                  -s.diff(U,t,2).subs(t,0)
                  -2*s.I*m*rate*s.diff(U,t).subs(t,0)
                  +(m*rate)**2*U.subs(t,0))==0

print(json.dumps({'status':'PASS','checks':8,'scope':'Exact Cartesian '
    'divergence, affine-core quotient and coordinate-rotation identities; '
    'no numerical pressure or NS evolution certification.'},indent=2))
