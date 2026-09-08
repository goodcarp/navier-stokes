#!/usr/bin/env python3
"""Check pressure polarization/support algebra; evaluate no profile integrals."""
from itertools import product
import sympy as sp

b, Om, c, A, nu=sp.symbols("b Om c A nu",real=True)
names=("S","W","P","v")
amp={"S":b,"W":Om,"P":c,"v":A}
region={"S":"core","W":"core","P":"outer","v":"outer"}
kind={"S":"pol","W":"azi","P":"pol","v":"azi"}
symbols={}


def pressure_bilinear(i,j,k):
    """Formal Pi(i,B(j,k)), using only the proved vanishing rules.

    Same-parity B outputs are poloidal and may be nonlocal. Mixed-parity
    B outputs are azimuthal, Leray leaves them unchanged, and they are local.
    """
    if region[j] != region[k]:
        return sp.Integer(0)  # Source itself vanishes before Leray.
    output_kind="pol" if kind[j]==kind[k] else "azi"
    if kind[i] != output_kind:
        return sp.Integer(0)
    if output_kind=="azi" and region[i] != region[j]:
        return sp.Integer(0)
    j,k=sorted((j,k),key=names.index)
    key=(i,j,k)
    if key not in symbols:
        symbols[key]=sp.Symbol(f"Pi_{i}_B{j}{k}")
    return symbols[key]


cubic=sp.expand(sum(-2*amp[i]*amp[j]*amp[k]*pressure_bilinear(i,j,k)
                    for i,j,k in product(names,repeat=3)))
expected=(
    b**3*(-2*pressure_bilinear("S","S","S"))
    +b*b*c*(-2*pressure_bilinear("P","S","S"))
    +b*c*c*(-2*pressure_bilinear("S","P","P"))
    +c**3*(-2*pressure_bilinear("P","P","P"))
    +Om*Om*b*(-2*pressure_bilinear("S","W","W")
               -4*pressure_bilinear("W","S","W"))
    +Om*Om*c*(-2*pressure_bilinear("P","W","W"))
    +A*A*b*(-2*pressure_bilinear("S","v","v"))
    +A*A*c*(-2*pressure_bilinear("P","v","v")
             -4*pressure_bilinear("v","P","v"))
)
assert sp.expand(cubic-expected)==0
assert len(sp.Poly(cubic,b,Om,c,A).terms())==8
assert pressure_bilinear("v","S","W")==0
assert pressure_bilinear("W","P","v")==0
assert pressure_bilinear("P","S","S")!=0
assert pressure_bilinear("S","P","P")!=0

visc=0
for i,j in product(names,repeat=2):
    if region[i]==region[j] and kind[i]==kind[j]:
        assert i==j
        visc += 2*nu*amp[i]*amp[j]*sp.Symbol(f"Pi_{i}_Delta{j}")
assert len(sp.Poly(visc,b,Om,c,A,nu).terms())==4
print("PASS: all ordered terms reproduce eight cubic and four viscous coefficients.")
print("PASS: mixed-swirl support cancellations do not discard poloidal Leray tails.")

# Axisymmetric pressure parity, with full cylindrical connection entries.
r=sp.symbols("r",positive=True)
ur,urr,urz,uzr,uzz,w,wr,wz=sp.symbols("ur urr urz uzr uzz w wr wz")
Gp=sp.Matrix([[urr,0,urz],[0,ur/r,0],[uzr,0,uzz]])
Gw=sp.Matrix([[0,-w/r,0],[wr,0,wz],[0,0,0]])
assert sp.trace(Gp*Gw)==0
print("PASS: poloidal/azimuthal pressure cross source vanishes identically.")

# Verify the two compact-core heat-tangent identities, without pressure integrals.
x,y,z=sp.symbols("x y z",real=True)
coords=(x,y,z)
X=sp.Matrix(coords)
radius=x*x+y*y+z*z
p,p1,p2,p3=sp.symbols("p p1 p2 p3")

def D(expr,xi):
    return (sp.diff(expr,xi)+2*xi*(p1*sp.diff(expr,p)
                +p2*sp.diff(expr,p1)+p3*sp.diff(expr,p2)))

def laplace(expr):
    return sp.expand(sum(D(D(expr,xi),xi) for xi in coords))

S0=sp.diag(-1,-1,2)
Q=X.cross(S0*X)
assert Q.applyfunc(lambda e:sum(sp.diff(e,xi,xi) for xi in coords))==sp.zeros(3,1)
eta=4*radius*p2+14*p1
assert (Q.applyfunc(lambda e:laplace(p*e))-eta*Q).applyfunc(sp.expand)==sp.zeros(3,1)
J=sp.Matrix([[0,-1,0],[1,0,0],[0,0,0]])
W=p*(J*X)
heat_swirl=(4*radius*p2+10*p1)*(J*X)
assert (W.applyfunc(laplace)-heat_swirl).applyfunc(sp.expand)==sp.zeros(3,1)
print("PASS: strain and swirl Laplacian identities supporting dS=dW=0.")

t300,t210,t120,t030,tWb,tWc,tvb,tvc=sp.symbols(
    "t300 t210 t120 t030 tWb tWc tvb tvc")
dS,dW,dP,dv=sp.symbols("dS dW dP dv")
mu,Theta,eps,CP,Cv=sp.symbols("mu Theta eps CP Cv")
H=(sp.Rational(38,7)+sp.Rational(2,5)*Theta-CP*mu*mu)/Cv
original=(t300*b**3+t210*b*b*c+t120*b*c*c+t030*c**3
          +Om*Om*(tWb*b+tWc*c)+A*A*(tvb*b+tvc*c)
          +nu*(dS*b*b+dW*Om*Om+dP*c*c+dv*A*A)+32*b**3)
normalized=sp.simplify(original.subs(Om*Om,b*b*Theta).subs(c,b*mu)
                       .subs(A*A,b*b*H).subs(nu,b*eps)/b**3)
target=(32+t300+t210*mu+t120*mu*mu+t030*mu**3
        +Theta*(tWb+tWc*mu)+H*(tvb+tvc*mu)
        +eps*(dS+dW*Theta+dP*mu*mu+dv*H))
assert sp.simplify(normalized-target)==0
polynomial=sp.Poly(sp.expand(normalized),mu,Theta,eps)
assert polynomial.degree(mu)<=3
assert polynomial.degree(Theta)<=1
assert polynomial.degree(eps)<=1
visc_chosen=sp.expand(target.subs({dS:0,dW:0})).coeff(eps)
assert sp.simplify(visc_chosen-(sp.Rational(38,7)*dv/Cv
             +(dP-CP*dv/Cv)*mu*mu+sp.Rational(2,5)*dv*Theta/Cv))==0
print("PASS: neutral normalized gate polynomial and the chosen-core simplification.")
print("All checks passed. No profile integral or remaining gate sign was evaluated.")
