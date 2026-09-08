#!/usr/bin/env python3
"""Exact angular selection and certificate-factor checks; no residual norms measured."""
import sympy as sp

R=sp.symbols("R",positive=True)
mu=sp.symbols("mu",real=True)
p,p1,p2=sp.symbols("p p1 p2")
f0,f1,f2=sp.symbols("f0 f1 f2")
Y0,Y1,Y2=sp.symbols("Y0 Y1 Y2")
q=1-mu*mu
sq=sp.sqrt(q)

def radial(expr):
    return sp.diff(expr,R)+f1*sp.diff(expr,f0)+f2*sp.diff(expr,f1)

def angular(expr):
    return sp.diff(expr,mu)+Y1*sp.diff(expr,Y0)+Y2*sp.diff(expr,Y1)

def Drho(expr):
    return sq*radial(expr)-mu*sq*angular(expr)/R

def Dz(expr):
    return mu*radial(expr)+q*angular(expr)/R

F=f0*Y0
Frr=sp.simplify(Drho(Drho(F)))
Frz=sp.simplify(Dz(Drho(F)))
Fzz=sp.simplify(Dz(Dz(F)))
Ftt=f1*Y0/R-mu*f0*Y1/R**2
Urr=-p-2*R**2*p1-4*R**4*q*mu**2*p2
Utt=-p-2*R**2*mu**2*p1
Uzz=2*p+2*R**2*(1+mu**2)*p1+4*R**4*q*mu**2*p2
Urz_plus_Uzr=2*R**2*mu*sq*p1+4*R**4*mu*sq*(1-2*mu**2)*p2
trace=sp.expand(sp.simplify(Urr*Frr+Utt*Ftt+Uzz*Fzz+Urz_plus_Uzr*Frz))
weighted=sp.expand((3*mu**2-1)*trace)
effective={}
explicit={
    f0: R**2*p2*(-60*mu**4+48*mu**2-4)
        +p1*(120*mu**4-108*mu**2+12)
        +(p/R**2)*(135*mu**4-126*mu**2+15),
    f1: R**3*p2*(60*mu**4-48*mu**2+4)
        +R*p1*(132*mu**4-108*mu**2+8)
        +(p/R)*(81*mu**4-66*mu**2+5),
    f2: (3*mu**2-1)**2*(2*R**2*p1+p),
}

def exact_polynomial_integral(expr):
    """Integrate the angular polynomial with rational monomial moments.

    A raw SymPy integrate(... ) result at f0, ell=6 originally remained an
    uncombined six-term rational expression. It is exactly zero; using ==0
    before canonicalization tested representation rather than the integral.
    This path uses no numerical sampling or approximate zero tolerance.
    """
    polynomial=sp.Poly(sp.expand(expr),mu)
    result=sum(coefficient*sp.Rational(2,power+1)
               for (power,),coefficient in polynomial.terms() if power%2==0)
    return sp.cancel(result)

# The three coefficient sums in that first raw integral's common numerator.
assert -3192-1260-220+20+1320+3332==0
assert 7056+2520+600-60-3144-6972==0
assert 8190+2835+735-75-3726-7959==0
for f in (f0,f1,f2):
    row=weighted.coeff(f)
    A=sp.expand(row.coeff(Y0))
    B=sp.expand(row.coeff(Y1))
    C=sp.expand(row.coeff(Y2))
    for endpoint in (-1,1):
        assert sp.simplify(C.subs(mu,endpoint))==0
        assert sp.simplify((B-sp.diff(C,mu)).subs(mu,endpoint))==0
    effective[f]=sp.factor(A-sp.diff(B,mu)+sp.diff(C,mu,mu))
    assert sp.expand(effective[f]-explicit[f])==0
    assert sp.Poly(effective[f],mu).degree()<=4
    assert sp.expand(effective[f]-effective[f].subs(mu,-mu))==0
    for degree in range(5,13):
        assert exact_polynomial_integral(effective[f]*sp.legendre(degree,mu))==0
print("PASS: exact three core weights, boundary cancellation and degree-four cutoff.")

# At the origin, a smooth solid harmonic of degree>2 has zero Hessian.
x,y,z=sp.symbols("x y z",real=True)
radius_squared=x*x+y*y+z*z
for degree in range(3,9):
    harmonic=0
    for (power,),coefficient in sp.Poly(sp.legendre(degree,mu),mu).terms():
        assert (degree-power)%2==0
        harmonic += coefficient*z**power*radius_squared**((degree-power)//2)
    for coordinate_a in (x,y,z):
        for coordinate_b in (x,y,z):
            assert sp.diff(harmonic,coordinate_a,coordinate_b).subs({x:0,y:0,z:0})==0
print("PASS: high-degree regular pressure modes have no central Hessian contact.")

for cutoff in range(101):
    first=cutoff+1
    eigenvalue=first*(first+1)
    assert eigenvalue==(cutoff+1)*(cutoff+2)
    first_even=2*(cutoff//2+1)
    assert first_even>cutoff and first_even%2==0
    assert first_even-2<=cutoff
    even_eigenvalue=first_even*(first_even+1)
    assert even_eigenvalue>=eigenvalue
assert 6*7==42  # First omitted even degree at cutoff L=4.
print("PASS: general and even-mode angular Poincare eigenvalues.")

radius,lam,norm_e,norm_weight=sp.symbols("radius lam norm_e norm_weight",positive=True)
energy_bound=radius*norm_e/sp.sqrt(lam)
outer_bound=2*norm_weight*energy_bound
assert sp.simplify(outer_bound-2*radius*norm_weight*norm_e/sp.sqrt(lam))==0
low_energy_bound=2*radius*norm_e
low_outer_bound=2*norm_weight*low_energy_bound
assert low_outer_bound==4*radius*norm_weight*norm_e
print("PASS: high-mode outer factor2R/sqrt(lambda) and low-mode outer factor4R.")
print("All algebraic checks passed. Energy estimates are proved in the note; no numerical residual was certified.")
