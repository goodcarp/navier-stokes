#!/usr/bin/env python3
"""Exact symbolic checks for cylindrical-NS-representation.md.

No pressure approximation, floating-point quadrature, NS evolution, or
interval certification is performed. This verifies independent algebraic
identities needed by a full cylindrical Fourier implementation.
"""
import json
import sympy as s

r = s.symbols('r', positive=True)
m, k = s.symbols('m k', real=True)
I = s.I
f = s.Function('f')(r)
p = s.Function('p')(r)
up = s.Function('up')(r)
um = s.Function('um')(r)
uz = s.Function('uz')(r)
passed = []


def zero(expr, label):
    answer = s.simplify(s.expand(expr))
    assert answer == 0, (label, answer)
    passed.append(label)


def lap(h, n=m):
    return s.diff(h, r, 2) + s.diff(h, r)/r - (n*n/r**2+k*k)*h


def grad(h):
    return (s.diff(h, r)-m*h/r, s.diff(h, r)+m*h/r, I*k*h)


def div(u):
    a, b, c = u
    return (s.diff(a, r)+(m+1)*a/r+s.diff(b, r)+(1-m)*b/r)/2+I*k*c


def vlap(u):
    return (lap(u[0], m+1), lap(u[1], m-1), lap(u[2], m))


def curl(u):
    a, b, c = u
    return (-k*a-I*(s.diff(c, r)-m*c/r),
            k*b+I*(s.diff(c, r)+m*c/r),
            (s.diff(a, r)+(m+1)*a/r-s.diff(b, r)-(1-m)*b/r)/(2*I))


U = (up, um, uz)
ur = (up+um)/2
ut = (up-um)/(2*I)
lr = lap(ur)-ur/r**2-2*I*m*ut/r**2
lt = lap(ut)-ut/r**2+2*I*m*ur/r**2
zero(lr+I*lt-lap(up, m+1), 'vector Laplacian: plus order m+1')
zero(lr-I*lt-lap(um, m-1), 'vector Laplacian: minus order m-1')
zero(div(grad(p))-lap(p), 'divergence of gradient')
for j, (left, right) in enumerate(zip(vlap(grad(p)), grad(lap(p)))):
    zero(left-right, f'diffusion-gradient intertwining component {j}')
zero(div(vlap(U))-lap(div(U)), 'divergence-diffusion intertwining')
zero(div(curl(U)), 'divergence of curl')
for j, term in enumerate(curl(grad(p))):
    zero(term, f'curl of gradient component {j}')
for j, (left, right, visc) in enumerate(zip(curl(curl(U)), grad(div(U)), vlap(U))):
    zero(left-right+visc, f'curl-curl identity component {j}')

# q stands for the conjugate pressure. This is the exact radial boundary
# identity behind G*=-D, including the helical 1/2 component weights.
q = s.Function('q')(r)
gradient_pair = r*((s.diff(q,r)-m*q/r)*up/2
                  +(s.diff(q,r)+m*q/r)*um/2-I*k*q*uz)
zero(gradient_pair+r*q*div(U)-s.diff(r*q*ur,r),
     'weighted gradient-divergence adjoint boundary identity')

# Each complex field is written in four independent real variables here;
# this independently checks the angular part of the Cartesian energy.
x, y, a, b = s.symbols('x y a b', real=True)
R, T = x+I*y, a+I*b
norm2 = lambda h: s.expand(h*s.conjugate(h))
zero(norm2(I*m*R-T)+norm2(I*m*T+R)
     -((m+1)**2*norm2(R+I*T)+(m-1)**2*norm2(R-I*T))/2,
     'angular gradient energy and shifted orders')
zero((norm2(R+I*T)+norm2(R-I*T))/2-norm2(R)-norm2(T),
     'helical energy weights')

# Differentiate a mode-b advected vector along an independent mode-a
# advector. This checks the basis curvature terms without Poisson algebra.
ab, bb = s.symbols('a_mode b_mode', real=True)
ar, at, az = s.symbols('ar at az')
br, bt, bz = s.symbols('br bt bz')
brr, btr, bzr, brz, btz, bzz = s.symbols('brr btr bzr brz btz bzz')
Nr = ar*brr+I*bb*at*br/r+az*brz-at*bt/r
Nt = ar*btr+I*bb*at*bt/r+az*btz+at*br/r
for sign in (1, -1):
    proposed = (ar*(brr+sign*I*btr)
                +I*(bb+sign)*at*(br+sign*I*bt)/r
                +az*(brz+sign*I*btz))
    zero(Nr+sign*I*Nt-proposed,
         f'nonlinear rotating-basis term sign {sign}')

# The streamfunction coefficient is unrestricted and may contain all
# radial and axial cutoff transitions. No flat-envelope assumption enters.
z = s.symbols('z', real=True)
psi = s.Function('Psi')(r, z)
sr, st = I*m*psi/r, -s.diff(psi,r)
zero(s.diff(sr,r)+sr/r+I*m*st/r,
     'full streamfunction Fourier seed divergence')
zero(sr+I*st+I*(s.diff(psi,r)-m*psi/r),
     'initial seed helical plus')
zero(sr-I*st-I*(s.diff(psi,r)+m*psi/r),
     'initial seed helical minus')

# Check the complete radial affine extension, without fixing its cutoff.
F = s.Function('F')
ss = r*r+z*z
s0 = s.symbols('s0', real=True)
Fp = s.Subs(s.diff(F(s0),s0),s0,ss)
Pr = -r*(F(ss)+2*z*z*Fp)
Pz = 2*z*(F(ss)+r*r*Fp)
zero(s.diff(Pr,r)+Pr/r+s.diff(Pz,z),
     'arbitrary radial affine core and pump divergence')

# A coordinate mapping must include its Jacobian in both operator and norm.
xi = s.symbols('xi', real=True)
RR = s.Function('R')(xi)
JJ = s.diff(RR,xi)
ff = s.Function('a')(xi)
mapped_divergence_form = s.diff(RR/JJ*s.diff(ff,xi),xi)/(RR*JJ)
mapped_chain_rule = (s.diff(ff,xi,2)/JJ**2
                     +(1/(RR*JJ)-s.diff(JJ,xi)/JJ**3)*s.diff(ff,xi))
zero(mapped_divergence_form-mapped_chain_rule,
     'mapped radial conservative operator')

# Factored regularity is checked for representative orders, including both
# helical orders for m=0 and m=4. Polynomial vanishing orders are exact.
for n in (0, 1, 3, 4, 5, 7, 8, 9):
    aa = 2+3*r*r-5*r**4+7*r**6
    zero(lap(r**n*aa,n)-r**n*(s.diff(aa,r,2)
         +(2*n+1)*s.diff(aa,r)/r-k*k*aa),
         f'factored regular scalar operator n={n}')
for mm in (0, 4, 8):
    polynomial = r**mm*(2+3*r*r-5*r**4)
    gp = s.diff(polynomial,r)-mm*polynomial/r
    gm = s.diff(polynomial,r)+mm*polynomial/r
    assert s.Poly(s.expand(gp/r**(mm+1)),r).degree() <= 2
    if mm:
        assert s.Poly(s.expand(gm/r**(mm-1)),r).degree() <= 4
    else:
        zero(gp-gm, 'axisymmetric pressure gradient ladders agree')
    passed.append(f'pressure ladder regularity for m={mm}')

# Free-space k=0 radial Green formula differentiated with independent
# cumulative integrals. Their endpoint derivatives fix the sign and 1/(2n).
n = s.symbols('n', positive=True)
A = s.Function('A')(r)
B = s.Function('B')(r)
P = (r**(-n)*A+r**n*B)/(2*n)
raw = -s.diff(P,r,2)-s.diff(P,r)/r+n*n*P/r**2
replacements = {
    s.diff(A,r,2):s.diff(r**(n+1)*f,r),
    s.diff(B,r,2):s.diff(-r**(1-n)*f,r),
    s.diff(A,r):r**(n+1)*f,
    s.diff(B,r):-r**(1-n)*f,
}
zero(raw.subs(replacements, simultaneous=True)-f,
     'zero-axial-frequency free-space Green normalization')

# Projection algebra follows from DG=L. Verify the associated finite
# dimensional weighted analogue, including a nontrivial helical weight.
H = s.diag(s.Rational(1,2),s.Rational(1,2),1)
G = s.Matrix([[1,0],[0,1],[1,2]])
D = -G.T*H
LL = D*G
PP = s.eye(3)-G*LL.inv()*D
zero((D*PP).norm()**2, 'compatible projection is divergence free')
zero((PP*PP-PP).norm()**2, 'compatible projection is idempotent')
zero((PP.T*H-H*PP).norm()**2, 'compatible projection is weighted self-adjoint')
W = s.Matrix(s.symbols('w0:3',real=True))
V = PP*W
C = W-V
zero((W.T*H*W-V.T*H*V-C.T*H*C)[0],
     'compatible projection Pythagorean identity')

print(json.dumps({'status':'PASS','checks':len(passed),'verified':passed,
    'scope':'Exact symbolic identities only; no discretization, NS evolution, '
            'whole-space truncation bound, or interval certificate.'},indent=2))
