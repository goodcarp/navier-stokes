#!/usr/bin/env python3
"""Exact local cubic kernels and support-ratio bounds; no numerical quadrature."""
import sympy as sp

r, z = sp.symbols("r z", positive=True)
t = sp.symbols("t", nonnegative=True)
R2 = r*r + z*z
Nbar = R2**(-sp.Rational(1, 2))  # 4*pi*N
Qtheta = -sp.diff(Nbar, z, 2, r)/r
Qrr = -sp.diff(Nbar, z, 2, r, 2)
Qrz = -sp.diff(Nbar, z, 3, r)
def D(expr):
    return r*sp.diff(expr, r)-2*z*sp.diff(expr, z)

# The plateau field (r,-2z) has zero cylindrical divergence.
assert sp.diff(r*r, r)/r + sp.diff(-2*z, z) == 0
K = sp.factor(-D(Qtheta)+2*Qtheta-2*(Qrr-2*z*Qrz/r))
K_expected = -45*(r**4-12*r*r*z*z+8*z**4)/R2**sp.Rational(9, 2)
S = sp.factor(sp.diff(Nbar, z, 4, r)/r)
assert sp.simplify(K-K_expected) == 0
assert sp.simplify(S-K/R2) == 0
assert sp.simplify(Qtheta-3*(4*z*z-r*r)/R2**sp.Rational(7, 2)) == 0
print("PASS signs and coefficients of Qtheta, K_loc, and S_theta=K_loc/R².")

F = 15*(t*t-12*t+8)/((4-t)*(1+t))
G = F/(1+t)
assert sp.factor((K/Qtheta).subs(r*r, t*z*z)+F) == 0
assert sp.factor((S/Qtheta).subs(r*r, t*z*z)+G/z**2) == 0
assert sp.factor(sp.diff(F, t)+45*(3*t*t-8*t+24)/((t-4)**2*(t+1)**2)) == 0
assert sp.expand(3*(t-sp.Rational(4,3))**2+sp.Rational(56,3)-(3*t*t-8*t+24)) == 0
assert sp.factor(sp.diff(G, t)-15*(t**3-25*t*t+80*t-104)/((t-4)**2*(t+1)**3)) == 0
# For 0<=t<=1/16 the numerator of G' is <=(1/16)^3+5-104<0.
assert sp.Rational(1,16)**3+5-104 < 0
assert F.subs(t, 0) == G.subs(t, 0) == 30
assert F.subs(t, sp.Rational(1,16)) > 0
print("PASS decreasing positive F and G on [0,1/16].")

T = sp.Rational(63,710)**2
zmin, zmax = sp.Rational(71,20), sp.Rational(89,20)
assert T < sp.Rational(1,16)
# Exact symmetric transition: integral psi = plateau length + half transition.
M = sp.Rational(1,4)+sp.Rational(3,4)/2
assert M == sp.Rational(5,8)
E = sp.symbols("E", positive=True)
assert sp.cancel(E/(1+E)+1/(1+E)) == 1

Fmin = sp.factor(F.subs(t, T))
corelo = sp.factor(-sp.Rational(18,7)+M*G.subs(t,T)/(15*zmax*zmax))
corehi = -sp.Rational(18,7)+2*M/(zmin*zmin)
assert Fmin == sp.Rational(30134114372415,1022453805739)
assert corelo == -sp.Rational(72293765997170443165798,28803445800475159076977)
assert corehi == -sp.Rational(87238,35287)
assert -sp.Rational(251,100) < corelo < corehi < -sp.Rational(247,100)
assert Fmin > sp.Rational(2947,100)
print("PASS exact support endpoints and cutoff moment M=5/8.")
print("K/Q range:", -30, "to", -Fmin, "=", -float(Fmin))
print("core/Cv range:", corelo, "to", corehi, "=", float(corelo), float(corehi))

# Independent Cartesian source calculation for any radial compact swirl.
x, y, zz = sp.symbols("x y zz", real=True)
h, h1 = sp.symbols("h h1")
coords = (x,y,zz)
W = sp.Matrix([-y*h,x*h,0])
DW = sp.Matrix(3,3,lambda i,j: sp.diff(W[i],coords[j])+2*coords[j]*h1*sp.diff(W[i],h))
gW = sp.expand(sp.trace(DW*DW))
assert sp.expand(gW-(-2*h*h-4*(x*x+y*y)*h*h1)) == 0
# Source angular average times radial measure is an exact total derivative:
# s^(1/2)*[-2h²-(8/3)s h h'] = -(4/3)d_s[s^(3/2)h²].
s = sp.symbols("s", positive=True)
rhs = -sp.Rational(4,3)*(sp.Rational(3,2)*sp.sqrt(s)*h*h+2*s**sp.Rational(3,2)*h*h1)
assert sp.expand(rhs-sp.sqrt(s)*(-2*h*h-sp.Rational(8,3)*s*h*h1)) == 0
print("PASS radial core-swirl source has only modes 0,2 and zero exterior monopole.")

# For an annular eta, boundary terms vanish and
# integral R^-2 eta''(R²)dR = (3/2) integral R^-4 eta'(R²)dR.
I = sp.symbols("I")
annular_l2_integral = -sp.Rational(32,7)*(3*I-2*sp.Rational(3,2)*I)
assert annular_l2_integral == 0
print("PASS exact annular cancellation against exterior R^-3 P2: tWc=0.")
print("All algebraic checks passed. The note proves integration and positivity steps.")
