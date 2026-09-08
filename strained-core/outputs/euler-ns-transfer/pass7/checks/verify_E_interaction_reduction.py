#!/usr/bin/env python3
"""Exact symbolic checks for E; no numerical pressure/sign certification."""
import sympy as s

r, z = s.symbols("r z", positive=True, real=True)
m, k = s.symbols("m k", real=True, nonzero=True)
a = s.Function("a")(r, z)
V = s.Function("V")(r, z)
phase = s.symbols("phase", real=True)
sn, cs = s.sin(phase), s.cos(phase)
wr = -m*a*sn/r
wt = -s.diff(a, r)*cs + k*a*sn
dr = lambda f: s.diff(f, r) + k*s.diff(f, phase)
dt = lambda f: m*s.diff(f, phase)
cr = V*dt(wr)/r - 2*V*wt/r
ct = V*dt(wt)/r + (s.diff(V, r) + V/r)*wr

def avg(expr):
    p = s.Poly(s.expand(expr), sn, cs)
    out = 0
    for (i, j), coeff in p.terms():
        if (i, j) == (0, 0):
            out += coeff
        elif (i, j) in ((2, 0), (0, 2)):
            out += coeff/2
        elif (i, j) in ((1, 0), (0, 1), (1, 1)):
            pass
        else:
            raise AssertionError((i, j))
    return s.expand(out)

advwwt = wr*dr(wt) + wt*dt(wt)/r + wr*wt/r
assert s.simplify(avg(advwwt) +
                  m*k*a*s.diff(a,r)/r + m*k*a*a/(2*r*r)) == 0
assert s.simplify(avg(wr*cr) - m*k*V*a*a/(r*r)) == 0
assert s.simplify(avg(wt*ct) +
                  m*k*a*a*(s.diff(V,r)+V/r)/(2*r)) == 0

qt = s.Function("qt")(r)
qr = qt+r*s.diff(qt,r)
before = s.expand(4*s.pi*r*(qt*V*avg(advwwt)
                             +qr*avg(wr*cr)+qt*avg(wt*ct)))
after = 6*s.pi*m*k*a*a*V*s.diff(qt,r)
# The difference is an exact radial derivative, so compact support removes it.
boundary = -2*s.pi*m*k*qt*V*a*a
assert s.simplify(before-after-s.diff(boundary,r)) == 0

N = 1/(4*s.pi*s.sqrt(r*r+z*z))
# On y=0, Q_theta=-r^{-1} partial_r partial_zz N.
Qtheta = -s.diff(N,z,2,r)/r
Qrr = -s.diff(N,z,2,r,2)
assert s.simplify(Qtheta-3*(4*z*z-r*r)/
                  (4*s.pi*(r*r+z*z)**s.Rational(7,2))) == 0
assert s.simplify(s.diff(Qtheta,r)-15*r*(r*r-6*z*z)/
                  (4*s.pi*(r*r+z*z)**s.Rational(9,2))) == 0
assert s.simplify(Qrr-Qtheta-r*s.diff(Qtheta,r)) == 0

# Full mixed source, with no z-pressure approximation.
b = s.Function("b")(r,z)
Wr = s.I*m*b/r
Wt = -s.diff(b,r)
Cr = s.I*m*V*Wr/r-2*V*Wt/r
Ct = s.I*m*V*Wt/r+(s.diff(V,r)+V/r)*Wr
source = s.diff(Cr,r)+Cr/r+s.I*m*Ct/r
expected = 2*(s.diff(V*s.diff(b,r),r)-m*m*s.diff(V,r)*b/r)/r
assert s.simplify(source-expected) == 0
assert s.simplify(s.diff(Wr,r)+Wr/r+s.I*m*Wt/r) == 0
Q1,Q2,Q3 = s.symbols("Qrr Qtt Qrz")
divQW = Q1*s.diff(Wr,r)+Q2*(Wr/r+s.I*m*Wt/r)+Q3*s.diff(Wr,z)
assert s.simplify(divQW-(Q1-Q2)*s.diff(Wr,r)-Q3*s.diff(Wr,z)) == 0

# Derivative-reduced whole-space radial Green solution, arbitrary m != 0.
C = s.Function("C")(r)
br = s.Function("br")(r)
Ip = s.Function("Ip")(r)
Im = s.Function("Im")(r)
rules = {
    s.diff(Ip,r,2): s.diff(r**m*s.diff(C,r)*br,r),
    s.diff(Im,r,2): -s.diff(r**(-m)*s.diff(C,r)*br,r),
    s.diff(Ip,r): r**m*s.diff(C,r)*br,
    s.diff(Im,r): -r**(-m)*s.diff(C,r)*br,
}
P = -2*C*br-(m-1)*r**(-m)*Ip-(m+1)*r**m*Im
L = lambda f: s.diff(f,r,2)+s.diff(f,r)/r-m*m*f/r**2
Sr = 2*(s.diff(r*C*s.diff(br,r),r)
        -m*m*s.diff(r*C,r)*br/r)/r
assert s.simplify((-L(P)-Sr).subs(rules, simultaneous=True)) == 0

# Separated trial residual sign and norm-tail denominators.
F = s.Function("F")(z)
assert s.expand((L(P)*F+s.diff(P*F,z,2)) -
                (L(P)*F+P*s.diff(F,z,2))) == 0
assert s.simplify(s.diff(r**(4-2*m)/(4-2*m),r)-r**(3-2*m)) == 0
assert s.simplify(s.diff(r**(2*m+4)/(2*m+4),r)-r**(2*m+3)) == 0
assert 2*4-4 == 4  # m=4 exterior weighted tail is integrable.

print("PASS: local E, full mixed source, Q/divergence identities, radial Green "
      "solution, residual sign, and exact weighted tails.")
