#!/usr/bin/env python3
"""
d1 -- the algebra of the doubling clock, for the Eulerian/Duhamel route.

Part A. SYMBOLIC verification of the three identities the whole argument rests on,
        derived from 3D Navier-Stokes in Cartesian coordinates with an axisymmetric
        no-swirl ansatz (nothing quoted, everything re-derived):
          (I1)  D_t u^r = -d_r p + nu (d_rr + (1/r) d_r - 1/r^2 + d_zz) u^r
          (I2)  D_t eta = nu (d_rr + (3/r) d_r + d_zz) eta ,  eta = omega^theta / r
          (I3)  D_t a   = -a^2 + beta + nu V ,   a = u^r/r,
                beta := -(1/r) d_r p ,  V := (1/r)(d_rr + (1/r)d_r - 1/r^2 + d_zz) u^r
        (I3) is (I1) divided by r, using D_t(1/r) = -a/r.  The a^2 compounding term is
        NOT a free gain: it is exactly the -(u^r)^2/r^2 produced by differentiating 1/r.

Part B. The three clock routes, with their exact constants.
        theta := a0 * T_d ; c2 := M T_d log Re_E -> (2/kappa) theta  (kappa = strain per e-fold).
          R0  frozen strain          theta = log(3/2)        [NOT a lower-bound route]
          R1  Taylor/Duhamel, |beta| <= C a0^2
          R2  Riccati,        beta   >= -q^2 a0^2
        plus the viscous correction (1+theta) e^{-eps theta} = 3/2, eps = lambda/a0.

Every number printed is computed here.  Numerics falsify, never prove.
"""
import json, math
import numpy as np
import sympy as sp

OUT = {}
LOG = []
def say(s=""):
    print(s); LOG.append(s)

# ======================================================================= PART A
say("="*78)
say("PART A -- symbolic re-derivation of the identities (sympy)")
say("="*78)

x, y, z, t, nu = sp.symbols('x y z t nu', real=True, positive=False)
r = sp.sqrt(x**2 + y**2)
ur = sp.Function('ur')(r, z, t)     # u^r(r,z,t)
uz = sp.Function('uz')(r, z, t)     # u^z(r,z,t)
pp = sp.Function('p')(r, z, t)

# Cartesian components of the axisymmetric no-swirl field
U = sp.Matrix([ur*x/r, ur*y/r, uz])
X = (x, y, z)

def lap(f):
    return sum(sp.diff(f, v, 2) for v in X)

# momentum residual  N = d_t U + (U.grad)U + grad p - nu Lap U
adv = sp.Matrix([sum(U[k]*sp.diff(U[i], X[k]) for k in range(3)) for i in range(3)])
gp  = sp.Matrix([sp.diff(pp, v) for v in X])
NS  = sp.Matrix([sp.diff(U[i], t) for i in range(3)]) + adv + gp - nu*sp.Matrix([lap(U[i]) for i in range(3)])

# project on e_r = (x/r, y/r, 0)
er = sp.Matrix([x/r, y/r, 0])
res_r = sp.simplify((NS.T*er)[0, 0])

R, Z, T = sp.symbols('R Z T', positive=True)
def to_rz(expr):
    """evaluate a (r,z)-symmetric expression at the point (x,y)=(R,0), z=Z, t=T"""
    e = expr.subs({x: R, y: 0, z: Z, t: T})
    return sp.simplify(sp.expand(e))

# the claimed identity (I1), same substitution
urF = sp.Function('ur'); uzF = sp.Function('uz'); pF = sp.Function('p')
UR = urF(R, Z, T); UZ = uzF(R, Z, T); P = pF(R, Z, T)
I1 = (sp.diff(UR, T) + UR*sp.diff(UR, R) + UZ*sp.diff(UR, Z)
      + sp.diff(P, R)
      - nu*(sp.diff(UR, R, 2) + sp.diff(UR, R)/R - UR/R**2 + sp.diff(UR, Z, 2)))
d1 = sp.simplify(to_rz(res_r) - I1)
say(f"(I1) residual  e_r . [NS] - claimed  =  {d1}")
OUT['I1_residual'] = str(d1)
assert d1 == 0, "I1 failed"

# vorticity: omega = curl U ; only the theta component survives
curl = sp.Matrix([sp.diff(U[2], y) - sp.diff(U[1], z),
                  sp.diff(U[0], z) - sp.diff(U[2], x),
                  sp.diff(U[1], x) - sp.diff(U[0], y)])
etheta = sp.Matrix([-y/r, x/r, 0])
om_th = sp.simplify(to_rz((curl.T*etheta)[0, 0]))
say(f"omega^theta = {sp.simplify(om_th)}   (expected d_z u^r - d_r u^z)")
OUT['omega_theta'] = str(sp.simplify(om_th))
assert sp.simplify(om_th - (sp.diff(UR, Z) - sp.diff(UZ, R))) == 0

# (I2): build the eta equation from the vorticity equation, symbolically in (r,z)
# vorticity eq for the scalar w = omega^theta :
#   d_t w + u^r d_r w + u^z d_z w = (u^r/r) w + nu (d_rr + (1/r)d_r - 1/r^2 + d_zz) w
# We VERIFY this from the Cartesian vorticity equation instead of quoting it.
Om = curl
vort = (sp.Matrix([sp.diff(Om[i], t) for i in range(3)])
        + sp.Matrix([sum(U[k]*sp.diff(Om[i], X[k]) for k in range(3)) for i in range(3)])
        - sp.Matrix([sum(Om[k]*sp.diff(U[i], X[k]) for k in range(3)) for i in range(3)])
        - nu*sp.Matrix([lap(Om[i]) for i in range(3)]))
res_th = sp.simplify(to_rz((vort.T*etheta)[0, 0]))
W = sp.diff(UR, Z) - sp.diff(UZ, R)
claim_w = (sp.diff(W, T) + UR*sp.diff(W, R) + UZ*sp.diff(W, Z) - (UR/R)*W
           - nu*(sp.diff(W, R, 2) + sp.diff(W, R)/R - W/R**2 + sp.diff(W, Z, 2)))
d2 = sp.simplify(sp.expand(res_th - claim_w))
say(f"(I2a) residual  e_theta . [vorticity eq] - claimed w-equation  =  {d2}")
OUT['I2a_residual'] = str(d2)
assert d2 == 0, "I2a failed"

# now eta = w/r : substitute w = r*eta and check the 5D Laplacian pops out
ETA = sp.Function('eta')(R, Z, T)
sub_w = {W: R*ETA}
lhs = (sp.diff(R*ETA, T) + UR*sp.diff(R*ETA, R) + UZ*sp.diff(R*ETA, Z) - (UR/R)*(R*ETA)
       - nu*(sp.diff(R*ETA, R, 2) + sp.diff(R*ETA, R)/R - (R*ETA)/R**2 + sp.diff(R*ETA, Z, 2)))
claim_eta = R*(sp.diff(ETA, T) + UR*sp.diff(ETA, R) + UZ*sp.diff(ETA, Z)
               - nu*(sp.diff(ETA, R, 2) + 3*sp.diff(ETA, R)/R + sp.diff(ETA, Z, 2)))
d3 = sp.simplify(sp.expand(lhs - claim_eta))
say(f"(I2b) residual  [w-eq with w=r eta] - r*[eta-eq with 5D Laplacian]  =  {d3}")
OUT['I2b_residual'] = str(d3)
assert d3 == 0, "I2b failed"

# (I3): D_t a = D_t(u^r)/r - a^2
A = UR/R
Dt = lambda f: sp.diff(f, T) + UR*sp.diff(f, R) + UZ*sp.diff(f, Z)
d4 = sp.simplify(Dt(A) - (Dt(UR)/R - A**2))
say(f"(I3a) residual  D_t(u^r/r) - [D_t(u^r)/r - a^2]  =  {d4}")
OUT['I3a_residual'] = str(d4)
assert d4 == 0, "I3a failed"

BETA = -sp.diff(P, R)/R
VISC = (nu/R)*(sp.diff(UR, R, 2) + sp.diff(UR, R)/R - UR/R**2 + sp.diff(UR, Z, 2))
d5 = sp.simplify(Dt(A) - (-A**2 + BETA + VISC))   # using I1
d5 = sp.simplify(d5.subs(sp.diff(UR, T),
        sp.solve(sp.Eq(I1, 0), sp.diff(UR, T))[0]))
say(f"(I3b) residual  D_t a - [-a^2 + beta + nu V]  (after using I1)  =  {d5}")
OUT['I3b_residual'] = str(d5)
assert d5 == 0, "I3b failed"
say("ALL SYMBOLIC IDENTITIES VERIFIED (residual identically 0).")
say()

# ======================================================================= PART B
say("="*78)
say("PART B -- the three clock routes and their exact constants")
say("="*78)
say("Notation: a0 = a(x0,0), theta = a0 T_d, and asymptotically")
say("          c2 := M T_d log Re_E = (2/kappa) theta   (kappa = strain per e-fold;")
say("          kappa = 1/2 for the bang-bang plateau, so c2 = 4 theta).")
say()

G = 1.5                                    # the (3/2) threshold

# ---- R0 : frozen strain (recorded, NOT a lower-bound route) -------------------
th0 = math.log(G)
say(f"R0 frozen strain a(t) == a0 :   theta = log(3/2) = {th0:.7f}   c2 = {4*th0:.7f}")
say("   *** R0 is NOT a valid route for an UPPER bound on T_d: a(t) is only constant if")
say("       beta == a^2 exactly.  beta < a^2 makes a decay, so exp(a0 t) OVERSTATES the")
say("       growth and log(3/2)/a0 UNDERSTATES T_d.  The frame's c2 = 4 log(3/2) = 1.62186")
say("       is the value of an idealisation, not of a bound.")
say()

# ---- R1 : Taylor / Duhamel with |beta| <= C a0^2 -----------------------------
# omega(tau) = omega0 + tau*a0*omega0 + int_0^tau (tau-s) beta omega ds
# |omega| <= (3/2)M up to the first crossing  =>  omega(tau) >= M(1 + x - (3/4)C x^2), x = a0 tau
def R1_theta(C):
    if C <= 0: return 0.5
    disc = 1.0 - 1.5*C
    if disc < 0: return None
    return (1.0 - math.sqrt(disc))/(1.5*C)
Cmax = 2.0/3.0
say("R1 Taylor/Duhamel, hypothesis  sup|beta| <= C a0^2  on the window:")
say("   need  1 + x - (3/4) C x^2 >= 3/2 ,  x = a0 tau .  Solvable iff C <= 2/3.")
say(f"   C_max = {Cmax:.7f}  (at C = C_max the only root is x = 2/(3C) = {2/(3*Cmax):.4f}, c2 = {4*2/(3*Cmax):.4f})")
rowsR1 = []
for C in (0.0, 0.05, 0.1, 0.2, 1/3, 0.5, 0.6, 2/3, 0.7, 1.0, 4.0):
    th = R1_theta(C)
    rowsR1.append(dict(C=C, theta=th, c2=(None if th is None else 4*th)))
    say(f"   C = {C:<6.4f}  theta = {'FAILS' if th is None else f'{th:.6f}'}"
        f"   c2 = {'--' if th is None else f'{4*th:.6f}'}")
OUT['R1'] = rowsR1
say()
say("   TRANSLATION INTO THE FRAME'S UNITS.  a0 = (M/2)L, L = log(R/rho0), so a0^2 = M^2L^2/4,")
say("   and |d_t^2 omega| = |beta||omega| <= C a0^2 (3/2)M = (3/8) C M^3 L^2.")
say(f"   C <= 2/3  <=>  |d_t^2 omega| <= {3/8*2/3:.4f} M^3 log^2(R/rho0)  on the window.")
say("   A bound |d_t^2 omega| <~ M^3 log^2 with an O(1) constant is therefore NOT enough:")
say("   at the constant 1 (C = 4) the O(tau^2) remainder is 1.5 x the O(tau) gain already at")
say("   x = 1/2, and no x reaches 3/2 at all; the route closes only below the constant 1/4,")
say("   i.e. C = 4 is 6 x the largest admissible C = 2/3.")
OUT['R1_second_derivative_threshold'] = 3/8*2/3

# remainder / gain at the R0 window, as a function of C
say()
say("   remainder/gain at x = theta:   (3/4) C x^2 / x = (3/4) C x")
for C in (0.25, 0.5, 2/3, 1.0, 2.0, 4.0):
    x = 0.5
    say(f"     C = {C:<5.3f}  at x = 1/2 :  remainder/gain = {0.75*C*x:.4f}")
say()

# ---- R2 : Riccati with beta >= -q^2 a0^2 -------------------------------------
# a' >= -a^2 - B, B = q^2 a0^2  =>  a(t) >= sqrt(B) tan(theta0 - sqrt(B) t), theta0=arctan(a0/sqrt B)
# int_0^t a >= log( cos(theta0 - sqrt(B) t) / cos theta0 )
def R2_theta(q):
    if q == 0.0:
        return 0.5                      # log(1+x) = log(3/2) -> x = 1/2
    th0_ = math.atan(1.0/q)
    c0 = math.cos(th0_)
    if G*c0 >= 1.0:
        return None
    thf = math.acos(G*c0)
    return (th0_ - thf)/q
qmax = math.sqrt(4.0/5.0)
say("R2 Riccati, hypothesis  beta + nu V >= -q^2 a0^2  along the trajectory:")
say("   a' >= -a^2 - B  =>  int_0^t a >= log[cos(th0 - sqrt(B) t)/cos th0],  th0 = arctan(1/q)")
say(f"   solvable iff cos th0 <= 2/3, i.e. q <= sqrt(4/5) = {qmax:.7f}")
rowsR2 = []
for q in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.894, qmax, 0.9, 1.0):
    th = R2_theta(q)
    rowsR2.append(dict(q=q, B_over_a0sq=q*q, theta=th, c2=(None if th is None else 4*th)))
    say(f"   q = {q:<6.4f} (B/a0^2 = {q*q:<7.5f})  theta = {'FAILS' if th is None else f'{th:.6f}'}"
        f"   c2 = {'--' if th is None else f'{4*th:.6f}'}")
OUT['R2'] = rowsR2
say()
say(f"   *** q = 0 (i.e. beta + nu V >= 0) gives theta = 1/2 exactly and  c2 = 2  exactly. ***")
say()

# numerical ODE cross-check of the R2 closed form
def ode_check(q, n=400001):
    B = q*q
    tmax = R2_theta(q)
    if tmax is None: return None
    ts = np.linspace(0.0, tmax, n)
    dt = ts[1]-ts[0]
    a = 1.0; I = 0.0                     # a0 = 1 units
    for k in range(n-1):                 # RK4 on a' = -a^2 - B, and I' = a
        def f(a_): return -a_*a_ - B
        k1=f(a); k2=f(a+0.5*dt*k1); k3=f(a+0.5*dt*k2); k4=f(a+dt*k3)
        an = a + dt/6*(k1+2*k2+2*k3+k4)
        I += dt/6*(a + 4*(a+0.5*dt*k1*0 + 0.5*(a+an-a)) + an) if False else 0.5*dt*(a+an)
        a = an
    return I
say("   ODE cross-check of the closed form (trapezoid on RK4, a0 = 1): int_0^theta a should be log(3/2)")
chk = []
for q in (0.0, 0.3, 0.6, 0.8):
    I = ode_check(q, 200001)
    chk.append(dict(q=q, I=I, target=math.log(G), relerr=abs(I-math.log(G))/math.log(G)))
    say(f"     q = {q:.2f}  int a = {I:.9f}  target {math.log(G):.9f}  rel {abs(I-math.log(G))/math.log(G):.2e}")
OUT['R2_ode_check'] = chk
assert all(c['relerr'] < 2e-8 for c in chk), "R2 closed form failed its own ODE check"
say()

# ---- viscous correction ------------------------------------------------------
# eta loses at most e^{-lambda t} along the trajectory (H3), so
#   |omega(X(t))| >= M (1 + a0 t) e^{-lambda t}   under beta + nuV >= 0
from scipy.optimize import brentq
def theta_visc(eps):
    f = lambda x: (1.0 + x)*math.exp(-eps*x) - G
    hi = 1.0
    while f(hi) < 0 and hi < 1e6: hi *= 2
    if f(hi) < 0: return None
    return brentq(f, 1e-12, hi, xtol=1e-14, rtol=1e-15)
say("Viscous correction (R2 + eta-loss e^{-lambda t}):  (1+theta) e^{-eps theta} = 3/2, eps = lambda/a0")
rowsV = []
for eps in (0.0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6):
    th = theta_visc(eps)
    rowsV.append(dict(eps=eps, theta=th, c2=(None if th is None else 4*th)))
    say(f"   eps = {eps:<5.3f}  theta = {'FAILS' if th is None else f'{th:.6f}'}"
        f"   c2 = {'--' if th is None else f'{4*th:.6f}'}")
OUT['viscous'] = rowsV
say("   eps = lambda/a0 with lambda = nu/r_*^2 = O(M) and a0 = (M/2)L  =>  eps = O(1/L),")
say("   so the viscous cost is O(1/log Re) in c2 and vanishes asymptotically.")
say(f"   d theta/d eps at eps = 0 = {(theta_visc(1e-6)-0.5)/1e-6:.6f}  (analytic 3/4 = 0.75)")
OUT['dtheta_deps'] = (theta_visc(1e-6)-0.5)/1e-6
say()

# ---- the c2 = 2 statement with the true offsets -------------------------------
say("Exact c2 for the plateau with the measured offsets (kappa = 1/2, a0 = (M/2)L + 0.216773 M,")
say("log Re_E = 2L - 0.7031660), R2 route with q = 0, no viscous loss:")
rows = []
for L in (5, 10, 20, 40, 80, 160, 1e4):
    a0 = 0.5*L + 0.216773
    Td = 0.5/a0
    c2 = Td*(2*L - 0.7031660)
    rows.append(dict(L=L, a0=a0, MTd=Td, logReE=2*L-0.703166, c2=c2))
    say(f"   L = {L:<8.6g}  a0/M = {a0:<10.5f}  M T_d = {Td:<10.6f}  log Re_E = {2*L-0.703166:<10.4f}  c2 = {c2:.6f}")
OUT['c2_with_offsets'] = rows
say("   c2 rises monotonically to 2 from below; c2 <= 2 for every L > 0 here.")
say()
say("HEADLINE:  the Duhamel/Taylor route needs a SIZE bound on beta (C <= 2/3) that a")
say("           generic 'M^3 log^2' estimate does not deliver; the Riccati route needs only")
say("           the SIGN (beta + nu V >= 0) and then gives c2 = 2 exactly.")

json.dump(OUT, open("d1_results.json", "w"), indent=1, default=str)
open("d1_log.txt", "w").write("\n".join(LOG) + "\n")
print("\n[d1 done]")
