#!/usr/bin/env python3
"""r1 -- refuter's own re-derivation of the kernel constants, Lemma 5.2's outer-region
constants, the shell total-variation density J, and the datum's behaviour at the origin.

Written from scratch (no import from hk2/).  Units rho0 = M = 1 unless stated.
"""
import json, math
import numpy as np
import sympy as sp
from scipy import integrate

RES = {}
def rec(k, v, msg=""):
    RES[k] = v
    print(f"  {k:45s} = {v}" + (f"   {msg}" if msg else ""))

print("=== A. kernel constants, sympy, my own route ===")
n = 5
# fundamental solution of -Delta in R^n:  1/((n-2)|S^{n-1}|) |x|^{2-n}
S4 = 2*sp.pi**sp.Rational(n,2)/sp.gamma(sp.Rational(n,2))
rec("S4_exact", str(sp.simplify(S4)))
rec("S4_float", float(S4))
x = sp.symbols('x1:6', real=True)
rho = sp.sqrt(sum(xi**2 for xi in x))
G = rho**(2-n)/((n-2)*S4)
rec("G5_exact", str(sp.simplify(G)))
assert sp.simplify(G - 1/(8*sp.pi**2*rho**3)) == 0
# check -Delta G = 0 off the origin and flux normalisation: -int_{S^4} dG/drho rho^4 = 1
lapG = sp.simplify(sum(sp.diff(G, xi, 2) for xi in x))
assert lapG == 0
Rr = sp.symbols('R', positive=True)
flux = sp.simplify(-S4*Rr**4*sp.diff(Rr**(2-n)/((n-2)*S4), Rr))
assert flux == 1
rec("flux_normalisation", str(flux))
K = -sp.diff(G, x[4])           # a = -d_z psi1, psi1 = G * eta  =>  a = K * eta
rec("K_exact", str(sp.simplify(K)))
CK = sp.Rational(3, 8)/sp.pi**2
assert sp.simplify(K - CK*x[4]/rho**5) == 0
rec("C_K_exact", str(CK)); rec("C_K", float(CK))
rec("C_K_times_S4", str(sp.simplify(CK*S4)))
gK = [sp.simplify(sp.diff(K, xi)) for xi in x]
gK2 = sp.simplify(sum(g**2 for g in gK)*rho**10/CK**2)
rec("|gradK|^2 rho^10 / C_K^2", str(sp.simplify(gK2)))
assert sp.simplify(gK2 - (1 + 15*x[4]**2/rho**2)) == 0
CgK = 4*CK
rec("C_gradK_exact", str(CgK)); rec("C_gradK", float(CgK))
rec("C_gradK_times_S4", str(sp.simplify(CgK*S4)))
# spherical mean of grad K on |w| = 1: components c*(delta_jz - 5 w_z w_j); <w_z^2> on S^4
th = sp.symbols('theta', real=True)
mz2 = sp.integrate(sp.cos(th)**2*sp.sin(th)**3, (th, 0, sp.pi))/sp.integrate(sp.sin(th)**3, (th, 0, sp.pi))
rec("<w_z^2>_S4", str(mz2)); assert mz2 == sp.Rational(1, 5)

print("\n=== A'. consistency of the kernel with the record's 'M/2 per e-fold' ===")
# a(0) for the bare plateau eta = -M sgn(z)/r on a shell of log-width L:
#   a(0) = int K(-x') eta(x') dx' = C_K * L * int_{S^4} |cos phi|/sin phi dOmega_4
#   int_{S^4} |cos|/sin sin^3 dphi dOmega_3 = |S^3| * 2 * int_0^{pi/2} cos sin^2 = 2 pi^2 * 2/3
ph = sp.symbols('phi', real=True)
ang = 2*sp.pi**2*2*sp.integrate(sp.cos(ph)*sp.sin(ph)**2, (ph, 0, sp.pi/2))
rate = sp.simplify(CK*ang)
rec("a(0) per e-fold, bare plateau (exact)", str(rate), "record (far-near lemma Consequence A): M/2")
assert rate == sp.Rational(1, 2)

print("\n=== B. Lemma 5.2 outer-region constants (my own derivation) ===")
# region rho' > 2 rho: |x-x'| >= rho' - rho >= rho'/2 ; density M rho'^2 J drho'
rp, rr = sp.symbols('rhop rho', positive=True)
for p in (4, 5):
    far = sp.integrate((2/rp)**p*rp**2, (rp, 2*rr, sp.oo))
    near = sp.integrate((2/rr)**p*rp**2, (rp, 0, rr/2))
    rec(f"Lambda_{p} outer (rho'>2rho) / (M J)", str(sp.simplify(far)))
    rec(f"Lambda_{p} inner (rho'<rho/2) / (M J)", str(sp.simplify(near)))

print("\n=== C. the shell total-variation density J, my own quadrature ===")
S3 = 2*math.pi**2
def J_ac(Hh, Hhp):
    # J_ac = |S^3| int_{-1}^1 sqrt(Hh^2 + (1-t^2) Hh'^2) (1-t^2) dt   (Theta == 1 shells)
    f = lambda t: math.sqrt(Hh(t)**2 + (1-t*t)*Hhp(t)**2)*(1-t*t)
    v, err = integrate.quad(f, -1, 1, limit=400, points=[0.0])
    return S3*v, S3*err
# (D-A): Hh(t) = sgn(t) h(phi)/sin(phi),  h = min(1, phi_ax/delta)
def DA_Hh(delta):
    def Hh(t):
        phi = math.acos(t); pax = min(phi, math.pi-phi); s = math.sin(phi)
        return math.copysign(min(1.0, pax/delta)/s, t) if t != 0 else 1.0
    def Hhp(t):     # d/dt of h(phi)/sin(phi), phi = acos t: dphi/dt = -1/s
        phi = math.acos(t); pax = min(phi, math.pi-phi); s = math.sin(phi); c = t
        h = min(1.0, pax/delta)
        hp = 0.0 if pax >= delta else (1.0/delta if phi < math.pi/2 else -1.0/delta)
        dHdphi = (hp*s - h*c)/s**2
        return math.copysign(1.0, t)*dHdphi*(-1.0/s)
    return Hh, Hhp
for dd in (7.5, 15.0):
    delta = math.radians(dd)
    Hh, Hhp = DA_Hh(delta)
    # split at the corners for the quadrature
    f = lambda t: math.sqrt(Hh(t)**2 + (1-t*t)*Hhp(t)**2)*(1-t*t)
    pts = [-math.cos(delta), 0.0, math.cos(delta)]
    v, err = integrate.quad(f, -1, 1, limit=400, points=pts)
    Jac = S3*v; Jj = 2*S3
    rec(f"J_ac_DA_delta{dd}", Jac, f"(quad err {S3*err:.1e}); hk2: 39.162729 / 38.308103")
    rec(f"J_jump_DA", Jj, "hk2: 39.478418 = 4 pi^2")
    rec(f"J_total_DA_delta{dd}", Jac+Jj)
# bare plateau: Hh = 1/sin, exact J_ac = 2 pi^2 (sqrt2 + asinh 1)
rec("J_ac_bare_exact", 2*math.pi**2*(math.sqrt(2)+math.asinh(1.0)))
# (D-B) angular profile  Hh(t) = tanh(s/sd)/s * tanh(t/w)
def DB_Hh(delta_deg, w):
    sd = math.sin(math.radians(delta_deg))
    def Hh(t):
        s = math.sqrt(max(1-t*t, 1e-300)); return math.tanh(s/sd)/s*math.tanh(t/w)
    def Hhp(t, h=1e-6):
        return (Hh(t+h)-Hh(t-h))/(2*h)
    return Hh, Hhp
for dd in (7.5, 15.0):
    Hh, Hhp = DB_Hh(dd, 0.20)
    v, err = J_ac(Hh, Hhp)
    rec(f"J_DB_delta{dd}_w0.20 (Theta==1 shells only)", v, "hk2: 65.625910 / 63.401831")

print("\n=== D. the TRUE shell density of (D-B): J(rho') including the radial ramps ===")
# eta = -M Hh(t) Theta(rho)/rho ; |grad eta| = M sqrt( (Hh (Theta/rho)')^2 + (1-t^2) Hh'^2 Theta^2/rho^4 )
# J(rho') := (1/(M rho'^2)) int_{S^4} |grad eta| rho'^4 dOmega_4 = rho'^2 |S^3| int sqrt(...) (1-t^2) dt
def Theta(rho, R, w0=0.25, w1frac=0.10):
    w1 = w1frac*R
    return 0.5*(math.tanh((rho-1.0)/w0) - math.tanh((rho-R)/w1))
def sech2(u):
    u = abs(u)
    return 0.0 if u > 350 else (2.0*math.exp(-u)/(1.0+math.exp(-2*u)))**2
def Thetap(rho, R, w0=0.25, w1frac=0.10):
    w1 = w1frac*R
    return 0.5*(sech2((rho-1.0)/w0)/w0 - sech2((rho-R)/w1)/w1)
def J_of_rho(rho, R, Hh, Hhp):
    T0 = Theta(rho, R); T1 = Thetap(rho, R)
    Gp = abs(T1/rho - T0/rho**2)          # |d/drho (Theta/rho)|
    f = lambda t: math.sqrt((Hh(t)*Gp)**2 + (1-t*t)*(Hhp(t)*T0/rho**2)**2)*(1-t*t)
    v, err = integrate.quad(f, -1, 1, limit=400, points=[0.0])
    return rho**2*S3*v
Hh, Hhp = DB_Hh(7.5, 0.20)
R = math.exp(10.0)
grid = np.concatenate([np.linspace(0.02, 3.0, 300), np.exp(np.linspace(math.log(3.0), math.log(1.5*R), 300))])
Jr = np.array([J_of_rho(float(x_), R, Hh, Hhp) for x_ in grid])
i = int(np.argmax(Jr))
rec("J_DB_bulk(rho'=e^5)", J_of_rho(math.exp(5.0), R, Hh, Hhp), "should equal the Theta==1 value")
rec("J_DB_sup_over_rho'", float(Jr[i]), f"attained at rho' = {grid[i]:.4f}")
rec("J_DB_at_rho'=1", J_of_rho(1.0, R, Hh, Hhp))
rec("J_DB_at_rho'=0.75", J_of_rho(0.75, R, Hh, Hhp))
rec("J_DB_at_rho'=0.5", J_of_rho(0.5, R, Hh, Hhp))
rec("J_DB_at_rho'=R", J_of_rho(R, R, Hh, Hhp))
rec("J_DB_at_rho'=0.9R", J_of_rho(0.9*R, R, Hh, Hhp))
rec("J_DB_sup/J_DB_bulk", float(Jr[i])/J_of_rho(math.exp(5.0), R, Hh, Hhp))
RES['J_DB_profile'] = [[float(g), float(j)] for g, j in zip(grid[::10], Jr[::10])]

print("\n=== E. the datum at the origin: Theta(0) and the 1/rho singularity of eta ===")
T0 = Theta(0.0, R)
rec("Theta(0) = (1/2)(1 - tanh(rho0/w0)) (L=10)", T0)
rec("Theta(0) exact form", 0.5*(1-math.tanh(4.0)))
# eta near the origin: eta ~ -M Theta(0) Hh(t)/rho ;  sup_t |Hh| = Hh at the layer edge ~ 1/sd
ts = np.linspace(-1+1e-9, 1-1e-9, 200001)
Hv = np.array([Hh(float(t)) for t in ts[::100]])
rec("sup_t |Hh(t)| (D-B, 7.5deg, w=.2)", float(np.max(np.abs(Hv))))
for rr_ in (1e-1, 1e-2, 1e-3, 1e-4):
    rec(f"sup|eta| on |x| = {rr_:g} rho0  (units M/rho0)", float(np.max(np.abs(Hv)))*Theta(rr_, R)/rr_)
# the l = 1 zonal coefficient of Hh (C_1^{3/2}(t) = 3t, weight (1-t^2), N_1 = 12/5)
C1 = lambda t: 3.0*t
N1 = integrate.quad(lambda t: C1(t)**2*(1-t*t), -1, 1)[0]
h1 = integrate.quad(lambda t: Hh(t)*C1(t)*(1-t*t), -1, 1, limit=400, points=[0.0])[0]/N1
rec("hhat_1 (my quadrature)", h1, "hk2 k4: 7.91e-1")
# the origin log:  p_1 = alpha rho log rho solves p'' + 4p'/rho - 4p/rho^2 = -hhat_1 G(rho) with
# G = -M Theta(0)/rho :  5 alpha / rho = hhat_1 M Theta(0)/rho  =>  alpha = hhat_1 M Theta(0)/5,
# psi1 ~ 3 alpha z log rho,  a = -d_z psi1 ~ -3 alpha log rho  -> +infinity as rho -> 0.
alpha = h1*T0/5.0
rec("alpha = hhat_1 Theta(0)/5", alpha)
rec("a(x) ~ -3 alpha log|x| + O(1): coefficient 3 alpha", 3*alpha)
rec("|grad a| ~ 3 alpha/|x|, ||Hess a|| ~ 3 alpha/|x|^2  =>  global K2 = +infinity for (D-B)", True)
# symbolic check of the particular solution
rho_s, al = sp.symbols('rho alpha', positive=True)
p1 = al*rho_s*sp.log(rho_s)
lhs = sp.simplify(sp.diff(p1, rho_s, 2) + 4*sp.diff(p1, rho_s)/rho_s - 4*p1/rho_s**2)
rec("ODE check: (rho log rho)'' + 4(.)'/rho - 4(.)/rho^2", str(lhs))
assert sp.simplify(lhs - 5*al/rho_s) == 0

json.dump(RES, open('r1_results.json', 'w'), indent=1)
print("\nWROTE r1_results.json")
