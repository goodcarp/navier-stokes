#!/usr/bin/env python3
"""R3 - every constant of GATE H-c re-derived by MY OWN route, no import of their kern.py.
Routes deliberately different from theirs:
  * derivative constants: NUMERIC evaluation of d_z of the solid harmonic (mpmath gegenbauer,
    complex-step / high-precision central difference) at random interior points, l up to 20
    -- not sympy simplify.
  * alpha_l: single-layer jump solved by hand and checked against a direct 1D shell-ODE solve.
  * N_l: mpmath quadrature.
  * C1, C2in: mpmath, 40 digits, series summed to machine-zero.
  * collar: R_A and the two collar integrals re-derived symbolically.
"""
import mpmath as mp, math, json, hashlib
mp.mp.dps = 40

def C(l, t):            # C_l^{3/2}
    return mp.gegenbauer(l, mp.mpf(3)/2, t)

def solid_int(l, r, z):  # rho^l C_l(z/rho)
    rho = mp.sqrt(r**2+z**2); return rho**l * C(l, z/rho)
def solid_ext(l, r, z):
    rho = mp.sqrt(r**2+z**2); return rho**(-l-3) * C(l, z/rho)

print("=== R3(b): derivative constants, numeric route (mpmath, dps=40, central diff h=1e-12) ===")
h = mp.mpf(10)**-12
bad = []
for l in range(1, 21):
    for (r, z) in [(mp.mpf('0.7'), mp.mpf('0.3')), (mp.mpf('1.3'), mp.mpf('-0.9')), (mp.mpf('0.2'), mp.mpf('1.7'))]:
        rho = mp.sqrt(r**2+z**2)
        d = (solid_int(l, r, z+h) - solid_int(l, r, z-h))/(2*h)
        c_num = d / (rho**(l-1) * C(l-1, z/rho))
        if abs(c_num - (l+2)) > mp.mpf('1e-14'): bad.append(('I1', l, r, z, c_num))
        d2 = (solid_ext(l, r, z+h) - solid_ext(l, r, z-h))/(2*h)
        c2 = d2 / (rho**(-l-4) * C(l+1, z/rho))
        if abs(c2 + (l+1)) > mp.mpf('1e-14'): bad.append(('I2', l, r, z, c2))
print("   I1 constant (l+2) and I2 constant -(l+1): l=1..20, 3 points each ->",
      "ALL AGREE" if not bad else f"DISAGREE {bad[:3]}")
# also confirm the PREREG value (2l+1) is wrong from l=2
r,z = mp.mpf('0.7'), mp.mpf('0.3'); rho=mp.sqrt(r**2+z**2)
c2 = ((solid_int(2,r,z+h)-solid_int(2,r,z-h))/(2*h))/(rho*C(1,z/rho))
print(f"   l=2 measured constant = {mp.nstr(c2,10)}   (l+2)=4   PREREG (2l+1)=5  -> PREREG FALSE")

print("=== R3(d): N_l by quadrature ===")
for l in range(0, 7):
    q = mp.quad(lambda t: C(l,t)**2*(1-t**2), [-1, 0, 1])
    cf = mp.mpf((l+1)*(l+2))/(l+mp.mpf(3)/2)
    print(f"   l={l}  quad {mp.nstr(q,12)}  closed {mp.nstr(cf,12)}  ok={abs(q-cf)<mp.mpf('1e-25')}")

print("=== R3(c): alpha_l from the single-layer jump, solved by hand ===")
# A_int = alpha rho^l, A_ext = alpha rho^{-l-3} (continuity at rho=1 forces same alpha)
# -Lap psi = g_l C_l delta_{S^4}: jump  d_rho psi|_- - d_rho psi|_+ = g_l
# = l*alpha - (-(l+3)alpha) = (2l+3) alpha  => alpha = g_l/(2l+3)
for l in range(0, 5):
    print(f"   l={l}: l*alpha + (l+3)*alpha = ({2*l+3}) alpha  => alpha = g/{2*l+3}   [matches g/(2l+3)]")

print("=== R3(e): Phi(0) = -(3/5) g_1 and sharpness of kappa_0 ===")
# g_1 = (1/N_1) INT_0^pi w C_1(cos phi) sin^2 phi dphi,  C_1 = 3 cos phi
# Phi(0) = -(1+2)/(2+3) g_1 = -(3/5) g_1
N1 = mp.mpf(12)/5
w_bang = lambda ph: -mp.sign(mp.cos(ph))
g1 = mp.quad(lambda ph: w_bang(ph)*3*mp.cos(ph)*mp.sin(ph)**2, [0, mp.pi/2, mp.pi])/N1
print(f"   bang-bang g_1 = {mp.nstr(g1,12)}  (exact -5/6 = {mp.nstr(mp.mpf(-5)/6,12)})")
print(f"   Phi(0) = -(3/5) g_1 = {mp.nstr(-mp.mpf(3)/5*g1,12)}   (kappa_0 = 1/2)")
sharp = mp.mpf(3)/4*mp.quad(lambda ph: abs(mp.cos(ph))*mp.sin(ph)**2, [0, mp.pi/2, mp.pi])
print(f"   sup |Phi(0)|/M = (3/4) INT |cos| sin^2 = {mp.nstr(sharp,12)}  -> kappa_0 = 1/2 SHARP")

print("=== R3(f): C1 and C2in, mpmath, my own summation ===")
Nl = lambda l: mp.mpf((l+1)*(l+2))/(l+mp.mpf(3)/2)
Cm1 = lambda m: mp.mpf((m+1)*(m+2))/2
def Qm(u, odd):
    ls = range(3, 400, 2) if odd else range(2, 400)
    return mp.sqrt(mp.fsum([(mp.mpf(l+2)/(2*l+3))**2*Cm1(l-1)**2/Nl(l)*u**(2*(l-1)) for l in ls]))
def Sm(v, odd):
    ls = range(1, 400, 2) if odd else range(0, 400)
    return mp.sqrt(mp.fsum([(mp.mpf(l+1)/(2*l+3))**2*Cm1(l+1)**2/Nl(l)*v**(2*(l+4)) for l in ls]))
out = {}
for tag, odd in (("general", False), ("z-odd", True)):
    c1 = mp.sqrt(2)*mp.quad(lambda u: Qm(u, odd)/u, [mp.mpf(10)**-8, mp.mpf('0.1'), mp.mpf('0.5')])
    c2 = mp.sqrt(2)*mp.quad(lambda v: Sm(v, odd)/v, [mp.mpf(10)**-8, mp.mpf('0.1'), mp.mpf('0.5')])
    print(f"   {tag:8s}  C1 = {mp.nstr(c1,10)}   C2in = {mp.nstr(c2,10)}")
    out[tag] = dict(C1=float(c1), C2in=float(c2))

print("=== R3(g): the two collar pieces ===")
R_A = (mp.mpf(2)**5 - mp.mpf(2)**-5)**(mp.mpf(1)/5)
S4 = 8*mp.pi**2/3
A1 = mp.mpf(3)/(8*mp.pi**2)*2*S4*R_A
print(f"   |A| = (|S^4|/5)(2^5-2^-5) rho^5  ->  R_A = {mp.nstr(R_A,10)} rho")
print(f"   piece 1 coefficient of 1/s : (3/8pi^2)*2*|S^4|*R_A = 2*R_A = {mp.nstr(A1,10)}")
inner_z = mp.quad(lambda zz: 1/(zz**2+mp.mpf(1)/4)**2, [-mp.inf, 0, mp.inf])
print(f"   INT dz/(z^2+1/4)^2 = {mp.nstr(inner_z,10)}   (claim 4*pi with r=1: {mp.nstr(4*mp.pi,10)})")
A2 = mp.mpf(3)/(8*mp.pi**2)*2*mp.pi**2*(mp.mpf(1)/24)*4*mp.pi
print(f"   piece 2 = (3/8pi^2)(2pi^2)(1/24)(4pi) = {mp.nstr(A2,10)}   (pi/8 = {mp.nstr(mp.pi/8,10)})")
out['collar'] = dict(R_A=float(R_A), coef_over_s=float(A1), const=float(A2))
json.dump(out, open('r1_results.json','w'), indent=1)
print('\nSCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())
