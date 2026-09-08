#!/usr/bin/env python3
"""r1 -- independent re-derivation and audit of the attempt's algebra.

A. L1 kernel, derived from 3D Biot-Savart from scratch (not quoted from the attempt).
B. sharp-plateau a(0,0) in closed form.
C. the FLOOR IDENTITY: exactly what  F(t) := a(t)(1+a0 t)/a0 >= 1  is equivalent to.
   -> tests the attempt's claim "(H2) is EQUIVALENT to F >= 1".
D. the a0 < 0 sign flip: what the attempt's 'reversed control' actually measures.
E. the Taylor-route constant: does 'M^3 log^2 with an O(1) constant' really mean C = 4?
F. the corollary's constants, recomputed with the attempt's OWN measured offset.
"""
import json, math
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

OUT = {}
def say(s=""):
    print(s, flush=True); LOG.append(s)
LOG = []

# ---------------------------------------------------------------- A. L1 kernel
say("="*78); say("A. L1 kernel a(0,0) = int K omega^theta dV, derived from 3D Biot-Savart"); say("="*78)
# on the axis at height z0, the azimuthal average of Biot-Savart gives
#   u^z(0,0,z0) = (1/4pi) int omega^theta(r,z) * r / (r^2+(z0-z)^2)^{3/2} dV
# so  d_z u^z(0) = (3/4pi) int omega^theta r z / rho^5 dV, and a = -(1/2) d_z u^z.
r, z, z0 = sp.symbols('r z z0', real=True)
kern = r/(r**2 + (z0-z)**2)**sp.Rational(3,2)
dk = sp.simplify(sp.diff(kern, z0).subs(z0, 0))
say(f"   d_z0 [ r/(r^2+(z0-z)^2)^{{3/2}} ] at z0=0  =  {sp.simplify(dk)}")
target = 3*r*z/(r**2+z**2)**sp.Rational(5,2)
say(f"   residual vs 3 r z / rho^5 : {sp.simplify(dk - target)}")
assert sp.simplify(dk-target) == 0
say("   => a(0,0) = -(1/2)(1/4pi) int omega^theta * 3 r z/rho^5 dV = int K omega^theta dV")
say("      with K = -(3/(8 pi)) r z / rho^5.   L1 CONFIRMED (independent derivation).")
OUT["A_L1_kernel_confirmed"] = True

# numeric confirmation on a smooth compact test field with an independent quadrature
say()
say("   numeric confirmation on a smooth compactly-supported test vorticity:")
def a_origin_kernel(om, n=400, rmax=6.0):
    # int K om dV, dV = 2 pi r dr dz, over z in (-rmax,rmax)
    rr = np.linspace(1e-6, rmax, n); zz = np.linspace(-rmax, rmax, 2*n)
    RR, ZZ = np.meshgrid(rr, zz, indexing='ij')
    rho = np.sqrt(RR**2+ZZ**2); rho = np.maximum(rho, 1e-12)
    K = -(3.0/(8*math.pi))*RR*ZZ/rho**5
    F = K*om(RR, ZZ)*2*math.pi*RR
    return np.trapz(np.trapz(F, zz, axis=1), rr)
def a_origin_direct(om, n=400, rmax=6.0, eps=1e-3):
    # a = -(1/2) d_z u^z(0,0,0), u^z on the axis by the exact azimuthal Biot-Savart
    rr = np.linspace(1e-6, rmax, n); zz = np.linspace(-rmax, rmax, 2*n)
    RR, ZZ = np.meshgrid(rr, zz, indexing='ij')
    def uz(zc):
        d = (RR**2 + (zc-ZZ)**2)**1.5
        F = om(RR, ZZ)*RR/d*2*math.pi*RR/(4*math.pi)
        return np.trapz(np.trapz(F, zz, axis=1), rr)
    return -0.5*(uz(eps)-uz(-eps))/(2*eps)
# odd-in-z smooth blob, vanishing near the axis and near the origin (so the log converges)
om_test = lambda R_, Z_: -np.exp(-((np.sqrt(R_**2+Z_**2)-2.0)**2)/0.5)*(R_**2/(R_**2+0.25))*np.tanh(Z_/0.7)
k1 = a_origin_kernel(om_test); k2 = a_origin_direct(om_test)
say(f"   kernel form   a(0,0) = {k1:.8f}")
say(f"   direct BS     a(0,0) = {k2:.8f}     rel diff = {abs(k1-k2)/abs(k1):.3e}")
OUT["A_numeric"] = dict(kernel=k1, direct=k2, rel=abs(k1-k2)/abs(k1))

# ---------------------------------------------------------- B. sharp plateau a(0)
say(); say("="*78); say("B. sharp plateau omega^theta = -M sgn z on rho0<rho<R : a(0,0) in closed form"); say("="*78)
p = sp.Symbol('p')
ang = sp.integrate(sp.sin(p)**2*sp.Abs(sp.cos(p)), (p, 0, sp.pi))
say(f"   int_0^pi sin^2 phi |cos phi| dphi = {ang}")
say(f"   a(0,0) = (3/4) M log(R/rho0) * {ang} = {sp.Rational(3,4)*ang} M log(R/rho0)")
OUT["B_a_origin_sharp"] = str(sp.Rational(3,4)*ang)
say("   => EXACTLY (M/2) log(R/rho0), with NO O(1) offset, for the sharp plateau at the ORIGIN.")
say("      (the frame's +0.216773 M is the offset at the material innermost SHELL, not the origin)")

# ---------------------------------------------------- C. what the floor really tests
say(); say("="*78); say("C. FLOOR IDENTITY  --  audit of 'H2 is EQUIVALENT to a(1+a0 t)/a0 >= 1'"); say("="*78)
say("   D_t a = -a^2 + B,  B := beta + nu V.  Put h = 1/a (a>0). Then h' = 1 - B/a^2, so")
say("       h(t) = 1/a0 + t - int_0^t B/a^2 ds ,")
say("   and  a(t) >= a0/(1+a0 t)  <=>  h(t) <= 1/a0 + t  <=>  int_0^t B/a^2 ds >= 0 .")
say("   Therefore  F(t) >= 1 for all t  <=>  (H2')  int_0^t B/a^2 ds >= 0 for all t,")
say("   which is STRICTLY WEAKER than the pointwise (H2) B >= 0.  The attempt's 'equivalent'")
say("   is a one-way implication.")
# explicit counterexample: B changes sign, floor never dips
a0 = 1.0
def Bfun(t):  # negative on [0.30,0.45], positive elsewhere; net integral positive at all t
    return 1.6 if t < 0.30 else (-0.9 if t < 0.45 else 1.6)
def rhs(t, y): return [-y[0]**2 + Bfun(t)]
sol = solve_ivp(rhs, [0, 1.0], [a0], max_step=1e-4, dense_output=True, rtol=1e-10, atol=1e-12)
ts = np.linspace(0, 1.0, 4001); aa = sol.sol(ts)[0]
F = aa*(1+a0*ts)/a0
say(f"   counterexample: B(t) = +1.6 / -0.9 / +1.6 on [0,.3)/[.3,.45)/[.45,1]")
say(f"     min_t B(t) = -0.9 < 0   (H2 FAILS)   but   min_t F(t) = {F.min():.6f} >= 1  (floor HOLDS)")
OUT["C_counterexample_minF"] = float(F.min())
assert F.min() >= 1.0 - 1e-9

# ------------------------------------------------- D. the a0<0 reversed 'control'
say(); say("="*78); say("D. the reversed control: what F = a(1+a0 t)/a0 does when a0 < 0"); say("="*78)
say("   For a0 < 0, dividing by a0 REVERSES the inequality:")
say("       F(t) >= 1  <=>  a(t) <= a0/(1+a0 t)  <=>  int_0^t B/a^2 ds <= 0 .")
say("   So with B >= 0 (i.e. (H2) HOLDING) and a0 < 0 the statistic necessarily reports F < 1.")
for B0 in (0.0, 2.0, 3.0):
    sol = solve_ivp(lambda t,y: [-y[0]**2 + B0*y[0]**2 if False else -y[0]**2 + B0], [0,0.8], [-1.0],
                    max_step=1e-4, rtol=1e-11, atol=1e-13, dense_output=True)
    ts = np.linspace(0,0.8,2001); aa = sol.sol(ts)[0]
    F = aa*(1-1.0*ts)/(-1.0)
    say(f"     a0 = -1, B == {B0:+.1f} (>=0, so (H2) HOLDS):  min_t F = {F.min():.4f}  -> statistic 'FIRES'")
    OUT.setdefault("D_reversed", []).append(dict(B=B0, minF=float(F.min())))
say("   => the reversed-datum control fires by SIGN CONVENTION, not because (H2) fails.")
say("      Their own d3 log records beta/a^2 in [+2.115,+3.251] on the reversed run, i.e. (H2)")
say("      HOLDS there, while the floor statistic reports 0.1130.  The control is NOT diagnostic")
say("      of (H2); it is a run where (H2) holds and the gate says 'violated'.")

# --------------------------------------------- E. the Taylor route constant C=4?
say(); say("="*78); say("E. Taylor route: does the brief's 'M^3 log^2 with an O(1) constant' mean C = 4?"); say("="*78)
say("   attempt's own translation:  |d_t^2 omega| = |beta||omega| <= C a0^2 * (3/2)M = (3/8) C M^3 L^2")
say("   (using a0 = (M/2)L and |omega| <= (3/2)M).  So constant kappa_2 := (3/8) C, i.e. C = (8/3) kappa_2.")
for k2 in (0.25, 1.0, 1.5):
    say(f"     kappa_2 = {k2:>4}  ->  C = {8*k2/3:.4f}   remainder/gain at x=1/2 = {(3/8)*(8*k2/3):.4f}")
say("   The attempt states 'at the brief's suggested O(1) constant (C = 4)'.  C = 4 corresponds to")
say("   kappa_2 = 3/2, NOT 1.  At kappa_2 = 1 the correct value is C = 8/3 and the remainder/gain")
say("   at x = 1/2 is 1.000, not the attempt's headline 1.5.")
xs = np.linspace(0, 3, 300001)
for C in (2/3, 8/3, 4.0):
    pol = 1 + xs - 0.75*C*xs**2
    say(f"     C = {C:.4f}:  max_x [1 + x - (3/4)C x^2] = {pol.max():.6f}   (needs >= 1.5)")
    OUT.setdefault("E", []).append(dict(C=C, maxpoly=float(pol.max())))
say("   The QUALITATIVE conclusion (the Taylor route fails at an O(1) constant) survives;")
say("   the headline number 1.5x does not track the attempt's own translation formula.")

# ------------------------------------------------------- F. corollary constants
say(); say("="*78); say("F. the corollary constant c2, recomputed with the attempt's OWN measured offset"); say("="*78)
cE = -0.7031660
def c2(L, theta, kappa, ca): return theta*(2*L + cE)/(kappa*L + ca)
say("   as published (kappa = 1/2, c_a = +0.216773, imported from the SHARP-plateau elliptic instrument):")
row = [c2(L, 0.5, 0.5, 0.216773) for L in (5,10,20,40,80,160,1e4)]
say("     L = 5/10/20/40/80/160/1e4 :  " + " / ".join(f"{v:.4f}" for v in row) + "   (rises to 2 from BELOW)")
say("   with the attempt's OWN d2 fit for the ADMISSIBLE datum it actually integrates")
say("   (delta = 7.5 deg, w = 0.20:  kappa = 0.48137, c = -0.19141):")
row2 = [c2(L, 0.5, 0.48137, -0.19141) for L in (5,10,20,40,80,160,1e4)]
say("     L = 5/10/20/40/80/160/1e4 :  " + " / ".join(f"{v:.4f}" for v in row2) + "   (falls to 2.077 from ABOVE)")
say("   with the w = 0.12 datum actually used in every run, using the run's own measured a0:")
for lab, L, a0m, lRe in (("N3", 2.0794, 0.65478, 8.021), ("N4", 2.7726, 0.99303, 9.407),
                          ("N5", 3.4657, 1.33209, 10.794)):
    say(f"     {lab}: L={L:.4f}  a0(tracer, measured) = {a0m:.5f}   formula (M/2)L+0.216773 = {0.5*L+0.216773:.5f}"
        f"   ratio = {a0m/(0.5*L+0.216773):.3f}")
    OUT.setdefault("F_a0", []).append(dict(lab=lab, L=L, a0_measured=a0m, a0_formula=0.5*L+0.216773))
OUT["F_published"] = row; OUT["F_own_constants"] = row2
say("   => 'c2 <= 2 at every finite L, approached from below' is an artefact of substituting an")
say("      offset measured on the SHARP (inadmissible) plateau for the one the attempt itself")
say("      measured on the admissible tapered datum.  With its own constant the limit is 2.077")
say("      and the approach is from ABOVE.")

# ----------------------------------------------------------- G. dimension check
say(); say("="*78); say("G. dimensions of the theorem's headline conclusion"); say("="*78)
say("   Theorem D concludes (1 + a0 tau) >= 3/2, i.e. tau >= 1/(2 a0), so  T_d <= 1/(2 a0).")
say("   The attempt writes 'M0 T_d <= 1/(2 a0)'.  [M0 T_d] = 1 (dimensionless), [1/(2a0)] = time.")
say("   Correct statement:  T_d <= 1/(2 a0), equivalently M0 T_d <= M0/(2 a0).  Typo, not a hole:")
say("   every downstream number (theta = a0 T_d <= 1/2, c2 = 2 theta/kappa) uses the correct form.")

json.dump(OUT, open("r1_results.json","w"), indent=1)
open("r1_log.txt","w").write("\n".join(LOG)+"\n")
print("\n[r1 done]")
