#!/usr/bin/env python3
"""C8 -- EXACT viscous first-order rate for the two data, symbolically.

For axisymmetric no-swirl,  D_t eta = nu Delta_5 eta,  Delta_5 = d_rr + (3/r) d_r + d_zz,  eta = omega^theta/r.
The relative first-order rate of the PEAK of |omega^theta| at a material point is
      (d/dt log |omega^theta|)|_{t=0} = a + nu (Delta_5 eta)/eta          (since omega = r eta, d log r/dt = a).

(A) BANG-BANG SHELL  eta = -M sgn(z)/r  on rho0 < |x| < R.  Away from the z = 0 sheet and the two edges
    eta is a smooth function of r alone, and Delta_5 (1/r) = 2/r^3 - 3/r^3 = -1/r^3, so
          nu (Delta_5 eta)/eta = -nu/r^2      EXACTLY, at every interior point.
    (Equivalently in omega form: nu(Delta - 1/r^2)omega^theta = -nu omega^theta/r^2 for locally constant omega^theta.)
    With the viscous floor rho0 = sqrt(nu/M): at the innermost radius nu/r^2 = M/sin^2(phi), an O(M)
    penalty independent of R -- while a = (M/2) log(R/rho0) grows.  Net growth once log(R/rho0) > 2/sin^2 phi.

(B) GAUSSIAN-CORE RING of width sigma at (r_k,z_k), sigma << r_k.  At the core centre
          Delta_5 eta / eta = -2/sigma^2 + (3/r_k)(d_r eta)/eta + ... = -2/sigma^2 + O(1/(sigma r_k)),
    so the penalty is 2 nu/sigma^2 = 2M (sigma_0/sigma)^2 with the floor sigma_0 = sqrt(nu/M):
    an O(M) penalty at the INNERMOST core, which the stack must beat with kappa(s) N per octave.
Both statements are verified symbolically below.
"""
import sympy as sp, numpy as np, json, hashlib
r, z, sig, rk, zk, nu, Mv, phi = sp.symbols('r z sigma r_k z_k nu M phi', positive=True)
L5 = lambda f: sp.diff(f, r, 2) + 3*sp.diff(f, r)/r + sp.diff(f, z, 2)

print("(A) shell:  eta = -M/r (upper half)")
eta_A = -Mv/r
print("    Delta_5 eta =", sp.simplify(L5(eta_A)), "    nu Delta_5 eta / eta =", sp.simplify(nu*L5(eta_A)/eta_A))
assert sp.simplify(nu*L5(eta_A)/eta_A + nu/r**2) == 0
print("    => EXACT relative viscous rate = -nu/r^2.")
# omega form cross-check
om = -Mv                          # locally constant omega^theta
lap_om = sp.simplify(sp.diff(om, r, 2) + sp.diff(om, r)/r + sp.diff(om, z, 2) - om/r**2)
print("    omega-form cross-check: (Delta - 1/r^2) omega^theta / omega^theta =", sp.simplify(lap_om/om), " (must be -1/r^2)")
assert sp.simplify(lap_om/om + 1/r**2) == 0

print("\n(B) Gaussian eta-ring: eta = -(M/r_k) exp(-((r-r_k)^2+(z-z_k)^2)/(2 sigma^2))")
eta_B = -(Mv/rk)*sp.exp(-((r-rk)**2 + (z-zk)**2)/(2*sig**2))
rate_B = sp.simplify((nu*L5(eta_B)/eta_B).subs({r: rk, z: zk}))
print("    nu Delta_5 eta/eta at the core centre =", rate_B, "  (must be -2 nu/sigma^2)")
assert sp.simplify(rate_B + 2*nu/sig**2) == 0

print("\nNUMBERS (viscous floor: the SMALLEST scale of the datum sits at sqrt(nu/M))")
ps = float(np.arccos(1/np.sqrt(3)))
print(f"  (A) shell, rho0 = sqrt(nu/M): penalty at the inner edge = nu/(rho0 sin phi)^2 = M/sin^2 phi;")
print(f"      at phi* : M/sin^2 phi* = {1/np.sin(ps)**2:.4f} M   (sin^2 phi* = 2/3)")
print(f"      stretching a = (M/2) log(R/rho0).  NET GROWTH iff log(R/rho0) > {2/np.sin(ps)**2:.4f}")
C6 = json.load(open('c6_results.json'))
C = C6['8']['E_over_M2R5']
Lstar = 2/np.sin(ps)**2
print(f"      i.e. log Re_E = 2 log(R/rho0) + (2/5)log(E/(M^2R^5)) > {2*Lstar + 0.4*np.log(C):.4f}"
      f"   => Re_E > {np.exp(2*Lstar + 0.4*np.log(C)):.1f}   -- an O(1) threshold.")
C2 = json.load(open('c2_results.json'))
print("  (B) ring stack, sigma_0 = sqrt(nu/M): penalty at the inner core = 2M; needs kappa(s) N > 2:")
for s in ['0.04','0.06','0.08','0.1','0.125']:
    k = C2['a_inner'][s]['slope']
    print(f"      s={float(s):6.4f}  kappa={k:.6e}  N > {2/k:7.1f}  => log Re_E > {2*np.log(2)*2/k:8.1f}"
          f"  (Re_E > 10^{2*np.log(2)*2/k/np.log(10):.0f})")
print("\n  => the shell is viscously NEUTRAL at O(1) Reynolds number; the discrete ring stack is not.")
print("     A ring's core has curvature 1/sigma^2 at ITS OWN peak; a plateau has none.")
print("SHA256", hashlib.sha256(open(__file__,'rb').read()).hexdigest())
