#!/usr/bin/env python3
"""C7 -- closed form for kappa(s), the per-octave strain increment of the ring stack, to O(s^4).

For a source ring at rho_k >> the evaluation point, the evaluation point may be taken at the origin, where
J5 = (pi/2)/rho'^5 exactly, so the ring's contribution to a is
      kappa = (3/(2 pi)) (pi/2) int int eta r'^3 (0 - z') rho'^{-5} dr' dz'
            = (3/4) (M/r_k) int int G(r',z') f(r',z') dr' dz' ,  f := r'^3 z' / rho'^5,
with G the unit-peak isotropic Gaussian of width sigma = s rho_k centred at (r_k,z_k) = rho_k(sin phi*, cos phi*).
Gaussian moments:  int G = 2 pi sigma^2 ,  int G delta_i delta_j = 2 pi sigma^4 delta_ij , so
      kappa = (3 pi/2)(M sigma^2/r_k) [ f + (sigma^2/2) Laplacian f + O(sigma^4) ] .
Everything below is exact symbolic algebra; the numbers are compared to the c2 quadrature.
"""
import sympy as sp, numpy as np, json, hashlib
r, z, s = sp.symbols('r z s', positive=True)
f = r**3*z/(r**2+z**2)**sp.Rational(5,2)
lap = sp.simplify(sp.diff(f, r, 2) + sp.diff(f, z, 2))
phi = sp.acos(1/sp.sqrt(3))
rk, zk = sp.sin(phi), sp.cos(phi)                      # rho_k = 1
f0 = sp.simplify(f.subs({r: rk, z: zk}))
l0 = sp.simplify(lap.subs({r: rk, z: zk}))
print("f  =", f)
print("Laplacian f =", sp.simplify(lap))
print("at phi* (rho=1):  f0 =", sp.nsimplify(f0), "=", float(f0), "  (= sin^3 cos = sin^2 cos * sin)")
print("                  lap f0 =", sp.nsimplify(sp.simplify(l0)), "=", float(l0))
print("ratio lap f0 / f0 =", sp.nsimplify(sp.simplify(l0/f0)), "=", float(l0/f0))
kap = (3*sp.pi/2)*(s**2/rk)*(f0 + (s**2/2)*l0)
kap = sp.simplify(sp.expand(kap))
lead = (3*sp.pi/2)*s**2*sp.sin(phi)**2*sp.cos(phi)
print("\nkappa(s)/M =", sp.nsimplify(sp.simplify(kap)))
print("leading    =", sp.nsimplify(sp.simplify(lead)), "= (pi/sqrt3) s^2 :", sp.simplify(lead - sp.pi/sp.sqrt(3)*s**2) == 0)
rat = sp.simplify(kap/lead)
print("kappa/leading =", sp.nsimplify(sp.expand(rat)))

print("\ncompare with the c2 quadrature (A_25, unnormalised):")
C2 = json.load(open('c2_results.json'))
rows = []
for sv in [0.04, 0.06, 0.08, 0.10, 0.125]:
    num = C2['kappa'][str(sv)]
    lead_v = float(np.pi/np.sqrt(3)*sv*sv)
    pred = float(kap.subs(s, sv))
    rows.append(dict(s=sv, quad=num, lead=lead_v, pred_O_s4=pred,
                     rat_quad=num/lead_v, rat_pred=pred/lead_v, rel=abs(num-pred)/num))
    print(f"  s={sv:6.4f}  quad {num:.8e}   O(s^4) formula {pred:.8e}   rel {abs(num-pred)/num:.2e}"
          f"    quad/lead {num/lead_v:.6f}   formula/lead {pred/lead_v:.6f}")
print("\nSHA256", hashlib.sha256(open(__file__,'rb').read()).hexdigest())
json.dump(dict(rows=rows, ratio=str(sp.nsimplify(sp.expand(rat)))), open('c7_results.json','w'), indent=1, default=float)
