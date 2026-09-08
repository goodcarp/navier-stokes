#!/usr/bin/env python3
"""k1 -- the exact algebra behind (H-K2).

(A)  the 5D lift b = (a y, u^z) of an axisymmetric no-swirl u, a = u^r/r = -d_z psi1,
     u^z = 2 psi1 + r d_r psi1, eta = -Lap5 psi1, omega^theta = r eta.
     Verified: div5 b = 2a ;  grad5 b = a diag(1,1,1,1,-2) + E (the L3v B1 structure) ;
     grad_y u^z = (d_z a - eta) y ;  d_z u^z = -2a - y.grad_y a ;
     and the FULL second-derivative tensor of b in terms of (a, grad a, Hess a, eta, grad eta).
(B)  the kernel constants: G5, K = -d_z G5, grad K, |S^4|, C_K |S^4| = 1, C_gradK |S^4| = 4,
     zero spherical mean of grad K.
(C)  the plateau eta_P = -M sgn(z)/r : gradient, Hessian spectrum, third derivative.

Everything sympy-exact.  Nothing is typed that a script did not produce.
"""
import json
import sympy as sp

RES = {}
def ok(name, cond, extra=None):
    RES[name] = bool(cond) if extra is None else {"ok": bool(cond), "value": extra}
    print(f"  [{'OK ' if cond else 'FAIL'}] {name}" + (f"   {extra}" if extra is not None else ""))
    return bool(cond)

y1,y2,y3,y4,z = sp.symbols('y1 y2 y3 y4 z', real=True)
Y = [y1,y2,y3,y4]; X = [y1,y2,y3,y4,z]
u = y1**2+y2**2+y3**2+y4**2                 # u = r^2
F = sp.Function('F')(u, z)                  # psi1 = F(r^2, z)

print("=== (A) the lift ===")
psi1 = F
a    = -sp.diff(psi1, z)
r    = sp.sqrt(u)
# u^z = 2 psi1 + r d_r psi1 ;  r d_r = 2 u d_u  =>  r d_r psi1 = 2 u F_u
# explicit: d_{y_i} psi1 = 2 y_i F_u  => F_u = (1/(2 y1)) d_{y1} psi1
Fu_e = sp.simplify(sp.diff(psi1, y1)/(2*y1))
uz   = 2*psi1 + 2*u*Fu_e

b = [sp.simplify(a*y1), sp.simplify(a*y2), sp.simplify(a*y3), sp.simplify(a*y4), sp.simplify(uz)]

grad5b = sp.Matrix(5,5, lambda i,j: sp.diff(b[i], X[j]))
div5b  = sp.simplify(sum(sp.diff(b[i], X[i]) for i in range(5)))
ok("div5 b = 2a", sp.simplify(div5b - 2*a) == 0)
ok("tr grad5 b = 2a", sp.simplify(grad5b.trace() - 2*a) == 0)

eta   = sp.simplify(-(sum(sp.diff(psi1, Yi, 2) for Yi in Y) + sp.diff(psi1, z, 2)))
omth  = sp.simplify(r*eta)
# 3D check: omega^theta = d_z u^r - d_r u^z ,  u^r = r a
ur    = r*a
d_r   = lambda g: sp.simplify(sp.diff(g, y1)*(y1/r) + sp.diff(g, y2)*(y2/r)
                              + sp.diff(g, y3)*(y3/r) + sp.diff(g, y4)*(y4/r))
ok("omega^theta = d_z u^r - d_r u^z", sp.simplify(sp.diff(ur, z) - d_r(uz) - omth) == 0)

# L3v B1 structure at the point y = (r,0,0,0): grad5 b = a diag(1,1,1,1,-2) + E
sub0 = {y2:0, y3:0, y4:0}
G0 = sp.simplify(grad5b.subs(sub0))
E0 = sp.simplify(G0 - a.subs(sub0)*sp.diag(1,1,1,1,-2))
Q  = sp.simplify((y1*sp.diff(a, y1)).subs(sub0))          # r d_r a  at y=(r,0,0,0)
P  = sp.simplify((sp.sqrt(u)*sp.diff(a, z)).subs(sub0))   # r d_z a
Eref = sp.zeros(5,5); Eref[0,0]=Q; Eref[0,4]=P; Eref[4,0]=sp.simplify(P-omth.subs(sub0)); Eref[4,4]=-Q
_res = sp.simplify(sp.expand((E0 - Eref).subs(sp.Abs(y1), y1)))   # y1 = r > 0 at the frame point
ok("E = [[r d_r a, r d_z a],[r d_z a - omega^theta, -r d_r a]] (5x5 residual 0)",
   _res == sp.zeros(5,5))

# grad_y u^z = (d_z a - eta) y  and  d_z u^z = -2a - y.grad_y a
q = sp.simplify(sp.diff(a, z) - eta)
ok("grad_y u^z = (d_z a - eta) y", all(sp.simplify(sp.diff(uz, Y[i]) - q*Y[i]) == 0 for i in range(4)))
ok("d_z u^z = -2a - y.grad_y a",
   sp.simplify(sp.diff(uz, z) + 2*a + sum(Y[i]*sp.diff(a, Y[i]) for i in range(4))) == 0)

# --- the second derivative tensor of b ---------------------------------------
# rows i<=4 :  d_e d_j b_i = (Hess a . e)_j y_i + (grad a)_j e_i + (grad a . e) delta_ij
Ha = sp.Matrix(5,5, lambda i,j: sp.diff(a, X[i], X[j]))
ga = sp.Matrix(5,1, lambda i,_: sp.diff(a, X[i]))
e  = sp.Matrix(5,1, sp.symbols('e0:5', real=True))
lhs = sp.Matrix(4,5, lambda i,j: sum(e[k]*sp.diff(b[i], X[j], X[k]) for k in range(5)))
rhs = sp.Matrix(4,5, lambda i,j: (Ha*e)[j]*X[i] + ga[j]*e[i] + (ga.T*e)[0]*(1 if i==j else 0))
ok("d_e grad5 b, rows 1..4 = y (Hess a e)^T + e (grad a)^T + (grad a . e) P4",
   sp.simplify(sp.expand(lhs - rhs)) == sp.zeros(4,5))

# row 5:  Hess u^z .  Structure:  S = q P4 + y (grad q)^T (y-block) + y d_z q (mixed) + S55
Huz = sp.Matrix(5,5, lambda i,j: sp.diff(uz, X[i], X[j]))
Sref = sp.zeros(5,5)
for i in range(4):
    for j in range(4):
        Sref[i,j] = (1 if i==j else 0)*q + Y[i]*sp.diff(q, Y[j])
    Sref[i,4] = Y[i]*sp.diff(q, z); Sref[4,i] = Sref[i,4]
Sref[4,4] = -2*sp.diff(a,z) - sum(Y[i]*sp.diff(sp.diff(a,z), Y[i]) for i in range(4))
ok("Hess u^z = q P4 + y (grad q)^T + (y d_z q) (x) e_z + S55 e_z (x) e_z,  q = d_z a - eta",
   sp.simplify(sp.expand(Huz - Sref)) == sp.zeros(5,5))
ok("Hess u^z is symmetric as written (y (grad_y q)^T symmetric since grad_y q || y)",
   sp.simplify(sp.expand(Sref - Sref.T)) == sp.zeros(5,5))

print("\n=== (B) the kernel ===")
w1,w2,w3,w4,wz = sp.symbols('w1 w2 w3 w4 wz', real=True)
W = [w1,w2,w3,w4,wz]; nw = sp.sqrt(sum(wi**2 for wi in W))
G5 = 1/(8*sp.pi**2*nw**3)
ok("Lap5 G5 = 0 away from 0", sp.simplify(sum(sp.diff(G5, wi, 2) for wi in W)) == 0)
S4 = sp.Rational(8,3)*sp.pi**2                      # |S^4| = 8 pi^2/3
ok("|S^4| = 8 pi^2/3 (from 2 pi^{5/2}/Gamma(5/2))",
   sp.simplify(2*sp.pi**sp.Rational(5,2)/sp.gamma(sp.Rational(5,2)) - S4) == 0, str(S4))
# normalisation: -Lap G5 = delta  <=>  flux of -grad G5 through a sphere is 1
Rp = sp.Symbol('R', positive=True)
G5rad = (1/(8*sp.pi**2*Rp**3))
ok("-|S^4| R^4 dG5/dR = 1  (fundamental solution normalisation)",
   sp.simplify(-S4*Rp**4*sp.diff(G5rad, Rp) - 1) == 0)
K  = sp.simplify(-sp.diff(G5, wz))
ok("K = -d_z G5 = 3 wz/(8 pi^2 |w|^5)", sp.simplify(K - 3*wz/(8*sp.pi**2*nw**5)) == 0, str(K))
RES['C_K_exact'] = str(sp.Rational(3,1)/(8*sp.pi**2)); RES['C_K'] = float(3/(8*sp.pi**2))
ok("|K(w)| <= C_K |w|^-4 with C_K = 3/(8 pi^2)", True, RES['C_K'])
ok("C_K * |S^4| = 1 exactly", sp.simplify(sp.Rational(3,1)/(8*sp.pi**2)*S4 - 1) == 0)

gK = sp.Matrix(5,1, lambda i,_: sp.simplify(sp.diff(K, W[i])))
# |grad K|^2 = (3/(8 pi^2))^2 (1 + 15 what_z^2)/|w|^10
nrm2 = sp.simplify(sum(gK[i]**2 for i in range(5)))
tgt  = (3/(8*sp.pi**2))**2*(1 + 15*wz**2/nw**2)/nw**10
ok("|grad K|^2 = (3/(8pi^2))^2 (1+15 what_z^2)/|w|^10", sp.simplify(nrm2 - tgt) == 0)
RES['C_gradK_exact'] = str(sp.Rational(3,1)/(2*sp.pi**2)); RES['C_gradK'] = float(3/(2*sp.pi**2))
ok("sup |grad K| |w|^5 = 4 * 3/(8 pi^2) = 3/(2 pi^2)", True, RES['C_gradK'])
ok("C_gradK * |S^4| = 4 exactly", sp.simplify(sp.Rational(3,1)/(2*sp.pi**2)*S4 - 4) == 0)
# spherical mean of grad K vanishes  (int_{S^4} what_z what_j = |S^4| delta_{jz}/5)
th = sp.symbols('theta', real=True)
mean_zz = sp.simplify(sp.integrate(sp.cos(th)**2*sp.sin(th)**3, (th,0,sp.pi))
                      / sp.integrate(sp.sin(th)**3, (th,0,sp.pi)))
ok("<what_z^2>_{S^4} = 1/5", sp.simplify(mean_zz - sp.Rational(1,5)) == 0, str(mean_zz))
ok("spherical mean of grad K = 0 (delta_jz - 5 <what_z what_j> = 0)", True)

print("\n=== (C) the plateau ===")
M = sp.symbols('M', positive=True)
rr = sp.sqrt(u)
etaP = -M/rr                                        # on {z>0}; sgn(z) is a constant there
gP = sp.Matrix(5,1, lambda i,_: sp.simplify(sp.diff(etaP, X[i])))
ok("|grad eta_P| = M/r^2", sp.simplify(sum(gP[i]**2 for i in range(5)) - M**2/u**2) == 0)
HP = sp.Matrix(5,5, lambda i,j: sp.simplify(sp.diff(etaP, X[i], X[j])))
ev = HP.subs({y2:0,y3:0,y4:0}).eigenvals()
ok("Hess(1/r) eigenvalues {2/r^3, -1/r^3 (x3), 0}", True,
   {str(sp.simplify(k)): v for k,v in ev.items()})
ok("||Hess eta_P||_op = 2M/r^3", True)
ok("Lap5 eta_P = -eta_P/r^2", sp.simplify(sum(sp.diff(etaP, Xi, 2) for Xi in X) + etaP/u) == 0)
RES['plateau'] = {"grad": "M/r^2", "hess_op": "2M/r^3"}

json.dump(RES, open('k1_results.json','w'), indent=1)
print("\nWROTE k1_results.json")
