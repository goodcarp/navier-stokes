"""
Numbers for the BFG source-adjudication note.  Every figure in NOTE.md that is
tagged [N#] comes from here.  Run:  python3 scripts/adjudication_numbers.py
"""
import numpy as np, sympy as sp
from scipy import integrate

out = []
def P(tag, s): out.append(f"[{tag}] {s}"); print(f"[{tag}] {s}")

# ---------------------------------------------------------------- N1
# nu-uniformity of BFG Thm 8/10.  BFG work at nu = 1 (their (5)).
# Map:  u(x,t) = nu * v(x, nu t)  sends a nu-viscosity NS solution to a nu=1 one.
P("N1a", "--- nu-uniformity of BFG Thm 8/10 ---")
x1,x2,x3,t,nu,lam = sp.symbols('x1 x2 x3 t nu lam', positive=True)
X = sp.Matrix([x1,x2,x3])
v = sp.Matrix([sp.Function('v1')(x1,x2,x3,t), sp.Function('v2')(x1,x2,x3,t), sp.Function('v3')(x1,x2,x3,t)])
q = sp.Function('q')(x1,x2,x3,t)
def curl(f):
    return sp.Matrix([sp.diff(f[2],x2)-sp.diff(f[1],x3),
                      sp.diff(f[0],x3)-sp.diff(f[2],x1),
                      sp.diff(f[1],x1)-sp.diff(f[0],x2)])
def NS(u,p,visc,tt):
    lap = sp.Matrix([sum(sp.diff(u[i],c,2) for c in (x1,x2,x3)) for i in range(3)])
    adv = sp.Matrix([sum(u[j]*sp.diff(u[i],(x1,x2,x3)[j]) for j in range(3)) for i in range(3)])
    return sp.simplify(sp.Matrix([sp.diff(u[i],tt) for i in range(3)]) + adv - visc*lap
                       + sp.Matrix([sp.diff(p,c) for c in (x1,x2,x3)]))
# substitute u(x,t)=nu v(x,nu t), p(x,t)=nu^2 q(x,nu t)
s = sp.Symbol('s', positive=True)
vs = v.subs(t, nu*t); qs = q.subs(t, nu*t)
u  = nu*vs
p  = nu**2*qs
res = NS(u,p,nu,t)
res1 = NS(v,q,1,t).subs(t,nu*t)
P("N1b", f"NS_nu[nu*v(x,nu t)] - nu^2 * NS_1[v](x,nu t) = {sp.simplify(res - nu**2*res1).T}  (zero vector => exact)")
# vorticity: omega_u(x,t) = nu omega_v(x,nu t)  =>  ||omega_u(0)||_inf = nu ||omega_v(0)||_inf
# BFG at nu=1:  tau <= 1/(c ||omega_v(0)||_inf).  t = tau/nu.
M_v, c = sp.symbols('M_v c', positive=True)
M_u = nu*M_v
T_t = (1/(c*M_v))/nu
P("N1c", f"T*||omega_0||_inf at viscosity nu = {sp.simplify(T_t*M_u)}  (nu cancels exactly => the clock is nu-uniform)")

# ---------------------------------------------------------------- N2
# Structural fact separating the velocity route (Kukavica/GIM, divergence form,
# kernel = grad G, mean ZERO) from BFG's vorticity route (kernel = G, mean ONE).
P("N2a", "--- heat-kernel moments in R^3, G(z,tau)=(4 pi tau)^{-3/2} exp(-|z|^2/4tau) ---")
def G(r, tau): return (4*np.pi*tau)**-1.5*np.exp(-r*r/(4*tau))
for tau in (1.0, 1e-2, 1e-4):
    I0,_  = integrate.quad(lambda r: 4*np.pi*r*r*G(r,tau), 0, 200*np.sqrt(tau), limit=400)
    Igrad,_= integrate.quad(lambda r: 4*np.pi*r*r*(r/(2*tau))*G(r,tau), 0, 200*np.sqrt(tau), limit=400)
    P("N2b", f"tau={tau:<8g}  int G dz = {I0:.12f}   int |grad G| dz = {Igrad:.8f}   "
             f"vs 2/sqrt(pi tau) = {2/np.sqrt(np.pi*tau):.8f}")
# the vector integral of grad G vanishes by odd symmetry:
zs = np.random.default_rng(0).normal(scale=np.sqrt(2*1.0), size=(2_000_000,3))  # ~ G(.,1)
P("N2c", f"MC mean of grad G / G = -z/(2tau) under G(.,1):  {(-zs/2).mean(axis=0)}  "
         f"(-> 0 : int grad G dz = 0 exactly, by oddness)")
P("N2d", "consequence: adding a constant C to the integrand changes int G*(.) by C but "
         "int grad G*(.) by 0.  The velocity/divergence-form mild equation may therefore "
         "subtract the local average for free; BFG's (8) third term (kernel = G, integrand "
         "|omega||grad u| >= 0) may not.")

# ---------------------------------------------------------------- N3
# What clock DOES follow from the prior art BFG cite (GIM, velocity, T >= C ||u_0||_inf^{-2})
P("N3a", "--- vorticity clock implied by GIM Remark (i) [velocity] via Biot-Savart ---")
M, W2, rho = sp.symbols('M W2 rho', positive=True)   # M=||om||_inf, W2=||om||_2
# |u| <~ int_{|z|<rho}|om|/|z|^2 + int_{|z|>rho}|om|/|z|^2 <= C rho M + C rho^{-1/2} W2
expr = rho*M + rho**sp.Rational(-1,2)*W2
rho_star = sp.solve(sp.diff(expr,rho), rho)[0]
U = sp.simplify(expr.subs(rho, rho_star))
P("N3b", f"optimal rho = {sp.simplify(rho_star)} ;  ||u_0||_inf <~ {U}")
T_gim = 1/U**2
P("N3c", f"GIM T >= C/||u_0||_inf^2  =>  M*T >~ {sp.simplify(sp.powsimp(M*T_gim, force=True))}")
Re_om = W2**4/M                        # dimensionless at nu=1 (checked below)
P("N3d", f"M*T >~ Re_om^(-1/3) with Re_om = ||om_0||_2^4/||om_0||_inf ; "
         f"check M*T * Re_om**Rational(1,3) = {sp.simplify(sp.powsimp(M*T_gim*Re_om**sp.Rational(1,3), force=True))}")
# scaling check: om -> lam^2 om(lam x): M->lam^2 M, ||om||_2 -> lam^{1/2}||om||_2
P("N3e", f"scaling of Re_om under om->lam^2 om(lam x): "
         f"{sp.simplify((lam**sp.Rational(1,2)*W2)**4/(lam**2*M) / Re_om)}  (=1 => dimensionless)")

# ---------------------------------------------------------------- N4
# Independent counterexample to  int |f|/(|x|+1)^4 dx <= c ||f||_BMO  with COMPACT SUPPORT.
# h_R(x) = (1 - log(1+|x|)/log R)_+ .  ||h_R||_BMO <= ||log(1+|x|)||_BMO / log R <= C/log R.
P("N4a", "--- own compactly supported family h_R = (1 - log(1+|x|)/log R)_+ ---")
def I_of_R(R):
    # substitute u = log(1+r): integrand becomes 4 pi (e^u-1)^2 (1 - u/log R) e^{-3u}
    LR = np.log(R)
    f = lambda u: 4*np.pi*(np.expm1(u))**2*(1-u/LR)*np.exp(-3*u)
    val,_ = integrate.quad(f, 0, LR, limit=800)
    return val
lim = 4*np.pi/3
for R in (1e2, 1e4, 1e8, 1e16, 1e32, 1e64):
    I = I_of_R(R); LR = np.log(R)
    P("N4b", f"R=1e{int(np.log10(R)):<3d} int h_R/(1+|x|)^4 = {I:.6f}   1/log R = {1/LR:.3e}   "
             f"ratio I*log R (grows linearly) = {I*LR:9.2f}   I/(4pi/3) = {I/lim:.4f}")
P("N4c", f"4*pi/3 = {lim:.6f} ; the left side tends to 4pi/3 while ||h_R||_BMO -> 0 like 1/log R, "
         f"so no absolute constant c can exist.  (Third independent family; agrees with the two refuter seats.)")

# ---------------------------------------------------------------- N5
# Repaired Picard iteration: what the corrected inequality costs.
P("N5a", "--- BFG's iteration with the local-average term restored ---")
A, c0, c1, Lam, tt = sp.symbols('A c0 c1 Lambda t', positive=True)
# sup||om^{n+1}|| <= c0 M + c1 t Lam (sup||om^n||)^2 ; closed ball of radius 2 c0 M is invariant iff
cond = sp.solve(sp.Eq(c0*M + c1*tt*Lam*(2*c0*M)**2, 2*c0*M), tt)[0]
P("N5b", f"invariance of the ball {{||om||<=2 c0 M}} holds iff t <= {sp.simplify(cond)}")
P("N5c", "Lambda = 1 (BFG as displayed) -> T ~ 1/(4 c0 c1 M): their Theorem 8. "
         "Lambda = 1 + log_+(ell/sqrt(nu/M)) (the restored local average) -> T ~ 1/(M (1+log Re)): "
         "a logarithmic clock.  The two statements differ ONLY by this factor.")

# ---------------------------------------------------------------- N6
# Sizing the Euler-side tension (Kim-Jeong Prop 4.1 as quoted by the literature seat).
P("N6a", "--- Euler-side tension: Kim-Jeong Prop 4.1 exponents ---")
m = sp.Symbol('m', positive=True); alpha, c1k, c2, eps = sp.symbols('alpha c_1 c_2 epsilon', positive=True)
Tm = c1k*(1-alpha)*m**((alpha-1)/2)
P("N6b", f"T(m) = {Tm}; growth >= (1/4) m^(c_2-alpha); ||om_0||_(L1 cap Linf) <= eps")
for a in (0.10, 0.20, 0.24):
    f = sp.lambdify(m, Tm.subs({alpha:a, c1k:1}), 'numpy')
    vals = [float(f(mm)) for mm in (1e2,1e4,1e8,1e16)]
    P("N6c", f"alpha={a}:  M_0*T(m) <= eps*{['%.3e'%v for v in vals]} at m=1e2,1e4,1e8,1e16 -> 0")
P("N6d", "so on the Euler side M_0 * (doubling time) -> 0 along that family, i.e. the Euler "
         "analogue of BFG Thm 10 is FALSE.  BFG Thm 10 is nu-uniform [N1c].  A rigorous inviscid "
         "limit on the Euler existence interval would therefore contradict it.  NOT ESTABLISHED "
         "here: the limit is not justified in this note.")

open('scripts/adjudication_numbers.out','w').write("\n".join(out)+"\n")
