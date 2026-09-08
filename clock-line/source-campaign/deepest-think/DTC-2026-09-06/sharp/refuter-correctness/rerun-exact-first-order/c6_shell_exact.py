#!/usr/bin/env python3
"""C6 -- the IDEAL (space-filling, extremal) endpoint: the bang-bang shell
   omega^theta = -M sgn(z) on rho0 < |x| < R,  ||omega_0||_inf = M,  a(0,0) = (M/2) log(R/rho0).
Energy by exact 5D zonal-harmonic separation; H_l exact rationals, matching solved in adaptive
precision, radial integrals in closed form.  Every mode's contribution must be >= 0 (it is that mode's
Dirichlet energy) -- checked.

Delta_5[f(rho) C_l^{3/2}(t)] = [f'' + (4/rho) f' - l(l+3) f/rho^2] C_l ;  dsigma_{S^4} = 2 pi^2 (1-t^2) dt.
eta = -M sgn(t)/(rho sqrt(1-t^2)) = (1/rho) sum_{l odd} H_l C_l(t),
   H_l = -2 M int_0^1 C_l sqrt(1-t^2) dt / N_l ,  N_l = 2(l+1)(l+2)/(2l+3).
A_l'' + (4/rho)A_l' - l(l+3)A_l/rho^2 = -H_l/rho ;  particular  P_l = H_l rho/((l+4)(l-1)) (l>1),
   P_1 = -(H_1/5) rho log rho.   Regions: a rho^l  |  P + b rho^l + c rho^{-l-3}  |  d rho^{-l-3}, C^1 matched.
E = (1/pi) int_{R^5} eta psi = 2 pi sum_l N_l H_l int_{rho0}^{R} rho^3 A_l(rho) drho .
"""
import sympy as sp, mpmath as mp, numpy as np, json, hashlib, time
t0 = time.time()
tt = sp.Symbol('t')
def N_l(l): return sp.Rational(2*(l+1)*(l+2), 2*l+3)
def mom(k): return sp.Rational(1,2)*sp.gamma(sp.Rational(k+1,2))*sp.gamma(sp.Rational(3,2))/sp.gamma(sp.Rational(k+1,2)+sp.Rational(3,2))
def h_coef(l):
    C = sp.Poly(sp.expand(sp.gegenbauer(l, sp.Rational(3,2), tt)), tt)
    return sp.nsimplify(sp.simplify(-2*sum(co*mom(k) for (k,), co in C.terms())/N_l(l)), rational=True)

LMAX = 201
LS = list(range(1, LMAX+1, 2))
H = {l: h_coef(l) for l in LS}
print("exact H_l:", [(l, H[l]) for l in LS[:5]])
# independent numeric check of H_l
mp.mp.dps = 30
for l in [1, 5, 21, 61]:
    f = sp.lambdify(tt, sp.gegenbauer(l, sp.Rational(3,2), tt), 'mpmath')
    num = -2*mp.quad(lambda x: mp.sqrt(1-x**2)*f(x), [0, 1])/mp.mpf(sp.Rational(N_l(l)).p)*mp.mpf(sp.Rational(N_l(l)).q)
    print(f"  H_{l}: exact {float(H[l]):+.14f}   quad {float(num):+.14f}")

def mode_integral(l, R, dps):
    """I_l = int_1^R rho^3 A_mid drho, exact matching solved at precision dps."""
    mp.mp.dps = dps
    Rm = mp.mpf(R); Hl = mp.mpf(sp.Rational(H[l]).p)/mp.mpf(sp.Rational(H[l]).q)
    if l == 1:
        P  = lambda r: -(Hl/5)*r*mp.log(r)
        dP = lambda r: -(Hl/5)*(mp.log(r) + 1)
    else:
        P  = lambda r: Hl*r/mp.mpf((l+4)*(l-1))
        dP = lambda r: Hl/mp.mpf((l+4)*(l-1))
    # unknowns a,b,c,d ;  rho0 = 1
    Amat = mp.matrix(4, 4); rhs = mp.matrix(4, 1)
    # a*1 - b*1 - c*1 = P(1)
    Amat[0,0], Amat[0,1], Amat[0,2], Amat[0,3] = 1, -1, -1, 0; rhs[0] = P(1)
    Amat[1,0], Amat[1,1], Amat[1,2], Amat[1,3] = l, -l, l+3, 0; rhs[1] = dP(1)
    # P(R) + b R^l + c R^{-l-3} - d R^{-l-3} = 0
    Amat[2,0], Amat[2,1], Amat[2,2], Amat[2,3] = 0, Rm**l, Rm**(-l-3), -Rm**(-l-3); rhs[2] = -P(Rm)
    Amat[3,0], Amat[3,1], Amat[3,2], Amat[3,3] = 0, l*Rm**(l-1), (-l-3)*Rm**(-l-4), (l+3)*Rm**(-l-4); rhs[3] = -dP(Rm)
    a, b, c, d = mp.lu_solve(Amat, rhs)
    # int_1^R rho^3 * (P + b rho^l + c rho^{-l-3})
    if l == 1:
        Ip = -(Hl/5)*(Rm**5*mp.log(Rm)/5 - (Rm**5 - 1)/25)
    else:
        Ip = Hl/mp.mpf((l+4)*(l-1))*(Rm**5 - 1)/5
    Ib = b*(Rm**(l+4) - 1)/(l+4)
    Ic = c*(mp.log(Rm) if l == 1 else (Rm**(1-l) - 1)/(1-l))
    return Ip + Ib + Ic

OUT = {}
for K in [8, 12, 16, 20]:
    R = mp.mpf(2)**K
    terms = []
    for l in LS:
        dps = 40 + int((2*l + 3)*K*0.30103) + 20
        terms.append(mode_integral(l, R, dps))
    mp.mp.dps = 40
    parts = np.array([float(2*mp.pi*mp.mpf(sp.Rational(N_l(l)).p)/mp.mpf(sp.Rational(N_l(l)).q)
                            * (mp.mpf(sp.Rational(H[l]).p)/mp.mpf(sp.Rational(H[l]).q)) * I) for l, I in zip(LS, terms)])
    neg = int(np.sum(parts < -1e-12*abs(parts).max()))
    E = float(parts.sum()); Rf = float(R); csum = np.cumsum(parts)
    print(f"\nR/rho0 = 2^{K} = {Rf:g}")
    print(f"   modes with NEGATIVE Dirichlet energy (must be 0): {neg}")
    print(f"   partial sums / E at l = 11,41,81,141,{LMAX}: " + ", ".join(f"{csum[i]/E:.7f}" for i in [5,20,40,70,len(LS)-1]))
    print(f"   last term/E = {parts[-1]/E:.3e}   (tail is positive and ~l^-3; truncation bias < {abs(parts[-1])*LMAX/2/E:.1e})")
    print(f"   E = {E:.10e}    E/(M^2 R^5) = {E/Rf**5:.8e}    E/(M^2 R^5 log(R/rho0)) = {E/(Rf**5*np.log(Rf)):.8e}")
    a = 0.5*np.log(Rf); nu = 1.0
    ReE = E**0.4/nu; c1 = (np.log(2)/a)*np.log(ReE)
    print(f"   a_inner = (M/2)log(R/rho0) = {a:.6f}   log Re_E = {np.log(ReE):.6f}  [2 log(R/rho0) = {2*np.log(Rf):.6f}]")
    print(f"   t_d = ln2/a = {np.log(2)/a:.8f}/M     c1 = t_d M log Re_E = {c1:.6f}     [4 ln 2 = {4*np.log(2):.6f}]")
    OUT[str(K)] = dict(E=E, E_over_M2R5=float(E/Rf**5), a=float(a), logRe=float(np.log(ReE)), c1=float(c1), neg=neg)

print("\nE(R) = C M^2 R^5 (1 + O(1/log)) : log Re_E = 2 log(R/rho0) + (2/5) log C + o(1),")
print("so c1 -> (ln2) * 2 log(R/rho0) / ((1/2) log(R/rho0)) = 4 ln 2 = %.10f exactly." % (4*np.log(2)))
print(f"\nelapsed {time.time()-t0:.1f}s  SHA256 {hashlib.sha256(open(__file__,'rb').read()).hexdigest()}")
json.dump(OUT, open('c6_results.json','w'), indent=1, default=float)
