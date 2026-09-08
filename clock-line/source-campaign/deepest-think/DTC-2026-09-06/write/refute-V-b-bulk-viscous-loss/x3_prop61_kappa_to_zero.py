#!/usr/bin/env python3
"""x3 -- REFUTER test of Proposition 6.1 (PROOF.md sec 6).

Prop 6.1 as written claims: over b_kappa(x) = beta kappa sin(x/kappa), ||b'|| = beta is fixed,
||b''|| = beta/kappa, "and (6.1) gives m_tau prop 1/kappa != 0.  Letting kappa -> 0 at fixed
beta makes m_tau arbitrarily large while ||grad u||_inf is fixed, so the supremum of the loss
over drifts of given ||grad u||_inf is +infinity."

The seat's own numerics (b4 S1) sweep kappa UPWARD (1 -> 16), i.e. ||b''|| -> 0, the direction
in which the claim is not at issue.  The claim needs kappa -> 0, exactly where the leading-
order-in-nu formula (6.1) loses its footing (its expansion parameter is sqrt(nu tau)/kappa).

Two independent instruments, both run DOWNWARD in kappa:
   (A) the seat's own Fourier-spectral PDE solver (copied verbatim from b4)
   (B) a Monte-Carlo of the backward SDE  dY = -b(Y) ds + sqrt(2 nu) dW, Y_0 = X(tau),
       loss = E[Y_tau] - x0   (representation-free of the PDE; independent code)
Pre-declared falsification rule X1: Prop 6.1 is REFUTED if loss*kappa is not monotone-ish
constant down to small kappa AND the loss itself turns over and decreases, i.e. if
sup_kappa |loss| is finite.
"""
import json, math
import numpy as np
R = {}

# ---------------- (A) the seat's solver, copied verbatim from b4_hessian_necessity.py ------
def solve(kappa, beta, nu, tau, x0, N=2048, nt=4000):
    Lx = 2*math.pi*kappa
    x = np.arange(N)*Lx/N
    k = 2*math.pi*np.fft.fftfreq(N, d=Lx/N)
    b = beta*kappa*np.sin(x/kappa)
    dt = tau/nt
    E1 = np.exp(-nu*k**2*dt); E2 = np.exp(-nu*k**2*dt/2)
    def NL(xih):
        xi_x = np.real(np.fft.ifft(1j*k*xih))
        return np.fft.fft(-b*xi_x - b)
    xih = np.zeros(N, dtype=complex)
    for _ in range(nt):
        a1 = NL(xih); a2 = NL(E2*(xih + 0.5*dt*a1)); a3 = NL(E2*xih + 0.5*dt*a2)
        a4 = NL(E1*xih + dt*E2*a3)
        xih = E1*xih + dt*(E1*a1 + 2*E2*(a2+a3) + a4)/6.0
    bf = lambda q: beta*kappa*math.sin(q/kappa)
    q = x0; h = tau/4000
    for _ in range(4000):
        k1=bf(q); k2=bf(q+0.5*h*k1); k3=bf(q+0.5*h*k2); k4=bf(q+h*k3)
        q = q + h*(k1+2*k2+2*k3+k4)/6.0
    qq = q % Lx
    xi_at = float(np.real(np.sum(xih*np.exp(1j*k*qq)))/N)
    return q + xi_at - x0, q

# ---------------- (B) independent backward-SDE Monte-Carlo -------------------------------
def mc_loss(kappa, beta, nu, tau, x0, npath=400000, nstep=4000, seed=11):
    """eta(X(tau),tau) = E[eta_0(Y_tau)] = E[Y_tau] with eta_0(x)=x ;  loss = E[Y_tau] - x0.
       Variance reduction by common random numbers against the b == 0 run, whose exact
       answer is E[Xtau + sqrt(2 nu) W_tau] = Xtau: estimator = mean(Y_tau - sqrt(2nu) W_tau)
       - x0, with the same W."""
    rng = np.random.default_rng(seed)
    bf = lambda q: beta*kappa*np.sin(q/kappa)
    # forward characteristic to get X(tau)
    q = x0; h = tau/8000
    for _ in range(8000):
        k1=bf(q); k2=bf(q+0.5*h*k1); k3=bf(q+0.5*h*k2); k4=bf(q+h*k3)
        q = q + h*(k1+2*k2+2*k3+k4)/6.0
    Xtau = float(q)
    dt = tau/nstep; sq = math.sqrt(2*nu*dt)
    Y = np.full(npath, Xtau); Wsum = np.zeros(npath)
    for _ in range(nstep):                    # dY = -b(Y) ds + sqrt(2 nu) dW  (b autonomous)
        dW = rng.standard_normal(npath)
        Y = Y - bf(Y)*dt + sq*dW
        Wsum += sq*dW
    est = Y - Wsum
    m = est.mean(); se = est.std(ddof=1)/math.sqrt(npath)
    return float(m + Xtau - Xtau - x0), float(se), Xtau

beta, tau, nu = 1.6, 0.5, 2e-3
print(f"beta={beta}  tau={tau}  nu={nu}   sqrt(2 nu tau)={math.sqrt(2*nu*tau):.5f}")
rows = []
for kappa in [16.0, 8.0, 4.0, 2.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625]:
    N, nt = 1024, 3000
    L_pde, q = solve(kappa, beta, nu, tau, 0.8*kappa, N=N, nt=nt)
    if kappa in (1.0, 0.125, 0.03125):
        ns = 2000
        L_mc, se, Xt = mc_loss(kappa, beta, nu, tau, 0.8*kappa, npath=60000, nstep=ns)
    else:
        ns, L_mc, se = 0, float('nan'), float('nan')
    rows.append(dict(kappa=kappa, K2=beta/kappa, N=N, nt=nt, nstep=ns,
                     loss_pde=L_pde, loss_pde_times_kappa=L_pde*kappa,
                     loss_mc=L_mc, mc_se=se, z_pde_vs_mc=(L_pde-L_mc)/se))
    print(f" kappa={kappa:9.6f}  ||b''||={beta/kappa:9.2f}   loss(PDE)={L_pde: .6e}"
          f"   loss*kappa={L_pde*kappa: .6e}   loss(MC)={L_mc: .6e} +- {se:.1e}  z={(L_pde-L_mc)/se:+.2f}")
R['kappa_down_sweep'] = rows
best = max(rows, key=lambda r: abs(r['loss_pde']))
R['sup_over_family'] = dict(argmax_kappa=best['kappa'], max_abs_loss=abs(best['loss_pde']),
                            loss_at_smallest_kappa=rows[-1]['loss_pde'],
                            ratio_max_to_smallest=abs(best['loss_pde'])/max(abs(rows[-1]['loss_pde']),1e-300))
print(f"\n sup over the family is attained at kappa = {best['kappa']} with |loss| = "
      f"{abs(best['loss_pde']):.6e};  at the smallest kappa tested the loss has fallen to "
      f"{rows[-1]['loss_pde']:.6e}")
# resolution control on the PDE solver at the two ends
ctrl=[]
for kappa in [1.0, 0.0625, 0.015625]:
    a,_=solve(kappa,beta,nu,tau,0.8*kappa,N=1024,nt=3000)
    b_,_=solve(kappa,beta,nu,tau,0.8*kappa,N=4096,nt=12000)
    ctrl.append(dict(kappa=kappa, loss_coarse=a, loss_fine=b_, rel=(abs(a-b_)/max(abs(b_),1e-300))))
    print(f"   resolution control kappa={kappa}: {a:.8e} vs {b_:.8e} rel {abs(a-b_)/max(abs(b_),1e-300):.2e}")
R['resolution_control']=ctrl
R['verdict_X1'] = ("REFUTED: sup is finite and attained at intermediate kappa"
                   if abs(rows[-1]['loss_pde']) < abs(best['loss_pde'])
                   else "not refuted: loss still growing")
print(" X1 verdict:", R['verdict_X1'])

json.dump(R, open('x3_results.json','w'), indent=1)
print("\nWROTE x3_results.json")
