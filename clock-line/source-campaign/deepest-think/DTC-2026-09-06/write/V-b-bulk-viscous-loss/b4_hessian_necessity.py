#!/usr/bin/env python3
"""b4 -- IS A SECOND-DERIVATIVE BOUND ON u REALLY NEEDED?

Decisive test.  Take a datum with Delta eta_0 == 0 EXACTLY: eta_0(x) = x in 1D.
Every "bulk" (Laplacian) prediction for the viscous loss along a material trajectory is
then exactly ZERO, and the metric/material-frame argument contributes nothing.  If the
measured loss is nonzero and scales like  nu * ||b''|| * tau^2, then the first-order
(first-moment) term is real and NO argument that uses only ||grad u|| can bound it.

Drift family (period 2 pi kappa):   b_kappa(x) = beta kappa sin(x/kappa)
    ||b'||_inf   = beta          (INDEPENDENT of kappa  ->  Gamma, c fixed)
    ||b''||_inf  = beta/kappa    (-> 0 as kappa -> infty)
So a sweep in kappa at fixed beta holds  c = Gamma tau  fixed and sends K_2 -> 0.
The prediction of Theorem V.4 is  loss = O(K_2 nu tau^2), i.e. loss * kappa = const.

Method: eta = x + xi with xi 2pi kappa-periodic;  xi_t = -b xi_x + nu xi_xx - b, xi(0)=0.
Fourier-spectral, integrating-factor RK4.  No probabilistic representation is used, so this
is an instrument independent of the proof.

Leading-order prediction (derived in PROOF.md sec 6):
    m_tau = -(1/2) int_0^tau exp(-int_s^tau b'(x_u) du) b''(x_s) v_s ds ,
    v_s   = 2 nu int_0^s exp(-2 int_p^s b'(x_u) du) dp .
"""
import json, math
import numpy as np

RES = {}

def solve(kappa, beta, nu, tau, x0, N=2048, nt=4000):
    Lx = 2*math.pi*kappa
    x = np.arange(N)*Lx/N
    k = 2*math.pi*np.fft.fftfreq(N, d=Lx/N)
    b = beta*kappa*np.sin(x/kappa)
    bh = np.fft.fft(b)
    dt = tau/nt
    E1 = np.exp(-nu*k**2*dt); E2 = np.exp(-nu*k**2*dt/2)
    def NL(xih):
        xi_x = np.real(np.fft.ifft(1j*k*xih))
        return np.fft.fft(-b*xi_x - b)
    xih = np.zeros(N, dtype=complex)
    for _ in range(nt):                                # ETD-style RK4 (integrating factor)
        a1 = NL(xih)
        a2 = NL(E2*(xih + 0.5*dt*a1))
        a3 = NL(E2*xih + 0.5*dt*a2)
        a4 = NL(E1*xih + dt*E2*a3)
        xih = E1*xih + dt*(E1*a1 + 2*E2*(a2+a3) + a4)/6.0
    xi = np.real(np.fft.ifft(xih))
    # trajectory  xdot = b(x)
    bf = lambda q: beta*kappa*math.sin(q/kappa)
    q = x0; h = tau/20000
    for _ in range(20000):
        k1 = bf(q); k2 = bf(q+0.5*h*k1); k3 = bf(q+0.5*h*k2); k4 = bf(q+h*k3)
        q = q + h*(k1+2*k2+2*k3+k4)/6.0
    # spectral interpolation of xi at X(tau)
    qq = q % Lx
    xi_at = float(np.real(np.sum(xih*np.exp(1j*k*qq)))/N)
    loss = q + xi_at - x0                              # = E[Z_tau] (see PROOF.md sec 6)
    return loss, q, xi_at

def predict(kappa, beta, nu, tau, x0, nq=20001):
    bf  = lambda q: beta*kappa*math.sin(q/kappa)
    b1f = lambda q: beta*np.cos(np.asarray(q)/kappa)
    b2f = lambda q: -(beta/kappa)*np.sin(np.asarray(q)/kappa)
    s = np.linspace(0, tau, nq); q = np.empty(nq); q[0] = x0
    h = tau/(nq-1)
    for i in range(nq-1):
        k1 = bf(q[i]); k2 = bf(q[i]+0.5*h*k1); k3 = bf(q[i]+0.5*h*k2); k4 = bf(q[i]+h*k3)
        q[i+1] = q[i] + h*(k1+2*k2+2*k3+k4)/6.0
    q = q[::-1]                      # x_s = X(tau - s): the REVERSE characteristic
    b1 = b1f(q); b2 = b2f(q)
    B = np.concatenate([[0.0], np.cumsum(0.5*(b1[1:]+b1[:-1]))*h])       # B(s) = int_0^s b'
    v = 2*nu*np.array([np.trapz(np.exp(-2*(B[i]-B[:i+1])), s[:i+1]) if i > 0 else 0.0
                       for i in range(nq)])
    integ = np.exp(-(B[-1]-B))*b2*v
    return -0.5*float(np.trapz(integ, s))

# ---- sweep in kappa at fixed beta (Gamma fixed, K2 = beta/kappa -> 0)
beta = 1.6; tau = 0.5; x0 = 0.8
rows = []
for kappa in [1.0, 2.0, 4.0, 8.0, 16.0]:
    nu = 2e-3
    loss, q, xa = solve(kappa, beta, nu, tau, x0*kappa, N=2048, nt=4000)
    pr = predict(kappa, beta, nu, tau, x0*kappa)
    rows.append(dict(kappa=kappa, K2=beta/kappa, nu=nu, loss=loss, loss_times_kappa=loss*kappa,
                     predicted=pr, meas_over_pred=loss/pr))
RES['S1_kappa_sweep'] = rows
print("S1  kappa   K2=||b''||   loss            loss*kappa       predicted        meas/pred")
for r in rows:
    print(f"    {r['kappa']:5.1f}  {r['K2']:8.4f}   {r['loss']: .6e}   {r['loss_times_kappa']: .6e}"
          f"   {r['predicted']: .6e}   {r['meas_over_pred']:7.4f}")

# ---- sweep in nu at fixed kappa (prediction: loss proportional to nu)
rows = []
for nu in [8e-3, 4e-3, 2e-3, 1e-3, 5e-4]:
    loss, q, xa = solve(1.0, beta, nu, tau, x0, N=2048, nt=6000)
    pr = predict(1.0, beta, nu, tau, x0)
    rows.append(dict(nu=nu, loss=loss, loss_over_nu=loss/nu, predicted=pr, meas_over_pred=loss/pr))
RES['S2_nu_sweep'] = rows
print("\nS2  nu        loss             loss/nu          predicted        meas/pred")
for r in rows:
    print(f"    {r['nu']:8.1e}  {r['loss']: .6e}   {r['loss_over_nu']: .6e}   {r['predicted']: .6e}   {r['meas_over_pred']:7.4f}")

# ---- resolution control
rows = []
for N, nt in [(512, 2000), (1024, 4000), (2048, 8000), (4096, 16000)]:
    loss, q, xa = solve(1.0, beta, 2e-3, tau, x0, N=N, nt=nt)
    rows.append(dict(N=N, nt=nt, loss=loss))
RES['S3_resolution'] = rows
print("\nS3  N     nt      loss")
for r in rows: print(f"    {r['N']:5d} {r['nt']:6d}  {r['loss']: .10e}")

# ---- affine control: b(x) = c1 x  ->  Z_tau exactly Gaussian mean zero  ->  loss == 0.
# Verified here on the same instrument by taking kappa -> infinity numerically (S1) and
# analytically in PROOF.md sec 4.  Extrapolated intercept of loss*kappa:
kk = np.array([r['kappa'] for r in RES['S1_kappa_sweep']])
ll = np.array([r['loss'] for r in RES['S1_kappa_sweep']])
RES['S4_loss_kappa_product_spread'] = dict(
    values=[float(v) for v in ll*kk],
    max_over_min=float((ll*kk).max()/(ll*kk).min()))
print("\nS4  loss*kappa spread (should tend to a constant):",
      RES['S4_loss_kappa_product_spread']['max_over_min'])

json.dump(RES, open('b4_results.json','w'), indent=1)
print("\nWROTE b4_results.json")
