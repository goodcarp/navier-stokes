"""
v2_flowmap_bounds.py -- gap-V-aronson.

The Lagrangian-frame coefficient bounds, and their sharpness.

Claim (Grönwall).  J solves  d_t J = (grad_5 b)(X(t),t) J,  J(0)=I.  If
Gamma := sup_{[0,tau]} || grad_5 b(.,t) ||_{op,inf} and Gamma*tau <= c then every
singular value of J lies in [e^{-c}, e^{c}], hence
    e^{-2c} I  <=  G = (J^T J)^{-1}  <=  e^{2c} I ,     e^{-nc} <= det J <= e^{nc}   (n=5).

Tested by direct numerical integration of d_t J = B(t) J for
 (i) random time-dependent B(t) with ||B||_op == Gamma exactly (worst case search),
 (ii) the physical field: the 5D lift of the uniform axisymmetric strain
      u = (a r, -2 a z), grad_5 b = diag(a,a,a,a,-2a), ||.||_op = 2a, div_5 = 2a.

Falsification rule (pre-declared): the claim is refuted if any run produces a
singular value outside [e^{-c}-1e-9, e^{c}+1e-9] or det J outside
[e^{-5c}-1e-9, e^{5c}+1e-9].
"""
import json, numpy as np
from scipy.integrate import solve_ivp

HERE = "~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/gaps/gap-V-aronson"
rng = np.random.default_rng(20260906)
n = 5
OUT = {}

def integrate(Bfun, tau, n=5):
    def rhs(t, Jf):
        J = Jf.reshape(n, n)
        return (Bfun(t) @ J).ravel()
    sol = solve_ivp(rhs, (0, tau), np.eye(n).ravel(), rtol=1e-12, atol=1e-14, dense_output=True)
    return sol.y[:, -1].reshape(n, n)

# ---- (i) random B(t) normalised to operator norm Gamma ----
rows = []
worst_hi = -np.inf; worst_lo = np.inf; worst_det_hi = -np.inf; worst_det_lo = np.inf
for trial in range(400):
    Gamma = float(rng.uniform(0.2, 5.0))
    tau = float(rng.uniform(0.05, 2.0))
    c = Gamma*tau
    K = 3
    As = []
    for k in range(K):
        A = rng.normal(size=(n, n))
        A = A/np.linalg.norm(A, 2)*Gamma          # operator norm exactly Gamma
        As.append(A)
    def Bfun(t, As=As, tau=tau, K=K):
        k = min(int(t/tau*K), K-1)
        return As[k]
    J = integrate(Bfun, tau)
    sv = np.linalg.svd(J, compute_uv=False)
    det = np.linalg.det(J)
    worst_hi = max(worst_hi, sv.max()/np.exp(c))
    worst_lo = min(worst_lo, sv.min()*np.exp(c))
    worst_det_hi = max(worst_det_hi, det/np.exp(n*c))
    worst_det_lo = min(worst_det_lo, det*np.exp(n*c))
    G = np.linalg.inv(J.T@J)
    ev = np.linalg.eigvalsh(G)
    rows.append(dict(Gamma=Gamma, tau=tau, c=c, smax=sv.max(), smin=sv.min(),
                     det=det, Gmax=ev.max(), Gmin=ev.min(),
                     ok=bool(sv.max() <= np.exp(c)+1e-9 and sv.min() >= np.exp(-c)-1e-9
                             and abs(det) <= np.exp(n*c)+1e-9 and abs(det) >= np.exp(-n*c)-1e-9
                             and ev.max() <= np.exp(2*c)+1e-9 and ev.min() >= np.exp(-2*c)-1e-9)))
n_ok = sum(r['ok'] for r in rows)
print("(i) random B(t), 400 trials:  all bounds hold in %d/400" % n_ok)
print("    worst  smax/e^c = %.9f   smin*e^c = %.9f   det/e^{5c} = %.9f   det*e^{5c} = %.9f"
      % (worst_hi, worst_lo, worst_det_hi, worst_det_lo))
OUT['i_n_ok'] = n_ok; OUT['i_trials'] = len(rows)
OUT['i_worst_smax_over_expc'] = worst_hi
OUT['i_worst_smin_times_expc'] = worst_lo
OUT['i_worst_det_over_exp5c'] = worst_det_hi
OUT['i_worst_det_times_exp5c'] = worst_det_lo

# ---- sharpness: constant B = diag(Gamma,...) gives smax = e^{c} exactly ----
Gamma, tau = 1.3, 0.7; c = Gamma*tau
B0 = np.diag([Gamma, -Gamma, Gamma, -Gamma, Gamma])
J = integrate(lambda t: B0, tau)
sv = np.linalg.svd(J, compute_uv=False)
print("    sharpness (constant diagonal B): smax = %.12f  vs e^c = %.12f ; smin = %.12f vs e^{-c} = %.12f"
      % (sv.max(), np.exp(c), sv.min(), np.exp(-c)))
OUT['sharp_smax'] = sv.max(); OUT['sharp_expc'] = float(np.exp(c))
OUT['sharp_smin'] = sv.min(); OUT['sharp_expmc'] = float(np.exp(-c))

# ---- (ii) the physical lift: uniform axisymmetric strain ----
# grad_5 b = diag(a,a,a,a,-2a).  Gamma = 2a.  det J = exp(int 2a) = the 5D volume factor.
print("\n(ii) 5D lift of u=(a r,-2 a z):  grad_5 b = diag(a,a,a,a,-2a)")
res2 = []
for a0 in (0.5, 1.0, 2.0):
    for tau in (0.2, 0.5, 0.8):
        Gamma = 2*a0; c = Gamma*tau
        J = integrate(lambda t: np.diag([a0, a0, a0, a0, -2*a0]), tau)
        sv = np.linalg.svd(J, compute_uv=False); det = np.linalg.det(J)
        lam = np.exp(a0*tau)                       # T_lambda strain factor
        res2.append(dict(a=a0, tau=tau, c=c, lam=lam, smax=sv.max(), smin=sv.min(),
                         det=det, det_pred=float(np.exp(2*a0*tau)),
                         smax_over_expc=sv.max()/np.exp(c), det_over_exp5c=det/np.exp(5*c)))
        print("    a=%.1f tau=%.1f  c=%.2f  lambda=%.4f  smax=%.6f (e^c=%.6f)  smin=%.6f  det=%.6f (e^{2 a tau}=%.6f)"
              % (a0, tau, c, lam, sv.max(), np.exp(c), sv.min(), det, np.exp(2*a0*tau)))
OUT['ii'] = res2
# the weight really is e^{int div_5 b} = e^{2 a tau}: check
OUT['ii_det_equals_exp_int_div5'] = max(abs(r['det']-r['det_pred']) for r in res2)
print("    max |det J - exp(int div_5 b)| = %.3e   (Liouville, so the weight IS the 5D volume factor)"
      % OUT['ii_det_equals_exp_int_div5'])
# and how far det J is from its worst-case envelope e^{5c}:
print("    det J / e^{5c} at c=%.2f: %.3e  -> the envelope e^{5c} is very lossy for the physical field"
      % (res2[-1]['c'], res2[-1]['det_over_exp5c']))

with open(HERE+"/v2_results.json", "w") as fh:
    json.dump(OUT, fh, indent=1, default=float)
print("\nwrote v2_results.json")
