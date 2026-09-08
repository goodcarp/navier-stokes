"""
v4_falsify.py -- gap-V-aronson.  Refuter runs for THEOREM V.1.

PRE-DECLARED FALSIFICATION RULES (fixed before any number was read):
 R1  (SDE, n=5) the theorem is REFUTED if the Monte-Carlo estimate of
     P(|Y_tau - x0| >= d), minus 3 standard errors, exceeds min(B1,B2) at any tested d.
 R2  (PDE, n=1) the theorem is REFUTED if the measured
     |eta_1(X(tau),tau) - eta_2(X(tau),tau)| / ||eta_1(.,0)-eta_2(.,0)||_inf
     exceeds the n=1 bound at any tested d.
 R3  the material frame is NECESSARY: the same PDE quantity measured at the FIXED
     Eulerian point x0 instead of at X(tau) must VIOLATE the bound (control must fire),
     otherwise the numerics are non-diagnostic.

Drift used (n=5): the 5D lift of a time-dependent axisymmetric strain plus a
genuinely nonlinear, compressible perturbation, so that div_5 b != 0 and b is not
linear:   b(x,t) = a(t)*(y, -2z) + kappa * g(x),  g bounded C^1.
Gamma is MEASURED (max operator norm of grad b over the sampled region), not assumed.
"""
import json, numpy as np
HERE = "~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/gaps/gap-V-aronson"
OUT = {}
rng = np.random.default_rng(1789)
n5 = 5

# ---------- bounds (same formulas as v3) ----------
_ES = np.linspace(1e-3, 1.4, 4001)
def B1(c, q, n): return 2*n*np.exp(-q*np.exp(-2*c)/(4*n))
def B2(c, q, n):
    k = c/(2*np.exp(2*c)*(np.exp(2*c)-1))
    logv = n*np.log(1+2/_ES) - (1-_ES*_ES/2)**2*k*q
    i = int(np.argmin(logv)); return float(np.exp(logv[i]))
def BOUND(c, q, n): return min(B1(c,q,n), B2(c,q,n))
def SHARP(c, q):    # not proved: exponent q*c/(2(e^{2c}-1)), attained by the pure strain
    return np.exp(-q*c/(2*(np.exp(2*c)-1)))

# =====================================================================
# TEST A -- n=5 SDE
# =====================================================================
a0, kap = 1.0, 0.35
def bfield(x, t):
    y = x[..., :4]; z = x[..., 4:5]
    lin = np.concatenate([a0*y, -2*a0*z], axis=-1)
    g = np.concatenate([np.sin(x[..., 1:2]+0.5*x[..., 4:5]),
                        np.cos(x[..., 2:3]),
                        np.sin(x[..., 3:4]-x[..., 0:1]),
                        np.cos(0.7*x[..., 4:5]),
                        np.sin(x[..., 0:1]+x[..., 2:3])], axis=-1)
    return lin + kap*g
def gradb(x, t, h=1e-5):
    J = np.zeros(x.shape[:-1]+(5,5))
    for j in range(5):
        e = np.zeros(5); e[j] = h
        J[..., :, j] = (bfield(x+e, t)-bfield(x-e, t))/(2*h)
    return J

tau = 0.4
x0 = np.array([0.9, 0.2, -0.3, 0.15, 0.6])
# deterministic forward trajectory
def rk4(x, t0, t1, N):
    h = (t1-t0)/N; t = t0
    for _ in range(N):
        k1 = bfield(x, t); k2 = bfield(x+h/2*k1, t+h/2)
        k3 = bfield(x+h/2*k2, t+h/2); k4 = bfield(x+h*k3, t+h)
        x = x + h/6*(k1+2*k2+2*k3+k4); t += h
    return x
Xtau = rk4(x0.copy(), 0.0, tau, 4000)
# Gamma measured along a tube around the trajectory
pts = []
for s in np.linspace(0, tau, 40):
    xs = rk4(x0.copy(), 0.0, s, 400) if s > 0 else x0.copy()
    pts.append(xs + 0.6*rng.normal(size=(200, 5)))
pts = np.concatenate(pts, axis=0)
Gam = float(np.max(np.linalg.norm(gradb(pts, 0.0), ord=2, axis=(-2, -1))))
c = Gam*tau
print("TEST A (n=5).  measured Gamma = %.5f   tau = %.3f   c = Gamma*tau = %.5f" % (Gam, tau, c))
print("               div_5 b at x0 = %.5f   (nonzero: the drift is compressible)"
      % np.trace(gradb(x0[None, :], 0.0)[0]))
OUT['A_Gamma'] = Gam; OUT['A_tau'] = tau; OUT['A_c'] = c
OUT['A_div5b_x0'] = float(np.trace(gradb(x0[None, :], 0.0)[0]))

def sim(nu, npath, nstep, seed):
    r = np.default_rng if False else np.random.default_rng(seed)
    h = tau/nstep
    Y = np.repeat(Xtau[None, :], npath, axis=0)
    t = 0.0
    for k in range(nstep):
        dW = r.normal(scale=np.sqrt(h), size=(npath, 5))
        Y = Y - bfield(Y, tau-t)*h + np.sqrt(2*nu)*dW
        t += h
    return Y

rowsA = []
for nu in (2e-3, 5e-4):
    Y = sim(nu, 250000, 600, 4242 + int(1e5*nu))
    Zn = np.linalg.norm(Y - x0[None, :], axis=1)
    print("  nu=%.1e  sqrt(nu*tau)=%.4e   mean|Z|=%.4e   rms|Z|=%.4e" %
          (nu, np.sqrt(nu*tau), Zn.mean(), np.sqrt((Zn**2).mean())))
    for frac in (4.0, 6.0, 8.0, 10.0, 12.0):
        d = frac*np.sqrt(nu*tau)
        q = d*d/(nu*tau)
        p = float((Zn >= d).mean()); se = np.sqrt(max(p, 1e-12)*(1-p)/len(Zn))
        b = BOUND(c, q, n5); sh = SHARP(c, q)
        viol = (p - 3*se) > b
        rowsA.append(dict(nu=nu, d_over_sqrt=frac, q=q, p_emp=p, se=se,
                          bound=b, sharp_rate=sh, violates=bool(viol)))
        print("     d/sqrt(nu tau)=%4.1f  q=%6.1f   P_emp=%.3e +-%.1e   BOUND=%.3e   (sharp-rate ref %.2e)   %s"
              % (frac, q, p, se, b, sh, "VIOLATION" if viol else "ok"))
OUT['A_rows'] = rowsA
OUT['A_any_violation'] = any(r['violates'] for r in rowsA)

# fitted empirical rate:  -log P = rate * q
fit = {}
for nu in (2e-3, 5e-4):
    sel = [r for r in rowsA if r['nu'] == nu and 1e-5 < r['p_emp'] < 0.2]
    if len(sel) >= 2:
        qq = np.array([r['q'] for r in sel]); pp = np.array([r['p_emp'] for r in sel])
        A = np.vstack([qq, np.ones_like(qq)]).T
        sl, ic = np.linalg.lstsq(A, -np.log(pp), rcond=None)[0]
        fit[str(nu)] = dict(rate=float(sl), intercept=float(ic))
        print("  fitted empirical rate at nu=%.1e :  %.5f   (proved B1 rate %.5f, B2 rate %.5f, sharp %.5f)"
              % (nu, sl, np.exp(-2*c)/(4*n5), c/(2*np.exp(2*c)*(np.exp(2*c)-1)), c/(2*(np.exp(2*c)-1))))
OUT['A_fit'] = fit
OUT['A_rate_B1'] = float(np.exp(-2*c)/(4*n5))
OUT['A_rate_B2'] = float(c/(2*np.exp(2*c)*(np.exp(2*c)-1)))
OUT['A_rate_sharp'] = float(c/(2*(np.exp(2*c)-1)))

# =====================================================================
# TEST B -- n=1 PDE, compressible drift, end-to-end
# =====================================================================
print("\nTEST B (n=1 PDE).  d_t eta + b(x) d_x eta = nu d_xx eta,  b' != 0.")
Gb = 1.6
bx = lambda x: Gb*np.sin(x)          # Lipschitz constant exactly Gb; b' = Gb cos x != 0
Nx, Lx = 24001, 30.0
xg = np.linspace(-Lx, Lx, Nx); hx = xg[1]-xg[0]
def solve(eta0, nu, tauB, cfl=0.20):
    dt = cfl*min(hx*hx/(2*nu), hx/max(abs(bx(xg)).max(), 1e-12))
    ns = int(np.ceil(tauB/dt)); dt = tauB/ns
    e = eta0.copy(); bb = bx(xg)
    for _ in range(ns):
        # upwind advection + centred diffusion (advective form, so no div term)
        dp = np.zeros_like(e); dm = np.zeros_like(e)
        dp[:-1] = (e[1:]-e[:-1])/hx; dm[1:] = (e[1:]-e[:-1])/hx
        adv = np.where(bb > 0, bb*dm, bb*dp)
        lap = np.zeros_like(e); lap[1:-1] = (e[2:]-2*e[1:-1]+e[:-2])/hx**2
        e = e - dt*adv + dt*nu*lap
        e[0] = e[1]; e[-1] = e[-2]
    return e, ns
tauB = 0.5; cB = Gb*tauB
x0B = 0.8
# forward trajectory of x0B
xt = x0B
NS = 20000; hstep = tauB/NS
for _ in range(NS):
    k1 = bx(xt); k2 = bx(xt+hstep/2*k1); k3 = bx(xt+hstep/2*k2); k4 = bx(xt+hstep*k3)
    xt = xt + hstep/6*(k1+2*k2+2*k3+k4)
XtauB = xt
print("  Gamma=%.2f tau=%.2f c=%.2f ;  x0=%.3f -> X(tau)=%.4f  (drift displacement %.4f)"
      % (Gb, tauB, cB, x0B, XtauB, XtauB-x0B))
OUT['B_c'] = cB; OUT['B_x0'] = x0B; OUT['B_Xtau'] = float(XtauB)
rowsB = []
for nu in (2e-3, 5e-4):
    for frac in (4.0, 6.0, 8.0):
        d = frac*np.sqrt(nu*tauB)
        base = np.exp(-((xg-x0B)/3.0)**2)*np.sin(xg)      # smooth reference datum
        pert = np.where(np.abs(xg-x0B) >= d, 1.0, 0.0)     # perturbation supported at distance >= d
        e1, ns = solve(base, nu, tauB); e2, _ = solve(base+pert, nu, tauB)
        dev_mat = float(abs(np.interp(XtauB, xg, e2-e1)))
        dev_eul = float(abs(np.interp(x0B,  xg, e2-e1)))
        q = d*d/(nu*tauB); bnd = BOUND(cB, q, 1)
        rowsB.append(dict(nu=nu, frac=frac, q=q, dev_material=dev_mat, dev_eulerian=dev_eul,
                          bound=bnd, violates=bool(dev_mat > bnd),
                          control_fires=bool(dev_eul > bnd)))
        print("   nu=%.0e d/sqrt(nu tau)=%4.1f q=%6.1f | material dev=%.3e  BOUND=%.3e  %s | Eulerian dev=%.3e %s"
              % (nu, frac, q, dev_mat, bnd, "VIOLATION" if dev_mat > bnd else "ok",
                 dev_eul, "(control FIRES)" if dev_eul > bnd else "(control silent)"))
OUT['B_rows'] = rowsB
OUT['B_any_violation'] = any(r['violates'] for r in rowsB)
OUT['B_control_fires'] = any(r['control_fires'] for r in rowsB)

print("\nSUMMARY  R1 violated: %s   R2 violated: %s   R3 control fires: %s"
      % (OUT['A_any_violation'], OUT['B_any_violation'], OUT['B_control_fires']))
with open(HERE+"/v4_results.json","w") as fh: json.dump(OUT, fh, indent=1, default=float)
print("wrote v4_results.json")
