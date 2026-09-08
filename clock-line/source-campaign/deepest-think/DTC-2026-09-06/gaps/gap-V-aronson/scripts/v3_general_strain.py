"""
v3 -- the two numerical refuters, run on a NON-uniform axisymmetric incompressible strain
(so that grad_alpha log J is NOT zero and the Lagrangian drift is genuinely present).

Field:   psi1 = -a0 z (1 + eps sin(k r))
         u^r = -r d_z psi1 = a0 r (1 + eps sin k r)          a = u^r/r
         u^z = 2 psi1 + r d_r psi1 = -2 a0 z (1+eps sin kr) - a0 z r eps k cos(kr)
         (div3 u = 0 exactly -- verified in v1 A10)

REFUTER R1 (locality).  Solve  d_t v + u.grad v = nu Delta5 v  with v0 a compactly supported
bump whose support is at distance >= d from the material label alpha0.  Read v at the material
point X(tau,alpha0).  Claim under test:
        |v(X(tau),tau)| <= ||v0||_inf * C exp( -d^2/(C' nu tau) ),   C' = 4 e^{2c}.
KILL: if the measured value exceeds ||v0||_inf * exp(-d^2/(4 e^{2c} nu tau)) by more than a
constant factor that grows with d, the claimed C' is wrong.

REFUTER R2 (first moment).  Monte-Carlo the reversed-time SDE
        dZ_s = -b(Z_s, tau-s) ds + sqrt(2 nu) dW_s ,  Z_0 = X(tau,alpha0)   (in R^5)
so that  eta(X(tau),tau) = E[eta_0(Z_tau)]  and  m := E[Z_tau] - alpha0  is exactly the first
moment of the Lagrangian kernel.  Two competing predictions:
   (A) Aronson class only:      |m| <~ c sqrt(nu tau)          [attained by layered media]
   (B) smooth flow:             |m| <~ ||grad^2 b||_inf nu tau^2
Common random numbers against the pure-strain control (whose m is EXACTLY 0, v1 A11) are used
so the estimator has small variance.
"""
import numpy as np, math, json
from scipy.special import ive

OUT = {}
def rec(k, v):
    OUT[k] = v
    print(f"{k} = {v}")

a0, eps, kk = 0.5, 0.35, 4.0
def a_of(r):      return a0*(1.0 + eps*np.sin(kk*r))
def ur_of(r, z):  return r*a_of(r)
def uz_of(r, z):  return -2*a0*z*(1.0+eps*np.sin(kk*r)) - a0*z*r*eps*kk*np.cos(kk*r)

def flow_rz(r, z, t0, t1, n=4000):
    """RK4 for (r,z) under (ur,uz); works for t1<t0 too."""
    h = (t1-t0)/n
    for _ in range(n):
        k1r, k1z = ur_of(r, z), uz_of(r, z)
        k2r, k2z = ur_of(r+.5*h*k1r, z+.5*h*k1z), uz_of(r+.5*h*k1r, z+.5*h*k1z)
        k3r, k3z = ur_of(r+.5*h*k2r, z+.5*h*k2z), uz_of(r+.5*h*k2r, z+.5*h*k2z)
        k4r, k4z = ur_of(r+h*k3r, z+h*k3z), uz_of(r+h*k3r, z+h*k3z)
        r = r + h*(k1r+2*k2r+2*k3r+k4r)/6.0
        z = z + h*(k1z+2*k2z+2*k3z+k4z)/6.0
    return r, z

def gradb_op_norm(r, z):
    """||grad5 b||_op = ||grad3 u||_op = max(|a|, ||2x2 block||)  (v1 A5)."""
    aa = a_of(r)
    dur_dr = a0*(1+eps*np.sin(kk*r)) + a0*r*eps*kk*np.cos(kk*r)
    dur_dz = 0.0*r
    duz_dr = -2*a0*z*eps*kk*np.cos(kk*r) - a0*z*(eps*kk*np.cos(kk*r) - r*eps*kk*kk*np.sin(kk*r))
    duz_dz = -2*a0*(1+eps*np.sin(kk*r)) - a0*r*eps*kk*np.cos(kk*r)
    B = np.array([[dur_dr, dur_dz], [duz_dr, duz_dz]])
    return max(abs(aa), np.linalg.norm(B, 2))

# ------------------------------------------------------------------ R1: locality PDE
def solve_rz(f0, nu, tau, rlo, rhi, zlo, zhi, Nr, Nz, cfl=0.20):
    r = np.linspace(rlo, rhi, Nr); z = np.linspace(zlo, zhi, Nz)
    hr = r[1]-r[0]; hz = z[1]-z[0]
    R, Z = np.meshgrid(r, z, indexing='ij')
    e = f0(R, Z).astype(float)
    UR = ur_of(R, Z); UZ = uz_of(R, Z)
    umax = max(abs(UR).max(), abs(UZ).max())
    dt = min(cfl*min(hr, hz)**2/(2*nu), 0.4*min(hr, hz)/max(umax, 1e-12))
    n = int(math.ceil(tau/dt)); dt = tau/n
    def L(e):
        er = np.zeros_like(e); ez = np.zeros_like(e)
        err = np.zeros_like(e); ezz = np.zeros_like(e)
        er[1:-1, :] = (e[2:, :]-e[:-2, :])/(2*hr)
        ez[:, 1:-1] = (e[:, 2:]-e[:, :-2])/(2*hz)
        err[1:-1, :] = (e[2:, :]-2*e[1:-1, :]+e[:-2, :])/hr**2
        ezz[:, 1:-1] = (e[:, 2:]-2*e[:, 1:-1]+e[:, :-2])/hz**2
        de = -UR*er - UZ*ez + nu*(err + 3.0/R*er + ezz)
        de[0, :] = de[-1, :] = 0; de[:, 0] = de[:, -1] = 0
        return de
    for _ in range(n):
        k1 = L(e); k2 = L(e+dt*k1); e = e + .5*dt*(k1+k2)
        e[0, :] = e[-1, :] = 0; e[:, 0] = e[:, -1] = 0
    return r, z, e, n

def bilerp(r, z, e, rq, zq):
    i = min(max(np.searchsorted(r, rq)-1, 0), len(r)-2)
    j = min(max(np.searchsorted(z, zq)-1, 0), len(z)-2)
    tr = (rq-r[i])/(r[i+1]-r[i]); tz = (zq-z[j])/(z[j+1]-z[j])
    return ((1-tr)*(1-tz)*e[i, j]+tr*(1-tz)*e[i+1, j]+(1-tr)*tz*e[i, j+1]+tr*tz*e[i+1, j+1])

def bump(rc, zc, s):
    def f(R, Z):
        q = ((R-rc)**2+(Z-zc)**2)/s**2
        return np.where(q < 1.0, np.exp(-1.0/np.maximum(1e-16, 1.0-q))*math.e, 0.0)
    return f

tau = 0.4; nu = 3.0e-3
r0, z0 = 1.0, 0.35
G = max(gradb_op_norm(rr, zz) for rr in np.linspace(0.35, 2.4, 400) for zz in [-1.2, -0.5, 0.0, 0.5, 1.2])
c = G*tau
rT, zT = flow_rz(r0, z0, 0.0, tau)
rec("R1_params", dict(a0=a0, eps=eps, k=kk, tau=tau, nu=nu, G=float(G), c=float(c),
                      alpha0=[r0, z0], Xtau=[float(rT), float(zT)]))

rows = []
for D in [0.10, 0.13, 0.16, 0.19, 0.22, 0.25]:
    s = 0.035
    d = D - s
    f0 = bump(r0+D, z0, s)
    r, z, e, n = solve_rz(f0, nu, tau, 0.30, 2.6, -1.35, 1.35, 641, 641)
    val = float(abs(bilerp(r, z, e, rT, zT)))
    ratio = d*d/(nu*tau)
    claim = math.exp(-ratio/(4*math.exp(2*c)))
    rows.append(dict(D=D, d=d, d2_over_nutau=ratio, measured=val, claim=claim,
                     measured_over_claim=val/claim, nsteps=n))
    print("   ", rows[-1])
rec("R1_locality_rows", rows)
# fitted exponential rate:  log|v| ~ -ratio/Cfit
xs = np.array([q['d2_over_nutau'] for q in rows]); ys = np.log(np.array([q['measured'] for q in rows]))
sl, ic = np.polyfit(xs, ys, 1)
rec("R1_fitted_Cprime", float(-1.0/sl))
rec("R1_claimed_Cprime_4e2c", float(4*math.exp(2*c)))
rec("R1_KILL_fires", bool(-1.0/sl > 4*math.exp(2*c)))
# grid convergence at one D
cg = []
for N in [401, 521, 641, 801]:
    r, z, e, n = solve_rz(bump(r0+0.16, z0, 0.035), nu, tau, 0.30, 2.6, -1.35, 1.35, N, N)
    cg.append(dict(N=N, val=float(abs(bilerp(r, z, e, rT, zT)))))
    print("   conv", cg[-1])
rec("R1_grid_convergence", cg)

# ------------------------------------------------------------------ R2: first moment MC
def b5(P, t):
    """P: (N,5).  returns drift b (5D lift of u)."""
    y = P[:, :4]; zz = P[:, 4]
    rr = np.sqrt(np.maximum((y*y).sum(1), 1e-300))
    urv = ur_of(rr, zz); uzv = uz_of(rr, zz)
    out = np.empty_like(P)
    out[:, :4] = (urv/rr)[:, None]*y
    out[:, 4] = uzv
    return out

def b5_pure(P, t):
    y = P[:, :4]; zz = P[:, 4]
    out = np.empty_like(P)
    out[:, :4] = a0*y
    out[:, 4] = -2*a0*zz
    return out

def mc_first_moment(nu, tau, nsteps, nsamp, seed):
    rng = np.random.default_rng(seed)
    # start points: the two flows' own X(tau,alpha0)
    rT1, zT1 = flow_rz(r0, z0, 0.0, tau)
    lam = math.exp(a0*tau); rT2, zT2 = r0*lam, z0/lam**2
    P1 = np.zeros((nsamp, 5)); P1[:, 0] = rT1; P1[:, 4] = zT1
    P2 = np.zeros((nsamp, 5)); P2[:, 0] = rT2; P2[:, 4] = zT2
    h = tau/nsteps; sq = math.sqrt(2*nu*h)
    for i in range(nsteps):
        dW = rng.standard_normal((nsamp, 5))*sq          # common random numbers
        s = i*h
        P1 = P1 - b5(P1, tau-s)*h + dW
        P2 = P2 - b5_pure(P2, tau-s)*h + dW
    # Lagrangian labels: for the pure strain E[Z2_tau] = alpha0 exactly.
    m1 = P1.mean(0); m2 = P2.mean(0)
    diff = P1 - P2
    md = diff.mean(0); sd = diff.std(0)/math.sqrt(nsamp)
    # alpha0 in R^5 is (r0,0,0,0,z0); the pure-strain control's exact mean is that point.
    alpha0 = np.array([r0, 0, 0, 0, z0])
    m = md + (m2 - alpha0)                               # = E[Z1] - alpha0
    return dict(m=m.tolist(), m_norm=float(np.linalg.norm(m)),
                se=float(np.linalg.norm(sd)),
                control_pure_m=(m2-alpha0).tolist(),
                control_pure_norm=float(np.linalg.norm(m2-alpha0)),
                se_control=float(np.linalg.norm(P2.std(0)/math.sqrt(nsamp))))

# ||grad^2 b||: finite differences of the 2x2 block / a
def grad2_scale():
    hh = 1e-4
    vals = []
    for rr in np.linspace(0.6, 1.6, 60):
        for zz in [-0.5, 0.0, 0.5]:
            f = lambda R, Zq: ur_of(R, Zq)
            d2 = (f(rr+hh, zz) - 2*f(rr, zz) + f(rr-hh, zz))/hh**2
            vals.append(abs(d2))
    return max(vals)

G2 = grad2_scale()
rec("R2_grad2b_scale", float(G2))
res = []
for nu_ in [3.0e-3, 1.0e-3, 3.0e-4]:
    out = mc_first_moment(nu_, tau, 400, 400000, seed=12345)
    predA = c*math.sqrt(nu_*tau)                       # Aronson-class bound
    predB = G2*nu_*tau*tau                             # smooth-flow prediction
    out.update(nu=nu_, predA_c_sqrt_nutau=predA, predB_grad2b_nu_tau2=predB,
               m_over_predA=out['m_norm']/predA, m_over_predB=out['m_norm']/predB,
               sqrt_nutau=math.sqrt(nu_*tau))
    res.append(out); print("   ", {k: out[k] for k in
        ('nu', 'm_norm', 'se', 'control_pure_norm', 'se_control',
         'predA_c_sqrt_nutau', 'predB_grad2b_nu_tau2', 'm_over_predA', 'm_over_predB')})
rec("R2_first_moment", res)

with open(__file__.replace('v3_general_strain.py', 'v3_results.json'), 'w') as fh:
    json.dump(OUT, fh, indent=1, default=str)
print("\nDONE v3")
