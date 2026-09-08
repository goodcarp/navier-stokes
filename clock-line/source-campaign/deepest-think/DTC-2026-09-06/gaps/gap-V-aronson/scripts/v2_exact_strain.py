"""
v2 -- the pure axisymmetric strain: an EXACTLY SOLVABLE instance of the Lagrangian reduction,
used (i) to exhibit the constants C, C' explicitly and (ii) to check the reduction itself
against a completely independent direct solve of the 5D-lift PDE in the (r,z) half-plane.

Pure strain    u^r = a0 r ,  u^z = -2 a0 z    (div3 u = 0).
Lift           b = (a0 y, -2 a0 z) on R^4 x R ,  div5 b = 4a0 - 2a0 = 2a0 = 2a   (compressible).
Flow           X5(t,alpha) = diag(e^{a0 t} I4, e^{-2a0 t}) alpha ,  A = that,  J = e^{2 a0 t}
               -> J is SPATIALLY CONSTANT, so grad_alpha log J = 0 and the Lagrangian equation
                  carries NO drift:
                     d_t etatilde = nu ( e^{-2a0 t} Delta_y + e^{+4a0 t} d_zz ) etatilde
               whose solution is convolution with the Gaussian of covariance
                     Sigma(t) = 2 nu diag( S1 I4 , S2 ),
                     S1 = int_0^t e^{-2a0 s} ds ,  S2 = int_0^t e^{+4a0 s} ds .
G := ||grad5 b||_op = 2 a0 ,  c := G tau .   Then  Sigma_max = 2 nu S2 <= 2 nu tau e^{2c}
and Sigma_min = 2 nu S1 >= 2 nu tau e^{-2c}: the whole c-dependence, explicitly.
"""
import numpy as np, json, math
from scipy.special import ive

OUT = {}
def rec(k, v):
    OUT[k] = v
    print(f"{k} = {v}")

# ---------------------------------------------------------------- exact Lagrangian kernel
def S12(a0, tau):
    S1 = (1.0 - math.exp(-2*a0*tau))/(2*a0)
    S2 = (math.exp(4*a0*tau) - 1.0)/(4*a0)
    return S1, S2

def lag_conv_axisym(f_rz, r0, z0, nu, a0, tau, rgrid, zgrid):
    """(G_Sigma * f)(alpha0) with alpha0 = (r0 e_1, z0) in R^5, f axisymmetric in R^4.
       4D isotropic part: variance s1 per component; 1D z part: variance s2.
       sphere average in R^4 of exp(u cos th) = 2 I_1(u)/u ."""
    S1, S2 = S12(a0, tau)
    s1 = 2*nu*S1
    s2 = 2*nu*S2
    R, Zz = np.meshgrid(rgrid, zgrid, indexing='ij')
    u = r0*R/s1
    # (2 pi s1)^{-2} exp(-(r0^2+R^2)/(2 s1)) * 2 I1(u)/u  ; use ive(1,u)=I1(u)e^{-u}
    with np.errstate(divide='ignore', invalid='ignore'):
        sph = np.where(u > 1e-12, 2*ive(1, u)/np.maximum(u, 1e-300), 1.0)
    ker4 = (2*np.pi*s1)**-2 * np.exp(-(r0**2 + R**2)/(2*s1) + u) * sph
    ker1 = (2*np.pi*s2)**-0.5 * np.exp(-(z0 - Zz)**2/(2*s2))
    w = 2*np.pi**2 * R**3                      # |S^3| r^3
    integ = f_rz(R, Zz) * ker4 * ker1 * w
    dr = rgrid[1]-rgrid[0]; dz = zgrid[1]-zgrid[0]
    return np.trapz(np.trapz(integ, dx=dz, axis=1), dx=dr)

# ---------------------------------------------------------------- direct (r,z) PDE solver
def solve_rz(f0, nu, a0, tau, rlo, rhi, zlo, zhi, Nr, Nz, cfl=0.20, report=None):
    """d_t eta + a0 r d_r eta - 2 a0 z d_z eta = nu (d_rr + (3/r) d_r + d_zz) eta,
       Dirichlet 0 on the box (data compactly supported well inside).  RK2, centred."""
    r = np.linspace(rlo, rhi, Nr); z = np.linspace(zlo, zhi, Nz)
    hr = r[1]-r[0]; hz = z[1]-z[0]
    R, Z = np.meshgrid(r, z, indexing='ij')
    e = f0(R, Z).astype(float)
    ur = a0*R; uz = -2*a0*Z
    umax = max(abs(ur).max(), abs(uz).max())
    dt = cfl*min(hr*hr, hz*hz)/(2*nu) if nu > 0 else 1e9
    dt = min(dt, 0.5*min(hr, hz)/max(umax, 1e-12))
    n = int(math.ceil(tau/dt)); dt = tau/n
    def L(e):
        de = np.zeros_like(e)
        er = np.zeros_like(e); ez = np.zeros_like(e)
        er[1:-1, :] = (e[2:, :] - e[:-2, :])/(2*hr)
        ez[:, 1:-1] = (e[:, 2:] - e[:, :-2])/(2*hz)
        err = np.zeros_like(e); ezz = np.zeros_like(e)
        err[1:-1, :] = (e[2:, :] - 2*e[1:-1, :] + e[:-2, :])/hr**2
        ezz[:, 1:-1] = (e[:, 2:] - 2*e[:, 1:-1] + e[:, :-2])/hz**2
        de = -ur*er - uz*ez + nu*(err + 3.0/R*er + ezz)
        de[0, :] = de[-1, :] = 0.0; de[:, 0] = de[:, -1] = 0.0
        return de
    for k in range(n):
        k1 = L(e); k2 = L(e + dt*k1)
        e = e + 0.5*dt*(k1 + k2)
        e[0, :] = e[-1, :] = 0.0; e[:, 0] = e[:, -1] = 0.0
    return r, z, e, n, dt

def bilerp(r, z, e, rq, zq):
    i = np.searchsorted(r, rq) - 1; j = np.searchsorted(z, zq) - 1
    i = min(max(i, 0), len(r)-2); j = min(max(j, 0), len(z)-2)
    tr = (rq-r[i])/(r[i+1]-r[i]); tz = (zq-z[j])/(z[j+1]-z[j])
    return ((1-tr)*(1-tz)*e[i, j] + tr*(1-tz)*e[i+1, j]
            + (1-tr)*tz*e[i, j+1] + tr*tz*e[i+1, j+1])

# ================================================================ PART 1: cross-check
a0 = 0.5; tau = 0.4; nu = 2.0e-3
G = 2*a0; c = G*tau
r0, z0 = 1.0, 0.3                                   # material label alpha0
lam = math.exp(a0*tau)
rT, zT = r0*lam, z0*lam**-2                         # X(tau, alpha0)
rec("PART1_params", dict(a0=a0, tau=tau, nu=nu, G=G, c=c, alpha0=[r0, z0], Xtau=[rT, zT]))
S1, S2 = S12(a0, tau)
rec("PART1_S1_S2_over_tau", [S1/tau, S2/tau])
rec("PART1_e_pm_2c", [math.exp(-2*c), math.exp(2*c)])
rec("PART1_S_within_e2c", bool(S1/tau >= math.exp(-2*c) - 1e-12 and S2/tau <= math.exp(2*c) + 1e-12))

# a compactly supported "feature" bump sitting at distance D from alpha0 at t=0
def make_bump(rc, zc, s):
    def f(R, Z):
        q = ((R-rc)**2 + (Z-zc)**2)/s**2
        return np.where(q < 1.0, np.exp(-1.0/np.maximum(1e-16, 1.0-q))*math.e, 0.0)
    return f

rows = []
for D in [0.10, 0.14, 0.18, 0.22]:
    s = 0.04
    f0 = make_bump(r0 + D, z0, s)
    d = D - s                                        # distance from alpha0 to supp v0
    ex = lag_conv_axisym(f0, r0, z0, nu, a0, tau,
                         np.linspace(0.2, 2.2, 1601), np.linspace(-1.2, 1.2, 1601))
    r, z, e, n, dt = solve_rz(f0, nu, a0, tau, 0.25, 2.6, -1.3, 1.3, 561, 561)
    pd = bilerp(r, z, e, rT, zT)
    rows.append(dict(D=D, d=d, exact=float(ex), pde=float(pd),
                     rel=float(abs(pd-ex)/max(abs(ex), 1e-300)), nsteps=n,
                     d2_over_nutau=d*d/(nu*tau)))
    print("   ", rows[-1])
rec("PART1_crosscheck", rows)

# grid refinement at one D
D = 0.14; s = 0.04; f0 = make_bump(r0+D, z0, s)
ex = lag_conv_axisym(f0, r0, z0, nu, a0, tau,
                     np.linspace(0.2, 2.2, 2001), np.linspace(-1.2, 1.2, 2001))
conv = []
for N in [281, 401, 561, 801]:
    r, z, e, n, dt = solve_rz(f0, nu, a0, tau, 0.25, 2.6, -1.3, 1.3, N, N)
    pd = bilerp(r, z, e, rT, zT)
    conv.append(dict(N=N, pde=float(pd), rel=float(abs(pd-ex)/abs(ex))))
    print("   conv", conv[-1])
rec("PART1_grid_convergence", dict(exact=float(ex), rows=conv))

# ================================================================ PART 2: the explicit bound
# exact tail of the Lagrangian Gaussian: |(G_Sigma * v0)(alpha0)| <= ||v0||_inf * P(|Xi| >= d),
# Xi ~ N(0, Sigma), Sigma = 2 nu diag(S1 I4, S2).  Since Sigma <= 2 nu tau e^{2c} I,
#     P(|Xi| >= d) <= P(|W| >= d),  W ~ N(0, 2 nu tau e^{2c} I5)
#                  = Q5( d / sqrt(2 nu tau e^{2c}) )   (chi_5 tail)
from scipy.stats import chi2
def chi5_tail(x):                                   # P(|W|>=x sigma) for W ~ N(0, sigma^2 I5)
    return chi2.sf(x*x, 5)
rows = []
for cc in [0.2, 0.4, 0.8, 1.2]:
    for ratio in [4.0, 9.0, 16.0, 25.0, 36.0, 49.0]:   # ratio = d^2/(nu tau)
        sig2 = 2*math.exp(2*cc)                         # Sigma_max/(nu tau)
        x = math.sqrt(ratio/sig2)
        tail = chi5_tail(x)
        # the claimed form C exp(-ratio/C'), C' = 4 e^{2c}:
        claim = math.exp(-ratio/(4*math.exp(2*cc)))
        rows.append(dict(c=cc, d2_over_nutau=ratio, chi5_tail=tail,
                         claim_exp=claim, C_needed=tail/claim))
rec("PART2_tail_vs_claim", rows)
rec("PART2_C_needed_max", max(r_['C_needed'] for r_ in rows))
rec("PART2_statement",
    "P(|Xi|>=d) <= C(c) exp(-d^2/(4 e^{2c} nu tau)) with C(c) = max over the sweep below")

with open(__file__.replace('v2_exact_strain.py', 'v2_results.json'), 'w') as fh:
    json.dump(OUT, fh, indent=1, default=str)
print("\nDONE v2")
