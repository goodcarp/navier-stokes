"""r2 - THE DECISIVE TEST.

prove-lagrangian sec 4(2) / status row T says the flow map is  T_{lambda(s)}  -- lambda depends on the
MATERIAL SHELL s -- up to relative O(c/L), and its closed form (4.1) is

    lambda(sigma,theta) = (1 - (1-sigma) kappa theta/2)^{-2},   sigma = log(rho/rho0)/L in [0,1],

which at the terminal time kappa*theta = 2(1-sqrt(2/3)) runs from lambda=3/2 (sigma=0) to 1 (sigma=1).

The gap-T seat's LEMMA T compares the flow map to a SINGLE T_lambda.  This script measures, with my
own quadrature (nothing imported from the seat):

  A. the sup relative displacement mu of Lambda against the BEST single T_lambdabar  (seat's (H2)),
     and the sup relative Jacobian deviation muJ against lambdabar^2 (seat's (H3)), as functions of L;
  B. the Lemma-T bound those measured constants produce, relative to a(0);
  C. the same two constants against the SHELL-LOCAL reference T_{lambda(rho)};
  D. the true a[eta0 o Lambda^-1](0) and its error against (i) the single-lambdabar reference
     (M/2) lambdabar L and (ii) the shell-resolved reference (M/2) int lambda dlog rho -- the one
     step P of prove-lagrangian actually consumes;
  E. the per-shell repair: LEMMA T applied on unit-log-width shells and summed.

Own quadrature.  With Lambda(x) = (lambda(rho) y, lambda(rho)^-2 z), eta transported,
   a = -(3/4) L int_0^1 dtau int_0^pi  cos(phi) sin^2(phi) omega(phi)
                  [1 + D(rho)(sin^2 phi - 2 cos^2 phi)] / g(phi;lambda)^5  dphi ,
   g = sqrt(lambda^2 sin^2 + lambda^-4 cos^2),  D = dlog(lambda)/dlog(rho),  rho = rho0 e^{L tau}.
(derived and symbolically checked in r1_algebra.py)
"""
import json, numpy as np

M, rho0 = 1.0, 1.0
KTH = 2*(1-np.sqrt(2.0/3.0))          # terminal kappa*theta of prove-lagrangian (4.1)

def lam_of_sigma(sig, kth=KTH):  return (1-(1-sig)*kth/2)**-2
def dlogl_dsig(sig, kth=KTH):    return -kth/(1-(1-sig)*kth/2)     # d log lambda / d sigma

def gl(a, b, n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5*(b-a)*x+0.5*(a+b), 0.5*(b-a)*w

def panels_gl(pts, n):
    xs, ws = [], []
    for a, b in zip(pts[:-1], pts[1:]):
        x, w = gl(a, b, n); xs.append(x); ws.append(w)
    return np.concatenate(xs), np.concatenate(ws)

def omega_bb(phi):   return -M*np.sign(np.cos(phi))
def omega_tap(phi, d=np.deg2rad(7.5)):
    pa = np.minimum(phi, np.pi-phi)
    return -M*np.sign(np.cos(phi))*np.minimum(1.0, pa/d)

def gfun(phi, lam):  return np.sqrt(lam**2*np.sin(phi)**2 + lam**-4*np.cos(phi)**2)

def a_shell_dep(L, omega, n_phi=200, n_tau=200, kth=KTH, use_D=True, taper=False):
    """a[eta0 o Lambda^-1](0) for the shell-dependent strain Lambda."""
    pts = [0.0, np.pi/2, np.pi]
    if taper:
        d = np.deg2rad(7.5); pts += [d, np.pi-d]
    ph, wph = panels_gl(np.array(sorted(set(pts))), n_phi)
    ta, wta = gl(0.0, 1.0, n_tau)
    PH, TA = np.meshgrid(ph, ta, indexing='ij'); WP, WT = np.meshgrid(wph, wta, indexing='ij')
    lam = lam_of_sigma(TA, kth)
    D = dlogl_dsig(TA, kth)/L if use_D else 0.0    # dlog lambda / dlog rho
    s2, c2 = np.sin(PH)**2, np.cos(PH)**2
    integ = np.cos(PH)*s2*omega(PH)*(1 + D*(s2-2*c2))/gfun(PH, lam)**5
    return -(3.0/4.0)*L*np.sum(integ*WP*WT)

def a_single(L, lam, omega, n_phi=400, taper=False):
    pts = [0.0, np.pi/2, np.pi]
    if taper:
        d = np.deg2rad(7.5); pts += [d, np.pi-d]
    ph, wph = panels_gl(np.array(sorted(set(pts))), n_phi)
    s2 = np.sin(ph)**2
    return -(3.0/4.0)*L*np.sum(np.cos(ph)*s2*omega(ph)*wph/gfun(ph, lam)**5)

def hyp_constants(L, lambar, kth=KTH, n=4001, m=2001):
    """sup over S of  |Lambda x - T_lambar x|/|T_lambar x|   and   |J_Lambda/lambar^2 - 1|."""
    sig = np.linspace(0, 1, n); ph = np.linspace(0, np.pi, m)
    S, P = np.meshgrid(sig, ph, indexing='ij')
    lam = lam_of_sigma(S, kth); D = dlogl_dsig(S, kth)/L
    s, c = np.sin(P), np.cos(P)
    # displacement of Lambda x vs T_lambar x, both applied to the same x = rho(s,c)
    dr = (lam-lambar)*s; dz = (lam**-2-lambar**-2)*c
    num = np.hypot(dr, dz); den = np.hypot(lambar*s, lambar**-2*c)
    mu = float(np.max(num/den))
    J = lam**2*(1 + D*(s**2-2*c**2))
    muJ = float(np.max(np.abs(J/lambar**2 - 1)))
    return mu, muJ

def hyp_constants_local(L, kth=KTH, n=4001, m=2001):
    """same, but against the SHELL-LOCAL reference T_{lambda(rho)}: mu = 0 identically."""
    sig = np.linspace(0, 1, n); ph = np.linspace(0, np.pi, m)
    S, P = np.meshgrid(sig, ph, indexing='ij')
    lam = lam_of_sigma(S, kth); D = dlogl_dsig(S, kth)/L
    s, c = np.sin(P), np.cos(P)
    J = lam**2*(1 + D*(s**2-2*c**2))
    return 0.0, float(np.max(np.abs(J/lam**2 - 1)))

def lemmaT_bound(lam, L, mu, muJ):
    if mu >= 1: return float('inf')
    return np.pi*lam*M*L*(3*mu*(1+muJ)/(2*(1-mu)**5) + 3*muJ/8.0)

# ---------------------------------------------------------------- A/B/C/D
lam_mean = None
rows = []
for L in (8.317766166719343, 20.0, 50.0, 100.0, 400.0):
    # shell-resolved reference  (M/2) int_0^L lambda dlog rho = (M/2) L int_0^1 lambda dsigma
    sg, wg = gl(0, 1, 400)
    lam_int = float(np.sum(lam_of_sigma(sg)*wg))
    lam_mean = lam_int
    a_ref_shell = M/2*L*lam_int
    a_true  = a_shell_dep(L, omega_bb)
    a_noD   = a_shell_dep(L, omega_bb, use_D=False)
    # best single lambdabar: minimise sup relative displacement
    grid = np.linspace(1.0, 1.5, 501)
    mus = [hyp_constants(L, lb, n=801, m=401)[0] for lb in grid]
    lb_star = float(grid[int(np.argmin(mus))]); mu_star = float(np.min(mus))
    mu_s, muJ_s = hyp_constants(L, lb_star)
    a_ref_single = M/2*lb_star*L
    mu_loc, muJ_loc = hyp_constants_local(L)
    rows.append(dict(
        L=L, lam_min=float(lam_of_sigma(0.0)), lam_max=float(lam_of_sigma(1.0)),
        lam_logmean=lam_int,
        lambar_star=lb_star, mu_single=mu_s, muJ_single=muJ_s,
        mu_local=mu_loc, muJ_local=muJ_loc,
        a_true=a_true, a_ref_shell=a_ref_shell, a_ref_single=a_ref_single,
        rel_err_shell_ref=abs(a_true-a_ref_shell)/a_true,
        rel_err_noD_vs_shellref=abs(a_noD-a_ref_shell)/a_true,
        bound_single_rel=lemmaT_bound(lb_star, L, mu_s, muJ_s)/a_true,
        bound_local_rel=lemmaT_bound(lam_int, L, mu_loc, muJ_loc)/a_true,
        rel_err_lam1=abs(a_true-M/2*1.0*L)/a_true,
        rel_err_lam15=abs(a_true-M/2*1.5*L)/a_true,
    ))

# ---------------------------------------------------------------- E per-shell repair
def per_shell_repair(L, n_shell):
    """split S into n_shell equal log-width pieces; on each, compare Lambda to T_{lambda(mid)}.
    Sum the LEMMA T bounds.  Returns (sum of bounds)/a_true and the measured per-shell (mu,muJ)."""
    edges = np.linspace(0, 1, n_shell+1)
    tot, mus, muJs = 0.0, [], []
    ph = np.linspace(0, np.pi, 2001)
    s, c = np.sin(ph), np.cos(ph)
    for i in range(n_shell):
        s0, s1 = edges[i], edges[i+1]
        smid = 0.5*(s0+s1); lb = lam_of_sigma(smid)
        sg = np.linspace(s0, s1, 401)
        SS, PP = np.meshgrid(sg, ph, indexing='ij')
        lam = lam_of_sigma(SS); D = dlogl_dsig(SS)/L
        ss, cc = np.sin(PP), np.cos(PP)
        num = np.hypot((lam-lb)*ss, (lam**-2-lb**-2)*cc); den = np.hypot(lb*ss, lb**-2*cc)
        mu = float(np.max(num/den))
        J = lam**2*(1+D*(ss**2-2*cc**2)); muJ = float(np.max(np.abs(J/lb**2-1)))
        mus.append(mu); muJs.append(muJ)
        tot += lemmaT_bound(lb, L*(s1-s0), mu, muJ)
    return tot, max(mus), max(muJs)

repair = []
for L in (8.317766166719343, 20.0, 50.0, 100.0, 400.0):
    a_true = a_shell_dep(L, omega_bb)
    for nsh in (max(1,int(round(L))),):     # unit-log-width shells
        tot, mu, muJ = per_shell_repair(L, nsh)
        repair.append(dict(L=L, n_shells=nsh, mu_max=mu, muJ_max=muJ,
                           sum_bound=tot, a_true=a_true, rel=tot/a_true))

# ---------------------------------------------------------------- convergence controls
conv = []
for L in (8.317766166719343, 100.0):
    v = [a_shell_dep(L, omega_bb, n_phi=n, n_tau=n) for n in (100, 200, 400)]
    conv.append(dict(L=L, vals=v, rel_200_400=abs(v[2]-v[1])/abs(v[2])))
# instrument cross-check against the estate: a_single at lambda=1 must be (M/2)L; taper -> kappa_delta
xchk = dict(
    a_single_lam1_over_ML=[float(a_single(L, 1.0, omega_bb)/(M*L)) for L in (8.317766166719343, 50.0)],
    kappa_delta=float(a_single(8.317766166719343, 1.0, omega_tap, taper=True)/(M*8.317766166719343)),
    a_single_lam_15_over_ML=float(a_single(8.317766166719343, 1.5, omega_bb)/(M*8.317766166719343)),
)

out = dict(kappa_theta_terminal=KTH, lam_logmean=lam_mean,
           sqrt_3_over_2=float(np.sqrt(1.5)), rows=rows, repair=repair,
           convergence=conv, cross_check=xchk)
json.dump(out, open('r2_results.json','w'), indent=1)

print(f"kappa*theta terminal = {KTH:.10f}   lambda(sigma=0)=3/2, lambda(sigma=1)=1")
print(f"log-mean of lambda over the shell = {lam_mean:.12f}   sqrt(3/2) = {np.sqrt(1.5):.12f}")
print("\ncross-check of my instrument:", xchk)
print("\nA/B/C/D  (bang-bang cap, terminal time):")
h = (f"{'L':>8s} {'lambar*':>8s} {'mu(H2)':>8s} {'muJ(H3)':>9s} {'muLOC':>6s} {'muJLOC':>8s} "
     f"{'relerr shellref':>16s} {'LemT rel bnd single':>20s} {'LemT rel bnd local':>19s}")
print(h)
for r in rows:
    print(f"{r['L']:8.3f} {r['lambar_star']:8.4f} {r['mu_single']:8.4f} {r['muJ_single']:9.4f} "
          f"{r['mu_local']:6.2f} {r['muJ_local']:8.5f} {r['rel_err_shell_ref']:16.3e} "
          f"{r['bound_single_rel']:20.4g} {r['bound_local_rel']:19.4g}")
print("\nrelative error if the single lambda is taken at the ENDS of the range prove-lagrangian names:")
for r in rows:
    print(f"  L={r['L']:8.3f}  lambda=1: {r['rel_err_lam1']:.4f}   lambda=3/2: {r['rel_err_lam15']:.4f}")
print("\nE  per-shell repair (unit log-width shells):")
for r in repair:
    print(f"  L={r['L']:8.3f} n={r['n_shells']:4d} mu_max={r['mu_max']:.5f} muJ_max={r['muJ_max']:.5f} "
          f"sum(bound)/a = {r['rel']:.5f}")
print("\nquadrature convergence:", [(c['L'], f"{c['rel_200_400']:.2e}") for c in conv])
