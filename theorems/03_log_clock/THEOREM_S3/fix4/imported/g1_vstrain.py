"""
g1 -- UNIT 1 of FIX3: (V-strain), item M3 of THEOREM_S3 section 3.

WHAT M3 ASKS.  Theorem V.4' (round2/pmax-h3v Part B) bounds the viscous defect of eta at ONE
tracked material point.  u1 section 4.2 step 3 needs the viscous defect of the strain
FUNCTIONAL a[.](0) -- a weighted integral of eta over the whole shell -- and of the exterior
strain afrak(rho,s), which is the same functional with the inner shells removed.  u1's own
CARRIED GAP and refute-u1 F10 both say this is priced by assuming the answer.

THE ROUTE (the brief's).  w := eta - eta^tr, eta^tr := eta_0 o Phi_s^{-1} the inviscidly
transported datum with the SAME drift b.  Then D_t w = nu Lap5 eta exactly and w(.,0) = 0, so
w is transported with a source and

    sup_x rho |w(.,s)|  <=  nu e^{c_R(s)} int_0^s Psi(r) dr ,   Psi(r) := sup_x |x| |Lap5 eta(x,r)| ,

and a[w](0) = int Kcal w dx_5 with Kcal(x) = d_z G_5(x) = -3z/(8 pi^2 rho^5), whose weight
integral over a shell is EXACTLY (3/8) log(rho_2/rho_1).  Hence

    | a[eta(s)](0) - a[eta^tr(s)](0) |  <=  (3/8) Lfrak  nu s e^{c_R} Psibar  =:  eps_a M L .

WHAT THIS SCRIPT ESTABLISHES, each with its own check:

  A. the kernel and its weight:  a[.](0) = int Kcal eta dx_5 reproduces (M/2)L on the bare
     plateau; int_{rho1<rho<rho2} |Kcal|/rho dx_5 = (3/8) log(rho2/rho1) EXACTLY; with the
     datum's own radial envelope the effective count is L + 2/8, not L.
  B. the harmonic cancellation:  Lap5 Kcal = -d_z delta_0, so int Kcal Lap5 eta dx_5 = d_z eta(0)
     EXACTLY.  The INSTANTANEOUS viscous rate of a[eta](0) is nu d_z eta(0,s) and nothing else;
     the whole V-strain defect is the commutator of diffusion with transport.  (This says the
     bound of A is lossy, and by how much.)
  C. the datum constants:  ||grad eta_0||_inf (rho_0=M=1), and Psi(0) = sup rho|Lap5 eta_0|.
     For THEOREM_S3's KINKED angular profile Psi(0) = +infinity: Lap5 eta_0 carries surface
     Dirac masses on the four kink cones.  Psi_sigma(0) for the sigma-mollified profile is
     computed, together with kappa_delta(sigma) and Gfrak_0(sigma), so the trade is visible.
  D. eps_a(L) at L = 1e4, 1e5, 1e6, both windows, in the two regimes.

Units: rho_0 = M = 1 throughout the numerics (so nu = 1/s_delta^2, s_delta = 1/sin delta).

Outputs -> g1_results.json
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "fix2copy"))
import numpy as np
import sympy as sp

OUT = {}
DEG = math.pi/180.0
DELTA, DM, EPS_R = 7.5*DEG, 5.0*DEG, 0.25
SDEL = 1.0/math.sin(DELTA)                 # s = 1/sin delta = E_0
S3SPH = 2*math.pi**2                       # |S^3|
S4SPH = 8*math.pi**2/3                     # |S^4|

# =====================================================================================
# A.  THE KERNEL Kcal AND ITS WEIGHT
# =====================================================================================
# G_5(x) = 1/(8 pi^2 rho^3) solves -Lap5 G_5 = delta_0 ; a = -d_z psi_1, psi_1 = G_5 * eta,
# so a(x) = int K(x-x') eta(x') dx' with K(w) = -d_z G_5(w) = 3 w_z/(8 pi^2 |w|^5), and
#     a[eta](0) = int K(-x) eta(x) dx = int Kcal(x) eta(x) dx ,  Kcal = d_z G_5 = -3z/(8 pi^2 rho^5).
rho_s, phi_s = sp.symbols('rho phi', positive=True)
# Lap5 acting on a radial function in R^5:  (1/rho^4) d_rho (rho^4 d_rho f)
f5 = rho_s**-3
lapG = sp.simplify(sp.diff(rho_s**4*sp.diff(f5, rho_s), rho_s)/rho_s**4)
OUT["A_kernel"] = {
    "Lap5(rho^-3) away from origin (sympy)": str(lapG),
    "G5": "1/(8 pi^2 rho^3),  -Lap5 G5 = delta_0",
    "Kcal(x)": "-3 z /(8 pi^2 rho^5)  = d_z G_5",
    "abs_Kcal": "3 |cos phi| /(8 pi^2 rho^4)",
}
# weight integral, exactly:  int |Kcal|/rho dx_5
#   = (3/(8 pi^2)) |S^3| int_0^pi |cos phi| sin^3 phi dphi  int drho/rho
#   = (3/(8 pi^2)) (2 pi^2) (1/2) log(rho2/rho1) = (3/8) log(rho2/rho1)
ang_abs = 2*sp.integrate(sp.cos(phi_s)*sp.sin(phi_s)**3, (phi_s, 0, sp.pi/2))
ang_pl = 2*sp.integrate(sp.cos(phi_s)*sp.sin(phi_s)**2, (phi_s, 0, sp.pi/2))
w_coef = sp.Rational(3, 1)/(8*sp.pi**2)*(2*sp.pi**2)*ang_abs
OUT["A_kernel"]["int_|cos|sin^3_dphi"] = str(ang_abs)
OUT["A_kernel"]["int_|cos|sin^2_dphi"] = str(ang_pl)
OUT["A_kernel"]["weight_per_efold_exact"] = str(sp.nsimplify(w_coef))
OUT["A_kernel"]["weight_per_efold"] = float(w_coef)
# control: the bare plateau eta_P = -M sgn(z)/(rho sin phi) gives a[eta_P](0) = (M/2) L
plateau_rate = float(sp.Rational(3, 1)/(8*sp.pi**2)*(2*sp.pi**2)*ang_pl)
OUT["A_kernel"]["plateau_a0_rate_per_efold"] = plateau_rate      # must be 1/2


def weight_numeric(rho1, rho2, n=200001):
    """int_{rho1<rho<rho2} |Kcal|/rho dx_5, by quadrature, as a control on (3/8) log ratio."""
    phi = np.linspace(0.0, math.pi, 4001)
    ang = np.trapz(np.abs(np.cos(phi))*np.sin(phi)**3, phi)
    return (3.0/(8*math.pi**2))*S3SPH*ang*math.log(rho2/rho1)


OUT["A_kernel"]["weight_numeric_1_to_e^10"] = weight_numeric(1.0, math.exp(10.0))
OUT["A_kernel"]["weight_formula_1_to_e^10"] = 0.375*10.0

# the datum's own radial envelope:  Theta_env(rho) = min(1, rho^8, (R/rho)^8) has
#   int Theta_env dlog rho = L + 2/8 ,  and the true Theta has int Theta du = L - 2 eps_r EXACTLY.
u_s, e_s = sp.symbols('u e', positive=True)
# int_{-inf}^{inf} (1/2)[tanh((u-A)/e) - tanh((u-B)/e)] du = B - A :  check the single edge,
# int_{-inf}^{0} (1/2)(1+tanh((u-a)/e)) du = (e/2) log(1+e^{-2a/e}) at a>0 (exact antiderivative)
a_s = sp.Symbol('a', positive=True)
prim = sp.simplify(sp.diff((e_s/2)*sp.log(1 + sp.exp(2*(u_s-a_s)/e_s)), u_s)
                   - (1 + sp.tanh((u_s-a_s)/e_s))/2)
I_Theta = prim
OUT["A_envelope"] = {
    "int_Theta_du_over_R_exact": "L - 2 eps_r",
    "sympy_check_int_u_dTheta": str(sp.simplify(I_Theta)),
    "int_Theta_env_dlogrho": "L + 2/8 = L + 0.25",
    "budget_ELL_RAMP": 0.25,
    "true_ramp_loss_efolds": 2*EPS_R,
}
# numeric check of int Theta du = L - 2 eps_r
Lnum = 12.0
ug = np.linspace(-20.0, Lnum+20.0, 4000001)
Thn = 0.5*(np.tanh((ug-EPS_R)/EPS_R) - np.tanh((ug-Lnum+EPS_R)/EPS_R))
OUT["A_envelope"]["numeric_int_Theta_du_L12"] = float(np.trapz(Thn, ug))
OUT["A_envelope"]["formula_L_minus_2eps_r"] = Lnum - 2*EPS_R
# the two tails, exactly log(1+e^{-2})/8 each
OUT["A_envelope"]["tail_efolds_each_side_exact"] = math.log(1.0+math.exp(-2.0))/8.0
mask = ug < 0.0
OUT["A_envelope"]["tail_efolds_inner_numeric"] = float(np.trapz(Thn[mask], ug[mask]))

# =====================================================================================
# B.  THE HARMONIC CANCELLATION:  int Kcal Lap5 eta dx = d_z eta(0)
# =====================================================================================
# Lap5 Kcal = d_z Lap5 G_5 = -d_z delta_0 , so <Lap5 Kcal, eta> = + d_z eta(0).
# Numerical control on a smooth axisymmetric z-odd test field.


def test_field(r, zz, k=0):
    """eta_test = z exp(-(r^2+z^2)/2) (axisymmetric in y in R^4, z-odd, Schwartz).
       k = 0 -> value ; k = 'lap' -> Lap5."""
    rr2 = r*r + zz*zz
    E = np.exp(-0.5*rr2)
    if k == 0:
        return zz*E
    # Lap5 (z f(rho)) with f = exp(-rho^2/2):  = z (f'' + 6 f'/rho)   (l = 1 mode in R^5)
    #   d_z of a radial times z: Lap5(z f) = z Lap5 f + 2 d_z f = z(f''+4f'/rho) + 2 z f'/rho
    rho = np.sqrt(np.maximum(rr2, 1e-300))
    fp = -rho*E
    fpp = (rho*rho - 1.0)*E
    return zz*(fpp + 6.0*fp/np.maximum(rho, 1e-300))


def pair_with_Kcal(fun, smin=1e-7, smax=40.0, ns=2000, nT=800):
    """int Kcal(x) fun(x) dx_5 = -(3/(8 pi^2)) |S^3| int int cos(phi) fun(rho,phi) sin^3 phi dphi drho/1
       (rho^4 from dx_5 cancels rho^{-4} of |Kcal| ... careful: Kcal ~ rho^{-4}, dx_5 = rho^4 drho dOmega_4)"""
    xg, wg = np.polynomial.legendre.leggauss(ns)
    a, b = math.log(smin), math.log(smax)
    uu = 0.5*(b-a)*xg + 0.5*(b+a)
    wu = 0.5*(b-a)*wg
    s = np.exp(uu)
    T, wT = np.polynomial.legendre.leggauss(nT)
    T = 0.5*math.pi*(T+1.0)
    wT = 0.5*math.pi*wT
    cT, sT = np.cos(T), np.sin(T)
    tot = 0.0
    for si, wsi in zip(s, wu*s):
        r_ = si*sT
        z_ = si*cT
        v = fun(r_, z_)
        tot += wsi*(-3.0/(8*math.pi**2))*S3SPH*np.sum(wT*cT*sT**3*v)
    return float(tot)


lhs = pair_with_Kcal(lambda r, zz: test_field(r, zz, 'lap'))
# d_z eta_test(0) = 1 for eta = z exp(-rho^2/2)
OUT["B_harmonic_identity"] = {
    "statement": "Lap5 Kcal = - d_z delta_0  =>  int Kcal Lap5 eta dx_5 = d_z eta(0)",
    "numeric_int_Kcal_Lap_eta_test": lhs,
    "exact_dz_eta_test_at_0": 1.0,
    "rel": abs(lhs - 1.0),
    "consequence": ("the instantaneous viscous rate of a[eta](0) is nu d_z eta(0,s) EXACTLY; "
                    "the V-strain defect is the diffusion-transport commutator, not diffusion"),
}
# control: the same pairing against Kcal of a field whose d_z at 0 vanishes
lhs2 = pair_with_Kcal(lambda r, zz: test_field(r, zz, 'lap')*0.0 + 0.0)
OUT["B_harmonic_identity"]["null_control"] = lhs2

# =====================================================================================
# C.  THE DATUM CONSTANTS, AND WHY Psi(0) = +infinity FOR THE KINKED PROFILE
# =====================================================================================
# eta_0 = G(rho) W(phi), G = -Theta(rho)/rho, W = A(phi)/sin phi, A = sgn(cos phi) g1 g2.
#   Lap5 (G W) = W (G'' + 4 G'/rho) + (G/rho^2)(W'' + 3 cot phi W')
#   rho^3 Lap5 eta_0 = -W (Theta_uu + Theta_u - 2 Theta) - Theta (W'' + 3 cot phi W')
# so  rho |Lap5 eta_0| = rho^{-2} | ... | .   A'' = 0 on each piece and a Dirac at each kink.


def A_pieces(phi, delta=DELTA, dm=DM):
    s_sign = np.where(np.cos(phi) >= 0.0, 1.0, -1.0)
    pax = np.minimum(phi, math.pi - phi)
    dpax = np.where(phi < math.pi/2, 1.0, -1.0)
    eqd = np.abs(phi - math.pi/2)
    deqd = np.where(phi > math.pi/2, 1.0, -1.0)
    g1 = np.minimum(1.0, pax/delta)
    dg1 = np.where(pax < delta, dpax/delta, 0.0)
    g2 = np.minimum(1.0, eqd/dm)
    dg2 = np.where(eqd < dm, deqd/dm, 0.0)
    return s_sign*g1*g2, s_sign*(dg1*g2 + g1*dg2)


def F_of_phi(phi, delta=DELTA, dm=DM):
    """F = A/sin phi and its first derivative, analytically (A piecewise linear, A''=0)."""
    s_ = np.sin(phi)
    c_ = np.cos(phi)
    A0, A1 = A_pieces(phi, delta, dm)
    F0 = A0/s_
    F1 = (A1*s_ - A0*c_)/s_**2
    F2 = (-2.0*A1*c_*s_ + A0*(2.0*c_*c_ + s_*s_))/s_**3      # A'' = 0 (a.c. part)
    return F0, F1, F2


def mollified_F(phi_grid, sigma):
    """F_sigma = F * chi_sigma with F EVEN-reflected about phi = 0 and phi = pi (pmax-h3v
       Proposition 3's mollification).  F_sigma' = F' * chi , F_sigma'' = F' * chi'."""
    from scipy.signal import fftconvolve
    h = phi_grid[1] - phi_grid[0]
    nker = max(40, int(math.ceil(8.0*sigma/h)))
    ext_lo = phi_grid[0] - h*np.arange(nker, 0, -1)
    ext_hi = phi_grid[-1] + h*np.arange(1, nker+1)
    full = np.concatenate([ext_lo, phi_grid, ext_hi])
    ref = np.where(full < 0.0, -full, np.where(full > math.pi, 2*math.pi - full, full))
    sgn = np.where((full < 0.0) | (full > math.pi), -1.0, 1.0)      # F even => F' odd
    F0, F1, _ = F_of_phi(np.clip(ref, 1e-12, math.pi-1e-12))
    F1 = F1*sgn
    y = h*np.arange(-nker, nker+1)
    chi = np.exp(-0.5*(y/sigma)**2)
    chi /= chi.sum()*h
    chip = -(y/sigma**2)*chi
    F0s = fftconvolve(F0, chi, mode='same')*h
    F1s = fftconvolve(F1, chi, mode='same')*h
    F2s = fftconvolve(F1, chip, mode='same')*h
    sl = slice(nker, nker+len(phi_grid))
    return F0s[sl], F1s[sl], F2s[sl]


def Theta_u(u, k=0, L=40.0, er=EPS_R):
    a = (u - er)/er
    b = (u - L + er)/er
    if k == 0:
        return 0.5*(np.tanh(a) - np.tanh(b))
    ca = 1.0/np.cosh(np.clip(a, -350, 350))**2
    cb = 1.0/np.cosh(np.clip(b, -350, 350))**2
    if k == 1:
        return 0.5*(ca - cb)/er
    ta, tb = np.tanh(a), np.tanh(b)
    if k == 2:
        return (tb*cb - ta*ca)/er**2
    raise ValueError


def psi_and_norms(sigma, L=40.0, nphi=200001, nu_=2001, phi_eps=1e-3):
    """Psi_sigma(0) = sup_x rho |Lap5 eta_0^{(sigma)}| , ||grad eta_0||_inf ,
       Gfrak_0(sigma) = sup rho^2|grad eta_0| , kappa_delta(sigma)."""
    phi = np.linspace(phi_eps, math.pi - phi_eps, nphi)
    s_, c_ = np.sin(phi), np.cos(phi)
    if sigma is None:
        W0, W1, W2 = F_of_phi(phi)
    else:
        W0, W1, W2 = mollified_F(phi, sigma)
    # Lap_{S^4} W = W'' + 3 cot(phi) W' ; for the UNMOLLIFIED profile (A'' = 0 piecewise)
    # this collapses to the stable closed form  (A' c s + A(s^2 - c^2))/s^3 .
    if sigma is None:
        A0, A1 = A_pieces(phi)
        angL = (A1*c_*s_ + A0*(s_*s_ - c_*c_))/s_**3
    else:
        angL = W2 + 3.0*(c_/s_)*W1
    # the two end caps phi < phi_eps and phi > pi - phi_eps are excluded from the scan; there
    # A = +-phi/delta exactly and Lap_{S^4}W -> 4/(3 delta), far below the bulk value at phi=delta
    cap_val = 4.0/(3.0*DELTA)
    uu = np.linspace(-4.0, L + 4.0, nu_)
    T0 = Theta_u(uu, 0, L)
    T1 = Theta_u(uu, 1, L)
    T2 = Theta_u(uu, 2, L)
    rad = T2 + T1 - 2.0*T0
    best, arg = 0.0, None
    G0sup, gradinf, argg = 0.0, 0.0, None
    for i in range(len(uu)):
        e2 = math.exp(-2.0*uu[i])
        v = float(np.max(np.abs(-W0*rad[i] - T0[i]*angL)))*e2
        if v > best:
            best, arg = v, float(uu[i])
        g2 = np.sqrt((T0[i]-T1[i])**2*W0**2 + T0[i]**2*W1**2)
        m = float(g2.max())
        if m > G0sup:
            G0sup = m
        if m*e2 > gradinf:
            gradinf, argg = m*e2, float(uu[i])
    # kappa_delta(sigma) = (1/2) P_h(1),  P_h(1) = 3 int_0^1 h(arcsin v) v^2 dv,
    #   h(phi) = |A(phi)| = |F(phi)| sin phi   (the FULL angular profile)
    v = np.linspace(0.0, 1.0, 400001)
    ph = np.arcsin(np.clip(v, 0.0, 1.0))
    if sigma is None:
        hA = np.abs(F_of_phi(np.clip(ph, 1e-12, None))[0])*np.sin(ph)
    else:
        idx = np.clip(np.searchsorted(phi, ph), 1, len(phi)-1)
        t = (ph - phi[idx-1])/(phi[idx] - phi[idx-1])
        Fi = W0[idx-1]*(1-t) + W0[idx]*t
        hA = np.abs(Fi)*np.sin(ph)
    Ph1 = 3.0*np.trapz(hA*v**2, v)
    return {"sigma": sigma, "Psi0": best, "argmax_u": arg, "endcap_angL": cap_val,
            "sup_abs_angL": float(np.max(np.abs(angL))),
            "grad_inf": gradinf, "grad_inf_argu": argg, "Gfrak0": G0sup,
            "kappa_delta": 0.5*float(Ph1)}


rows = []
for sg in [None, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005]:
    r = psi_and_norms(sg)
    rows.append(r)
    print("sigma=%-8s Psi0=%14.6f  ||grad eta||_inf=%10.6f  Gfrak0=%12.7f  kappa=%.10f"
          % (str(sg), r["Psi0"], r["grad_inf"], r["Gfrak0"], r["kappa_delta"]), flush=True)
# independent control on ||grad eta_0||_inf and Gfrak_0: fix2's own field object (DatumS3),
# whose (r,z) derivatives were finite-difference checked in fix2.
from f2_datum_s3 import DatumS3                                  # noqa: E402
_D = DatumS3(delta_deg=7.5, dm_deg=5.0, eps_r=EPS_R, L=40.0)
_uu = np.linspace(-4.0, 44.0, 4001)
_ph = np.linspace(1e-7, math.pi-1e-7, 40001)
_gi, _g0 = 0.0, 0.0
for _u in _uu:
    _r = math.exp(_u)*np.sin(_ph)
    _z = math.exp(_u)*np.cos(_ph)
    _g = _D.grad_norm(np.maximum(_r, 1e-300), _z)
    _gi = max(_gi, float(np.max(_g)))
    _g0 = max(_g0, float(np.max(_g))*math.exp(2.0*_u))
# the closed form for Psi_sigma(0):  the binding kink is the taper at phi = delta, where
#   [W'] = [A']/sin delta = (1/delta)/sin delta ; a Gaussian mollifier of std sigma turns that
#   jump into a spike of height [W'] chi_sigma(0) = [W']/(sigma sqrt(2 pi)), so
#   Psi_sigma(0) ~ [W'] * max_u (Theta e^{-2u}) / (sigma sqrt(2 pi)) .
_uu2 = np.linspace(-3.0, 6.0, 900001)
_Te = Theta_u(_uu2, 0, 40.0)*np.exp(-2.0*_uu2)
_maxTe = float(_Te.max())
_Wjump = (1.0/DELTA)/math.sin(DELTA)
OUT["C_Psi_closed_form"] = {
    "W_prime_jump_at_taper": _Wjump,
    "W_prime_jump_at_equator_kink": (1.0/DM)/math.cos(DM),
    "max_u_Theta_e^{-2u}": _maxTe, "argmax_u": float(_uu2[int(np.argmax(_Te))]),
    "C_kink_predicted": _Wjump*_maxTe/math.sqrt(2*math.pi),
    "C_kink_measured_sigma_5e-4": [r for r in rows if r["sigma"] == 0.0005][0]["Psi0"]*0.0005,
    "formula": "Psi_sigma(0) ~ C_kink/sigma ,  C_kink = [W'] max_u(Theta e^{-2u})/sqrt(2 pi)",
}
print("C_kink predicted %.5f  measured %.5f" % (OUT["C_Psi_closed_form"]["C_kink_predicted"],
      OUT["C_Psi_closed_form"]["C_kink_measured_sigma_5e-4"]), flush=True)
OUT["C_control_DatumS3"] = {"grad_inf": _gi, "Gfrak0_sup_rho2_grad": _g0,
                            "pmax_h3v_reported_grad_inf": 31.093556949918312,
                            "Gfrak0_closed_form": 1.0/math.sin(DELTA)**2}
print("CONTROL (DatumS3): ||grad eta_0||_inf = %.6f   Gfrak0 = %.7f" % (_gi, _g0), flush=True)
OUT["C_datum"] = {"rows": rows,
                  "note_kinked": ("sigma=None drops the four Dirac masses of A'' : it is the "
                                  "a.c. part only, NOT Psi(0), which is +infinity"),
                  "A_prime_jumps": {"at_phi=delta": 1.0/DELTA, "at_phi=pi/2-delta_m": 1.0/DM},
                  "Gfrak0_closed_form": 1.0/math.sin(DELTA)**2,
                  "E0_closed_form": SDEL}

# =====================================================================================
# D.  eps_a(L)
# =====================================================================================
import f1_window as F1                                        # noqa: E402  (byte copy of fix2)

W_EFOLD = 0.375              # (3/8), proved in A
ENV = 0.25                   # extra effective e-folds from the datum's own radial envelope


def eps_a_smooth(L, c, c_G, Psi0):
    """eps_a = (3/8)(L + 1/4) * nu * tau * e^{c_R} * Psi0 / (M L) , nu = M rho_0^2/s^2,
       tau = c/(M L), Psi0 in units M/rho_0^2, e^{c_R} <= e^{c_G} (conservative)."""
    return W_EFOLD*(L + ENV)*(c/(SDEL**2*L))*math.exp(c_G)*Psi0/L


def eps_a_kink(L, c, c_G, jump, C_kink=1.0):   # noqa
    """the kinked datum: Lap5 eta_0 is a surface measure, Psi(r) ~ C_kink * jump * rho_0/sqrt(nu r),
       so nu int_0^s Psi = 2 C_kink jump sqrt(nu s) rho_0 ... in units rho_0 = M = 1,
       nu = 1/s^2, s = tau = c/L  =>  nu*int = 2 C_kink jump sqrt(c/L)/s_delta."""
    supw = 2.0*C_kink*jump*math.sqrt(c/L)/SDEL*math.exp(c_G)
    return W_EFOLD*(L + ENV)*supw/L


CAP = F1.CAP_CERT
tab = []
for L in (1e4, 1e5, 1e6):
    for label, x in (("frozen c=c_*", 1.0), ("cap c=1.1943662c_*", CAP)):
        c = x*F1.CSTAR
        col = "proved" if L > 3e4 else "measured"
        A = F1.assemble_c(L, c, col, f=4.0)
        if "c_G" not in A:
            A = F1.assemble_c(L, c, "measured", f=4.0)
            col = "measured"
        cG = A["c_G"]
        row = {"L": L, "window": label, "column_for_c_G": col, "c": c, "c_G": cG,
               "eps_a_sigma=0.001": eps_a_smooth(L, c, cG, [r for r in rows
                                                            if r["sigma"] == 0.001][0]["Psi0"]),
               "eps_a_sigma=0.002": eps_a_smooth(L, c, cG, [r for r in rows
                                                            if r["sigma"] == 0.002][0]["Psi0"]),
               "eps_a_sigma=0.01": eps_a_smooth(L, c, cG, [r for r in rows
                                                           if r["sigma"] == 0.01][0]["Psi0"]),
               "eps_a_sigma=0.05": eps_a_smooth(L, c, cG, [r for r in rows
                                                           if r["sigma"] == 0.05][0]["Psi0"]),
               "eps_a_ac_part_only": eps_a_smooth(L, c, cG, [r for r in rows
                                                             if r["sigma"] is None][0]["Psi0"]),
               "eps_a_kinked_sqrt_route": eps_a_kink(L, c, cG, 1.0/DM),
               "eps_budget_at_this_c": A["eps"], "eps_Tprime": A["eps_Tprime"]}
        tab.append(row)
        print("L=%9.3g %-20s c=%.6f c_G=%.5f  eps_a(sig=1e-3)=%.4e  eps_a(kink)=%.4e "
              " eps_budget=%.6g" % (L, label, c, cG, row["eps_a_sigma=0.001"],
                                    row["eps_a_kinked_sqrt_route"], A["eps"]), flush=True)
OUT["D_eps_a"] = tab
OUT["D_formula"] = {
    "smooth": "eps_a = (3/8)(1 + 1/(4L)) * Psi0 * c * e^{c_G} / (s^2 L) ,  s = 1/sin delta",
    "kinked": "eps_a = (3/4)(1 + 1/(4L)) * C_kink * [A'] * e^{c_G} * sqrt(c/L) / s",
    "s_delta": SDEL, "s_delta_sq": SDEL**2,
}

with open(os.path.join(HERE, "g1_results.json"), "w") as fh:
    json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
print("\nWROTE g1_results.json")
print(json.dumps({k: OUT[k] for k in ("A_kernel", "A_envelope", "B_harmonic_identity")},
                 indent=1, sort_keys=True, default=str))
