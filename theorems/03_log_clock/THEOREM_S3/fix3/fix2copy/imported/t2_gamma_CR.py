"""
t2 -- (Gamma-off)'s constant C'' as a function of L (U2's fixed point, re-implemented here),
      and the velocity-remainder constant  C_R  with U2's proved r|grad a| bound in place of
      U1's unproved collar-gradient placeholder G_collar.

TWO REPLACEMENTS, both forced by the round-2 seats:

 (1) C'' replaces C'.  ASSEMBLY / U1 carry hypothesis (Gamma-off) with C' = 151.15, the
     exact-plateau constant of write/L3v-and-gamma-bound + its refuter, and C_1 = 1.
     U2 PROVES (Gamma-off) on the whole of R^5 for the true field, with C_1 = 0 and
        C'' = 2 Chat_a + (3 pi^2/16) e^{p c_G} Gfrak_0 + lambda ,
     a fixed point in (c_G, p) because c_G = (lambda + C''/L) c.  This file solves that fixed
     point at each L; it is NOT a constant.  The control below reproduces U2's own
     C''(10^4) = 582.24 (sigma_* = 1/2) and 801.05 (sigma_* = 0, the whole of R^5) from
     U2's datum norms, by an instrument written here.

 (2) G_collar is discharged.  U1 sec.3.2 bounds the z-component of the velocity remainder with
     G_1 <= 1/2 + G_far + G_inner + G_collar and carries (R-z): "G_collar an absolute constant,
     NOT PROVED".  U2 sec.6 proves  r|grad a| <= (3 pi^2/16) e^{p c_G} Gfrak_0 M sin phi
     on all of R^5, so  G_1 <= Ghat sin phi  with no decomposition and no unproved piece.
     (R-z) is therefore discharged, at the price of Ghat ~ 10^3 instead of ~2.2 + G_collar.

 (3) The angular supremum.  refute-u1 F3 keeps the geometry U1 drops (|R_y| <= C_a(phi) M rho
     sin phi, so the collar's 1/sin phi cancels; and the exact zeta-integral
     int_0^|z| log(rho/rho_zeta) dzeta = |z| - r arctan(|z|/r), so there is no log(1/sin phi)).
     U2 sec.7(b) supplies the axis-safe collar |a_collar| <= 2 R_A e^{c_G} E_0 M with NO
     angular dependence, PROVED, needing no taper-transport statement.  Together the
     supremum of C_R over ALL phi is finite, so hypothesis (A-cone) is DISCHARGED.

Everything is expressed against M (the datum scale), not M_s: the collar branch and the
r|grad a| bound are absolute in M, while the far/near and (P2) pieces carry lam_om = M_s/M.

Outputs -> t2_results.json.  Imported as a module by t3_budget.py.
"""
import json, math
import numpy as np

# ------------------------------------------------------------------ exact kernel constants
C_K = 3.0/(8.0*math.pi**2)                 # |K(w)| <= C_K |w|^-4 ,  K = -d_z G_5
S4 = 8.0*math.pi**2/3.0                    # |S^4| = 2 pi^{5/2}/Gamma(5/2) = 8 pi^2/3
CK_S4 = C_K*S4                             # = 1 exactly
RIESZ = math.pi**4/2.0                     # int |x-x'|^-4 |x'|^-2 dx' = (pi^4/2)/|x| in R^5
COEF_G = C_K*RIESZ                         # = 3 pi^2/16
R_A = (2.0**5 - 2.0**-5)**0.2              # annulus rho/2<|x'|<2rho has the volume of B(0,R_A rho)
TWO_RA = 2.0*R_A

# far/near kernel lemma, z-odd constants.  FOREIGN: rebuild/far-near-kernel-lemma sec.2-3,
# re-derived independently by s3close/round2/u2 (u3) to 0.29199853 / 0.01475427 / 3.9992184.
C_FAR = 0.29199853
C_INNER = 0.01475427

def riesz_check(nxi=400001, nc=20001):
    """int_{R^5} |x-x'|^-4 |x'|^-2 dx' at |x| = 1.  In polar coordinates centred on x,
       w = x - x', dw = xi^4 dxi dOmega_4, dOmega_4 = |S^3| sin^3(alpha) dalpha = 2 pi^2 sin^3 a da,
       and |x'|^2 = 1 - 2 xi cos a + xi^2, so the xi^4 cancels the |w|^-4 exactly:
           I = 2 pi^2 int_0^inf dxi int_{-1}^{1} (1-c^2)/(1 - 2 xi c + xi^2) dc .
       No singularity survives.  Compare with the Gamma-function value pi^4/2."""
    c = np.linspace(-1.0, 1.0, nc)
    wc = 1.0 - c*c
    xi = np.concatenate([np.linspace(0.0, 40.0, nxi//2),
                         40.0 + np.tan(np.linspace(0.0, math.atan(1e7), nxi//2))])
    inner = np.array([np.trapz(wc/(1.0 - 2.0*x*c + x*x), c) for x in xi])
    return 2.0*math.pi**2*float(np.trapz(inner, xi))

# ------------------------------------------------------------------ (Gamma-off) fixed point
def chat_a(sigma_star, lam, c_G, E0):
    """U2 (7.3): |a(x,s)| <= a(0,s) + Chat_a M ,  in units of M."""
    br1 = (TWO_RA*lam/sigma_star + math.pi*lam/8.0) if sigma_star > 0 else float("inf")
    br2 = TWO_RA*math.exp(c_G)*E0
    return (C_FAR + C_INNER)*lam + min(br1, br2)

def gamma_off(L, c, lam, E0, G0, sigma_star=0.0, p_cap=3.0, damp=0.5, tol=1e-12, itmax=4000):
    """Solve  C'' = 2 Chat_a + COEF_G e^{p c_G} G0 + lam ,  c_G = (lam + C''/L) c ,
              A = (lam/2) L + Chat_a ,  Gam_rad = max(A, A + 2 Ghat + lam/2, 2 Chat_a + 2 Ghat + lam/2) ,
              c_R = Gam_rad c / L ,  p = min(p_cap, 1 + 2 c_R/c_G) .
       Returns None when the iteration runs away (no fixed point at that L)."""
    c_G, p = lam*c, 2.0
    Cpp = None
    for _ in range(itmax):
        Ca = chat_a(sigma_star, lam, c_G, E0)
        if p*c_G > 500.0:
            return None
        Ghat = COEF_G*math.exp(p*c_G)*G0
        Cpp_new = 2.0*Ca + Ghat + lam
        c_G_new = (lam + Cpp_new/L)*c
        A = 0.5*lam*L + Ca
        Gam_rad = max(A, A + 2.0*Ghat + lam/2.0, 2.0*Ca + 2.0*Ghat + lam/2.0)
        c_R = Gam_rad*c/L
        p_new = min(p_cap, 1.0 + 2.0*c_R/c_G_new)
        if c_G_new > 60.0 or not math.isfinite(c_G_new):
            return None
        d = max(abs(c_G_new - c_G), abs(p_new - p))
        c_G = c_G + damp*(c_G_new - c_G)
        p = p + damp*(p_new - p)
        Cpp = Cpp_new
        if d < tol:
            Ca = chat_a(sigma_star, lam, c_G, E0)
            Ghat = COEF_G*math.exp(p*c_G)*G0
            Cpp = 2.0*Ca + Ghat + lam
            return {"C_pp": Cpp, "c_G": c_G, "p": p, "Ghat": Ghat, "Chat_a": Ca,
                    "c_R": Gam_rad*c/L, "Gam_rad": Gam_rad, "L": L}
    return None

def L_gamma_star(c, lam, E0, G0, sigma_star=0.0, lo=1.0, hi=1e12):
    """the smallest L at which the (Gamma-off) fixed point exists"""
    if gamma_off(hi, c, lam, E0, G0, sigma_star) is None:
        return None
    for _ in range(200):
        mid = math.sqrt(lo*hi)
        if gamma_off(mid, c, lam, E0, G0, sigma_star) is None:
            lo = mid
        else:
            hi = mid
        if hi/lo < 1.0 + 1e-9:
            break
    return hi

# ------------------------------------------------------------------ C_R over the angle
def C_R_of_phi(phi, lam_om, c_G, E0, Ghat, ca_flat=None, grad_sinphi=True, nq=20001):
    """|R(X,s)| <= C_R(phi) M |X| , rho = 1, with the geometry kept (refute-u1 F3) and
       U2's axis-safe collar and proved r|grad a|:

         R_y  : Chat_a(sin phi) * sin phi
         term A (2 int_0^|z| [a - frak_a(rho_zeta)] ) : 2 int_0^|z| Chat_a(sin phi_zeta) dzeta
         term B (2 int_0^|z| [frak_a(rho_zeta) - frak_a(rho)]) :
                 lam_om * ( |z| - r arctan(|z|/r) )                    [ |d frak_a/dlog rho| <= M_s/2 ]
         term C (int_0^|z| r d_r a) : Ghat * r * log((|z|+1)/r)        [ r|grad a| <= Ghat M sin phi ]
                 or Ghat * |z|                                          [ flat, for a measured Ghat ]
    """
    s = math.sin(phi)
    z = abs(math.cos(phi))
    ca = (lambda sp: ca_flat) if ca_flat is not None else \
         (lambda sp: chat_a(sp, lam_om, c_G, E0))
    Ry = ca(s)*s
    if z <= 0.0:
        return Ry
    zg = np.linspace(0.0, z, nq)
    rz = np.sqrt(s*s + zg*zg)
    sp = s/rz
    if ca_flat is not None:
        cav = np.full_like(zg, ca_flat)
    else:
        br1 = TWO_RA*lam_om/sp + math.pi*lam_om/8.0
        br2 = TWO_RA*math.exp(c_G)*E0
        cav = (C_FAR + C_INNER)*lam_om + np.minimum(br1, br2)
    tA = 2.0*float(np.trapz(cav, zg))
    tB = lam_om*(z - s*math.atan(z/s))
    tC = Ghat*(s*math.log((z + 1.0)/s) if grad_sinphi else z)
    return Ry + tA + tB + tC

def C_R_sup(lam_om, c_G, E0, Ghat, ca_flat=None, grad_sinphi=True, phi_lo_deg=0.0, n=1200):
    lo = max(phi_lo_deg*math.pi/180.0, 1e-4)
    phis = np.linspace(lo, math.pi/2.0, n)
    vals = [C_R_of_phi(p, lam_om, c_G, E0, Ghat, ca_flat, grad_sinphi, nq=4001) for p in phis]
    k = int(np.argmax(vals))
    # refine
    a = phis[max(k-1, 0)]
    b = phis[min(k+1, n-1)]
    ph2 = np.linspace(a, b, 201)
    v2 = [C_R_of_phi(p, lam_om, c_G, E0, Ghat, ca_flat, grad_sinphi, nq=20001) for p in ph2]
    k2 = int(np.argmax(v2))
    return float(v2[k2]), float(ph2[k2]*180/math.pi)

# ------------------------------------------------------------------ the radial-ramp e-folds
def ell_ramp(eps_r=0.25, half=60.0, n=2000001):
    """e-folds lost at ONE tanh edge of Theta against the sharp indicator:
       int (1_{u<0} - (1/2)(1 - tanh((u+eps_r)/eps_r))) du  -- the outer edge of the shell.
       (The inner edge lies inside the label cut ell_loss and is not charged twice.)"""
    u = np.linspace(-half, half, n)
    Th = 0.5*(1.0 - np.tanh((u + eps_r)/eps_r))
    return float(np.trapz(np.where(u < 0.0, 1.0, 0.0) - Th, u))

# ------------------------------------------------------------------ main
if __name__ == "__main__":
    T1 = json.load(open("t1_results.json"))
    OUT = {"kernel_constants": {
        "C_K_3_over_8pi2": C_K, "S4": S4, "C_K_times_S4": CK_S4,
        "riesz_pi4_over_2": RIESZ, "riesz_quadrature": riesz_check(),
        "coef_3pi2_over_16": COEF_G, "R_A": R_A, "two_R_A": TWO_RA,
        "far_near_zodd_quoted": {"C_far": C_FAR, "C_inner": C_INNER}}}
    OUT["kernel_constants"]["riesz_relerr"] = abs(
        OUT["kernel_constants"]["riesz_quadrature"] - RIESZ)/RIESZ

    # ---- CONTROL: reproduce U2's own fixed point from U2's own datum norms -------------
    cU2 = 2.0*math.log(1.5)
    ctl = {}
    for sig, tag in [(0.5, "sigma*=1/2"), (0.0, "sigma*=0 (R^5)")]:
        r = gamma_off(1e4, cU2, 1.5, 7.66060196, 20.9197070, sigma_star=sig)
        ctl[tag] = r
    ctl["u2_reported"] = {"sigma*=1/2": {"C_pp": 582.24, "c_G": 1.2636, "p": 2.1068,
                                         "L_star": 4997.6},
                          "sigma*=0 (R^5)": {"C_pp": 801.05, "c_G": 1.2814, "p": 2.1097,
                                             "L_star": 5313.8}}
    ctl["relerr_C_pp"] = {k: abs(ctl[k]["C_pp"] - ctl["u2_reported"][k]["C_pp"])
                          / ctl["u2_reported"][k]["C_pp"] for k in ["sigma*=1/2", "sigma*=0 (R^5)"]}
    ctl["L_gamma_star_u2_datum_sigma0"] = L_gamma_star(cU2, 1.5, 7.66060196, 20.9197070, 0.0)
    ctl["L_gamma_star_u2_datum_sigma_half"] = L_gamma_star(cU2, 1.5, 7.66060196, 20.9197070, 0.5)
    ctl["Chat_a_lam1_sigma_half_at_c_G_eq_c"] = chat_a(0.5, 1.0, 1.0*cU2, 7.66060196)
    ctl["assembly_C_a_30deg"] = 8.697888
    OUT["control_reproduce_u2"] = ctl

    # ---- PRODUCTION: the theorem's own datum -------------------------------------------
    C = T1["clock"]
    CSTAR = C["c_star"]
    LAM_OM = 1.5
    E0 = T1["datum_norms"]["THEOREM_datum"]["E0_sup_rho_eta"]
    G0 = T1["datum_norms"]["THEOREM_datum"]["G0_sup_rho2_grad_eta"]
    OUT["theorem_datum_inputs"] = {"c_star": CSTAR, "E0": E0, "G0": G0, "lam_om": LAM_OM,
                                   "lam_max": C["lam_max_apriori"], "kappa": C["kappa_delta"]}
    OUT["ell_ramp_outer_edge"] = ell_ramp()

    tab = {}
    for L in [1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6, 1e7, 1e9]:
        r = gamma_off(L, CSTAR, LAM_OM, E0, G0, sigma_star=0.0)
        tab["L=%g" % L] = r
    OUT["gamma_off_fixed_point"] = tab
    OUT["L_gamma_star_theorem_datum"] = L_gamma_star(CSTAR, LAM_OM, E0, G0, 0.0)
    OUT["L_gamma_star_theorem_datum_sigma_half"] = L_gamma_star(CSTAR, LAM_OM, E0, G0, 0.5)
    OUT["L_gamma_star_if_G0_were_DC"] = L_gamma_star(CSTAR, LAM_OM, 7.66060196, 20.9197070, 0.0)

    # ---- C_R at a representative large L ------------------------------------------------
    ref = gamma_off(1e6, CSTAR, LAM_OM, E0, G0, sigma_star=0.0)
    cG, Gh = ref["c_G"], ref["Ghat"]
    G_MEAS = 0.98142          # MEASURED: u2 sec.9, max r|grad a|, lam=3/2, C^1-perturbed map
    CA_MEAS = 0.069           # MEASURED: R9, material-point strain offset
    cr = {}
    for tag, kw in [("proved_supALLphi", dict(ca_flat=None, grad_sinphi=True, phi_lo_deg=0.0)),
                    ("proved_cone_phi>=delta", dict(ca_flat=None, grad_sinphi=True,
                                                    phi_lo_deg=7.5)),
                    ("proved_Ca_measured_supALLphi", dict(ca_flat=CA_MEAS, grad_sinphi=True,
                                                          phi_lo_deg=0.0)),
                    ("measured_supALLphi", dict(ca_flat=CA_MEAS, grad_sinphi=False,
                                                phi_lo_deg=0.0))]:
        gh = G_MEAS if tag.startswith("measured") else Gh
        v, arg = C_R_sup(LAM_OM, cG, E0, gh, **kw)
        cr[tag] = {"C_R": v, "argmax_phi_deg": arg, "Ghat_used": gh, "c_G": cG, "L": 1e6}
    cr["C_R_at_phi0_30deg_proved"] = C_R_of_phi(math.pi/6, LAM_OM, cG, E0, Gh)
    cr["Chat_a_sup_all_phi_proved"] = chat_a(0.0, LAM_OM, cG, E0)
    cr["Chat_a_at_phi0_proved"] = chat_a(0.5, LAM_OM, cG, E0)
    cr["Chat_a_lam1_phi0_control_vs_assembly_8.697888"] = chat_a(0.5, 1.0, cG, E0)
    OUT["C_R"] = cr

    with open("t2_results.json", "w") as fh:
        json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
    print(json.dumps(OUT, indent=1, default=str))
