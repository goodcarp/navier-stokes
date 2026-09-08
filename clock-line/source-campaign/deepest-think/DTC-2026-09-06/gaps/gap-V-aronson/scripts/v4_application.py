"""
v4 -- every number quoted in NOTE.md sections 5-7: the application of Theorem V to the
prove-lagrangian datum, and the resulting eta-loss on the innermost material shell over the
doubling window.   Pure arithmetic from the campaign's own constants; nothing is fitted.

Frame (prove-lagrangian NOTE, sections 1 and 4):
   omega^theta_0 = -M sgn(z) h_delta(phi) on rho0 < rho < R ;  eta = omega^theta/r
   L := log(R/rho0) ;  viscous floor rho0 delta = sqrt(nu/M)  =>  nu = M (rho0 delta)^2
   strain a = (M/2) log(R/rho) + O(M) ;  ||grad u||_op = 2a = M L (1+O(1/L))
   window: theta := M L T_d ;  accelerated 2(1-sqrt(2/3))/kappa = 0.73432 ; conservative
           log(3/2)/kappa = 0.81093    (kappa = 1/2)
   => c := ||grad u||_inf tau = theta  (exactly, to leading order)
"""
import math, json
OUT = {}
def rec(k, v):
    OUT[k] = v
    print(f"{k} = {v}")

# ---------------- window constants straight out of prove-lagrangian ----------------
kappa = 0.5
theta_acc  = 2*(1-math.sqrt(2/3))/kappa
theta_cons = math.log(1.5)/kappa
rec("theta_accelerated", theta_acc)
rec("theta_conservative", theta_cons)
c = theta_cons                                    # the larger, i.e. the conservative window
rec("c_window", c)
rec("e2c", math.exp(2*c)); rec("Cprime_4e2c", 4*math.exp(2*c))

deg = math.pi/180
for delta_deg, phi0_deg, rho_over_rho0 in [(7.5, 30, 1.0), (7.5, 30, 2.0), (7.5, 62.83, 2.0),
                                           (5.0, 30, 2.0), (15.0, 30, 2.0)]:
    delta = delta_deg*deg; phi0 = phi0_deg*deg; q = rho_over_rho0
    # geometry in units of rho0 ; the viscous floor sets rho0 delta = sqrt(nu/M) = 1 length unit
    r_over_rho0 = q*math.sin(phi0)
    # distances (in units of rho0) from alpha0=(q,phi0) to the boundary of the region where the
    # datum coincides with the pure untapered plateau  {rho0<rho<R, delta<phi<pi-delta}
    d_inner = q - 1.0
    d_taper = 2*q*math.sin(abs(phi0-delta)/2)
    d_equat = q*math.cos(phi0) if phi0 < math.pi/2 else 0.0
    d = min(x for x in [d_inner, d_taper, d_equat] if x > 0) if min(d_inner, d_taper, d_equat) > 0 else 0.0
    tag = f"delta={delta_deg}deg phi0={phi0_deg}deg rho={q}rho0"
    row = dict(d_inner=d_inner, d_taper=d_taper, d_equator=d_equat, d=d,
               r_over_rho0=r_over_rho0, delta=delta)
    if d > 0:
        # d^2/(nu tau) = (d/rho0)^2 * L / (delta^2 c)     [since nu tau = (rho0 delta)^2 c / L]
        row["d2_over_nutau_coeff_times_L"] = d*d/(delta*delta*c)
        row["exponent_coeff_times_L"] = d*d/(delta*delta*c*4*math.exp(2*c))
        for L in [5, 10, 20, 50, 100]:
            row[f"bound_exp_at_L{L}"] = math.exp(-row["exponent_coeff_times_L"]*L)
    # local smooth loss:   sqrt(nu tau)/r = delta sqrt(c/L) / (r/rho0)   ;  nu tau/r^2 = delta^2 c/(L (r/rho0)^2)
    row["ell_sqrt_coeff_over_sqrtL"] = delta*math.sqrt(c)/r_over_rho0
    row["ell_lin_coeff_over_L"] = delta*delta*c/r_over_rho0**2
    for L in [5, 10, 20, 50, 100]:
        row[f"ell_sqrt_at_L{L}"] = row["ell_sqrt_coeff_over_sqrtL"]/math.sqrt(L)
        row[f"ell_lin_at_L{L}"]  = row["ell_lin_coeff_over_L"]/L
    rec("CASE " + tag, row)

# ------- the brief's own conservative choice d = rho0 delta  -------
row = dict()
row["d2_over_nutau"] = "L/c"
row["exponent"] = "L/(4 e^{2c} c)"
row["coeff"] = 1.0/(4*math.exp(2*c)*c)
for L in [5, 10, 20, 50, 100, 200]:
    row[f"bound_at_L{L}"] = math.exp(-row["coeff"]*L)
rec("BRIEF_d_equals_rho0_delta", row)

# ------- amplitude prefactor:  N/|eta_0(alpha0)| -------
for delta_deg, phi0_deg, q in [(7.5, 30, 2.0), (7.5, 30, 1.0)]:
    delta = delta_deg*deg; phi0 = phi0_deg*deg
    N_over = 1.0/(1.0*math.sin(delta))              # ||eta_0||_inf = M/(rho0 sin delta)
    eta0 = 1.0/(q*math.sin(phi0))                   # |eta_0(alpha0)| = M/(q rho0 sin phi0)
    rec(f"prefactor_2N_over_eta0 delta={delta_deg} phi0={phi0_deg} rho={q}", 2*N_over/eta0)

# ------- translation of an eta-loss ell into a shift of c2 -------
# omega^theta(X(tau)) = M exp(int a) (1-ell);  reaching 3/2 needs int a >= log(3/2) - log(1-ell)
# theta shifts by -log(1-ell) ~ ell ;  c2 = 2 theta / kappa  => Delta c2 = 2 ell / kappa = 4 ell
rows = []
for L in [5, 10, 20, 50, 100]:
    ell_sqrt = 7.5*deg*math.sqrt(c)/(2*math.sin(30*deg))/math.sqrt(L)
    ell_lin  = (7.5*deg)**2*c/(2*math.sin(30*deg))**2/L
    c2acc = 2*theta_acc          # c2 = theta * (log Re_E / L) = 2 theta = 8(1-sqrt(2/3))
    rows.append(dict(L=L, ell_sqrt=ell_sqrt, dc2_sqrt=2*ell_sqrt,
                     ell_lin=ell_lin, dc2_lin=2*ell_lin,
                     dc2_move_to_2rho0=math.log(2.0)*c2acc/L))
    print("   ", rows[-1])
rec("c2_shift_table", rows)
rec("c2_accelerated_2theta", 2*theta_acc)      # = 8(1-sqrt(2/3))
rec("c2_check_8_1_minus_sqrt23", 8*(1-math.sqrt(2/3)))
rec("c2_conservative_4log32", 4*math.log(1.5))

# ------- the divergence term, sized -------
# in the 5D Lagrangian frame the ONLY trace of div5 b = 2a is the drift nu grad_alpha log J . g^{-1}
#   |grad_alpha log J| <= 2 tau e^{c} ||grad a||_inf,  ||grad a||_inf ~ M/rho0 for the plateau
# scale-free drift parameter  beta_bar := sqrt(tau/nu) * |b_L|
rows = []
for L in [5, 10, 20, 50, 100]:
    delta = 7.5*deg
    # units M=rho0=1 : nu = delta^2, tau = c/L, grad a ~ 1
    nu_ = delta**2; tau_ = c/L
    gradlogJ = 2*tau_*math.exp(c)*1.0
    bL = nu_*gradlogJ*math.exp(2*c)
    beta_bar = math.sqrt(tau_/nu_)*bL
    rows.append(dict(L=L, nu=nu_, tau=tau_, grad_logJ=gradlogJ, bL=bL, beta_bar=beta_bar))
    print("   ", rows[-1])
rec("divergence_drift_table", rows)

with open(__file__.replace('v4_application.py', 'v4_results.json'), 'w') as fh:
    json.dump(OUT, fh, indent=1, default=str)
print("\nDONE v4")
