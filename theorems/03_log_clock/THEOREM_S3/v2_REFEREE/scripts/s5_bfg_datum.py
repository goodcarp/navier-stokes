"""
s5 -- THE BFG HYPOTHESES, re-checked independently for the sigma-mollified datum, and the
      log Re_E dictionary arithmetic.

BFG (arXiv:1704.05546v4) Theorem 10 needs: omega_0 in L^2 cap L^inf, and a unique mild
solution in C_w([0,T], L^inf).  THEOREM_S3_v2 sec.4 asserts, for the mollified profile:
  sup|A_sigma|/N_sigma = 1.0 exactly ;  z-oddness ;  ||omega_0||_2^2/(M^2 R^3) = 1.4381082627 ;
  int_{-1}^1 w(t)^2 dt = 1.8239293735 <= 2 ;  decay rho^{-8} / rho^{8} .
All are recomputed here from the independent instrument of s1.
Outputs -> s5_results.json
"""
import json, math, os
import numpy as np
from scipy.integrate import quad
import s1_profile_indep as P

HERE = os.path.dirname(os.path.abspath(__file__))
DEG = math.pi/180.0
N_SIG = 1.0124508490
EPS_R = 0.25

def Theta(u, L, er=EPS_R):
    return 0.5*(math.tanh((u-er)/er) - math.tanh((u-L+er)/er))

if __name__ == "__main__":
    R = {"N_sigma_used": N_SIG}
    brk = [7.5*DEG, math.pi/2 - 5.0*DEG, math.pi/2, math.pi/2 + 5.0*DEG, math.pi - 7.5*DEG]

    # (i) sup |A_sigma|/N   and z-oddness
    xs = np.linspace(1e-7, math.pi-1e-7, 40001)
    a = np.array([P.As(x) for x in xs])/N_SIG
    R["sup_abs_A_over_N"] = float(np.abs(a).max())
    odd = np.array([abs(P.As(math.pi-x) + P.As(x)) for x in np.linspace(0.05, math.pi/2, 400)])
    R["max_z_oddness_residual"] = float(odd.max()/N_SIG)

    # (ii) int_{-1}^1 w(t)^2 dt = int_0^pi (A_sigma/N)^2 sin phi dphi
    edges = [0.0] + brk + [math.pi]
    I2 = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        v, _ = quad(lambda p: (P.As(p)/N_SIG)**2*math.sin(p), lo, hi, limit=200,
                    epsabs=1e-12, epsrel=1e-11)
        I2 += v
    R["int_w2_dt"] = I2
    R["fix4_int_w2_dt"] = 1.8239293735
    R["parseval_hypothesis_satisfied"] = bool(I2 <= 2.0)

    # (iii) ||omega_0||_2^2/(M^2 R^3) = 2 pi (int_0^pi (A/N)^2 sin) * e^{-3L} int Theta^2 e^{3u} du
    rad = {}
    for L in (10.0, 20.0, 40.0):
        v, _ = quad(lambda u: Theta(u, L)**2*math.exp(3.0*(u-L)), -60.0, L+8.0,
                    limit=400, epsabs=1e-14, epsrel=1e-12)
        rad["L=%g" % L] = v
    R["radial_factor"] = rad
    R["omega0_L2sq_over_M2R3"] = {k: 2.0*math.pi*I2*v for k, v in rad.items()}
    R["fix4_omega0_L2sq_over_M2R3"] = 1.4381082627

    # kinked control
    def Ak(p):
        pax = min(p, math.pi-p); eq = abs(p-math.pi/2)
        return math.copysign(1.0, math.cos(p))*min(1.0, pax/(7.5*DEG))*min(1.0, eq/(5.0*DEG))
    I2k = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        v, _ = quad(lambda p: Ak(p)**2*math.sin(p), lo, hi, limit=200, epsabs=1e-13,
                    epsrel=1e-12)
        I2k += v
    R["int_w2_dt_kinked"] = I2k
    R["fix4_int_w2_dt_kinked"] = 1.8751740891
    R["omega0_L2sq_over_M2R3_kinked_L40"] = 2.0*math.pi*I2k*rad["L=40"]
    R["fix4_omega0_L2sq_kinked"] = 1.4785130340

    # (iv) the log Re_E dictionary
    s_rec = 1.0/math.sin(7.5*DEG)
    R["log_Re_E_dictionary"] = {
        "s": s_rec, "two_log_s": 2.0*math.log(s_rec),
        "C_E_quoted": 0.0551906693,
        "shift_recomputed": 2.0*math.log(s_rec) + 0.4*math.log(0.0551906693),
        "shift_quoted": 2.9135781820,
        "logLambda_star_from_quoted_Lstar":
            2*1424610.4953 + 2.0*math.log(s_rec) + 0.4*math.log(0.0551906693),
        "logLambda_star_quoted": 2849223.9041}

    # (v) the tail exponents claimed in sec.1.1
    R["tail_check"] = {
        "Theta_inner_exact_form": "rho^8/(e^2 + rho^8) at eps_r = 1/4, since 2/eps_r = 8",
        "Theta(u=-3,L=40)": Theta(-3.0, 40.0),
        "rho8_over_e2_plus_rho8 at rho=e^{-3}": math.exp(-24.0)/(math.e**2 + math.exp(-24.0)),
        "Theta(u=L+3,L=40)": Theta(43.0, 40.0),
        "predicted_outer e^{-8(u-L+eps_r)}": math.exp(-8.0*(43.0-40.0+0.25))}

    json.dump(R, open(os.path.join(HERE, "s5_results.json"), "w"), indent=1, default=str)
    print(json.dumps(R, indent=1, default=str))
