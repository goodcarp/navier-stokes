"""
g5 -- UNIT 1, the propagation half: Psi(s) from Psi(0).

The weighted-maximum-principle family of pmax-h3v is P1 (weight rho^1 on |eta|) and P2
(weight rho^2 on |grad eta|).  The (V-strain) route needs one more member, on the HESSIAN,
with weight rho^1:

    Psi_2(s) := sup_x |x| |grad^2_5 eta(x,s)|_F ,      Psi(s) := sup_x |x| |Lap5 eta| <= sqrt(5) Psi_2 .

From Lemma D of pmax-h3v (differentiate the equation twice),

    D_t (d_i d_j eta) = nu Lap5 (d_i d_j eta) - (d_i d_j b).grad eta
                        - (d_i b).grad(d_j eta) - (d_j b).grad(d_i eta)

so, by the regularised Kato inequality of pmax-h3v Lemma C applied to the Hessian and Lemma B
with k = 1 (whose zeroth-order coefficient k(k-3) = -2 <= 0 at n = 5),

    D_t ( R_eps |H| )  <=  (Gamma_rad + 2 Gamma) R_eps|H|  +  R_eps K_2 |grad eta|
                           + nu [ Lap - (2/R_eps^2) x.grad - q ] (R_eps|H|) ,   q >= 0 .

The source is the only new object.  It is NOT bounded by K_2 N_2/rho alone (that blows up at
the origin) nor by K_2 rho ||grad eta||_inf alone (that blows up at infinity); the product of
the two weighted bounds is what closes it:

    rho K_2 |grad eta|  <=  K_2 min( rho G_inf , N_2/rho )  <=  K_2 sqrt( G_inf N_2 ) ,
    G_inf(s) := ||grad eta(.,s)||_inf <= e^{c_Gamma(s)} G_inf(0) ,        (max principle, k=0)
    N_2(s)   := sup rho^2|grad eta|   <= e^{p c_G} Gfrak_0 M .            (P2)

Hence, with c_R = int Gamma_rad and c_Gamma = int Gamma over [0,s],

    Psi_2(s) <= e^{c_R + 2 c_Gamma} [ Psi_2(0) + s K_2 sqrt( G_inf N_2 ) ] .

This script evaluates the second term against the first, at the window's own constants, so the
size of the propagation is on the record.  Units rho_0 = M = 1.

Outputs -> g5_results.json
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "fix2copy"))
import math                                                        # noqa: F811
import f1_window as F1                                             # noqa: E402

DEG = math.pi/180.0
DELTA = 7.5*DEG
SDEL = 1.0/math.sin(DELTA)
GFRAK0 = 1.0/math.sin(DELTA)**2                # PROVED (fix2 sec.4(c))
GINF0 = 20.3339                                # g1 section C, two instruments
K2HAT = 161.7735                               # hk2, PROVED over the window (t3_budget import)
P_MAX = 3.0                                    # p = 1 + 2 c_R/c_G <= 3 (pmax-h3v P2')

OUT = {"inputs": {"Gfrak_0": GFRAK0, "grad_eta_inf_0": GINF0, "K2hat": K2HAT,
                  "sqrt(Ginf*Gfrak0)": math.sqrt(GINF0*GFRAK0), "p_max": P_MAX}}
rows = []
for L in (1e4, 1e5, 1e6):
    for x, lab in ((1.0, "frozen c=c_*"), (F1.CAP_CERT, "cap c=1.1943662c_*")):
        c = x*F1.CSTAR
        col = "proved" if L > 3e4 else "measured"
        A = F1.assemble_c(L, c, col, f=4.0)
        if "c_G" not in A:
            A = F1.assemble_c(L, c, "measured", f=4.0)
            col = "measured"
        cG = A["c_G"]
        # conservative: c_R <= c_Gamma <= c_G, p <= 3
        Ginf = math.exp(cG)*GINF0
        N2 = math.exp(P_MAX*cG)*GFRAK0
        src = K2HAT*math.sqrt(Ginf*N2)
        tau = c/L                                     # M = 1
        growth = math.exp(3.0*cG)                     # e^{c_R + 2 c_Gamma} <= e^{3 c_G}
        rows.append({"L": L, "window": lab, "column": col, "c": c, "c_G": cG,
                     "G_inf(s)": Ginf, "N_2(s)": N2,
                     "source_K2_sqrt(Ginf N2)": src, "tau": tau,
                     "tau*source": tau*src, "growth_factor_e^{3c_G}": growth,
                     "Psi2_0_sigma=0.002": 4130.11/math.sqrt(5.0),
                     "propagation_over_Psi2_0":
                         tau*src/(4130.11/math.sqrt(5.0))})
        print("L=%8.3g %-20s c_G=%.5f  tau*K2*sqrt(Ginf N2) = %.6g   (Psi_2(0) at sigma=0.002 "
              "= %.1f)  ratio %.3e" % (L, lab, cG, tau*src, 4130.11/math.sqrt(5.0),
                                       tau*src/(4130.11/math.sqrt(5.0))))
OUT["rows"] = rows
OUT["conclusion"] = ("the propagation term tau K_2 sqrt(G_inf N_2) is at most "
                     "%.3g in units M/rho_0^2 at L = 1e4 and falls like 1/L, against a datum "
                     "value Psi_2(0) of order 1e3 for a sigma = 0.002 mollification: the "
                     "propagation is not the binding half of (H-Delta-eta); the datum is."
                     % max(r["tau*source"] for r in rows))
with open(os.path.join(HERE, "g5_results.json"), "w") as fh:
    json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
print("\nWROTE g5_results.json")
