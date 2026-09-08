"""
p5 -- (V-strain), item M3, for the mollified datum: eps_a is finite and priced.

fix3 sec.1 PROVED, for any datum:
  * Lemma 1 (the kernel weight): int_{rho_1<|x|<rho_2} |Kcal|/|x| dx_5 = (3/8) log(rho_2/rho_1),
    exactly, with Kcal = d_z G_5 ;
  * Lemma 2 (transport with a source): sup |x||w(s)| <= nu e^{c_R} int_0^s Psi ,  w = eta - eta^tr ;
  * the theorem of sec.1.4:
        eps_a = (3/8)(1 + 1/(4L)) * Psibar * c * e^{c_R} / (s^2 L) ,   s = 1/sin delta ,
    conditional on (H-Delta-eta), whose only input is  Psibar = sup rho |Lap_5 eta| ;
  * the propagation (sec.1.7):  Psi_2(s) <= e^{c_R + 2 c_Gamma}[Psi_2(0) + s K_2 sqrt(G_inf N_2)],
    Psi <= sqrt(5) Psi_2 .

For THEOREM_S3's KINKED datum Psi(0) = +infinity: Lap_5 eta_0 carries a Dirac mass on each of
the four kink cones.  For the sigma-mollified datum of p1 it is finite, Psi_sigma(0) = C_kink/sigma,
and M3 closes.  This file computes eps_a(L) at the SELF-CONSISTENT window of p4 and compares it
with the budget's own eps(L).  e^{c_R} is majorised by e^{c_G} throughout (conservative, as in
fix3 g1 section D).

Outputs -> p5_results.json
"""
import json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
src = open(os.path.join(HERE, "p4_budget.py")).read().split('if __name__ == "__main__":')[0]
NS = {"__file__": os.path.join(HERE, "p4_budget.py"), "__name__": "p4mod"}
exec(compile(src, "p4_budget.py", "exec"), NS)
os.chdir(HERE)
P1 = NS["P1"]; Datum = NS["Datum"]; bind = NS["bind"]
assemble_c = NS["assemble_c"]; self_consistent = NS["self_consistent"]
T3 = NS["T3"]

DEG = math.pi/180.0
DELTA, DM = 7.5*DEG, 5.0*DEG
SDEL = 1.0/math.sin(DELTA)
W_EFOLD = 0.375
ENV = 0.25

def eps_a_smooth(L, c, c_G, Psi0):
    return W_EFOLD*(L + ENV)*(c/(SDEL**2*L))*math.exp(c_G)*Psi0/L

def eps_a_kink_sqrt(L, c, c_G, jump, C=1.0):
    supw = 2.0*C*jump*math.sqrt(c/L)/SDEL*math.exp(c_G)
    return W_EFOLD*(L + ENV)*supw/L


if __name__ == "__main__":
    SG = float(sys.argv[1])
    SHIFT = float(sys.argv[2]) if len(sys.argv) > 2 else 2.9135782
    P1RES = json.load(open(os.path.join(HERE, "p1_results.json")))
    row = P1RES["rows"]["sigma=%g" % SG]
    D = Datum(SG, SHIFT)
    E0, G0 = row["E0"], row["Gfrak0"]
    CAP = 1.0 + row["eps_cap_certified"]
    bind(D, E0, G0)
    OUT = {"sigma": SG, "Psi_sigma_0": row["Psi0"], "C_kink": row["C_kink"],
           "kinked_Psi_0": "infinity (Dirac masses on the four kink cones; the a.c. part is "
                           "%.4f)" % P1RES["rows"]["kinked"]["Psi0"],
           "formula": "eps_a = (3/8)(1 + 1/(4L)) Psi_sigma(0) c e^{c_G}/(s^2 L), s = 1/sin delta",
           "s_delta": SDEL}
    tab = []
    for L in (1e4, 1e5, 1e6, 1e7):
        for col in ("proved", "measured"):
            s = self_consistent(D, L, col, E0, G0, CAP, fs=NS["FS_F4"])
            if not s["self_consistent"]:
                tab.append({"L": L, "column": col, "eps": float("inf")})
                continue
            A = assemble_c(D, L, s["c"], col, E0, G0, f=s["best_f"])
            ea = eps_a_smooth(L, s["c"], A["c_G"], row["Psi0"])
            ek = eps_a_kink_sqrt(L, s["c"], A["c_G"], 1.0/DM)
            tab.append({"L": L, "column": col, "c": s["c"], "c_over_c_star": s["c_over_c_star"],
                        "c_G": A["c_G"], "eps_budget": A["eps"], "eps_a_mollified": ea,
                        "eps_a_kinked_sqrt_route": ek,
                        "eps_a_over_eps": ea/A["eps"] if np.isfinite(A["eps"]) else None,
                        "eps_with_eps_a_added":
                            assemble_c(D, L, s["c"], col, E0, G0, f=s["best_f"],
                                       eps_a_extra=ea)["eps"]})
            print("L=%.0e %-8s c/c_*=%.6f c_G=%.5f eps=%.6g eps_a=%.4e (%.4f%% of eps)"
                  % (L, col, s["c_over_c_star"], A["c_G"], A["eps"], ea,
                     100.0*ea/A["eps"] if np.isfinite(A["eps"]) else float('nan')), flush=True)
    OUT["table"] = tab
    # ---- Psi_2(0) = sup rho |grad^2 eta_0|_F , the datum input of fix3's propagation lemma,
    # computed on fix2's DatumS3 with the angular factor replaced by the mollified profile.
    sys.path.insert(0, os.path.join(HERE, "imported"))
    from f2_datum_s3 import DatumS3                                   # noqa: E402

    class DatumMoll(DatumS3):
        def __init__(self, prof, **kw):
            DatumS3.__init__(self, **kw)
            self.p = prof

        def W(self, phi, k=0):
            ph = np.abs(np.asarray(phi, dtype=float))
            return np.interp(ph, self.p.phi, (self.p.w0, self.p.w1, self.p.w2)[k])

    DM_ = DatumMoll(D.prof, delta_deg=7.5, dm_deg=5.0, eps_r=0.25, L=40.0)
    ph_ = np.linspace(1e-4, math.pi-1e-4, 40001)
    best2, bestF = 0.0, None
    for u_ in np.linspace(-4.0, 6.0, 1001):
        rho_ = math.exp(u_)
        r_ = np.maximum(rho_*np.sin(ph_), 1e-300)
        z_ = rho_*np.cos(ph_)
        hf = DM_.hess_F(r_, z_)
        v_ = float(np.max(hf))*rho_
        if v_ > best2:
            best2, bestF = v_, float(u_)
    OUT["Psi_2_0_sup_rho_hessF"] = best2
    OUT["Psi_2_0_argmax_u"] = bestF
    OUT["Psi_over_sqrt5_Psi2"] = row["Psi0"]/(math.sqrt(5.0)*best2)
    print("Psi_2(0) = sup rho |grad^2 eta_0|_F = %.4f at u = %s ; Psi_sigma(0)/(sqrt5 Psi_2(0))"
          " = %.4f" % (best2, bestF, OUT["Psi_over_sqrt5_Psi2"]), flush=True)

    # the propagation half of (H-Delta-eta), fix3 sec.1.7
    G_INF = row["grad_eta0_inf"]
    K2 = T3.RECORD["K2_proved_window"]
    prop = {}
    for L in (1e4, 1e6):
        c = D.c_star
        tau = c/L
        prop["L=%g" % L] = {"tau_K2_sqrt(G_inf N2)": tau*K2*math.sqrt(G_INF*G0),
                            "Psi_2_0": best2,
                            "relative": tau*K2*math.sqrt(G_INF*G0)/best2}
    OUT["propagation"] = prop
    OUT["grad_eta0_inf"] = G_INF
    json.dump(OUT, open(os.path.join(HERE, "p5_results.json"), "w"), indent=1,
              sort_keys=True, default=str)
    print(json.dumps({k: v for k, v in OUT.items() if k != "table"}, indent=1, default=str))
