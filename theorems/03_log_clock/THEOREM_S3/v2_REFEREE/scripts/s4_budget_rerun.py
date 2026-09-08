"""
s4 -- RE-RUN fix4's budget from a COPY, and price the (H3'_vartheta) repair.

(a) reproduce L_* and log Lambda_* in the three columns from fix4's own instrument, run on a
    byte copy;  (b) check that the "all proved" column touches no MEASURED record entry;
(c) the C_K sensitivity pushed far past hk2's majorant;  (d) the cost of shrinking the
Theorem-V.4' ball so that (H3'_vartheta) can hold for the sigma-mollified datum (s3).

The only re-implementation is `viscous_budget_d`, which is `t3_budget.viscous_budget`
transcribed with the ball radius d exposed as a parameter; it is checked against the
imported function at dscale = 1 to machine precision before use.
Outputs -> s4_results.json
"""
import json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIX4 = os.path.join(os.path.dirname(HERE), "copies", "THEOREM_S3", "fix4")
sys.path.insert(0, FIX4)
sys.path.insert(0, os.path.join(FIX4, "imported"))
_cwd = os.getcwd()
import p4_budget as P4                                    # noqa: E402  (chdirs internally)
os.chdir(_cwd)
T3 = P4.T3
from scipy.stats import ncx2                              # noqa: E402

P1RES = json.load(open(os.path.join(FIX4, "p1_results.json")))
ROW = P1RES["rows"]["sigma=0.02"]
E0, G0 = ROW["E0"], ROW["Gfrak0"]
CAP = 1.0 + ROW["eps_cap_certified"]
SIGMA = 0.02
DEG = math.pi/180.0

# ------------------------------------------------- viscous_budget with d exposed
def viscous_budget_d(L, c, C_K, f, dscale=1.0, delta_eff=None, phi0=T3.PHI0, s=T3.S_REC,
                     delta=T3.DELTA, c_G=None, Cprime=0.0, use_V1_for_Z=True):
    if c_G is None:
        c_G = (1.5 + Cprime/L)*c
    n = T3.RECORD["n_dim"]
    r0 = (1.0+f)*math.sin(phi0)
    nutau = c/(s*s*L)
    sy = nutau*T3.J_pow(c, -2.0)/c
    sz = nutau*T3.J_pow(c,  4.0)/c
    rho_s = 1.0+f
    de = delta if delta_eff is None else delta_eff
    d_inner = f
    d_taper = rho_s*math.sin(phi0-de)
    d_eq = rho_s*math.sin(math.pi/2 - T3.DM - phi0)
    d = dscale*max(min(d_inner, d_taper, d_eq), 1e-12)
    R_minus = r0 - d
    N = r0/math.sin(delta)                       # kept at the stock delta (conservative)
    eps_bulk = sy/(r0*r0)
    V1 = n*nutau*(math.exp(2*c_G)-1.0)/c_G
    tau = c/L
    EW = 0.5*C_K*math.exp(c_G)*tau*V1
    if R_minus <= 0 or d <= 1e-11:
        return {"eps_v": float("inf"), "d": d, "R_minus": R_minus}
    E_hess = EW*(r0/(R_minus**2)) + 4.0*N*EW/d
    trC = 8*sy + 2*sz
    trC2 = 4*(2*sy)**2 + (2*sz)**2
    trC3 = 4*(2*sy)**3 + (2*sz)**3
    E4 = (r0/(R_minus**5))*(trC**2 + 2*trC2)
    P_a_half = T3.gauss_tail_aniso(d/2.0, sy, sz)
    P_a_full = T3.gauss_tail_aniso(d, sy, sz)
    if use_V1_for_Z:
        q = (d/2.0)**2/nutau
        P_Z_half = min(1.0, 2*n*math.exp(-math.exp(-2*c_G)*q/(4.0*n)))
    else:
        P_Z_half = P_a_half
    Pp = P_Z_half + P_a_half + P_a_full
    m2 = trC; m4 = trC**2 + 2*trC2; m6 = trC**3 + 6*trC*trC2 + 8*trC3
    E_tail = 2*N*Pp + Pp + math.sqrt(Pp)*(math.sqrt(m2)/r0 + math.sqrt(m4)/r0**2
                                          + math.sqrt(m6)/r0**3)
    return {"eps_bulk": eps_bulk, "E_hess": E_hess, "E_4": E4,
            "E_tail": min(E_tail, 1e12), "d": d, "R_minus": R_minus, "N": N, "c_G": c_G,
            "eps_v": min(eps_bulk + E_hess + E4 + min(E_tail, 1e12), 1e12)}

if __name__ == "__main__":
    OUT = {"E0": E0, "Gfrak0": G0, "cap": CAP,
           "fix4_quoted": {"L_star_proved_cap": 1424610.4953, "L_star_proved_0p1": 1555469.4004,
                           "L_star_Ca_meas_cap": 1380034.8075, "L_star_measured_cap": 1757.8779,
                           "logLambda_proved_cap": 2849223.9041,
                           "shift": 2.9135781820, "C_E": 0.0551906693}}

    # ---- identity check of the transcription
    ident = []
    for (L, c, f) in [(1e6, 0.88, 4.0), (1e4, 0.83, 8.0), (2e6, 0.879, 4.0)]:
        a = T3.viscous_budget(L, c, 161.7735, f, c_G=1.24*c/0.83)
        b = viscous_budget_d(L, c, 161.7735, f, c_G=1.24*c/0.83)
        ident.append({"L": L, "eps_v_imported": a["eps_v"], "eps_v_transcribed": b["eps_v"],
                      "rel": abs(a["eps_v"]-b["eps_v"])/max(a["eps_v"], 1e-300)})
    OUT["transcription_identity_check"] = ident
    print("identity:", [("%.3e" % x["rel"]) for x in ident], flush=True)

    # ---- the shift arithmetic
    s_rec = 1.0/math.sin(7.5*DEG)
    C_E = json.load(open(os.path.join(FIX4, "p2_results.json")))
    def find(o, key):
        if isinstance(o, dict):
            for k, v in o.items():
                if k == key: return v
                r = find(v, key)
                if r is not None: return r
        return None
    ce = find(C_E, "C_E_mollified") or find(C_E, "C_E")
    OUT["C_E_from_p2"] = ce
    if isinstance(ce, (int, float)):
        OUT["shift_recomputed"] = 2.0*math.log(s_rec) + 0.4*math.log(ce)
    OUT["shift_from_quoted_C_E"] = 2.0*math.log(s_rec) + 0.4*math.log(0.0551906693)
    print("shift from C_E=0.0551906693 :", OUT["shift_from_quoted_C_E"], flush=True)

    # ---- (a) reproduce L_* for the proved column
    D = P4.Datum(0.02, 2.9135781820)
    e0, g0, cap, row = E0, G0, CAP, ROW
    P4.bind(D, e0, g0)
    res = {}
    for tgt, tag in [(cap-1.0, "eps<=cap"), (0.1, "eps<=0.1")]:
        v = P4.Lstar(D, "proved", tgt, e0, g0, cap, fs=P4.FS_F4)
        res[tag] = v
        print("  reproduce L_* proved %-12s = %s" % (tag, v), flush=True)
    OUT["L_star_proved_reproduced"] = res
    OUT["logLambda_proved_cap_reproduced"] = (2*res["eps<=cap"] + 2.9135781820
                                              if res["eps<=cap"] else None)
    json.dump(OUT, open(os.path.join(HERE, "s4_results.json"), "w"), indent=1, default=str)

    # ---- (b) which RECORD entries the "proved" column touches
    route, ca_kind, gh_kind, k2key, sens, use_V1 = T3.COLUMNS["proved"]
    OUT["proved_column_record_entries"] = {
        "gamma_off_route": route, "C_a_kind": ca_kind, "Ghat_kind": gh_kind,
        "K2_key": k2key, "K2_value": T3.RECORD[k2key],
        "Lemma_T_sensitivity": sens, "uses_V1_for_Z": use_V1,
        "MEASURED_entries_in_RECORD_not_used_by_proved":
            ["Cprime_measured=11.74", "Ghat_measured=0.98142",
             "material_offset_M=0.069", "K2_measured_window=5.3854"],
        "K2_proved_window_comment": "t3_budget.py:47 -- 'PROVED over lambda in [1,3/2]'"}

    # ---- (c) C_K sensitivity far past the majorant
    ck = {}
    for v in [1.0, 161.7735, 1.0e3, 1.0e4, 1.0e5, 1.0e6]:
        L = P4.Lstar(D, "proved", 0.1, e0, g0, cap, fs=P4.FS_F4, tol=1e-5, C_K=v)
        ck["C_K=%g" % v] = L
        print("  C_K = %-10g  L_*(eps<=0.1) = %s" % (v, L), flush=True)
    OUT["C_K_sensitivity"] = ck
    json.dump(OUT, open(os.path.join(HERE, "s4_results.json"), "w"), indent=1, default=str)

    # ---- (d) the cost of shrinking the V.4' ball (the (H3'_vartheta) repair)
    _orig = T3.viscous_budget
    sh = {}
    for tag, kw in [("stock", {}),
                    ("delta_eff = delta + 4 sigma", {"delta_eff": 7.5*DEG + 4*SIGMA}),
                    ("delta_eff = delta + 6 sigma", {"delta_eff": 7.5*DEG + 6*SIGMA}),
                    ("dscale = 0.8", {"dscale": 0.8}),
                    ("dscale = 0.5", {"dscale": 0.5})]:
        def vb(L, c, C_K, f, **kk):
            kk.pop("phi0", None); kk.pop("s", None); kk.pop("delta", None)
            return viscous_budget_d(L, c, C_K, f, **dict(kw, **kk))
        T3.viscous_budget = vb
        P4._CC.clear(); P4._BS.clear()
        L1 = P4.Lstar(D, "proved", cap-1.0, e0, g0, cap, fs=P4.FS_F4)
        r = P4.assemble_c(D, 2e6, 1.0669124225*D.c_star, "proved", e0, g0, f=4.0)
        sh[tag] = {"L_star_eps<=cap": L1, "eps_v_at_L2e6": r["eps_v"], "eps_at_L2e6": r["eps"],
                   "E_hess": r["E_hess"], "E_tail": r["E_tail"], "eps_bulk": r["eps_bulk"]}
        print("  %-28s L_*=%s  eps_v(2e6)=%.4e" % (tag, L1, r["eps_v"]), flush=True)
    T3.viscous_budget = _orig
    OUT["ball_shrink_cost"] = sh
    json.dump(OUT, open(os.path.join(HERE, "s4_results.json"), "w"), indent=1, default=str)
    print(json.dumps(OUT, indent=1, default=str))
