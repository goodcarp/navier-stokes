"""
q4 -- what (H3-V)'s repair costs: THEOREM_S3's budget re-run with vartheta carried.

L-14, INVERSE CONVENTION, DECLARED.  The object under test here is THEOREM_S3's budget
instrument itself, so it is IMPORTED, not re-implemented: `copies/t3_budget.py` and
`copies/t2_gamma_CR.py` are byte copies (hashes in `copies/IMPORTED_SHA256SUMS.txt`;
`copies/u5_budget.py` is there too and reproduces the sha256 4c2cdbc32cb0866a... that
THEOREM_S3 sec.4 quotes for its own provenance).  Nothing in them is edited.  What is
re-implemented here is the CLAIM being tested: the extra error terms that Theorem V.4'
carries when (H3) is weakened to (H3'_vartheta).

THE EXTRA TERMS (derived in PROOF.md sec.B3, from Steps 1 and 2 of Theorem V.4):

  put  Q := (1+vartheta)/(1-vartheta) ,  vt := vartheta/(1-vartheta) ,  and

  eps_v'  =  eps_bulk * [ 1 + vt (9 + 2 sigma_z/sigma_y) ]        <- Step 2's Hessian contraction
                                                                      and the normalisation
           +  Q * ( E_hess + E_4 + E_tail )                       <- every V.4 remainder: each has
                                                                     a numerator (1+vartheta) and
                                                                     the denominator |eta_0(x_*)|
                                                                     >= (1-vartheta) M/r_*

  and every one of these reduces to THEOREM_S3's own eps_v at vartheta = 0.

vartheta(f, d) is read from q3_results.json.  A grid-plus-polish search returns a lower bound
on a supremum, so the budget is run at SAFETY * vartheta with SAFETY = 2, and the SAFETY = 1
and SAFETY = 4 columns are reported beside it.

Outputs -> q4_results.json
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
COPIES = os.path.join(HERE, "copies")
sys.path.insert(0, COPIES)
_cwd = os.getcwd()
os.chdir(COPIES)
import t3_budget as T3                      # noqa: E402   (module-level json.load needs this cwd)
import t2_gamma_CR as T2                    # noqa: E402
os.chdir(_cwd)

import numpy as np                          # noqa: E402

THETA = json.load(open(os.path.join(HERE, "q3_results.json")))
OUT = {"provenance": {
    "imported": "copies/t3_budget.py, copies/t2_gamma_CR.py, copies/t1_results.json",
    "hashes": open(os.path.join(COPIES, "IMPORTED_SHA256SUMS.txt")).read().strip().split("\n"),
    "edited": "nothing; the vartheta terms are added on top in this file",
}}

CSTAR = T3.CSTAR
FS = T3.FS


def vartheta_of(f, safety):
    row = THETA["f_grid_d_from_assembly_rule"].get("f=%g" % f)
    if row is None:
        return None
    return safety * float(row["vartheta"])


def assemble_theta(L, c, column, f, safety, vartheta_override=None):
    """THEOREM_S3.assemble with the (H3'_vartheta) terms added to eps_v."""
    K = T3.column_constants(L, c, column)
    if K is None:
        return {"L": L, "column": column, "f": f, "eps": float("inf"),
                "gamma_off_exists": False}
    bs = T3.bootstrap(L, c, K["C_pp"], K["C_R"], K["ca_sup"], C1=0.0, lam_mx=T3.LAM_MAX)
    vb = T3.viscous_budget(L, c, K["C_K"], f, c_G=K["c_G"], use_V1_for_Z=K["use_V1"])
    vt_raw = vartheta_of(f, safety) if vartheta_override is None else vartheta_override
    if vt_raw is None or vt_raw >= 1.0:
        eps_v = float("inf")                       # (H3'_vartheta) is vacuous at this f
        pieces = {"vartheta": vt_raw, "admissible": False}
    else:
        Q = (1.0 + vt_raw) / (1.0 - vt_raw)
        vt = vt_raw / (1.0 - vt_raw)
        ratio = vb["sigma_z"] / vb["sigma_y"]
        extra_bulk = vb["eps_bulk"] * vt * (9.0 + 2.0 * ratio)
        eps_v = (vb["eps_bulk"] + extra_bulk
                 + Q * (vb["E_hess"] + vb["E_4"] + vb["E_tail"]))
        pieces = {"vartheta": vt_raw, "admissible": True, "Q": Q, "vt": vt,
                  "sigma_z_over_sigma_y": ratio, "extra_bulk": extra_bulk,
                  "eps_v_theta0": vb["eps_v"], "eps_v_theta": eps_v,
                  "E_hess": vb["E_hess"], "E_4": vb["E_4"], "E_tail": vb["E_tail"],
                  "eps_bulk": vb["eps_bulk"]}
    epsT = T3.eps_Tprime(bs["mu"], bs["muJ"], T3.R_H, K["sens"])
    ell = T3.ell_loss(f, min(bs["mu"], 0.999) if np.isfinite(bs["mu"]) else 0.999)
    eps_a = ell / L + epsT + K["ca_traj"] / (T3.KAPPA * L)
    eps_delta = 1.0 - 2.0 * T3.KAPPA
    if (not np.isfinite(eps_v)) or eps_v >= 1.0 or eps_a >= 1.0 or not np.isfinite(eps_a):
        eps = float("inf")
    else:
        eps = ((1.0 + math.log(1.0 / (1.0 - eps_v)) / math.log(1.5))
               / ((1.0 - eps_a) * (1.0 - eps_delta)) - 1.0)
    out = {"L": L, "c": c, "column": column, "f": f, "gamma_off_exists": True,
           "eps_ell": ell / L, "eps_Tprime": epsT, "eps_a": eps_a, "eps_v": eps_v,
           "eps": eps, "d": vb["d"], "R_minus": vb["R_minus"], "N": vb["N"]}
    out.update(pieces)
    return out


def eps_at(L, col, safety):
    best, arg = float("inf"), None
    for f in FS:
        r = assemble_theta(L, CSTAR, col, f, safety)
        if r["eps"] < best:
            best, arg = r["eps"], f
    return best, arg


def Lstar(col, target, safety):
    lo, hi = 2.0, 1e12
    if eps_at(hi, col, safety)[0] > target:
        return None
    for _ in range(300):
        mid = math.sqrt(lo * hi)
        if eps_at(mid, col, safety)[0] <= target:
            hi = mid
        else:
            lo = mid
        if hi / lo < 1.0000001:
            break
    return hi


# --------------------------------------------------------------- control: reproduce THEOREM_S3
ctrl = {}
for col in ["proved", "proved_Ca_measured", "measured"]:
    ctrl[col] = {"L_star_eps<=0.5": T3.Lstar_for(col, 0.5), "L_star_eps<=0.1": T3.Lstar_for(col, 0.1)}
OUT["control_reproduce_THEOREM_S3"] = ctrl
OUT["FOREIGN_THEOREM_S3_L_star"] = {
    "proved_eps<=0.5": 311511.9, "proved_Ca_measured_eps<=0.5": 302704.3,
    "measured_eps<=0.5": 422.4, "proved_eps<=0.1": 1056242.9,
    "origin": "round2/THEOREM_S3/THEOREM_S3.md sec.4.5",
}
OUT["SHIFT"] = {"logReE_shift_theorem_datum": T3.SHIFT}

# --------------------------------------------------------------- vartheta admissibility
adm = {}
for f in FS:
    for safety in [1.0, 2.0, 4.0]:
        v = vartheta_of(f, safety)
        adm["f=%g,safety=%g" % (f, safety)] = {"vartheta": v, "admissible": (v is not None and v < 1.0)}
OUT["admissibility"] = adm

# --------------------------------------------------------------- the term table at L = 3e5, 1e6
tt = []
for L in [1e5, 3e5, 1e6]:
    for f in FS:
        r = assemble_theta(L, CSTAR, "proved", f, 2.0)
        tt.append({k: r.get(k) for k in ["L", "f", "vartheta", "admissible", "eps_bulk",
                                         "extra_bulk", "E_hess", "E_4", "E_tail",
                                         "eps_v_theta0", "eps_v", "eps_a", "eps"]})
OUT["term_table_proved_safety2"] = tt

# --------------------------------------------------------------- L_* with vartheta
res = {}
for col in ["proved", "proved_Ca_measured", "measured"]:
    for target in [0.5, 0.1]:
        for safety in [1.0, 2.0, 4.0]:
            k = "%s_eps<=%g_safety=%g" % (col, target, safety)
            Ls = Lstar(col, target, safety)
            e, bf = (eps_at(Ls, col, safety) if Ls else (float("inf"), None))
            res[k] = {"L_star": Ls, "best_f": bf,
                      "log_Lambda_star": (2 * Ls + T3.SHIFT) if Ls else None}
OUT["L_star_with_vartheta"] = res

# --------------------------------------------------------------- the comparison the brief wants
base = ctrl["proved"]["L_star_eps<=0.5"]
OUT["headline"] = {
    "L_star_theta0_proved_eps<=0.5": base,
    "L_star_theta_proved_eps<=0.5_safety2": res["proved_eps<=0.5_safety=2"]["L_star"],
    "ratio": (res["proved_eps<=0.5_safety=2"]["L_star"] / base) if base else None,
    "best_f_theta0": T3.eps_at(base, "proved")[1] if base else None,
    "best_f_theta": res["proved_eps<=0.5_safety=2"]["best_f"],
    "log_Lambda_star_theta0": (2 * base + T3.SHIFT) if base else None,
    "log_Lambda_star_theta": res["proved_eps<=0.5_safety=2"]["log_Lambda_star"],
}

# --------------------------------------------------------------- eps_v at the chosen f, both ways
ev = {}
for L in [1e5, 3e5, 1e6]:
    f0 = OUT["headline"]["best_f_theta"]
    r0 = assemble_theta(L, CSTAR, "proved", f0, 0.0, vartheta_override=0.0)
    r1 = assemble_theta(L, CSTAR, "proved", f0, 2.0)
    ev["L=%g" % L] = {"f": f0, "eps_v_theta0": r0["eps_v"], "eps_v_theta": r1["eps_v"],
                      "ratio": r1["eps_v"] / r0["eps_v"] if r0["eps_v"] else None,
                      "vartheta": r1.get("vartheta")}
OUT["eps_v_comparison_at_best_f"] = ev

print(json.dumps({k: OUT[k] for k in ["control_reproduce_THEOREM_S3", "headline",
                                      "eps_v_comparison_at_best_f"]}, indent=1, default=str))
print("\nadmissibility (safety = 2):")
for f in FS:
    a = adm["f=%g,safety=2" % f]
    print("  f=%-6g vartheta=%-12.4g  admissible=%s" % (f, a["vartheta"], a["admissible"]))
print("\nL_* table:")
for k, v in sorted(res.items()):
    print("  %-42s L_* = %-14s best_f = %s" % (k, v["L_star"], v["best_f"]))
json.dump(OUT, open(os.path.join(HERE, "q4_results.json"), "w"), indent=1, sort_keys=True,
          default=str)
print("\nwrote q4_results.json")
