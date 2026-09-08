"""
r2b -- same machinery as r2, at hk2 PROOF.md's corrected window-sup constants
        K2 <= 161.7735 M/rho0 (proved, lambda = 3/2 end of the window) and 5.3854 (measured window sup);
        the assembly's a8 used 66.6622 / 3.0202, which hk2 section 14 item 7 withdrew.

r2 -- is the assembly's L_* the threshold of the theorem it states?

ASSEMBLY.md section 1.3 states:  for every L >= L_*,  T_d <= 2 log(3/2) (1+eps) / (M L),
proved by contradiction on the window  tau = c/(M L),  c = 2 log(3/2)(1+eps).
The bootstrap (2.1)-(2.2) must therefore be run on the window of length c, and the
closing condition is

        c  >=  c_*(1 + eps(L, c))        with  c_* = log(3/2)/kappa_delta,

i.e. a SELF-CONSISTENT window.  a2_budget.py computes eps(L, c_*) only and defines
L_* by eps(L, c_*) <= 1/2 (or 0.1).  Because eps(L, c) grows like e^{p c} with
p c_* ~ 34, eps(L, c_*) <= 1/2 does not imply that any window closes.

This script re-uses the assembly's OWN machinery (a copy of a2_budget.py, executed
as a4_subwindows.py does) and computes, per column,

   F(L) = min over c in [c_*, c_max] and f in F_SET of  [ c - c_*(1 + eps(L, c)) ]

and bisects for the smallest L with F(L) >= 0.  That is the L_* of the theorem as
stated.  It also evaluates eps at the assembly's own L_* on the window c = 1.5 c_*
to show directly that the window the headline needs does not close there.

Everything below is computed here; nothing is typed from the record.
"""
import json, math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
COPIES = os.path.join(HERE, "copies")
os.chdir(COPIES)
# a2 reads a1_results.json from cwd; the re-run of a1 in copies/ writes it
while not os.path.exists(os.path.join(COPIES, "a1_results.json")):
    time.sleep(5)
time.sleep(2)
src = open("a2_budget.py").read().split('if __name__ == "__main__":')[0]
G = {"__name__": "r2"}
exec(src, G)
assemble = G["assemble"]; bootstrap = G["bootstrap"]; COLUMNS = G["COLUMNS"]
RECORD = G["RECORD"]; KAPPA = G["KAPPA"]; C_a_proved = G["C_a_proved"]; PHI0 = G["PHI0"]
ell_loss = G["ell_loss"]
CSTAR = math.log(1.5)/KAPPA
EPS_DELTA = 1.0 - 2.0*KAPPA

OUT = {"c_star": CSTAR, "eps_delta": EPS_DELTA}

def eps_N(L, N, column, c, f, C_K=1.0):
    """eps on a window of length c (in theta units) with N sub-windows for the map error
    (N = 1 is the assembly's a2; N > 1 is its a4 restart pricing), viscous budget on the
    whole window.  Same formulas as a2.assemble / a4.eps_N, with c free."""
    cp_key, ca_kind, use_V1, sens = COLUMNS[column]
    Cprime = RECORD[cp_key]
    C_a = C_a_proved(math.sin(PHI0)) if ca_kind == "far_near" else RECORD["material_offset_M"]
    bs = bootstrap(L, c/N, Cprime, C_a, sens=sens)
    if not bs["closed"]:
        return float("inf")
    vb = G["viscous_budget"](L, c, C_K, f, c_G=(1.5 + Cprime/L)*c, Cprime=Cprime,
                             use_V1_for_Z=use_V1)
    eps_a = ell_loss(f)/L + bs["eps_Tprime"] + 1.5*C_a/(KAPPA*L)
    eps_v = vb["eps_v"]
    if (not np.isfinite(eps_v)) or eps_v >= 1.0 or eps_a >= 1.0 or not np.isfinite(eps_a):
        return float("inf")
    return (1.0 + math.log(1.0/(1.0-eps_v))/math.log(1.5))/((1.0-eps_a)*(1.0-EPS_DELTA)) - 1.0

# consistency with the assembly's own assemble() at c = c_*
chk = {}
for col in ["proved", "measured"]:
    for L in [160.0, 1e15, 1e18]:
        a = assemble(L, CSTAR, col, f=1.0)["eps"]
        b = eps_N(L, 1, col, CSTAR, 1.0)
        chk["%s_L=%g" % (col, L)] = {"assemble": a, "eps_N": b}
OUT["consistency_with_a2_assemble"] = chk

F_SET = [0.25, 1.0, 4.0]
def slack(L, N, column, C_K=1.0, ratios=None):
    """max over c of  c - c_*(1+eps(L,c)) ; >= 0 means a self-consistent window exists.
    Returns (best_slack, best_c_over_cstar, best_eps)."""
    if ratios is None:
        ratios = np.concatenate([np.linspace(1.0, 1.1, 21), np.linspace(1.1, 1.5, 17)[1:],
                                 np.linspace(1.5, 3.0, 7)[1:]])
    best = (-float("inf"), None, None)
    for r in ratios:
        c = CSTAR*r
        e = min(eps_N(L, N, column, c, f, C_K) for f in F_SET)
        if not np.isfinite(e):
            continue
        s = c - CSTAR*(1.0+e)
        if s > best[0]:
            best = (s, r, e)
    return best

def Lstar_selfconsistent(N, column, C_K=1.0, lo=1.0, hi=1e40):
    if slack(hi, N, column, C_K)[0] < 0:
        return None, None
    for _ in range(80):
        mid = math.sqrt(lo*hi)
        s = slack(mid, N, column, C_K)
        if s[0] >= 0:
            hi = mid; last = s
        else:
            lo = mid
        if hi/lo < 1.001:
            break
    return hi, last


sc = {}
for col, N, C_K in [("proved", 16, 161.7735), ("measured", 16, 161.7735), ("measured", 16, 5.3854), ("proved", 1, 161.7735)]:
    t0 = time.time()
    Ls, info = Lstar_selfconsistent(N, col, C_K)
    key = "%s_N=%d_CK=%g" % (col, N, C_K)
    sc[key] = {"L_star_selfconsistent": Ls, "c_over_cstar_at_threshold": info[1] if info else None,
               "eps_at_threshold": info[2] if info else None, "seconds": time.time()-t0}
    print(key, sc[key], flush=True)
# E_hess ratio at fixed (L, c, f): linear in C_K
vbf = G["viscous_budget"]
sc["E_hess_ratio_161.7735_over_66.6622"] = vbf(640.0, CSTAR, 161.7735, 1.0, Cprime=119.33)["E_hess"]/vbf(640.0, CSTAR, 66.6622, 1.0, Cprime=119.33)["E_hess"]
OUT["L_star_selfconsistent_hk2_corrected"] = sc
with open(os.path.join(HERE, "r2b_results.json"), "w") as fh:
    json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
print(json.dumps(sc, indent=1, default=str))
