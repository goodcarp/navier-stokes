"""
u7 -- the angular sup.

mu := sup_{x in S} |Phi-Lambda|/|Lambda| is a supremum over the WHOLE shell, and the drive
constant C_a(phi) = 0.291999 + 0.014754 + 3.999218/sin(phi) + pi/8 of the far/near kernel
lemma diverges like 1/sin(phi) on the axis.  ASSEMBLY sec.2.4 evaluates it at phi = phi_0 = 30
degrees ("at the worst angle over the window") and so does the headline budget here.  That is
the sup over the TRAJECTORY, not over the shell.

This file prices the honest alternative: the majorant on the cone {delta <= phi <= pi-delta}
with C_a(delta), delta = 7.5 deg -- i.e. the whole shell except the taper cone, on which the
datum vanishes at the axis but the kernel bound does not.
"""
import json, math, os
import numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
src = open("u5_budget.py").read().split('if __name__ == "__main__":')[0]
G = {"__name__": "u7"}
exec(src, G)
assemble = G["assemble"]; COLUMNS = G["COLUMNS"]; CSTAR = G["CSTAR"]
C_a_proved = G["C_a_proved"]; RECORD = G["RECORD"]; LOGREE_SHIFT = G["LOGREE_SHIFT"]
U2 = G["U2"]
DEG = math.pi/180.0
FS = [0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0]

def C_R_at(phi, G_collar=0.0):
    return (3.0*C_a_proved(math.sin(phi)) + math.log(1.0/math.sin(phi))
            + 0.5 + RECORD_G_far + RECORD_G_inner + G_collar)
RECORD_G_far = U2["RECORD"]["G_far"]
RECORD_G_inner = U2["RECORD"]["G_inner"]

OUT = {}
for tag, phi in [("phi0=30deg", 30.0*DEG), ("phi=delta=7.5deg", 7.5*DEG),
                 ("phi=10deg", 10.0*DEG), ("phi=90deg", 90.0*DEG)]:
    OUT[tag] = {"C_a": C_a_proved(math.sin(phi)), "C_R": C_R_at(phi)}

# the drive is (C_R + 2 log lam_max); the near-field deficit term is lam_om C_a/(kappa L).
# Re-run the budget with both taken at phi = delta.
GC = C_R_at(7.5*DEG) - C_R_at(30.0*DEG)      # extra C_R, fed through the G_collar slot
def eps_at(L, col, gc):
    return min(assemble(L, CSTAR, col, f=f, G_collar=gc)["eps"] for f in FS)
def Lstar(col, target, gc):
    lo, hi = 2.0, 1e18
    if eps_at(hi, col, gc) > target:
        return None
    for _ in range(60):
        mid = math.sqrt(lo*hi)
        if eps_at(mid, col, gc) <= target:
            hi = mid
        else:
            lo = mid
        if hi/lo < 1.0000001:
            break
    return hi

OUT["extra_C_R_from_phi0_to_delta"] = GC
for col in ["proved", "measured"]:
    OUT["%s_L_star_cone_delta" % col] = Lstar(col, 0.5, GC)
    OUT["%s_L_star_phi0" % col] = Lstar(col, 0.5, 0.0)
OUT["logLambda_star_cone_delta_proved"] = (2*OUT["proved_L_star_cone_delta"] + LOGREE_SHIFT
                                           if OUT["proved_L_star_cone_delta"] else None)

with open("u7_results.json", "w") as fh:
    json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
print(json.dumps(OUT, indent=1, sort_keys=True, default=str))
