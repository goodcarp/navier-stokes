#!/usr/bin/env python3
"""Closing arithmetic: relate the measured invariant Q' = T32*M*log(R/s*) to log Re_E.
Every constant here is computed from the recorded runs, none typed."""
import json, math, os, hashlib, statistics as st
H = os.path.dirname(os.path.abspath(__file__))
L = lambda t: json.load(open(os.path.join(H, "results_%s.json" % t)))["runs"]
main = sorted(L("main") + L("n7"), key=lambda r: r["N"])
rs = L("resweep")

lr = [r["ell"] / r["R"] for r in main]
print("ell/R over the primary sweep : %s   mean %.5f  spread %.2e"
      % (["%.4f" % x for x in lr], st.mean(lr[1:]), max(lr[1:]) - min(lr[1:])))
c_ell = st.mean(lr[1:])

vlim = [r for r in rs if r["Re0"] <= 16 and r["T32"]]
k = [r["s_star_at_T32"] / math.sqrt(r["nu"]) for r in vlim]
print("viscosity-limited branch (Re0 <= 16): s*/sqrt(nu/M) = %s  mean %.3f"
      % (["%.2f" % x for x in k], st.mean(k)))
kbar = st.mean(k)

# log(R/s*) = (1/2) log(M R^2/nu) - log kbar ;  M R^2/nu = Re_E / (ell/R)^2
# log(R/s*) = (1/2)log(M R^2/nu) - log kbar,  and  M R^2/nu = Re_E / (ell/R)^2
# so  log(R/s*) = (1/2) log Re_E - log(ell/R) - log kbar
off = math.log(kbar) + math.log(c_ell)
print("=> log(R/s*) = (1/2) log Re_E - %.4f     [offset = log(s*/sqrt(nu/M)) + log(ell/R)]" % off)
for r in vlim:
    pred = 0.5 * r["logReE"] - off
    print("   Re0=%-5g measured log(R/s*) = %.4f   predicted %.4f   diff %.4f"
          % (r["Re0"], math.log(r["R"] / r["s_star_at_T32"]), pred,
             math.log(r["R"] / r["s_star_at_T32"]) - pred))

Qp = [r["T32"] * math.log(r["R"] / r["s_star_at_T32"]) for r in main if r["T32"]]
print("\nPRIMARY invariant Q' = T32*M*log(R/s*) : mean %.4f  sd %.4f  range [%.4f, %.4f]"
      % (st.mean(Qp), st.pstdev(Qp), min(Qp), max(Qp)))
Qv = [r["T32"] * math.log(r["R"] / r["s_star_at_T32"]) for r in vlim]
print("VISCOSITY-LIMITED  Q' : %s  mean %.4f" % (["%.3f" % q for q in Qv], st.mean(Qv)))
print("\n=> in the viscosity-limited regime the measured law reads")
print("   T32 * M * [ (1/2) log Re_E - %.3f ] = %.3f   (mean over Re0 <= 16 at N = 5)"
      % (off, st.mean(Qv)))
print("   i.e.  T32 = %.3f / ( M ( log Re_E - %.3f ) )" % (2 * st.mean(Qv), 2 * off))
P = [r["T32"] * r["logReE"] for r in main if r["T32"]]
print("\nP = T32*M*log Re_E on the primary sweep: %s" % ["%.3f" % x for x in P])
print("   monotone decreasing: %s ; asymptote 2*Q' = %.3f" %
      (all(P[i] > P[i + 1] for i in range(len(P) - 1)), 2 * st.mean(Qp)))
print("SCRIPT-SHA256", hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest())
