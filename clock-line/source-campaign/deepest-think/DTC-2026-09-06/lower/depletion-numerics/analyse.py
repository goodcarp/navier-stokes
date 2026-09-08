#!/usr/bin/env python3
"""Analysis for lower/depletion-numerics. Applies the PREREG.md gates verbatim."""
import sys, os, json, math, glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SHARP = os.path.join(os.path.dirname(os.path.dirname(HERE)), "sharp", "viscous-numerics")
LEVELS = ["T32", "T2", "T4", "T8"]
CT = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0]


def load(fn):
    p = os.path.join(HERE, fn)
    return json.load(open(p))["runs"] if os.path.exists(p) else []


def allruns():
    out = []
    for fn in ("dep_main.json", "dep_n7.json", "dep_n7c.json", "dep_ctrl.json",
               "dep_conv.json", "dep_conv6.json", "dep_n2.json"):
        out += load(fn)
    # fall back to partials only for runs with no final file
    have = {(r["N"], round(r["h"], 5), r["lam"], r["sign"], r["label"]) for r in out}
    for p in sorted(glob.glob(os.path.join(HERE, "dep_*_partial.json"))):
        r = json.load(open(p))["runs"][0]
        k = (r["N"], round(r["h"], 5), r["lam"], r["sign"], r["label"])
        if k not in have:
            r["stop_reason"] = "PARTIAL-" + r["stop_reason"]
            out.append(r); have.add(k)
    return out


def Cat(r, t):
    h = r["hist"]
    if t > h["t"][-1]:
        return None
    return float(np.interp(t, h["t"], h["C"]))


def primary(runs):
    return {r["N"]: r for r in runs
            if r["label"] in ("main",) and abs(r["h"] - 0.125) < 1e-9 and r["sign"] > 0
            and r["lam"] == 1.5}


def main():
    runs = allruns()
    P = primary(runs)
    Ns = sorted(P)
    print("=== PRIMARY SWEEP  datum A, Re0=100, h=1/8, box 1.5R ===")
    print(f"{'N':>2} {'grid':>10} {'logReE':>7} {'T32':>8} {'T2':>8} {'T4':>8} {'T8':>8} "
          f"{'T4/T32':>7} {'C(0)':>6} {'stop':>22} {'t_end':>6} {'peak':>6}")
    for N in Ns:
        r = P[N]; c = r["cross"]
        g = lambda k: c[k]["t"] if k in c else None
        f = lambda v: ("%8.4f" % v) if v is not None else "      --"
        ratio = (g("T4") / g("T32")) if g("T4") and g("T32") else None
        print(f"{N:>2} {str(r['grid']):>10} {r['logReE']:7.3f} "
              f"{f(g('T32'))} {f(g('T2'))} {f(g('T4'))} {f(g('T8'))} "
              f"{(('%7.3f'%ratio) if ratio else '     --')} {r['hist']['C'][0]:6.3f} "
              f"{r['stop_reason']:>22} {r['t_final']:6.3f} {r['peak_ommax']:6.3f}")

    # box contamination flag (D5)
    print("\n--- crossing diagnostics: s*, s*/L (CONTAMINATED if s*/L > 0.5) ---")
    for N in Ns:
        r = P[N]
        for k in LEVELS:
            if k in r["cross"]:
                cc = r["cross"][k]
                flag = "  CONTAMINATED" if cc["sL"] > 0.5 else ""
                print(f"  N={N} {k}: t={cc['t']:.4f} s*={cc['s_star']:.3f} "
                      f"s*/L={cc['sL']:.3f} a*={cc['a_star']:.4f} C={cc['C']:.4f}{flag}")

    # C(t) table
    print("\n--- coherence count C(t) = max_x a / ||omega||_inf ---")
    print("  N " + "".join(f"{('t=%.2f'%t):>9}" for t in CT) + "   first t with C<1.5")
    firstC = {}
    for N in Ns:
        r = P[N]; h = r["hist"]
        row = "".join((f"{Cat(r,t):9.4f}" if Cat(r, t) is not None else "       --") for t in CT)
        Cv = np.array(h["C"]); tv = np.array(h["t"])
        idx = np.where(Cv < 1.5)[0]
        tf = float(tv[idx[0]]) if len(idx) else None
        firstC[N] = tf
        print(f" {N:>2} {row}   {('%.4f'%tf) if tf is not None else 'never'}")

    # a_max(t) itself, and where it and the vorticity maximum sit
    print("\n--- a_max(t) = max_x a  (the NUMERATOR of C), and the radii "
          "s_amax / s_ommax  (supplementary) ---")
    print("  N " + "".join(f"{('t=%.2f'%t):>9}" for t in CT)
          + "   | s of a_max at t=0,1 | s of omega_max at t=0, 0.5, 1, 1.5")
    for N in Ns:
        r = P[N]; h = r["hist"]
        row = "".join((f"{float(np.interp(t,h['t'],h['amax'])):9.4f}"
                       if t <= h["t"][-1] else "       --") for t in CT)
        sa = [float(np.interp(t, h["t"], h["s_amax"])) if t <= h["t"][-1] else float("nan")
              for t in (0.0, 1.0)]
        so = [float(np.interp(t, h["t"], h["s_star"])) if t <= h["t"][-1] else float("nan")
              for t in (0.0, 0.5, 1.0, 1.5)]
        print(f" {N:>2} {row}   | {sa[0]:6.3f} {sa[1]:6.3f}     | "
              + " ".join(f"{v:7.3f}" for v in so) + f"   (R={r['R']:g})")

    # marginal doubling times (supplementary diagnostic, not gated)
    print("\n--- marginal times per doubling (supplementary, NOT gated) ---")
    print("  N   T32      T2       T4-T2    T8-T4     (T4-T2)/T2   log(R/s*) at T2   at T4")
    for N in Ns:
        c = P[N]["cross"]
        t32 = c.get("T32", {}).get("t"); t2 = c.get("T2", {}).get("t")
        t4 = c.get("T4", {}).get("t"); t8 = c.get("T8", {}).get("t")
        d2 = (t4 - t2) if (t4 and t2) else None
        d3 = (t8 - t4) if (t8 and t4) else None
        lr2 = math.log(P[N]["R"] / c["T2"]["s_star"]) if "T2" in c else None
        lr4 = math.log(P[N]["R"] / c["T4"]["s_star"]) if "T4" in c else None
        fm = lambda v: ("%8.4f" % v) if v is not None else "      --"
        print(f" {N:>2} {fm(t32)} {fm(t2)} {fm(d2)} {fm(d3)}   "
              f"{fm(d2/t2 if d2 and t2 else None)}   {fm(lr2)}   {fm(lr4)}")

    # fits
    print("\n--- fitted exponents  log(T_k*M) = alpha - beta log N  (over N with a crossing,"
          " excluding CONTAMINATED) ---")
    fits = {}
    for k in LEVELS:
        xs, ys, used = [], [], []
        for N in Ns:
            c = P[N]["cross"]
            if k in c and c[k]["sL"] <= 0.5:
                xs.append(math.log(N)); ys.append(math.log(c[k]["t"])); used.append(N)
        if len(xs) >= 2:
            A = np.vstack([np.ones(len(xs)), -np.array(xs)]).T
            coef, *_ = np.linalg.lstsq(A, np.array(ys), rcond=None)
            fits[k] = (coef[0], coef[1], used)
            print(f"  {k:>4}: alpha={coef[0]:+.4f}  beta={coef[1]:+.4f}   N used = {used}")
        else:
            fits[k] = None
            print(f"  {k:>4}: fewer than 2 usable crossings ({used})")

    # 1/T_k linear in log Re_E ?  (the "T = c/(M(log Re_E - b))" form; supplementary)
    print("\n--- fit  1/(T_k M) = A_k log Re_E + B_k   =>  T_k = 1/(M(A_k log Re_E + B_k))"
          "  (supplementary, NOT gated) ---")
    for k in LEVELS:
        xs = [P[N]["logReE"] for N in Ns
              if k in P[N]["cross"] and P[N]["cross"][k]["sL"] <= 0.5]
        ys = [1.0 / P[N]["cross"][k]["t"] for N in Ns
              if k in P[N]["cross"] and P[N]["cross"][k]["sL"] <= 0.5]
        nn = [N for N in Ns if k in P[N]["cross"] and P[N]["cross"][k]["sL"] <= 0.5]
        if len(xs) >= 2:
            A = np.vstack([np.array(xs), np.ones(len(xs))]).T
            cf, *_ = np.linalg.lstsq(A, np.array(ys), rcond=None)
            res = np.array(ys) - A @ cf
            print(f"  {k:>4}: A={cf[0]:+.5f}  B={cf[1]:+.5f}  =>  1/A={1/cf[0]:.4f}, "
                  f"offset b=-B/A={-cf[1]/cf[0]:+.4f}   max|resid|/y={np.max(abs(res/np.array(ys))):.4f}"
                  f"   N={nn}")
        else:
            print(f"  {k:>4}: too few points {nn}")
    # the same fit read as a MEAN STRETCHING RATE:  ln(lambda_k)/T_k = alpha_k log Re_E + beta_k
    print("  read as a mean stretching rate  <a>_k := ln(lambda_k)/T_k "
          "= alpha_k (log Re_E - b_k):")
    lv = dict(T32=1.5, T2=2.0, T4=4.0, T8=8.0)
    for k in LEVELS:
        xs = [P[N]["logReE"] for N in Ns
              if k in P[N]["cross"] and P[N]["cross"][k]["sL"] <= 0.5]
        ys = [math.log(lv[k]) / P[N]["cross"][k]["t"] for N in Ns
              if k in P[N]["cross"] and P[N]["cross"][k]["sL"] <= 0.5]
        if len(xs) >= 2:
            A = np.vstack([np.array(xs), np.ones(len(xs))]).T
            cf, *_ = np.linalg.lstsq(A, np.array(ys), rcond=None)
            print(f"    {k:>4}: alpha={cf[0]:.5f}  b={-cf[1]/cf[0]:+.4f}   "
                  f"(depletion would need alpha to FALL sharply with the level)")

    # variation factors
    print("\n--- variation factors V_k = max_N T_k / min_N T_k (uncontaminated crossings) ---")
    V = {}
    for k in LEVELS:
        vals = [P[N]["cross"][k]["t"] for N in Ns
                if k in P[N]["cross"] and P[N]["cross"][k]["sL"] <= 0.5]
        nn = [N for N in Ns if k in P[N]["cross"] and P[N]["cross"][k]["sL"] <= 0.5]
        if len(vals) >= 2:
            V[k] = max(vals) / min(vals)
            print(f"  V_{k} = {V[k]:.4f}   over N={nn}  (values {[round(v,4) for v in vals]})")
        else:
            V[k] = None
            print(f"  V_{k} = -- (needs >=2, has {nn})")

    # the sharp seat's invariant, extended to every level (supplementary, NOT gated)
    print("\n--- Q_k = T_k * M * log(R/s*_k)  and mean log-growth rate  ln(k)/T_k "
          "(supplementary, NOT gated) ---")
    print("  N " + "".join(f"{('Q_'+k):>9}" for k in LEVELS)
          + "   |" + "".join(f"{('rate_'+k):>10}" for k in LEVELS))
    for N in Ns:
        c = P[N]["cross"]; R = P[N]["R"]
        lv = dict(T32=1.5, T2=2.0, T4=4.0, T8=8.0)
        qs, rs_ = [], []
        for k in LEVELS:
            if k in c:
                qs.append("%9.4f" % (c[k]["t"] * math.log(R / c[k]["s_star"])))
                rs_.append("%10.4f" % (math.log(lv[k]) / c[k]["t"]))
            else:
                qs.append("       --"); rs_.append("        --")
        print(f" {N:>2} " + "".join(qs) + "   |" + "".join(rs_))

    # frozen-initial-strain overhead (supplementary, NOT gated)
    print("\n--- overhead over the FROZEN-INITIAL-STRAIN clock  T_k^frozen = ln(k)/a_max(0),"
          "  a_max(0) = C(0)*M  (supplementary, NOT gated) ---")
    print("  N   a_max(0)   T32/T32^fr   T2/T2^fr   T4/T4^fr   T8/T8^fr")
    for N in Ns:
        c = P[N]["cross"]; a0 = P[N]["hist"]["C"][0] * P[N]["om0"]
        lv = dict(T32=1.5, T2=2.0, T4=4.0, T8=8.0)
        cells = []
        for k in LEVELS:
            if k in c:
                cells.append("%10.4f" % (c[k]["t"] / (math.log(lv[k]) / a0)))
            else:
                cells.append("        --")
        print(f" {N:>2} {a0:9.4f} " + " ".join(cells))
    # a_max(0) vs octave count
    if len(Ns) >= 2:
        x = np.array([math.log(P[N]["R"]) for N in Ns])
        y = np.array([P[N]["hist"]["C"][0] * P[N]["om0"] for N in Ns])
        A = np.vstack([x, np.ones(len(x))]).T
        cf, *_ = np.linalg.lstsq(A, y, rcond=None)
        print(f"  fit a_max(0) = {cf[0]:.4f} * log(R/rho0) + {cf[1]:+.4f}"
              f"   (first-order theory: 0.5 * log(R/rho0) + 0.2168)")

    # common-set variation factors (the like-for-like reading of the KILL clause)
    common = [N for N in Ns
              if all(k in P[N]["cross"] and P[N]["cross"][k]["sL"] <= 0.5 for k in ("T32", "T4"))]
    Vc = {}
    print(f"\n--- COMMON-SET variation factors, N = {common} (both T32 and T4 clean) ---")
    for k in ("T32", "T2", "T4"):
        vals = [P[N]["cross"][k]["t"] for N in common if k in P[N]["cross"]]
        Vc[k] = (max(vals) / min(vals)) if len(vals) >= 2 else None
        print(f"  Vc_{k} = {('%.4f' % Vc[k]) if Vc[k] else '--'}   values {[round(v,4) for v in vals]}")

    # PRE-REGISTERED GATE
    print("\n=== PRE-REGISTERED GATE (PREREG.md §4, verbatim) ===")
    reach4 = [N for N in Ns if "T4" in P[N]["cross"]]
    reach4_clean = [N for N in Ns if "T4" in P[N]["cross"] and P[N]["cross"]["T4"]["sL"] <= 0.5]
    S1 = (V["T4"] is not None) and V["T4"] < 1.5
    S2 = (V["T32"] is not None) and V["T32"] > 2
    S3 = all(firstC.get(N) is not None and firstC[N] <= 1.5 for N in Ns)
    print(f"  (S1) V_T4 < 1.5      : {V['T4']} -> {S1}")
    print(f"  (S2) V_T32 > 2       : {V['T32']} -> {S2}")
    print(f"  (S3) C<1.5 by t=1.5  : {firstC} -> {S3}")
    kill = None
    if V["T4"] and V["T32"]:
        k_own = abs(V["T4"] / V["T32"] - 1.0)
        print(f"  KILL, own-set    |V_T4/V_T32 - 1| = |{V['T4']:.4f}/{V['T32']:.4f} - 1| "
              f"= {k_own:.4f}  <= 0.30 -> {k_own <= 0.30}")
    if Vc.get("T4") and Vc.get("T32"):
        k_com = abs(Vc["T4"] / Vc["T32"] - 1.0)
        kill = k_com <= 0.30
        print(f"  KILL, common-set |V_T4/V_T32 - 1| = |{Vc['T4']:.4f}/{Vc['T32']:.4f} - 1| "
              f"= {k_com:.4f}  <= 0.30 -> {kill}   [the like-for-like reading; GATE USES THIS]")
        if Vc["T4"] > Vc["T32"]:
            print("    note: V_T4 EXCEEDS V_T32 -- T4 varies MORE with N than T32 does, "
                  "which is the opposite of depletion in the same direction as the KILL clause")
    missing = [N for N in Ns if N not in reach4]
    verdict = None
    if S1 and S2 and S3:
        verdict = "DEPLETION SURVIVES"
    elif kill:
        verdict = "DEPLETION KILLED"
    else:
        verdict = "INCONCLUSIVE"
    if missing:
        verdict += f"  (reach caveat: N={missing} never reached 4 M0)"
    if set(reach4) != set(reach4_clean):
        verdict += f"  (box caveat: T4 CONTAMINATED at N={sorted(set(reach4)-set(reach4_clean))})"
    print(f"\n  VERDICT: {verdict}")

    # controls
    print("\n=== CONTROLS ===")
    frozen = {}
    for fn in ("results_main.json", "results_n7.json"):
        p = os.path.join(SHARP, fn)
        if os.path.exists(p):
            for r in json.load(open(p))["runs"]:
                if r["label"] == "main" and abs(r["h"] - 0.125) < 1e-9:
                    frozen[r["N"]] = r
    print("  D1 continuity with the frozen sharp-seat run (must be < 1%):")
    for N in Ns:
        if N in frozen:
            for k, kk in (("T32", "T32"), ("T2", "T2")):
                a = P[N]["cross"].get(k, {}).get("t"); b = frozen[N][kk]
                if a and b:
                    print(f"    N={N} {k}: mine={a:.12f} frozen={b:.12f} "
                          f"rel={abs(a-b)/b:.3e}")
    print("  D2 maximum principle sup|eta| non-increasing (rel overshoot, must be <= 5e-3):")
    for r in runs:
        h = r["hist"]; e = np.array(h["etamax"])
        over = float(np.max(np.maximum.accumulate(e) - e[0]) / e[0])
        run_over = float(np.max(np.maximum(0.0, np.maximum.accumulate(e) - e)) / e[0])
        print(f"    {r['label']:>8} N={r['N']} h=1/{1/r['h']:.0f} lam={r['lam']:g} "
              f"sign={r['sign']:+.0f}: overshoot above eta0 = {over:.3e}, "
              f"final ratio = {r['etamax_ratio']:.4f}")
    ctl = [r for r in runs if r["label"] == "D3-sign"]
    if ctl:
        c0 = ctl[0]["hist"]
        print("  D3-derived TRUST HORIZON: ||omega||_inf of the SIGN-REVERSED datum at each"
              " primary crossing time (the anti-signal; 1.000 = no growth):")
        for N in Ns:
            for k in LEVELS:
                if k in P[N]["cross"]:
                    tt = P[N]["cross"][k]["t"]
                    v = float(np.interp(tt, c0["t"], c0["ommax"])) if tt <= c0["t"][-1] else float("nan")
                    tag = "  <-- SIGN-CONTROL DEGRADED" if v >= 1.2 else ""
                    print(f"    N={N} {k} at t={tt:.4f}: reversed-sign ||omega|| = {v:.4f}{tag}")
    for r in runs:
        if r["label"] == "D3-sign":
            print(f"  D3 sign control N={r['N']}: peak ommax = {r['peak_ommax']:.4f} "
                  f"(must stay < 1.5), crossings = {list(r['cross'])}, "
                  f"final om = {r['ommax_final']:.4f}, stop={r['stop_reason']}")
    print("  D4/D5 resolution and box:")
    base4 = P.get(4, {}).get("cross", {}).get("T4", {}).get("t")
    base6 = P.get(6, {}).get("cross", {}).get("T4", {}).get("t")
    base7 = P.get(7, {}).get("cross", {}).get("T4", {}).get("t")
    for r in runs:
        if r["label"] in ("D4-res", "D5-box"):
            t4 = r["cross"].get("T4", {}).get("t")
            t32 = r["cross"].get("T32", {}).get("t")
            base = {4: base4, 6: base6, 7: base7}.get(r["N"])
            rel = (abs(t4 - base) / base) if (t4 and base) else None
            print(f"    {r['label']} N={r['N']} h=1/{1/r['h']:.0f} lam={r['lam']:g}: "
                  f"T32={t32} T4={t4} base T4={base} rel={('%.4f'%rel) if rel else '--'} "
                  f"stop={r['stop_reason']}")

    json.dump(dict(V={k: V[k] for k in V}, Vc=Vc, common_N=common, fits={k: (fits[k][:2] if fits[k] else None) for k in fits},
                   firstC=firstC, verdict=verdict,
                   crossings={N: {k: P[N]["cross"][k] for k in P[N]["cross"]} for N in Ns},
                   Ct={N: {str(t): Cat(P[N], t) for t in CT} for N in Ns}),
              open(os.path.join(HERE, "gates.json"), "w"), indent=1, default=float)
    print("\nwrote gates.json")


if __name__ == "__main__":
    main()
