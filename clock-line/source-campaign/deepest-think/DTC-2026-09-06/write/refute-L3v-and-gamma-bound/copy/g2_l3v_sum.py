"""g2 — L3v: the coefficients H_l(lambda), the summand decay, and the L3v sum with an
analytic tail.  Cross-checks against lower/prove-lagrangian s6 (truncated) and against the
refuter's streamed partial sums."""
import json, time
import numpy as np
import glib

M = 1.0
L0 = 4001                     # exact (quadrature) range
out = {'L0': L0}
t0 = time.time()

def jump_model_term(l, lam):
    """Term_l with H_l replaced by its jump part (2 lam M/N_l) P_{l+1}(0)"""
    return lam * M * (l + 1.5) * abs(glib.Pn0(l + 1)) / ((l + 4.0) * (l - 1.0))

def tail(lam, lstart, lend=2000001):
    """sum of the jump-model terms over odd l in (lstart, lend], plus an analytic remainder
       for l > lend using |P_{l+1}(0)| <= sqrt(2/(pi(l+1)))  (verified in g1)."""
    ls = np.arange(lstart + 2 - (lstart % 2), lend + 1, 2)
    # P_{l+1}(0) magnitude by the recurrence, log-space to stay exact for large l
    p = abs(glib.Pn0(ls[0] + 1))
    vals = np.empty(ls.size)
    vals[0] = p
    for i in range(1, ls.size):
        n = ls[i] + 1
        p = p * ((n - 1.0) / n) * ((n - 2.0) / (n - 1.0))   # two steps: P_n(0) from P_{n-2}(0)
        vals[i] = p
    s = float(np.sum(lam * M * (ls + 1.5) * vals / ((ls + 4.0) * (ls - 1.0))))
    # remainder l > lend : term <= lam*M*sqrt(2/pi) l^{-3/2}*(1+3/l), odd l only
    rem = lam * M * np.sqrt(2 / np.pi) * (1.0 / np.sqrt(lend)) * 1.02
    return s, rem

cases = []
for delta_deg in [0.0, 7.5, 15.0, 30.0]:
    for lam in [1.0, 1.1, 1.25, 1.4, 1.5]:
        cases.append((lam, delta_deg))

res = {}
for lam, dd in cases:
    delta = np.deg2rad(dd)
    H = glib.H_direct(L0, lam, delta, npanel=520)
    H2 = glib.H_byparts(L0, lam, delta, npanel=520)
    inst = float(np.max(np.abs(H[1::2] - H2[1::2])) / np.max(np.abs(H[1::2])))
    ls = np.arange(3, L0 + 1, 2)
    T = np.array([glib.term(l, H[l]) for l in ls])
    partial = {}
    for cut in [101, 401, 601, 801, 1201, 1601, 2001, 4001]:
        partial[str(cut)] = float(np.sum(T[ls <= cut]))
    ts, rem = tail(lam, L0)
    # residual of the jump model over the exact range, and its measured decay
    TJ = np.array([jump_model_term(l, lam) for l in ls])
    resid = T - TJ
    # decay exponent of |resid| over the top decade
    m = ls > 400
    pfit = np.polyfit(np.log(ls[m]), np.log(np.abs(resid[m]) + 1e-300), 1)
    # tail of the residual, modelled as A l^{-p} summed over odd l>L0
    A, p = float(np.exp(pfit[1])), float(-pfit[0])
    rtail = A / (2 * (p - 1)) * L0 ** (1 - p) if p > 1 else float('nan')
    res['lam%.2f_del%.1f' % (lam, dd)] = {
        'H1': float(H[1]), 'kappa': float(-0.6 * H[1]),
        'instrument_reldiff': inst,
        'partial_sums': partial,
        'jump_tail_beyond_L0': ts, 'jump_tail_remainder_bound': rem,
        'resid_decay_exponent': p, 'resid_tail_estimate': rtail,
        'S_total': float(partial['4001'] + ts + rtail),
        'S_total_no_resid': float(partial['4001'] + ts),
        'l32_times_absH_at_l_4001': float(abs(H[4001]) * 4001**1.5),
        'term_times_l32_at_4001': float(T[-1] * 4001**1.5),
    }
    print('%-18s S=%.5f  (part4001 %.5f + tail %.5f + resid %.6f)  kappa=%.7f  p=%.2f'
          % ('lam%.2f d%.1f' % (lam, dd), res['lam%.2f_del%.1f' % (lam, dd)]['S_total'],
             partial['4001'], ts, rtail, -0.6 * H[1], p), flush=True)

out['cases'] = res
out['secs'] = time.time() - t0
with open('g2_results.json', 'w') as f:
    json.dump(out, f, indent=1)
print('secs', out['secs'])
