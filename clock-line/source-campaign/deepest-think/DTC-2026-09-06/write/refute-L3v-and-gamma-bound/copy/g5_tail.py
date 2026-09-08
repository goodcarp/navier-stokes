"""g5 — the analytic tail of the L3v sum.

Term_l = |H_l| ||C_l^{3/2}||_inf / (l(l+3)-4).  Write H_l = H_l^J + R_l with the JUMP part
H_l^J = (2 lam M/N_l) P_{l+1}(0)  (g2 measures |R_l| decaying two orders faster).  Then
Term_l^J = lam M (l+3/2) |P_{l+1}(0)| / ((l+4)(l-1)),  and with l+1 = 2m,
|P_{2m}(0)| = binom(2m,m)/4^m, which obeys the Wallis bounds
        1/sqrt(pi(m+1/2))  <=  binom(2m,m)/4^m  <=  1/sqrt(pi m)          (checked here).
The tail is summed exactly to l = 2*LEND_m-1 by the two-term recurrence and bounded beyond it.
"""
import json
import numpy as np

out = {}

# ---- Wallis bounds, checked ---------------------------------------------------------
LEND_m = 2000000
m = np.arange(1, LEND_m + 1)
logc = np.cumsum(np.log((2 * m - 1) / (2 * m)))
c = np.exp(logc)                      # c[i] = binom(2m,m)/4^m for m = i+1
lo = 1 / np.sqrt(np.pi * (m + 0.5))
hi = 1 / np.sqrt(np.pi * m)
out['wallis_lower_ok'] = bool(np.all(c >= lo * (1 - 1e-13)))
out['wallis_upper_ok'] = bool(np.all(c <= hi * (1 + 1e-13)))
out['wallis_worst_lower_ratio'] = float(np.min(c / lo))
out['wallis_worst_upper_ratio'] = float(np.max(c / hi))

L0 = 4001

def jump_tail(lam, L0=L0):
    """sum over odd l > L0 of lam*(l+3/2)|P_{l+1}(0)|/((l+4)(l-1)), exact to l = 2*LEND_m-1,
       plus an upper bound on the remainder."""
    m0 = (L0 + 1) // 2                    # l = L0  <->  m = m0 = 2001
    mm = m[m0:]                           # m = m0+1, ...   i.e. l = 2m-1 > L0
    cc = c[m0:]
    l = 2.0 * mm - 1.0
    s = float(lam * np.sum((l + 1.5) * cc / ((l + 4.0) * (l - 1.0))))
    mr = float(LEND_m)
    # term(m) = lam (2m+0.5) c_m /((2m+3)(2m-2)) <= lam (1/(2m-2)) (pi m)^{-1/2}
    # sum_{m>mr} <= lam pi^{-1/2} int_mr^inf dx/((2x-2) sqrt(x)) = lam pi^{-1/2} * (1/ sqrt(mr-1)) * atanh-ish
    rem = float(lam / np.sqrt(np.pi) * (1.0 / np.sqrt(mr - 1.0)) * 1.05)
    return s, rem

res = {}
d = json.load(open('g2_results.json'))
for key, v in d['cases'].items():
    lam = float(key.split('_')[0][3:])
    jt, rem = jump_tail(lam)
    p4 = v['partial_sums']['4001']
    rt = v['resid_tail_estimate']
    res[key] = {'partial_4001': p4, 'jump_tail': jt, 'jump_tail_rem_bound': rem,
                'resid_tail_est': rt, 'S': p4 + jt + rt,
                'S_lo': p4 + jt - rem - abs(rt), 'S_hi': p4 + jt + rem + 3 * abs(rt),
                'kappa': v['kappa'], 'resid_p': v['resid_decay_exponent']}
    print('%-18s S = %.4f  [%.4f, %.4f]   (part %.5f + tail %.5f + resid %.6f)'
          % (key, res[key]['S'], res[key]['S_lo'], res[key]['S_hi'], p4, jt, rt))
out['L3v_sums'] = res
out['sup_over_lambda'] = {}
for dd in ['0.0', '7.5', '15.0', '30.0']:
    vals = [res[k]['S'] for k in res if k.endswith('del' + dd)]
    out['sup_over_lambda'][dd] = max(vals)
print(json.dumps(out['sup_over_lambda'], indent=1))
print('wallis', out['wallis_lower_ok'], out['wallis_upper_ok'],
      out['wallis_worst_lower_ratio'], out['wallis_worst_upper_ratio'])
with open('g5_results.json', 'w') as f:
    json.dump(out, f, indent=1)
