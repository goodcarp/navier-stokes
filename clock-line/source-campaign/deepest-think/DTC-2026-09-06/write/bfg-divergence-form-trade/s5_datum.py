#!/usr/bin/env python3
"""
s5_datum.py -- what each clock says for the campaign's own datum.

For the bang-bang shell (the delta -> 0 endpoint of the mollified delta-tapered
plateau) the estate has, computed exactly in sharp/exact-first-order (c6_shell_exact):

    E_shell = 0.172403978 M^2 R^5 ,   log Re_E = 2 L - 0.7031660 ,  L = log(R/rho_0),

and its frozen-frame doubling window is  M t_d = ln2 / (a_inner/M) = 2 ln 2 / L
(a_inner = (M/2) L).  Those two numbers are QUOTED from that seat (with its file name);
everything else on this page is recomputed here.

The delta-taper shifts the Reynolds number by 2 log(1/delta) (quoted from
gaps/refute-gap-T-lipschitz: log Re_E = 2L + 2log(1/delta) - 0.7031659); using the
shell value is the choice most favourable to the power clock, so the conclusion below
is a lower bound on how bad the power clock is.
"""
import json
import mpmath as mp
mp.mp.dps = 40

s4 = json.load(open('s4_results.json'))
c_div = mp.mpf(repr(s4['c_div']))
out = {'c_div': float(c_div), 'quoted': {
    'log_Re_E_offset': -0.7031660,
    'source': 'sharp/exact-first-order/NOTE.md sec.3 (c6_shell_exact.py)'}}

print("QUOTED from sharp/exact-first-order (not re-derived here):")
print("   log Re_E = 2L - 0.7031660        E_shell = 0.172403978 M^2 R^5")
print("   frozen-frame doubling window     M t_d = 2 ln2 / L")
print()
hdr = "%-9s %-14s %-16s %-16s %-16s %-12s" % ("L", "log Re_E", "Re_E", "(b) c_div/Re_E", "(a) 1/(1+logRe)", "M t_d model")
print(hdr); print("-"*len(hdr))
rows = []
for L in [mp.mpf('8.3178'), mp.mpf(20), mp.mpf(50), mp.mpf(100)]:
    lR = 2*L - mp.mpf('0.7031660')
    ReE = mp.e**lR
    b = c_div/ReE
    a = 1/(1+lR)                      # this is (a)/c_1
    md = 2*mp.log(2)/L
    rows.append([float(L), float(lR), float(ReE), float(b), float(a), float(md), float(a/b)])
    print("%-9.4f %-14.6f %-16.6e %-16.6e %-16.6e %-12.6f" % tuple(rows[-1][:6]))
out['rows'] = rows
print()
print("ratio (a)/(b) = c_1 * Re_E/(1+log Re_E):")
for r in rows:
    print("   L = %-8.4f   (a)/(b) = %.6e * c_1" % (r[0], r[6]))
print()
print("ratio of the model window to the proved power window,  (M t_d)/(b):")
for r in rows:
    print("   L = %-8.4f   %.6e" % (r[0], r[5]/r[3]))
out['model_over_power'] = [r[5]/r[3] for r in rows]
print()
print("ratio of the model window to the proved LOG window,  (M t_d)/((a)*c_1):")
for r in rows:
    print("   L = %-8.4f   %.6f / c_1" % (r[0], r[5]/r[4]))
out['model_over_log'] = [r[5]/r[4] for r in rows]

json.dump(out, open('s5_results.json','w'), indent=1)
print("\ns5: done")
