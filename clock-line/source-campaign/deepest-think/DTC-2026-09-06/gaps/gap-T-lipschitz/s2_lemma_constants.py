"""s2 - the explicit constants of LEMMA T and its corollaries."""
import json, numpy as np
pi = np.pi
def Cfull(mu, muJ): return pi*(3*mu*(1+muJ)/(2*(1-mu)**5) + 3*muJ/8)
def Ceq(mus):       return pi*(3*(1+mus)/(2*(1-mus)**5) + 3.0/8)   # coefficient of mu when muJ = mu <= mus

out = {}
out['C0_small_mu']        = 15*pi/8                 # muJ = mu -> 0
out['C0_small_mu_value']  = float(15*pi/8)
out['C0_relative_form']   = float(15*pi/4)          # divided by a_T = (M/2) lambda L
out['gradK_const']        = 3/(2*pi**2)
out['K_const']            = 3/(8*pi**2)
out['pi3']                = pi**3
tab = {}
for mus in (0.0, 0.01, 0.05, 0.10, 0.20, 0.25, 0.50):
    tab[f'{mus:.2f}'] = dict(C_mu=float(Ceq(mus)), C_rel=float(2*Ceq(mus)))
out['C_table_muJ_eq_mu'] = tab
out['split'] = dict(term_I_small_mu=float(3*pi/2), term_II_small_mu=float(3*pi/8))
json.dump(out, open('s2_results.json','w'), indent=1)
print(f"C0 (muJ=mu->0)            = 15*pi/8 = {15*pi/8:.7f}")
print(f"C0 relative (/ (M/2)lamL) = 15*pi/4 = {15*pi/4:.7f}")
print(" mu*    C(mu*)    C_rel(mu*)")
for k,v in tab.items(): print(f" {k}   {v['C_mu']:9.4f}  {v['C_rel']:9.4f}")
print(f"term I coefficient 3pi/2 = {3*pi/2:.7f}, term II 3pi/8 = {3*pi/8:.7f}")
