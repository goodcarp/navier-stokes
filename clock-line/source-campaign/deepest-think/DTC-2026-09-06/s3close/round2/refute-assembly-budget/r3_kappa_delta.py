"""
r3 -- kappa_delta from a closed form, to adjudicate the assembly's correction of
      prove-lagrangian's kappa_delta(5 deg) = 0.4999917  ->  0.4999171 (section 1.2) /
      0.4999169 (section 0.1).

kappa_delta = (1/2) P_h(1),  P_h(1) = 3 int_0^1 h(arcsin v) v^2 dv  (A = 0, B = 1 at lambda = 1),
h = min(1, phi/delta).  Split at v = sin(delta):

   kappa_delta = 1/2 - (3/2) int_0^{sin delta} (1 - arcsin(v)/delta) v^2 dv .

Evaluated in 40-digit mpmath, and the small-delta law 1 - 2 kappa = delta^3/4 + O(delta^5)
derived by series (sympy).  Independent of a1's scipy quad on the (Av^2+B)^{-5/2} form.
"""
import json, math, os
import mpmath as mp
import sympy as sp

mp.mp.dps = 40
OUT = {}
def kappa(delta_rad):
    d = mp.mpf(delta_rad)
    w = mp.sin(d)
    I = mp.quad(lambda v: (1 - mp.asin(v)/d)*v*v, [0, w])
    return mp.mpf(1)/2 - mp.mpf(3)/2*I

deg = mp.pi/180
tab = {}
for dd in [3, 5, 7.5, 10, 15, 30]:
    k = kappa(dd*deg)
    tab["%g" % dd] = {"kappa": mp.nstr(k, 16), "1-2kappa": mp.nstr(1-2*k, 10),
                      "delta^3/4": mp.nstr((dd*deg)**3/4, 10)}
OUT["kappa_table"] = tab
OUT["record_prove_lagrangian_quoted"] = {"5": 0.4999917, "7.5": 0.4997212, "15": 0.4978077, "30": 0.4836252}
OUT["assembly_quoted"] = {"5_sec1.2": 0.4999171, "5_sec0.1": 0.4999169, "7.5": 0.4997212305210886}
k5 = kappa(5*deg)
OUT["kappa5_16digits"] = mp.nstr(k5, 16)
OUT["kappa5_rounded7"] = float(mp.nstr(k5, 7))
OUT["one_minus_2kappa5"] = float(1-2*k5)
OUT["record_value_implies_1-2kappa"] = 1-2*0.4999917
# series law
d, v = sp.symbols('delta v', positive=True)
I = sp.integrate((1 - sp.asin(v)/d)*v**2, (v, 0, sp.sin(d)))
ser = sp.series(1 - 2*(sp.Rational(1, 2) - sp.Rational(3, 2)*I), d, 0, 6).removeO()
OUT["series_1_minus_2kappa"] = str(sp.simplify(ser))
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "r3_results.json"), "w") as fh:
    json.dump(OUT, fh, indent=1, sort_keys=True)
print(json.dumps(OUT, indent=1, sort_keys=True))
