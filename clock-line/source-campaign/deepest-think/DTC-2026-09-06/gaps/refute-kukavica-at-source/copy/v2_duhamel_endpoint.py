"""
V3 of v1_verify_numbers.py integrates s over (0,t) directly; the grad-kernel Duhamel integrand
has a 1/sqrt(t-s) endpoint singularity there, so plain Gauss-Legendre in s stalls at ~4 digits
(8.3715 vs 8pi/3 = 8.37758).  This script removes the singularity with s = t - sigma^2 and
also reports the closed forms, so the sqrt(t)-vs-t dichotomy is pinned exactly.

Closed forms used (verified numerically below):
  int_{R^3} (|z|+a)^{-4} dz = 4 pi int_0^inf r^2 (r+a)^{-4} dr = 4 pi / (3a)
  => int_0^t int (|z|+sqrt(t-s))^{-4} dz ds        = (8 pi / 3) sqrt(t)
  => int_0^t int sqrt(t-s)(|z|+sqrt(t-s))^{-4} dz ds = (4 pi / 3) t
"""
import json, math
import numpy as np

x, w = np.polynomial.legendre.leggauss(600)
_C = {}
def inner(a, n=1200):
    if n not in _C:
        xx, ww = np.polynomial.legendre.leggauss(n)
        s = 0.5*(xx+1.0); ww = 0.5*ww
        _C[n] = (s/(1.0-s), ww/(1.0-s)**2)
    r, wr = _C[n]
    return 4.0*math.pi*np.sum((wr*r*r)[None, :]/(r[None, :]+np.atleast_1d(a)[:, None])**4, axis=1)

res = {}
# check the radial closed form 4 pi /(3a)
res['radial_closed_form_check'] = [
    dict(a=a, numeric=float(inner(a)[0]), closed=4*math.pi/(3*a)) for a in (0.1, 1.0, 3.0)
]
rows = []
for t in (0.25, 1.0, 4.0):
    sig = 0.5*math.sqrt(t)*(x+1.0); ws = 0.5*math.sqrt(t)*w   # s = t - sigma^2, sigma in (0, sqrt t)
    inn = inner(sig)
    A = float(np.sum(ws*2.0*sig*inn))          # grad-kernel Duhamel  (2 sigma dsigma = -ds)
    B = float(np.sum(ws*2.0*sig*sig*inn))      # G-kernel Duhamel     (extra sqrt(t-s) = sigma)
    rows.append(dict(t=t,
                     gradG_integral=A, gradG_over_sqrt_t=A/math.sqrt(t),
                     G_integral=B, G_over_t=B/t,
                     closed_8pi_over_3=8*math.pi/3, closed_4pi_over_3=4*math.pi/3))
res['duhamel_desingularised'] = rows
print(json.dumps(res, indent=1))
