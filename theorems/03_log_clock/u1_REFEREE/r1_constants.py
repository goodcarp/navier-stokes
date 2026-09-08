"""Independent (from-scratch) re-derivation of the angular strain constants.
No import of any u*.py.  Derivation:
   a_ref = (3M/4) INT dlog rho INT_0^pi h(phi)|cos phi| sin^2 phi / g(phi;lam)^5 dphi
   g^2 = lam^2 sin^2 phi + lam^-4 cos^2 phi
   v = sin phi on [0,pi/2], doubled (h symmetric about pi/2)
   => a_ref = (M/2) P_h(lam),  P_h(lam) = 3 INT_0^1 h(arcsin v) v^2 (A v^2 + B)^{-5/2} dv
"""
import math
import mpmath as mp
mp.mp.dps = 30

DEG = mp.pi/180

def h_axis(phi, delta):
    pax = min(phi, mp.pi - phi)
    return min(mp.mpf(1), pax/delta)

def h_eq(phi, dm):
    return min(mp.mpf(1), abs(phi - mp.pi/2)/dm)

def P(lam, hfun):
    A = lam**2 - lam**(-4); B = lam**(-4)
    f = lambda v: hfun(mp.asin(v))*v**2*(A*v**2 + B)**mp.mpf(-2.5)
    return 3*mp.quad(f, [0, mp.mpf(1)])

d75 = mp.mpf('7.5')*DEG
dm5 = mp.mpf(5)*DEG
h75   = lambda p: h_axis(p, d75)
hfull = lambda p: h_axis(p, d75)*h_eq(p, dm5)
h1    = lambda p: mp.mpf(1)

print("P_1(1)                 =", mp.nstr(P(1, h1), 18))
print("P_1(1.5)/1.5           =", mp.nstr(P(mp.mpf('1.5'), h1)/mp.mpf('1.5'), 18), " (should be 1: P_1(lam)=lam)")
print("P_1(1.8377407156)/lam  =", mp.nstr(P(mp.mpf('1.8377407156083092'), h1)/mp.mpf('1.8377407156083092'), 18))
k75   = P(1, h75)/2
kfull = P(1, hfull)/2
print("kappa_delta(7.5, axis only) =", mp.nstr(k75, 18))
print("kappa_delta(7.5 & dm=5 FULL)=", mp.nstr(kfull, 18))
print("relative difference          =", mp.nstr((k75-kfull)/kfull, 8))
cs_axis = mp.log(mp.mpf(3)/2)/k75
cs_full = mp.log(mp.mpf(3)/2)/kfull
print("c_* (axis only)  =", mp.nstr(cs_axis, 12))
print("c_* (full datum) =", mp.nstr(cs_full, 12))
print("c_2 = 2log(3/2)/kappa  axis =", mp.nstr(2*mp.log(mp.mpf(3)/2)/k75, 12),
      " full =", mp.nstr(2*mp.log(mp.mpf(3)/2)/kfull, 12))

LM = mp.exp(3*cs_axis/4)
print("lam_max = exp(3c/4)  =", mp.nstr(LM, 18))

# r_h on [1, lam] by a fine scan + refinement
def r_h(hfun, lam_hi, n=241):
    best = None
    for i in range(n+1):
        lam = 1 + (lam_hi-1)*mp.mpf(i)/n
        val = P(lam, hfun)/lam
        if best is None or val < best[0]:
            best = (val, lam)
    return best

for nm, hf in [("h_axis 7.5", h75), ("h_full 7.5&5", hfull)]:
    v15, l15 = r_h(hf, mp.mpf('1.5'), 61)
    vlm, llm = r_h(hf, LM, 61)
    print("r_h[%s] on [1,1.5]      = %s at lam=%s" % (nm, mp.nstr(v15,12), mp.nstr(l15,6)))
    print("r_h[%s] on [1,lam_max]  = %s at lam=%s" % (nm, mp.nstr(vlm,12), mp.nstr(llm,6)))

# monotonicity of P_h on [1, 2.2]
prev = None; minc = None; turn = None
for i in range(0, 241):
    lam = 1 + mp.mpf('0.005')*i
    val = P(lam, h75)
    if prev is not None:
        inc = val - prev
        if minc is None or inc < minc: minc = inc
        if inc < 0 and turn is None: turn = lam
    prev = val
print("P_h75 min increment over [1,2.2] step .005 =", mp.nstr(minc, 8), " first decrease at lam =", turn)
