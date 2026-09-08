#!/usr/bin/env python3
"""refute-u2 / r2 -- independent recomputation of the datum constants and the kappa table.
Written from the definitions in u2/PROOF.md sec.1 and ASSEMBLY sec.1.2-1.3; no u2 code imported."""
import math, numpy as np, mpmath as mp
SD = math.sin(math.radians(7.5)); W = 0.20; EPS = 0.25
def Th_rTh(lr, L):
    a1=(lr-EPS)/EPS; a2=(lr-L+EPS)/EPS
    return 0.5*(np.tanh(a1)-np.tanh(a2)), 0.5*(np.cosh(a1)**-2-np.cosh(a2)**-2)/EPS
def hh(p):  s=np.sin(p); c=np.cos(p); return np.tanh(s/SD)*np.tanh(c/W)/s
def hhp(p):
    s=np.sin(p); c=np.cos(p); T=np.tanh(s/SD); Z=np.tanh(c/W)
    return ((np.cosh(s/SD)**-2)*(c/SD)*Z + T*(-(np.cosh(c/W)**-2)*(s/W)))/s - T*Z*c/s**2
ph = np.linspace(1e-9, math.pi-1e-9, 4000001); h = hh(ph); hp = hhp(ph)
print("E0 = sup|h| =", np.abs(h).max(), " tanh(1/w)/sin d =", math.tanh(1/W)/SD)
print("G0 = bulk sup sqrt(h^2+h'^2) =", np.hypot(h,hp).max(), " at phi =",
      math.degrees(ph[int(np.argmax(np.hypot(h,hp)))]), "deg   [u2 grid value 20.919707]")
for L in (10.0, 40.0):
    lr = np.linspace(-6, L+6, 4000001); Th,_ = Th_rTh(lr, L)
    k = int(np.argmax(Th*np.exp(-lr)))
    print(f"L={L}: ||eta_0||_inf = {(Th*np.exp(-lr))[k]*np.abs(h).max():.6f} M/rho0 at rho={math.exp(lr[k]):.5f}"
          f"   [claimed M/(rho0 sin delta) = {1/SD:.6f}]")
# ASSEMBLY 1.1 angular family
d = math.radians(7.5); dm = math.radians(5.0)
p2 = np.linspace(1e-9, math.pi/2, 4000001)
hA = np.minimum(1, p2/d)*np.minimum(1, np.abs(p2-math.pi/2)/dm)/np.sin(p2)
hAp = np.gradient(hA, p2)
print("ASSEMBLY 1.1 angular family: G0 =", np.hypot(hA,hAp)[10:-10].max(),
      " at phi =", math.degrees(p2[10+int(np.argmax(np.hypot(hA,hAp)[10:-10]))]), "deg; E0 =", hA.max())
# kappa table
mp.mp.dps = 30
D = mp.pi*mp.mpf('7.5')/180; DM = mp.pi*mp.mpf(5)/180
SDm = mp.sin(D); Wm = mp.mpf('0.20')
def kap(f, pts): return mp.mpf(3)/2*mp.quad(lambda p: f(p)*mp.sin(p)**2*mp.cos(p), pts)
T = {'bang-bang': (lambda p: mp.mpf(1), [0, mp.pi/2]),
     'sharp taper 7.5': (lambda p: min(mp.mpf(1), p/D), [0, D, mp.pi/2]),
     'tanh taper 7.5': (lambda p: mp.tanh(mp.sin(p)/SDm), [0, D, mp.pi/2]),
     'linear equat 5': (lambda p: min(mp.mpf(1), (mp.pi/2-p)/DM), [0, mp.pi/2-DM, mp.pi/2]),
     'tanh equat 0.20': (lambda p: mp.tanh(mp.cos(p)/Wm), [0, mp.pi/2]),
     'ASSEMBLY 1.1': (lambda p: min(mp.mpf(1), p/D)*min(mp.mpf(1), (mp.pi/2-p)/DM), [0, D, mp.pi/2-DM, mp.pi/2]),
     'campaign (D-C)': (lambda p: mp.tanh(mp.sin(p)/SDm)*mp.tanh(mp.cos(p)/Wm), [0, D, mp.pi/2])}
for n,(f,pts) in T.items():
    k = kap(f, pts)
    print(f"{n:18s} kappa={mp.nstr(k,12):16s} 1-2k={mp.nstr(1-2*k,8):14s} floor={mp.nstr((1-2*k)/(2*k),8)}")
# P1/P2 sharpness on the reference strain for (D-C)
for lam in (1.25, 1.5):
    s = np.sin(ph); c = np.cos(ph)
    scale = np.sqrt(lam**2*s**2 + lam**-4*c**2)
    p1 = (np.abs(h)*scale).max()
    A = h; B = -hp
    gr = (A*s + B*c)/lam; gz = (A*c - B*s)*lam**2
    p2v = ((lam**2*s**2 + lam**-4*c**2)*np.hypot(gr, gz)).max()
    print(f"lam={lam}: P1 exact {p1:.6f} vs bound {lam*np.abs(h).max():.6f} loss {lam*np.abs(h).max()/p1:.4f}"
          f" | P2 exact {p2v:.6f} vs bound {lam**4*np.hypot(h,hp).max():.6f} loss {lam**4*np.hypot(h,hp).max()/p2v:.4f}")
