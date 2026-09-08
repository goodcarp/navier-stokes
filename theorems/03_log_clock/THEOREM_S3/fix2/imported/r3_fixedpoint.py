"""R3: independent re-derivation of the window fixed point, L_* and C''."""
import math, numpy as np
C1, C2 = 0.291998532549, 0.0147542715827
RA2 = 2*(32-1/32)**0.2
E0, G0 = 7.66060196237754, 20.91971219785737
RIESZ = 3*math.pi**2/16
c = 2*math.log(1.5)
def sexp(x):
    return math.exp(x) if x < 700 else float('inf')
def Chat_a(lam, sig, cG):
    A = math.inf if sig <= 0 else (RA2/sig + math.pi/8)*lam
    B = RA2*sexp(cG)*E0
    return (C1+C2)*lam + min(A, B)
def F(Gbar, L, lam, sig, mode='proved', mu=0.0):
    cG = Gbar*c/L
    if cG > 700: return float('inf'), None
    Ca = Chat_a(lam, sig, cG)
    a0 = 0.5*lam*L
    p = 1.0
    for _ in range(2000):
        fac = sexp(3*cG) if mode=='p3' else (((1+mu)*lam)**2*sexp(cG) if mode=='map' else sexp(p*cG))
        if not np.isfinite(fac): return float('inf'), None
        Gh = RIESZ*fac*G0
        A = a0 + Ca
        Grad = max(A, A+2*Gh+lam/2, 2*Ca+2*Gh+lam/2)
        pn = 1 + 2*min(Grad, Gbar)/Gbar
        if abs(pn-p) < 1e-14: p = pn; break
        p = pn
    fac = sexp(3*cG) if mode=='p3' else (((1+mu)*lam)**2*sexp(cG) if mode=='map' else sexp(p*cG))
    if not np.isfinite(fac): return float('inf'), None
    Gh = RIESZ*fac*G0
    Cpp = 2*Ca + lam + Gh
    return lam*L + Cpp, dict(Cpp=Cpp, cG=cG, p=p, Ghat=Gh, Ca=Ca)
def smallest_fp(L, lam, sig, mode='proved', mu=0.0):
    base = lam*L
    for g in np.geomspace(base*(1+1e-14), base*1e10, 3000):
        f,_ = F(g, L, lam, sig, mode, mu)
        if f <= g:
            a, b = base, g
            for _ in range(400):
                m = 0.5*(a+b)
                f2,_ = F(m, L, lam, sig, mode, mu)
                if f2 <= m: b = m
                else: a = m
            f3, info = F(b, L, lam, sig, mode, mu)
            info['Gbar'] = b; return info
    return None
def Lstar(lam, sig, mode='proved', mu=0.0):
    lo, hi = 10.0, 1e9
    if smallest_fp(hi, lam, sig, mode, mu) is None: return None
    for _ in range(120):
        m = math.sqrt(lo*hi)
        if smallest_fp(m, lam, sig, mode, mu) is None: lo = m
        else: hi = m
        if hi/lo < 1+1e-12: break
    return hi
def show(tag, lam, sig, mode='proved', mu=0.0):
    Ls = Lstar(lam, sig, mode, mu); i = smallest_fp(1e4, lam, sig, mode, mu)
    print(f"{tag:34s} L*={Ls if Ls is None else round(Ls,1)}  "
          f"C''(1e4)={i['Cpp']:.2f} c_G={i['cG']:.4f} p={i['p']:.4f} Chat_a={i['Ca']:.4f}")
for lam in (1.0,1.25,1.5): show(f"proved sig=1/2 lam={lam}", lam, 0.5)
for lam in (1.0,1.25,1.5): show(f"proved sig=0   lam={lam}", lam, 0.0)
show("crude p3 sig=1/2 lam=1.5", 1.5, 0.5, 'p3')
for mu in (0.0,0.02,0.05): show(f"map mu={mu} sig=1/2 lam=1.5", 1.5, 0.5, 'map', mu)
s7 = math.sin(math.radians(7.5))
show("proved sig=sin7.5 lam=1.5", 1.5, s7)
show("map0.05 sig=sin7.5 lam=1.5", 1.5, s7, 'map', 0.05)
show("map0.05 sig=0 lam=1.5", 1.5, 0.0, 'map', 0.05)
print("--- 8.1 frozen c_G = lam*c ---")
for lam in (1.0,1.25,1.5):
    cG = lam*c
    for sig,tag in ((0.5,'half'),(0.0,'axis')):
        Ca = Chat_a(lam,sig,cG)
        for nm,fac in (('p3',math.exp(3*cG)),('p2',math.exp(2*cG)),('map0.05',((1.05)*lam)**2*math.exp(cG))):
            Gh = RIESZ*fac*G0
            print(f"  lam={lam} {tag:4s} {nm:8s} Chat_a={Ca:.4f} Ghat={Gh:.2f} C''={2*Ca+lam+Gh:.2f}")
print("--- boundary fixed point at L* (lam=1.5) ---")
for sig,tag in ((0.5,'half'),(0.0,'axis')):
    Ls = Lstar(1.5,sig); i = smallest_fp(Ls*(1+1e-9),1.5,sig)
    print(f"  {tag}: L*={Ls:.1f} C''={i['Cpp']:.2f} c_G={i['cG']:.4f} p={i['p']:.4f}")
print("--- crossover sin phi ---")
for lam,cG in ((1.5,1.5*c),(1.0,1.0*c)):
    print(f"  lam={lam}: sin phi = {RA2*lam/(RA2*math.exp(cG)*E0 - math.pi*lam/8):.8f} "
          f"= {math.degrees(math.asin(RA2*lam/(RA2*math.exp(cG)*E0 - math.pi*lam/8))):.4f} deg")
