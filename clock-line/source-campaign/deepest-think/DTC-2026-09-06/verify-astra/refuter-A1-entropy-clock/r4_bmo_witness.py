#!/usr/bin/env python3
"""R4: numeric witness that the SEMINORM-ONLY weighted-BMO step used in
Bradshaw-Farhat-Grujic (arXiv:1704.05546v4) Thm 8's sketch --
   int |f(z)| /(|z|+1)^4 dz  <=  c ||f||_BMO
-- is FALSE with a universal c, and that the deficit is EXACTLY a
logarithm of a scale ratio: i.e. the quantity Astra's (5)-(6) restores.
Witness family: f_R(z) = log_+(R/|z|)  (truncated log), which is the
saturating example for  |grad u| ~ ||omega||_inf * log(scale ratio)."""
import mpmath as mp
mp.mp.dps=25

def W(R):
    # int_{R^3} log_+(R/|z|)/(1+|z|)^4 dz = 4 pi int_0^R r^2 log(R/r)/(1+r)^4 dr
    return 4*mp.pi*mp.quad(lambda r: r**2*mp.log(R/r)/(1+r)**4,[0,1,R])

def mean_ball(R,rho):
    # centered ball B(0,rho): average of f_R
    if rho<=0: return mp.mpf(0)
    num=mp.quad(lambda r: r**2*(mp.log(R/r) if r<R else mp.mpf(0)),[0,min(rho,R),rho] if rho>R else [0,rho])
    return 3*num/rho**3

def osc_ball(R,rho):
    m=mean_ball(R,rho)
    num=mp.quad(lambda r: r**2*abs((mp.log(R/r) if r<R else mp.mpf(0))-m),[0,min(rho,R),rho] if rho>R else [0,rho])
    return 3*num/rho**3

print("R        weighted-int W(R)     max centered mean-osc over rho in [1e-6,1e6]   W(R)/log R")
for k in [1,2,3,4,6,8,10]:
    R=mp.mpf(10)**k
    w=W(R)
    rhos=[mp.mpf(10)**e for e in range(-6,7)]
    o=max(osc_ball(R,rho) for rho in rhos)
    print(f"1e{k:<3d}  {float(w):>16.6f}   {float(o):>14.6f}                    {float(w/mp.log(R)):>10.6f}")
print()
print("READING: the mean-oscillation (BMO-seminorm proxy) stays O(1) uniformly in R,")
print("while the weighted integral grows LINEARLY in log R. Hence no universal c can")
print("make  int |f|/(1+|z|)^4 <= c ||f||_BMO.  The correct inequality carries the")
print("local mean |f_B| (BFG's own parenthetical), and bounding |f_B| for f=grad u")
print("needs a LENGTH, which ||omega||_inf alone does not supply. The cheapest length")
print("is Astra's ell = E^{1/5} M^{-2/5}; at the heat scale sqrt(nu(t-s)) < ell the")
print("mean costs exactly C*M*log(ell/sqrt(nu(t-s))) -- Astra's (5)-(6)+(14) logarithm.")
