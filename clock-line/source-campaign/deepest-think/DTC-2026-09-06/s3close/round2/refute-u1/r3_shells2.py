"""Independent integration of (5.1) with the shift d treated CONTINUOUSLY
(cumulative trapezoid + linear interpolation at sigma+d), so the answer is not
quantised by the grid spacing.  Grid = node values on [0,1].
"""
import math, json
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import CubicSpline

DEG = math.pi/180.0
delta = 7.5*DEG
def h_axis(phi):
    pax = min(phi, math.pi-phi); return min(1.0, pax/delta)
def P_exact(lam):
    A = lam**2 - lam**-4; B = lam**-4
    f = lambda v: h_axis(math.asin(min(v,1.0)))*v*v*(A*v*v+B)**-2.5
    val,_ = quad(f,0.0,1.0,limit=400,epsabs=1e-13,epsrel=1e-13,points=[math.sin(delta)])
    return 3.0*val
lg = np.linspace(1.0,2.6,801); Psp = CubicSpline(lg, np.array([P_exact(l) for l in lg]))

KAPPA = P_exact(1.0)/2.0
CSTAR = math.log(1.5)/KAPPA
LAMMAX = math.exp(0.75*CSTAR)
ELL = math.log(2.0)+2.0*math.log(LAMMAX)

def run(L, N):
    d = ELL/L
    s = np.linspace(0.0,1.0,N+1); ds = 1.0/N
    def rhs(th, x):
        Pv = Psp(np.clip(np.exp(x),1.0,2.6))
        # cumulative integral from sigma to 1, trapezoid, at the nodes
        seg = 0.5*(Pv[:-1]+Pv[1:])*ds
        tail = np.concatenate([np.cumsum(seg[::-1])[::-1],[0.0]])   # tail[i] = INT_{s_i}^1
        # value at s_i + d by linear interpolation of the antiderivative
        sd = np.clip(s+d, 0.0, 1.0)
        val = np.interp(sd, s, tail)
        return 0.5*val
    ev = lambda th,x: x[0]-math.log(1.5); ev.terminal=True; ev.direction=1
    sol = solve_ivp(rhs,(0.0,2.0),np.zeros(N+1),rtol=1e-12,atol=1e-14,events=ev,max_step=0.005)
    return float(sol.t_events[0][0]), d

print("kappa=%.16f c_*=%.12f lam_max=%.12f ell=%.12f"%(KAPPA,CSTAR,LAMMAX,ELL))
U1 = {40:0.7816769,160:0.7479297,640:0.7399101}
for L in [40.0,160.0,640.0]:
    row=[]
    for N in [200,400,800,1600,3200]:
        th,d = run(L,N); row.append((N,th))
    # Richardson-ish extrapolation from the last two
    ext = row[-1][1] + (row[-1][1]-row[-2][1])
    print("L=%5g d=%.7f  "%(L,d) + "  ".join("N=%d:%.8f"%r for r in row) +
          "  | extrap~%.8f  u1(N=200)=%.7f"%(ext, U1[int(L)]))
    # comparators
    s0=0.0
    thC = math.log(1.5)/(KAPPA*(1.0-d))
    print("        theta_C (sigma0=0) = %.8f"%thC)
