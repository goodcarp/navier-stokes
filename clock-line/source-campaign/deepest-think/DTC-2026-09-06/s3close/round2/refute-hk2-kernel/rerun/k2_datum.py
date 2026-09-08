#!/usr/bin/env python3
"""k2 -- the datum's derivative fields, verified; the shell total-variation density; the
geometry of N_tau over the window.  Units rho0 = M = 1."""
import json, math
import numpy as np
from hk2lib import DatumA, DatumB, geometry, THETA_MAX, lam_of_theta

RES = {}
def say(s): print(s)

# ---------------------------------------------------------------- 1. FD verification of DatumB
say("=== 1. DatumB analytic derivatives vs central finite differences ===")
DB = DatumB(delta_deg=7.5, w=0.20, L=10.0)
pts = [(0.5,0.866),(0.75,0.3849),(0.9,0.2),(0.3,1.2),(1.5,2.0),(0.2,0.05)]
worst = 0.0; rows=[]
f = lambda R,Z: DB.eta_rz(np.array([R]),np.array([Z]))[0][0]
def d1(g, h):   # 4th order central
    return (-g(2*h)+8*g(h)-8*g(-h)+g(-2*h))/(12*h)
def d2(g, h):
    return (-g(2*h)+16*g(h)-30*g(0.0)+16*g(-h)-g(-2*h))/(12*h*h)
for (r,z) in pts:
    e,er,ez,err,erz,ezz = [v[0] for v in DB.eta_rz(np.array([r]), np.array([z]))]
    rho = math.hypot(r,z); scale = abs(e)
    h = 1e-3
    fr  = d1(lambda u: f(r+u,z), h);      fz  = d1(lambda u: f(r,z+u), h)
    frr = d2(lambda u: f(r+u,z), h);      fzz = d2(lambda u: f(r,z+u), h)
    frz = d1(lambda u: d1(lambda v: f(r+u,z+v), h), h)
    ana = [er,ez,err,erz,ezz]; num = [fr,fz,frr,frz,fzz]
    sc  = [scale/rho, scale/rho, scale/rho**2, scale/rho**2, scale/rho**2]
    rel = max(abs(a-n)/max(abs(n), 1e-3*s_) for a,n,s_ in zip(ana,num,sc))
    worst = max(worst, rel); rows.append(dict(r=r,z=z,rel=rel))
    say(f"   (r,z)=({r},{z})  worst scaled err = {rel:.3e}")
RES['DB_fd_worst_rel'] = worst
assert worst < 1e-6, worst

# ---------------------------------------------------------------- 2. plateau identities
say("\n=== 2. plateau values (exact, used as the near-ball sups for DatumA) ===")
RES['plateau'] = {"grad": "M/r^2", "hessF": "sqrt(7) M/r^3", "hessOp": "2 M/r^3",
                  "sqrt7": math.sqrt(7)}
say(f"   |grad eta_P| = M/r^2 ;  ||Hess||_F = sqrt(7) M/r^3 = {math.sqrt(7):.7f} M/r^3 ;"
    f"  ||Hess||_op = 2M/r^3")

# ---------------------------------------------------------------- 3. shell TV density J(rho)
say("\n=== 3. shell total-variation density:  |D eta|(shell) = M rho'^2 J drho'  ===")
# J = int_{S^4} sqrt(Hh^2+Hh'^2) dOmega_4  +  jump contributions.
# dOmega_4 = sin^3 phi dphi |S^3| , |S^3| = 2 pi^2.
S3 = 2*math.pi**2; S4 = 8*math.pi**2/3
def J_of_A(delta_deg, n=200001):
    D = DatumA(delta_deg=delta_deg)
    phi = np.linspace(1e-9, math.pi-1e-9, n)
    Hh, Hhp = D.Hh_of_phi(phi)
    integ = np.sqrt(Hh**2 + Hhp**2)*np.sin(phi)**3
    ac = S3*np.trapz(integ, phi)                      # absolutely continuous part
    # equatorial jump: [eta] = 2 M h(pi/2)/r' = 2M/r' on {z=0}; its shell measure per drho'
    #   int_{rho'<|y|<rho'+d} 2M/|y| d^4y = 2M |S^3| rho'^2 drho'  -> J_jump = 2 |S^3|
    jj = 2*S3
    return ac, jj, ac+jj
for dd in (7.5, 15.0):
    ac, jj, tot = J_of_A(dd)
    RES[f'J_A_delta{dd}'] = dict(ac=ac, jump=jj, total=tot)
    say(f"   (D-A) delta={dd:4.1f} deg :  J_ac = {ac:.6f}   J_jump = {jj:.6f}   J = {tot:.6f}")
def J_of_B(DB, n=200001):
    t = np.linspace(-1+1e-12, 1-1e-12, n)
    H0 = DB.Hh(t,0); H1 = DB.Hh(t,1)
    s = np.sqrt(1-t*t)
    # eta = G(rho) Hh(t);  |grad eta| at Theta=1 : radial M Hh/rho^2 , angular (1-t^2)^{1/2}|dHh/dt|... 
    # in (rho,phi): d_phi Hh = -s * dHh/dt ; |grad| = (M/rho^2) sqrt(Hh^2 + s^2 (dHh/dt)^2)
    integ = np.sqrt(H0**2 + (s*H1)**2)*s**2          # sin^3 phi dphi = s^2 dt (since dphi = -dt/s)
    return S3*np.trapz(integ, t)
for dd, ww in ((7.5,0.20),(15.0,0.20)):
    DBx = DatumB(delta_deg=dd, w=ww, L=10.0)
    jb = J_of_B(DBx); RES[f'J_B_delta{dd}_w{ww}'] = jb
    say(f"   (D-B) delta={dd:4.1f} deg, w={ww} :  J = {jb:.6f}   (no jump: profile is smooth)")
# pure plateau reference (no taper, no jump-mollification): Hh = 1/sin, Hh' = -cos/sin^2
phi = np.linspace(1e-9, math.pi-1e-9, 400001); s=np.sin(phi)
Jp = S3*np.trapz(np.sqrt(1+np.cos(phi)**2)/s**2*s**3, phi)
RES['J_plateau_ac'] = Jp
say(f"   bare plateau (no taper): J_ac = {Jp:.6f} = 2 pi^2 (sqrt2 + arcsinh 1) = "
    f"{2*math.pi**2*(math.sqrt(2)+math.asinh(1)):.6f}")

# ---------------------------------------------------------------- 4. geometry over the window
say("\n=== 4. N_tau geometry: the tracked point and its distances over lam in [1,3/2] ===")
tab = []
for lam in (1.0, 1.1, 1.25, 1.5):
    g = geometry(phi0_deg=30.0, f=0.0, delta_deg=7.5, lam=lam)
    tab.append({k: float(v) for k,v in g.items()})
    say(f"   lam={lam:4.2f}  (r,z)=({g['r']:.5f},{g['z']:.5f})  rho={g['rho']:.5f}  "
        f"phi={g['phi_deg']:.3f} deg   d_eq={g['d_eq']:.5f}  d_taper={g['d_taper']:.5f}  "
        f"cone={g['taper_cone_deg']:.3f} deg")
RES['geometry_phi0_30'] = tab
# admissible phi0 range: require d_eq(lam) >= d_min for all lam<=3/2.  d_eq = lam^-2 cos phi0 (f=0)
dmins = {}
for dmin in (0.20, 0.25, 0.30, 0.3827):
    # lam^-2 cos phi0 >= dmin  at lam = 3/2  ->  cos phi0 >= (9/4) dmin
    c = 2.25*dmin
    dmins[dmin] = math.degrees(math.acos(min(1.0, c))) if c <= 1 else None
RES['phi0_max_for_dmin'] = {str(k): v for k,v in dmins.items()}
say("   admissible phi0 (f=0): d_eq(lam) = lam^-2 cos phi0 >= d_min for all lam<=3/2 requires")
for k,v in dmins.items():
    say(f"      d_min = {k:.4f}  ->  phi0 <= {v:.3f} deg" if v is not None else f"      d_min={k}: impossible")
# distance-shrink factor: dist(T_lam x, T_lam S) >= sigma_min(T_lam) dist(x,S), sigma_min = lam^-2
RES['sigma_min_T_lam_at_3over2'] = (1.5)**-2
say(f"   sigma_min(T_lam) at lam=3/2 : {1.5**-2:.6f} = 4/9  (every material distance shrinks by"
    f" at most this)")

# ---------------------------------------------------------------- 5. local sups on B(x,d)
say("\n=== 5. local sups over B(x,d), DatumB, delta=7.5deg w=0.20 L=10, at the tracked point ===")
def sups_on_ball(D, r0, z0, d, n=241):
    # sample the 2D (r,z) disc of radius d about (r0,z0) -- B(x,d) in R^5 meets the half-plane
    # in exactly this disc, and the fields are axisymmetric, so the sup is attained there.
    g = np.linspace(-d, d, n)
    RR, ZZ = np.meshgrid(r0+g, z0+g, indexing='ij')
    m = (RR-r0)**2 + (ZZ-z0)**2 <= d*d
    RR = np.where(m, RR, r0); ZZ = np.where(m, ZZ, z0)
    RR = np.maximum(RR, 1e-6)
    gn = D.grad_norm(RR, ZZ); hf = D.hess_F(RR, ZZ)
    # boundary ring
    th = np.linspace(0, 2*math.pi, 4001)
    gb = D.grad_norm(np.maximum(r0+d*np.cos(th),1e-6), z0+d*np.sin(th))
    return float(np.max(gn)), float(np.max(hf)), float(np.max(gb))
rows5 = []
for lam in (1.0, 1.25, 1.5):
    g = geometry(30.0, 0.0, 7.5, lam)
    d = 0.9*min(g['d_eq'], g['d_taper'])
    s1, s2, sb = sups_on_ball(DB, g['r'], g['z'], d)
    loc = DB.eta_rz(np.array([g['r']]), np.array([g['z']]))
    rows5.append(dict(lam=lam, d=d, sup_grad=s1, sup_hessF=s2, sup_grad_bdry=sb,
                      eta=float(loc[0][0]), eta_r=float(loc[1][0]), eta_z=float(loc[2][0])))
    say(f"   lam={lam:4.2f}  d={d:.5f}  sup_B|grad eta|={s1:.5f}  sup_B||Hess||_F={s2:.5f}  "
        f"sup_dB|grad eta|={sb:.5f}   |eta(x)|={abs(loc[0][0]):.5f}  |grad eta(x)|="
        f"{math.hypot(loc[1][0],loc[2][0]):.5f}")
RES['local_sups_DB'] = rows5

json.dump(RES, open('k2_results.json','w'), indent=1)
say("\nWROTE k2_results.json")
