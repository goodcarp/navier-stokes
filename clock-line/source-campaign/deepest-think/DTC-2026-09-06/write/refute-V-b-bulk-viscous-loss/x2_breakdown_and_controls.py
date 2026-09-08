#!/usr/bin/env python3
"""x2 -- (a) term-by-term breakdown of Theorem V.4's E_tail and E_hess;
        (b) isolate WHICH convention change turns the seat's sec-5 tail table into
            the theorem's own E_tail (threshold d -> d/2 ; 1-D projection -> |Z| in R^5 ;
            bare P -> the 2N*P + P + sum sqrt(P) combination);
        (c) Monte-Carlo control on my own tail quadrature (it must be able to fail);
        (d) mutation control on the seat's check_constants.py."""
import json, math, subprocess, os, shutil, tempfile, re
import numpy as np
from scipy.stats import chi2
from scipy.special import roots_hermitenorm

R = {}
TH = 4*(1-math.sqrt(2/3)); I2 = 4/5 - 16*math.sqrt(6)/135; I4 = -4/7 + 27*math.sqrt(6)/28
C  = math.sqrt(1.5)*TH; n = 5; sd = math.sin(math.radians(7.5))
_hx, _hw = roots_hermitenorm(400); _hw = _hw/_hw.sum()

def P_norm_gt(a, sy, sz):
    if a <= 0: return 1.0
    g = math.sqrt(2*sz)*_hx
    rem = a*a - g*g
    p = np.where(rem <= 0, 1.0, chi2.sf(np.maximum(rem, 0.0)/(2*sy), df=4))
    return float((_hw*p).sum())

# ---------- (c) Monte-Carlo control on P_norm_gt (pre-declared: refuted if any |z| > 4) ----
rng = np.random.default_rng(7)
mc = []
for (sy, sz, a) in [(1.0,3.5,2.0),(1.0,3.5,4.0),(0.05,0.18,0.5),(0.005,0.018,0.25),(2.0,0.5,3.0)]:
    NN = 4_000_000
    Zy = math.sqrt(2*sy)*rng.standard_normal((NN,4)); Zz = math.sqrt(2*sz)*rng.standard_normal(NN)
    p_mc = float(((Zy**2).sum(1)+Zz**2 > a*a).mean())
    p_q  = P_norm_gt(a, sy, sz)
    se = math.sqrt(max(p_mc*(1-p_mc),1e-12)/NN)
    mc.append(dict(sy=sy, sz=sz, a=a, p_mc=p_mc, p_quad=p_q, z=(p_q-p_mc)/se))
R['C_mc_control'] = mc
print("(c) MC control on the tail quadrature (|z| must be <= 4):")
for m in mc: print(f"    sy={m['sy']:6.3f} sz={m['sz']:6.3f} a={m['a']:4.2f}  quad={m['p_quad']:.6e}"
                   f"  mc={m['p_mc']:.6e}  z={m['z']:+.2f}")
assert all(abs(m['z'])<=4 for m in mc), "tail quadrature control FIRED"

# ---------- (a)/(b) breakdown ----------
def pieces(L, f=0.0, phi0_deg=30.0, K2hat=1.0, d=None, Rminus=None):
    M=1.0; rho0=1.0; nu=(rho0*sd)**2*M; tau=TH/(M*L)
    r0=(1+f)*math.sin(math.radians(phi0_deg))
    sy=nu*I2/(M*L); sz=nu*I4/(M*L); sC=sy/r0**2; N=r0/sd
    if d is None: d = rho0*sd
    Rm = (r0-d) if Rminus is None else Rminus
    K2=K2hat*M/rho0; V1=n*nu*tau*(math.exp(2*C)-1)/C; EW=0.5*K2*math.exp(C)*tau*V1
    trC=8*sy+2*sz; trC2=4*(2*sy)**2+(2*sz)**2
    m2=trC; m4=trC**2+2*trC2
    g=math.sqrt(2*sz)*_hx; a2=2*sy
    m6=float((_hw*(a2**3*192+3*a2**2*24*g*g+3*a2*4*g**4+g**6)).sum())
    Ph=P_norm_gt(d/2,sy,sz); Pf=P_norm_gt(d,sy,sz); P=2*Ph+Pf
    t_2NP=2*N*P; t_P=P
    t_sq=[r0**(-k)*math.sqrt([m2,m4,m6][k-1])*math.sqrt(P) for k in (1,2,3)]
    return dict(L=L,f=f,r0=r0,d=d,Rm=Rm,N=N,sy=sy,sz=sz,s_C=sC,
                Eh1=(r0/Rm**2)*EW, Eh2=(4*N*EW/d if d>0 else float('inf')),
                Eh1_Rm_is_r0=EW/r0,
                E4=(r0/Rm**5)*m4, P_half=Ph, P_full=Pf, Pcal=P,
                tail_2NP=t_2NP, tail_P=t_P, tail_sqrt=t_sq, E_tail=t_2NP+t_P+sum(t_sq))

print("\n(a) E_tail term-by-term, d = rho0 sin delta, f=0, phi0=30deg")
print(" L      2N*P        P         k=1 sqrt     k=2 sqrt     k=3 sqrt     E_tail      eps_bulk")
rows=[]
for L in [10,20,40,80,160,320,640,1280]:
    p=pieces(L,K2hat=0.0); rows.append(p)
    print(f"{L:5d}  {p['tail_2NP']:.3e}  {p['tail_P']:.3e}  {p['tail_sqrt'][0]:.3e}  "
          f"{p['tail_sqrt'][1]:.3e}  {p['tail_sqrt'][2]:.3e}  {p['E_tail']:.3e}  {p['s_C']:.3e}")
R['E_tail_breakdown'] = [{k:v for k,v in p.items()} for p in rows]

print("\n(b) which convention change costs what, at d = rho0 sin delta (numbers are probabilities)")
print(" L    seat 1-D thr=d   1-D thr=d/2    |Z|>d (R^5)   |Z|>d/2 (R^5)")
phi0=math.radians(30.0)
rows=[]
for L in [10,20,40,80,160]:
    var_n=2*(I2*math.sin(phi0)**2+I4*math.cos(phi0)**2)/L
    a=0.5*math.exp(-1.0/(2*var_n)); b=0.5*math.exp(-0.25/(2*var_n))
    sy=sd**2*I2/L; sz=sd**2*I4/L
    c_=P_norm_gt(sd,sy,sz); d_=P_norm_gt(sd/2,sy,sz)
    rows.append(dict(L=L, seat_1D_d=a, oneD_dhalf=b, R5_d=c_, R5_dhalf=d_))
    print(f"{L:5d}   {a:.4e}      {b:.4e}    {c_:.4e}    {d_:.4e}")
R['convention_table']=rows

print("\n(a2) E_hess: term1 (R_-=r0), term1 (R_-=r0-d), term2, ratio term2/term1")
print(" L      Eh1(R=r0)    Eh1(R=r0-d)     Eh2          Eh2/Eh1     seat's sized value")
rows=[]
for L in [10,40,160,640]:
    p=pieces(L,K2hat=1.0)
    seat=0.5*math.exp(C)*n*(math.exp(2*C)-1)/C*TH*0.5/L*p['s_C']
    rows.append(dict(L=L, Eh1_r0=p['Eh1_Rm_is_r0'], Eh1=p['Eh1'], Eh2=p['Eh2'],
                     ratio=p['Eh2']/p['Eh1'], seat=seat, full_over_seat=(p['Eh1']+p['Eh2'])/seat,
                     eh1r0_over_seat=p['Eh1_Rm_is_r0']/seat))
    print(f"{L:5d}  {p['Eh1_Rm_is_r0']:.4e}  {p['Eh1']:.4e}  {p['Eh2']:.4e}  "
          f"{p['Eh2']/p['Eh1']:8.2f}   {seat:.4e}")
R['E_hess_breakdown']=rows
print("   Eh1(R_-=r0)/seat's value  =", [round(r['eh1r0_over_seat'],6) for r in rows],
      " (= theta_max/I2 =", round(TH/I2,6), ")")

# ---------- (d) mutation control on the seat's own gate ----------
src = "~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/write/V-b-bulk-viscous-loss"
tmp = tempfile.mkdtemp()
for fn in os.listdir(src):
    if fn.endswith('.json') or fn == 'check_constants.py': shutil.copy(os.path.join(src,fn), tmp)
base = subprocess.run(['python3', os.path.join(tmp,'check_constants.py')], capture_output=True, text=True)
caught = missed = 0; missed_keys=[]
for jf in ['b1_results.json','b2_results.json','b3_results.json','b4_results.json','b5_results.json']:
    orig = json.load(open(os.path.join(tmp,jf)))
    def walk(o, path=()):
        if isinstance(o, dict):
            for k,v in o.items(): yield from walk(v, path+(k,))
        elif isinstance(o, list):
            for i,v in enumerate(o): yield from walk(v, path+(i,))
        elif isinstance(o, (int,float)) and not isinstance(o,bool): yield path
    paths = list(walk(orig))
    step = max(1, len(paths)//60)
    for p in paths[::step]:
        d = json.loads(json.dumps(orig))
        cur = d
        for t in p[:-1]: cur = cur[t]
        cur[p[-1]] = cur[p[-1]]*1.5 + 0.37
        json.dump(d, open(os.path.join(tmp,jf),'w'))
        r = subprocess.run(['python3', os.path.join(tmp,'check_constants.py')], capture_output=True, text=True)
        if r.returncode != 0: caught += 1
        else: missed += 1; missed_keys.append(jf+":"+"/".join(str(t) for t in p))
    json.dump(orig, open(os.path.join(tmp,jf),'w'))
R['gate_mutation'] = dict(baseline_rc=base.returncode, caught=caught, missed=missed,
                          missed_fraction=missed/max(caught+missed,1), missed_examples=missed_keys[:25])
print(f"\n(d) mutation of the seat's own gate: baseline rc={base.returncode}, "
      f"caught {caught}, missed {missed} ({100*missed/max(caught+missed,1):.1f}%)")
for k in missed_keys[:15]: print("    missed:", k)

json.dump(R, open('x2_results.json','w'), indent=1, default=str)
print("\nWROTE x2_results.json")
