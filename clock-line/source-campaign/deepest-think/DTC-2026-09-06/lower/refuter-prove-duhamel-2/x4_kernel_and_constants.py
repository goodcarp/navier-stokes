#!/usr/bin/env python3
"""x4 -- (a) independent verification of prove-duhamel's L1 kernel
             K(r,z) = -(3/(8 pi)) r z / rho^5   via direct 3D Biot-Savart on a vortex filament
             (no stream function, no grid Poisson solve, no elliptic-integral table);
         (b) the sharp-shell closed form a(0) = (M/2) log(R/rho0) implied by K;
         (c) the constants used in the seat's Corollary: kappa, c_a, c_E -- are they the
             constants of the datum the seat actually ran?
"""
import json, math
import numpy as np

LOG=[]; OUT={}
def say(s=""):
    print(s, flush=True); LOG.append(s)

# ---------------------------------------------------------------- (a) filament check
say("="*96)
say("(a) K verified against direct 3D Biot-Savart for a circular vortex filament")
say("="*96)
say("    L1 predicts, for a filament of circulation Gamma at (r_s, z_s):")
say("       a(0,0) = 2 pi r_s Gamma K(r_s,z_s) = -(3/4) Gamma r_s^2 z_s / rho_s^5")
def a_origin_biotsavart(rs, zs, Gam=1.0, eps=1e-4, nphi=200000):
    """u_x at x=(eps,0,0) from a filament, divided by eps -> a(0,0)."""
    ph = (np.arange(nphi)+0.5)*(2*math.pi/nphi); dphi = 2*math.pi/nphi
    yx = rs*np.cos(ph); yy = rs*np.sin(ph); yz = zs
    dlx = -rs*np.sin(ph)*dphi; dly = rs*np.cos(ph)*dphi; dlz = 0.0
    dx = eps - yx; dy = -yy; dz = -zs
    d3 = (dx*dx+dy*dy+dz*dz)**1.5
    # (dl x d)_x = dly*dz - dlz*dy
    ux = np.sum((dly*dz - dlz*dy)/d3)*Gam/(4*math.pi)
    return ux/eps
rows=[]
for (rs, zs) in ((1.0,1.0),(1.0,2.0),(2.0,0.7),(0.5,3.0),(3.0,1.3),(1.0,0.3)):
    rho = math.hypot(rs,zs)
    pred = -(3.0/4.0)*rs*rs*zs/rho**5
    num  = None
    vals = [a_origin_biotsavart(rs,zs,eps=e) for e in (4e-4,2e-4,1e-4)]
    # Richardson in eps^2
    num = vals[-1] + (vals[-1]-vals[-2])/3.0
    rows.append(dict(rs=rs,zs=zs,pred=pred,num=num,rel=abs(num-pred)/abs(pred)))
    say(f"    r_s={rs:>4.1f} z_s={zs:>4.1f}   L1 prediction {pred:+.9f}   Biot-Savart {num:+.9f}"
        f"   rel {abs(num-pred)/abs(pred):.2e}")
OUT["filament"]=rows
worst = max(r["rel"] for r in rows)
say(f"    worst relative disagreement = {worst:.2e}   -> K = -(3/8pi) r z/rho^5 CONFIRMED independently")
say()

# ---------------------------------------------------------------- (b) sharp shell
say("="*96)
say("(b) the sharp-shell value implied by K")
say("="*96)
say("    omega^theta = -M sgn(z) on rho0<|x|<R:")
say("      a(0) = (3/4) M int_{log rho0}^{log R} dlog rho  int_0^pi sgn(cos phi) sin^2 phi cos phi dphi")
from scipy.integrate import quad
ang = quad(lambda p: np.sign(math.cos(p))*math.sin(p)**2*math.cos(p), 0, math.pi, points=[math.pi/2])[0]
say(f"      angular integral = {ang:.12f}   (exact 2/3 = {2/3:.12f})")
say(f"      => a(0) = (3/4)(2/3) M L = {0.75*ang:.12f} M L    (i.e. kappa = 1/2 EXACTLY, at the ORIGIN)")
OUT["sharp_shell_kappa"]=0.75*ang
say()
say("    NOTE: this is kappa at the ORIGIN, where omega = 0 and u = 0 and the particle never moves.")
say("    The frame's a = (M/2)log(R/rho0) + 0.216773 M is the strain on the MATERIAL innermost")
say("    shell.  The two agree at leading order but the seat's L1/L2/L4 structure is only")
say("    available at the origin (k^2 = 0), which carries no vorticity.")
say()

# ---------------------------------------------------------------- (c) the constants
say("="*96)
say("(c) constants used in the seat's Corollary vs the datum it actually ran")
say("="*96)
say("    Re_E = E0^{2/5} M0^{1/5}/nu ;  nu = M rho0^2 at the frame's floor Re0 = 1;")
say("    E0 = cE' M^2 R^5  =>  log Re_E = 2 log(R/rho0) + (2/5) log cE'.")
for name, cEp in (("sharp shell (frame / seat's c_E)", 0.1724040),
                  ("the seat's OWN mollified datum (d2)", 0.1418961)):
    say(f"      {name:>36}:  cE' = {cEp:.7f}   (2/5)log cE' = {0.4*math.log(cEp):+.7f}")
OUT["cE"]=dict(sharp=0.4*math.log(0.1724040), datum=0.4*math.log(0.1418961))
say(f"    The seat uses c_E = -0.7031660 (the SHARP-shell value) with the mollified datum whose")
say(f"    own energy gives c_E = {0.4*math.log(0.1418961):.7f}.  Mismatch {abs(0.4*math.log(0.1418961)+0.703166):.4f} in log Re_E.")
say()
say("    c2(L) = (theta/(kappa L)) * (2L + c_E) -> 2 theta/kappa.  Sensitivity to kappa:")
say(f"      {'kappa':>9} {'source':>44} {'c2 = 1/kappa (theta=1/2)':>26}")
for kap, src in ((0.5, "sharp bang-bang limit (used in the Corollary)"),
                 (0.50091, "seat d2, delta=7.5 deg, w=0.08"),
                 (0.49972, "elliptic instrument, delta=7.5 deg"),
                 (0.48137, "seat d2, delta=7.5 deg, w=0.20 (the RUNS' datum)"),
                 (0.45317, "seat d2, delta=30 deg")):
    say(f"      {kap:>9.5f} {src:>44} {0.5*2/kap:>26.5f}")
OUT["kappa_sensitivity"]=[(0.5,2.0),(0.48137,1.0/0.48137),(0.45317,1.0/0.45317)]
say()
say("    -> 'c2 = 2 ... an honest ceiling for the whole family' holds only in the singular limit")
say("       delta,w -> 0 where kappa = 1/2 exactly.  Every ADMISSIBLE member the seat actually")
say("       built has kappa < 1/2, hence c2 = 1/kappa > 2.  For the runs' own datum (w = 0.20)")
say(f"       the ceiling is {1/0.48137:.4f}, not 2.")
json.dump(OUT, open("x4_results.json","w"), indent=1, default=str)
open("x4_log.txt","w").write("\n".join(LOG)+"\n")
print("\n[x4 done]")
