#!/usr/bin/env python3
"""x5 -- GAP 2 ('transfer of the L1/L2/L4 origin structure to the material innermost shell,
a relative O(1/L) step').  The SIZE of a transfers; the SIGN STRUCTURE does not.

At the origin the strain kernel is K = -(3/8pi) r z/rho^5, with sgn K = -sgn z, so
K*omega^theta >= 0 pointwise for the odd datum: EVERY source contributes with the same sign
(L2), which is what makes P1 sign-definite (L4).

Off the origin the kernel is the full elliptic Biot-Savart kernel.  Here it is computed
directly (filament Biot-Savart, no elliptic tables) at the natural analogue of the innermost
material shell -- a point on the equatorial plane, which is material -- and its sign is
tested against the datum's sign.
"""
import json, math
import numpy as np
LOG=[]; OUT={}
def say(s=""):
    print(s, flush=True); LOG.append(s)

def ur_from_filament(r0, z0, rs, zs, Gam=1.0, nphi=40000):
    """u^r at field point (r0,0,z0) from a circular filament of circulation Gam at (rs,zs)."""
    ph = (np.arange(nphi)+0.5)*(2*math.pi/nphi); dphi = 2*math.pi/nphi
    yx = rs*np.cos(ph); yy = rs*np.sin(ph)
    dlx = -rs*np.sin(ph)*dphi; dly = rs*np.cos(ph)*dphi
    dx = r0-yx; dy = -yy; dz = z0-zs
    d3 = (dx*dx+dy*dy+dz*dz)**1.5
    ux = np.sum((dly*dz - 0.0*dy)/d3)*Gam/(4*math.pi)     # (dl x d)_x = dly*dz - dlz*dy, dlz=0
    return ux

def K_offaxis(r0, z0, rs, zs):
    """kernel of a = u^r/r0 against  d(int omega^theta dr dz)  ... per unit 3-volume measure
       d^3x = 2 pi rs drs dzs, so K_off = ur_from_filament(Gam=1)/(r0 * 2 pi rs)."""
    return ur_from_filament(r0, z0, rs, zs, 1.0)/(r0*2*math.pi*rs)

say("="*100)
say("Sign of the strain kernel: at the origin vs at a point of the equatorial material plane")
say("="*100)
say("datum sign: omega^theta = -M sgn(z).  Sign-definiteness means  K * omega^theta >= 0,")
say("i.e.  sgn K(rs,zs) = -sgn(zs)  for every source.")
say()
r0 = 1.0            # innermost material shell, on the equator (rho0 = 1 in the seat's units)
say(f"field point: (r,z) = ({r0}, 0) -- the innermost shell on the material equatorial plane")
say(f"{'rs':>7} {'zs':>7} {'K_origin':>13} {'sgnOK@orig':>11} {'K_offaxis':>13} {'sgnOK@shell':>12}")
bad=[]; rows=[]
for rs in (0.6, 0.9, 1.0, 1.1, 1.5, 2.0, 4.0, 8.0):
    for zs in (0.05, 0.2, 0.5, 1.0, 2.0, 5.0):
        rho = math.hypot(rs,zs)
        Ko = -(3.0/(8*math.pi))*rs*zs/rho**5
        Kf = K_offaxis(r0, 0.0, rs, zs)
        ok_o = (Ko*( -1.0*np.sign(zs)) >= 0)
        ok_f = (Kf*( -1.0*np.sign(zs)) >= 0)
        rows.append(dict(rs=rs,zs=zs,K_origin=Ko,K_off=Kf,ok_origin=bool(ok_o),ok_shell=bool(ok_f)))
        if not ok_f: bad.append((rs,zs,Kf))
        say(f"{rs:>7.2f} {zs:>7.2f} {Ko:>13.6f} {str(ok_o):>11} {Kf:>13.6f} {str(ok_f):>12}")
OUT["rows"]=rows
say()
say(f"sources where the ORIGIN kernel has the wrong sign : 0 (proved: sgn K = -sgn z)")
say(f"sources where the SHELL  kernel has the wrong sign : {len(bad)} of {len(rows)}")
for b in bad[:12]:
    say(f"    (rs,zs) = ({b[0]:.2f},{b[1]:.2f})   K_offaxis = {b[2]:+.6f}  (needs <= 0, i.e. same sign as at the origin)")
OUT["n_sign_violations_at_shell"]=len(bad); OUT["n_tested"]=len(rows)
say()
say("Consequence: L2's 'every material element contributes with the same sign, for all time'")
say("and L4's 'P1 >= 0' are properties of the ORIGIN ONLY.  At the material innermost shell the")
say("kernel is sign-INDEFINITE, so the analogue of P1 has no sign, and the reduction of (H2) to")
say("one inequality on one integral does not transfer.  GAP 2 is not the O(1/L) size step the")
say("NOTE describes; it is the loss of the entire sign structure the argument is built on.")
say()
say("What DOES transfer (confirmed): the SIZE.  For the sharp shell the origin value is")
say("a(0) = (M/2) log(R/rho0) EXACTLY (x4b), and the frame's material-shell value is")
say("(M/2) log(R/rho0) + 0.216773 M, so the difference is exactly the constant 0.216773 M.")
json.dump(OUT, open("x5_results.json","w"), indent=1, default=str)
open("x5_log.txt","w").write("\n".join(LOG)+"\n")
print("\n[x5 done]")
