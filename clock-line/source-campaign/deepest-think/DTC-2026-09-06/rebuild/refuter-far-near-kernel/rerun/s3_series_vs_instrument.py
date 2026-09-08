#!/usr/bin/env python3
"""S3 - KILL K1: does the exact shell/zonal series reproduce the validated Biot-Savart instrument?
Also KILL K6: the exact closed form for the inner-edge plateau constant c_edge(phi) vs +0.216773.
Instrument: /Users/.../toys/g2_selfstretch_biot_savart.py (calibrated 1e-5 on the axis identity)."""
import numpy as np, math, importlib.util, json, hashlib, sys, time
sys.path.insert(0, '.')
from kern import geg_table, Nl, bangbang_g, a_series, cesaro
G2 = "~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/toys/g2_selfstretch_biot_savart.py"
spec = importlib.util.spec_from_file_location("g2", G2); g2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(g2)

t0 = time.time(); LMAX = 601
g = bangbang_g(LMAX)
print("g_l (=h_l) check vs S1 exact rationals:")
exact = {1:-5/6, 3:3/40, 5:-247/1680, 7:1513/40320, 9:-2773/42240, 11:217075/9225216, 13:-3561403/92252160}
worst = 0.0
for l, v in exact.items():
    d = abs(g[l]-v); worst = max(worst, d); print(f"   l={l:2d}  quad {g[l]:+.12f}  exact {v:+.12f}  diff {d:.2e}")
print(f"   worst |diff| = {worst:.2e}")
# Parseval:  SUM g_l^2 N_l = INT_{-1}^1 w^2 dt = 2 M^2
ls = np.arange(LMAX+1); P = float(np.sum(g**2*Nl(ls)))
print(f"Parseval  SUM_l g_l^2 N_l = {P:.6f}  (exact 2.000000; truncation at L={LMAX} loses the tail)")

rho0, R = 1.0, 4096.0
print("\n--- K1: exact series vs g2 Biot-Savart instrument (three test points) ---")
print(f"{'rho':>7} {'phi/pi':>8} {'a_series(L=601)':>16} {'a_series(Cesaro)':>17} {'a_instrument':>13} {'diff':>11}")
rows = []; k1 = False
for rho, phid in [(8.0, 1/6), (128.0, 0.45), (512.0, 1/3)]:
    phi = phid*math.pi
    terms = a_series(rho, phi, rho0, R, g, LMAX)
    Pl, Cs = cesaro(terms)
    r_, z_ = rho*math.sin(phi), rho*math.cos(phi)
    ai = g2.strain_at(r_, z_, rho0, R)
    d = Cs - ai
    rows.append(dict(rho=rho, phi_over_pi=phid, a_series=Pl, a_cesaro=Cs, a_instrument=ai, diff=d))
    print(f"{rho:7.1f} {phid:8.4f} {Pl:16.8f} {Cs:17.8f} {ai:13.8f} {d:+11.2e}")
    if abs(d) > 1e-3: k1 = True
print("KILL K1", "FIRED" if k1 else "did not fire", "(threshold 1e-3 absolute)")

print("\n--- axis identity control: a(0,0) must be (M/2) log(R/rho0) ---")
ax = g2.strain_at(0.0, 0.0, rho0, R); print(f"   instrument {ax:.8f}  exact {0.5*math.log(R/rho0):.8f}  rel {abs(ax-0.5*math.log(R/rho0))/(0.5*math.log(R/rho0)):.2e}")

print("\n--- K6: the plateau constant.  c_edge(phi) = -SUM_{l odd>=3} ((l+2)g_l/((2l+3)(l-1))) C_{l-1}(cos phi) ---")
def c_edge(phi, LM=LMAX):
    t = math.cos(phi); C = geg_table(np.array([t]), LM+1)[:, 0]
    ls = np.arange(3, LM+1, 2)
    return cesaro(-((ls+2)*g[ls]/((2*ls+3)*(ls-1)))*C[ls-1])
for phid in [5, 10, 15, 30, 45, 60, 90]:
    Pl, Cs = c_edge(math.radians(phid))
    print(f"   phi={phid:3d} deg   partial {Pl:+.6f}   Cesaro {Cs:+.6f}")
p10, c10 = c_edge(math.radians(10.0))
print(f"\n   c_edge(10 deg) = {c10:+.6f}   refuter's measured value +0.216773   diff {c10-0.216773:+.2e}")
k6 = abs(c10-0.216773) > 2e-3
print("KILL K6", "FIRED" if k6 else "did not fire", "(threshold 2e-3)")
# finite-(R/rho0) version, exactly the refuter's geometry: rho = rho0(1+1e-7), R = 64,256,1024,4096
print("\n   finite-R check at the refuter's own geometry (rho = rho0+, phi = 10 deg):")
for Rv in [64.0, 256.0, 1024.0, 4096.0]:
    terms = a_series(1.0000001, math.radians(10.0), 1.0, Rv, g, LMAX)
    Pl, Cs = cesaro(terms)
    print(f"      R={Rv:7.0f}  a_series(Cesaro) {Cs:.6f}   (M/2)logR {0.5*math.log(Rv):.6f}   offset {Cs-0.5*math.log(Rv):+.6f}")
json.dump(dict(rows=rows, c_edge_10deg=c10, parseval=P, worst_g_err=worst), open('s3_results.json','w'), indent=1)
print(f"\nelapsed {time.time()-t0:.1f}s  SCRIPT-SHA256 {hashlib.sha256(open(__file__,'rb').read()).hexdigest()}")
