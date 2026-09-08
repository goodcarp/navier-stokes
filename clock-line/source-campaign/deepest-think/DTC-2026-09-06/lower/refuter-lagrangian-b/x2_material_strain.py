#!/usr/bin/env python3
'''x2 driver -- see gegen5d.py for the instrument.'''
import numpy as np, json, sys
from gegen5d import Field, h1, htap
M=1.0
OUT={}
# ---------------- VALIDATION at lam = 1 against the attempt / refuter-correctness
val={}
for Rv in [64.0, 256.0, 1024.0, 4096.0]:
    L=np.log(Rv); F=Field(1.0,h1,L)
    a0,_ = F.a(1e-4, np.cos(np.deg2rad(10.0)))
    am,pp= F.a(1.0,   np.cos(np.deg2rad(10.0)))
    am4,_= F.a(1.0,   np.cos(np.deg2rad(10.0)), lmax=401)
    val['R=%d'%round(Rv)]={'a0_minus_halfL':a0-0.5*L,'material_offset':am-0.5*L,
                           'material_offset_lmax401':am4-0.5*L}
OUT['validation_lam1']=val
print("VALIDATION lam=1:"); print(json.dumps(val,indent=1)); sys.stdout.flush()

# ---------------- PROBE: strain at the tracked material point vs lam
phi0=np.deg2rad(30.0); probe={}
for Rv in [256.0, 4096.0, 65536.0]:
    L=np.log(Rv); row={}
    for lam in [1.0,1.1,1.25,1.4,1.5]:
        F=Field(lam,h1,L)
        r=lam*np.sin(phi0); z=lam**-2*np.cos(phi0)
        rho_ev=float(np.hypot(r,z)); t_ev=float(z/rho_ev)
        am,_=F.a(rho_ev,t_ev); a0,_=F.a(1e-4,0.5)
        row['lam=%g'%lam]={'rho_ev':rho_ev,'phi_ev_deg':float(np.rad2deg(np.arccos(t_ev))),
            'a_origin':a0,'a_origin_minus_half_lam_L':a0-0.5*lam*L,
            'a_material':am,'a_material_minus_half_lam_L':am-0.5*lam*L,
            'a_material_over_half_lam_L':am/(0.5*lam*L)}
    probe['R=%d'%round(Rv)]=row
    print("PROBE R=%d"%round(Rv)); print(json.dumps(row,indent=1)); sys.stdout.flush()
OUT['material_probe']=probe
json.dump(OUT,open(__file__.replace('.py','_results.json'),'w'),indent=1,default=str)
print("DONE")
