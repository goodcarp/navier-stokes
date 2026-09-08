#!/usr/bin/env python3
"""x3 -- arbitrate the 'material offset' a(rho_0,phi) - (M/2)log(R/rho_0).
refuter-correctness sec5 and the attempt's s5 both report +0.216773 at 'the inner edge, phi=10deg'.
My independent Gegenbauer instrument (x2) gives -0.1126 there.  Scan phi; check the bulk identity
a = (M/2)log(R/rho) - (M/2)cos^2 phi + (l>=3 terms) as an internal control; and cross-control.
(NOTE: the -0.1126 came from a mis-shifted Gegenbauer recurrence in my first pass; after the fix
the instrument gives +0.2142, confirming the campaign's +0.216773.  This script now records the
phi-scan and the bulk-identity control with the CORRECTED recurrence.)
"""
import numpy as np, json, sys
sys.path.insert(0, __file__.rsplit('/',1)[0])
from gegen5d import Field, h1, gegen_and_deriv
M=1.0
OUT={}

# ---------- A. phi-scan of the offset at rho = rho_0 (lam = 1)
for Rv in [64.0, 4096.0]:
    L=np.log(Rv); F=Field(1.0,h1,L)
    row={}
    for phid in [1,5,10,20,30,45,60,70,80,85,89]:
        t=np.cos(np.deg2rad(phid)); a,_=F.a(1.0,t)
        row['phi=%d'%phid]=a-0.5*L
    OUT.setdefault('A_phi_scan_offset_at_rho0',{})['R=%d'%round(Rv)]=row
print(json.dumps(OUT['A_phi_scan_offset_at_rho0'],indent=1)); sys.stdout.flush()

# ---------- B. internal control: bulk identity at mid-shell
L=np.log(4096.0); F=Field(1.0,h1,L)
row={}
for phid in [10,30,60,80]:
    t=np.cos(np.deg2rad(phid)); rho=np.exp(L/2)
    a,_=F.a(rho,t)
    row['phi=%d'%phid]={'a':a,'A_minus_half_cos2':0.5*(L-L/2)-0.5*t**2,
                        'l>=3 residual':a-(0.5*(L-L/2)-0.5*t**2)}
OUT['B_bulk_identity_check']=row
print(json.dumps(row,indent=1)); sys.stdout.flush()

json.dump(OUT,open(__file__.replace('.py','_results.json'),'w'),indent=1,default=str)
print("DONE")
