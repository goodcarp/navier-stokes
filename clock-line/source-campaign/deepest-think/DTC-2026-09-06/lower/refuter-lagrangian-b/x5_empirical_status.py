#!/usr/bin/env python3
"""x5 -- what the only viscous NS data set says about the ACCELERATED constant.
Inputs are the measured numbers reported in
  .../DTC-2026-09-06/sharp/viscous-numerics/NOTE.md
(Q' = T32 * M * log(R/s*) = 0.988 +- 0.018 over N = 2..7; t=0 strain law kappa = 0.437/0.415/0.411
at N = 3/5/6 for datum A), transcribed here and combined with the two model predictions:
   frozen strain     Q'_frozen = log(3/2)/kappa
   accelerated (4.1) Q'_accel  = 2(1-sqrt(2/3))/kappa
No new physics; this is an arithmetic comparison, done in a script so the numbers are not typed.
"""
import numpy as np, json
Qmeas, Qsd = 0.988, 0.018
kappas = {'N=3':0.437,'N=5':0.415,'N=6':0.411}
lg = np.log(1.5); acc = 2*(1-np.sqrt(2/3))
OUT={'measured_Qprime':Qmeas,'sd':Qsd,'log32':lg,'accel_theta_kappa':acc}
rows={}
for k,v in kappas.items():
    qf, qa = lg/v, acc/v
    rows[k]={'kappa':v,'Qprime_frozen':qf,'Qprime_accel':qa,
             'z_frozen':(qf-Qmeas)/Qsd,'z_accel':(qa-Qmeas)/Qsd}
OUT['per_kappa']=rows
OUT['c2_frozen']=4*lg; OUT['c2_accel']=8*(1-np.sqrt(2/3))
OUT['claimed_improvement_fraction']=1-OUT['c2_accel']/OUT['c2_frozen']
OUT['model_shortfall_at_kappa_0.411']=1-rows['N=6']['Qprime_accel']/Qmeas
# what c2 the measurement itself implies, if T32*M*log(R/s*) -> const and log Re_E = 2 log(R/s*)+c
OUT['implied_c2_from_data_if_2L']=2*Qmeas
print(json.dumps(OUT,indent=1))
json.dump(OUT,open(__file__.replace('.py','_results.json'),'w'),indent=1,default=str)
