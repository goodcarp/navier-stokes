"""r4 - is the seat's CONTROL C1 actually a diffeomorphism?

NOTE sec 5 and s4_controls.py both assert:
   "This is a genuine C^1 SO(4)-equivariant diffeomorphism (1 + beta' >= 1 - mu_beta > 0)".
beta(phi) = mu_b m(phi/phi_c) m((pi-phi)/phi_c) cos(phi),  m(t)=t^2(3-2t) on [0,1] then 1.
The chain rule puts terms of size mu_b/phi_c into beta', so the QUOTED inequality cannot be right.
What actually matters for (H1) is (a) min(1+beta') > 0, (b) phi + beta maps [0,pi] onto [0,pi].
Both are tested here on a grid refined inside the ramp, at every (mu_b, phi_c) the seat ran.
Also reported: the true min of 1+beta' vs the seat's claimed lower bound 1-mu_b.
"""
import json, numpy as np
def m(t):  return np.where(t < 1, t**2*(3-2*t), 1.0)
def dm(t): return np.where(t < 1, 6*t*(1-t), 0.0)
def beta(p, mub, pc):  return mub*m(p/pc)*m((np.pi-p)/pc)*np.cos(p)
def dbeta(p, mub, pc):
    A, B = m(p/pc), m((np.pi-p)/pc)
    dA, dB = dm(p/pc)/pc, -dm((np.pi-p)/pc)/pc
    return mub*((dA*B + A*dB)*np.cos(p) - A*B*np.sin(p))

rows = []
for mub in (0.9, 0.5, 0.3):
    for pc in (1e-1, 1e-2, 1e-4, 1e-8, 1e-12):
        # grid: dense inside both ramps plus uniform bulk
        ramp = np.linspace(0, min(pc, np.pi/2), 200001)
        p = np.unique(np.concatenate([ramp, np.pi-ramp, np.linspace(0, np.pi, 200001)]))
        d = dbeta(p, mub, pc); b = beta(p, mub, pc); img = p + b
        rows.append(dict(mu_beta=mub, phi_c=pc,
                         min_1_plus_dbeta=float(np.min(1+d)),
                         max_1_plus_dbeta=float(np.max(1+d)),
                         seat_claimed_lower_bound=1-mub,
                         claim_holds=bool(np.min(1+d) >= 1-mub-1e-12),
                         monotone=bool(np.all(np.diff(img) >= -1e-15)),
                         image_min=float(img.min()), image_max=float(img.max()),
                         onto=bool(abs(img.min()) < 1e-12 and abs(img.max()-np.pi) < 1e-12),
                         sup_rel_disp=float(np.max(2*np.abs(np.sin(b/2)))),
                         mu_target=float(2*np.sin(mub/2))))
out = dict(rows=rows,
           any_claim_violation=bool(any(not r['claim_holds'] for r in rows)),
           all_monotone=bool(all(r['monotone'] for r in rows)),
           all_onto=bool(all(r['onto'] for r in rows)))
json.dump(out, open('r4_results.json','w'), indent=1)
hdr = ('%5s %8s %12s %12s %12s %9s %9s %6s %10s %8s' %
       ('mu_b','phi_c','min(1+bp)','max(1+bp)','seat >=','claim ok','monotone','onto','sup disp','mu'))
print(hdr)
for r in rows:
    print(f"{r['mu_beta']:5.2f} {r['phi_c']:8.0e} {r['min_1_plus_dbeta']:12.5g} {r['max_1_plus_dbeta']:12.5g} "
          f"{r['seat_claimed_lower_bound']:12.3f} {str(r['claim_holds']):>9s} {str(r['monotone']):>9s} "
          f"{str(r['onto']):>6s} {r['sup_rel_disp']:10.6f} {r['mu_target']:8.6f}")
print(f"\nseat's stated bound 1+beta' >= 1-mu_beta violated somewhere: {out['any_claim_violation']}")
print(f"map monotone everywhere: {out['all_monotone']};  onto [0,pi]: {out['all_onto']}")
