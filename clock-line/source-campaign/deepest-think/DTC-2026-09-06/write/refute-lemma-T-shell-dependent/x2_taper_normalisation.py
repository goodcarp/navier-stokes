"""x2 - the target seat's 'unresolved cross-check' (its Corrections item 4 / gap 5) is a
NORMALISATION error, not a discrepancy.  prove-lagrangian s2_taper.py line 13 defines
    Phi_h(lam) := P_h(lam)/P_h(1) ,
and the SYNTHESIS line 47 tabulates Phi_{h_delta}(3/2).  The seat compared its UNNORMALISED
P_{h_delta}(3/2) against that ratio.  Divide by P_h(1) and the 3.4% 'gap' at 30 deg closes.
"""
import json, mpmath as mp
mp.mp.dps = 30
def g2(p,l): return l**2*mp.sin(p)**2 + l**-4*mp.cos(p)**2
def P_h(l,deg):
    l = mp.mpf(l); d = mp.mpf(deg)*mp.pi/180
    f = lambda p: 3*min(mp.mpf(1), p/d)*mp.cos(p)*mp.sin(p)**2/g2(p,l)**mp.mpf(2.5)
    return mp.quad(f, [0, min(d, mp.pi/2), mp.pi/2])
SYNTH = {'3':1.4980,'5':1.4914,'7.5':1.4736,'10':1.4444,'15':1.3597,'20':1.2584,'30':1.0674}
SEAT  = {'3':1.497990,'5':1.491146,'7.5':1.472733,'10':1.442517,'15':1.353721,'20':1.245539,'30':1.032403}
rows=[]
for k,synth in SYNTH.items():
    P1 = P_h(1.0,float(k)); P15 = P_h(1.5,float(k))
    ratio = float(P15/P1)
    rows.append(dict(delta_deg=k, P_h_1=float(P1), P_h_1p5=float(P15),
                     seat_P_h_1p5=SEAT[k], my_minus_seat=abs(float(P15)-SEAT[k]),
                     ratio_Phi_h=ratio, synthesis_Phi_h=synth,
                     resid_unnormalised=abs(float(P15)-synth)/synth,
                     resid_normalised=abs(ratio-synth)/synth))
out=dict(rows=rows,
         max_resid_unnormalised=max(r['resid_unnormalised'] for r in rows),
         max_resid_normalised=max(r['resid_normalised'] for r in rows),
         max_my_minus_seat=max(r['my_minus_seat'] for r in rows))
json.dump(out, open('x2_results.json','w'), indent=1)
print("delta   P_h(1)      P_h(3/2)   seat P_h(3/2)  ratio=Phi_h  synthesis  |unnorm|   |norm|")
for r in rows:
    print("%5s  %.7f  %.7f  %.7f   %.6f   %.4f   %.2e  %.2e" %
          (r['delta_deg'], r['P_h_1'], r['P_h_1p5'], r['seat_P_h_1p5'],
           r['ratio_Phi_h'], r['synthesis_Phi_h'], r['resid_unnormalised'], r['resid_normalised']))
print("\nmax relative residual  UNNORMALISED (the seat's comparison): %.4e" % out['max_resid_unnormalised'])
print("max relative residual  NORMALISED   (the correct comparison): %.4e" % out['max_resid_normalised'])
print("my P_h(3/2) vs the seat's, max abs diff: %.2e  (the seat's own numbers are right)"
      % out['max_my_minus_seat'])
