#!/usr/bin/env python3
"""x7 -- the c2 inflation that follows from the CORRECTED inset, so the headline
comparison is apples-to-apples.  c2 = 8(1-sqrt(2/3)) * L/(L - log(1+f)) (the seat's and
gap-V's own formula)."""
import json, math
R={}
c2inf = lambda L,f: 8*(1-math.sqrt(2/3))*L/(L-math.log(1+f))
seat_f={10:0.25972,20:0.19228,40:0.14189,80:0.10441,160:0.07663}     # PROOF.md sec 5
gapV_f={10:0.9131,20:0.6590,40:0.4694,80:0.3382,160:0.2437}          # gap-V sec 5(ii)
thm_f ={10:0.5076972822670149,20:0.37636092251960546,40:0.2752757189297324,
        80:0.20323505876469117,160:0.1500621905476369}               # x4, K2hat=0, vs 1/L
thmK1 ={10:3.590,20:3.590,40:3.590,80:3.590,160:3.590}               # x5, K2hat=1 optimum
print("  L    c2(seat f)  c2(theorem f, K2hat=0)  c2(theorem f, K2hat=1)  c2(gap-V f)   asymptote")
rows=[]
for L in [10,20,40,80,160]:
    a,b,c,d = c2inf(L,seat_f[L]), c2inf(L,thm_f[L]), c2inf(L,thmK1[L]), c2inf(L,gapV_f[L])
    rows.append(dict(L=L,c2_seat=a,c2_theorem_K0=b,c2_theorem_K1=c,c2_gapV=d))
    print(f"{L:5d}   {a:.5f}      {b:.5f}                 {c:.5f}                 {d:.5f}     {8*(1-math.sqrt(2/3)):.7f}")
R['c2_table']=rows
R['asymptote']=8*(1-math.sqrt(2/3))
json.dump(R,open('x7_results.json','w'),indent=1)
print("WROTE x7_results.json")
