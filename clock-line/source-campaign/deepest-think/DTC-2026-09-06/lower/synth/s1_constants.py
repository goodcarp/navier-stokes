import math, json
out={}
# clock constants
out['4log(3/2)']=4*math.log(1.5)
out['8(1-sqrt(2/3))']=8*(1-math.sqrt(2/3))
out['ratio_accel_over_frozen']=2*(1-math.sqrt(2/3))/math.log(1.5)
out['theta_frozen']=math.log(1.5)
out['theta_accel']=2*(1-math.sqrt(2/3))
out['theta_riccati']=0.5
# c2 = 2*theta/kappa with kappa=1/2
for name,th in [('frozen',math.log(1.5)),('accel',2*(1-math.sqrt(2/3))),('riccati',0.5)]:
    out['c2_'+name+'_kappa_half']=2*th/0.5
# c2 at the datum's measured kappa values (riccati theta=1/2)
for k in [0.5,0.49972,0.48137,0.45317,0.437,0.415,0.411,0.400]:
    out['c2_riccati_kappa_%.5f'%k]=2*0.5/k
    out['Qprime_frozen_kappa_%.5f'%k]=math.log(1.5)/k
    out['Qprime_accel_kappa_%.5f'%k]=2*(1-math.sqrt(2/3))/k
# z-scores vs viscous measurement Q'=0.988 +- 0.018
for k in [0.437,0.415,0.411]:
    out['z_frozen_k%.3f'%k]=(math.log(1.5)/k-0.988)/0.018
    out['z_accel_k%.3f'%k]=(2*(1-math.sqrt(2/3))/k-0.988)/0.018
# advertised gain vs shortfall
out['advertised_gain_pct']=100*(1-8*(1-math.sqrt(2/3))/(4*math.log(1.5)))
out['shortfall_pct_k0.411']=100*(1-(2*(1-math.sqrt(2/3))/0.411)/0.988)
# log Re_E bookkeeping
out['c_E_sharp']=(2/5)*math.log(0.172403978)
out['c_E_mollified']=(2/5)*math.log(0.1418961)
out['c_E_mismatch']=out['c_E_sharp']-out['c_E_mollified']
# corollary: effective loglog / log coefficients
def kappa_eff(s,K=1.0,c1=1.0):
    m=c1/(s*(1+max(0.0,math.log(K*s**(-0.2)))))
    return m*s*math.log(1/s)/c1
for e in [6,12,30,100]:
    out['kappa_eff_1e-%d'%e]=kappa_eff(10.0**-e)
print(json.dumps(out,indent=1,sort_keys=True))
json.dump(out,open('s1_results.json','w'),indent=1,sort_keys=True)
