#!/usr/bin/env python3
"""Outward 1D range certificate for the full mixed coefficient E.

The trial pressure is the exact horizontal Green solution P(r)F(z),
including both radial tails. The whole-space error is bounded by an
angular-mode energy estimate. This is not a PDE time integration.
"""
import sys,json,argparse,time
from pathlib import Path
from fractions import Fraction as F
here=Path(__file__).resolve().parent
for p in [here.parent/'pass5',here.parents[2]/'outputs/euler-ns-transfer/pass5/certificate',here.parents[1]/'pass5/certificate']:
    if (p/'interval_local_integrals.py').exists():sys.path.insert(0,str(p));break
from interval_local_integrals import iv,interval,cutoff_jets,exact_endpoints

def mesh(points,n):
    result=[]; total=points[-1]-points[0]
    for a,b in zip(points[:-1],points[1:]):
        count=max(1,int(n*(b-a)/total)+1)
        result.extend((a+(b-a)*i/count,a+(b-a)*(i+1)/count) for i in range(count))
    return result

def pack(x):return dict(display=str(x),exact=exact_endpoints(x))

def run(n=2048,dps=35):
    start=time.time();iv.dps=dps
    m,k=4,20;a=interval(F(63,200));r0=interval(F(21,100));delta=interval(F(7,50))
    rpoints=list(map(F,['.07','.14','.1575','.21','.28','.315','.35']))
    cells=[]
    for lo,hi in mesh(rpoints,n):
        r=interval(lo,hi);dr=interval(hi-lo)
        e,e1=cutoff_jets(((r-r0)/delta)**2,1);er=2*(r-r0)/delta**2*e1
        C,C1=cutoff_jets((r/a)**2,1);Cp=2*r/a**2*C1
        co,si=iv.cos(k*r),iv.sin(k*r)
        gl=[r**m*Cp*e*co,r**m*Cp*e*si]
        gh=[r**(-m)*Cp*e*co,r**(-m)*Cp*e*si]
        cells.append(dict(r=r,dr=dr,e=e,er=er,C=C,co=co,si=si,gl=gl,gh=gh,lo=lo,hi=hi))
    # Prefix/suffix integrals enclose actual oscillatory Green moments.
    prefix=[[iv.mpf(0),iv.mpf(0)]]
    for c in cells:prefix.append([prefix[-1][j]+c['dr']*c['gl'][j] for j in (0,1)])
    suffix=[None]*(len(cells)+1);suffix[-1]=[iv.mpf(0),iv.mpf(0)]
    for i in range(len(cells)-1,-1,-1):
        c=cells[i];suffix[i]=[suffix[i+1][j]+c['dr']*c['gh'][j] for j in (0,1)]
    totals={key:iv.mpf(0) for key in ['Pnorm','trial_A','trial_B','d_A','d_B','local_R']}
    for i,c in enumerate(cells):
        r,dr,e,er,C=(c[j] for j in ['r','dr','e','er','C'])
        partial=interval(0,c['hi']-c['lo'])
        left=[prefix[i][j]+partial*c['gl'][j] for j in (0,1)]
        right=[suffix[i+1][j]+partial*c['gh'][j] for j in (0,1)]
        P=[-2*C*e*c[phase]-(m-1)*r**(-m)*left[j]-(m+1)*r**m*right[j] for j,phase in enumerate(['co','si'])]
        P2=P[0]**2+P[1]**2;Pa=iv.sqrt(P2)
        H=(er-e/r)**2+k*k*e*e
        vals=dict(Pnorm=r**3*P2,trial_A=r**2*Pa*iv.sqrt(H),trial_B=r*e*Pa,
                  d_A=r**5*H,d_B=r**3*e*e,local_R=r*r*e*e*C)
        for key,value in vals.items():totals[key]+=dr*value
    cin=[-(m+1)*v for v in suffix[0]];cout=[-(m-1)*v for v in prefix[-1]]
    inner=(cin[0]**2+cin[1]**2)*interval(rpoints[0])**(2*m+4)/(2*m+4)
    outer=(cout[0]**2+cout[1]**2)*interval(rpoints[-1])**(4-2*m)/(2*m-4)
    totals['Pnorm']+=inner+outer
    axial={key:iv.mpf(0) for key in ['Y1','Y2','Z1','Z2','Fpp']}
    zpoints=list(map(F,['3.4','3.55','3.7','3.775','4','4.225','4.3','4.45','4.6']))
    for lo,hi in mesh(zpoints,n):
        z=interval(lo,hi);dz=2*interval(hi-lo) # both symmetric packets
        y=(z-4)/interval(F(3,5));x=(z-4)/interval(F(9,20))
        q,q1,q2=cutoff_jets(y*y,2);v,v1,v2=cutoff_jets(x*x,2)
        qz=2*y*q1/interval(F(3,5));vz=2*x*v1/interval(F(9,20))
        qzz=(2*q1+4*y*y*q2)/interval(F(3,5))**2
        vzz=(2*v1+4*x*x*v2)/interval(F(9,20))**2
        Fpp=qzz*v+2*qz*vz+q*vzz
        vals=dict(Y1=q*q*v/z**7,Y2=abs(q*qz)*v/z**6,Z1=q*q/z**14,Z2=qz*qz/z**12,Fpp=Fpp*Fpp)
        for key,value in vals.items():axial[key]+=dz*value
    # Analytic kernel inequalities on r<=.35 and |z|>=3.4.
    # ||r d|| uses the triangle inequality between q and q_z terms.
    nd=m*iv.sqrt(iv.pi)*(45/(2*iv.pi)*iv.sqrt(totals['d_A']*axial['Z1'])+15/iv.pi*iv.sqrt(totals['d_B']*axial['Z2']))
    ne=iv.sqrt(iv.pi*totals['Pnorm']*axial['Fpp'])
    app=m*(45*totals['trial_A']*axial['Y1']+30*totals['trial_B']*axial['Y2'])
    tmax=interval(F(7,68)**2)
    K=(6-tmax)/(1+tmax)**5
    local_gain=interval(F(45,2)*m*k)*K*totals['local_R']*axial['Y1']
    error=interval(F(2,m*m))*nd*ne
    Eupper=app+error-local_gain
    threshold=interval(F(5093339,210000))
    # Preserve validity even if a sharper run eventually returns a negative bound.
    scaled=iv.mpf([max((780*Eupper).a,(1020*Eupper).a),max((780*Eupper).b,(1020*Eupper).b)])
    compatibility=4+scaled
    below_one_sixtieth=bool(Eupper.b<interval(F(1,60)).a)
    result=dict(status='PASS' if below_one_sixtieth and compatibility.b<threshold.a else 'INCONCLUSIVE',panels_requested=n,dps=dps,
        E_upper_strictly_below_one_sixtieth=below_one_sixtieth,
        radial={key:pack(val) for key,val in totals.items()},axial={key:pack(val) for key,val in axial.items()},
        inner_tail=pack(inner),outer_tail=pack(outer),d_norm_upper_expression=pack(nd),residual_norm_upper_expression=pack(ne),
        trial_pairing_absolute_upper_expression=pack(app),pressure_error_upper_expression=pack(error),
        local_negative_gain_lower_expression=pack(local_gain),E_upper_expression=pack(Eupper),
        compatibility_upper_expression=pack(compatibility),compatibility_threshold=str(F(5093339,210000)),
        seconds=time.time()-start,
        scope='Certified upper expressions from outward range integrals plus analytic inequalities; no two-sided E enclosure, no time evolution or return theorem.')
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--panels',type=int,default=2048);ap.add_argument('--dps',type=int,default=35);ap.add_argument('--output');args=ap.parse_args()
    result=run(args.panels,args.dps)
    if args.output:Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['status','panels_requested','d_norm_upper_expression','residual_norm_upper_expression','trial_pairing_absolute_upper_expression','pressure_error_upper_expression','local_negative_gain_lower_expression','E_upper_expression','compatibility_upper_expression','seconds']},indent=2),flush=True)
