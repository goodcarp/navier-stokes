#!/usr/bin/env python3
"""New leading-envelope initial Gtt enclosure with a new full-pressure error.

No NS time integration. Radial Green moments/norms use outward range sums,
both exact tails, and the unchanged axial F'' interval from pass7.
"""
import argparse,json,os,sys,time
from pathlib import Path
from fractions import Fraction as F
here=Path(__file__).resolve().parent
def locate(candidates,name):
    for candidate in candidates:
        if (candidate/name).is_file():
            return candidate
    raise FileNotFoundError("Missing sibling dependency: "+name)
p5dir=locate([here.parent/"pass5",here.parents[1]/"pass5"/"certificate"],
             "interval_local_integrals.py")
p8dir=locate([here.parent/"pass8",here.parents[1]/"pass8"/"checks"],
             "certify_fixed_circle_acceleration.py")
p7dir=locate([here.parent/"pass7",here.parents[1]/"pass7"/"checks"],
             "outer-E-certificate.json")
sys.path.insert(0,str(p5dir))
from interval_local_integrals import iv,interval,cutoff_jets,exact_endpoints
sys.path.insert(0,str(p8dir))
from certify_fixed_circle_acceleration import gjets

def endpoint(x):
    return (-1 if x["sign"] else 1)*F(int(x["mantissa"]))*F(2)**x["exponent"]
def recover(x):
    return interval(endpoint(x["exact"]["lower"]),endpoint(x["exact"]["upper"]))
def pack(x):
    return dict(display=str(x),exact=exact_endpoints(x))
def mesh(points,n):
    length=points[-1]-points[0]
    for a,b in zip(points[:-1],points[1:]):
        count=max(1,int(n*(b-a)/length)+1)
        for j in range(count):
            yield a+(b-a)*j/count,a+(b-a)*(j+1)/count

def seed_jets_at_root(r):
    width=interval(F(1,10))
    y=(r-interval(F(7,50)))/width
    x=y*y
    p,d1,d2=cutoff_jets(x,2)
    t=(4*x-1)/3
    assert t.a>0 and t.b<1
    z1=(-t**-2-(1-t)**-2)*interval(F(4,3))
    z2=(2*t**-3-2*(1-t)**-3)*interval(F(4,3))**2
    z3=(-6*t**-4-6*(1-t)**-4)*interval(F(4,3))**3
    d3=p*(1-p)*((1-6*p+6*p*p)*z1**3+3*(1-2*p)*z1*z2+z3)
    return (p,2*y*d1/width,(2*d1+4*y*y*d2)/width**2,
            (12*y*d2+8*y**3*d3)/width**3)

def run(n=2048):
    started=time.time();iv.dps=40
    old=json.loads((p8dir/"fixed-circle-acceleration-certificate.json").read_text())
    lead=json.loads((p8dir/"leading-envelope-certificate.json").read_text())
    assert old["status"]==lead["status"]=="PASS"
    assert lead["parameters"]["k"]==-20 and lead["parameters"]["m"]==4
    assert lead["parameters"]["radial_center"]=="7/50"
    assert lead["parameters"]["radial_width"]=="1/10"
    oldE=json.loads((p7dir/"outer-E-certificate.json").read_text())
    # Axial profiles are unchanged, but the radial norm below is NEW.
    Fpp=recover(oldE["axial"]["Fpp"])
    rlo,rhi=map(F,old["unique_global_radial_maximum_interval"])
    rpoints=[F(1,25),F(13,200),F(9,100),F(7,50),F(63,400),
             F(19,100),rlo,rhi,F(6,25)]
    av=interval(F(63,200));width=interval(F(1,10))
    m,k=4,-20
    cells=[]
    for lo,hi in mesh(rpoints,n):
        r=interval(lo,hi);dr=interval(hi-lo)
        e=cutoff_jets(((r-interval(F(7,50)))/width)**2,0)[0]
        C,C1=cutoff_jets((r/av)**2,1)
        Cp=2*r*C1/av**2
        co,si=iv.cos(k*r),iv.sin(k*r)
        gl=[r**m*Cp*e*co,r**m*Cp*e*si]
        gr=[r**(-m)*Cp*e*co,r**(-m)*Cp*e*si]
        cells.append(dict(lo=lo,hi=hi,r=r,dr=dr,e=e,C=C,co=co,si=si,gl=gl,gr=gr))
    prefix=[[iv.mpf(0),iv.mpf(0)]]
    for cell in cells:
        prefix.append([prefix[-1][j]+cell["dr"]*cell["gl"][j] for j in (0,1)])
    suffix=[None]*(len(cells)+1);suffix[-1]=[iv.mpf(0),iv.mpf(0)]
    for i in range(len(cells)-1,-1,-1):
        cell=cells[i]
        suffix[i]=[suffix[i+1][j]+cell["dr"]*cell["gr"][j] for j in (0,1)]
    radial_norm=iv.mpf(0)
    for i,cell in enumerate(cells):
        r=cell["r"];e=cell["e"];C=cell["C"]
        part=interval(0,cell["hi"]-cell["lo"])
        left=[prefix[i][j]+part*cell["gl"][j] for j in (0,1)]
        right=[suffix[i+1][j]+part*cell["gr"][j] for j in (0,1)]
        pp=[-2*C*e*cell[phase]-(m-1)*r**(-m)*left[j]-
            (m+1)*r**m*right[j] for j,phase in enumerate(["co","si"])]
        radial_norm+=cell["dr"]*r**3*(pp[0]**2+pp[1]**2)
    cin=[-(m+1)*v for v in suffix[0]]
    cout=[-(m-1)*v for v in prefix[-1]]
    inner=(cin[0]**2+cin[1]**2)*interval(rpoints[0])**(2*m+4)/(2*m+4)
    outer=(cout[0]**2+cout[1]**2)*interval(rpoints[-1])**(4-2*m)/(2*m-4)
    radial_norm+=inner+outer
    assert radial_norm.a>0 and Fpp.a>0

    # Enclose the Green derivative at the exact critical radius.
    rootcell=[i for i,cell in enumerate(cells) if cell["lo"]==rlo and cell["hi"]==rhi]
    assert len(rootcell)==1
    i=rootcell[0];r=interval(rlo,rhi)
    e,er,err,errr=seed_jets_at_root(r)
    C,C1=cutoff_jets((r/av)**2,1);Cp=2*r*C1/av**2
    co,si=iv.cos(k*r),iv.sin(k*r);part=interval(0,rhi-rlo)
    left=[prefix[i][j]+part*r**m*Cp*e*trig for j,trig in enumerate([co,si])]
    right=[suffix[i+1][j]+part*r**(-m)*Cp*e*trig for j,trig in enumerate([co,si])]
    Pr=-2*C*(er*co-k*e*si)+m*(m-1)*r**(-m-1)*left[0]-m*(m+1)*r**(m-1)*right[0]
    Pi=-2*C*(er*si+k*e*co)+m*(m-1)*r**(-m-1)*left[1]-m*(m+1)*r**(m-1)*right[1]

    # NEW full-pressure error from harmonic slab, m=4, h=9/40.
    # pi cancels between ne²=pi*radial_norm*Fpp and the slab evaluation.
    slab_factor=interval(F(35,4))*r**6/interval(F(9,20))**9*interval(F(10001,10000))
    error_pressure=iv.sqrt(radial_norm*Fpp*slab_factor/(m*m))
    residual_norm=iv.sqrt(iv.pi*radial_norm*Fpp)

    G=gjets(r,4)
    c,nu=interval(F(7,5)),interval(F(1,1000))
    L2=c*c*r*r*G[2]-2*c*nu*r*G[3]+nu*nu*(G[4]-2*G[3]/r+3*G[2]/r**2)
    Dre=(k*k*r+m*m/r)*e-r*err-er
    Dim=k*(e+2*r*er)
    product=(co*Dre+si*Dim)*Pr-(co*Dim-si*Dre)*Pi
    local=m*m*e*e*G[2]/(2*r*r)-C*(er*er+k*k*e*e+r*er*(err+k*k*e)-m*m*e*er/r)
    S_trial=local+product/2
    S_error=iv.sqrt(Dre*Dre+Dim*Dim)*error_pressure/2
    S=S_trial+iv.mpf([-S_error.b,S_error.b])
    F1=2*e*er;F2=2*(er*er+e*err)
    Jnu=(6*(er*err+e*errr)+6*e*err/r-4*k*k*e*er-2*k*k*e*e/r-
         (4*m*m+8)*e*er/r**2+(2*m*m+4)*e*e/r**3)
    pump=-2*c*(2*F1+r*F2)
    non_A=interval(F(m*k,2))*(pump+nu*Jnu)
    amp=interval(900,1020);lam2=interval(F(1,64),F(1,16))
    Gtt=amp*L2+lam2*amp*S+lam2*non_A
    Gtt_trial=amp*L2+lam2*amp*S_trial+lam2*non_A
    Tunit=interval(F(m*k,2))*(F1+e*e/r)
    Tprime=interval(F(m*k,2))*(F2+F1/r-e*e/r**2)
    Gt=nu*amp*G[2]+lam2*Tunit
    Gtr=amp*(-c*r*G[2]+nu*(G[3]-G[2]/r))+lam2*Tprime
    radial_correction=Gtr**2/(-amp*G[2])
    tracked_Gtt=Gtt+radial_correction
    def exact_lo(x):return endpoint(exact_endpoints(x)["lower"])
    def exact_hi(x):return endpoint(exact_endpoints(x)["upper"])
    predicates=dict(Gt_above_11=exact_lo(Gt)>11,
        Gtt_below_minus_13000=exact_hi(Gtt)<-13000,
        tracked_Gtt_below_minus_12000=exact_hi(tracked_Gtt)<-12000,
        NEW_pressure_error_expression_below_13_over_10=exact_hi(error_pressure)<F(13,10))
    result=dict(status="PASS" if all(predicates.values()) else "INCONCLUSIVE",
        exact_rational_predicates=predicates,panels=n,dps=iv.dps,
        dependencies={name:os.path.relpath(path,here) for name,path in {
            "interval_engine":p5dir/"interval_local_integrals.py",
            "unit_mean_root_and_jets":p8dir/"fixed-circle-acceleration-certificate.json",
            "unit_mean_jet_implementation":p8dir/"certify_fixed_circle_acceleration.py",
            "leading_family":p8dir/"leading-envelope-certificate.json",
            "unchanged_axial_Fpp":p7dir/"outer-E-certificate.json"}.items()},
        parameters=lead["parameters"],root_interval=[str(rlo),str(rhi)],
        radial_pressure_weighted_square=pack(radial_norm),
        inner_pressure_tail=pack(inner),outer_pressure_tail=pack(outer),
        unchanged_axial_Fpp_square=pack(Fpp),
        NEW_weighted_residual_norm=pack(residual_norm),
        NEW_full_pressure_radial_error_bound=pack(error_pressure),
        Green_gradient_real=pack(Pr),Green_gradient_imag=pack(Pi),
        seed_jets=[pack(x) for x in [e,er,err,errr]],
        unit_mean_jets=[pack(x) for x in G],
        unit_transport_diffusion=pack(L2),
        S_local=pack(local),S_pressure_trial=pack(product/2),
        S_trial=pack(S_trial),S_pressure_error=pack(S_error),S_actual=pack(S),
        non_A_pump_viscosity_coefficient=pack(non_A),
        Gtt_horizontal_trial_diagnostic=pack(Gtt_trial),
        Gtt_actual_full_pressure=pack(Gtt),Gt_actual=pack(Gt),
        moving_radial_correction=pack(radial_correction),
        moving_radial_value_Gtt_actual=pack(tracked_Gtt),
        pressure_error_interval_semantics="The pressure-error interval encloses an upper-bound expression. Only its upper endpoint bounds the actual complex pressure-gradient error; its lower endpoint is not a lower bound on that error.",
        seconds=time.time()-started,
        scope="Actual initial critical-circle jet enclosure, with a newly recomputed leading-envelope pressure residual and both radial tails. Independent A/lambda intervals enlarge the result. No finite-time NS evolution or global-maximum acceleration assertion.")
    return result

if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--panels",type=int,default=2048);ap.add_argument("--output")
    args=ap.parse_args();result=run(args.panels)
    if args.output:Path(args.output).write_text(json.dumps(result,indent=2)+"\n")
    summary={k:result[k] for k in ["status","panels","NEW_weighted_residual_norm",
      "NEW_full_pressure_radial_error_bound","S_local","S_pressure_trial","S_actual",
      "Gtt_actual_full_pressure","Gt_actual","moving_radial_value_Gtt_actual","seconds"]}
    print(json.dumps(summary,indent=2),flush=True)
