#!/usr/bin/env python3
"""Exploratory radial-annular pump: exact pressure functional cutoff at l=4.

Imports the audited pass4 evaluator without editing that archived package.
Only initial derivatives are evaluated; quadrature is not interval certified.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import sys,json,argparse,time
from pathlib import Path
here=Path(__file__).resolve().parent
for parent in (here.parent,here.parent.parent):
    for path in (parent/'pass4',parent/'pass4'/'experiments'):
        if (path/'evaluate_pressure_gate.py').exists(): sys.path.insert(0,str(path))
from evaluate_pressure_gate import Evaluator,cutoff,profile_functions
import numpy as np
import sympy as sp
from scipy.special import eval_legendre


def radial_pump_function(inner=2.5,plateau_end=5.,outer=6.):
    r,z=sp.symbols('r z',positive=True)
    f=sp.Function('f'); s=r*r+z*z
    eta=-(1-f(s/inner**2))*f(sp.Rational(1,4)+sp.Rational(3,4)*(s-plateau_end**2)/(outer**2-plateau_end**2))
    field=(-r*(eta+z*sp.diff(eta,z)),sp.Integer(0),2*z*eta+r*z*sp.diff(eta,r))
    modules={'f':lambda x:cutoff(x,0)}
    for n in range(1,4): modules['f'+str(n)]=lambda x,n=n:cutoff(x,n)
    def lower(expr):
        replacements={}
        for atom in expr.atoms(sp.Subs):
            if isinstance(atom.expr,sp.Derivative) and atom.expr.expr.func==f:
                n=sum(k for _,k in atom.expr.variable_count)
                replacements[atom]=sp.Function('f'+str(n))(atom.point[0])
        return expr.xreplace(replacements)
    exprs=[]
    for u in field:
        exprs.extend([u,sp.diff(u,r),sp.diff(u,z),sp.diff(u,r,2),sp.diff(u,r,z),sp.diff(u,z,2)])
    return sp.lambdify((r,z),[lower(e) for e in exprs],modules=[modules,'numpy'],cse=True)


class RadialEvaluator(Evaluator):
    def __init__(self,nr=1200,nmu=1024,lmax=4,inner=2.5,plateau_end=5.,outer=6.):
        start=time.time()
        if inner/2<=1 or inner>=3.55 or plateau_end<=4.47 or outer<=plateau_end:
            raise ValueError('Pump must avoid core and equal -S0x on all swirl support.')
        if lmax<4 or lmax%2: raise ValueError('Use an even cutoff at least four.')
        self.R=(np.arange(nr)+.5)*(outer+.25)/nr
        self.mu,self.weights=np.polynomial.legendre.leggauss(nmu)
        self.RR=self.R[:,None]; self.MM=self.mu[None,:]
        self.sn=np.sqrt(1-self.MM**2)
        self.r=self.RR*self.sn; self.z=self.RR*self.MM
        self.ells=np.arange(0,lmax+1,2)
        self.PL=np.array([eval_legendre(int(l),self.mu) for l in self.ells])
        self.PL1=np.array([np.zeros_like(self.mu) if l==0 else l*(self.mu*self.PL[j]-eval_legendre(int(l-1),self.mu))/(self.mu**2-1) for j,l in enumerate(self.ells)])
        self.PL2=(2*self.mu*self.PL1-self.ells[:,None]*(self.ells[:,None]+1)*self.PL)/(1-self.mu**2)
        self.project=((2*self.ells+1)[:,None]*self.weights[None,:]*self.PL/2).T
        funcs=profile_functions(.7,1.,4.,.45)
        funcs[2]=radial_pump_function(inner,plateau_end,outer)
        self.fields=np.array([np.array([np.broadcast_to(a,self.r.shape) for a in fn(self.r,self.z)]).reshape(3,6,nr,nmu) for fn in funcs])
        self.metadata=dict(nr=nr,nmu=nmu,lmax=lmax,pump='radial annular',inner_transition_end=inner,plateau_end=plateau_end,outer_support=outer)
        self.pressures=[self.pressure_source_hzz(self.source(self.combine(np.eye(4)[j])),j==0,j==1) for j in range(4)]
        # These two values vanish analytically by radial-profile pressure
        # invariance; keep their quadrature values as diagnostics only.
        self.Cp_numeric=self.outer_stress(2)
        self.dP_numeric=self.outer_stress(2,gradient=True)
        self.Cp=0.; self.dP=0.
        self.Cv=self.outer_stress(3)
        self.dv=self.outer_stress(3,gradient=True)
        self.setup_seconds=time.time()-start

    def pressure_hessian(self,g):
        result=super().pressure_hessian(g)
        self.last_pr=result[-2]; self.last_pz=result[-1]
        return result

    def evaluate(self,**kwargs):
        row=super().evaluate(**kwargs)
        c=row['c']; A=row['A']
        uouter=c*self.fields[2,:,0]+A*self.fields[3,:,0]
        gradp=np.array([self.last_pr,np.zeros_like(self.last_pr),self.last_pz])
        pressure_term=2*self.integrate_outer_tensor(np.einsum('irm,jrm->ijrm',uouter,gradp))
        row.update(outer_pressure_part=float(pressure_term),outer_local_part=row['outer_part']-float(pressure_term),Cp_numeric=self.Cp_numeric,dP_numeric=self.dP_numeric,
                   angular_scope='Both scalar functionals have exact l<=4 pressure cutoff for this radial meridional class. They may disagree numerically through retained-mode and local quadrature errors.')
        return row


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--nr',type=int,default=1200)
    ap.add_argument('--nmu',type=int,default=1024)
    ap.add_argument('--lmax',type=int,default=4)
    ap.add_argument('--c',type=float,nargs='+',default=[0.,1.,2.])
    args=vars(ap.parse_args()); cs=args.pop('c')
    ev=RadialEvaluator(**args)
    print(json.dumps(dict(metadata=ev.metadata,pressures=ev.pressures,Cv=ev.Cv,dv=ev.dv,Cp_numeric=ev.Cp_numeric,dP_numeric=ev.dP_numeric,setup_seconds=ev.setup_seconds)),flush=True)
    for c in cs: print(json.dumps(ev.evaluate(b=1.,omega=1.,c=c)),flush=True)
