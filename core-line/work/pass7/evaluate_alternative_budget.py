"""Compact-support quadrature of alternative E bounds; exploratory only."""
import sys,json
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from evaluate_outer_E import cutoff,fields

def evaluate(nr=500,nz=1000):
    dr=.28/nr; dz=1.2/nz
    r=(.07+(np.arange(nr)+.5)*dr)[:,None]
    z=(3.4+(np.arange(nz)+.5)*dz)[None,:]
    m,k=4.,20.
    x=(r-.21)**2/.14**2
    e=cutoff(x); er=cutoff(x,1)*2*(r-.21)/.14**2
    q=cutoff((z-4)**2/.6**2)
    qz=cutoff((z-4)**2/.6**2,1)*2*(z-4)/.6**2
    a=e*q; ar=er*q; az=e*qz
    vz=cutoff((z-4)**2/.45**2)
    vzz=cutoff((z-4)**2/.45**2,1)*2*(z-4)/.45**2
    B=cutoff(r*r/.315**2)*vz
    Br=cutoff(r*r/.315**2,1)*2*r/.315**2*vz
    Bz=cutoff(r*r/.315**2)*vzz
    Vr=B+r*Br
    R2=r*r+z*z
    Qt=3*(4*z*z-r*r)/(4*np.pi*R2**3.5)
    Qtr=15*r*(r*r-6*z*z)/(4*np.pi*R2**4.5)
    Qrr=Qt+r*Qtr
    Qrz=-15*r*z*(4*z*z-3*r*r)/(4*np.pi*R2**4.5)
    loc,d,_=fields(r,z)
    integ=lambda f: float(2*np.pi*np.sum(r*f)*dr*dz)
    energy=integ(m*m*a*a/r**2+ar*ar+k*k*a*a)
    Cw=integ(Qrr*m*m*a*a/r**2+Qt*(ar*ar+k*k*a*a))
    nQ=np.sqrt(integ((Qrr*Qrr+Qrz*Qrz)*m*m*a*a/r**2+Qt*Qt*(ar*ar+k*k*a*a)))
    nG=np.sqrt(integ(B*B*(ar*ar+k*k*a*a)+Vr*Vr*m*m*a*a/r**2))
    nH=np.sqrt(integ((1+m*m)*Br*Br*a*a+(Bz*a+B*az)**2))
    nd=np.sqrt(integ(r*r*abs(d)**2))
    Eloc=float(2*loc.sum()*dr*dz)
    Eapp=float(-8*np.pi*m*k*np.sum((Qrr-Qt)*B*a*a)*dr*dz)
    return dict(nr=nr,nz=nz,Cw=Cw,seed_energy=energy,E_local=Eloc,
                norm_Qw=nQ,norm_r_d=nd,norm_w_grad_v=nG,norm_H=nH,
                Epressure_bound_direct=4*nQ*nG,
                Epressure_bound_angular_projection=4*nd*nG/m,
                Epressure_app=Eapp,Epressure_app_remainder_bound=4*nd*nH/m,
                status='EXPLORATORY midpoint quadrature; no rigorous enclosure')

if __name__=='__main__':
    print(json.dumps(evaluate()),flush=True)
