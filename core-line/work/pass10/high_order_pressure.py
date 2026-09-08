"""Sixth-order radial collocation diagnostic, not a compatible NS discretization.

Every angular pressure mode here is even. Reflect the scalar stencil across
the axis, solve the exact cylindrical ODE, and impose the decaying Bessel
Robin condition at the last node (outside the compact initial sources).
The zero-zero mode uses an outer gauge and reports the computed outer flux;
it does not silently repair source compatibility. The associated numerical
velocity is not asserted exactly divergence free.
"""
import numpy as np
from scipy.fft import fft,ifft
from scipy.linalg import solve_banded
from scipy.special import kve
from scipy.sparse import csr_matrix
from full_initial_evolution import Cylinder,differentiation


def reflected_derivative(x,order,width=7,parity=1):
    n=len(x);half=width//2
    rows=[];cols=[];vals=[]
    for i in range(n):
        start=min(i-half,n-width)
        js=np.arange(start,start+width)
        indices=np.where(js<0,-js-1,js)
        nodes=np.where(js<0,-x[indices],x[indices])
        scale=max(abs(nodes-x[i]));yy=(nodes-x[i])/scale
        rhs=np.zeros(width);rhs[order]=1 if order==1 else 2
        w=np.linalg.solve(np.array([yy**k for k in range(width)]),rhs)/scale**order
        w*=np.where(js<0,parity,1)
        rows.extend([i]*width);cols.extend(indices);vals.extend(w)
    return csr_matrix((vals,(rows,cols)),shape=(n,n))


class HighOrderCylinder(Cylinder):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.Dr=differentiation(self.r,1,width=7)
        self.Drr=differentiation(self.r,2,width=7)
        self.eDr=reflected_derivative(self.r,1)
        self.eDrr=reflected_derivative(self.r,2)
        self.oDr=reflected_derivative(self.r,1,parity=-1)
        self.oDrr=reflected_derivative(self.r,2,parity=-1)

    def pressure_gradient(self,p,m):
        assert m%2==0
        return np.stack((self.eDr@p,1j*m*p/self.r[:,None],self.dz(p)))

    def gradient(self,u,m):
        assert m%2==0
        rr=self.r[:,None]
        return np.stack((np.stack((self.oDr@u[0],self.oDr@u[1],self.eDr@u[2])),
            np.stack(((1j*m*u[0]-u[1])/rr,(1j*m*u[1]+u[0])/rr,1j*m*u[2]/rr)),
            np.stack([self.dz(v) for v in u])),axis=1)

    def vector_lap(self,u,m):
        rr=self.r[:,None]
        out=np.stack((self.oDrr@u[0]+(self.oDr@u[0])/rr,
                      self.oDrr@u[1]+(self.oDr@u[1])/rr,
                      self.eDrr@u[2]+(self.eDr@u[2])/rr))
        out+=self.dz(u,2)-m*m*u/rr**2
        out[0]-=(u[0]+2j*m*u[1])/rr**2
        out[1]-=(u[1]-2j*m*u[0])/rr**2
        return out

    def poisson(self,source,m,tag=''):
        assert m%2==0
        n=self.nr;band=6
        op=-self.eDrr-self.eDr.multiply(1/self.r[:,None])
        op=op.tolil();op[-1,:]=self.eDr[-1,:]
        op=op.tocoo()
        ab=np.zeros((2*band+1,n))
        for i,j,a in zip(op.row,op.col,op.data):ab[band+i-j,j]+=a
        ab[band,:-1]+=m*m/self.r[:-1]**2
        rhs=fft(source,axis=1,norm='forward');out=np.empty_like(rhs)
        radius=self.r[-1]
        for j,k in enumerate(abs(self.kz)):
            a=ab.copy();v=rhs[:,j].copy();v[-1]=0
            a[band,:-1]+=k*k
            if k==0 and m==0:
                for col in range(max(0,n-1-band),n):a[band+n-1-col,col]=0
                a[band,-1]=1
            else:
                kap=k*kve(abs(m-1),k*radius)/kve(m,k*radius)+m/radius if k else m/radius
                a[band,-1]+=kap
            out[:,j]=solve_banded((band,band),a,v,check_finite=False)
        if m==0:
            moment=np.sum(self.vol*rhs[:,0])
            flux=radius*complex((self.eDr[-1]@out[:,0])[0])
            self.repairs.append(dict(tag=tag,source_repaired=False,
                raw_midpoint_mean_source_radial_integral=[float(moment.real),float(moment.imag)],
                computed_outer_zero_mode_radial_flux=[float(flux.real),float(flux.imag)]))
        return ifft(out,axis=1,norm='forward')
