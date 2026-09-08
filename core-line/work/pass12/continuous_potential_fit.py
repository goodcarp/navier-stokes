#!/usr/bin/env python3
"""Continuous L2 potential projection of a specified MAC interpolant.

The target is the parity-preserving cubic radial interpolant of the full
MAC increment, retaining axial Fourier phases. This changes the fitting
metric, not the fluid data or the actual NS equation. It is still an
unvalidated whole-space H7 reconstruction.
"""
import argparse,json,time,hashlib
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.sparse import csr_matrix
from scipy.linalg import cholesky_banded,cho_solve_banded
from scipy.fft import fft,ifft
import compact_potential_reconstruction as base
from observe_compact_reconstruction import gauss_spans


def interpolation_matrix(nodes,target,R,odd=True):
    """Cubic interpolant with radial parity, zero at R; no extra axis jet fit."""
    nodes=np.asarray(nodes);n=len(nodes)
    x=np.r_[-R,-nodes[::-1],nodes,R]
    I=np.eye(n)
    y=np.vstack((np.zeros((1,n)),(-1 if odd else 1)*I[::-1],I,np.zeros((1,n))))
    return CubicSpline(x,y,bc_type='not-a-knot',extrapolate=False)(target)


def fit_mode(g,delta_hat,j,intervals=96):
    m=int(g.modes[j]);basis=base.RadialBasis(g.R,intervals,g.mapping,m)
    breaks=np.unique(np.r_[np.sqrt(np.unique(basis.knots)),g.rf[1:-1],g.r])
    r,w=gauss_spans(breaks,m+18)
    B,Dr,Lap,mB=basis.matrices(r);wr=w*r
    Dsp,Lsp,Msp=map(csr_matrix,(Dr,Lap,mB))
    # The continuous toroidal/poloidal cross term is zero by horizontal IBP.
    H=(Dsp.T@Dsp.multiply(wr[:,None])+Msp.T@Msp.multiply(wr[:,None])).toarray()
    Q=(Lsp.T@Lsp.multiply(wr[:,None])).toarray()
    Sr=interpolation_matrix(g.rf[1:-1],r,g.R,odd=True)
    St=interpolation_matrix(g.r,r,g.R,odd=True)
    Sz=interpolation_matrix(g.r,r,g.R,odd=False)
    dR=Dsp.T@(wr[:,None]*Sr);mR=Msp.T@(wr[:,None]*Sr)
    mT=Msp.T@(wr[:,None]*St);dT=Dsp.T@(wr[:,None]*St)
    lZ=Lsp.T@(wr[:,None]*Sz)
    Ur,Ut,Uz=delta_hat[g.slr],delta_hat[g.slt],delta_hat[g.slz]
    rhsP=-1j*g.kz*(dR@Ur)-g.kz*(mT@Ut)-lZ@Uz
    rhsT=-1j*(mR@Ur)-dT@Ut
    hb=base.lower_banded(H,9);qb=base.lower_banded(Q,9)
    P=np.empty_like(rhsP);T=np.empty_like(rhsT);n=basis.n;backward=0.
    def solve(band,rhs):
        scale=1/np.sqrt(band[0]);eq=band.copy()
        for d in range(len(eq)):
            eq[d,:n-d]*=scale[d:]*scale[:n-d]
        fac=cholesky_banded(eq,lower=True,check_finite=False)
        return scale[:,None]*cho_solve_banded((fac,True),scale[:,None]*rhs,check_finite=False)
    T[:]=solve(hb,rhsT)
    for col,k in enumerate(g.kz):
        P[:,col]=solve(qb+k*k*hb,rhsP[:,col:col+1])[:,0]
        if col in (0,1,g.nz//2,g.nz-1):
            residual=(Q+k*k*H)@P[:,col]-rhsP[:,col]
            backward=max(backward,float(np.linalg.norm(residual)/(np.linalg.norm(rhsP[:,col])+1.e-300)))
    Bc,Dc,Lc,mc=basis.matrices(g.r);Bf,Df,Lf,mf=basis.matrices(g.rf[1:-1])
    fitted=np.concatenate((1j*g.kz*(Df@P)+1j*mf@T,-g.kz*(mc@P)-Dc@T,-Lc@P))
    return basis,P,T,fitted,backward


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--checkpoint',required=True)
    ap.add_argument('--intervals',type=int,default=96);ap.add_argument('--output',required=True)
    a=ap.parse_args()
    base.fit_mode=fit_mode
    result=base.run(a)
    result['fit_metric']='Continuous L2 of the specified radial parity cubic interpolant before axial localization.'
    result['driver_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['radial_quadrature']='Gauss order m+18 on union of potential and target knots; polynomial Gram and RHS.'
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
