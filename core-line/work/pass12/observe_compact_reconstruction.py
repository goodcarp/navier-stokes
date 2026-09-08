#!/usr/bin/env python3
"""Continuous-form quadrature diagnostics for compact potential endpoints.

No interval bounds. Radial potential Gram integrands are piecewise
polynomials integrated by adequate Gauss order; analytic base cross terms
and axial trapezoidal quadrature still require convergence/error bounds.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.fft import fft,ifft
from scipy.optimize import minimize_scalar
from numpy.polynomial.legendre import leggauss
from compact_potential_reconstruction import RadialBasis,axial_window,evaluate_increment,fields


def gauss_spans(breaks,order):
    x,w=leggauss(order);breaks=np.unique(breaks)
    a=breaks[:-1,None];b=breaks[1:,None]
    return ((a+b)/2+(b-a)/2*x).ravel(),((b-a)/2*w).ravel()


def synth(coeff,k,z_origin,L,nz):
    # Continuous Fourier series, with exact signed-frequency placement.
    ix=np.rint(k*L/(2*np.pi)).astype(int)
    if 2*np.max(abs(ix))>=nz:
        raise ValueError('Output grid must retain all input signed frequencies.')
    z=(np.arange(nz)-nz//2)*L/nz
    target=np.zeros((len(coeff),nz),complex)
    target[:,ix%nz]=coeff*np.exp(1j*k*(z[0]-z_origin))
    return z,ifft(target,axis=-1,norm='forward')


def diagnostics(filename,nz=4099,gauss=36):
    data=np.load(filename);p=json.loads(str(data['parameters']))
    k=data['k'];origin=float(data['z_origin']);intervals=int(data['intervals']);L=p['L']
    b0=RadialBasis(p['R'],intervals,p['mapping'],0)
    bs0=b0.spline(np.array([0.]),nu=1)[0]
    at_zero=np.exp(-1j*k*origin)
    omega=1-2*float(np.real(bs0@data['T_0']@at_zero))
    b=1-2*float(np.real(bs0@(data['P_0']*(1j*k))@at_zero))
    def mean_G(r):
        u=fields(np.array([r]),np.array([4.]))['u'][0]
        inc=evaluate_increment(b0,data['P_0'],data['T_0'],k,origin,L,[r],[4.])
        return float(r*np.real(u[1,0]+inc[1,0,0]))
    opt=minimize_scalar(lambda r:-mean_G(r),bounds=(.18,.3),method='bounded',options={'xatol':1.e-13})
    modal=[];cross4=0.;initialK=0.
    for m in range(4,4*p['J']+1,4):
        basis=RadialBasis(p['R'],intervals,p['mapping'],m)
        if gauss < m+18:
            raise ValueError('Insufficient Gauss order for exact radial Gram polynomial degrees.')
        breaks=np.sqrt(np.unique(basis.knots))
        r,w=gauss_spans(breaks,gauss)
        B,Dr,Lap,mB=basis.matrices(r)
        wr=w*r
        H=Dr.T@(wr[:,None]*Dr)+mB.T@(wr[:,None]*mB)
        Q=Lap.T@(wr[:,None]*Lap)
        z,P=synth(data[f'P_{m}'],k,origin,L,nz)
        _,Pz=synth(data[f'P_{m}']*(1j*k),k,origin,L,nz)
        _,T=synth(data[f'T_{m}'],k,origin,L,nz)
        q,qp=axial_window(z,outer=L/2)
        Pz=q*Pz+qp*P;P=q*P;T=q*T
        quad=lambda A,G:float(np.real(np.sum(A.conj()*(G@A))))
        # K=1/2 norm² and a positive angular mode has conjugate factor two.
        energy=2*np.pi*L/nz*(quad(T,H)+quad(Pz,H)+quad(P,Q))
        modal.append({'mode':m,'increment_energy':energy})
        if m==4:
            # Only the original m=4 seed contributes to the base cross term.
            seed_breaks=np.unique(np.r_[.04,.09,.19,.24,breaks[(breaks>.04)&(breaks<.24)]])
            rc,wc=gauss_spans(seed_breaks,gauss)
            axial=(abs(z-4)<.6)|(abs(z+4)<.6)
            zz=z[axial];Pc=P[:,axial];Pzc=Pz[:,axial];Tc=T[:,axial]
            cross=base_sq=0.
            for start in range(0,len(rc),40):
                rr=rc[start:start+40];ww=wc[start:start+40]*rr
                Bc,Dc,Lc,mc=basis.matrices(rr)
                inc=np.stack((Dc@Pzc+1j*mc@Tc,1j*mc@Pzc-Dc@Tc,-Lc@Pc))
                base=fields(rr[:,None],zz[None,:])['u'][4]
                cross+=float(np.real(np.sum(ww[None,:,None]*base.conj()*inc)))
                base_sq+=float(np.sum(ww[None,:,None]*abs(base)**2))
            cross4=4*np.pi*L/nz*cross
            initialK=2*np.pi*L/nz*base_sq
    return {'scope':'Floating continuous-form endpoint diagnostics; no NS residual or error certificate.',
       'file':str(filename),'nz_quadrature':nz,'gauss_order_per_radial_span':gauss,
       'core_Omega':omega,'core_b':b,'core_beta':b/omega,'core_margin_b_minus_Omega':b-omega,
       'mean_G_fixed':mean_G(.21547557533348247),'mean_G_selected':mean_G(opt.x),
       'mean_G_selected_radius':float(opt.x),'initial_K_quadrature':initialK,
       'increment_modal_energies':modal,'initial_increment_cross':cross4,
       'endpoint_K_quadrature':initialK+sum(a['increment_energy'] for a in modal)+cross4}


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('input');ap.add_argument('--nz',type=int,default=4099)
    ap.add_argument('--gauss',type=int,default=36);ap.add_argument('--output',required=True)
    a=ap.parse_args();res=diagnostics(a.input,a.nz,a.gauss)
    Path(a.output).write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
