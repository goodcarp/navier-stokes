#!/usr/bin/env python3
"""Nonlinear planar slice diagnostic, NOT full compact 3D Navier--Stokes.

Evolves z-independent horizontal vorticity in the prescribed affine strain
(c x,c y,-2c z). All resolved angular interactions and planar pressure are
retained by vorticity/Biot--Savart. Axial endcaps and their pressure are absent.
Finite-volume radial / dealiased Fourier angular / CN-Heun time discretization.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import argparse, json, sys, time
from pathlib import Path
import numpy as np
from scipy.fft import rfft, irfft
from scipy.linalg.lapack import zgttrf, zgttrs
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
dependency = next(p for p in [HERE.parents[1] / 'outputs/euler-ns-transfer/pass4/experiments',
                             HERE.parents[1] / 'pass4/experiments']
                  if (p / 'evaluate_pressure_gate.py').is_file())
sys.path.insert(0, str(dependency))
from evaluate_pressure_gate import cutoff


class Tri:
    def __init__(self, lo, diag, hi):
        self.factors = zgttrf(np.asarray(lo, complex), np.asarray(diag, complex), np.asarray(hi, complex))
        assert self.factors[-1] == 0
    def solve(self, rhs):
        sol, info = zgttrs(*self.factors[:-1], np.asarray(rhs, complex))
        assert info == 0
        return sol


class Slice:
    def __init__(self, nr=768, nphi=48, R=.8, dt=1e-6, nu=.001, c=1.4, A=950., lam=.25):
        assert nphi % 6 == 0
        self.nr, self.nphi, self.R, self.dt = nr, nphi, R, dt
        self.nu, self.c, self.A, self.lam = nu, c, A, lam
        self.h = h = R/nr
        self.r = r = (np.arange(nr)+.5)*h
        self.rf = np.arange(nr+1)*h
        self.nmode = nphi//2+1
        self.keep = nphi//3
        self.m = 4*np.arange(self.nmode)
        self.poisson, self.step, self.L = [], [], []
        for m in self.m:
            lo = self.rf[1:-1]/(r[1:]*h*h)
            hi = self.rf[1:-1]/(r[:-1]*h*h)
            diag = -(self.rf[:-1]+self.rf[1:])/(r*h*h)-m*m/r**2
            # Zero outer vorticity at the boundary face; no inner radial flux.
            lapdiag = diag.copy(); lapdiag[-1] -= self.rf[-1]/(r[-1]*h*h)
            # Conservative affine radial advection: -div(c r q).
            Ll = nu*lo+c*self.rf[1:-1]**2/(2*r[1:]*h)
            Lu = nu*hi-c*self.rf[1:-1]**2/(2*r[:-1]*h)
            Ld = nu*lapdiag-c*(self.rf[1:]**2-self.rf[:-1]**2)/(2*r*h)
            # Zero interpolated outer vorticity also removes its advective flux.
            Ld[-1] += c*self.rf[-1]**2/(2*r[-1]*h)
            self.L.append((Ll,Ld,Lu))
            self.step.append(Tri(-dt*Ll/2,1-dt*Ld/2,-dt*Lu/2))
            if m:
                alpha = (1-m*h/(2*R))/(1+m*h/(2*R))
                d = diag.copy(); d[-1] += self.rf[-1]/(r[-1]*h*h)*alpha
                self.poisson.append(Tri(-lo,-d,-hi))
            else:
                self.poisson.append(None)
        x=(r/.315)**2
        C=cutoff(x); Cr=cutoff(x,1)*2*r/.315**2
        a=(r-.14)/.1
        e=cutoff(a*a)
        ep=cutoff(a*a,1)*2*a/.1
        epp=cutoff(a*a,2)*(2*a/.1)**2+cutoff(a*a,1)*2/.1**2
        k=-20.; m=4.
        b=e*np.exp(1j*k*r)
        bp=(ep+1j*k*e)*np.exp(1j*k*r)
        bpp=(epp+2j*k*ep-k*k*e)*np.exp(1j*k*r)
        q=np.zeros((nr,self.nmode),complex)
        q[:,0]=A*(2*C+r*Cr)
        q[:,1]=-.5*lam*(bpp+bp/r-m*m*b/r**2)
        self.q=q
        self.initial_psi = .5*lam*b
        self.initial_theta = -.5*lam*bp
        self.initial_radial = .5*lam*1j*m*b/r

    def velocities(self,q):
        psi=np.zeros_like(q)
        for j in range(1,self.keep):
            psi[:,j]=self.poisson[j].solve(q[:,j])
        faces=np.zeros((self.nr+1,self.nmode),complex)
        faces[1:-1]=.5*(psi[:-1]+psi[1:])
        for j in range(1,self.keep):
            faces[-1,j]=psi[-1,j]/(1+self.m[j]*self.h/(2*self.R))
        ur=np.zeros_like(faces)
        ur[1:]=1j*self.m[None,:]*faces[1:]/self.rf[1:,None]
        vt=-(faces[1:]-faces[:-1])/self.h
        # Cell-centered mean swirl from the integrated vorticity.
        circulation=self.h*np.cumsum(self.r*q[:,0].real)
        meanfaces=np.zeros(self.nr+1); meanfaces[1:]=circulation/self.rf[1:]
        vt[:,0]=.5*(meanfaces[:-1]+meanfaces[1:])
        return psi,ur,vt

    def nonlinear(self,q):
        _,ur,vt=self.velocities(q)
        qf=np.zeros((self.nr+1,self.nmode),complex)
        qf[1:-1]=.5*(q[:-1]+q[1:])
        rf=rfft(irfft(ur,n=self.nphi,axis=1,norm='forward')*
                 irfft(qf,n=self.nphi,axis=1,norm='forward'),axis=1,norm='forward')
        tf=rfft(irfft(vt,n=self.nphi,axis=1,norm='forward')*
                 irfft(q,n=self.nphi,axis=1,norm='forward'),axis=1,norm='forward')
        result=-(self.rf[1:,None]*rf[1:]-self.rf[:-1,None]*rf[:-1])/(self.h*self.r[:,None])
        result-=1j*self.m[None,:]*tf/self.r[:,None]
        result[:,self.keep:]=0
        return result

    def linear(self,q):
        out=np.zeros_like(q)
        for j in range(self.keep):
            lo,di,hi=self.L[j]
            out[:,j]=di*q[:,j]
            out[1:,j]+=lo*q[:-1,j]
            out[:-1,j]+=hi*q[1:,j]
        return out

    def solve_step(self,rhs):
        out=np.zeros_like(rhs)
        for j in range(self.keep):
            out[:,j]=self.step[j].solve(rhs[:,j])
        return out

    def advance(self):
        dt=self.dt; q=self.q
        base=q+.5*dt*self.linear(q)
        n0=self.nonlinear(q)
        predictor=self.solve_step(base+dt*n0)
        self.q=self.solve_step(base+.5*dt*(n0+self.nonlinear(predictor)))

    def observe(self,t):
        q=self.q; psi,ur,vt=self.velocities(q)
        uc=1j*self.m[None,:]*psi/self.r[:,None]
        G=self.r*vt[:,0].real
        interp=CubicSpline(self.r,G)
        ri=brentq(interp.derivative(),.18,.28)
        fluct=2*np.pi*self.h*np.sum(self.r[:,None]*(abs(uc[:,1:])**2+abs(vt[:,1:])**2))
        mean=np.pi*self.h*np.sum(self.r*abs(vt[:,0])**2)
        diss=2*np.pi*self.h*np.sum(self.r*(abs(q[:,0])**2+2*np.sum(abs(q[:,1:])**2,axis=1)))
        covariance=2*np.real(np.sum(uc[:,1:]*np.conj(vt[:,1:]),axis=1))
        torque=-CubicSpline(self.r,self.r**2*covariance).derivative()(ri)/ri
        initial_errors=None
        if t==0:
            mask=(self.r>.03)&(self.r<.35)
            initial_errors=dict(psi_max=float(np.max(abs(psi[mask,1]-self.initial_psi[mask]))),
                radial_velocity_max=float(np.max(abs(uc[mask,1]-self.initial_radial[mask]))),
                tangential_velocity_max=float(np.max(abs(vt[mask,1]-self.initial_theta[mask]))))
        # The stored name means this spline critical branch in (.18,.28).
        # Energy integrals stop at R; irrotational exterior tails are excluded.
        return dict(t=float(t),mean_maximum=float(interp(ri)),mean_maximum_radius=float(ri),
            mean_fixed_initial_receiver=float(interp(.2154755753335)),
            fluctuation_energy_per_axial_length=float(fluct),
            material_slice_weighted_fluctuation_energy=float(np.exp(-2*self.c*t)*fluct),
            total_horizontal_energy_per_axial_length=float(mean+fluct),
            dissipation_integrand=float(diss),torque_at_radial_maximum=float(torque),
            mean_total_circulation_over_2pi=float(self.h*np.sum(self.r*q[:,0].real)),
            highest_retained_mode_vorticity_L2=float(np.sqrt(self.h*np.sum(self.r*abs(q[:,self.keep-1])**2))),
            initial_Biot_Savart_errors=initial_errors)


def run(nr=768,nphi=48,R=.8,dt=1e-6,nu=.001,c=1.4,A=950.,lam=.25,T=.002,outputs=40,output=None):
    start=time.time(); obj=Slice(nr,nphi,R,dt,nu,c,A,lam)
    steps=round(T/dt); assert abs(steps*dt-T)<1e-12
    stride=max(1,steps//outputs); records=[obj.observe(0)]
    for step in range(1,steps+1):
        obj.advance()
        if step%stride==0 or step==steps:
            obs=obj.observe(step*dt)
            assert np.isfinite(obs['total_horizontal_energy_per_axial_length'])
            records.append(obs)
    tt=np.array([x['t'] for x in records]); EE=np.array([x['total_horizontal_energy_per_axial_length'] for x in records])
    ZZ=np.array([x['dissipation_integrand'] for x in records])
    energy_balance=float(EE[-1]-EE[0]+nu*np.trapz(ZZ,tt))
    result=dict(status='EXPLORATORY finite nonlinear planar model, not compact 3D NS and not an interval certificate',
        parameters=dict(nr=nr,nphi=nphi,R=R,dt=dt,nu=nu,c=c,A=A,lam=lam,T=T),
        missing_full_3D_effects=['axial envelopes and pressure','evolving exterior meridional strain','central core and older field','finite amplitude 3D stability'],
        records=records,total_energy_balance_residual=energy_balance,
        total_energy_balance_relative_residual=energy_balance/EE[0],seconds=time.time()-start)
    text=json.dumps(result,indent=2)+'\n'
    if output:Path(output).write_text(text)
    else:print(text)
    print(json.dumps(dict(output=output,seconds=result['seconds'],first=records[0],last=records[-1],relative_energy_balance=result['total_energy_balance_relative_residual'])),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    for name,typ,default in [('nr',int,768),('nphi',int,48),('R',float,.8),('dt',float,1e-6),('nu',float,.001),('c',float,1.4),('A',float,950.),('lam',float,.25),('T',float,.002),('outputs',int,40)]:
        ap.add_argument('--'+name,type=typ,default=default)
    ap.add_argument('--output')
    run(**vars(ap.parse_args()))
