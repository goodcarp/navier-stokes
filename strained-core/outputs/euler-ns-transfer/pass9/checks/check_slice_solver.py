#!/usr/bin/env python3
"""Independent analytic-solution and conservation tests for slice discretization.

Passing confirms selected numerical properties, not a continuum error bound.
"""
import json
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'experiments'))
from evolve_circular_slice import Slice

def check(nr):
    h=.8/nr; dt=2e-6; T=.001; a=.15; c=1.4; nu=.001
    s=Slice(nr=nr,nphi=24,dt=dt,lam=0,A=0)
    r=s.r; m=4
    psi=r**m*np.exp(-(r/a)**2)
    source=4*r**m*((m+1)/a**2-r*r/a**4)*np.exp(-(r/a)**2)
    q=np.zeros_like(s.q);q[:,1]=source
    pp,ur,vt=s.velocities(q)
    mask=r<.45
    psierror=float(np.max(abs(pp[mask,1]-psi[mask]))/np.max(psi))
    div=(s.rf[1:,None]*ur[1:]-s.rf[:-1,None]*ur[:-1])/(h*r[:,None])+1j*s.m[None,:]*vt/r[:,None]
    diverror=float(np.max(abs(div)))
    # Exact axisymmetric Gaussian evolution under viscosity and affine dilation.
    K=100.;s.q[:]=0;s.q[:,0]=K*np.exp(-(r/a)**2)
    for _ in range(round(T/dt)):s.advance()
    b2=a*a*np.exp(2*c*T)+(2*nu/c)*np.expm1(2*c*T)
    exact=K*a*a/b2*np.exp(-r*r/b2)
    gaussianerror=float(np.max(abs(s.q[mask,0]-exact[mask]))/K)
    # The integral of the nonlinear vorticity flux vanishes by periodicity
    # and the zero outer vorticity face. Use the actual nonaxisymmetric seed.
    s=Slice(nr=nr,nphi=24,dt=dt,A=950,lam=.25)
    NN=s.nonlinear(s.q)
    circulation=float(abs(s.h*np.sum(s.r*NN[:,0])))
    return dict(nr=nr,poisson_relative_max_error=psierror,discrete_divergence_max=diverror,
                Gaussian_relative_max_error=gaussianerror,nonlinear_circulation_error=circulation)

results=[check(n) for n in (192,384,768)]
for row in results:
    assert row['discrete_divergence_max']<1e-10
    assert row['nonlinear_circulation_error']<1e-9
for a,b in zip(results,results[1:]):
    assert 3.5<a['poisson_relative_max_error']/b['poisson_relative_max_error']<4.5
    assert 3.5<a['Gaussian_relative_max_error']/b['Gaussian_relative_max_error']<4.5
print(json.dumps(dict(status='PASS',scope='Analytic Gaussian evolution, known Poisson solution, discrete incompressibility and circulation; no PDE certification',results=results),indent=2))
