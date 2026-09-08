"""Independent exterior harmonic-energy and compatible-MAC checks.

The Green identities are analytic; numerical Bessel quadrature and matrix
tests below are floating-point checks, not whole-space NS validation.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import kve
import sympy as sp


def dtn(m,k,R):
    m=abs(int(m));k=abs(float(k))
    if not k:
        if not m:raise ValueError('The zero mode requires q=0 for finite radial energy.')
        return m/R
    x=k*R
    return k*(kve(abs(m-1),x)+kve(m+1,x))/(2*kve(m,x))


def exterior_jet(r,m,k,R):
    kap=dtn(m,k,R)
    if k:
        x=abs(k)*r;x0=abs(k)*R
        f=-np.exp(x0-x)*kve(m,x)/(kap*kve(m,x0))
        fp=np.exp(x0-x)*abs(k)*(kve(abs(m-1),x)+kve(m+1,x))/(2*kap*kve(m,x0))
    else:
        f=-R/m*(R/r)**m
        fp=(R/r)**(m+1)
    fpp=-fp/r+(m*m/r**2+k*k)*f
    return f,fp,fpp


def continuous_check(m,k,R=2.7):
    kap=dtn(m,k,R)
    def energy(r):
        f,fp,_=exterior_jet(r,m,k,R)
        return r*(fp*fp+(m*m/r**2+k*k)*f*f)
    def dissipation(r):
        f,fp,fpp=exterior_jet(r,m,k,R)
        diagonal=fpp*fpp+(fp/r-m*m*f/r**2)**2+k**4*f*f
        cross=m*m*(fp/r-f/r**2)**2+k*k*fp*fp+(m*k*f/r)**2
        return r*(diagonal+2*cross)
    mass_quad=quad(energy,R,np.inf,epsabs=2.e-11,epsrel=2.e-11)[0]
    stiff_quad=quad(dissipation,R,np.inf,epsabs=2.e-11,epsrel=2.e-11)[0]
    mass=R/kap
    stiff=1+2*R*(m*m/R**2+k*k)/kap+m*m/(R*kap)**2
    me=abs(mass_quad/mass-1);se=abs(stiff_quad/stiff-1)
    assert me<2.e-10 and se<2.e-10,(m,k,me,se)
    return {'m':m,'k':k,'mass':mass,'full_stiffness':stiff,
            'mass_quadrature_relative_error':me,'stiffness_quadrature_relative_error':se}


def compatible_matrix_check(m,k,n=7,R=2.7):
    rf=R*np.linspace(0,1,n+1)**1.25
    rc=(rf[1:]+rf[:-1])/2
    vol=(rf[1:]**2-rf[:-1]**2)/2
    delta=R-rc[-1];kap=dtn(m,k,R)
    W=np.r_[rf[1:-1]*np.diff(rc),R*(delta+1/kap)]
    masses=np.r_[W,vol,vol]
    D=np.zeros((n,3*n),complex)
    for i in range(n):
        D[i,i]=rf[i+1]/vol[i]
        if i:D[i,i-1]=-rf[i]/vol[i]
        D[i,n+i]=1j*m/rc[i]
        D[i,2*n+i]=1j*k
    G=-D.conj().T*vol[None,:]/masses[:,None]
    P=np.eye(3*n)-G@np.linalg.solve(D@G,D)
    MH=np.diag(masses)
    div=np.max(abs(D@P))
    idem=np.max(abs(P@P-P))
    orth=np.max(abs(P.conj().T@MH-MH@P))
    robin=-kap/(1+kap*delta)
    assert abs(G[n-1,n-1]-robin)<1.e-13
    assert np.count_nonzero(abs(G[n-1])>1.e-14)==1
    assert div<2.e-13 and idem<2.e-13 and orth<2.e-13
    # Random pressure gradient is completely removed by this same projector.
    rng=np.random.default_rng(423)
    pressure=rng.normal(size=n)+1j*rng.normal(size=n)
    gradient_leak=np.max(abs(P@(G@pressure)))
    assert gradient_leak<2.e-13
    return {'m':m,'k':k,'outer_mass':W[-1],'boundary_gradient':robin,
            'divergence':float(div),'idempotence':float(idem),
            'kinetic_self_adjoint_defect':float(orth),
            'gradient_removal_error':float(gradient_leak)}


def exact_boundary_check():
    R,kap,m,k=sp.symbols('R kap m k',positive=True,real=True)
    u=sp.Matrix([1,-sp.I*m/(R*kap),-sp.I*k/kap])
    ur=sp.Matrix([-1/R-(m*m/R**2+k*k)/kap,
                  sp.I*m*(1/R+1/(R*R*kap)),sp.I*k])
    d=-R*(u.conjugate().T*ur)[0]
    assert sp.simplify(d-(1+2*R*(m*m/R**2+k*k)/kap+m*m/(R*kap)**2))==0
    assert sp.simplify(d.subs({k:0,kap:m/R})-2*(m+1))==0


if __name__=='__main__':
    exact_boundary_check()
    pairs=[(0,.7),(0,2.3),(1,0),(4,0),(4,.7),(8,0),(8,2.3)]
    result={'status':'PASS','scope':'Exact boundary algebra and floating-point exterior/matrix tests; no evolved-field or NS residual bound',
            'continuous':[continuous_check(*p) for p in pairs],
            'discrete':[compatible_matrix_check(*p) for p in pairs]}
    Path(__file__).with_name('exterior-harmonic-mass-check.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
