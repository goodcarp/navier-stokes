"""Independent manufactured checks for the common-k projection reference.

Uses analytic Gaussian Hankel pairs, real-space differentiation, exact
projection algebra, weighted energy and forward/inverse checks.
Finite quadrature errors reported here are not rigorous enclosures.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.special import jn_zeros
from shared_hankel_projection import (symbols,project,forward_helical,
    inverse_helical,kinetic_inner)


def gauss(n,upper):
    x,w=np.polynomial.legendre.leggauss(n)
    return upper*(x+1)/2,upper*w/2


def gaussian_hat(m,alpha,k):
    return k**m/(2*alpha)**(m+1)*np.exp(-k*k/(4*alpha))


def exact_algebra():
    k,b=sp.symbols('k b',real=True)
    H=sp.diag(sp.Rational(1,2),sp.Rational(1,2),1)
    for sign in (1,-1):
        G=sp.Matrix([-k,sign*k,sp.I*b])
        D=sp.Matrix([[k/2,-sign*k/2,sp.I*b]])
        P=sp.eye(3)+G*D/(k*k+b*b)
        assert sp.simplify(D*G)[0]==-k*k-b*b
        assert sp.simplify(D+G.conjugate().T*H)==sp.zeros(1,3)
        assert sp.simplify(D*P)==sp.zeros(1,3)
        assert sp.simplify(P*P-P)==sp.zeros(3)
        assert sp.simplify(P.conjugate().T*H-H*P)==sp.zeros(3)
    # Independent toroidal/poloidal mass and viscous stiffness symbols.
    tor=sp.Matrix([sp.I*k,sp.I*k,0])
    pol=sp.Matrix([-sp.I*b*k,sp.I*b*k,k*k])
    assert sp.simplify(tor.conjugate().T*H*pol)[0]==0
    assert sp.simplify(tor.conjugate().T*H*tor)[0]==k*k
    assert sp.simplify(pol.conjugate().T*H*pol)[0]==k*k*(k*k+b*b)
    assert sp.expand((k*k+b*b)*(tor.conjugate().T*H*tor)[0])==k**4+b*b*k*k
    assert sp.expand((k*k+b*b)*(pol.conjugate().T*H*pol)[0])==k**6+2*b*b*k**4+b**4*k*k


def manufactured(m,signed=True,n=256):
    r,wr=gauss(n,12.)
    k,wk=gauss(n,20.)
    kap=1.3;alpha=.8;beta=1.1
    # W=curl(psi e_z)+grad(phi), with phi=r^m exp(-alpha*r^2).
    phi=r**m*np.exp(-alpha*r*r)
    psi=r**m*np.exp(-beta*r*r)
    phipr=(m/r-2*alpha*r)*phi
    psipr=(m/r-2*beta*r)*psi
    tor=np.array([1j*(m*psi/r-psipr),1j*(m*psi/r+psipr),np.zeros_like(r)])
    grad=np.array([phipr-m*phi/r,phipr+m*phi/r,1j*kap*phi])
    W=tor+grad
    What=forward_helical(W,r,wr,k,m,signed)
    G,D,q2=symbols(k,kap,m,signed)
    phihat=gaussian_hat(m,alpha,k)
    psihat=gaussian_hat(m,beta,k)
    torhat=np.array([1j*k*psihat,1j*k*psihat,np.zeros_like(k)])
    if not signed and m==0:torhat[1]*=-1
    exact=torhat+G*phihat
    transform_error=float(np.max(abs(What-exact))/max(1.,np.max(abs(exact))))
    projected=project(What,k,kap,m,signed)
    projection_error=float(np.max(abs(projected-torhat))/max(1.,np.max(abs(torhat))))
    recovered=inverse_helical(projected,k,wk,r,m,signed)
    recovery_error=float(np.max(abs(recovered-tor))/max(1.,np.max(abs(tor))))
    recovered_r=inverse_helical(projected,k,wk,r,m,signed,derivative=1)
    divergence=.5*(recovered_r[0]+(m+1)*recovered[0]/r
                      +recovered_r[1]+(1-m)*recovered[1]/r)+1j*kap*recovered[2]
    divergence_error=float(np.max(abs(divergence))/max(1.,np.max(abs(recovered))))
    reprojected=inverse_helical(project(forward_helical(recovered,r,wr,k,m,signed),
                                      k,kap,m,signed),k,wk,r,m,signed)
    roundtrip_error=float(np.max(abs(reprojected-recovered))/max(1.,np.max(abs(recovered))))
    physical_norm=np.sum(wr*r*(.5*abs(W[0])**2+.5*abs(W[1])**2+abs(W[2])**2))
    spectral_norm=kinetic_inner(What,What,wk,k).real
    parseval_error=float(abs(spectral_norm-physical_norm)/physical_norm)
    for error in (transform_error,projection_error,recovery_error,divergence_error,
                  roundtrip_error,parseval_error):
        assert error<2.e-10,(m,signed,error)
    return dict(m=m,signed_orders=signed,nr=n,nk=n,
                analytic_transform_relative_error=transform_error,
                gradient_removal_relative_error=projection_error,
                solenoidal_recovery_relative_error=recovery_error,
                physical_divergence_relative_error=divergence_error,
                projected_roundtrip_relative_error=roundtrip_error,
                Parseval_relative_error=parseval_error)


def random_projection():
    rng=np.random.default_rng(731)
    k,wk=gauss(97,30.);kap=2.3
    U=rng.normal(size=(3,len(k)))+1j*rng.normal(size=(3,len(k)))
    P=project(U,k,kap);R=U-P
    again=project(P,k,kap)
    G,D,q2=symbols(k,kap)
    div=float(np.max(abs(np.sum(D*P,axis=0))))
    idem=float(np.max(abs(again-P)))
    E=kinetic_inner(U,U,wk,k).real
    defect=E-kinetic_inner(P,P,wk,k).real-kinetic_inner(R,R,wk,k).real
    assert div<5.e-14 and idem<5.e-15 and abs(defect/E)<1.e-14
    return dict(divergence=div,idempotence=idem,orthogonal_energy_relative_defect=float(defect/E))


def incompatible_grid_demo():
    # Exact transforms of one smooth gradient, sampled at three distinct
    # Bessel-zero grids, cannot be projected as if their indices shared k.
    m=4;kap=1.3;alpha=.8;N=40;radius=12.
    grids=[jn_zeros(n,N)/radius for n in (m+1,m-1,m)]
    bad=np.array([-grids[0]*gaussian_hat(m,alpha,grids[0]),
                   grids[1]*gaussian_hat(m,alpha,grids[1]),
                   1j*kap*gaussian_hat(m,alpha,grids[2])])
    leftover=project(bad,grids[2],kap)
    relative=float(np.linalg.norm(leftover)/np.linalg.norm(bad))
    assert relative>.01
    return dict(naive_mismatched_grid_gradient_leak=relative,
                scope='Deliberately incorrect indexwise mixing; demonstrates why shared k or a compatible non-diagonal map is necessary.')


if __name__=='__main__':
    exact_algebra()
    tests=[manufactured(m,signed) for m in (0,4,8) for signed in (True,False)]
    out=dict(status='PASS',scope='Exact operator and toroidal/poloidal block algebra; Gaussian quadrature tests; no nonlinear evolution or error enclosure',
        manufactured=tests,random_projection=random_projection(),
        incompatible_grid_demo=incompatible_grid_demo())
    Path(__file__).with_name('shared-hankel-projection-check.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
