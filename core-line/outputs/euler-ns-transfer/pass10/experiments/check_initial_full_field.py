"""Independent structural and Cartesian-difference checks for the initial field.

Finite-difference checks are numerical consistency tests, not enclosures.
The divergence and mode-source identities are also checked symbolically.
"""
import json
from pathlib import Path
import numpy as np
import sympy as sp
from initial_full_field import fields,reconstruct,cartesian_velocity,cutoff,_quadratic_composition,_axial_cutoff


def exact_checks():
    r,z,m,k=sp.symbols('r z m k',real=True,nonzero=True)
    f=sp.Function('f'); b=sp.Function('b');q=sp.Function('q')
    s=r*r+z*z
    Phi=f(s)
    ur=-r*(Phi+z*sp.diff(Phi,z))
    uz=2*z*Phi+r*z*sp.diff(Phi,r)
    assert sp.simplify(sp.diff(ur,r)+ur/r+sp.diff(uz,z))==0
    wr=sp.I*m*b(r)*q(z)/r;wt=-sp.diff(b(r),r)*q(z)
    assert sp.simplify(sp.diff(wr,r)+(wr+sp.I*m*wt)/r)==0
    V=sp.Function('V')(r)
    Gv=sp.Matrix([[0,-V/r,0],[sp.diff(V,r),0,0],[0,0,0]])
    w=sp.Matrix([wr,wt,0])
    Gw=sp.Matrix.hstack(w.diff(r),sp.Matrix([(sp.I*m*wr-wt)/r,
                         (sp.I*m*wt+wr)/r,0]),w.diff(z))
    expected=(2/r)*(sp.diff(V*sp.diff(b(r),r),r)-m*m*sp.diff(V,r)*b(r)/r)*q(z)
    assert sp.simplify(2*sp.trace(Gv*Gw)-expected)==0


def rotate_vector(u,theta):
    co=np.cos(theta);si=np.sin(theta)
    return np.stack((co*u[0]-si*u[1],si*u[0]+co*u[1],u[2]))


def points():
    # Deliberate core, annular joins, seed-envelope joins, both axial joins,
    # axis, exterior, and full-affine neighborhoods; not plateau-only tests.
    rr=np.array([0.,0.,0.,.1,.6,.3,1.6,2.2,3.1,5.3,5.8,
        .2154755753335,.205,.229,.23,.22,.08,.06,.15,.29,.1,.215,.4,6.2])
    zz=np.array([0.,4.,5.4,.1,.2,.7,.2,.3,.7,.3,.1,
        4.,3.65,4.31,4.52,-3.62,-4.38,3.52,4.,-4.25,3.58,-4.55,4.,.3])
    th=np.linspace(.137,5.931,len(rr))
    return rr*np.cos(th),rr*np.sin(th),zz,rr,th


def run():
    exact_checks()
    x,y,z,r,theta=points()
    data=fields(r,z)
    uc=reconstruct(data['u'],theta)
    G=reconstruct(data['gradient'],theta)
    Q=np.zeros((3,3,len(r)));Q[0,0]=Q[1,1]=np.cos(theta)
    Q[1,0]=np.sin(theta);Q[0,1]=-np.sin(theta);Q[2,2]=1
    Gcart=np.einsum('ikn,kln,jln->ijn',Q,G,Q)
    axes=np.eye(3)
    def velocity_at(offset):return cartesian_velocity(x+offset[0],y+offset[1],z+offset[2])
    first_errors=[];second_errors=[]
    lap_exact=rotate_vector(reconstruct(data['laplacian'],theta),theta)
    source_exact=reconstruct(data['source'],theta)
    div_convection_errors=[]
    for h in (2.e-5,1.e-5):
        numerical_gradient=np.stack([(velocity_at(axes[j]*h)-velocity_at(-axes[j]*h))/(2*h)
                                      for j in range(3)],axis=1)
        first_errors.append(float(np.max(abs(numerical_gradient-Gcart))/max(1.,np.max(abs(Gcart)))))
        numerical_lap=sum((velocity_at(axes[j]*h)-2*rotate_vector(uc,theta)+velocity_at(-axes[j]*h))/h**2
                          for j in range(3))
        second_errors.append(float(np.max(abs(numerical_lap-lap_exact))/max(1.,np.max(abs(lap_exact)))))
        div_a=np.zeros_like(r)
        for j in range(3):
            values=[]
            for sign in (1,-1):
                xx=x+sign*h*axes[j,0];yy=y+sign*h*axes[j,1];zz=z+sign*h*axes[j,2]
                rr=np.hypot(xx,yy);tt=np.arctan2(yy,xx)
                a=rotate_vector(reconstruct(fields(rr,zz)['convection'],tt),tt)
                values.append(a[j])
            div_a+=(values[0]-values[1])/(2*h)
        div_convection_errors.append(float(np.max(abs(div_a-source_exact))/max(1.,np.max(abs(source_exact)))))
    assert first_errors[-1]<2.e-6
    assert second_errors[-1]<2.e-5
    assert div_convection_errors[-1]<4.e-6
    assert first_errors[-1]<first_errors[0]/3.5
    assert div_convection_errors[-1]<div_convection_errors[0]/3.5
    divergence=max(float(np.max(abs(v))) for v in data['divergence'].values())
    assert divergence<1.e-10
    source_direct=np.einsum('ij...,ji...->...',G,G)
    source_reconstruction=float(np.max(abs(source_exact-source_direct))/max(1.,np.max(abs(source_direct))))
    assert source_reconstruction<2.e-15
    convection_direct=np.einsum('ij...,j...->i...',G,uc)
    convection_reconstruction=float(np.max(abs(reconstruct(data['convection'],theta)-convection_direct))/max(1.,np.max(abs(convection_direct))))
    assert convection_reconstruction<2.e-15

    # Independent explicit full mixed source from the audited pressure equation.
    safe=np.where(r==0,1.,r);m=4;k=-20;A=1020;lam=.25
    e,er,err=_quadratic_composition(r,7/50,1/10,2)
    phase=np.exp(1j*k*r);b=e*phase
    bp=(er+1j*k*e)*phase;bpp=(err+2j*k*er-k*k*e)*phase
    C=cutoff((r/(63/200))**2)
    Cr=2*r/(63/200)**2*cutoff((r/(63/200))**2,1)
    V=r*C;Vp=C+r*Cr
    q=_axial_cutoff(z,3/5)[0]*_axial_cutoff(z,9/20)[0]
    mixed=(lam*A/safe)*(Vp*bp+V*bpp-m*m*Vp*b/safe)*q
    mixed_error=float(np.max(abs(data['source'][4]-mixed))/max(1.,np.max(abs(mixed))))
    assert mixed_error<2.e-14
    assert np.max(abs(fields(r,z,A=0)['source'][4]))<1.e-10
    B1=bp/safe-b/safe**2;B2=bp/safe-m*m*b/safe**2
    qs=_axial_cutoff(z,3/5)[0]
    seed0=lam**2*qs**2*(m*m*abs(B1)**2-np.real(B2*np.conj(bpp)))
    seed8=-(lam**2/2)*qs**2*(m*m*B1**2+B2*bpp)
    seed0_from_full=data['source'][0]-fields(r,z,lam=0)['source'][0]
    assert np.max(abs(seed0-seed0_from_full))/max(1.,np.max(abs(seed0)))<2.e-13
    assert np.max(abs(seed8-data['source'][8]))/max(1.,np.max(abs(seed8)))<2.e-14

    affine=fields(np.array([0.,.1]),np.array([0.,.1]))
    assert np.max(abs(affine['source'][0]-4))<1.e-14
    assert np.max(abs(affine['laplacian'][0]))<1.e-14
    outer_axis=fields(0.,4.)
    assert abs(outer_axis['source'][0]-(6*1.4**2-2*1020**2))<1.e-8
    for group in data.values():
        for value in group.values():assert np.all(np.isfinite(value))
    return dict(status='PASS',scope='Analytic identities plus floating-point consistency tests, not interval or PDE certification',
        exact_checks=['radial strain divergence','arbitrary envelope divergence','full mixed swirl/seed pressure source'],
        numerical_divergence_max=divergence,
        Cartesian_gradient_relative_errors=first_errors,
        Cartesian_laplacian_relative_errors=second_errors,
        Cartesian_divergence_of_convection_relative_errors=div_convection_errors,
        Fourier_source_reconstruction_relative_error=source_reconstruction,
        Fourier_convection_reconstruction_relative_error=convection_reconstruction,
        full_mixed_source_relative_error=mixed_error)


if __name__=='__main__':
    result=run()
    Path(__file__).with_name('initial-full-field-check.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
