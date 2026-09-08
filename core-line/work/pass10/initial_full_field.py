"""Analytic initial full-field evaluator, using genuine Fourier coefficients.

u = S_(chi+c eta) + W_chi + A v + lam w_leading, with defaults
A=1020, lam=1/4, c=7/5. No pressure-neutral retuning is implicit.

fields(r,z) broadcasts r,z and returns:
  jets[m]: (3,6,*shape), derivatives (value,r,z,rr,rz,zz), m=0,4;
  u[m]: (3,*shape), gradient[m]: (3,3,*shape), m=0,4;
  divergence[m], laplacian[m], m=0,4;
  source[m] = (-Delta p)_m, convection[m] = ((u.grad)u)_m, m=0,4,8.
Rows of gradient are output components and columns are differentiation
directions, in the orthonormal cylindrical basis (r,theta,z).

All positive modes are Fourier coefficients, NOT real-wave amplitudes:
f(theta)=f[0]+2 Re(sum_(m>0) f[m] exp(im theta)). Negative coefficients
are conjugates. The axial/radial cutoffs and all their transitions are kept.
This is floating-point evaluation of analytic formulas, not an interval
certificate or an NS time integrator. Only numpy is required.
"""
from math import comb
import numpy as np

DERIVATIVES = ('value', 'r', 'z', 'rr', 'rz', 'zz')
COMPONENTS = ('r', 'theta', 'z')


def cutoff(x, n=0):
    """The unchanged C-infinity cutoff and argument derivatives through order 3.

    chi=1 for x<=1/4, chi=0 for x>=1. On the transition, t=(4x-1)/3,
    chi=logistic(1/t-1/(1-t)). The endpoint guard discards only derivatives
    below floating-point underflow, avoiding 0*infinity intermediates.
    """
    if n not in (0, 1, 2, 3):
        raise ValueError('Only cutoff derivatives 0,1,2,3 are implemented.')
    x=np.asarray(x,dtype=float)
    t=(4*x-1)/3
    inside=(t>0)&(t<1)
    a=np.clip(t,1.e-7,1-1.e-7)
    h=1/a-1/(1-a)
    b=np.exp(-np.abs(h))
    f=np.where(h>=0,1/(1+b),b/(1+b))
    if n==0:
        return np.where(t<=0,1.,np.where(t>=1,0.,f))
    # f*(1-f) via the symmetric form avoids losing small left-tail jets.
    logistic_first=b/(1+b)**2
    h1=-a**-2-(1-a)**-2
    h2=2*a**-3-2*(1-a)**-3
    h3=-6*a**-4-6*(1-a)**-4
    if n==1:
        out=logistic_first*h1*(4/3)
    elif n==2:
        out=logistic_first*((1-2*f)*h1*h1+h2)*(4/3)**2
    else:
        out=logistic_first*((1-6*f+6*f*f)*h1**3
                           +3*(1-2*f)*h1*h2+h3)*(4/3)**3
    return np.where(inside,out,0.)


def eta(s,n=0):
    """Radial annular profile, derivatives with respect to s=|x|^2.

    Exactly equal to -(1-chi(4s/25))*chi((3s-64)/44), since the two
    transition supports are disjoint. Zero at R<=5/4 and R>=6;
    equals -1 throughout 5/2<=R<=5.
    """
    return (4/25)**n*cutoff(4*np.asarray(s)/25,n)\
           -(3/44)**n*cutoff((3*np.asarray(s)-64)/44,n)


def _product(a,b):
    """Product of scalar two-variable jets in the declared six-entry order."""
    return np.stack((a[0]*b[0], a[1]*b[0]+a[0]*b[1],
        a[2]*b[0]+a[0]*b[2],
        a[3]*b[0]+2*a[1]*b[1]+a[0]*b[3],
        a[4]*b[0]+a[1]*b[2]+a[2]*b[1]+a[0]*b[4],
        a[5]*b[0]+2*a[2]*b[2]+a[0]*b[5]))


def _radial_jet(r,z,f0,f1,f2):
    return np.stack((f0,2*r*f1,2*z*f1,2*f1+4*r*r*f2,
                     4*r*z*f2,2*f1+4*z*z*f2))


def _quadratic_composition(x,center,width,order=3):
    """Derivatives in x of chi(((x-center)/width)^2)."""
    y=(x-center)/width
    argument=y*y
    out=[cutoff(argument)]
    if order>=1:out.append(2*y/width*cutoff(argument,1))
    if order>=2:out.append((2*cutoff(argument,1)+4*y*y*cutoff(argument,2))/width**2)
    if order>=3:out.append((12*y*cutoff(argument,2)+8*y**3*cutoff(argument,3))/width**3)
    return out


def _axial_cutoff(z,width):
    plus=_quadratic_composition(z,4.,width,2)
    minus=_quadratic_composition(z,-4.,width,2)
    return [a+b for a,b in zip(plus,minus)]


def _gradient(jets,m,r):
    u=jets[:,0]; ur=jets[:,1]; uz=jets[:,2]
    safe=np.where(r==0,1.,r)
    angular=np.stack(((1j*m*u[0]-u[1])/safe,
                      (1j*m*u[1]+u[0])/safe,1j*m*u[2]/safe))
    if m==0:
        angular[0]=np.where(r==0,-ur[1],angular[0])
        angular[1]=np.where(r==0,ur[0],angular[1])
    # Every m=4 jet vanishes in a fixed axis neighborhood for this datum.
    return np.stack((ur,angular,uz),axis=1)


def _laplacian(jets,m,r):
    safe=np.where(r==0,1.,r)
    u=jets[:,0]
    out=jets[:,3]+jets[:,1]/safe+jets[:,5]-m*m*u/safe**2
    out[0]-=(u[0]+2j*m*u[1])/safe**2
    out[1]-=(u[1]-2j*m*u[0])/safe**2
    if m==0:
        out[0]=np.where(r==0,0.,out[0])
        out[1]=np.where(r==0,0.,out[1])
        out[2]=np.where(r==0,2*jets[2,3]+jets[2,5],out[2])
    return out


def fields(r,z,A=1020.,lam=.25,c=1.4):
    """Evaluate the full datum and its exact analytic modal differential terms."""
    r,z=np.broadcast_arrays(np.asarray(r,dtype=float),np.asarray(z,dtype=float))
    if np.any(r<0):raise ValueError('Cylindrical radius r must be nonnegative.')
    zero=np.zeros_like(r); one=np.ones_like(r)
    Rj=np.stack((r,one,zero,zero,zero,zero))
    Zj=np.stack((z,zero,one,zero,zero,zero))
    R2=_product(Rj,Rj); Z2=_product(Zj,Zj); Sj=R2+Z2
    s=r*r+z*z
    ps=[cutoff(s,n) for n in range(4)]
    ph=[ps[n]+c*eta(s,n) for n in range(4)]
    Phi=_radial_jet(r,z,*ph[:3]); Phi1=_radial_jet(r,z,*ph[1:])
    Psi=_radial_jet(r,z,*ps[:3]); Psi1=_radial_jet(r,z,*ps[1:])
    meridional_r=-_product(Rj,Phi+2*_product(Z2,Phi1))
    meridional_z=2*_product(Zj,Phi+_product(R2,Phi1))
    core_swirl=_product(Rj,Psi+(2/3)*_product(Sj,Psi1))

    # Unit axisymmetric outer swirl v=r chi(r^2/a^2) q_v(z), a=63/200.
    f,fr,frr=_quadratic_composition(r,0.,63/200,2)
    fjet=np.stack((f,fr,zero,frr,zero,zero))
    q,qz,qzz=_axial_cutoff(z,9/20)
    qjet=np.stack((q,zero,qz,zero,zero,qzz))
    outer_swirl=_product(_product(Rj,fjet),qjet)
    mean=np.stack((meridional_r,core_swirl+A*outer_swirl,meridional_z)).astype(complex)

    # Leading seed curl[e(r)q_s(z) cos(m theta+k r) e_z].
    m=4; k=-20.
    e=_quadratic_composition(r,7/50,1/10,3)
    phase=np.exp(1j*k*r)
    b=[phase*sum(comb(n,j)*(1j*k)**(n-j)*e[j] for j in range(n+1))
       for n in range(4)]
    q,qz,qzz=_axial_cutoff(z,3/5)
    safe=np.where(r==0,1.,r)
    br=b[1]/safe-b[0]/safe**2
    brr=b[2]/safe-2*b[1]/safe**2+2*b[0]/safe**3
    seed_r=1j*m*np.stack((b[0]*q/safe,br*q,b[0]*qz/safe,
                          brr*q,br*qz,b[0]*qzz/safe))
    seed_theta=-np.stack((b[1]*q,b[2]*q,b[1]*qz,b[3]*q,b[2]*qz,b[1]*qzz))
    seed_z=np.zeros_like(seed_r)
    seed=(lam/2)*np.stack((seed_r,seed_theta,seed_z))
    jets={0:mean,4:seed}
    velocity={m:J[:,0] for m,J in jets.items()}
    gradients={m:_gradient(J,m,r) for m,J in jets.items()}
    divergence={m:np.einsum('ii...->...',G) for m,G in gradients.items()}
    laplacian={m:_laplacian(J,m,r) for m,J in jets.items()}
    G0,G4=gradients[0],gradients[4]
    u0,u4=velocity[0],velocity[4]
    trace_product=lambda G,H:np.einsum('ij...,ji...->...',G,H)
    matvec=lambda G,v:np.einsum('ij...,j...->i...',G,v)
    source={0:np.real(trace_product(G0,G0)+2*trace_product(G4,G4.conj())),
            4:2*trace_product(G0,G4),8:trace_product(G4,G4)}
    convection={0:np.real(matvec(G0,u0)+matvec(G4,u4.conj())+matvec(G4.conj(),u4)),
                4:matvec(G0,u4)+matvec(G4,u0),8:matvec(G4,u4)}
    return dict(jets=jets,u=velocity,gradient=gradients,divergence=divergence,
                laplacian=laplacian,source=source,convection=convection)


def reconstruct(coefficients,theta):
    """Reconstruct a real cylindrical scalar/vector/tensor from modal coefficients.

    theta must broadcast against the coefficient arrays (after their leading
    component dimensions). For an extra angular axis, add a trailing singleton
    spatial dimension to r,z when calling fields.
    """
    result=np.real(coefficients[0]).copy()
    for m,value in coefficients.items():
        if m:result=result+2*np.real(value*np.exp(1j*m*np.asarray(theta)))
    return result


def cartesian_velocity(x,y,z,**kwargs):
    """Convenience view for independent Cartesian difference checks."""
    x,y,z=np.broadcast_arrays(np.asarray(x),np.asarray(y),np.asarray(z))
    r=np.hypot(x,y); theta=np.arctan2(y,x)
    u=reconstruct(fields(r,z,**kwargs)['u'],theta)
    co=np.cos(theta); si=np.sin(theta)
    return np.stack((co*u[0]-si*u[1],si*u[0]+co*u[1],u[2]))
