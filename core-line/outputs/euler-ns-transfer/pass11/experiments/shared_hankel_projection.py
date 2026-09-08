"""Common-radial-frequency helical Hankel projection reference.

Continuous transform convention H_n f(k)=integral_0^infinity r J_n(kr) f(r) dr,
inverse f(r)=integral_0^infinity k J_n(kr) H_n f(k) dk.

The projection is exact algebra on a common k. Quadrature matrices below
approximate the transforms; they are NOT exact finite-dimensional transforms.
Finite unwindowed sums of point-frequency Bessel waves are not whole-space L2
reconstructions. This module is an independent operator/roundtrip benchmark,
not an NS time integrator.
"""
import numpy as np
from scipy.special import jv,jvp


def helical_orders(m,signed=True):
    orders=(m+1,m-1,m)
    return orders if signed else tuple(abs(n) for n in orders)


def symbols(k,kappa,m=0,signed=True):
    """Return G=(3,*shape), D=(3,*shape), q2 for nonnegative m.

    With signed Bessel orders, G=(-k,k,i*kappa), D=(k/2,-k/2,i*kappa).
    If all transform orders are replaced by their absolute values, m=0
    instead has G=(-k,-k,i*kappa), D=(k/2,k/2,i*kappa).
    """
    if m<0 and not signed:
        raise ValueError('Absolute-order branch is implemented for m>=0 only.')
    k,kappa=np.broadcast_arrays(np.asarray(k,float),np.asarray(kappa,float))
    sign=-1 if (not signed and m==0) else 1
    G=np.stack((-k,sign*k,1j*kappa))
    D=np.stack((k/2,-sign*k/2,1j*kappa))
    return G,D,k*k+kappa*kappa


def project(U,k,kappa,m=0,signed=True):
    """Kinetic-orthogonal projection of helical spectra (Uplus,Uminus,Uz).

    U must have leading dimension 3. The zero (k,kappa) mode is left unchanged;
    constant whole-space modes are not L2 and are not represented by the
    Gaussian quadrature tests.
    """
    G,D,q2=symbols(k,kappa,m,signed)
    div=np.sum(D*U,axis=0)
    return U+G*div/np.where(q2==0,1.,q2)


def forward(values,r,weights_r,k,order):
    return jv(order,np.outer(k,r))@(np.asarray(weights_r)*r*values)


def inverse(spectrum,k,weights_k,r,order,derivative=0):
    matrix=jv(order,np.outer(r,k)) if derivative==0 else jvp(order,np.outer(r,k),derivative)*k[None,:]**derivative
    return matrix@(np.asarray(weights_k)*k*spectrum)


def forward_helical(values,r,weights_r,k,m,signed=True):
    return np.array([forward(value,r,weights_r,k,n)
                     for value,n in zip(values,helical_orders(m,signed))])


def inverse_helical(spectra,k,weights_k,r,m,signed=True,derivative=0):
    return np.array([inverse(value,k,weights_k,r,n,derivative)
                     for value,n in zip(spectra,helical_orders(m,signed))])


def kinetic_inner(U,V,weights_k,k):
    """Continuous-Plancherel metric approximated on the common k quadrature."""
    return np.sum(weights_k*k*(.5*np.conj(U[0])*V[0]
                              +.5*np.conj(U[1])*V[1]+np.conj(U[2])*V[2]))
