"""Shared: Gegenbauer C_l^{3/2} by recurrence, bang-bang shell coefficients g_l = h_l,
and the EXACT term-by-term e-fold integral of the shell operator (the far/near split made explicit).

Shell representation (S1):  a(x) = INT_{rho0}^{R} Phi[w(rho' .)](x/rho') drho'/rho'.
For the scale-invariant bang-bang cap w = -M sgn(z) every shell has the same w, so with
u = rho/rho' in (rho/R, rho/rho0):
   a(rho,phi) = INT_{rho/R}^{rho/rho0} Phi(u,phi) du/u
              = -SUM_{l odd} ((l+2) g_l/(2l+3)) A_l C_{l-1}^{3/2}(t)
                +SUM_{l odd} ((l+1) g_l/(2l+3)) B_l C_{l+1}^{3/2}(t)
   A_1 = log(R/rho);  A_l = (1-(rho/R)^{l-1})/(l-1), l>=3;  B_l = (1-(rho0/rho)^{l+4})/(l+4).
The l=1 interior term is the whole log: -(3 g_1/5) log(R/rho) = (M/2) log(R/rho)  (kappa_0 = 1/2).
"""
import numpy as np, math

def geg_table(t, LMAX):
    """C_l^{3/2}(t) for l=0..LMAX; t array. Recurrence l C_l = (2l+1) t C_{l-1} - (l+1) C_{l-2}."""
    t = np.atleast_1d(np.asarray(t, float))
    C = np.empty((LMAX+1, t.size))
    C[0] = 1.0
    if LMAX >= 1: C[1] = 3.0*t
    for l in range(2, LMAX+1):
        C[l] = ((2*l+1)*t*C[l-1] - (l+1)*C[l-2])/l
    return C

def Nl(l): return (l+1)*(l+2)/(l+1.5)

def bangbang_g(LMAX, nquad=4000):
    """g_l = h_l for w = -M sgn(z), M=1:  g_l N_l = -2 INT_0^{pi/2} C_l^{3/2}(cos th) sin^2 th dth."""
    x, wq = np.polynomial.legendre.leggauss(nquad)
    th = 0.5*(math.pi/2)*(x+1.0); wq = 0.5*(math.pi/2)*wq
    C = geg_table(np.cos(th), LMAX)
    I = C @ (wq*np.sin(th)**2)
    ls = np.arange(LMAX+1)
    g = -2.0*I/Nl(ls)
    g[ls % 2 == 0] = 0.0          # even l vanish by z-parity
    return g

def a_series(rho, phi, rho0, R, g, LMAX):
    """exact matched-shell series for a(rho,phi); rho0 <= rho <= R."""
    t = math.cos(phi)
    C = geg_table(np.array([t]), LMAX+1)[:, 0]
    ls = np.arange(1, LMAX+1, 2)
    A = np.where(ls == 1, math.log(R/rho), (1.0 - (rho/R)**np.maximum(ls-1, 1))/np.maximum(ls-1, 1))
    B = (1.0 - (rho0/rho)**(ls+4))/(ls+4)
    inner = -((ls+2)*g[ls]/(2*ls+3)) * A * C[ls-1]
    outer =  ((ls+1)*g[ls]/(2*ls+3)) * B * C[ls+1]
    return inner + outer                       # per-l terms (caller sums / accelerates)

def cesaro(terms, ntail=None):
    """partial sums + Cesaro average of the last ntail partial sums (kills the l^-2 oscillation)."""
    P = np.cumsum(terms)
    if ntail is None: ntail = max(4, len(P)//4)
    return float(P[-1]), float(P[-ntail:].mean())
