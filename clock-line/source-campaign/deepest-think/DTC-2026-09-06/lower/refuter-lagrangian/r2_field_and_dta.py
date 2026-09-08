#!/usr/bin/env python3
"""r2: my own 5D Gegenbauer solve for the bang-bang plateau, used to test three things the
attempt asserts:
  (a) a(0,0) = (M/2) L  and the material-point offset a(rho0,10deg) - (M/2)L  (~0.2168 claimed)
  (b) whether a(x,0) >= 0 pointwise  (the attempt's LEMMA 2 is only proved AT THE ORIGIN)
  (c) D_t a(0)|_{t=0} for the TRUE Euler dynamics, against the attempt's model value +M^2 L^2/8.
      Uses  d_t a(0) = (3/8pi^2) \int eta [ u.grad_5 W + 2 a W ] dy_5,   W = (-z)|y|^-5,
      where the second term is div_5 u = 2a  (the 5D lift is NOT divergence free -- dropping it
      flips the sign of the answer).
All quadrature is mine; nothing is imported from the attempt.
"""
import numpy as np, json
from numpy.polynomial.legendre import leggauss

M = 1.0

def gegen(lmax, t):
    """C_l^{3/2}(t) and dC/dt for l=0..lmax, shape (lmax+1, len(t))"""
    C = np.zeros((lmax+2, t.size)); D = np.zeros((lmax+2, t.size))
    C[0] = 1.0; C[1] = 3*t; D[0] = 0.0; D[1] = 3.0
    for n in range(1, lmax+1):
        C[n+1] = (2*(n+1.5-1)*t*C[n] - (n+2*1.5-2)*C[n-1])/(n+1)
        D[n+1] = (2*(n+0.5)*(t*D[n]+C[n]) - (n+1)*D[n-1])/(n+1)
    return C[:lmax+1], D[:lmax+1]

def norm_l(l):  # \int_-1^1 (1-t^2) [C_l^{3/2}]^2 dt
    return (l+1)*(l+2)/(l+1.5)

def coeffs(lmax, ngl=20000):
    """h_l for f(t) = rho*eta = -M sgn(t)/sqrt(1-t^2) (bang-bang plateau)."""
    x, w = leggauss(ngl)
    C, _ = gegen(lmax, x)
    f = -M*np.sign(x)/np.sqrt(1-x**2)
    h = np.zeros(lmax+1)
    for l in range(1, lmax+1, 2):
        h[l] = np.sum(w*f*C[l]*(1-x**2))/norm_l(l)
    return h

def radial(l, h_l, rho0, R, rho):
    """psi_l(rho), psi_l'(rho) for the four-region C^1 matched solution.
    Source h_l rho^-1 C_l on (rho0,R); Delta_5(rho^g C_l) = [g(g+3)-l(l+3)] rho^{g-2} C_l."""
    gp, gm = float(l), -(l+3.0)
    if l == 1:
        P  = lambda r: -h_l*r*np.log(r)/5.0
        Pp = lambda r: -h_l*(np.log(r)+1.0)/5.0
    else:
        d = l*(l+3.0)-4.0
        P  = lambda r: h_l*r/d
        Pp = lambda r: h_l/d*np.ones_like(r)
    # region I  rho<rho0 : A rho^gp ; II: P + B rho^gp + Cc rho^gm ; III: D rho^gm
    # unknowns A,B,Cc,D from C^1 at rho0 and R
        # build 4x4
    r0, r1 = rho0, R
    Mx = np.zeros((4,4)); rhs = np.zeros(4)
    # at r0: A r0^gp - B r0^gp - Cc r0^gm = P(r0)
    Mx[0] = [ r0**gp, -r0**gp, -r0**gm, 0]; rhs[0] = P(np.array([r0]))[0] if l==1 else P(r0)
    Mx[1] = [ gp*r0**(gp-1), -gp*r0**(gp-1), -gm*r0**(gm-1), 0]
    rhs[1] = (Pp(np.array([r0]))[0] if l==1 else Pp(np.array([r0]))[0])
    # at R: B R^gp + Cc R^gm + P(R) - D R^gm = 0
    Mx[2] = [0, r1**gp, r1**gm, -r1**gm]; rhs[2] = -(P(np.array([r1]))[0] if l==1 else P(r1))
    Mx[3] = [0, gp*r1**(gp-1), gm*r1**(gm-1), -gm*r1**(gm-1)]
    rhs[3] = -(Pp(np.array([r1]))[0])
    A,B,Cc,D = np.linalg.solve(Mx, rhs)
    rho = np.atleast_1d(rho)
    psi = np.zeros_like(rho); dpsi = np.zeros_like(rho)
    i1 = rho < r0; i2 = (rho>=r0)&(rho<=r1); i3 = rho>r1
    psi[i1]  = A*rho[i1]**gp;            dpsi[i1] = A*gp*rho[i1]**(gp-1)
    psi[i2]  = P(rho[i2]) + B*rho[i2]**gp + Cc*rho[i2]**gm
    dpsi[i2] = Pp(rho[i2]) + B*gp*rho[i2]**(gp-1) + Cc*gm*rho[i2]**(gm-1)
    psi[i3]  = D*rho[i3]**gm;            dpsi[i3] = D*gm*rho[i3]**(gm-1)
    return psi, dpsi, (A,B,Cc,D)

class Field:
    def __init__(self, L, lmax=401, rho0=1.0):
        self.rho0, self.R, self.L, self.lmax = rho0, rho0*np.exp(L), L, lmax
        self.h = coeffs(lmax)
        self.alpha1 = None
        self.sol = {l: radial(l, self.h[l], self.rho0, self.R, np.array([1.0]))[2]
                    for l in range(1, lmax+1, 2)}
        self.alpha1 = self.sol[1][0]
    def a_origin(self):
        return -3.0*self.alpha1                      # psi -> alpha1 rho C_1 = 3 alpha1 z
    def eval(self, rho, t):
        """psi, dpsi/drho, dpsi/dt summed over modes"""
        rho = np.atleast_1d(np.asarray(rho, float)); t = np.atleast_1d(np.asarray(t, float))
        C, Dc = gegen(self.lmax, t)
        psi = np.zeros((rho.size, t.size)); dr = np.zeros_like(psi); dt = np.zeros_like(psi)
        for l in range(1, self.lmax+1, 2):
            p, dp, _ = radial(l, self.h[l], self.rho0, self.R, rho)
            psi += np.outer(p, C[l]); dr += np.outer(dp, C[l]); dt += np.outer(p, Dc[l])
        return psi, dr, dt
    def a_uz(self, rho, t):
        psi, dr, dt = self.eval(rho, t)
        RH = rho[:,None]; T = t[None,:]
        # a = -d_z psi = -( t dpsi/drho + (1-t^2)/rho dpsi/dt )
        a = -(T*dr + (1-T**2)/RH*dt)
        # d_r psi = sinphi dpsi/drho - (t sinphi/rho) dpsi/dt ; r = rho sinphi
        s = np.sqrt(np.maximum(1-T**2, 0.0))
        dpsi_dr = s*dr - (T*s/RH)*dt
        uz = 2*psi + (RH*s)*dpsi_dr
        ur = (RH*s)*a
        return a, ur, uz, psi

out = {}
for L in (np.log(64), np.log(256), np.log(1024), np.log(4096)):
    F = Field(L, lmax=401)
    a0 = F.a_origin()
    out.setdefault('a_origin', []).append([float(np.exp(L)), float(a0), float(a0-0.5*L)])
# material offset at rho=rho0, phi=10deg
F = Field(np.log(4096), lmax=801)
t10 = np.cos(np.deg2rad(10.0))
a,_,_,_ = F.a_uz(np.array([1.0]), np.array([t10]))
out['material_offset_rho0_phi10'] = float(a[0,0] - 0.5*F.L)

# (b) is a >= 0 everywhere?  scan the shell and outside
rr = np.exp(np.linspace(-0.5, F.L+0.5, 240))
tt = np.cos(np.linspace(0.001, np.pi-0.001, 181))
A,_,_,_ = F.a_uz(rr, tt)
imin = np.unravel_index(np.argmin(A), A.shape)
out['min_a_over_domain'] = {'a_min': float(A[imin]),
    'rho/rho0': float(rr[imin[0]]), 'phi_deg': float(np.rad2deg(np.arccos(tt[imin[1]]))),
    'a_origin_for_scale': float(F.a_origin()), 'frac_of_domain_with_a_negative': float((A<0).mean())}

# (c) D_t a(0)|_0 for true Euler, all modes, at several L
def dta_origin(F, nr=1200, nphi=1200):
    lg = np.linspace(0.0, F.L, nr); rho = np.exp(lg)
    phi = np.linspace(1e-7, np.pi-1e-7, nphi); t = np.cos(phi)
    a, ur, uz, _ = F.a_uz(rho, t)
    RH = rho[:,None]; T = t[None,:]; S = np.sqrt(np.maximum(1-T**2,0))
    r = RH*S; z = RH*T
    eta = -M*np.sign(T)/r
    W  = -z*RH**-5
    dWdr = 5*r*z*RH**-7
    dWdz = -RH**-5 + 5*z**2*RH**-7
    integ = eta*(ur*dWdr + uz*dWdz + 2*a*W)
    # dy_5 = 2pi^2 r^3 dr dz ; (r,z)=(rho,phi): dr dz = rho drho dphi ; drho = rho dlogrho
    dens = integ*(r**3)*(RH**2)
    val = np.trapz(np.trapz(dens, phi, axis=1), lg)
    return 0.75*val          # (3/8pi^2)*2pi^2 = 3/4
res = []
for Lv in (np.log(256), np.log(1024), np.log(4096), np.log(65536)):
    Fv = Field(Lv, lmax=301)
    d = dta_origin(Fv)
    res.append([float(np.exp(Lv)), float(Lv), float(d), float(d/(Lv**2)), 0.125])
out['dta_origin'] = res
print(json.dumps(out, indent=1))
json.dump(out, open('r2_results.json','w'), indent=1)
