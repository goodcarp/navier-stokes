#!/usr/bin/env python3
"""Exact algebra for actual initial covariance/torque jets; no time integration."""
import sympy as s

r,z = s.symbols("r z", positive=True, real=True)
m,k,lam,A,c,nu = s.symbols("m k lam A c nu", real=True)
ph = s.symbols("ph", real=True)
V = s.Function("V")(r)
pp = s.Function("pp")(r,z)
G = s.Function("G")(r)
wr = -lam*m*s.sin(ph)/r
wt = lam*k*s.sin(ph)
dr = lambda f: s.diff(f,r)+k*s.diff(f,ph)
dt = lambda f: m*s.diff(f,ph)

def avg(expr):
    poly = s.Poly(s.expand(expr),s.sin(ph),s.cos(ph))
    val = 0
    for (i,j),coef in poly.terms():
        if (i,j)==(0,0):
            val += coef
        elif (i,j) in ((2,0),(0,2)):
            val += coef/2
        elif (i+j) in (1,3) or (i,j)==(1,1):
            pass
        else:
            raise AssertionError((i,j))
    return s.expand(val)

R = -lam**2*m*k/(2*r)
assert s.simplify(avg(wr*wt)-R)==0
Rrr = lam**2*m**2/(2*r**2)
Rtt = lam**2*k**2/2
assert s.simplify(avg(wr**2)-Rrr)==0
assert s.simplify(avg(wt**2)-Rtt)==0
assert avg(wr**2*wt)==0 and avg(wt**3)==0
crossgrad = avg(dr(wr)*dr(wt)+dt(wr)*dt(wt)/r**2)
assert s.simplify(crossgrad +
       lam**2*m*k*(k*k+m*m/(r*r))/(2*r))==0
Lap0=lambda f:s.diff(f,r,2)+s.diff(f,r)/r
viscR=nu*(Lap0(R)-2*R/r**2-2*crossgrad)
assert s.simplify(viscR+nu*R*(2*k*k+(2*m*m+1)/r**2))==0

# Direct vector Laplacian gives the same covariance viscosity.
scalarlap=lambda f:dr(dr(f))+dr(f)/r+dt(dt(f))/r**2
Lwr=scalarlap(wr)-wr/r**2-2*dt(wt)/r**2
Lwt=scalarlap(wt)-wt/r**2+2*dt(wr)/r**2
assert s.simplify(avg(nu*(wt*Lwr+wr*Lwt))-viscR)==0

# Reynolds equation local terms, retaining mean advection.
Rtime = (-c*r*s.diff(R,r)-Rrr*(s.diff(V,r)+V/r)
         -R*c+(2*V*Rtt-c*r*R)/r+viscR)
explicit = (-c*R+lam**2*(k*k*V/r-m*m*(s.diff(V,r)+V/r)/(2*r*r))
            -nu*R*(2*k*k+(2*m*m+1)/r**2))
assert s.simplify(Rtime-explicit)==0

T0=lam**2*m*k/(2*r)
assert s.simplify(-s.diff(r*r*R,r)/r-T0)==0
Tnonp=-s.diff(r*r*explicit,r)/r
Texpected=(-c*T0+nu*T0*((2*m*m+1)/r**2-2*k*k)
  +lam**2*(m*m*(s.diff(V,r,2)+s.diff(V,r)/r-V/r**2)/(2*r)
          -k*k*(s.diff(V,r)+V/r)))
assert s.simplify(Tnonp-Texpected)==0

# Collapse both pressure covariance fluxes with the FULL Poisson equation.
eps=s.exp(-s.I*k*r)
Rp=-eps*(s.I*k*s.diff(pp,r)+m*m*pp/r**2)
Zp=-eps*s.I*k*s.diff(pp,z)
Tp=-s.diff(r*r*Rp,r)/r-r*s.diff(Zp,z)
source=-(s.diff(pp,r,2)+s.diff(pp,r)/r+s.diff(pp,z,2)-m*m*pp/r**2)
collapsed=eps*((k*k*r+m*m/r+s.I*k)*s.diff(pp,r)-s.I*k*r*source)
assert s.simplify(Tp-collapsed)==0

# The local source contribution on the plateau has real part 2 k² V'.
source_plateau=2*(s.I*k*s.diff(V,r)-k*k*V-m*m*s.diff(V,r)/r)/r
v0,v1=s.symbols("v0 v1",real=True)
local_source_term=(-s.I*k*r*source_plateau).subs({s.diff(V,r):v1,V:v0})
assert s.simplify(s.expand(local_source_term).as_real_imag()[0]-2*k*k*v1)==0

# Fixed-circle second jet at an initial mean maximum.
Dstar=lambda f:s.diff(f,r,2)-s.diff(f,r)/r
L0=lambda f:-c*r*s.diff(f,r)+nu*Dstar(f)
L2=s.expand(L0(L0(G)))
expectedL2=(c*c*(r*r*s.diff(G,r,2)+r*s.diff(G,r))
 -2*c*nu*r*s.diff(G,r,3)
 +nu*nu*(s.diff(G,r,4)-2*s.diff(G,r,3)/r
         +3*s.diff(G,r,2)/r**2-3*s.diff(G,r)/r**3))
assert s.simplify(L2-expectedL2)==0
assert s.simplify(L0(T0)-c*T0-3*nu*T0/r**2)==0
VasG=G/r
shear=lam**2*(m*m*(s.diff(VasG,r,2)+s.diff(VasG,r)/r-VasG/r**2)/(2*r)
              -k*k*(s.diff(VasG,r)+VasG/r))
assert s.simplify(shear-lam**2*(m*m*s.diff(G,r,2)/(2*r*r)
                   -(k*k/r+m*m/(2*r**3))*s.diff(G,r)))==0

# The exact mean-shear production in global fluctuation energy.
C=s.Function("C")(r)
e=s.Function("e")(r)
q=s.Function("q")(z)
qv=s.Function("qv")(z)
Rfull=-lam**2*m*k*e**2*q**2/(2*r)
Vfull=A*r*C*qv
production=-2*s.pi*r*Rfull*(s.diff(Vfull,r)-Vfull/r)
assert s.simplify(production -
 s.pi*lam**2*m*k*A*r*s.diff(C,r)*e**2*q**2*qv)==0
print("PASS: covariance production and viscosity, pressure/axial flux collapse, "
      "fixed-circle Gtt, and exact global fluctuation-energy production.")
