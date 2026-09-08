#!/usr/bin/env python3
"""Independent exact general-envelope initial NS jet algebra; no time solver."""
import sympy as s

r,z=s.symbols("r z",positive=True,real=True)
m,k,lam,A,c,nu=s.symbols("m k lam A c nu",real=True)
phase=s.symbols("phase",real=True)
e=s.Function("e")(r)
V=s.Function("V")(r)
g=s.Function("g")(r)
p=s.Function("p")(r,z)
ep=s.diff(e,r); epp=s.diff(e,r,2); eppp=s.diff(e,r,3)
sn,co=s.sin(phase),s.cos(phase)
wr=-lam*m*e*sn/r
wt=-lam*ep*co+lam*k*e*sn
dr=lambda f:s.diff(f,r)+k*s.diff(f,phase)
dt=lambda f:m*s.diff(f,phase)

def avg(expr):
    poly=s.Poly(s.expand(expr),sn,co)
    out=0
    for (i,j),coef in poly.terms():
        if (i,j)==(0,0):out+=coef
        elif (i,j) in ((2,0),(0,2)):out+=coef/2
        elif (i+j) in (1,3) or (i,j)==(1,1):pass
        else:raise AssertionError((i,j))
    return s.expand(out)

B0=lam**2*m*k/2
F=e**2
R=-B0*F/r
T=B0*(s.diff(F,r)+F/r)
assert s.simplify(avg(wr*wt)-R)==0
assert s.simplify(-s.diff(r*r*R,r)/r-T)==0
assert avg(wr**2*wt)==avg(wt**3)==0
Rrr=lam**2*m*m*e*e/(2*r*r)
Rtt=lam**2*(ep*ep+k*k*e*e)/2
assert s.simplify(avg(wr**2)-Rrr)==0
assert s.simplify(avg(wt**2)-Rtt)==0

# Complete vector-Laplacian covariance derivative.
crossgrad=avg(dr(wr)*dr(wt)+dt(wr)*dt(wt)/r**2)
cross_expected=B0/r*(e*epp-2*ep**2+2*e*ep/r-(k*k+m*m/r**2)*e**2)
assert s.simplify(crossgrad-cross_expected)==0
Lap0=lambda f:s.diff(f,r,2)+s.diff(f,r)/r
visc_from_cov=nu*(Lap0(R)-2*R/r**2-2*crossgrad)
Hnu=-4*e*epp+2*ep**2-2*e*ep/r+2*k*k*e**2+(2*m*m+1)*e**2/r**2
assert s.simplify(visc_from_cov-nu*B0*Hnu/r)==0
scalarlap=lambda f:dr(dr(f))+dr(f)/r+dt(dt(f))/r**2
Lwr=scalarlap(wr)-wr/r**2-2*dt(wt)/r**2
Lwt=scalarlap(wt)-wt/r**2+2*dt(wr)/r**2
assert s.simplify(avg(nu*(wt*Lwr+wr*Lwt))-visc_from_cov)==0

# Mean pump and shear production, keeping partial (not material) R_t.
pump=-c*r*s.diff(R,r)-2*c*R
assert s.simplify(pump-c*T)==0
shear=-Rrr*(s.diff(V,r)+V/r)+2*V*Rtt/r
shear_expected=lam**2*(V*(ep*ep+k*k*e*e)/r-
                     m*m*e*e*(s.diff(V,r)+V/r)/(2*r*r))
assert s.simplify(shear-shear_expected)==0

# General-envelope modal source directly from mixed advection.
b=e*s.exp(s.I*k*r)
Wr=s.I*m*b/r; Wt=-s.diff(b,r)
Cr=s.I*m*V*Wr/r-2*V*Wt/r
Ct=s.I*m*V*Wt/r+(s.diff(V,r)+V/r)*Wr
source=s.diff(Cr,r)+Cr/r+s.I*m*Ct/r
source_expected=2*(s.diff(V*s.diff(b,r),r)-m*m*s.diff(V,r)*b/r)/r
assert s.simplify(source-source_expected)==0

# Both pressure covariance fluxes reduce to p_r and the exact source.
phasefactor=s.exp(-s.I*k*r)
BB=ep-s.I*k*e
Rp=phasefactor*(BB*s.diff(p,r)-m*m*e*p/r**2)
Zp=phasefactor*BB*s.diff(p,z)
Tp=-s.diff(r*r*Rp,r)/r-r*s.diff(Zp,z)
ss=-(s.diff(p,r,2)+s.diff(p,r)/r+s.diff(p,z,2)-m*m*p/r**2)
DD=(k*k*r+m*m/r)*e-r*epp-ep+s.I*k*(e+2*r*ep)
assert s.simplify(Tp-phasefactor*(DD*s.diff(p,r)+r*BB*ss))==0

# Real explicit source term at the mean maximum V=r C, V'=-C.
E0,E1,E2,C0=s.symbols("E0 E1 E2 C0",real=True)
source_real=2*(s.diff(V,r)*ep+V*(epp-k*k*e)-m*m*s.diff(V,r)*e/r)/r
source_imag=2*k*(s.diff(V,r)*e+2*V*ep)/r
source_part=(r*(ep*source_real+k*e*source_imag)/2)
source_expected=(s.diff(V,r)*(ep**2+k*k*e**2-m*m*e*ep/r)
                 +V*ep*(epp+k*k*e))
assert s.simplify(source_part-source_expected)==0
at_max={s.diff(V,r):-C0,V:r*C0}
source_max=s.expand(source_part.subs(at_max))
shear_torque=-s.diff(r*r*shear_expected,r)/(r*lam**2)
shear_as_g=s.simplify(shear_torque.subs(V,g/r).doit())
shear_max=s.simplify(shear_as_g.xreplace({s.diff(g,r):s.Integer(0)}))
assert s.simplify(shear_max-
 (m*m*e**2*s.diff(g,r,2)/(2*r*r)-2*g*ep*(epp+k*k*e)/r))==0
g2=s.symbols("g2",real=True)
# Keep g'' as its own local jet when replacing g by r²C0.
combined=s.expand(shear_max.xreplace({g:r*r*C0,s.diff(g,r,2):g2})+source_max)
expected_combined=m*m*e**2*g2/(2*r*r)-C0*(ep**2+k*k*e**2+r*ep*(epp+k*k*e)-m*m*e*ep/r)
assert s.simplify(combined-expected_combined)==0

# Add mean transport/diffusion of T to the covariance torque derivative.
Dstar=lambda f:s.diff(f,r,2)-s.diff(f,r)/r
pump_total=-c*r*s.diff(T,r)-s.diff(r*r*pump,r)/r
assert s.simplify(pump_total+2*c*B0*(2*s.diff(F,r)+r*s.diff(F,r,2)))==0
visc_total=nu*Dstar(T)-s.diff(r*r*visc_from_cov,r)/r
Jnu=(6*(ep*epp+e*eppp)+6*e*epp/r-4*k*k*e*ep-2*k*k*e**2/r
      -(4*m*m+8)*e*ep/r**2+(2*m*m+4)*e**2/r**3)
assert s.simplify(visc_total-nu*B0*Jnu)==0

# Exact flat-envelope limit reproduces the earlier audited formula.
flat={e:1,ep:0,epp:0,eppp:0}
assert s.simplify(DD.subs(flat)-(k*k*r+m*m/r+s.I*k))==0
assert s.simplify(pump_total.subs(flat))==0
assert s.simplify(Jnu.subs(flat)-((2*m*m+4)/r**3-2*k*k/r))==0
print("PASS: arbitrary-envelope covariance/viscosity, full modal source and "
      "pressure-flux cancellation, maximum-circle shear/source combination, "
      "pump and viscous transport corrections, and exact flat-envelope limit.")
