#!/usr/bin/env python3
"""Exact horizontal D kernels, conservative signs, and energy prefactors."""
import sympy as s
from fractions import Fraction as F
r,z,R,t=s.symbols('r z R t',positive=True)
aj,bj=s.symbols('aj bj',real=True)
R2=r*r+z*z
N=1/s.sqrt(R2)  # Common 1/(4*pi) omitted from Q and J throughout.
Qrr=-s.diff(N,z,2,r,2)
Qtt=-s.diff(N,z,2,r)/r
Qrz=-s.diff(N,z,3,r)
Pdot=lambda f:r*s.diff(f,r)-2*z*s.diff(f,z)
QP_r=r*Qrr-2*z*Qrz
assert s.simplify(-2*s.diff(QP_r,r)-Pdot(Qrr)+2*Qrr+3*Pdot(Qrr))==0
assert s.simplify(-2*QP_r/r-Pdot(Qtt)+2*Qtt+3*Pdot(Qtt))==0

J=(s.Rational(16,7)/R2**s.Rational(3,2))*s.legendre(2,z/s.sqrt(R2))
J+=(s.Rational(36,7)/R2**s.Rational(3,2)+aj*R2**2+bj/R2**s.Rational(5,2))*s.legendre(4,z/s.sqrt(R2))
q=35*t*t-35*t+4
p=231*t**3-336*t*t+119*t-6
h=(21*t*t-14*t+1)/(5*t-1)
j=52920*t**3-77840*t*t+28280*t-1504
expected={
 'rr':(-s.Rational(18,7)+s.Rational(5,8)*p/(R*R*q),
       -j/(28*q)+aj*R**7*(7*t-3)/q+s.Rational(5,4)*bj*p/(R*R*q)),
 'tt':(-s.Rational(18,7)+s.Rational(5,8)*h/(R*R),
       s.Rational(24,7)+(5*bj/(4*R*R)-10)*h+aj*R**7)}
def convert(f):
    return s.factor(s.powdenest(s.simplify(f).subs({r:R*s.sqrt(1-t),z:R*s.sqrt(t)}),force=True))
for name,Q,Jhess in [('rr',Qrr,s.diff(J,r,2)),('tt',Qtt,s.diff(J,r)/r)]:
    core=-s.Rational(18,7)*Q+s.diff(Q,z,2)/24
    pump=-3*Pdot(Q)-2*Jhess
    assert s.simplify(convert(core/Q)-expected[name][0])==0
    assert s.simplify(convert(pump/Q)-expected[name][1])==0

x=s.symbols('x',nonnegative=True)
assert s.expand(q.subs(t,1-x)-(4-35*x+35*x*x))==0
assert s.expand(p.subs(t,1-x)-(8-140*x+357*x*x-231*x**3))==0
assert s.expand(j.subs(t,1-x)-(1856-31360*x+80920*x*x-52920*x**3))==0
assert s.simplify(s.diff(h,t)-3*(35*t*t-14*t+3)/(5*t-1)**2)==0

xmax=F(1,90); Rmin2=F(17,5)**2
assert F(49,4673)<xmax
qmin=4-35*xmax; pmax=8+357*xmax*xmax
jmin=1856-31360*xmax-52920*xmax**3
assert qmin>0 and jmin>0 and 7*(1-xmax)-3>0
pump_rr_upper=-jmin/112+F(25,4)*pmax/(Rmin2*qmin)
core_rr_upper=-F(18,7)+F(5,8)*pmax/(Rmin2*qmin)
assert pump_rr_upper<-12 and core_rr_upper<-F(12,5)
hm=s.factor(h.subs(t,s.Rational(89,90)))
assert hm>s.Rational(19,10)
assert F(25,4)/Rmin2-10<0
pump_tt_upper=F(24,7)+(F(25,4)/Rmin2-10)*F(19,10)
core_tt_upper=-F(18,7)+F(5,4)/Rmin2
assert pump_tt_upper<-12 and core_tt_upper<-F(12,5)
assert F(12,5)+F(7,5)*12==F(96,5)

# Exact horizontal Laplacian decomposition used in the separated energy.
th=s.symbols('theta',real=True);mm,kk=s.symbols('m k',real=True)
ee=s.Function('e')(r);qq=s.Function('q')(z);phase=mm*th+kk*r
potential=ee*qq*s.cos(phase)
lap_h=s.diff(potential,r,2)+s.diff(potential,r)/r+s.diff(potential,th,2)/r**2
expected_lap=qq*((s.diff(ee,r,2)+s.diff(ee,r)/r-(kk*kk+mm*mm/r**2)*ee)*s.cos(phase)
                -kk*(2*s.diff(ee,r)+ee/r)*s.sin(phase))
assert s.simplify(lap_h-expected_lap)==0
viscosity_upper=F(6,1000)*F(5,17)**5*(F(12,5)*120000+F(6764,75)*F(160,3))
assert viscosity_upper<4
assert F(5093339,210000)>4
print('PASS: both local advection contractions reduce to -3 P.dot.grad Q.')
print('PASS: exact full-pressure horizontal core/pump kernel ratios.')
print('PASS: conservative horizontal bounds core<-12/5 and pump<-12.')
print('PASS: separated horizontal Laplacian and rational viscosity prefactors.')
print('Coarse viscous upper bound if certified R1<120000:',viscosity_upper)
print('Conclusion with interval certificate: D_nu < 4-(96/5) C_w <4.')
