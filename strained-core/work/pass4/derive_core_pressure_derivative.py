#!/usr/bin/env python3
"""Exact radial reduction of the core pressure derivative; no time simulation."""
import sympy as sp
r,z,s,t=sp.symbols('r z s t',real=True)
p=sp.symbols('p0:4')
f=sp.symbols('f0:4')
h0=sp.symbols('h0 h0p h0pp')
h2=sp.symbols('h2 h2p h2pp')
h4=sp.symbols('h4 h4p h4pp')
a0=sp.symbols('a0 a0p a0pp')
a2=sp.symbols('a2 a2p a2pp')
families=[p,f,h0,h2,h4,a0,a2]

def D(expr,v):
 return sp.diff(expr,v)+2*v*sum(q[i+1]*sp.diff(expr,q[i]) for q in families for i in range(len(q)-1))

def avg(expr,weighted=False):
 poly=sp.Poly(sp.expand(expr),r,z)
 out=0
 for (nr,nz),coef in poly.terms():
  assert nr%2==0 and nz%2==0,(nr,nz)
  out+=coef*s**((nr+nz)//2)*(1-t)**(nr//2)*t**(nz//2)
 if weighted: out*=3*t-1
 poly=sp.Poly(sp.expand(out),t)
 return sp.factor(sum(coef/sp.Integer(2*n[0]+1) for n,coef in poly.terms()))

def trace_poly(ar,az,vr,vz):
 return D(ar,r)*D(vr,r)+(ar/r)*(vr/r)+D(az,z)*D(vz,z)+D(ar,z)*D(vz,r)+D(az,r)*D(vr,z)

def trace_swirl(aw,vw):
 return -(aw/r)*D(vw,r)-(vw/r)*D(aw,r)

sr=-r*(p[0]+2*z*z*p[1]); sz=2*z*(p[0]+r*r*p[1]); w=r*f[0]
assert sp.expand(D(sr,r)+sr/r+D(sz,z))==0
R2=r*r+z*z
H2=(3*z*z-R2)/2
H4=(35*z**4-30*R2*z*z+3*R2**2)/8
pS=h0[0]+h2[0]*H2+h4[0]*H4
pW=a0[0]+a2[0]*H2
cr=sp.expand(-(sr*D(sr,r)+sz*D(sr,z)))
cz=sp.expand(-(sr*D(sz,r)+sz*D(sz,z)))
swr= r*f[0]**2; swz=sp.Integer(0)
cvw=sp.expand(-(sr*D(w,r)+sz*D(w,z)+sr*w/r))

# Cubic pressure derivative sources: b^3, b*Omega^2.
g300_local=2*trace_poly(sr,sz,cr,cz)
g300_pressure=-2*trace_poly(sr,sz,D(pS,r),D(pS,z))
gWb_local=2*trace_poly(sr,sz,swr,swz)+2*trace_swirl(w,cvw)
gWb_pressure=-2*trace_poly(sr,sz,D(pW,r),D(pW,z))

if __name__=='__main__':
 for name,expr in [('A300local',g300_local),('A300pressure',g300_pressure),('AWblocal',gWb_local),('AWbpressure',gWb_pressure)]:
  out=avg(sp.expand(expr),True)
  print(name,'=',sp.collect(sp.expand(out),list(h0)+list(h2)+list(h4)+list(a0)+list(a2)))
 print('cvw/r=',sp.factor(cvw/r))
 print('g300(0)=',sp.expand(sp.expand(g300_local+g300_pressure).subs({r:0,z:0,p[0]:1,p[1]:0,p[2]:0,p[3]:0})))
 print('gWb(0)=',sp.expand(sp.expand(gWb_local+gWb_pressure).subs({r:0,z:0,p[0]:1,p[1]:0,p[2]:0,p[3]:0,f[0]:1,f[1]:0,f[2]:0,f[3]:0})))

# Exact compact-support polynomial proxies (Euler cubic coefficients only).
# psi=(1-s)^n on [0,1], zero outside; not flat at the center and not C-infinity
# at the outer endpoint. They are used only to test cutoff dependence of the
# continuous inviscid cubic functionals, never the flat-core viscosity claim.
def weighted_integral(poly,power):
 return sp.expand(sum(c*s**(k[0]+power+1)/(k[0]+power+1) for k,c in sp.Poly(sp.expand(poly),s).terms()))

def inside_pressures(ps):
 ps1=sp.diff(ps,s)
 fs=sp.expand(ps+sp.Rational(2,3)*s*ps1)
 H0p=-ps**2-sp.Rational(4,5)*s*ps*ps1+sp.Rational(4,15)*s*s*ps1**2
 H2=-sp.Rational(2,7)*ps**2+sp.Rational(4,21)*s**(-sp.Rational(5,2))*weighted_integral(ps1**2,sp.Rational(7,2))
 H4=(sp.Rational(48,35)*weighted_integral(ps*ps1,sp.Rational(7,2))+sp.Rational(16,15)*weighted_integral(ps1**2,sp.Rational(9,2)))*s**(-sp.Rational(9,2))+sp.Rational(32,21)*sp.integrate(ps1**2,(s,s,1))
 A0p=fs**2/3
 A2=-sp.Rational(1,3)*s**(-sp.Rational(5,2))*weighted_integral(fs**2,sp.Rational(3,2))
 return tuple(sp.expand(e) for e in (fs,H0p,H2,H4,A0p,A2))

def proxy(n):
 ps=(1-s)**n
 fs,H0p,H2,H4,A0p,A2=inside_pressures(ps)
 # Check every scalar Poisson ODE independently of coefficient evaluation.
 ps1=sp.diff(ps,s); ps2=sp.diff(ps,s,2)
 GS0=6*ps**2+16*s*ps*ps1+sp.Rational(16,5)*s*s*ps*ps2-sp.Rational(8,15)*s*s*ps1**2-sp.Rational(32,15)*s**3*ps1*ps2
 GS2=sp.Rational(8,21)*(21*ps*ps1+6*s*ps*ps2+2*s*ps1**2-4*s*s*ps1*ps2)
 GS4=-sp.Rational(64,35)*(3*ps*ps2-13*ps1**2-2*s*ps1*ps2)
 GW0=-2*fs**2-sp.Rational(8,3)*s*fs*sp.diff(fs,s)
 GW2=sp.Rational(8,3)*fs*sp.diff(fs,s)
 for val in [4*s*sp.diff(H0p,s)+6*H0p+GS0,
             4*s*sp.diff(H2,s,2)+14*sp.diff(H2,s)+GS2,
             4*s*sp.diff(H4,s,2)+22*sp.diff(H4,s)+GS4,
             4*s*sp.diff(A0p,s)+6*A0p+GW0,
             4*s*sp.diff(A2,s,2)+14*sp.diff(A2,s)+GW2]:
  assert sp.expand(val)==0
 subs={p[i]:sp.diff(ps,s,i) for i in range(4)}
 subs.update({f[i]:sp.diff(fs,s,i) for i in range(4)})
 subs[h0[2]]=sp.diff(H0p,s)
 subs[a0[2]]=sp.diff(A0p,s)
 for syms,val in [(h2,H2),(h4,H4),(a2,A2)]:
  subs.update({syms[i]:sp.diff(val,s,i) for i in range(3)})
 Q300=sp.expand(avg(sp.expand(g300_local+g300_pressure),True).subs(subs)/(2*s))
 QWb=sp.expand(avg(sp.expand(gWb_local+gWb_pressure),True).subs(subs)/(2*s))
 T300=sp.Rational(20,7)+sp.integrate(Q300,(s,0,1))
 TWb=sp.Rational(52,15)+sp.integrate(QWb,(s,0,1))
 return T300,TWb

if __name__=='__main__':
 for n in [4,5,6]:
  vals=proxy(n)
  print('Exact polynomial proxy n =',n,': t300,tWb =',vals,'decimal =',tuple(float(v) for v in vals))

# Verify the general radial Green formulas and their reduced angular source,
# with cumulative integrals treated by their exact fundamental derivatives.
M,N,V,T=sp.symbols('M N V T')
def radial_D(e):
 return (sp.diff(e,s)+sum(p[i+1]*sp.diff(e,p[i]) for i in range(3))
         +s**sp.Rational(7,2)*p[1]**2*sp.diff(e,M)
         +s**sp.Rational(7,2)*p[0]*p[1]*sp.diff(e,N)
         +s**sp.Rational(9,2)*p[1]**2*sp.diff(e,V)
         -p[1]**2*sp.diff(e,T))
H0p=-p[0]**2-sp.Rational(4,5)*s*p[0]*p[1]+sp.Rational(4,15)*s*s*p[1]**2
H2=-sp.Rational(2,7)*p[0]**2+sp.Rational(4,21)*s**(-sp.Rational(5,2))*M
H4=(sp.Rational(48,35)*N+sp.Rational(16,15)*V)*s**(-sp.Rational(9,2))+sp.Rational(32,21)*T
GS2=sp.Rational(8,21)*(21*p[0]*p[1]+6*s*p[0]*p[2]+2*s*p[1]**2-4*s*s*p[1]*p[2])
GS4=-sp.Rational(64,35)*(3*p[0]*p[2]-13*p[1]**2-2*s*p[1]*p[2])
assert sp.expand(4*s*radial_D(radial_D(H2))+14*radial_D(H2)+GS2)==0
assert sp.expand(4*s*radial_D(radial_D(H4))+22*radial_D(H4)+GS4)==0
moment_subs={h0[2]:radial_D(H0p)}
for syms,val in [(h2,H2),(h4,H4)]:
 moment_subs.update({syms[0]:val,syms[1]:radial_D(val),syms[2]:radial_D(radial_D(val))})
actual=sp.expand((avg(sp.expand(g300_local+g300_pressure),True)/(2*s)).subs(moment_subs))
local=(-sp.Rational(824,35)*p[0]**2*p[1]-sp.Rational(3664,245)*s*p[0]**2*p[2]
       -sp.Rational(64,35)*s*s*p[0]**2*p[3]-sp.Rational(2848,105)*s*p[0]*p[1]**2
       -sp.Rational(1216,105)*s*s*p[0]*p[1]*p[2]
       +sp.Rational(128,105)*s**3*p[0]*p[1]*p[3]+sp.Rational(128,105)*s**3*p[0]*p[2]**2
       +sp.Rational(544,105)*s*s*p[1]**3-sp.Rational(64,105)*s**3*p[1]**2*p[2])
nested=(-sp.Rational(128,735)*s**(-sp.Rational(5,2))*(3*p[1]-2*s*p[2])*M
        +sp.Rational(3072,245)*s**(-sp.Rational(5,2))*p[2]*N
        +sp.Rational(1024,105)*s**(-sp.Rational(5,2))*p[2]*V
        -sp.Rational(256,735)*(63*p[0]+96*s*p[1]+20*s*s*p[2])*T)
assert sp.expand(actual-local-nested)==0

# Integration-by-parts reductions used in equation (2). I,J denote the
# integrals in the note; R=int Q*p'^2 and Z=int s^(-7/2)*p'*K.
I,J,R,Z=sp.symbols('I J R Z')
int_local=(sp.Rational(1728,245)*I+sp.Rational(64,5)*J+sp.Rational(200,49))
int_M=-sp.Rational(256,735)*J
int_T=-sp.Rational(256,735)*(7*R+56*I+20*J)
int_NV=sp.Rational(3072,245)*(-I-sp.Rational(7,9)*J+sp.Rational(5,2)*Z)
closed=(sp.Rational(340,49)-sp.Rational(18368,735)*I-sp.Rational(448,105)*J
        -sp.Rational(256,105)*R+sp.Rational(1536,49)*Z)
assert sp.expand(sp.Rational(20,7)+int_local+int_M+int_T+int_NV-closed)==0

def proxy_closed(n):
 ps=(1-s)**n; ps1=sp.diff(ps,s)
 Iv=sp.integrate(s*ps*ps1**2,(s,0,1)); Jv=sp.integrate(s*s*ps1**3,(s,0,1))
 Qv=weighted_integral(ps,sp.Integer(0))
 Kv=weighted_integral(ps*ps1,sp.Rational(7,2))+sp.Rational(7,9)*weighted_integral(ps1**2,sp.Rational(9,2))
 Rv=sp.integrate(sp.expand(Qv*ps1**2),(s,0,1))
 Zv=sp.integrate(sp.expand(s**(-sp.Rational(7,2))*ps1*Kv),(s,0,1))
 return (closed.subs({I:Iv,J:Jv,R:Rv,Z:Zv}),
         sp.Rational(244,105)-sp.Rational(64,15)*Iv-sp.Rational(64,45)*Jv)

if __name__=='__main__':
 # Independent exact rational targets from the full cubic-source calculation.
 expected={4:(sp.Rational(51414292252,11712375675),sp.Rational(9188,4725)),
           5:(sp.Rational(16529009732,3720401685),sp.Rational(23908,12285))}
 for n,targets in expected.items():
  assert all(sp.simplify(a-b)==0 for a,b in zip(proxy_closed(n),targets))
 print('PASS: general radial Poisson ODEs, full angular cubic-source reduction,')
 print('      exact integration-by-parts arithmetic, and independent rational proxy checks.')
 print('No PDE evolution, interval sign bound, or inherited profile closure was verified.')
