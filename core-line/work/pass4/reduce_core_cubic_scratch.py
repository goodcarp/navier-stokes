import importlib.util,sympy as q
spec=importlib.util.spec_from_file_location('m','work/pass4/derive_core_pressure_derivative.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
s=m.s;p=m.p
M,N,V,T=q.symbols('M N V T')
# M=int tau^(7/2)p'^2; N=int tau^(7/2)p p'; V=int tau^(9/2)p'^2; T=int_s^infty p'^2
D=lambda e:q.diff(e,s)+sum(p[i+1]*q.diff(e,p[i]) for i in range(3))+s**q.Rational(7,2)*p[1]**2*q.diff(e,M)+s**q.Rational(7,2)*p[0]*p[1]*q.diff(e,N)+s**q.Rational(9,2)*p[1]**2*q.diff(e,V)-p[1]**2*q.diff(e,T)
H0p=-p[0]**2-q.Rational(4,5)*s*p[0]*p[1]+q.Rational(4,15)*s*s*p[1]**2
H2=-q.Rational(2,7)*p[0]**2+q.Rational(4,21)*s**(-q.Rational(5,2))*M
H4=q.Rational(48,35)*s**(-q.Rational(9,2))*N+q.Rational(16,15)*s**(-q.Rational(9,2))*V+q.Rational(32,21)*T
subs={m.h0[2]:D(H0p)}
for syms,val in [(m.h2,H2),(m.h4,H4)]: subs.update({syms[0]:val,syms[1]:D(val),syms[2]:D(D(val))})
A=q.expand(m.avg(q.expand(m.g300_local+m.g300_pressure),True)/(2*s))
A=q.expand(A.subs(subs))
for sym in [M,N,V,T]:print(sym,':',q.factor(A.coeff(sym)))
print('LOCAL',q.expand(A.subs({M:0,N:0,V:0,T:0})))
# ODE tests.
GS4=-q.Rational(64,35)*(3*p[0]*p[2]-13*p[1]**2-2*s*p[1]*p[2])
GS2=q.Rational(8,21)*(21*p[0]*p[1]+6*s*p[0]*p[2]+2*s*p[1]**2-4*s*s*p[1]*p[2])
assert q.expand(4*s*D(D(H4))+22*D(H4)+GS4)==0
assert q.expand(4*s*D(D(H2))+14*D(H2)+GS2)==0
print('GENERAL H2,H4 ODE PASS')
