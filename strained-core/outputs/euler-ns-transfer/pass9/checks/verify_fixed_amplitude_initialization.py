#!/usr/bin/env python3
"""Exact rational and off-neutral core algebra for the A=1020 datum."""
from fractions import Fraction as F
from pathlib import Path
import json
import sympy as s

X,L,q,T,d=s.symbols('X L q T d',real=True)
b1=2+d; O1=2
b2=-4*b1-(T-32)/2
O2=2*b1+4
beta2=s.expand(b2-O2-2*b1*O1+2*O1**2)
assert s.simplify(beta2+T/2+10*d)==0
Tup=102-s.Rational(23199,1000)*X+q*(7-s.Rational(96,5)*L+340)
lower=s.expand(-Tup/2-10*(X+q*L-s.Rational(204,35))/2)
expected=-s.Rational(153,7)+s.Rational(13199,2000)*X-q*s.Rational(347,2)+q*s.Rational(23,5)*L
assert s.simplify(lower-expected)==0
A=F(1020); x=A*A*F(113,20000000)
beta1=(x-F(204,35))/2
tupper=102-F(23199,1000)*x+F(347,16)
beta2lower=-F(153,7)+F(13199,2000)*x-F(347,32)
selected_G=F(1950,16)-F(19,1000)*A
assert beta1==F(17391,700000)>F(2484,100000)
assert tupper==-F(634112687,50000000)<-F(1268,100)
assert beta2lower==F(4264878809,700000000)>F(609,100)
assert selected_G==F(20499,200)>F(10249,100)
result=dict(status='PASS',selected_parameters=dict(A='1020',lam='1/4',nu='1/1000'),
    family_lam_interval=['1/8','1/4'],X_lower=str(x),beta_prime_lower=str(beta1),
    complete_pressure_T_upper=str(tupper),beta_second_lower=str(beta2lower),
    selected_mean_derivative_lower=str(selected_G),
    uniform_mean_derivative_lower='8871/800',uniform_fractional_fluctuation_energy_derivative_lower='12853/2205',
    scope='Actual initial identities and inequalities for an explicit amplitude; no full 3D time evolution, duration, inherited return or blowup.',
    dependencies=['pass5 all-amplitude T(H)<102-(23199/1000)H and C_v lower bound',
                  'pass8 leading torque, energy, D_L and E_L bounds',
                  'pass9 full off-neutral core identity and independent-amplitude mean-jet formula'])
here=Path(__file__).resolve().parent
(here/'fixed-amplitude-initialization-certificate.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
