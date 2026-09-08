from fractions import Fraction as F
import sympy as s

c=F(7,5); m=4; k=20; nu=F(1,1000)
rlo=F(63,400); rhi=F(63,250)
Slo=2*m*780*F(1,2)/rhi
Shi=2*m*1020/rlo
assert Slo==F(260000,21) and Shi==F(1088000,21)
assert c*k==28
assert Slo-c*k==F(259412,21)
assert Shi-c*k==F(1087412,21)
assert F(k)/(Shi-c*k)==F(105,271853)>F(1,2600)
assert F(k)/(Slo-c*k)==F(105,64853)<F(1,600)

D=lambda t,S,r:nu*((k*k+(m/r)**2)*t+k*S*t*t+S*S*t**3/3)
assert D(F(1,600),Shi,rlo)==F(22357,2551500)<F(1,100)
assert D(F(9,1000),Shi,rlo)==F(20550697,27562500)<1
assert D(F(3,100),Slo,rhi)==F(536479,330750)>1
assert F(k)*rlo/(2*m)==F(63,160)
assert F(k)*rhi/(2*m)==F(63,100)

# Gradient at Gamma_r=0 and conversion to a cylindrical frame rotating at Omega.
cc,Om,kr,kt=s.symbols('c Omega kr kt',real=True)
A=s.Matrix([[cc,-Om],[-Om,cc]])
J=s.Matrix([[0,-1],[1,0]])
kap=s.Matrix([kr,kt])
assert -A.T*kap-Om*J*kap==s.Matrix([-cc*kr+2*Om*kt,-cc*kt])

# Pressure coefficient in the affine Kelvin polarization preserves kappa.a=0.
q=s.Matrix(s.symbols('q1:4')); a=s.Matrix(s.symbols('a1:4'))
AA=s.Matrix(3,3,s.symbols('A0:9'))
adot=-AA*a+2*q*(q.dot(AA*a))/q.dot(q)
qdot=-AA.T*q
assert s.simplify(qdot.dot(a)+q.dot(adot))==0
assert F(290649,8750)-F(53997,2000)==F(435297,70000)
print('PASS: initial shear, mean gradient, frozen winding/damping, and Kelvin pressure constraint')
