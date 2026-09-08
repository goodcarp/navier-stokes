"""Independent cylindrical pressure-source and full feedback-jet checks."""
import sympy as s


def check(name,expression):
    value=s.simplify(expression)
    assert value==0,(name,value)
    print(f"PASS: {name}")


R2,Z2,radius2,mu2=s.symbols('R2 Z2 radius2 mu2',real=True)
p,p1,p2=s.symbols('p p1 p2',real=True)
# For the unit poloidal field:
# u_r=-r(psi+2z²psi'), u_z=2z(psi+r²psi').
# Direct cylindrical strain/curvature source, independent of the Cartesian checker.
ur_r=-(p+2*(R2+Z2)*p1+4*R2*Z2*p2)
ur_over_r=-(p+2*Z2*p1)
uz_z=2*(p+(R2+2*Z2)*p1+2*R2*Z2*p2)
cross=-16*R2*Z2*(3*p1+2*Z2*p2)*(2*p1+R2*p2)
g=s.expand(ur_r**2+ur_over_r**2+uz_z**2+cross)
check('poloidal divergence',ur_r+ur_over_r+uz_z)
g_sphere=s.expand(g.subs({R2:radius2*(1-mu2),Z2:radius2*mu2}))
weighted=s.Poly(s.expand((3*mu2-1)*g_sphere),mu2)
angular=sum(coef/s.Integer(2*power[0]+1) for power,coef in weighted.terms())
expected=s.Rational(16,105)*radius2*(-4*radius2**2*p1*p2+6*radius2*p*p2+2*radius2*p1**2+21*p*p1)
check('independent cylindrical strain angular average',angular-expected)

I=s.symbols('I',real=True)
pv=s.Rational(8,105)*(-4*(-I)+6*(s.Rational(1,2)-I)+2*I+21*(-s.Rational(1,2)))
check('cutoff-dependent integral cancels',pv+s.Rational(4,7))
check('full axial strain pressure including contact',-2+pv+s.Rational(18,7))

b,Om,bp,Omp,pzz,pzzp=s.symbols('b Om bp Omp pzz pzzp',real=True)
L=s.Matrix([[-b,-Om,0],[Om,-b,0],[0,0,2*b]])
Lt=s.Matrix([[-bp,-Omp,0],[Omp,-bp,0],[0,0,2*bp]])
source_time=2*s.trace(L*Lt)
tuned={bp:2*b*b,Omp:2*b*Om}
check('tuned pressure-source time derivative at origin',source_time.subs(tuned)-(24*b**3-8*b*Om**2))
check('tuned contact term',(-source_time/3).subs(tuned)-(-(24*b**3-8*b*Om**2)/3))

# Exact beta' numerator, whose value is zero on the tuned boundary.
eta_b,eta_w,eta_bp,eta_wp=s.symbols('eta_b eta_w eta_bp eta_wp',real=True)
beta=b/Om
beta_dot=(-4*b*b-pzz/2+eta_b-beta*eta_w)/Om
deriv=s.diff(beta_dot,b)*bp+s.diff(beta_dot,Om)*Omp+s.diff(beta_dot,pzz)*pzzp+s.diff(beta_dot,eta_b)*eta_bp+s.diff(beta_dot,eta_w)*eta_wp
tuning={**tuned,pzz:-8*b*b,eta_b:0,eta_w:0,eta_bp:0,eta_wp:0}
check('full tuned normalized-strain second derivative',deriv.subs(tuning)+(pzzp+32*b**3)/(2*Om))

CFp,Dp=s.symbols('CFp Dp',real=True)
CF=s.Rational(38,7)*b*b+s.Rational(2,5)*Om*Om
Fprime=CFp/Om**2-4*b*CF/Om**2
pressure_prime=-s.Rational(72,7)*b**3+s.Rational(8,5)*b*Om**2-CFp+Dp
beta_second=-(pressure_prime+32*b**3)/(2*Om)
check('exterior cone derivative needs pressure-profile defect',beta_second-(Om*Fprime/2-Dp/(2*Om)))
print('All checks passed. The full cubic pressure integral remains unevaluated; no sign for beta second derivative is asserted.')
