"""Exact algebra for the same-solution radius iteration audit."""
import sympy as s


def check(name, expression):
    value=s.simplify(expression)
    assert value==0,(name,value)
    print(f"PASS: {name}")


a=s.symbols('a',positive=True)
L0,L1,L2,U0,nu=s.symbols('L0 L1 L2 U0 nu',positive=True)
Re0=L0*U0/nu
Re1=Re0*(L0/L1)**a
Re2=Re1*(L1/L2)**a
check('exact two-step Reynolds telescoping',Re2-Re0*(L0/L2)**a)
check('exact matched velocity power',nu*Re2/L2-U0*(L0/L2)**(1+a))
check('velocity time exponent exceeds one half', (1+a)/(2+a)-s.Rational(1,2)-a/(2*(2+a)))

t,K,Q,Omega,D=s.symbols('t K Q Omega D',positive=True)
beta=(K+Q/(Omega**2*D))/2
q=1-beta*t*t/(a+2)
central_vorticity_factor=1+K*t*t/2
normalized=s.series(q**(a+2)*central_vorticity_factor,t,0,3).removeO()
check('normalized central rotation loses potential-energy fraction',normalized-(1-Q*t*t/(2*Omega**2*D)))

b11,b12,b21,b22=s.symbols('b11 b12 b21 b22',real=True)
B=s.Matrix([[b11,b12],[b21,b22]])
w=(b21-b12)/2
sigma=-s.trace(B)
check('two-by-two skew square identity',(B**2)[1,0]/2-(B**2)[0,1]/2-s.trace(B)*w)
eta_w,eta_s,pzz=s.symbols('eta_w eta_s pzz',real=True)
wdot=sigma*w+eta_w
sdot=-sigma*sigma-pzz+eta_s
ratio_dot=(sdot*w-sigma*wdot)/w**2
expected=(-2*sigma*sigma-pzz+eta_s-(sigma/w)*eta_w)/w
check('exact strain-to-rotation evolution',ratio_dot-expected)

R,W=s.symbols('R W',positive=True)
ratio_rhs=(-2*R**2*W**2-pzz+eta_s-R*eta_w)/W
check('lower ratio boundary condition',ratio_rhs.subs(R,0)-(-pzz+eta_s)/W)

q_scale=s.symbols('q_scale',positive=True)
sigma_e,w_e=s.symbols('sigma_e w_e',nonzero=True)
check('isotropic normalization cannot change strain ratio',
      (q_scale**(a+2)*sigma_e)/(q_scale**(a+2)*w_e)-sigma_e/w_e)
check('matched RMS normalization exponent',q_scale**(a+1)*q_scale-q_scale**(a+2))

print('All identities passed. A scalar radius rule is not a proved profile-class return.')
