from fractions import Fraction as F
import sympy as s
r,c,Ar,C,Cr=s.symbols('r c A C C_r',positive=True)
wr,wt=s.symbols('wr wt',real=True)
grad=s.Matrix([[c,-Ar*C,0],[Ar*(C+r*Cr),c,0],[0,0,-2*c]])
w=s.Matrix([wr,wt,0])
assert s.expand((w.T*grad*w)[0]-(c*(wr**2+wt**2)+Ar*r*Cr*wr*wt))==0
assert F(63,400)>F(7,50)
assert (F(7,25)/F(63,200))**2>F(5,8)
assert 2*2*F(9,20)/2==F(9,10)
work=4*20*780*F(63,800)*F(9,10)
energy=F(1,2)*F(6764,75)*F(12,5)
assert work/energy==F(552825,13528)
rate=F(14,5)+work/energy
assert rate==F(2953517,67640)>F(4366,100)
print('PASS: full initial strain contraction and rational fractional fluctuation-energy drain; no positive-time decay rate inferred.')
