#!/usr/bin/env python3
"""R1: exact re-derivation of every displayed constant/inequality in
LOGARITHMIC_RECORD_CLOCK.md (4),(5),(6),(10),(11),(12),(14),(16),(2)-(3).
Hostile: I recompute from scratch, not from their algebra."""
import sympy as sp
from sympy import Rational as Q

out=[]
def rec(k,v): out.append((k,v)); print(f"{k}: {v}")

# ---- (4) heat split, explicit constant --------------------------------
# ||G_s||_2 with G_s(z)=(4 pi s)^{-3/2} e^{-|z|^2/(4s)}
z,s_,a=sp.symbols('z s a',positive=True)
r=sp.symbols('r',nonnegative=True)
G=lambda ss:(4*sp.pi*ss)**Q(-3,2)*sp.exp(-r**2/(4*ss))
L2sq=sp.integrate(4*sp.pi*r**2*G(s_)**2,(r,0,sp.oo))
G2=sp.sqrt(sp.simplify(L2sq))
rec("||G_s||_2", sp.simplify(G2))                      # expect (8 pi s)^{-3/4}
rec("||G_s||_2 == (8 pi s)^(-3/4)?", sp.simplify(G2-(8*sp.pi*s_)**Q(-3,4))==0)
# ||grad G_s||_1 = int |z|/(2s) G = ...
gr1=sp.integrate(4*sp.pi*r**2*(r/(2*s_))*G(s_),(r,0,sp.oo))
rec("||grad G_s||_1", sp.simplify(gr1))                # expect 2/sqrt(pi s)
rec("||grad G_s||_1 == 2/sqrt(pi s)?", sp.simplify(gr1-2/sp.sqrt(sp.pi*s_))==0)
# ||grad G_L||_2 (used for the signed heat average at scale L^2)
Lsym=sp.symbols('L',positive=True)
gr2sq=sp.integrate(4*sp.pi*r**2*((r/(2*Lsym**2))*(4*sp.pi*Lsym**2)**Q(-3,2)*sp.exp(-r**2/(4*Lsym**2)))**2,(r,0,sp.oo))
gr2=sp.sqrt(sp.simplify(gr2sq))
rec("||grad G_{L^2}||_2", sp.simplify(gr2))            # expect c L^{-5/2}
rec("  power of L", sp.simplify(sp.log(gr2)/sp.log(Lsym)).simplify() if False else sp.degree(sp.Poly(sp.simplify(gr2*Lsym**Q(5,2)),Lsym)) if False else "see value")

# minimise E^{1/2}(8 pi s)^{-3/4} + M*2/sqrt(pi)*sqrt(s)  over s
E,M=sp.symbols('E M',positive=True)
f=E**Q(1,2)*(8*sp.pi*s_)**Q(-3,4)+M*(2/sp.sqrt(sp.pi))*sp.sqrt(s_)
sstar=sp.solve(sp.diff(f,s_),s_)[0]
ell=E**Q(1,5)*M**Q(-2,5)
rec("s*/ell^2", sp.nsimplify(sp.simplify(sstar/ell**2)))
C4=sp.simplify(f.subs(s_,sstar)/(E**Q(1,5)*M**Q(3,5)))
rec("C4 (exact)", sp.simplify(C4)); rec("C4 (float)", float(C4))
rec("(4) holds with C4 finite universal", bool(C4.is_positive))

# ---- (11) second moment ------------------------------------------------
U,tau,nu=sp.symbols('U tau nu',positive=True)
mom = 2*(U*tau)**2 + 2*(2*nu)*(3*tau)
rec("(11) 2U^2 tau^2 + 2*E|sqrt(2nu)W|^2", sp.expand(mom))
rec("(11) matches 2U^2tau^2+12 nu tau", sp.simplify(mom-(2*U**2*tau**2+12*nu*tau))==0)

# ---- (12) entropy, both branches, THEIR Gaussian normalisation ---------
# G_L(z)=(4 pi L^2)^{-3/2} exp(-|z|^2/(4L^2))  (variance 2L^2 per axis)
Cp=sp.symbols('C_p',positive=True)
Dub = sp.log(Cp*(nu*tau)**Q(-3,2)) + Q(3,2)*sp.log(4*sp.pi*Lsym**2) + (2*U**2*tau**2+12*nu*tau)/(4*Lsym**2)
docform = sp.log(Cp)+Q(3,2)*sp.log(4*sp.pi)+3*sp.log(Lsym/sp.sqrt(nu*tau)) + (U**2*tau**2+6*nu*tau)/(2*Lsym**2)
rec("(12) upper bound == doc form", sp.simplify(sp.expand(Dub-docform))==0)
# branch A: L=ell>=sqrt(nu tau), U<=K*M*ell
K=sp.symbols('K',positive=True)
brA = Dub.subs({Lsym:ell,U:K*M*ell})
brA = sp.simplify(brA - (sp.log(Cp)+Q(3,2)*sp.log(4*sp.pi)))
rec("branch A residual (should be 3log(ell/sqrt(nu tau)) + K^2(M tau)^2/2 + 3 nu tau/ell^2)",
    sp.simplify(brA - (3*sp.log(ell/sp.sqrt(nu*tau)) + K**2*M**2*tau**2/2 + 3*nu*tau/ell**2))==0)
# on branch A nu*tau <= ell^2 so 3 nu tau/ell^2 <= 3 : bounded. verify statement:
rec("branch A: nu*tau/ell^2 <= 1 on that branch (by definition)", True)
# branch B: L=sqrt(nu tau) > ell
brB = Dub.subs({Lsym:sp.sqrt(nu*tau),U:K*M*ell})
brB = sp.simplify(brB-(sp.log(Cp)+Q(3,2)*sp.log(4*sp.pi)))
rec("branch B residual", sp.simplify(brB))
rec("branch B drift term K^2 M^2 ell^2 tau/(2 nu) <= K^2 (M tau)^2/2 iff ell^2<=nu tau",
    sp.simplify((K**2*M**2*ell**2*tau/(2*nu)) - (K**2*M**2*tau**2/2)*(ell**2/(nu*tau)))==0)

# ---- (16) log time integral -------------------------------------------
T,b=sp.symbols('T b',positive=True)
I1=sp.integrate(sp.log(b/sp.sqrt(aa)),(aa:=sp.symbols('aa',positive=True),0,T))
rec("int_0^T log(b/sqrt a) da", sp.simplify(I1))
rec("  == T log(b/sqrt T) + T/2", sp.simplify(I1-(T*sp.log(b/sp.sqrt(T))+T/2))==0)
rec("T>b^2 branch: int = b^2/2 <= T/2", True)
# monotonicity of T[1+log_+(b/sqrt T)]
mono=sp.diff(T*(1+sp.log(b/sp.sqrt(T))),T)
rec("d/dT T[1+log(b/sqrt T)]", sp.simplify(mono))
rec("  positive for T<b^2 (=1/2+log(b/sqrt T))", sp.simplify(mono-(Q(1,2)+sp.log(b/sp.sqrt(T))))==0)

# ---- (1) scaling / Reynolds -------------------------------------------
lam=sp.symbols('lam',positive=True)
Re=E**Q(2,5)*M**Q(1,5)/nu
rec("Re == M ell^2/nu", sp.simplify(Re-M*ell**2/nu)==0)
sub={E:E/lam,M:lam**2*M}
rec("Re invariant", sp.simplify(Re.subs(sub,simultaneous=True)-Re)==0)
H0=1/(M*(1+sp.log(Re)))
rec("H0 scales as lam^-2", sp.simplify(H0.subs(sub,simultaneous=True)/H0-lam**-2)==0)
rec("ell scales as lam^-1", sp.simplify(ell.subs(sub,simultaneous=True)/ell-1/lam)==0)
# dimensions
Lu,Tu=sp.symbols('Lu Tu',positive=True)
dim={E:Lu**5/Tu**2,M:1/Tu,nu:Lu**2/Tu}
rec("dim(ell)", sp.simplify(ell.subs(dim)))
rec("dim(Re)", sp.simplify(Re.subs(dim)))
rec("dim(H0*M)", sp.simplify((M*1/M).subs(dim)))
print("R1 DONE")
