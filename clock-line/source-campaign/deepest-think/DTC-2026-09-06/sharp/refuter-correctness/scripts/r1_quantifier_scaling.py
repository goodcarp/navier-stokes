"""
R1. Quantifier / scaling audit of the candidate TWO-SIDED DOUBLING CLOCK.

Candidate as written:
    c1/(M (1+log_+ Re_E))  <=  inf_data T_d  <=  c2/(M log Re_E)

with E = ||u||_2^2, M = ||omega||_inf, ell = E^{1/5} M^{-2/5}, Re_E = M ell^2 / nu,
T_d = first time ||omega||_inf reaches (3/2)||omega_0||_inf.

The middle term is a single number (an infimum over ALL smooth finite-energy data);
the two outer terms depend on a datum through M and Re_E.  This script settles what
the literal statement asserts, using the exact NS scaling symmetry.

Scaling: u_lambda(x,t) = lambda u(lambda x, lambda^2 t) solves NS with the SAME nu.

Everything below is computed, not asserted: E and M are evaluated by quadrature on an
explicit smooth compactly-decaying divergence-free field and on its rescalings.
"""
import numpy as np, json, hashlib, sys, time

t0=time.time()
out={}

# ---------------------------------------------------------------- explicit datum
# u = curl A,  A = (0, 0, x*y*exp(-|x|^2))  -> smooth, divergence free, finite energy.
# A3 = x y e^{-r2};  u = (d_y A3, -d_x A3, 0)
def field(X,Y,Z,lam=1.0):
    # u_lam(x) = lam * u(lam x)   (time slice t=0)
    x,y,z=lam*X,lam*Y,lam*Z
    r2=x*x+y*y+z*z; e=np.exp(-r2)
    A3=x*y*e
    dA3_dx=(y-2*x*x*y)*e
    dA3_dy=(x-2*y*y*x)*e
    dA3_dz=(-2*z*x*y)*e
    u1=dA3_dy; u2=-dA3_dx; u3=0.0*x
    # vorticity = curl u ; A3 only -> omega = curl curl A = grad(div A) - lap A = -lap A (div A = dA3/dz)
    # do it directly by finite difference-free analytic curl:
    # u1 = dA3/dy, u2 = -dA3/dx, u3 = 0
    # w1 = d_y u3 - d_z u2 = d_z dA3/dx
    # w2 = d_z u1 - d_x u3 = d_z dA3/dy
    # w3 = d_x u2 - d_y u1 = -(d_xx + d_yy) A3
    d2_zx=(-2*z)*(y-2*x*x*y)*e
    d2_zy=(-2*z)*(x-2*y*y*x)*e
    d2_xx=(-4*x*y-2*x*y+4*x*x*x*y)*e      # d/dx[(y-2x^2 y)e] = (-4xy)e + (y-2x^2y)(-2x)e
    d2_yy=(-4*x*y-2*x*y+4*y*y*y*x)*e
    w1=d2_zx; w2=d2_zy; w3=-(d2_xx+d2_yy)
    return (lam*u1,lam*u2,lam*u3),(lam*lam*w1,lam*lam*w2,lam*lam*w3)

def EM(lam,n=161,L=4.0):
    g=np.linspace(-L,L,n); h=g[1]-g[0]
    X,Y,Z=np.meshgrid(g,g,g,indexing='ij')
    u,w=field(X,Y,Z,lam)
    E=(u[0]**2+u[1]**2+u[2]**2).sum()*h**3
    M=np.sqrt(w[0]**2+w[1]**2+w[2]**2).max()
    return E,M

# divergence control (must be ~0)
def divcheck(lam=1.0,h=1e-5):
    P=np.array([0.31,-0.47,0.62])
    def u(p):
        return np.array(field(np.array([p[0]]),np.array([p[1]]),np.array([p[2]]),lam)[0]).ravel()
    d=0.0
    for i in range(3):
        e=np.zeros(3); e[i]=h
        d+=(u(P+e)[i]-u(P-e)[i])/(2*h)
    return d

print("CONTROL  div u at a generic point (must be ~0):",f"{divcheck():.3e}")
out['div_control']=float(divcheck())

lams=[0.5,1.0,2.0,3.0]
rows=[]
E1,M1=EM(1.0)
nu=1.0
for lam in lams:
    E,M=EM(lam)
    ell=E**0.2*M**-0.4
    Re=M*ell**2/nu
    rows.append(dict(lam=lam,E=E,M=M,ell=ell,Re=Re,
                     E_over_pred=E/(E1/lam), M_over_pred=M/(M1*lam**2)))
print("\nSCALING TABLE  (nu fixed; predicted E ~ lam^-1, M ~ lam^+2, Re_E invariant)")
print(f"{'lam':>5} {'E':>14} {'E*lam/E1':>10} {'M':>14} {'M/(lam^2 M1)':>13} {'Re_E':>14}")
for r in rows:
    print(f"{r['lam']:>5} {r['E']:>14.8e} {r['E_over_pred']:>10.7f} {r['M']:>14.8e} {r['M_over_pred']:>13.9f} {r['Re']:>14.10f}")
out['scaling_rows']=rows
Re_spread=max(r['Re'] for r in rows)/min(r['Re'] for r in rows)-1
print(f"\nRe_E relative spread over lambda in [0.5,3]: {Re_spread:.3e}   (exact invariance would be 0; residual is quadrature)")
out['Re_spread']=float(Re_spread)

# ------------------------------------------------- consequence for the literal claim
print("""
CONSEQUENCE.  M T_d and Re_E are BOTH scale invariant (M T_d because T_d ~ lam^-2).
So for any datum realising a pair (M, Re_E, T_d) the whole one-parameter family
{lambda} realises (lam^2 M, Re_E, lam^-2 T_d).  Hence at ANY fixed value of Re_E,

        inf over data of T_d  =  0        (let lambda -> infinity),

and the infimum over the whole admissible class is 0 a fortiori.""")
for lam in [1,10,100,1000]:
    print(f"   lam={lam:>5}:  M -> {lam**2:>9} M,  Re_E unchanged,  T_d -> {1.0/lam**2:.1e} T_d")
print("""
Therefore, read literally (M and Re_E free symbols, inf_data a number):
  * the LEFT inequality c1/(M(1+log_+Re_E)) <= inf_data T_d = 0 is FALSE for every c1>0;
  * the RIGHT inequality 0 <= c2/(M log Re_E) is TRUE but vacuous.
The statement is not false-because-of-fluid-mechanics; it is ill-quantified.
The only reading under which both halves say something is the per-Reynolds one:
    T(Re) := inf { M T_d : smooth finite-energy data with Re_E <= Re },
    claim  c1/(1+log Re) <= T(Re) <= c2/log Re.
Note T(Re) is non-increasing in Re, which the two-sided form requires.""")

# ------------------------------------------- arithmetic tying E_shell to log Re_E
C=0.172403978
print(f"\nCHECK of the exact-first-order seat's shell arithmetic:")
print(f"  Re_E = C^(2/5) M R^2/nu with C={C}:  C^(2/5) = {C**0.4:.9f},  log C^(2/5) = {np.log(C**0.4):+.7f}")
print(f"  quoted offset -0.7031660 ; computed {np.log(C**0.4):+.7f} ; abs diff {abs(np.log(C**0.4)+0.7031660):.2e}")
out['shell_offset']=float(np.log(C**0.4))
out['shell_offset_quoted']=-0.7031660

# c1/c2 book-keeping: doubling vs 3/2
print(f"  a=(M/4)(log Re_E+0.703): time to factor 2   = 4 ln2 /(M log Re_E) = {4*np.log(2):.7f}/(M log Re_E)")
print(f"                           time to factor 3/2 = 4 ln1.5/(M log Re_E) = {4*np.log(1.5):.7f}/(M log Re_E)")
print("  NOTE: the seat labels 4 ln 2 = 2.7726 as 'c1'.  In the candidate statement c1 is the")
print("  constant of the LOWER bound (Astra) and the construction supplies c2.  The shell number")
print("  is a c2, and for the (3/2) threshold in T_d it is 4 ln(3/2) = %.7f, not 4 ln 2." % (4*np.log(1.5)))
out['c2_shell_32']=float(4*np.log(1.5)); out['c2_shell_double']=float(4*np.log(2))

json.dump(out,open('r1_results.json','w'),indent=1)
print(f"\nelapsed {time.time()-t0:.1f}s  SHA256 {hashlib.sha256(open(__file__,'rb').read()).hexdigest()}")
