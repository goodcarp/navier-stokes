#!/usr/bin/env python3
"""G3 — exact: how much positive strain integral ('transport charge') does a steadily translating axisymmetric
structure deposit at a FIXED axis point as it passes?  Claim: for u^z(0,z,t) = U f(z - U t) with f -> 0 at +-inf,
   ∫_{-inf}^{inf} a(0,z,t) dt = 0   and   ∫ a_+(0,z,t) dt = (1/2) * (total decrease of f) = (1/2) max f  (unimodal f),
using a(0,z) = -(1/2) d_z u^z(0,z) (G1 AXIS-STRETCH). Worked exactly on Hill's spherical vortex.
Hill's vortex (lab frame, translation speed U, radius R, omega^theta = A r inside): derived here from the Stokes
stream function, not quoted. Controls: a wrong-frame profile (co-moving instead of lab) must give a DIFFERENT toll;
the identity must fail for a profile with f(+inf) != f(-inf).
"""
import sympy as sp, sys, hashlib
r,z,R,U,A,zeta,t,z0 = sp.symbols('r z R U A zeta t z0', positive=True)
# co-moving-frame interior stream function of Hill's vortex: psi = (A/10) r^2 (R^2 - r^2 - z^2)
psi_in = sp.Rational(1,10)*A*r**2*(R**2 - r**2 - z**2)
ur = -sp.diff(psi_in,z)/r; uz = sp.diff(psi_in,r)/r
# it must carry vorticity omega^theta = A r  (E^2 psi = -omega r with the Stokes operator)
om = sp.diff(ur,z) - sp.diff(uz,r)
chk_vort = sp.simplify(om - A*r)
# the sphere is a streamline (psi=0 on r^2+z^2=R^2)
chk_sphere = sp.simplify(psi_in.subs(r**2, R**2 - z**2))
# matching to exterior potential flow past a sphere (co-moving): psi_out = -(U/2) r^2 (1 - R^3/(r^2+z^2)^(3/2)); tangential velocity continuous fixes U = 2 A R^2/15
psi_out = -sp.Rational(1,2)*U*r**2*(1 - R**3/(r**2+z**2)**sp.Rational(3,2))
uz_in_axis  = sp.simplify(sp.cancel(uz)).subs(r,0)                     # co-moving axis velocity inside
uz_out_axis = sp.limit(sp.simplify(sp.cancel(sp.diff(psi_out,r)/r)), r, 0) # co-moving axis velocity outside (z>0)
uz_in_eq  = sp.simplify(uz.subs({r:R, z:0}))                          # co-moving tangential (axial) velocity at the equator, inside
uz_out_eq = sp.simplify((sp.diff(psi_out,r)/r).subs({r:R, z:0}))       # same, outside
U_match = sp.solve(sp.Eq(uz_in_eq, uz_out_eq), U)[0]
chk_U = sp.simplify(U_match - 2*A*R**2/15)
# lab-frame axis profile f(zeta) = u^z_lab(0,zeta)/U, zeta = z - U t (vortex moving in +z)
f_in  = sp.simplify((uz_in_axis + U)/U).subs(A, 15*U/(2*R**2))
f_out = sp.simplify((uz_out_axis + U)/U)
f = sp.Piecewise((f_in, sp.Abs(zeta) < R), (f_out.subs(z, sp.Abs(zeta)), True)).subs(z, zeta)
f_center = sp.simplify(f_in.subs(z,0)); f_edge_in = sp.simplify(f_in.subs(z,R)); f_edge_out = sp.simplify(f_out.subs(z,R))
# strain at the fixed axis point z0 as the structure passes: a(t) = -(1/2) d_z u^z = -(U/2) f'(z0 - U t)
fz = sp.diff(f_in, z)  # inside |zeta|<R ; outside use f_out
# total signed toll: (1/2) [f(-inf) - f(+inf)]  by substitution zeta = z0 - U t, dt = -d zeta / U
toll_signed = sp.Rational(1,2)*(sp.limit(f_out.subs(z,-zeta).subs(zeta,-sp.oo) if False else f_out.subs(z, sp.oo), z, sp.oo) - f_out.subs(z, sp.oo))
toll_signed = sp.simplify(toll_signed)
# positive toll: (1/2) * ∫ (f')_- d zeta = (1/2) * total decrease of f.  f rises from 0 at -inf to f_center at 0 then falls to 0: decrease = f_center.
# compute directly: on zeta in (0,R): f_in decreasing? and on zeta>R: f_out = (R/zeta)^3 decreasing.
dec_in  = sp.simplify(f_in.subs(z,0) - f_in.subs(z,R))
dec_out = sp.simplify(f_out.subs(z,R) - sp.limit(f_out, z, sp.oo))
toll_pos = sp.simplify(sp.Rational(1,2)*(dec_in + dec_out))
mono_in  = sp.simplify(sp.diff(f_in,z))   # must be <= 0 on (0,R)
mono_out = sp.simplify(sp.diff(f_out,z))  # must be <= 0 on (R, inf)
print('Hill vortex derived: omega^theta - A r =', chk_vort, '| psi on sphere =', chk_sphere, '| U - 2AR^2/15 =', chk_U)
print('lab-frame axis profile: f(0) =', f_center, ' f(R-) =', f_edge_in, ' f(R+) =', f_edge_out, ' f(zeta>R) =', f_out)
print("f' inside (0,R):", mono_in, "  f' outside:", mono_out)
print('signed toll  ∫ a dt          =', toll_signed)
print('positive toll ∫ a_+ dt        =', toll_pos, '  (= (1/2) max f = ', sp.Rational(1,2)*f_center, ')')
# controls
ctrl_frame = sp.simplify(sp.Rational(1,2)*(uz_in_axis/U).subs(A,15*U/(2*R**2)).subs(z,0))  # co-moving 'toll' would be (1/2)(3/2)=3/4 != 5/4
g = 1 + sp.tanh(zeta)/3  # profile with unequal tails (2/3 at -inf, 4/3 at +inf): signed toll must be nonzero
ctrl_tail = sp.Rational(1,2)*(sp.limit(g, zeta, -sp.oo) - sp.limit(g, zeta, sp.oo))
ok = (chk_vort==0) and (chk_sphere==0) and (chk_U==0) and (toll_signed==0) and (toll_pos==sp.Rational(5,4)) and (f_edge_in==f_edge_out==1) \
     and (ctrl_frame != toll_pos) and (ctrl_tail != 0) and (mono_in.subs(z,R/2) < 0) and (mono_out.subs(z,2*R) < 0)
print('CTRL co-moving-frame toll (must differ):', ctrl_frame, '| CTRL unequal-tail signed toll (must be nonzero):', ctrl_tail)
print('RESULT', 'PASS' if ok else 'FAIL', '| SCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())
sys.exit(0 if ok else 1)
