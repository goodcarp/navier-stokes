#!/usr/bin/env python3
"""G1 — exact (sympy) verification of the axisymmetric-NS identities this run leans on.
Exit 0 iff every check is identically zero. Controls: two deliberately wrong forms must be NONZERO.
Conventions: cylindrical (r,theta,z); u = u^r e_r + u^theta e_theta + u^z e_z; axisymmetric.
Stokes stream function psi: u^r = -psi_z / r, u^z = psi_r / r (divergence-free by construction).
q = u^theta/r, Gamma = r u^theta, eta = omega^theta/r, a = u^r/r, L5 = d_rr + (3/r)d_r + d_zz.
"""
import sympy as sp, sys, hashlib
r,z,t,nu,alpha = sp.symbols('r z t nu alpha', positive=True)
psi = sp.Function('psi')(r,z,t); v = sp.Function('v')(r,z,t); p = sp.Function('p')(r,z,t)
ur = -sp.diff(psi,z)/r; uz = sp.diff(psi,r)/r
Dt  = lambda f: sp.diff(f,t) + ur*sp.diff(f,r) + uz*sp.diff(f,z)
Lap = lambda f: sp.diff(f,r,2) + sp.diff(f,r)/r + sp.diff(f,z,2)
L5  = lambda f: sp.diff(f,r,2) + 3*sp.diff(f,r)/r + sp.diff(f,z,2)
Er  = Dt(ur) - v**2/r + sp.diff(p,r) - nu*(Lap(ur) - ur/r**2)
Eth = Dt(v) + ur*v/r - nu*(Lap(v) - v/r**2)
Ez  = Dt(uz) + sp.diff(p,z) - nu*Lap(uz)
q = v/r; Gam = r*v; a = ur/r
om = sp.diff(ur,z) - sp.diff(uz,r); eta = om/r
checks = {}
checks['DIV-FREE']      = sp.diff(r*ur,r)/r + sp.diff(uz,z)
checks['SWIRL-q  : D_t q - nu L5 q + 2 a q  == Eth/r'] = Dt(q) - nu*L5(q) + 2*a*q - Eth/r
checks['SWIRL-Gam: D_t Gam - nu(Lap - 2/r d_r)Gam == r Eth'] = Dt(Gam) - nu*(Lap(Gam) - 2*sp.diff(Gam,r)/r) - r*Eth
lhs = r*(Dt(eta) - nu*L5(eta) - sp.diff(q**2,z))
rhs = sp.diff(Er,z) - sp.diff(Ez,r)
checks['ETA-EQ   : r(D_t eta - nu L5 eta - d_z q^2) == d_z Er - d_r Ez'] = lhs - rhs
# weighted swirl functional identities (pointwise divergence forms), F generic
F = sp.Function('F'); G = sp.Symbol('G')  # F(Gamma)
Fg = F(Gam)
transport = sp.diff(ur*r**alpha*Fg, r) + sp.diff(uz*r**alpha*Fg, z) - r**alpha*(ur*sp.diff(Fg,r)+uz*sp.diff(Fg,z) + (alpha-1)*a*Fg)
checks['WEIGHT-TRANSPORT: div(b r^a F) == r^a[b.grad F + (a-1) a F]'] = transport
# viscous identity: generic power F = Gamma^m with symbolic m (monomials span analytic F)
m = sp.Symbol('m', positive=True)
Fm = Gam**m
Fm_p = m*Gam**(m-1); Fm_pp = m*(m-1)*Gam**(m-2)
visc_lhs = r**alpha * Fm_p * (Lap(Gam) - 2*sp.diff(Gam,r)/r)
visc_rhs = (sp.diff(r**alpha*Fm_p*sp.diff(Gam,r), r) + sp.diff(r**alpha*Fm_p*sp.diff(Gam,z), z)
            - r**alpha*Fm_pp*(sp.diff(Gam,r)**2 + sp.diff(Gam,z)**2)
            - sp.diff((alpha+1)*r**(alpha-1)*Fm, r) + (alpha**2-1)*r**(alpha-2)*Fm)
checks['WEIGHT-VISC: r^a F\'(G)(Lap-2/r d_r)G == div(...) - r^a F\'\'|grad G|^2 - d_r((a+1) r^{a-1} F) + (a^2-1) r^{a-2} F'] = visc_lhs - visc_rhs
# axis facts with the smooth parametrization u^r = r A(s,z), u^z = W(s,z), u^theta = r Q(s,z), s = r^2,
# taken as Taylor polynomials in s to order 2 (sufficient for every r->0 limit below)
A0,A1,A2,W0,W1,W2,Q0,Q1,Q2 = [sp.Function(n)(z) for n in ['A0','A1','A2','W0','W1','W2','Q0','Q1','Q2']]
sA = r**2
urA = r*(A0 + A1*sA + A2*sA**2); uzA = W0 + W1*sA + W2*sA**2; vA = r*(Q0 + Q1*sA + Q2*sA**2)
divA = sp.diff(r*urA,r)/r + sp.diff(uzA,z)
checks['AXIS-DIV: div|_{r=0} == 2 A(0,z) + W_z(0,z)'] = sp.limit(sp.simplify(divA), r, 0) - (2*A0 + sp.diff(W0,z))
omz = sp.diff(r*vA, r)/r
checks['AXIS-OMZ: omega^z(0,z) == 2 Q(0,z)'] = sp.limit(sp.simplify(omz), r, 0) - 2*Q0
qA = vA/r
checks['AXIS-L5q: (L5 q)(0,z) == 8 Q_s(0,z) + Q_zz(0,z)'] = sp.limit(sp.simplify(L5(qA)), r, 0) - (8*Q1 + sp.diff(Q0,z,2))
# a(0,z) = -(1/2) d_z u^z(0,z) is a CONSEQUENCE of div-free (not of the Taylor parametrization): check it
# with the stream function psi = r^2 (P0(z) + P1(z) r^2), divergence-free by construction
P0,P1 = sp.Function('P0')(z), sp.Function('P1')(z)
psiT = r**2*(P0 + P1*r**2)
urT = -sp.diff(psiT,z)/r; uzT = sp.diff(psiT,r)/r
checks['AXIS-STRETCH: a(0,z) == -(1/2) d_z u^z(0,z) for div-free fields'] = sp.limit(sp.simplify(urT/r), r, 0) + sp.Rational(1,2)*sp.limit(sp.simplify(sp.diff(uzT,z)), r, 0)
# controls that MUST fail
controls = {}
controls['CTRL-1 (wrong sign in swirl eq: +2aq -> -2aq) must be NONZERO'] = Dt(q) - nu*L5(q) - 2*a*q - Eth/r
controls['CTRL-2 (L5 replaced by Lap in eta eq) must be NONZERO'] = r*(Dt(eta) - nu*Lap(eta) - sp.diff(q**2,z)) - rhs
ok = True
for k,e in checks.items():
    val = sp.simplify(sp.expand(e))
    good = (val == 0)
    ok &= good
    print(('PASS ' if good else 'FAIL ') + k + ('' if good else f'  -> residual: {val}'))
for k,e in controls.items():
    val = sp.simplify(sp.expand(e))
    good = (val != 0)
    ok &= good
    print(('PASS ' if good else 'FAIL ') + k + ('' if good else '  -> control did NOT fire'))
print('SCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())
sys.exit(0 if ok else 1)
