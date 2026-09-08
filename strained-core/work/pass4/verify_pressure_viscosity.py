#!/usr/bin/env python3
"""Independent second-jet viscosity coefficient for affine-core initial data.

Finite symbolic identities are assertions. Numerical comparisons are exploratory
quadrature diagnostics, NOT certified accuracy or PDE/sign verification.
The callable viscosity_from_second_jets(ev, amplitudes) can reuse a resolved
Evaluator without rebuilding profiles or solving pressure.
"""

import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'experiments'))
import numpy as np
import sympy as sp
from scipy.integrate import simpson


def symbolic_checks():
    x, y, z, r, xi = sp.symbols('x y z r xi', real=True)
    aa, bb, cc, dd, ee, ff = sp.symbols('aa bb cc dd ee ff', real=True)
    F = aa*z + bb*xi*z + cc*z**3
    W = dd + ee*xi + ff*z**2
    Fz = sp.diff(F, z)
    uz = 2*F + 2*xi*sp.diff(F, xi)
    cart = sp.Matrix([-x*Fz-y*W, -y*Fz+x*W, uz]).subs(xi, x*x+y*y)
    variables = (x, y, z)
    G = cart.jacobian(variables)
    assert sp.simplify(sp.trace(G)) == 0
    S_cart = sum(sp.trace(G.diff(v)**2) for v in variables)
    g = sp.trace(G**2)
    lapg = sum(sp.diff(g, v, 2) for v in variables)
    lapG = sum((G.diff(v, 2) for v in variables), sp.zeros(3))
    assert sp.simplify(2*sp.trace(G*lapG)-lapg+2*S_cart) == 0

    v = (-r*Fz).subs(xi, r*r)
    w = (r*W).subs(xi, r*r)
    q = uz.subs(xi, r*r)
    U = sp.Matrix([[sp.diff(v,r), -w/r, sp.diff(v,z)],
                   [sp.diff(w,r), v/r, sp.diff(w,z)],
                   [sp.diff(q,r), 0, sp.diff(q,z)]])
    J = sp.Matrix([[0,-1,0],[1,0,0],[0,0,0]])
    angular = sp.trace((J*U-U*J)**2)/r**2
    S_cyl = sp.trace(U.diff(r)**2)+sp.trace(U.diff(z)**2)+angular
    at_theta_zero = S_cart.subs({x:r,y:0})
    assert sp.simplify(S_cyl-at_theta_zero) == 0
    assert sp.simplify(angular) != 0, 'Test must detect omitted basis rotation.'
    print(json.dumps(dict(symbolic_checks='PASS',
        assertions=['divergence-free generic six-parameter axisymmetric polynomial',
                    'h_nu = Laplacian(g) - 2 S2',
                    'Cartesian S2 equals cylindrical derivatives plus commutator',
                    'angular commutator is nonzero in this test'])), flush=True)


def viscosity_from_second_jets(ev, amplitudes):
    """Return coefficient of nu in p_zz'(0), via -(4/5) int (S2)_2/R.

    Requires ev's compact profiles, exact affine inner ball R<.5, radial grid,
    second component derivatives, and Legendre degree two. No pressure solve.
    """
    U = ev.combine(amplitudes)
    G = ev.gradient(U)
    v, w, q = U[:,0]
    vr, wr, qr = U[:,1]
    vz, wz, qz = U[:,2]
    vrr, wrr, qrr = U[:,3]
    vrz, wrz, qrz = U[:,4]
    vzz, wzz, qzz = U[:,5]
    r = ev.r
    zero = np.zeros_like(v)
    Gr = np.array([[vrr, -wr/r+w/r**2, vrz],
                   [wrr, vr/r-v/r**2, wrz],
                   [qrr, zero, qrz]])
    Gz = np.array([[vrz, -wz/r, vzz],
                   [wrz, vz/r, wzz],
                   [qrz, zero, qzz]])
    J = np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,0.]])
    comm = np.einsum('ij,jkrm->ikrm',J,G)-np.einsum('ijrm,jk->ikrm',G,J)
    S2 = (np.einsum('ijrm,jirm->rm',Gr,Gr)
          +np.einsum('ijrm,jirm->rm',Gz,Gz)
          +np.einsum('ijrm,jirm->rm',comm,comm)/r**2)
    # This is exact for every allowed combination of ev's profiles. It prevents
    # subtraction noise near the axis from masquerading as an inner source.
    S2[ev.R < .5] = 0.
    degree2 = int(np.nonzero(ev.ells == 2)[0][0])
    source_l2 = S2 @ ev.project[:,degree2]
    return float(-4/5*simpson(source_l2/ev.R,x=ev.R))


def numerical_comparison(ev):
    names = ['strain_core','rotation_core','meridional_outer','swirl_outer']
    coefficients = [viscosity_from_second_jets(ev,np.eye(4)[j]) for j in range(4)]
    reference = [0.,0.,ev.dP,ev.dv]
    rows = []
    for name, value, other in zip(names,coefficients,reference):
        rows.append(dict(profile=name,second_jets=value,stress_or_exact=other,
                         absolute_difference=value-other,
                         relative_difference=None if other == 0 else (value-other)/abs(other)))
    amplitudes = np.array([1.,.8,.3,2.])
    direct = viscosity_from_second_jets(ev,amplitudes)
    from_parts = float(np.dot(amplitudes**2,coefficients))
    # This assertion checks quadratic additivity on the same grid, not accuracy.
    assert np.isclose(direct,from_parts,rtol=1e-9,atol=1e-8)
    return dict(metadata=ev.metadata,rows=rows,
                mixed_amplitudes=amplitudes.tolist(),mixed_direct=direct,
                mixed_from_second_jet_parts=from_parts,
                mixed_stress_reference=float(ev.dP*amplitudes[2]**2+ev.dv*amplitudes[3]**2),
                status='Exploratory same-grid cross-check; no quadrature error bound or certified sign.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--symbolic-only',action='store_true')
    parser.add_argument('--nr',type=int,default=400)
    parser.add_argument('--nmu',type=int,default=512)
    args = parser.parse_args()
    symbolic_checks()
    if not args.symbolic_only:
        from evaluate_pressure_gate import Evaluator
        # A broader moderate geometry limits angular undersampling in this audit.
        # Lmax=2 is sufficient: this checker does not reconstruct pressure.
        ev = Evaluator(nr=args.nr,nmu=args.nmu,lmax=2,
                       sigma=1.,height=.7,distance=2.,swirl_ratio=.4)
        print(json.dumps(numerical_comparison(ev)),flush=True)


if __name__ == '__main__':
    main()
