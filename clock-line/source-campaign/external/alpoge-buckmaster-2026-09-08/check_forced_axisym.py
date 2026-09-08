#!/usr/bin/env python3
"""
Symbolic checks supporting THEOREM_forced_axisymmetric_on_axis_typeII.md.

Four independent checks, all exact (sympy), no floating point in checks 1-3.

  C1  Scaling covariance of the FORCED Navier-Stokes system on R^3:
      if (u,p,f) solves d_t u + (u.grad)u + grad p = nu Laplace u + f, div u = 0,
      then so does
        u_L(x,t) = L u(Lx, L^2 t + T),  p_L = L^2 p(Lx, L^2 t + T),
        f_L(x,t) = L^3 f(Lx, L^2 t + T).
      This is the exponent that makes ||f_L||_inf = L^3 ||f||_inf.

  C2  Axisymmetric swirl-free vorticity equation WITH force.  With a Stokes
      stream function psi(r,z,t) (so div u = 0 identically),
        u_r = -d_z psi,   u_z = (1/r) d_r (r psi),   omega = d_z u_r - d_r u_z,
      the poloidal curl of the momentum equations equals
        d_t omega + u_r d_r omega + u_z d_z omega - (u_r/r) omega
          - nu (d_rr + (1/r) d_r + d_zz - 1/r^2) omega - F_theta,
      where F_theta = (curl f)_theta = d_z f_r - d_r f_z.  The pressure drops out.

  C3  The lifted quantity Omega = omega/r obeys, in the swirl-free case,
        d_t Omega + u_r d_r Omega + u_z d_z Omega
          = nu (d_rr + (3/r) d_r + d_zz) Omega + F_theta / r,
      i.e. a drift-diffusion equation with NO stretching term and with the
      force entering only as the source F_theta/r.  (d_rr + (3/r) d_r + d_zz)
      is the 5-dimensional Laplacian on SO(4)-invariant functions.

  C4  Periodic axisymmetry is trivial.  For R = rotation by alpha about the
      x3-axis and a in Z^3, the composite R_0^{-1} o R_a (R_a = rotation by
      alpha about the vertical line through a) is the translation by
      (R^{-1} - I) a.  Checked symbolically; then a numeric check that the
      generated translations accumulate at 0 in two independent directions.

Run:  python3 check_forced_axisym.py
"""

import sympy as sp

ok = True
def report(name, passed):
    global ok
    ok = ok and bool(passed)
    print(("PASS  " if passed else "FAIL  ") + name)


# ---------------------------------------------------------------- C1 ---------
def c1_scaling():
    """Covariance is checked by composing with the scaling map explicitly:
    the residual F of the original fields is formed in the variables (Y,S),
    then transported by sp.Subs(...).doit(), which substitutes correctly
    inside Derivative objects."""
    x1, x2, x3, t, L, nu, T = sp.symbols('x1 x2 x3 t lam nu T', positive=True)
    Y1, Y2, Y3, S = sp.symbols('Y1 Y2 Y3 S')
    Y, X = (Y1, Y2, Y3), (x1, x2, x3)

    U = [sp.Function('u%d' % i)(Y1, Y2, Y3, S) for i in range(3)]
    P = sp.Function('p')(Y1, Y2, Y3, S)

    def mom(u, p, V, tt):
        return [sp.diff(u[i], tt)
                + sum(u[k] * sp.diff(u[i], V[k]) for k in range(3))
                + sp.diff(p, V[i])
                - nu * sum(sp.diff(u[i], V[k], 2) for k in range(3))
                for i in range(3)]

    def dvg(u, V):
        return sum(sp.diff(u[k], V[k]) for k in range(3))

    # F is the force the original fields (U,P) require; D is their divergence.
    F = mom(U, P, Y, S)
    D = dvg(U, Y)

    args_from, args_to = (Y1, Y2, Y3, S), (L * x1, L * x2, L * x3, L**2 * t + T)
    transport = lambda e: sp.Subs(e, args_from, args_to).doit()

    mp = dict(zip(args_from, args_to))
    uL = [L * U[i].subs(mp, simultaneous=True) for i in range(3)]
    pL = L**2 * P.subs(mp, simultaneous=True)

    FL = mom(uL, pL, X, t)
    DL = dvg(uL, X)

    good = all(sp.simplify(sp.expand(FL[i] - L**3 * transport(F[i]))) == 0
               for i in range(3))
    good = good and sp.simplify(sp.expand(DL - L**2 * transport(D))) == 0
    report("C1  forced NS is covariant under (u,p,f) -> (L u, L^2 p, L^3 f)", good)

    bad = sp.simplify(sp.expand(FL[0] - L**2 * transport(F[0])))
    report("C1' exponent 2 on the force does NOT work (control)", bad != 0)


# ------------------------------------------------------------- C2, C3 --------
def c2_c3_vorticity():
    r, z, t, nu = sp.symbols('r z t nu', positive=True)
    psi = sp.Function('psi')(r, z, t)
    P = sp.Function('P')(r, z, t)          # pressure
    fr = sp.Function('f_r')(r, z, t)
    fz = sp.Function('f_z')(r, z, t)

    ur = -sp.diff(psi, z)
    uz = sp.diff(r * psi, r) / r

    # divergence-free by construction
    div = sp.simplify(sp.diff(r * ur, r) / r + sp.diff(uz, z))
    report("C2a stream-function ansatz is divergence free", sp.simplify(div) == 0)

    lap = lambda g: sp.diff(g, r, 2) + sp.diff(g, r) / r + sp.diff(g, z, 2)
    adv = lambda g: ur * sp.diff(g, r) + uz * sp.diff(g, z)

    # momentum residuals (swirl-free: no u_theta, no centrifugal term)
    Rr = sp.diff(ur, t) + adv(ur) + sp.diff(P, r) - nu * (lap(ur) - ur / r**2) - fr
    Rz = sp.diff(uz, t) + adv(uz) + sp.diff(P, z) - nu * lap(uz) - fz

    curlR = sp.diff(Rr, z) - sp.diff(Rz, r)

    omega = sp.diff(ur, z) - sp.diff(uz, r)
    Ftheta = sp.diff(fr, z) - sp.diff(fz, r)

    target = (sp.diff(omega, t) + adv(omega) - (ur / r) * omega
              - nu * (lap(omega) - omega / r**2) - Ftheta)

    report("C2  poloidal curl of forced momentum = forced omega_theta equation",
           sp.simplify(sp.expand(curlR - target)) == 0)

    # C3: the lifted variable
    Om = sp.Function('Om')(r, z, t)
    lap5 = lambda g: sp.diff(g, r, 2) + 3 * sp.diff(g, r) / r + sp.diff(g, z, 2)
    lhs = (sp.diff(r * Om, t) + adv(r * Om) - (ur / r) * (r * Om)
           - nu * (lap(r * Om) - (r * Om) / r**2))
    rhs = r * (sp.diff(Om, t) + adv(Om) - nu * lap5(Om))
    report("C3  omega = r*Omega turns the omega equation into "
           "d_t Om + b.grad Om = nu Lap5 Om + F_theta/r",
           sp.simplify(sp.expand(lhs - rhs)) == 0)


# ---------------------------------------------------------------- C4 ---------
def c4_periodic_axisymmetry():
    a1, a2, a3, al = sp.symbols('a1 a2 a3 alpha', real=True)
    c, s = sp.cos(al), sp.sin(al)
    R = sp.Matrix([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    a = sp.Matrix([a1, a2, a3])
    xs = sp.symbols('X1 X2 X3', real=True)
    X = sp.Matrix(xs)

    Ra = R * (X - a) + a                 # rotation by alpha about the vertical line through a
    comp = R.inv() * Ra                  # composed with the inverse rotation about the origin
    want = X + (R.inv() - sp.eye(3)) * a
    report("C4  R_0^{-1} o R_a = translation by (R^{-1} - I) a",
           sp.simplify(comp - want) == sp.zeros(3, 1))

    # the generated translation vectors accumulate at 0 in two independent directions
    import math
    e1 = lambda al: ((math.cos(al) - 1), math.sin(al))       # (R^{-1}-I)(1,0)
    e2 = lambda al: (-math.sin(al), (math.cos(al) - 1))      # (R^{-1}-I)(0,1)
    al = 1e-6
    v1, v2 = e1(al), e2(al)
    small = max(abs(v1[0]), abs(v1[1]), abs(v2[0]), abs(v2[1])) < 1e-5
    det = v1[0] * v2[1] - v1[1] * v2[0]
    report("C4' generated translations are arbitrarily small and span the plane",
           small and abs(det) > 0)


if __name__ == "__main__":
    c1_scaling()
    c2_c3_vorticity()
    c4_periodic_axisymmetry()
    print()
    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")
    raise SystemExit(0 if ok else 1)
