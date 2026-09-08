# Strict local feedback without exact neutral tuning

The neutral construction is useful for identifying the second derivative, but favorable local feedback does not require exact tuning to that boundary. The same component bounds also certify a datum whose every velocity amplitude is rational. This is a robustness statement about initial data and short-time solutions, not a regeneration theorem.

Keep the radial profile, core strain and rotation, annular strength c=7/5 and viscosity nu=1/1000 from the main certificate. Let the exterior swirl amplitude A vary, and put H=A² C_v and e=H−204/35. Initially b=Omega=1. The exact central equations give

    b'=2+e/2,  Omega'=2,  beta'=e/2,
    S:=b'−2b²=e/2,
    Omega''=8+e.

The flat-core viscous curvature terms and their first derivatives still vanish. Differentiating the actual central equations, without evolving the initial radial ansatz, gives

    beta''=−T(H)/2−5e,
    S'=−T(H)/2−4e,
    T(H):=p_zz'(0)+32.

The full pressure reconstruction is linear in H once the meridional field is fixed. The same certified coefficients therefore imply for every H>0

    T(H)<102−(23199/1000)H.

Writing delta=290649/17500, it follows for e>=0 that

    beta''>delta+(13199/2000)e,
    S'>delta+(15199/2000)e.

The initial logarithmic growth bound for the specified outer swirl stress divided by Omega² is independent of A and remains greater than four. Meridional acceleration depends on A, but the initial azimuthal equation has no pressure component and uses the same initial meridional velocity. This independence is asserted only at the initial instant.

## A direct rational-amplitude example

Choose A=1100. The already computed outward interval endpoints give

    113/20000000 < C_v < 279/40000000,
    13673/2000 < H < 33759/4000.

In particular H>204/35. Hence this fully explicit compact datum has all of

    beta'(0)>0,  beta''(0)>0,
    b'(0)−2b(0)²>0,  (b'−2b²)'(0)>0,
    (log(C/Omega²))'(0)>4.

The associated solution retains the strict inequalities on some positive interval by smooth local existence and continuity. The master theorem keeps its exactly neutral amplitude because that cleanly separates the second-order pressure test. The rational-amplitude example is a separate initial datum.

The inequalities are also stable under sufficiently small smooth, axisymmetric, odd, divergence-free perturbations measured in H10, with an exterior stress weight fixed as in the finite-time corollary. Point evaluation of the required jets and the initial pressure/time-derivative functionals are continuous in that topology. This gives an unspecified open neighborhood within the symmetry class. It does not give a numerical perturbation radius, an invariant neighborhood under evolution, or a uniform sequence of inherited states.

`certificate/verify_robust_initial_data.py` checks the exact jet algebra, stored interval endpoint comparisons and resulting rational margins. Its analytic dependencies are the main reconstruction and flat-core identities.
