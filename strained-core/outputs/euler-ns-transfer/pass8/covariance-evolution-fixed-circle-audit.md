# Independent audit of the actual fixed-circle acceleration certificate

Accepted: certify_fixed_circle_acceleration.py uses the actual full-pressure
second-jet formula derived in covariance-evolution.md. Its saved outward
certificate proves
\[
 \boxed{G_{tt}(0,r_*,\pm4)<-250000}
\]
uniformly for the actual Pass7 family \(3/4\le\lambda\le1\).
The recorded upper expression is below \(-293333\).
This is an initial derivative of the actual smooth NS solution, not a
finite-time evolution estimate.

## 1. Cutoff jets and the maximizing radius

The derivatives of the logistic composition through fourth order are
correct. In particular the fourth derivative is
\[
 p(1-p)\left[
 (1-14p+36p^2-24p^3)z_1^4
 +6(1-6p+6p^2)z_1^2z_2
 +3(1-2p)z_2^2+4(1-2p)z_1z_3+z_4\right].
\]
The code includes every \(4/3\) chain-rule factor in the derivatives
of the transition exponent. It uses this formula only in a compact
interior transition interval, so no endpoint derivative singularity
is evaluated.

For \(g(r)=r^2\chi(r^2/a^2)\), the displayed coefficients
\((24,156,112,16)\) in \(g''''\) are also correct.
The standalone audit checker independently generates these derivatives.

The earlier Pass6 argument puts every positive global radial maximum
in \((a/2,4a/5)\). The new range sums prove \(g'>0\) below .2,
\(g'<0\) above .23, and \(g''<0\) on \([.2,.23]\).
Endpoint signs in
\[
 0.2154755753334<r_*<0.2154755753336
\]
then isolate its unique root and prove the uniqueness of the global
radial maximum. Multiplication by \(A_\lambda>0\) does not change
this radius. The flat seed region contains that entire interval.

## 2. Radial pressure and the second-jet upper bound

The Green moment limits are correct: \(C_0'=0\) below \(a/2\)
and above \(a\). The seed is flat at the maximizing radius.
The small partial-cell terms enclose the unknown exact endpoint
inside its isolating interval. The exact derivative is
\[
 P'=-40iC_0e^{20ir}+12r^{-5}I_+-20r^3I_-.
\]
Its real/imaginary implementation and the phase-rotated multiplication
by \(K=(k^2r+m^2/r)+ik\) have the correct signs.

At a true maximum \(V_0'=-C_0\). Thus the full source is
\[
 e^{-ikr}s_m=2C_0(m^2/r^2-k^2-ik/r),\qquad
 \tfrac12\operatorname{Re}[-ikr e^{-ikr}s_m]=-k^2C_0.
\]
Together with the independently proved whole-space
\(|p_{m,r}-P'|<5\), the actual seed coefficient is at most
\[
 S_+=\frac{m^2g''}{2r^2}
 +\frac12\operatorname{Re}(e^{-ikr}KP')
 -k^2C_0+\frac52|K|.
\]
This matches the code. It is an upper bound; its interval lower
endpoint does not bound the actual coefficient from below.

The exact formula has the form
\[
 G_{tt}=A_\lambda L_0^2g+\lambda^2 A_\lambda S
 +\nu\mathcal T_0[(2m^2+4)/r^2-2k^2],\qquad S\le S_+.
\]
The saved ranges prove \(L_0^2g<0\), \(S_+<0\), and the last bracket
negative. Since \(A_\lambda>780\) and \(\lambda^2\ge9/16\), all
inequality directions in
\[
 G_{tt}<780\left[L_0^2g+\frac9{16}S_+\right]
\]
are valid. Discarding the negative viscous torque term is conservative.
The fourth-jet expression used for \(L_0^2g\) is exactly
\[
 c^2r^2g''-2c\nu r g'''+
 \nu^2(g''''-2g'''/r+3g''/r^2)
\]
at the true critical radius; terms proportional to \(g'\) vanish
exactly rather than being numerically presumed zero.

## 3. Radial recentering is a distinct, computable correction

For each fixed \(z=\pm4\), \(G_{rr}(0,r_*,z)=A_\lambda g''(r_*)<0\).
The implicit-function theorem supplies a smooth local radial
maximizing branch \(r(t,z)\). Its initial speed is
\[
 \dot r(0,z)=-\frac{G_{tr}}{A_\lambda g''}
 =cr-\nu(g'''/g''-1/r)
       +\frac{\mathcal T_0}{rA_\lambda g''},
\]
where
\[
 G_{tr}=A_\lambda[-crg''+\nu(g'''-g''/r)]-\mathcal T_0/r.
\]
For the recentered value \(H(t,z)=G(t,r(t,z),z)\),
\[
 \boxed{H_{tt}(0,z)=G_{tt}(0,r_*,z)
                  -\frac{G_{tr}^2}{A_\lambda g''}.}
\]
The correction is positive. It may not be omitted when discussing
the radially tracked value. The extended certificate correctly
adds its outward upper bound, and proves
\[
 H_{tt}(0,\pm4)<-240000,\qquad
 85<H_t(0,\pm4)<186,\qquad
 \dot r(0,\pm4)>1/5.
\]
Treating \(A_\lambda\) and \(\lambda^2\) as independent intervals
enlarges the enclosure and is valid despite their exact retuning
relation.

## 4. Why this does not determine the global maximum's acceleration

Initially the exterior mean is constant in \(z\) across each axial
plateau \(|z\mp4|\le9/40\). Its global maximizing set therefore
contains a whole interval of maximizing circles in each packet.
The axial Hessian vanishes there; a nondegenerate two-dimensional
maximum theorem does not apply.

The first mean derivative is the same on the interior of each
plateau, but the full pressure in the second derivative may vary
with \(z\). The present \(|p_{m,r}-P'|<5\) estimate is centered on
the two planes \(z=\pm4\); it does not give the same uniform error
up to the axial plateau edges. A global maximizing value can select
a different axial location. Flat cutoff edges also require care in
any second-order envelope expansion; existence of a classical
second derivative for the global supremum is not asserted.

Even for the tracked radial branch at \(z=4\), negative initial
acceleration and positive initial velocity do not prove a cessation
time near \(G_t/|G_{tt}|\). That ratio is only a local linearized
scale without control of subsequent derivatives. Nothing here
rules out later recovery of torque or increased global mean maximum.

The audit checker independently verifies the logistic and profile
jets, pressure-source simplification, phase multiplication, moving
radial-branch identities, and exact saved dyadic inequalities.
No additional profile quadrature was required for this audit.
