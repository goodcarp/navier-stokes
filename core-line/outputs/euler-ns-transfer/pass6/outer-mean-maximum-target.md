# An outer chiral seed that increases the actual mean angular-momentum maximum

There is an explicit compact nonaxisymmetric perturbation for which the **initial azimuthal-mean maximum** increases, while the central pressure is retuned exactly to the earlier neutral value. A conservative sufficient seed-amplitude interval is `[3/4,1]` in the vector-potential normalization below. Whether the full initial pressure-derivative feedback inequality survives at these amplitudes is not established here; that is the remaining compatibility test.

This is an initial derivative and short-time mean-maximum statement for ordinary unforced NS, with viscosity `nu=1/1000`. It does not supply persistent torque, a return, or a singularity.

## 1. Base datum and an outer-supported seed

Retain the pass 5 radial meridional datum `S_Phi`, core rotation `W_psi`, and unit outer swirl `v`. Its parameters are

\[
a=\frac{63}{200},\quad h=\frac9{20},\quad d=4,
\quad H=\frac{204}{35},\quad C_v=-\Pi(v,v)>0.
\]

The positive exterior mean angular momentum is

\[
\Gamma_{\rm out}(r,z)=A r^2\chi(r^2/a^2)
\big[\chi((z-d)^2/h^2)+\chi((z+d)^2/h^2)\big].
\]

Set `r0=21/100`, `delta=7/50`, and `h_s=3/5`. Define

\[
e(r)=\chi((r-r_0)^2/\delta^2),\qquad
q(z)=\chi((z-4)^2/h_s^2)+\chi((z+4)^2/h_s^2),
\]
\[
E(r,z)=e(r)q(z),\qquad
w=\nabla\times\big[E(r,z)\cos(4\theta+20r)e_z\big].
\tag{1}
\]

The seed is supported in `7/100<=r<=7/20`, `17/5<=|z|<=23/5`. Its support stays away from the axis and the core, lies in `5/2<R<5`, and satisfies `r/|z|<=7/68<1/4`. It is smooth, compact, solenoidal, and has zero azimuthal mean. Its vector potential is even under full spatial inversion and invariant under rotation by `pi/2`, so the velocity is odd and `C4` equivariant. No axisymmetry is asserted.

The radial envelope equals one on `[7/50,7/25]=[0.14,0.28]`; the axial envelope equals one near `z=+4` and `z=-4`. The single-mode computation in `nonaxisymmetric-torque-budget.md` gives

\[
\mathcal T_w=\frac{40}{r}\partial_r(rE^2),
\qquad \mathcal T_w=\frac{40}{r}\quad\hbox{where }E=1.
\tag{2}
\]

For the actual seed `lambda w`, the initial mean torque is `lambda^2 T_w`. Linear seed/base cross terms average to zero.

## 2. All relevant maximizing circles lie in the seed plateau

Write `xi=r^2/a^2` and `f(xi)=xi chi(xi)`. Every positive maximum of `f` occurs in

\[
\frac14<\xi_*<\frac58.
\tag{3}
\]

Below `1/4`, `f=xi` is strictly increasing, and `f'(1/4)=1`. To justify the other endpoint without a numerical root search, write the transition argument `t=(4xi-1)/3` and `p=chi(xi)`. For `xi>=5/8` and `xi<1`, `t>=1/2`, `1-p>=1/2`, and `t^(-2)+(1-t)^(-2)>=8`. Thus

\[
f'=p\left[1-\frac43\xi(1-p)
\{t^{-2}+(1-t)^{-2}\}\right]
\le p(1-\tfrac{16}{3}\xi)<0.
\]

This proves (3), without assuming uniqueness of the maximum. In particular,

\[
\frac a2<r_*<a\sqrt{5/8}<\frac45a<\frac7{25},
\qquad \frac a2>\frac7{50}.
\]

At every maximizing circle `(r_*,z=+4)` or `(r_*,z=-4)`, the seed envelope is flat and equals one. These circles are also axial maxima of the exterior base profile. Therefore

\[
\boxed{\mathcal T_w(r_*,\pm4)>
\frac{40}{(4/5)a}=\frac{10000}{63}.}
\tag{4}
\]

## 3. Exact initial viscous cost

Let `G=overline(Gamma)` and `Delta_*=partial_rr-r^(-1)partial_r+partial_zz`. At a maximizing circle, mean meridional advection vanishes because `G_r=G_z=0`. The axial cutoff is flat there. The mean equation consequently gives the exact insertion identity

\[
\boxed{\partial_tG(0,r_*,\pm4)
=\lambda^2\mathcal T_w(r_*,\pm4)
+4\nu A\xi_*[2\chi'(\xi_*)+\xi_*\chi''(\xi_*)].}
\tag{5}
\]

The second term is nonpositive because it is `nu` times the radial second derivative at a maximum. Direct interval ranges for the unchanged cutoff certify

\[
|\chi'|<4,\qquad |\chi''|<32.
\tag{6}
\]

Using `xi_*<5/8` gives the conservative lower bound

\[
4\nu A\xi_*[2\chi'+\xi_*\chi'']\ge-70\nu A.
\tag{7}
\]

Thus an exact necessary-and-sufficient sign test at the stated circle is `lambda^2 T_w > -nu Delta_* G_base`; a sufficient uniform inequality is

\[
\boxed{\lambda^2\frac{10000}{63}>70\nu A.}
\tag{8}
\]

The actual mean-zero seed has no direct mean viscous term, and the azimuthal pressure torque averages to zero. These observations do not discard its pressure effect on the separate central feedback functional.

## 4. Exact neutral retuning is feasible throughout a finite amplitude interval

Define the seed pressure coefficient by its actual whole-space value

\[
C_w=-\Pi(w,w)=\int Q_{ij}w_iw_j\,dx,
\qquad Q_{ij}=-\partial_{zzij}(4\pi|x|)^{-1}.
\]

The seed is horizontal and separated from the origin. On its cone support, `Q_rr,Q_thetatheta` are strictly positive, so `C_w>0`. An inexpensive analytic upper bound is also available:

\[
\boxed{C_w<\frac32.}
\tag{9}
\]

Here are the details, retaining the full tensor. Angular integration of the seed gives exactly

\[
C_w=\pi\int r\left[
Q_{rr}\frac{16E^2}{r^2}
+Q_{\theta\theta}(E_r^2+400E^2)\right]dr\,dz.
\]

On the support, both tensor entries are at most `3/(pi R^5)<=3/[pi(17/5)^5]`. Since the axial bumps are disjoint, `integral q^2 dz<=12/5`. Also

\[
16\int\frac{e^2}{r}dr\le16\log5<\frac{80}{3},\qquad
400\int re^2dr\le\frac{588}{25}.
\]

The radial bump has total variation two. Equation (6) gives `||e_r||_infinity<=400/7`, hence

\[
\int r e_r^2dr\le\frac7{20}\frac{400}{7}\,2=40.
\]

Together these imply

\[
C_w<3(\tfrac5{17})^5\frac{12}{5}
\left(\frac{80}{3}+\frac{588}{25}+40\right)<\frac32.
\]

All seed/base bilinear pressure terms at the origin vanish by the nonzero Fourier mode. Accordingly define the **named exact amplitude**, not a rounded approximation,

\[
\boxed{A_\lambda=\sqrt{\frac{H-\lambda^2C_w}{C_v}},\qquad
u_{0,\lambda}=S_\Phi+W_\psi+A_\lambda v+\lambda w.}
\tag{10}
\]

For `3/4<=lambda<=1`, (9) makes the radicand positive. The earlier outward certificate

\[
\frac{113}{20000000}<C_v<\frac{279}{40000000}
\]

implies the convenient bounds `780<A_lambda<1020`. Equation (10) preserves `p_zz(0)=-8` exactly, as its total exterior contribution is `-A_lambda^2 C_v-lambda^2 C_w=-H`. The seed vanishes near the origin, and its odd/C4 symmetries preserve the affine central matrix and the corresponding initial neutral first-jet identities.

The positive global mean maximum is on the outer packets: the core positive mean angular momentum is at most one, while the exterior value at `xi=1/4,z=4` exceeds `780 a^2/4>19`. The mean-zero seed does not change the initial mean beyond the explicitly retuned `A_lambda`.

## 5. A finite initial mean-maximum gain, and the remaining test

Combining (4), (7), `nu=1/1000`, and `A_lambda<1020`, every `lambda in[3/4,1]` satisfies

\[
\boxed{\partial_tG(0,r_*,\pm4)>
\frac9{16}\frac{10000}{63}-\frac{357}{5}
=\frac{626}{35}>0.}
\tag{11}
\]

At `lambda=1`, the stronger lower bound is `27509/315>87`. Smooth local existence therefore gives `sup G(t)>sup G(0)` for all sufficiently small positive times: evaluate `G(t)` at any one of these initially maximizing circles. No differentiability of the maximizer as a function of time is needed. This conclusion concerns the azimuthal-mean maximum; it does not assert the same derivative bound for the maximum of the fully angle-dependent `Gamma`.

The neutral central pressure in (10) does **not** establish `partial_t p_zz(0)+32<0` for this new family. The explicit small-perturbation continuity radius from the other pass 6 note is not known to contain the finite amplitudes `[3/4,1]`. The next test is the complete actual pressure-derivative functional for `u_{0,lambda}` at an amplitude satisfying (8), with its nonlocal acceleration retained. Passing both signs would still leave torque persistence and return on the inherited solution unproved.

`verify_outer_mean_maximum.py` recomputes the cutoff-jet range certificate using 64 rational panels, then checks the rational support, pressure-budget, retuning, and torque inequalities. Its exact dyadic cutoff bounds are saved in `outer-mean-maximum-certificate.json`. Its trust base is the displayed analysis and `mpmath.iv`; no full pressure-derivative quadrature or PDE simulation is performed.
