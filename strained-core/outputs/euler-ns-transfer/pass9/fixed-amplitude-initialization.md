# An explicit amplitude for the next full-field evolution

The next numerical or analytic evolution can start from the concrete field

\[
 u_0=M+1020v+\tfrac14 w_L,\qquad \nu=1/1000,
\]

with the unchanged compact profiles from pass8. Every scalar parameter is
now rational. The outer amplitude no longer depends on evaluating a pressure
integral to very high precision. This is a new choice of initial amplitude,
not an evolved or regenerated endpoint. The pressure-neutral pass8 family
and its certificates remain unchanged historical results.

Exact pressure neutrality is unnecessary for the retained-gain target. The
explicit amplitude gives an initially strict core cone, while preserving
the receiving and energy gains. The calculation below covers the same
field with any seed amplitude `1/8 <= lambda <= 1/4`; the selected next
evolution uses `lambda=1/4`.

## Full initial core identities without neutrality

Put `X=1020^2 C_v`, `L=C_L`, `q=lambda^2`, and

\[
 H_0=204/35,\qquad d=(X+qL-H_0)/2.
\]

The exact initial pressure is `p_zz=-76/35-X-qL`. Symmetry and the
unchanged affine initial core give

\[
 b=\Omega=1,\quad b'=2+d,\quad\Omega'=2,\quad\beta'=d,
 \qquad\beta=b/\Omega.
\]

The local NS acceleration is harmonic in this affine-core neighborhood:
the pressure Laplacian is spatially constant there. Thus both initial
viscous curvature functionals and their first time derivatives vanish.
Writing `T=p_zz'(0)+32`, the full off-neutral identity is

\[
 \boxed{\beta''(0)=-T/2-10d.}
\]

Using the old neutral formula `-T/2` here would be incorrect.

The already certified pressure normalization satisfies
`C_v>113/20000000`, so

\[
 X>293913/50000>H_0,\qquad
 \boxed{\beta'(0)>17391/700000>0.02484.}
\]

This is an actual first-order cone margin for the fixed-amplitude field.

## Complete pressure and the second core derivative

The all-amplitude base-pressure estimate and the exact angular-mode
selection give, before imposing any neutrality relation,

\[
 T<102-\frac{23199}{1000}X+q(D_L+1020E_L).
\]

The unchanged leading-seed certificates give
`D_L<7-(96/5)L`, `|E_L|<1/3`, and `L>0`. Hence

\[
 T<-\frac{634112687}{50000000}<-12.68.
\]

To estimate `beta''`, keep the shared `X,L` dependence instead of taking
unrelated upper and lower bounds on `T` and `d`:

\[
 \beta''(0)>
 -\frac{153}{7}+\frac{13199}{2000}X
 -\frac q2(7+1020/3)+\frac{23}{5}qL.
\]

Substituting `X>293913/50000`, `q<=1/16`, and discarding the positive
last term yields

\[
 \boxed{\beta''(0)>4264878809/700000000>6.09.}
\]

All pressure contributions, including viscosity and the full nonlocal
swirl/seed interaction, remain in these bounds.

## Receiving and supplying fields

The pass8 torque and energy bounds used only `900<A<=1020`, not pressure
neutrality. They therefore apply to this amplitude. Uniformly over the
stated lambda interval,

\[
 G_t(0,r_*,4)>8871/800,\qquad K'(0)/K(0)>12853/2205.
\]

For the selected `lambda=1/4`, the first bound improves directly to

\[
 \boxed{G_t(0,r_*,4)>1950/16-(19/1000)1020
                    =20499/200>102.49.}
\]

The pass9 whole-pressure initial mean-jet formula also needs no core
neutrality: at the initial mean maximum the mean-acceleration term is
multiplied by the zero initial mean gradient. Its independent amplitude
interval includes 1020. Accordingly its deceleration bounds still apply
to this explicitly initialized field; favorable first derivatives do not
give indefinite acceleration or an actual turning time.

## What this removes, and what remains

The H4 validation bridge can now compare a reconstructed approximation
with one initial field whose amplitude is exactly representable. Cutoff,
spatial reconstruction, domain and numerical errors still need bounds.
Pressure normalization is still needed for the mathematical estimates,
but uncertainty in an implicitly defined amplitude is no longer an
initial-data error source.

The finite planar pulse is a guide to selecting a test interval. It remains
necessary to evolve the compact full three-dimensional field, including
the axial pressure response and older field, enclose its full residual,
and carry the same actual endpoint into a usable successor class. This
amplitude choice supplies none of those missing evolution or return claims.
