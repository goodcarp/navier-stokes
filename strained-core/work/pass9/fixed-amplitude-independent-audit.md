# Independent audit of the exact fixed-amplitude variant

Accepted: setting the outer-swirl coefficient to the exact rational
\(A=1020\), with the unchanged leading seed and
\(1/8\le\lambda\le1/4\), defines a computationally explicit initial
family that no longer needs an integral-defined amplitude.
It preserves the previously proved positive initial receiver and
fluctuation-energy gains and gives positive off-neutral core-ratio
velocity and acceleration. This is a separate explicit variant,
not a regenerated endpoint or a validated positive-time stage.

Write \(L=C_{\rm lead}>0\), \(X=A^2C_v\), and \(H_0=204/35\).
The seed, core and pump geometry is unchanged. The actual compact
field is
\[
 u_{0,\lambda}=M+1020v+\lambda w_{\rm lead}.
\]
The central affine neighborhood remains exactly
\(\operatorname{diag}(-1,-1,2)+J\), and odd/fourfold symmetry
remains valid. Its full initial pressure is now
\[
 h=p_{zz}(0)=-8-(X+\lambda^2L-H_0).
\]
Thus define
\[
 \delta=\frac{X+\lambda^2L-H_0}{2}.
\]
The certified lower bound on \(C_v\) implies
\[
 X>\frac{293913}{50000},\qquad
 \boxed{\delta>\frac{17391}{700000}>0.024844.}
\]
No exact pressure neutrality is claimed.

## Off-neutral core identities

The initial affine core still gives
\(d_b=d_\Omega=d_{b1}=d_{\Omega1}=0\).
Indeed \(\Delta u_0=0\) in its neighborhood, the local pressure
Laplacian is constant, and
\(\Delta N(u_0)=\nu\Delta^2u_0-\Delta(u_0\cdot\nabla u_0)
 +\nabla\operatorname{tr}(\nabla u_0\nabla u_0)=0\) there.
The nonlocal harmonic part of the pressure need not be zero for
these Laplacian identities to hold.

Consequently,
\[
 b_1=2+\delta,\quad \Omega_1=2,\quad
 b_2=-4(2+\delta)-F/2,\quad \Omega_2=8+2\delta,
 \quad F=p_{zz}'(0).
\]
The full quotient formula gives exactly
\[
 \boxed{\beta'(0)=\delta,\qquad
 \beta''(0)=-\frac{F+32}{2}-10\delta.}
\]
This explicitly pays the off-neutral correction rather than applying
the old neutral shortcut.

## Full pressure bound at fixed amplitude

The general pressure estimate was derived for arbitrary positive
outer amplitude, before the later neutrality substitution:
\[
 T=F+32<102-\frac{23199}{1000}X
              +\lambda^2(D_L+A E_L).
\]
The leading-envelope bounds remain
\[
 D_L<7-\frac{96}{5}L,\qquad |E_L|<\frac13 .
\]
They depend on the unchanged unit profiles, support geometry,
whole-space pressure projection and viscosity, not on choosing
an implicit neutral \(A_\lambda\).

Combining the estimates with the exact off-neutral quotient yields
\[
 \boxed{\beta''(0)>
 -\frac{153}{7}+\frac{13199}{2000}X
 -\frac{\lambda^2}{2}(7+A/3)+\frac{23}{5}\lambda^2L.}
\]
Here the coefficient of \(X\) and the coefficient of \(L\) are both
positive. Substituting \(A=1020\), the certified lower bound for
\(X\), and \(\lambda^2\le1/16\), then dropping only the positive
\(L\) term, gives
\[
 \boxed{\beta''(0)>
 \frac{4264878809}{700000000}>6.09268.}
\]
The same ingredients separately give
\[
 \boxed{T<-\frac{634112687}{50000000}<-12.68225.}
\]
The algebra was independently checked symbolically, including every
coefficient in the off-neutral correction.

## Which earlier bounds extend

The mean derivative estimate uses the strict unit torque bound
\(>1950\), \(g''(r_*)>-19\), and \(A\le1020\). It therefore remains
strict at this exact endpoint:
\[
 G_t(0,r_*,\pm4)>\frac{8871}{800}>11.
\]
The energy estimate uses \(A\ge900\), the unchanged \(I_r,I_z,R_0,R_1\)
and viscosity bounds. Hence
\[
 K'(0)/K(0)>\frac{12853}{2205}>5.82
\]
also remains valid. The positive outer initial mean maximum still
dominates the core maximum.

The new general-envelope initial-jet certificate explicitly relaxed
the amplitude and seed coefficient to the independent closed intervals
\(A\in[900,1020]\), \(\lambda^2\in[1/64,1/16]\).
Its actual full-pressure fixed-circle and fixed-height radial-branch
second-derivative bounds therefore include \(A=1020\).
Those outer jet identities did not use central pressure neutrality.

This choice removes the initial *amplitude quadrature* error from
the full-field validation problem. Numerical representation,
spatial and temporal residuals, and any remaining profile errors
still require bounds. No amplitude is retuned during evolution.
Positive initial core-ratio derivatives do not establish a
retained endpoint cone or a duration by themselves.
