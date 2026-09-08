# Full initial field and pressure source in cylindrical Fourier modes

`initial_full_field.py` evaluates the analytic datum

```
u0 = S_Phi + W_psi + A v + lambda w_L,
Phi(s)=chi(s)+(7/5) eta(s),       psi(s)=chi(s),       s=r^2+z^2,
A=1020,                         lambda=1/4.
```

The displayed defaults are the requested fixed amplitudes. The evaluator
does **not** substitute the Pass8 pressure-neutral square-root tuning. All
cutoff transitions are retained. Outputs are floating-point evaluations of
analytic formulas, not interval enclosures or a pressure/evolution solve.

## Profiles on their entire support

The unchanged cutoff is one for `x<=1/4`, zero for `x>=1`, and on the
transition, with `t=(4x-1)/3`,

```
chi(x) = 1/[1+exp(-1/t+1/(1-t))].
```

The radial annular profile is

```
eta(s) = chi(4s/25)-chi((3s-64)/44)
       = -[1-chi(4s/25)] chi((3s-64)/44).
```

The equality is exact: the first cutoff has ended before the second leaves
its plateau. Thus `eta=0` for `R<=5/4` and `R>=6`, and `eta=-1` on
`5/2<=R<=5`, with the previously fixed inner and outer joins. Differentiating
the difference is particularly simple and retains every join derivative.

Write the axisymmetric components as `(a,W,d)`. Directly from the audited
radial vector-potential construction,

```
a = -r[Phi+2z^2 Phi'],
d =  2z[Phi+r^2 Phi'],
W = r[psi+(2s/3)psi'] + A V,
V = r chi(r^2/(63/200)^2) q_v(z),
q_v(z) = chi(((z-4)/(9/20))^2)+chi(((z+4)/(9/20))^2).
```

The core swirl is `r[psi+(2s/3)psi']`; replacing it by `r psi` would change
the field and its pressure. The module keeps the vector-potential version.
Every `Phi` and `psi` derivative here is with respect to `s`.

For the leading seed put

```
m=4,    k=-20,
e(r)=chi(((r-7/50)/(1/10))^2),
q_s(z)=chi(((z-4)/(3/5))^2)+chi(((z+4)/(3/5))^2),
b(r)=e(r) exp(ikr).
```

The genuine positive Fourier coefficient is

```
u_4 = (lambda/2) q_s (im b/r, -b', 0).                         (1)
```

The full real field is `u_0+u_4 exp(4i theta)+conj(u_4) exp(-4i theta)`.
Both derivatives of the nonconstant radial envelope and derivatives of its
axial cutoff are present in the returned first and second jets.

## Covariant gradient, source, and nonlinear acceleration

For any coefficient `(a_m,b_m,d_m)`, the cylindrical orthonormal gradient
has **output component as row, differentiation direction as column**:

```
G_m = [[a_m,r,  (im a_m-b_m)/r, a_m,z],
       [b_m,r,  (im b_m+a_m)/r, b_m,z],
       [d_m,r,          im d_m/r, d_m,z]].                       (2)
```

The angular-basis terms in the second column are essential. With
`G_-4=conj(G_4)`, pressure uses `tr(G^2)`, rather than a Frobenius norm:

```
g_0 = tr(G_0^2)+2 tr(G_4 conj(G_4)),
g_4 = 2 tr(G_0 G_4),
g_8 = tr(G_4^2),                  -Delta p_m = g_m.             (3)
```

These are the complete nonnegative Fourier modes of the initial source.
The analogous complete nonlinear acceleration coefficients are

```
B_0 = G_0 u_0+G_4 conj(u_4)+conj(G_4) u_4,
B_4 = G_0 u_4+G_4 u_0,
B_8 = G_4 u_4,                    B=(u.grad)u.                 (4)
```

The evaluator forms (2)--(4) before applying any support cancellation, so
every mean/core/outer/seed interaction is represented. The mode-zero
axisymmetric contribution can independently be written

```
tr(G_0^2) = a_r^2+(a/r)^2+d_z^2+2 a_z d_r-2 W W_r/r.          (5)
```

For this particular datum, the seed lives wholly within the annular pump's
affine region and is disjoint from the core. Hence its mixed source with
the meridional field and core swirl vanishes; this is a checked consequence,
not an implementation assumption. Globally the remaining mixed coefficient
also equals

```
g_4 = (lambda A/r)[(V0 b')'-(m^2/r)V0' b] q_s q_v,
V0(r)=r chi(r^2/(63/200)^2).                                  (6)
```

The `lambda/2` in (1) and the factor two in the mixed source have both been
included in (6). It is valid through the seed's radial and axial joins.
Absence of explicit `q_s'` in (6) does not license deleting axial velocity
derivatives from the gradient or nonlinear acceleration.

For another useful direct check, with
`B1=b'/r-b/r^2` and `B2=b'/r-m^2 b/r^2`, the seed-only sources are

```
(g_0)_seed = lambda^2 q_s^2 [m^2 |B1|^2-Re(B2 conj(b''))],
 g_8        = -(lambda^2/2) q_s^2 [m^2 B1^2+B2 b''].
```

No pressure solve or harmonic truncation is hidden in these expressions.
Although the initial source stops at mode eight, a subsequent nonlinear
time step can generate higher modes.

## Divergence, vector viscosity, and axis limits

The meridional formula obeys `a_r+a/r+d_z=0` identically; (1) has
`(u_4,r)_r+u_4,r/r+im u_4,theta/r=0` identically. These two identities are
checked symbolically, with arbitrary differentiable profiles.

The vector Laplacian is reconstructed from the six scalar jet entries via

```
L_m = partial_rr+r^-1 partial_r+partial_zz-m^2/r^2,
(Delta u)_r = L_m u_r-u_r/r^2-2im u_theta/r^2,
(Delta u)_theta = L_m u_theta-u_theta/r^2+2im u_r/r^2,
(Delta u)_z = L_m u_z.                                       (7)
```

At `r=0`, the code takes the smooth limits explicitly. The seed vanishes
on a fixed axis neighborhood; for the mean, `a/r -> a_r` and
`W/r -> W_r`. Its two horizontal vector-Laplacian components vanish on the
axis, and the axial one is `2d_rr+d_zz`. The affine core yields

```
G(0)=[[-1,-1,0],[1,-1,0],[0,0,2]],  g(0)=4,  Delta u(0)=0.
```

Given a separately solved full pressure, the initial actual NS derivative
is `u_t=nu Delta u-B-grad p`. The returned convection is positive `B`, not
its negative.

## Python interface and checks

```
from initial_full_field import fields, reconstruct
data = fields(r, z, A=1020., lam=.25, c=1.4)
u_cyl = reconstruct(data['u'], theta)
g = reconstruct(data['source'], theta)
```

`r,z` are NumPy-broadcast inputs. Components are `(r,theta,z)`;
`data['jets'][m]` has shape `(3,6,*broadcast_shape)` and derivative order
`(value,r,z,rr,rz,zz)`. Other fields and their mode sets are documented in the
module header. `cartesian_velocity` is provided for independent checks.

`check_initial_full_field.py` passed exact symbolic divergence and mixed
source identities. It also checked Cartesian finite differences at the
core, both annular joins, the seed's radial transitions, both axial joins,
the axis and exterior. Halving the step gave second-order convergence:
relative errors at `h=1e-5` were `1.84e-7` for the full Cartesian gradient,
`1.44e-6` for the vector Laplacian, and `2.31e-7` for
`div[(u.grad)u]` against (3). Modal reconstructions agree with direct matrix
contractions to about `1e-16`; evaluated divergence was below `2.9e-14`.
These are numerical consistency checks, not rigorous floating-point or PDE
error bounds. The recorded results are in `initial-full-field-check.json`.
