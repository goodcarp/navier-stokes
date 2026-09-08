# Independent whole-space initial references

For the requested full compact datum
`u0=M+1020 v+(1/4)w_L`, `c=7/5`, `nu=1/1000`, the independent quadrature
references are

| Quantity | Numerical reference |
|---|---:|
| Initial core `p_zz(0)` | `-8.709035202638` |
| Initial `beta'`, `beta=b/Omega` | `0.354517601319` |
| Full kinetic energy `E=||u0||_2^2/2` | `81037.42658035` |
| Full kinetic-energy derivative `E'` | `-4611.62709024` |
| Actual nonaxisymmetric energy `K` | `6.443771766575` |
| Actual nonaxisymmetric derivative `K'` | `447.4261166273` |
| `K'/K` | `69.4354382550` |

These are high-accuracy numerical references, **not interval enclosures**.
They integrate the compact whole-space initial field, including both axial
endcaps. No finite-cylinder or periodic-axial pressure solver is used. They
are suitable checks of initialization and normalization, not certified
positive-time evolution.

The script is `initial_whole_space_integrals.py`; all runs are stored in
`initial-whole-space-integrals.json`. It uses the independently checked
`initial_full_field.py` and composite Gauss quadrature split at every cutoff
plateau/support join.

## Two independent initial pressure calculations

Only the angular-zero pressure mode contributes `p_zz` at the origin.
With `N=1/(4pi |x|)` and `g=tr((grad u)^2)`, the source formula is

```
p_zz(0)=-g(0)/3 + PV integral (3z^2-|x|^2)/(4pi |x|^5) g(x) dx.
```

Here `g(0)=4`, so the contact term is `-4/3`; omitting it gives the wrong
reference. The radial field `M=S_Phi+W_psi` is integrated in spherical
coordinates over `0<=R<=6`:

```
(p_zz)_M=-4/3 +(1/2) integral_0^6 dR/R
                  integral_-1^1 (3mu^2-1) g_M(R,mu) dmu.
```

In the affine ball `g_M=4`, whose angular principal-value contribution is
exactly zero. The code subtracts this constant there before numerical
integration. It retains all core and annular transition contributions.
At order 72 it obtains `-2.171428571428968`, agreeing with the independent
exact radial identity `-76/35=-2.171428571428571…`.

The remaining source is integrated in cylindrical coordinates over both
compact packets. By evenness in `z`, its pressure contribution is

```
integral_(z>0) integral_(r>0)
  r(2z^2-r^2)/(r^2+z^2)^(5/2) [g_full,0-g_M,0] dr dz.
```

Separately, the nonsingular stress-kernel calculation uses

```
Q_rr=[12(r^2+z^2)^2-105r^2z^2]/[4pi(r^2+z^2)^(9/2)],
Q_tt=3(4z^2-r^2)/[4pi(r^2+z^2)^(7/2)],
Cv=integral Q_tt V^2 dx,
CL=integral [Q_rr average_theta(w_L,r^2)
                  +Q_tt average_theta(w_L,theta^2)] dx,
p_zz=-76/35-A^2 Cv-lambda^2 CL.
```

The unit coefficients are numerically

```
Cv = 6.27156152642478e-6,
CL = 0.202784305876279.
```

At order 72 the source and stress results differ by `1.36e-13`.
The exact pressure identity and the two differently conditioned integrations
provide a useful cross-check; agreement does not itself certify that many
digits. The already certified fixed-amplitude initialization independently
gives `beta'>17391/700000`, equivalently
`p_zz<-8-2(17391/700000)≈-8.04968857`. The reference is consistent with it.

## Actual three-dimensional fluctuation energy and its derivative

Let `U_4` and `G_4` be genuine positive Fourier coefficients, including the
`lambda/2` factor. The tensor quadrature uses

```
K = 2pi integral_all_z integral r |U_4|^2 dr dz,
D_seed = 4pi integral_all_z integral r |G_4|_F^2 dr dz,
production = -4pi integral_all_z integral r Re(conj(U_4).G_0 U_4) dr dz,
K' = production - nu D_seed.
```

Its results are

```
production = 504.989536310509,
D_seed = 57563.4196832157,
viscous loss = 57.5634196832157.
```

For an independent calculation, define the complete radial/axial moments

```
R0=integral r[e'^2+(k^2+m^2/r^2)e^2] dr,
R1=integral r[(e''+e'/r-(k^2+m^2/r^2)e)^2
                       +k^2(2e'+e/r)^2] dr,
Ir=-integral r C'(r)e^2 dr,
Z0=integral q_s^2 dz,   Z1=integral (q_s')^2 dz,
Iz=integral q_s^2 q_v dz.
```

With `m=4`, `k=-20`, the exact separated identities are

```
K = (pi lambda^2/2) R0 Z0,
D_seed = pi lambda^2(R1 Z0+R0 Z1),
K' = pi lambda^2[80 A Ir Iz-c R0 Z0-nu(R1 Z0+R0 Z1)].
```

They agree with the tensor computation to rounding at order 72. In
particular, the full `q_s'` viscous contribution is retained. This is the
actual compact three-dimensional initial derivative, not the unweighted
planar-slice energy derivative or a later-time model prediction.

## Total energy and numerical convergence

The radial meridional field is pointwise orthogonal to the azimuthal fields,
and the core swirl and outer swirl have disjoint supports. Their gradient
Frobenius cross terms vanish as well. Angular averaging removes the mean/
seed cross terms. Thus total energy and dissipation decompose exactly as

```
E = E_M+A^2 E_v+K,
D = D_M+A^2 D_v+D_seed,        E'=-nu D.
```

The last equality is the actual unforced whole-space energy identity. The
quadrature values are `E_M≈77623.4542737116`,
`E_v≈0.00327521004889325`, `D_M≈2879977.45284969`, and
`D_v≈1.60907940956397`. The large annular mean field dominates the total
energy; it must be included when comparing with a full-field pilot.

| Gauss order per piece | `p_zz` from source | `K'` from tensor |
|---:|---:|---:|
| 32 | `-8.7090368908755` | `447.4260827962` |
| 48 | `-8.7090352036234` | `447.4261166236` |
| 72 | `-8.7090352026380` | `447.4261166273` |

The stress pressure was already stable at the displayed precision at order
32. The changes from order 48 to 72 were about `2.7e-9` in total energy and
`5.4e-7` in total energy derivative. These convergence comparisons are
observations, not rigorous quadrature-error estimates. A periodic or
finite-radius pressure pilot can additionally differ because of its actual
boundary problem, even if its local derivative implementation is correct.

This bounded reference calculation does not compute `p_zz'`; that quantity
requires the full initial pressure response in the acceleration, which the
separate pressure-solver work is addressing.
