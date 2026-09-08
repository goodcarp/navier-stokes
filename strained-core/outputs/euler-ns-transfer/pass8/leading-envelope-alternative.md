# A leading envelope that gains both mean angular momentum and fluctuation energy

2026-09-08. There is an explicit smooth compact replacement for the trailing outer seed that passes **two actual initial NS tests simultaneously**: the mean angular-momentum maximum increases, and the nonaxisymmetric fluctuation gains kinetic energy after paying both pump work and viscosity. Central pressure neutrality can be retuned exactly. The new full central pressure-derivative gate is a separate test; no later return or singularity is asserted here.

## 1. Why positive torque need not drain the seed

For the horizontal divergence-free chiral potential

`w=curl[e(r)q(z) cos(m theta+kr) e_z]`,

the exact initial covariance and azimuthal-mean torque are

`R=overline(w_r w_theta)=-mk e^2q^2/(2r)`,

`T=-(1/r)partial_r(r^2 R)=mk q^2/(2r) partial_r(r e^2)`.

The phase/envelope derivative cross term averages to zero in R; a nonconstant radial envelope is nevertheless retained in T. With `mk<0`, R is positive. When `e_r/e<-1/(2r)`, the radial stress flux decreases strongly enough that T is also positive. On the decreasing background angular velocity, positive R extracts energy from the mean. Thus positive mean torque and positive fluctuation-energy extraction are compatible; the older flat-envelope seed did not exploit this combination.

## 2. An explicit compact datum

Use the unchanged smooth cutoff chi from the existing certificates, and define

`e_new(r)=chi(((r-7/50)/(1/10))^2)`,

`q(z)=chi(((z-4)/(3/5))^2)+chi(((z+4)/(3/5))^2)`,

`w_new=curl[e_new(r)q(z) cos(4theta-20r) e_z]`.       (1)

It is smooth, compact, solenoidal, odd under inversion, and C4-equivariant. Its support satisfies

`1/25<=r<=6/25`, `17/5<=|z|<=23/5`.

It avoids the axis and the core, lies in the actual annular pump's affine plateau `5/2<|x|<5`, and satisfies `r/|z|<=6/85<1/4`. The seed has zero azimuthal mean. Its radial plateau is `[9/100,19/100]`; the maximizing circle is deliberately in its decreasing transition.

Let `M=S_Phi+W_psi`, retain the unit axisymmetric outer swirl v, and set

`C_new=-Pi(w_new,w_new)>0`,

`A_lambda=sqrt[(204/35-lambda^2 C_new)/C_v]`,

`u_(0,lambda)=M+A_lambda v+lambda w_new`,

`1/8<=lambda<=1/4`, `nu=1/1000`.       (2)

Equation (2) denotes one actual initial field. It does not mean that the outer swirl is retuned again during its later evolution.

## 3. Certified torque at the unchanged mean maximum

The positive outer mean remains `A_lambda r^2 C(r)q_v(z)`, because the new seed has zero mean. Its unique maximizing radius has the already certified enclosure

`1077377876667/5000000000000 < r_* < 269344469167/1250000000000`.

Outward-rounded cutoff jets on this interval give

`0.64645152592248 < e_new(r_*) < 0.64645152593035`,

`-39.29483882084 < e_new'(r_*) < -39.29483881870`,

and the actual unit torque satisfies

`T_new(r_*)=-40/r_* [e_new^2+2r_* e_new e_new'] > 1950`.       (3)

The stored sharper interval is approximately `[1954.599521999,1954.599522140]`. It is a range enclosure, not a sampled root evaluation. The axial seed cutoff is one on the entire initial axial maximum plateau `|z∓4|<=9/40`.

At every such initially maximizing circle, the mean-advection terms vanish, the axial mean second derivative vanishes, and the mean equation gives

`partial_t overline(Gamma)(0)=lambda^2 T_new+nu A_lambda [r^2C(r)]''_(r=r_*)`.

The existing root certificate gives the last unit derivative greater than `-19`. The amplitude bounds below give `A_lambda<1020`, so all amplitudes in (2) satisfy

`partial_t overline(Gamma)(0)>1950/64-(19/1000)1020`

`                                      =8871/800>11`.       (4)

The outer positive maximum exceeds the core maximum, as before. Evaluating the evolving mean at any one initial maximizing circle proves `sup overline(Gamma)(t)>sup overline(Gamma)(0)` for sufficiently small positive t. No differentiability or acceleration of the moving global maximizer is asserted.

## 4. Pressure budget and exact neutrality

Define

`R0=integral r[(e_new')^2+(400+16/r^2)e_new^2]dr`,

`Z0=integral q^2 dz`.

The established `|chi'|<4`, total variation two of e_new, and `log6<9/5` give

`R0 <= (6/25)*80*2+400[(6/25)^2-(1/25)^2]/2+16*(9/5)=392/5`,

`Z0<=12/5`, `||w_new||_2^2=pi R0 Z0`.

The horizontal tensor Q is positive on the seed cone and its largest horizontal eigenvalue is at most `3/[pi(17/5)^5]`. Therefore

`0<C_new<=3(5/17)^5 R0 Z0<=1764000/1419857<3/2`.       (5)

Using the existing certified bounds for C_v and `lambda^2<=1/16` yields

`900<A_lambda<1020`.

Thus the square root in (2) is positive, and its exact definition preserves `p_zz(0)=-8`. Angular selection cancels all base/seed bilinear pressure terms. The unchanged central affine neighborhood and C4 symmetry preserve `b(0)=Omega_core(0)=1`, `b'(0)=Omega_core'(0)=2`, and `beta'(0)=0`. This statement does not determine beta'' for the new seed.

## 5. Actual net fluctuation-energy gain, including viscosity

Let `K=||u-P0u||_2^2/2`, where P0 is azimuthal averaging. Initially `u-P0u=lambda w_new`. The full NS energy identity gives exactly

`K'(0)/(pi lambda^2)`

` =80 A_lambda I_r I_z-c R0 Z0-nu(R1 Z0+R0 Z1)`,       (6)

where

`I_r=-integral r C'(r)e_new(r)^2 dr`,

`I_z=integral q(z)^2q_v(z)dz`, `Z1=integral(q')^2dz`,

`R1=integral r[(e_new''+e_new'/r-(400+16/r^2)e_new)^2`

`                         +400(2e_new'+e_new/r)^2]dr`.

The pressure makes no contribution to this integrated solenoidal fluctuation energy. The fluctuation self-advection also integrates to zero. The positive first term is the exact work taken from the decreasing mean angular velocity. The other two terms retain both pump dilation and the complete three-dimensional viscous cost; the horizontal Laplacian identity used for R1 includes all cylindrical angular-basis derivatives.

New outward-rounded one-dimensional range quadratures certify

`I_r>0.02196178344685>1/50`,

`R1<175257.31003420<200000`.

The common axial plateaus give `I_z>=9/10`; the unchanged axial cutoff gives `Z1<=160/3`. Substitution into (6), with `A_lambda>900`, proves

`K'(0)/(pi lambda^2)`

` >80*900*(1/50)*(9/10)-(7/5)*(392/5)*(12/5)`

`       -(1/1000)[200000*(12/5)+(392/5)*(160/3)]`

` =205648/375>548`.       (7)

Since `K(0)<=pi lambda^2 (392/5)(12/5)/2`, it follows that

`K'(0)/K(0)>12853/2205>5`.       (8)

This is actual initial growth of nonaxisymmetric kinetic energy despite viscosity. Total kinetic energy still dissipates; the seed draws its gain from the mean. Smooth local existence makes both K and the mean maximum increase on some positive interval for each datum in (2), without supplying a numerical duration.

## Scope and the next test

This candidate corrects the specific initial energy-depletion problem of the previous trailing, flat-envelope seed. Its positive covariance and declining spatial envelope are both essential. The steep envelope is not a small WKB correction, so the earlier single-plane-wave Kelvin solution cannot be silently reused as its exact local evolution.

The [separate full-pressure calculation](leading-full-pressure-feedback.md) now proves the favorable central pressure derivative for this same new field (2). Still, the actual inherited envelope, stress divergence, pressure, and fluctuation energy must persist and return together. No subsequent rebuilding of the seed is implicit in this result.

`certify_leading_envelope.py` reproduces the interval jets, one-dimensional integrals, pressure budget, and rational margins. `leading-envelope-certificate.json` stores the exact dyadic enclosures and its existing root/cutoff-certificate dependencies. No PDE time integration or numerical pressure solve is used for (3)--(8).
