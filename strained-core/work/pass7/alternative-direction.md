# PASS7: chirality, core orientation, and a certifiable pressure bound

This concerns the actual PASS6 outer seed and the same exactly retuned datum. Its finite mean-maximum gain is already proved. The aim here is to make the **full initial pressure-derivative** test manageable without replacing any part of the initial field.

Write

`u=S_Phi+sigma W_psi+A v+lambda w_(m,k)`, `sigma in {+1,-1}`,

`w_(m,k)=curl[f e_z]`, `f=e(r)q(z) cos(m theta+kr)`,

with the unchanged envelope, `m=4`, `k=20`, and exact pressure tuning

`A^2 C_v=204/35-lambda^2 C_w`.

Both choices of core-rotation sign have the same radial base pressure certificate and the same initial neutrality. The meaningful central ratio in the negative-rotation case is `b/|Omega|`.

## 1. The carrier reflection does not freely select the dangerous sign

For the physical reflection `R(x,y,z)=(x,-y,z)`, push-forward of velocity gives

`R_* S_Phi=S_Phi`, `R_* v=-v`, `R_* W_psi=-W_psi`,

`R_* w_(m,k)=-w_(m,-k)`.

The last minus sign matters: the vector potential is axial in the curl formula, and a spatial reflection reverses curl orientation. It also follows directly from the cylindrical components

`w_r=-m e q sin(m theta+kr)/r`,

`w_theta=-e_r q cos(m theta+kr)+k e q sin(m theta+kr)`.

The scalar `zz` pressure pairing at the origin and the NS projection commute with this reflection. For the actual interaction `L(U;w)` defined in `outer-pressure-compatibility-test.md`, therefore,

`E(k):=L(v;w_k)` is odd in `k`,

`F(k):=L(W_psi;w_k)` is odd in `k`,

`L(S_Phi;w_k)`, `C_w(k)`, and `2nu Pi(w_k,Delta w_k)` are even in `k`.

Changing the seed's sign or its constant angular phase does not change these quadratic seed interactions. Taking `(-m,-k)` produces the identical scalar potential. Changing `k` alone changes the useful mean torque's sign, since that torque is proportional to `mk`. Changing both the outer swirl sign `A` and `k` makes the negative mean minimum grow in absolute value, but leaves `A E(k)` unchanged. Thus these orientation symmetries do **not** independently choose the sign of `A E` while preserving the same signed maximum-growth objective.

For the **actual separated outer seed**, there is a stronger cancellation: `F=L(W_psi;w)=0` exactly. The local mixed advections vanish by disjoint support. The remaining exterior-pressure acceleration pairs an axisymmetric azimuthal core with a scalar Hessian; angular averaging selects the axisymmetric pressure, whose Hessian contraction with `grad W_psi` vanishes. Thus changing the core's rotation sign also cannot modify this pressure test. The general odd-parity statement alone would not have established this zero.

## 2. The local part of E already has the favorable chirality

Let `N=(4pi|x|)^(-1)`, `Q=-D^2 partial_zz N`, and let `U=v_theta>=0` denote the unit axisymmetric outer swirl speed. On the seed support,

`Q_theta=3(4z^2-r^2)/(4pi(r^2+z^2)^(7/2))`,

`partial_r Q_theta=15r(r^2-6z^2)/(4pi(r^2+z^2)^(9/2))<0`,

`Q_rr-Q_theta=r partial_r Q_theta`.

After integration by parts, the cubic convective pressure term is `-int (u dot grad Q):(u tensor u)`. Its coefficient containing one `v` and two `w` factors is exactly

`E_local=6pi m k int e(r)^2 q(z)^2 U(r,z) partial_r Q_theta(r,z) dr dz < 0`.       (1)

Indeed the cylindrical cubic integrand contributing this coefficient is `-6 U partial_r Q_theta w_r w_theta`; its angular average uses `overline(w_r w_theta)=-mk e^2 q^2/(2r)`. Thus the sign compatible with the positive outer mean-maximum gain already favors the local central-pressure contribution. The remaining nonlocal pressure piece, rather than the convective piece, requires enclosure.

## 3. Exact whole-space projection bounds use no seed second derivatives

Let `p_mix=2Pi_field(v,w)`, where `-Delta Pi_field(v,w)=tr(grad v grad w)`. Since both fields are solenoidal,

`tr(grad v grad w)=div[(w dot grad)v]`.

Put `G=(w dot grad)v`. Then

`p_mix=2 N*div G`, `grad p_mix=-2 P_grad G`,

where `P_grad=I-P` is the orthogonal L2 gradient projection. The part pairing `Qv` with the seed-only pressure integrates to zero in the azimuthal variable. Hence the complete nonlocal coefficient is

`E_pressure=2 int Qw dot grad p_mix=-2 int d p_mix`,

`d=div(Qw)`.

It follows that

`|E_pressure| <= 4 ||P_grad(Qw)||_2 ||G||_2`

`                 <= (4/m) ||r d||_2 ||G||_2`.       (2)

The last step uses the exact nonzero angular mode: `||p/r||_2<=||grad p||_2/m`. It is a whole-space estimate, with no artificial boundary condition. In particular, `||grad p_mix||_2<=2||G||_2` uses derivatives of the unit swirl `v`, not second derivatives of the oscillatory seed.

For `B=U/r`, angular integration gives

`||G||_2^2=pi int r [B^2(e_r^2+k^2 e^2)q^2+U_r^2 m^2 e^2 q^2/r^2] dr dz`.

An optional explicit pressure trial is `p_app=-2B f`. Exactly,

`G=grad(Bf)+H`,

`H_r=-B_r f`, `H_theta=B_r partial_theta f`, `H_z=-partial_z(Bf)`.

This yields the computable leading pairing

`E_app=-4pi m k int (Q_rr-Q_theta) B e^2q^2 dr dz`,

and error at most `(4/m)||r d||_2||H||_2`. It is an exact splitting, not a claim that the trial approximates pressure well.

Compact-support midpoint quadrature in `evaluate_alternative_budget.py` gives exploratory values

`E_local ~ -0.00477244`, `||r d||_2 ~ 0.00333918`, `||G||_2 ~ 14.6846`.

Consequently the bound (2) is about `0.04903`; it is too broad for the desired `0.02` pressure allowance. The direct `4||Qw||||G||` bound is about `0.967` and is substantially worse. The simple `-2Bf` trial does not adequately fix this. These diagnostic integrals are not used as rigorous enclosures.

A stronger **explicit adjoint residual certificate** can improve (2). For any smooth finite-energy mode-m trial `eta`, let `h=d+Delta eta`. The exact identity is

`E_pressure=4 int grad eta dot G+4 int grad(N*h) dot G`,

so

`|E_pressure-4 int grad eta dot G| <= (4/m)||r h||_2||G||_2`.       (3)

Only the displayed trial integral and residual norm require rigorous quadrature. A pressure-side trial has the analogous energy-error formula. This provides a concrete route to certify the interaction without treating agreement of finite-box pressure solves as a whole-space proof.

## 4. The requested norm certificate

For the actual seed, write `Qdiff=Q_rr-Q_theta`. Its exact complex mode gives

`|r d|^2=m^2[(k Qdiff e q)^2+(Qdiff(e_r-e/r)q+Q_rz e q_z)^2]`,

`Qdiff=15r^2(r^2-6z^2)/(4pi(r^2+z^2)^(9/2))`,

`Q_rz=-15rz(4z^2-3r^2)/(4pi(r^2+z^2)^(9/2))`.

The full spatial norm squared is `2pi int_upper r |r d|^2 dr dz`: one factor `pi` comes from the real angular mode, and the factor two accounts for the reflected packet. `certify_outer_adjoint_norm.py` evaluates this expression by outward-rounded tensor-box range quadrature using the existing audited cutoff-jet engine. Its JSON output preserves exact dyadic bounds; floating-point adaptivity only chooses boxes and does not determine enclosure validity. The 18,000-box result proves

`||r d||_2^2 <= 9923203144385424509210371461965 * 2^(-119)`,

and therefore `||r d||_2<773/200000=0.003865`. It does not reach the originally requested `0.0036` target; no such claim is made. The parent certificate already succeeds with a weaker bound, so further refinement is unnecessary.

## 5. Frequency changes have a pressure floor and a viscous cost

For a fixed envelope, write

`C_w=m^2 C_r+k^2 C_k+C_e`,

where `C_r=pi int Q_rr e^2q^2/r`, `C_k=pi int r Q_theta e^2q^2`, and `C_e=pi int r Q_theta e_r^2q^2`. These are positive. Fixing the plateau torque strength `K=lambda^2 m k>0` gives

`lambda^2 C_w=K[C_r/rho+rho C_k+C_e/(m^2 rho)]`, `rho=k/m`.

Its exact minimum over positive real carrier ratios is

`2K sqrt[C_k(C_r+C_e/m^2)]`,

at `rho=sqrt[(C_r+C_e/m^2)/C_k]`. Thus taking arbitrarily large radial frequency at fixed torque is not a way to make the pressure budget disappear. For a narrow radial band where the horizontal entries of Q are similar, the balanced carrier is approximately `k=m/r`; the present `k=20`, `m=4`, `r~0.21` already lies near this simple balance.

There is also a strictly adverse viscous term. Harmonicity of `Q` on the seed support and integration by parts give exactly

`2nu Pi(w,Delta w)=2nu int Q_ij partial_l w_i partial_l w_j >0`.

The horizontal block of Q is positive there. At fixed torque the angular derivatives grow with frequency, so the inviscid envelope optimization is not the full optimization for ordinary NS.

## 6. A separate chiral compensation parameter, if E needs adjustment

There is an exact same-datum alternative to changing the primary carrier's sign. Add a second compact mode `m=12`, supported away from every initial mean-maximizing circle, while retaining the primary mode `m=4`. Both are odd and C4-equivariant, and a rotation by `pi/4` negates both seeds. Distinct modes remove every base/seed/seed cross term; this rotation also removes every pure cubic seed term. Neutral pressure tuning and the full pressure correction therefore split into the sum of the two quadratic seed contributions.

The auxiliary chirality can be reversed without changing the primary torque at its maximizing circles. If it too is separated from the core, `F_aux=0`, and its odd correction `A E_aux` can be made nonpositive by choosing its carrier sign. Its even meridional and viscous terms and its positive `C_aux` retuning cost remain. This is a concrete extra parameter, not an automatic favorable result: one must check that those even costs are smaller than the favorable odd term and that the remaining exact pressure budget is positive. An auxiliary seed is unnecessary if the present one passes the requested pressure certificate.

`verify_alternative_direction.py` checks the reflection/component identities, kernel and local-cubic identities, frequency optimization, and the nonresonant mode selection. It does not replace the norm enclosure or the full coefficient certificate.
