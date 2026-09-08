# Independent audit of the compact fourfold initial family

**Accepted.** The explicit seed, pressure retuning, central symmetry, norm-form parameter range, and strict initial receiver torque in rotational-receiver-torque.md and three-dimensional-initial-gate.md are mutually consistent. No substantive correction is needed. The independent continuity-constant audit in three-dimensional-gate-audit.md remains an explicit dependency.

## Seed support and symmetry

For \(L=1/4\), the envelope support satisfies
\[
\frac{11L}{16}\le r\le\frac{13L}{16},\qquad |z|\le\frac L{16}.
\]
Hence
\[
\frac{121}{256}\le\frac{r^2+z^2}{L^2}
\le\frac{170}{256}.
\]
Both endpoints lie strictly inside the receiver's decreasing cutoff collar \((1/4,1)\). The seed is also inside the original affine ball of radius \(1/2\), while remaining a positive distance from the axis and origin. Flat cutoff endpoints allow its zero extension to be smooth in Cartesian coordinates.

The cylindrical curl of \(a(r,z)\cos(4\theta+kr)e_z\) is exactly
\[
w_r=-4a\sin(4\theta+kr)/r,\quad
w_\theta=-a_r\cos(4\theta+kr)+ka\sin(4\theta+kr),\quad w_z=0.
\]
It is divergence free and has zero cylindrical component means. Under inversion, \(r\) is unchanged, \(z\) changes sign and \(\theta\) gains \(\pi\); the even envelope and even angular frequency make the vector potential even, hence its curl odd. A quarter-turn changes the scalar phase by \(2\pi\), giving correct vector equivariance under \(C_4\). Its nonzero angular modes cannot be canceled by an axisymmetric base, so every nonzero seed amplitude gives nonaxisymmetric data.

## Exact cross-pressure cancellation and neutral tuning

Let \(\mathcal T_\vartheta\) be the action of a rotation about the \(z\)-axis on vector fields. Rotational covariance of the Newton operator gives
\[
\Pi(\mathcal T_\vartheta a,\mathcal T_\vartheta d)=\Pi(a,d)
\]
for the scalar \(zz\) Hessian functional at the origin. Every axisymmetric base field \(U\) is fixed by this action, whereas averaging \(\mathcal T_\vartheta w\) over a full rotation gives zero. By bilinearity,
\[
\Pi(U,w)=\Pi\!\left(U,\frac1{2\pi}
                     \int_0^{2\pi}\mathcal T_\vartheta w\,d\vartheta\right)=0.
\]
This argument applies to the overlapping affine core as well as the remote axisymmetric swirl. It does not wrongly assert that the bilinear source vanishes pointwise.

The remaining seed pressure is \(\varepsilon^2J_w\). Since the exterior swirl contributes \(-A^2C_v\), the positive-sign retuning
\[
A_\varepsilon^2=(H_0+\varepsilon^2J_w)/C_v
\]
is correct and gives \(p[u_\varepsilon]_{zz}(0)=-8\) exactly.

Oddness fixes the central velocity at zero. A real matrix commuting with a quarter-turn has a transverse block \(aI+\Omega J\), no transverse–axial entries, and an axial entry \(d\). Incompressibility imposes \(d=-2a\), so its form is \(bD+\Omega J\). A symmetric pressure Hessian commuting with the same rotation has equal transverse diagonal entries and no mixed entries. Uniqueness preserves these symmetries.

The seed and retuning do not change the initially affine neighborhood. There \(\Delta u_0=0\) and \(\Delta p_0=-\operatorname{tr}[(D+J)^2]\) is constant, so \(\Delta u_t(0)=0\) too. The initial viscous curvature terms and their first time derivatives vanish. Consequently the central identities \(b'=\Omega'=2\), \(\beta'=0\), and
\(\beta''=-(p_{zz}'+32)/2\) remain exact for this fourfold family. Global axisymmetry is not needed for this particular central calculation.

## The stated nonzero parameter interval is sufficient

With \(W=\|w\|_{10}\), the established bilinear bound gives \(|J_w|\le128W^2\). The second restriction in \(\varepsilon_0\) ensures
\[
|\varepsilon^2J_w|\le H_0/2,
\]
so \(A_\varepsilon\) is real and strictly positive. Rationalizing the square root gives
\[
|A_\varepsilon-A_0|
\le\frac{128A_0W^2}{H_0}\varepsilon^2.
\]
The retuned datum therefore satisfies
\[
\|u_\varepsilon-u_*\|_{10}
\le W|\varepsilon|+D\varepsilon^2.
\]
The final two restrictions in \(\varepsilon_0\) allocate at most \(d_*/2\) to each term. Thus the *entire retuned field*, rather than the seed alone, lies in the audited continuity neighborhood. The pressure-derivative half-margin and
\(\beta_\varepsilon''(0)\ge290649/35000\) follow. All norms are finite and all denominators positive. No decimal size or practically useful lower bound on \(\varepsilon_0\) is implied.

## Receiver torque, its strict sign, and its limitation

The compact test field \(\phi_L=\chi_LJx\) is divergence free. Integration by parts in the actual NS equation removes pressure and yields the stated receiver identity. Its nonlinear contraction is
\[
Q_L(u)=\frac2{L^2}\int\chi'(|x|^2/L^2)
                         (u\cdot x)(u\cdot Jx)\,dx.
\]
For cylindrical angular mean \(U\) and mean-zero fluctuation \(v\), both cross terms vanish after angular integration. The viscous test field \(\Delta\phi_L\) is axisymmetric and azimuthal, so its pairing with \(v\) vanishes too.

The seed covariance is exactly
\(\overline{w_rw_\theta}=-mk a^2/(2r)\), including the sign of the radial phase derivative. Incorporating the cylindrical volume factor gives
\[
Q_L(w)=-\frac{2\pi mk}{L^2}
 \int\chi'((r^2+z^2)/L^2)r^2a^2\,dr\,dz>0.
\]
The positivity is strict because \(mk>0\), the envelope is nonzero on an open set, and its whole support lies in the strictly decreasing collar.

The axisymmetric base is exactly \(Dx+Jx\) throughout the receiver. Its viscous pairing is zero, its normalized initial rotation is one, and its receiver derivative is two. Thus
\[
\Omega_L'(0)=2+\varepsilon^2Q_L(w)/(2I_L)>2
\]
for every allowed nonzero parameter. Changing only the remote swirl amplitude does not affect this localized identity.

The geometric inequality behind
\(|Q_L(v)|\le C_\chi\|v\|_2^2\) is valid because \(x\perp Jx\) and \(|Jx|\le|x|\). Its integrated energy cost concerns only the fluctuation's receiver contribution; mean transport and viscosity are not assigned to it. The notes correctly distinguish this initial transfer from sustained transfer, whole-field return, or repeatability. They also explicitly retain the seed's added initial energy and the absence of an initial fluctuation torque inside its central hole.

The two supplied symbolic checkers were independently read and rerun successfully. This audit leaves the earlier certificates unchanged.
