# Independent audit of the three-dimensional pressure-gate continuity estimate

The proposed constants are valid in the precise derivative-sum Sobolev norm used in pass 3. The estimate controls the full whole-space initial pressure derivative for arbitrary smooth divergence-free three-dimensional data. No symmetry or compact-support assumption is needed for the norm inequality itself. Its application to the neutral strain-ratio gate must also retain the neutral pressure condition, or separately control the off-neutral correction.

## Norm convention and the pressure bilinear form

Use

\[
\|a\|_{H^s}^2=\sum_{|\alpha|\le s}\sum_{j=1}^3
\|\partial^\alpha a_j\|_{L^2(\mathbb R^3)}^2.
\]

For a scalar field use the same norm without the component sum. Pass 3 proves the scalar embedding `||h||_infinity<=||h||_H2`; its unitary-Fourier constant is in fact at most `1/(2 sqrt(pi))<1`. The multiplier

\[
R_{zz}:=\partial_{zz}(-\Delta)^{-1},\qquad
\widehat{R_{zz}h}(\xi)=-\frac{\xi_z^2}{|\xi|^2}\widehat h(\xi)
\]

is a contraction on every one of these Sobolev norms. Defining its value at frequency zero is immaterial to this `L^2` multiplier. Therefore

\[
\Pi(a,b):=\left[R_{zz}\operatorname{tr}(\nabla a\nabla b)\right](0)
\]

is a bounded bilinear form once its source is in `H2`. This definition includes the entire distributional Hessian, including the contact term; it is not a principal-value-only approximation. The source is symmetric under interchanging `a,b` by the matrix trace identity.

For each multiindex `|alpha|<=2`, the differentiated source has nine component products and total Leibniz weight `2^|alpha|<=4`. Every such product has derivative orders `1+|beta|` and `1+|alpha-beta|`. Their sum is at most four, so one factor has derivative order at most two and can be put in `L-infinity` using its `H4` norm. The other factor has at most three derivatives and can be put in `L2`. Consequently

\[
\|\partial^\alpha\operatorname{tr}(\nabla a\nabla b)\|_{L^2}
\le36\|a\|_{H^4}\|b\|_{H^4}.
\]

There are `binomial(5,3)=10` output multiindices, yielding

\[
\boxed{|\Pi(a,b)|\le36\sqrt{10}\|a\|_{H^4}\|b\|_{H^4}
<128\|a\|_{H^4}\|b\|_{H^4}.}
\tag{1}
\]

Indeed, `(36 sqrt(10))^2=12960<16384=128^2`. No missing vector-component factor is required: the nine source products have already been counted, and each component norm is bounded by the full vector norm.

## Leray nonlinearity and telescoping the actual derivative

For `B(a,b)=P_Leray(a dot grad b)`, the pass 3 estimate is

\[
\|B(a,b)\|_{H^4}\le512\|a\|_{H^6}\|b\|_{H^6},
\qquad \|\Delta a\|_{H^4}\le3\|a\|_{H^6}.
\tag{2}
\]

The first constant exceeds the elementary sufficient value `3*2^4*sqrt(binomial(7,3))=48 sqrt(35)<512`. The Leray multiplier is an `L2` contraction commuting with derivatives. If the symmetric definition of `B` from the coefficient-expansion note is used instead, the same bound holds and `B(u,u)` is unchanged.

For the full initial pressure derivative define

\[
F(u)=2\Pi(u,\nu\Delta u-B(u,u)).
\]

For smooth divergence-free data this equals `partial_t p_zz(0,0)` under the actual unforced NS acceleration. Put `e=u-v`, `d=||e||_H10`, `X=||u||_H10`, and `Y=||v||_H10`. The viscosity difference is

\[
2\nu[\Pi(e,\Delta u)+\Pi(v,\Delta e)],
\]

and the nonlinear difference is minus twice

\[
\Pi(e,B(u,u))+\Pi(v,B(e,u))+\Pi(v,B(v,e)).
\]

These are exact telescoping identities; no symmetry of the unsymmetrized `B` is assumed. Equations (1) and (2) imply precisely

\[
\boxed{|F(u)-F(v)|\le256d\left[
3\nu(X+Y)+512(X^2+XY+Y^2)\right].}
\tag{3}
\]

Thus, around a fixed datum `u_*` of norm `X0`, restricting `d<=1` allows `X,Y<=Z=X0+1` and gives the requested Lipschitz constant

\[
\boxed{C=256[6\nu Z+1536Z^2],\qquad
|F(u)-F(u_*)|\le C\|u-u_*\|_{H^{10}}.}
\tag{4}
\]

The derivative requirements actually close with `H6` at this step, but using the existing `H10` norm is valid and keeps the same topology as the local-time estimates.

## What the norm neighborhood preserves

The pass 5 neutral datum has the certified scalar margin

\[
F(u_*)+32\le-\delta,\qquad
\delta=\frac{290649}{8750}>33.
\]

Therefore any divergence-free datum with

\[
\|u-u_*\|_{H^{10}}\le
\min\left\{1,\frac{\delta}{2C}\right\}
\tag{5}
\]

satisfies `F(u)+32<=-delta/2`. This is an explicit **norm-form** neighborhood. Without numerical enclosures of `X0` and the chosen seed's `H10` norm, (5) is not a certified decimal perturbation amplitude or a useful-size asymmetry statement.

For an unretuned insertion `u=u_*+epsilon w`, (5) follows from `|epsilon| ||w||_H10` satisfying its right side. But this alone need not retain the neutral condition `p_zz(0)=-8`, even when `w` is supported away from the unchanged affine core. An `m=4` perturbation has zero linear cross contribution to the scalar pressure Hessian at the axis by Fourier orthogonality, while its quadratic self contribution can change it by `epsilon^2 Pi(w,w)`.

If the seed/base geometry supplies this orthogonality, exact retuning of the exterior swirl can use

\[
A_\epsilon^2=A_0^2+
\frac{\epsilon^2\Pi(w,w)}{C_v}.
\tag{6}
\]

The positive radicand must be justified. The norm distance in (5) is then the distance of the **full retuned datum**, bounded by

\[
|\epsilon|\,\|w\|_{H^{10}}
+|A_\epsilon-A_0|\,\|v_{\rm outer}\|_{H^{10}}.
\tag{7}
\]

The sign in (6) follows from `p_zz[u_*]=-8` and the exterior contribution `-A^2 C_v`. Alternatively, one can keep the unretuned datum and include the appropriate off-neutral terms in the actual central quotient derivatives. Equations (3)--(5) by themselves retain the pressure-derivative scalar margin, not every central dynamical identity used at the neutral endpoint.

Discrete rotational symmetry and oddness, if imposed on the chosen compact seed, may preserve the required central matrix form by uniqueness. Those symmetry facts and any pressure retuning are separate from this general functional-continuity estimate. Nothing here shows persistence of nonaxisymmetric torque or a returning state.
