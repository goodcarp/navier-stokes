# Independent audit of the twelve-coefficient pressure expansion

**Verdict:** the polarization factors, eight cubic terms, four displayed viscous terms, and neutral substitution in Section 4 of `full-pressure-feedback-gate.md` are correct. One locality justification should be explicit. For the particular radial core profiles already selected, two viscous coefficients vanish exactly, and the outer-swirl viscous coefficient is strictly positive. Thus there are ten potentially nonzero coefficients for those profiles, although the twelve-slot formula is valid more generally.

No profile integral has been numerically evaluated in this audit.

## 1. Normalizations and the nonlocal expansion

Use the source definitions
\[
 B(a,b)=\frac12P_L[(a\cdot\nabla)b+(b\cdot\nabla)a],\qquad
 \Pi(a,b)=\partial_{zz}\big[N*\operatorname{tr}(\nabla a\nabla b)\big](0).
\]
The trace makes \(\Pi\) symmetric; there is no missing factor of two in its definition. The usual pressure is \(p[u]=p[u,u]\). Hence
\[
 p_{zz}'=2\Pi(u_0,\nu\Delta u_0-B(u_0,u_0))
\]
has the correct polarization. The distributional contact term is included in \(\Pi\); these coefficient definitions must continue to use that full distribution, not only the integral away from zero.

For \(u_0=bS+\Omega W+cP+Av\), strict disjointness of the core and outer supports makes their cross **advection sources** vanish before Leray projection. Thus \(B(S,P)=B(S,v)=B(W,P)=B(W,v)=0\) is valid despite the nonlocal projection.

The other generated fields need different treatment:

* \(B(S,S),B(W,W),B(P,P),B(v,v)\) are poloidal and generally have nonlocal pressure tails. The cross-region terms such as \(\Pi(P,B(S,S))\) must be retained; the table correctly retains them.
* \(B(S,W)\) and \(B(P,v)\) are azimuthal. An axisymmetric azimuthal vector field is divergence free, so Leray acts as the identity on these mixed sources. They therefore retain compact support in the core and outer region, respectively. This justifies the less immediate cancellations
  \[
  \Pi(v,B(S,W))=\Pi(W,B(P,v))=0.
  \]
  Without this support statement, discarding the prospective \(b\Omega A\) and \(c\Omega A\) terms merely because the original fields were disjoint would be unjustified.

For a poloidal field \(a=a_r e_r+a_z e_z\) and azimuthal field \(w=w_\theta e_\theta\), their cylindrical gradient matrices have the forms
\[
 \nabla a=\begin{pmatrix}a_{r,r}&0&a_{r,z}\\0&a_r/r&0\\a_{z,r}&0&a_{z,z}\end{pmatrix},\qquad
 \nabla w=\begin{pmatrix}0&-w_\theta/r&0\\w_{\theta,r}&0&w_{\theta,z}\\0&0&0\end{pmatrix}.
\]
Their product has zero trace identically. This proves every pressure parity cancellation used in the expansion, even if the poloidal field has a nonlocal tail. The apparent coordinate singularity is resolved by smoothness on the axis.

Expanding with these rules gives exactly the eight cubic coefficients in the table. In particular,
\[
 t_{Wb}=-2\Pi(S,B(W,W))-4\Pi(W,B(S,W)),
\]
and
\[
 t_{vc}=-2\Pi(P,B(v,v))-4\Pi(v,B(P,v)).
\]
The factors four arise from the two ordered mixed terms in \(B(u_0,u_0)\). Both terms in \(t_{vc}\) are required; pressure parity does not discard either one.

## 2. Viscous terms: the table is valid but can be sharpened

The Laplacian preserves poloidal/azimuthal character and the support of a smooth compact initial field. Therefore all cross terms in \(2\nu\Pi(u_0,\Delta u_0)\) vanish by parity or core/outer disjointness. Its four diagonal entries are exactly
\[
 \nu\big(d_Sb^2+d_W\Omega^2+d_Pc^2+d_vA^2\big),
 \qquad d_j=2\Pi(j,\Delta j).
\]

**For the chosen radial vector-potential core, \(d_S=d_W=0\).** This is an exact identity and needs no profile integration beyond the already established universal core pressure formulas.

For the unit strain matrix \(S_0=\operatorname{diag}(-1,-1,2)\), write
\[
 S_\psi=\operatorname{curl}\left[-\frac{\psi(|x|^2)}3x\times(S_0x)\right].
\]
The quadratic vector polynomial \(x\times(S_0x)\) is harmonic and homogeneous of degree two. Consequently
\[
 \tag{1} \Delta S_\psi=S_\eta,\qquad
 \eta(s)=4s\psi''(s)+14\psi'(s).
\]
Because \(\psi=1\) on a neighborhood of zero, \(\eta=0\) there. The perturbed cutoff \(\psi+t\eta\) is still smooth, compact, and identically one near zero. The universal pressure identity therefore gives
\[
 \Pi(S_\psi+t\Delta S_\psi,S_\psi+t\Delta S_\psi)=-18/7
\]
for every such parameter \(t\). Differentiation at zero proves \(d_S=0\). This is a profile variation used to evaluate a quadratic functional, not an assumption that the NS evolution follows that variation.

Similarly, the unit rotating profile is \(W=h(s)Jx\), with \(h=\psi+2s\psi'/3\). Direct differentiation gives
\[
 \tag{2} \Delta W=(4s h''+10h')Jx.
\]
The new radial coefficient is zero near the origin. The universal pressure formula \(p_{zz}[h(s)Jx](0)=2/5\) for compact radial \(h=1\) near zero is unchanged along \(W+t\Delta W\). Hence \(d_W=0\). Arbitrary compact affine extensions need not share these cancellations; they apply to the actual radial profiles specified in `affine-core-pressure.md`.

**Also \(d_v>0\) for the selected outer swirl.** The support is separated from the origin, so pressure integration by parts gives, with \(Q_{ij}=-\partial_{zzij}N\),
\[
 \Pi(v,\Delta v)=-\int Q_{ij}v_i\Delta v_j
                =\int Q_{ij}\partial_kv_i\partial_kv_j.
\]
The last equality retains the tensor and uses \(\Delta Q_{ij}=0\) on the support. Every Cartesian derivative of \(v\) is horizontal, and the horizontal restriction of \(Q\) is positive definite on the support cone. Since the compact swirl is nonzero,
\[
 \tag{3} d_v=2\int Q_{ij}\partial_kv_i\partial_kv_j>0.
\]
This agrees with the negative viscous derivative of the outer coefficient \(C_v=-p_{zz}[v](0)\). No sign for \(d_P\) follows from this argument because its meridional derivatives can have vertical components where \(Q\) is not positive definite.

Accordingly, the sentence that no signs or values of any twelve coefficients have yet been established should be qualified: the chosen radial core gives two exact zeros, and the selected outer swirl gives one strict positive sign.

## 3. Explicit normalized polynomial and admissible domain

Let \(b>0\), \(\theta=\Omega/b\), \(\mu=c/b\), \(\varepsilon=\nu/b\), and define
\[
 H(\mu,\theta^2)=\frac{38/7+(2/5)\theta^2-C_P\mu^2}{C_v}.
\]
Neutral tuning is exactly \(A^2=b^2H\). It requires \(H\ge0\); for the chosen \(C_P<0,C_v>0\), this holds strictly for every real \(\mu,\theta\). The intended positively rotating, forward-pumping regime has \(\theta>0,\mu>0,\varepsilon\ge0\), plus the separately derived pump threshold if that additional property is required.

The full pressure gate \(p_{zz}'+32b^3<0\) becomes
\[
 \tag{4}
 \begin{aligned}
 \mathscr G={}&32+t_{300}+t_{210}\mu+t_{120}\mu^2+t_{030}\mu^3\\
 &+\theta^2(t_{Wb}+t_{Wc}\mu)
   +H(t_{vb}+t_{vc}\mu)\\
 &+\varepsilon\,[d_S+d_W\theta^2+d_P\mu^2+d_vH]<0.
 \end{aligned}
\]
For the selected profiles one can set \(d_S=d_W=0\). This is a polynomial of degree at most three in \(\mu\), linear in \(\theta^2\), and linear in \(\varepsilon\), with mixed products permitted. The neutral substitution introduces no square-root dependence on \(A\), since all its surviving terms are even in \(A\).

For reference, the tuned combinations multiplying the six cubic monomials \(1,\mu,\mu^2,\mu^3,\theta^2,\mu\theta^2\) are
\[
 \begin{aligned}
 &t_{300}+\frac{38}{7C_v}t_{vb},\qquad
 t_{210}+\frac{38}{7C_v}t_{vc},\\
 &t_{120}-\frac{C_P}{C_v}t_{vb},\qquad
 t_{030}-\frac{C_P}{C_v}t_{vc},\\
 &t_{Wb}+\frac{2}{5C_v}t_{vb},\qquad
 t_{Wc}+\frac{2}{5C_v}t_{vc}.
 \end{aligned}
\]
The additive \(32\) belongs to the gate, not to a pressure coefficient. The three viscous combinations multiplying \(\varepsilon\), \(\varepsilon\mu^2\), and \(\varepsilon\theta^2\) are
\[
 d_S+\frac{38}{7C_v}d_v,\qquad
 d_P-\frac{C_P}{C_v}d_v,\qquad
 d_W+\frac{2}{5C_v}d_v.
\]
Although the first and third are strictly positive for the chosen profiles, the middle one remains unevaluated. This does not establish a sign for the full viscous contribution or for \(\mathscr G\).

The accompanying checker expands all ordered bilinear terms under the proved support/parity rules, checks the coefficient factors and neutral polynomial, and verifies the two Laplacian identities used for the exact core zeros. It does not evaluate any nonlocal profile integral or test the sign of the remaining gate.
