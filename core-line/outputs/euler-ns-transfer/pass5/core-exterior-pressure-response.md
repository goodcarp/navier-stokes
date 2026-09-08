# Exact affine-core response to exterior harmonic pressure

For the radial affine core used in pass 4, the contribution to the initial pressure derivative from the exterior pressure acceleration is
\[
\boxed{\mathcal C(-\nabla p_{\rm out})
 =\frac{36}{7}bH_2-\frac85 bM_\psi H_4
 =\frac{18}{7}b\,p_{{\rm out},zz}(0)
    -\frac{bM_\psi}{15}p_{{\rm out},zzzz}(0).}
\tag{1}
\]
Here \(p_{\rm out}\) is axisymmetric and harmonic on a neighborhood of the **entire core support**,
\[
p_{\rm out}(r,\mu)=\sum_{\ell\ge0}H_\ell r^\ell P_\ell(\mu),
\qquad \mu=z/r,\qquad M_\psi=\int_0^\infty\psi(s)\,ds.
\tag{2}
\]
The degree-four term survives radial integration. Degrees above four vanish exactly, and the constant mode contributes nothing. Formula (1) requires only one cutoff moment; it removes the core spatial quadrature for this exterior-pressure part of the acceleration. It does not evaluate the remaining parts of the full feedback gate or certify a numerical sign.

## Definitions and the contact term

Let \(S_0=\operatorname{diag}(-1,-1,2)\), \(Jx=(-y,x,0)\), and let \(\psi\) be smooth, compactly supported in the squared-radius variable, and equal to one near zero. Write
\[
\begin{aligned}
S_\psi(x)&=\left(\psi+\frac{2|x|^2}{3}\psi'\right)S_0x
               -\frac23\psi'(x\cdot S_0x)x,\\
W_\psi(x)&=\left(\psi+\frac{2|x|^2}{3}\psi'\right)Jx,
\qquad u_c=bS_\psi+\Omega W_\psi.
\end{aligned}
\]
Every cutoff derivative is with respect to its argument \(s=|x|^2\). With \(N(x)=1/(4\pi|x|)\), define the bilinear pressure and core functional by
\[
\Pi(a,w)=\partial_{zz}\!\left[N*
                   \operatorname{tr}(\nabla a\nabla w)\right](0),
\qquad \mathcal C(w)=2\Pi(u_c,w).
\]
For \(w=-\nabla\phi\), put
\[
h_c=-2\operatorname{tr}(\nabla u_cD^2\phi).
\]
The exact distributional formula is
\[
\mathcal C(-\nabla\phi)
 =-\frac13h_c(0)+
   \frac12\int_0^\infty\frac{dr}{r}
       \int_{-1}^1(3\mu^2-1)h_c(r,\mu)\,d\mu.
\tag{3}
\]
The azimuthal velocity \(W_\psi\) has zero contraction with the Hessian of any axisymmetric scalar. Hence only \(bS_\psi\) enters (3). This cancellation retains the full exterior pressure; it is not a local pressure approximation.

## Exact angular and radial evaluation

Set
\[
T_\ell=\operatorname{tr}\!\left(\nabla S_\psi\,
                  D^2[r^\ell P_\ell(\mu)]\right),
\quad
I_\ell(r)=\int_{-1}^1(3\mu^2-1)T_\ell(r,\mu)\,d\mu.
\]
The degree-four angular selection identity proved in
work/pass4/angular-residual-certificate.md gives \(I_\ell=0\) for every \(\ell>4\), and odd degrees vanish by parity. Direct substitution of \(f(r)=r^\ell\) into its three exact weights gives
\[
\begin{aligned}
I_0(r)&=0,\\
I_2(r)&=\frac{32}{35}r^4\psi''(r^2)
                   +\frac{16}{5}r^2\psi'(r^2),\\
I_4(r)&=\frac{64}{7}r^6\psi''(r^2)
                   +\frac{1536}{35}r^4\psi'(r^2)
                   +\frac{144}{5}r^2\psi(r^2).
\end{aligned}
\tag{4}
\]
The companion checker independently obtains (4) from the Cartesian gradient and solid-harmonic Hessians, rather than assuming the angular weights.

Compact support and \(\psi(0)=1\) give, by \(s=r^2\) and integration by parts,
\[
\begin{array}{c|c}
\text{radial integral}&\text{value}\\ \hline
\int_0^\infty r^3\psi''(r^2)\,dr&1/2\\
\int_0^\infty r\psi'(r^2)\,dr&-1/2\\
\int_0^\infty r^5\psi''(r^2)\,dr&M_\psi\\
\int_0^\infty r^3\psi'(r^2)\,dr&-M_\psi/2\\
\int_0^\infty r\psi(r^2)\,dr&M_\psi/2
\end{array}
\tag{5}
\]
Consequently,
\[
\int_0^\infty I_2(r)\frac{dr}{r}=-\frac87,\qquad
\int_0^\infty I_4(r)\frac{dr}{r}=\frac85M_\psi.
\tag{6}
\]
Because \(h_c=-2b\sum_\ell H_\ell T_\ell\), the principal-value contribution is
\[
\frac87bH_2-\frac85bM_\psi H_4.
\]
The quadratic harmonic is \(r^2P_2(\mu)=z^2-(x^2+y^2)/2\), whose Hessian is \(S_0\). Thus \(h_c(0)=-12bH_2\), and the contact contribution in (3) is **\(4bH_2\)**. Adding it proves (1). Smooth solid harmonics of degree greater than two have zero central Hessian, so there are no omitted higher-degree contact terms.

For a smooth harmonic function on a neighborhood of the compact support, the regular solid-harmonic expansion and its required derivatives converge uniformly there. The mode-by-mode calculation therefore extends to the whole exterior pressure. If the exterior source is separated from the core, this harmonicity and convergence follow directly from the Newton kernel.

## Computing the two exterior multipoles

For a compact exterior pressure source \(g_{\rm out}\), use the conventional coefficients
\[
g_{{\rm out},\ell}(\rho)=
\frac{2\ell+1}{2}\int_{-1}^1
                g_{\rm out}(\rho,\mu)P_\ell(\mu)\,d\mu.
\]
The regular radial solution of \(-\Delta p_{\rm out}=g_{\rm out}\) inside the source-free core is
\[
H_\ell=\frac1{2\ell+1}
          \int_0^\infty \rho^{1-\ell}
                       g_{{\rm out},\ell}(\rho)\,d\rho.
\tag{7}
\]
There is no singularity in this formula because the exterior source avoids the origin. In particular,
\[
H_2=\frac15\int\rho^{-1}g_{{\rm out},2}(\rho)\,d\rho,
\qquad
H_4=\frac19\int\rho^{-3}g_{{\rm out},4}(\rho)\,d\rho.
\]
These multipoles still require exact evaluation or rigorous enclosure. If neutral tuning has been proved exactly, it can determine \(H_2\) through the initial central Hessian; it does not determine \(H_4\).

The dimensions are consistent: \(M_\psi\) has dimensions of length squared, and \(H_4M_\psi\) has the same dimensions as \(H_2\). For a nonnegative cutoff, \(M_\psi>0\), so the sign of the degree-four response is opposite to \(bH_4\). Arbitrary smooth cutoffs need not have a positive moment.

This is an exact identity for the initial functional. It supplies no later-time preservation of a radial core, no repeated amplification stage, and no Navier–Stokes singularity theorem.

The file verify_core_exterior_pressure_response.py checks the Cartesian angular contractions, the contact terms, all five radial moments, and both forms of (1) using exact symbolic arithmetic. Its additional compact polynomial profiles verify the integrations algebraically; they are not substitutes for the smooth cutoff assumed in the proof.
