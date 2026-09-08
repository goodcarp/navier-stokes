# Independent audit of the initial periodic-image tail

Reviewed `verify_periodic_image_tail.py`, `initial_periodic_images.py`, `initial-periodic-images.md`, the exact stored pass6 cutoff-derivative endpoints, and the pass8 seed-energy bounds. No retained-image quadrature was rerun. Independent symbolic differentiation of the Newton kernel and parsing of the exact cutoff endpoints passed.

**Verdict:** the omitted-image estimate and its coarse energy bound are valid for the unchanged initial datum \(A=1020,\lambda=1/4,c=7/5\). They do not enclose the retained-image quadratures or validate a positive-time periodic-to-whole-space correction. One wording correction was made with the parent's authorization: seed/mean orthogonality holds after angular integration, not pointwise in angle.

## Kernel and image-count factors

Write \(q=r^2+z^2\). Independent differentiation of
\(4\pi\partial_{zz}N=(2z^2-r^2)q^{-5/2}\) gives

\[
4\pi\partial_{zzrr}N=(105r^2z^2-12q^2)q^{-9/2},
\]
\[
4\pi(\partial_r\partial_{zz}N)/r=3(r^2-4z^2)q^{-7/2},
\]
\[
4\pi\partial_{zzzz}N=(105z^4-90qz^2+9q^2)q^{-9/2},
\]
\[
4\pi\partial_{rzzz}N=15rz(4z^2-3r^2)q^{-9/2}.
\]

These agree with the certificate's dimensionless matrix. The script's `Qrr,Qtt,Qzz,Qrz` are the negatives of these Hessian entries, and its final subtraction restores the correct positive stress-kernel contraction. The two principal-minor identities and positive first pivots prove both \(24I-H\) and \(24I+H\) are positive semidefinite on \(0\le t\le1\). The azimuthal eigenvalue lies in \([-12,3]\). Hence the operator norm bound \(24/(4\pi R^5)\) is valid; the constant is attained on the axis.

For \(E=\tfrac12\int|u_0|^2\), one image contributes at most
\(12E/[\pi(nL-6)^5]\). Including both signs gives \(24E/\pi\) times the scalar sum. The omitted indices are \(n>N\), and

\[
\sum_{n>N}(nL-6)^{-5}
\le\int_N^\infty(Lx-6)^{-5}\,dx
=\frac1{4L(LN-6)^4}.
\]

Thus the stated tail \(6E/[\pi L(LN-6)^4]\) has the correct factors of two, the correct starting index, and the correct power of \(L\). The assumed \(LN>6\) suffices for the tail inequality; the initial disjoint-copy interpretation separately uses \(L>12\). No image contact term is present.

## Energy upper bound

For the actual radial meridional extension, spherical averaging followed by integration by parts gives

\[
E(S_\Phi)=\frac{8\pi}{15}\int_0^\infty s^{7/2}\Phi'(s)^2\,ds.
\]

The unintegrated angular-average speed squared is
\(2s\Phi^2+(8/5)s^2\Phi\Phi'+(8/15)s^3\Phi'^2\). In the energy integral the first two terms cancel because \(\int s^{5/2}\Phi\Phi'=-(5/4)\int s^{3/2}\Phi^2\). Compact support and the factor at zero justify both boundary terms.

The three derivative supports of \(\Phi\) are disjoint: \(s\in[1/4,1]\), \([25/16,25/4]\), and \([25,36]\). On them the derivative maxima are bounded by \(4\), \(16c/25\), and \(3c/11\), and the total variations are \(1,c,c\). The respective maximum radial weights are \(1,(5/2)^7,6^7\). This proves the exact weighted-derivative bound used in the checker. The stored pass6 first-derivative lower endpoint is greater than \(-4\) and its upper endpoint is nonpositive, as needed; monotonicity and the endpoint values give these total variations.

The core-swirl speed coefficient is bounded by \(1+8/3=11/3\), giving \(E(W_\psi)\le484\pi/135\). The current outer swirl uses radius \(a=63/200\), and its two axial supports have total length \(9/5\), so
\(E(v)\le\pi(a^4/4)(9/5)=9\pi a^4/20\). This radius is the present datum's radius, not an earlier profile's radius.

For the seed, \(\|w_L\|_2^2=\pi R_0Z_0\), and the two-packet bound is already included in \(Z_0\le12/5\). With \(\lambda=1/4\), its energy is at most \((\pi/32)(392/5)(12/5)\). No additional factor of two is needed. Meridional and swirl components are pointwise orthogonal, the core and outer swirl supports are disjoint, and all seed/mean cross terms vanish after angular integration. Therefore the energy sum used in the checker is valid. With \(\pi<22/7\), it is strictly below

\[
\frac{2015718896182553}{7560000000}<270000.
\]

Combining this with \(\pi>3\) gives exactly the stated rational tail bound \(540000/[L(LN-6)^4]\).

## Full source, stress, and both packets

The radial `M` contribution is integrated over the full sphere, with weight \(2\pi R^2\,dR\,d\mu\). The outer correction integrates only the positive axial packet with weight \(4\pi r\,dr\,dz\); the factor two for the negative packet is valid because the paired kernel and the relevant source/diagonal stresses are even in \(z\). The meridional radial-axial stress has the corresponding odd parity and is already included in the full-sphere part.

The evaluator stores genuine positive Fourier coefficients. Consequently the zero-mode seed covariance is \(2|u_4|^2\), as the image script uses, and the zero-mode source has \(2\operatorname{tr}(G_4\overline{G_4})\). The omitted radial-azimuthal and axial-azimuthal stress components contract with zero entries of the axis-centered kernel. Mean/seed linear terms have nonzero angular frequency and vanish in this axis pressure functional. The outer meridional/azimuthal cross stress also has only those zero-kernel components. Subtracting the radial `M` source in the cylindrical correction therefore avoids double counting without losing a contributing term.

Finally, \(g=\partial_i\partial_j(u_i u_j)\) and compact support prove \(\int g=0\) exactly. The paired source-kernel constant at the source origin is \(1/[\pi(nL)^3]\), so its removal is analytically valid. It is an exact cancellation built into a floating-point quadrature, not a certificate of that quadrature. Source/stress agreement and quadrature refinement remain useful numerical diagnostics only.

The energy paragraph now says: “Angular averaging, pointwise orthogonality of meridional and azimuthal components, and disjoint core/outer swirl supports remove the cross-energy terms.” No mathematical correction to the certificate was found.
