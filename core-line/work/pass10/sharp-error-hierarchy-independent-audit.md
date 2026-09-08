# Independent audit of the ordered-tensor error hierarchy

Reviewed `sharp-error-hierarchy.md` and `verify_sharp_error_hierarchy.py`. The checker completed with PASS, 105 exact algebraic/combinatorial checks. I independently derived the linear contractions, nonlinear scaling, and axial Poincare lower bound. No correction is required. This is a conditional analytic error estimate, not a certificate for an evolved approximate field or a useful endpoint.

## Linear tensor contractions

For the unweighted ordered tensor \(T=\nabla^j e\), differentiation of \(v\cdot\nabla e\) gives \(j\) terms with exactly one derivative on \(v\). Each acts as \(\nabla v\), or its transpose, on one derivative index of \(T\). In the full Frobenius pairing, the antisymmetric part is skew and contributes zero. Each symmetric part contributes at most \(S\|T\|_2^2\). The term \((\nabla^j e)\cdot\nabla v\) acts on the velocity-component index and contributes one additional \(S\|T\|_2^2\). Thus the diagonal is exactly bounded by \((j+1)S\); no component-counting factor is needed.

For each fixed derivative subset, the other tensor products are contractions of two Frobenius tensors. Cauchy–Schwarz in the contracted index bounds the output Frobenius norm by the product of the input norms with constant one. If the lower error derivative has order \(h<j\), the two linear families contribute coefficients

\[
\binom jh+\binom j{h-1}=\binom{j+1}h,
\]

where the second term is absent at \(h=0\). After the stated \(\ell^j\) weighting, their coefficient is
\(\binom{j+1}h\ell^{j-h}M_{j-h+1}\). This reproduces every entry of the displayed lower-triangular matrix. In particular the last row is \((\ell^4M_5,5\ell^3M_4,10\ell^2M_3,10\ell M_2,5S)\). Pressure is removable from these energy pairings because every ordered derivative of the error is still divergence free. The dropped viscous contribution has the favorable sign.

## Embedding, nonlinear constant, and norm scaling

The possible subtlety is whether applying scalar H2 embedding to a tensor loses a multiplicity factor. It does not. For a tensor of derivative order \(a\le2\), and a fixed total multiindex \(\gamma\), the summed scalar-H2 coefficient is

\[
\sum_{\substack{|\alpha|=a\\\alpha\le\gamma}}
\frac{a!}{\alpha!}.
\]

The ordered tensor norm at total order \(|\gamma|\) has coefficient \(|\gamma|!/\gamma!\), equal to the same sum with additional factors \((|\gamma|-a)!/(\gamma-\alpha)!\), all at least one. Hence the scalar-H2 sum is bounded by the available ordered norm. Summing vector components introduces no additional constant either.

The scalar derivative-sum H2 embedding constant can independently be bounded below one without a numerical computation. Its Fourier weight is at least one on \(|\xi|\le1\) and at least \(|\xi|^4/2\) on \(|\xi|\ge1\). The squared Fourier Cauchy–Schwarz constant is therefore at most

\[
(2\pi)^{-3}4\pi\left(\frac13+2\right)
=\frac7{6\pi^2}<1.
\]

Rescaling \(x=\ell y\) thus gives the claimed
\(\|\nabla^a e\|_{\infty,F}\le\ell^{-a-3/2}E_\ell\) for \(a\le2\). Every self-commutator has derivative orders \(a,b\) with \(a+b=j+1\), \(a\le2\), \(b\le4\). Its scaled bound is

\[
\ell^j\ell^{-a-3/2}E_\ell\ell^{-b}E_\ell
=\ell^{-5/2}E_\ell^2.
\]

There are binomial weights totaling \(2^j-1\) after cancellation of the top transport term. This proves \(q=(0,1,3,7,15)\); in particular the nonlinear L2 row is zero. The scalar reduction uses \(\|q\|_2=\sqrt{284}<17\), and the linear quadratic form is bounded by the maximum eigenvalue of the symmetric part of the full displayed matrix. No omitted higher-derivative term can be replaced merely by the strain norm.

The norm conversion in the note is also correct: ordered derivative multiplicities through degree four range from one to twelve. These weights define a fully rotation-invariant norm; they are different from, and consistent with, the horizontal-rotation norm used in the separate co-rotating note. At general \(\ell\),
\(\|e\|_{H^4_{\rm der}}\le\max(1,\ell^{-4})E_\ell\). The velocity and gradient pointwise conversion factors are respectively \(\ell^{-3/2}\) and \(\ell^{-5/2}\). If a direct C2 geometry error is needed, the same proof gives \(\ell^{-7/2}E_\ell\).

The cooperative comparison system is valid on nonnegative coordinates. Its squared-radius nonlinearity is coordinatewise nondecreasing there, and all entries of the upper-bound matrix are nonnegative. The stated L2 fluctuation-energy conversion, mean-value conversion, and central coefficient bounds consequently follow. A successor class still needs the complete stated spatial/amplitude rescaling and a margin larger than the converted error.

## Uniform-H5 baseline obstruction

The radial plateau \(r\le63/400\) has exact outer swirl \(Arq_v(z)\). Its two disjoint axial intervals each have length \(d=9/10\), and their plateau lengths sum to \(d\). For every derivative order used, the smooth cutoff and its derivatives vanish at the interval endpoints, so the Dirichlet Poincare inequality can be applied repeatedly to each bump. Angular averaging eliminates the nonzero-mode seed cross term; poloidal components are orthogonal and the compact core swirl is disjoint. Thus retaining only the pure \(z\)-derivative gives precisely

\[
\|u_0\|_{H^j_{\rm der}}^2\ge
\frac{A^2r_0^4\pi^{2j+1}}{2d^{2j-1}},\quad j=4,5.
\]

The resulting rational lower bounds, \(843075135/64>3600^2\) and \(2341875375/16>12000^2\), and exponent comparisons are correct. With only H4 initialization error \(e_0\), the reverse triangle inequality gives \(V_5\ge\|v(0)\|_{H^4}>3600-e_0\), as claimed. The argument applies to a uniform supremum used as the growth coefficient. It does not infer a lower bound on a time-dependent H5 integral from an initial value, and the note explicitly preserves that distinction.

The new hierarchy avoids charging all higher derivatives to a common exponential coefficient, but whether it actually improves the present certificate enough requires enclosures of the full field's \(S,M_2,\ldots,M_5\), all five residual components, and the comparison solution. No favorable duration follows from the algebra alone.
