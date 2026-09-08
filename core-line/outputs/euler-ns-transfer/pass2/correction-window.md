# Correction depth versus the viscous frequency ceiling

**Status: a conditional compatibility lemma and two quantitative rejection tests, not a Navier–Stokes construction.** The proposed core scales do leave an algebraic all-order correction window, provided a new viscous stage supplies uniform normalized estimates for its **entire** force residual. The Euler paper's terminal correction estimate alone does not supply this hypothesis. Its activation and phase-independent mean terms are separate obligations.

## 1. What can actually be imported from the paper

References below are to the printed page and equation numbers of [Alpöge–Buckmaster, *Blowup for the Euler equations with smooth forcing*](https://cims.nyu.edu/~tristanb/euler.pdf), locally `work/sources/euler.pdf`.

* **Finite inverse-frequency corrections:** Proposition 6.1, pp. 23–25, (6.2)–(6.6), gives a level-\(n\) contribution bounded by \(a_nY_*\delta^n\), after checking every quotient in (6.3). At output order \(k\) and depth \(J\), the input order allocation is \(r_*=k+4J+9\), \(r_F=r_*+4\). The constants are computed by (6.4), not assumed uniformly bounded in depth.
* **Small amplitude is an independent condition:** the fast circulation-square row of (6.3), p. 24, explicitly requires \(T\sigma\varpi_*z_0Y_*/\delta\le1\). The paper says this is an amplitude condition, not an inverse-frequency gain. Consequently even the schematic substitution \(\delta=M^{-1}\) below needs a new argument for a new viscous stage.
* **Three force channels:** Proposition 6.2, p. 27, (6.9), bounds activation, the terminal sum \(S_{\rm rem}=\sum_{g=J+1}^{2J+4}s_g\delta^g\), and phase-independent means separately. The means in (6.7), p. 26, have quadratic amplitudes and slow derivatives. Increasing terminal depth does not multiply every mean by \(\delta^{J+1}\). The separate exponents in (12.15), p. 104, and the all-order summability conclusion (12.17), confirm this distinction.
* **Inherited constants are not free:** (10.26)–(10.28), p. 75, control a recursive majorant with multiplier \(\mu_j=1024(J_j+1)\). The displayed bound for \(\log G\) contains \(\prod_{j=2}^{m-1}\mu_j\); the published super-separated \(\log N_m\) pays this product. That estimate cannot simply be inserted into a geometric core schedule with \(\log(1/\rho_m)=O(m)\). Failure of this upper majorant to fit is not a lower bound on actual constants or a general impossibility theorem.

The paper's fixed off-axis inverse and Euler amplitude equations also do not automatically provide the corresponding estimates for a shrinking three-dimensional viscous core.

## 2. Dimensional bookkeeping for the proposed core

Use the provisional exponents in `outputs/euler-ns-transfer/NEXT_TARGET.md`:
\[
 L=L_0\rho,\qquad \mathrm{Re}=\mathrm{Re}_0\rho^{-1/4},\qquad
 U=\frac{\nu\mathrm{Re}_0}{L_0}\rho^{-5/4},\qquad
 M=M_0\rho^{-1/16}.
\]
Thus the turnover time and force unit are
\[
 \tau=L/U=\tau_0\rho^{9/4},\qquad U^2/L=F_0\rho^{-7/2},
\]
while the carrier diffusion per turnover is
\[
 \frac{M^2}{\mathrm{Re}}=\frac{M_0^2}{\mathrm{Re}_0}\rho^{1/8}\longrightarrow0.
\]

Take each stage's \(\rho_m\) to be a fixed scale for its coordinate change \(y=(x-x_m)/L_m\), \(s=(t-t_m)/\tau_m\). For a vector residual **after the pressure has been selected**, write
\[
 f_m(x,t)=\frac{U_m^2}{L_m}\mathcal R_m(y,s).
\]
This identity must include all perturbation/base interactions and joins. A moving or deforming coordinate frame adds chain-rule terms; controlling those terms is part of the hypothesis, not a consequence of this formula. A continuously changing scale cannot be substituted without paying its dilation and time-reparametrization defects.

Here is one conservative, explicit candidate estimate for a new finite-stage lemma:
\[
 \tag{H}
 \|\partial_y^\gamma\partial_s^j\mathcal R_m\|_\infty
 \le A_{J,K}M_m^{k+2j-(J+1)},\qquad k=|\gamma|,\quad k+j\le K.
\]
The factor \(M^{2j}\) budgets up to two carrier derivatives per normalized time derivative. It is an assumed estimate, not a universal rule for the nonlinear evolution. Fixed units, combinatorial factors, and polynomial-in-depth harmonic costs can be included in \(A_{J,K}\). Crucially, \(A_{J,K}\) must be uniform over stages and the normalized inherited states admitted by the proposed stage lemma.

The exact power of \(\rho\) in the physical mixed derivative bound is then
\[
 \tag{1}
 \|\partial_x^\gamma\partial_t^j f_m\|_\infty
 \le \widetilde A_{J,K}
 \rho_m^{\{J+1-56-17k-38j\}/16}.
\]
We absorb \(F_0L_0^{-k}\tau_0^{-j}M_0^{k+2j-J-1}\) into the order constants. The worst case at total order \(K\) is \(k=0,j=K\). In particular,
\[
 \tag{2} J(K)=38K+87
 \quad\Longrightarrow\quad
 \|\partial_x^\gamma\partial_t^j f_m\|_\infty
 \le \widetilde A_K\rho_m^2\quad(k+j\le K).
\]
This deliberately reserves one power of \(\rho\) to pay order constants and another for summability. At this depth, the paper's analogous derivative allocation would require
\[
 \tag{3} r_F=K+4J(K)+13=153K+361.
\]
This is a significant input-regularity obligation; it is not evidence that such normalized jets are available.

## 3. The finite correction recursion permits an explicit order majorant

For comparison with (6.4), fix the working derivative order and let \(E\ge1\) dominate \(a_0\), \(H\ge1\) its Leibniz constant, and \(D\ge1\) its inhomogeneous solution constant. The source recursion is
\[
 b_n=a_n+H(a_{n-1}+a_{n-2}),\qquad
 s_n=H\left[6a_{n-1}+3a_{n-2}
 +3\!\sum_{i+j=n-1}a_i(a_j+b_j)
 +3\!\sum_{i+j=n-2}a_i(a_j+b_j)\right],
 \quad a_n=4Ds_n.
\]
Zero indices below zero have value zero. Set \(B=256DEH^2\), and let \(C_n\) be the Catalan numbers. Induction using \(\sum_{i+j=n-1}C_iC_j=C_n\) gives
\[
 \tag{4} a_n\le EC_nB^n\le E(4B)^n.
\]
Indeed \(b_n\le3HEC_nB^n\); the two convolutions and linear terms are bounded by
\(H[9E+6E^2(1+3H)]C_nB^{n-1}\).
The needed inequality is \(4DH[9+6E(1+3H)]\le132DEH^2\le B\).
The same upper sequence bounds the terminal sources when actual coefficients above \(J\) are set to zero. Therefore
\[
 \tag{5}
 S_{\rm rem}\le 2E(4B\delta)^{J+1}
 \qquad\text{if }4B\delta\le\tfrac12.
\]
This is a majorant of the **displayed scalar recursion**, not a new viscous PDE estimate. It exhibits both the terminal power and its depth-dependent constant. It also explains why one should test the effective parameter \(4B/M\), not merely \(1/M\).

As an additional quantitative hypothesis, suppose normalized input jets and propagator bounds give
\[
 \log E_r+\log H_r+\log D_r\le C(r+1)\log(r+2),
\]
and the other evaluation constants obey comparable bounds. This is compatible with controlled fixed Gevrey profiles and order-uniform normalized coefficient bounds; it is not established for the candidate evolved backgrounds. With \(r=O(J+K)\), (5) leads to the conservative bound
\[
 \tag{6} \log\widetilde A_K\le C'(K+1)^2\log(K+2).
\]
The factor quadratic in \(K\) retains the repeated correction constants. A claim of merely factorial final constants would need a sharper proof.

## 4. Conditional all-order compatibility lemma

Let \(\rho_m=2^{-m}\) and suppose (H) holds for each admitted finite order, with finite stage-independent constants \(\widetilde A_K\) and \(B_K\). Replace these constants by their running maxima over orders up to \(K\). Then there is a sequence \(K_m\to\infty\) for which
\[
 \tag{7}
 \log\widetilde A_{K_m}\le m\log2,
 \qquad M_m\ge8B_{K_m},
 \qquad J_m=38K_m+87.
\]
For example, admit order \(K\) only once both finite thresholds hold, and additionally require \(K\le m\). Every fixed order is eventually admitted. It follows from (2) that
\[
 \sum_m\|\partial_x^\gamma\partial_t^j f_m\|_\infty<\infty
 \quad\text{for every fixed }\gamma,j,
\]
provided the finitely many increments preceding admission have their own finite derivative bounds.

For a smooth extension at a common accumulation time, each fixed old increment must also be known to extend smoothly through that time, or to have bounded derivatives of all orders up to it (the next time derivative then gives the one-sided limits). Smooth flat joins and compact temporal support strictly before the accumulation time are one sufficient alternative. The diagonal tail estimate alone says nothing about uncontrolled high derivatives of a finite old prefix. Under these prefix and joining hypotheses, uniform derivative summability gives a smooth total force. Compact spatial support in a common set, or appropriate weighted spatial estimates, must be checked separately for the desired global admissibility statement.

Under (6) and comparable bounds for \(B_K\), the explicit schedule
\[
 \tag{8} K_m=\lfloor m^{1/3}\rfloor,\qquad J_m=38K_m+87
\]
satisfies (7) for all sufficiently large \(m\), since \(m^{2/3}\log m=o(m)\). Any finitely many early stages require separate admission or a larger starting scale index. This is an asymptotic compatibility example, not a numerically practical construction or a uniform bound proved by the paper.

Even if active harmonics relevant to amplification reach \(n\le C J_m^\chi\), for a fixed \(\chi\), their largest diffusion-to-strain ratio obeys
\[
 \tag{9} \frac{(nM_m)^2}{\mathrm{Re}_m}
 \lesssim J_m^{2\chi}2^{-m/8}\longrightarrow0.
\]
Polynomial depth-dependent carrier costs therefore leave room in this model. Smooth profiles have infinitely many Fourier modes; (9) is a conditional finite-band test, not an assertion that all modes of an arbitrary smooth profile must grow or evade diffusion.

## 5. A Reynolds-number loss can close the window

If a correction inversion actually costs
\[
 \delta_{\rm eff}\sim C\,\mathrm{Re}^{p}M^{-s},\qquad s>0,
\]
then the proposed \(M\asymp\mathrm{Re}^{1/4}\) gives a vanishing correction parameter only when
\[
 \tag{10} p<s/4.
\]
Equality is a constant-dependent boundary; a power estimate with \(p>s/4\) does not certify contraction at the proposed scales. More generally, the correction floor \(M\gg\mathrm{Re}^{p/s}\) and viscous ceiling \(M\ll\mathrm{Re}^{1/2}\) have a polynomial gap only if
\[
 \tag{11} p/s<1/2.
\]
At equality the required constants and margins decide; beyond equality this particular perturbative frequency-window test fails. It does not prove that every possible construction fails.

For the fixed exponents, put \(\theta=s/16-p/4\). Replacing \(M^{-(J+1)}\) in (H) by \(\delta_{\rm eff}^{J+1}\), while retaining the same derivative costs, changes the raw exponent to
\[
 \theta(J+1)-\frac72-\frac{17k}{16}-\frac{19j}{8}.
\]
Thus \(\theta>0\) still permits linear depth in the admitted order, with a larger coefficient as the gap shrinks. If one treats a residual viscous term perturbatively, \(M^{-1}+M^2/\mathrm{Re}\asymp\rho^{1/16}\) would algebraically retain the original gap; actually cancelling that viscous residual and controlling the resulting inverses remain PDE obligations.

## 6. An unsuppressed mean is a separate obstruction

The previous lemma needs the **full** residual bound (H). Source (6.9) does not obtain this from terminal depth: it separately pays activation with a small seed and means with quadratic amplitude and slow rates. The source's published scale system successfully pays those channels in (12.15); this success cannot be assumed at the proposed shrinking-core scales.

For a concrete conditional falsifier, suppose the uncancelled, pressure-projected mean of a packet has a fixed nonzero normalized profile and normalized size \(\epsilon_m^2\), with \(\epsilon_m=\rho_m^c\) for one fixed finite \(c\). If it is assigned to the force, its physical magnitude is of order \(\rho_m^{2c-7/2}\). For a spatial derivative of order \(k\) whose profile derivative is nonzero, the magnitude at a corresponding point is
\[
 \tag{12} \rho_m^{2c-7/2-k}.
\]
Fixed nonzero compact smooth profiles have nonzero derivatives of arbitrarily high order. Some such \(k>2c-7/2\) therefore gives an unbounded derivative along the shrinking cores. Suppose also that the packet times accumulate at the proposed terminal time and their centers stay in a fixed compact set. With time-disjoint packets, no cancellation by the other packets, and no competing same-stage cancellation, a convergent subsequence of those centers then gives a finite spacetime accumulation point where the force cannot extend smoothly. Centers escaping to spatial infinity would not imply this local nonsmoothness conclusion. An analogous temporal derivative scales as \(\rho_m^{2c-7/2-9j/4}\). Without a nonzero profile/lower bound or without the no-cancellation hypothesis, these powers only show that this upper-bound method fails to certify smoothness; they are not a nonsmoothness theorem.

One sufficient tail requirement for a slow normalized mean at total order \(K\) is a bound of size \(A_K\rho^{11/2+(9/4)K}\) on its normalized mixed derivatives, with \(A_K\le\rho^{-1}\); the physical derivatives are then at most \(\rho\). This is faster than every fixed power once \(K\to\infty\). Alternatively, means may cancel exactly after pressure selection or be incorporated into the evolved background and corrected there. For a localized oscillatory packet, the slow Reynolds-stress divergence can drive a useful mean flow; it is not legitimate to count that effect as evolved strain and simultaneously dispose of it as a small external force.

Small activation data can pay an activation residual. It does not, by itself, pay an order-one mean left after amplification. Likewise, adding more oscillatory correction levels does not prove that a low-grade phase-independent mean disappears.

## 7. Decision for the next stage lemma

The algebra does **not** reject polynomially growing Reynolds number plus \(M=\rho^{-1/16}\): linear-in-order depth with a sufficiently slow all-order schedule can fit below the diffusion ceiling. To use that window, a candidate stage must provide (i) a contracting correction inverse, including its Reynolds-number losses; (ii) stage-uniform normalized jet constants through the allocated orders; and (iii) a full residual estimate covering means, activation, joins, and inherited interactions. Failure of a proposed stage on any of these explicit tests is a reason to reject that stage design before building an infinite cascade.

`verify_correction_window.py` checks the exact exponents, the scalar Catalan majorant, representative finite terminal sums, and the asymptotic polynomial/exponential comparison. It does not test existence, stability, regeneration, force cancellation, or any imported PDE code.
