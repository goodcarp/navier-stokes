# A return can survive while its amplification disappears

Status: exact auxiliary-model lemma and analytic deductions, not a Navier–Stokes construction and not a novelty claim.

## 1. The phase variables

At a fixed material label, retain only the highest fast-phase derivatives of the axisymmetric equations in the source paper's coordinates
\(y=(z,r^2/2)\), with phase \(s=Np\cdot A(y,t)\) and \(\zeta=\nabla_y(p\cdot A)\).
Write \(\kappa=\zeta_1^2/r^2+\zeta_2^2\) and \(q=\zeta_1^2+r^2\zeta_2^2=r^2\kappa\).
For mean-zero periodic phase fields, let
\[
 W=(\gamma,\eta)^T,\qquad \eta=\frac{\sigma}{N}\partial_s^{-1}\xi.
\]
Here \(\gamma\) is circulation perturbation, \(\xi\) reduced azimuthal vorticity, and the inverse derivative acts only on mean-zero phase functions. With the older background held in the coefficients, the principal system is
\[
 \partial_tW=B(t)W+d_\nu(t)\partial_s^2W,
 \quad B(t)=\begin{pmatrix}0&-d_{\rm old}/(\sigma\kappa)\\ \sigma c\zeta_1&0\end{pmatrix},
 \quad d_\nu=\nu N^2q.
\]
This representation matches the paper's pair when \(\gamma=T F\) and \(\xi=\Omega F'\): then \(\eta=\widehat\Omega F\). It does not assume \(F''=-F\). The actual source profile is smooth periodic and linear near zero; treating it as a single sine would be unjustified.

This is a principal system. Material-label derivatives, varying coefficients across a packet, cutoffs, nonlinear interactions and the different lower-order terms in the two viscous equations are omitted here. Their control is a separate PDE obligation.

Heat evolution also changes the carrier itself: its exact linear cell \(F(s)=s\) and the normalization \(F'(0)=1\) are not preserved. The source's pointwise center and inherited-core identities must therefore be rebuilt; the upper bounds below do not give a lower bound for a growing physical gradient.

## 2. Exact periodic-profile lemma

Let \(B(t)\) be a continuous real 2-by-2 matrix independent of phase, and let \(d_\nu(t)\ge0\) be continuous. Let \(\Phi(t,s)\) be the fundamental matrix of \(R'=B(t)R\), and put \(D(t,s)=\int_s^t d_\nu(\tau)\,d\tau\). For smooth periodic data on the circle of length \(2\pi\),
\[
 W(t)=\Phi(t,s)e^{D(t,s)\partial_\theta^2}W(s).
\]
Proof: matrix multiplication acts only on components; the heat operator acts only on phase. They commute. Differentiation gives the equation and initial data. Equivalently, the nth Fourier coefficient is
\[
 W_n(t)=e^{-n^2D(t,s)}\Phi(t,s)W_n(s).
\]
This does not require matrices \(B(t_1),B(t_2)\) to commute: all their time ordering is retained in \(\Phi\).

For phase-mean-zero data and every integer \(j\ge0\), Parseval gives
\[
 \|\partial_\theta^jW(t)\|_{L^2_\theta}
 \le \|\Phi(t,s)\|e^{-D(t,s)}\|\partial_\theta^jW(s)\|_{L^2_\theta}.
\]
For \(j\ge1\) the same estimate holds without a mean-zero assumption, since differentiation removes the zero mode. The corresponding statement holds in full periodic Sobolev norms for mean-zero data, and in derivative seminorms otherwise, so this is not a single-harmonic artifact.

For the source paper's controlled propagator bound (3.5),
\[
 \|\Phi(t,s)\|\le C_R\rho(t)/\rho(s),
\]
the net upper bound becomes \(C_R[\rho(t)/\rho(s)]e^{-D(t,s)}\). A large Euler gain by itself is therefore insufficient: its logarithm must be compared with integrated phase diffusion.

With an inhomogeneous residual \(h\), the exact extra term is
\[
 \int_s^t \Phi(t,\tau)e^{D(t,\tau)\partial_\theta^2}h(\tau)\,d\tau.
\]
Thus an application cannot simply ignore injection or other corrections; it must bound this term in the norm used to certify a growing gradient.

## 3. The ratio blind spot

For one Fourier harmonic, write
\[
 T'=a(t)\Omega-d(t)T,\qquad
 \Omega'=b(t)T-d(t)\Omega.
\]
Where \(T\ne0\), the ratio \(v=\Omega/T\) satisfies
\[
 v'=b-av^2,
 \qquad (\log|T|)'=av-d.
\]
The common damping vanishes from the first equation, not the second. With \(b=0\) and \(\Omega=0\) after a return, one has \(T'=-dT\), rather than the constant-amplitude Euler hold.

An exact control makes this distinction concrete. On \(0\le t\le\pi/4\), use \(a=1,b=-1,d=2\) and \((T,\Omega)(0)=(1,1)\). Then
\[
 (T,\Omega)=e^{-2t}(\cos t+\sin t,\cos t-\sin t).
\]
The same perfect reset \(\Omega(\pi/4)=0\) occurs with or without damping. Yet the retained T amplitude is \(\sqrt2e^{-\pi/2}<1\), whereas without damping it is \(\sqrt2>1\). The strict inequality follows exactly from \(e^\pi>1+\pi>2\); its decimal value is about 0.294. A test that counts return success without measuring retained amplitude accepts a contracting reset pulse.

With different damping \(d_1,d_2\), the quotient becomes
\[
 v'=b+(d_1-d_2)v-av^2.
\]
Consequently even the quotient invariance must not be applied to lower-order diffusion errors without estimating their difference.

## 4. What would make this useful in a PDE construction

One needs a stage lemma controlling the entire solution, with a lower bound on retained amplification after paying diffusion and all corrections. A candidate condition is a positive, uniform net logarithmic margin
\[
 \int_{I_m} g_m(t)\,dt-\nu N_m^2\int_{I_m}q_m(t)\,dt\ge\lambda_*>0,
\]
plus analogous control on transition and holding windows, and a bound on the residual Duhamel term. Here \(g_m\) must be a proved component/cone growth rate, not an upper bound for a propagator norm. This expression is a design requirement for this principal-mode mechanism, not a proven sufficient condition for Navier–Stokes blowup. The source's return, global error and geometric bounds must also be rebuilt. Smooth forcing is permitted in Clay C/D; the obligation is the smoothness and admissibility of the total force, not the absence of a control.

Source: Alpöge–Buckmaster, *Blowup for the Euler equations with smooth forcing*, equations (3.3)–(3.5), (7.12)–(7.17), and Proposition 12.3, downloaded September 8, 2026 from https://cims.nyu.edu/~tristanb/euler.pdf. The periodic heat calculation and damped-pulse control above are derived here.
