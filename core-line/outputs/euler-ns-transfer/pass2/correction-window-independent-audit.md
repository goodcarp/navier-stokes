# Independent audit: correction-window.md

2026-09-08. Reviewed `work/pass2/correction-window.md` and its verifier, checked the displayed source recursion against the Euler manuscript's printed pp. 24–26, and ran `python3 work/pass2/verify_correction_window.py`. No original file was modified by this auditor.

**Verdict:** the exact exponents, scalar majorant and conditional diagonal argument are sound. They do not establish the assumed full residual estimate (H), uniform normalized inherited backgrounds, or a regenerative viscous stage. One missing spatial-compactness hypothesis in the mean nonsmoothness test was reported to the author and has been added.

## Exact derivative and depth count

With force unit `rho^(-7/2)`, fixed-stage spatial scale `rho`, time scale `rho^(9/4)`, carrier `M=rho^(-1/16)`, and assumed normalized derivative cost `M^(k+2j-J-1)`, the physical exponent is exactly

`(J+1-56-17k-38j)/16`.

For total order `k+j<=K`, its smallest value occurs at `k=0,j=K`. The proposed `J=38K+87` gives exponent at least `2`. The source-style derivative allowance `K+4J+13=153K+361` is also exact. The statement correctly treats `M^(2j)` as a hypothesis and uses fixed stage coordinates; moving/deforming-frame derivatives are not silently omitted.

## Scalar Catalan bound

The source recursion (6.4) has only previously available coefficients on the right. If `a_n<=E C_n B^n`, then the companion coefficients obey `b_n<=3HE C_n B^n`. Both convolution sums are bounded with the Catalan convolution identity, giving the sufficient condition

`4DH[9+6E(1+3H)] <=132DEH^2 <= B=256DEH^2`.

This proves the stated coefficient majorant by induction. Setting actual coefficients above depth `J` to zero only decreases the nonnegative terminal-source recursion. Summing its geometric bound when `4B delta<=1/2` gives the claimed terminal estimate. The result controls the stated scalar recursion only: it is not a reason to omit any of the source's smallness rows, amplitude conditions, means or activation terms.

## Order schedule and hidden quantifiers

The arbitrary-constant diagonal choice works **only because** the assumed constants `A_K,B_K` are finite and uniform in the stage and all admitted normalized inherited states. After taking running maxima and imposing `K<=m`, every fixed order has a finite admission threshold. The resulting tails obey `||D_x^gamma D_t^j f_m||_infty<=2^(-m)` at each fixed derivative order.

Under the additional `log A_K=O(K^2 log K)` growth assumption and comparable bounds for `B_K`, `K_m=floor(m^(1/3))` eventually pays both the force and contraction constants. Indeed its logarithmic cost is `O(m^(2/3)log m)=o(m)`, whereas `log M_m=(m/16)log2+O(1)`. These are asymptotic assertions with finitely many early stages handled separately, as the text states. The numerical examples do not measure actual PDE constants.

The harmonic test `J_m^(2chi)2^(-m/8)->0` is correct for each fixed polynomial overhead. It does not cover every Fourier mode of a general smooth carrier; the note expressly restricts the claim.

## Endpoint extension and the old prefix

The final wording includes the necessary extra hypothesis: each fixed early increment must independently extend smoothly to the accumulation time, or have uniformly bounded derivatives of every order up to that time. In the latter case the next time derivative gives one-sided limits; this is enough to obtain smoothness up to the endpoint. Flat compact temporal support is another sufficient option.

Without that condition, diagonal estimates through increasing orders do not prove a smooth total force. For example, one old increment `f_0(x,t)=chi(x)(T-t)^(K_0+1/2)` and all later increments zero has harmless tails at every admitted order but fails to be `C^(K_0+1)` at `T`. Its derivatives through order `K_0` are bounded. The current prefix/join conditions correctly exclude this example. Smoothness at internal joins must likewise come from full flatness or proved matching, not only the finite order currently budgeted.

## Mean obstruction and global location

The mean exponent `2c-7/2-k` is correct for a fixed nonzero normalized mean profile of size `rho^(2c)`. A nonzero compact smooth spatial profile has nonzero derivatives at arbitrarily high orders. Under no cancellation, time-disjoint packets, and a finite spatial accumulation point, some high derivative then diverges at points approaching `(x_*,T)`.

The first draft omitted the last spatial condition. Packet centers escaping to infinity could have unbounded global derivative norms while the force remained smooth near every finite spatial point through `T`. The revised note explicitly requires bounded centers and accumulating times and explains this distinction. Its separate common-support/weighted-admissibility requirements remain necessary.

The slow-mean sufficient exponent `11/2+(9/4)K` leaves two powers of `rho` before paying `A_K<=rho^(-1)`, hence one summable power afterward. That calculation is consistent.

## Verification

The standard-library verifier completed successfully: exact derivative/depth exponents, the finite scalar recursion and Catalan bounds, terminal remainder, Reynolds-loss thresholds, illustrative order-constant budgets, harmonic overhead and fixed-power mean obstruction all passed. No Lean code, PDE simulation, or actual finite-stage force residual was verified.
