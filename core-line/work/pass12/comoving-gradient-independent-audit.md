# Independent review of the comoving-gradient lemmas

Scope: bounded internal mathematical review of
`comoving-gradient-decay.md`, including its Section 6 extension. This is
not external peer review or a formal PDE verification.

The independent assessor read the positive-time proof and checked:

- The cross-term absorption leaves the cutoff coefficient
  `4n(1-d)+24d`; `C_n=max(4n,24)` is valid.
- The drift terms combine to `2L(1-d^2)`.
- Translating with the actual drift removes dependence on the cylinder
  center without changing the diffusion tensor.
- The global polynomial-tail factor, the nonzero-origin drift variant,
  and the integral `16*pi^2*R^5/105` are correct.
- The classical/approximation premise and the initial-trace qualification
  prevent the original weak-solution regularity overclaim.

The assessor separately read and accepted Section 6 under its explicit
representation and `C_b^1` datum hypotheses, checking the reversed-time
SDE sign, `J'=-Db J`, comparison with the noiseless flow, the `exp(3Lt)`
weighted-square estimate, Doob moment constants, weighted norm-value
initial trace, and `c_(5,10)<10/3`. The review agrees that convergence of
weighted norm values is distinct from convergence of weighted
differences. No defect was found in this bounded review.

The parent independently checked the same Bernstein constants and the
Section 6 Itô/flow/weight calculation. The companion symbolic checker
passes 26 exact finite algebra checks; its output explicitly excludes
existence, general solution regularity, NS bootstrap, and S3 certification.

The note and checker are frozen after this review. Their use remains
conditional on their stated solution, global drift, and initial-data
premises; no source theorem or research goal is declared complete.
