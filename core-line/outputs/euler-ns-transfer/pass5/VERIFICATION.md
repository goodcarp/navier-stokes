# Verification and limitations

The certificate combines exact analytic reductions with inclusion-preserving range quadrature using mpmath.iv at 35 decimal digits. Rational panel partitions, endpoint cutoff estimates, cumulative-integral interval propagation, and exact dyadic endpoint comparison are documented in the mathematical notes. Approximate priorities choose which panels to subdivide; they do not determine the enclosure error.

`python3 certificate/run_checks.py` runs the finite symbolic identities, support/kernel bounds, cutoff consistency tests, central-jet and receiver arithmetic, and exact certificate aggregation. CHECK_RESULTS.txt records the packaged run. These checks supplement the analytic proofs and audits; they do not mechanically prove those entire arguments.

`python3 certificate/certify_full_initial_gate.py --recompute` independently regenerates the four integral enclosures used by the assembler: the combined radial meridional cubic, core strain/swirl coefficient, positive outer swirl stress and its viscous coefficient. FULL_INITIAL_GATE_CERTIFICATE.json records the recomputation and all exact endpoints. Omitting --recompute only aggregates stored component bounds.

The analytical dependencies include the earlier radial pressure/cubic identity, extended to the present sign-changing radial profile; the degree-four cutoff for this scalar functional; the core and annular adjoint kernels; support and viscosity cancellations; and local smooth NS theory. Separate audits check the full decomposition, cumulative quadrature, and finite-time deductions. mpmath's interval elementary functions and Python's exact integer/rational operations remain trusted software. No Lean or other formal proof-kernel replay has been performed.

The stated inequality concerns one actual initial-time pressure derivative. Local existence and strict margins also give a positive interval of the specified improving observables. There was no NS time-stepping simulation. No numerical lower bound for a useful stage duration, invariant profile class, inherited-state return, indefinite amplification, assembled all-order smooth forcing, or singular solution is certified.

The exploratory cylindrical and coarse radial pressure evaluations are not inputs to this certificate. Earlier pass2, pass3 and pass4 packages are unchanged. File hashes identify the delivered artifacts and do not enlarge the mathematical claim.
