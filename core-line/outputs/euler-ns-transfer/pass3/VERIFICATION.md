# Verification and scope of the third pass

Run `python3 checks/run_checks.py` from any directory. SymPy is required; no other external Python dependency is used. CHECK_RESULTS.txt records the packaged run.

1. `verify_quantitative_core_stage.py` checks the explicit Sobolev counting constants, derivative-norm Fourier comparison, time-derivative bound scaling, and finite remainder/gain inequalities. Its sample scalar inputs are model norm values, not numerically evaluated profile integrals.
2. `verify_rotational_receiver.py` checks the local projected acceleration, radial harmonic moment examples, exact receiver scaling and uniform gain/radius arithmetic.
3. `verify_material_core_transfer.py` checks material deformation/area identities, the rigid-rotation/Laplacian commutator and the harmonic initial-jet relations used in the circulation argument.
4. `verify_receiver_flux.py` checks the receiver divergence, nonlinear angular-momentum flux, viscous weight Laplacian, harmonic quadratic coefficient and matched energy exponent.

The analytical arguments retain the full Leray projection/pressure and viscosity. Standard local smooth NS existence is used as a known theorem, with an explicit high-order energy estimate giving the interval used here. An independent agent reviewed the Sobolev constants and derivative ladder. The rotational projection and root-proposed exact scale matching were derived and checked independently. The material-loop and fixed-receiver calculations were compared to ensure their different observables were not conflated. Independent review of the flux/terminal-profile calculations is recorded in the accompanying audit.

These are analytical finite local results with symbolic algebra checks. There was no numerical fluid simulation or formal kernel replay. No reusable terminal profile, inherited-state invariant class, full return map or infinite singular construction has been verified. A successful verifier exit does not certify any of those missing claims.

External primary sources are recorded in SOURCES.json and receiver-source-check.md. Their use is restricted to standard local theory and comparison of known mechanisms; neither selected external cascade/reconnection source proves the missing ordinary-NS transfer.
