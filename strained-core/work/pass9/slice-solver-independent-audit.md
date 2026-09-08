# Independent audit of the nonlinear circular-slice diagnostic

Reviewed `evolve_circular_slice.py` and `check_slice_solver.py` against

```
zeta_t + div_h[(c x+v) zeta] = nu Delta_h zeta,
curl_h v=zeta,     div_h v=0.
```

I found no sign, factor-of-two, or Fourier-dealiasing error in the implemented
interior discretization. The result is a finite-domain numerical diagnostic,
not a certified solution of the compact three-dimensional problem.

## Equation, pressure and Fourier normalization

* The `m=4j` array correctly represents the fourfold angular symmetry of the
  initial datum. The forward-normalized real FFT treats positive modes as
  Fourier coefficients, so the initial `q[:,1]` is one half of the physical
  complex amplitude. The mean needs no such factor.
* `psi=(-Delta_m)^-1 zeta`, `v_r=im psi/r`, `v_theta=-psi_r` have the correct
  vorticity sign. Recovering these velocities at every stage retains planar
  incompressible pressure and the generated mean; it is not a passive scalar
  model. Axial pressure and the exterior meridional field remain absent.
* Keeping exactly `j<nphi/3` de-aliases every quadratic interaction affecting
  the retained modes. The strict inequality matters; it is used consistently
  in all linear solves, velocities and output projections.
* The finite-volume pump matrix represents `-r^-1 partial_r(c r^2 zeta)`.
  It contains both the radial transport and the `-2c zeta` stretching term.
  There is no omitted stretching source for this prescribed affine slice.
* The nonlinear radial and angular fluxes use mutually compatible face and
  cell velocities. Their discrete divergence cancels. The two Heun
  evaluations with the same Crank--Nicolson linear solve give a second-order
  additive time step; they do not solve a fully implicit nonlinear equation.

## Energy and circulation

With the script's Fourier coefficients, the continuum quantities are

```
K_fluct = 2pi integral_0^R r sum_(j>0)(|v_r,j|^2+|v_theta,j|^2) dr,
K_mean  = pi integral_0^R r |v_theta,0|^2 dr,
Z = 2pi integral_0^R r(|zeta_0|^2+2 sum_(j>0)|zeta_j|^2) dr.
```

The implemented factors agree. On the whole plane, with sufficiently
decaying velocity, the per-unit-axial-length horizontal energy obeys
`K'=-nu Z`: the pump work `-c integral|v|^2` is canceled by the area dilation
in the energy integral. Multiplication by `e^-2ct` gives a distinct weighted
observable and must be labeled accordingly. It is not the energy law for
this unweighted slice or for the compact three-dimensional field.

The mean swirl recovered from
`r v_theta,0(r)=integral_0^r rho zeta_0(rho) drho` has the right normalization.
The conservative nonlinear term preserves that integral, up to rounding.

## Boundaries and diagnostic limits

1. The pressure/Biot--Savart boundary for each nonzero mode is the centered
   approximation of the exterior harmonic condition
   `psi_r(R)=-m psi(R)/R`. Its ghost coefficient
   `(1-mh/(2R))/(1+mh/(2R))` and boundary-face interpolation are consistent.
   The vorticity boundary is instead `zeta(R)=0`. Together these define a
   truncated model; no exact whole-plane boundary-error enclosure is present.
2. The reported energy sums stop at `R`. For the harmonic nonzero-mode
   continuation, the exterior instantaneous energy is exactly
   `2pi sum_(m>0) m |psi_m(R)|^2`. Finite-radius energy balance also has
   boundary terms. Consequently `Delta K+nu integral Z` is a diagnostic
   residual, not an exact finite-domain conservation law or error bound.
3. A nonzero residual mean circulation produces a `1/r` mean velocity outside
   the grid and logarithmically divergent whole-plane energy. The intended
   compact initial swirl has zero total circulation; numerical residual
   circulation is a representation error to remove or enclose in any future
   whole-space reconstruction. It must not be silently treated as a harmless
   exact whole-plane field.
4. The spline/bracket solve tracks a critical point inside `(.18,.28)`.
   Without an all-radius sign/maximum check, `mean_maximum` is a tracked
   critical branch, not a proved global radial or spacetime maximum. Likewise,
   the reported Reynolds-stress torque excludes pump and viscosity in the
   total mean-angular-momentum time derivative.
5. Central radial advection and Fourier projection do not supply a positivity
   or nonlinear stability certificate. Resolution/time-step comparisons are
   necessary numerical controls; the analytic tests alone do not establish
   convergence of this particular nonlinear run.

The parent has acknowledged the boundary, energy and tracked-branch scope
and is labeling the numerical report accordingly. These qualifications do
not reverse the observed finite pulse, but delimit what the observation
establishes.

## Tests reproduced

Ran `python3 work/pass9/check_slice_solver.py` independently. It passed:
the known `r^4 exp(-(r/a)^2)` Poisson profile and Gaussian affine-diffusion
evolution converge at second order over 192, 384 and 768 radial cells;
discrete divergence is below `1.8e-15`, and nonlinear circulation error below
`1.8e-14` in those runs. These test the selected identities and consistency,
not a continuum enclosure or three-dimensional evolution.
