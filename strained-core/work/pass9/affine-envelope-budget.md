# A finite-time energy pulse for the whole leading envelope

The exact linear rotating-affine model below gives a finite-time gain for the
entire new radial envelope, including its pressure projection, shear, pump and
viscosity. At the carrier-unwinding time it has **more than 15% energy gain**.
It also has a finite positive shear-work budget. Neither statement proves that
the actual cylindrical angular-momentum torque stays positive for that long.

## 1. Exact model and its energy normalization

Let `r*=0.2154755753…`, `h=4/r*`, `Omega=A chi((r*/.315)^2)`,
`gamma=2 Omega`, and `S=gamma h`. Take `c=7/5`, `nu=1/1000` and

```
A0 = [[c,-Omega,0],[-Omega,c,0],[0,0,-2c]],
A(t)=R(Omega t) A0 R(Omega t)^T.
```

The trace-free, symmetric matrix field `U(x,t)=A(t)x` solves unforced NS with a
quadratic pressure, since `A'+A^2` is symmetric. It has infinite energy. Its
gradient matches the compact datum's pump/swirl gradient at the selected
maximum at `t=0`; prescribing its later evolution is a reduced-model choice.

In coordinates rotating at `Omega`, the horizontal transport velocity is
`(cx, cy-gamma x)`. For a horizontal, axially independent perturbation, the
linearized vorticity satisfies

```
zeta_t + (cx,cy-gamma x).grad zeta = -2c zeta + nu Delta zeta.       (1)
```

The rotation of the basis, and thus the Coriolis contribution in the velocity
equation, have been included. Inverting vorticity below includes the full
linear incompressible pressure polarization. An arbitrary radial envelope
does not have the single-wave cancellation of perturbation self-advection:
(1) is an exact **linearized** equation, not the nonlinear compact NS flow.

Set `X=e^-ct x`, `Y=e^-ct y`, and use the unitary Fourier transform in `X`.
Initially the complex streamfunction is

```
psi0(X) exp(ihY),
psi0(X)=e_new(r*+X) exp(-iaX),       a=20,
e_new(r)=chi(((r-7/50)/(1/10))^2).
```

The real part is the physical perturbation. Its amplitude may be multiplied
by any `lambda`; all energies and fluxes below then acquire `lambda^2`.
For `zeta_hat0=(xi^2+h^2) psi_hat0`, define

```
K=xi+St,   Q=K^2+h^2,
I_c(xi,t)=integral_0^t e^(-2cs)[(xi+Ss)^2+h^2] ds,
D=exp(-nu I_c),
Psi(X,t)=(2pi)^(-1/2) integral zeta_hat0 D/Q exp(iKX) dxi.
```

The exact complex horizontal velocity is
`e^-ct (ih Psi,-Psi_X) exp(ihY)`. The energy used here is

```
E(t)=1/2 integral_X average_Y |w|^2 dX
    =e^-2ct/4 integral |zeta_hat0|^2 D^2/Q dxi.                    (2)
```

This is an energy integral in comoving radial labels with an angular-label
average. The physical energy in one expanding horizontal period is
`(2pi/h)e^(2ct) E(t)`. Consequently its identity lacks the `-2cE` term in the
identity for (2). Neither normalization is the compact three-dimensional
fluctuation energy. The parent nonlinear slice's per-unit-axial-length energy
should not be identified with (2).

## 2. Energy and local torque require different information

Parseval gives the comoving integrated covariance

```
Rint(t)=integral_X average_Y(w_x w_y) dX
       =-h e^-2ct/2 integral K |zeta_hat0|^2 D^2/Q^2 dxi,
E'=-2cE+gamma Rint-(nu/2)e^-4ct integral |zeta_hat0|^2 D^2 dxi.    (3)
```

Every initially leading frequency `xi<0` contributes positive stress only
until `xi+St=0`. A nonconstant compact envelope has a broad spectrum, so its
whole stress need not reverse when the carrier alone unwinds.

The local covariance and flat momentum-flux divergence retain interference:

```
R(X,t)=-(h/2)e^-2ct Im(conj(Psi) Psi_X),
-partial_x R=(h/2)e^-3ct Im(conj(Psi) Psi_XX).                    (4)
```

An increase in (2) does not fix the sign of (4). In particular, (4) is not the
cylindrical angular-momentum flux divergence
`-(1/r)partial_r(r^2 R_rtheta)-partial_z(r R_ztheta)`. Its curvature term,
axial endcaps, evolved meridional flow, finite-amplitude interactions, and
central pressure response are absent here. They cannot be restored by
multiplying (4) by the initial radius at positive time.

## 3. A certified finite-time gain for this entire envelope

Take `tau=a/S`. Put `eta=xi+a` and `H=h^2`. The Fourier power of the real
envelope is even in `eta`; translating its center affects only Fourier phase.
Let `Aj=integral |partial_r^j e_new|^2 dr`, for `j=0,1,2`. At zero viscosity
and with the pump factor in (2) removed,

```
E0=[A1+(a^2+H)A0]/4,
E_inv(tau)-E0
  =(a^2/4) integral |e_hat(eta)|^2
       (5eta^2+a^2+H)/(eta^2+H) deta.                            (5)
```

Thus the complete spectrum, including initially trailing components, has
strictly positive gain at this time. Equation (5) follows by pairing `eta`
and `-eta` in `[((eta-a)^2+H)^2/(eta^2+H)]`; no narrow-band approximation is
used.

For this seed the established cutoff bounds and range quadrature give

```
1/10 <= A0 <= 1/5,     A1 <= 160,     A2 < 800000.
```

The last enclosure produced by `verify_affine_envelope_budget.py` is
`480868.5363 < A2 < 641080.3545`. Only the weaker upper bound is needed.
The coarse certified maximum-location and amplitude bounds imply

```
1000/63 <= h <= 1600/63,    Omega >= 450,
tau <= 7/5000,             E0 < 93.
```

With these bounds (5) gives an additive inviscid gain greater than `16`.
To include viscosity, use `I_c<=I_0` and `exp(-2nu I_c)>=1-2nu I_0`. At
`t=tau`,

```
I_0=tau(eta^2-a eta+H+a^2/3).
```

The even part of `[((eta-a)^2+H)^2 I_0/tau]/(eta^2+H)` equals

```
eta^4 +(2H+31a^2/3)eta^2 + H^2+7Ha^2/3+7a^4
       + a^4(a^2-16H)/(3(eta^2+H)).                              (6)
```

The remainder is negative throughout the allowed `H` range. The viscosity
loss in the energy before multiplying by `e^-2c tau` is therefore at most

```
(nu tau/2)[A2+(2H+31a^2/3)A1+(H^2+7Ha^2/3+7a^4)A0]
 <= 2062985492/1406514375 < 3/2.
```

Finally `e^-2c tau>249/250`, so

```
E(tau) >= (249/250)[E0+16-3/2]
        > E0+27927869/1984500
        > (23/20) E0.                                          (7)
```

The checker verifies all algebra by symbolic identities and all displayed
rational inequalities exactly. Its only new numerical enclosure is a
one-dimensional outward-rounded range integral. Equation (7) is a genuine
finite-time result for (1), not a Taylor coefficient or a fitted simulation.
It establishes net gain at `tau`; it does not claim `E'>0` throughout the
interval, or positive local torque there.

## 4. Finite impulse and eventual exhaustion

Write `Z0=(1/4) integral |zeta_hat0|^2 dxi`. From (2),

```
E(t) <= e^-2ct Z0/h^2.                                         (8)
```

More specifically, positive integrated shear work has the finite budget

```
integral_0^infinity max(Rint(t),0) dt
 <= (1/(4gamma)) integral_(xi<0) |zeta_hat0|^2
                  [1/h^2-1/(xi^2+h^2)] dxi
 <= Z0/(gamma h^2).                                            (9)
```

For each `xi<0`, discard the factors `e^-2ct D^2<=1` and integrate the
positive kernel until `K=0`; `dK/dt=gamma h`. Then use
`max(integral f,0)<=integral max(f,0)`. This proves (9) even when the local
envelope has interference and the integrated energy has multiple pulses.
Changing `lambda` changes both available work and initial energy by
`lambda^2`; it does not provide unlimited work.

At `c=0`, completing the square gives

```
I_0=t(xi+St/2)^2+h^2t+S^2t^3/12,
E(t) <= (Z0/h^2) exp(-2nu h^2t-nu S^2t^3/6).                   (10)
```

For `c>0` the exact moments are

```
J0=(1-e^-2ct)/(2c),
J1=[1-e^-2ct(1+2ct)]/(4c^2),
J2=[1-e^-2ct(1+2ct+2c^2t^2)]/(4c^3),
I_c=(xi^2+h^2)J0+2xi S J1+S^2J2
   =J0(xi+S J1/J0)^2+h^2J0+S^2(J2-J1^2/J0).
```

A simple finite-time bound is
`I_c>=e^-2ct[h^2t+S^2t^3/12]`. The cubic viscous exponent in (10) must not
be extrapolated to arbitrarily large time with `c>0`: dilation makes `I_c`
converge for each fixed `xi`. Even then (8) forces the energy (2) to decay.
In the inviscid unpumped model, the Schwartz initial spectrum gives
`E(t)~Z0/(S^2t^2)` and `Rint(t)~-2h Z0/(S^3t^3)` as `t` grows. Splitting
`|xi|<St/2` from the rapidly decreasing tail justifies these expansions.
Thus that exact model ultimately returns energy to the mean; a compact
envelope postpones exhaustion but cannot regenerate its leading spectrum.

## 5. Diagnostic consequence and remaining lemma

The carrier swing time alone is not an obstruction to finite-time energy
gain: (7) defeats that inference without WKB assumptions. Conversely,
energy gain is not evidence for a persistent receiver torque. A next actual
stage must control the **time integral of the spatial flux divergence** at
an inherited receiver together with the central pressure gate, and must
replace spent leading spectral weight through the actual nonlinear flow.
Equations (4) and (9) specify distinct quantities to test.

`affine-envelope-spectrum.py` evaluates (2)--(4) by Fourier quadrature, with
results in the two `affine-envelope-spectrum-*.json` files. Doubling the
domain at fixed grid spacing agrees closely for the listed values, but is
not a rigorous quadrature enclosure. Its flat-force timing is not a timing
prediction for the parent's nonlinear circular slice or for full 3D NS.

No return map, compact three-dimensional persistence interval, or blowup
claim follows from this reduced model.
