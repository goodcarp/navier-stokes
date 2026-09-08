# A certified initial pressure-feedback inequality for ordinary Navier–Stokes

This is a computer-assisted **initial-data and local-time result**, using exact analytical reductions and outward-rounded interval arithmetic. It does not prove a singularity, a returning profile, or repeated amplification. The interval-arithmetic implementation and the displayed analytical identities form its trust base; no formal proof-kernel replay is claimed.

## Explicit compact smooth datum

Let chi(s)=1 for s<=1/4 and chi(s)=0 for s>=1. On the transition define

    chi(s)=1/[1+exp(-1/t+1/(1-t))],  t=(4s-1)/3.

This is the same C-infinity cutoff used in pass4. Let s=|x|^2 and set

    psi(s)=chi(s),
    eta(s)=-(1-chi(4s/25))*chi(1/4+3(s-25)/44),
    Phi(s)=psi(s)+(7/5)eta(s).

For D=diag(-1,-1,2), define the radial straining extension

    S_f(x)=curl[-f(|x|^2)(x cross Dx)/3].

It is smooth, compact and divergence free. Where f is constant it equals f Dx. The pump S_eta is zero for |x|<=5/4, equals -Dx on 5/2<=|x|<=5, and vanishes for |x|>=6. It is deliberately different from the previous thin cylindrical pump.

The core rotation is

    W_psi=[psi(s)+(2s/3)psi'(s)](-y,x,0).

In cylindrical coordinates (r,theta,z), the unit outer swirl is

    v=v_theta e_theta,
    v_theta=r chi(r^2/(63/200)^2)
      *[chi((z-4)^2/(9/20)^2)+chi((z+4)^2/(9/20)^2)].

This defines a smooth Cartesian vector field, including on the axis. Its support is contained in r<=63/200, 71/20<=|z|<=89/20. It is disjoint from the core and lies entirely in the pump's -Dx plateau.

Use the decaying whole-space pressure convention p=N*tr((grad u)^2), N=1/(4pi|x|). Put

    C_v=-p_zz[v](0)>0,    H=204/35,    A=sqrt(H/C_v),
    u_0=S_Phi+W_psi+A v,    nu=1/1000.

The definition of A uses the actual integral C_v, not a rounded numerical approximation. The interval enclosure establishes its strict positivity. This is one fixed compact smooth datum for unforced ordinary Navier–Stokes on R3.

## Certified statement

Let u be its local smooth NS solution and let p be its pressure. At the origin define

    b(t)=partial_z u_z(t,0)/2,
    Omega(t)=omega_z(t,0)/2,
    beta(t)=b(t)/Omega(t).

Initially b=Omega=1 and

    p_zz(0,0)=-8,   b'(0)=Omega'(0)=2,   beta'(0)=0.

The full initial derivative satisfies the certified upper bound

    T := partial_t p_zz(0,0)+32
       <= -290649/8750 < -33.

Consequently

    beta''(0) >= 290649/17500 > 16.

The initial outer swirl-pressure ratio also has positive derivative, with logarithmic derivative strictly greater than 4. These are statements about the actual initial acceleration and pressure of the coupled NS datum. No pressure control is prescribed independently of the velocity.

## Exact reconstruction of the full derivative

The radial meridional part m=S_Phi absorbs every core/pump meridional interaction. Its cubic pressure derivative is the functional C300(Phi) in [the radial integral audit](streaming-radial-audit.md), derived in [pass4](../pass4/core-pressure-derivative.md). The radial core-swirl coefficient is tWb(psi). The annular response to the radial core swirl is exactly zero.

The two swirl supports are disjoint, so no core/outer swirl cross term survives. Axisymmetric poloidal/azimuthal pressure cross terms vanish identically. The exact meridional and radial core-swirl viscosity coefficients are zero. The remaining viscous term is nu A^2 d_v, where

    d_v=2 integral Q_ij partial_k v_i partial_k v_j,
    Q_ij=-partial_zzij N.

Write k_core and k_pump for the exact outer-swirl kernels normalized by Q_theta. The complete identity is

    T=32+C300(Phi)+tWb(psi)
      +H * average_[Q_theta v_theta^2](k_core+(7/5)k_pump)
      +(1/1000)H d_v/C_v.

The average is with respect to the positive measure Q_theta v_theta^2 dx, normalized by C_v. The [local kernel calculation](local-swirl-feedback-kernel.md), [adjoint pressure calculation](swirl-adjoint-kernel.md), and [complete reconstruction audit](full-gate-certificate-audit.md) account for the meridional pressure response as well as swirl transport. The full pressure functional has an exact angular cutoff at degree four for this radial meridional class; no unevaluated high-angular remainder is discarded.

## Bounds used in the certificate

On the actual outer support, the analytical kernel estimates give

    k_core <= -12/5,    k_pump <= -31/2.

These bounds require no pressure quadrature. The pump bound uses only its specified join radii, -1<=eta<=0 and the packet support. The cutoff moment used in k_core is exactly integral psi(s) ds=5/8.

Outward-rounded interval integration of the unchanged smooth cutoff gives

| Quantity | Enclosure used |
|---|---|
| C300(Phi) | [36.0778592591, 65.9877263169], hence <66 |
| tWb(psi) | [2.5300883924, 2.5970809811], hence <4 |
| C_v | [0.0000056503282502, 0.0000069746389609], strictly positive |
| d_v | [0.0020053494699, 0.0050869205474] |
| d_v/C_v | <901 |

The displayed decimal intervals are rounded outward for readability. The certificate uses the exact dyadic endpoints in FULL_INITIAL_GATE_CERTIFICATE.json. For C300, the quadrature carries intervals for both cumulative integrals Q and K across every radial panel; it does not replace them by independently sampled values.

Combining the coarse rational bounds gives exactly

    T <= 32+66+4+(204/35)[-12/5-(7/5)(31/2)+901/1000]
       = -290649/8750.

The swirl-ratio estimate is

    (log(C_theta/Omega^2))'(0)
      >= (2381/357)(7/5)-4-901/1000 >4.

The flat initial core has vanishing viscous curvature terms and their first time derivatives. Differentiating the actual central equations therefore gives beta''(0)=-T/2, as established independently in pass4. This uses the complete pressure derivative, not a differentiated radial-profile ansatz.

## Viscosity one and local consequences

Set U(t,x)=1000 u(1000t,x) and P(t,x)=1000000 p(1000t,x). Then U solves ordinary unforced NS with viscosity one. Its central b,Omega and initial velocity amplitudes are multiplied by 1000, its pressure-feedback gate by 1000^3, and the favorable sign is preserved.

Local smooth existence and continuity turn the strict derivative margins into a positive, possibly very short time interval of improving observables. [The finite-time corollary](finite-feedback-corollary.md) states precisely what follows. No duration comparable to a blowup time is obtained here. The radial profile is imposed only at the initial time and is not assumed to remain radial during evolution.

Reproducibility: run `python3 certificate/certify_full_initial_gate.py --recompute` from the packaged pass directory. This recomputes all four interval integrals used above and combines their exact endpoints with the analytically proved rational kernel bounds. The symbolic checkers and verification note cover the separate identities.

The remaining main problem is same-solution regeneration: retaining the favorable feedback and supplied swirl through a usable normalized terminal state, then iterating compatible stages. The present certificate settles an initial compatibility test only.

[A separate robustness corollary](robust-initial-data.md) uses the same coefficient bounds to obtain strictly favorable first and second initial derivatives with the direct rational amplitude A=1100. The neutral datum above remains the subject of the displayed second-derivative theorem.
