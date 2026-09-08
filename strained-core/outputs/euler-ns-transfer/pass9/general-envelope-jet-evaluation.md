# A new full-pressure initial-jet enclosure for the leading packet

The general-envelope formula was evaluated using the actual Pass8 leading
profile, its exact maximizing-radius enclosure, a newly integrated radial
Green solution, and a newly bounded whole-space pressure error.
Uniformly for \(1/8\le\lambda\le1/4\) and its named exact retuning,
\[
 \boxed{11<G_t(0,r_*,\pm4)<106,\qquad
 -101000<G_{tt}(0,r_*,\pm4)<-13000.}
\]
For the radially maximizing branch at each fixed axial center,
\[
 \boxed{-101000<
 \frac{d^2}{dt^2}G(t,r(t,\pm4),\pm4)\big|_{t=0}
 <-12000.}
\]
The outward certificate gives the sharper upper endpoints
\(-13501.0611\) and \(-12125.3276\), respectively.
These concern the new leading family. They are not imported from the
old trailing-family deceleration result.

## Inputs and a newly enclosed nonlocal pressure

The exact mean radial critical point remains
\[
 r_*\in[0.2154755753334,\,0.2154755753336].
\]
The smooth cutoff derivatives through order three of the new radial
envelope and through order four of the mean profile are range-enclosed
at that interval.

The source is exactly separable as \(S_{\rm lead}(r)F_z(z)\).
The trial \(P_{\rm lead}(r)F_z(z)\) uses the complete radial Green
formula with the new center \(7/50\), width \(1/10\), and carrier
\(k=-20\). Outward prefix/suffix moment ranges on 2048 requested
radial panels enclose \(P_{\rm lead}\), \(P_{\rm lead}'(r_*)\),
and \(\int r^3|P_{\rm lead}|^2\,dr\).
The exact regular \(r^4\) tail below \(r=.04\) and decaying \(r^{-4}\)
tail above \(r=.24\) are both included.

The unchanged axial profiles allow reuse of the previously certified
\(\int|F_z''|^2\,dz\). The radial norm is newly computed:
the old trailing residual norm and its numerical \(<5\) pressure
bound are not reused.

The new weighted residual satisfies
\[
 \|r(P_{\rm lead}F_z'')_{\rm real}\|_2<8.556120.
\]
Applying the exact harmonic-slab estimate about \(z=\pm4\) gives
\[
 \boxed{|p_{m,r}(r_*,\pm4)-P_{\rm lead}'(r_*)|
       <1.298391<13/10.}
\]
Its computed interval encloses an upper-bound expression; the lower
endpoint does not bound the actual pressure error from below.
The \(\pi\) factor cancels between the physical residual norm and
the slab evaluation inequality, giving the implemented squared bound
\[
 \frac1{16}
 \left(\int r^3|P_{\rm lead}|^2\,dr\right)
 \left(\int|F_z''|^2\,dz\right)
 \frac{35r_*^6}{4(9/20)^9}\frac{10001}{10000}.
\]

The explicit local part of \(\mathcal S\) in the general-envelope
formula is approximately \(-10714.66\). The trial pressure contributes
between \(9495.05\) and \(9545.93\). The full error enclosure leaves
\[
 -1543.287<\mathcal S<-845.058.
\]
This substantial cancellation illustrates why the pressure contribution
must be retained.

## The actual maximizing-branch correction

For the new non-flat envelope,
\[
 G_{tr}=A[-crg''+\nu(g'''-g''/r)]+\lambda^2\mathcal T_{\rm unit}' ,
\]
\[
 \mathcal T_{\rm unit}'=
 \frac{mk}{2}\left[(e^2)''+(e^2)'/r-e^2/r^2\right].
\]
The fixed-height radial branch adds the positive correction
\(-G_{tr}^2/(Ag'')\) to the fixed-circle second derivative.
Its complete interval is included in the second displayed bound.
Relaxing the exact relation between \(A_\lambda\) and \(\lambda^2\)
to the independent intervals \(A\in[900,1020]\),
\(\lambda^2\in[1/64,1/16]\) only enlarges these enclosures.

This result does not give a turning time by dividing the first
derivative by the second. It does not establish the acceleration
of the unrestricted global mean maximum, whose initial axial
plateau is degenerate. The actual second derivative can change
after insertion. The packet's previously certified positive
initial fluctuation-energy gain is fully compatible with the
negative initial acceleration of this mean receiver.

The certificate script records all dependency paths, interval
inputs, both tails, and exact dyadic endpoints. Its final PASS
requires rational comparisons proving \(G_t>11\),
\(G_{tt}<-13000\), the tracked derivative \(<-12000\), and the
new pressure-error bound expression \(<13/10\).
No positive-time NS integration or formal proof-kernel replay is
performed.
