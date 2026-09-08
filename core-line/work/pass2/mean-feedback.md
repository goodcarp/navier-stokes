# The packet mean must be part of the coupled stage

This is an exact algebraic diagnostic for a proposed localized-wave construction. It is not an approximate solution estimate or a regeneration theorem.

## A divergence-free packet has a slow stress of leading size

Let chi(x,y,z) be a fixed smooth envelope, M>0, and phi a global phase shift. Define the actual divergence-free field

    w_M^phi = curl[(chi/M) sin(My+phi) e_z]
            = chi cos(My+phi) e_x
              +(sin(My+phi)/M)(chi_y e_x-chi_x e_y).

Average over phi in [0,2pi]. Writing a=(chi_y,-chi_x,0), the quadratic stress is exactly

    <w_M^phi tensor w_M^phi>_phi
       = (chi^2/2) e_x tensor e_x + (1/(2M^2)) a tensor a.

The phase average commutes with physical differentiation and with the linear Leray projection P. Since each field is divergence free,

    <P[(w_M^phi dot grad)w_M^phi]>_phi
       = P div[(chi^2/2) e_x tensor e_x] + O(M^-2)

in every fixed norm in which the indicated fixed-envelope derivatives and projection are bounded. The exact remainder is P div(a tensor a)/(2M^2); no unspecified carrier term is hidden in this expression. Norm constants depend on the envelope. The displayed leading vector before projection is (partial_x(chi^2/2),0,0). Its curl is

    (0, partial_z partial_x(chi^2/2), -partial_y partial_x(chi^2/2)).

For a generic localized envelope this is nonzero, so pressure cannot remove the leading slow force. It does not acquire an inverse-M factor merely because the carrier oscillates rapidly. By contrast, a spatially homogeneous single plane wave has constant mean stress and no stress divergence. Its leading self-advection cancellation alone therefore supplies no self-generated mean strain.

This phase average is an average of explicitly defined velocity fields over an auxiliary shift parameter. It is not claimed to be a literal pointwise spatial average or an actual NS solution family. If a full family of actual solutions is used, its mean obeys the usual exact Reynolds-stress equation with its actual covariance; the covariance then has to be evolved, not prescribed by this initial ansatz.

## The useful force is also the term that must be accounted for

For physical packet amplitude epsilon U and envelope length L, this leading mean acceleration is of size epsilon^2 U^2/L. Over one turnover L/U it gives a velocity change of scale epsilon^2 U if the stress remains at that scale. Small epsilon alone therefore does not give an order-one relative feedback on that interval. The compact-packet pressure lemma in active-core-pressure.md gives an exact example of the favorable strain sign, without using a phase average.

There are three mathematically different options for a proposed mean term:

1. It cancels after pressure selection, in which case it does not drive the desired mean flow.
2. It is included in the actual slow velocity evolution, so its feedback must be controlled through growth and return.
3. It is left in the external force, which must meet the all-order physical smoothness requirements in correction-window.md.

Increasing the depth of phase-mean-zero corrections does not by itself execute option 2 or 3. A useful design should place the intended strain-producing stress in the coupled velocity dynamics and obtain an all-order estimate for the remaining error. Counting the same term as beneficial evolved feedback and as an omitted small residual would be inconsistent.

## Why simply adding exact homogeneous waves does not close this problem

For two waves b sin(k dot x), c sin(l dot x), with k dot b=l dot c=0, the cross-advection creates frequencies k+l and k-l with respective vector coefficients

    [(b dot l)c+(c dot k)b]/2,
    [(b dot l)c-(c dot k)b]/2.

At each nonzero frequency, the coefficient must be projected onto that frequency's orthogonal plane; a zero-frequency sine term vanishes. In general the projected coefficients are nonzero and the sum of independently evolved waves is not a solution. Degenerate cases exist: with zero affine background gradient, k=e_x, l=e_y and both polarizations parallel to e_z, the two heat-decaying waves coexist exactly and have zero nonlinear interaction. Thus a universal prohibition on every noncollinear two-wave solution would be false. These degenerate examples supply no nonlinear transfer.

The distinction between a spatially valid solution and wave ODEs valid only along a chosen trajectory is also the central issue in Le Dizes–Leblanc's [2006 primary note](https://www.irphe.fr/~ledizes/web/JFM2006a.pdf). We use that source as a warning about the proposed bootstrapping method, not as a blanket classification theorem for all multifrequency flows. Localized, mutually interacting waves require a full PDE argument with the pressure and new frequencies retained.
