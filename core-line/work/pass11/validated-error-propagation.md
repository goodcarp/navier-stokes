# Rational slab propagation for the sharp error hierarchy

`validated_error_hierarchy.py` now implements a validated comparison-system
propagator for the five-component inequality proved in pass10. It accepts
supplied bounds; it does **not** certify that those bounds describe an NS
approximation. The accompanying examples and checks are manufactured
comparison systems, not data from the selected fluid evolution.

## 1. Inputs and the inequality being propagated

For a fixed verification length \(\ell>0\), the error components are
\(E_j=\ell^j\|\nabla^j(u-v)\|_2\), \(j=0,\ldots,4\), using ordered
Cartesian derivative tensors. On each entire supplied time slab, the caller
must have upper bounds

\[
 \|\operatorname{sym}\nabla v\|_{\infty,op}\le S,\qquad
 \|\nabla^j v\|_{\infty,Frob}\le M_j\quad(2\le j\le5),
 \qquad\ell^j\|\nabla^j\mathcal R\|_2\le\delta_j.
\]

Here \(v\) is the complete smooth, divergence-free whole-space
approximation and \(\mathcal R\) its full projected NS residual. Supremum
bounds from only selected spatial regions or time samples do not meet this
premise. The coefficients are held at their supplied nonnegative upper
bounds on the whole slab.

The comparison inequality is

\[
 D^+E\le BE+cq\|E\|_2^2+\delta,
 \quad q=(0,1,3,7,15)^T,\quad c=\ell^{-5/2},                 \tag{1}
\]

with \(B_{jj}=(j+1)S\),
\(B_{ji}=\binom{j+1}{i}\ell^{j-i}M_{j-i+1}\) for \(i<j\), and zeros
above the diagonal. Ordering derivatives from high to low would instead
make this matrix upper triangular. Initial componentwise error bounds are
required. There is no initial error reset at later slabs.

All inputs are integers or exact fraction/decimal strings. JSON decimal
tokens are parsed as strings before conversion to rational numbers. The
Python API rejects binary floats. An irrational \(c\) is bounded above
using an integer-square-root construction satisfying exactly
\(c_{\rm up}^2\ell^5\ge1\). All following arithmetic uses integers and
`Fraction`, including every acceptance comparison.

## 2. Finite-path linear flow without an off-diagonal exponential

Let \(d_j=B_{jj}\). For each increasing index path
\(i=i_0<i_1<\cdots<i_k=j\), multiply its off-diagonal entries

\[
 w=\prod_{a=1}^k B_{i_a i_{a-1}}.
\]

Its contribution to \((e^{Bh})_{ji}\) is \(w\) times the convolution of
the \(k+1\) scalar functions \(e^{d_{i_a}t}\). Summing all increasing
paths is exact: it follows by recursive variation of constants in the
triangular ODE. There are finitely many paths, with at most four
off-diagonal edges. For the time integral \(J(h)=\int_0^h e^{Bt}dt\),
prepend one zero rate to the same convolution.

For nonnegative rates \(d_0,\ldots,d_k\), that kernel has the positive
series

\[
 K(h)=h^k\sum_{n=0}^\infty
 \frac{h_n(d_0h,\ldots,d_kh)}{(n+k)!},                       \tag{2}
\]

where \(h_n\) is the complete homogeneous polynomial of degree \(n\).
The implementation computes its coefficients from the generating product
\(\prod_a(1-(d_a h)t)^{-1}\) using a finite recurrence. There are no
alternating divided differences and no near-equal-rate subtraction.

If \(x\ge\max_a d_ah\), then
\(h_n(d_0h,\ldots,d_kh)\le\binom{n+k}{k}x^n\). Thus the remainder
after degree \(N\) is bounded by

\[
 \frac{h^k}{k!}\sum_{n>N}\frac{x^n}{n!}.
\]

When \(x/(N+2)<1\), the latter exponential tail is at most

\[
 \frac{x^{N+1}}{(N+1)!}
       \frac1{1-x/(N+2)}.                                  \tag{3}
\]

All quantities in (2)–(3) are rational. The program increases \(N\) until
the tail is below its requested tolerance or reports a series-budget
failure. A common choice \(x=5Sh\) works for every path. Crucially, the
large coefficients \(M_j\) appear only in finite path products; none is
inserted into this exponential argument. If \(S=0\), all path kernels are
exactly \(h^k/k!\), and the linear flow is a finite polynomial even when
an off-diagonal coefficient is enormous.

This produces entrywise rational enclosures for
\(H=e^{Bh}\) and \(J=\int_0^h e^{Bt}dt\).

## 3. Close one nonlinear slab by an exact quadratic inequality

From the linear enclosures form componentwise upper bounds

\[
 a\ge H E_{\rm in}+J\delta,\qquad b\ge c_{\rm up}Jq.
\]

Find a rational \(K\ge0\) such that

\[
 \|a+bK\|_2^2\le K.                                       \tag{4}
\]

This gives a bound for the entire slab, not just its endpoint. Indeed, all
linear solutions involved are componentwise nondecreasing because the
matrix and sources are nonnegative. The endpoint vectors \(a,b\) therefore
dominate their intermediate-time versions. On the box of continuous paths
\(0\le z(t)\le a+bK\), the Volterra map for the comparison ODE satisfies

\[
 e^{Bt}E_{\rm in}
  +\int_0^t e^{B(t-s)}[\delta+c_{\rm up}q\|z(s)\|^2]ds
 \le a+bK.
\]

It preserves that box. Standard Picard iteration converges there by the
bounded local Lipschitz estimate (the Volterra iterates have factorial
convergence bounds; one need not assert a one-step contraction). It produces
a comparison solution on the full slab. The componentwise Dini comparison
then bounds the error by that solution and hence by \(a+bK\).

To solve (4), let

\[
 A=\|a\|^2,\quad P=a\cdot b,\quad C=\|b\|^2,
 \quad t_0=1-2P,\quad D=t_0^2-4AC.
\]

If \(A=0\), choose \(K=0\). Otherwise, sufficient conditions are
\(t_0>0\) and \(D\ge0\), since the entirely rational choice

\[
 K=\frac{2A}{t_0}
\]

satisfies
\(CK^2-t_0K+A=-AD/t_0^2\le0\). If \(C=0\), the linear root is used.
Optional rational bisection between \(A\) and the feasible value above
approaches the smaller root while always retaining a feasible upper
candidate. No floating-point square root decides acceptance.

The implementation rounds \(a,b\) upward to dyadic rationals before
checking (4). It rounds the componentwise endpoint upward again afterward,
so denominator sizes do not grow unchecked across many slabs. The stored
`closure_K` certifies the unrounded \(a+bK\); the separately stored
whole-slab component bound and its squared norm include the final upward
rounding. Confusing those two radii would lose the exact closure check.

## 4. Inheritance, subdivision and failure

The next slab starts with the previous certified componentwise endpoint.
A supplied `jump_error` is added once at an original slab boundary if the
approximation changes representation there. It must bound
\(\ell^j\|\nabla^j(v_{\rm new}-v_{\rm old})\|_2\). Automatic bisection
does not apply that jump again.

If a slab fails the quadratic test or exhausts the exponential-series
budget, the program can bisect it, reuse the same valid coefficient upper
bounds, and propagate the first half's endpoint into the second half.
If the allowed depth or minimum duration is reached, it returns only the
validated prefix and labels the remaining interval unverified. Failure of
this conservative box test is not a blowup result; one manufactured test
deliberately fails a single-slab test while its exact solution remains
finite, then succeeds after subdivision.

The verification length is fixed within a `propagate` call. The helper
`reweight_endpoint` provides the exact conversion
\(E_j^{\rm new}\le(\ell_{\rm new}/\ell_{\rm old})^j E_j^{\rm old}\)
if a later invocation deliberately changes that norm. Such a change does
not erase prior errors.

## 5. Interface and checks

Run from the directory containing the scripts:

```text
python3 validated_error_hierarchy.py input.json --output certificate.json
python3 verify_validated_error_hierarchy.py
```

The input structure is:

```json
{
  "ell": "1/10",
  "initial_error": ["0", "0", "0", "0", "0"],
  "slabs": [
    {
      "duration": "1/1000000",
      "S": "1",
      "M": ["1", "1", "1", "1"],
      "residual": ["1e-20", "1e-20", "1e-20", "1e-20", "1e-20"],
      "jump_error": ["0", "0", "0", "0", "0"]
    }
  ]
}
```

This schema example is synthetic. `M` lists orders two through five, and
`residual` and `initial_error` list the scaled tensor orders zero through
four. Optional controls set rational rounding precision, exponential-series
budget, radius-root bisections, adaptive slab bisections and minimum slab
duration. Every reported numerical bound is an exact rational string.

The checker currently passes 44 checks, including:

- Exact nilpotent linear flow versus an independently computed matrix
  polynomial, including all path multiplicities and integrated forcing.
- A nonzero-diagonal, nontrivial path versus independent elementary
  exponential formulas, enclosing both its homogeneous and forced flow.
- The exact nonlinear manufactured solution \(z=q\alpha(t)\),
  \(\alpha'=284\alpha^2\), and enclosure after inherited subdivision.
- An off-diagonal coefficient \(10^{40}\) with \(S=0\), whose manufactured
  high-order endpoint remains below \(1.1\times10^{-20}\). The
  exponential argument is exactly zero.
- Outward rounding, nonsquare \(\ell\), representation jumps, zero-error
  preservation, invalid inputs (including direct Python constructors) and
  both explicit failure mechanisms.

`manufactured_hierarchy_example.json` reproduces the large off-diagonal
test. None of its values are bounds for the selected NS datum.

To use this for the actual project, the missing inputs are still the
whole-slab bounds for the full approximation and its residual, initial and
representation errors, and the endpoint observable margins. Once those
inputs are certified, this program supplies the comparison propagation;
the pass10 pointwise, energy and inherited H4 conversions then apply.
