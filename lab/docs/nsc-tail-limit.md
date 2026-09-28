# Infinite depth, finite response, and the missing domain condition

The [tail calculation](../scripts/check_nsc_tail_limit.py) tests the proposed
infinite hierarchy on the existing smooth spatial motif. Its
[record](../results/nsc-5-tail-limit.json) gives an exact endpoint criterion
and a response-error formula, rather than interpreting decimal digits.

## A decimal expansion is not a physical signature

An infinite sequence of contributions can sum to a finite number. It need not
produce a nonterminating decimal: `sum 2^(-n)=1`. Conversely, `1/3` has an
infinite decimal expansion without specifying a nested physical geometry, and
it terminates in base three. Fractal structure must be established through
the operator/geometry and its scale relation, not the notation of its result.

A non-double-counting definition is `delta_n=O[Gamma_n]-O[Gamma_(n-1)]`, where
Gamma_n is the depth-n response and O is a specified physical observable.
If `|delta_n|<=C q^n` with `0<q<1` is derived, the omitted tail after depth N
is bounded by `C q^(N+1)/(1-q)`. Neither positivity of delta_n nor q=1/Omega
is automatic. The initial observable must also be computed; placing 0.95 in
front of an otherwise unspecified series would impose the desired baseline.

Even the rate of convergence needs proof. With Omega=4, b=2 and K=2, the scalar
recursion becomes `Gamma_(n+1)=2-1/Gamma_n`, starting at Gamma_1=2. Exactly,

\[
\Gamma_n=\frac{n+1}{n}\longrightarrow1,
\qquad |\Gamma_n-1|=\frac1n.
\]

The limit exists, but convergence is algebraic despite Omega>1. Thus a factor
1/Omega in the recursion does not alone prove corrections of order Omega^(-n).

## The same smooth motif supplies a sharp criterion

Use the existing a=1, R=2, kappa=1 smooth-motif lattice with m=16 intervals,
and the declared scaling `H_n=Omega^n H`, `B_n=Omega^n B`. Interleave its node
and edge degrees of freedom and remove hopping signs by a unitary gauge.
This gives a scalar Jacobi path with positive coefficients

\[
a_{2mn+2i}=\Omega^n\alpha_i,\quad
a_{2mn+2i+1}=\Omega^n\beta_i,
\quad \alpha_i=1/h-w_i/2,\quad\beta_i=1/h+w_i/2.
\]

The last beta links to the next motif. Define

\[
q=\prod_i\alpha_i/\beta_i\simeq0.04463272338.
\]

At zero energy the independent even- and odd-site endpoint solutions acquire
motif factors `(-1)^m q` and `(-1)^m/(Omega q)`. Their square summability gives

\[
\boxed{\Omega q\le1:\ \text{limit point};\qquad
\Omega q>1:\ \text{limit circle}.}
\]

In the first case the endpoint supplies no extra self-adjoint boundary datum.
In the second, the operator needs a self-adjoint boundary condition at infinite
depth. A domain-selection rule supplied by the full recursive action would
have to determine it. Equality belongs to the first case. The second formal
solution need not satisfy the regular left endpoint; it still determines the
classification at infinity.

This uses the Weyl endpoint alternative and real-energy criterion in
[Teschl, Sections 2.4 and 2.6](https://www.mat.univie.ac.at/~gerald/ftp/book-jac/jacop.pdf).
It applies to the declared unweighted direct-sum lattice, not every possible
Hilbert measure or the full continuum nested-space operator. Its critical
Omega is approximately 22.4051; refinement changes that lattice value.

## Bound the remote boundary's effect

For `G00(z)=<e0,(z-J)^(-1)e0>` and y=Im(z)>0, generate

\[
p_{-1}=0,\ p_0=1,\qquad
p_{j+1}=\frac{zp_j-a_{j-1}p_{j-1}}{a_j}.
\]

All self-adjoint completions of the first M sites lie in a Weyl disk of radius

\[
R_M=\frac{1}{2y\sum_{j=0}^{M-1}|p_j|^2}.
\]

Consequently any two such completions differ in G00 by at most `2 R_M`.
The bound also admits a passive terminal self-energy with Im Sigma<=0. The
runner verifies the disk against finite real terminal potentials, a passive
terminal response, direct matrix inversion and the original block recursion.

For z=0.3+0.4i, the numerical disk diameters behave as follows:

| Omega | 8 motifs | 64 motifs | Endpoint |
|---:|---:|---:|---|
| 1.3 | about 1.77e-20 | about 5.81e-159 | unique endpoint completion |
| 4 | about 2.04e-13 | about 3.17e-97 | unique endpoint completion |
| 50 | about 0.00504 | about 0.00504 | an infinite-depth boundary condition remains |

The last finite Dirichlet sequence can appear numerically settled while other
permitted completions remain separated. A stable-looking list of digits is
therefore not sufficient to establish the unique infinite response.

The radius formula is exact; displayed values use high-precision arithmetic,
not certified interval rounding. They concern the frozen binary64 lattice
coefficients. Spatial-discretization errors, uncertainty in physical inputs,
and the limit toward real spectral poles are separate. No physical constant
is established to 159 decimal places.

## Turning this into an observational prediction

The bound is for one projected Green function. A dark density requires the
common action's metric variation and a state; a dark fraction requires a
cosmological solution. Their sensitivity to the operator, discretization and
parameter errors must be propagated before comparison with experimental
uncertainty. An arbitrary `0.95 + tail` would insert the desired answer.

The [plateau conditions](nsc-plateau-conditions.md) separately test cosmic-time
evolution. Depth convergence and a late-time attracting ratio are distinct
requirements that the same completed model must satisfy.

Reproduce: `python3 scripts/check_nsc_tail_limit.py --check`.
