# Angular pairing removes the inverse-square raw source-current term

For the original equal-weight angular pair, the raw changed coincidence
N,beta integrands after the leading phase/density cancellation are
O(e^-3). Consequently the omitted **paired** raw source tail in the
[cutoff bridge](nsc-source-cutoff-bridge.md) is O(Lambda^-2).
This improves its previous valid O(Lambda^-1) estimate. No numeric
coefficient at160/320 is supplied here.

Inverse-cubic is the first **allowed** power. This result does not assert
that its complete coefficient is nonzero: an individual term such as
`s*(f'_s)^2` cannot establish that before all density and interference
terms are combined. No failed optimization or NON-EXISTENCE conclusion
is involved.

## The source-frequency recurrence has an angular involution

Use the same E=s*e expansion, with e>0, as the cutoff bridge. Set
`D_s=Pi_s-Pi_{-s}=s*S3`. Real geometry gives
`V(-ell)=conjugate(V(ell))`, and D_s anticommutes with V. Thus

\[
L_0(-\ell)D_s\overline A=D_s\overline{L_0(\ell)A},\qquad
\boxed{A_j(-\ell)=(-1)^jD_s\overline{A_j(\ell)}}.
\]

The boxed rule follows inductively from the finite minor recurrence and
major transport equation, including their initial data. To check the
initial affine-vacuum convention, write its minor/major ratio v and
t=1/e. The background Riccati equation is

\[
v_\rho=\frac{i}{a}V_{-s,s}
+\frac{2i}{a^2t}v-\frac{i}{a}V_{s,-s}v^2.
\]

It is invariant under
`(ell,t,v)->(-ell,-t,-conjugate(v))`. Its finite formal solution with
zero principal minor coefficient therefore has the same involution.
Choosing the positive real major normalization
`(1+|v|²)^(-1/2)` preserves it. This fixes a mathematical vacuum-column
phase; no original physical source column is renormalized. The owned
affine approximation error and the exponential coherent-source remainder
remain separate and do not alter these finite asymptotic coefficients.

## The complete current coefficients, including interference

Let A=sum_j e^-j A_j and use the actual input-frequency current

\[
J=\operatorname{Im}(A^\dagger A_z)-se A^\dagger A.
\]

The coefficient at e^-p contains **both** the differentiated terms with
j+k=p and the density terms with j+k=p+1. Under the boxed rule, its
angular parity is `(-1)^(p+1)`. The same is true with S3 in the current.
The bounded N insertion A†VA has that parity too, since
`D_s V(-ell) D_s=-conjugate(V(ell))`. Therefore each complete raw N or
beta coefficient at order e^-p has parity `(-1)^(p+1)`.

The changed leading term vanishes because the intrinsic fields agree on
Sigma. The e^-1 current cancels by
`delta n2_s=s*f'_s`, and the bounded order-e^-1 N term is unchanged on
Sigma. The remaining e^-2 coefficient is odd in ell. It cancels in the
actual equal-multiplicity angular pair with the same source-energy grid,
weights and cutoff. Both source signs are retained explicitly. No
symmetry of finite-energy occupations or reflection phases is assumed;
their exponentially small ultraviolet remainder still requires its bound.
The ell=0 radius response remains identically zero with its baseline kept.

## Remainder order and the unresolved number

Use the finite expansion through **M=4** to obtain an envelope H2 error
O(e^-4), hence a differentiated-kernel/current error O(e^-3). M=2 does
not give this bound, and the M=3 current remainder alone would still be
O(e^-2). Combining the sufficient remainder with the angular cancellation
proves the claimed O(e^-3) paired integrand and O(Lambda^-2) tail.

An actual directed coefficient still needs the initial expansion error,
H2 norms of L0 A4, the lower coefficient and observable norms, and the
owned potential-derivative propagation integrals. The same quantities for
the retarded variation are needed for a tangent bound. Smoothness gives
finite constants for each fixed admissible history, not a numerical
constant uniform over the whole unrestricted history class. The required
numbers are not inferred from the frozen-background source certificate.

The [executable check](../scripts/check_nsc_source_tail_angular_parity.py)
verifies the generator/Riccati involutions and complete observable parity
at orders1,2,3. It also checks the original paired source grids, weights
and ledger. It does not evaluate a new physical source or field.

```sh
python scripts/check_nsc_source_tail_angular_parity.py --check
```
