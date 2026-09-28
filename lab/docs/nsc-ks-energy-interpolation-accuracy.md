# Interpolation component of the retained-source matter error

This record combines the certified geometric profile bounds, the reviewed
row-norm energy-derivative estimate and the nodal polynomial of each actual
stored operator table. It propagates this one error component through all
sixty retained families and their separate negative-energy partners. It runs
no field or source evolution.

For the existing one-amplitude control, write `delta r = alpha*d`,
`d=chi(s)*s*w_direction(z)`. The radius record bounds the physical profile
`W=|alpha*w_direction|` and its first axial derivative `W_z`, with `U=0`.
Let `Delta rho=rho_up-1`, `sigma` be the fixed normal support, and let
`a_min,r_min` be the record's exact rational lower bounds. Then

\[
D\le\Delta\rho/a_{\min}^2,\qquad
K_0\le {\Delta\rho\over a_{\min}}(|m|+|\ell|/r_{\min}),
\]
\[
K_1\le {\Delta\rho|\ell|\sigma W_z\over a_{\min}r_{\min}^2},\qquad
J_0\le {\Delta\rho|\ell|\sigma W\over a_{\min}|\alpha|r_{\min}^2},
\]
\[
J_1\le {\Delta\rho|\ell|\sigma W_z\over a_{\min}|\alpha|r_{\min}^2}
       (1+2\sigma W/r_{\min}).
\]

These follow directly from differentiating `ell/r`: the mixed derivative
contains `d_z/r^2 - 2*d*r_z/r^3`. Only absolute integrals are used. This is
an enclosure on the owned preparation slab, not an inferred duration.
The source and geometry are bound to the nonzero control's saved history.

The actual-node remainder factor multiplies the nth energy-derivative
majorant, with n equal to the node count. Physical columns use their original
weighted `A_up` norms and the required square-root-of-two conversion from
row norms. The carrier term is included in `F_z`. The finite-matter error
owner then propagates these field remainders into N and beta. Positive and
negative source covariances have separately bounded matrix norms; their error
uppers are added, without doubling a physical source contribution.

This certifies the ideal interpolation remainder at the stored binary nodes.
Errors in the numerical nodal evolution, interpolation arithmetic, upstream
preparation, retained source quadrature, the changed-history infinite tail
and the between-node constraint residual remain separate. A component PASS
is not a physical incoming-gate PASS.
