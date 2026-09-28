# Reciprocal-radius remainder for continuous validation

The finite-history radius potential can be validated by a short polynomial
plus an exact geometric-series remainder. This is a numerical enclosure of
the existing reciprocal radius, with no new action or matter term.

Write `x=delta r/r_ref`, with sup bounds `|x|<=eta<1`, `|x_z|<=eta1`,
and `|x_zz|<=eta2`. For order p and n=p+1,

\[
\frac{-x}{1+x}=\sum_{j=1}^{p}(-x)^j+R_p(x),\qquad
|R_p^{(j)}(x)|\le f_j\quad (j=0,1,2),
\]

where the following positive rational expressions bound the absolute tail
and its derivatives:

\[
f_0=\frac{\eta^n}{1-\eta},\quad
f_1=\frac{n\eta^{n-1}}{1-\eta}+\frac{\eta^n}{(1-\eta)^2},
\]
\[
f_2=\frac{n(n-1)\eta^{n-2}}{1-\eta}
+\frac{2n\eta^{n-1}}{(1-\eta)^2}+\frac{2\eta^n}{(1-\eta)^3}.
\]

The operator remainder bounds through two axial derivatives are
`|ell|/(a_min r_ref,min)` times `(f0, f1 eta1, f2 eta1^2+f1 eta2)`.
For the owned compatible family,
`eta_j <= (sigma W_j+sigma^3 U_j/6)/r_ref,min`. All operations in the
bound are exact rational arithmetic on the supplied constants.

## Axial profile constants

The helper also bounds the analytic `LocalAxialFunction` through order two.
Its binary Chebyshev coordinate-map coefficients are treated as exact and
the argument interval includes the entire cutoff support. Conversion to
power coefficients uses exact rational recurrence. Absolute coefficient
sums bound the polynomial and its two derivatives; this can be conservative
at high degree.

For the normalized transition `eta(u)=expit(1/u-1/(1-u))`, symmetry reduces
the estimate to `0<u<=1/2`. Set v=1/u>=2. Its logistic product is at most
`exp(2-v)`, the first inner derivative is at most `v^2+4`, and the second
has absolute value at most `2v^3`. Thus `|eta'|<=8`. For the second
derivative, bounding each term in
`exp(2-v)(v^4+2v^3+8v^2+16)` at its maximum and using `e>8/3` gives
`36+20.25+32+16<105`. These bounds apply on both transitions, and the
derivatives extend by zero at the flat endpoints. The product rule then
gives the stored profile bounds with the actual transition width.

These constants bound the represented analytic profile and the reciprocal
series tail. Floating profile evaluation, polynomial/FFT arithmetic and
the continuous time residual still need separate bounds. This helper alone
does not populate the physical gate's field-error budget.
