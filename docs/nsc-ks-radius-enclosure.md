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


## Mixed geometry and continuous coupling inputs

The same history has $N=1$, $\beta_K=0$, $a=a_{\rm ref}(T)$ and
$r=r_{\rm ref}(T)+\chi(s)[sw(z)+s^3U(z)/6]$. Directed mixed Taylor jets
through total order9 use the recurrence $\rho_T=-a(\rho)$ to compose the
existing background rho-jets into canonical-time jets. The increment
$\rho(T_0+\tau)-\rho(T_0)$ has exact zero constant, including interval
basepoints; independent interval subtraction cannot establish that zero.
Coefficients divide derivatives by $t!z!$. The preparation slab is
$\rho\in[1,33/32]$, and unresolved plateau boundaries require subdivision.
Positivity is enclosed rather than supplied by a floor.

For whole-cell axial bounds, combine the time coefficients with the same
authenticated global profile norms:

$$[\chi(s)s]_t\,w^{(j)}(z)/j!
+[\chi(s)s^3/6]_t\,U^{(j)}(z)/j!.$$

All needed w/U orders and their original amplitudes must be present, the
analytic profile identities must agree, and the numerical cell must contain
the compact support. These boxes bound actual mixed derivatives, including
transition terms such as $\partial_T^4\partial_z^5r$; independent higher
geometric slots are not introduced.

The radius-coupling recorder supplies the previously missing continuous
integrals for the same reference/difference field equation. With
$r=r_0+\delta r$,

$$M=-\frac{i\ell}{a}\frac{\delta r}{r_0r}\sigma_2,
\qquad M_z=-\frac{i\ell}{a}\frac{\delta r_z}{r^2}\sigma_2.$$

For $-\sigma\le s\le0$, $|d\rho|=a|ds|$ and global $W_j,U_j$, define
$J_j=\sigma^2W_j/2+\sigma^4U_j/24$. Then

$$\int\|M\||d\rho|\le\frac{|\ell|J_0}{r_{0,\min}r_{\min}},\qquad
\int\|M_z\||d\rho|\le\frac{|\ell|J_1}{r_{\min}^2}.$$

The covariance-growth recorder also uses full nonlinear
$q=1/r_g-1/r_{\rm ref}$ Fourier moments through order9. Its comparison is

$$\binom{M_0}{M_1}_{\Sigma}\le
e^{2K_0}\begin{pmatrix}1&0\\K_1&1\end{pmatrix}
\binom{M_0(D_{\rm up})+\int M_0(R)}{M_1(D_{\rm up})+\int M_1(R)}.$$

The initial and true-defect moments stay missing. For ell=0 the comparison
is identity and the radius response is zero, with baseline retained.
Geometry/coupling enclosures alone fill no complete UV or physical gate.
