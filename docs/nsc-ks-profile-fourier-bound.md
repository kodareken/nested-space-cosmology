# Certified Fourier data for the compact axial profiles

The reciprocal-radius expansion separates each potential term into a
time coefficient and `f(z)=w(z)^p U(z)^q`. This owner bounds Fourier
coefficients of those actual compact profiles, independently of source
energy. It does not certify the changed state's source-energy UV tail.

The Fourier/aliasing argument is classical; see the treatment of aliasing
in [Trefethen and Weideman, SIAM Review 56 (2014)](https://doi.org/10.1137/130932132).
Here the cutoff is smooth but not globally analytic, so the bound uses a
finite derivative norm instead of claiming exponential convergence.

For period L and `B_p >= integral |f^(p)| dz`, integration by parts gives

\[
|\widehat f_k|\le C_p |k|^{-p},\qquad
C_p=\frac{B_p}{L}\left(\frac{L}{2\pi}\right)^p.
\]

The ball-arithmetic DFT of M exact-grid samples encloses the discrete
coefficient. Its difference from the continuum coefficient at `|k|<=K<M/2`
is at most

\[
A=\frac{2C_p p}{(p-1)(M-K)^p}.
\]

This follows from the exact alias identity, `|k+lM|>=|l|(M-K)`, and
`sum_(l>=1) l^-p <= p/(p-1)`. Each retained coefficient receives that
additional real/imaginary enclosure. For derivative order `j<p-1`, the
continuous omitted-band supremum is bounded by

\[
T_j=\frac{2C_p(2\pi/L)^j K^{j-p+1}}{p-j-1}.
\]

`B_p` is obtained by summing cell widths times interval-jet derivative
bounds over the entire compact support. The flat-endpoint Cauchy bounds
are included; cells that cannot be classified are subdivided. There is no
sample-max assumption or empirical tail fit. Profiles vanish with their
derivatives near the numerical period's endpoints, so the integration-by-
parts boundary terms vanish.

The current record prepares the four powers of w in the saved amplitude
.001, U=0 control. It uses p=8, 256 cells per transition (plus the inner
partition and required subdivisions), M=16384 and K=2048, at 120-bit
precision. It is capped at 120 CPU seconds and performs no field or source
solve. Coefficient balls are stored as exact dyadic midpoint/radius pairs.
Replay checks hashes, the alias/tail inequalities and real-profile
conjugacy from the stored arrays.

These data feed the continuous residual enclosure. The time coefficient
remainders, reciprocal-radius remainder, complete source coverage and
physical constraint residual still require their own bounds. The local
physical gate remains OPEN.
