# Continuum embedding of the local incoming interval

The canonical KS equation has principal part
`S3/a(rho)^2 * partial_z`. Its mass and radius-dependent terms are zeroth
order. Consequently each Dirac column has characteristic speed magnitude
`1/a^2`, independently of its source-energy label. On the owned slab,
`a>=4/5` gives the exact rational travel bound

\[
D\le (\rho_{\rm up}-1)/(4/5)^2.
\]

The bound reuses the established background lower bound. The history must
have positive radius; no radius floor or changed physical geometry is used.
The numerical cell contains the compact axial support with padding greater
than2D on each side. Its periodic continuation and the original compact
history therefore agree on the whole backward cone of I, with a strict
neighborhood margin. Uniqueness and finite propagation for the same Dirac
operator identify their continuum propagators on that cone.

This applies to the covariance as well as each column: the two upstream
arguments in `U C_up U^dagger` both lie in the backward cones. Supplying the
same upstream covariance retains its correlations. It is not necessary to
freeze the incoming covariance or to impose a physical periodic boundary.
The strict margin also permits local axial differentiation on I.

For the owned translation-invariant upstream source, exact Bloch resolution
has `k_n=q+2*pi*n/L`, cell-normalized basis `exp(i*k_n*z)`, and reconstruction
measure `dq/(2*pi) sum_n`. Partitioning the original continuous k integral
into these intervals gives this measure directly; there is no additional
`1/L` in the physical kernel. The coefficients remain the original
`A_up(E=-k) C_src(E=-k) A_up^dagger`. Intrafiber source coherences are retained.
This is an exact representation of the same source, not a new preparation.

The embedding is a continuum statement. A finite Fourier differentiation
matrix is nonlocal and needs its own numerical error bound against that
continuum equation. The completed fine-field proof supplies one such bound
for its saved four-energy control; it is not a bound for all other operator
tables or for an incomplete energy quadrature.

The rational margin verifier rejects insufficient padding, a foreign
preparation slab, or unestablished radius positivity. Its result closes only
the geometric domain-of-dependence condition, without a global parent/child
matching or metric-evolution requirement. The local physical gate remains
OPEN.
