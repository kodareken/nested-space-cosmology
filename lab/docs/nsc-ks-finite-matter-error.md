# Numerical matter error for a fixed finite source table

The incoming fields and actual axial derivatives enter the same coherent
contractions as the existing constraint owner. This helper propagates their
numerical errors without changing the archived covariance or source weights.

For the archived source matrix C it bounds
`||C||_2 <= sqrt(||C||_1 ||C||_infinity)` with interval arithmetic, including
all off-diagonal coherences. It therefore need not silently assume that a
rounded matrix has operator norm exactly one. Its difference from the exact
physical source law is a separate preparation/source error.

If f,z bound the reconstructed weighted Frobenius norms of F,Fz and e,ez
bound their errors, with `gamma>=||C||_2`, then

\[
D=\gamma(2fe+e^2),\qquad J=\gamma(ze+fe_z+ee_z),
\]
\[
\epsilon_N\le\mu[\sqrt{m^2+\ell^2/r_{min}^2}\,D+J/a_{min}],
\qquad\epsilon_\beta\le\mu J.
\]

The potential norm uses the orthogonal Pauli coefficients. Quadrature weights
are already in the columns and the original signed-family multiplicity mu
enters once. Current and reference errors must both be counted for a
matter difference unless a correlated estimate is proved.

Endpoint field norms are bounded from the enclosed Fourier coefficients of
the saved trajectory. The axial norm retains the actual envelope mode
`2*pi*k/L-E`; it is not replaced by the source-energy label. These bounds
describe numerical field error for a finite source table. They neither
close spectral coverage nor certify a physical incoming gate.
