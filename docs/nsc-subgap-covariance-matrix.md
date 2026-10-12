# Full matrix identity for the inherited subgap source

The source-accuracy audit found that four stress projections cannot generally
reconstruct the full canonical2x2 covariance. In particular, the massless
stress vertices omit sigma1. This owner instead uses the complete probe basis
`(I,sigma1,sigma2,sigma3)` and retains all matrix components.

In the owned positive subgap interval, the resolved columns have the form
`(R*v,w,0)`. The original horizon covariance is

\[
C_H=\begin{pmatrix}f&-is\\is&1-f\end{pmatrix}.
\]

Therefore the full real-frequency canonical source covariance is exactly

\[
C=f|R|^2vv^\dagger+(1-f)ww^\dagger+K+K^\dagger,
\qquad K=-isRvw^\dagger.
\]

The numerical Gram matrix is not replaced with I, and the reflection modulus
is not forced to one. The source's existing f and s are used unchanged. No
adiabatic reference is subtracted from this full matrix, and it is not
returned as an evolved incoming state.

For any2x2 matrix, all entries are recovered by
`C=(1/2)sum_i Tr(C*sigma_i)*sigma_i`, including sigma0=I. The diagonal
terms use real frequencies. The coherence admits the same analytic
continuation already used by the scalar subgap owner:

\[
K(z)=-i s(z)R(z)v(z)w(\bar z)^\dagger.
\]

The independently propagated dual at conjugate frequency is required.
No complex-frequency source covariance or holomorphic continuation of the
piecewise incoming occupation is invented. For a real polynomial weight,
its integrated coherence contributes `J+J^dagger` to the matrix moment.
Taking twice the entrywise real part of J would lose imaginary off-diagonal
coherence and is rejected by the tests.

This extends the existing bilinear identity, not the physical source. A
future source-moment evaluator can use these full matrices while the metric
history still acts through the actual Dirac propagator from the same fixed
upstream preparation. Freezing this source representation on the incoming
surface under a changed history would remain incorrect.

The record checks all304 archived positive subgap rows across38 signed
angular families against the original three-column source sandwich, using
independent archived unsewn fundamental columns. No field or source solver
is run. The identity residual is an algebra/replay check, not an accuracy
bound on those input modes or on an energy integral. Validated continuation,
quadrature/contour remainders and the upstream transformation remain OPEN.
