# Linear high-energy contraction at fixed axial transfer

This connection reuses the physical retarded projector-pair identity and
endpoint coefficients in the [fixed-transfer owner](nsc-incoming-fixed-transfer.md).
It concerns linear response about the reference history at fixed finite
Fourier transfer. It does not certify the finite changed history or supply
a numerical tail constant on the interval I.

Set `E_o=E+omega/2`, `E_i=E-omega/2`. At the fixed intrinsic incoming
surface the pure-radius response has

\[
Q=E^{-2}Q_2+E^{-3}Q_3+O(E^{-4}),\qquad
Q_2=-\frac{\ell a^2\widehat w}{4r^2}\sigma_1,
\quad P_0=\operatorname{diag}(1,0),
\quad P_1=\frac a2(m\sigma_1-\ell\sigma_2/r).
\]

The exact relation `P_o Q+Q P_i=Q` at order `E^-3` imposes

\[
P_0Q_3+Q_3P_0-Q_3=-\{P_1,Q_2\},\qquad
(Q_3)_{diag}=\frac{m\ell a^3\widehat w}{4r^2}\sigma_3.
\]

The two off-diagonal entries of `Q3` remain undetermined by this equation.
No complex-conjugacy relation at a single transfer is imposed on them.
The midpoint offsets first affect this projector relation at order `E^-4`,
since their first correction to `P1/E` is of order `E^-2` and Q starts
at `E^-2`.

The normal vertex is `H=H0+E H1`, with
`H0=-m sigma1+ell sigma2/r` and `H1=-sigma3/a`. Consequently

\[
\operatorname{tr}(H_1Q_2)=0,\qquad
\operatorname{tr}(H_0Q_2+H_1Q_3)=0,\qquad
\operatorname{tr}Q_2=\operatorname{tr}Q_3=0.
\]

These zeros hold separately for each angular sign. Keeping only Q2 would
produce the false normal coefficient `m ell a^2 what/(2 r^2)`; the forced
diagonal Q3 cancels it. Multiplicities and the action's overall sign do not
alter a zero. The shift vertex is proportional to the midpoint momentum
times identity, so its corresponding coefficients vanish as well.

With the existing fixed-transfer `O(E^-4)` pair remainder, the at-most-linear
vertices leave an `O(E^-3)` contracted linear vacuum remainder. The full
coherent source correction is still part of the source and requires its
existing exponential estimate. This coefficient identity neither makes the
constants uniform in axial transfer nor justifies Fourier inversion of the
bound. The previously owned large-transfer linear estimate and a quantitative
bound on the finite-history remainder remain distinct requirements.
