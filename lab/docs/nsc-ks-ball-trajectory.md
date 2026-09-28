# Ball representation of the saved field polynomial

`BallFourierSegment` reads the anchored trajectory's binary endpoints and
six dense correction vectors exactly. It converts the degree-seven time
polynomial and its spatial nodal data to ball-arithmetic Fourier
coefficients. The homogeneous reference contributes only to the zero mode.
The original quadrature column weights enter once.

The numerical period and Fourier origin belong to the reconstruction. The
envelope is reconstructed; the rapid source carrier is not Fourier
differentiated. Spatial derivatives multiply each mode by its actual
`2*pi*i*k/L`. Derivatives in the time variable are derivatives with respect
to the accepted step fraction. Affine restrictions preserve the original
polynomial and retain exact dyadic substep coordinates.

Power coefficients can be converted to Bernstein coefficients using

\[
b_k=\sum_{j=0}^{k}\frac{\binom{k}{j}}{\binom{n}{j}}a_j.
\]

This standard polynomial identity, together with the nonnegative
Bernstein partition of unity, supplies a continuous-in-time norm bound.
The implementation uses the DFT and directed balls from
[python-flint](https://python-flint.readthedocs.io/en/latest/acb.html).
The optional dependency is pinned in `requirements-validation.txt`; the
ordinary evolution environment is unchanged.

No field equation is solved in this conversion. An enclosure of the saved
polynomial becomes a solution error bound only after the continuum PDE
residual, initial error and geometric domain conditions are established.
Exact dyadic midpoint/radius data are used for stored bounds; decimal
displays are not substituted for them.
