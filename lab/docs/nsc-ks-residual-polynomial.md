# Polynomial part of the continuous KS residual

For a saved step of signed width h and fraction xi, the normalized residual
is `S=X_xi-h L_g X`. Hence integrating its norm over xi gives the same
quantity as integrating the physical rho residual over absolute d-rho.

The algebra combines the ball Fourier field polynomial with time polynomial
enclosures of the existing operator. For each envelope mode k and source
energy E, its diagonal rate contains the actual `2*pi*k/L-E`; the off-diagonal
rates retain both mass and the signed angular contribution. The radius
correction uses the same `i sigma2` coupling. The nonzero-wave test detects
discarding the axial envelope derivative.

The resulting time polynomial is represented in a common Bernstein degree.
For each spatial position, its row norm over the entire time cell is bounded
by the largest coefficient row norm. A spatial Taylor enclosure includes a
full-cell derivative remainder; global Fourier coefficient moments provide
the field part of that remainder. Geometry profiles supply their own
interval derivatives.

This bounds the polynomial residual part continuously in time and space.
Scalar time-coefficient remainders, the reciprocal-radius remainder and any
profile Fourier truncation error must still be added. Their omission cannot
yield a complete field error or a physical local gate certificate. No
reference stress, new occupation or force is supplied by this module.
