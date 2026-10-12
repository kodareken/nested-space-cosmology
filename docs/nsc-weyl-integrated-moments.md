# Momentum-integrated auxiliary defect on a preparation-time box

The operator-error moments use the nuclear norm of each2x2 Fourier
coefficient, summed in the axial mode and integrated in canonical k with
`dk/(2*pi)`. For any twice differentiable periodic matrix symbol,

\[
\sum_l\|a_l(k)\|_1
\le\sup_z\|a(z,k)\|_1+{L^2\over12}\sup_z\|a_{zz}(z,k)\|_1.
\]

The zero mode is bounded by the first supremum. Integration by parts twice
bounds each nonzero mode by the second divided by `omega_l^2`; their sum is
`L^2/12`. The nuclear norm of a2x2 matrix is bounded directly through
`sqrt(||A||_F^2+2|det A|)`.

The [uniform axial geometry](nsc-ks-uniform-axial-geometry.md) supplies the
required whole-cell bounds. A pointwise geometry box cannot be promoted to
this Fourier norm. The kinetic/time contribution integrates its `|k|^-5`
coefficient analytically. Small and large transfer use the
[weighted Fourier sums](nsc-radius-transfer-sums.md); each sign of k and the
measure enter once. Source/channel multiplicity and the rho integral remain
separate.

## Optional bounded auxiliary at low momentum

For error analysis only, define

\[
Q_{\rm aux}=P_0+\chi_K(k)(P_1+P_2+P_3+P_4),
\]

where `chi_K=1-plateau(|k|;K/4,K/2)` is the owned smooth plateau complement.
It is zero through K/4 and one from K/2. The actual physical source and the
physical reference subtraction are unchanged. Reconstruction must retain
`Q_aux-P_ref`; it is not a new action term or a free counterforce.

For `|k|>=K`, the cutoff is constant at the expansion point. All small-transfer
shifted arguments have magnitude at least K/2, where it is also one. The
existing fourth-order small-transfer identity is therefore unchanged. In
large transfer, a shifted argument can approach zero. The cutoff suppresses
higher auxiliary orders there, and their global norm can instead be bounded
by `(4/K)^(j+1)` times the regular B_j norm on `mu in [0,4/K]`.
No cutoff derivative is discarded in the claimed high-k domain. The low-k
residual and the finite reconstruction correction are separate required
terms if this auxiliary is used for a full covariance computation.

The pilot compares the original P0..P4 auxiliary and this bounded variant on
one time box, with the actual saved profile and original group14 parameters.
It encloses the entire axial cell. The split K is canonical momentum, not an
input-source energy cutoff. Whole-time coverage, initial covariance error,
low-k residual and the source/subtraction matching are not supplied. The
local incoming gate remains OPEN.
