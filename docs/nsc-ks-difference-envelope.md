# Reference and difference coordinates for the prepared incoming field

The source law remains `C_Sigma[g] = F[g] C_src F[g]^dagger` with fixed
upstream columns and `delta C_src = 0`. This implementation reuses the
[KS envelope equation](nsc-ks-source-envelope.md), source matrix, chart,
radius family and retarded variation. Its purpose is to avoid differentiating
roundoff in the homogeneous background when the history response is small.

Write `X=A+D`, where `A` is independent of `z` and obeys the owned reference
radial generator `A'=G_ref A`, initialized with the same `A_up`. In the same
ODE solve evolve

\[
D'=L_g D+\Delta G_g A,\qquad
Y'=L_gY+\delta G_g(A+D),\qquad D_{up}=Y_{up}=0,
\]

with

\[
\Delta G_g=-\frac{i\ell}{a}\frac{\delta r}{r_{ref}(r_{ref}+\delta r)}S_2,
\qquad \delta G_g=-\frac{i\ell}{a}\frac{\partial_\alpha r_g}{r_g^2}S_2.
\]

The reciprocal difference is evaluated by this algebraically identical
expression to reduce cancellation. The Fourier derivative of the spatially
constant `A` is exactly zero in the declared discrete operator, so adding
the first two equations gives the original equation for `X`. There is no
empirical stress correction. The tangent includes the full evolved state.

At the incoming surface reconstruct

\[
F=e^{-iEz}(A+D),\qquad F_z=e^{-iEz}(D_z-iE(A+D)),
\]

and the analogous retarded tangents. The implementation stores `A`, `D`,
and `D_z` as well as the standard prepared-state arrays, with the same
fixed-preparation digest and sampled-history binding. It accepts the same
physical interval `I=S(1)+[.12,.18]` and requires the upstream slice to be
outside the normal support, positive radius, and the inherited causal buffer.

The tests compare the sum equation with the original operator, the zero
history identity, nonzero evolution in both representations, and finite
differences of the field, its axial derivative and coherent covariance.
These are implementation checks. Integration, spectral, full-source and
between-node error bounds remain separate requirements for the local gate.
No metric timestep or new action term is introduced.
