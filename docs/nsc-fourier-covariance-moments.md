# Covariance transport with the two moments used by the constraints

For the periodically padded coefficient problem, expand a covariance-error
symbol as `D(z,k)=sum_l exp(i*omega_l*z) D_l(k)`. Define

\[
M_j(D)=\sum_l\int_{\mathbb R}{dk\over2\pi}|k|^j\|D_l(k)\|_1,
\qquad j=0,1,
\]

where the norm is the nuclear norm of each2x2 matrix coefficient. It is not
an asserted nuclear norm of the whole Weyl operator. This distinction avoids
using the false rule that a symbol L1 norm automatically bounds an operator
trace norm.

The owned Weyl kernel directly gives uniform coincidence bounds

\[
\|D(z,z)\|_1\le M_0(D),\qquad
\|J_D(z)\|_1\le M_1(D),\qquad
J_D={\partial_z-\partial_{z'}\over2i}D(z,z')\big|_{z'=z}.
\]

The derivative of the Weyl midpoint cancels in the difference, leaving the
factor k. Thus these are precisely the zeroth and first moments needed by
the existing N,beta insertion. Source energy has not been substituted for
outgoing momentum.

## Fourier-diagonal evolution and spatial mixing

Split the same Hamiltonian into a Hermitian Fourier-diagonal part H0(rho,k)
and the actual spatially varying potential delta V. The reference mass and
radius terms can remain in H0. For each l,k, the diagonal evolution acts as

\[
D_l(k)\longmapsto U_0(k+\omega_l/2)D_l(k)U_0(k-\omega_l/2)^\dagger.
\]

Both matrices are unitary, so the coefficient's nuclear norm is unchanged.
This remains true for the time-dependent reference Hamiltonian; it does not
require its matrices at different times to commute.

The [exact Weyl shifts](nsc-weyl-commutator-remainder.md) give the remaining
commutator as a sum of `delta V_r D_{l-r}(k-omega_r/2)` and the right-multiplied
opposite shift. Submultiplicativity of the matrix nuclear norm and translation
of the k integral imply

\[
\dot M_0\le2 V_0 M_0+R_0,\qquad
\dot M_1\le2 V_0 M_1+V_1 M_0+R_1,
\]

with `V_j=sum_r |omega_r|^j ||delta V_r||_op` and
`R_j=M_j(R_rho)`. Absolute preparation time is used. The first-moment
coefficient is V1: two commutator sides combine with the half-transfer
`|omega_r|/2`. It is not V1/2. A compact shifted-kernel test detects that
missing factor.

For `K_j=integral V_j |d rho|` and
`A_j=M_j(D_up)+integral R_j |d rho|`, the positive comparison system yields

\[
M_0(D_\Sigma)\le e^{2K_0}A_0,\qquad
M_1(D_\Sigma)\le e^{2K_0}(A_1+K_1A_0).
\]

These conditional estimates retain the initial error and the true operator
defect. Neither can be replaced with zero when unknown. They use the same
Dirac evolution as the Bloch/Sobolev helper, with additional absolute Fourier
summability as the input condition.

For the pure-radius history, `delta V_rho=(ell/a) q(z) sigma2`, where
`q=1/r_g-1/r_ref`. Certified Fourier-l1 moments of q give
`K_j <= (rho_up-1)|ell| A_j(q)/a_min`. This is the original rho clock and
potential. It neither changes the source nor supplies a new action term.

The existing Pauli vertex then gives

\[
|\Delta N|\le\mu\left[\sqrt{m^2+\ell^2/r^2}\,M_0+M_1/a\right],\qquad
|\Delta\beta|\le\mu M_1.
\]

A [continuum local embedding](nsc-ks-local-embedding.md) must be established
before identifying the periodic problem with the physical interval I.
Numerical operator error and source accuracy remain separate. In particular,
this helper is not a bound for a truncated input-energy source just because
it is a bound for the full auxiliary error `D=C-P`. The actual subtraction,
its initial covariance and any spectral split must be accounted for exactly.
The physical local gate remains OPEN.
