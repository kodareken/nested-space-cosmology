# Error on the local backward characteristic cone

This is a more local alternative to the global-period H2 estimate. For the
same canonical KS radius family, the principal matrix is diagonal and its
speeds have magnitude `1/a^2` in rho. Diagonal source-energy terms are purely
imaginary. They therefore do not amplify the l2 norm of a spin row along
its characteristic, including all original quadrature-weighted source
columns.

Let J(rho) be the backward cone of I, with radius enlarged by the remaining
characteristic travel distance. Histories and initial fields must agree on
this cone; no physical periodic identification is assumed. Supremum errors
over J(rho) have no incoming boundary term because the domain shrinks along
the extremal characteristic speeds toward I.

For rowwise field error e0 and first axial derivative error e1, the two
characteristics are coupled only by the existing off-diagonal potential.
The norm of that coupling is

\[
c_0=\sqrt{m^2+\ell^2/r_g^2}/a,
\qquad c_1=|\ell r_{g,z}|/(a r_g^2).
\]

If K0,K1 bound their absolute rho integrals on the cone, and R0,R1 bound
the integrals of the reconstructed PDE residual sup norms, variation of
constants and the triangular comparison system give

\[
E_0\le e^{K_0}(e_0+R_0),\qquad
E_1\le e^{K_0}[e_1+R_1+K_1(e_0+R_0)].
\]

The diagonal phase cancels in the row-norm equation; no exponential in the
source-energy label is introduced. The first differentiated equation
includes the required `B_z e` term. The owned global K1 bound can be reused
as a conservative cone bound. For example,
`K0 <= DeltaT*(m+|ell|/r_min)` is valid when the same positive radius bound
holds throughout the cone.

Convert two row bounds to the weighted Frobenius norms by sqrt(2):

\[
\epsilon_F\le\sqrt2 E_0,\qquad
\epsilon_{F_z}\le\sqrt2(E_1+E_{max}E_0).
\]

The existing coherent matter-error inequality then bounds N and beta.
This uses the actual axial derivative of the carrier and envelope, not an
outgoing-momentum substitution. Source/preparation and arithmetic errors
remain separate inputs. A finite grid's apparent causal support is not
used: the residual is tested against the continuum equation on the cone.

The implementation uses [python-flint ball arithmetic](https://python-flint.readthedocs.io/en/latest/arb.html)
for the exponential and exports directed dyadic upper bounds. Inputs must
already be continuous residual bounds. Samples, solver tolerances and
resolution differences do not satisfy that condition. The helper alone
does not certify a physical local gate.
