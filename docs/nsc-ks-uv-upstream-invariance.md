# Upstream scalar normalization and the incoming UV coefficient

This result concerns the **formal occupied-vacuum coefficient difference on the
incoming surface**. It does not rescale the physical source, prove the actual
upstream column error, or give a numerical C4/C_M. The implementation is
`src/recursive_horizons/nsc_ks_uv_upstream_invariance.py`.

## Lower coefficients from the actual operator

Use the existing L0, potential V and source sign s:

\[
L_0=\partial_\rho-S_3a^{-2}\partial_z-iV/a,\quad
V=-mS_1+\ell S_2/r,\quad T_s=\partial_\rho-sa^{-2}\partial_z.
\]

Write A0=e_s, A1_major=iq, A1_minor=l=a(m-is ell/r)/2, and
A2_major=R+ih. Then

\[
Q=m^2+\ell^2/r^2,\quad |l|^2=a^2Q/4,\quad
n_2=2R+q^2+|l|^2.
\]

The exact density/current continuity identity implies

\[
T_s n_2=-2s a^{-2}\partial_z|l|^2
=s\ell^2 r_z/r^3=s T_s q_z.
\]

Thus k=n2-s q_z is transported upstream data. **It is not set to zero.**
The two histories have the same k because they share preparation and the same
T_s characteristics (a depends only on rho and is unchanged).

For J_I=Im(A† A_z)-s e A†A and J_S3 defined with S3 inserted, the three lowest
coefficients (e^1, e^0, e^-1) are

| Bilinear | e^1 | e^0 | e^-1 |
|---|---:|---:|---:|
| J_I | -s | 0 | -s k |
| N bracket = A†VA + J_S3/a | -1/a | 0 | -a Q/2-k/a |

These follow by direct Pauli substitution, retaining q_z. No outgoing momentum
is replaced by its source-energy label. The action sign is checked against
`source_column_matter`: N is minus the N bracket, while beta is plus J_I,
with the common multiplicity/measure. The raw beta momentum vertex is -I and
the Gaussian action minus is applied once. This explicit convention, not a
shorthand sign in an older coefficient description, owns the new identity.

On Sigma, the declared radius history chi(s_normal)*(s_normal*w+s_normal^3 U/6)
vanishes. Therefore r_g=r_ref and a_g=a_ref. All three lower **differences**
above vanish on Sigma. The N e^-1 difference need not vanish elsewhere.

## What an unknown upstream constant changes

Within one formal branch, the recurrence fixes the minor coefficient
algebraically and transports the major coefficient. A change in the higher
upstream major constants can be represented, through order four, by

\[
A'(e)=f(e) A(e),\qquad
f(e)=1+c_2 e^{-2}+c_3 e^{-3}+c_4 e^{-4}+O(e^{-5}).
\]

The c_j are common to the two histories and independent of rho,z. This uses
z-homogeneous upstream data, not an arbitrary spatially varying reweighting.
L0 and the fast projector term commute with this scalar. More explicitly,
the triangular recurrence gives delta A2=c2 A0,
delta A3=c2 A1+c3 A0, and delta A4=c2 A2+c3 A1+c4 A0. At the upstream surface,
c2,c3,c4 can be selected successively to match the changed major data because
A0_major=1 and the histories agree in an upstream neighborhood. Every higher
coefficient must be transformed consistently.

A Hermitian bilinear with at most one carrier derivative is multiplied by
|f|². Its cubic inverse-energy coefficient changes by

\[
J_3' - J_3 = 2\operatorname{Re}(c_2)J_1
+2\operatorname{Re}(c_3)J_0
+(|c_2|^2+2\operatorname{Re}(c_4))J_{-1}.
\]

Taking history-minus-reference **on Sigma** makes this zero for both owned
N and beta vertices, since their lower differences vanish there. Consequently
these common scalar upstream constants do not determine the formal C4 difference
on Sigma. The fixed subtraction is not modified. This is not a claim that every
upstream state choice, other branch, coherence, or remainder is irrelevant.

## Why varying R alone gives the wrong test

In the massless case J3 contains h3_z+R q_z-q R_z-s n4 plus fixed minor terms.
A real upstream shift R->R+c also induces h3_z->h3_z+c q_z. With the common
fourth-order density boundary fixed, continuity gives
n4->n4+2c(n2-n2_up)=n4+2cs q_z. Their contributions to J3 are
`c*q_z + c*q_z - 2*c*q_z = 0`. Treating h3 and n4 as independent while varying
R leaves a spurious uncancelled term. The corresponding remaining N shift is
proportional to a Q; its difference vanishes when the incoming radii agree.

Negative controls retain a nonzero cubic shift when the radii do not match.
A spatially varying factor also fails the commutation condition:
`[L0,z] psi = -S3 psi/a^2`. Those cases are outside this theorem.

## Normalization example and remaining work

If the *full* upstream major is chosen real positive as
(1+|v|²)^(-1/2), with v=l_up/e+O(e^-2), algebra gives
Re(A2_up)=-|l_up|²/2 and Im(A2_up)=0. The existing
`vacuum_upstream_conditions()` describes that formal convention. Unit norm
alone does not fix the phase: multiplying by exp(i alpha/e²) preserves the
norm but changes Im(A2_up) to alpha. Neither observation establishes physical
source matching or its error. The cancellation theorem above does not require
n2_up=0 or use this normalized-gauge example.

Actual C4 still needs the constructed transported coefficients and their
contracted difference on I. C_M needs the controlled higher-order remainder;
the whole-slab minor-A3 kernel bound still retains upstream information.
Physical source and full incoming-gate certification remain OPEN.
