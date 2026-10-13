# Opposite-angular negative-energy map of the KS difference incoming state

The new inhomogeneous source-fixed KS incoming state is mapped to its actual
opposite-angular negative-energy partner. The original source, scales and
real radius history are unchanged. This is a signed-family identity, not a
factor-two folding and not a physical gate or source-accuracy claim.

## Reuse

| Input | Owner |
|---|---|
| Joint reference/difference incoming state | `KSDifferenceIncoming` |
| Actual-radius envelope PDE | `nsc_ks_source_envelope` |
| Homogeneous generator `G_E` | `ks_generator` |
| Opposite-angular antiunitary map | `ReferenceSourcePanel.negative_partner` |
| Negative-energy seed/source law | `PairedHorizonSeedMap.at_radius` for `E<0` |
| Weights once, full `C_src` coherences | `source_column_matter` / `FixedSourcePreparation` |

No existing hash-pinned file is edited. No `Gamma_rest`, metric step, source
generator or field campaign is introduced.

## Envelope law and antiunitary map

On a real radius history the envelope unknown of `F=e^{-iEz}X` obeys

\[
\partial_\rho X
=\frac{S_3}{a^2}\partial_z X
+\frac{i}{a}\Bigl(-m S_1+\frac{\ell}{r}S_2-\frac{E}{a}S_3\Bigr)X.
\]

The `z` derivative is part of the PDE. The same operator is mapped by

\[
(E,\ell,X)\;\longmapsto\;\bigl(-E,-\ell,S_3\overline{X}\bigr).
\]

`S_3K` acts on the two-component spinor. Because `S_3 S_1 S_3=-S_1`,
`S_3 S_2 S_3=-S_2` and \(\overline{S_2}=-S_2\), the full right-hand side
including \(\partial_z X\) and the actual-radius \(S_2\) term is
antiunitary. The split form that uses `ks_generator` on \(r_{\rm ref}\) plus
the algebraic radius correction is the same operator.

## Incoming state

`negative_angular_partner(prepared, negative_source)` requires a
`KSDifferenceIncoming` and an original `FixedSourcePreparation` for the
negative source. It checks

- `energies = -` (positive energies),
- identical column weights, applied once,
- the owned source-law diagnostic
  `C_neg -(I-\overline{C_{\rm src}})` without replacing `C_neg` by a
  rounded complement.

The map sends `A_up`, `A_ref`, `D`, `D_z` and the retarded envelope
tangents through `S_3K`. It then reconstructs

\[
F=e^{-i(-E)z}(A+D),\qquad
F_z=e^{-i(-E)z}\bigl(D_z-i(-E)(A+D)\bigr)
\]

and the analogous tangents, so the negative carrier is consistent with the
difference-state identities. The negative preparation digest is rebound
from the supplied source, mapped `A_up`, the same mass, \(-\ell\), and the
same `rho` endpoints. The sampled real history binding is copied, not
re-evolved.

Envelope differences are mapped, not set to zero. The result is a second
signed family with `energy_folding_factor=1`.

## Source law versus evolved covariance

The relation `C_neg=I-\overline{C_{\rm src}}` is the matching horizon source
law already used by `negative_partner` and by `at_radius` for `E<0`. It is
recorded as `source_complement_residual`. The supplied covariance is stored
as given; a roundoff-scale departure from the exact complement is not
overwritten.

After inhomogeneous evolution one must not assume the per-energy isometry
`F(E,z)F(E,z)^\dagger=I`. The incoming covariance remains

\[
C_\Sigma=F\,C_{\rm src}\,F^\dagger
\]

with weights applied once. The negative kernel is the transformed field with
the supplied `C_negative`. It is not `I-C_{\rm positive}(z)`. A
non-isometric point-column control distinguishes those two matrices.

## API

- `s3_conjugate(field)` — `S_3K` on the second-to-last spinor axis
- `actual_radius_envelope_rhs(X, Xz, rho, energies, mass, angular, radius)`
- `source_complement_residual(C_pos, C_neg)`
- `negative_angular_partner(prepared, negative_source)`

Physical source error and the local incoming gate remain OPEN.


## Massive paired-spinor information boundary

The earlier full-spinor PG identification uses the constant unitary
$W=\eta_1\otimes(I+i\sigma_2)/\sqrt2$ to give

$$H_T^{(4)}=-iI_\eta\otimes(v\partial_\rho+v'/2)
-NmI_\eta\otimes\sigma_1
-N(\lambda/r)\eta_3\otimes\sigma_3,
\qquad v=(N/q_{\rm PG})\sigma_2-\beta I.$$

Its radial signed-frequency involution is
$\mathfrak C_T=(\eta_1\otimes\sigma_3)K$, with
$B_TH_T^{(4)*}B_T^\dagger=-H_T^{(4)}$ as a differential-operator identity.
The constant change preserves the transparent seam and adds no connection.
It supplies neither a charged four-dimensional state identification nor
the needed global paired-state preparation.

Groups1–12 ($m=0,\lambda>0$) admit $\sigma_1K$ in T's two-component
basis; groups13,23 ($m>0,\lambda=0$) admit $\sigma_3K$. Mixed groups14–22
and24–32 require the angular-paired map. A same-block constant unitary
would need $B\sigma_i+\sigma_iB=0$ for all three Pauli matrices, whose
only solution is $B=0$. A varying axis instead carries the connection
$-ivB'B^\dagger$ and cannot be presented as the unchanged operator.

Stored mixed-group two-component data provide only $X$ in
$\mathcal C_+=\left(\begin{smallmatrix}X&Z\\Z^\dagger&Y\end{smallmatrix}\right)$.
There are16 real Hermitian paired coefficients and4 observed ones, giving
a12-dimensional linear observation kernel. Algebraic witnesses
$\operatorname{diag}(X,0)$ and $\operatorname{diag}(X,I)$ share the same
observed block and radial signed pairing but differ by Frobenius norm
$\sqrt2$ in the unseen marginal. Neither is adopted as a physical state.
C1b-M remains OPEN: both angular marginals, cross block, Cauchy/CAR
normalization and full-complement projection must come from the actual
horizon/infinity-to-paired-PG mode map. The signed KS identity above
does not erase this earlier information boundary.
