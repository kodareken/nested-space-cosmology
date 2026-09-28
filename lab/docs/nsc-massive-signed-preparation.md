# Massive signed-spinor map and paired-state sufficiency gate

This result applies the existing paired-angular Dirac representation to T's
current basis. It supplies an exact operator map and a hard information
certificate for the mixed compact/angular groups. It does not assign a
spatial covariance to any massive group.

**C1b-M remains OPEN.** The minimum delivered here is the explicitly permitted
certificate that the stored two-component data do not specify the required
paired state. The LLL preparation and PDF v0.25.0 are unchanged.

## 1. Exact map to the existing T operator

Import the paired radial matrices from
[the full-spinor boundary owner](nsc-chiral-boundary.md). The angular label
$\eta$ is an existing eigenvalue-sign label, not a room or a new particle
copy. Its principal, angular and mass matrices are
$I\otimes\sigma_2$, $\eta_3\otimes\sigma_1$ and $I\otimes\sigma_3$.

The constant unitary

$$
W=\eta_1\otimes\frac{I+i\sigma_2}{\sqrt2}
$$

maps these coefficients to $I\otimes\sigma_2$,
$-\eta_3\otimes\sigma_3$ and $-I\otimes\sigma_1$. Thus the same PG
operator in T's paired basis is

$$
H_T^{(4)}=-iI_\eta\otimes(v\partial_\rho+v'/2)
-Nm\,I_\eta\otimes\sigma_1
-N\frac{\lambda}{r}\eta_3\otimes\sigma_3,
\qquad v=\frac{N}{q_{\mathrm{PG}}}\sigma_2-\beta I.
$$

The $\eta=+$ block is exactly the stored T expression. The opposite angular
block has the opposite sign of $\lambda/r$. The map is constant, so it adds
no basis-connection term. The fixed transparent seam is preserved; the old
reflecting spatial walls are not imported into this PG domain.

The radial signed-frequency involution is

$$
\mathfrak C_T=B_TK,\qquad B_T=\eta_1\otimes\sigma_3.
$$

Here $K$ is complex conjugation. As a differential-operator identity,
$B_TH_T^{(4)*}B_T^\dagger=-H_T^{(4)}$. For its Fourier symbol, this becomes
$B_TH_T^{(4)}(p)^*B_T^\dagger=-H_T^{(4)}(-p)$, including the shift term.
The T normal-current metric is preserved. Pullback through the existing
same-seam trace map gives $\eta_1\otimes\sigma_3$ in the stored seed basis.

This is an exact radial operator/domain identification. A complete charged
four-dimensional conjugation additionally needs the angular/gauge-state
identification; it is not inferred from this matrix identity alone.

## 2. Which groups close in two components?

The 32 massive retained groups separate by their actual operator coefficients:

| Groups | Coefficients | Radial signed-frequency map |
|---|---|---|
| 1–12 | $m=0$, $\lambda>0$ | $\sigma_1K$ in T; $-iK$ in the seed basis |
| 13, 23 | $m>0$, $\lambda=0$ | $\sigma_3K$ in both bases |
| 14–22, 24–32 | $m>0$, $\lambda>0$ | Paired map $\eta_1\otimes\sigma_3K$ |

For the mixed groups, a constant two-component involution on the same T
operator would have to reverse all three Pauli matrices. The principal
coefficient fixes the $\sigma_2$ condition; two distinct radii separate the
constant compact mass from $\lambda/r$ and fix the other two conditions.
The exact linear system $B\sigma_i+\sigma_iB=0$ for $i=1,2,3$ has only
$B=0$, which is not unitary. This is a local same-block obstruction, not a
nonexistence theorem for the four-component or full NSC theory.

A pointwise axis flip can reverse the multiplication term at a single radius.
When that axis varies with $r$, its derivative carries the required
basis-connection contribution $-ivB'B^\dagger$. Such a connection is part of
a position-dependent change of frame; it cannot be discarded while claiming
an intertwiner of the unchanged T operator. The record evaluates this term
on the existing geometry without inserting it as a new physical source.

The fourteen simple maps are operator maps, not completed global PG states.
Their scattering/Cauchy reconstruction and physical preparation still need
evaluation. None receives the LLL flow kernel.

## 3. Hard certificate for the stored mixed-group data

For each mixed group at positive recorded frequency, the full paired seed
covariance has the form

$$
\mathcal C_+(E)=\begin{pmatrix}X(E)&Z(E)\\Z(E)^\dagger&Y(E)\end{pmatrix}.
$$

Only $X(E)$ is in the serialized two-component block. Even granting the
charge-related completion condition

$$
\mathcal C_-(-E)=I-B_T\mathcal C_+(E)^*B_T^\dagger.
$$

as an additional constraint on the control matrices, it relates opposite
angular and frequency labels; it does not supply the
missing positive-frequency angular marginal $Y(E)$. There are 16 real
Hermitian coefficients in the paired matrix and four in the observed block.
The resulting linear observation kernel has dimension 12. This count is not
a claim of twelve free physical parameters: CAR and the physical preparation
can impose further conditions.

The executable certificate therefore uses a stronger finite witness. Keeping
the actual stored $X$ as control data, the two algebraic completions
$\mathrm{diag}(X,0)$ and $\mathrm{diag}(X,I)$ have exactly the same observed
block. Both obey CAR within the inherited numerical tolerance and admit the
same signed-frequency pairing. The unobserved marginal, and the negative-
frequency block in the original angular sector, differ by Frobenius norm
$\sqrt2$. This remains true even with angular cross-correlation set to zero.

**Neither completion is adopted as a state.** They are information witnesses,
not horizon-prepared solutions, occupations to fit, or spatial covariances.
They establish insufficiency of the stored marginal plus CAR and the radial
pairing. They do not establish insufficiency of the declared action or of a
future reconstruction from the existing horizon/incoming preparation.
An operator involution alone does not prove invariance of the physical state;
the actual preparation must establish its appropriate charge/gauge relation.

## 4. Exact data still required

The required owner is the horizon/infinity-to-paired-PG mode map, in the
current and gauge conventions of T. It must produce both angular marginals
$X,Y$ and their cross block $Z$ from the already declared $C_H$ and incoming
occupation. If a state symmetry makes $Z=0$, that must follow from this
identification. Degeneracy, identical-copy counting and reference momentum
pairing do not supply it. Unfolding a nonzero angular pair changes its matrix
description, not the already counted number of modes.

That same global map must supply the Cauchy/CAR normalization and projection
onto $J$ and its full complement. Without it, a list of seed matrices cannot
produce the requested $C_{\mathrm{PG}}$ or $G^{<,>,K}$. No massive spatial
state is manufactured by this record. C1b-full, transmitting branch jets,
the boundary remainder and Z4 remain open; the physical history, stress and
metric timestep remain unassigned.

The [result record](../results/development/nsc-massive-signed-preparation.json)
contains exact identities, per-group residuals, the mixed-group witnesses and
authenticated input hashes. Its focused verifier runs no historical source
generator.

```sh
python3 scripts/derive_nsc_massive_signed_preparation.py --check
python3 -m pytest -q tests/test_nsc_massive_signed_preparation.py
```
