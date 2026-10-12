# Exact KS envelope of the approved local prepared Dirac law

This owner is an economical re-representation of the same source-fixed
incoming Dirac law already used by the PG CF4/exact-phase path. The metric
family, source matrix, Pauli matrices, chart and plateau are reused. No new
action, force, Γ_rest, normalization or metric step is introduced. The
independent adaptive integrator exists to compare that continuum law with
stored PG CF4 arrays; it does not copy the PG discretization.

## Reuse

| Input | Owner |
|---|---|
| Radial KS generator `G_E` | `ks_generator` |
| Canonical map `R=inverse/r` | `canonical_map` |
| PG potential | `CommonTimeBulkSplit.local_potential` |
| KS frame `sqrt(-v_±)` | `nsc_nonlinear_ks_source` (identity only) |
| Radius family | `LocalIncomingFamily`, `CompatibleIncomingMetric` |
| Normal plateau `χ(s)` | `CompatibleRadiusDirection` |
| Source fibers/weights | `FixedSourcePreparation` |
| Covariance product | `restriction_covariance_tangent` |

Archived amplitudes at the actual `rho_up` remain caller data. This module
never loads a 7601-grid, the full retained source, or a physical producer.

## Envelope law

The approved local field of a source column with energy label `E` is

\[
F_E(T,z)=e^{-iEz}X_E,\qquad
i\partial_T X
=\bigl[-mS_1+\ell/r\,S_2-E/a\,S_3\bigr]X
-\frac{i}{a}S_3\partial_z X,
\]

with `a=a_{\mathrm{ref}}(T)`. There is no extra `r_z` term: the evolved
unknown is the half-density variable `u=r\sqrt{a}\,\psi`. To avoid inverting
`T`, integrate in `rho` from `rho_up≥1.03` down to `rho_σ=1`:

\[
\partial_\rho X
=\frac{S_3}{a^2}\partial_z X
+G_E^{\mathrm{ref}}(\rho)\,X
+\frac{i\ell}{a}\Bigl(\frac1{r_g}-\frac1{r_{\mathrm{ref}}}\Bigr)S_2 X.
\]

`G_E^{\mathrm{ref}}` is exactly `ks_generator`. The remaining terms are the
periodic envelope derivative and the actual-radius correction of the same
compatible family

\[
r_g=r_{\mathrm{ref}}
+\sum_j\alpha_j\,\chi(s)\bigl(s\,w_j(z)+s^3 U_j(z)/6\bigr),
\qquad s=T(\rho)-T(1).
\]

`w` and `U` are sampled once on the fixed `z`-grid. The RHS never rebuilds
the chart over `z`. Parameter tangents of the same linear operator are

\[
\partial_\rho Y
=LY
-\frac{i\ell}{a}\frac{\delta r}{r_g^2}S_2 X,
\qquad Y(\rho_{\mathrm{up}})=0
\]

by the fixed upstream preparation. `tangents='all'` or `'zero'` selects those
directions; the nonlinear envelope is always evolved.

## Carrier reconstruction

After the radial solve, the trigonometric polynomial of `X` and `∂_z X` is
evaluated at the target nodes from the Fourier coefficients. Then

\[
F=e^{-iEz}X,\qquad
F_z=e^{-iEz}(\partial_z X-iEX),
\]

and likewise for `Y`. The rapid carrier is never FFT-differentiated. Source
frequency is a column label of `G_E` and of this phase; it does not choose
the axial mesh.

## Domain, padding and preparation

The computational `z`-grid is supplied by the caller, uniform, periodic, and
centered at `S(1)+.15`. Default length is `.4`. The fixed physical interval
is `I=S(1)+[.12,.18]`. Axial support must be passed explicitly; the usual
window is `S(1)+[.09,.21]`. Targets must lie in `I` and in the computational
domain.

The continuum characteristic length is

\[
D=\int_1^{\rho_{\mathrm{up}}}\frac{d\rho}{a_{\mathrm{ref}}^2}.
\]

Padding beyond the active axial support must exceed `2D` (conservative).
Normal support must lie outside the initial `rho`, `a>0` on the trapped
interval, and the radius lower estimate must stay positive. Owned Chebyshev
functions supply a global coefficient bound; arbitrary callbacks supply only
a sampled estimate, explicitly distinguished in diagnostics. The periodic envelope is numerical padding only. It is not a
physical periodic `Π` boundary condition. Finite Fourier-grid causality and
truncation remain OPEN.

The fixed-preparation digest contains the source, `A_{\mathrm{up}}`, `m`,
`ℓ`, `rho_up` and `rho_σ`. It does not contain the numerical `z` mesh, so
the same physical preparation can be refined in Fourier resolution.
`require_history` binds the evaluated `w/U` samples, amplitudes, windows,
grid, rho endpoints and solver options. Matching those samples is not
callback identity.

Caller `rtol`, `atol` and `max_step` are required, stored exactly, and grant
no physical confidence. Covariance uses the original source matrix with
weights applied once.

The DOP853 solver retains only its last accepted field state. It uses the
same stepping implementation as solve_ivp and records the actual step count,
without retaining a full-field history whose memory grows with every step.

## Frame identity

`R=inverse_{\mathrm{trace}}/r` depends only on the fixed `β/a`; the radius
cancels. The current-frame map of the owned PG potential is the KS
potential,

\[
U^\dagger M_{\mathrm{PG}} U=-mS_1+\ell/r\,S_2,
\]

with residual at roundoff (`~2.5\times10^{-16}` on the probe). This is the
same generator already checked by `local_frame_checks`.

## Returned object

`KSEnvelopeIncoming` carries target `z`, `F`, `dF`, `F_z`, `dF_z`, weighted
properties, the original source matrix/weights/labels, `m/ℓ`, upstream
columns/`rho_up`, the immutable numerical profile/preparation bindings, and
diagnostics. It is not a PG `HistoryBinding` and it does not fabricate
`node_fields`.

Physical source error, continuum error, integration error and Fourier
truncation remain OPEN. Parent comparison of `α=0` and `α=.001` against
stored 7601 arrays is outside this owner.

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_ks_source_envelope.py
```
