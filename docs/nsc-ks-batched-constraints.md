# Bounded-memory batched KS constraint assembly

Numerical adapter that evaluates the owned [KS local incoming split](nsc-ks-local-constraints.md)
on disjoint [reference source panels](nsc-ks-source-inventory.md) without
allocating one dense `C_src` across all energies. The local incoming gate
remains OPEN. No archived science is rerun.

Existing one-state-per-angular-family assembly still holds for a single
complete fiber run. The retained inventory has many disjoint panels and
energy batches; this owner consumes those batches.

## Reuse

| Input | Owner |
|---|---|
| Complete three-column energy fibers | `ReferenceSourcePanel` / `panel.batches` |
| Homogeneous upstream continuation | `ReferenceSourcePanel.fixed_upstream` |
| Joint reference/difference incoming state | `KSDifferenceIncoming` |
| Stable coherent matter difference | `source_fixed_matter_difference` |
| Channel multiplicity, ell=0 documentation | `_channel_record` / `_family_matter` |
| Opposite-angular negative source | `negative_angular_partner` / `map_negative_batch` |
| Local/reference geometry, D[0]=0 | `surface_geometry_response` |
| Error budget and residual assessment | `local_error_budget`, `assess_local_residual` |

No frozen incoming C0, new source cutoff, stress refit, `Gamma_rest`, or
metric step is accepted. LLL columns are not invented.

## Equation

\[
\mathcal E_B[g]
=
\mathcal E_B[g_{\rm ref}]
+\Delta\mathcal E_B^{\rm local+reference}[j_\Sigma g]
+\sum_{\mathrm{batches}}
\bigl(G_B[C_\Sigma[g]]-G_B[C_{\rm ref}]\bigr),
\qquad B=N,\beta.
\]

Each batch keeps its own block of `C_src`, weights once, and complete
three-column coherences. Because the retained source is block-diagonal
across energies, the contraction of a concatenated fiber equals the sum
of the batch contractions. Original ledger multiplicity

\[
\mu=\frac{\texttt{copy\_count}\times\texttt{degeneracy}}{n_{\rm actual\ angular\ signs}}
\]

is applied once per signed family, not multiplied by the number of batches.
Baseline and geometry are added once at finalization.

Positive and negative energies are both explicit. The negative family is
the supplied negative source mapped by the owned antiunitary identity
`(E,ell,X) -> (-E,-ell,S_3\overline{X})`. A positive result is never
multiplied by two, and `F(E,z)F(E,z)^\dagger=I` is not assumed.

## API

`bind_upstream_batch(panel, *, rho_up, rtol, atol, max_step) -> KSUpstreamBatch`

Calls `panel.fixed_upstream` once. The returned object binds panel identity,
original row range, `FixedSourcePreparation`, `A_up`, mass, actual angular
sign, energy sign, and `rho_up`. It contains no trial history.

`map_negative_batch(batch, negative_source) -> KSUpstreamBatch`

Applies the owned signed map to `A_up` and rebinds the digest from the
caller-supplied negative source. Same original rows, opposite energy sign.

`KSConstraintAccumulator(provider, baseline_gradient, coefficients, interval, error_budget, coverage=None)`

`add(prepared, batch, channel, batch_id=None)` accepts a genuine
`KSDifferenceIncoming`, the bound upstream batch, and the original channel
ledger. It checks source and `A_up`/preparation digests, common target `z`
and direction order, and rejects duplicate or overlapping original panel
rows per energy sign. Distinct computational meshes are allowed.

`finalize()` adds baseline and `surface_geometry_response` once, then
returns raw N,beta residuals, the full retarded Jacobian, per-batch matter
corrections, and an explicit OPEN scope.

`assemble_batched_incoming(entries, provider, baseline_gradient, coefficients, interval, error_budget, coverage=None)`

`entries` is a sequence of `(prepared, batch, channel)` or
`(id, prepared, batch, channel)`.

`signed_incoming_pair(prepared, batch, negative_source)` maps one evolved
positive batch to its explicit negative partner.

`declare_radius_invariant_channel(channel)` reuses the exact zero response
for the ledger's ell=0 groups 0,13,23. In this pure-radius history family,
the canonical KS generator depends on radius only through `ell/r`, and the
intrinsic incoming matter vertex is unchanged. Therefore both the matter
change and its history tangent vanish at every source energy. Each group's
baseline allocation remains in the baseline supplied to the accumulator.
This declaration neither constructs LLL columns nor certifies the baseline's
source accuracy. A second analytic declaration or a numerical batch for an
already declared group is rejected, in either insertion order.

Numerical batches still accept only the retained real-field groups 1–32;
LLL columns are not available through this API. Missing coverage and
error-budget components stay `None`/`OPEN`.

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_ks_batched_constraints.py
```

## Limitations

This is residual/Jacobian assembly for a finite disjoint source table. It
does not certify full spectral coverage, changed-history tails, preparation
field error, or a physical constraint root. Continuation, integration and
source-approximation errors remain the existing OPEN bounds.
