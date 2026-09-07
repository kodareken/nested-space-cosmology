# Compact result archive

This directory contains the exact 58-artifact dependency closure that ends at
`NSC-2-ZETA1-RECURSION-MAP`. The machine-readable order, hashes, categories,
dependencies, and paper claim links are in [manifest.json](manifest.json).

The categories are intentionally different:

| Category | Meaning |
|---|---|
| `imported_benchmark_reproduction` | Established literature reproduced or bound as a constraint; not project novelty |
| `repository_derived_exact_identity` | An exact consequence of the stated repository premises |
| `repository_derived_numerical_result` | A finite numerical calculation with its recorded nonclaims |
| `diagnostic_nonpass` | A tested truncation or interpretation that did not pass |
| `superseded_candidate` | A useful intermediate value later removed from physical consideration |
| `current_frontier` | The latest derived equation and the next unresolved calculation |

There is exactly one current frontier:

- [nsc-2-zeta1-recursion-map.json](nsc-2-zeta1-recursion-map.json) derives the
  unitary parent/child dilation and the energy-resolved recursive tail. It does
  not select \(\Omega\), solve the mode-resolved tail, or derive physical
  \(\zeta\).

For the scientific interpretation, including the determinant roots invalidated
by anomaly compensation, read [the current-result document](../docs/current-result.md).
For isolated recomputation without rewriting these files, run `make reproduce`.
