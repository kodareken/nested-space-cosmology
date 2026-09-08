# Compact result archive

This directory contains the original 58-artifact v0.1.0 dependency closure that
ends at `NSC-2-ZETA1-RECURSION-MAP`, preserved byte-for-byte, plus the 59th
follow-up `NSC-2-ZETA1-UNIT-CLOSURE-CHECK`. The machine-readable order, hashes,
categories, dependencies, and paper claim links are in
[manifest.json](manifest.json).

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

- [nsc-2-zeta1-recursion-map.json](nsc-2-zeta1-recursion-map.json) derives a
  unitary parent/child dilation and the energy-resolved recursive tail. It does
  not select \(\Omega\), solve the mode-resolved tail, or derive physical
  \(\zeta\). Dilation alone does not derive \(\zeta=\Omega^2\).

The 59th record is a follow-up diagnostic, not a second frontier:

- [nsc-2-zeta1-unit-closure-check.json](nsc-2-zeta1-unit-closure-check.json)
  corrects the additive-gap units to \(q_j=\lambda_j/\zeta+\mu^2\) and proves
  that the inherited subtracted derivative is strictly negative on this finite,
  fixed-geometry family. It does not select a physical scale.

For the scientific interpretation, including the determinant roots invalidated
by anomaly compensation, read [the current-result document](../docs/current-result.md).
For isolated recomputation without rewriting these files, run `make reproduce`.
