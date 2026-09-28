# Continuous constraint enclosure on 257 nodes

The 2026-09-28 Mac integration writes the reviewed successor to
`results/development/nsc-ks-continuous-constraint-enclosure-v2.json`.
It retains schema compatibility and records revision 2 and the original hash.
The original record remains byte-for-byte preserved under its earlier name.


The final local incoming certificate needs a continuous supremum of both raw
N,beta residuals on `I=S(1)+[.12,.18]`. A denser sample is not that supremum.
This successor owns the immutable 257-node / 256-cell contract and the
remaining assembly arithmetic that v4 left explicit.

The candidate is the live profile
`0b0e4cedfdb695d342c7f7347740631e916bdaee3a977df56808c725b25dacb0`. The grid
owner is `LocalIncomingFamily.collocation_nodes(257)`. The residual parts
come from the source-fixed evaluator chain
`KSHistoryEvaluator` / `KSCutoffBridgeEvaluator` /
`assemble_residual(baseline, geometry, matter, edge)`.

## Between-node contract

Production input is an outward nodal `|E_B|` enclosure at every verification
node and a directed derivative upper bound on every cell, independently for
`B=N,beta`. The Lipschitz intersection is the existing
[`nsc_ks_continuous_constraint_bound`](nsc-ks-continuous-constraint-bound.md)
primitive, restricted to 257 nodes. The owner rejects missing cells, a
reordered or non-increasing grid, endpoint movement outside a derivative
bound, a profile / state-law / source / dependency-hash mismatch, NaN or
negative values, and partial coverage. Dependency hashes are checked against
file bytes. Outward binary64 conversion uses the exact binary rational, so a
value such as `1/15` is not rounded below the bound. Samples may be stored as
diagnostics. They never replace the derivative majorant.

The production between-node component stays `null` until those directed
inputs exist. The historical geometry interpolation remainder is not the
full residual remainder.

## Remaining arithmetic

The remaining arithmetic scope is matrix products, contractions, signed-family
summation and baseline+geometry+matter addition. N and beta stay independent.
Phase contraction arithmetic is already enclosed by `phase_value` and is not
counted again. The exact-rational `gamma_n` term is applied last. Complete
120-family coverage cannot be inferred from a supplied array or a caller flag.
Until the original weights, signs, degeneracies and coherences are
authenticated from the source inventory, the production arithmetic component
stays `null` and `certificate_use=false` even when synthetic arrays are
present. Observed float drift is not a bound. Synthetic exact-rational
controls remain diagnostic.

## Record

`scripts/derive_nsc_ks_continuous_constraint_enclosure.py` writes
`results/development/nsc-ks-continuous-constraint-enclosure.json`. The
record is OPEN, lists the concrete missing arrays, and may contain a
manufactured Lipschitz / exact-rational control only under
`certificate_use=false`. It does not issue EXISTENCE or NON-EXISTENCE.

```sh
PYTHONPATH=src python3 -m pytest -q \
  tests/test_nsc_ks_continuous_constraint_bound.py \
  tests/test_nsc_ks_continuous_constraint_enclosure.py \
  tests/test_nsc_ks_assembly_arithmetic.py
python3 scripts/derive_nsc_ks_continuous_constraint_enclosure.py --check
```
