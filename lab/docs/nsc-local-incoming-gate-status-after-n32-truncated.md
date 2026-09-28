# Best measured n=32 truncated residual

The honest small step reached
`(3.3176810525573086e-07, 5.478125919396171e-06)` with prediction minus
measure `(6.56009485570778e-10, 4.076400187454782e-12)`. A later
coefficient step of norm about `0.014` improved beta and missed N.
Both stay far above `3e-11`. The gate is OPEN.

```sh
python3 scripts/derive_nsc_local_incoming_gate_status_after_n32_truncated.py --check
```
