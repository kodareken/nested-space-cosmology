# Deciding bounds after assembled Chebyshev n=32

The n=32 assembled leftover gap is order `1e-3` above `3e-11`. Geometry
between-node `1e-13` cannot move that gap. UV, full between-node, field and
low/subgap stay `None`. No deciding epsilon spend.

```sh
python3 scripts/derive_nsc_ks_deciding_bounds_after_assembled_chebyshev_n32.py --record
python3 scripts/derive_nsc_ks_deciding_bounds_after_assembled_chebyshev_n32.py --check
```
