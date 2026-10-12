# Assembled Chebyshev n=16 Gauss–Newton leftover stall

The iterate4 assembled Jacobian in the same Chebyshev n=16 basis is
recomputed against the next5 unclipped rectangular direction. The leftover
cannot reach `3e-11`. Clipped prediction is recorded only. No Dirac.

```sh
python scripts/derive_nsc_ks_assembled_chebyshev_gn_stall.py --check
```
