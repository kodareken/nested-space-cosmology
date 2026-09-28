# n=32 retarded evolve is not authorized

The n=32 geometry leftover equals the current residual: the 64-column
geometry linearization is rank/conditioning deficient
(`condition ~ 2.36e16`). Truncated SVD leftovers stay about `10^{-4}`.
No 60-family evolve is started.

```sh
python scripts/derive_nsc_ks_n32_evolve_refusal.py --check
```
