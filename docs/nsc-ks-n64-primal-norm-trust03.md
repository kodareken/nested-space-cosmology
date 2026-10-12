# Trust-region shrink of the measured ball direction

Parent is `6534e55b…`. The rejected step `aa7e1982…` fixes the N
quadratic coefficient of its own direction, `1.50329526e6`. Linear N
gain meets that term at norm `8.61998730e-10`.

This step is the same direction at `0.3` times that crossing. The
predicted N gain is `3.33` times the quadratic term there, and the
predicted beta gain stays above `1e-10`. Nearby 10–20% perturbations
of this direction lose the N margin, so they are not evolved. Rank 74
and above stay unevolved.

```sh
python3 scripts/derive_nsc_ks_n64_primal_norm_trust03.py --run --cpu-budget 14400
python3 scripts/derive_nsc_ks_n64_primal_norm_trust03.py --check
```
