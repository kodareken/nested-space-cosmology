# Assembled Chebyshev n=32 finite-image stall

The zero-padded iterate4 history was evaluated with a fresh assembled
Jacobian: reused n=16 matter columns plus newly evolved high-mode
directions. The rectangular condition is about `1.05e14`. The unclipped
and clipped leftovers equal the measured residual
`(0.004107193193, 0.005790788098)` and stay far above `3e-11`. No linear
step is available, so no clipped trial is authorized. Named stall:
`finite_wu_n32_assembled_image_cannot_reach_existence`. This is finite-image
evidence, not scoped NON-EXISTENCE of the predeclared class.

```sh
python3 scripts/derive_nsc_ks_assembled_chebyshev_n32_stall.py --record
python3 scripts/derive_nsc_ks_assembled_chebyshev_n32_stall.py --check
```
