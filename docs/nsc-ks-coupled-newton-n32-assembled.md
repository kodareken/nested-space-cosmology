# Assembled Chebyshev n=32 Jacobian at iterate4

Zero-pads the accepted iterate4 history `e15982e2…` from 16 to 32 Chebyshev
coefficients per function. Reuses the already-evolved n=16 matter columns by
index remap into the 64-direction layout. Evolves only the new high-mode
directions through the retained-source path, with the low-mode history frozen
as one amplitude-1 direction so the base radius stays the same physical `g`.

No slot-transfer. No geometry-only map. No invented matter columns. Do not
evolve `2d2588c3…` or `949a1881…`. The physical local incoming gate remains
OPEN until a later closed certificate.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python3 scripts/derive_nsc_ks_coupled_newton_n32_assembled.py --run --cpu-budget 1200 --max-new 60
python3 scripts/derive_nsc_ks_coupled_newton_n32_assembled.py --check
```
