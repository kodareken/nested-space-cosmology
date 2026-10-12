# Primal-block remeasure of the iterate6 radius

The iterate6 coefficients are padded with zero modes 32 through 63. The
radius polynomial is unchanged. The evolution uses 128 retarded directions
and `evolve_primal_controlled`, so those directions do not enter the DOP853
error norm. No coefficient step is taken. The comparison with the published
iterate6 residual is a solver change, not a Newton acceptance.

```sh
python3 scripts/derive_nsc_ks_n64_primal_norm_zeropad.py --run --cpu-budget 14400
python3 scripts/derive_nsc_ks_n64_primal_norm_zeropad.py --check
```
