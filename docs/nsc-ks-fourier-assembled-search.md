# Assembled Fourier-on-I EXISTENCE search

Fourier-on-I coefficients act on the iterate4 assembled residual
(baseline + geometry + matter + edge). The retarded Jacobian is the
iterate4 assembled map transferred by a slot change of basis. No new
matter columns. No Dirac unless the leftover can reach `3e-11`.

```sh
python scripts/derive_nsc_ks_fourier_assembled_search.py --check
```
