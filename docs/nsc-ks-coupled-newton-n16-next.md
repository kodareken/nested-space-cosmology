# First n=16 Newton proposal from the zero-pad residual

After the n=8 rectangular search stalls, the accepted history is evaluated at
16 Chebyshev coefficients with 32 retarded directions. This owner forms the
first n=16 rectangular Gauss-Newton proposal from that measured residual.

The output is a linear prediction. Fresh evolution of the selected n=16
candidate is required before the residual can be accepted. The physical local
incoming gate remains OPEN.

```sh
python scripts/derive_nsc_ks_coupled_newton_n16_next.py --check
```
