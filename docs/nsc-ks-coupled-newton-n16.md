# n=16 zero-pad evaluation of the accepted n=8 history

The first rectangular n=8 trial is accepted as a solver step. Its next n=8
Newton proposal does not improve the measured residual, so that candidate is
not evolved. This owner raises the Chebyshev count to 16 by padding the
accepted coefficients with zeros, then evolves all sixty families with 32
retarded directions on the 16-solve plus 33-verification union (47 nodes).

The same frozen upstream archive is restored. Operators from any n=8 history,
including group14, are not reused. The output is a measured residual of the
same accepted `g` at higher resolution, not a new Newton step and not a
physical root. The local incoming gate remains OPEN.

```sh
python scripts/derive_nsc_ks_coupled_newton_n16.py --run
python scripts/derive_nsc_ks_coupled_newton_n16.py --check
```
