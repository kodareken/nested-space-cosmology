# Accepted first Newton trial and next rectangular proposal

The [first rectangular trial](nsc-ks-coupled-newton-trial.md) is accepted as a
solver step only. Its measured N,beta maxima on the complete 23-node union are
below the previous history, the radius bound stays positive, and all sixty
operator channels were evolved for that `g`. The linear prediction is not the
accepted residual.

This owner then forms the next full-state Newton/Gauss-Newton proposal from
that measured residual and retarded Jacobian. The selected rectangular
candidate still requires a fresh evolution. The physical local incoming gate
remains OPEN: source, field, tail, phase-quadrature and between-node errors
stay `None`.

```sh
python scripts/derive_nsc_ks_coupled_newton_next.py --check
```
