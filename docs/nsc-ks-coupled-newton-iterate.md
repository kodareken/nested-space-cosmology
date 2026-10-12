# Next rectangular Newton-trial evolution

This owner evaluates the rectangular candidate selected by the
[accepted first trial and next proposal](nsc-ks-coupled-newton-next.md). It
reuses the first-trial evaluator, the same frozen upstream archive, grid64,
degree32, twenty-three nodes and the 192-node phase. Operators from the
previous history, including the first trial and the group14 pilot, are not
reused.

The output is a measured residual. The new linear prediction is not that
residual. The physical local incoming gate remains OPEN.

```sh
python scripts/derive_nsc_ks_coupled_newton_iterate.py --run
python scripts/derive_nsc_ks_coupled_newton_iterate.py --check
```
