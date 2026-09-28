# Next clipped n=16 rectangular trial

This owner evolves the next clipped n=16 rectangular candidate selected by
[the accepted-trial proposal](nsc-ks-coupled-newton-n16-damped-next.md). The
unclipped full-scale identities and the stalled n=8 candidate are not used.

The same frozen upstream archive is restored. Operators from any other `g`
are not reused. The evaluation keeps 32 retarded directions, the 47-node
union, and the 192-node phase. The output is a measured residual. The local
incoming gate remains OPEN.

```sh
python scripts/derive_nsc_ks_coupled_newton_n16_damped_iterate.py --run
python scripts/derive_nsc_ks_coupled_newton_n16_damped_iterate.py --check
```
