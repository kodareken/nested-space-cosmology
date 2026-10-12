# Clipped n=16 rectangular trial

This owner evolves the clipped n=16 rectangular candidate selected by
[the damped proposal](nsc-ks-coupled-newton-n16-damped.md). The unclipped
full-scale identity `2d2588c3…` is not used. The stalled n=8 candidate
`949a1881…` is not used.

The same frozen upstream archive is restored. Operators from any other `g`,
including the n=16 zero-pad history, are not reused. The evaluation keeps 32
retarded directions, the 47-node solve/verification union, and the 192-node
phase. The output is a measured residual. The linear prediction is not this
residual, and the physical local incoming gate remains OPEN.

```sh
python scripts/derive_nsc_ks_coupled_newton_n16_damped_trial.py --run
python scripts/derive_nsc_ks_coupled_newton_n16_damped_trial.py --check
```
