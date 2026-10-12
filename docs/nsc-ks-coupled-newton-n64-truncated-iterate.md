# Truncated n=64 step from measured columns

Iterate6 is the best measured n=32 history. Modes 32 through 63 of the same
`(w, U)` class are already measured there. This owner takes one equilibrated
truncated step, clips it to the unit ball, damps it by the requested scale,
and evolves that candidate from the same upstream source. A linear prediction
is not a residual. The unclipped rank-96 step is not evolved.

```sh
python3 scripts/derive_nsc_ks_coupled_newton_n64_truncated_iterate.py --run \
  --parent results/development/nsc-ks-coupled-newton-n32-truncated-iterate6.json \
  --kept-modes 80 --scale 1 --name iterate7 --cpu-budget 7200
```
