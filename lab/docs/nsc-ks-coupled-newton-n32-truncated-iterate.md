# Next truncated n=32 step from a measured Jacobian

Each iterate reads a completed assembled residual, takes one
column-equilibrated truncated step whose norm stays at most 1, and
evolves that candidate from the same upstream source. The linear
prediction is not the residual. A step is kept only when the measured
residual improves on both components.

`--name` must match `[a-z0-9-]+`. `--campaign-best` makes the prediction
beat that measured residual instead of the immediate parent. Iterate6 is
the repair from iterate5 that beats iterate3, the campaign best at that
time. `--propose` rewrites only a proposal whose candidate identity and
coefficients already match. `--check` replays one saved family from its
operator payload and reassembles the residual. It does not evolve.

```sh
python3 scripts/derive_nsc_ks_coupled_newton_n32_truncated_iterate.py --run \
  --parent results/development/nsc-ks-coupled-newton-n32-truncated-trial.json \
  --cutoff 5e6 --scale 1 --name iterate2 --cpu-budget 7200
```
