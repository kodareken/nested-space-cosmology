# CURRENT-HISTORY nodal phase-value enclosure v2

This is the successor of the [n=8 coupled-history phase record](nsc-ks-coupled-phase-accuracy.md) for the live local history `0b0e4ced…` (`results/development/nsc-ks-gate-history-lm-broyden.json`). The history is n=32 in `w` and `U`. The collocation is the current-family Chebyshev set on `I` (129 or 257). The upstream slice is the original archive binary64 `rho_up=1.0300000000000002`.

It reuses `directed_source_phase_z_enclosure` and the production unit-angular contraction `nsc_ks_value_evaluator.edge_change`: the physical one-direction `(w,U)` metric, not the 64 coefficient tangents. Numerical angular-square weights are `signed_family_angular_square_weights` of the original inventory (120 signed families; 107952 each). The true N,beta factor uses the exact integer 107952 and the enclosed background `a(1)`, and also records numerical `a`/weight/float-contraction error. Production 192 and 384 Gauss gradients are compared to that **true-gradient ball**, including midpoint displacement, not only the ball radius.

The component target is `1e-12` per N and beta at each requested node. Between-node phase remainder, source/field error, finite tail and the physical local gate stay `None` / `OPEN`. No source or field evolution is performed.

On the recorded 3-node current-history pilot, the inherited enclosure radius is about `1.08e-12` and `1.20e-12` in N at the endpoints of `I`, so 192 and 384 both miss `1e-12` there. The middle node is inside `1e-12`. 384 does not shrink that miss: production sits at the true-ball midpoint and the limiter is the directed remainder, not Gauss bias. Numerical `a`/weight/float contraction is about `1e-18`. Mean cost is about `0.7` CPU seconds per node, so 129 nodes are about `90` CPU seconds. That is an extrapolation, not a 129-node certificate.

The n=8 owner hardcoded 21 nodes and a 60 CPU-second cap. This owner accepts `--history`, `--nodes {129,257}`, `--record`/`--check`, `--output-prefix`, and `--cpu-budget`. Node proofs are exclusive, resumable files. `--check` ignores CPU fields. Default `--selection pilot` is nodes `0`, middle and last. All 129/257 nodes are refused until that pilot cost is recorded (`--allow-all`).

```sh
.venv/validation/bin/python scripts/derive_nsc_ks_current_phase_accuracy_v2.py \
  --history results/development/nsc-ks-gate-history-lm-broyden.json \
  --nodes 129 --selection pilot --cpu-budget 120 --output-prefix pilot --record
.venv/validation/bin/python -m pytest -q tests/test_nsc_current_phase_accuracy_v2.py
```
