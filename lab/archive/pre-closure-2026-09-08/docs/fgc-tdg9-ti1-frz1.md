# FGC-1-TDG9-TI1-FRZ1: same-state SSPRK3-tableau counterfactual freeze

`FGC-1-TDG9-TI1-FRZ1` is the direct, still-unexecuted successor to sealed
LOC2 binder commit `de7fb1dd52f85d6f618bd91e5e5b20953baa463f`. LOC2 proved that
all ten selected RK4-2049 temporal failures contain independently failing
endpoint-state and width-scaled-RHS components. TI1 asks one narrower question:

> At the exact same restored RK4-2049 predecessor state, does replacing only
> the RK4 explicit time tableau with SSPRK3 clear the ten complete-C failures?

The three predecessor checkpoints, accepted binary64 state and time, 2049-point
grid, SBP4 spatial RHS, regular-centre projector, transaction, tracers, temporal
ledger, retry widths, owned rows, channels, exact `p >= 3/2` contraction rule,
primary/v2 localizers, root-proof limits, and candidate ceilings remain fixed.
The only changed input to the seven-path shadow compositor is the explicit time
tableau selector.

The runtime constant named
`second_order_diagonal_norm_SBP_plus_SSPRK3` is used only to select its SSPRK3
tableau. It does **not** change the inherited RK4-2049 SBP4 spatial operator.
The experiment is therefore named
`same_state_SSPRK3_tableau_counterfactual_on_RK4_2049_SBP4_operator`. It is not
a production SSPRK3 comparison, independent-method agreement, convergence
study, accepted evolution, continuation, calibration, or physics result.

The bounded shadow may construct exactly three TDG6 compositors: seven shadow
proposals and 28 SSPRK3 stage/endpoint RHS records per retry, for maxima of 21
proposals and 84 records. It may publish only canonical `manifest.json` and
`terminal.json` under
`runs/fgc-2-sf1/tdg9-ti1/ssprk3-tableau-on-rk4-2049-retries-3-5`.
It may not invoke temporal admission, accept or commit a fine path, acquire a
campaign writer, append a journal, publish a checkpoint, add a fourth width, or
raise a localization budget.

Every completed occurrence records the original complete failure, the shadow
complete failure/pass/inconclusive booleans, both ownership classes, and the
exact transition matrix. A changed ownership label is not clearance. All ten
complete-C contractions must pass or be exact zero for the all-clear terminal.
Any unresolved threshold or bounded localization makes the result
inconclusive. No terminal selects a successor remedy or authorizes state
advance.

Ordinary verification reads only the compact prospective certificate. The
status command authenticates the sealed authority and reports an absent or
complete two-leaf output without executing a shadow.

```bash
make fgc-tdg9-ti1-frz1
make verify-fgc-tdg9-ti1-prelaunch
make status-fgc-tdg9-ti1
make run-fgc-tdg9-ti1
```
