# FGC-1-PRO13-FRZ1: method-owned calibration-ladder restart freeze

**Status:** prospective, machine-reproduced protocol freeze. It advances no
trajectory and authorizes no calibration or candidate execution.

```text
PROTO13_frozen = true
PROTO13_method_owned_ladders_frozen = true
PROTO13_restart_manifest_frozen = true
PROTO13_successor_runtime_implemented = false
PROTO13_fresh_GR0_dynamic_calibration_authorized = false
fresh_GR0_dynamic_calibration_completed = false
GR0_case_eligible = false
FGCQR_holdout_execution_authorized = false
```

## Question owned by this artifact

RSP2/PREF14 established one narrow numerical fact: for the amplitude-`3`
SSPRK3 control, the radial-momentum constraint order is
`1.9817436463819333` on `8193 -> 16385`, and every finite complete-constraint
order on that pair is at least `1.8428133958883794`. The earlier
`4097 -> 8193` value of `1.4990947384363291` remains public. The finer result
therefore identifies the CAL9 miss as pre-asymptotic on the tested ladder; it
does not itself constitute a fresh collapse calibration.

PROTO13 answers the next prospective design question:

> What exact, immutable restart problem may test whether amplitude `3`
> produces the required GR-0 calibration geometry through `t=32`, without
> weakening either method's convergence burden or manufacturing a fresh
> initial state after seeing the RSP2 result?

Its answer is a six-member restart manifest and two method-owned ladders. The
answer is frozen before the successor runtime and before the new output
namespace exist.

## Frozen numerical experiment

All six members restart from their complete accepted state at
`t = 23/16`. The primary ladder is:

```text
RK4/D4-2: 2049 -> 4097 -> 8193
```

The independent comparator ladder is:

```text
SSPRK3/D2-1: 4097 -> 8193 -> 16385
```

Each method must independently pass its own complete constraint and spatial
error admissions. A resolution from one method cannot fill a missing
resolution in the other. Cross-method physical observables are compared only
through Richardson or interval enclosures at common physical radii; unequal
finest grids are never treated as collocated samples.

The frozen continuation has:

- branch `GR-0` and amplitude `3` only;
- restart time `23/16` and final coordinate time `32`;
- first new common event at `3/2` and common spacing `1/16`;
- recoverable checkpoints every `1/4`;
- the unchanged minimum finest-pair constraint order `p >= 3/2`;
- the unchanged maximum nested-tail ratio `1/4`;
- eight consecutive qualified trapped common events;
- the unchanged four-times-combined-error trapped-sign margin;
- the complete inherited SRC4 source, REF1 equations, constraints, health,
  boundary, scale, transaction, CFL, dissipation, retry, and orientation
  contracts.

The six restart members are not reconstructed approximately. Their `u`, `p`,
`q`, tracer positions, tracer proper times, event proper-time histories, event
field histories, stage counters, transaction counters, retry ledgers, and
causal debits are restored from the immutable PROTO12 and RSP2 checkpoints.
The freeze binds both the primary state hash and a hash of the complete restart
payload for every member.

## Immutable lineage and temporal provenance

The artifact consumes checkpoint commit
`ef5c78eb5c981883ea47dacb7a968bb660eff4dd`. It verifies the exact tracked
PROTO12, CAL9/PREF13, RSP2 configuration, and RSP2/PREF14 result blobs at that
commit and also verifies the two ignored restart checkpoint hashes.

Before generating the canonical result, it proves both new roots are absent:

```text
runs/fgc-2-sf1/proto13/calibration
runs/fgc-2-sf1/proto13/holdout
```

The freeze creates neither path. This absence is temporal prelaunch evidence.
After an authorized successor creates the calibration root, later result
verification must bind this immutable freeze commit rather than incorrectly
rerun the historical absence observation.

## Adversarial controls

The reproducer attacks and rejects mutations to:

- the primary ladder;
- the comparator ladder;
- the `3/2` threshold;
- the restart time;
- a candidate-authorization claim;
- one restart state hash;
- and the new namespace.

It also rejects a changed predecessor blob, checkpoint, state array, history,
counter, coordinate time, source/CFL retry ledger, or incomplete restart
payload. These checks make the freeze a machine contract rather than a prose
intention.

## Claim boundary

PROTO13 contains no new trajectory. It does not establish that amplitude `3`
forms a trapped interval, does not make a GR-0 case eligible, and does not
authorize SGB-L or FGC-QR. It changes no physical equation, coupling, initial
profile, source evaluator, constraint, health rule, boundary rule, threshold,
or outcome definition.

The successor runtime owner is **FGC-1-HLT11-MON11**. That artifact must
restore and re-admit the six states, implement recoverable continuation, prove
the new calibration namespace is unused relative to this immutable freeze,
and separately authorize exactly one GR-0 calibration. Only a later post-run
certificate may classify the calibration.

No retained-EFT, defocusing, singularity-resolution, child-domain,
dark-sector, varying-light-speed, or broader gradient claim follows.

## Reproduction

Before the successor namespace exists:

```bash
make fgc-pro13-frz1
```

The canonical result is `results/fgc-1-pro13-frz1.json`. The two ignored
checkpoint files named in the configuration must be present with their bound
hashes; the command advances no state and creates no run directory.
