# FGC-1-CAL10-PREF15: terminal PROTO13 temporal-gate result

`FGC-1-CAL10-PREF15` binds and independently reconstructs the completed
PROTO13 GR-0 continuation. It is an outcome-neutral numerical calibration
result. It is not an FGC-QR trajectory, a collapse-mechanism result, a
physical obstruction, or evidence against Finite Gradient Closure.

## Frozen outcome

The authorization commit is
`682e6cd91419853a4e5bf132c163326f9cc6a029`. The campaign restored the six
method-owned states at `t=23/16`, advanced forty new common events, and reached
the first event at which the inherited causal-past temporal spectrum was
available:

```text
terminal event:       63
terminal time:        63/16 = 3.9375
history samples:      64
normal-flow tracers:  48
fields per tracer:    6
```

The run terminated normally with
`calibration_failed_no_eligible_GR0_case`. The typed stop was
`causal_past_temporal_spectral_admission`. Every member had zero source
retries. Ordinary adaptive CFL reductions occurred and remain separately
recorded; the stop was not CFL-retry exhaustion.

The public decision record is:

```text
PROTO13_campaign_terminated_normally = true
PROTO13_terminal_checkpoint_recomputed = true
PROTO13_terminal_spatial_and_constraint_admissions_passed = true
PROTO13_temporal_admission_failed = true
PROTO13_temporal_failure_cause_derived = false
GR0_case_eligible = false
temporal-gate diagnosis may begin = true
The terminal checkpoint may not resume = true
```

## What passed at the terminal event

Both complete method-owned common-event admissions still passed:

```text
RK4/D4-2:      2049 -> 4097 -> 8193
SSPRK3/D2-1:   4097 -> 8193 -> 16385
```

That means their reference-balanced physical constraints and spatial spectra
passed on their own frozen ladders. The trapped observable was evaluated on
the exact common 2,049-node physical mesh with separate Richardson intervals.
Its sign remained negative, so no trapped event occurred. This is merely the
state of the GR-0 calibration at `t=63/16`; the run stopped far before the
frozen final time `t=32`.

## What failed

The first 64-sample temporal assessment failed for both RK4 and SSPRK3. The
failure is not one Boolean hiding an otherwise passing calculation. Both
internal layers failed:

1. at least one individual windowed temporal power budget failed on each
   method and resolution ladder; and
2. the raw temporal field- and derivative-tail ratios did not all contract
   below the unchanged `1/4` threshold across adjacent spatial resolutions.

The serialized adjacent-ratio counts are:

| Method | Field ratios at or above `1/4` | Derivative ratios at or above `1/4` | Total comparisons per kind |
|---|---:|---:|---:|
| RK4 | 510 | 511 | 576 |
| SSPRK3 | 453 | 453 | 576 |

PREF15 restores all six terminal states and all six 64-by-48-by-6 tracer
histories from the atomic checkpoint. It independently recomputes the two
common-event admissions, the common-node trapped assessment, both temporal
admissions, and every serialized temporal ratio without calling the campaign
runner's terminal classifier. The external event log and result must equal the
checkpoint-embedded copies exactly.

## Interpretation boundary

This result says that PROTO13's frozen GR-0 calibration contract rejected its
own trajectory at the first temporal gate. It does **not** yet say why.

In particular, PREF15 does not decide whether the failed normalized temporal
budgets are dominated by:

- binary64-scale deviations divided by tiny total power;
- the short 64-sample window or its taper;
- a physical transient that is resolved rather than disappearing with spatial
  refinement;
- tracer/proper-time interpolation;
- insufficient temporal sampling;
- or a genuine continuum regularity problem.

Those alternatives make different predictions and require a separately frozen
diagnostic. Changing the `1/4` threshold after seeing the result, discarding
the temporal gate, resuming the terminal checkpoint, or opening SGB-L/FGC-QR
is not authorized.

## Successor permission

Only a bounded temporal-gate diagnosis may now be designed. It must retain the
raw terminal histories and separate absolute signal scale, window response,
temporal sampling, spatial-resolution convergence, and arithmetic enclosure.
It must prospectively state which observation distinguishes a broken
instrument from an admissible physical temporal spectrum before it evaluates
any replacement admission rule.

The GR-0 calibration remains incomplete and no case is eligible. SGB-L,
FGC-QR, COL1, DEF1, retained-EFT evolution, transition, singularity resolution,
child-domain, dark-sector, and varying-local-light-speed claims remain closed.

Reproduce the terminal binder when the ignored raw bundle is present with:

```bash
make fgc-cal10-pref15
```

A clean clone may verify the compact hash-bound certificate without the large
ignored checkpoint. A partial raw bundle always fails closed.

## Arithmetic and temporal provenance

The terminal states and histories were advanced under the authorization's
pinned CPython `3.14.3`, NumPy `2.5.1`, Accelerate, Darwin `arm64` runtime.
Exact recomputation of the floating common-event and temporal records therefore
belongs to that runtime partition. The independent analytic OpenBLAS partition
may verify the compact canonical certificate, source hashes, classifications,
and nonclaims, but it must not demand byte-identical regeneration of
Accelerate-owned reductions whose last bits are backend dependent.

Likewise, PRO13's two output roots were proved absent before launch. The
authorized calibration root now exists, so successor-time verification reads
the exact PRO13 result blob from immutable authorization commit `682e6cd…`
and requires it to equal the current tracked record. It does not falsely rerun
the historical absence observation after the state of the filesystem changed.
This changes neither the prelaunch fact nor any scientific gate.
