# FGC-1-RSP2-PREF14: checkpoint-recomputed finer-ladder result

**Authors:** Douglas Ek & ChatGPT 5.6 Sol

**Status:** machine-reproduced post-run numerical result; RSP2's target and
complete constraint admission pass, while calibration and candidate physics
remain unopened

```text
RSP2_study_terminated_normally = true
RSP2_terminal_checkpoint_recomputed = true
RSP2_target_order_cleared = true
RSP2_complete_constraint_admission_passed = true
RSP2_target_preasymptotic_on_tested_ladder = true
PROTO13_design_may_begin = true
PROTO13_frozen = false
fresh_GR0_dynamic_calibration_completed = false
FGCQR_holdout_execution_authorized = false
```

## Result

RSP2 was frozen at immutable commit
`81412c1b6b86648a289e865d8ab40b9bbab26370` before its raw namespace
existed. It then advanced exactly one amplitude-`3`, SSPRK3/D2-1 member on
`16385` points to the unchanged endpoint `t=23/16`. The run completed in
`8078.001129082986` wall seconds with `1771` accepted steps, `7084` accepted
stages, zero source retries, zero CFL retries, no typed stop, and terminal
state hash
`91e64d0864390cb4ac63cfdbdbbcb8ca5ad6c2f47cdd4d6741e05ef560adad0d`.

The two immutable CAL9 predecessor states and the new endpoint give the raw
owned-domain radial-momentum norms

```text
4097:  0.0060603980047813765
8193:  0.0021440219481765766
16385: 0.000542836914675905
```

The frozen binary64 enclosures, used only for order classification, are

```text
4097:  2.296474121976644e-9
8193:  3.5915945773012936e-9
16385: 6.4437699620611966e-9
```

After the declared subtraction, the effective norms are

```text
4097:  0.0060603957083072545
8193:  0.0021440183565819993
16385: 0.0005428304709059429
```

The independently reconstructed adjacent-grid orders are therefore

```text
4097 -> 8193:  1.4990947384363291
8193 -> 16385: 1.9817436463819333
```

The first value remains below the strict `3/2` gate by
`0.0009052615636708783`. The genuinely finer pair passes it by
`0.48174364638193334`. The CAL9 miss was therefore pre-asymptotic on this
tested ladder; it is not rounded, discarded, or reclassified.

Complete constraint admission also passes. The four finite finest-pair orders
are:

| component | `8193 -> 16385` order |
|---|---:|
| Hamiltonian | `2.1246424153634873` |
| radial momentum | `1.9817436463819333` |
| gauge-t | `1.8516167425060845` |
| gauge-r | `1.8428133958883794` |

The remaining six reduction-constraint pairs are enclosed at zero by the
declared accumulated-roundoff allowance. The minimum finite order is
`1.8428133958883794`, monotone refinement passes, both absolute guards pass,
and there is no non-target obstruction.

## Independent checkpoint reduction

PREF14 does not trust only the runner's terminal JSON and does not call the
runner's endpoint classifier. Its reproducer:

1. verifies the immutable authorization commit and every implementation blob
   named in the launch manifest;
2. binds the raw manifest, 25-record canonical event log, terminal checkpoint,
   and terminal result by SHA-256;
3. requires exactly one initial event, 23 ordered accepted `1/16` boundaries,
   and one terminal assessment;
4. requires the checkpoint's embedded event log and terminal result to equal
   the external files exactly;
5. restores the CAL9 `4097` and `8193` arrays and the RSP2 `16385` arrays from
   their checkpoints and verifies all three array-content hashes;
6. recomputes the unchanged full and owned-domain constraint residuals,
   roundoff enclosures, guards, statuses, both adjacent target orders, and the
   complete admission; and
7. requires exact serialized equality with both the terminal event and the
   external study result.

Threshold, grid, time, target-component, claim-promotion, accepted-stage, and
one-bit state mutations all fail closed. The ignored raw files remain outside
Git; the tracked certificate contains their hashes and the compact evidence
needed to audit the conclusion.

## What this permits

The result removes one numerical prerequisite and authorizes the **design** of
FGC-2-SF1-PROTO13. PROTO13 must be prospectively frozen before another
trajectory is read. Under the roadmap it uses method-owned resolution ladders:

```text
RK4:    2049 -> 4097 -> 8193
SSPRK3: 4097 -> 8193 -> 16385
```

Amplitude `3` is the only prospectively designated calibration case. Each
method must pass its own convergence and error budget, and cross-method
agreement must compare extrapolated or interval-enclosed observables rather
than pretending unequal grids are identical.

## Scientific boundary

RSP2 is not the fresh PROTO13 calibration, does not make a GR-0 case eligible,
and does not establish collapse or a trapped interval. It reads no SGB-L or
FGC-QR state, tests no regulator activation, computes no affine-null
defocusing margin, authorizes no holdout, changes neither EFT1's `3/12` result
nor RUN1's `7/8` stop, and derives no singularity resolution, child domain,
dark sector, varying local speed of light, or verdict on the general gradient
programme.

The precise conclusion is:

> On the frozen amplitude-`3` GR-0 SSPRK3 constraint study, the late
> radial-momentum order loss clears on the `8193 -> 16385` pair and every
> complete constraint-admission component passes. The next calibration
> protocol may be designed prospectively.

## Reproduction

With the immutable ignored raw bundle present:

```bash
python3 scripts/reproduce_fgc_rsp2_pref14.py
python3 scripts/reproduce_fgc_rsp2_pref14.py --check
```

The canonical configuration is
[`configs/fgc/fgc-1-rsp2-pref14.toml`](../configs/fgc/fgc-1-rsp2-pref14.toml),
and the compact result is
[`results/fgc-1-rsp2-pref14.json`](../results/fgc-1-rsp2-pref14.json).
