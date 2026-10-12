# FGC-1-PRO19-CFL1-FRZ1 — exact CFL-contract continuation

## Result

The first authenticated PRO19 production attempt advanced `RK4-2049` through
five finite accepted macro-steps and then stopped fail-closed before publishing
a proposed CFL-retry record. The terminal status was:

```text
invalid_implementation_or_nonconverged_run
cfl_rejection cap reduction differs
```

This is an implementation-contract diagnosis, not a GR-0 or physical result.
No SGB-L or FGC-QR state was opened, no common event was completed, and no
calibration, trapping, activation, defocusing, or transition claim was earned.

## Exact diagnosis

The outer controller requested:

```text
0x1.aaa9612df9ba5p-8
```

TDG7 lawfully selected the slightly smaller lattice-safe macro width actually
attempted:

```text
0x1.aaa9612df9800p-8
```

The numerical attempt correctly proposed half of the actual macro width as the
next CFL cap:

```text
0x1.aaa9612df9800p-9
```

The durable coordinator incorrectly serialized the outer requested cap as the
attempted cap. Its schema therefore rejected a record that mixed two different
quantities. The corrected coordinator serializes the exact executed-plan macro
width and still requires the successor to be its exact binary half. The CFL
threshold, event target, physical input, retry cap of 32, and every scientific
stop remain unchanged.

## Preserved recovery boundary

The failed proposal restored its accepted boundary and published no rejection,
terminal, orphan, or staging suffix. The exact safe ancestor is:

```text
original authority commit: c11ba422ce49ddd6f6175de1a9d2da675b97f2de
original plan SHA-256:      f3b365978aeb2ae6063742807a3161bf5d6aa0c8392cd283fdec2a79d3e46060
checkpoint generation:     6
checkpoint SHA-256:         6277b3be4ffe1bd79b858e2b4770c244b88b9971cbff67ffb1a40d40f14e2e63
journal sequence:           5
journal tip SHA-256:        bf92cbbab28fc55259f1d6cd1249fa73d5e740f067e4e4204e0dd3d1cfb854e5
member:                     RK4-2049
accepted time:              0x1.78554de5a30e0p+0
```

The other five members remain at `23/16`. All source, CFL, and TDG6 retry
counters remain zero at the checkpoint.

## One-use continuation contract

This static freeze does not claim that its future corrected Git image is
already authenticated and does not authorize execution by itself. The
successor is not allowed to relabel a new Git commit as the campaign's
original authority. A post-commit preflight must authenticate the corrected
committed image, then bind that image to the original plan and the exact
generation-6 checkpoint.
Only the existing same-host, verified-dead writer takeover may retire the stale
lease. Manual lock deletion, checkpoint rollback, historical GEN0 reimport,
alternate-plan adoption, or candidate execution is forbidden.

The continuation preflight must pass without creating output or advancing
state before the explicit resume command can run.

## Claim boundary

Passing this static gate means only that the exact one-use continuation
contract is frozen. Execution becomes authorized only if the later committed
image and the unchanged generation-6 store both pass the independent
continuation preflight. It does not mean the event will complete, the
calibration is eligible, or any physical mechanism has passed or failed.
