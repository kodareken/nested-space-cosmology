# FGC-1-TDG8-RCV3-PREF2: outcome-neutral terminal binder

## Result

`FGC-1-TDG8-RCV3-PREF2` independently binds the completed REC1 attempt as a
finite numerical/instrumental terminal. It reconstructs the exact recovery
edge and every later durable transition without importing the recovery
authority, event runner, retry controller, campaign store, or evolution
decision code.

The authenticated terminal is:

```text
anchor:       generation 9 / sequence 10
recovery:     exact predicted generation 10 / sequence 12
terminal:     generation 29 / sequence 50
event:        23, target 3/2
member:       RK4-2049
owner:        temporal
reason:       temporal_retry_exhausted
retry count:  22
disposition:  invalid_terminal
```

No accepted post-restart macro-step was added. The affected descriptor remains
`77847126...`; its decoded physical `u,p,q` state remains `3cd9f576...`. The
five other members are byte-for-byte unchanged from the generation-nine
anchor. The terminal checkpoint preserves the generation-28 member map
exactly.

## Independent reconstruction

The explicit binder performs one read-only binding operation over two
independently stabilized terminal-store snapshots. Each snapshot uses two
non-following low-level tree scans, for four scans total. The binder proves:

- REC1 authority commit `98f2e4f...` and its exact config, compact result,
  authority module, and runner bytes;
- the independent PREF1 scanner/decoder at commit `46abf807...`;
- the non-following 115-leaf, 25,734,805-byte terminal manifest
  `46b57bc3...`;
- content addresses, raw hashes, parent links, and journal tips for all 30
  checkpoints and 51 journal records;
- the exact prospective recovery checkpoint `6dc263e7...`;
- 19 later rejection/cursor-transition pairs, followed by sequence 49's
  rejection and sequence 50's terminal lock;
- every TDG6 rejection's 18-channel inventory, numerical classification,
  accepted-state preservation, unchanged accumulated debit, and retry count;
- all 17 state descriptors, all 11 NPZ payloads, and all 99 finite little-
  endian `float64` arrays with `allow_pickle=False`;
- the exact terminal checkpoint `13eb3ffc...`, rejection `928cfd1b...`,
  terminal journal tip `29aee787...`, and terminal lock `298ec9bd...`.

The final attempted macro-step is `0x1.aaa9600000000p-30`; the rejected retry
step is `0x1.aaa9600000000p-31`. The record explicitly classifies the stop as
numerical, not physical.

## Interpretation

PREF2 records what happened; it does not decide how to repair it. Temporal
retry exhaustion before an accepted post-restart step means the present
numerical admission/runtime path could not progress this member under its
frozen rules. It does not test GR collapse, FGC-QR activation, trappedness, or
the Raychaudhuri target.

Therefore all of these remain false:

- common-event completion;
- GR-0 trapping calibration;
- SGB-L or FGC-QR execution;
- candidate-action success or rejection;
- mechanism or physical result;
- retained-EFT evolution or physical-transition authorization.

No successor remedy is part of PREF2. Any successor must be designed and
frozen prospectively under a separate artifact.

## Verification boundary

Ordinary verification is compact and never opens the ignored store:

```bash
make fgc-tdg8-rcv3-pref2
```

The explicit independent live proof is:

```bash
make verify-fgc-tdg8-rcv3-pref2
```

The live proof is read-only. It neither replays a proposal nor mutates a
checkpoint, cursor, journal, payload, lock, or accepted state.
