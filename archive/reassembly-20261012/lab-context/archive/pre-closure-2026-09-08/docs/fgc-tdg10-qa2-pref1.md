# FGC-1-TDG10-QA2-PREF1 — independent two-width nonpass terminal binder

## Current status

`FGC-1-TDG10-QA2-PREF1` independently binds the completed QA2 retry-4/5
SSPRK3-on-inherited-SBP4 exact complete-C robustness diagnostic. The compact
result is `results/fgc-1-tdg10-qa2-pref1.json`, SHA-256
`08f88f02966ea5936e417a3a6a2c583e11783341a3d76f4ef5dc821e992691f6`,
with class
`independently_bound_retry4_5_ssprk3_sbp4_exact_complete_C_two_width_nonpass_terminal`.

The historical `FGC-1-TDG10-QA2-FRZ1` certificate remains unchanged and
truthful: at freeze it was unexecuted, measured zero shadows, recorded no
channel outcome, set `width_robustness_passed=false`, and advanced no state.
The later diagnostic run does not rewrite that prelaunch history.

## Authenticated raw terminal

The committed QA2 authority is
`615144fe72b402835acd6dd89d28f89cbe690cce`, whose parent is the sealed
QA1-PREF1 commit `a7915925a23ceaff6a491b8d08c572a70d242985`. The gitignored
raw namespace contains exactly two canonical leaves:

| Leaf | SHA-256 |
|---|---|
| `manifest.json` | `8638ac8caf2fb4adc0ff419e2147080e78e13e4cb3d5640f38fd841a261fd1eb` |
| `terminal.json` | `252810f4a388ab87e382c33c7cdd7fa576be67a317752828a89c247def49867c` |

The runner terminal class is
`completed_one_or_more_width_all_channel_nonpass_no_state_advance`.
Retry 4 and retry 5 each constructed seven shadow paths, seven proposals, and
28 SSPRK3 stage-plus-endpoint records. Both widths report the same three
nonadmitted channels: `u:alpha`, `u:lambda`, and `u:R`. Across both widths the
raw work budget is therefore 14 paths, 14 proposals, and 56 records.

The historical campaign store remained exactly 115 leaves with snapshot
SHA-256
`5917d70de3080dcc6b7abeb7fdb7cc3370c6140cbeb37bde0c142d6a2d36a445`
before and after the diagnostic. No diagnostic endpoint, payload, accepted
state, cursor, journal, checkpoint, temporal admission, or campaign mutation
was published.

The independent binder authenticates those raw bytes and the 16-blob QA2
authority delta, restores the retry-4 generation-10 and retry-5 generation-11
predecessors, and reconstructs 14 shadow paths, 14 proposals, 56 SSPRK3
stage-plus-endpoint records, and 36 channel assessments without importing QA2
runner or authority decision code. Both exact-rational localization routes
agree for every channel. Independent reduction matches the raw terminal:
retry 4 and retry 5 each fail exactly `u:alpha`, `u:lambda`, and `u:R`; both
width classes are `one_or_more_channel_nonpass`; and
`width_robustness_passed=false`.

The frozen outcome is narrow: this exact SSPRK3-on-inherited-SBP4 plus
exact-complete-C remedy is rejected as a general remedy on the tested retry-4
and retry-5 neighborhood. No successor remedy is selected.

## Compact boundary

The one-time live construction authenticated the two-leaf raw namespace and
sealed store and left both unchanged. Ordinary verification consumes only the
tracked compact config/result in the current repository routing bundle and
remains raw-, store-, shadow-, Git-, and QA2-runner-blind:

```bash
make fgc-tdg10-qa2-pref1
make verify-fgc-tdg10-qa2-pref1
python scripts/check_repo.py --only-tdg10-qa2-pref1
```

These ordinary routes never fall through to `--live`, launch QA2, or access
`runs/` or the campaign store.

## Nonclaims

QA1-PREF1 licenses QA2 method design only. It never licenses RA1, diagnostic
endpoint transplantation, or old-member adoption. This independently bound
nonpass does not earn a production SSPRK3 method or comparator,
independent-method agreement, state advance, common event, GR-0 calibration,
SGB-L or FGC-QR execution, candidate execution, mechanism, retained-EFT,
transition, physics, publication, release, or push.
