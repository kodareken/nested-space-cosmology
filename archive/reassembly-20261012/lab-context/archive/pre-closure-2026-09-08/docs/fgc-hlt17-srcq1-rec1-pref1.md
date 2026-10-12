# FGC-1-HLT17-SRCQ1-REC1-PREF1 — independent post-result binder

This document owns the independent no-write binder for the completed
`FGC-1-HLT17-SRCQ1-REC1` remaining-member qualification. The tracked compact
certificate is complete. It is not a campaign, event execution, calibration,
or physical claim.

[REC1-FRZ1](fgc-hlt17-srcq1-rec1-frz1.md) remains the sealed one-run
authority. PREF1 authenticates that immutable commit. It does not import
REC1 runner, authority classification, or publication decisions, and it
does not rewrite Attempt2, REC1 raw, or any namespace.

## Status

Completed. The one no-write live reproduction independently returned the
classification below in 379.107645666998 seconds with 1104265216-byte peak
RSS. It wrote nothing and has not been rerun. The tracked config SHA-256 is
`a0681084...d1c63ead`; the canonical compact result SHA-256 is
`ea69bc1f...69b34e75` (5325 bytes).

If and only if the independent predicates close, the classification is:

```text
independently_bound_all_six_source_c1r1_qualified_no_state_advance
```

That class licenses only a separately frozen PRO20 first-event authority.
It does not license event execution, calibration, SGB-L health, candidate
execution, mechanism, or physics.

## Authority and freeze image

| Item | Identity |
|---|---|
| Authority commit | `be6a1792a30c39fc48ea7435368cbc8c17718d4c` |
| Parent A2 | `1fc17e5e2327cb3d8776386fc12d92a48af07f9c` |
| Exact freeze delta | `configs/fgc/fgc-1-hlt17-srcq1-rec1-frz1.toml`, `docs/fgc-hlt17-srcq1-rec1-frz1.md`, `mk/current-foundation.mk` |
| Freeze config SHA-256 | `6a11a65d…86989525` |
| Source closure SHA-256 | `48e507ca…fd258ce0` (248 files at A2) |
| Origin capture | `2824af06…2a83fe72` |

PREF1 authenticates that historical commit. It does not require later
development `HEAD` to remain at `be6a179`.

## Raw claims treated as facts to reproduce

The one REC1 raw terminal exists only in the canonical checkout:

```text
runs/fgc-2-sf1/hlt17-srcq1-rec1-source-origin-qualification/qualification.json
```

| Raw fact | Bound identity |
|---|---|
| Leaf SHA-256 | `ec9cd18b48ebfae93508182444ed5b689c566a7bc4cadcdc6a002d3550b17f16` |
| Leaf bytes | 105375 |
| Directory payload SHA-256 | `567d14095cc8631ed04715b332a4a35968d5a43dee1983bf062dc7647aa823ea` |
| Halt | none (`halted_outcome` is JSON null) |
| RK4-2049/4097/8193 | referenced Attempt2 prepared predecessors, not remeasured |
| SSPRK3-4097 | one independently evidenced CFL overlay, then `prepared_overlay` |
| SSPRK3-8193 / 16385 | `prepared_fresh` |
| State / endpoint / store | not advanced |
| Recorded wall / RSS | `371.6832666249975` s, `1133412352` bytes |

Those labels are never sufficient. PREF1 independently reconstructs the
source/C1R1 paths and compares complete facts. Outcome, disposition, and
admission text from the runner are compared only after independent
establishment.

Attempt2 remains the sealed predecessor:

- sealed compact `78604913…` at `efb53b83…`
- metadata-linked compact `ba656a63…` at `4891feab…`
- raw leaf `889a8d17…`, 45589 bytes, directory `7b06deb6…`

RK receipts are bound from Attempt2 raw, then confirmed as REC1
references:

| Member | Attempt2 C1R1 receipt |
|---|---|
| RK4-2049 | `25e6e1de…4211d8cb` |
| RK4-4097 | `005577c7…1bc83523` |
| RK4-8193 | `f3461dbd…9fc936f2` |

## Independent reconstruction

Live reconstruction, when later invoked, does the following and nothing
else:

1. Authenticate authority `be6a179` / parent A2 / three-path freeze
   delta, freeze config, live environment, and A2 source closure.
2. Authenticate both Attempt2 compact identities and Attempt2 raw.
3. Authenticate the REC1 raw leaf and one-leaf directory payload.
4. Inventory `runs/` excluding the REC1 namespace. The bound baseline is
   300 files, 150437927 bytes, sha256sum-line digest
   `77aa79f0300bd6866fcb491fe38aabf813a0064cbc7d511ef067364b4998ed88`.
5. Confirm Planck hashes against the collected-data manifest when the
   optional FITS files are present; absence is allowed, identity is not
   retuned.
6. Confirm the PRO20 event namespace is absent and that no endpoint,
   accepted state, or campaign store is present.
7. Rebuild the byte-bound origin/runtime from pro20, H17, and C1R1
   owners. Derive the six requested caps from retained
   `cfl_maximum * grid_spacing / inherited_previous_speed_upper`.
8. Skip RK remeasurement.
9. Evaluate the bound source once per newly owned SSP member.
   Source mutation checks, diagnostic normalization, cap derivation, and
   C1R1 assessment-to-decision reduction are implemented by PREF1 rather
   than imported from the SRCQ1/REC1 runner record builders.
10. For SSPRK3-4097, the frozen initial cap must independently yield
    canonical positive `observed_ratio > maximum_ratio` with kind
    `cfl_retry` (`0x1.000000003e80dp-3` over `0x1.0000000000000p-3`).
    Absorb that overlay only on an isolated in-memory checkpoint, follow
    only the returned `next_plan` cap `0x1.aaaaaa9908000p-10`, then
    reconstruct all eighteen C1R1 decisions, debits, and receipts.
11. Reconstruct SSPRK3-8193 and SSPRK3-16385 as fresh prepares.
12. Compare independent receipts, assessments, families, eighteen-channel
    facts, overlay counts, and plans to raw. Do not classify from runner
    labels.
13. Serialize no endpoint. Write no file.

Planck and historical-run hashing is streamed from stable no-follow file
identities; the live binder does not materialize the multi-gigabyte FITS image
as one Python byte string.

Recorded raw wall and RSS are authenticated as resource claims inside the
freeze ceilings (3600 s, 4 GiB). Live wall is not required to equal the
recorded wall.

## Compact verification

Default verification is raw-, runs-, source-, Git-, and shadow-blind:

```text
python -I -B scripts/reproduce_fgc_hlt17_srcq1_rec1_pref1.py
python -I -B scripts/reproduce_fgc_hlt17_srcq1_rec1_pref1.py --check
```

Those commands read only tracked compact config/result bytes. They
do not reconstruct, capture origin, or inspect `runs/` or Git. While those
files are absent or altered, `--check` fails closed.

## Completed live command history

The coordinator ran the following no-write reconstruction once from the
canonical checkout that holds the REC1 raw terminal:

```text
/opt/homebrew/Caskroom/miniconda/base/bin/python3 -I -B \
  scripts/reproduce_fgc_hlt17_srcq1_rec1_pref1.py --live
```

`--live` printed the compact result and config to stdout and wrote nothing.
Do not run it again. The reviewed exact bytes are now tracked at:

```text
configs/fgc/fgc-1-hlt17-srcq1-rec1-pref1.toml
results/fgc-1-hlt17-srcq1-rec1-pref1.json
```

Ordinary `--check` validates those tracked bytes only.

## Owned paths

| Path | Role |
|---|---|
| `src/recursive_horizons/fgc/evolution/hlt17_srcq1_rec1_pref1_binder.py` | Independent binder |
| `scripts/reproduce_fgc_hlt17_srcq1_rec1_pref1.py` | Default check / explicit no-write live print |
| `tests/test_fgc_hlt17_srcq1_rec1_pref1_binder.py` | Synthetic and adversarial tests |
| this document | Owner map and nonclaim boundary |

Lower-level mathematical/runtime imports from pro20, H17, and C1R1 are
allowed. Reduction, comparison, and classification are owned here.

## Nonclaims

PREF1 does not establish GR-0 eligibility, calibration completion, SGB-L
health, DEF1, ROB1, a production method, a candidate, a mechanism, or a
physical result. It does not authorize PRO20 event execution. Attempt2
remains `completed_terminal_reports_cfl_evidence_incomplete_no_state_advance`
and is not reclassified. REC1 runner labels are not scientific proof.
