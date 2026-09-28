# FGC-1-HLT17-SRCQ1-REC1 — Attempt2 recovery qualification

This document owns the committed `FGC-1-HLT17-SRCQ1-REC1-FRZ1` freeze for the
prospective recovery harness. It is not a completed qualification, campaign,
compact result, or physical claim. It authorizes exactly one remaining-member
measurement only from the clean committed freeze image and only when the
read-only status target returns `safe_to_run=true`. The accepted
[PLAN](../PLAN.md) owns the SF1 finish line.
[SRCQ1](fgc-hlt17-srcq1-frz1.md), [PRO20-EV1](fgc-pro20-ev1-frz1.md),
[HLT17-MON17](fgc-hlt17-mon17.md), [SRCB1](fgc-pro20-source-binding.md),
[SRCF1](fgc-pro20-source-factory.md), and [C1R1](fgc-tdg11-c1r1.md) remain
the retained origin, binding, factory, and C1R1/HLT17 contracts.
AGENTS.md stays unchanged creative context.

Attempt1 and Attempt2 implementation, config, compact results, and the Attempt2
raw terminal are sealed. This recovery does not edit them and does not
reinterpret Attempt2 as an independently reproduced CFL stop.

Implementation A is retained in `4da705f...`; the live-checkout test adaptation
is `e379ed4...`; authority implementation A2 is
`1fc17e5e2327cb3d8776386fc12d92a48af07f9c`. Freeze B adds only this owner,
the exact config, and current Make routing as the direct successor of A2.

## Why REC1 exists

Attempt2 authority is consumed. The three RK4 members prepared all eighteen
channels. SSPRK3-4097 completed origin source evaluation and the runner
reported `cfl_retry_required`, but the raw terminal omitted
`result.evidence` observed/max ratio. Attempt2 is bound as
`completed_terminal_reports_cfl_evidence_incomplete_no_state_advance`.
Runner labels alone never prove a CFL stop.

## Frozen Attempt2 identities

These are sealed predecessor facts, not a new measurement. REC1 authenticates
**both** compact identities. It does not replace the sealed hash with the
later metadata-linked hash, edit the compact result, or alter any scientific
field. Future authority independently binds the historical Git blob and the
exact one-field metadata extension.

| Identity | Value |
|---|---|
| Compact path | `results/fgc-1-hlt17-srcq1-attempt2.json` |
| Sealed historical compact | commit `efb53b839a36a7f7f52018a44ef4eeaee42c1752`, SHA-256 `7860491360fe058d6c5019f780103e6bf13b2c243d3c1d59ec3f74109636ab5f`, no `derivation_document` |
| Metadata-linked compact | commit `4891feab58922e449a06585a848eb6a11df994ea`, SHA-256 `ba656a63343a87b389e388e39e142e8cfc94194c96c33326f6da0963f4504c48`, `derivation_document=docs/fgc-hlt17-srcq1-frz1.md` |
| Raw leaf | `runs/fgc-2-sf1/hlt17-srcq1-source-origin-qualification/qualification.json` |
| Raw leaf SHA-256 | `889a8d17551a1eea6ed56f43f602fefb3f61ab78637bbdd0c995905e22259d89` |
| Directory payload SHA-256 | `7b06deb68c989734f8dd669f2f1a8f4531aa90b82cc0e9ff34330db363123009` |
| Raw leaf bytes | 45589 |
| Origin capture | `2824af06c674336d479ec53faa4269f02fb598055860ec22f2ad84602a83fe72` |
| Attempt2 authority | `f3a31dc6fce35baa90d305a61e89cee9706fed38` |
| Attempt2 implementation | `7123aa385b3f4165186acd241f532b58ebbb3a9b` / `e801532e...` |

The metadata-linked blob is the sealed scientific payload plus exactly one
added field. A `ba656a…` digest is accepted only with that exact
`derivation_document`. Altered scientific content or a different metadata
field is refused. The sealed `786049…` identity remains explicitly bound.

State-safe flags that REC1 authenticates and must preserve:

- `accepted_state_advanced=false`
- `campaign_store_written=false`
- `endpoint_adopted_or_serialized=false`
- `pro20_namespace_absent=true`
- `retry_authorized=false`

RK4-2049, RK4-4097, and RK4-8193 remain predecessor evidence. REC1 does not
remeasure their member qualification.

## Owned implementation

| Path | Role |
|---|---|
| `src/recursive_horizons/fgc/evolution/hlt17_srcq1_rec1.py` | Deterministic no-write recovery core |
| `scripts/qualify_fgc_hlt17_srcq1_rec1.py` | Guarded CLI; default plan/status only |
| `tests/test_fgc_hlt17_srcq1_rec1.py` | Small synthetic/mocked objects and temporary output |
| this document | Prospective map, resource recommendations, freeze-B seam |

The core authenticates Attempt2 hashes and state-safe flags; rebuilds the
byte-bound origin/runtime under a later freeze; skips RK members; resumes at
SSPRK3-4097 in canonical remaining order; evaluates the bound source once per
newly owned member; prepares at the already frozen SRCQ1 prospective initial
cap; retains the complete `HLT17BridgeResult` evidence and attempted/next
plans; independently requires evidence kind `cfl_retry`, canonical finite
positive ratios, and `observed_ratio>maximum_ratio` before any CFL stop;
requires bridge `sink_writes` to be the integer zero; absorbs source/CFL
overlays only in an isolated in-memory checkpoint; follows only returned
`overlay.next_plan`; accounts source and CFL owners independently through
the frozen 32/32 caps, with lawful exhaustion at owner count 33 and mixed
source/CFL sequences inside the total bounded loop; retains full per-attempt
bridge evidence when evidence validation or resource checking stops a
member; and stops typed on source/CFL exhaustion, temporal, inconclusive,
resource, or invalid outcomes.

If SSPRK3-4097 prepares, the walk continues to SSPRK3-8193 then
SSPRK3-16385. The result combines predecessor RK evidence references with
newly qualified SSP records and full fingerprints/resources.

There is no store or campaign publication, no accepted-state advance, and no
corrected or shadow endpoint serialization. Atomic terminal writing of the
qualification record may use neutral `evidence_io`. An uncertain publication
is a typed stop. Old `run_fgc_gr0_calibration_v*` runners are not imported.
Attempt2 authority code is not imported.

## API and state contract

Public core APIs:

- `plan_status()` — identities, frozen caps, Attempt2 pins, recommended
  ceilings, seam text; no I/O
- `authenticate_attempt2_predecessor` — exact compact/raw hashes and
  state-safe flags
- `require_independent_cfl_ratio` — `observed_ratio>maximum_ratio`; runner
  labels are never sufficient
- `retain_bridge_result` — complete endpoint-free `HLT17BridgeResult` record
- `isolate_working_checkpoint` — overlay absorption is forbidden on the live
  member checkpoint
- `qualify_srcq1_rec1` — no-write remaining-member walk
- `qualify_and_publish` — authority seam, then exclusive directory publish
- `require_authority_delta` / `load_coordinator_authority_delta_validator`

Initial caps remain the frozen SRCQ1 formula/hex. Subsequent caps and plans
come only from persisted HLT17 overlay objects. No post-result fitting,
manual halving, grid/method change, or SSPRK3-on-SBP4 hybrid.

Prepared dispositions after overlay follow-through are `prepared_fresh`,
`prepared_overlay`, or `prepared_successor`. A CFL classification requires
kind `cfl_retry`, canonical finite positive `observed_ratio_hex` /
`maximum_ratio_hex`, and `observed_ratio>maximum_ratio`. Missing ratio
evidence is `inconclusive`, not a CFL stop. Nonzero or malformed
`sink_writes` is `invalid`. Source and CFL overlay counts are independent;
there is no shared 32-overlay ceiling.

## Resource ceilings — coordinator freeze-B values

These remain the SRCQ1 recommended guards. They are not scientific
thresholds and were not fitted to a measurement.

Already frozen and not retunable here:

| Control | Value | Source |
|---|---|---|
| Source retry cap | 32 | HLT17 / PRO20-EV1 |
| CFL retry cap | 32 | HLT17 / PRO20-EV1 |
| Temporal retries per macro step | 32 | TDG6 |
| Per-member enclosure | production N and 8N/16N | HLT17 production table |
| Initial requested caps | SRCQ1 formula/hex | SRCQ1 freeze B |

**Recommended conservative process ceilings for freeze B to pin or replace.**

| Ceiling | Recommended default |
|---|---:|
| Process RSS | 4 GiB |
| Per-member wall | 900 s |
| Total wall | 3600 s |

`coordinator_freeze_required=true` and `scientific_threshold_fitted=false`
are recorded on every result. Freeze B confirms these values unchanged.

## Freeze-B committed image

The tracked config is
`configs/fgc/fgc-1-hlt17-srcq1-rec1-frz1.toml`. It binds:

- implementation parent `1fc17e5e2327cb3d8776386fc12d92a48af07f9c`;
- core `faf15958...`, runner `fd594702...`, and authority `294ae85e...`;
- source closure `48e507ca674c7d164b905bb47d9c6778902886979ddb3dbe906e7d5cfd258ce0`;
- the exact CPython 3.14.3 / NumPy 2.5.1 / Accelerate arithmetic image;
- both Attempt2 compact identities, raw leaf/directory identity and exact
  six-member source-configuration map;
- origin capture `2824af06...`, unchanged caps, independent `32/32` retry
  owners and the existing production enclosure/resource ceilings;
- absent REC1 and PRO20 event namespaces, no active scientific runner, and
  false state/store/endpoint/campaign authority flags.

The config is deterministically emitted from A2 and is not an authority by
itself. The validator authenticates the config, exact direct-successor Git
delta, live bytes, environment, raw evidence and process boundary together.

## Authority seam and one real run

Implementation A does not accept a boolean skip/force flag. Real CLI
execution requires both of:

```text
--authority-commit <exact 40- or 64-digit lowercase commit>
--output-directory <absent namespace>
```

Default CLI (no arguments) prints `plan_status()` and does not capture the
historical origin, authenticate Attempt2 from disk beyond the pinned
constants, call the physical source, or write. Incomplete real-run arguments
refuse. An existing output namespace refuses. A missing parent directory is a
typed absent stop, not a mkdir workaround.

Authority/delta validation is implemented in A2 and activated only by the
committed freeze-B delta:

```text
recursive_horizons.fgc.evolution.hlt17_srcq1_rec1_auth1.validate_authority_delta
```

The real CLI path stops at this seam before origin capture and factory
construction on any mismatch. A validator must return a mapping, never
`True`/`False`. The
receipt must contain exact matching `authority_commit` and
`implementation_sha256` plus an exact lowercase `source_closure_sha256`.
Missing or malformed `source_closure_sha256` refuses before physical work.
There is no unauthenticated `0`*64 source-closure fallback.

After a valid authority receipt, typed repository-input or origin-capture
failures publish one atomic terminal with
`physical_source_qualification_attempted=false`. A `Pro20OriginError` is an
invalid premise, not a source outcome. A `Pro20SourceFactoryError` from the
actual builder after `qualify_srcq1_rec1` begins is a typed source/factory
stop with `physical_source_qualification_attempted=true`. That flag means
the qualification invocation began; it does not claim source execution.
Unexpected programmer exceptions are not classified as science. Publication
errors retain the already tracked attempt flag and never force it true
merely because a terminal became visible. CLI reports
`physical_source_qualification_attempted` from that explicit flag only.

### Freeze B pins

1. Authority implementation A2 and `hlt17_srcq1_rec1.py` SHA-256.
2. Exact three-path Git delta from A2 (`evidence_io.InspectDelta`).
3. Final RSS/wall numbers, replacing or confirming the recommended table.
4. Exact evolution environment image (Python/NumPy/BLAS/platform).
5. Exact `source_closure_sha256` on the authority receipt. Matching it is not
   authentication of the physical source. The field is required; the CLI does
   not invent a placeholder digest.
6. Exact six-member `prospective_requested_cap_hex_by_member` identical to the
   already frozen SRCQ1 mapping. Do not tune a cap after seeing source/C1R1
   evidence. Subsequent overlay plans come only from HLT17 overlay objects.
7. Re-authenticated Attempt2 compact/raw hashes, both compact identities, and
   state-safe flags. Authority independently binds the sealed Git blob
   `786049…` and the one-field metadata extension `ba656a…`.
8. The `validate_authority_delta` mapping receipt, including matching
   `authority_commit`, `implementation_sha256`, `source_closure_sha256`, and
   the absent REC1 output namespace.

### Prelaunch and one real invocation, no informal retry

Before execution:

```text
make fgc-hlt17-srcq1-rec1-frz1
make verify-fgc-hlt17-srcq1-rec1-prelaunch
make status-fgc-hlt17-srcq1-rec1
```

The last command must return `safe_to_run=true` from the clean authority HEAD.

The single authorized invocation is routed by:

```text
make run-fgc-hlt17-srcq1-rec1
```

Its underlying command is:

```text
/opt/homebrew/Caskroom/miniconda/base/bin/python3 -I -B \
  scripts/qualify_fgc_hlt17_srcq1_rec1.py \
  --authority-commit <freeze-B-commit> \
  --output-directory <new-absent-directory>
```

Run it once. Do not remeasure RK members, retune caps, reduce grids,
substitute the SSPRK3-on-SBP4 hybrid, reuse a namespace, or retry
informally. A typed source/CFL/temporal/resource/inconclusive stop is the
result.

Recommended output namespace:

```text
runs/fgc-2-sf1/hlt17-srcq1-rec1-source-origin-qualification
```

## Typed recovery branches

| Branch | Independent requirement | Absorb overlay? | Continue members? |
|---|---|---|---|
| `prepared` | C1R1 eighteen-channel admission; `prepared_fresh`/`prepared_overlay`/`prepared_successor` | no | next remaining SSP member |
| `cfl` / `cfl_retry_required` | kind `cfl_retry`, finite positive observed/max hex, and `observed_ratio>maximum_ratio` | isolated checkpoint only | follow `next_plan` until that owner exhausts at count 33 or lattice exhaustion |
| `cfl` / `cfl_retry_exhausted` | same independent CFL predicate | isolated checkpoint only | stop |
| `source` / `source_retry_required` | wired `source_retry` evidence | isolated checkpoint only | follow `next_plan` until exhaustion |
| `source` / `source_retry_exhausted` | same source evidence | isolated checkpoint only | stop |
| `temporal` | HLT17 temporal disposition; not a source/CFL overlay | no source/CFL absorb | stop |
| `inconclusive` | missing CFL/source evidence, order-inconclusive channels, or uncertain publication | no | stop |
| `resource` | RSS/wall ceiling | no | stop |
| `invalid` | identity, isolation, cap, authority, nonzero `sink_writes`, owner count past 33 without exhaustion, or `observed_ratio<=maximum_ratio` | no | stop |

A runner label of `cfl_retry_required` without observed/max evidence is
`inconclusive`. It is not a CFL classification.

## Remaining independent PREF binder obligations

REC1 does not close these. They remain independent binder duties:

- Reconstruct numerical/state evidence without trusting runner labels or
  this recovery's outcome field.
- Independently re-hash implementation A, Attempt2 compact/raw terminals,
  static inputs, origin capture, source-configuration images, and the
  qualification payload.
- Independently recompute CFL from retained `observed_ratio_hex` and
  `maximum_ratio_hex`; do not accept a CFL stop from disposition text.
- Confirm RK members were referenced, not remeasured.
- Authenticate the physical source/import/environment image; SRCB1 and the
  factory still leave `physical_source_authenticated=false` locally.
- Environment/compiler gate as a PRO20 authority input, not this harness.
- Independent check that no endpoint was adopted and no campaign store was
  written.
- First-event live runner and Make target only after committed PRO20
  authority.
- Historical source leaves, Planck inputs, Attempt1/Attempt2 results, and
  pre-existing results remain intact and are not rewritten by REC1.

## Nonclaims

This harness does not establish GR-0 eligibility, calibration, SGB-L health,
DEF1, ROB1, a physical claim, or PROTO19 write authority. Status and ordinary
compact verification must not launch source evolution. Attempt2 is not
reclassified as a scientific or resource nonpass.
