# FGC-1-HLT17-SRCQ1 — physical-source/origin/resource qualification

This is implementation A of the prospective SRCQ1 harness. It is not a
completed certificate, not permission to run the physical six-member
measurement, and not a campaign or compact result. The accepted
[PLAN](../PLAN.md) owns the SF1 finish line.
[PRO20-EV1](fgc-pro20-ev1-frz1.md), [HLT17-MON17](fgc-hlt17-mon17.md),
[SRCB1](fgc-pro20-source-binding.md), and [SRCF1](fgc-pro20-source-factory.md)
remain the retained origin, binding, factory, and C1R1/HLT17 contracts.
AGENTS.md stays unchanged creative context.

A later exact-delta authority freeze B must pin this implementation before
the one real measurement. Freeze B, not this document, authorizes that run.

## Frozen authority B

The implementation image is commit
`7123aa385b3f4165186acd241f532b58ebbb3a9b`; harness SHA-256 is
`e801532e01061122b03bd2e5b8b1359926a65445b676d695fcdc6da371784bab`.
The complete 244-file Recursive Horizons source/runner closure plus pinned
Python/NumPy binaries has SHA-256
`616fc0c740d6b00cf2983c4e19b576ad25441904e81a2853229f59055d57f3ab`.
The only authority-delta paths are this document and
`configs/fgc/fgc-1-hlt17-srcq1-frz1.toml`.

The freeze pins the retained origin capture `2824af06...`, exact CPython3.14.3
/NumPy2.5.1/Accelerate/Darwin-arm64 image,4GiB RSS,900s/member,3600s total,
and all six prospective initial caps. The output namespace is
`runs/fgc-2-sf1/hlt17-srcq1-source-origin-qualification`. It must be absent.

## Owned implementation

| Path | Role |
|---|---|
| `src/recursive_horizons/fgc/evolution/hlt17_srcq1.py` | Deterministic no-write qualification core |
| `scripts/qualify_fgc_hlt17_srcq1.py` | Guarded CLI; default plan/status only |
| `tests/test_fgc_hlt17_srcq1.py` | Small synthetic/mocked objects and temporary output |
| this document | Prospective map, resource recommendations, freeze-B seam |

The core revalidates a retained origin capture; reads the four static factory
inputs through `evidence_io.read_regular_file` (no-follow) and their fixed
hashes; builds six fresh static GR-0 shells only when the real factory path
is invoked; restores the generation-1 origin; checks the six method-owned
grids, operators, inherited counters, and TDG6 sibling; evaluates each bound
source once at the physical origin with source/raw/kinetic/causal/constraint
diagnostics and no input/monitor/causal/tracer mutation; then prepares fixed
C1R1 HLT17 shadows at a prospectively extracted cap without commit, overlay
absorption, or endpoint adoption. Members run one at a time. Typed
source/CFL/temporal/resource/inconclusive/invalid outcomes are distinct.

There is no store or campaign publication, no accepted-state advance, and no
corrected or shadow endpoint serialization. Atomic terminal writing of the
qualification record may use neutral `evidence_io`. An uncertain publication
is a typed stop. Old `run_fgc_gr0_calibration_v*` runners are not imported.

## API and state contract

Public core APIs:

- `plan_status()` — identities, recommended ceilings, seam text; no I/O
- `derive_prospective_requested_caps(origin,runtime,expected_cap_hex_by_member)`
  — independently recomputes each initial request as
  `cfl_maximum*grid_spacing/inherited_previous_speed_upper`. HLT15 is
  `FRESH_READY` and contains no plan, so no plan is fabricated. The authority
  mapping must match this derivation; no remaining-interval fallback or
  post-result retuning is accepted
- `read_static_factory_inputs` / `static_input_identities`
- `evaluate_bound_source_once`
- `qualify_physical_source_origin` — no-write six-member walk
- `qualify_and_publish` — authority seam, then exclusive directory publish
- `require_authority_delta` / `load_coordinator_authority_delta_validator`

Result identities always include implementation, environment, static-input,
origin, and source-configuration hashes; wall/RSS facts; all eighteen
channel decisions/debits when a C1R1 family exists; C1R1 actual identity
versus the IMP1 reference wire; and byte fingerprints before/after each
member. C1R1 is
`tdg11_c1r1_integer_exponent_direct_ring_v1`. IMP1 remains
`FGC-1-TDG11-IMP1` / `tdg11_imp1_exact_bernstein_with_dual_fallback_v1`.
They are not the same execution.

## Resource ceilings — coordinator freeze-B values

These are resource guards from already preserved provisional controls. They
are not scientific thresholds and were not fitted to a measurement.

Already frozen and not retunable here:

| Control | Value | Source |
|---|---|---|
| Source retry cap | 32 | HLT17 / PRO20-EV1 |
| CFL retry cap | 32 | HLT17 / PRO20-EV1 |
| Temporal retries per macro step | 32 | TDG6 |
| Per-member enclosure | production N and 8N/16N | HLT17 production table |

**Recommended conservative process ceilings for freeze B to pin or replace.**
They are active fail-closed defaults in implementation A so an unbounded
run is not implied, but they are not a measured physical-RHS budget.

| Ceiling | Recommended default | Rationale |
|---|---:|---|
| Process RSS | 4 GiB | Preserved baseline full-size controls peaked near580MiB; this keeps substantial non-scientific process margin. |
| Per-member wall | 900 s | Baseline largest-member exact assessment was about229s; source work and C1R1 retain conservative margin. |
| Total wall | 3600 s | Six members remain bounded to one hour total; not a scientific threshold. |

`coordinator_freeze_required=true` and `scientific_threshold_fitted=false`
are recorded on every result. **Final numbers are a freeze-B decision.**

## Authority seam and one real run

Implementation A does not accept a boolean skip/force flag. Real CLI
execution requires both of:

```text
--authority-commit <exact 40- or 64-digit lowercase commit>
--output-directory <absent namespace>
```

Default CLI (no arguments) prints `plan_status()` and does not capture the
historical origin, call the physical source, or write. Incomplete real-run
arguments refuse. An existing output namespace refuses. A missing parent
directory is a typed absent stop, not a mkdir workaround.

Authority/delta validation is coordinator-owned freeze B:

```text
recursive_horizons.fgc.evolution.hlt17_srcq1_auth1.validate_authority_delta
```

That module is not part of implementation A. Until freeze B provides it, the
real CLI path stops at the seam before origin capture and factory
construction. A validator must return a mapping, never `True`/`False`.

### Freeze B must pin

1. Implementation A commit and `hlt17_srcq1.py` SHA-256.
2. Exact git delta from A to the freeze-B commit (`evidence_io.InspectDelta`).
3. Final RSS/wall numbers, replacing or confirming the recommended table.
4. Exact evolution environment image (Python/NumPy/BLAS/platform).
5. Source-closure digest that the factory stores as an unauthenticated
   external reference; matching it is not authentication. The freeze-B
   validator receipt may supply `source_closure_sha256`; otherwise the CLI
   records the explicit unauthenticated placeholder `0`*64.
6. Exact six-member `prospective_requested_cap_hex_by_member`. Each value is
   recomputed from retained causal `previous_speed_upper`, frozen
   `cfl_maximum=1/8`, and method-owned grid spacing. The HLT15 cursor contains
   no plan and remains unmodified. TDG7 records the separately selected
   executable width at or below the request. Do not tune a cap after seeing
   source/C1R1 evidence.
7. The `validate_authority_delta` mapping receipt.

### One real invocation, no informal retry

After freeze B is committed and `hlt17_srcq1_auth1.py` exists:

```text
/opt/homebrew/Caskroom/miniconda/base/bin/python3 -I -B \
  scripts/qualify_fgc_hlt17_srcq1.py \
  --authority-commit <freeze-B-commit> \
  --output-directory <new-absent-directory>
```

Run it once. Do not retune caps, reduce grids, substitute the SSPRK3-on-SBP4
hybrid, rename IMP1 execution, reuse a namespace, or retry informally. A
typed source/CFL/temporal/resource/inconclusive stop is the result.

## Remaining independent PREF binder obligations

SRCQ1 does not close these. They remain independent binder duties:

- Reconstruct numerical/state evidence without trusting runner labels or
  this qualification's outcome field (`FGC-1-PRO20-EV1-PREF1` and any later
  SRCQ1-PREF).
- Independently re-hash implementation A, static inputs, origin capture,
  source-configuration images, and the qualification payload.
- Authenticate the physical source/import/environment image; SRCB1 and the
  factory still leave `physical_source_authenticated=false` locally.
- Environment/compiler gate as a PRO20 authority input, not this harness.
- Independent check that no endpoint was adopted and no campaign store was
  written.
- First-event live runner and Make target only after committed PRO20
  authority.
- Historical source leaves, Planck inputs, and pre-existing results remain
  intact and are not rewritten by SRCQ1.

## Nonclaims

This harness does not establish GR-0 eligibility, calibration, SGB-L health,
DEF1, ROB1, a physical claim, or PROTO19 write authority. Status and ordinary
compact verification must not launch source evolution.

## FGC-1-HLT17-SRCQ1-ATTEMPT1 terminal status

The one run under authority `271ab71b...` stopped during implementation
validation before the first member qualification. The harness dereferenced the
C1R1 ledger's encoded `inherited_snapshot` string as though it were a TDG6
object and raised `AttributeError`. Static-factory construction may already
have called its existing time-zero source checks; no accepted state advanced,
no endpoint or campaign store was written, and no qualification terminal was
published. The output and PRO20 namespaces remain absent.

This is `invalid_implementation_no_terminal_no_state_advance`, not a source,
resource, temporal, model or physics nonpass. The authority is consumed and
does not authorize retry. A corrected implementation and new prospective
freeze are required before another measurement.

Implementation `7123aa3` decodes the retained TDG6 sibling through its owning
C1R1/IMP1 ledger codec and adds a regression using the real encoded snapshot
type. The new freeze above supersedes the consumed Attempt1 authority; it does
not reinterpret or erase Attempt1.

### FGC-1-HLT17-SRCQ1-ATTEMPT2 terminal status

The corrected run under authority `f3a31dc...` atomically published one
45,589-byte terminal (`889a8d17...`; directory payload `7b06deb6...`). The
three RK4 members completed source evaluation and C1R1 preparation with all18
channels reporting `resolved_order_pass`. SSPRK3-4097 completed its origin
source evaluation and the runner reported `cfl_retry_required`; accepted state
and endpoints remained unchanged. Total wall was206.55s and peak RSS981.7MB,
within the frozen resource envelope.

The raw terminal omitted the stop's observed/max CFL-ratio evidence. Therefore
Attempt2 is bound as
`completed_terminal_reports_cfl_evidence_incomplete_no_state_advance`, not an
independently reproduced CFL classification or scientific/resource nonpass.
The authority is consumed and authorizes no retry. Any recovery must retain
full source/CFL evidence and receive a new prospective freeze.
