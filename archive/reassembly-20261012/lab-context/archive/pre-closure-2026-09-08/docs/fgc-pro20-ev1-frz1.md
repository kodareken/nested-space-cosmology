# FGC-1-PRO20-EV1-FRZ1 — fresh first event prospective freeze

This document owns the corrected committed prospective first-event authority.
The earlier `7266e3a...` freeze stopped during prelaunch because a parent-only
test still required the config to be absent; it constructed no source or
namespace and authorizes no retry. This refreeze pins corrected implementation
parent `6e24358319940d31a4b914663d698e88e6c73905`, the
exact freeze delta, environment, source/origin/history inputs and operational
ceilings. The Make status/run routes are live only on that exact clean authority
commit. No event has been executed or result earned. The accepted PLAN.md and
research roadmap own the SF1 finish line.
The selected numerical object remains TDG11-IMP1: exact accumulation-corrected
comparison, unchanged p>=3/2 on all eighteen channels, retained full debit,
and only the original binary64 fine endpoint eligible for adoption.

## Fixed event and state boundary

- New protocol: `FGC-2-SF1-PROTO19`.
- New namespace: `runs/fgc-2-sf1/pro20-event1/calibration`.
- Only origin: synchronized HLT15 common event23 at `t=23/16`.
- Only event authorized by this PRO20 freeze: common event24 at
  `t=3/2`, one member at a time in the existing canonical order.
- RK4 sizes remain 2049/4097/8193; SSPRK3 sizes remain 4097/8193/16385, each
  on its own inherited spatial operator. The diagnostic SSPRK3-on-SBP4
  hybrid is not substituted.

RCV3, QA and MSEL endpoints are not possible origins. No old campaign store
is rewritten. Source/CFL/temporal retry caps remain32 and all applicable
stage/source/boundary controls remain required. A synchronized first event
does not supply GR-0 eligibility, the eight trapped events, SGB-L health,
DEF1, ROB1 or a physical claim. This freeze binds the applicable common-event
diagnostic owners; HLT16's synchronization helper alone is not a calibration
or observable certificate.

## Authenticated historical origin path

The source directory is `runs/fgc-2-sf1/proto17/calibration`. It now contains
the preserved HLT16 history in addition to the original eight HLT15 leaves.
Do not call the sealed PREF27 prelaunch-only exact-eight-leaf checker on the
expanded directory, and do not weaken that historical checker.

Read-only inspection established this route:

1. Authenticate the eight original leaves from the immutable PREF27 compact
   record and their exact raw hashes.
2. Authenticate the complete preserved expanded source tree with repairs
   disabled, then select `HLT16CampaignStore.authenticated_checkpoint_at_generation(1)`.
   This is the closed import bridge, not a later evolved checkpoint.
3. Read its six committed descriptor/payload pairs through `load_state`.
   Each has metadata generation0 and PROTO17 predecessor provenance.
4. Independently compare every original `physical_arrays` hash (u, p, q,
   tracer positions/proper times and both event histories), full original
   cursor, input hash, monitor, causal state and counters against the original
   HLT15 state object. All six comparisons passed in read-only inspection.
5. Bind the additional grid/label/static-template identity when the exact
   GR-0 shell is reconstructed. Do not invoke an old calibration runner chain.

The selected bridge checkpoint content ID is
`23f4a1ac64af5e1463862bf67ca240980921611a4a77d31050b0f92c81d44895`.
The complete proto17 source subtree has56 leaves and6054860 bytes, with
cleanup-manifest digest
`8e0b1f3ba70a4de7061d322ed041ff0328ae4e6ca7ed78da79421883f47ad777`.
That digest uses the `CLEANUP-BASELINE-TREE-v1` domain and repository-relative
paths; it is different from PREF27's eight-leaf projection digest.

| Origin identity | Content ID | Raw file SHA-256 |
|---|---|---|
| HLT15 checkpoint0 | `7c059dcc2197a9f57baf8012f4113c240917090171ee06cc639e9e8734603d2d` | `9bcbc68ca2f168d9f7ae87c37346d9ebf90aeaf2568f26546441634734511b3b` |
| HLT15 receipt | `78ae89a9f74193658003cf51f8567c4088db80fd6fb86867c2d581a0c16d5baf` | `7e99fb5fb70b6d2d590232d21222e9aa510b823dfdb3342aa86489b21d4d2f4a` |

PREF27's original eight-leaf projection digest is
`95a1bb96570216b8007e32a0cc05f4c31f420052f1d90c1b7134fa27ce958585`.
The tracked PREF27 compact file raw SHA-256 is
`102b385d058066d668e9e8541c974c736c583f6d5c4afaf9b99e95a16aabeb83`.
Capture authenticates those compact bytes, not only the eight-leaf tree hash.
The generation-zero receipt carries both `receipt_sha256` (its own content ID
`78ae89a9…`) and `checkpoint_sha256` (the referenced HLT15 checkpoint
`7c059dcc…`). Content IDs are path/schema-owned: a receipt is identified by
`receipt_sha256`, a checkpoint by `checkpoint_sha256`, an HLT16 descriptor by
`descriptor_sha256`, and an HLT15 state object by its raw canonical bytes.
Content IDs, raw-file hashes, cursor-chain IDs, full cursor-byte hashes,
physical u/p/q hashes and descriptor addresses are not interchangeable.

### Read-only capture layer

`src/recursive_horizons/fgc/evolution/pro20_origin.py` implements that
coordinator-verified route as a read-only historical-origin capture. It is
not a runner, not an execution freeze, not a compact certificate, and not
HLT17, source-runtime, resource, or environment qualification. Ordinary
focused tests in `tests/test_fgc_pro20_origin.py` use synthetic temporary
fixtures, mocks, and pure comparison helpers; they do not open the closed
historical store or launch a source/RHS.

Public capture pins:

- selected method `TDG11-IMP1`;
- the actual generation-1 bridge content ID above;
- original PREF27 compact/raw eight-leaf identities, including the distinct
  checkpoint/receipt content IDs and raw-file hashes;
- the complete safe `runs/fgc-2-sf1/proto17` inventory, leaf count, byte
  count, and `CLEANUP-BASELINE-TREE-v1` digest;
- canonical six-member order;
- common event 23 at `t=23/16` (`0x1.7000000000000p+0`);
- original physical arrays `u`, `p`, `q`, tracer positions/proper times, and
  both event histories as captured little-endian C bytes;
- the full original PROTO17 cursor bytes, input hash, monitor state, causal
  state, and inherited step/transaction/source/CFL counters.

Returned records are frozen. Derived labels are rebuilt from retained PREF27
compact/original checkpoint/state bytes and the actual generation-1
descriptor/payload bytes; a coordinated dataclass replacement plus recomputed
`captured_bytes` is not a verified origin. `revalidate_pro20_origin_capture`
performs that in-memory recheck without live files or success labels. The
capture retains HLT16 descriptor JSON, payload NPZ, and decoded
grid/label/physical arrays so a later factory can restore the bridge snapshot
without rereading a mutable store. Compact bytes and the complete source
inventory are compared before and after member reads. The caller-supplied
repository root must already be an absolute canonical path; a symlink root is
refused rather than silently resolved.

The capture does not build time-zero GR-0 shells, does not call PREF27's
exact-eight prelaunch checker against the expanded root, and does not repair,
publish, recover, or write the closed historical store. A later separately
qualified physical-source factory and new HLT17 rebase may consume this
slice; this module does not perform those steps.

### Inherited counters must survive

All members are at the same exact time, `0x1.7000000000000p+0`. The imported
TDG6 ledger is the HLT15 zero temporal ledger at that time, but the physical
step/transaction and source/CFL counters already contain history:

| Member | Step index | Transaction serial | Source retries | CFL retries |
|---|---:|---:|---:|---:|
| RK4-2049 | 342 | 1710 | 0 | 217 |
| RK4-4097 | 673 | 3365 | 0 | 443 |
| RK4-8193 | 976 | 4880 | 0 | 161 |
| SSPRK3-4097 | 631 | 2524 | 0 | 358 |
| SSPRK3-8193 | 987 | 3948 | 0 | 187 |
| SSPRK3-16385 | 1771 | 7084 | 0 | 0 |

The new IMP1 origin indices use those actual indices. Its frozen TDG6 sibling
must remain byte-identical. New per-macro retry counts start at zero; inherited
lifetime source/CFL counts do not. The new campaign identity is a separately
authenticated rebase, not permission to alter physical/tracer/history bytes.

The in-memory HLT17 member/codec and cursor/bridge joint checkpoint contract
is now consumed by a separately named production persistence sibling. That
sibling is not a boolean flip of the synthetic HLT17/PROTO19 store.

## Production persistence sibling

`evolution/pro20_ev1_protocol.py` and `evolution/pro20_ev1_store.py` own the
first production-shaped journal/checkpoint/store slice. They do not modify
`protocol_v19.py` or `hlt17_campaign_store.py`. Those synthetic owners still
require `store_kind=synthetic_temporary`, false authority flags, the
RK4-9/SSPRK3-9 cohort, and they still refuse any `/runs/` production path.

This sibling is a different owner because it:

- validates the exact six live members and their immutable C1R1 bundles;
- binds an externally supplied immutable authority receipt into a production
  root manifest and writer capability;
- publishes the entire `runs/fgc-2-sf1/pro20-event1` container atomically,
  with the seeded `calibration` store already complete inside it, under a
  caller-supplied canonical repository root;
- appends one validated member transition at a time, in canonical cohort
  order, with exclusive writer sessions;
- authenticates unique generations read-only;
- publishes `terminal.lock` for both a typed scientific terminal and a
  successful six-member `first_event_complete` event;
- refuses automatic takeover, stale-writer recovery, or in-place repair.

Flipping `production_write_authorized` on the synthetic schema would not
create those contracts. The synthetic module would still be a temporary
two-member qualification wire. The production module would still lack a
receipt type, live-member order, atomic namespace publication, and the extra
nonclaim flags.

### External root authority and generation-level nonclaims

This slice does not issue scientific or execution authority. It accepts only a
fully validated external receipt type
(`FGC-1-PRO20-EV1-authority-receipt-v1`). Construction of that receipt in
tests is not live authorization. The root manifest binds:

- protocol `FGC-2-SF1-PROTO19`;
- artifact `FGC-1-PRO20-EV1-FRZ1`;
- exact authority, implementation, config, source, origin, and environment
  SHA-256 identities;
- the production namespace;
- the six-member order;
- `t=23/16` to `t=3/2`;
- every nonclaim flag false.

Runner labels are not evidence. There is no `runner_sha256` identity in the
receipt or root. Generation journal, checkpoint, and manifest documents keep
every authority/claim boolean false even after the root has bound a receipt.
Calibration eligibility, SGB-L, DEF1, candidate, mechanism, and physics remain
false. A successful six-member event is terminally locked even when the inner
numerical disposition is `first_event_complete`.

The production protocol reuses PROTO19's already-tested transition reducers,
not its synthetic root or publisher. Every reused private helper is named in
`PINNED_PROTOCOL_V19_PRIVATE_DEPENDENCIES`; the later implementation/source
closure binds their bytes. Production parsing independently authenticates the
new root, journal, checkpoint and manifest schemas and refuses promoted
generation, session or terminal nonclaim fields.

### Publication and crash semantics

All publication is exclusive, no-replace, and no-follow through the neutral
evidence I/O helper. The seed store is one exclusive directory. Later
generations are exclusive directories inside that store. A free OS lock is not
authority to take over an unclosed session. Post-publication uncertainty stays
visible and poisonous on the writer that experienced it. Complete visible
bytes may be read as evidence; they do not authorize resume, repair, or
automatic recovery.

Focused tests are pure/synthetic and use canonical temporary repository roots
with the same existing `runs/fgc-2-sf1` ancestor. They build the six-member
production seed once, use one bounded RK4-2049 transition, and test completion
at the pure cohort/terminal-lock layer; they do not numerically walk six
production members to the target. They do not open the real
`runs/fgc-2-sf1/pro20-event1/calibration` tree, add an authority TOML,
status/run CLI, live target, compact result, or physical-source call.

## One-shot runner and authority implementation parent

`pro20_ev1_runtime.py` owns the one-shot state machine, while
`pro20_ev1_authority.py` owns the prospective Git/environment/resource gate.
The thin `run_fgc_pro20_ev1.py` and freeze-config reproducer are now bound by
the committed config and exact-delta authority. Their default/status paths
construct no origin or source and publish no namespace. The live Make run
target first requires the full prelaunch/status gate and cannot accept a force
or receipt bypass.

The runner:

- validates authority before origin or source construction;
- binds the captured origin, selected generation-one bridge, source closure
  and all six source-configuration identities before seed publication;
- derives each fresh cap as
  `min(remaining, cfl_maximum*grid_spacing/previous_speed_upper)` and uses only
  cursor-owned plans on pending retries;
- intercepts prepare- and admission-time temporal exhaustion, publishing the
  unchanged accepted bundle with non-executable ledger evidence;
- publishes post-seed wall/RSS/disk/namespace/attempt stops as typed resource
  terminals and closes the writer;
- leaves an unclosed forensic session on unexpected or post-publication
  uncertainty; and
- has no receipt/force bypass or automatic resume/takeover path.

Operational ceilings are 86400 seconds total wall, 4 GiB RSS, 16 GiB namespace
bytes, 32 GiB minimum free disk and 4096 published attempts. The approximate
224 accepted-step count derived from the initial caps is planning evidence,
not an admission or scientific gate.

The authority implementation authenticates a clean exact one-parent freeze,
the complete declared freeze delta, implementation/runner/authority/config
hashes, source closure and environment, REC1-PREF1's separate-authority
license, the selected HLT15 generation-one capture, static inputs, Planck and
historical-runs baselines, process isolation, resources and absent output or
staging paths. It reads the historical origin but never constructs the
physical source. The later freeze must add the config and truthful Make,
catalog, checker and public-frontier routing atomically.

## Remaining implementation and freeze obligations

The origin capture is only verified input. The persistence sibling is only a
store/protocol owner. Neither authorizes PRO20 state advancement, a
first-event run, or a physical result.

`PRO20EV1CampaignStore.publish_seed_store` is the authority's mechanical
publication primitive, not a supported direct live entrypoint. The catalog
names `run-fgc-pro20-ev1` as the sole authorized state-changing route; status
and ordinary verification remain read-only.

The physical-source factory must use the preserved GR-0 library construction
and exact static inputs, not new imports of `run_fgc_gr0_calibration_v*`.
Authenticate the complete source/configuration closure and the pinned
evolution environment independently of a synthetic callable fingerprint.
Qualify actual source/operator integration and all production sizes under
prospectively declared resource limits. Synthetic exponential size timings
are not physical-RHS or total campaign wall-clock evidence.

Only after the coherent source image, full HLT17 prelaunch audit and exact
committed-image status return `safe_to_run=true` may the explicit first-event
run begin. This freeze leaves these seams explicit:

- one authorized live first-event run;
- independent PREF1 binder.

Status and ordinary compact verification must not launch source evolution.
The resulting independent PREF1 must reconstruct the numerical/state evidence
without trusting runner classifications or publication decisions. Every
historical source leaf, Planck input and pre-existing result remains intact.
