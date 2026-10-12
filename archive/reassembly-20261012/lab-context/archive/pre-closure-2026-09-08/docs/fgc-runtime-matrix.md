# FGC runtime matrix — arithmetic provenance, not physics

The repository contains two immutable numerical histories that were generated
with NumPy `2.5.1` but different Python/BLAS builds. Their last bits are not
interchangeable:

| Partition | Frozen runtime | Ownership |
|---|---|---|
| Analytic, pre-RSP2, and TDG4 theorem certificate chain | CPython `3.14.6`, OpenBLAS `0.3.33`, Darwin `arm64` | Exact regeneration of the established analytic/calibration certificates, TDG4's no-history rational freeze controls, and its post-freeze sampling-identifiability theorem binder |
| RSP2 raw endpoint, PROTO13 restart, CAL10/PREF15 reconstruction, and TDG1/PREF16, TDG2/PREF17, plus TDG3/PREF18 actual-history diagnosis | CPython `3.14.3`, Accelerate, Darwin `arm64` | Exact replay of the immutable RSP2/PREF14 endpoint, every continuation from it, the terminal common-event assessment, and all three read-only six-history-ladder diagnoses |
| PDF rendering | CPython `3.12.13`, ReportLab `4.4.9` | Deterministic document dependency only |

This split does not change an equation, threshold, state, or conclusion. It
prevents one arithmetic backend from being treated as if it had generated the
bytes owned by another. The PROTO13 runner records its complete observed
Python, NumPy, BLAS, and platform metadata and refuses drift before evaluating
the source or advancing state.

The prospective RSRC1 per-member seed/attempt child belongs to the evolution partition.
Its focused controls are synthetic/store-blind, but any later physical child
must use the exact evolution runtime and a separately frozen environment. The
child imports no PRO20 production store; a parent alone may eventually own a
new store. No live target or authority currently exists.

`scripts/run_fgc_test_partition.py` routes each test module to the runtime that
owns its immutable floating evidence. Exact and synthetic tests that do not
consume the RSP2 raw endpoint remain in the analytic partition. PREF14 raw
replay, HLT11 reproduction, the PROTO13 restart tests, CAL10/PREF15's
checkpoint recomputation, and TDG1/PREF16's terminal-history binder use the
evolution partition. A partially present raw bundle is an error.

Process-global settings have their own test boundary. The sealed UHYP1/MD1
reproducers disable Python's integer-to-decimal digit cap when imported, while
the unmodified TDG11-IMP1 enclosure tests explicitly prove operation under a
finite cap. The partition runner executes that exact test module in a fresh
interpreter from the same Python installation with `int_max_str_digits=4300`.
It passes every selected test ID, preserves child failure output and exit
status, and reports in-process and isolated case counts separately. No test
is skipped or weakened, no historical source/result hash is changed, and the
legacy process retains its own setting. This is test-process isolation, not a
change to scientific arithmetic, IMP1's rational resource limits or a campaign
environment.

The partial SGB1-CTL1 annular/source algebra, source and whole-cell admission,
initial-health, cone/continuity, matched-family principal, point-local
symmetrizer, bounded synthetic runtime, and finite no-trap proof instruments
belong to the analytic partition. The prospective SOL1 orbit-local
Picard--Lindelof continuation is also analytic: its frozen nominal result is
interval-inconclusive before the first tube and opens no health or execution.
SOL1-FRZ1 adds compact Git/source/raw/`runs/`/reconstruction-blind verification;
independent SOL1-PREF1 binds the same analytic scoped outcome, and its ordinary
compact verification is likewise Git/source/raw/`runs`/reconstruction blind.
DEF1-STAB1 error/margin, local geometric
sensitivity, exact conversion/provider qualification, and its nonpromoting
freeze-contract payload belong there as well. They use exact
rationals/intervals, explicit outward binary64 storage and small declared
synthetic arrays, consume no ignored campaign arrays, and have no live
runner. DEF1-STAB1-FRZ1 adds a Git/raw/`runs/`/source/trajectory-blind compact
instrument freeze; ordinary verification is `make fgc-def1-stab1-frz1`.
Independent STAB1-PREF1 is the same analytic partition: Git/raw/`runs/`/
source/trajectory-blind compact verification is `make fgc-def1-stab1-pref1`.
`verify-fgc-wave0-algebra` tests the bounded instruments. These routes do
not establish aggregate SGB-L health, authenticate trajectory-error radii,
or complete holdout. PREF1 sets `DEF1_error_map_passed` only as
`map_readiness_only`.

FGC-1-TDG1-FRZ1's synthetic counterexample and checkpoint-hash-only freeze
stay in the analytic partition because they load no raw history array;
FGC-1-TDG1-PREF16's now-bound actual-history diagnosis belongs to the
evolution partition. FGC-1-TDG2-FRZ1 likewise stays in the analytic partition:
its Decimal/binary64 and interpolation controls are synthetic, and it verifies
only the immutable checkpoint hash without loading a history. The later TDG2
actual-history binder must belong to the evolution partition.
FGC-1-TDG2-PREF17 now supplies that binder: it loads the immutable Accelerate-
owned arrays read-only, advances no state, and remains in the evolution
partition.
FGC-1-TDG3-FRZ1 is again analytic: it reads only PREF17's compact tracked
diagnosis, may hash but never loads the raw checkpoint, and exercises its two
native-grid estimators only on synthetic histories. FGC-1-TDG3-PREF18 is the
evolution-owned actual-history binder: it restores the immutable arrays
read-only, applies both frozen estimators to all `576` ladders, advances no
state, and defines no replacement gate. FGC-1-TDG4-FRZ1 returns to the
analytic partition: it reads only PREF18's compact tracked record and checks
exact rational row-space, smooth support-geometry, uniform-alias, and
regularity-bound controls. It loads no campaign array and does not define the
replacement gate whose premises it is designed to constrain.
FGC-1-TDG4-PREF19 remains analytic. It independently binds the frozen theorem,
the immutable PROTO4/PROTO13 temporal rule, and CAL10's compact stop record
without reading a history array. It retires that sampled rule only as a
continuum admission and authorizes prospective replacement-gate design; it
does not define the replacement or authorize any trajectory.
FGC-1-TDG5-FRZ1 and FGC-1-TDG5-PREF20 remain analytic. They read only
compact tracked results and immutable source blobs and execute exact rational
cubic-Hermite, restriction, extremum, Butcher-order, conditional Richardson,
engine-shape, and mutation controls. Neither reads a raw history or checkpoint,
allocates a production state, or performs a trajectory step. PREF20 authorizes
implementation and synthetic validation only.

FGC-1-TDG5-IMP1 is the first NumPy-backed TDG5 runtime instrument. It creates
isolated synthetic GR-0 shadow transactions for both methods, advances one
coarse full step and two fine half steps through the unchanged engine, and
commits only the fine endpoint after revalidating the real state and ledgers.
It also evaluates complete owned-state cubic differences with separated
binary64 and Bernstein certification debits. This is synthetic implementation
evidence, not the production campaign compositor: it freezes no threshold,
defines no replacement temporal admission, advances no production state, and
opens no PROTO14 or candidate branch.

FGC-1-TDG6-FRZ1 remains in the analytic partition. Its design module and
reproducer use exact rational arithmetic to freeze the three-level
`p >= 3/2` contraction test, method stage accounting, channel classes, retry
limits, and public debit semantics. It reads only compact tracked TDG5 evidence
and immutable source blobs. It loads no raw campaign bundle, allocates no
production state, advances no trajectory, and authorizes no production
compositor. The later TDG6 binder and runtime implementation remain separate
artifacts.

FGC-1-TDG6-PREF21 also remains in the analytic partition. Its independent
standard-library theorem module derives direct and composed quarter maps,
exact interval-order and channel/debit controls, and parses immutable runtime
source syntax. It loads no NumPy campaign arrays or checkpoint and advances no
state. The result authorizes only implementation and synthetic qualification
of the still-absent production compositor; a future runtime implementation
will belong to the pinned evolution/restart partition.

FGC-1-TDG6-IMP2 remains in the analytic partition because its NumPy runtime
and certificate execute only nine-node synthetic controls: the seven-path
RK4/SSPRK3 compositor, normal-flow tracer previews, complete continuous channel
admission, fine-only atomic commit, temporal retry, and checkpoint extension.
It reads no campaign trajectory or checkpoint and advances no production
state. A future PROTO14 integration and trajectory remain separately frozen
evolution/restart-runtime actions.

FGC-1-PRO14-FRZ1 is split at its actual arithmetic boundary. The pure
`protocol_v14` validator and its mutation tests remain analytic because they
load no campaign array. The canonical freeze reproducer and
`test_fgc_pro14_frz1_reproduction` belong to the evolution/restart partition:
they restore the six ignored Accelerate-owned payloads to verify their exact
state, tracer, history, monitor, causal, and counter hashes, although they
advance no state. The compact tracked certificate, namespace absence, exact
restart manifest, and no-trajectory boundary remain portable audit evidence
when the complete ignored bundle is absent from a clean clone. Future
HLT12/MON12 runtime preflight and every PROTO14 trajectory step also belong to
the pinned evolution/restart partition.

FGC-1-HLT12-MON12 belongs to the evolution/restart partition. Its canonical
preflight restores the six Accelerate-owned `t=23/16` payloads, binds the
PROTO14 successor member and runner, initializes the six zero TDG6 ledgers,
and attacks checkpoint, retry-ownership, namespace, and nonpromotion
boundaries without advancing a production state. The small synthetic
`test_fgc_proto14_runtime` module loads no raw checkpoint and remains eligible
for the analytic partition; the HLT12 reproducer and PROTO14 runner tests are
runtime-bound. The later production calibration is deliberately absent from
ordinary verification and may run only from the immutable HLT12 authorization
commit under this same pinned Accelerate runtime.

FGC-1-CAL11-PREF22 also belongs to the evolution/restart partition. Its
canonical terminal binder consumes the all-or-none PROTO14 raw bundle only to
replay the checkpoint-owned result: common event `23` is complete, while event
`24` aborts at the pre-shadow runtime validator. The binder confirms zero
accepted post-restart macro steps for all six members and rolls the terminal
checkpoint back to event `23`; it creates no output namespace and runs no
replacement trajectory. The compact certificate remains portable when the
whole raw bundle is absent, while any partial raw bundle fails closed.

Once that complete PROTO14 bundle exists, the evolution partition first runs
CAL11/PREF22 in `--check` mode. Only after that independent terminal binding
passes does the partition seal the three known prelaunch-only test IDs listed
in the runtime matrix. It asserts that every listed ID exists. The contract
proves that no other test was filtered. This preserves historical
authorization-time evidence
rather than asking prelaunch artifacts to become postlaunch reproducers. All
other evolution tests, including CAL11's raw-terminal binder, remain selected.
The analytic partition remains unchanged and never reads the PROTO14 raw
bundle.

FGC-1-TDG7-FRZ1 belongs to the analytic partition. Its prospective **stage-safe
exact-binary64 shared lattice** consumes compact CAL11 lineage only: `G=8Q`
makes the actual RK4/SSPRK3 `c={0,1/2,1}` stage times and endpoints exact. The
historical independently rounded one/two/four boundary arrays are retained as
evidence of the failure and differ from the selected `8Q` boundaries. The
diagnosis neither loads the raw bundle nor invokes the campaign runtime,
mutates PROTO14, or resumes its checkpoint. The independent TDG7 binder is
authorized but incomplete; runtime repair and fresh calibration remain outside
this analytic freeze.

FGC-1-TDG7-PREF23 also belongs to the analytic partition. It independently
binds the immutable TDG7/CAL11 lineage and historical numerical-engine source,
reconstructs the `G=8Q` witness, and checks the actual RK4/SSPRK3 stage
abscissae without loading a raw campaign bundle or invoking a runner. It
authorizes only later runtime-repair implementation plus synthetic
qualification. The repair, PROTO14 mutation/resume, a new protocol/namespace,
fresh calibration, candidate execution, and physical claims remain outside the
analytic binder.

FGC-1-TDG7-IMP3 belongs to the synthetic runtime-instrument partition. It
implements the authorized shared-lattice adapter and validates its RK4/SSPRK3
paths, typed pre-work stops, fine-only commit, and TDG6 retry lineage without
loading CAL11/PROTO14 arrays or creating a runs namespace. The normal minimum
step is TDG6-owned; TDG7's corresponding stop is defensive-only for forged
evidence. A stateless adapter cannot own the discarded-retry distinction, so
only a future protocol-freeze design with a durable cursor is authorized.

FGC-1-HLT13-MON13 remains in the synthetic runtime-instrument partition: its
cursor/journal/checkpoint qualification consumes compact lineage and temporary
synthetic stores only. PROTO16/PRO16-FRZ1 remain in the analytic compact
partition: they freeze the missing successful `COMMON_EVENT_COMMIT` contract
and future external GenesisSpec authority, but implement neither transition nor
production runner and create neither PROTO16 namespace. Their later HLT14
implementation/preflight cannot begin: PREF24 statically finds that PROTO16
does not freeze complete successor cursor/high-water and external
generation-zero construction rules to the deterministic standard. At that
PREF24 boundary HLT14 was unauthorized; only a prospective PROTO17 freeze
could repair the contract.

PROTO17/PRO17-FRZ1 remains in the analytic compact partition as well. It adds
only a standard-library pure reference construction for the predecessor-only
common-event and external-genesis mapping; it does not use raw state arrays,
write a namespace, or qualify a production runner. At that freeze boundary,
HLT14 remained a later runtime/synthetic-qualification task; the completed
HLT14 result is recorded separately below.

PRO17-PREF25 is also analytic and compact. Its separate theorem may read the
sealed PROTO17 blobs from Git and construct synthetic objects in memory, but it
imports neither PROTO17 implementation module and does not inspect `runs/`.
Its immutable, schema, hash, recovery, and mutation checks authorized only the
later HLT14 implementation plus synthetic qualification. The four raw-input,
launch-authentication, payload-replay, and namespace/foreign-store requirements
remained uncompleted and are still open after the synthetic HLT14 result.

**FGC-1-HLT14-MON14** is now complete in the synthetic runtime-instrument
partition. It uses temporary production-shaped Git authority, real-NPZ bundle,
and store fixtures to reproduce PREF25's three sealed hashes and qualify the
generation-zero, common-event, recovery, semantic, reuse, and foreign-store
implementation controls across `72=4+6+18+44` unique cases. It is not a
production launch: temporary raw hashes depend on the fixture inputs, and its
exception-injection coverage is not a power-loss durability proof. The actual
raw bundles, launch-time Git/blob/live-import authority, historical-payload
replay, and production namespace reuse/foreign-store rejection were
`required=true, completed=false` at the HLT14 boundary. The only successor authorization was separate
production-prelaunch design; no namespace, pretrajectory, calibration,
eligibility, SGB-L, FGC-QR, DEF1, EFT, transition, or physical result follows.

**FGC-1-PRO18-PREF26** is a separate, deliberately raw-dependent local binder,
not part of the portable aggregate. Its focused target reads only the two
declared PROTO12/RSP2 checkpoints and journals, checks exact byte/hash and
bounded archive identities, reconstructs six restart payloads, and replays
their legacy semantics. It completes only HLT14 duties one and three. It does
not inspect the future PROTO17 roots or run the HLT14 runtime. The portable
default verifies the compact PREF26 artifact statically; deliberate raw
reproduction remains behind `make verify-fgc-pro18-pref26`.

**FGC-1-PRO18-AUTH1** is content construction in the evolution runtime
partition. It constructs a canonical commit-independent **calibration-only**
authority payload binding compact PREF26 evidence, the complete
GR-0/amplitude-3 `GenesisSpec`, the production source closure including the
standard-library-only HLT15 bootstrap/authority module, the exact evolution
environment, raw-source configuration, and a deterministic fresh campaign ID.
At that historical AUTH1 content boundary, it did not authenticate an external
commit/Git/blob or live imports, touch a future root, authorize or execute
HLT15, materialize a namespace/generation zero, or advance a state.

**FGC-1-HLT15-GEN1** is a one-time, same-invocation evolution-runtime launch
operation, not a continuation. It performs the externally supplied AUTH1
tuple, committed-byte, source/import, environment, and raw-reimport checks
before the first calibration-root observation, then atomically materializes the
six-member generation-zero store. It advances no state. **FGC-1-PRO18-PREF27**
is the sealed read-only binder for that eight-leaf store and its temporary
no-adoption controls. Actual admission rejects stale sibling HLT15 staging
residue; the binder records three corresponding no-staging-residue controls.
A complete ignored store may be absent in a clean clone; a partial presence
fails closed. Ordinary verification consumes compact
PREF27 evidence and performs no root re-observation, raw reimport, GEN0
materialization, or temporary race/reuse control. The explicit focused target
alone reruns those bounded read-only/temporary checks.

Complete absence in a clean clone skips only the raw-bound
partition; `scripts/check_repo.py` still verifies the tracked compact
certificates and hashes. The analytic OpenBLAS runtime may inspect the compact
PREF15 record, but it must not demand bitwise regeneration of the Accelerate-
owned common-event reductions.

The executable locations in the Makefile are local defaults and may be
overridden with `PYTHON`, `EVOLUTION_PYTHON`, and `PAPER_PYTHON`. An override
does not relax the fingerprint: a mismatching runtime fails before tests or
evolution.

GitHub runs `make verify-portable-clean-clone` on Linux. That gate installs the
pinned NumPy and ReportLab dependencies, exercises the clean-clone test surface,
verifies all compact certificates and public data, and rebuilds the paper. It
does not claim byte reproduction of either local macOS arithmetic history. The
ignored RSP2/PROTO13 raw bundles are absent there, so the runtime-bound modules
are skipped and remain covered by the hash-bound local canonical gate.

**FGC-1-PRO18-FRZ1** is a static overlay, not an additional runtime partition.
It targets unchanged **FGC-2-SF1-PROTO17** and freezes the real
two-source-container/fresh-output prelaunch chain
`PREF26 -> FGC-1-PRO18-AUTH1 immutable authority -> HLT15-GEN1 -> PREF27`.
It opens no raw archive/history and creates no namespace. At the PRO18 freeze
boundary, the four HLT14 production duties were `required=true,
completed=false`; only PREF26 implementation/design was authorized. No pretrajectory, calibration, GR-0
eligibility, SGB-L, FGC-QR, DEF1, EFT, transition, or physical result follows.

**FGC-1-PRO18-AUTH1** now constructs the canonical commit-independent
**calibration-only** authority payload in the evolution runtime partition. It
binds compact PREF26 evidence, the complete GR-0/amplitude-3 `GenesisSpec`,
production source closure including the standard-library-only HLT15
bootstrap/authority module, exact evolution environment, raw-source
configuration, and deterministic fresh campaign ID. At that historical
content boundary it did not authenticate external commit/Git/blob or live
imports, authorize or execute HLT15, observe a future root, materialize
namespace/generation zero, or advance state. The completed HLT15/PREF27 status
is recorded above; neither artifact advances the state.

**FGC-1-HLT16-MON16 / FGC-2-SF1-PROTO18** is qualified in the pinned
Accelerate evolution/restart partition. It covers finite little-endian payload
encoding, descriptor/journal/checkpoint publication, exclusive bootstrap and
writer ownership, source/CFL/TDG6 retry separation, deterministic recovery,
status, and the bounded first-event runner. The first authenticated attempt
then completed five finite accepted `RK4-2049` macro-steps and stopped invalid
before publishing a CFL-retry record: the coordinator serialized the outer cap
instead of TDG7's executed lattice-safe width. **FGC-1-PRO19-CFL1-FRZ1**
freezes one exact continuation from clean nonterminal generation `6`,
checkpoint `6277b3be...`; the corrected image and unchanged store still require
post-commit authentication. No common event, GR-0 calibration, candidate
execution, activation, trapping, DEF1, or physical claim follows.

**FGC-1-PRO19-SID1-FRZ1** is a read-only exact pre-recovery authority in the
same evolution/restart partition. It authenticates the real generation-`7`
checkpoint plus its uncheckpointed sequence-`7` TDG6 rejection, recursively
decodes the stored TDG6 evidence in memory, and freezes one deterministic
sequence-`8` cursor/ledger transition and complete generation-`8` checkpoint.
It takes no writer lease, publishes no recovery, reruns no rejected proposal,
and leaves the accepted finite descriptor/payload/`u,p,q` state unchanged. It
therefore supplies no common event, calibration, candidate, activation,
trapping, DEF1, retained-EFT, transition, or physical result.

**FGC-1-PRO19-SID2-PREF1** also belongs to the pinned evolution/restart
partition because its one-time live audit reopens the persisted NumPy payload.
It independently binds the exact sequence-`8` cursor transition and complete
generation-`8` checkpoint published after immutable SID1 commit `c83ef958...`,
then proves that the six descriptors and affected `u,p,q` arrays did not
change. The explicit SID2 post-recovery target owns that temporal live-store
observation. Ordinary aggregate verification consumes only the compact SID2
record after the mutable campaign advances and must not reinterpret lawful
later leaves as drift in the historical generation-`8` boundary. SID2 is a
recovery-integrity binder, not trajectory-resume or physics authority.

**FGC-1-PRO19-SID3-AUTH1** is the historical committed static execution
authority consumed by the now-terminal GR-0 resume attempt. It binds the exact
committed execution closure owned by immutable implementation commit
`4d7133c...` and the exact compact generation-`8` inputs without reading or
mutating the live campaign store.
Shared analytic dependencies may contain dormant candidate-capable definitions
and default parameter records; SID3 nevertheless opens no candidate-only
module or evaluator, runtime configuration, initial or evolution state, branch,
or outcome. Its artifact-local static decision remains
`safe_to_resume_trajectory=false`; the separately guarded live-store preflight
later authenticated generation `8`, the interpreter/runtime image, and the
restart premises for exactly one bounded writer attempt. That attempt is now
bound by PREF28 below. No common event, GR-0 calibration, SGB-L, FGC-QR,
activation, trapping, DEF1, retained-EFT, transition, mechanism, or physical
result exists.

SID3's six closure, factory, authority, bootstrap, runner, and focused-audit
test modules belong to the pinned evolution/restart partition. They
authenticate the exact NumPy package entry point, native core extension, and
distribution metadata used by the resume image. The analytic runtime has the
same NumPy version but a different package origin and metadata inventory, so
running those identity tests there is a partition error rather than an
alternative arithmetic check. This routing change weakens no identity test and
changes no equation, threshold, state, or claim.

**FGC-1-PRO19-PREF28** remains in the pinned evolution/restart partition. Its
explicit post-attempt audit reads the immutable terminal campaign, validates
the generation-eight to generation-nine append-only edge, and confirms that
all accepted member identities are unchanged. Ordinary verification consumes
only its compact record and does not reopen the ignored run store. PREF28 binds
an invalid runtime terminal; it neither diagnoses the unpersisted exception
message nor supplies a common event, calibration, candidate, or physical
result.

**FGC-1-TDG8-RCV3-FRZ1** is a compact, store-blind recovery-fork freeze. It
preserves the invalid generation-`10` source terminal, freezes the exact
generation-`9`/sequence-`10` 50-leaf prefix, and binds the sole permitted
projection-only bootstrap image. The runtime and bootstrap bytes are
authorized only to install that exact prefix under the fixed external wrapper;
ordinary verification never invokes them and never repeats the one-time
destination-absence observation.

**FGC-1-TDG8-RCV3-PREF1** now independently authenticates that installed
projection. It reads the unchanged terminal source, destination, and external
receipt without following links; proves exact 59-to-50-leaf equivalence;
reconstructs checkpoints 0--9 and journals 0--10; safely decodes and hash-binds
all 11 payloads; and verifies all six current member identities plus the exact
depth-two pending retry. This is restart-state authentication only. Bounded
execution remains unauthorized, no common event or GR-0 calibration is
complete, candidate branches remain closed, and no physical result exists.

**FGC-1-TDG8-RCV3-AUTH1** consumes that immutable binder and prospectively
authorizes only the existing persisted GR-0 event-`23` continuation. Its
compact verifier is store-blind; status reauthenticates the external receipt
and live descendant read-only; the explicit runner has no generation-zero,
bootstrap, SGB-L, FGC-QR, calibration-promotion, or direct-PDE alternative.
AUTH1 itself executes nothing and earns no event, calibration, or physics.

The first AUTH1 attempt later appended only two nonaccepted records: the TDG6
rejection at journal sequence `11` and its cursor-only retry transition at
sequence `12`. It then failed before generation-`10` checkpoint publication
because the live in-memory evidence still used tuples where the persisted JSON
schema requires lists. The recovered record bytes themselves validate after
JSON normalization, so this was a schema-boundary integration defect rather
than a state, equation, constraint, or physical stop. Commit `ee888156...`
normalizes the evidence at that boundary.

**FGC-1-TDG8-RCV3-REC1-AUTH1** is the prospective two-phase successor. It binds
the exact generation-`9` checkpoint, sequence-`11`/`12` records, stale lease,
committed fix, and unique predicted metadata-only generation-`10` checkpoint.
The explicit runner must first publish or recognize exactly that checkpoint;
any different live boundary stops invalid before PDE work. Only after exact
reconciliation may the existing engine continue the same GR-0 event `23`.
Compact verification remains store-blind. REC1 itself publishes no recovery,
advances no accepted state, completes no event or calibration, opens no
candidate branch, and earns no physical result.

The authorized REC1 run has now terminated and is bound by
**FGC-1-TDG8-RCV3-PREF2**. The independent binder reconstructs the exact
metadata-only generation-`10` recovery and every subsequent record through
generation `29` / sequence `50`. No accepted post-restart macro-step exists:
20 TDG6 rejections advance only temporal-retry metadata through retry count
`22`, after which the run stops `temporal_retry_exhausted`. All accepted state
identities and accumulated debit remain unchanged. This is a finite
numerical/instrumental terminal, not GR-0 calibration or candidate physics;
PREF2 freezes no successor remedy. Ordinary verification is compact and
store-blind, while the explicit historical binder is read-only.

**FGC-1-TDG9-AR1-AUTH1** is a prospective exact-arithmetic diagnostic
authority, not another campaign continuation. It freezes authenticated
deterministic replays of only RCV3 retry widths `3`, `4`, and `5` under the
exact CPython/NumPy/Darwin arm64 environment and full
GR-0/amplitude-3/RK4-2049/event-23 selection. Prelaunch requires the evaluator
commit at `HEAD` and exactly the declared 17-path delta; execution requires the
supplied authority commit as the evaluator commit's single direct successor at
clean `HEAD`, exact ancestry, no replacement refs, and the same committed
delta, thereby locking unlisted transitive modules.
That runtime commit is recorded in the raw manifest rather than self-bound in
the tracked authority bytes. The legacy TDG6 evidence must reproduce before
classification, and all eighteen channels pass through two independent
exact-dyadic Hermite evaluators. `exact_zero` is not an order pass; because its
exact discrete difference vanishes, it maps to the legacy-enclosure-owner
aggregate. The runner may write only a separate ignored diagnostic
manifest/result and must leave the campaign store byte-identical. AUTH1 records
no replay outcome, does not replace TDG6, and authorizes no event completion,
calibration, candidate branch, mechanism, retained-EFT, transition, or physical
claim.
**FGC-1-TDG9-AR1-PREF1** independently binds the completed two-leaf AR1
terminal without importing the runner or persisted-retry wrapper. It repeats
all `108` exact evaluator calls from generic TDG6 proposal reconstruction and
finds ten sufficient contraction failures. Thus the legacy enclosure is not
the sole owner on the three frozen samples; continuum order, event completion,
calibration, candidate dynamics, and physics remain untested.

**FGC-1-TDG9-LOC1-FRZ1** is the still-unexecuted read-only successor. It fixes
PREF1's ten failures, exactly `122,640` base D01/D12 cubics split into `V`, `S`,
and `C` (`367,920` component cubics), separate `490,560` per-family and
`1,471,680` aggregate candidate ceilings, two exact localization routes,
co-maximizer multiplicity, ownership/cancellation and Hermite-input/ULP
reporting, and the owned-row region map. It authorizes no
fourth width, resource escalation, Stage-2 source decomposition, comparator,
campaign mutation, candidate, mechanism, or physics.

LOC1's attempted diagnostic published no raw namespace: the v1 independent
root route created `368` false interior candidates in `184` tiny cubics and
stopped on `independent_disagreement`, even though both routes shared the same
two co-maximizers and exact maximum. **FGC-1-TDG9-LOC1-ERR1** binds that
instrument diagnosis. **FGC-1-TDG9-LOC2-FRZ1** prospectively changes only the
independent root route to exact endpoint deflation, exact-square evaluation,
scale-invariant bounded nonsquare isolation, and exact Sturm-count checking.
Its coefficient-depth audit is premise-only. **FGC-1-TDG9-LOC2-PREF2** now
independently binds the completed two-leaf result against the sealed 115-leaf
store and immutable authority commit. All 60 primary/v2 component-level
localizations agree; all ten selected failures reduce to independent endpoint-
state and width-scaled-RHS failures rather than mixed-only cancellation. The
result selects no temporal remedy, and all campaign continuation, calibration,
candidate, mechanism, and physical claims remain closed.

**FGC-1-TDG9-TI1-FRZ1** stopped invalid before its first shadow proposal when
the replay fingerprint rejected a legitimate finite binary64 monitor value.
It is preserved as instrument history, not reclassified as a tableau result.
**FGC-1-TDG9-TI2-FRZ1** prospectively repairs only the fingerprint encoding:
finite built-in floats use exact binary64 hexadecimal identity while the
scientific serializer stays strict. All state, grid, SBP4 RHS, projector,
transaction, tracers, ledger, widths, thresholds, exact routes, budgets, and
nonclaims remain unchanged. TI2 is still not production SSPRK3 or an
independent-method result.

**FGC-1-TDG9-TI2-PREF1** independently binds the completed tableau-only
shadow. All 60 exact primary/v2 route reconstructions agree: four of the ten
complete failures clear, six persist, two of the persistent cases become
endpoint-state-owned, and four retain independent value/slope failure. The
campaign store and PDE state remain unchanged. This finite-dimensional mixed
result selects no temporal remedy and does not constitute production SSPRK3,
independent-method agreement, event completion, calibration, candidate,
mechanism, or physical evidence.

**FGC-1-TDG9-AC1-FRZ1** prospectively authorizes one retry-3
SSPRK3-on-inherited-SBP4 read-only all-18-channel production-classifier
audit. TI2/PREF1 bound ten preselected complete-C occurrences, sampled
retry-3 `u:alpha` and `u:R`, and did not serialize all-18
`continuous_admission`. The compact certificate contains no audit outcome.
Prelaunch restores retry 3, constructs zero proposals, and leaves the
115-leaf store unchanged. All-pass would satisfy only this frozen state/width
against the production TDG6 all-of classifier.

**FGC-1-TDG9-AC1-PREF1** independently binds the completed AC1 terminal. It
reclassifies all 18 exact intervals without importing AC1 runner decision
code or reconstructing the seven SSPRK3 proposals. Seventeen channels admit;
only `u:R` is `order_inconclusive`. The campaign store and PDE state remain
unchanged. This licenses only a separately frozen diagnosis of `u:R` or a
sharper discriminator and does not constitute a production method, selected
remedy, state advance, calibration, candidate, mechanism, or physical
evidence.

**FGC-1-TDG9-UR1-FRZ1** prospectively authorizes one retry-3
SSPRK3-on-inherited-SBP4 read-only `u:R` envelope-owner diagnostic. AC1-PREF1
found 17/18 channels admitted; only `u:R` is `order_inconclusive`. The compact
certificate binds retry 3, the inherited RK4-2049/SBP4 state/operator, the TI2
`u:R` row hash `b4a48...`, the AC1 D01/D12 intervals, exact Fraction lower
reconstruction, and a closed three-value lower-owner enum. It restores one
fingerprint, executes zero shadow proposals, and contains no owner outcome.
Seven proposals and 28 records are reserved only for a later explicit RUN1.
The freeze authorizes only that later read-only diagnostic. It does not earn a
production method, select a remedy, advance state, admit a retry, calibrate
GR-0, open SGB-L/FGC-QR, or earn mechanism or physical evidence.

**FGC-1-TDG9-UR1-PREF1** independently binds the completed UR1 diagnostic. It
reconstructs the retry-3 `u:R` D01/D12 zero-lower owners without importing UR1
runner decision code or reconstructing the seven SSPRK3 proposals. Both zero
lowers are owned by
`coefficient_plus_arithmetic_debit_clips_positive_raw_maximum`. Both D01/D12
lowers are zero under the current global aggregation; raw maxima are strictly
positive; the current globally aggregated coefficient-construction debit
alone is sufficient on this frozen state/width with exact construction/raw
ratios D01 `654019454794891/25865933815808` (~25.285) and D12
`163504863789045/5012309316608` (~32.621); the arithmetic debit is negligible
and arithmetic alone is not sufficient. UR1 does not prove that raw and
construction maxima are co-located, and it does not prove that
candidatewise/local pairing repairs admission. Bernstein is recorded/excluded;
the row hash matches TI2; the intervals match AC1; and the production
classifier remains `order_inconclusive`. The campaign store and PDE state
remain unchanged. This licenses only prospective design of a separately frozen
envelope/admission change and does not constitute a production method or
comparator, selected remedy, state/retry admission, event, GR-0 calibration,
SGB-L/FGC-QR opening, candidate, mechanism, EFT, transition, or physical
evidence. FGC-QR remains unopened.

**FGC-1-TDG10-QA1-FRZ1** prospectively authorized one retry-3
SSPRK3-on-inherited-SBP4 exact complete-C all-18-channel qualification. **At
freeze it was unexecuted:** zero shadows had been measured, no channel
outcome existed, and no state was advanced. Retry 3 restores generation 9 /
journal sequence 10 only; sequence 11 remains the historical retry-3
rejection; generation 10 is not an executable retry-3 predecessor. The freeze
authorized only a later diagnostic RUN1 and remains the historical compact
freeze. A later gitignored RUN1 now reports raw class
`completed_exact_all_channel_pass_no_state_advance` under authority
`aeaf0498...`, with all 18 raw channel summaries `resolved_order_pass`, 7
shadows/28 records, a nonaccepted diagnostic endpoint, and an unchanged
115-leaf store. **FGC-1-TDG10-QA1-PREF1** independently binds that completed
terminal: compact SHA-256
`3a107bcede479df4326aac5e194beefa1042860a8be498d7dc5d0519fdbc8603`, class
`independently_bound_retry3_ssprk3_sbp4_exact_complete_C_all_channel_pass_terminal`.
Independent replay reconstructs 7 shadow paths / 7 proposals / 28 SSPRK3
stage+endpoint records. All 18 channels independently classify
`resolved_order_pass`; `failed_channels=[]`. Both exact rational localizer
routes agree. The tightest channel is `u:R`; exact ratio `8*U12^2/L01^2` is
approximately `0.9931888298284983`, leaving about `0.006811170171501636`
below the pass boundary. Four raw hashes and the 115-leaf store remain
unchanged. The diagnostic endpoint remains nonaccepted. Ordinary verification
consumes only the compact certificate and is raw/store/shadow/Git blind.
`--live` is explicit one-time construction history. QA1 licenses only
QA2/method design, never old-member adoption or an RA1 transplant of the
already-known diagnostic endpoint. RA1 does not exist and is not authorized
by PREF1 itself. It does not earn a production SSPRK3 comparator, two-method
agreement, state advance, common event, GR-0 calibration, SGB-L/FGC-QR, or
mechanism or physical evidence.

**FGC-1-TDG10-QA2-FRZ1** prospectively authorizes one retry-4/5
SSPRK3-on-inherited-SBP4 exact complete-C two-width robustness qualification.
**At freeze it is unexecuted:** zero shadows have been measured, no channel
outcome exists, `width_robustness_passed` is false, and no state is advanced.
Retry 4 restores generation 10 / journal sequence 12 only; sequence 13 remains
the historical rejection. Retry 5 restores generation 11 / journal sequence 14
only; sequence 15 remains the historical rejection. QA1 licenses this
method-design freeze, never old-member adoption. The later RUN1 is
diagnostic/read-only: no temporal admission, no campaign write, and no
diagnostic endpoint serialization. Both-width all-18 pass would license
design of a truthful production method/admission owner; any nonpass rejects
this exact SSPRK3+SBP4+exact-C remedy as a general production remedy on the
tested retry neighborhood. It does not earn a production method, two-method
agreement, common event, GR-0 calibration, SGB-L/FGC-QR, candidate, mechanism,
or physical evidence.

The later explicitly authorized QA2 RUN1 produced a two-leaf raw terminal.
**FGC-1-TDG10-QA2-PREF1 now independently binds it** with compact SHA-256
`08f88f02966ea5936e417a3a6a2c583e11783341a3d76f4ef5dc821e992691f6`
and class
`independently_bound_retry4_5_ssprk3_sbp4_exact_complete_C_two_width_nonpass_terminal`.
The binder authenticates raw hashes `8638ac8c...` and `252810f4...`, then
independently reconstructs 14 paths, 14 proposals, 56 SSPRK3
stage-plus-endpoint records, and all 36 channel decisions without importing
QA2 decision code. Both exact-rational routes agree. Retry 4 and retry 5 each
classify exactly `u:alpha`, `u:lambda`, and `u:R` as nonadmitted; the sealed
115-leaf store remains unchanged. The result rejects only this exact
SSPRK3-on-inherited-SBP4 plus exact-C remedy on the tested retry neighborhood
and selects no successor. Ordinary PREF1 verification is compact-only and
never launches QA2 or opens the raw namespace or campaign store. QA1 licenses
QA2 method design only, never RA1 or old-member adoption.

Ordinary analytic verification still treats the completed UR1-PREF1 compact
binder as the owner of one **test-partition filtering rule**. This is runtime
routing, not selection of a scientific or production successor. QA2-PREF1
rejects the tested exact-C remedy on retries 4/5 and selects no successor. The
QA1 freeze remains the historical unexecuted compact freeze. Before the
analytic partition removes any listed case, it runs
`scripts/reproduce_fgc_tdg9_ur1_pref1.py` in its check-only, store-blind form:
no `--live`, `--run`, or status path. If that compact check fails, filtering
does not occur.

The nine prelaunch-only observations listed under
`partition.tdg9_ur1_pref1_postlaunch.sealed_prelaunch_test_ids` are then sealed
out of analytic ordinary verification. They preserve historical live
projection, copied-store recovery, owner-token absence, and output-absence
premises that were true at authorization time. They are not current-store
reproducers, and this is not a blanket seal of every successor-state mismatch.

Two later individual reconstructions listed under
`evolution_restart_reconstruction_test_ids` are not sealed. They remain live
runtime recomputes: the analytic partition proves both IDs exist, removes
exactly those two cases, and the evolution/restart partition appends exactly
those two under `EVOLUTION_PYTHON` when the complete local raw bundle exists.
They reconstruct a sealed LOC1 tuple and a retry-3 D01 family. Compact
siblings in the same modules stay analytic. When the raw bundle is absent from
a clean clone, evolution remains skipped and those two reconstructions do not
run.

Pure arithmetic, mutation, and scientific gates in the same files remain
active in the analytic partition. No equation, threshold, campaign state, or
physical claim is changed by this routing.

## TDG11 selection and independent binding ownership

[FGC-1-TDG11-MSEL1-FRZ1](fgc-tdg11-msel1-frz1.md) keeps compact verification
Git-, raw-, store-, shadow-, and environment-blind. Its explicit live diagnostic
requires the frozen CPython 3.14.3 / NumPy 2.5.1 / Accelerate image, with the
Python and NumPy extension hashes and OS image recorded prospectively.
The Homebrew analytic runtime is not a substitute for that live image.

It restores generations 9, 10, and 11 for retries 3, 4, and 5. The synchronized
HLT15 `t=23/16` origin remains a separate future production start. The fixed
maximum is 42 proposals, 210 stage/endpoint records, and 42 separately counted
accepted-source prechecks; the six static templates have their separate
factory prechecks. Only a two-leaf diagnostic namespace may be published.
An unbound raw selection grants no IMP1 or production permission.

[MSEL1-PREF1](fgc-tdg11-msel1-pref1.md) independently reconstructs the
completed original-RK4 families and reproduces all three compensated
arithmetic premises. Actual accounting is 21 completed proposals, 105
stage/endpoint records, 24 accepted-source prechecks, and 129 shadow RHS
calls, plus the six static shells. The completed run did not use the full
prospective ceiling. It serialized no endpoint and changed no campaign store.

The binder's live operation uses the same exact evolution image. Its ordinary
compact route is raw-, store-, shadow-, environment-observation-, and
Git-blind: it verifies tracked source bindings and reconstructs the original
raw terminal commitment from saved independent records. The PREF1-specific
localizer adapter ID remains visible; only that ID is mapped when comparing
the frozen raw wire format. The bound selection opens a separate IMP1 owner,
not production evolution or a physical result.

## TDG11-IMP1 implementation and synthetic qualification

[IMP1](fgc-tdg11-imp1.md) retains original binary64 RK4/SBP4 and
SSPRK3/SBP2 proposals. Comparison uses exact accumulation-corrected Hermite
cubics, exact Bernstein sufficient proofs and the inherited dual-exact
localizers on unresolved comparisons. The order law stays `p>=3/2`.
Its public debit uses the retained Bernstein finest-pair upper plus the full
fine accumulation bound; any extra enclosure slack is separate from the
intersection-gate debit and is not called roundoff or PDE error.

The new in-memory ledger freezes the complete inherited TDG6 snapshot and
tracks IMP1 time, counters, exact debt and rejection history separately.
Only the actual binary64 four-substep endpoint can be adopted. Source/health
stops, stale boundaries, invalid retries and late tracer failures do not
produce a partial real commit. A new checkpoint codec is not itself a
durable campaign cursor or proof of historical-origin authentication.

Explicit `--qualify --write-result` construction uses the prospectively
pinned CPython 3.14.3 / NumPy 2.5.1 / Accelerate environment before and after
synthetic controls. It verifies both integrators with independent 9-point
reference reconstruction and 129/2049-point size controls; no physical source,
campaign store, raw namespace or scientific endpoint is consumed or written.
The deterministic report commitment is `a3170dd4...`; timings are observations,
not gates or part of that commitment.

Ordinary `make fgc-tdg11-imp1` checks compact/config/source commitments without
Git, environment observation or shadow execution. Focused tests supplement
the SF1 foundation and canonical verification.

`FGC-1-HLT17-SRCQ1-REC1-PREF1` bound the physical origin/source and C1R1
preparation for all six GR-0 members without writing an endpoint or campaign
state. That licensed only the separate PRO20 freeze. `FGC-1-PRO20-EV1-PREF1`
now independently binds the closed production terminal. Ordinary
`make fgc-pro20-ev1-pref1` is compact-only.

The `PRO20-EV1` persistence sibling implemented the exact six live-member
production-shaped root, journal, checkpoint, generation manifest,
writer-session and terminal-lock schemas. The corrected freeze owned the
one-shot run. Independent PREF1 now binds that closed terminal as a typed
`resource_exhausted` stop with partial accepted state and no common event.
Ordinary compact verification is Git/raw/store/source/Planck/replay-blind.
The live run target is tombstoned. All generation/session/terminal claim
flags remain false for calibration, SGB/DEF, holdout, candidate, mechanism
and physics.
