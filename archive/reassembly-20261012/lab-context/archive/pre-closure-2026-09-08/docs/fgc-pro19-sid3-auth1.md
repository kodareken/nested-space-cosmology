# FGC-1-PRO19-SID3-AUTH1: committed generation-eight resume image

**Status:** premise-only committed-image authority. SID3 freezes the exact
GR-0 execution-source and immutable-input image that may be presented to a
later read-only live preflight. The compact artifact does not inspect the
mutable campaign store and does not authorize trajectory resume by itself.

## Two-commit authority

The implementation commit is fixed exactly as:

```text
A = 4d7133c0737524712220f77ee5ae4d9c3550c75a
```

Commit A owns the bootstrap, execution-closure verifier, exact module
permission map, resume authority, static factory, runner, and tests. The later
authority commit C is intentionally not embedded in A or in this document:
the operator supplies C as a full commit identifier, and the bootstrap
requires the current clean `HEAD` to equal C. A must be an ancestor of C, and
the complete `A..C` changed-path set must equal the closure's declared
non-importable authority delta. Every closure-bound Python or input file must
still have A's Git blob and exact live bytes at C. This permits the prospective
config, closure, result, and publication updates without relabelling C as the
implementation source or admitting later source drift.

The tracked compact bundle is:

- `configs/fgc/fgc-1-pro19-sid3-auth1.toml`;
- `configs/fgc/fgc-1-pro19-sid3-execution-closure.json`;
- `results/fgc-1-pro19-sid3-auth1.json`; and
- this owner document.

`make fgc-pro19-sid3-auth1` reconstructs and checks only that compact
authority plus its bounded source tests. The separate
`make verify-fgc-pro19-sid3-prelaunch` command additionally runs the focused
repository check. Neither target opens
`runs/fgc-2-sf1/proto17/calibration`, acquires a writer, or executes a PDE
proposal.

## Exact execution closure

The closure binds **75 files** at A: **67 Python files** and **8 static input
files**. Sixty-six Python files form the complete A-owned
module-to-repository-origin permission map; the remaining Python file is the
committed-stdin bootstrap itself. The import graph uses exactly four synthetic
namespaces:

```text
recursive_horizons
recursive_horizons.fgc
recursive_horizons.fgc.evolution
scripts
```

Synthetic namespaces prevent repository package initializers from expanding
the graph. The 66 module names and their conventional `src`-before-root
origins must equal A's complete approved map, not merely contain a mandatory
subset. The generic closure separately binds every permitted file's Git blob
identifier and SHA-256, its no-follow live bytes, the isolated interpreter
identity, and the exact NumPy package, extension, and distribution-metadata
origins. Symlinks, hard-link aliases, untracked or ignored source shadows,
duplicate origins, path escapes, dirty state, later source drift, and any
observed repository module outside this map fail closed.

The eight immutable inputs are:

```text
configs/fgc/fgc-1-cal9-run1.toml
configs/fgc/fgc-1-pro19-sid2-pref1.toml
configs/fgc/fgc-1-rsp2-run1.toml
results/fgc-1-hlt10-mon10.json
results/fgc-1-pro18-auth1.json
results/fgc-1-pro19-frz1.json
results/fgc-1-pro19-sid2-pref1.json
results/fgc-1-rsp2-frz1.json
```

The isolated runtime is the recorded CPython `3.14.3`/NumPy `2.5.1` Darwin
`arm64` environment. An equivalent-looking interpreter, package install,
extension, metadata file, or import origin is not interchangeable with the
recorded image.

## Bound generation-eight ancestor

SID3 composes SID2's metadata-only recovery with the original GR-0 amplitude-3
event-23 plan. Its exact recovered ancestor is:

```text
checkpoint generation:     8
checkpoint SHA-256:         6d2b6e37e7c7cfd4ce64504353d886099003f959ea2aaee1a14fc778ce14fa46
journal sequence:           8
journal tip SHA-256:        949b29e46081bf5a42c34461b3894419db86aeca6c733b7a2636c91c6cb2ee4a
member:                     RK4-2049
member accepted time:       0x1.78554de5a30e0p+0
member cursor SHA-256:      6d594999f6aaebb5d56f0fd6ff12499e2b22a47e3daf74ac8b13b15d273fa177
member mode:                RETRY_PENDING
pending owner:              temporal
pending cap:                0x1.aaa9612df9000p-10
```

Generation eight is nonterminal. Its recovery advanced no accepted physical
field and did not replay the rejected proposal. The original authorization
commit remains `c11ba422ce49ddd6f6175de1a9d2da675b97f2de`, and the original
plan SHA-256 remains
`f3b365978aeb2ae6063742807a3161bf5d6aa0c8392cd283fdec2a79d3e46060`.

## Shared candidate-capable definitions

The exact GR-0 import graph contains shared analytic implementation modules,
not 66 files whose definitions are all GR-0-only. In particular,
`vectorized_source._coerce_action` is generic code that can construct an
`ActionParameters` record labelled `ModelID.FGC_QR` when it is called with the
corresponding frozen record shape. Other shared modules may likewise define
or instantiate default parameter records used across branches.

That source-level capability is disclosed rather than relabelled as candidate
execution. SID3's static factory and runner bind GR-0, amplitude `3`, the
PROTO12 GR-0 operator, the six declared members, and the original progression
plan. Candidate-only module and script prefixes are forbidden. No candidate
runtime configuration, initial or evolved state, branch selection, evaluator
entrypoint, or outcome is opened. Loading shared definitions is not evidence
that every callable they contain ran, and SID3 does not claim that every
permitted file contains only GR-0 definitions.

## Transient live preflight and explicit resume

The compact result deliberately records:

```text
safe_to_resume_trajectory = false
trajectory_resume_authorized_by_artifact_alone = false
```

Only the isolated committed-stdin bootstrap may turn that static premise into
a transient live decision. It loads its own bytes from A, resolves the current
clean authority commit C, verifies the complete A-to-C execution closure,
installs the guarded `src`-before-root importer, authenticates the exact
generation-eight ancestor, reconstructs all six current GR-0 members in
memory, and checks every runtime premise without taking a writer or creating
output.

Run that read-only operation explicitly with:

```bash
make preflight-fgc-pro19-sid3-resume
```

A passing response may report `safe_to_resume_trajectory=true` only for that
observation. It is not a durable artifact claim and expires if the commit,
worktree, interpreter, imports, store, writer state, or campaign lineage
changes.

The state-advancing command is deliberately separate and absent from every
ordinary test, check, compact verification, and foundation dependency:

```bash
make resume-fgc-pro19-sid3-event1
```

Both Make targets obtain the bootstrap with `git show A:...`, pipe it directly
to the pinned evolution interpreter as `python -I -B -`, clear the ambient
environment except for a fixed executable search path and `LC_ALL=C.UTF-8`, disable
optional Git locks, filesystem monitoring, and the untracked cache, and pass
the full current `HEAD` as C. Resume does not rely on an earlier preflight: the
guarded runner repeats the complete read-only preflight immediately before it
may acquire its writer lease and mutate the campaign.

## Exact nonclaims

SID3 does not, by compact verification alone:

- inspect, authenticate, repair, lock, or mutate the live campaign store;
- acquire a writer lease, execute or replay a PDE proposal, reimport GEN0,
  advance a state, complete a common event, or resume a trajectory;
- make `safe_to_resume_trajectory` true or predict that a later live preflight
  will pass;
- complete the GR-0 calibration, select an eligible case, or classify the
  eventual event outcome;
- authorize or execute SGB-L, FGC-QR, DEF1, ROB1, or any candidate branch;
- open a candidate runtime configuration, initial state, evolved state,
  branch-selection path, evaluation entrypoint, or outcome;
- establish retained-EFT validity, collapse, trapping, activation,
  defocusing, a transition, singularity resolution, a child universe, or any
  other physical result;
- prove that shared candidate-capable definitions are absent from the GR-0
  source graph; or
- protect against an actively malicious same-account process that can rewrite
  the repository, executable environment, or campaign while the operation is
  running. That stronger claim requires an external or operating-system trust
  root not supplied here.

The only positive compact claim is that one exact prospective committed GR-0
resume image and its immutable generation-eight inputs are frozen. Actual
permission to resume exists only inside a passing guarded live preflight and
the runner's immediate pre-mutation recheck.
