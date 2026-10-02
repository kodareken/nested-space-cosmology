# Stage-4 source family

The question is whether a change in regional source shape or imbalance
moves contraction, dispersion, or retention off the separated baseline.
The executable owner is
[nsc_discovery_family.py](../src/recursive_horizons/nsc_discovery_family.py).
Checks live in
[test_nsc_discovery_family.py](../tests/test_nsc_discovery_family.py).
The dense radius solve remains
[nsc_nested_parent_child.py](../src/recursive_horizons/nsc_nested_parent_child.py).
The march from the prepared slice to the handoff is
`nsc_discovery_prediction.propagate_baseline`. Later stations are the
existing episode continuation. This note does not add a third time stepper.

## Family

Imbalance \(a\) is one of \(-1,-1/2,0,1/2,1\). The six CAR weights are

\[
c=(0.5+0.25a,\ 0.5+0.25a,\ 0.5,\ 0.5,\ 0.5-0.25a,\ 0.5-0.25a).
\]

Their sum is \(3\) at every \(a\). The audited occupations are the
endpoint \(a=+1\). The uniform member is \(a=0\). Every listed vector
lies in \((0,1]\), and the covariance eigenvalues of an orthonormal
frame with those weights stay in \([0,1]\).

Child half-width \(w\) is the scale \(0.8\), \(1\), or \(1.2\) times the
baseline half-width \(1\). The source supports are the central interval
\([2-w,2+w]\) and the outer annuli \([0,2-w]\) and \([2+w,4]\). Clocks
stay at \(x=1,2,3\), and the geometric child window stays \((1,3)\).
The source windows are not a new nested interval.

Packets use the owned bump and lobe preparation. Outer columns are
projected off the child and Löwdin-orthonormalized. That step is
recorded. The child columns are not included in it. At unit scale the
child pair is the owned original middle pair, unchanged. At the other
scales the child lobes are an explicit declared frame. If that pair is
not already orthonormal, an explicit two-column Löwdin is recorded
before the child is frozen. The mean source current is measured and is
not deleted. Gram defect and fine-grid tails are measured.
`exact_spatial_support` stays false.

The chart coefficients \(A\), \(C_W\), flux, and \(\kappa\) are the
locked audited values. The geometry basis for \(n_f=128\) and \(n_f=256\)
is the frozen portable \(W\), checked as a complete orthonormal frame.
The observer is the explicit rank-6 reference, stored separately from
the separated source.

## Initial data and handoff

Each evolved member calls `initial_state` on its own columns and
weights. The solve is the source homotopy of the owned projected radius
constraint. The saved \(T=0.3\) radius is not copied. The full shift
residual maximum, RMS and mean, and the mean current stay in the initial record.

At \(T=0\), \(\Phi\) equals those prepared columns by value. The arrays
are copies, not aliases. The episode pin check treats that value
equality as an alias (`state Phi is aliased to the source columns`).
This driver therefore does not hand the \(T=0\) slice to that runner.
`propagate_baseline` advances the prepared state, and the handoff is
that output: the evolved field, the stage-integrated child clock at
\(x=2\), endpoint-trapezoid clocks at \(x=1\) and \(x=3\), and the
source identity. The continuation sets `initial_state_called` false so
the episode does not solve again.

The default non-production duration is \(0.01\). Coordinate time
\(0.3\) requires `--production`. That call belongs to the root executor
after this interface is frozen.

## Six cases, one pool

The admissible catalog has all \(5\times 3\) members. The exploration
batch evolves six, each on one thread, at \(n_f=128\):

| Case | \(a\) | width scale | Role |
|---|---:|---:|---|
| baseline | \(+1\) | \(1\) | audited occupations, nominal lobes |
| uniform | \(0\) | \(1\) | imbalance removed |
| endpoint | \(-1\) | \(1\) | opposite endpoint |
| narrow | \(+1\) | \(0.8\) | shape only |
| wide | \(+1\) | \(1.2\) | shape only |
| crossed | \(-1\) | \(0.8\) | endpoint imbalance and narrow shape |

\(a=\pm 1/2\) remain in the catalog and are not given evolution slots in
this bounded batch. The corner \(a=-1\), width scale \(1.2\) is in the
column catalog. At \(n_f=64\) and \(n_f=128\) the owned positive
homotopy stalls above its radius-correction floor. That member is not
given an evolution slot, and no other radius is substituted for it. The `--prepare --all-members` option admits all 15 catalog members into the same
six-worker pool and aggregate budget. Default preparation admits the six cases
above; the saved catalog identifies the remaining nine. Dense-solve failures
are retained with their original blocker, elapsed CPU and an unresolved fate.
At the short test resolution, other shapes can also fail the homotopy; a
constructor failure does not imply non-existence. Each successful handoff is
committed before the next source is attempted. Confirmation is a separate directory at \(n_f=256\).
Its successor confirmation admits the five predeclared candidates below,
even when their qualitative signs agree. Up to one additional distinct-sign
or failed-preparation member can fill the sixth source slot. Candidates are
deduplicated by imbalance and width; the baseline appears once.

| Imbalance | Width | Admission purpose |
|---:|---:|---|
| +1 | 0.8 | quantitative narrow-width sensitivity |
| +1 | 1 | nominal width and baseline comparator |
| +1 | 1.2 | quantitative wide-width sensitivity |
| 0 | 1 | uniform-imbalance control |
| -1 | 1.2 | NF256 retry of unresolved initial preparation |

Each selection retains the predecessor case, reported fate, original
preparation blocker when present, and selection reasons. Admission schedules
an experiment; it does not force a physical conclusion or assert that the
higher-resolution solve will succeed. The bounded default is five independent
source preparations at NF256, not a demand for different qualitative signs.

`--frozen-width-controls` optionally adds frozen-geometry branches for the
three +1 width cases. Each branch copies its own source's unchanged propagated
handoff, including source weights, field, metric, momenta and clocks. It does
not borrow the nominal-width metric, solve again, or reset the field. IDs end
in `_frozen_geometry` and retain the coupled parent ID. These are three control
branches in addition to the source-candidate cap of six; the default therefore
has five sources and at most eight continuation cases. All cases share one
six-worker episode executor and the aggregate budget. A failed source has no
frozen substitute branch.

The aggregate CPU budget is \(21600\) seconds. Measured preparation
CPU is subtracted before the continuation pool opens. Preparation includes the
dense solves, failed attempts, measurement consumers and the serial handoff
evolution. The record separately reports initial-solve and handoff CPU.
Continuation counts accumulated child CPU and coordinator assessment CPU
across resumed calls. The six-hour
cap covers preparation and the child CPU together. The serial prefix checks
the budget between members; it cannot interrupt an individual dense solve or
`propagate_baseline` call mid-flight. Intermediate prefix states are not
checkpointed by that existing propagator. The pool is the
episode runner, opened once. Forecast and checkpoint chunks stay at
\(64\,\mathrm{MiB}\). The rate-sample forecast does not include the dense
solves and is not a stop.

## Regime observable

Classification reads actual checkpoints, with each member's own initial sample.
It stores child proper-length change, packet-width change, retained fractions
and their rank, and a 3-by-3 source-packet/window probability matrix. Rows of
that matrix are the left, child and right source pairs; columns are their
three declared coordinate support windows. The same tagged child columns
2 and 3 also report retention inside the fixed spatial child window `(1,3)`,
with the same total tagged-pair probability denominator. Own preparation-window
retention is retained under a separate window label. This prevents a changed
width and a changed measurement window from being mistaken for the same effect.
Differences give measured packet
transfer without treating computational ancestry as independent energy parcels.

The metric supplies areal-radius change. The tidal consumer supplies actual
Ricci scalar and radial/angular tide norms and their growth, with no use of
chi as a curvature substitute. Those signs join contraction, dispersion,
retention rank and transfer signs in the comparison. Geometry and probability
positivity are retained. Cases ending at different coordinate times from the
baseline remain unresolved. If a consumer raises, the record retains its
exception text and missing measurements. A missing curvature consumer or
nonfinite geometry also leaves the pattern unresolved. These finite-grid
sign comparisons carry no continuum or renewal certificate. Raw quantitative
differences in proper length, packet width, both retention readouts, areal radius
and actual curvature norms are also reported against the baseline at the same
coordinate time. Frozen branches additionally report differences against their
own coupled parent, without applying a universal percentage gate.

A scalar fate class does not apply a universal \(1\%\) gate. A small
change and a large change with the same signs are the same class. The
reported discriminating observable is the child-source retained
fraction, together with the left-minus-right occupation gap, which
equals \(a\) for this weight formula.

## Reproduction

```sh
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_family.py -q
```

The tests stay at or below coordinate time `0.01`. Confirmation preparation
reads the immutable predecessor and creates a new successor. Run these from
the repository root; `scripts/lab.py` changes the scientific working directory
to `lab`, so relative evidence paths below start with `results/`:

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_family.py \
  --confirm --exploration results/development/nsc-discovery-family-v1 \
  --output results/development/nsc-discovery-family-v2 --frozen-width-controls
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_family.py \
  --materialize --output results/development/nsc-discovery-family-v2 \
  --production --duration 0.3
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_family.py \
  --run --output results/development/nsc-discovery-family-v2 --workers 6 --backend fft
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_family.py \
  --check --output results/development/nsc-discovery-family-v2
```

For a fresh NF128 exploration, `--prepare --output <new successor>` remains
available; it accepts `--all-members` and `--frozen-width-controls`.
Materialization and continuation reject writes to the immutable family-v1
path. Read-only `--check` can still verify that predecessor.

Successors store producer hashes for the local Python import closure. A file
gets a historical commit pin only when the immutable Git blob authenticates
its exact bytes. Uncommitted producer bytes are explicitly marked as working
tree pins. Commit the frozen producer before production to obtain historical
replay identity. The checker authenticates committed historical bytes through
Git or the existing source-history resolver, independently of current verifier
hashes; it never heals a scientific JSON record.

Family-v1 itself has no producer identity envelope. The checker says so. When
`observed-run-binding.json` is present, it separately authenticates that
external envelope's historical producer blobs and artifact hashes. That
binding remains a post-run coordinator observation, not a retroactive pre-run
attestation. The current verifier's identity is reported separately.
