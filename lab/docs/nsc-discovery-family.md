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
It admits up to five members whose measured sign pattern differs from the
baseline, together with an independently prepared baseline comparator at the
same confirmation resolution. An unresolved pattern is not regenerated.

The aggregate CPU budget is \(21600\) seconds. Measured preparation
CPU is subtracted before the continuation pool opens. Preparation includes the
dense solves, failed attempts, measurement consumers and the serial handoff
evolution. The record separately reports initial-solve and handoff CPU.
Continuation counts accumulated child CPU across resumed calls. The six-hour
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
three declared coordinate support windows. Differences give measured packet
transfer without treating computational ancestry as independent energy parcels.

The metric supplies areal-radius change. The tidal consumer supplies actual
Ricci scalar and radial/angular tide norms and their growth, with no use of
chi as a curvature substitute. Those signs join contraction, dispersion,
retention rank and transfer signs in the comparison. Geometry and probability
positivity are retained. Cases ending at different coordinate times from the
baseline remain unresolved. If a consumer raises, the record retains its
exception text and missing measurements. A missing curvature consumer or
nonfinite geometry also leaves the pattern unresolved. These finite-grid
sign comparisons carry no continuum or renewal certificate.

A scalar fate class does not apply a universal \(1\%\) gate. A small
change and a large change with the same signs are the same class. The
reported discriminating observable is the child-source retained
fraction, together with the left-minus-right occupation gap, which
equals \(a\) for this weight formula.

## Reproduction

```sh
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_family.py -q
```

The tests stay at or below coordinate time `0.01`. After freeze, the
root executor materializes the handoff and continues it:

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_family.py \
  --prepare --output lab/results/development/nsc-discovery-family-v1
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_family.py \
  --materialize --output lab/results/development/nsc-discovery-family-v1 --production --duration 0.3
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_family.py \
  --run --output lab/results/development/nsc-discovery-family-v1
```

`--confirm --output <confirmation> --exploration <exploration>` writes
the selected \(n_f=256\) list in a different directory. It does not
evolve those members until a later `--materialize --production`.

`--check --output lab/results/development/nsc-discovery-family-v1` validates
the specification and immutable episode chunks. These commands create a new
successor; existing committed chunks are preserved. The implementation tests
do not create repository evidence or run a production trajectory.
