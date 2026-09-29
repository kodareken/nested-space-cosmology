# Resumable bounds for the middle source window

The [reviewed correction](nsc-vacuum-source-correction.md) encloses one original
group14/mid24_1 row. The [uniform expansion](nsc-vacuum-source-remainder.md)
covers the same family's archived columns only for absolute energy at least 32.
This owner schedules the correction across the remaining original angular-+1
window, $16\le E<32$: 48 positive rows and their 48 archived negative
partners. It does not replace those proofs and it does not start that campaign
by being imported.

Commit `d8ccbf5` is the reviewed single-row correction and continuous
source-insertion method. The row mathematics below is that method. This
resumable driver has not had a separate independent review.

## What each completed row contains

The positive energy is evolved as the correction

\[e_y=A_E(y)\,e-r_M(y),\]

with the same order-8 vacuum expansion, DOP853 settings, degree-12 defect
integral and start $y_0=-18$ as the reviewed row. $A_E$ is the skew Bloch
generator. $r_M$ is the exact finite-expansion defect. The numerical
correction starts at zero, and that value is not treated as the true horizon
state. The horizon-to-start remainder

\[\|e(y_0)\|\le C_M^{\mathrm{prep}}(y_0)/E^M\]

is required to be finite and strictly positive, and it is counted once inside
the validated vacuum Bloch error.

Thermal and coherent occupations stay on the original exponential law. They are
added outside the vacuum curve, with the vacuum expansion constant set to zero
in that term so the curve is not counted twice. The bound is finite for every
row in this window.

The archived comparison uses the corrected Bloch endpoint on the positive
columns. The negative partner uses the complemented endpoint $(n_x,-n_y,-n_z)$
against that partner's own archived columns and covariance. Copying the positive
scalar bound, or forgetting the complement, is rejected. Source columns,
weights, occupations and the original preparation are not rewritten.

## Checkpoints, resume and aggregation

Each completed positive row is one exclusively published directory

`rows/row-group14__mid24_1-NNNNN/`

containing canonical `record.json` and `witness.npz`. Publication is atomic and
does not replace an existing directory. The record binds the panel, row, both
energy signs, both source digests, both preparation digests, the numeric
settings, and the source/dependency hashes of this owner, the reviewed
correction driver, and the inputs those owners read.

Resume authenticates every saved witness before it is allowed to capture.
Authentication replays the stored curve and does not call the solver. `--check`
does only that replay. `--row-budget` is the maximum number of new positive
rows solved in that invocation; zero revalidates saved rows and solves nothing.
The CPU budget is process time for those new solves. Each success is published
before the next row starts. A timeout discards the unfinished row. A saved row
is never solved again.

The aggregate lists completed rows only. Missing rows contribute an identity,
not a bound. Incomplete coverage stays `OPEN`. Completing this window still
leaves `physical_upstream_budget_component` null: lower energies, the other
positive angular family, quadrature, response contraction and the rest of the
local gate are not supplied by these 96 signed comparisons.

## Commands

From the repository root, with the validation environment:

```sh
.venv/validation/bin/python scripts/lab.py -m pytest -q tests/test_nsc_middle_source_coverage.py
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_middle_source_coverage.py --check --output <checkpoint-directory>
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_middle_source_coverage.py --resume --row-budget 1 --cpu-budget 60 --output <checkpoint-directory>
```

`--check` writes nothing when the directory is absent and never solves.
`--resume` creates the directory and will not run unless both budgets are
present. The default output path is
`results/development/nsc-middle-source-coverage-v1`. No production campaign
has been published there.

The tests solve at most two original positive rows, in a temporary directory.
They then replay and resume those saved witnesses without another solve, and
they refuse a changed or missing payload, changed settings, a changed source,
a duplicated identity, and a negative comparison that skips the complement.

## Limitations

This is a local effective-source bound on one declared slice. It is not an
inheritance identity, not a derivation of an infinite nest, and not an N/beta
residual. The saved source law is an input. Quadrature weights are not applied
to these unweighted operator distances, and the distances are not contracted
through the response kernel. Other groups, group14's other positive angular
family, and energies outside $16\le E<32$ remain out of this window. A
partial checkpoint directory must not be overwritten in place; a missing
payload or a rejected identity stops that resume instead of inventing the row.
The physical gate remains OPEN.

## Observed bounded pilot

On 29 September 2026, the two-row temporary-directory pilot covered energies
16.019251120011916 and 16.10108577611476 with their actual negative partners.
Each witness used 484 cells. The largest unweighted covariance-error upper was
1.0027417326013562e-11. Capture plus validation took about 8.50 CPU seconds;
saved-witness replay took about 2.88 seconds and captured zero new rows. The
remaining 46 positive rows were explicitly missing. Parent integration ran
these tests together with the spinor tests: ten passed in 25.13 seconds.
These measurements are pilot evidence, not full-window coverage.
