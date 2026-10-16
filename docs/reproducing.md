# Reproducing the root project

There is one package, `nested-space-cosmology`, importing `recursive_horizons`
from root `src`. Use the pinned versions in [pyproject.toml](../pyproject.toml).
Scientific conclusions retain the domains in their owning notes.
The current paper and test dependencies use Matplotlib 3.11.2; the scientific
stack includes python-flint 0.9.0 and Numba 0.68.0.

Saved evidence bytes, input identities, array shapes, discrete counts and
symbolic identities retain exact checks. Recomputed floating-point arithmetic
uses the numerical allowance of its validator. The clock boundary checker
uses relative `2e-11` / absolute `8e-12`; the rectangle diagnostic comparison
uses relative `2e-12` / absolute `2e-13`. Compact-source checks propagate
binary64 coordinate rounding through the source polynomial and its derivatives.
Source-frame checks also propagate Hermite reconstruction roundoff through
derived constraint residuals; cancellation is assessed against its operands.
Reference-geometry derivative allowances account for the actual spline knot
spacing, while retained coordinates and copied inputs remain exact where declared.
Restored forcing grids allow `2e-13` absolute drift in the independently rebuilt
wall coordinate. Their saved affine node spacing and fixed source endpoint
still use the original eight-epsilon check; saved arrays are never moved.
Restored reference splines evaluate a saved endpoint just outside their live
domain at the live endpoint, within the already admitted eight-epsilon bound.
Interior queries are unchanged; larger excursions and nonfinite values fail.
Departure-history constraint replay propagates rounding through the analytic
vacuum, terminal subtraction, centred stencils and moving-grid coefficients.
Only `E_Ubar`, `E_x` and their numeric summary leaves receive those additional
operand-dependent allowances. Mass, source constraints, counts, masks and
other diagnostics keep their existing checks. The eight-unit elementary-function
assumption is a replay policy, not a certified SciPy accuracy bound.
The compatible, background, SSPRK3 and midpoint-complete validators all
forward these transient bounds to their departure-history comparisons.
Midpoint profile brackets are checked against their own floating query and
knot coordinates. A changed integer bracket must be adjacent and correspond
to one unique common knot within sixteen epsilon units of both coordinates;
crossing counts must exactly follow those brackets. Saved-cut coordinates
and metadata remain exact. This handles a one-ULP knot-side difference seen
in GitHub's PDF authentication, without a general integer tolerance or any
change to saved data. The standalone midpoint check passes with 1,418,330
charged points, two final radial stages and no source solve or evolution.
Profile-fit replay similarly propagates rounding through its polynomial basis
conversion and differentiated residuals. Adaptive inverse iterations and their
work totals are checked for internal consistency and algorithmic bounds in each
run; floating-point stopping decisions need not match across platforms. Historical
budget receipts use the saved costs, while replay charges its actual work.
These are replay allowances, not physical uncertainty or convergence bounds.
One-ULP conditioning controls and deliberate larger corruptions test the
acceptance and rejection sides.

The [validation successor record](storage/validation-source-successors.json)
identifies the original Git bytes and the reviewed current checker bytes.
Old result manifests remain intact, current validation code runs, and newly
created records bind current code. Unknown source changes and all data-hash
changes still fail. Verified original Git bytes are cached within a process;
current files, registry entries and checker hashes are checked on every call.
The shared synthetic-test budgets count only CPU time
inside their tests, retaining cumulative point and active-CPU limits.

When a publication replay fails on another platform, the checker reports
the rejected numeric leaves or separates wall-coordinate drift from internal
forcing-grid drift. These diagnostics read the failed comparison's operands;
they do not rerun a producer or change acceptance. Investigate those values
before modifying a replay allowance.
Conditioned comparisons retain the actual rejected subtree for that report,
so accepted conditioned differences cannot obscure a later failing value.

The Verify workflow lets an active full-suite run finish while subsequent
commits are saved. It keeps the latest pending run in the same concurrency
group; intermediate pending runs can be replaced, as documented by
[GitHub's concurrency rules](https://docs.github.com/en/actions/how-tos/write-workflows/choose-when-workflows-run/control-workflow-concurrency).
A successful older run is evidence for that commit, not for a newer queued one.

## Current work

The [shadow-radial comparison](nsc-atlas-finite-amplitude.md#separating-radial-quadrature-from-reconstructed-field-transport)
has exact controls for its local identity, wall conversion, differentiated
panel-error sum and moving-interface signs:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_mass_transport_identity.py -q -s
```

These controls share 20 active process-CPU seconds and take no numerical
field samples. They do not evaluate a retained cap or authorize a new
trajectory. The numerical shadow-flow/error evaluator is still to be built.
Its selected reference uses the smooth analytic potential; the affine panel
component and cap composition are now implemented for manufactured controls.

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_radial_reference.py -q -s
```

These ten manufactured controls compare the potential and its small implementation
gap with 120-digit arithmetic, check exact moving affine maps and an independent
fixed-radius variation, and exercise the born-panel and budget failure paths.
They share 200,000 charged values / 20 active CPU seconds and use no retained
cap. The returned ODE tolerances and Taylor remainder estimates are separate
quantities; neither asserts a physical or uniform integration-error bound.
Development uses 9,356 charged values / 0.03125 active CPU seconds.

The read-only field adapter supplies fixed-radius and radial derivatives
of the same reconstructed fields, with explicit owner-cell traces:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_transport_fields.py -q -s
```

Its 27 manufactured controls pass. Cumulative development uses
18,526 charged values / 0.609375 active CPU seconds of the
200,000-point / 20-second allowance. No retained cap or radial solve is used.
Native samples exclude axis jets and exact formula-switch queries. Reference
samples provide explicitly named small-radius limits, and a separate method
returns regular axis coefficients without claiming axis field jets.

The cap composer retains the producer's native panels and splits only the
reference at the fixed formula switch. Local map errors and their derivatives
are summed independently of the final mass/rate subtraction:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_radial_comparison.py -q -s
```

Seven controls check an exact constant-field cap, a nonconstant reference
against independent fixed-radius variation, moving switch coincidences,
frozen partitions and invalid native evidence. Cumulative development uses
472,439 charged values / 6.625 active CPU seconds of the 1,000,000-point /
40-second allowance. A native born panel exactly on a formula switch remains
outside the existing tangent's domain. No retained cap has been evaluated;
the lapse, residual/jump integrals and integration-error control remain open.

The next component integrates mass, lapse and the three field-residual channels
on one declared smooth piece, with no input from the native mass tangent:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_transport_reference.py -q -s
```

Five manufactured checks pass using 9,828 charged values / 0.046875 active CPU
seconds of a 200,000-point / 20-second allowance. Exact polynomial integrals
and a 60-digit independent mass integral check the states and attenuation.
The complete interface assembly remains to be implemented. The expanded seven
transport-algebra controls also check wall/centre coincidences and separate
interface velocities; that run uses zero field samples / 0.8125 CPU seconds.
Neither component evaluates the retained cap or provides an error enclosure.

The full cap assembler now connects these components and every ordered
interface event, while keeping the affine and residual-ODE reference solutions
independent:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_transport_balance.py -q -s
```

Nine manufactured controls pass. Cumulative development, including failed
and corrected attempts, uses 890,431 charged values / 12.328125 active CPU
seconds of a 2,000,000-point / 60-second allowance. One coupled-rate test uses
the original RHS and discrete ledger; the new channel sum differs by about
`8.2e-16`. The full-balance wall remains off the native formula switches, and
the centre lapse is independently checked against the frozen field. Reported
reference mass/rate/jump discrepancies still require a declared accuracy
budget before interpreting a retained-state result. That next fixed terminal
comparison is specified in the calculation owner; it has not run yet.

The interface algebra uses no numerical fields or root search:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_interface_identity.py -q -s
```

All three exact controls pass in 1.71875 active CPU seconds within their
20-second allowance. They establish primitive continuity, formula-switch
orders and the necessity of the moving potential-jump term.

The discrete cap mass derivative and departure-flux ledger are implemented
in a diagnostic companion to the unchanged producer. Run their controls with:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_mass_derivatives.py tests/test_nsc_atlas_cap_mass_tangent.py tests/test_nsc_atlas_cap_mass_balance.py -q -s
```

The [mathematical owner](nsc-atlas-finite-amplitude.md#discrete-mass-derivative-and-departure-flux-balance)
states the branch-specific derivative and telescoping identity. The controls
have passed. The [completed four-case result](current-result.md#discrete-cap-mass-balance-16-october-2026)
accounts for the upwind contribution and retains the remaining discrepancy.
The [calculation record](prior-art-and-open-claim.md#discrete-cap-mass-balance-16-october-2026)
preserves the failed rounding check and its bounded diagnosis.

Full payloads belong in ignored `.build/results`. Creation refuses to
overwrite an existing output. Reconstruct a missing payload, or validate
an existing one, using the frozen producer and input bindings:

```sh
python scripts/run_python.py scripts/diagnose_nsc_atlas_cap_mass_balance.py --output .build/results/nsc-atlas-cap-mass-balance-v1.json
python scripts/run_python.py scripts/diagnose_nsc_atlas_cap_mass_balance.py --output .build/results/nsc-atlas-cap-mass-balance-v1.json --check
python scripts/run_python.py scripts/summarize_result.py --output results/summaries/nsc-atlas-cap-mass-balance-v1.json --check
```

The last command checks hashes and structure without a scientific replay.
The compact summary contains the exact regeneration argument vector.

The [first prediction-register entry](prediction-register.md#nsc-p001--perturbative-cap-and-moving-wall)
binds the closed perturbative cap–wall relation and four computed
coefficients to their original evidence. Authenticate its projections without
repeating a calculation:

```sh
python scripts/run_python.py scripts/check_prediction_register.py
```

The [terminal cap mass-functional derivative](current-result.md#terminal-cap-mass-functional-derivative-16-october-2026)
uses the [analyzer](../scripts/analyze_nsc_atlas_cap_mass_direction.py) and
[recorder](../scripts/diagnose_nsc_atlas_cap_mass_direction.py) on the authenticated
terminal cap states. Its 26 mathematical, orchestration and provenance controls
pass with cumulative usage of 171,957 points and 1.34375 process CPU seconds,
including the initial epsilon-tolerance failure and its correction:

```sh
python scripts/run_python.py -m pytest tests/test_analyze_nsc_atlas_cap_mass_direction.py tests/test_nsc_atlas_cap_mass_direction_orchestration.py tests/test_diagnose_nsc_atlas_cap_mass_direction.py -q
```

The [compact summary](../results/summaries/nsc-atlas-cap-mass-direction-v1.json)
uses 1,711,610 points and 6.6875 process CPU seconds including creation
validation. One standalone check passes with 855,805 points and 4.078125
CPU seconds. The 520,080-byte record has SHA256
`d0fec7048ab1cb476e15672ccf334c2c8e628c059c1fdec709e702c9dafd6f74`.
The original producer was frozen at `fc71709`. Full data remains local in
ignored `.build/results/`. A clean checkout can regenerate it when needed:

```sh
python scripts/run_python.py scripts/diagnose_nsc_atlas_cap_mass_direction.py --output .build/results/nsc-atlas-cap-mass-direction-v1.json
```

Creation refuses to overwrite an existing record. Append `--check` to check
an existing raw record. For an inexpensive summary/input authentication that
does not repeat the numerical diagnostic, use:

```sh
python scripts/run_python.py scripts/summarize_result.py --output results/summaries/nsc-atlas-cap-mass-direction-v1.json --check
python scripts/run_python.py scripts/check_result_storage.py
```

The first command reports whether the raw file is present; its absence on
a clean checkout is expected. [Storage policy](nsc-storage.md#current-calculation-outputs-16-october-2026)
keeps future arrays local and their summaries, hashes and commands in Git.

Each invocation pins the prior cap audit, verifies its seven producer bindings
within the twelve audit bindings, and authenticates the 104 inherited bindings
once. Each analysis uses 20 cap metrics and four cap RHS evaluations at the
terminal cut, holding the directions fixed for backward probes. Creation
performs one analysis and one validation under a shared 60-CPU-second,
four-million-point allowance; a check performs one analysis under 30 CPU
seconds/two million points. No source, full stage or trajectory is replayed.
The three control suites share 30 CPU seconds/300,000 points. Fresh creation
requires a new output file without `--check`; first failures leave an ignored
`.build` receipt without retry.

The finest order-eight finite directional residual is `-0.6612504470`, compared
with `-0.6407403055` from saved history expressed with the same rates; its
order-eight minus order-four change is `1.09497e-5`. These are qualified
backward functional diagnostics, without a convergence-order claim. At zero
amplitude the normalized residual is not a physical gravitational defect.

The [selected cap mass and wall-curve audit](current-result.md#selected-cap-mass-and-wall-curve-audit-16-october-2026)
uses the [recorder](../scripts/diagnose_nsc_atlas_cap_consistency.py) and
[analyzer](../scripts/analyze_nsc_atlas_cap_consistency.py) on the fixed
operator-comparison histories. Its 34 manufactured and provenance controls
pass on their first attempt with 62,119 shared points and 0.34375 process
CPU seconds:

```sh
python scripts/run_python.py -m pytest tests/test_analyze_nsc_atlas_cap_consistency.py tests/test_nsc_atlas_cap_consistency_orchestration.py tests/test_diagnose_nsc_atlas_cap_consistency.py -q
```

The [completed record](../results/development/nsc-atlas-cap-consistency-v1.json)
uses 4,393,156 points and 15.359375 process CPU seconds including creation
validation. One standalone check passes with 2,196,578 points and 8.453125
CPU seconds. The 820,093-byte record has SHA256
`79d1e02917ba7c1a7d3326cd80612b0467f19b0b5b7b1e69b6e1e2f1f3ae9728`.
To authenticate it:

```sh
python scripts/run_python.py scripts/diagnose_nsc_atlas_cap_consistency.py --output results/development/nsc-atlas-cap-consistency-v1.json --check
```

Each invocation authenticates the 104 inherited bindings once and checks
seven audit source/control bindings. Each analysis reconstructs 64 cap
metrics and 12 selected exterior prefixes and hybrid wall maps, comparing
quadrature orders 4/8 and radial/time widths one/two at cuts 63, 68 and 128.
Creation performs one analysis and one validation under a shared 120-CPU-second,
eight-million-point allowance; a standalone check performs one analysis under
60 CPU seconds/four million points. Neither invokes a source solve, full
coupled stage, cap RHS or trajectory replay. The three control suites share
30 CPU seconds/one million points. A fresh creation requires a new output
file without `--check`; the first failure leaves an ignored `.build` receipt
with its cut, order and phase when available, without retry or completed output.

Pulse residuals remain strongly sensitive to radial derivative width despite
negligible order-4/8 changes. At the final finite cut the temporal mass
residual is about `-8.32e-5`, while radial mass and lapse residuals are below
`1.7e-9` and `9.3e-9`. The actual wall-curve acceleration also differs from
the hybrid prediction. The terminal mass-functional calculation above now
tests the cap transport connection; neither audit supplies a new wall law or
a new evolution.

The [near-wall saved-row audit](current-result.md#saved-near-wall-row-audit-16-october-2026)
uses the [bounded diagnostic](../scripts/diagnose_nsc_atlas_boundary_rows.py)
and [eight synthetic controls](../tests/test_diagnose_nsc_atlas_boundary_rows.py).
Controls pass with 60,460 aggregate points/0.15625 CPU seconds. Creation
including validation uses 97,350 points/1.59375 CPU seconds; its one
standalone check uses 48,675 points/1.5625 CPU seconds. Each authenticates
the 104 input bindings once and evaluates three fixed cuts, without a
stage or trajectory. Both production allowances are 30 CPU seconds and
150,000 points. To authenticate the retained result:

```sh
python scripts/run_python.py scripts/diagnose_nsc_atlas_boundary_rows.py --output results/development/nsc-atlas-boundary-rows-audit-v1.json --check
```

Fresh reproduction selects a new output file without `--check`; existing
files are protected. The result reports stencil sensitivity and leaves
the origin of the boundary shape unresolved. It does not impose a wall
derivative; the selected cap audit above now evaluates the independent
mass and wall-curve consistency questions.

The [same-grid exterior-operator comparison](prior-art-and-open-claim.md#fixed-grid-exterior-operator-comparison-16-october-2026)
uses the [quadratic derivative](../src/recursive_horizons/nsc_atlas_departure_gradient.py),
[explicit coupled stage](../src/recursive_horizons/nsc_atlas_gradient_evolution.py)
and [callback SSPRK3 driver](../src/recursive_horizons/nsc_atlas_stage_ssprk3.py).
Its [recorder](../scripts/record_nsc_atlas_operator_comparison.py) and
[analyzer](../scripts/analyze_nsc_atlas_operator_comparison.py) reuse the
authenticated 64-cell baseline and unchanged source/gauge cache. Focused
controls are:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_departure_gradient.py tests/test_nsc_atlas_gradient_evolution.py tests/test_nsc_atlas_stage_ssprk3.py tests/test_record_nsc_atlas_operator_comparison.py tests/test_analyze_nsc_atlas_operator_comparison.py -q
```

All 89 controls pass before physical creation, reporting 2,047,112 points
and 3.84375 process CPU seconds across their admitted shared allowances.
The [completed comparison](current-result.md#fixed-grid-exterior-operator-comparison-16-october-2026)
uses 55,398,966 total points/256.25 process CPU seconds. Its one standalone
check passes with 4,040,624 points; the [manifest](../results/development/nsc-atlas-operator-comparison-v1/manifest.json)
binds 104 paths. To authenticate the retained result:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_operator_comparison.py --output results/development/nsc-atlas-operator-comparison-v1 --check
```

A fresh reproduction uses a new output directory without `--check`.
The numerical comparison uses 64 exterior/cap
cells, 512 time steps and both original amplitudes. It shares 450 CPU
seconds/80 million points including creation validation; its standalone
check has 60 CPU seconds/eight million points and no trajectory replay.
Old baseline arrays are referenced rather than duplicated. Same-grid cap
comparisons retain total physical cell integrals without adjacent summation.

The [completed fixed-input spatial comparison](current-result.md#fixed-input-spatial-refinement-16-october-2026)
uses the [recorder](../scripts/record_nsc_atlas_spatial_refinement.py) and
[saved-state analyzer](../scripts/analyze_nsc_atlas_spatial_refinement.py).
It adds zero and finite-amplitude 64-cell windows to the retained 32-cell
comparison, holding the 512 SSPRK3 steps, 1,025 exact cache rows, source
snapshot and computational profile fixed. Run its focused controls or
authenticate the [saved record](../results/development/nsc-atlas-spatial-refinement-v1/spatial.json):

```sh
python scripts/run_python.py -m pytest tests/test_analyze_nsc_atlas_spatial_refinement.py tests/test_record_nsc_atlas_spatial_refinement.py -q
python scripts/run_python.py scripts/record_nsc_atlas_spatial_refinement.py --output results/development/nsc-atlas-spatial-refinement-v1 --check
```

The original creation and validation use 60,087,926 points and 256.8125
process CPU seconds. One independent standalone check passes with
6,684,624 authentication points. The [manifest](../results/development/nsc-atlas-spatial-refinement-v1/manifest.json)
binds 92 files. Creation requires a fresh directory and shares 450 process
CPU seconds and 80,000,000 points across authentication, both trajectories, saved diagnostics
and creation validation. A standalone check has 60 CPU seconds/8,000,000
points; it reconstructs saved physical fields, final stages and diagnostics,
without source integration, cache reconstruction or trajectory replay. The
two new control suites together have 60 CPU seconds/2,000,000 points.
Failures leave an ignored `.build` receipt and no completed canonical
record. Comparisons sum adjacent total physical cap integrals and distinguish
native from restricted derivative stencils, retaining the nonzero source
constraint. Two spatial grids do not establish a convergence order or error
enclosure. The paired wall-scalar spatial difference is 0.774698% of the
fine response peak, while the paired cap `Iq` difference is about 7% of
the old 32-cell paired peak at common normalized cells. Wall curvature and
quiet-profile support remain unresolved. A fresh numerical repeat uses the same producer with
`--output` naming a fresh directory and without `--check`; the completed
record is never overwritten.

The [incoming wall-constraint audit](current-result.md#incoming-wall-curvature-and-transported-constraint-16-october-2026)
uses pinned histories, four one-sided spatial stencils and causal birth
interpolation. Authenticate the completed [record](../results/development/nsc-atlas-wall-constraint-audit-v1.json)
without replaying its source or trajectory:

```sh
python scripts/run_python.py -m pytest tests/test_diagnose_nsc_atlas_wall_constraint.py -q
python scripts/run_python.py scripts/diagnose_nsc_atlas_wall_constraint.py --output results/development/nsc-atlas-wall-constraint-audit-v1.json --check
```

Creation requires a fresh output path; it never replaces the completed
record. Creation and standalone checking each have a shared allowance of
30 CPU seconds and 2,000,000 points. Creation uses 1,391,719 points and
1.046875 CPU seconds; one independent check passes with the same point
count. Both use 129 analytic vacuum evaluations and no source, profile,
cache, coupled-stage or time-evolution replay. The 19 controls share
30 CPU seconds/200,000 points and pass with 71,773 points and 0.0625 CPU
seconds. All four wall subtractions worsen the disjoint-subset RMS;
spatial/transport accuracy remains unresolved, with no junction-failure
or error-enclosure claim.

The [constraint transport audit](current-result.md#saved-constraint-transport-and-sampling-16-october-2026)
uses only pinned saved arrays, separating the inherited source residual,
diagnostic time-stencil sensitivity and incoming-characteristic origin:

```sh
python scripts/run_python.py -m pytest tests/test_diagnose_nsc_atlas_constraint_transport.py -q
python scripts/run_python.py scripts/diagnose_nsc_atlas_constraint_transport.py --output results/development/nsc-atlas-constraint-transport-audit-v1.json
python scripts/run_python.py scripts/diagnose_nsc_atlas_constraint_transport.py --output results/development/nsc-atlas-constraint-transport-audit-v1.json --check
```

Its 21 controls pass with aggregate 72,316 points and 0.359375 CPU seconds.
The completed audit uses 689,270 points and 0.921875 CPU seconds. Creation
and a standalone check each have 30 CPU seconds/1,000,000 points; controls
share 30 CPU seconds/300,000 points. Creation requires a fresh path and
checking repeats only saved arithmetic and hash authentication, with no
source/cache/profile/stage or time-evolution replay. The producer, controls
and input bindings remain fixed with the result.
One standalone check passes with 689,270 points.

The current cap representation carries an exact characteristic background and
evolves departures. Run its 72 manufactured, endpoint and independent radial controls:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_characteristic_background.py tests/test_nsc_atlas_background_cap.py tests/test_nsc_atlas_background_profile.py tests/test_nsc_atlas_background_evolution.py tests/test_nsc_atlas_background_endpoint.py -q
```

The profile suite has one 30-CPU-second/20,000-point allowance; the cap
suite has 30 CPU seconds and 100,000 points for positive partial quadrature
and an independent narrow-pulse order-4/8/16 comparison. The input adapter
has 30 CPU seconds/20,000 points; the coupled stage has 30 CPU
seconds/100,000 points. The endpoint suite has 30 CPU seconds/20,000 points;
it accepts the exact terminal profile label and rejects an excursion beyond it.
The
[owner](nsc-atlas-finite-amplitude.md#carrying-a-characteristic-background-and-evolving-its-correction)
specifies the background's finite domain, centre limits and free incoming
correction. These tests do not replay the prepared physical source.

The [six-case recorder](../scripts/record_nsc_atlas_background_evolution.py)
keeps the raw cap departures and reconstructed physical integrals separately.
Its controls check profile fitting versus withheld rows, physical reconstruction,
paired comparisons, shared accounting, immutable output and failed-attempt receipts:

```sh
python scripts/run_python.py -m pytest tests/test_record_nsc_atlas_background_evolution.py -q
```

The completed [first physical comparison](current-result.md#characteristic-background-coupled-evolution-16-october-2026) uses 228.9375 CPU seconds and 43,551,872 points. Its shared allowance is 600 CPU seconds and
100 million evaluated points, including authentication. It reuses the bound
source cache and extracts the computational profile once. Its [calculation
record](prior-art-and-open-claim.md#coupled-characteristic-background-comparison-16-october-2026)
fixes six cases and stopping conditions. Creation uses a new output directory;
`--check` authenticates saved data without replaying time evolution:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_background_evolution.py --output results/development/nsc-atlas-background-evolution-v1
python scripts/run_python.py scripts/record_nsc_atlas_background_evolution.py --output results/development/nsc-atlas-background-evolution-v1 --check
```

The [four-panel comparison](figures/nsc-atlas-background-evolution.png) is
drawn only from saved arrays and authenticated input hashes. Its original
render uses 1.15625 CPU seconds and 36,051 arithmetic points. To reproduce
it into a fresh local output (the renderer never overwrites an existing file):

```sh
python scripts/run_python.py scripts/plot_nsc_atlas_background_evolution.py --output .build/nsc-atlas-background-evolution.png
```

The [scalar-time diagnostic](../scripts/diagnose_nsc_atlas_background_time.py)
uses the retained reference profile and four frozen stages. Its 35 controls
exercise independent quadrature and chain identities, hash-only input
authentication and one shared allowance. It performs no time evolution:

```sh
python scripts/run_python.py -m pytest tests/test_diagnose_nsc_atlas_background_time.py -q
python scripts/run_python.py scripts/diagnose_nsc_atlas_background_time.py --output results/development/nsc-atlas-background-time-audit-v1.json
python scripts/run_python.py scripts/diagnose_nsc_atlas_background_time.py --output results/development/nsc-atlas-background-time-audit-v1.json --check
```

Each audit has 30 CPU seconds/150,000 points; controls have 30 CPU
seconds/100,000 points. Creation requires a new output path. Checking
repeats the same bounded diagnostic and excludes elapsed CPU time from
the numerical equality comparison.

The [saved audit](../results/development/nsc-atlas-background-time-audit-v1.json)
completes with 57,152 points and 1.0625 CPU seconds; its standalone check
passes. It exposes substantial cancellation and does not isolate the
accumulated scalar time error. The subsequent [local accepted-step comparison](prior-art-and-open-claim.md#actual-scalar-step-balance-16-october-2026)
uses identical saved starting states and the existing source/gauge cache:

```sh
python scripts/run_python.py -m pytest tests/test_diagnose_nsc_atlas_background_step.py -q
python scripts/run_python.py scripts/diagnose_nsc_atlas_background_step.py --output results/development/nsc-atlas-background-step-audit-v1.json
python scripts/run_python.py scripts/diagnose_nsc_atlas_background_step.py --output results/development/nsc-atlas-background-step-audit-v1.json --check
```

Its 22 controls use 5,312 points and 0.125 CPU seconds. The physical
diagnostic uses 138,891 points and 1.546875 CPU seconds, within its
30-second/250,000-point allowance. It takes six local SSPRK2 steps with
twelve coupled stages and saves an exact five-part scalar balance; it
replays no complete trajectory. Creation is exclusive and `--check`
repeats only those same local steps. The [outcome](current-result.md#actual-scalar-step-balance-16-october-2026)
retains the distinction between local increments and full-history error.

The [SSPRK3 operator](../src/recursive_horizons/nsc_atlas_background_ssprk3.py)
uses the existing coupled stage and exact supplied midpoint cache keys.
Its 20 controls take 27,607 points and 0.5625 CPU seconds, plus two
intentional budget-exhaustion charges. The 17 local-recorder controls use
3,808 points and 0.140625 CPU seconds:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_background_ssprk3.py tests/test_diagnose_nsc_atlas_ssprk3_step.py -q
python scripts/run_python.py scripts/diagnose_nsc_atlas_ssprk3_step.py --output results/development/nsc-atlas-ssprk3-step-audit-v1.json
python scripts/run_python.py scripts/diagnose_nsc_atlas_ssprk3_step.py --output results/development/nsc-atlas-ssprk3-step-audit-v1.json --check
```

The [saved local result](../results/development/nsc-atlas-ssprk3-step-audit-v1.json)
uses six stages, 70,118 points and 1.1875 CPU seconds, within its
30-second/150,000-point allowance. It reads the pinned earlier local
SSPRK2 record without replaying those steps. The weighted five-part
ledger, common starting state and stage keys are verified independently.
Creation requires a new path; `--check` repeats only the six new stages.
The [full original-window comparison](prior-art-and-open-claim.md#fixed-input-ssprk3-coupled-window-16-october-2026)
uses four fixed 128/256-step cases and no new source. Its recorder's 27
controls pass with 37,821 points and 1.796875 CPU seconds:

```sh
python scripts/run_python.py -m pytest tests/test_record_nsc_atlas_ssprk3_evolution.py -q
python scripts/run_python.py scripts/record_nsc_atlas_ssprk3_evolution.py --output results/development/nsc-atlas-ssprk3-evolution-v1
python scripts/run_python.py scripts/record_nsc_atlas_ssprk3_evolution.py --output results/development/nsc-atlas-ssprk3-evolution-v1 --check
```

Creation has one 300-CPU-second/50-million-point allowance for four
evolutions, diagnostics and authentication together. The standalone
60-second/5-million-point check reconstructs saved fields, four final
stages and diagnostics without a trajectory or source replay. Creation
requires a fresh directory; failures retain their first context under
ignored `.build`. The fixed model and 80 bound paths remain unchanged.

The [completed four-case comparison](current-result.md#fixed-input-ssprk3-coupled-window-16-october-2026)
uses 130.734375 CPU seconds and 24,939,592 points, including validation.
Standalone checking passes with 956,575 points. Its scalar time difference
remains 21.54% of the fine response peak; closer cross-method agreement
does not replace a temporal error estimate. The saved `diagnostics`
contain `finite_minus_zero`, `paired_time`, `paired_versus_ssprk2`,
constraint histories and pre-support controls. The old SSPRK2 paired
histories supply the pointwise comparison `(4*paired_512-paired_256)/3`;
compare this array with the new paired 256-step history without a solve.

The [midpoint-complete comparison](prior-art-and-open-claim.md#midpoint-complete-ssprk3-comparison-16-october-2026)
adds only the missing cache labels and two 512-step cases. The original
source snapshot, profile and spatial operators stay fixed. Its 17 cache
and 19 recorder controls pass with 540,731 points and 2.890625 CPU seconds,
plus an intentional 513-point exhaustion attempt:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_refined_cache.py tests/test_record_nsc_atlas_ssprk3_refinement.py -q
python scripts/run_python.py scripts/record_nsc_atlas_ssprk3_refinement.py --output results/development/nsc-atlas-ssprk3-refinement-v1
python scripts/run_python.py scripts/record_nsc_atlas_ssprk3_refinement.py --output results/development/nsc-atlas-ssprk3-refinement-v1 --check
```

Creation shares 300 CPU seconds/50 million points across cache evaluation,
gauge-order comparison, evolution and validation. A standalone check has
60 CPU seconds/5 million points; it reproduces the midpoint cache and two
final radial stages, with no trajectory or source-ODE replay. The 86
input/source bindings remain frozen. Creation requires a new directory.

The [completed 512-step pair](current-result.md#midpoint-complete-ssprk3-comparison-16-october-2026)
uses 174.90625 CPU seconds and 33,604,170 points, including validation.
The standalone check uses 1,418,330 points. Its cache construction uses
575,112 points, including both gauge orders; all original 513 rows are
preserved. The paired scalar temporal difference is `5.84100e-7`, or
1.67616% of the fine peak. `diagnostics` retains full signed comparison
arrays and maximum locations, quiet channels and profile-knot brackets.

The [three-resolution figure](figures/nsc-atlas-ssprk3-refinement.png)
uses only authenticated saved arrays. Its original render uses 21,801
arithmetic points and 1.09375 CPU seconds within a 30-second/200,000-point
allowance. Reproduce it at a fresh path:

```sh
python scripts/run_python.py scripts/plot_nsc_atlas_ssprk3_refinement.py --output .build/nsc-atlas-ssprk3-refinement.png
```

The [compatible transport audit](../scripts/diagnose_nsc_atlas_compatible_transport.py)
uses the fixed [compatible result](../results/development/nsc-atlas-compatible-evolution-v1/compatible.json)
without replaying its source or time evolution. Its 25 controls exercise
the exact incoming cap pullback, a closed polynomial prefix oracle,
derivative/integration separation, changed-input rejection and the shared
allowance:

```sh
python scripts/run_python.py -m pytest tests/test_diagnose_nsc_atlas_compatible_transport.py -q
```

Authenticate only the pinned input payload, manifest and its 51 source
bindings, without reconstructing a scientific state:

```sh
python scripts/run_python.py -c "import sys; sys.path.insert(0, 'scripts'); import diagnose_nsc_atlas_compatible_transport as audit; print(audit.authenticate()[2])"
```

Reproduce the saved-state audit into a file that does not exist:

```sh
python scripts/run_python.py scripts/diagnose_nsc_atlas_compatible_transport.py --output .build/compatible-transport-audit-review.json
```

Omitting `--output` prints the JSON without saving it. Each invocation has
one cumulative 30-process-CPU-second/100,000-point allowance, including
the before/after input hash checks, two zero-amplitude cap pullbacks and
four frozen finite-amplitude stages. Each stage reconstructs its cap once,
then crosses two departure derivative rules with two prefix quadratures;
it does not integrate a new trajectory. The
[completed audit](../results/development/nsc-atlas-compatible-transport-audit-v1.json)
uses 48,910 points and 0.859375 CPU seconds. Its input identity, audit
source hashes, runtime and counts are saved in the JSON. The
[outcome](current-result.md#compatible-transport-audit-16-october-2026)
distinguishes cap-profile smearing and derivative sensitivity from
quadrature sensitivity, without treating the postprocessed alternatives
as evolved solutions or error bounds.

The compatible successor evolves exact-vacuum geometric departures and
conserved cap integrals. Its single common scalar wall value comes from
the cap gradient primitive. Run its manufactured and analytic controls:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_vacuum_departure.py tests/test_nsc_atlas_reduced_cap.py tests/test_nsc_atlas_compatible_evolution.py tests/test_nsc_atlas_departure_diagnostics.py -q
```

The 94 controls include actual incoming-wave propagation, moving vacuum
preservation, finite-amplitude wall conversion and independent pure-null
constraint violations. Each regional/control suite shares its declared
30-CPU-second/20,000-point budget. The
[owner](nsc-atlas-finite-amplitude.md#conservative-reduced-variables-and-exact-vacuum-departures)
specifies the integral weights, regular-axis parity and background
subtraction. These controls do not replay the prepared physical source.

The [successor recorder](../scripts/record_nsc_atlas_compatible_evolution.py)
authenticates and reuses the first coupled run's source and reference caches.
Its [completed five-case comparison](../results/development/nsc-atlas-compatible-evolution-v1/compatible.json)
has the same grids, time steps and amplitude, with one shared
180-CPU-second/8,000,000-point allowance including inherited authentication
and final checks. Run the recorder controls and authenticate the canonical
result without a physical replay:

```sh
python scripts/run_python.py -m pytest tests/test_record_nsc_atlas_compatible_evolution.py -q
python scripts/run_python.py scripts/record_nsc_atlas_compatible_evolution.py --check --output results/development/nsc-atlas-compatible-evolution-v1
```

A fresh reproduction requires a directory that does not exist:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_compatible_evolution.py --output .build/compatible-evolution-review
```

Checking reconstructs raw-to-physical saved fields, final cap constraints,
wall/prefix stages and independent diagnostics. It performs no source,
reference or time evolution. The raw finite-difference and analytic-vacuum
departure diagnostics remain distinct. Fine cap integrals are restricted
by summing adjacent cells. The first failure is saved under the ignored
`.build/compatible-evolution-failures` area with no completed result
directory, automatic retry or change of numerical settings.

The saved creation and validation total 64.40625 process CPU seconds and
5,371,737 charged points. Standalone authentication uses 1,202,991 points;
it reuses the inherited source/reference caches and evolves no trajectory.
All five cases retain 129 slices. Exact vacuum and zero-outgoing recovery,
paired spatial responses, raw finite-case time changes and remaining cap
transport/constraint limits are reported in the
[current result](current-result.md#compatible-coupled-evolution-16-october-2026).

The simultaneous exterior, wall and polar-cap stage uses a source snapshot,
ordered exterior integrals, the constrained regular cap and the existing
hybrid wall law. Its fixed reference wall label evolves the actual wall
lapse through the acceleration equation. Run the analytic and manufactured
controls without a new physical source integration:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_characteristic_trace.py tests/test_nsc_atlas_exterior_characteristic_slice.py tests/test_nsc_atlas_polar_cap.py tests/test_nsc_atlas_reference_wall_gauge.py tests/test_nsc_atlas_coupled_evolution.py -q
```

The [scientific owner](nsc-atlas-finite-amplitude.md#a-simultaneous-exterior-wall-and-cap-stage)
specifies the coordinate, axis, characteristic and averaging conventions.
The first-order spatial scheme and fixed-step SSPRK2 retain measured
residuals; they do not supply a continuum error certificate.

The [coupled recorder](../scripts/record_nsc_atlas_coupled_evolution.py)
completes five declared cases, with two spatial resolutions and one
separate time refinement. The independent saved-history diagnostics and
original-rectangle overlap have their own manufactured controls:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_coupled_diagnostics.py tests/test_record_nsc_atlas_coupled_evolution.py -q
python scripts/run_python.py scripts/record_nsc_atlas_coupled_evolution.py --check --output results/development/nsc-atlas-coupled-evolution-v1
```

For a fresh reproduction, use a directory that does not exist:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_coupled_evolution.py --output .build/coupled-evolution-review
```

Creation uses one cumulative 180-process-CPU-second/8,000,000-evaluation
budget including authentication, with the original bounded source children.
The record separates calculation and verification charges. It retains
129 accepted slices per case, exact-index cached source and reference
inputs, final radial/wall data and independent comparisons. Source hashes
are checked before and after calculation. A first failure is saved under
the ignored `.build/coupled-evolution-failures` area, without creating a
completed result directory or changing the declared grid.

Standalone checking performs no source integration, reference quadrature
or time evolution. It authenticates the original rectangle, reconstructs
the source interpolant and five final radial-cap/prefix states, and
recomputes the saved-array diagnostics. The creation-only
[result](../results/development/nsc-atlas-coupled-evolution-v1/coupled.json)
retains unresolved spatial/interface and cap field-compatibility errors;
see the [current outcome](current-result.md#first-coupled-exterior-wall-and-cap-window-16-october-2026).

The [nonlinear rectangle marcher](../src/recursive_horizons/nsc_atlas_nonlinear_rectangle.py)
assembles the exact source and quiet bottom in the original affine null
rectangle. Its independent nodal diagnostics compute pure null constraints
and geometric mass without using the solved mixed equations as derivatives.
The [four-grid recorder](../scripts/record_nsc_atlas_nonlinear_rectangle.py)
compares zero amplitude and half the source-only bound on 16/32 cells per
direction. Run its implementation controls and authenticate the saved result:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_nonlinear_rectangle.py tests/test_record_nsc_atlas_nonlinear_rectangle.py -q
python scripts/run_python.py scripts/record_nsc_atlas_nonlinear_rectangle.py --check --output results/development/nsc-atlas-nonlinear-rectangle-v1
```

For a fresh reproduction, select a directory that does not exist:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_nonlinear_rectangle.py --output .build/nonlinear-rectangle-review
```

Creation has one 120-process-CPU-second, 60,000-evaluation budget, including
initialization and all four marches. Source children retain their own
30-second/20,000-evaluation limits and debit the parent. Boundary queries
finish before marching. Failure records its phase and any failed cell under
`.build/nonlinear-rectangle-failures`, without writing a completed result
directory. Authentication checks source/output hashes, settings, boundaries,
counts and recomputed diagnostics with zero source integrations or cell
solves. The [owner](nsc-atlas-finite-amplitude.md#original-exterior-rectangle-and-independent-diagnostics)
and [admission](prior-art-and-open-claim.md#nonlinear-exterior-rectangle-16-october-2026)
define the finite scope; sampled differences do not certify a continuum
error bound or the later moving-wall and cap evolution.

The [exact source integrator](../src/recursive_horizons/nsc_atlas_nonlinear_source.py)
retains the original compact preparation and quiet corner. Its amplitude
control comes from a continuous source-only bound. The
[characteristic-cell solver](../src/recursive_horizons/nsc_atlas_characteristic_cell.py)
advances one positive-radius mixed cell with a bounded Newton iteration.
Their independent controls run without a prepared-source regional march:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_nonlinear_source.py tests/test_nsc_atlas_characteristic_cell.py -q -s
```

The source tolerance test prints its bounded comparison receipt. Source
queries retain the shared CPU/evaluation budget and reject extrapolation.
Cell diagnostics report the nonlinear residual, correction and iteration
counts; they do not certify the regional Einstein constraints or continuum
error. Both implementations and their physical domains are specified in
the [same nonlinear owner](nsc-atlas-finite-amplitude.md).

The [normalized bulk operator](../src/recursive_horizons/nsc_atlas_nonlinear_bulk.py)
and [hybrid wall map](../src/recursive_horizons/nsc_atlas_nonlinear_boundary.py)
connect the exact regional equations to a cap whose lapse is determined by
its central clock. Their [owner](nsc-atlas-finite-amplitude.md) derives the
regular axis system and identifies the source, moving-grid and constraint
requirements for a coupled evolution. Run the independent canonical-action,
vacuum, coordinate and mixed-chart controls together:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_nonlinear_bulk.py tests/test_nsc_atlas_nonlinear_boundary.py tests/test_nsc_atlas_hybrid_boundary.py -q
```

These pointwise helpers do not evolve a spacetime. In particular, bulk mass
rates labelled `on_constraints` require the null constraints; first jets
alone cannot measure those residuals. The hybrid map requires a constrained
untrapped cap and a positive incoming radial sum `H`.

The [exact boundary and central-clock helpers](../src/recursive_horizons/nsc_atlas_nonlinear_boundary.py)
evaluate the same-action nonlinear matching algebra without evolving a
spacetime. Their [owning derivation](nsc-atlas-finite-amplitude.md)
specifies compatible gauge inputs and the remaining coupled problem.
Run their independent vacuum, canonical-action and coordinate controls:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_nonlinear_boundary.py -q
```

These tests create no physical evolution record and select no amplitude
interval for the full model.

The [profile-moment successor](nsc-atlas-profile-moment.md) compares the
pre-arrival mass identity against the original radial evaluator, then
compares local lifts on 16/32 active cells with 40 withheld queries:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_profile_moment.py --output .build/profile-moment-review
python scripts/run_python.py scripts/record_nsc_atlas_profile_moment.py --check --output results/development/nsc-atlas-profile-moment-v1
python scripts/run_python.py scripts/record_nsc_atlas_local_lift.py --output .build/local-lift-review
python scripts/run_python.py scripts/record_nsc_atlas_local_lift.py --check --output results/development/nsc-atlas-local-lift-v1
```

Each complete physical comparison has one shared 120-process-CPU-second,
30-million-kernel budget. The mass pair uses moment orders 4/8 and the
original radial/ray order 16 reference; the local lifts share 81 direct
moment-order-8 boundary jets. Both retain the same 129-node field and
257-knot geometry. Creation refuses an existing output directory. The
local-lift `--check` authenticates its predecessor and reconstructs only
the saved polynomials; it performs no physical solve. Numerical outcomes
are in the [current result](current-result.md#profile-moment-and-local-lift-16-october-2026).

The [conservative evolution implementation](nsc-atlas-scalar-evolution.md)
adds coefficient caching and the forced fixed-coordinate scalar equation.
Its manufactured tests do not themselves supply a physical waveform.

The first forcing-cache comparison shares 65 direct slices between 17/33
training times and 32 withheld times. The latter include radial midpoints
to distinguish time-only from combined time/radial discrepancies:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_forcing_cache.py --output .build/forcing-cache-review
python scripts/run_python.py scripts/record_nsc_atlas_forcing_cache.py --check --output results/development/nsc-atlas-forcing-cache-v1
```

Its 120-process-CPU-second, 30-million-kernel budget includes preparation.
Authentication reconstructs the interpolants from retained slices with no
physical fallback. A cached evaluation inside the scalar solver performs
no new field or metric integration.

The boundary-component diagnosis uses only the authenticated 81 saved
jets and reconstructs eight polynomial fits:

```sh
python scripts/run_python.py scripts/assess_nsc_atlas_lift_components.py --output .build/lift-components-review.json
python scripts/run_python.py scripts/assess_nsc_atlas_lift_components.py --check --output results/development/nsc-atlas-lift-components-v1.json
```

Its 30-process-CPU-second budget excludes physical solves entirely. The
record binds its producer and input, and `--check` repeats the inexpensive
component arithmetic and withheld comparisons.

The clock refinement reuses those outgoing cap jets, recomputes wall clocks
with 16/32 panels, and keeps the same boundary knots and holdouts:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_lift_clock.py --output .build/lift-clock-review
python scripts/run_python.py scripts/record_nsc_atlas_lift_clock.py --check --output results/development/nsc-atlas-lift-clock-v1
```

One field preparation and two wall solves share 120 process-CPU seconds
and 30 million kernels. The check restores saved first jets and clock
polynomials without new physical solves. The 32-panel result meets the
original sampled boundary target; the [current result](current-result.md#clock-resolved-scalar-boundary-16-october-2026)
records the numerical comparison and subsequent propagation status.

The first exterior scalar propagation restores these authenticated inputs
and runs four cases separating cache time, correction grid and ODE tolerance:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_scalar_evolution.py --output .build/scalar-evolution-review
python scripts/run_python.py scripts/record_nsc_atlas_scalar_evolution.py --check --output results/development/nsc-atlas-scalar-evolution-v1
```

The complete comparison shares 120 process-CPU seconds and 30 million
kernels, with at most 6,000 RHS calls per solve. It performs no new
leading-field, wall, forcing-slice or boundary solve. Authentication restores
the coefficient polynomials and reconstructs 68 discrete wall rates from
retained interior arrays, without replaying propagation. The
[result](current-result.md#exterior-scalar-correction-16-october-2026) records
separate sensitivities and the limit of the formal coefficient.

The central reconstruction reuses the earlier cap clock and current
exterior/wall records, and compares new reflected incoming-ray integrals
at radial orders 16/32 with ray order 16:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_central_scalar.py --output .build/central-scalar-review
python scripts/run_python.py scripts/record_nsc_atlas_central_scalar.py --check --output results/development/nsc-atlas-central-scalar-v1
```

The five common interior arrival cuts and one leading-field preparation
share 120 process-CPU seconds and 30 million kernels. It performs no new
wall, cap-clock or exterior evolution. The check authenticates saved inputs,
quadrature inventory and reconstruction algebra without repeating incoming
integration. The [central result](current-result.md#central-scalar-coefficient-16-october-2026)
retains all terms, the quiet control and separate numerical sensitivities.
Its independent incoming-ray successor retains that saved radial-32/ray-16
case and computes only radial-32/ray-32 at the same five cuts:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_central_ray.py --output .build/central-ray-review
python scripts/run_python.py scripts/record_nsc_atlas_central_ray.py --check --output results/development/nsc-atlas-central-ray-v1
```

One fresh leading preparation and five incoming integrals share a new
120-process-CPU-second, 30-million-kernel budget. Prior costs are recorded
separately. Its check authenticates the reused coarse side and new
reconstruction algebra without physical replay.

The [scalar-input successor](nsc-atlas-scalar-forcing.md) compares a 9/17-knot
Hermite boundary lift against 16 withheld queries, or cumulative exterior
forcing with radial orders 2/4 at nine times:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_scalar_inputs.py --comparison lift --output .build/scalar-input-lift-review
python scripts/run_python.py scripts/record_nsc_atlas_scalar_inputs.py --comparison forcing --output .build/scalar-input-forcing-review
```

Each complete pair shares a 120-process-CPU-second, 30-million-kernel
budget, including the one field and one wall preparation. Both hold 129
field nodes and 257 geometry knots fixed. Existing wall settings and cap
radial/ray orders 16 remain unchanged. The boundary fits share their 33
unique direct queries; withheld points never become fitted knots in that
comparison. The interpolation target comes from the saved field sensitivity
and is a numerical diagnostic. Exceeding it is recorded as an unresolved
discretization, rather than silently accepting the fit.

Creation refuses an existing output directory. Use `--check --output PATH`
to authenticate a saved `inputs.json` and its manifest without a physical
solve. These inputs do not yet propagate the scalar correction.

The saved pairs are under `results/development/nsc-atlas-scalar-inputs-v1`,
in `lift` and `forcing`. Both use current source bindings. The
[result entry](current-result.md#scalar-lift-and-exterior-forcing-16-october-2026)
reports the unresolved coarse interpolation and the separate exterior
quadrature/transport discrepancies.

The subsequent saved-data assessment changes local Hermite stencil size
while keeping the original training and withheld samples fixed:

```sh
python scripts/run_python.py scripts/assess_nsc_atlas_scalar_lift.py --output .build/scalar-lift-assessment-review
python scripts/run_python.py scripts/assess_nsc_atlas_scalar_lift.py --check --output results/development/nsc-atlas-scalar-lift-assessment-v1
```

This calculation uses no field or wall solve. Its separate creation-only
record binds the predecessor evidence and the assessment producer; `--check`
repeats only the inexpensive interpolation arithmetic. The results motivate
additional samples near the active pulse edges, rather than higher degree
on the unchanged knots.

The [scalar backscatter calculation](nsc-atlas-scalar-metric-response.md)
evaluates the outgoing correction and displaced exterior boundary derivative
from the same prepared field, cap metric and wall solution:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_scalar_boundary.py --comparison field --output .build/scalar-boundary-review
python scripts/run_python.py scripts/record_nsc_atlas_scalar_boundary.py --check --output results/development/nsc-atlas-scalar-boundary-v1/field
```

Its five comparisons are `field` (129/257), `radial` (8/16), `ray` (8/16),
`radial_fine` (16/32) and `ray_fine` (16/32). The finer radial/ray controls
hold 129 field nodes fixed; the other quadrature controls use 257. Geometry
knots remain 257, and cap radial/ray baseline orders are 16/16. Upstream
wall settings remain radial order 2, field-derivative order 32, 8 clock
panels, source order 64 and ODE tolerance $2\times10^{-9}$. The cap field
keeps derivative order 2, centre quadrature 8 and threshold $10^{-4}$.

Each pair shares one 120-process-CPU-second, 30-million-kernel budget,
including field/wall preparation and both independent integral coordinates.
The producer preserves existing output directories and binds source bytes
before and after the calculation. All five saved records use current
sources; `--check` performs zero physical solves. Nine samples, signed
contributions and domain margins document a boundary input, not a boundary
interpolant or propagated waveform. The [current result](current-result.md#scalar-backscatter-boundary-16-october-2026)
states all measured numerical sensitivities.

The [smooth cap calculation](nsc-atlas-cap-response.md) reconstructs the
interior metric and central null-arrival correction from the same prepared
field and wall solution. Its two clock coefficients keep proper wall time
and exterior null label fixed, respectively. Create a comparison, then
authenticate its saved settings, source bindings and bytes:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_cap_response.py --comparison field --output .build/cap-review
python scripts/run_python.py scripts/record_nsc_atlas_cap_response.py --check --output .build/cap-review
```

Each complete comparison shares one 120 process-CPU-second and
30-million-scalar-kernel budget across field preparation, wall evolution,
cap reconstruction, rays and diagnostics. Authentication performs no
physical solve. The baseline uses field/geometry 257/257, cap radial/ray
orders 16/16 and 16 cap time panels; upstream wall radial order 2 and 8
time panels retain the existing field-derivative order 32, source order 64
and ODE tolerance $2\times10^{-9}$.

| Comparison | Changed setting | Fixed field nodes |
| --- | --- | --- |
| `field` | 129 to 257 field nodes | Varied |
| `radial` | Cap radial order 8 to 16 | 257 |
| `ray` | Ray order 8 to 16 | 257 |
| `time` | Cap time panels 8 to 16 | 257 |
| `geometry` | Geometry knots 129 to 257 | 129 |
| `derivative` | Continuum second-derivative order 2 to 4 | 129 |
| `centre_order` | Small-radius quadrature order 8 to 16 | 129 |
| `centre_radius` | Small-radius threshold $10^{-4}$ to $10^{-3}$ | 129 |
| `ray_fine` | Ray order 16 to 32 | 129 |
| `radial_fine` | Cap radial order 16 to 32 | 129 |

The field-129 overrides are explicit in each record's `comparison_base`;
all other settings stay fixed within each pair. The default derivative
order is 2, small-radius order 8 and threshold $10^{-4}$. These comparisons
measure numerical sensitivity of formal metric and null-clock coefficients.
They do not compute the higher-order scalar amplitude or a finite-amplitude
error enclosure. Signed ray integrals and extrema refer to the quadrature
samples, with no claim of continuum extrema.

In the `time` comparison, the prepared field, wall and radial/ray rules are
identical. Only the wall-time coefficient's integration panels change. The
second solution reuses the first solution's nine reference-ray integrals;
`ray_nodes` counts fresh work and separate reused-ray counts preserve that
distinction. The shared budget still includes the first ray calculation
and both wall-time integrations.

The retained `field`, `radial`, `ray`, `centre_radius` and `ray_fine` records
keep the original source context `19983a7b4730d1fdc8e59ae4ba2f3efd560913a6`.
Their bytes are unchanged after the declaration checker gained a
$10^{-13}$ absolute tolerance for geometry-endpoint roundoff. Authenticate
those five records with an explicit context, for example:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_cap_response.py --check --source-context 19983a7b4730d1fdc8e59ae4ba2f3efd560913a6 --output results/development/nsc-atlas-cap-response-v1/field
```

This reads the pinned source bytes from a full Git checkout, then checks
the saved record with the current structural validator. It does not run
the historical producer. Missing Git objects fail explicitly; there is no
substitution of current sources. The `geometry`, `derivative` and
`centre_order` records similarly retain context
`5e58c9fd2866f8adcc107ee22f13f389e7f15c77`; pass that full hash to
`--source-context` for their authentication. The `time` and `radial_fine`
records use current sources. The publication checker selects all three
contexts by explicit comparison name, and authentication reports the
selected context and zero cap solves.

The [moving-wall calculation](nsc-atlas-wall-evolution.md) extends the same
prepared field to its retained wedge and evaluates two mass paths, wall
displacement and the exterior null clock. Create or authenticate a comparison:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_wall_evolution.py --comparison field --output .build/wall-review
python scripts/run_python.py scripts/record_nsc_atlas_wall_evolution.py --check --output .build/wall-review
```

Its independent comparison names are `field`, `radial`, `derivative`, `time`,
`geometry`, `source` and `ode`. Each two-setting comparison includes its field
solves and queries in one 120 process-CPU-second, 30-million-kernel budget.
The saved summaries are formal quadratic coefficients and numerical
refinement diagnostics. The smooth-cap calculation above supplies the next
metric and central-clock connection.

The independent radial/mixed overlap uses the same field on the old rectangle
and a different gravitational integration path:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_wall_overlap.py --output .build/wall-overlap-review
python scripts/run_python.py scripts/record_nsc_atlas_wall_overlap.py --check --output results/development/nsc-atlas-wall-overlap-v1
```

It records nine common cuts, mass and lapse differences, the mixed solver's
own constraint diagnostics, and the cumulative compute budget. Authentication
performs no field, metric or wall evolution.

The [smooth exterior metric](nsc-atlas-bulk-response.md) has separate
finite comparisons for metric grid, leading field, geometry interpolation,
field derivative quadrature, source quadrature and time tolerance. For
example, create a fresh calculation or authenticate an existing one with:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_bulk_response.py --comparison metric --output .build/bulk-review
python scripts/run_python.py scripts/record_nsc_atlas_bulk_response.py --check --output .build/bulk-review
```

`--check` checks the recorded sources, runtime declaration, configuration and
output integrity; it does not repeat the metric solve. The six comparison
names are `metric`, `field`, `geometry`, `derivative`, `source` and `time`.
Each complete two-solve comparison has its own cumulative 120 process-CPU
second and 30-million-scalar-kernel ceiling. The extended prepared-field
interface leaves the original v8 evidence untouched; replay that earlier
record explicitly with `scripts/replay_history.py atlas-v8`.

The six retained bulk-v1 records were produced at `3470152e1766`, before
the wedge API was added. Authenticate them with their original recorder,
or explicitly repeat all six original numerical studies:

```sh
python scripts/run_python.py scripts/replay_history.py check-bulk-v1
python scripts/run_python.py scripts/replay_history.py bulk-v1
```

The first command checks retained bytes and declarations without a metric
solve. The second requires the recorded runtime and uses each original
comparison's compute budget. Both select the named source context explicitly;
they never replace a failing current-source check with historical code.

```sh
python -m pip install -e '.[paper,dev]'
python scripts/check_lab_snapshot.py
python scripts/check_publication.py
python scripts/check_claim_bindings.py
python -m pytest tests/
python scripts/build_paper.py --check
```

`make verify` runs current checks, every retained test and the deterministic
paper check. CI uses the same pinned stack on Ubuntu and Windows. Historical
replay is an explicit separate operation, selected through `make history-check`
or the commands below.

The history-evaluator test records CPU time as a diagnostic. Set
`NSC_ENFORCE_TEST_CPU_BUDGETS=1` to enforce its original 30-second reference
on a comparable, idle host, running that module serially. Selecting only the
timing test does not run the work it measures. Operational calculation budgets
remain enforced; parallel pytest workers do not provide comparable timings.
Both this control and the geometry replay carry the `cpu_budget` marker.
The geometry replay retains its producer's 20-second cap. For a parallel run,
`-m 'not cpu_budget'` optionally excludes those two whole tests, including the
geometry replay's numerical comparisons; the default serial suite includes them.
The [14 October Win11 follow-up](win11-verification-20261014.md) lists the
11 requested reruns. The [subsequent receipt update](receipt-rebinding-20261014.md)
records the metadata regeneration that resolved its four binding failures.

`check_lab_snapshot.py` authenticates the current paths, hashes and sizes in the
v2 [retained-file receipts](lab-snapshot.json), plus the declared object archive
and unchanged original pins in [source history](../archive/source-history/manifest.json).
The [history compaction record](repository-consolidation.md#explicit-historical-replay)
describes the sparse authentication transport retained after pruning. It
preserves the bytes required by current checks, without keeping full parent
ancestry. Complete old checkouts and campaigns require the Zenodo archive.
`check_publication.py` checks retained result identities, publication routing
and explicit DOI/hash receipts for archived evidence. It does not download or
verify a remote archive. Current checks do not require the retired reassembly
manifest or migration scripts; their carriers are cited at
[Zenodo](https://doi.org/10.5281/zenodo.23090015).

Immutable records can bind earlier source bytes. Authentication then uses the
record's explicitly declared Git context, preserving its original numerical
claims and hashes. A pass authenticates that recorded calculation; it does not
identify its original source with the current implementation. Tests of active
APIs use this root package. Original-context authentication and replay are
never a silent fallback from a failed current-source check.

The threshold18 import has two normalized creator path labels. Its metadata
checker authenticates the exact old and current bytes and the single string
replacement, while numerical replay uses current local files. It records the
original and current tree hashes separately; any other change is rejected.

The bounded prepared-history receipt update has a separate exact
[ledger](../results/metadata/prepared-history/receipt-rebind.json). Verify its
current carriers with
`python scripts/run_python.py scripts/rebind_nsc_prepared_history_receipts.py --verify`.
The original execution signatures and numerical members remain fixed.

When changing the manuscript's figures, the retained figure producers are:

```sh
python scripts/build_geometry_overview.py
python scripts/build_characteristic_initialization_figure.py
python scripts/build_source_mass_figure.py
```

Run the paper check after figure changes and inspect the result visually.

The compact analytic successor authenticates the original v7/v8 evidence and
the proof of the fixed-model preparation inverse. Its checker records no new
waveform or numerical inverse:

```sh
python scripts/record_nsc_atlas_observer_inference.py --check
```

The [inference owner](nsc-atlas-observer-inference.md) supplies the mathematical
argument. The checker authenticates that argument's bytes and declaration;
passing the checker is not a machine proof. Creation requires a fresh output
directory, for example `--output .build/observer-inference-review`, and refuses
an existing directory. Its source bindings cover the proof, producer and tests.
Its pinned predecessor manifests and all stored output hashes are checked
without importing historical producers or substituting current source bytes.

The separate [scalar wall check](nsc-atlas-scalar-wall-consistency.md) records
the canonical probe's singular surface stress and regular tangential flux.
It authenticates the fixed preparation and includes exact symbolic controls;
it does not rerun the waveform or evolve the geometry:

```sh
python scripts/record_nsc_atlas_scalar_wall_consistency.py --check
python scripts/run_python.py scripts/record_nsc_atlas_scalar_wall_consistency.py --check-algebra
python scripts/run_python.py -m pytest tests/test_nsc_atlas_scalar_wall_consistency.py
```

Its owner, producer and tests are source-bound. As with the inverse, creating
another record requires a fresh output directory. Independent analytic review
and exact-algebra checks are distinct from provenance authentication. Bulk
backreaction, normal-force balance and amplitude bounds remain open.

The [characteristic initialization](nsc-atlas-characteristic-initialization.md)
evaluates one borrowed null constraint for the original compact source.
Its corner geometry and affine gauge are fixed; it runs no field or coupled
Einstein evolution:

```sh
python scripts/record_nsc_atlas_characteristic_initialization.py --check
python scripts/run_python.py scripts/record_nsc_atlas_characteristic_initialization.py --check-numerics
python scripts/run_python.py -m pytest tests/test_nsc_atlas_characteristic_initialization.py
```

Authentication checks the source-bound owner, implementation, tests and
prerequisites. Numerical replay separately checks exact identities and
the finite one-dimensional quadratures in the recorded budget. Creation
requires a fresh output directory. Refinement agreement is a numerical
diagnostic, not a rigorous enclosure or a finite-amplitude metric error bound.
That intrinsic constraint alone leaves the transverse geometry open.

The [source-leg mass successor](nsc-atlas-source-mass.md) closes the
transverse radius equation on the same initial null surface. It compares
the geometric mass using the full frame transformation, mixed geometry
and improved-stress transport:

```sh
python scripts/record_nsc_atlas_source_mass.py --check
python scripts/run_python.py scripts/record_nsc_atlas_source_mass.py --check-numerics
python scripts/run_python.py -m pytest tests/test_nsc_atlas_source_mass.py
```

Its numerical replay checks exact identities, independent transverse-field
integrals and separate mass/geometry refinements within the recorded finite
budget. Ordinary authentication performs no numerical integration. The
geometric mass is invariant under coordinate changes, but generally differs
between conformal frames during the input. The negative post-input
coefficient is relative to the fixed reference; it does not select a finite
amplitude or establish a vacuum future. Later wall matching, bulk evolution
and the gravitational correction to the observer signal remain open.

The [normal-force successor](nsc-atlas-scalar-normal-force.md) checks the
remaining regular scalar stress jump on the fixed wall. Its ordinary
check authenticates the new declaration and original scalar-wall evidence;
explicit algebra replay independently tests the projection and force sign:

```sh
python scripts/record_nsc_atlas_scalar_normal_force.py --check
python scripts/run_python.py scripts/record_nsc_atlas_scalar_normal_force.py --check-algebra
python scripts/run_python.py -m pytest tests/test_nsc_atlas_scalar_normal_force.py
```

Creation requires a fresh output directory and stays inside the declared
60-process-CPU-second budget. It runs no waveform, pressure history or
coupled metric evolution. Its necessary mean-curvature condition does not
set the wall displacement to zero; individual one-sided pressures and
finite-amplitude validity remain uncomputed.

The [wall-mass closure](nsc-atlas-wall-mass-closure.md) expresses the
exterior mass coefficient using the already fixed leading wall trace.
It checks both null fluxes, full frame conversion and original-frame stress
transport analytically, without a numerical mass curve:

```sh
python scripts/record_nsc_atlas_wall_mass_closure.py --check
python scripts/run_python.py scripts/record_nsc_atlas_wall_mass_closure.py --check-algebra
python scripts/run_python.py -m pytest tests/test_nsc_atlas_wall_mass_closure.py
```

Ordinary checking authenticates stored evidence; explicit algebra replay
checks the identities and omission controls. Creation requires a new output
directory within the 60-process-CPU-second budget. A numerical evaluation
needs its own quadrature, field/geometry error control and admission. The
existing waveform samples do not certify this new integral. Coupled cap/wall
geometry and the corrected observer signal were left to the successor below.

The [cap/wall response](nsc-atlas-cap-wall-response.md) derives the
vacuum-subtracted cap mass, wall-displacement functional and compatible
temporal/normal conditions. It records the common model specification with
the two-sided result:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_cap_wall_response.py --check
python scripts/run_python.py scripts/record_nsc_atlas_cap_wall_response.py --check-algebra
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_wall_response.py tests/test_record_nsc_atlas_cap_wall_response.py tests/test_nsc_cosmological_atlas.py
```

Ordinary checking authenticates the stored identities, current sources and
unchanged exterior evidence. Explicit algebra replay also repeats the exact
stress, frame and junction calculations. Creating a record uses a fresh
directory and a 60-process-CPU-second limit. These commands integrate no
mass, trajectory or bulk metric. `field_mass_ingredients()` exposes the
first-derivative integrands and endpoints; `wall_response_coefficients()`
exposes the regular displacement equation. The corrected observer signal
requires the next bulk radius/lapse and field-response calculation.

Atlas v8 remains the finite exterior-to-interior waveform result. Its original
replay is explicitly selected below. The fixed input, numerical budget,
refinement diagnostics and independent analytic tail bound remain in the
[response owner](nsc-atlas-prepared-response.md). Refinement diagnostics are not
a rigorous numerical enclosure. The [v8 guide](nsc-cosmological-atlas-v8.md)
describes the outputs and causal limits.

The paper builder renders local figures, authenticates their bytes and compares
two builds in check mode. Visual inspection is still required for a manuscript
change. The software release version and the manuscript revision are independent.

Revision `working-2026-10-16.5` has passed visual inspection of the cover,
contents, abstract, new Section 13.18, quantitative table, conclusion and
evidence appendix. Its 37-page PDF is 2,810,596 bytes with SHA256
`5d28ba600c527daa00e8497203faf9aa001b17d0fc283c90191bed093fb95aaf`.
Two clean builds match the tracked PDF byte for byte. The current paper
and catalog bindings pass, and the source has 32 claims with 112 resolved
evidence links. The saved near-wall audit is documented in its research
owner; it does not change the manuscript's fixed-grid numerical result.

The canonical PDF uses the pinned Zopfli compressor for ReportLab streams, so
standard zlib and zlib-ng installations produce the same encoded document.
This controls compression only; fonts, images, equations and layout retain their
bound inputs. The builder checks its rendering dependencies and keeps strict
byte equality. A failed comparison saves diagnostic PDFs and runtime details
under ignored `.build/paper-check`; the paper workflow preserves those files.

## Preserved editions and source contexts

Routine authentication uses the retained receipts and explicitly pinned Git
objects. Historical publication tests run the original result contracts against
baseline `04a0cb9378977772243ebf7c9be3c74182b030f2`, without restoring deleted
producers into the active package. The original result manifest and numerical
claims remain unchanged. [Preservation and source history](repository-consolidation.md)
describes these receipts; retired migration and result carriers are cited at
[Zenodo](https://doi.org/10.5281/zenodo.23090015).

Scientific recomputation is separate. `replay_history.py` selects an original
commit and materializes its source layout under ignored `.build/replay`.
Its unchanged checkers enforce the original evidence hashes and comparison
policies. A replay requiring an archived carrier needs its exact hash-bound
bytes; missing evidence or source objects fail explicitly. Edited current
owners are never substituted for historical source bytes.

### Atlas v1–v6 and supporting editions

The v1–v6 manifests record Python 3.14.8 and Matplotlib 3.11.2. Their original
image checks require byte equality in the recorded Windows environment;
matching package versions alone does not guarantee identical images on another
machine. [requirements-history.txt](../requirements-history.txt) preserves
these scientific pins and supplies the replay tools. Pillow 12.3.0 is a replay
dependency; Pillow and pytest versions were not recorded in the atlas manifests.
For a deliberate original-context replay, use an isolated environment without
installing the current project package there. After history compaction, restore
the exact archived context from Zenodo before these full-context commands;
the sparse authentication store alone does not supply the complete old layout:

```powershell
py -3.14 -m venv .venv/history
.\.venv\history\Scripts\python.exe -m pip install -r requirements-history.txt
.\.venv\history\Scripts\python.exe scripts/replay_history.py prepare
.\.venv\history\Scripts\python.exe scripts/replay_history.py all
```

`prepare` selects the baseline context without scientific recomputation. `all`
recomputes the six atlas bundles within their recorded finite domains and
budgets, authenticates the finite and OPEN companions, and runs the original
publication tests. It does not restart their full research campaigns. The
individual selectors are `atlas`, `finite`, `companion` and `tests`.

With Make installed, the optional historical command is:

```sh
make history-check HISTORY_PYTHON=.venv/history/Scripts/python.exe
```

`HISTORY_PYTHON` defaults to `PYTHON`; it never selects another interpreter after
a failure. `make verify` checks the current project and does not run historical
replay. Current CI runs the retained suite on Linux and Windows with Matplotlib
3.11.2; the former separate historical rendering jobs are retired.

### Atlas v7: original root context

Atlas v7 uses commit `7fee894ad5a4da2ba378f724e45e926401f82670` under
`.build/replay/7fee894ad5a4`:

```sh
python scripts/replay_history.py prepare-v7
python scripts/replay_history.py atlas-v7
```

`prepare-v7` authenticates the original evidence and sources without a field
solve or renderer. `atlas-v7` runs the unchanged original checker, including
exact figure comparisons. That manifest requires Windows, Python 3.14.8,
NumPy 2.5.1, SciPy 1.17.1, SymPy 1.14.0 and Matplotlib 3.11.1. This immutable
historical requirement differs from the current 3.11.2 pin; the current
validation environment is not an exact-image replay environment for v7.

### Atlas v8: original prepared-response context

Atlas v8 uses commit `1154baf290d901e73058392e143d0acb37d22951` under
`.build/replay/1154baf290d9`:

```sh
python scripts/replay_history.py prepare-v8
python scripts/replay_history.py atlas-v8
```

`prepare-v8` authenticates the original source bindings and retained evidence
without running a field solve. `atlas-v8` runs the unchanged original producer
within its finite numerical budget, including exact image comparisons. It
requires the same recorded Windows/Python 3.14.8/Matplotlib 3.11.1 profile as v7.
Use these full replay selectors only with that original environment; the
current renderer does not promise historical PNG identity.

Neither v7 nor v8 is added to legacy `all`. Changed preserved outputs, missing
source objects and incompatible runtimes fail explicitly. Original image
comparisons retain byte equality; numerical comparisons retain their recorded
tolerances. On compatible x86-64 Windows hosts, the established replay setting
is `$env:OPENBLAS_CORETYPE='Haswell'`. CPU dispatch and BLAS configuration matter
in addition to version pins, and a later passing run does not excuse a mismatch.

The retired `make reproduce`, `make reproduce-development` and
`make reproduce-exact` routes stop with guidance. A new campaign needs a named
question, fresh admission and a finite budget. Its producer must run in an
explicitly selected source context; pinned hashes are never changed to fit the
current package.

## Windows runtime

The existing validation environment uses Python 3.14.8 and the project's pinned
scientific stack. If its moved Windows redirector cannot locate `pyvenv.cfg`, use
`python scripts/run_python.py <script-or-python-arguments>` with the base interpreter;
for example, `python scripts/run_python.py -m pytest tests/`. Do not install a
second project package or restore the old directory. Normal installations use
their activated interpreter. Root `scripts/lab.py` is only a compatibility launcher:
it now runs in this root and does not switch scientific packages.

## Scope

The main [paper](../paper/nested-space-cosmology.md) is a living atlas manuscript.
The [paper guide](papers.md) routes archived and supporting editions. Spherical
symmetry, specified states and fixed parameters are assumptions of individual
tests. No new constants, propagation or recurrence campaign was admitted during
repository reassembly. New calculations use the [existing admission rule](prior-art-and-open-claim.md).
