# Reproducing the root project

There is one package, `nested-space-cosmology`, importing `recursive_horizons`
from root `src`. Use the pinned versions in [pyproject.toml](../pyproject.toml).
Scientific conclusions retain the domains in their owning notes.
The current paper and test dependencies use Matplotlib 3.11.2; the scientific
stack includes python-flint 0.9.0 and Numba 0.68.0.

## Current work

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
