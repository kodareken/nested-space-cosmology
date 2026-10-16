# Prior art and admission of new research

## Nonlinear exterior rectangle, 16 October 2026

Assemble the exact source and implicit characteristic cells into the original
exterior rectangle. This is the first regional evolution of the nonlinear
equations for the fixed preparation. The source and quiet bottom determine
this rectangle without future wall data: a compatible future timelike wall
has increasing `Ubar` and lies beyond its quiet-corner value `Ubar_b` while
the chart remains regular. Continuation of that wall and the cap remains the
next coupled problem.

| Admission item | Fixed decision |
| --- | --- |
| Borrowed and connected | The exact source, normalized mixed equations and constraint transport in the [finite-amplitude owner](nsc-atlas-finite-amplitude.md), together with the tested implicit-midpoint cell. No new action, preparation or free boundary data. |
| Missing connection | Propagation of the prescribed finite-amplitude source across the rectangle, with independently differenced null constraints and a measured discretization change. |
| Fixed domain | `Ubar=CJ0*(U-U0)`, `CJ0=-3.2/U0`; `0<=Ubar<=CJ0*(Ua-U0)`, `.4<=x<=x_last`. Use `source_geometry(64)` once and retain its original support fronts. |
| Fixed comparisons | Four uniform grids: zero amplitude and `alpha=alpha_star/2`, each with 16 and 32 cells per direction. The normalized scalar remains present at zero amplitude. Compute raw state/constraint changes and finite-minus-zero differences separately. |
| Input accuracy | Source DOP853 `rtol=2e-10,atol=1e-12`, default maximum step `source_length/32`, split at the original support edges. Sample both sources and quiet bottoms on the 33-node grids before marching; coarse boundaries are their even-index subsets. |
| Cell settings | Existing absolute and relative tolerances `1e-12`, reference scales `(1,1,1)`, at most 12 Newton updates and 8 backtracking trials, maximum scaled condition number `1e12`. Fixed state guard `r in [.2,.8]`, `log(C) in [-2,2]`, `c in [-4,4]`. These are engineering guards, not an invariant-tube theorem. |
| Cumulative budget | One parent budget of 120 process CPU seconds and 60,000 charged evaluations, including reference geometry, source construction, queries and all four grids. Each source uses a child of the existing Budget with at most 30 CPU seconds and 20,000 evaluations, debiting the parent for every charge. Neither budget resets during the calculation. Record each grid's RHS/Jacobian counts and source charges. |
| Independent controls | Before the original-source run, use manufactured polynomial constraint violations and exact shifted Schwarzschild grids; check a zero-field march against analytic geometry. Manufactured implementation tests have an aggregate 30 CPU-second/20,000-cell-charge ceiling. In the physical run compare the zero-amplitude metric with exact Schwarzschild, retain quiet-field recovery, and compare source transverse derivatives with the integrated first jets, including the nonzero post-pulse scalar derivative. |
| Diagnostics | Compute first and pure second derivatives from the nodal fields, never from the solved mixed equation. Retain raw constraints and derivative-defined mass. Compare common interior nodes under refinement and report all-interior, front-near and away-from-front summaries; boundary derivative defects are separate. Small cell residuals are not constraint convergence. |
| Outcome | A source-bound numerical rectangle and measured refinement, or a precise failed cell, tube or budget diagnosis. A useful result supplies data for subsequent wall-born cells; it does not establish cap evolution, a central observation, stability or a uniform finite-amplitude error enclosure. |
| Stop or redirect | Stop on the first budget, iteration or tube failure with its cell indices, coordinates, input corners and cumulative accounting. Do not enlarge the tube, alter the amplitude, subdivide or restart silently. Preserve a useful diagnosis before designing any changed calculation. |

The recorder is creation-only. It binds its producer and inputs before and
after calculation, saves the four nodal arrays and exact boundaries, and
authenticates them without reintegrating a source or solving a cell. A failed
run leaves no completed result directory; its diagnostic is kept in the
ignored build area until a concise scientific diagnosis is incorporated here.
Two grid sizes measure a difference; they do not by themselves establish an
asymptotic order or a rigorous error bound. No new leading-field solve is
included merely to supply another comparator.

**Outcome:** all four grids complete, with 2,560 successful cells, 23,237
charged evaluations and 5.71875 process CPU seconds including validation.
The 31 new implementation/authentication tests and independent reviews
pass. The [saved result](../results/development/nsc-atlas-nonlinear-rectangle-v1/rectangle.json)
authenticates without a source integration or cell solve. Vacuum radius and
lapse errors decrease about fourfold. On common interior nodes the raw
finite-amplitude `E_x` maximum decreases from `0.2074424` to `0.0591010`;
the same finite-difference diagnostic on the exact source decreases from
`0.2086466` to `0.0594446`. Its analytic source constraint remains below
`1.06e-15`. The original source's steep derivatives therefore require
explicit treatment in extracting outgoing jets. Retain the exact source
and interpolate/differentiate the propagated increment for the next
rectangle-to-wall connection, carrying its measured input sensitivity.
There is no additional rectangle scan or amplitude change in this outcome.

## Exact source and characteristic-cell implementation, 16 October 2026

Use the [source first-jet system](nsc-atlas-finite-amplitude.md#completing-the-source-first-jets)
and the implemented normalized mixed operator to construct the first
nonlinear characteristic input. Preserve the original `C3` reduced-field
pulse, source geometry, quiet Schwarzschild bottom and corner mass. Integrate
the scaled Jordan radius changes and the three transverse jets in the affine
source gauge; compute profile derivatives analytically. A bounded dense
evaluation supplies first jets without differencing saved samples.

For initialization controls, use the continuous profile derivative bounds
and source-only continuation threshold derived in the owner, with
`r_min=1/5`, `v_max=3/2`, `F_min=1/2`. The deterministic nonzero test value
is half that threshold. It selects an initial-leg test, not a full-spacetime
amplitude range or a physical calibration. Check zero-amplitude recovery,
the transformed null constraint, quiet intervals, corner mass, source
orientation, bounded interpolation and numerical refinement. Each test has
30 process CPU seconds; each source construction has at most 20,000 counted
reference/profile evaluations, including its initialization when applicable.
No new physical trajectory record is implied by a passing unit control.

Implement a local implicit-midpoint characteristic cell using the same bulk
operator. Its nonlinear iteration has an explicit finite count, positive
radius/lapse guards and a residual check. Independent Minkowski and analytic
vacuum cell refinements precede any prepared-source regional march. A failed
iteration or non-refining result is diagnosed before changing resolution.
The first prepared-source regional calculation requires its own specified
finite rectangle and shared budget. Wall-born cells, the cap and the central
readout remain the intended continuation of this same solver.

**Implementation outcome:** 38 source controls pass in 1.71875 process CPU
seconds, and 35 local-cell controls pass in 0.59375 process CPU seconds.
Independent reviews confirm both implementations and their source-gauge
connection. The source's half-threshold control is
`alpha=0.00012983127142551443`; its two-tolerance comparison uses 3,063
shared evaluations including geometry and queries. Its maximum radius and
mass changes are `1.44e-14` and `1.71e-13`. These changes measure tolerance
sensitivity, while the reported analytical constraint residual checks the
frame identity. Curved exact-corner cell errors approach the expected
quartic local scaling under two balanced halvings. No prepared-source
regional march or full-spacetime amplitude bound is supplied by these tests.

## Exact regional operators, 16 October 2026

The [finite-amplitude owner](nsc-atlas-finite-amplitude.md) supplies the
Einstein-frame action, mixed equations and null constraints. Connect the
implemented wall traces to regional evolution by implementing the normalized
positive-radius bulk operator and deriving a regular polar-areal cap system.
Keep the same fixed cap/exterior potentials, source and material wall. The
bulk helper accepts finite first jets with positive radius and lapse; it
does not reconstruct missing pure second derivatives or solve a gauge.
Return constraint-surface mass rates with their residual qualification.
Preserve the continuous zero-amplitude scalar limit and evaluate the
potential increment without subtracting nearly equal vacuum terms.

Independent canonical-field conversion, exact vacuum geometries, constraint
residuals and coordinate rescalings supply the controls. Each test has a
30-process-CPU-second ceiling. No physical evolution, amplitude choice or
new source record is admitted here. The cap derivation must identify regular
axis variables, the central clock, incoming characteristic data and its
conversion to the existing wall map before implementation. An inconsistency
in constraints, characteristic count or gauge is resolved in that connection
before a coupled physical run is designed.

The independently checked polar-cap connection now provides that boundary
count: reconstruct the cap lapse from the central clock, prescribe only
the compatible exterior lapse trace, and require `A_cap>0` and
`H=R_plus+S*r>0`. Angular matching determines the wall speed; one incoming
cap scalar characteristic completes transmission. Implement this mixed-chart
boundary map alongside the existing two-null-chart map. Check its vacuum
branch, canonical Lorentz-frame transmission, common field derivative,
exterior coordinate rescaling and the noncircular constraint implication.
The same 30-process-CPU-second test ceiling and no-physical-evolution scope
apply. The actual cap fields and their radial constraint reconstruction
remain caller inputs, not extra freely chosen boundary functions.

**Implementation outcome:** 38 bulk controls and 32 hybrid-wall controls,
together with the 32 existing boundary/clock controls, pass in one combined
run: 102 tests, 0.57 seconds of pytest-reported duration and 0.734375 process
CPU seconds. Independent review checks the equations, implementation,
trapped Schwarzschild signs, canonical boosts and constraint qualification.
No physical evolution is run. The next executable connection is exact
source initialization and a characteristic update inside the existing
rectangle, followed by wall-born cells and the regular cap march; it uses
the same source and model, with its own finite run budget.

## Exact boundary-operator implementation, 16 October 2026

The [finite-amplitude owner](nsc-atlas-finite-amplitude.md) supplies the
same-action scalar, Israel and normal-force laws. The next implementation
is their pointwise outgoing-to-incoming trace map, with fixed `b=2`,
`h_minus=1`, `h_plus=0`. Use the invariant normalized field `c`, so all
algebra has a regular `alpha=0` limit and never divides by `beta_plus`.
The helper accepts compatible positive wall lapse traces and their
derivatives; it does not independently choose a second central-clock gauge.
A separate regular-centre readout converts arbitrary positive null-label
rates to Jordan proper time and the physical normalized scalar.

Add these pure helpers and independent tests to the existing package.
Check exact vacuum recovery, the reference zero of `beta_plus`, canonical
scalar and force identities at nonzero amplitude, null-coordinate
rescaling, central-clock conversion and malformed domains. Each test has
a 30-process-CPU-second ceiling. The domain is finite real input with
positive Einstein radius and lapse, nonnegative amplitude parameter and
compatible timelike rates; the coupled caller owns the physical domain
and chart compatibility. No PDE solve, physical amplitude selection,
new preparation or evolution record is included. An inconsistency in
frame, orientation, gauge or matching is resolved in the specific identity
before using this map as a nonlinear boundary solver.

**Implementation outcome:** the pure boundary and central-clock helpers
pass 32 independent controls in 0.28 seconds of pytest-reported duration.
They retain normalized scalar information at zero amplitude, handle the
reference zero of `beta_plus`, and recover canonical-action matching and
observer quantities under coordinate changes. No physical solve or
evolution record is produced. The next implementation connection is
compatible regional gauge and bulk evolution with regular central data;
the pointwise operator does not supply that coupled solution map.

## Central scalar reconstruction admission, 16 October 2026

The existing central-clock cap equation and regular reflection determine
the central cubic scalar coefficient by incoming-ray quadrature once the
exterior wall value and its material derivative are known. Two independent
derivations give the same relation; a second cap PDE solve is unnecessary
for this particular readout. The missing connection is numerical evaluation
of the incoming clock derivative and the reflected cap forcing integral.

Keep the same action, reduced-field preparation, reference geometry, quiet
anchors and finite no-return domain. Use the full regular cap field and
radial metric on incoming rays: the pre-arrival profile-moment simplification
does not apply where the reflected profile is nonzero. Obtain clock
derivatives from the existing first-jet optical identity, not differences
of stored clock samples. Report the inherited continuum/arrival-spline rate
discrepancy. Evaluate the Jordan central field at fixed central proper time;
no second arrival shift is added to that already pulled-back observable.

Implementation and independent analytic/manufactured tests are admitted
while the exterior propagation is resolved. Tests have a 30-process-CPU-second
ceiling each and cover the optical derivative, incoming integration signs,
regular reflection, clock pullback, zero source and cubic scaling. Physical
central evaluation requires the exterior result and a separately specified
shared budget and quadrature comparison. A centre cancellation failure or
domain exit requires diagnosis. This connection does not supply a physical
preparation law, finite-amplitude remainder or detector quantum state.

**First central comparison admitted:** after the implementation and
independent controls are checked, compare incoming radial orders 16/32
at ray order 16 and five retained emission cuts
$t=1/8,3/8,1/2,5/8,7/8$. Their central arrival labels lie strictly inside
the retained window. Use one shared 129-node leading field and 257-knot
geometry with the existing cap derivative/centre settings. Authenticate and
reuse the second solution of the cap `radial_fine` record (radial 32,
ray 16, cap clock 16), the new clock-32 wall/boundary record and the
fine-cache/fine-grid/tighter-tolerance exterior readout. Check all common
coordinates, leading profile and material wall data before reuse. The cap
wall clock $T_1$ and material radius data are independent of the exterior
clock panel refinement. Recompute the incoming travel integral alongside
the forcing at each new radial setting, then reconstruct the fixed-$V$
conversion using the new $V_1$. Compare the fresh travel integral with its
saved radial-32/ray-16 value; the fine case should reproduce it. No new
wall, cap-clock or scalar evolution is needed.

The complete pair, including its single leading-field preparation and all
new incoming integrals, shares 120 process-CPU seconds and 30 million scalar
kernels. Keep all physical inputs fixed. Record the optical time integral,
clock derivative, reflected forcing and their contributions to the central
coefficient. Propagate the already measured exterior cache/grid/tolerance
differences through the same reconstruction algebra without another field
solve. Report quiet recovery and numerical centre behaviour. The static/null
frame cancellation is an algebraic check, not an independently evolved
second solution. Refine ray order or other inputs only if the outcome
requires it; a budget failure or non-refining result gets a diagnosis first.

**First central outcome:** the five-cut radial 16/32 pair takes 27.9375
process-CPU seconds and 14,498,108 scalar kernels. The maximum central
coefficient magnitude is 34,413.8691, and its maximum radial change is
$1.09781\times10^{-5}$. The fine incoming travel integral reproduces the
authenticated earlier calculation within $1.43\times10^{-14}$. Quiet
recovery is below $3.49\times10^{-137}$. This is a large formal coefficient
relative to the sampled leading-field peak 42.2137; it does not fix an
admissible finite amplitude.

The incoming forcing and optical time derivative are new integrals, so
their ray-resolution sensitivity must also be measured. Authenticate and
reuse the radial-32/ray-16 central result, then evaluate radial-32/ray-32
at the same five cuts. One fresh 129-node field with 257 geometry knots
and the unchanged derivative/centre settings supplies the new rays; the
old wall, cap clock, outgoing boundary and exterior response remain reused.
The complete new calculation, including preparation, has 120 process-CPU
seconds and 30 million kernels. No new wall/cap-clock/scalar evolution is
admitted. Compare the incoming forcing, optical derivative, travel
integral, every central contribution and total coefficient. This separate
ray comparison retains the same physical inputs; unresolved drift or a
budget failure requires diagnosis before another resolution increase.

**Incoming-ray outcome:** the new radial-32/ray-32 side takes 50.046875
process-CPU seconds and 26,350,673 scalar kernels. The maximum central
coefficient change is $5.28497\times10^{-5}$, with separate optical-time,
forcing and travel-integral changes $1.33819\times10^{-6}$,
$3.60495\times10^{-7}$ and $5.20715\times10^{-10}$. Quiet recovery
persists. The saved coarse side and its previous cost remain explicit.
The [record](../results/development/nsc-atlas-central-ray-v1/ray.json)
completes the two admitted incoming quadrature comparisons; no further
resolution increase is required to report these finite sensitivities.

The next mathematical connection is [finite-amplitude control](nsc-atlas-finite-amplitude.md).
Its independently reviewed exact action reduction supplies the nonlinear
source, bulk and wall equations and a source-only positivity argument.
Analytic investigation of the coupled interface/centre estimate is useful
before another physical calculation. A source-leg bound or small truncated
coefficient does not establish that full-domain estimate. No nonlinear
physical run or larger numerical budget is admitted by this paragraph.

**Reference-collar outcome:** the exact normal-geodesic and Jacobi equations
give positive lapse and bounded scalar coefficients on compatible patches
of half-width $\ell/20$ about the original wall interval. Rational
inequalities close the reference scalar energy estimate with coefficient
$119/8$. Its outer flux remains explicit, and homogeneous transmission is
required even when the field represents a residual. This is an analytic
reference result, not a new physical evolution or a sourced-metric bound.

The next bounded analytic tasks are to close the first-commuted collar
energy using those constants and to derive a forced regular-centre
readout estimate for the same reference cap. Keep the same action, source
regularity and finite domain. Retain supplied outer traces, source norms
and patch compatibility; no whole-domain inverse follows from omitting
them. Independent derivations and exact arithmetic are sufficient for this
step. Any supporting algebraic check has a 30-process-CPU-second ceiling;
no physical solve or source refinement is included. A derivative loss or
unprovided boundary datum redirects the estimate at that specific term.

**Scalar-estimate and constraint outcome:** independent derivations close
the first-commuted reference collar energy with coefficient $1549/8$,
and give both a characteristic and an energy-trace estimate at the regular
centre. Initial energy, forcing and outer incoming traces remain explicit.
The exact bulk equations propagate $rE_U$ along constant-$U$ rays and
$rE_V$ along constant-$V$ rays. Their wall audit identifies the remaining
incoming constraint mode and the independent normal-force equation that
removes it when full Israel and scalar matching are imposed. This is an
analytic identity, not a solved nonlinear interface evolution. The next
step selects a noncircular independent boundary system and compatible
metric estimates under the same finite analytic scope; no physical run is
admitted by these results.

**Boundary-selection outcome:** independently imposing the incoming
exterior null constraint, temporal junction and normal force makes the
angular mismatch obey a homogeneous first-order equation. Initial angular
matching therefore propagates, and the cap's incoming constraint vanishes.
The derivation retains residuals until this conclusion; it does not use
the constraint-dependent temporal shortcut. This fixes a consistent
constraint prescription for sufficiently regular solutions. Construction
and bounds for the coupled boundary solution map remain open.

## Conservative scalar evolution admission, 16 October 2026

The [evolution connection](nsc-atlas-scalar-evolution.md) reuses the already
derived fixed-coordinate equation and quiet anchors. Its missing numerical
inputs are a resolved wall lift and a checked interpolation of the exterior
forcing. Implementation and manufactured controls may proceed while the
local-lift comparison is evaluated; numerical propagation of the physical
input follows its measured accuracy and the first cache comparison.

The first cache comparison holds the prepared field at 129 nodes, geometry
at 257 knots, radial Gauss order 2, field derivative order 32, source order
64 and 129 radial cache nodes. Compare 17/33 uniformly spaced training times,
including both source edges, with the 32 fine-mesh midpoints withheld from
both fits: 65 unique direct forcing slices. Each withheld slice includes
the 129 radial knots and their 128 midpoints, staying within the 513-point
field interface. Report time errors at the radial knots and combined
time/radial errors at midpoints separately. Both fits share the original
field preparation and use the same source state. Track separate per-block
and per-coordinate discrepancies, exact source/wall/quiet identities and
the analytic source derivative. Preserve the post-input source correction.

The finite budget is 120 process-CPU seconds and 30 million scalar kernels
for the complete first comparison, including its preparation and all
reference slices. At most 129 correction nodes and 6,000 correction RHS
calls are allowed in the subsequent first propagation pilot, sharing a
fresh complete 120-second/30-million-kernel budget with any inputs it
recomputes. Each manufactured test has a 30-process-CPU-second ceiling.
No field or metric integration occurs inside the correction RHS.

Independent controls include the zero input, nonconstant exact source and
lift, a nonlinear moving-wall manufactured solution with nonzero flux and
potential, source subtraction, time normalization, fixed-$U$ clipping and
wall material derivatives. The physical comparison must report input-cache
effects separately from correction-grid and ODE-tolerance effects. Sampled
input differences are not uniform bounds. Budget exhaustion, unresolved
cache structure or an incompatible source/quiet identity requires a
diagnosis before another pilot. A formal wall coefficient does not complete
the cap detector waveform or finite-amplitude control.

**First cache outcome:** 65 direct slices cost 22.84375 process-CPU seconds
and 4,431,777 scalar kernels. Time-knot errors refine from 0.06551 to
0.007002 in $F$ and 7.91622 to 0.090853 in $G$; source-pullback error
refines from $1.206\times10^{-4}$ to $1.415\times10^{-5}$. Radial-midpoint
errors are recorded separately as combined time/radial sensitivities.
This measures the first cache's cost and accuracy. The next scalar pilot
must compare the two caches on one correction grid, then change only the
correction grid at the finer cache. It must retain these sampled input
effects and the unresolved uniform-error distinction.

**First restored-data propagation admitted:** after the clock-refined lift
meets its existing sampled boundary target, authenticate the saved forcing
slices and clock boundary jets. Reconstruct the original 257-knot reference
geometry and both interpolants with a shared fresh budget; perform no new
field, cap or wall solve. Run three correction evolutions: 17/33 cache times
at 129 correction nodes, and the 33-time cache at 65 nodes. Hold the
clock-32 lift, relative ODE tolerance $2\times10^{-10}$, absolute tolerance
$10^{-12}$ and maximum normalized step $1/64$ fixed. All three solves and
their readouts share 120 process-CPU seconds and 30 million scalar kernels;
each solve permits at most 6,000 RHS calls. Check retained reference
coordinates against the reconstructed geometry before propagation.

Report wall values and material derivatives at 17 cuts, source/quiet/wall
identities, actual cache and grid differences, and their size relative to
the computed correction. The existing conditional estimate gives sampled
input-sensitivity indicators when supplied with withheld differences, not
uniform bounds; do not compare $F,G$ directly with the differently scaled
boundary target. This first pilot computes an exterior formal coefficient.
The cap detector signal, ODE-tolerance comparison where consequential, and
finite-amplitude remainder retain their separate next questions.

## Profile moment and local lift admission, 16 October 2026

The [profile-moment successor](nsc-atlas-profile-moment.md) checks the exact
pre-arrival mass identity, then supplies more boundary samples where the
saved interpolation comparison failed.

| Admission item | Fixed decision |
| --- | --- |
| Retained mathematics | Quiet regular cap, the full radial mass constraint and flux, the derived first-derivative profile moment, outgoing scalar forcing, material wall pullback and Hermite interpolation. |
| Missing connection | A cheaper numerical cap-mass evaluation agreeing with the original radial integral; then a resolved local boundary function and exactly integrated lift for scalar propagation. |
| Fixed model and state | Same conformal scalar, h=1, b=2, m=0.4, original prescribed reduced-field input, quiet anchors and finite wedge. No finite alpha, new source or fitted physical parameter. |
| Moment comparison | Field 129, geometry 257, unchanged wall/cap settings, moment Gauss 4/8 and outgoing ray 16. Share one field/wall preparation across both settings. Compare mass at wall and two interior ray fractions at up to nine cuts (at most 27 points), and J_minus, its material derivative, b and b_U against nine original radial BoundaryJet queries. Split moment cells at existing profile edges, field time cuts and arrival spline knots. |
| Local lift comparison | Proceed after the moment discrepancy is explained and small enough for the intended input accuracy. Hold the selected moment/ray settings fixed. Compare active 16/32 cells with four cells in each outer source block, using five-node local Hermite stencils. Fine 41 training knots plus 40 withheld midpoints give 81 unique queries. Track value/slope/lift differences and floating continuity defects; do not reuse fitted points as holdouts. |
| Budget | At most 120 process-CPU seconds and 30 million scalar kernels per complete pilot or two-setting comparison, including preparation, wall evolution, moment prefixes and independent reference queries. Tests of analytic/manufactured identities use at most 30 process-CPU seconds each. Temporary arrays remain in .build; consequential records retain settings, source bytes and costs. |
| Independent controls | Exact polynomial moments, radial mass and flux identities, joint moving-ray derivatives, original radial metric at interior and wall points, zero source, finite-domain/budget rejection, polynomial reproduction, C1 assembly and exact lift differentiation across changing stencils. |
| Stop or redirect | Diagnose non-refining moment/radial differences before adopting the faster evaluator. A dense-field/continuum-jet discrepancy is not removed by increasing moment order. If local fit errors miss the saved field-based target, locate the failing cells and assess new samples or field resolution; do not raise polynomial degree indefinitely on unchanged data. Neither interpolation agreement nor metric coefficients establish a finite-amplitude solution. |

The likely failures are a numerical profile/derivative mismatch in the exact
mass identity, loss of a source edge or quiet integration constant, and a
local polynomial whose values or slopes fail near a stencil transition.
Independent radial queries, split cumulative integration and genuinely
withheld cell checks expose these failures before scalar evolution.

**Moment outcome:** the complete 4/8 pair and original radial reference take
12.25 process-CPU seconds and 2,637,474 scalar kernels. The maximum mass
quadrature change is $1.07\times10^{-14}$; the independent radial difference
is $2.50\times10^{-9}$. For the boundary value those differences are
$2.67\times10^{-15}$ and $1.52\times10^{-10}$, respectively. The separate
radial/dense-field representation discrepancy is not removed by moment
refinement, but its sampled boundary effect is over 1,500 times below the
existing interpolation target. The full wall-$U$ width times the slope
discrepancy is below $1.8\times10^{-9}$. This admits the already specified
local lift pair at moment order 8, holding all physical inputs fixed.

**Local-lift outcome:** the 81-query pair takes 17.9375 process-CPU seconds
and 3,737,156 scalar kernels. Withheld value errors are
$3.95548\times10^{-5}$ and $2.90055\times10^{-5}$, both above the original
target. The maximum lift difference is $2.47058\times10^{-8}$; assembly
and transport roundoff are below $5.33\times10^{-15}$. The fine fit now
fails most strongly through the pulse interior rather than its edges.
Before another physical lift calculation, use the saved 81 first jets to
decompose the same boundary into backscatter, material-radius and clock
pullback terms, reconstructing the same fits and withheld comparisons.
Budget: 30 process-CPU seconds, zero physical solves or new scalar kernels.
Check value and derivative sum closure. In particular, inspect the inherited
eight-panel clock's value/continuum-slope consistency; the earlier nine-cut
comparisons sampled only its panel boundaries. Do not silently accept the
unresolved interpolant or increase its degree again.

**Saved-component diagnosis:** zero additional physical solves isolate the
fine-fit error in the clock pullback. Maximum withheld value discrepancies
are $5.99\times10^{-10}$ for backscatter, $1.67\times10^{-8}$ for material
radius and $2.901\times10^{-5}$ for the clock. Direct value and derivative
sums close to floating roundoff. The source inspection confirms that clock
values use panel-polynomial integrals while their supplied slopes use
continuum rates at the query. This identifies the term requiring refinement,
without by itself bounding its error against the continuum clock.

**Clock refinement admitted:** compare wall clock panels 16/32 with the
same field 129, geometry 257, original wall settings and fixed fine lift
(41 training knots, 40 withheld queries). Reuse the authenticated outgoing
$J_-,(J_-)_s$ from the local-lift record: these depend only on the unchanged
leading field and geometry. Recompute each wall and its boundary pullback,
checking coordinates and leading profile consistency. Share one field
preparation and one 120-process-CPU-second/30-million-kernel budget across
both wall solutions and all 81 queries per solution. Compare the retained
eight-panel result, withheld value/slope errors, lift differences and the
clock's polynomial/continuum slope discrepancy. Success must resolve the
original target without changing physical data. A plateau or budget failure
requires a new diagnosis; no further degree increase is assumed.

**Clock outcome:** the complete 16/32 pair takes 35.234375 process-CPU
seconds and 5,002,056 scalar kernels. Withheld boundary-value errors fall
from $1.73136\times10^{-6}$ to $5.74819\times10^{-8}$; the finer result
passes the unchanged target. The maximum lift difference is
$3.24913\times10^{-10}$. Continuum versus panel-polynomial clock slope
differences decrease independently, and their boundary identity closes
below $1.92\times10^{-12}$. The 81 outgoing cap jets are reused with no
new cap queries. This admits the saved-input propagation specified above.

**First propagation attempt:** the 17-time-cache/129-node case exhausts
its 6,000-RHS allowance at relative tolerance $2\times10^{-10}$ before
producing a result. No completed record is retained. Inspect the saved
forcing's radial interpolation crossings and solver split inventory
without a new field or scalar evolution. Diagnose whether nonsmooth
coefficient derivatives make that tolerance disproportionate to measured
input errors; any revised calculation needs its stated controls and a
separate tolerance comparison. The physical model and saved inputs stay fixed.

**Cost diagnosis and revised propagation:** the 129-node grid has 107
existing evolution segments, but crosses 2,609 distinct radial cache knots;
2,189 crossings after the quiet quarter have nonzero coefficient-slope
jumps. Splitting at every such knot would require at least 41,760 DOP853
RHS evaluations before any adaptive refinement, versus a 1,712-call minimum
for the current splits. The saved-input diagnosis uses no propagation.
Retain the cache and existing splits. Revise the first three cache/grid
comparisons to relative tolerance $2\times10^{-6}$ and absolute tolerance
$10^{-8}$, and add one fine-cache/129-node evolution at $2\times10^{-7}$
and $10^{-9}$ solely to measure ODE sensitivity. Keep maximum step $1/64$
and the 6,000-call ceiling per solve. All four solves and readouts still
share 120 process-CPU seconds and 30 million kernels. Report the actual
tolerance effect separately; its smallness is a measured outcome, not an
assumption from comparing the requested tolerance with input indicators.

**Propagation outcome:** all four cases complete in 5.9375 process-CPU
seconds with 950,926 scalar kernels and no new physical input solves.
The sampled wall peak is $0.2679995064$. Maximum cache/grid/ODE value
changes are $1.75501\times10^{-3}$, $6.98421\times10^{-6}$ and
$3.72805\times10^{-9}$; the derivative changes are recorded independently.
Source-value, source-rate and wall-flux defects vanish. This supplies the
formal exterior coefficient and material derivative for central
reconstruction. The 17-cut sensitivities are not uniform error bounds.

## Scalar lift and forcing admission, 16 October 2026

The sampled outgoing boundary is now known. The
[next input calculation](nsc-atlas-scalar-forcing.md) constructs its lift
and the cumulative exterior forcing needed by the scalar evolution.

| Admission item | Fixed decision |
| --- | --- |
| Retained mathematics | The derived forced scalar equation, quiet outgoing cap flux, first-jet mass constraints, displaced-wall pullback, Hermite interpolation and the existing clipped Volterra operator. |
| Missing connection | A differentiable wall input and its exactly integrated lift; exterior mass/lapse and scalar forcing on a full radial slice with one cumulative integration. |
| Fixed model and state | Original conformal scalar, h=1, b=2, m=0.4, prescribed reduced-field input, original source gauge and quiet anchors. No new source, fitted coupling or finite alpha. |
| Finite domain | The original retained exterior wedge and wall x in [0.45,0.55], with cap queries confined to the already certified pre-return outgoing rays. |
| Boundary implementation | First-jet analytic J_minus and wall derivatives; Hermite interpolation in U_W with an exact antiderivative. Initial 9/17 knots, source edges retained, 16 additional withheld midpoints, 33 unique locations total. Field/geometry 129/257, cap radial/ray 16. Record b, width-scaled slope and lift discrepancies; refine only failing panels after diagnosis. |
| Forcing implementation | One source trace and cumulative radial constraint per V across original cells and requested coordinates. Nine wall-time slices, radial orders 2/4, field derivative order 32. Check source compatibility, quiet data, wall flux and the independent timelike Jordan mass. |
| Budget | At most 120 process-CPU seconds and 30 million scalar kernels per complete pilot or two-setting comparison, including preparation and wall solve. One initial pilot per connection; refinements respond to the measured discrepancy. Temporary arrays stay in .build; source-bound compact records retain consequential inputs, outputs and costs. |
| Independent checks | Symbolic first-jet identities, nonconstant manufactured ray derivatives, exact polynomial lift differentiation, cumulative polynomial integrals, source-pullback derivative, zero source, clipping and independent wall mass transport. Withheld samples and separate refinements test numerical sensitivity. |
| Stop or redirect | Diagnose domain exit, slope/value inconsistency, failed refinement or budget exhaustion before another calculation. Select interpolation tolerance below measured field-input sensitivity, with an explicit roundoff floor; a sampled tolerance is not a uniform theorem. Scalar evolution and its time cache require the measured cost and accuracy of these inputs first. |

The likely failures are inconsistent boundary slopes, quadrature across
unresolved structure, and repeated field queries exhausting the shared
budget. Withheld interpolation checks, split cumulative quadrature and
shared counters identify those failures without changing the physical model.

**First interpolation outcome:** the 9/17-knot pair completes in 17.05
process-CPU seconds and 4.89 million scalar kernels. Its maximum withheld
boundary error decreases from 0.555 to 0.0611, above the diagnostic target
$2.28\times10^{-7}$ derived from the saved field sensitivity. The largest
first-jet identity defect is $2.88\times10^{-13}$; source-edge splitting
and higher interpolation order therefore merit investigation before another
physical evaluation. A read-only postprocessing pilot may use four, five,
seven or nine local Hermite training nodes from the existing 17-knot mesh,
with stencils confined to the three source-edge blocks. Keep the same 16
withheld queries out of every fit. Budget: 30 process-CPU seconds, zero field
or wall solves. Report value and width-scaled slope errors per block;
non-improvement or a noise plateau redirects to new independent samples and
the analytically derived pre-arrival mass moment, rather than accepting the
unresolved cubic fit.

**Postprocessing outcome:** the admitted four stencils require 0.078125
process-CPU seconds, zero new field/wall solves and zero scalar kernels.
Their maximum value errors are 0.0123173, 0.00956095, 0.0849837 and
0.765149 for four, five, seven and nine requested nodes. None meets the
original target. The post-pulse block reaches $3.01\times10^{-10}$ with
its five available nodes; the unresolved active-pulse errors peak near
its source edges. The
[assessment](../results/development/nsc-atlas-scalar-lift-assessment-v1/assessment.json)
therefore redirects to additional local samples and a numerical comparison
of the derived mass-moment formula. Further degree increases on the same
training data are not supported by this outcome.

**First forcing outcome:** the cumulative 2/4 radial pair completes in
9.52 process-CPU seconds and 1.20 million scalar kernels. The maximum
changes in $F,G$ are $7.20\times10^{-14}$ and $4.45\times10^{-10}$;
the separate timelike mass discrepancy is $6.50\times10^{-7}$ at both
orders. Further radial-order increases would not address that larger
field/transport difference. The
[saved inputs](current-result.md#scalar-lift-and-exterior-forcing-16-october-2026)
retain the source pullback after the leading pulse, quiet anchors and costs.
The known forcing can next be cached in time with a measured interpolation
error; it has not yet been used to evolve the correction.

## Scalar backscatter boundary admission, 16 October 2026

The [scalar-metric equations](nsc-atlas-scalar-metric-response.md) determine
the outgoing correction at the wall before the exterior correction is
evolved. This calculation evaluates that missing boundary input from the
same prepared field and computed cap metric.

| Admission item | Fixed decision |
| --- | --- |
| Retained mathematics | Forced cap equation in central-clock null coordinates; the no-return certificate; regular radial Einstein constraint; weighted null transmission and the displaced-wall pullback. |
| Missing connection | The numerical outgoing trace J_minus and its induced exterior areal derivative c1_r, with backscatter and displacement contributions kept separate. |
| Fixed model and state | Original conformal scalar, h=1, b=2, m=0.4, prescribed reduced-field input and quiet anchors; no new source, fitted parameter or finite alpha. |
| Finite domain | Nine reference wall cuts, x in [0.45,0.55]. For each cut integrate q from q_a to q_W on the outgoing ray u=q_W-2 atanh(x); T=(u+q)/2 and r=tanh((q-u)/2). These points remain before central arrival and within the existing cap. |
| Implementation | A boundary adapter reuses the same prepared field, cap radial metric and wall displacement/exterior clock. Split quadrature at the existing profile edges. Compute the same integral in q and r, with the correct Jacobian, and compare the general and quiet-profile potential formulas. No new cap clock solve is needed. |
| Budget | At most 120 process-CPU seconds and 30 million scalar kernels per complete pilot or two-setting comparison, including preparation and wall solve. Field nodes 129/257; cap radial and outgoing-ray orders 8/16/32. At most one initial pilot; subsequent field, radial and ray pairs vary one setting. A further 16/32 pair is justified only by an unresolved discrepancy. Transient pilots stay under .build. |
| Independent checks | Manufactured forced-null integration, zero source, cubic scaling, q/r reparameterization, potential identity, wall flux reconstruction and regularity across beta_plus=0. Numerical sensitivity is recorded separately from a rigorous error enclosure. |
| Stop or redirect | Diagnose a domain exit, disagreement of the independent formulas, non-refining residual or exhausted budget before another run. Nine boundary samples do not define a controlled interpolant or a propagated scalar waveform; those are subsequent connections. |

The likely numerical difficulties are missing a pulse edge, losing a sign
in the outgoing flux, and allowing repeated metric queries to exceed the
shared budget. Split integration, two coordinate evaluations, signed
contributions and cumulative counters expose those failures directly.

The first field pair completed in 34.61 process-CPU seconds and 8.45 million
scalar kernels. The field 129/257 change in the exterior wall derivative
coefficient is $1.73\times10^{-6}$. Radial order 8/16 changes it by
$3.99\times10^{-7}$, while outgoing-ray order 8/16 changes it by
$3.26\times10^{-3}$. The latter is much larger than the field sensitivity;
the radial contribution remains comparable in the outgoing trace itself.
These observed discrepancies justify the admitted separate 16/32 radial
and ray comparisons at 129 field nodes, without changing the shared budget.

**Outcome:** all five pairs complete in 19.06--34.61 process-CPU seconds,
with at most 8.45 million scalar kernels. The separate 16/32 radial and ray
changes in the wall derivative fall below $1.17\times10^{-10}$ and
$1.31\times10^{-10}$. The [result record](current-result.md#scalar-backscatter-boundary-16-october-2026)
retains the nonzero outgoing coefficient, signed contributions, coordinate
comparison and remaining field sensitivity. This closes the sampled
boundary-input calculation. A controlled interpolant and the forced scalar
evolution require their own discretization and admission; no further
boundary refinement is justified by the current discrepancies.

## Scalar response to the computed metric, 15 October 2026

The [cap and clock calculation](nsc-atlas-cap-response.md) now completes the
formal quadratic metric connection. The next field coefficient requires
the coupled cap/exterior scalar equation on that geometry. The
[scalar-metric owner](nsc-atlas-scalar-metric-response.md) derives that
problem before a numerical evolution is designed.

| Admission item | Fixed decision |
| --- | --- |
| Retained mathematics | Conformal covariance of the scalar equation; spherical reduced-field equations; the computed Einstein mass/lapse coefficients, material displacement and clock embeddings; the original scalar and metric junction laws. |
| Missing connection | Order-epsilon-alpha bulk forcing, compatible transmission at the displaced wall, the source-coordinate pullback and the waveform at fixed central proper time. |
| Fixed continuation | Same action, couplings, quiet data and finite reference geometry. For this candidate continuation keep the prescribed original reduced-field input exactly chi(U0,V)=epsilon g(V). This fixes its higher-order input; it does not select a physical preparation mechanism. |
| Domain | Original exterior wedge and regular cap; first fixed-clock waveform comparisons lie on a compact subinterval inside the certified central arrival interval. A moved finite-amplitude terminal arrival needs its own domain analysis. |
| Current work and budget | Derive the linear inhomogeneous transmission problem and check consequential identities symbolically or on manufactured fields, at most 30 process-CPU seconds per check. No scalar evolution or new physical parameter scan is performed in this increment. Numerical discretization, boundary treatment and its finite compute budget will be recorded here before a coupled solve. |
| Independent checks | Conformal reduced-field invariance; expansion of the exact spherical wave operators; radial-constraint cancellations; wall-normal normalization and transmission; regular-centre cancellation; fixed-clock pullback without counting the null shift twice. |
| Stop or redirect | A sign, frame or source-gauge mismatch is resolved at its equation. The leading no-return wall condition cannot be reused if the changed cap metric produces backscatter. These equations do not yet supply a scalar waveform or finite-amplitude remainder. |

The likely failures are confusing the conformally transformed scalar with
the canonical Einstein scalar, dropping the source or wall pullback, and
imposing the leading outgoing condition on the correction. The separate
frame definitions, expanded material normals and coupled transmission
equations make those errors directly checkable.

**Analytic outcome:** fourteen independent symbolic and manufactured-field
checks pass in 3.3125 process-CPU seconds total. They verify the reduced
operators, curvature and radial-constraint cancellations, source pullback,
central null coordinates, waveform reconstruction and weighted transmission.
A regular static manufactured field has nonzero cap scattering potential
$K_1=4r^2/135+4r^4/189+O(r^6)$. Within the original no-return interval,
the outgoing correction at the wall is a known integral of this type of
bulk forcing, independent of the unknown incoming correction. This permits
a sequential transmission calculation at this order. No physical scalar
correction has been evolved.

## Smooth cap and central-clock admission, 15 October 2026

The wall now has computed mass, displacement and exterior-clock coefficients.
The next connection reconstructs the cap metric from the same transmitted
scalar and calculates central proper-time arrivals. The
[cap owner](nsc-atlas-cap-response.md) states the equations, regular centre,
clock convention and finite domain.

| Admission item | Fixed decision |
| --- | --- |
| Retained mathematics | Spherical Einstein radial constraints and mass flux; the existing conformal-frame map, transformed fixed cap vacuum, regular reflected scalar solution and prepared wall trace; the checked two-sided mass/displacement laws and exterior clock. |
| Missing connection | The smooth cap mass/lapse and its central proper-clock normalization, an independent interior check of the wall mass, and the corrected null-arrival map at fixed proper wall time and at fixed exterior null label. |
| Fixed model and state | Original conformal scalar, compact reduced-field preparation and quiet cap, h=1, b=2, m=0.4; same action and source gauge, formal alpha=kappa4 epsilon^2/ell^2, no fitted source or finite amplitude. |
| Finite domain | The reference cap under the original wall interval x in [0.45,0.55], closed after the last wall event by its inward null ray to the centre. Profile labels end at the existing last arrival q_b. Reuse the no-return certificate for quiet earlier labels; retain the reflected term and transmitted tail. |
| Implementation | New cap field and response modules reuse the existing dense field and wall solution without changing their source-bound records. Profile derivatives come from the continuum integral equation. Integrate regular-centre radial constraints, wall time and inward null rays; keep the two clock comparisons distinct. |
| Budget | At most 120 process-CPU seconds and 30 million scalar-kernel values per complete pilot or two-setting comparison, including field/wall solves and queries. Field nodes 129/257, geometry knots 129/257, at most 64 Gauss nodes per radial/ray panel, 32 wall-time panels and 129 displayed cuts. Initial pilots determine practical quadrature settings before compact evidence is frozen; transient arrays remain under .build. |
| Independent checks | Analytic manufactured fields; scalar equation and mass conservation; wall derivatives; radial cap mass against the existing timelike integral including its potential and frame term; zero source; Jordan centre clock; independent wall proper-normalization expansion; fixed-s to fixed-V conversion; separate numerical refinements. |
| Stop or redirect | Diagnose centre cancellation, missing/out-of-domain profile labels, a mass-path defect that fails to refine, incompatible clock normalization or budget exhaustion before another run. A metric ray correction does not establish the higher-order scalar amplitude, finite-amplitude remainder, formation or stability. |

The three likely difficulties are cancellation near the regular centre,
quadrature that misses moving pulse structure, and conflating the wall's
proper time with the external label or central clock. Continuum derivatives,
split quadratures, independent mass transport and explicit clock conversion
make those failures detectable while keeping the physical model fixed.

The first time-panel pair exhausted the shared budget while evaluating the
same arrival rays twice. Changing only the wall-time panels changes the
$T_1$ integral, while the leading cap field and reference rays stay fixed.
The revised producer reuses those rays after checking their identical inputs,
records the saved work separately, and keeps the original budget. The
observed radial-order 8/16 arrival difference of $6.30\times10^{-3}$ also
motivates a 16/32 radial comparison with 129 field nodes held fixed. This
tests the interior ray metric independently of the wall mass-path residual;
it uses the same already admitted orders and compute ceilings.

**Outcome, 15 October 2026:** ten separately varied comparisons now connect
the smooth cap mass/lapse, moving-wall cap time and central null arrival.
The primary mass-path discrepancy is $2.52\times10^{-9}$. The finer
16/32 radial and ray controls change the arrival coefficient by less than
$2\times10^{-9}$ each at fixed 129-node field resolution; field 129/257
changes it by $6.42\times10^{-7}$, while cap time panels 8/16 change it
by $1.92\times10^{-5}$. Each complete pair used 35.98--56.31 process-CPU
seconds and fewer than 24.94 million scalar kernels. The
[current result](current-result.md#smooth-cap-and-central-arrival-clock-15-october-2026)
retains the settings and all ten evidence bundles. The next connection is
the higher-order scalar response on this changing geometry, followed by
finite-amplitude control; the null-clock coefficients alone do not supply it.

## Exterior wedge and moving-wall calculation admission, 15 October 2026

The rectangle has supplied a checked gravitational input. The next connection
uses the same leading scalar on the full retained wedge to evaluate the radial
mass and lapse at the wall, compare the independent timelike mass functional,
and calculate the displacement and exterior clock correction. The cap's smooth
metric and corrected central observation follow afterward.

| Admission item | Fixed decision |
| --- | --- |
| Retained mathematics | The [bulk radial mass/lapse laws](nsc-atlas-bulk-response.md), [two-sided mass and displacement functionals](nsc-atlas-cap-wall-response.md), original prepared-field Volterra equation and signed wall normalization. No new action, source, coupling or amplitude. |
| Missing connection | Controlled evaluation beyond U_a up to U_W(V), including the moving lower endpoint of the transverse derivative; the resulting metric input to the material wall's proper-time evolution and exterior null clock. |
| Fixed data | Original compact C3 reduced-field preparation, quiet cap/bottom/corner, h=1, b=2, m=0.4, same Jordan source gauge and canonical conformal scalar. All response quantities are formal coefficients of alpha=kappa4 epsilon^2/ell^2. |
| Domain | Original no-return wedge U_0<=U<=U_W(V), V_a<=V<=V_b, wall radius in [0.45,0.55]. Evaluate the exterior kernel only on this retained domain. For U>U_a, start chi_U at the inverse wall label L(U), with chi_U=0 there. |
| Implementation | Expose the existing dense field solution additively; paired field evaluation and split quadrature supply first derivatives. Integrate radial mass/lapse and the independent proper-time masses, then the quiet-anchored displacement and exterior clock. Reuse source data and field values; never differentiate stored waveform samples. |
| Compute budget | At most 120 process-CPU seconds and 30 million scalar-kernel evaluations for each complete pilot or two-setting comparison, including its field solves and queries. Start with 129/257 field nodes, temporal Gauss orders 32/64 and per-cell radial orders 2/4. At most 129 displayed time cuts; dense grids and exploratory output stay under .build. Diagnose exhaustion before redesigning a run. |
| Independent checks | Rectangle overlap; exact quiet/source/wall limits; nonzero physical p_U despite chi_U=0 at the wall; inverse-wall bracketing; radial versus timelike Jordan mass; both displacement equations; angular and clock normalization; regular beta_+=0 crossing; field, radial, temporal, source and geometry refinements varied separately; zero source and coordinate/unit controls. |
| Useful output | A finite numerical wall response and exterior clock coefficient with measured refinement differences and mass-path defects. Preserve the six original rectangle records under explicit commit-pinned replay when their source API evolves. |
| Stop or redirect | Diagnose an unprovided datum, field query outside its domain, mass-path disagreement that does not refine, singular normalization or exhausted budget. A formal coefficient alone supplies neither a finite-amplitude remainder nor the unsolved central-observer correction. |

Three likely failures are crossing the moving domain while evaluating a kernel,
losing a frame or cap-vacuum term in the mass comparison, and confusing the
exterior null coordinate with the central observer's clock. Domain checks,
independent mass paths and separately named observer maps make each failure
detectable without changing the physical inputs.

**Outcome, 15 October 2026:** seven separately varied comparisons now evaluate
the wall masses, displacement and exterior clock. The 129/257 field comparison
reduces the maximum mass-path discrepancy from 6.49e-7 to 1.63e-7; the clock's
8/16-panel comparison changes its ratio coefficient by at most 2.41e-8.
Every comparison completed within 120 process-CPU seconds and 30 million
scalar-kernel values. A separate radial-versus-mixed rectangle overlap checks
both mass and lapse with a distinct quadrature. The
[current result](current-result.md#moving-wall-and-exterior-clock-15-october-2026)
records all eight evidence bundles. The smooth cap and central clock are the
next physical connection; these coefficients do not establish a finite-amplitude
remainder.

## Smooth exterior metric propagation admission, 15 October 2026

The next step carries the same prepared scalar and its already fixed
characteristic gravitational data into the smooth exterior. The full target
remains the sourced cap/wall/exterior configuration and corrected observer
response. First solve the characteristic rectangle
$U_0\leq U\leq U_a$, $V_a\leq V\leq V_b$: it provides the transverse
initial leg needed to reach the emerging wall strip, rather than assuming
metric data in the excised cap. The next wall strip must use the derived
material matching and coordinate embedding.

| Admission item | Fixed decision |
| --- | --- |
| Borrowed and retained | The spherical double-null Einstein equations and smooth conformal-frame map already derived in the [cap/wall owner](nsc-atlas-cap-wall-response.md), the [source-leg constraint](nsc-atlas-characteristic-initialization.md), the [transverse mass data](nsc-atlas-source-mass.md), and the existing prepared-field Volterra solver. |
| Missing connection | Propagate the radius and lapse coefficients, compare their geometric mass with independent scalar flux, and supply the resulting characteristic data for matching to the moving wall. |
| Fixed model and state | Same conformal scalar, canonical normalization, h=1, b=2, m=0.4, original compact source and quiet bottom/cap. Use the original null coordinates, fixed Jordan source lapse and quiet corner; transform source data to the Einstein frame without changing the affine coordinate silently. No physical finite amplitude is selected. |
| Finite domain | Original characteristic rectangle through V_b, inside the retained exterior wedge; any later extension U>U_a requires explicit compatible wall data. |
| Implementation | Expose bounded leading field values and first derivatives from the existing dense solution and its integral equations. Solve the linear radius/lapse equations with the supplied characteristic data; use independent mass transport and null constraints as diagnostics. |
| Compute budget | At most 120 process-CPU seconds per complete pilot or refinement comparison, 30 million scalar-kernel values, and fixed metric grids no larger than 257 U nodes. Begin at 33/65 nodes; admit a finer comparison only if its discrepancy identifies a resolvable numerical error. No parameter survey, new source fit, or arbitrary waveform. |
| Checks and useful outputs | Zero source; original source/bottom values; both null constraints and mixed equations; geometric mass versus independent scalar flux; field, geometry, quadrature and metric-grid refinement separately; unit/orientation conventions; cumulative budgets. Save compact coefficients and diagnostics, not dense trial grids or failed-run debris. |
| Stopping or redirection | Diagnose a missing boundary datum, incompatible frame/gauge data, constraint drift that does not refine, or budget exhaustion before designing another run. A solved rectangle supplies initial data for the wall strip; it does not replace the full finite-domain objective. |

**Outcome, 15 October 2026:** the [bulk result](nsc-atlas-bulk-response.md)
propagates the formal radius/lapse response across the initial rectangle
and supplies the next transverse characteristic leg. Six separately varied
comparisons distinguish metric, scalar-field, geometry, derivative
quadrature, source quadrature and time-tolerance effects. The 33/65
metric comparison reduces the largest sampled mass-flux defect from
$1.65\times10^{-6}$ to $4.31\times10^{-7}$, with both integrated null
constraints refining. Each comparison used 17.16--21.91 process-CPU
seconds and fewer than 10.45 million scalar-kernel evaluations. The
[current-result entry](current-result.md#smooth-exterior-metric-response-15-october-2026)
records the coefficients and all six evidence links. The next domain is
the moving-wall strip, followed by the cap and observer correction; the
formal rectangle result does not provide a finite-amplitude remainder.

## Cap mass and moving-wall response admission, 15 October 2026

The accepted programme first closes the cap-side input and the coupled wall
equations of the existing conformal reference. The exterior mass functional
is already available. A new numerical waveform or a cosmological fit would
not resolve this missing connection.

| Admission item | Fixed decision |
| --- | --- |
| Borrowed | The existing full conformal Hilbert stress and scalar/material junction laws of Aviles--Maeda--Martinez, 1910.07534v2; Hayward's general spherical mass transport, gr-qc/9710089v2; the complete mass frame transformation of Faraoni--Vitagliano, 1401.1189v1; the already checked exterior functional and normal-force identity. Exact equations and reading boundaries remain in the existing source cards. |
| Missing | Does the existing transmitted field fix the cap's vacuum-subtracted mass and a compatible quadratic areal-radius response of the material wall? Which bulk or initial data remain necessary afterward? |
| Fixed model | Canonical real massless scalar, xi=1/6, no direct scalar-wall interaction, scalar-independent vacuum stress, constant material tension; h=H ell=1, K=k ell=2, m=M/ell=0.4; original C3 preparation, quiet cap/bottom, same corner and common cap-to-exterior normal. |
| Domain and order | Original no-return wall interval 0.45 <= R/ell <= 0.55 in the isolated reference. Formal order alpha=kappa_4 epsilon^2/ell^2, with perturbed proper wall time anchored at the quiet starting event. Keep cap vacuum mass h^2 R^3/2 and its displacement correction. No formation, late-time, finite-amplitude existence or glued-parent future claim. |
| Decision and unlock | Derive cap transport directly with the Jordan improved stress and independently using the smooth cap frame and transformed vacuum potential. Combine both one-sided masses with the angular junction; test temporal and normal compatibility, identify integration constants and any missing metric data. A closed response permits a separately controlled finite evaluation; an obstruction names the missing datum or incompatibility. |
| Budget | One analytic derivation and independent review; exact symbolic verification within 60 process-CPU seconds per complete certificate. No waveform re-solve, mass quadrature, parameter sweep or bulk PDE solve in this admission. |
| Checks | Vacuum subtraction on the actual displaced radius; both null/scalar derivatives; the cap curvature contribution; transformed vacuum potential; quiet anchors; orientation and length scaling; zero source; regularity at beta_plus=0; differentiation of angular matching with tangential flux; normal Gauss identity; omission controls. |
| Stop or redirect | Stop on an unprovided integration datum, inconsistent material/scalar jumps, singular response coefficient in the admitted regular domain, mismatched independent transport, or exhausted budget. A formal wall coefficient does not establish a full bulk solution or a finite-amplitude error bound. |

The likely obstacles are losing the cap vacuum/displacement terms, mistaking
an angular equation for all junction conditions, and promoting a formal wall
response to a complete spacetime. Independent frame/stress derivations,
explicit temporal/normal checks and the remaining-data statement address them.

**Outcome, 15 October 2026:** the
[two-sided derivation](nsc-atlas-cap-wall-response.md) closes the cap
vacuum-subtracted mass and quiet-anchored displacement functionals. The
saved certificate has 26 exact zero identities and six omission controls;
creation used 1.546875 process-CPU seconds with no numerical integration.
Independent original-frame and frame-transform calculations agree, including
the temporal and normal conditions. The next missing connection is the
smooth bulk radius/lapse solution and wall embedding, followed by observer
clock and field corrections. No additional independent cap mass constant or
wall-pressure history was needed for this formal boundary reduction.

This file prevents Nested-Space Cosmology from spending compute on an already known component and then presenting the reproduction as project novelty. Import a result with its assumptions and attribution. The research contribution must close a previously missing connection, distinguish candidate exteriors, or establish why a proposed connection is unavailable in a named domain.

The research direction adopted on 10 October 2026 is exterior first: build a
finite mathematical map of possible parent regions, children and descendants;
then investigate its consequences for an observer inside a region. The
[cosmological atlas](nsc-cosmological-atlas.md) compares geometric,
exterior-response and connected-mouth candidates. Each test fixes its own
theory, state and domain. No candidate's action is compulsory for the others.

The labeled theory is [THEORY.md](../THEORY.md). The live checkpoint is
[docs/current-result.md](current-result.md). Primary literature is cited in
the [technical paper](../paper/nested-space-cosmology.md).

## Imported results — do not pursue as research targets

| Subject | Imported result | Project use |
|---|---|---|
| Electron wave/particle behavior | Individual electron detections build an interference pattern. | Experimental constraint on any particle ontology. |
| Quantum state | QFT describes particles as field excitations and quantum mechanics assigns amplitudes to their states. PBR constrains purely epistemic wavefunctions under preparation independence. | Distinguish localized energy from the amplitude of the same field state. |
| Bell nonlocality | Loophole-free Bell tests exclude local-hidden-variable explanations. | A gradient medium cannot be an ordinary local three-dimensional fluid; its many-body state must be nonlocal/contextual and no-signalling. |
| Stable field particles | Solitons and Skyrmions are localized finite-energy field configurations; collective quantization can give fermionic baryons. | Construction pattern for stationary charge sectors. |
| Nuclear energy | Nuclear masses and fusion/fission \(Q\)-values are differences of complete bound-state energies; BPS/near-BPS Skyrme models calculate nuclear binding. | Required low-energy projection and calibration/holdout source. |
| Self-gravitating nuclear matter | Full Einstein–BPS calculations produce compact-star mass-radius branches from nuclear-fixed parameters. | Gravity benchmark, not a discovery target. |
| Trapped collapse to expansion | Exact regular black-universe solutions contain a trapped interior, positive complete-null-Raychaudhuri \(Q\), a finite radius minimum, and an expanding cosmological interior. | Limiting solution and transition benchmark. |
| Charged-fermion self-sourcing | [Maldacena, Milekhin and Popov](https://arxiv.org/html/1807.04726v3) derive a semiclassical Einstein–Maxwell throat supported by charged-fermion vacuum energy. | Reuse their source and matching equations; the NSC gauge, light-mass, state and coefficient identification is the task. This is not an NSC-computed solution or generic child formation. |
| Outside geometry replacing local phantom matter | An explicit RS2 effective-brane construction produces black-universe solutions with a positive-kinetic local scalar; the exotic null response is carried by the projected five-dimensional Weyl term. | Structural realization of \(\Gamma_{\mathrm{outside}}\); bulk completion and stability remain open. |
| Matter in another room acting as dark matter | Linearized two-brane gravity makes matter on the other wall gravitate locally as shadow matter and predicts \(25\%\) weaker light deflection at equal Newtonian mass. | Realization of outside particles acting on local fabric, plus a fixed lensing/dynamics discriminator. |
| Local particles, gauge fields, Higgs, and gravity from one spectrum | The Chamseddine–Connes spectral action produces the Standard Model coupled to Einstein plus Weyl gravity from one Dirac operator and imposes high-scale coupling relations. | Imported local-room action; its finite Dirac/Yukawa spectrum remains input unless the recursive fixed point determines it. |
| Local common scale anomaly | The regulated fermionic determinant and its induced spectral action have a derived local Higgs–dilaton form including curvature, gauge, Higgs, and dilaton-gradient terms. | Import the common local Weyl mode; derive only the new relative parent/child sheet mode and its off-diagonal spectrum. |
| Spectral interaction between two geometries | A doubled Dirac operator with two four-metrics and one off-diagonal field derives a canonical relative-metric interaction and an explicit leading FLRW action. | Import the leading interaction; the full exponential relative kernel is a heat-kernel diagnostic, not the complete physical anomaly. |
| Spectral resolution diagnostics | Full Euclidean spectral-action calculations exhibit ultraviolet suppression and local-decoupling features in particular kernels. | Derive the physical causal response before interpreting these features as a propagation wall. |
| Baby universes and inheritance | Black-hole offspring universes and inherited/mutated constants have prior literature, including cosmological natural selection. | Prior hypothesis class; not observational proof. |
| Stable NEC violation | Beyond-Horndeski/DHOST examples can violate the null condition without the simple phantom kinetic instability. | Candidate mathematical source for the transition boundary. |
| Reciprocal energy slopes | Dualities and energy functions invariant under scale inversion are known; the BPS core scales as compression plus volume energy. | Fixed-point scaffold for the parent/child map. |
| Exact effective-action flow | Wetterich’s exact flow equation for the effective action with changing resolution. | Imported mathematics; the new task is the recursive boundary condition. |
| Fermionic Gaussian determinant | Klich’s Fock-space trace equals a one-particle determinant (`arXiv:cond-mat/0209642`). | Implemented for the actual geometry vertex; not a discovery of the identity. |
| Influence functional for gravity | Martín and Verdaguer’s stochastic semiclassical gravity (`arXiv:gr-qc/9904021`). | Executable connection to this Dirac geometry; not a discovery of the influence approach. |
| Dirac heat coefficients | Vassilevich’s heat-kernel expansion (`arXiv:hep-th/0306138`). | Independently reconstructed on the present operator; not new-to-world coefficients. |

Primary anchors include the Tonomura single-electron experiment, PBR (`arXiv:1111.3328`), loophole-free Bell tests (`PRL 115, 250401` / `250402`), the BPS-Skyrme nuclear papers (`arXiv:1309.0820`, `arXiv:1312.2960`), the full gravitating BPS calculation (`arXiv:1503.03095`), the exact black-universe solution (`arXiv:gr-qc/0611022`), and Smolin’s cosmological natural-selection proposal (`arXiv:gr-qc/9404011`). The direct outside-geometry bridge is Bronnikov and Donskoy’s RS2 black-universe construction (`arXiv:0910.4930`, DOI `10.1134/S0202289310010068`). Garriga and Tanaka’s two-wall result (`arXiv:hep-th/9911055`, DOI `10.1103/PhysRevLett.84.2778`) supplies the shadow-matter projection. The local spectral-room authority is Chamseddine and Connes (`arXiv:hep-th/9606001`, DOI `10.1007/s002200050126`). Kurkov and Lizzi’s local Higgs–dilaton anomaly (`arXiv:1210.2663`, DOI `10.1142/S0217732312502033`) supplies the common local scale mode but not independent parent/child sheet resolutions. Bochniak and Sitarz derive the doubled-FLRW spectral action (`arXiv:2012.06401`, DOI `10.1103/PhysRevD.103.044041`) and the general flat-background interaction through third order (`arXiv:2201.03839`, DOI `10.1088/1475-7516/2022/04/055`). The spectral propagation-wall authority is Kurkov, Lizzi, and Vassilevich (`arXiv:1312.2235`, DOI `10.1016/j.physletb.2014.02.053`). The fermionic scale-anomaly origin of the bosonic spectral action is Andrianov and Kurkov (`arXiv:1001.2036`). The exact effective-action flow is Wetterich, DOI `10.1016/0370-2693(93)90726-X`. The Weyl-anomaly identity used for the invariant partition function is `arXiv:1106.3263`, DOI `10.1007/JHEP10(2011)001`. That paper distinguishes its normalization scale from the spectral cutoff; it is not a licence to interchange raw finite determinants with a proper-time Hessian. The scalar recursive square-root response is established surface Green-function mathematics, for example Nemec, Tomanek and Cuniberti, Appendix A.3 (`arXiv:0711.1088`). The BFK gluing formula (`arXiv:math/0304347`) assumes compactness, a product collar, and spectral conditions; those assumptions are not a proof for the present noncompact warped throat. Loss of ellipticity for Dirac Hamiltonians at horizons is documented by Finster and Röken (`arXiv:1512.00761`). Floquet/periodic confinement is established ODE spectral theory; the periodic-throat gap is a project application, not a new-to-world mechanism. Klich’s fermionic determinant identity (`arXiv:cond-mat/0209642`) and Martín–Verdaguer influence methods (`arXiv:gr-qc/9904021`) are imported; the project contribution is their executable connection to this geometry, not those identities.

## What is established and what is not

The latest source/state calculation also reuses:

- [Flachi–Moss–Toms](https://arxiv.org/html/hep-th/0103138v2) for massless interval determinants, conformal reduction and Wilson-line mechanisms. The project evaluates the finite curved interaction while retaining the distinct bulk determinant and boundary cocycle.
- [Dunne](https://arxiv.org/html/0711.1178v1) and [Kirsten–McKane](https://arxiv.org/html/math-ph/0403050v1) for general boundary determinants. Their formulas identify the conditional holonomy minimum of the current positive radial partners.
- [Krishnan–Raju](https://arxiv.org/html/1605.01603v1) for the Einstein/GHY metric variation, used in the compact boundary momentum balance.
- [Montero–Van Riet–Venken](https://arxiv.org/html/1910.01648v4) for charged de Sitter geometry, used to test the common source coefficients before importing a near-extremal seed.

These mechanisms and formulas retain their prior attribution. The curved
source, coefficient and domain matches are scoped project calculations,
not a claim that the generic Casimir, Wilson-line or gluing mechanisms are new.

Observed electron interference proves wave behavior and localized detections. It does not select a unique ontology in which a wavefunction is literally the material medium. QFT removes the classical particle/wave dichotomy by making a particle a quantum state of fields, but the wavefunction remains the amplitude assigned to a state, not an experimentally isolated second substance.

Exact black-universe solutions show that a black-hole exterior and expanding interior can belong to one regular four-geometry. They do not show that our observed universe is such an interior. A curvature singularity or an infinite limit alone also does not encode a child domain; the transition map and its conserved charges must be supplied.

The spectral action unifies known local particles, gauge fields, Higgs, and gravity as spectral data of one Dirac operator. Its finite Dirac/Yukawa matrix still contains measured fermion masses and mixings as input in the imported construction. Nested-Space completion is stronger: the combined functional flow and recursive fixed point must determine that finite spectrum with the same \(\Theta\). Until it does, the spectral action is an imported unifying scaffold, not the completed derivation of all particle masses from the room.

## Longer-term universal-action conjecture

One retained candidate programme seeks to construct and test one object

$$
(S_{\mathrm{one}},\Theta,\mathcal{T})
$$

such that the same fields, operator coefficients, boundary kernel, and scale map jointly produce:

1. finite stationary electron, proton, and neutron charge sectors;
2. quantum amplitudes and detector couplings of those same configurations;
3. nuclear masses and reaction \(Q\)-values without inserted nucleon masses;
4. a self-gravitating branch using those same derived matter coefficients;
5. a stable trapped-to-defocusing-to-expanding transition; and
6. a non-tunable dimensionless prediction identifying our cosmology with the child side of that transition.

The individual ingredients have prior work. A completed and independently tested invariant closure would be a candidate project contribution; a priority claim would still require comparison with that literature. This is a longer-term conjecture, not the sole project-level target or an admission prerequisite for the exterior atlas. Particle and nuclear closure must not block a finite geometric, causal or inheritance question.

## Compute admission rule

A new scientific calculation must identify at least one of:

- a coefficient shared by two previously separate scales;
- a stationary configuration of a declared candidate action that closes a missing connection;
- a prediction whose datum was not used to determine \(\Theta\);
- a parent/child attachment or response with a compatible source, state and conservation law, including a scoped stability result; or
- a dimensionless discriminator against ordinary cosmology; or
- a new operator, domain or state compatibility condition that changes which physical calculation is admissible; or
- a causal, inheritance or observer map that distinguishes candidate exteriors, or a demonstrated observational degeneracy with a named possible discriminator.

Before execution, its compact record names:

1. **Borrowed result:** version, equation/theorem, assumptions and existing implementation or evidence.
2. **Missing connection:** the question the source did not answer and why it matters to the atlas.
3. **Decision and unlock:** which outcomes change a connection or enable a subsequent question.
4. **Finite domain:** regions, observer, state, time/scale window and cutoff. An infinite limit may be imported as a theorem; it is not an unbounded numerical work queue.
5. **Budget and stop:** explicit compute bound and outcomes that stop or redirect this branch.
6. **Verification and attribution:** independent comparison, error statement and focused prior-art search before claiming novelty.

The [programme dependency graph](nsc-cosmological-atlas.md) makes item 3
operational: identify prerequisites by name, the unresolved input or compatibility
condition, and which subsequent calculations become justified under each outcome.
An equation may already be available while its physical inputs are not. Prefer
the smallest connection that unlocks useful consequences or excludes a branch.
This research graph is acyclic; physical recurrence is a separate model question.

Evaluating imported formulas to assemble the atlas and checking adapters are
permitted integration work, labelled as such. They do not become scientific
discoveries because tests pass. Hold a candidate's theory and coefficients
fixed during perturbation tests; changing them starts a different candidate.
Do not insert a target dark-sector fraction, an inherited scale ratio, or a
desired antimatter/time interpretation to select the exterior. Compare with
local observations only after the candidate has supplied an observable.

An incompatible junction or a null response is useful only within the stated
assumptions. It retires that attempted connection, not every nested geometry.
Search coverage and unread sources stay explicit: a bounded literature search
does not establish worldwide priority.

Standalone repetitions of interference, a known soliton, a published nuclear table, a black-universe metric, or a generic test suite are not project results.

The [atlas v2 application](nsc-cosmological-atlas.md) records
the first two decisions under this rule: import Brill–Lindquist spatial data
to supply a common parent, and stop the constant-tension/comoving-radiation
attachment after its scoped incompatibility check. Neither is claimed as
new-to-world mathematics. Moving-boundary dynamics and matching to the
nonspherical parent require their own source/domain admission.

The [v3 application](nsc-cosmological-atlas.md) closes that
motion-only question through imported conservation, and imports a distinct
signed de Sitter wall as a finite control. Published N-body gluing now supplies
a separate parent family with spherical collars; recreating its existence
proof or elliptic solver is not the next research target. Its missing connection
was complete wall-body Cauchy data, the chosen Schwarzschild side, and a source
regularity class covered by the theorem. Neither matching a mass parameter nor
passing a local shell check closes that interface.

The [v4 application](nsc-atlas-connection-slice.md) supplies that
interface in one model: a regular cap, fixed wall, correct vacuum bridge and
AE end with bulk and distributional constraints checked. The gluing theorem
acts on smooth vacuum ends only; identical open overlaps preserve the sourced
body. This is an application of published mathematics plus local constraint
matching, not a new existence theorem. A finite observer domain demonstrates
geometric degeneracy for assemblies sharing the retained vacuum initial data.
Global shell evolution, formation, stability and an observable that breaks
that equivalence remain separate admission questions. No global smoothness
requirement is imposed on the wall merely to use the end-only theorem.

The [v5 application](nsc-atlas-causal-classification.md) borrows the
isolated wall's causal classification and null coordinates, locates the v4
slice, and records access for named interior, wall and parent observers. It
distinguishes early geometric paths from excluded new emissions after the
slice, and the parent's escape horizon from the central de Sitter horizon.
This closes a project integration question, not a new causal-classification
theorem. The reference's singular past and the glued parent's unknown global
future remain explicit. A response calculation must next specify the field,
state and wall transmission law, or justify another communicating attachment;
an available null path alone does not supply an observable imprint.

The [v6 application](nsc-atlas-field-transfer.md) specifies a
classical conformal scalar, zero earlier interior data and one incoming
characteristic on the exterior wall side. Borrowed distributional matching
and conformal propagation give a conditional central waveform. This closes
the state-to-observer interface for the declared input, not its remote
preparation or a global quantum state. The shell cancellation does not remove
smooth Schwarzschild scattering. A later persistent-response comparison must
use its own field equation and wall law; a global exterior calculation needs
admissible initial data and the correct causal branch. A waveform alone cannot
identify geometry while its incident state remains freely adjustable.

## Finite exterior preparation admission, 12 October 2026

The [preparation owner](nsc-atlas-field-preparation.md) closes a bounded part of
that v6 gap. This is an application of borrowed field equations and standard
integral-equation estimates; no worldwide-priority claim is made.

| Admission item | Fixed decision |
| --- | --- |
| Borrowed | The v5 isolated signed wall and null clocks; v6 conformal reduced-field matching and regular-centre propagation; Dolan–Ottewill 1106.4318v1 Eq. 11 scalar potential; standard contraction principle |
| Missing | Independent admissible data on the retained white-hole-side geometry, rather than an arbitrary wall trace |
| Prerequisites | `geometry-baseline`, `boundary-classification`, `wall-transmission`; keep the right-parent causal exclusion |
| Decision | Does a finite characteristic strip away from the singularity define a unique transmission problem? |
| Unlocks | `parent-preparation` supplies a conditional local preparation; `wall-trace` can then be evaluated, enabling `observer-response` and a later `observer-inference` comparison |
| Domain | Fixed H ell=1, kappa ell=2, M/ell=0.4; real classical massless l=0 conformal scalar; wall R/ell=0.45–0.55; exterior R/ell>=0.4; zero earlier cap data |
| Budget | Existing adapters only, at most 513 samples and 60 CPU seconds per scientific record; analytic transfer bound, no PDE solve, parameter sweep, backreaction or formation simulation |
| Stop or redirect | Stop if the retained-side or no-return conditions fail, the corner data conflict, or the contraction bound fails to certify this strip. Report the missing condition; do not retune the theory. A failed bound alone would not prove nonexistence. |
| Independent check | Recover the scalar potential from -R_UV/R; differentiate the wedge integral; bound its area analytically independent of clock quadrature; compare existing 32/64 clock evaluations and unit rescaling |

The analytic operator norm is at most 0.138821, which certifies a unique finite
reduced solution. A numerical waveform from those independent data remains
unevaluated. A global physical source for the characteristic preparation,
formation of the geometry and quantum-state selection remain open. Any proposed
solver must identify the specific observer uncertainty it resolves and obtain a
separate finite admission here before execution.

## Prepared scalar response admission, 12 October 2026

This separate admission succeeds the finite preparation certificate above.
The [prepared-response owner](nsc-atlas-prepared-response.md) supplies the exact
constant-extension representation and an independent analytic positive-tail
bound before numerical propagation. The earlier certificate remains a record
of what had been established without evaluating a waveform.

| Admission item | Fixed decision |
| --- | --- |
| Borrowed | The v7 retained IV strip and source, v6 conformal wall matching and regular-centre response; standard Volterra integration, SciPy DOP853, derivative-informed cubic Hermite interpolation, PCHIP and clipped composite trapezoidal integration; Bessel-series manufactured control |
| Missing | The independently prepared exterior data have a certified transfer operator, but its central waveform has not been evaluated or separated from numerical error |
| Prerequisites | Keep the v7 full retained-domain and no-return certificates, compatible zero bottom and cap data, fixed positive C3 compact input, and right-parent causal exclusion; the finite IV input is not a drive from P |
| Decision | What central signal does this independent preparation produce, and does its post-input tail agree with the analytic positive bound under refinement? |
| Unlocks | Close `wall-trace` and `observer-response` for this preparation; support a subsequent source-fixed `observer-inference` comparison without asserting that one waveform uniquely identifies an exterior |
| Fixed domain | H ell=1, kappa ell=2, M/ell=0.4; real classical massless l=0 conformal scalar, no additional wall coupling, linear probe; wall R/ell=0.45–0.55, exterior R/ell>=0.4; the existing finite arrival window only |
| Representation | Constant U extension of the exact retained-wedge operator for storage; evolve z=chi-g, integrate only to min(U,U_W(V)), and evaluate J=chi_V from the right-hand side rather than differentiating saved output |
| Budget | At most 257 direct geometry knots per grid and 300 per adapter call; spatial grids 129/257/513 and 257 output points; at most 6,000 right-hand-side calls per solve, 10 million kernel evaluations and 60 CPU seconds for the complete recorded experiment; no sweep, new general PDE framework or theory retuning |
| Acceptance | Final spatial and geometry waveform differences at every output point at most 1e-6+1e-4 times the absolute fine value; temporal tolerances tighter than this threshold; post-input J refinement and independent reference-quadrature differences below 10% of the analytic positive lower bound |
| Independent check | Positive analytic bracket (1-q_*) times B(V)<=J(V)<=B(V), using an exact transverse kernel integral and Gauss 32/64 double quadrature at three fixed post-input points; zero-potential and constant-potential Bessel controls; unit scaling, linearity, compatibility and finite-clock rejection |
| Stop or redirect | Stop at the budget, on failed retained/no-return assumptions, unresolved tail/error separation, or disagreement with the independent bracket. Report the unresolved interface rather than changing parameters to rescue a signal. Other atlas routes remain open. |

The post-input derivative trace has the exact positive lower bound
(1-q_*) K_min^2 (U_a-U_0)^2 G/2, approximately 6.53246e-5 epsilon in the
recorded geometry, where G=64 epsilon (V_b-V_a)/315. This analytic check is
independent of the numerical field solution. Its decimal value uses the
existing geometric quadrature diagnostics, not interval arithmetic. The
result will report a finite classical waveform and numerical error evidence;
it will not claim a global parent source, late-time persistence, energy flux,
unique exterior inference or worldwide novelty of scalar scattering.

### First bounded attempt and event-aligned integration clarification

The initial 129-node prototype stopped at its 6,000-right-hand-side-call cap
near normalized source coordinate s=0.575; no completed waveform was accepted.
The [response owner](nsc-atlas-prepared-response.md#first-execution-stop-and-algorithm-clarification)
records the failure. The partial spatial integration cell changes at known
moving-wall/node crossings, so the ODE driver's steps need to respect those
discrete events.

Before the redirected attempt, split integration at the roots U_W(V)=U_i
of the existing monotone Hermite wall interpolant and at both input support
edges. Use a bracketed scalar root finder and carry the correction state
continuously between DOP853 segments, initially trying each segment's full
length. The source, theory, spatial integral, grids, tolerances, output window
and acceptance tests remain fixed. The 6,000-call cap is the sum across all
segments of one solve; the 10-million-kernel and 60-CPU-second caps still apply
to the complete new record. This is a recorded algorithm correction after a
bounded failure, not permission to retune the physical model or reset a cap
within a solve. Stop again if these limits or the original accuracy checks fail.

## Observer preparation inverse admission, 12 October 2026

The v8 signal now permits an inference question before any additional wave
evolution. This admission tests the fixed-model inverse, rather than fitting
another signal. The original compact positive pulse remains unchanged evidence.

| Admission item | Fixed decision |
| --- | --- |
| Borrowed | The v7 retained wedge, positive Schwarzschild kernel, no-return wall law, analytic norm bound q_*<1 and monotone central clock; v8 evaluated response; standard contraction/Neumann-resolvent mathematics |
| Missing | Whether a complete compatible central trace uniquely determines its characteristic preparation in the fixed model, and which preparation restrictions are needed before using that trace to compare geometries |
| Prerequisites | Same scalar theory, wall, retained domain, zero bottom and cap; complete observer history from the initial endpoint and its identified proper-time map; preserve the causal exclusion of the right parent |
| Decision | Prove or reject a bijection between C2 source and wall traces vanishing on a fixed initial open interval, with direct PDE, boundary and corner checks; give the associated central-trace inverse and its limited stability statement |
| Unlocks | A preparation-constrained response-range comparison or a source-linked independent channel. An unrestricted inverse does not justify comparing freely fitted waveform families as geometry evidence. |
| Domain | Original wall R/ell in [0.45,0.55], exterior floor 0.4, H ell=1, kappa ell=2, M/ell=0.4, classical massless conformal l=0 probe and the original finite central clock interval |
| Budget | Proof and focused primary-source review only; one independent analytic review, elementary derivative/corner checks and evidence packaging. No wave solve, kernel quadrature, inverse iteration, parameter scan or new floating-point scientific estimate. |
| Independent check | Differentiate the wall-prescribed integral to recover f_UV=-Kf and f_U at the wall; recover source and bottom; check the L-switch at U_a; compare with the independent forward-map resolvent identity |
| Stop or redirect | Stop any broader claim if the inverse needs points outside the retained wedge, fails corner regularity, or leaves the permitted source class. Without a second independently compatible geometry, report only fixed-model invertibility and conditional range-overlap criteria. |
| Claim boundary | A scoped application of established integral-equation methods, not a new general theorem, numerical source reconstruction, physical source generator, universal degeneracy, geometric reconstruction or priority claim |

The three foreseeable failures are a boundary/domain mismatch, an inverse
that leaves the declared preparation class, and confusing one trace with a
controlled family of source-response measurements. Direct differentiation,
an explicit source-class test, and the focused inverse-wave literature
comparison address them. A failing check changes the scientific conclusion;
it does not authorize retuning the model or launching a larger simulation.

## Conformal probe stress and shell admission, 13 October 2026

The fixed-model inverse leaves freely chosen preparation as an inference
limitation. Before introducing a source-selection law or another measurement
channel, check whether the existing transmitted probe has compatible stress
at its wall. This is a finite compatibility application of published
scalar-tensor junction conditions, not a new gravity theory.

| Admission item | Fixed decision |
| --- | --- |
| Borrowed | Aviles, Maeda and Martinez, [1910.07534v2](https://arxiv.org/html/1910.07534v2), Eqs. (2.10)-(2.14), (3.42)-(3.46), and the non-null distributional derivation; the existing signed pure-tension wall and conformal scalar transmission; the v7/v8 smooth compact preparation |
| Missing | Does the order-amplitude-squared scalar stress require an additional singular surface source or tangential energy exchange that the fixed wall action does not supply? |
| Prerequisites | Same four-dimensional Jordan-frame geometry, normal from cap to exterior, continuous scalar, piecewise sufficient differentiability, fixed wall k and xi=1/6; use the original C3 preparation rather than infer full four-dimensional stress regularity for every C2 inverse input |
| Action and normalization | Fix the scalar kinetic coefficient to one, zero scalar potential, f(phi)=1/(2 kappa_4)-xi phi^2/2 with kappa_4=8 pi G, scalar-independent background vacuum stress, and material wall L=-sigma; k=4 pi G sigma is distinct from kappa_4. No direct scalar-wall interaction. |
| Perturbative meaning | Zero scalar background, phi=epsilon phi_1; evaluate its quadratic Hilbert stress on the unchanged background. Include xi G_mu_nu phi^2 exactly once. Assume f>0 for the weak-field interpretation; this is not a quantitative backreaction bound or solved Einstein data. |
| Decision and unlock | Separate the singular scalar surface coefficient, smooth tangential flux jump and remaining bulk Einstein residual. Compatibility permits a subsequent preparation-energy or sourced-geometry question; incompatibility identifies the missing surface law without adjusting the present wall. |
| Finite domain | The existing regular timelike wall segment R/ell in [0.45,0.55], finite retained exterior and regular cap; local jump identities with H ell=1, k ell=2, M/ell=0.4 fixed. No formation, global energy, late-time or stability conclusion. |
| Budget | One analytic derivation, one independent derivation/review, and exact symbolic checks within 60 CPU seconds. No field re-solve, numeric derivative of the saved waveform, parameter sweep, new PDE solver or backreaction integration. |
| Independent check | Compare the distributional improved-stress coefficient with the imported variational scalar-tensor junction; compute tangential flux using one-sided regular traces and the common-normal Gauss-Codazzi conventions. Retain zero-field, normalization and orientation checks. |
| Stop or redirect | Stop on an undefined distributional product, insufficient regularity, an indispensable unspecified wall coupling, or a failed identity. Report any bulk or normal-force condition left unresolved; shell compatibility alone cannot establish a self-consistent radiating geometry. |

The likely failures are confusing two stress-accounting conventions, using
insufficient regularity at the interface or centre, and treating a surface
cancellation as a solved bulk response. The action convention, one-sided
trace calculation and separate bulk residual are the checks against these
failures. A new amplitude estimate or coupled solution needs a later admission.

## Characteristic gravitational initialization admission, 13 October 2026

The probe's surface compatibility leaves its gravitational initial data
unspecified. A generic finite stress bound would not close that gap: smoothness
on the compact wedge already supplies boundedness, and the amplitude remains
free. First ask which null constraint correction the fixed preparation forces.

| Admission item | Fixed decision |
| --- | --- |
| Borrowed | Spherical null Raychaudhuri/Einstein constraints and intersecting-characteristic data, including Madler-Gannouji-Gallo [2503.24162v3](https://arxiv.org/html/2503.24162v3), Eqs. (16), (29a), (86)-(92); the existing canonical conformal stress; Aviles-Maeda-Martinez [1910.07534v2](https://arxiv.org/html/1910.07534v2), Eqs. (2.11)-(2.14), (2.21)-(2.24), for an independent smooth-vacuum frame check |
| Missing | Compatible gravitational data on the admitted source and zero bottom legs; leaving both the radius and null normalization unperturbed has not been justified for a nonzero pulse |
| Fixed prerequisites | Original v7 C3 reduced-field preparation, v8 leading field, xi=1/6, canonical positive kinetic coefficient, zero scalar background, unchanged corner Schwarzschild geometry and both corner expansions; bottom scalar data zero. Preserve the right-parent causal exclusion. |
| Gauge and order | Hold the null lapse at its reference value along the source leg, equivalently keep the corner-anchored background areal-radius affine parameter affine there. Solve the order-amplitude-squared radius constraint correction with zero corner correction and derivative; do not also prescribe the original radius all along the leg. |
| Decision | Determine the required initial radius and affine radial-derivative correction and whether it vanishes after the compact scalar trace ends. Give an analytic sign or obstruction where possible and evaluate the coefficient for the recorded preparation. No amplitude is selected to make the answer small. |
| Unlocks | A constraint-consistent input for a later coupled exterior/wall problem, or a precise incompatibility requiring a changed initial-data choice. A source-leg result alone does not determine a mass shift, corrected wall, interior waveform or stability. |
| Finite domain | Existing U=U0, Va<=V<=Vb in retained sector IV, R/ell in [0.4,0.55], and zero bottom Va, U0<=U<=Ua. No new future domain, global source or full bulk evolution. |
| Budget | One-dimensional constraint quadratures and exact algebra only; at most 129 displayed source samples, Gauss orders 32/64 split at the fixed support edges, at most 20000 scalar integrand evaluations and 60 CPU seconds per complete evidence creation. Reuse 32/64 geometric adapters. No PDE solver, field re-solve, parameter sweep or bulk derivative campaign. |
| Independent checks | Derive the constraint from null Raychaudhuri and the full improved stress; compare direct stress integration with integration by parts and the smooth-vacuum Einstein-frame focusing law, including the change in affine parameter. Check zero input, corner/bottom compatibility, support edges, unit scaling and rejection of changed theory/domain. Numerical refinement is a diagnostic, not a rigorous enclosure. |
| Acceptance | Exact constraint and frame identities agree; independent 32/64 evaluations of the dimensionless correction coefficients agree within 1e-8 absolute plus 1e-8 relative at the displayed points; the last coefficient and any sign statement have an independent analytic check. Report the remaining characteristic and wall equations explicitly. |
| Stop or redirect | Stop on inconsistent corner data, unanchored affine normalization, a singular frame map, disagreement of independent constraints, exhausted budget or unresolved quadrature. If gravitational initialization needs an unspecified physical datum, expose it instead of silently tuning it. |

The likely failures are overspecifying the radius and lapse together,
mistaking a gauge-labelled expansion for a mass or observer signal, and
using the vacuum frame transformation across the material wall without its
transformed sources. The anchored characteristic data, restricted output and
exterior-only independent check address these three risks. This is an
application of borrowed constraints to the existing preparation, without a
claim of worldwide mathematical priority.

## Invariant source-leg mass admission, 13 October 2026

The initialized null radius is intrinsic data. The next missing connection is
an invariant spherical gravitational response to that same fixed preparation.
The geometric mass also needs a transverse radius derivative; its corner datum
is already fixed, and the leading transverse scalar is determined by the
existing characteristic field equation. Check closure before numerical work.

| Admission item | Fixed decision |
| --- | --- |
| Borrowed | An--Lim [2005.04090v1](https://arxiv.org/html/2005.04090v1), mass definition (1.3), mixed radial equation (2.4), mass transport (2.14)-(2.15); Madler--Gannouji--Gallo [2503.24162v3](https://arxiv.org/html/2503.24162v3), geometric mass (21) and transverse transport (84c)/(91); the already imported improved conformal stress and smooth-exterior frame map. |
| Missing | The Jordan-frame geometric mass coefficient and its connection to the initialized source leg. Minimal-matter transport cannot be inserted unchanged into the Jordan-frame equation. |
| Fixed prerequisites | Same original C3 reduced-field source, zero bottom trace, canonical conformal scalar with xi=1/6 and zero potential, no additional wall interaction, fixed corner radius, both expansions and lapse, and fixed source-leg lapse. Reuse the sealed initialization evidence and leading field equation. |
| Decision | Establish whether the mixed equation closes with those data; calculate the order-amplitude-squared geometric mass response along the same leg, including whether a response persists after the compact trace ends. Report its frame, coordinate invariance and sign with their actual assumptions. |
| Unlocks | An invariant exterior consequence of the admitted preparation and compatible gravitational data; a concrete input for later sourced wall matching. This alone gives no wall motion, central-observer correction, stability or asymptotic energy. |
| Finite domain | Original U=U0, Va<=V<=Vb in retained sector IV, with the original corner and zero bottom. No new future domain or finite amplitude selection. |
| Budget | Analytic adaptation of borrowed equations and one-dimensional quadratures only. At most 129 displayed samples, Gauss orders 32/64 split at source support, at most 100000 counted integrand evaluations including reused geometry/constraint quadratures, and 60 process-CPU seconds per complete evidence creation. No PDE solver, wave re-solve, parameter scan or bulk evolution. |
| Independent checks | Compare geometric mass from mixed-radius transport with borrowed Einstein-frame mass flux and the full frame transformation, including transverse scalar derivatives. Independently check Jordan improved-stress mass transport. Test zero source, fixed corner, support edges, coordinate/unit rescaling, changed theory/data rejection and source/evidence tampering. |
| Acceptance | Exact symbolic identities close with declared data; independent quadrature/geometry refinements and mass representations agree within 1e-8 absolute plus 1e-8 relative at displayed points. Any post-input sign has an analytic check. Numerical agreement is not an enclosure or a finite-amplitude remainder bound. |
| Stop or redirect | Expose any unprovided transverse datum, incompatible corner, singular frame map, normalization mismatch, unresolved refinement or exhausted budget. Do not tune a missing datum, import the exterior-region monotonicity assumption into sector IV, or use a historical Dirac/Bondi result to close this scalar calculation. |

The failure forecast is an incorrect matter/frame normalization, a hidden
transverse-data requirement, or promotion of a geometric mass coefficient to
global energy or completed evolution. Independent metric and stress checks,
the explicit corner, and the finite source-leg scope address those risks.
Known mass laws are borrowed; their application here makes no worldwide
priority claim.

## Scalar normal-force admission, 13 October 2026

The prepared pulse has compatible singular surface stress and tangential
flux, and its initial null geometry now has an invariant mass response.
Before transporting that geometry to the wall, check the remaining local
normal-force condition. This is an application of existing junction and
stress laws, not a new general junction theorem.

| Admission item | Fixed decision |
| --- | --- |
| Borrowed | Full nonminimal Hilbert stress and scalar/material jumps from Aviles--Maeda--Martinez [1910.07534v2](https://arxiv.org/html/1910.07534v2), (2.11)-(2.14), (3.42)-(3.46); normal Gauss constraint and shell evolution identity from Lobo--Crawford [gr-qc/0507063v2](https://arxiv.org/html/gr-qc/0507063v2), (5)-(6). The existing scalar-wall owner fixes conventions. |
| Missing | The quadratic regular scalar normal-stress jump for this conformal probe and constant-tension wall, and its role in the unsolved geometric normal-force condition. Surface and tangential cancellation do not supply it. |
| Fixed prerequisites | Original C3 preparation and no-return window; canonical real massless scalar, xi=1/6, zero scalar background/potential; unchanged scalar-independent vacuum terms and material tension. Continuous induced metric and scalar; sufficient one-sided regularity; common cap-to-exterior normal; [n phi]=k phi, [K_ab]=-k gamma_ab, [K]=-3k. |
| Decision | Project the full stress onto the normal and reduce its jump analytically on the background. Determine whether the prepared field supplies a nonzero quadratic normal forcing, an automatic compatibility identity, or an unprovided datum. Keep any metric and wall correction explicitly unsolved. |
| Unlocks | A checked interface condition for the later sourced geometry/wall initial-boundary problem, or a precise obstruction. If cancellation is automatic, borrow the identity and proceed to transport rather than generating pressure curves. |
| Domain and order | Existing wall R/ell in [0.45,0.55] and original observer/field domain; formal order-amplitude-squared scalar forcing evaluated on the reference geometry. No finite amplitude, new future domain or assertion of an exact unbackreacted spacetime. |
| Budget | Analytic reduction and exact symbolic verification only, at most 60 process-CPU seconds per complete evidence creation. No waveform rerun, numerical differentiation, parameter scan, pressure curve, bulk solver or wall trajectory integration. |
| Independent check | Compare direct Gaussian-normal projection with one-sided product/jump identities and the borrowed Gauss/material-shell bookkeeping. Retain the cap vacuum term and the nonminimal Einstein-tensor term exactly once. Check normal-and-side reversal, zero amplitude, dimensional homogeneity and omission/alteration controls. Do not use the unknown sourced Einstein equations to force cancellation. |
| Stop or redirect | Stop after a proved identity or a named nonzero forcing; do not turn an automatic condition into a new simulation. Stop if undefined products, unprovided traces or incompatible source terms remain. Individual pressure curves or coupled transport require a separate admission. |

The likely failures are a normal/sign convention mismatch, omission of the
cap's vacuum Einstein tensor, or promotion of local compatibility to a solved
trajectory. Oriented one-sided checks, explicit vacuum terms and an unsolved
geometric remainder distinguish these failures. The bounded reading of the
two imported nonminimal-junction papers did not find this evaluated jump;
that observation does not certify worldwide priority.

Outcome: the [normal-force calculation](nsc-atlas-scalar-normal-force.md)
proves zero additional regular scalar normal-stress jump at the selected
coupling. Gauss/Israel and normal distributional conservation fix the plus
force sign in the project's conventions; the opposite sign printed in
Lobo--Crawford (6) is not copied. The resulting geometric condition still
contains unsolved metric and wall-displacement corrections. This branch
stops at the identity. A promising next bridge is exterior-side mass
transport along the existing wall traces, before a full coupled evolution;
it requires a separate admission and an independent transport check.

## Exterior wall-mass closure admission, 13 October 2026

The initialized mass and local wall conditions do not yet connect the field
to gravity on the moving boundary. Test whether a one-sided exterior mass
functional closes from the already available leading traces before selecting
a bulk evolution method. This admission is analytic; a numerical mass curve
needs its own quadrature and error decision after closure.

| Admission item | Fixed decision |
| --- | --- |
| Borrowed | Both null mass laws in An--Lim [2005.04090v1](https://arxiv.org/html/2005.04090v1), (2.14)-(2.15); the smooth-exterior field/frame map of Aviles--Maeda--Martinez [1910.07534v2](https://arxiv.org/html/1910.07534v2), (2.21)-(2.24); full geometric mass transformation in Faraoni--Vitagliano [1401.1189v1](https://arxiv.org/html/1401.1189v1), (7); general-matter mass law in Hayward [gr-qc/9710089v2](https://arxiv.org/html/gr-qc/9710089v2), (2.9d). Reuse the existing source-mass normalization and reading boundaries. |
| Missing | A quadratic Jordan geometric mass coefficient on the exterior limit of the reference wall, with both null fluxes and the full frame correction accounted for. |
| Fixed prerequisites | Same canonical real massless conformal scalar, original C3 preparation, zero bottom/cap data, no-return condition and fixed material/vacuum parameters. Existing chi_W, J=d chi_W/dV and chi_U=0 on the wall; the outward normal is (-Udot,Vdot). Keep the quiet corner mass and the right-parent causal exclusion. |
| Decision | Establish whether the leading scalar trace, background clocks and radial derivatives close the exterior-side mass along the timelike reference wall. Identify any unavailable datum, higher derivative or unsolved metric correction that prevents this. |
| Unlocks | A separately admitted finite mass evaluation and an exterior input to later coupled wall/cap matching, or the smallest missing dynamical equation. This alone supplies no cap mass, wall motion, material-shell mass, corrected observer signal or amplitude bound. |
| Domain and order | Original wall R/ell in [0.45,0.55], retained exterior sector IV, formal quadratic scalar order and fixed preparation/corner. The smooth exterior frame map is not continued across the material wall. |
| Budget | One analytic closure and independent exact-algebra check, at most 60 process-CPU seconds per complete evidence creation. No numerical mass curve, waveform rerun, PDE solver, parameter sweep or enlarged future domain. |
| Independent check | Compare both Einstein-frame null fluxes plus the full frame map with Hayward's Jordan-frame transport retaining the improved stress. Check quiet anchor, zero field, null reparametrization, unit scaling and regularity at beta_+=0. Omitted transverse/flux/frame terms must fail meaningful controls. |
| Stop or redirect | Stop with the exact obstruction if an unavailable datum or metric correction remains, if independent routes disagree, or if the proposed expression is singular at the regular beta_+=0 crossing. Do not alter the state or sources to close it. Numerical evaluation needs an explicit new budget and independent refinement; v8 waveform convergence alone does not certify a quadratic mass integral. |

The three likely failures are a missing second null flux, incorrect frame
or timelike-path accounting, and confusion between one-sided mass and a
coupled wall solution. Two independent mass routes, regular wall variables
and an explicit unsolved cap/wall remainder distinguish those failures.
