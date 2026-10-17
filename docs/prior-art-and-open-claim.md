# Prior art and admission of new research

## Cap transport and radial-quadrature separation, 16 October 2026

The completed mass ledger below leaves a definite remainder. Before another
physical calculation, connect it to the radial constraints and reconstructed
field equations using the same full direction as the ledger. The mathematical
inputs are the existing mass-transport identity, variation of constants for
a scalar linear ODE, and differentiation of moving interfaces.

This step admits exact symbolic controls only: the general local identity,
moving-interface integration, the wall conversion, and affine panel-error
telescoping with its full derivative. The controls share 20 active process-CPU
seconds, with no numerical field samples, retained-state evaluations, evolution,
or grid scans. Independent derivation and omission controls check the signs
and the terms that would otherwise disappear in a wall-only comparison.

**Outcome:** four exact controls pass, with zero field samples and 0.640625
active CPU seconds. The [mathematical owner](nsc-atlas-finite-amplitude.md#separating-radial-quadrature-from-reconstructed-field-transport)
now specifies the shadow radial problem, full-direction field residuals,
moving-interface signs and differentiated panel-error sum. Grok independently
confirmed the corrected full-direction balance. No retained cap was evaluated.

The shadow radial solution is a mathematical comparison using the same
reconstructed fields and centre lapse. It is not a replacement cap or a new
physical candidate. Its numerical evaluation and any corrected evolution
remain unadmitted until the exact bridge, all formula-switch interfaces,
one-sided branches, finite evaluation budget, and independent error check
are specified. In particular, defining a quadrature derivative by subtracting
the two total derivatives cannot then serve as an independent closure check.

The likely failures are omission of moving-interface terms, reuse of the
producer's attenuation in an exact-flow error sum, and circular verification
of the quadrature derivative. Exact omission controls address these before
any further physical run.

The next implementation component is a read-only field-jet adapter for this
identity. It reuses the existing branch-selected reconstruction and actual
background/potential partials. Two tangent columns carry the full stored
direction at fixed physical radius and the radial query derivative. Callers
specify the owner cell for every query and trace; no side is inferred at an
interface. Axis jets and samples exactly on either implementation switch
remain outside this component's domain.

Manufactured controls are admitted under one cumulative 200,000-point,
20-active-CPU-second allowance. They cover regular polynomial oracles,
simultaneous radius/clock/moment motion, one-sided limiter openings,
material trace derivatives, source-domain guards, immutable input data and
budget failures. No retained cap or shadow radial solve is admitted by
these controls. Failure stops the component for diagnosis before any
larger calculation. The main risks are confusing fixed-radius and material
derivatives, silently selecting a trace side, and replacing the implemented
small-radius or potential formula with a different expression.

**Field-adapter outcome:** 23 manufactured controls pass. All development
executions plus the final suite use 11,632 charged values / 0.21875 active
CPU seconds. Independent Grok source review found no blocking defect in
the two columns, branch selection, domain guards or copied inputs. Both
field columns have nonconstant independent difference controls; coefficient
jets have regular polynomial oracles. This component supplies field inputs
only. The retained-state shadow flow, switch/interface inventory and error
decomposition have not been evaluated.

The interface inventory now receives exact algebra controls, under a separate
20-active-CPU-second allowance with zero field samples. Prove continuity of
the primitive-derived scalar at every reconstruction face, and retain the
nonzero differences between the implemented potential and small-radius
formula branches. Derive threshold velocities from the full fixed-radius
field derivative. These controls admit no root search, retained cap or
radial integration; unresolved tangencies or coincident branch motion must
be identified before that numerical comparison is designed.

**Interface outcome:** all three exact controls pass with zero field samples
and 1.71875 active CPU seconds. The mathematical owner records continuity,
nonzero switch gaps and full-direction root speeds. An initial Grok review
misstated the small-radius leading order and treated sub-roundoff gaps as
zero; direct series algebra corrected both, and the follow-up review agreed.
There is no retained-state root inventory or shadow-flow result yet.

**Selected reference, after rethinking:** a literal integration of the
producer's truncated potential requires threshold-root isolation. Simply
dropping its tiny jumps would instead need bounds on root velocities,
including nearly tangent crossings. Use the smooth analytic potential of
the declared model for the reference radial flow, while leaving every native
quadrature node and its actual branch derivative unchanged. The panel error
then includes both radial quadrature and potential implementation error;
an explicit `R_x*(b_star-b_Q)` wall term keeps the balance exact. This avoids
root-location state and does not change the candidate or its evolution.
The literal-branch interface identities above remain valid comparisons,
but a potential-root inventory is no longer on the selected execution path.

The full-direction identity and field adapter are complete. Active work is
the smooth-potential reference and affine radial panel map; next are cap
interface limits, composition and independent error estimates. This component
uses one caller-owned budget and returns only its map, derivative and work
count. Exact polynomial/constant-coefficient solutions, independent fixed-r
variation, moving endpoints and born-panel controls are admitted under a
cumulative 200,000-point / 20-active-CPU-second allowance. No retained cap,
trajectory, root search or parameter scan is included. A failed map identity,
unresolved interface or exhausted allowance stops the comparison locally.

**Reference component outcome:** ten manufactured controls pass, using
9,356 charged values / 0.03125 active CPU seconds across their development
executions. The wall-conversion control was extended to independent native
and reference densities; its full run and corrected symbolic-form assertion
use 1.0 additional active CPU second, with zero field samples. Grok source
review found no concrete blocker in the potential tails or moving map.

The existing field-adapter allowance also covers the reference's named
small-radius limits and regular axis coefficients. All 27 current controls
pass across the complete run and one corrected axis-roundoff assertion.
Cumulative field-adapter development is 18,526 charged values / 0.609375
active CPU seconds. The axis check keeps five exact zero coefficients and
allows 16 eps relative roundoff in `a_r(0)` reconstructed from the spline.
Both switch traces have independent 60-digit cubic-profile checks. No cap
composition or retained-state comparison has been run.

The next component composes the native and smooth-reference affine maps
over the same frozen cap partition. Reuse the native cut selector, preserve
its backward owner choices, and split only the reference at the fixed
small-radius formula switch. Keep any first-order born pieces. Construct
the error derivative from local map differences, numerical entering masses
and independently accumulated reference attenuation; the difference of the
two final mass derivatives is a check only.

Manufactured cap controls share 1,000,000 charged points / 40 active CPU
seconds. They compare native primal/tangent results with the existing
producer, use an exact zero-coupling constant-field mass, check moving
interfaces and switch coincidences, and compare reference composition with
independent fixed-radius variation. No retained terminal cap, source fit,
trajectory or grid scan is admitted here. Failed ownership, telescoping,
primal/tangent comparison or an unresolved endpoint stops the component.

**Composition outcome:** seven manufactured controls pass. Cumulative
development uses 472,439 charged values / 6.625 active CPU seconds, including
the first failed assertion and its focused correction. That assertion tried
to resolve an order-eight map error already at roundoff; the nonconstant
control now uses order four to exercise a nonzero error derivative. No
equation or comparison tolerance was changed. Grok independently confirmed
the affine chain rule, local-error derivative and downstream weighting.
The recurrence bound treats the computed maps as inputs; reference integration
and coefficient errors remain outside it. Native born panels exactly on a
formula switch retain the existing tangent's rejection, while reference-only
born pieces have named side limits. No retained cap or new evolution was run.

Before the residual-integral implementation, extend the existing exact
transport controls to interfaces coincident at the wall or centre. Their
backward ordering determines which traces survive. An interior cut reaching
the wall still has its moving jump in the limit; selecting jumps solely by
the instantaneous inequality `0 < r < R` would erase that term. These exact
polynomial controls use the existing 20-active-CPU-second algebra allowance
and zero numerical field samples. They admit no physical cap evaluation.

**Endpoint outcome:** all seven transport-algebra controls pass with zero
field samples and 0.8125 additional active CPU seconds (2.453125 cumulative).
They reject a merged velocity for coincident cuts and the wrong limiting wall
trace. Grok independently confirmed the ordered-interface treatment.

The next numerical component integrates one declared smooth reference piece:
mass, lapse, three field-residual channels and attenuation of an incoming
jump channel. Its inputs are explicit field jets, the same centre lapse and
one caller-owned budget. It does not discover interfaces or assemble the
final cap balance. Manufactured polynomial and constant-field controls,
independent fixed-radius variation, born-piece and finite-domain guards
share 200,000 charged values / 20 active CPU seconds. A sampled metric margin
is reported as sampled, not a proof between solver queries. No retained cap,
trajectory or scan is admitted; a failed identity stops this component.

**Smooth-piece outcome:** all five controls pass on their first run, using
9,828 charged values / 0.046875 active CPU seconds. Grok's independent source
review confirmed the weighted residuals, regular centre start, attenuation,
born-piece identity and shared budget. The component reports only sampled
metric margins and requested solver tolerances. Full ordered-interface
assembly and a controlled retained-state comparison remain outstanding.

The full-cap assembler now has a defined bounded scope: reuse frozen native
panels and their reference switch pieces, propagate mass/lapse/residual
states, and add every ordered interior interface event, including born pieces
at coincident radii. Authenticate the centre lapse against the frozen field.
Construct the final comparison from the independent local-map error derivative,
wall conversion, bulk channels and jump channel; compare it with the original
native tangent balance. The two derivative constructions remain independent.

Manufactured full caps and exact piecewise-density controls share 2,000,000
charged values / 60 active CPU seconds. They cover both amplitudes, nonzero
interface motion, separate coincident velocities, wall and axis events,
unchanged inputs and deliberate centre/partition corruption. The full-balance
domain retains the native wall-sample switch guards. No retained state, new
trajectory or scan is admitted by these controls. Missing jumps, mismatched
wall traces or independent identity failure stop the assembler before a
physical calculation. Numerical solver tolerances remain distinct from an
error estimate; controlled retained-state evaluation will need its own budget.

**Full-cap outcome:** all nine manufactured controls pass. Development uses
890,431 charged values / 12.328125 active CPU seconds, including the first
small-radius constant-field rejection, its focused correction, and the final
mass/jump assertions. Decimal moment/rate rounding selected different tied
limiter representatives in the two arithmetic paths. The exact constant-field
oracle now uses dyadic data; the common-branch guard remains unchanged. The
coupled-rate control recovers the original ledger balance `-0.12288686492081752`
with assembled defect `-8.187894806610529e-16`. Grok independently confirmed
the formula and required explicit reference mass/rate/jump agreement checks.

### Next fixed terminal comparison

After the recorder and its manufactured orchestration checks are frozen,
evaluate only the original zero/finite terminal cut 128 at 64 cells and native
quadrature order eight. Reuse the saved profile, wall clocks and incoming
trace, authenticate the existing mass-ledger summary, and reconstruct two
metrics, two RHS directions and two tangents per analysis. Compare each with
the saved order-eight mass ledger before interpreting the new channels.

Use exactly two reference tolerance pairs, `(rtol,atol)=(2e-10,2e-12)` and
`(2e-12,2e-14)`, on those same states and directions. This is a fixed numerical
accuracy comparison, with no time step, spatial grid scan or source fit.
Per analysis there are four full reference balances. Creation plus independent
reconstruction share **48,000,000 charged points / 300 process CPU seconds**;
a later standalone check has **24,000,000 / 150 seconds**. Failure records
the phase, counts, scalar defects and consumed allowance before any redesign.

The intended resolution is `1e-6` in the normalized mass balance, comfortably
below the retained `0.1`-scale remainder. Require reference mass disagreement
at most `1e-8`, its wall-density effect `abs(R_x*a_W*delta_n)` at most `1e-7`,
reference-rate disagreement at most `2e-7`, and jump attenuation disagreement
at most `1e-7`. Native map telescopes retain their derived roundoff bounds.
Each channel's change between the two tolerance settings must be at most
`5e-7`, and the full balance defect at most `1e-6`. Report all checks and
sampled metric margins. These are declared diagnostic accuracy checks, not a
uniform integration-error enclosure or proof of a positive metric between
samples. A failure is not permission to loosen its limit without diagnosis.

Keep full piece/channel data under ignored `.build/results`. The tracked
summary contains the scalar channels, accuracy comparisons, source/input
identities, raw hash and regeneration command. This comparison diagnoses the
existing state; it does not establish a corrected evolution or closure.

The [terminal recorder](../scripts/diagnose_nsc_atlas_cap_transport_balance.py)
implements this fixed inventory. Its 21 manufactured orchestration and
authentication controls pass; cumulative development work is **195,757 charged
values / 1.765625 active CPU seconds** within 250,000 / 20. They exercise
the individual accuracy gates, native-ledger mismatch, tolerance sensitivity,
creation-only output, hash-only input authentication, producer changes and
failure receipts. Grok identified the raw-output boundary as the remaining
recorder risk. The recorder now requires ignored `.build/results` output and
its `--summary` option uses the existing summary tool with the explicit scalar
selection. Controls verify that tracked full-output paths are rejected.

The recorder reconstructs its scientific payload with absolute `1e-8` and
relative `2e-12` arithmetic tolerance. Saved scalar summaries are also checked
against their detailed balances and must independently pass the original
accuracy gates. Schema, inputs, producers, fixed call inventory and tolerance
declarations remain exact. Adaptive integration counts are retained as bounded
work receipts; they are not required to be identical across platforms.

The first attempt with producer `4ce93f9` stops on the zero case before any
reference flow: `native tangent partition/limiter differs`. It uses
1,667,473 points / 2.921875 CPU seconds and one metric/RHS/tangent. The exact
branch guard remains in force. A single prefix-only diagnosis may reconstruct
that zero metric/RHS and compare the adapter's branches to the saved tangent,
with at most 500,000 points / 30 CPU seconds, carried inside the original
48-million-point / 300-second allowance. It performs no reference integration,
time step or grid comparison. The candidate cause is component-first versus
direction-first rounding in branch selection; that must be established before
changing the adapter or repeating the retained comparison.

**Diagnosis:** 27 minus-characteristic branch choices differ at zero-valued
secants. The 128 panel endpoints and owners agree. Rebuilding in the native
130-component basis reproduces every saved branch; contracting input rates
before reconstruction changes cancellation roundoff. The prefix uses
199,435 points / 2.515625 CPU seconds. Grok independently confirms the cause
and the component-first correction. The adapter now contracts reconstruction
jets only after branch selection, preserving the native tangent and exact guard.

All 44 affected manufactured checks pass, including the previously failing
decimal constant cap. Cumulative full-cap controls are 1,327,411 points /
18.25 active CPU seconds; field controls are 53,753 / 0.734375; radial-composer
controls are 691,326 / 9.40625. Each remains within its original allowance.
The retained comparison may resume after freezing the corrected producer,
carrying **1,866,908 points / 5.4375 CPU seconds** already consumed. Its
physical states, reference tolerances and acceptance limits are unchanged.

**Retained outcome:** producer `273af3e` passes creation and reconstruction
on both cases at both tolerances. The largest balance discrepancy is
`5.76473e-10`; the largest channel tolerance change is `3.56406e-10`. Radial
map error is negligible at the declared diagnostic resolution, while the
bulk field residuals and ordered interfaces recover the original imbalance.
See the [full interpretation and compact evidence](current-result.md#terminal-cap-transport-attribution-16-october-2026).
Creation plus reconstruction uses 23,516,078 points / 185.40625 CPU seconds.
With the failed attempt and prefix diagnosis carried, the total is 25,382,986 /
190.84375 within the original allowance. This diagnoses the fixed state; it
does not close the cap-wall evolution.
The standalone reconstruction passes with 11,758,039 points / 93.265625 CPU
seconds within its 24-million-point / 150-second allowance.

### Next analytical connection: cell fluxes and reconstructed transport

Use the existing finite-volume updates, native mass-component weights and
the proved transport identity. Derive the discrete integration-by-parts
relation between one interior face's central/upwind fluxes and the bulk
Pi/Phi residuals plus ordered jumps. Keep the field-value/potential term,
background time dependence, moving cells and physical endpoint terms. The
small finite field-value channel is about `3.26e-5`, not zero at the current
diagnostic resolution.

Exact symbolic and manufactured algebraic controls may use 20 active CPU
seconds, with no retained-field evaluations, reference solves, new trajectories
or grid scans. First establish whether the current flux/reconstruction admits
the required common mass balance and which algebraic terms prevent it.
Upwind dissipation and reconstruction defects must be distinguished before
selecting a corrected discretization. A scalar mass projection or deletion
of one channel does not provide the missing field evolution. This analytical
step supplies a decision for the next implementation, not another physical run.

**Analytical outcome:** the [mass adjoint and projection identity](nsc-atlas-finite-amplitude.md#variational-mass-response-and-flux-projection)
are proved for the exact primitive-conforming model. The reaction appears
as a nonlocal tail in the q-moment derivative; it converts the original
field equations into a common mass-gradient system. Native moment weights
need not equal conjugate face traces, and the difference explains why
upwind work alone has no guaranteed mass sign. The reference functional
must not be confused with the earlier diagnostic reference when a native
small-radius field approximation violates the exact primitive relation.

The zero-backreaction geometric energy calculation includes moving cell
widths, the wall primitive, incoming wall data and centre reflection. Its
complete semidiscrete identity has explicit nonnegative numerical losses.
Grok independently reviewed the adjoint, reaction and gradient derivations
and the two-cell sign counterexample. Fifteen exact controls pass in 39
development executions: **zero field samples / 2.984375 active CPU seconds**.
No retained physical state, radial ODE or trajectory was evaluated.

**Decision:** develop a discretization from the common mass gradient, with
the geometric energy result as its zero-backreaction control. First resolve
the boundary work and regular reconstruction, then background and moving
geometry compatibility. A central telescoping identity by itself is not a
consistent evolution, and the existing reaction cannot be counted twice.
The next compute admission must specify that construction and its controls
before any retained-field comparison or evolution. The previous physical
candidate, wall law and observed transport defect remain unchanged.

### Regular-centre geometric mass control

The next component is a fixed-radius, zero-backreaction control of the
derived mass-gradient route. Use continuous piecewise linear `chi,u` in
`z=atanh(r)`, with `chi(0)=u(0)=0`, exact element mass/stiffness integrals,
and the same incoming `w_plus=Pi+Phi` at the wall. The mathematical inputs
are the established energy/mass identity, polynomial integration, a symmetric
weak derivative, the trace inequality and a finite tridiagonal solve.

At a declared `0<R<1`, implement one instantaneous weak RHS and its mass
gradient, physical boundary flux, regular centre limits and boundary trace
discrepancy. Include the wall quadratic term in the stiffness; derive its
positivity and incoming-flux estimate. Distinguish the weak boundary flux
from the reconstructed last-element derivative. This component is a reference
control in the existing package, not a second cap evolution or a replacement
for the retained background/departure producer.

Manufactured controls may use at most **200,000 coefficient evaluations /
20 active CPU seconds**, with two to seven fixed nodes. Check exact linear
travelling-in-time data, independent element integrals and matrix gradients,
the mass/flux identity, centre regularity, nonzero boundary discrepancy,
small-radius arithmetic, immutable inputs and invalid/budget cases. No saved
field, radial ODE, time integrator, mesh scan or physical fit is admitted.
If the identity or regularity fails, diagnose it before another construction.
A pass supplies a regular reference and a specific finite-backreaction
extension question; it does not establish that extension or the evolving wall.

**Outcome:** the [regular reference](nsc-atlas-finite-amplitude.md#regular-centre-fixed-cap-mass-gradient-control)
has an exact positive quadratic mass, the required weak incoming flux and
an explicit completed-square input bound. Twenty-six controls pass in 97
development executions: **10,600 charged coefficient/manufactured values /
0.53125 active CPU seconds**. The reconstructed boundary derivative is
reported separately. Grok independently confirmed the weak equations,
incoming conversion, positivity and centre limits. The completed source
review confirms element indexing, gradient pairing, gap remainder and weak
flux handling. Its small-radius boundary concern is corrected using the
stable wall stiffness action, with a tiny-radius incoming-residual check;
the centre is also checked on kinked data. No physical trajectory is run.

### Fixed-cap finite-backreaction mass gradient

Extend the regular geometric reference to the same smooth conformal potential
and exact radial constraints at nonnegative alpha. The field representation
is continuous linear `chi,k` in `z=atanh(r)`, with centre values zero and
`k=A0 p`; at nonzero alpha, `k` is not the scalar time derivative. The exact
primitive supplies `f=chi/r`, `q=chi_z/A0`, and `p=k/A0`. Reuse the positive
geometric gap for the regular first-element formulas.

Differentiate the radial mass equation with respect to all free nodal
coefficients, holding radius and alpha fixed. Check those derivatives by
independent attenuation-weighted quadrature, and recover the existing exact
quadratic reference at alpha zero. Construct one weak RHS from the resulting
gradient and the original incoming wall condition. Report the difference
between weak and reconstructed momentum and normal derivatives at the wall.
An algebraic mass balance alone does not establish boundary consistency.

Before radial integration, bound the exact metric factor over the whole
declared interval using the linear field slope bounds and the positive
affine mass equation. Evaluate the sufficient bound with Arb at 160 bits,
using exact binary input values and an exact dyadic lower endpoint. A failed
sufficient bound stops this reference; it does not prove the geometry invalid.
This bound covers the exact radial solution for the declared fields, not
the numerical ODE error or a subsequent time evolution.

Manufactured controls share **2,000,000 charged coefficient/derivative values
/ 60 active CPU seconds**. Use two to seven fixed nodes, alpha zero and a
declared small nonzero alpha, quiet data and regular nonconstant fields.
Two fixed radial tolerances and independent 32-node per-element gradient
quadrature check the calculation; no retained state, mesh scan or physical
history is involved. Require zero-alpha mass/gradient recovery within `2e-9`,
nonzero-alpha forward/adjoint gradient agreement and tolerance change within
`2e-8`, and mass/weak-flux algebra within `1e-10` for these manufactured
scales. Include a directional mass check, exact centre limits, potential and
unit conventions, immutable inputs, rejected bound/domain cases, and shared
budget exhaustion. No physical trajectory or change to the coupled wall law
is admitted. A failure produces a specific diagnosis before more computation.

**Manufactured outcome:** all 19 controls pass in 71 development executions,
using **1,704,332 coefficient/derivative values / 3.859375 active CPU seconds**.
The finite-alpha gradient routes differ by `3.75699e-13`; the absolute weak mass-flux
residual is `7.80626e-18`. The exact-model metric certificate has lower value
`0.7837849351318801` for that prescribed case. Independent momentum and normal
wall discrepancies remain nonzero and are reported in the [owner](nsc-atlas-finite-amplitude.md#finite-backreaction-fixed-cap-variational-reference).
The component rejects a deliberately corrupted gradient before the weak RHS.
The completed mathematical review confirms the equations and bound, and
identifies a missing inherited centre-clock constant in the first prototype.
The lapse now starts at `-log cosh(sqrt(alpha/6) f_C)`; an independent lapse
integral and rate-normalization control would reject the former choice.
The mass functional itself is unchanged. Completed Grok source review found
no unresolved blocker in this fixed-radius reference. Its unread boundaries
were the unchanged metric helper after the lapse edit, the background-cap
owner and the shared-budget test; it did not execute tests.
All earlier prototype work is included in the cumulative receipt; no retained
field or physical time evolution has been evaluated.

### Exact Galerkin wall-momentum projection bound

Use the finite-cap variational derivative to identify `V` as the continuous
linear L2 projection of `g=Bp` in the geometric coordinate. Bound its wall
trace by the interpolation remainder and the strictly dominant rows of the
centre-reduced mass matrix. Differentiate `B/A0` with the same mass and lapse
constraints, retaining the inherited centre lapse. Use explicit field norm
bounds and Arb arithmetic to enclose the resulting spatial projection error.

This bound covers the exact Galerkin operator at the declared fixed state;
numerical radial/gradient errors, the normal-wall discrepancy and moving
geometry remain separate. It must vanish exactly at alpha zero and when the
momentum field vanishes. Reuse the sufficient positive-metric certificate,
and stop if its hypotheses or the resulting finite enclosures are unresolved.

Exact algebra and manufactured enclosure controls may use **100,000 charged
coefficient values / 20 active CPU seconds**, with no radial solve, retained
data, evolution or mesh scan. Check the local interpolation constant,
centre/interior/wall matrix margins, constrained derivatives of `B/A0`,
independent high-precision enclosure values, zero limits, invalid inputs and
budget stops. A bound too broad for a later physical calculation is reported
as such; it does not justify replacing the missing numerical error enclosure.

**Outcome:** thirteen controls pass. The exact manufactured momentum bound
is `0.0005850596031365239`, sharpened from `0.021594667452231934` by using
each affine element's own normal-field bound. No mesh or physical input was
changed. Independent Grok derivation confirms the projection constant and
constraint-based norm estimate. One symbolic test initially used an
indefinite logarithmic primitive with the wrong complex branch for the real
interval; checking the real primitive by its derivative and centre value
resolves that test without changing the bound.

Development has 27 test executions and 8,872 charged values. Recovered CPU
receipts total 1.53125 seconds; one intervening single-test output was not
retained across context transfer, so its CPU time is unavailable. The final
focused replay passes with 16 values / 0.25 seconds. Do not represent the
partial CPU sum as the full development cost. This calculation stops at
the exact spatial projection bound; radial numerical error, normal trace
compatibility and moving-radius transport remain separate work.

Completed Grok source review found no correctness blocker, including the
elementwise Phi sharpening and directed endpoints. It reused the previously
reviewed metric/input helpers and did not rederive the factor three; the
separate mathematical review covers that factor. The high-precision test
checks evaluation of the proved formula, not a sampled field maximum.

### Moving-cap geometric transport identity

Use the exact zero-alpha conforming mass functional to test the next proposed
connection: a uniformly moving geometric mesh with L2-projected coordinate
advection. Keep the original incoming boundary law and inherited clocks.
Derive the matrix shape derivative and compare it with the physical Reynolds
wall-density term, including the explicit geometric wall potential. This
tests a proposed discretization, not the physical nested configuration.

Exact symbolic controls have **20 active CPU seconds**, no field samples,
radial solves, mesh scans or evolution. Use a general element matrix identity,
a globally affine primitive and one fixed two-element kinked primitive.
Check both wall-speed signs and the static limit. If projection leaves a
nonzero transport remainder, retain it and redirect the construction; do
not subtract it by imposing the desired mass law. Grok independently derives
the same moving-domain comparison. Finite-alpha motion requires its own
shape derivative after this reference is understood.

**Outcome:** seven exact controls pass in ten executions with zero field
samples. The first suite measured 0.578125 test-body CPU seconds (shared
matrix-fixture construction was outside that timer); the three final
timelike-speed controls measured 0.265625 seconds with fixture construction
included. Their pytest wall times were 1.30 and 0.79 seconds. The timer now
includes shared algebra setup. The [derived remainder](nsc-atlas-finite-amplitude.md#moving-mesh-projection-remainder)
is nonzero: for the fixed two-element primitive of amplitude gamma and
timelike `Z_T=1/2`, it is exactly `3 gamma^2/7`. The incoming boundary data
can match the reconstructed derivative at that instant, so that discrepancy
does not explain this counterexample. The proposed projected-advection rule
is not adopted. No finite-alpha solve or physical trajectory was run.
Completed independent Grok derivation confirms the full remainder and exact
two-element coefficient. The next connection is the moving-domain mass shape
derivative and a slope-compatible transport rule, not another grid scan.

### Dilation-compatible regular cap space

The P1 counterexample redirects the moving-cap construction to a field space
closed under coordinate dilation. Test global odd polynomials in `xi=z/Z`,
with powers `1,3,...,2m-1`; oddness preserves the regular centre and the
Euler derivative `xi d_xi` stays in the same space. Keep the action,
potential, incoming boundary law and inherited clock unchanged.

Derive the exact moving-domain mass shape identity for arbitrary admissible
nonnegative alpha, then check its geometric reference limit by independent
symbolic differentiation of the quadratic mass. Exact matrix/algebra
controls have **20 active CPU seconds**, zero field samples and no radial
or time integration. Use at most three odd basis functions, general
coefficient symbols, both velocity signs and the fixed-wall limit. Check
the regular centre, Gram and stiffness matrices, full weak mass rate and
the distinction between the weak and reconstructed wall derivatives.

Grok independently checks the shape argument. Stop if that argument omits
a field, clock or boundary term. This calculation may establish a compatible
moving-domain reference space; it does not supply nonlinear numerical
gradients, their error bounds or strong wall transmission. The P1 projection
constant is not transferred to this different space.

**Exact-control outcome:** six controls pass in seven executions, using
zero field samples and **1.296875 active CPU seconds**, including algebra
fixture construction. An unexpanded symbolic comparison failed in the
first run; expansion proves the zero difference and the focused replay
passes. The [shape identity](nsc-atlas-finite-amplitude.md#dilation-compatible-regular-field-space)
supplies exactly the moving wall-density contribution without modifying the
mass or source law. It holds for the exact constrained mass and gradients;
the independent quadratic control checks its alpha-zero limit. Numerical
gradients and strong material-interface transmission remain open.
Completed independent Grok review confirms the shape identity, clock
normalization and gradient pairing. Before numerical implementation, the
global derivative bound and basis conditioning must be handled explicitly;
the old P1 helper cannot certify the new polynomial coefficients as if they
were nodal samples.

### Numerical regular-polynomial cap reference

Implement one moving-cap weak RHS in the proven odd-polynomial space.
Reuse the existing radial mass/sensitivity integration and independent
adjoint quadrature through a common engine; keep the original P1 wrapper
and its public results. A Grok worker owns that extraction while the root
owns the new field evaluator and its controls. No production evolution is
replaced.

Use one to three odd modes in the dimensionless geometric coordinate, with
the same potential and centre lapse. Certify whole-cap derivative norms
from the polynomial coefficients, then apply the existing radial metric
inequality. Check the Gram solve's conditioning and residual. Independently
integrate the radius-shape sensitivity on fixed `xi=z/Z`, including the
Jacobian derivative; compare it with the exact shape identity and a centred
radius perturbation. Do not set it equal to the desired mass law.

New manufactured controls share **2,000,000 coefficient/derivative values /
60 active CPU seconds**. One separate P1 regression may use **500,000 values /
15 seconds** after the extraction, and one projection-bound regression may
use **10,000 values / 10 seconds** for the reused metric helper. These are
new refactor checks, not a reopening of the earlier fixed-cap allowance.
Use the declared nonconstant three-mode state at `Z=.5`, `alpha=.002`,
`chi=(.03,.02,-.005)`, `k=(.01,.015,-.002)`, `T_x=.7` and timelike
`R_x=+/-.04`; include quiet/zero-alpha limits and fixed-wall recovery.
Two fixed radial tolerances, a joint coefficient/radius directional check,
an independent high-precision Gram solve and deliberately corrupted gradient
or shape values test the new connection. Gradient and shape agreement must
be within `2e-8`, manufactured mass-rate balance within `1e-10`, and exact
zero-alpha mass recovery within `2e-9`.

Stop on a failed domain, inaccurate solve, inconsistent derivative or
exhausted budget and identify its source before more computation. Main risks
are a missing shape Jacobian, importing P1 bounds into the polynomial space,
and confusing the weak boundary flux with strong material transmission.
Separate shape integration, certified global slopes and explicit momentum,
normal and scalar-tangent discrepancies address these risks. No retained
state, new waveform, mesh scan or time trajectory is admitted.

**Manufactured outcome:** all 26 monomial-reference controls pass in 54 executions, using
**1,252,570 charged values / 2.84375 active CPU seconds**. The separate extraction regression
uses **490,356 / 1.15625** for 19 P1 controls and **4,620 / 0.5** for 13
projection controls, within their separate new allowances. P1 values remain
unchanged. No existing budget is silently reset.

The polynomial shape identity residual is `1.37469e-15`, independent
gradient discrepancy `2.91746e-14`, and weak moving mass-balance residual
`7.05545e-17`. Normal and scalar-tangent discrepancies remain nonzero in
the manufactured state, as recorded in the [owner](nsc-atlas-finite-amplitude.md#independent-numerical-shape-derivative).
The initial scalar comparison used the weak characteristic pair. Following
the existing hybrid interface through its reconstructed outgoing field gives
the actual scalar discrepancy `-0.0560200` and interface mass-rate discrepancy
`-0.000880323`. Both trace errors are now retained in the exact decomposition;
the earlier weak-pair difference is named separately. Five affected controls
pass after that interpretation correction, without changing fields or rates.
One additional check constructs local exterior jets and calls the existing
hybrid boundary helper to verify the actual transmitted pair and rates.
Grok confirms the shape equation and global slope bound mathematically.
The delegated extraction terminated without changing its owned file; the
root completed the extraction and ran the regressions. Completed Grok source
review of `a36f16a` found no unresolved correctness blocker. The subsequent
26-control run adds 578,573 charged values / 1.328125 CPU seconds for the
review follow-ups and monomial regression. No physical coupled trajectory
has been run.

### Conditioned regular cap basis

The numerical three-mode reference verifies the equations but cannot yet
represent a resolved prepared pulse. Before increasing its field dimension,
replace the monomial coordinates with odd Legendre polynomials on `xi in
[0,1]`. This is a basis change in the same regular polynomial space; the
action, field, source, clocks and wall laws stay fixed. Derive the mass and
stiffness matrices, wall/centre evaluation and dilation operator together.

Evaluate `P_(2j+1)(xi)/xi` and its derivative in `xi^2` by a stable recurrence,
so the normal field retains its centre cancellation. Prove the global
derivative bound from the Legendre equation and derivative expansion, then
reuse the existing radial metric inequality with these actual bounds.
The old P1 projection constant remains inapplicable.

Basis-only exact algebra and arithmetic controls share **100,000 charged
values / 20 active CPU seconds**. Admit at most sixteen odd modes. Compare
the matrices against exact integration for three modes; check the centre,
wall, dilation and slope-bound identities; compare one fixed sixteen-mode
evaluation with high-precision polynomial values and derivatives at five
declared points. This is a basis arithmetic check, with no cap radial solve,
physical samples, mesh scan or time evolution. Stop on failed identities,
unresolved cancellation or lost precision before integrating the new basis
with the cap solver. Grok independently checks the derivation.

**Basis outcome:** twelve controls pass in thirteen executions with
**34,140 charged values / 1.390625 active CPU seconds**. The fixed sixteen-mode check differs from
85-digit values by at most `1.95399e-14` in the reported raw basis quantities.
Its diagonal Gram matrix has exact condition number 21. No radial solve is
included in those controls.

### Orthogonal-basis cap integration

Extend the existing polynomial adapter with an explicit Legendre option;
preserve its monomial default. Transform the derivative bound, centre clock,
Gram pairing, wall evaluation and dilation operator consistently. Reuse the
same radial sensitivity engine and shape equation. No second gravitational
implementation is introduced.

Orthogonal-adapter controls have **2,500,000 charged values / 60 active CPU
seconds**. Compare one fixed three-mode field in both bases at nonzero alpha;
then check one fixed sixteen-mode profile at zero and nonzero alpha, with
exact quadratic gradients for the former and independent forward/adjoint
and shape checks for the latter. These are manufactured operator checks,
not an evolving pulse or a mesh scan. Test known zero limits and reject
unsupported basis/count inputs. Existing monomial regressions use their
remaining numerical-reference allowance above, without resetting it.

The completed Grok source review of `a36f16a` found no unresolved correctness
blocker. Follow up its useful checks here: make `dZ/dx` explicit in the output
name, exercise the actual timelike `T_x B` threshold, check both Gram solves
against high precision, and force the residual-rejection branch. An
orthogonal basis removes monomial conditioning growth; it does not itself
bound field truncation or close material transmission.

**Adapter outcome:** five controls pass with **1,835,038 charged values /
0.6875 active CPU seconds**. The same three-mode field agrees in both bases.
The fixed sixteen-mode profile recovers the exact quadratic limit and
passes independent finite-alpha gradient and shape checks; the condition
number is 21. No physical trajectory or truncation claim is made. Grok
confirms the basis mathematics. Completed independent Grok source review
found no correctness defect in the extension. It read the basis controls
but did not execute tests or re-read the shared radial solver body.

### Fixed-profile projection diagnosis

The next missing connection is representation of the already supplied cap
background. Use the existing `CharacteristicProfile` without changing its
knots or source data. On one fixed slice, project its reflected primitive
and canonical density onto sixteen odd Legendre modes. Split at every
mapped knot: 19-point Gauss integration is algebraically exact for the
coefficient integrands, and 32-point Gauss for the squared residuals.
Compute the residual L2 norms, derivative seminorms, centre values and wall
first jets directly. Their finite arithmetic discrepancies are diagnostic;
they do not enclose interpolation error or finite-alpha evolution.

The mathematical check reuses the existing independent Grok reviewer.
Implementation controls have **2,000,000 charged values / 30 active CPU
seconds** for exact polynomial reproduction, an independently integrated
piecewise profile, centre/wall identities, quiet limits and domain/budget
rejection. The single input diagnosis has **4,000,000 charged values / 30
CPU seconds**, using the existing operator-comparison profile and the
zero-case terminal cut 128. Authenticate its existing fixed hashes; no old
trajectory is replayed. Store only compact norms, traces, hashes and the
command. No cap radial solve, new preparation, time integration or scan is
admitted. Stop on representation error that prevents use of this field
space; identify the missing resolution mechanism before further evolution.

**Projection outcome:** eight controls pass with 234,742 charged values /
0.90625 CPU seconds. The single authenticated input diagnosis uses
2,328,407 / 2.296875; no radial solve or time step is performed. Relative
primitive/density errors are 6.79%/14.60%, and derivative errors are
40.88%/56.45%; large wall errors and a spurious centre response exclude
using the current sixteen global modes for this narrow background pulse.
The [owning derivation](nsc-atlas-finite-amplitude.md#fixed-profile-representation-diagnosis)
records norms, trace errors and fixed hashes. Retaining the background
directly is the next candidate representation; its explicit time-dependent
mass balance must be derived before another physical evolution. The current
diagnosis is complete. Independent Grok source review found no correctness
blocker; its nonzero sixteen-mode exact-polynomial control was added and
passed within the same development allowance.

### Affine-background mass connection

The fixed pulse diagnosis requires retaining the supplied characteristic
background and representing its departure. Before implementing that RHS,
derive the exact explicit-time and moving-radius contributions to the same
mass functional. Check whether subtracting projected background rates from
the polynomial weak equations preserves the moving mass law. Admit one
independent Grok derivation and exact symbolic controls of the general
coefficient contraction and one polynomial background counterexample,
**10,000 charged algebra values / 10 CPU seconds**. No field samples,
radial integrations, source fits or physical evolution. If a remainder
survives, expose it and identify what additional compatibility is needed;
do not silently install the naive background-subtracted update.

**Affine outcome:** the [exact mass connection](nsc-atlas-finite-amplitude.md#mass-balance-with-a-retained-characteristic-background)
contains an explicit-time remainder. Three independently integrated
symbolic controls pass with 6,000 charged algebra values / 0.21875 CPU
seconds, zero field samples. A one-mode quartic-background example has
remainder `64 Z^6/5` for an incompatible imposed incoming value; the
compatible incoming limit has zero remainder and zero departure rates.
The moving shape identity survives because the background stays in physical
coordinates. No background-subtracted physical update is installed.
Independent Grok derivation confirms the full mass split, shape identity,
geometric-limit remainder and counterexample.

### Retained-background jets and local metric domain

To evaluate the derived explicit-time remainder on the narrow supplied
profile, retain its original Hermite input arrays with the existing
`CharacteristicProfile`. Add a bounded adapter for stable `f,p,h` and their
T/z derivatives, using exact polynomial moments near the centre and the
existing endpoint jets elsewhere. The action, source data and physical
preparation stay fixed. Exact-C2 Hermite mathematics refers to the original
binary input values; numerical interpolation-coefficient rounding remains
part of the numerical error, not a new physical datum.

Derive a local positive-metric certificate from the original profile knots.
Exact-binary Arb Bernstein coefficients bound the C2 quintic and its first
two derivatives. Combine these with the regular polynomial departure bounds,
then integrate the positive radial source upper bound interval by interval.
Bound `n/r` by the endpoints of each affine cumulative envelope. This
replaces an overly pessimistic global slope estimate; it is not a grid scan.

Jet and derivative controls have **300,000 charged values / 20 CPU seconds**;
certificate controls have **500,000 / 30 seconds**. Use independent
high-precision polynomial jets and exact radial envelopes, regular-centre
and quiet limits, inward finite-domain checks and failure before integration.
An independently reviewed fixed-profile certificate may use **100,000 /
10 seconds**, without a radial solve or time evolution. If the certificate
fails, report that sufficient-bound failure before proposing further work.

**Jet/control outcome:** fourteen jet controls pass in fifteen executions,
using 9,430 charged values / 0.34375 CPU seconds. The first exact-equality
assertion exposed ordinary coefficient-conversion rounding; the explicit
2e-14 check passes. Interior-knot tests cover both evaluation branches.
Seven local-certificate controls pass with 7,213 / 0.203125. The delegated
certificate implementation run ended without edits or a final report;
the parent implemented the independently derived recurrence. Mathematical
and source reviews are complete. The fixed-profile certificate passes with
78,388 charged values / 2.15625 CPU seconds, giving the outward metric lower
bound `0.6843674647287019` on 65 existing intervals. No radial solve or time
step is used for this certificate.

### Explicit-time retained-background radial reference

Extend the existing polynomial reference with a retained-background option,
reusing its radial mass/sensitivity engine, shape check and wall accounting.
Add the exact explicit-time source derived above and split at existing
profile knots. Project the background time jets with algebraically exact
19-point split Gauss integration. Report the explicit-time remainder and
its identity discrepancy separately; do not absorb it into a zero balance
residual or impose it as a new physical source.

Controls have **4,000,000 charged values / 60 CPU seconds**. Compare a
polynomial background against the same full field in the existing adapter,
including coefficient/shape/time chain rules, and check quiet and boundary
limits. Run the affected monomial controls within their remaining recorded
allowance; one sixteen-mode unchanged-path regression has a separate
**1,000,000 / 30-second** allowance because the adapter changes. After the
local metric certificate and source review pass, a separately admitted
fixed-profile calculation may evaluate the actual explicit-time remainder.
No physical time evolution is admitted by these operator controls.

**Reference-control outcome:** seven new controls pass with 378,959 charged
values / 1.0 CPU seconds. Independent Grok source review confirms the time
source, shape, sensitivity indexing, projection subtraction and separate
remainder. The 26 affected monomial controls pass with another 578,573 /
1.375 (cumulative 1,831,143 / 4.21875, eighty executions). The separately
admitted sixteen-mode regression passes with 837,921 / 0.25. The original
Legendre-adapter allowance is unchanged. No physical time step is run.

### Fixed retained-pulse time-remainder comparison

After the local certificate and background-adapter reviews pass, evaluate
one fixed retained background with zero departure coefficients and sixteen
modes. Use the existing zero-case terminal geometry, with at most one
explicit inward radius ulp to respect the exact finite profile endpoint.
Compare alpha zero against the finite alpha already recorded in the same
operator-comparison payload. Use the reconstructed incoming trace from that
same background, `T_x=1` and `R_x=0`. These are two instantaneous operator
references, not a reconstructed finite evolution or a new physical source.

The shared allowance is **20,000,000 charged values / 120 CPU seconds** for
one authentication pass and the two references. The geometric limit must
recover zero departure rates and its analytic time flux within `2e-7`;
that is a numerical acceptance threshold, not a uniform error enclosure.
Keep the existing `2e-8` forward/adjoint and independent shape checks. Require
the reported mass remainder identity discrepancy below `2e-10`. Report the
actual finite-alpha remainder and both wall trace discrepancies without
retuning them. Print compact values, input/producer hashes and execution
counts; no raw file, grid scan or time step. On failed domain, budget or
identity checks, report the consumed work and diagnose before another run.

**Fixed-input outcome:** both references complete within their shared
allowance: 15,264,426 charged values / 7.984375 CPU seconds. Alpha zero
recovers the retained pulse with maximum departure rate `5.50e-10` and
weak mass discrepancy `-1.39928e-10`. At the recorded finite alpha, the
weak discrepancy is `-0.02337491766289285`; the explicit-time remainder
accounts for it to `3.44e-15`. The actual hybrid mass discrepancy is
`-0.023386447174947225`, and the normal trace gap is `-0.00258401679`.
Forward/adjoint and independent shape checks pass in both cases. The
[scientific receipt](nsc-atlas-finite-amplitude.md#fixed-retained-pulse-mass-remainder-diagnosis)
retains values, input/producer hashes and the one-ulp domain adjustment.
No scan, time step or raw output file is generated. This diagnoses the
stated affine coupling; a compatible evolving construction remains open.

### Moving tangent-enrichment regularity check

The finite remainder motivates testing a space containing the background
time tangents `a=chi0_T` and `b=k0_T`. At a frozen slice, their exact inclusion
eliminates the projection remainder by the chain rule, with actual basis
values at the wall. Affine departure rates still subtract the represented
background rates. This does not establish a time-dependent finite space.

Before an enriched radial implementation, check the regularity of its
moving basis. Use the exact manufactured C2 profile
`F(u)=(u-kappa)_+^3`, `0<T<kappa`, `ell=kappa-T` inside `(0,Z)`.
For positive z its reflected primitive is `(z-ell)_+^3`; its second time
tangent is `b=6(z-ell)_+`. Adding `q*b` to the primitive introduces a moving
gradient jump. Derive the actual cap quadratic mass by piecewise integration,
compare its total derivative to ordinary pointwise differentiation, and
check the H1 difference quotient. A fixed rational three-function Gram
calculation can independently measure the missing third time tangent.

Allowance: **50,000 charged algebra values / 15 CPU seconds**, zero field
samples, radial solves or time steps. Retain the original profile and
existing numerical results unchanged. A nonzero moving-interface term
rules out simply feeding this enriched basis to the smooth explicit-time
source; derive interface conditions or choose a different field space
before a physical step. Independent Grok review checks the regularity and
interface sign. No grid or mode scan is admitted.

**Regularity-control outcome:** three controls pass across four executions,
using 24,000 charged algebra values / 0.765625 CPU seconds. One assertion
initially compared algebraically equivalent unsimplified expressions; its
focused corrected repeat passes. Exact integration gives the missing
moving-interface contribution `9 q^2`; the H1 gradient quotient diverges
as `36 q^2/epsilon`, and the frozen three-function L2 projection error of
`b_T` is `3/2`. The [derivation](nsc-atlas-finite-amplitude.md#tangent-enrichment-and-its-moving-knot-obstruction)
separates frozen tangent closure from a valid evolving space. This excludes
naive insertion of the moving enrichment into the smooth sensitivity
source. Independent Grok review confirms the full mass derivative, positive
interface sign, H1 obstruction and distinction from the unchanged C2
background. It corrects the earlier smooth-pairing-only enrichment proposal.

### Internal wavefront evolution and mass compatibility

Derive the missing internal-interface law from the existing scalar equations,
with continuous scalar, mass and lapse and no singular scalar source. At
`r_j(T)`, use the actual `B=A exp(delta)` and right-minus-left jumps. Check
the Rankine-Hugoniot residuals against the existing density/flux contribution
`v[rho]+T_x[J]`. This must distinguish satisfying the scalar equations from
merely cancelling one mass residual. A nonzero derivative jump at the
reference background speed need not be a solution at finite coupling.

Use exact algebra to check the general residual identity, both characteristic
orientations and a mass-cancelling but scalar-incompatible counterexample.
Then construct the exact geometric-limit evolution of the previous ramp
initial datum with zero initial velocity: it splits into two fronts. Verify
continuity, both weak scalar equations and the full cap mass, including
its wall term and moving-wall derivative, on the interval before either
front reaches the axis or wall. This is a manufactured geometric-limit
control, not a new source for the retained physical candidate.

Allowance: **50,000 charged algebra values / 20 CPU seconds**, no field
samples, radial integration, grid scan or numerical time stepping. A failed
identity is diagnosed before further computation. Independent Grok review
checks the nonlinear interface derivation and its finite-domain assumptions.
The actual finite-alpha field evolution and material-wall transmission
remain required before claiming a coupled physical configuration.

**Interface-algebra outcome:** five controls pass across six executions,
using 36,000 charged algebra values / 0.390625 CPU seconds. One initial
comparison was a structural equality of factored/expanded expressions;
its corrected difference vanishes on the focused repeat. The two-front
solution obeys both jump equations, the full geometric cap mass and its
moving-wall derivative. A separate counterexample proves that mass
cancellation alone does not enforce scalar compatibility. Independent Grok
review confirms the scalar jump and mass pairing, with the explicit
qualification that `B=A0` can occur exceptionally at finite alpha.

### Characteristic interface position and amplitude

The internal jump law must also advance its amplitude. Derive it from the
existing reduced scalar equations `p_T=(Bq)_r-S`, `q_T=(Bp)_r`, where
`p=r Pi`, `q=chi_r` and the source is continuous across a source-free front.
The radial constraints give a continuous `B_r` even when the derivative
fields jump. With orientation `sigma=+1` for the outgoing `p-q` front and
`-1` for the incoming `p+q` front, test the proposed laws
`r_T=sigma B`, `j_T=-sigma B_r j`.

Add a bounded stage component in the existing package using the same
normalized potential, constrained mass and inherited lapse. It supplies
position/amplitude rates and the physical jump polarization, not a separate
metric or an evolved cap. Use independent radial-constraint cancellation,
the exact geometric-limit characteristic and amplitude, orientation/clock
controls, quiet limits and invalid-domain rejection. Keep each front
strictly between axis and wall; collisions and endpoint interactions remain
separate boundary calculations.

Stage controls have **50,000 charged values / 15 CPU seconds**. No radial
solve, grid search or numerical trajectory is admitted. Independent Grok
mathematical and source reviews precede adopting this component. The field
solution must supply its metric trace; this interface ODE alone does not
replace the bulk scalar evolution or material-wall transmission.

**Amplitude-stage outcome:** five stage controls pass across six executions
with 404 charged values / 0.015625 CPU seconds. The initial run used 340
values with CPU time below timer resolution; the focused guard repeat used
64 / 0.015625 after scalar-real orientation validation was tightened.
Independent Grok mathematical review confirms both amplitude orientations,
the cancellation giving continuous `B_r`, the actual physical jump factors
and the exact geometric curves. There is no extra `B_T` in the rate, and
`j B` is not generally constant in a time-dependent metric. The common
bulk still has to maintain the complementary-family continuity. Source
review finds no implementation blocker. The suggested time-dependent
`j B` invariant in that review is excluded by the mathematical derivation;
the code does not use it. No coupled time trajectory has been run.

### Fixed physical-coordinate field space

The route review separates the original smooth pulse from the manufactured
wavefront problem. Its C2 input has zero derivative jumps; the source-free
front amplitude law preserves zero. Building more front machinery cannot
by itself fix the original bulk projection remainder.

Compare moving tangent enrichment, total fields on a fixed local polynomial
space in physical z, and a full nonlinear characteristic integral solver.
The next core candidate is the fixed local space: represent total chi and k,
not a moving background plus departures. Use the existing profile-derived
knots once, and a C1 quintic element space that can represent both initial
quintic chi and quartic C1 k exactly. Keep the wall within one cut cell for
the first finite patch and require a nonsingular active Gram matrix.
This avoids adding front state or a new physical preparation, and reuses
the constrained mass/gradient and wall laws. The direct characteristic
route remains an alternative if the local-space compatibility fails.

First verify the exact shape and weak mass identity for functions fixed in
physical z, including a non-dilation-closed test space. Here `N_Z=A0 rho_W`
has no coefficient-advection term. The Gram matrix depends on the clipped
domain, but field coefficient rates are not derivatives of the moments
`M d`; inserting their `M_Z d` term would be a different update. Check both
statements against an independently integrated cap mass. Also establish
exact C1 quintic reconstruction of arbitrary local quintics and the
quartic initial momentum. Respect axis parity and avoid interface/basis
birth events; strong wall traces and field error remain separate checks.

These exact controls have **60,000 charged algebra values / 20 CPU seconds**,
zero field samples or numerical evolution. Independent Grok review checks
the fixed-basis derivation. Only a successful identity and representation
check admits design of the numerical local-space adapter; this record does
not admit a grid scan or a coupled physical trajectory.

**Fixed-space outcome:** all four controls pass once with 20,000 charged
algebra values / 0.265625 CPU seconds. The independently integrated mass
obeys `N_Z=A0 rho_W` without dilation, the weak coefficient rates give the
same wall flux, and an incorrect moment-derivative insertion produces
`-Z_x k_W^2/2`. The Hermite-plus-two-bubbles form exactly spans each local
quintic, including quartic C1 initial momentum. Independent Grok review
confirms the algebra and local initial representation. Axis parity must
be enforced separately; a merely C2 profile with an axis exactly at a knot
requires its existing regularity qualification. The [next core route](nsc-atlas-finite-amplitude.md#fixed-physical-coordinate-bulk-space)
is a bounded local basis/initial transfer and mass-gradient adapter. No
coupled physical evolution is claimed or admitted by these controls.

### Local fixed-basis implementation and initial transfer

Implement the selected total-field space with at most 128 fixed cells.
Use C1 quintic endpoint Hermite functions plus two interior bubbles on
ordinary cells, and an odd quintic on the axis cell. The three axis-cell
degrees of freedom are its scaled axis derivative and the next node value
and scaled derivative. Shared node slopes use a fixed neighboring cell
scale. The ordering gives half-bandwidth five for the Gram matrix.

Retain exact rational knot coordinates and polynomial shapes as the
mathematical inputs. Map original binary Hermite labels at the selected T
using exact differences. Put the final basis endpoint one equal interval
beyond the wall, making the initial cut fraction one half. Form each initial
chi/k polynomial from the original exact-binary Hermite data, using only
source pieces that intersect the physical domain. The last polynomial
extension defines basis support, not new source data. Check exact C1
matching and axis parity before converting the initial degrees of freedom
to floating values. Report conversion error separately from later numerical
evaluation/ODE error. Preserve the existing source/profile implementations.

For numerical operators, assemble the degree-ten Gram integral with six
Gauss points per full or clipped cell, normalize its diagonal, and use
banded Cholesky with a finite condition check. Keep the last-cell cut in
`[1/4,3/4]`; reject rank/domain failures rather than remeshing implicitly.

Basis and transfer controls have **2,000,000 charged values / 30 CPU
seconds**. Check local endpoint identities, C1 continuity, axis parity,
independent exact polynomial integrals, Gram solves, clipped-domain shape,
exact source-polynomial transfer, quiet data, input immutability and domain
rejection. Grok independently checks the construction and implementation.
No retained-profile calculation or physical evolution is admitted by these
manufactured controls; any later fixed-input check receives its own budget.

**Initial control receipt:** all ten original implementation controls pass
with **76,132 charged values / 0.03125 active CPU seconds**. The exact Gram
comparison uses independent polynomial integration. Cell indices are zero
based: cell 0 is the special axis cell, and only cells 1 through N-1 have
bubbles. Thus N=1 has three degrees of freedom and N=2 has seven. Three
additional focused controls check that bookkeeping and the independent
Bernstein evaluator below. All three pass with another **49,308 charged
values / 0.125 active CPU seconds**: the combined thirteen-control receipt
is **125,440 values / 0.15625 active CPU seconds**, within the original budget.
Independent Grok source review finds no blocking formula error. Its wall-trace
and solve-residual suggestions add direct axis/wall comparisons and a banded
backward-residual guard. The two affected controls pass with another 71,810
values / 0.125 active CPU seconds, giving **197,250 values / 0.28125 active
CPU seconds across fifteen executions of thirteen controls**. A reviewer's
initial indexing objection was withdrawn after distinguishing cell and node
numbering; its bubble shorthand was corrected against the exact reconstruction.

### One retained fixed-basis transfer check

Authenticate the existing operator payload and use its already selected
terminal zero-case slice, including the previously declared one-ulp inward
wall coordinate. Transfer its original binary Hermite data exactly into
the local total-field space, then round coefficients once. Independently
construct degree-five Bernstein controls from the original endpoint jets
and evaluate them using exact de Casteljau arithmetic. On each physical
piece, six distinct rational coordinates establish polynomial equality for
chi, k and their first z derivatives. The cut-cell checks remain inside the
wall; its extended basis support does not authorize source extrapolation.
Compare the exact axis and physical-wall traces directly as well.

This single calculation has **3,000,000 charged values / 30 CPU seconds**:
one authentication, one transfer, one independent comparison and one
clipped Gram factor. Record coefficient-conversion bounds, observed
floating discrepancies and diagonal-scaled conditioning. The predetermined
numerical consistency tolerance for each channel is
`5e-12 * max(1, reference_peak)`; sampled floating agreement is not a uniform
error enclosure. Gram admission retains `condition * eps <= 1e-8`.
No radial solve, time evolution, resolution scan or source refit is admitted.
Grok reviews the implementation and independent comparison. Failure reports
its phase and consumed work before a new calculation is considered.

Success establishes the retained initial representation, not a finite-alpha
solution. Next design the fixed-space mass/gradient adapter and its metric
domain estimate; the old original-target bound does not automatically cover
rounded coefficients or future states. Preserve compact evidence and hashes.

**Outcome:** one retained check passes with **936,018 charged values /
2.484375 CPU seconds**, one transfer/factor and zero radial solves or time
steps. There are 65 cells and 259 coefficients per field. All 260 polynomial
identities are exact; direct axis/wall comparisons also pass. The largest
sampled floating discrepancy is `1.2505552149377763e-12` in k_z, within the
declared tolerance. The normalized Gram condition `30121278.38782701` passes
its current-point limit but requires continued numerical accuracy checks.
The [compact receipt](../results/summaries/nsc-atlas-cap-fixed-transfer-v1.json)
retains all input/source hashes and conversion bounds. Independent Grok
review finds no sign or domain blocker for this calculation. The result
admits design of the metric-domain and mass-gradient adapter, not evolution.

### Fixed local-field metric-domain certificate

Bound the exact rational local polynomials whose degrees of freedom are
the supplied binary floating coefficients. Use Bernstein control hulls on
each existing physical cell, clipping only the final cell at the current
wall. The axis cell is odd, so derivative bounds control chi/r and k/A0
without division by zero. On later cells bound chi, k and chi_z locally.
Use the existing nonnegative potential and cumulative radial source
envelope for `n_r+a*n=b`: `0<=n<=U`, hence `A>=A0-2*alpha*U/r`.
The calculation certifies these exact fields; it excludes floating field
evaluation, radial ODE error and future evolved states.

Manufactured controls have **500,000 charged values / 20 active CPU
seconds**. Independently check exact power-to-Bernstein identities, clipped
domains, the geometric-limit metric floor, an analytically integrated
linear-field mass below the envelope, input immutability and invalid-domain
rejection. Independent Grok review examines the certificate mathematics
and the later adapter design. Failure of a sufficient bound leaves the
physical metric sign undecided. No retained-input certificate, radial solve,
new partition, grid scan or time evolution is admitted by these controls.

**Initial controls:** all five pass with **17,600 charged values**; active
CPU time is below the reported timer granularity (0.0 seconds). The exact
Bernstein expansion is checked independently; the alpha=0 linear-field
mass follows from its integrated energy identity. Source review and the
single retained-input certificate remain pending.

After source review, one retained application has **1,000,000 charged
values / 15 CPU seconds**. Authenticate the same saved terminal slice,
regenerate its exact transfer once (full coefficient arrays are not saved),
and certify the resulting binary coefficients at the recorded finite alpha.
There is one certificate, no Gram factor, comparison replay, radial solve
or time step. Bind input/source hashes in a compact receipt. A failed
positivity bound reports the phase/work and stops; it does not trigger a
subdivision scan. Success permits the numerical mass-adapter design to use
this declared initial state; numerical evaluation errors remain separate.

**Retained outcome:** independent Grok source review finds no blocker.
The single application passes with **278,272 charged values / 2.296875 CPU
seconds**, including one regenerated transfer. Its 65-cell certificate
gives `A>=0.6870371535642128` and `n(R)<=21.971840781627296` for the exact
rational-basis fields with binary coefficient inputs. The certificate itself
uses 83,136 values. The [compact receipt](../results/summaries/nsc-atlas-cap-fixed-domain-v1.json)
binds sources and the retained selection. There is no radial solve or time
evolution. Numerical field evaluation and ODE errors still need control in
the mass adapter.

### Fixed local-field mass-gradient adapter

Reuse the existing radial forward-sensitivity and independent adjoint
quadrature engine for the total fixed-space chi/k coefficients. Integrate
each existing physical cell through an element coordinate with Jacobian
`h*A0`; clip only the final cell. Use stable odd-axis reduced-field formulas,
the new exact-field metric certificate and the banded Gram factor. Keep
the actual wall traces distinct from the mass-gradient-projected wall value.
The shape contribution is the boundary density `N_Z=A0*rho_W`; there is no
affine-background time derivative, dilation or moment-derivative insertion.

Manufactured adapter controls have **5,000,000 charged values / 60 active
CPU seconds**. Compare alpha=0 mass and gradients with independently
integrated local polynomial energy, check the corresponding coefficient
rates, and test a finite-alpha coefficient directional derivative and a
fixed-coefficient wall-shape derivative. Also check quiet data, the common
centre lapse, moving-wall weak mass identity, explicit strong/weak trace
differences, input preservation and domain rejection. Independent Grok
mathematical/source review supplies a separate check. These controls admit
no retained-input mass solve, scan, trajectory or publication claim; such
a calculation requires its own budget after the adapter checks succeed.

The first five controls use **1,832,284 charged values / 1.6875 active CPU
seconds**. Four pass. The alpha=0 independently integrated mass and gradients
agree, but the coefficient-rate check fails by `3.40552204e-9`: the clipped
Gram inverse amplifies the radial gradient error at requested tolerances
`rtol=2e-11, atol=2e-13`. Keep the rate assertion unchanged. One focused
repeat tightens the radial tolerances to `2e-13, 2e-15`, within the original
control budget. Propagated forward/adjoint gradient error also needs review
before any retained calculation. This is a diagnosed accuracy correction,
not a resolution scan or a new physical input.

The focused zero-coupling repeat passes unchanged with **376,974 values /
0.296875 CPU seconds**. Three affected controls at the tighter defaults use
another **2,251,384 values / 2.015625 CPU seconds**: quiet/domain controls
pass, and finite coefficient/shape mass differences pass, but the new
propagated canonical-coefficient rate difference exceeds `2e-8`. The total
is **4,460,642 values / 4.0 active CPU seconds**. Before changing accuracy
policy, use at most the remaining **539,358 values / 56 CPU seconds** for
one manufactured base state and its independently integrated Gram norm.
Print the gradient, coefficient-rate, physical L2 and wall discrepancies.
No repeated directional/shape pair or retained-input solve is admitted.

The single diagnostic uses **380,646 values / 0.328125 CPU seconds**.
The raw gradient difference is `7.70217e-16`, propagated canonical rate
coefficient difference `7.44041e-8`, and coefficient rate peak `93205.5378`.
Its independently integrated physical L2 difference is `1.01493e-12`;
the wall-velocity difference is `2.31570e-15`. Thus a fixed absolute bound
on representation-dependent coefficient infinity norms is not the right
field accuracy criterion. Keep those values reported, and compare physical
L2 and wall discrepancies instead. The original zero-coupling rate identity
and independent mass/shape assertions stay unchanged.

The development total is now **4,841,288 values / 4.328125 CPU seconds**.
A specifically admitted **1,000,000-value / 15-active-CPU-second follow-up**
checks one stable Gram norm against independently integrated energy and one
finite base-state propagated comparison. Compute `||L^T S v||^2` from the
scaled Gram factor `M=S L L^T S`, avoiding cancellation in `v^T M v`.
Both physical L2 rate discrepancies and the wall-velocity discrepancy must
be below `2e-10`; this is still a numerical consistency criterion, not a
uniform error bound. Do not repeat the finite-difference pairs or scan a
resolution/tolerance family. Source review must cover the new norm method.

The new norm control passes with **385,454 values / 0.3125 CPU seconds**.
It agrees with independently integrated energy and gives rate L2 gaps
`2.11376e-16` (chi) and `1.01493e-12` (k), with wall-velocity gap
`2.31570e-15`. Aggregate development work is **5,226,742 values / 4.640625
active CPU seconds**, including the diagnosed failures and one standalone
base-state diagnostic. No additional tolerance repeat is needed.
Independent Grok source and mathematical reviews confirm the element map,
delta preprocessing, weak rates, projected differences and stable Gram norm.
The accidental test ordering of coefficient-infinity and L2 values was
removed: there is no general ordering between those norms. Physical and
wall accuracy assertions remain unchanged. Coordinate rounding remains a
numerical error question; review prose does not supply its enclosure.

### One retained fixed-space mass reference

After adapter and norm review, evaluate one instantaneous reference using
the same authenticated terminal profile, exact initial transfer and recorded
finite alpha. Keep `T_x=1`, `R_x=0` and the previously specified incoming
value `(p_background+h_background)/R` from the same supplied source; this is
an initial trace comparison, not a new preparation or an exterior evolution.
The scope is one 65-cell, 518-sensitivity constrained mass solve, the existing
adjoint comparison, and the fixed-space weak coefficient update.

Budget: **25,000,000 charged values / 120 CPU seconds**, including
authentication, transfer, certificate, Gram and field work. Use the checked
defaults `rtol=2e-13, atol=2e-15`; no scan or automatic retry is admitted.
Require finite output, positive metric domain, weak mass-balance discrepancy
below `2e-8`, and forward/adjoint projected L2 rate and wall-velocity gaps
below `2e-8`. These are numerical consistency checks. Record the actual
strong/weak wall differences independently without forcing them to pass.
Print compact results, work, hashes and failure phase; preserve no raw
coefficient history. A failure diagnoses the specific condition before
another calculation is designed. Success precedes coupled wall evolution
and its numerical/error analysis; it is not that evolution.

**First outcome:** the run stops at **25,000,008 charged values / 7.515625
CPU seconds** in `_radial_variations` at its adjoint coefficient evaluation
(line 117), before completing the independent gradient comparison. The
forward loop has integrated every cell before that code is reached. No
retained mass/RHS is accepted from this interrupted reference.

The work diagnosis is explicit. For D=259, an ordinary adjoint sample
charges at most `140+160+(2D+48)+12+2D = 1396` values: basis evaluation,
reduced fields, coefficient vector, potential and adjoint accumulation.
All `65*32` adjoint points therefore cost at most **2,903,680** values.
The two wall evaluations, four banded solves, two norms and final scalar
checks add below 90,000. The completed forward prefix was already below
25,000,008, so replaying the unchanged reference needs below **27,993,688**
values under this accounting. After preserving the failed attempt, admit
one **30,000,000-value / 120-CPU-second** completion replay with identical
fields, mesh, source, alpha and radial tolerances. This changes the finite
work allowance, not the calculation. Preserve both receipts; no automatic
retry or further increase follows if that one replay fails.

**Completion outcome:** the failed attempt and its producer are preserved
in `e5d7cb55c35844bae20d0af73eee0aac85cdb5d6`. One unchanged-input replay
with the admitted larger allowance passes using **25,463,040 values /
7.734375 CPU seconds**. Aggregate retained work is **50,463,048 values /
15.25 CPU seconds** across those two attempts. The completed reference has
17,741 field sites and 518 sensitivity channels, at the original 65 cells.
Rescaled mass is `17.633757131085453`; wall lapse is `0.01342227558686664`.
The weak balance is `-8.68175e-15`, momentum gap `3.35819e-14`, normal
gap `-2.45602e-14`, and projected physical L2 gradient gaps are
`3.52426e-13` (chi) / `8.16598e-10` (k). All declared consistency checks
pass. The [compact result](../results/summaries/nsc-atlas-cap-fixed-reference-v1.json)
preserves scalar outputs and twelve source hashes, without coefficient
histories. It does not establish a material-wall evolution or numerical
error enclosure. The next connection is admissible initial junction data
and a common finite wall/exterior/field update using the existing candidate.

### Fixed-space cap and moving hybrid compatibility

Before replacing the cap inside the existing coupled stage, derive the
interface residuals of that replacement. The canonical primitive is
`chi=r*f`; the exterior invariant field is instead `c=r*f*sinhc(q_f)`.
Assume their values match at the instant under study. Let p,h denote the
actual reduced cap traces, P=V_W/B the mass-gradient momentum, and
Q=R*w_plus-P its weak normal trace. The old hybrid uses the actual outgoing
combination `(p-h)/R` and the exterior-determined incoming characteristic.
Its matched Pi,Phi need not equal either pair away from zero trace gaps.

Derive independently the transported common-field discrepancy and the
fixed-radius temporal mass discrepancy, retaining the actual wall boost,
frame factor and radial density. Also check what changes if a weak outgoing
combination is supplied instead. Do not identify weak mass closure with
material-wall closure or silently replace the actual normal derivative.

Exact symbolic controls have **40,000 charged algebra values / 20 active
CPU seconds**, zero field samples, radial solves or time steps. They check
the moving trace chain rule, exterior frame conversion, strong and matched
mass products, and a counterexample to sufficiency of weak balance alone.
Grok independently checks these residuals and maps the existing initial
junction/source/gauge requirements. The result selects the next coupled
connection; it does not authorize a prescribed-wall substitute or a scan.

**Exact outcome:** five controls pass with **36,000 charged algebra values /
0.296875 active CPU seconds**, without field samples or numerical evolution.
For `dp=p-P, dh=h-Q` and the actual outgoing characteristic, the canonical
scalar tangent discrepancy is `(-A*dp/H+H*dh)/(2R)`. The rescaled fixed-radius
mass discrepancy is `-A^2*exp(delta)*(P*dh+Q*dp+dp*dh)/2`. Sending a weak
outgoing pair instead leaves the moving primitive discrepancy `v*dh/R`.
The exact counterexample shows why weak balance by itself cannot prove
material-wall compatibility. Grok independently confirms the trace formulas
and identifies the mass normalization that must accompany them. The metric
identity gives `m=r^3/2+alpha*N`, hence the physical residual is `alpha*E_N`;
on a moving wall its vacuum density `3R^2/2` must also be subtracted. The
remaining 4,000 algebra values cover that normalization control only.

That sixth control passes with **4,000 values / 0.1875 active CPU seconds**.
All six controls together use **40,000 / 0.484375**, exhausting the original
algebra allowance without a repeated suite. The physical and rescaled mass
residuals now have an explicit common dictionary, including the moving
vacuum contribution.

### Fixed cap initial junction audit

The retained cap reference used a zero-case profile at finite amplitude,
with diagnostic `T_x=1,R_x=0`. Check its common value against the existing
stored terminal walls before trying a coupled step. Authenticate the
existing payload and its bindings, and bind the completed compact cap
receipt. Transform its canonical field by `c=R*f*sinhc(sqrt(alpha/6)*f)`;
compare radius, time, amplitude and scalar value separately. A discrepancy
is recorded, not removed by choosing a new exterior endpoint.

Use at most **10,000 charged values / 15 CPU seconds**, one authentication,
two terminal comparisons, and the existing quiet source/gauge corner as a
separate possible launch event. No field projection, radial solve, prefix
solve, time step or grid scan is included. If the terminal values differ,
do not call the hybrid on a fabricated common event. The quiet corner may
be checked algebraically against the vacuum junction using its stored
source/gauge row; this establishes initial compatibility only, not a
nonzero-time finite-amplitude solution. Save scalar differences and hashes.

**Outcome:** one audit passes with **2,180 charged values / 2.1875 CPU
seconds**. Neither terminal wall supplies the cap's common data. Relative
to the stored zero case, its invariant scalar differs by `-1.339906988e-5`
and the amplitude is different. Relative to the stored finite case,
radius differs by `7.671842997e-4`, canonical time by `1.280312330e-3`, and
invariant scalar by `-4.181379753e-5`. No terminal hybrid is evaluated.
The existing quiet corner at `x=.4,R=.45,T=0` passes its source, shared
value, incoming jets and gauge-rate checks; the largest residual is
`9.325873407e-15`. Its material rates are `R_x=1.084575413` and
`T_x=1.840737097`. The [compact audit](../results/summaries/nsc-atlas-cap-initial-junction-v1.json)
binds all inputs. The next finite window starts from this quiet event,
using one common field value and explicit scalar/mass residuals.

### One fixed-cap hybrid stage

Connect the checked fixed-cap operator to the existing hybrid through one
callback after the radial constraint solve. The callback receives immutable
scalar wall data, constructs the common exterior value, obtains the exterior
first jets and returns the hybrid's actual `T_x,R_x,w_plus`. These inputs
are mutually exclusive with diagnostic prescribed rates. There is one
radial solve and one exterior callback per stage; no circular rerun.

Keep the fixed-coordinate boost law and derive the physical lapse rate
from `K_plus`. Report the scalar tangent discrepancy, the rescaled temporal
mass discrepancies against both actual and matched traces, and the physical
mass factors including vacuum motion. Do not impose a second scalar value,
replace strong traces, or claim that weak balance makes the connection exact.

Implementation controls use one small fixed space, a quiet vacuum event and
manufactured nonzero coefficients/outgoing jets. A shared **1,500,000-value /
45-active-CPU-second** allowance covers callback/direct equivalence, common
value and derivative identities, vacuum recovery, immutable inputs and
invalid-port rejection. The existing alpha-zero energy control may run once
as an independent regression within that allowance. No retained nonzero
reference, preparation fit, time step or resolution scan is admitted here.
Any failure is diagnosed before changing an assertion or repeating work.

**Implementation outcome:** four controls pass once with **375,431 charged
values / 0.375 active CPU seconds**. The zero-coupling energy oracle is
included directly on the small fixed space; the older suite was not rerun.
The nonzero control differentiates the moving scalar trace independently
without another radial solve. It also checks the physical mass subtraction,
the strong outgoing channel and the residual that remains when it differs
from the weak mass pair. No assertions were loosened and no retry occurred.
The original cap receipt's twelve producer hashes are now authenticated
against their pinned `8f6c73c` Git bytes, so later adapter changes remain
distinct from the earlier retained calculation. A hash-only check passes;
the initial-event numerical audit was not repeated for this provenance edit.

### Fixed-cap and exterior common stage

Connect the instantaneous hybrid adapter to the existing exterior departure
prefix and its one shared transverse derivative. Preserve the source's
values/first jets and the evolved physical wall lapse; only the reference
wall label and boost are gauge inputs. The fixed cap stores total polynomial
coefficients, not the old finite-volume `Ip,Iq` integrals.

The exterior endpoint scalar is the cap value. Distinguish its actual
coefficient tangent from the hybrid's requested tangent, and report both
wall kinematic residuals and the independent numerical-minus-required
incoming first jets. Radius conversion roundoff is reported separately.
If the source supplies its own second radius jet, recompute the pure source
constraint from that jet; absence of curvature is an unavailable diagnostic,
not a zero residual. No mixed equation is reused as an independent check.

Reuse `vacuum_slice`, `departure_slice`, the departure-gradient operator
and the hybrid adapter. Frozen historical stage producers stay unchanged;
an extracted orchestration helper is checked against their actual output.
Focused controls have **2,000,000 charged values / 60 active CPU seconds**:
quiet vacuum recovery, unchanged exterior prefixes and source anchoring,
one derivative shared by prefix and coordinate transport, nonzero common
scalar/jet residual identities, and input/domain failures. No retained
physical trajectory, resolution scan or time step is included. The cap's
clipped Galerkin stability limit must not be replaced by the exterior
characteristic Courant number. A finite update needs its own justified
time-operator and error allowance after this connection is checked.

**Stage outcome:** four controls pass with **176,442 charged values /
0.203125 active CPU seconds**. The extracted exterior reconstruction and
prefix outputs agree bit for bit with the frozen stage. One focused repeat
after naming the requested and actual tangents distinctly passes with
**73,806 / 0.0625**; aggregate development work is **250,248 / 0.265625**.
The actual endpoint defect is the hybrid kinematic defect plus the cap's
scalar transport defect. All source/endpoint values stay fixed and the
numerical incoming jets are measured rather than replaced.

### Frozen cap time-operator bound

Before choosing an explicit time step, derive the actual linear generator
in the independently checked `alpha=0`, fixed-radius, zero-incoming limit.
The mass matrices use the same fixed basis. For `S=K-beta*e*e^T` and
`beta=A0/R`, the endpoint trace inequality supplies positivity of S.
Check dissipation in `H=diag(S,M)` and derive an energy-norm generator bound
from exact rational M,K,e, not from the scaled Gram condition alone.

At most **80,000 charged algebra values / 20 active CPU seconds** cover
exact small rational matrix controls of the energy identity, positivity
bound, norm bound and cubic Runge--Kutta polynomial. A rigorous scalar
remainder estimate may then bound one constant-operator step and its
repetition. No physical fields, source preparation, trajectory, step-size
scan or finite-amplitude stability claim is included. This determines the
missing time-operator requirement; moving geometry and nonlinear feedback
need their own bounds before a physical interval can be certified.

**Exact outcome:** four matrix controls pass once with **60,000 charged
algebra values / 6.5 active CPU seconds**. The constant cap generator is
dissipative in its actual energy metric; its squared weighted norm is
bounded by `2tr(M^-1 K)+(e^T M^-1 e)^2`. The one-cell odd-quintic control
at `Z=1/2` gives the exact rational value `3444` for that upper bound.
The controls verify the cubic step polynomial and an exact rational
energy-amplification inequality for one predetermined manufactured step.
There is no trajectory or finite-amplitude stability result. The coupled
stage's Grok source review also completes without an actionable defect.

Use the remaining algebra allowance to check the moving-cut energy term in
that same zero-coupling limit. Differentiate S and M exactly, keeping the
wall trace and its derivative, and compare the resulting energy rate with
the cap's moving-density identity. One additional control is limited to
**16,000 algebra values** within the original 80,000 / 20 total. The
result may bound cap energy and a continuous time residual for prescribed
geometry and incoming data; it cannot replace the exterior feedback or
finite-amplitude estimate.

**Moving-cut outcome:** the additional control passes once with **16,000
algebra values / 1.484375 active CPU seconds**, bringing this calculation to
**76,000 / 7.984375**. It verifies `S_Z=vv^T+2A0*ee^T`, `M_Z=ee^T`, the
source flux and the single moving-domain density. Their energy inequality
and continuous-residual estimate apply to the zero-coupling cap with the
same prescribed geometry and incoming data. They do not close the exterior
or bound a finite-amplitude trajectory.
Independent Grok review confirms these displayed identities and the
continuous-residual estimate under the same prescribed-input assumptions.

### Fixed-coefficient time update and continuous reconstruction

The common stage has a different state from the historical `Ip/Iq`
integrator. Reuse its exact-time validator, state-combination operation and
SSPRK(3,3) coefficients in a small adapter for the fixed cap coefficients.
Keep the authenticated historical producers intact. Supply one additional
RHS evaluation at the accepted endpoint, then reconstruct the continuous
cubic Hermite path and evaluate its derivative minus the same stage RHS.
An interior residual is a declared-point diagnostic; the previous linear
energy estimate requires a bound throughout the interval and identical
prescribed geometry/input. Finite coupling needs its own combined estimate.

Admit **3,000,000 charged values / 60 active CPU seconds** for known cubic
and exponential ODE controls, endpoint-slope and callback-failure checks,
and one manufactured nonzero common step with its midpoint residual.
The manufactured source is `c=(x-0.4)^3` with its analytic derivative,
`alpha=0.002`, four exterior cells, one odd-quintic cap cell and one
`h=1/1024` step. This is an integration control, not the retained physical
preparation, a resolution selection, a scan or a physical prediction.
Keep all stage arrays in memory. Record accounting and any failure before
designing a follow-up. Independent Grok review checks time weights, endpoint
slopes, the continuous derivative and the limits of the residual claim.

The physical-source constructor will separately reuse the original 129
Hermite knots, not interpolate the 1025-row evaluation cache. Before
advancing a quiet state to the source front, check the complete saved
source, including gravitational interpolation, rather than only its
analytically masked scalar. That input audit does not yet admit a retained
physical trajectory.

Admit one fixed-input quiet-front audit, with **50,000 charged values / 20
CPU seconds**: authenticate the existing 104 source bindings, reconstruct
only the original finite-amplitude 129-knot `SourceSnapshot`, inspect the
exact binary Hermite coefficients of the intervals up to the lower scalar
front, and query the last pre-front knot, the midpoint from it to the front,
and the front. Evaluate the reference wall at the knot and front. An exact
polynomial check supplies the independent radial-interpolation control.
Keep source curvature residuals and any nonzero pre-front gravitational
data visible. Save only a compact summary and source hashes. No source ODE,
profile fit, spatial mesh selection or physical trajectory is authorized
by this audit; a nonzero result redirects the starting-data choice.

The first audit stopped at the independent polynomial comparison, before
wall evaluation: a result-relative tolerance is unsuitable when the Hermite
polynomial cancels nearly to zero. Its precise CPU consumption was not
printed; the enforced upper bounds are **12,423 charged values / 20 CPU
seconds**. Retain that failed attempt in the final compact receipt. One
corrected comparison uses `128*eps` times the absolute polynomial-coefficient
scale for each derivative, with **37,577 values / 20 CPU seconds** remaining
for the completion. Thus the two-attempt ceilings are 50,000 values / 40 CPU
seconds. The source coefficients, knots and queries do not change.

**Time-control outcome:** four controls pass with 267,485 values / 0.28125
active CPU seconds. Grok source review finds no adapter bug and requests
derived-scalar and physical-residual coverage. One focused repeat adds the
moving-trace derivative/difference comparison, conformal scalar residual
and physical temporal mass identities; it passes with 263,520 / 0.25.
Total **531,005 / 0.53125**, with no retained physical trajectory.

**Source-audit outcome:** completion passes with **13,260 values / 2.34375
CPU seconds**, so total work including the failed-attempt upper bound is
at most 25,683 values / 22.34375 CPU seconds. The crossing Hermite interval
is nonzero before the analytically masked scalar starts. Its actual
constraint residual and the earlier quiet value/first-jet knot are recorded
in the [scientific owner](nsc-atlas-finite-amplitude.md#retained-source-at-the-first-scalar-front).
Grok review confirms the polynomial and interpretation; transverse source
splines remain outside this audit.

### First retained-source fixed-cap interval

Use the original source interval 35, from `0.42228630694544744` to
`0.42292305857246026`, split at its scalar front. Start with quiet cap
coefficients and departures at the checked earlier knot. Keep the original
129-knot source, its finite amplitude, full curvature residual and the
existing 64-cell quadratic exterior operator. This advances the retained
preparation rather than another manufactured waveform.

Derive a two-cell trial space from the reference endpoint clocks: for
`Z0=atanh(R0)`, `dZ=Z1-Z0`, `dT=T1-T0`, set
`a=max(4*dZ,2*dT)` and knots `[0,Z0-a,Z0+a]`. The initial cut is one half;
the reference final cut is at most five eighths. The last cell contains
the reference incoming layer because `a>=dT`. These are geometric scale
and domain choices, not spatial accuracy or finite-amplitude causality
certificates. All actual stage radii must remain in the existing cut domain.

First integrate exact rational basis products for M,K,e and evaluate
`L_*^2=2tr(M^-1 K)+(e^T M^-1 e)^2` at the initial cut. Two independent
small controls have **60,000 values / 15 active CPU seconds**. Use the
integer ceiling L and `T_max=ceil(max(reference T_x at endpoints))` to
choose each time-piece count `ceil(16*L*T_max*piece_length)`. Guard actual
positive `T_x<=T_max`. This preserves the checked initial frozen time scale;
it is not a stability theorem for the moving coupled system. No step-size
or spatial scan is admitted.

One numerical interval attempt has **600,000,000 charged values / 180 CPU
seconds**, at most 1024 accepted steps and five common stage calls per
step (including one interior continuous residual). Stop on an exceeded
budget, domain, clock guard, or failed independent identity; save a compact
diagnosis before considering another attempt. No retry is hidden. Keep
stage arrays in memory and save only aggregate residuals, physical endpoint
values, inputs, hashes and execution counts. Cap residuals use physical L2
norms, with the derived wall-scalar residual separate; source, scalar,
mass and wall defects remain explicit. This is the first retained-source
numerical trajectory of this fixed-space operator, not a continuum error
enclosure or a new observational prediction.

Reuse the existing independent pure-U constraint algebra on centre nodes
28..36 and wall-interior nodes 56..63 of that same 64-cell exterior, with
the quadratic radial curvature of the stored departures. This is a selected
spatial residual, not a full history derivative or a constraint enclosure.

**Attempt outcome:** the derived initial bound is `L=12106`, clock ceiling
3, and piece counts 236+135. The attempt stopped at its CPU check after
**341/371 accepted steps**, **1,708 stage calls**, **235,384,613 values /
180.015625 CPU seconds**. The last accepted label is
`0.42287166568847995`. The compact receipt retains the failure, accepted
history hash and maxima, including attempted intermediate stages. There
is no completed right endpoint and no retry is admitted here.

The observed scalar tangent gap is much larger than the sampled time
residual, while weak mass balance remains near roundoff. Diagnose the
finite-space boundary response before another evolution. For smooth
prescribed geometry in the zero-coupling limit, derive the onset orders
of the same Galerkin ODE under `C(t)=c_m*t^m+O(t^(m+1))`, with zero scalar
state. Admit **30,000 algebra values / 15 active CPU seconds** for two
exact controls: an independent cubic-drive Taylor recurrence and the
moving scalar-tangent conversion. No field or source replay, mode/step
scan, or weak-flux modification. The actual hybrid drive's leading power
and coefficient require their own derivation; this diagnosis must not
assume they equal the prescribed source profile.

**Algebra outcome:** both onset controls pass once with **24,000 values /
0.40625 active CPU seconds**. For a nonzero finite-order drive and a fixed
finite basis, the strong incoming trace is delayed by one power of time;
the displayed conversion to `E_f,E_c` is exact in the zero-coupling limit.
Grok review confirms those statements. This does not establish the actual
hybrid drive's coefficient. The next action is its one-sided expansion,
not another numerical trajectory or resolution scan.

### Actual hybrid onset and incoming cap reference

Connect the preceding finite-space onset result to the existing C3 source
and hybrid, rather than identifying the drive with the source by assumption.
Expand at the mathematical support root in the original source-radius
coordinate. Use the scalar prefix, fixed-space zero initial state and actual
hybrid outgoing port. Separate this asymptotic statement from binary source
masking and from a bound over the attempted finite interval.

Derive the early incoming alpha-zero cap directly from its wave equation,
including its wall characteristics and moving mass identity. Identify the
existing prepared-field value/first-derivative interface that can supply its
incoming history. This comparison must not prescribe a second independent
waveform or use a fitted terminal profile. Finite-alpha extension requires
its own characteristic and explicit-time mass terms.

Three new controls have **40,000 charged algebra/numeric values / 20 active
CPU seconds**: source and scalar-prefix onset, the actual hybrid leading
coefficient, and the incoming reference including its integrated mass and
production hybrid. No source solve, trajectory, resolution scan or replay
is included. Grok independently reviews the equations and maps the existing
prepared-field interface. An algebraic contradiction redirects the derivation
before any evolution is attempted.

**Outcome:** three controls pass once using **38,000 charged values /
1.28125 active CPU seconds**, with no source or field replay. The hybrid
drive has leading coefficient `2*A0*V_plus/H*64*kappa^4`, and the scalar
tangent gap is `-64*kappa^4*t^3+O(t^4)`. An integral inequality supplies
the exterior O(t^4) order under the regular finite-equation assumptions;
the source bump alone would not establish it. The alpha-zero incoming
cap closes the strong boundary relation and its quartic-layer mass
identity. These are local algebraic results, not a finite-window error
bound or an evolving finite-alpha construction. The [owner](nsc-atlas-finite-amplitude.md#actual-hybrid-onset-from-the-retained-source)
keeps the floating support edge and prepared-interface coordinate conversion
explicit. The next work derives compatible finite-alpha characteristic
and mass terms before implementing an incoming-history representation.

### Finite-coupling characteristic transport and wall port

Use the existing reduced scalar PDE and radial constraints to carry the
incoming field without projecting its boundary drive into the fixed
primitive space. In z=atanh(r), derive the conservative equations for
`E_plus=A0*(p+q)` and `E_minus=A0*(p-q)`, where here `q=chi_r`.
Retain their common finite-alpha source, the primitive compatibility
constraint, the regular-axis reflection and the original central clock.

Implement a bounded smooth-point operator in the existing package and a
shared-value wall port that calls the existing hybrid. The exterior callback
supplies its jets on that same field value. Report the radial mass density
and temporal mass target as targets; no regional mass derivative has yet
been calculated by this interface. Reuse the stable potential-reaction
quotient and analytic constrained metric derivative. Do not fit a source
history, change the finite candidate or replay the stopped interval.

Exact PDE/primitive identities and manufactured point/port controls share
**100,000 charged values / 30 active CPU seconds**. Check both families,
the alpha-zero incoming solution, an independent finite-alpha polar PDE,
clock/coordinate conversion, primitive inconsistency reporting, common
wall value and scalar tangent, metric/domain failures and caller budgets.
No source solve, radial ODE, spatial grid scan or time evolution is included.
Independent Grok review checks the characteristic equations and maps existing
history APIs. A failure redirects the implementation before a trajectory.

**Outcome:** all six controls pass in seven executions with **27,948 values /
0.234375 active CPU seconds**. The original five use 25,680 values;
the stable-frame control and affected port repeat add 2,268. Independent
70-digit finite-alpha polar arithmetic and 80-digit frame arithmetic agree.
The point operator and shared-value port reuse the potential, radial
constraint identities and hybrid, and report primitive discrepancy.
Their mass outputs are conditional targets, not a verified radial mass
derivative. Grok independently confirms the PDE signs and wall conversion;
its small-frame cancellation concern is addressed in the new port. The
existing hybrid is unchanged. The next connection is the retained incoming
history and its explicit-time constrained mass derivative, without a
source refit or scan.

### Retained incoming history in the fixed physical cap space

Extend `fixed_cap_reference()` with the existing `RetainedCapBackground`.
Keep the fixed C1 quintic basis for departures and retain the background
primitive and density directly. Reuse the radial forward/adjoint engine,
including an explicit-time sensitivity channel, and the established affine
mass remainder. The wall callback receives total fields. Preserve the
existing no-background path and the actual inherited central lapse.

The all-radius certificate combines the original Hermite control hulls with
the exact fixed-cell polynomial hulls on the union of their existing cuts.
Whole-cell bounds on subpanels are conservative; no subdivision search or
resolution campaign is admitted. Radial integration and background-time
projection split at existing fixed cells, mapped history knots and the
small-z formula switch. The C2 background leaves the radial source
continuous across its moving knots; no added tangent basis is introduced.

Controls share **6,000,000 charged values / 60 active CPU seconds**. Compare
an affine polynomial history with the identical total fixed polynomial,
check the explicit-time gradient independently, retain a nonrepresentable
incoming quartic front without projection, check the actual remainder and
wall trace, recover quiet/cancelling fields, and reject unresolved domain,
metric or budget conditions. Include the existing fixed-reference controls
once to cover the unchanged path. Independent Grok reviews the certificate
composition and mass derivative. No retained source solve or trajectory is
part of these controls; any failed identity is diagnosed before a rerun.

**Outcome:** all four new and six existing controls pass across 15
development executions, with **4,569,376 values / 5.28125 active CPU
seconds**. The first execution uses 629,879 new-control values plus
3,020,160 existing-control values (1.28125 + 2.4375 CPU seconds). Three
failures identify a one-ULP Runge--Kutta endpoint overrun. After restoring
only such endpoints, the three targeted controls use 441,108 / 0.796875;
one passes. The remaining checks diagnose coefficient amplification in
the clipped basis and a cancellation below mass resolution. Physical
Gram/wall norms and explicit unresolved-mass handling give the final two
passes with 478,229 / 0.765625. No equation, mesh, radial tolerance or source
is tuned. The [owner](nsc-atlas-finite-amplitude.md#incoming-history-with-fixed-departures-and-an-actual-mass-derivative)
records the calculated incoming mass, actual wall residuals and nonzero
finite-coupling remainder. Grok independently confirms the mass identity;
the metric bound covers exact input interpolants, excluding numerical
evaluation and evolution errors. No physical source replay occurs.

### Preserving the physical incoming onset during history transfer

Before importing the prepared profile, check whether its existing C2
quintic representation preserves the quartic onset established above.
The first three zero jets alone permit a cubic interpolation term. Derive
that coefficient, test an exact higher-order perturbation of the quartic,
and construct the unique first-interval polynomial that retains the known
quartic coefficient while matching the next value and first two jets.

This admits symbolic controls only, sharing **20 active CPU seconds and
zero field samples**. Check endpoint interpolation, onset derivatives and
uniqueness independently by solving the coefficient system. Keep the
prepared source's continuum-integral second jet distinct from numerical
interpolant derivatives. No prepared source solve, knot scan or trajectory
is included. A failure of onset order or coefficient redirects the adapter
before any physical calculation.

**Outcome:** three symbolic controls pass with zero field samples and
**0.21875 active CPU seconds**. Ordinary quintic interpolation can add a
cubic term to a quartic physical onset. A sextic removes that term but
does not in general retain the specified quartic coefficient. The unique
degree-seven first active piece preserves both while matching the next
three jets. The [owner](nsc-atlas-finite-amplitude.md#preserving-the-source-onset-when-retaining-its-history)
records the construction and its one-sided front derivative. Implementation
is supported by independent Grok algebra of the same endpoint systems and
source-to-arrival coefficient. It
must extend the same profile, centre moments and Bernstein input controls
together before the prepared history is imported. No physical data were
evaluated or changed.

The next implementation extends `CharacteristicProfile` with an optional
one-sided first active degree-seven piece and retains its source-fixed
quartic coefficient in `RetainedCapBackground`. All preceding supplied
jets must be exactly zero. Extend both centre moment routines by the
actual polynomial degree and route both metric certificates through the
same original-input piece controls. The existing quintic path remains
unchanged. Six Gauss points still suffice for the fixed-basis time
projection: its highest product degree becomes eleven.

Manufactured controls share **2,000,000 charged values / 30 active CPU
seconds**. Use independent degree-seven polynomials, one-sided front and
C2 joins, exact polynomial moment integrals across a quiet front, the
small-z switch, positive-radius constrained mass and time derivatives,
both original-input certificates, input-copy and invalid-domain checks.
Run affected existing profile/jet/bound controls once within this allowance.
Grok independently checks the Bernstein and moment formulas. No prepared
source solve, numerical evolution or scan is included. Unresolved front
coefficients, lost centre terms or inconsistent certificates stop this
component before the physical history is imported.

**Implementation outcome:** 73 distinct controls pass in 75 development
executions, with **787,030 charged values / 2.0625 active CPU seconds**.
The new controls use 769,126 values initially and seven on focused repair;
the existing profile, jet, bound and mass-partial suites use 861, 9,405,
7,213 and 418 respectively. Two assertions were corrected without changing
the numerical equations or tolerances: third jets need not match at a C2
join, and the underflow fixture must avoid triggering an earlier overflow.
Independent Grok algebra confirms the controls, derivative hulls and full
moment ranges. Exact quintic total transfer now rejects onset histories;
the retained-background route preserves them directly. No prepared source
or physical evolution was run.

The prepared-history adapter now has a concrete input route: use the union
of existing prepared-field cuts, arrival-spline knots and named source
edges, with at most 513 unique labels. Query `CapLeadingField.profile`
once with `second=True`; use its original value and first two jets, and
the source/arrival formula for the first active quartic coefficient.
Keep the mathematically quiet jets zero only under explicit numerical
quiet tolerances (zero by default), recording the original discrepancies.
Retain the binary front's offset from its exact input-derived support
root. Neither a quiet tolerance nor a profile coefficient fits observations.

Adapter controls use **100,000 charged values / 20 active CPU seconds**
on manufactured prepared-field interfaces, without a source solve. They
check native knots, source/arrival normalization, derivative provenance,
quiet-discrepancy rejection/recording, first-piece onset, finite domain,
input immutability and shared budget accounting. A physical replay and
its interpolation-error comparison still require their own fixed-input
calculation; the adapter alone supplies no such error bound.

**Adapter outcome:** five new and 28 existing wall-profile controls pass
once, with **8,596 charged values / 0.09375 CPU seconds** (467 / 0.015625
for the new adapter). They verify the actual API connection on manufactured
input, source normalization, derivative provenance, tolerance accounting,
copied inputs, finite-domain rejection and a shared prepared/cap budget
counter. No physical source was solved in these controls.

### One prepared physical history transfer, 17 October 2026

Run the existing fixed preparation once at its 257-field-node,
257-geometry-knot setting, using the existing order-four differentiated
wall integral for F''. Transfer the native knot union through the checked
adapter with **zero quiet tolerances**. Compare F,F',F'' at the midpoint of
each existing interval against direct calls to the same prepared-field
interface. This is a finite interpolation comparison, not a continuum
error enclosure or a resolution scan. Retain both absolute maxima and
their global peak-relative values, the onset coefficient, inverse bracket,
input hashes and work count. Do not save the full field or history arrays.

The complete source solve, transfer and comparisons share **10,000,000
charged values / 60 active CPU seconds**. Stop at the first domain,
quiet-jet, arithmetic or budget failure and save its phase and measured
reason. No automatic retry, altered source, fitted coefficient, cap time
trajectory or finite-amplitude claim is part of this calculation. Native
knot capacity, numerical quiet-front leakage and unresolved interpolation
error are the three specific risks. Its outcome selects the next actual
source/cap connection rather than another grid campaign.

**First outcome:** the source completes and the exact-zero transfer check
stops at `[0, 5.37477524314984e-40, 1.5572687669723015e-39]` for the quiet
value/first/second jets. The immutable [receipt](../results/summaries/nsc-atlas-prepared-history-v1.json)
records **1,052,512 charged values / 2.03125 CPU seconds**, one source
integration and no retry. No interpolation comparison was reached.

Three source-only endpoint arithmetic queries identify the sensitivity:
the nominal front gives zeta `-1.0000000000000016`, while the next V float
gives `-0.9999999999999956`. That neighbour has `g=6.22302e-57`,
`g'=7.38294e-41`, and `g''=6.56929e-25`. The nominal front differs from
the exact support root of the binary endpoints by `1.11022e-16`. These
are arithmetic evaluations, without another field solve.

The correction uses the known causal branch. Before the compact source
starts, `g=0`; the zero-initial Volterra equation for `z=chi-g` is
homogeneous with a bounded kernel. Uniqueness gives `z=0`, and the C3
source has zero jets through order three at the front. `CapLeadingField`
therefore returns that known zero through the declared source-front
arrival, rather than only through the earlier preparation corner. The
floating front and numerical arrival map remain distinct from their
continuum values. This does not erase a computed active-source residual.

Three targeted fake-prepared controls pass in 1.93 seconds of pytest wall
time, covering exact quiet/front jets, the still-active subsequent source,
and the independent zero-kernel second-jet chain. Those existing fixtures
did not print a separate aggregate charge counter; no physical source was
solved. The new physical attempt keeps **zero quiet tolerances**, the same
257/257 source, and the same 10-million-value / 60-second ceiling. It is
admitted only after this equation and arithmetic diagnosis; write a new
`v2` receipt and preserve `v1`. No source shape, amplitude or grid changes.

**Second outcome:** the [v2 receipt](../results/summaries/nsc-atlas-prepared-history-v2.json)
completes the transfer and all 392 native-interval midpoint comparisons.
It retains 393 knots, including onset index 87 at q=`-0.5049233867280902`,
with A=`429711.59455406887` and arrival derivative `1.0289023850141723`.
All three quiet discrepancies and tolerances are zero. Maximum absolute
F,F',F'' comparison errors are `5.96343e-12`, `1.82137e-10`, `9.60804e-5`;
global peak-relative errors are `5.99031e-12`, `7.69983e-12`, `8.03301e-8`.
The largest F'' discrepancy occurs at q=`-0.5040732796583627`, in the
first active interval. The existing inverse bracket is at most `4.44089e-16`.

The second run uses **1,362,327 values / 2.3125 CPU seconds**; the two
physical attempts total **2,414,839 values / 4.34375 CPU seconds**. There
are no cap time steps or raw histories. Both receipts' eleven producer
bindings were checked against their frozen Git commits without a replay;
array hashes require regeneration to verify their full numerical contents.
Independent Grok review confirms the quiet-branch argument. This known
branch does not measure possible leakage in a later numerical evolution;
such solver residuals must still be inspected separately. The next
connection uses the retained history with the actual evolving cap clock
and fixed-space departures, keeping its explicit-time mass derivative.

### Retained history in the common stage, 17 October 2026

Pass the retained background through the existing common stage and hybrid;
its time is the physical wall T, including the state's clock departure.
The hybrid still determines T_x, and the cap mass rate keeps T_x N_T.
Check a cubic history against the identical total polynomial represented
in the fixed basis, including nonzero clock departure. Independently
difference the constrained mass along the returned state/radius/clock
direction. Then take one manufactured common SSPRK3 step and evaluate
its continuous midpoint residual using the same history at every stage.

The two new controls share **8,000,000 values / 60 active CPU seconds**;
the existing stage and time controls retain their own unchanged limits.
No prepared-source solve, scan or physical trajectory is admitted here.
Stop and diagnose the first mismatch. A frozen reference clock, omitted
N_T, or double-added background would appear in the total-polynomial and
directional-derivative comparisons. Success admits the next physical
departure calculation, not a finite-interval error or stability claim.

**Outcome:** both new controls and eight existing stage/time controls
pass, in twelve development executions using **2,102,999 charged values /
3.28125 active CPU seconds**. The first comparison initially supplied
physical-z polynomial coefficients to a local-coordinate transfer; the
test now includes the exact cell scale. A subsequent raw extension
coefficient gap of `8.54453e-9` motivated comparison of the actual cap
field L2 norm and wall trace, each below the original `2e-9` allowance.
Neither the production equations nor the radial tolerances changed.

The independent mass-direction difference is `2.2976403795706864e-5`,
against the returned `2.2976377177617093e-5`, a gap of `2.66181e-11`.
The explicit clock contribution is `5.145624455085506e-6`. The single
manufactured step completes; sampled midpoint residual maxima for exterior,
wall and cap groups are `1.72367e-10`, `1.77314e-11`, `4.62579e-7`, in their
respective component units. These are diagnostics, not uniform errors.
Grok independently confirms the clock chain rule and single N_T term;
its review did not rederive the hybrid geometry or inspect the test bodies.

### Exact quiet-prefix integration cuts, 17 October 2026

The physical retained history has onset index 87. Its first 87 pieces
are identically zero by their original endpoint jets and the constructor
contract. Interior quiet labels 1..86 therefore need not split the radial
ODE or time-projection integrals. Keep label 0, the onset label, every
active label, fixed basis boundaries and the centre formula switch.
Original inputs, profile evaluation and domain checks remain intact.

Grok's independent review identifies a separate certificate issue: a
closed-interval hull can include the full active piece even at point
contact with its front. Coarsening those certificate intervals can weaken
the positivity bound. Keep all original certificate cuts and hulls.

One manufactured history with a zero prefix and degree-seven onset compares
the reduced integration partition with the original partition. Compare
both projections with exact rational polynomial integrals, and compare
mass, lapse, clock derivative and the remainder identity. Include one
point-contact hull check. Allow **3,000,000 values / 60 active CPU seconds**,
with no physical source solve or scan. A missing onset cut, changed
certificate, or rounding gap beyond the declared numerical tolerances
stops this connection for diagnosis.

**Outcome:** the single control passes once with **386,961 charged values /
0.609375 active CPU seconds**. Its radial partition goes from five to
four panels; the four-interval metric certificate is exactly unchanged.
Both time projections agree with independent rational polynomial integrals,
and the point-contact control retains the full onset hull. No prepared
source or coupled trajectory is rerun for this change.

### Source-front clock alignment, 17 October 2026

Before a physical coupled window, compare the retained history's declared
front label with `T+atanh(R)` from the existing reference wall gauge at the
original source-front x. Use one gauge query and the original order-64
source geometry under **10,000 values / 10 CPU seconds**. Compare q_V by
the independent chain `q_x/V_x` as well. This is a clock-interface check,
not a source solve, fit, scan or error enclosure. Any discrepancy remains
visible before choosing initial data for the numerical departure window.

The first query finds q differs by `1.0664592301894058`, the negative of
the prepared history's initial T. The gauge explicitly sets T=0 at the
quiet anchor, whereas the history has `Ta=-1.0664592301893958`. Pass this
inherited additive offset into the common stage; its derivative is zero.
The first handwritten V_x expression also had an erroneous factor four
from the source's normalized bump coordinate; use the actual derivative
of `V=(1-x/.8)*exp(x/.8)/U0`. That was a diagnostic expression, not solver
code. The query used 1,412 values / 0.015625 CPU seconds.

One further gauge query checks the corrected labels and chain rule. One
manufactured common-state control shifts all history labels and its clock
origin by the same constant, checking invariant fields, rates and mass.
This added control uses the existing 8-million-value / 60-second common
history allowance. No physical integration or source parameter changes.

**Outcome:** the corrected front labels differ by `1.0103029524088925e-14`;
the independent q_V is `1.0289023849503227` versus the retained
`1.0289023850141723`, a difference of `-6.384959228000753e-11`.
The second gauge query uses 1,412 values, with measured CPU below the
timer increment (reported 0); both queries total 2,824 values / 0.015625
CPU seconds. These floating discrepancies are recorded, not fitted away.

The clock-translation control passes with **272,778 values / 0.40625 CPU
seconds**. Including this control, the common-stage/time development total
is **2,375,777 values / 3.6875 CPU seconds** across thirteen executions of
eleven distinct controls. The common stage accepts the explicit constant
`background_clock_offset`; a nonzero offset without a history is rejected.

### Quiet source data before the retained-history window, 17 October 2026

Before skipping the old numerical quiet prefix, inspect the actual retained
source interpolant. Its scalar bump is analytic and masked to zero before
the front, but a metric Hermite interval can straddle that front. Load the
original authenticated 129-row source once and query exactly three declared
events: the old source-cell left edge, the penultimate quiet time node of
the old fixed schedule, and the front. Compare h and its derivatives,
scalar jets, source constraint and transverse jets against the quiet
vacuum values. Use **20,000 values / 30 CPU seconds**, with no field or
source ODE solve. This determines whether a shortened quiet initialization
is compatible with the numerical input; discrepancies are not snapped away.

**Outcome:** three queries use 135 charged values / 0.015625 CPU seconds.
The scalar value and first two derivatives are zero at all three events.
At the proposed shortened start, `h=1.4337521e-14`,
`h_x=6.2991614e-14`, `h_xx=-5.2705143e-7` and the source constraint is
`-6.8427801e-11`. At the front the constraint is `-6.9010356e-11`.
The transverse scalar interpolation is also nonzero, `3.2382790e-11` at
the shortened start. The common exterior RHS consumes the source values
and first x derivatives; these transverse snapshot jets and the second
derivative constraint are diagnostics. That distinction does not establish
equivalence to an evolved quiet prefix. Retain the original left endpoint
and its original steps in the next history-based comparison.

### Derived residuals with the evolving history clock, 17 October 2026

The wall scalar and constrained cap mass depend on the background clock as
well as the fixed departure coefficients and moving radius. Differentiate
those actual functions: include the background `k delta T` in `delta chi`,
the total radial trace, `N_T delta T`, and the cap vacuum mass derivative.
Check the contraction against two independently recomputed perturbed slices
in an arbitrary direction changing coefficients, clock and radius together,
once with the manufactured cubic history and once without it. Use the
existing common-stage control allowance, with at most **2 million additional
values / 60 active CPU seconds** for these two controls; no physical source
solve or trajectory. Require scalar and physical mass differences below
`2e-10`, normalized mass below `2e-8`, and clock contributions large enough
that omission fails. Stop at the first mismatch before physical execution.

**Outcome:** both controls pass with 630,439 values / 0.796875 CPU seconds.
The largest scalar, normalized-mass and physical-mass differences are
`3.14315e-13`, `4.24353e-15` and `5.68635e-13`. Four existing time controls
also pass with 270,183 values / 0.265625 CPU seconds. Four orchestration
controls use substitute producers to check success/failure accounting,
preserved clock origin and refusal to overwrite an existing receipt.

### One full fixed-interval evolution with retained history, 17 October 2026

Use the existing onset runner with its retained-history mode. The new
representation addresses the derived onset lag; it is not a refinement of
the earlier stopped polynomial-only attempt. Keep the original 129 source
rows, alpha, both endpoints, two-cell seven-dimensional cap departure
space, 64 exterior cells, frozen initial operator and 371-step schedule.
Retain the entire quiet prefix because its source interpolation is not
identically vacuum. No source or geometry parameter is fitted.

Regenerate the same 257/257 prepared leading field once, retain its 393
native history knots and quartic onset with zero quiet tolerances, and use
the inherited `Ta` at every common stage. This history is a representation
of the cap background; the physical finite-coupling departures and clock
are evolved. The original shared budget includes preparation, radial
constraints, forward/adjoint gradients, stages and midpoint residuals.
Allow **2 billion charged values / 900 process-CPU seconds**, with the
preparation subphase capped at 10 million / 60 seconds from the original
clock. The larger allowance includes the history jets and split radial
integration absent from the previous run. It does not authorize a retry.

Keep the independently checked scalar and mass chain rules, affine mass
remainder, all wall defects, source constraint, selected exterior pure-null
residuals, metric floor and actual time residuals. Preserve initial vacuum
recovery; record tiny pre-front scalar leakage arising from the declared
history/gauge front discrepancy rather than silently zeroing it. The
frozen alpha-zero timescale is not a coupled stability certificate.

The decision is whether the prepared history removes the leading scalar
boundary mismatch while preserving the actual constrained mass calculation,
and which residual then limits closure. Save a compact creation-only receipt,
source/input hashes and a streamed accepted-state hash; no raw history.
Freeze the producer before execution. A completed finite interval alone is
not a continuum or observational prediction. Stop at the first domain,
metric, gradient, algebra or budget failure and diagnose it before another
calculation. Likely failure modes are unresolved radial conditioning,
persistent cap/interface mismatch, and source or time residuals masking the
comparison; the retained diagnostics distinguish them without a grid scan.

**Outcome:** the interval completes all 371 steps / 1,855 common stages,
using 1,176,614,645 values / 619.4375 CPU seconds, including one source
preparation. The 18,513-byte [receipt](../results/summaries/nsc-atlas-fixed-history-onset-v1.json)
binds 70 producer files at `2551c9f`, the original input/geometry/schedule,
and the same prepared-history input hash as v2. All bindings and stage
accounting were checked without replay. The maximum conformal wall-flux
gap is `5.4098290e-8`, versus `7.0767723e-4` in the earlier 341-step
attempt; the new full interval includes that shorter domain. The same
derived scalar time residual is `5.3880648e-13`. This is a representation
comparison, not a convergence estimate or an observational prediction.
The affine weak-balance remainder is retained: weak balance reaches
`5.2841575e-12`, while its remainder identity is within `2.3940593e-20`.
The source constraint maximum is `3.5440928e-7` at the final source knot.
No automatic refinement or further evolution follows from this result.

### One-sided source curvature at the completed endpoint, 17 October 2026

`SourceSnapshot` is C1 in its metric correction and explicitly uses the
right Hermite interval for curvature at an interior knot. The completed
window ends at such a knot. Before interpreting its largest source
constraint as an interior defect, evaluate both polynomial second
derivatives at that same endpoint from the authenticated original source
rows. Use the same continuous value and first derivative, source profile,
amplitude and frame. Compare the actual two constraints with
`alpha*(h_xx-Q)/sqrt(F)` and their jump. At most **20,000 values / 30 CPU
seconds**, one endpoint with two analytic polynomial jets, no nearby-grid
query, source solve, field solve or evolution. The decision is whether the
reported maximum belongs to the interval's left limit or the next source
piece; neither side is substituted to make a constraint vanish. Stop after
this attribution or the first failed frame identity.

**Outcome:** the two actual curvatures at source knot 36 are
`h_xx(left)=-1.1354839067362777e-6` and
`h_xx(right)=0.0027239436776577454`, with common ODE target
`-5.8245777685680185e-6`. Their direct constraints are
`6.087902982207517e-10` and `3.5440928280303467e-7`; the latter reproduces
the saved endpoint exactly. Both agree with the curvature-defect formula
within `1e-15`. This uses 134 values and less than one process-clock tick
(reported CPU time 0.0 seconds). The large saved maximum belongs to the
next interpolation piece at the endpoint, not the left limit of the
completed interval. Both are retained; no constraint or RHS is overwritten.

### Scalar boundary channels before further evolution, 17 October 2026

Use the existing strong/weak factorization to split `E_f` into the cap's
weak-velocity discrepancy and its incoming-characteristic discrepancy.
Prove the split by independently differentiating the moving field and
subtracting the hybrid's matched pair. Include controls with each channel
alone nonzero. Add these diagnostics to the existing hybrid return without
changing its rates, source or wall law. At most **6,000 new symbolic values /
20 CPU seconds**, plus the existing quiet and nonzero hybrid controls under
their original 1.5-million-value / 45-second allowance. No source solve or
trajectory. A failure of the independent chain rule stops the change.

The completed receipt keeps aggregate and endpoint residuals, not the final
coefficient state, so it cannot numerically apportion these new channels
without regeneration. Future onset receipts will retain one bounded final
state, its hash, wall rates and endpoint diagnostics; this small reusable
input is distinct from a full time history. Preserve the completed v1 bytes.
This diagnostic change includes no reproduction or new evolution.

The first exact split passes (6,000 algebra values / 0.203125 CPU seconds).
The quiet and nonzero numerical controls pass in three executions using
201,657 values / 0.484375 CPU seconds. The first nonzero test incorrectly
assumed both channels were large; its velocity channel is only `4.1094e-13`.
That size assertion was removed, while the independently differentiated
identity and symbolic single-channel omission controls remain. No equation
or physical tolerance changed.

One further analytic bridge derives the velocity channel from the actual
mass sensitivity: differentiate the mass constraint at fixed chi, T and R,
apply its integrating factor, then convert dr to A0 dz. Check the variation
coefficients and an exact nonconstant polynomial projection with **8,000
additional algebra values / 20 CPU seconds**. This tests whether the ideal
velocity discrepancy is the endpoint projection error of `(B/A0-1)k`, and
vanishes at alpha zero in the declared lapse gauge. It does not estimate
the numerical gradient error or evaluate the completed physical state.

**Outcome:** the additional exact control passes with 8,000 algebra values /
0.1875 CPU seconds. The two new identities together use 14,000 / 0.390625.
The mass-sensitivity integrating factor gives the claimed projection formula
and its alpha-zero cancellation. No numerical attribution of the physical
endpoint follows without its coefficient state and gradient diagnostics.

Independent Grok review confirms the fixed-radius density variation and
the channel factors. Its unread lapse-stepper dependency was checked in
`_radial_variations`: the implemented radial RHS is
`jacobian*[rho,a,partials-a*sensitivity]`, with the declared centre lapse.
No separate lapse factor or moving-radius derivative was omitted.

### Physical-history endpoint velocity limit control, 17 October 2026

The exact alpha-zero identity supplies a cheap, discriminating next check
before any trajectory reproduction. Rebuild the unchanged 257/257 prepared
history once, require its input hash to match v2, and use the completed
receipt's exact two-cell basis, endpoint radius and stored background time.
Set cap departures to zero and alpha to zero for this explicit limit
control. This is not a reconstruction of the completed finite-alpha state.
Use one fixed radial-reference event with T_x=1, R_x=0 and incoming value
zero; these port choices do not enter its weak primitive velocity.

Compare the returned `v_W-B p_W` with the proved value zero and compare
the forward/adjoint projected wall velocities. Independently contract the
mass-gradient discrepancy against the background's Gauss-exact polynomial
time projection. This can reveal radial-gradient/projection error on the
actual prepared history without attributing the completed junction gap to
it. Preserve the source and original result. Use **10 million charged values /
60 process-CPU seconds**, including preparation, basis and the one radial
event, with one common clock and no retries. No grid scan, time evolution
or finite-alpha endpoint reconstruction. Freeze the code before execution;
stop after this one diagnostic or the first domain/budget failure. A
nonzero velocity is a numerical defect of this control, not an additional
physical force; its size determines whether the radial calculation needs
attention before any further physical evolution.

**Outcome:** the frozen `025abbb` control completes with 1,696,469 charged
values / 2.546875 CPU seconds, one preparation and no time steps. Its
velocity discrepancy is `2.10085855e-15`. Independent forward and adjoint
contractions are `2.10098798e-15` and `-3.23501807e-18`; the contraction
identity differs by `1.29427e-19`. The prepared input identity matches both
original receipts. This local limit check gives no evidence that a change
to radial mass-gradient assembly is needed.

### Leading endpoint transfer and common-port comparison, 17 October 2026

The completed physical trajectory lacks its final coefficients. Use an
explicit leading reference event to check the source-to-cap transfer before
considering another evolution. Rebuild the same 257/257 prepared field once;
authenticate the original source geometry and retained-history input hash.
At the original right endpoint compare the direct prepared wall value and
first derivative, the inherited arrival label/rate, and the retained
profile evaluated with the anchored reference clock. Convert V derivatives
with the analytic V(x) Jacobian. No stored waveform is differentiated.

Sample the same prepared field at the 63 original exterior interior nodes
using `WedgeField`; its retained-side guard must pass. Set the reference
metric, alpha zero and zero cap departures, and evaluate one unchanged
64-cell exterior/common-cap stage. Report both scalar channels, the direct
versus prefix derivative and all clock/profile differences. This defines
a new reference state; it is not a reconstruction of the actual finite-alpha
endpoint. Geometry is not adjusted to make the transfers agree.

Grok review confirms V_x, scalar normalization, the W/U conversion and
budget accounting. Its proposed clock objection is resolved by reading
`RetainedCapBackground.reduced_jets(T,z)`: the common stage supplies the
central time T and the background evaluates the two feet T+z and T-z.
The script explicitly checks the quiet reflected foot and equality of all
five original null-domain endpoints. Label/rate differences remain reported,
and partial-profile versus material-x derivatives are labelled separately.

The decision is whether a mismatch is already present in this one leading
event and where it enters. Use **10 million charged values / 60 process-CPU
seconds**, including preparation, queries, basis and the one stage, on one
shared clock. The preparation child is charged exactly once. Stop on the
first failed domain, input identity or budget check, or after this comparison.
No time steps, grid scan, source retuning or change to the evolution operator.
Possible failure modes are incompatible clock labels, an exterior query on
the excised side, and overinterpreting the reference event as the missing
physical state; the separate differences, retained-side guard and explicit
state definition keep those distinctions observable.

Before interpreting its cap rate, independently reduce the alpha-zero
mass density to wave energy plus its boundary derivative. A single exact
control, **6,000 algebra values / 20 CPU seconds**, checks that identity
and its integrated first variation. For a retained solution of the linear
wave equation with zero departures, it predicts
`d_x=-T_x*A0*R*d_plus*M^{-1}e`. This identifies a boundary-force comparison,
not a numerical claim about the previously evolved endpoint.

**Exact-control outcome:** the wave-energy density and first-variation
identities pass with 6,000 algebra values / 0.203125 CPU seconds. The local
comparison now reports the density-rate difference from the predicted
incoming boundary forcing, leaving numerical gradient/projection error visible.

**Numerical outcome:** the frozen `f154ecc` reference comparison completes
with 1,708,414 values / 2.5625 CPU seconds, one preparation and no time
steps. The history-minus-direct prepared material rate is `6.76460751e-9`;
direct-minus-exterior-prefix is `-3.22726177e-12`. The common conformal
gap is `6.76138468e-9`. Label and label-rate differences are `2.54241e-14`
and `2.14395e-11`. The canonical-density rate has L2 norm `6.52894299e-7`;
its difference from the proved incoming boundary force is `4.03972e-13`.
This identifies a local history-transfer contribution without assigning
the earlier evolved gap entirely to it.

### Local dense-polynomial derivative audit, 17 October 2026

The earlier midpoint comparison is sampled evidence, not a bound everywhere
inside a history interval. At this endpoint its retained F' differs from
the prepared interface F' by about `1.88e-9`. The interface value comes from
DOP853 dense output, whereas its first jet recomputes the field RHS. Before
changing either preparation or transfer, separate those two numerical
objects analytically at the same label.

The existing endpoint checker gains an explicit `--audit-dense-piece`
mode. It keeps the same source, clock, geometry and stage, reconstructs the
same input identity once, and retains only the one relevant dense polynomial
and one history piece. Differentiate the pinned SciPy DOP853 nested
polynomial by the product rule and compare against independently expanded
power-polynomial derivatives: **4,000 algebra values / 10 CPU seconds**.
Read the installed dense-output and segment-selection code; record the
SciPy version and its source-file hash with the audit.

Then measure `retained F' - dense F'` and `dense F' - RHS-derived F'`, with
the same analytic source derivative and arrival Jacobian. Their sum must
equal the original local transfer difference. This is an interpolation/ODE
defect diagnosis; the dense derivative does not replace the physical first
jet. The complete one-event audit retains the **10-million-value / 60-second**
shared budget. No grid scan, new trajectory, tolerance tuning or change to
the field source. Stop after this split, or on an unexpected dense type,
out-of-domain polynomial query, changed input identity or exhausted budget.

The independent power-polynomial control passes with 1,336 counted values;
active CPU time is below the reported clock tick (`0.0` seconds). Grok's
read-only review independently proposes the same three-jet split. Its
printed clock-rate difference was corrected from the stored operands to
`2.14395e-11`; its assumption that the queried piece is the first onset
piece is not used. The audit records the actual interval index and inputs.

**Audit outcome:** the frozen `ac3d84d` calculation uses 1,710,662 values /
2.46875 CPU seconds and no time steps. The actual piece is interval 88,
after the special onset interval 87. Retained-minus-dense F' is
`1.88328538e-9`; dense-minus-equation F' is `-3.06454081e-14`. The split
identity and dense-value reconstruction are exactly zero in this evaluation.
The original preparation and common-stage outputs reproduce v1. This points
to the ordinary quintic history interpolation, not the ODE's local dense
derivative, as the measured contribution at this point.

### Degree-seven interior-jet enrichment, 17 October 2026

For a quintic endpoint Hermite interpolant the degree-six remainder is a
multiple of `(q-a)^3(q-b)^3`. Its first derivative vanishes at the midpoint;
the previous midpoint-only comparison can miss its leading derivative error.
Check this identity and the following enrichment exactly, using no physical
field queries and at most **20 process-CPU seconds** for the two symbolic
controls.

For xi=(q-a)/h and psi=xi^3*(1-xi)^3, add
`psi*(A+B*(xi-eta))` to the ordinary quintic. At an interior label eta,
`A=DeltaF/psi(eta)` and
`B=(h*DeltaFprime-psi_prime(eta)*A)/psi(eta)` preserve both endpoint
value/first/second jets and match the supplied interior value/first jet.
At the midpoint this is `64*psi*(DeltaF+h*DeltaFprime*(xi-1/2))`.
It reproduces every polynomial through degree seven. Preserve the special
quartic-onset piece; no physical source, coupling or knot grid is changed.
Store the actual binary interior query label so original-input certificates
and runtime interpolation constrain the same mathematical point.

The implementation must carry optional interior jets through the profile,
retained input identity and both original-input metric certificates. Existing
degree-seven moments and six-point background projection can be reused;
the quintic total-field transfer must reject enriched input explicitly.
Keep the existing constructor's default behavior and historical evidence.
Numerical controls and the first enriched physical comparison require their
specific budget and decision after this algebra and consumer review.

**Exact outcome:** both new symbolic controls pass with zero physical field
samples / 0.125 CPU seconds. The degree-six derivative blind spot and the
general interior-label degree-seven reconstruction are verified independently
of the production interpolation code.

Implementation verification admits **2 million counted values / 60 active
CPU seconds** for the new independent runtime controls: rational degree-seven
reproduction at a non-midpoint interior label, unchanged quiet/onset pieces,
C2 joins, input copies and rejection, original-input Arb reconstruction,
third moments, one independent alpha-zero mass integral and both metric
certificates. A zero-endpoint interior bubble must be visible to each
certificate, and the quintic total-field transfer must reject it. The
prepared-adapter synthetic control uses its existing 100,000-value / 20-second
budget and checks one additional profile query without source integration.
Run affected existing controls under their unchanged per-module budgets.
Stop and diagnose any failure; do not use physical evolution as a unit test.

**Implementation controls:** six new and 105 affected existing controls pass
in 112 development executions, with 2,428,045 counted values / 4.015625
active CPU seconds. One new forward-velocity assertion was too strict for
the existing ODE absolute error (`6.04190e-13` versus `3e-13`). The corrected
control tests the independent weighted-adjoint contraction against the exact
zero identity, retaining the forward/adjoint discrepancy separately. No
production tolerance or operator changed. A mistyped test filename caused
one collection failure before any calculation. Forty source/storage controls
also pass after updating the two existing implementation-successor hashes;
all original receipts remain unchanged.

### First enriched leading endpoint comparison, 17 October 2026

Run the existing one-event checker once with `--enrich-history` and
`--audit-dense-piece`, after freezing its source. The original 257/257
prepared field, source geometry, 393 knots, onset coefficient, reference
endpoint, basis, exterior cells and zero-coupling reference state stay fixed.
Query one value/first pair at the binary midpoint of each ordinary active
history interval; these are representation data, not new source parameters.
Authenticate the original four nodal arrays against the old history hash
and record the enriched identity separately. Keep the old comparison files.

Use **10 million charged values / 60 process-CPU seconds**, including one
preparation, interior queries, certificate and one common stage. The local
decision limits, in the existing dimensionless x normalization, are
`1e-10` for history material-rate and conformal wall-gap differences,
`1e-11` for prepared/prefix rate and density-force L2 differences, and
`1e-12` for the scalar velocity channel. These are diagnostic decision
limits, not continuum error bounds. A miss stops with its measured values;
do not retune to pass. No grid scan, extra time trajectory or change to the
physical model is admitted by this comparison.

**Outcome:** the frozen `34c4422` enriched reference passes every predeclared
limit, with 1,866,626 values / 2.75 CPU seconds. The same prepared source
and 393 knots now include 304 ordinary-piece interior pairs. The original
nodal hash remains `a8dcfb845ab24cf1c97392d5156081517d4685587a4448f7a6a47ec8e8076066`;
the enriched hash is `11867d7614f2289025985d6f35000756002ef9c2f17c09a88466cc11a2517e04`.
History material-rate difference is `3.25911e-13`, conformal gap
`-4.70422e-13`, prepared/prefix gap `-8.02150e-13` and density-force L2
difference `5.81967e-13`. Direct prepared values, clocks and exterior sample
hash are unchanged. Grok's static review found no functional inconsistency
in the enrichment or original-input certificates. This is still one leading
reference event, not the finite-alpha trajectory.

### Enriched finite-coupling onset interval, 17 October 2026

The local transfer correction now has exact controls and a successful
physical-source endpoint comparison. Run the original finite-coupling
onset interval once with that checked enriched history: unchanged original
129 source knots, alpha, initial quiet state, clock offset, two-cell basis,
64 exterior cells, x start/front/end and the original 236+135 time steps.
Do not skip the quiet prefix. Authenticate the committed enriched comparison
and require both its old nodal and new enriched input hashes before evolution.

The existing **two-billion-value / 900 process-CPU-second** shared budget
includes one preparation, certificates, 1,855 stages and final diagnostics.
The preparation child retains its **10-million-value / 60-second** limit
from the original run start. No retry, source retuning, time/grid scan or
new amplitude is admitted. Retain both scalar channels, forward/adjoint
wall velocity discrepancy, source constraints and explicit-time mass
remainder. Save one small final coefficient/exterior state and its hash,
never a raw trajectory.

The decision is whether the original interval completes and its largest
sampled conformal scalar gap is at most `1e-10`. Completion and this
scalar decision are distinct fields; neither proves the remaining
constraints or a continuum error estimate. A failed decision or budget
stop requires its channel/cost diagnosis before another calculation.
First run the bounded orchestration tests with fake preparation and
transfer: verify enrichment forwarding, original-input receipt and single
charging on success/failure, plus early rejection of an invalid mode.

**Outcome:** frozen `42c7741` completes all 371 steps / 1,855 stages with
1,244,922,621 charged values / 660.703125 CPU seconds. One preparation,
no retries, and the original input/geometry/schedule are retained. The
maximum conformal scalar gap `2.7542985253647956e-13` passes `1e-10`;
the physical temporal-mass residual is `7.752513817634473e-21`.
Incoming and weak-velocity channel maxima are `5.692138312796339e-13`
and `2.0064215331527202e-14`. The selected exterior constraint remains
`2.0281018462806577e-11`, and the source curvature diagnostics are unchanged.
Grok review found no material orchestration bug; parent inspection also
confirmed the four SSPRK/endpoint abscissae and independent midpoint call.
All 70 producer bindings and the saved final-state hash authenticate.
The 28,312-byte result SHA256 is
`f6b5995960483d9285a4226e2c7ae1be4e6bd08ad53fdef2a30855edae7be62f`.
This supports the transfer diagnosis and ends this scalar-interpolation
comparison. Continuous source/exterior error bounds and the coupled
finite-domain estimate remain the next work.

### Whole-cell source representation bound, 17 October 2026

After the enriched scalar comparison, the source interpolation error remains
an input to the coupled validity estimate. `exterior_departure_stage` and
`departure_slice` consume `r, log_C, c, r_x, log_C_x, c_x`; the saved
transverse splines are not inputs to that evolution. Bound these six source
quantities on the original first active cell before changing any source or
running another trajectory.

Use the actual binary coefficients of `SourceSnapshot`'s cubic, its original
cell endpoints and compact pulse. Its left h and h' are exactly zero. The
borrowed mathematical tools are the integral remainder for a second-order
ODE, interval arithmetic and the mean-value theorem. Enclose the cubic's
own residual against `h''=Q(x,h,h')`; do not replace its second derivative
with Q. Split only at the already declared source front, enclosing the
analytic compact support and any tiny front-mask difference explicitly.
With cell length l and derivative bounds Lh,Lw, prove the common positive
radius tube and require `theta=Lh*l^2/2+Lw*l<1`. Then the radial errors are
bounded by `l^2*D/[2*(1-theta)]` and `l*D/(1-theta)`, where D encloses the
whole-cell residual. Convert these errors to the six used frame quantities.

The certificate concerns the exact polynomials with binary inputs and
analytic source formulas. Floating evaluation and the coupled exterior
evolution still need their own error account. Use the current operators and
original source; no new source integration, stage, grid, amplitude or time
evolution is admitted. Develop the helper and independent algebra, exact
manufactured ODE and original-input enclosure controls within **150,000
charged values / 15 active process-CPU seconds**. After those pass and the
source is frozen, allow one **100,000-value / 30 process-CPU-second** physical
cell certificate. Save finite outward bounds, input/source hashes and the
two panel summaries. If the tube or theta cannot be certified, retain that
diagnosis before designing another calculation.

**Development:** six controls pass in eleven executions, plus one manufactured
arithmetic diagnostic: 120,144 recorded values / 0.40625 active CPU seconds.
The first four controls passed; two interval-box controls exposed positive
denominators widened through zero by repeated ball powers. Monotone endpoint
powers preserve the known positive floors. The two repaired controls and
independent high-precision frame derivative control then passed. Charging
was moved before the interval calculations so failures are counted too;
two affected controls passed again. A conservative further 8,192-value
allowance covers the initially uncharged failing frame work, giving a
128,336-value total ceiling. No physical source, original-cell certificate
or time evolution was run in this stage.

**First physical certificate:** frozen `ac98a98` certifies the exact binary
cubic on its closed source cell, using 10,372 values / 2.4375 CPU seconds.
All 59 producer bindings authenticate. The v1 receipt is 13,697 bytes,
SHA256 `1f84dbf883857cc02465e1910692d21d3fb9581c7a69d1120acee6d59588e5d4`.
It bounds radius error by `7.86156e-16` and its first jet by `2.48453e-12`.
A subsequent exact-rational audit verifies its stored inverse inequalities
and identifies a right-endpoint distinction: the next PPoly owner's value
differs from this cubic's exact left limit by `-2.08603e-30`, and its first
jet by `-2.85265e-28`. That audit uses 1,024 rational-control values /
0.15625 CPU seconds. The v1 certificate remains valid for its declared cubic.

Extend the six-input bound to include the actual next-piece value/first jet
selected at the final knot. Add the exact rational jump magnitudes to the
two radial error bounds, and certify this endpoint in the common tube.
One manufactured endpoint-ownership control fits the remaining original
development allowance. After it passes and code is frozen, admit one v2
certificate under the same **100,000-value / 30-second** limit. No source,
grid, time interval or trajectory changes; preserve v1. The scalar front
mask remains a separate explicitly bounded input difference.

The endpoint-ownership control passes with 11,240 further values and active
CPU below timer resolution. The combined development ceiling is 139,576
values, including the earlier failed-work allowance. Grok's review of the
first certificate found no material soundness bug; the final-knot correction
is a separate follow-up to that review.

**Endpoint-inclusive outcome:** frozen `191a05b` v2 completes in 10,372
values / 2.5625 CPU seconds, with no integration or time step. All 59 source
bindings authenticate. The previous residual, source envelopes and theta
are identical; exact right-owner jumps are now included. The displayed
frame bounds remain `7.86156e-16` for r and `2.48453e-12` for r_x. The
14,751-byte receipt SHA256 is
`707bb8d76359edaf2a26228657a219933ebd164edcf89eb87742f9e66bec2bdd`.
Grok's follow-up finds the endpoint correction sound for all six inputs.
Floating evaluation and coupled propagation are still separate tasks; this
whole-cell representation certificate does not admit another evolution run.

### Connected scalar characteristic inverse, 17 October 2026

The source-cell certificate now supplies first-jet input errors. The missing
connection is their propagation through the exterior, two-way material wall
and cap. Derive a continuum comparison in those connected regions before
another evolution. Use the canonical reduced scalar on both sides, the
existing mixed equations, characteristic integration, the mean-value theorem
and a finite nonnegative comparison matrix. Keep the scalar-dependent wall
term: canonical scalar normal derivatives are not continuous at finite alpha.

Fix a positive-radius exterior and a cap collar with a quiet inner edge,
positive A, lapse and null rates, and a timelike expanding wall. Bound the
scalar inverse conditional on compatible geometry and coefficient bounds;
metric, clock and moving-domain differences must appear explicitly as forcing,
not disappear into a claim about the full coupled problem. The instantaneous
source/prefix/transmission dependencies must be closed before applying a
time-integrated estimate. A frozen finite-dimensional field matrix is not a
substitute for this continuum argument.

Admit one analytic derivation and independent Grok review, with exact symbolic,
rational comparison and manufactured transport controls within **100,000
charged values / 30 active CPU seconds**. No physical source readout, new
trajectory, grid, mode or amplitude scan is admitted. Save the equations and
small code/tests in their existing owners; no new numerical history or JSON
receipt is needed. Stop with the specific feedback or regularity obstruction
if the comparison cannot close. A successful scalar inverse leaves the metric
feedback bound and the uniform continuous residuals as subsequent work.

**Development outcome:** the connected continuum scalar inequalities have
an acyclic instantaneous graph, `B^5=0`; the exact positive inverse precedes
the short-time bound. Six controls pass in seven executions using 88,000
values / 0.59375 active CPU seconds. The first run passed five controls and
caught an omitted proper-time/source-time conversion in the manufactured
hybrid comparison; that test was corrected and passed. The equations and
production hybrid were unchanged. The canonical nonlinear normal jump is
nonzero in the finite-alpha control. A same-radius cap constraint estimate
exhibits the remaining metric coupling without new scalar derivatives.
No physical coefficient bound or full Einstein-wall inverse is asserted.

Independent Grok review confirms the wall algebra, scalar inverse and radial
variation, and identifies the required source-frame conversion. Add its full
Jacobian and conservative rational enclosure; this uses four of the existing
six error bounds and requires no second source derivative. One further
symbolic/75-digit comparison control passes with 11,024 values / 0.203125
active CPU seconds. The total is seven controls in eight executions,
**99,024 values / 0.796875 CPU seconds**, within the original allowance.

After freezing these sources, admit one canonical-source application within
**4,096 charged values / 5 process-CPU seconds**. Authenticate the existing
v2 whole-cell receipt and use its exact outward dyadics, original alpha and
proved radius/conformal/first-jet tube. Bound the common canonical source
segment algebraically: `r_min=F_min*rJ_min` since `sqrt(F_min)>=F_min`,
`u_max=g_max/r_min`, `cx_max=gx_max` and
`rx_max=vJ_max/F_min+alpha*g_max*gx_max/(6*r_min)`. The source mask only
removes the corresponding scalar values/jets, so these amplitude bounds
cover both representations. Convert its four source errors with the checked
Jacobian bound and save one compact hash-bound summary. No original field
readout, solver, trajectory, subdivision or scan is needed. Stop if the
authenticated tube does not justify these inequalities.

The first command at `9f0ec9c` stops during source authentication, before
constructing the scientific budget or evaluating the conversion: the input
JSON was incorrectly listed in the implementation-only source graph. Keep
the receipt's separate exact SHA256 authentication and remove it from that
code-only list. No physical result file or scientific evaluation was produced;
freeze the corrected runner before the one admitted conversion.

**Application outcome:** frozen `d365a13` authenticates the original v2 input
and all 13 producer bindings, then certifies canonical source value error
`<=4.658e-38` and first-jet error `<=2.576e-33`. Its 1,024 charged values
finish below process-CPU timer resolution; no source readout, integration or
time step occurs. A separate 128-value rational audit verifies those outward
decimal displays. The 6,437-byte summary SHA256 is
`6817ecccf46d2eb0620ea081d1ffd1cc67a4d319431a9b1e9cdd80b41dfe3b9a`.
Grok's final review confirms the complete source differential, common-segment
tube and rational conversion. This supplies canonical forcing inputs; the
full coupled inverse and floating implementation error are still open.

### Metric feedback through integrated exterior jets, 17 October 2026

The scalar comparison leaves cap-mass feedback through the exterior metric.
Keep the same fixed null-label gauge W(x) for the exact and approximate
solutions. Use spatial L1 errors of the exterior transverse jets, whose
moving-domain derivative contains the wall input multiplied by W_x and
integrated in time. A boundary change at a single endpoint is not a change
throughout the exterior. Bound the radius, canonical-scalar and lapse prefixes
explicitly from those integrated errors and the existing source inputs.

Combine these prefix bounds with the same-radius cap constraints and the
fixed-domain cap characteristics. Determine the instantaneous dependency
graph before imposing a time-smallness condition. Preserve the nonlinear
canonical wall jump, cap lapse anchor, metric/clock variations, and all
boundary and primitive residuals. The estimate is conditional on a positive
common coefficient tube and bounded reference jets; uniform physical bounds
and a continuous numerical residual still need their own certificate.

Admit analytic derivation and independent Grok review, plus exact algebra and
manufactured moving-interval controls within **30,000 charged values / 15
active CPU seconds**. No source readout, evolution, mesh change or scan is
admitted. Retain the conditional comparison, its required norms and any
specific unclosed feedback; do not claim a physical validity interval before
the coefficient tube, residual and continuation margins have been checked.

**Development:** the explicit exterior prefix and moving-integral bounds
give a mixed-norm graph with `B^2=0`; cap metric and clock feedback enter the
time block. Four exact/manufactured controls pass with 29,000 values / 0.0625
active CPU seconds. The scalar lapse term is independently checked through
70-digit canonical conversion against the production bulk equation and is
`-alpha*f_U*f_x`, without an extra C factor. Initial transverse jets,
source-frame conversion, primitive/constraint residuals, reference derivative
bounds and separate pointwise tube checks remain explicit requirements.
Independent Grok review confirms the actual prefix inequalities, cap metric
identity and nilpotent coupled graph. A suggested extra C in the lapse term
was rejected by the source equations and the independent conversion control;
the reviewer confirms the canonical convention. The common tube, compatible
initialization, uniform reference jets and all residual enclosures remain
the inputs required to evaluate this conditional bound on the original cell.

### Uniform exact reference gauge, 17 October 2026

The coupled mixed-norm estimate now needs uniform physical coefficients.
First enclose its existing analytic reference gauge on the original first
source cell, using the exact corner `a=4/5,x0=2/5,R0=9/20`. The numerical
gauge and original source remain unchanged. Their binary normalization,
quadrature, inversion and floating evaluation errors remain separate inputs.

Admit exactly two certified endpoint integrals, from `9/20` to the fixed
radii `474/1000` and `477/1000`, and interval evaluation of the algebraic
wall/clock and vacuum-wedge coefficients. Use 160-bit balls, requested
absolute/relative quadrature tolerance `2^-80`, degree limit 64, evaluation
limit 512 and depth limit 20 per integral, with one shared hard allowance
of **100,000 charged values / 20 process-CPU seconds**. Each integrand
callback charges 64 values; other interval work charges 2,048. No source
query, waveform reconstruction, trajectory, root search or scan is admitted.
Stop if the two endpoint inequalities, real-path regularity or positive
wall coefficients are not proved. Save the compact dyadic enclosures and
authenticated original cell/producer hashes, not quadrature histories.

Independent exact junction, normalization and high-precision rapidity
controls have a **20,000-value / 15-active-CPU-second** allowance. Three
controls pass with 16,000 values / 0.09375 active CPU seconds. Grok reviews
the coefficient formulas and monotonic bracket independently. The likely
failures are an overly wide radius enclosure, an uncertified square-root
branch in complex integration, or conflating this exact reference with the
floating gauge and sourced tube. Fixed bracket inequalities, the analytic
branch flag, and a separate initial/gauge error ledger address these cases.

**Outcome:** frozen producer `81b0170` certifies both fixed bracket
inequalities with 5,248 charged values and process CPU below timer resolution.
The exact reference wall stays strictly inside `(0.474,0.477)` throughout
the original cell, with positive wall/clock coefficients. The 13,822-byte
[summary](../results/summaries/nsc-atlas-reference-gauge-bound-v1.json)
contains 17 authenticated producer bindings. Independent Grok review and a
separate 256-value rational audit confirm the formulas, inequalities and
outward displays. No source query, numerical root search, time step or retry
was performed. Numerical gauge errors and the sourced common tube remain
separate work; this result does not assert the coupled error threshold.

### Compatible exact initial anchor, 17 October 2026

The reference enclosure leaves the numerical gauge/initial comparison open.
The existing state stores departures, so specify its continuous reconstruction
using the exact reference and the unchanged binary departure coefficients.
At the original initial event all departure coefficients vanish. Check the
analytic source is before its exact support front and the entire cap lies in
the retained history's exact quiet prefix. Then the reference vacuum exterior,
de Sitter cap and signed wall supply compatible initial data. The recorded
floating wall is compared with that construction, not declared identical.

Admit one fixed-seed inverse-residual integral and two clock integrals, at
the radius already saved in the enriched-onset receipt. Reuse the proved
whole-cell bracket for the mean-value inverse and clock derivative bounds.
Use 160-bit balls, `2^-80` requested integration tolerances, degree limit 64,
evaluation limit 512 and depth limit 20 per integral; one shared allowance is
**100,000 charged values / 20 process-CPU seconds**. Each callback charges
64 values and algebra/metadata work charges 3,072. Authenticate the previous
cell, reference, onset and source-geometry inputs. No new numerical gauge
root, preparation, trajectory or scan is admitted. Stop if the initial
correction leaves the bracket, either quiet margin fails, or the corrected
wall leaves the original fixed-basis collar.

Independent manufactured inverse, exact junction/rate, quiet-profile and
moving-coordinate controls have **20,000 values / 15 active CPU seconds**.
Grok reviews the initial compatibility and later residual interpretation.
The three likely failures are circular localization of the seed, a corrected
clock crossing the retained support edge, or reusing a floating stage defect
as the exact reconstruction's residual. The prior full-cell bracket, explicit
quiet-margin calculation and common-coordinate chain rule address these.
Future whole-cell residual bounds remain required. No assertion of uniform
libm accuracy or monotonicity of the floating quadrature is made.

**Outcome:** four controls pass with 15,004 values / 0.40625 active CPU
seconds. Independent Grok code review confirms the initial enclosure and
scope. Frozen producer `33445f7` then completes the three-integral application
with 7,872 values and CPU below timer resolution. The saved physical wall
radius and relative clock differ from the exact initial reference by at most
`1.125e-16` and `1.244e-16`, respectively. The entire cap remains quiet with
margin above `0.0014530219858557`; the source is strictly before its exact
front and the original fixed-basis cut remains valid. The 12,790-byte
[receipt](../results/summaries/nsc-atlas-reference-initial-bound-v1.json)
retains the numerical source-curvature defect. A separate 512-value rational
audit verifies the displays, strict inequalities, 43 producer bindings and
four input hashes. Future sourced residuals and continuation remain open.

### Actual source forcing for the coupled tube, 17 October 2026

Compare the prescribed source with the exact zero-field vacuum reference.
The needed inhomogeneous data are the actual source values and first jets,
not the previously certified interpolation differences. Reuse the existing
uniform analytic pulse bounds G0, G1 and G2 and the zero-left-value/first-jet
source constraint. Bound its radial correction by integrating a uniform
acceleration envelope, verify the original Jordan/frame continuation margins,
then convert to the six invariant source components, canonical X and Q, and
the radius-product input to the exterior primitive.

Admit exact rational evaluation only, with outward rounding to a fixed
128-bit dyadic denominator. The bounds cover `0<=alpha<=alpha_max`, where
alpha_max is the unchanged recorded coupling. One physical application has
a **4,096-value / 5-process-CPU-second** allowance and reads the authenticated
source-cell receipt. Independent symbolic frame identities, zero-source and
zero-coupling limits, manufactured 80-digit frame values and rejection of an
unclosed source tube have **30,000 values / 15 active CPU seconds**.
No source integration, new waveform, trajectory, root search or scan is
admitted. Stop if the source continuation margins fail. Save only rational
bounds, input/producer hashes and their outward display values.

The three likely failures are substituting representation error for physical
forcing, losing a metric/field cross term in the frame change, or assuming
the vacuum speed bounds also enclose the sourced solution. Actual boundary
differences, independent frame identities and a later sourced-tube check
address them. The inner quiet-core condition is the wall-to-core travel bound
`sup(beta*T_x)*duration < inf(Z_wall-z_star)`; the initial whole-cap history
margin cannot be used to require the wall itself to remain quiet. Only
`E_minus=0` enters from that core; a zero outgoing-to-core trace follows from
non-arrival. A closed a priori bound still needs a local construction and
continuation argument before it establishes a sourced spacetime.

**Outcome:** four controls pass with 21,192 values / 0.28125 active CPU
seconds, and Grok confirms the complete source/frame inequalities and
outward dyadic rounding. Frozen producer `7025f41` completes the physical
application with 2,048 values and CPU below timer resolution. The full
coupling-interval bounds include `abs(X)<2.018e-7`, `abs(Q)<0.003501`,
`abs(r-x)<7.404e-16` and `abs(L)<4.940e-18`. The 6,365-byte
[certificate](../results/summaries/nsc-atlas-source-forcing-bound-v1.json)
retains 19 producer hashes. A separate 512-value rational audit verifies
the displays, ceilings, input hash and continuation margins. No source
query, integration, root search, time step or retry was performed.

### One conditional coupled tube image, 17 October 2026

The exact initial event, reference gauge and actual source forcing are now
certified. Evaluate a single predeclared continuum bootstrap box on the same
cell and coupling interval. The fixed bounds are `J_r,J_L,abs(deltaT)<=1e-6`,
`J_X<=1e-7`, cap `abs(E_plus)<=0.004`, `abs(E_minus)<=1e-8`,
`abs(chi)<=1e-5`, pointwise exterior `abs(delta r_U)<=1e-5`,
`abs(delta L_U)<=0.001`, `abs(X_U)<=1e-7`, and `abs(r_x)<=1.01`.
The core edge is the original first positive cap knot. There is no search
over boxes, source parameters, grids or time intervals.

Use the existing mixed-norm equations, explicit finite differences of the
wall maps and the cap radial constraints. Keep the sourced wall radius in
the A0 upper bound and the complete `-2*beta_bar/R` curvature contribution.
Check both characteristic margins, cap A, wall H, expansion, the longitudinal
derivative bound and inward travel to the quiet core. The physical application
has **20,000 charged values / 10 process-CPU seconds**, with a fixed charge
of 8,192 for one 160-bit interval evaluation. Authenticate the three compact
certificates and original collar metadata. Save all image components, margins,
coefficient bounds and hashes, including an unclosed result. Stop before
designing any replacement if an image or positivity check fails.

Four independent potential, zero/finite-coupling wall and manufactured
envelope controls have **40,000 values / 20 active CPU seconds**. No waveform,
source integration, numerical root, trajectory or refinement is admitted.
Grok reviews the inequalities and implementation. The three likely failures
are omitting a wall coefficient, substituting an L1 norm for a pointwise jet,
or reading strict bootstrap improvement as existence. Direct wall controls,
separate suprema and an explicit conditional result address them. A local
construction and derivative/continuation estimate remain required even if
every component improves strictly.

The first control invocation's output was unavailable after context transfer;
the pytest cache recorded its nodes without a failure. One bounded repeat
to recover an authoritative result passes all four controls with 26,192
values and active CPU below timer resolution. Charging the unavailable
invocation at its entire allowance gives a conservative cumulative upper
limit of 66,192 values / 20 active CPU seconds. This does not repeat or
admit a physical application.

The derivative-continuation follow-up is analytic: differentiate the exterior
prefixes, cap transport, radial constraints and wall traces within the same
positive tube. A separate Grok reviewer checks the dependency order. Three
symbolic controls have a shared **30,000-value / 20 active CPU-second**
allowance (8,000 values each). They independently differentiate the canonical
quotients/kernel, commute the cap PDE, and differentiate its radial metric
constraint. No physical coefficient evaluation or trajectory is added.
An unresolved derivative cycle or use of an unproved incoming constraint
stops the proposed continuation argument before an existence claim.
All three symbolic controls pass with 24,000 values / 0.328125 active CPU
seconds. The transverse-jet system closes as a triangular linear ODE without
substituting an incoming null constraint.
Independent Grok review confirms the conditional closure and identifies the
remaining local construction and constraint-propagation proof. Its boundary
constraint terminology is read against the owner: wall-born exterior rays
need the incoming U constraint, while the outgoing x constraint comes from
the exact source. The review's stronger regularity exclusion is not used;
the proof needs the stated finite source derivatives, not a claim about the
maximum possible regularity of the resulting fields.

Grok's completed implementation review confirms the conservative wall,
prefix, cap and travel inequalities. Frozen producer `706fe56` evaluates
the original box once: all ten components improve strictly and all seven
margins are positive, with 8,192 charged values and CPU below timer
resolution. No source query, integration, time step or retry occurs.
One separate exact-rational receipt audit is admitted with **512 charged
values / 5 CPU seconds**, checking every enclosure and display, box ratios,
input identities, producer hashes and the stated scope. It performs no
physical calculation or replay.
The audit passes with 512 values / 0.09375 CPU seconds: 24 producer bindings,
four input identities and all 64 interval records authenticate. The
[29,878-byte receipt](../results/summaries/nsc-atlas-coupled-tube-bound-v1.json)
has SHA-256 `192caf213cd4b67961f815203687a87748a9668c4a97d2d77f9ee2f5d78a8dab`.
The largest image ratio is approximately `0.515020212`; the exact upper
endpoint is below `0.516`. The next mathematical requirement is a compatible
local construction, using the conditional derivative estimates above.

### Local construction through the existing null wall map, 17 October 2026

The next investigation keeps the same action, source and quiet physical
initial data. Use the existing two-double-null boundary operator to define
one frozen-forcing characteristic map, with wall lapse choices treated as
coordinate gauges and the central observer clock recovered afterward.
The immediate question is whether the map preserves field primitives and
common wall values before taking a nonlinear fixed point. A separate Grok
worker reviews this route while the original worker investigates compatible
iteration in the hybrid polar chart.

Three exact symbolic controls share **36,000 charged values / 30 active
CPU seconds**, charging 10,000 each. They check a fixed source strip with
wall-born points, a spacelike initial curve with wall-born points, and the
six incoming wall jets against their common tangent rates. They evaluate
manufactured polynomials and identities only. No physical source, new
trajectory, coefficient scan or approximate solution is admitted. A complete
local-existence claim still requires an invariant complete function space,
controlled moving-domain extension, contraction and corner compatibility;
the map identities alone do not establish those requirements.
All three exact controls pass with 30,000 values / 0.203125 active CPU seconds.
The [owner](nsc-atlas-finite-amplitude.md#a-primitive-preserving-local-construction-map)
now gives the outgoing prefix, wall ODE, incoming transport and primitive
reconstruction explicitly. A nonlinear fixed point and the compatible chart
conversions remain unproved.

A fourth symbolic control uses 4,000 of the remaining 6,000 admitted values.
It records the obstruction to treating a small uniform/L1 iterate error as
a small time-derivative error, even under a common Lipschitz bound. Only this
new control is run; the three completed controls are not repeated. The
separate Grok review supports the double-null route and its constraint/gauge
logic, while leaving a moving-domain contraction proof to be supplied.
The new control passes with 4,000 values / 0.046875 CPU seconds, for a
cumulative 34,000 values / 0.25 seconds. The owner's fixed-point constraint
implication and exterior relabeling are explicit. Claims of contraction
based only on unspecified small-time constants or on strong convergence of
uncontrolled derivative defects are not adopted.
The owner now specifies a coordinate-clamping extension that covers the
zero-width cap corner, the exact exponentially weighted prefix inequality,
and the scaled moving-wall norm. Invariant Lipschitz and contraction bounds
are explicitly retained as proof obligations with their required order of
constant selection. No additional numerical control or physical evaluation
is run for this analytic follow-up.

The next analytic follow-up defines M as a weighted derivative seminorm,
retains the actual nonzero source in a finite `D(lambda)`, and derives the
derivative and difference bounds with lambda chosen before M and h. Moving
entry times and output projections retain their displacement terms. A
larger cap strip supplies the initial data beyond the terminal outgoing
wall label; its unit weight and short characteristic lengths avoid a
corner extension factor. These derivations are submitted to Grok review.
No source evaluation, new trajectory, scan or numerical control is run.
The remaining connection is chart coverage and continuation on the original
physical cell, rather than a replacement of it by a smaller quiet interval.

Independent review confirms the weighted norm estimates on the small wall
wedge and the conditional native cap coordinate, scalar and lapse first-jet
bridge. A second differentiated Volterra/transport hierarchy now identifies
how the physical wall derivatives would be generated from the same outgoing
data; its review is pending. Three exact controls are admitted with
**24,000 charged values / 30 active CPU seconds** (8,000 each): nonlinear
chain rules and lower-curve derivatives, differentiated wall-entry data,
and the polar/null metric and scalar identities. They use a manufactured
smooth solution and symbolic chart variables only. No physical source,
coefficient bound, trajectory or new tube evaluation is admitted. Their
outcome checks those identities, not full-domain existence or chart coverage.
All three controls pass with 24,000 values / 1.90625 active CPU seconds.
No prior control or retained physical calculation is repeated.
The local-map review also confirms the expanded cap strip: both integration
lengths are at most rho*h, unit cap weight avoids a long-diagonal penalty,
and its clamp and entry estimates close without an inverse-width factor.
This resolves local domain construction; physical-slice coverage is reviewed
separately from it.
The second derivative review confirms the triangular hierarchy and its
corner compatibility. It produces C3 fields and, through the constraints,
C4 radius without asking for an independently prescribed higher incoming
trace. The next review concerns actual physical-slice coverage, including
radii newly enclosed by the expanding wall.
That coverage review confirms the conditional inclusion after a wall-entry
split for new radii. It adds an explicit sufficient inner-end margin for
those points. This is an analytic domain estimate only; all numerical
allowances and retained receipts remain unchanged.

### Continuation on the original coupled cell, 17 October 2026

Assemble the reviewed local null construction, the strict conditional tube,
native derivative hierarchy, quiet-core overlap and physical-slice inclusion
on the original cell and coupling interval. The active question is whether
these exclude every finite-time obstruction, not whether a smaller quiet
window exists. One Grok review audits the joint chart/tube first-exit
argument; a second checks uniform higher derivatives, terminal continuity,
native exterior gauges and restarted Cauchy data.

The owner now retains the actual curved restart line, its negative slope
bound and both short characteristic lengths. It also writes the normalized
radius-gradient flow and its Jacobian equation, distinguishing local
invertibility from coverage of the native domain. This is analytic work:
no physical evaluation, scan, new trajectory or further symbolic test is
admitted by this step. A surviving obstruction is to be recorded at its
specific equation or domain connection before another calculation is designed.
The two reviews confirm uniform terminal traces and the completed
continuation assembly after adding the subinterval and curved-restart
estimates. The original onset cell and full admitted coupling interval
are retained. The result is regional existence with the ten certified
departure bounds, not a new numerical history or an observational test.

### Parameter variation of the onset solution, 17 October 2026

The next analytic question is whether the constructed onset solution has
two controlled amplitude derivatives in a topology sufficient for its wall
geometry and clocks. Keep the same source cell, action, original C3 source,
exact quiet anchor and fixed native gauges. At alpha=0 the geometry is the
vacuum reference but the normalized scalar is the nonzero prepared probe;
the whole normalized state is not zero.

Derive first and second variations of the primitive-preserving transport,
including the moving lower endpoint, inverse wall-entry time and evaluation
of traces on the moving wall. The clipped extension used for the existence
contraction is Lipschitz and is not a differentiable parameter chart.
The available spatial derivative budget must determine the variation norm;
smooth finite-dimensional bulk coefficients alone do not establish smooth
dependence of this moving-domain solution. Independent Grok review targets
those exact shape terms and the needed difference-quotient estimate.

This admits analytic derivation and review only: no physical evaluation,
trajectory, scan or new symbolic run. Record a proved variation statement
or its precise remaining estimate before designing a numerical remainder
calculation. The onset core remains exactly quiet throughout the certified
cell; its strict travel margin does not assert immediate central arrival
just beyond the endpoint. A central waveform requires a later-domain argument.

The [whole-cap comparison](nsc-atlas-finite-amplitude.md#a-common-comparison-domain-for-complete-cap-slices)
is now supplied by the already proved polar chart. The fixed source label
is common; the final proper time is not. The positive interval from a
fixed quiet radius to the actual wall covers all old and newly enclosed
radii. Independent Grok review confirms its coverage and Jacobian rules.
The [local first-variation inverse](nsc-atlas-finite-amplitude.md#a-first-variation-inverse-on-a-full-positive-width-prefix)
retains the full spatial prefix and also passes review, including primitive
and entry compatibility. Its constants are selected in the order spatial
weight, wall norm scale, then time. A spatially truncated argument and a
parameter theorem inferred from smooth coefficients alone were not adopted.

The [complete comparison operator](nsc-atlas-finite-amplitude.md#the-cap-and-wall-equations-in-the-common-source-clock)
now retains the fixed-W wall/clock equations and the cap's radial mass/lapse
variation with its radius-dependent vacuum mass. Separate Grok derivations
confirm characteristic signs, the single-gauge lapse law, the radial
attenuation bound and the vacuum limit. At alpha zero the normalized probe
recovers the quadratic mass and lapse sources; it is not silently removed.

The [coupled inverse](nsc-atlas-finite-amplitude.md#coupled-linear-inverse-on-the-complete-comparison-domain)
uses that radial map, solves all wall-state terms in the wall ODE, and
retains instantaneous cap feedback into the exterior incoming trace.
Choose cap norm scale, exterior weight, then time. The corrected block
estimate is `N_new<=N/2`, including `kappa K/2`, and the polar primitive
propagates through its own transport and wall tangent identity. The native
chart's compatibility bracket is not used in its place. Grok's follow-up
review confirms these corrections. This is a local continuous linear
inverse, not yet a nonlinear parameter-derivative theorem.

The [spatial derivative estimate](nsc-atlas-finite-amplitude.md#first-spatial-derivatives-of-the-coupled-linear-solution)
uses the order outgoing exterior derivative, cap gradients, incoming
exterior derivative. The scalar primitive supplies f_1,eta directly.
Actual C1 initial data and corner compatibility are retained at every
restart; a quiet cap anchor does not make the whole exterior prefix or
its variations zero. Grok's final review confirms the estimate after
correcting its intermediate confusion between `p_1,Ubar` and the mixed
derivative `p_1,x=q_1,Ubar`. The owner gives the actual old/wall-born
corner equation, with no independent normal derivative of the source trace.

The [actual regularity inputs](nsc-atlas-finite-amplitude.md#source-data-and-regularity-in-the-fixed-comparison-charts)
now close independently. The source ODE is parameter smooth into spatial
C3 and has the same exact quiet initial data. The native rapidity equation
raises eta to C3; the null rates are therefore C3 and the exterior
relabeling C4, preserving a C3 fixed-gauge lapse. Grok's initial lower
derivative count was corrected by this explicit bootstrap. The signs of
the rate exponents were checked against the owned wall map. The polar
flow supplies C2 comparison first jets and uniform spatial moduli.

The [finite-difference argument](nsc-atlas-finite-amplitude.md#finite-amplitude-differences-and-the-first-parameter-derivative)
now proves Lipschitz dependence and convergence of first parameter
quotients in continuous normalized jets. Positive endpoint denominators
replace an unjustified convex-tube assumption. The cap speed difference
is retained, and interpolation includes the actual stability constant.
Common-slab induction carries the inherited quotient errors. Separate
Grok review confirms this limit and the uniform square-root quotient
rate, which gives a first-order remainder of order alpha^(3/2). Its
constant is finite and unevaluated; it does not identify a saved numerical
coefficient or certify the saved history.

Next establish convergence of first parameter quotients with one spatial
derivative. First use the proved C0 stability and uniform highest spatial
moduli to obtain continuity of second spatial derivatives on the fixed
comparison domains; spatial equicontinuity alone is not that conclusion.
Then control the differentiated quotient equations in their triangular
order, with actual inherited data and corner identities. A second
parameter equation may use that first spatial variation in a continuous
norm. Its quotient convergence and uniform bound must be proved separately;
C3 source regularity does not license a second parameter derivative into
C1 jets. These steps remain analytic derivation and review only: no new
grid, physical evaluation, symbolic run or numerical Taylor constant is
admitted until the corresponding estimate has been established.

**Parameter-estimate outcome:** the [first quotient converges in C1](nsc-atlas-finite-amplitude.md#convergence-of-the-first-quotients-with-one-coordinate-derivative)
and the [second quotient converges in C0](nsc-atlas-finite-amplitude.md#second-amplitude-derivative-and-quadratic-onset-remainder)
on the original onset cell, uniformly in the admitted amplitude interval.
An explicit C2 edge extension and second-difference estimate supply the
needed parameter continuity of the highest spatial derivatives. The radial
constraints close the cap gradient coefficients without derivatives of
its first jets. The second quotient has two copies of the speed/field-gradient
variation; the unknown second derivative remains in the principal operator.

Two independent Grok workers reviewed the respective arguments and their
final written formulas. They confirm the derivative budget, actual initial
and wall-corner data, finite-slab induction and output-specific constant.
The resulting Taylor remainder is at most `M_2,O * alpha^2 / 2`, where
M_2,O is finite and has not been numerically evaluated. No new physical
values, symbolic controls, trajectory or scan were run. This result uses
normalized fields and does not assert C1 regularity of all second varied
jets or twice-alpha differentiability of `sqrt(alpha)*f` at zero.

The next numerical question must name an observable and match the exact
source, gauge and coefficient to its implementation before selecting a
finite evaluation budget. A saved numerical coefficient is not identified
by the analytic Taylor theorem alone. Later central arrival also needs
an extended-domain estimate.


### Connecting the existing wall coefficients, 17 October 2026

The first concrete output is the Einstein wall-radius coefficient and its
polar-slice central-clock coefficient at a fixed original source label.
The existing wall and cap solvers already supply the corresponding
coefficients at fixed Jordan proper wall time. Their exterior embedding
coefficient supplies the required change of event. The source affine
label and reference wall radius must remain distinct.

Reuse the existing cap response module and solution object. Add a scalar
pullback and a wall-only query that does not evaluate a central-arrival
ray. The mathematical inputs are inverse-function differentiation,
first-order conformal conversion and differentiation of the owned source
ODE. Check the source constraint against the earlier analytic source
functional; a Hermite approximation is not a new physical input.

Manufactured controls are admitted under one cumulative allowance of
20,000 charged scalar/algebra values and 20 active process-CPU seconds.
They use an exactly invertible nonzero label perturbation, an independently
differentiated conformal radius, a known clock integral, quiet and label
rescaling limits, input guards and a ray evaluator that fails if called.
An exact source-ODE control checks the cancellation of the reduced-field
terms. No retained physical field, source integration, trajectory, grid
scan or numerical Taylor constant is admitted. Authentication checks may
read existing records and original Git bytes without running producers.

Independent Grok investigation corrected its initial incomplete search:
the nonzero first clock coefficient already exists in `CapSolution`;
`ReferenceWallGauge` alone only supplies its zero-amplitude chart factor.
The implementation will retain the complete event shift and distinguish
it from the existing central-arrival shift. A failed sign, frame, source
or query-count check stops this component for diagnosis before any
physical evaluation is designed.


**Interface outcome:** 13 manufactured controls pass with 13,380 charged
values / 0.328125 active CPU seconds. The scalar pullback and wall-only
query are implemented in the existing cap module; the original arrival
sample retains its arithmetic and schema. Forty-two source-authentication
controls pass. The retained field comparison authenticates under its
explicit `19983a7b4730d1fdc8e59ae4ba2f3efd560913a6` context with zero cap
solves. An initial current-context check correctly rejected the historical
owner-document mismatch; its published context was then used, without
weakening that check. Original cap-module bytes also authenticate against
the registered successor. Separate Grok source reviews confirm the written
identities and the implemented query, including the original arrival path.
No physical coefficient has been evaluated.

The next analytic check is whether the alpha-zero geometric subsystem
closes over the common leading scalar, without its next normalized
coefficient, and whether the old mass, wall and clock functionals satisfy
that subsystem with the same quiet constants. Separate cap/radial and
wall/exterior Grok investigations are assigned those explicit equations.
This admits reading and derivation only; no numerical coefficient,
trajectory or error-constant evaluation is introduced.

**Invariant identification outcome:** the [continuum wall and clock
functionals](nsc-atlas-finite-amplitude.md#identifying-the-continuum-wall-and-clock-coefficients)
are the actual first geometric derivatives of the constructed onset
family, after the fixed-source event shift. The exact source identity
`(r_E r_E,Ubar)'=-F/4` gives the earlier source mass and fractional lapse
coefficients with their quiet constants. The first gravitational sources
use only the common leading scalar. Improved-stress mass transport,
the angular junction, the exterior lapse and embedding equation, and
the cap radial constraints and central clock form an anchored triangular
system. Ordinary ODE uniqueness identifies its outputs without computing
the next normalized scalar or reconstructing all perturbed null coordinates.

The cap vacuum subtraction retains `3 r^2 r_1/2`. Its exact normal
conversion and proper-time conversion reproduce both conformal terms
in the old clock rate. The old absolute cap time is translated by its
inherited, alpha-independent quiet-anchor offset before comparing leading
fields; this translation leaves the first clock coefficient unchanged.
Independent Grok source review confirms the source/radial formulas and
the complete invariant-output argument. This resolves the earlier full
fixed-W reconstruction caveat for these two outputs, rather than assuming
the missing coordinate fields. No new numerical values, symbolic runs,
physical trajectories or scans were used.

The next concrete numerical target is the Einstein wall-radius and
anchored polar-slice clock first coefficient at the original right source
endpoint. First bound the required continuum leading-field values and
derivatives, reference geometry and functional integration errors. A
Volterra residual estimate must include the actual source and moving
limit; agreement between two sampled grids is not a replacement. A
read-only method investigation compares reuse of the existing dense
field with a direct short-onset integral enclosure. No physical evaluation
or Taylor-constant computation is admitted by that investigation.

## Discrete cap mass balance, 16 October 2026

**Outcome:** the [completed ledger](current-result.md#discrete-cap-mass-balance-16-october-2026)
passes creation validation and one standalone check. The finite order-eight
balance is `-2.12697432799`, including upwind work `-2.01783226502` and
remaining balance `-0.10914206297`. The corrected source freeze is `8d575f1`.
The original failed attempt plus successful creation/validation uses
11,275,091 points / 8.75 CPU seconds; the separately admitted diagnostic and
standalone check retain their own receipts. Full data stays local; the
[summary](../results/summaries/nsc-atlas-cap-mass-balance-v1.json) binds its
hash and regeneration command. No new evolution follows from this ledger.
The admission and development record below retain their original sequence.

The completed directional check below leaves
the cap mass functional and its discrete transport as the unresolved connection.
Differentiate the implemented radial quadrature, then use finite-volume
summation by parts to separate its terms. The mathematical inputs are the
chain rule, polynomial integration and a finite telescoping sum; the physical
candidate, preparation and wall law stay those of the existing calculation.

| Item | Fixed discrete balance calculation |
| --- | --- |
| Inputs and domain | The same authenticated terminal cut 128, two amplitudes, 64 departure cells, fixed stored profile and orders four/eight. Reuse the cap-consistency input reader and stored hybrid clock/wall rates. No profile fitting, new trajectory, grid scan or physical source. |
| Derivative | Reconstruct four base metrics and four frozen cap RHS vectors. A diagnostic companion in the existing package reconstructs the mass-edge recurrence, checks its primal face masses against the unchanged producer, and differentiates its actual quadrature, small-radius formulas and potential polynomial. Select a common backward-direction limiter and moving-cut branch before component weights. Include any born zero-width panel; an unresolved branch or domain condition stops that case. |
| Balance | Compare a forward tangent with reverse contraction, then telescope the departure fluxes. Retain geometry/background, both endpoints, interior faces and reaction separately; only the centre q-flux vanishes. Split interior central and upwind work, without assuming the latter has a sign. Compare the total with the same moving-wall mass target and with the already saved directional estimates. Large cancelling terms do not individually identify a faulty equation. |
| Controls | Exact polynomial mass and summation-by-parts oracles, smooth manufactured differences, C2 knot crossings, quiet extension, actual small-radius/potential derivatives, limiter openings from zero, moving-cut collisions, nonzero centre p-flux, incoming-wall sensitivity, immutable inputs and budget failures. These controls share at most 1,000,000 charged points and 30 active process-CPU seconds; no retained physical calculation is a control. |
| Physical allowance | After controls and independent review, freeze producer/control sources. One creation including validation is limited to 12,000,000 charged points and 90 process-CPU seconds; one standalone check to 6,000,000 points and 45 seconds. Charge field samples and derivative components, not just the four case count. Four edge tapes have at most 6,192 primary/edge field sites at the present panel bound; full gradient components and original metric/RHS work are additional charges. Raw data stays in ignored local storage with a tracked compact summary and hashes. |
| Decision | A failed primal, tangent/reverse or telescoping identity identifies an implementation error and stops the physical evaluation. A verified ledger must account for the total mass defect; attribute a term only against its own mathematical counterpart. Unresolved cancellation, branch sensitivity or budget exhaustion produces a diagnosis, not another scan or an altered wall law. |

Independent Grok reviews confirmed the edge recurrence and the background
moment derivatives. Review corrections retain one-sided limiter openings,
distinguish a positive profile margin from a boundary-limited zero case,
and distinguish term size from attribution of the complete defect.

Both orientations of a born panel are now controlled, including a sliver
owned by the cell to the left of the original face. Exact integration
samples at either the small-radius or potential formula threshold are
excluded conservatively; selecting an infinitesimal branch there would
need a separate argument about the two implemented formulas. The current
54 controls pass in 55 executions with 62,213 charged points / 1.578125
active CPU seconds, below the cumulative control allowance. The analyzer
and recorder receive a final independent source review before creation.

The review found no demonstrated blocker; its unread summary layout,
SHA256 and budget implementation were checked directly. Sources were
frozen at `bc70162`. The first physical attempt stopped in the zero case,
order four, at the ledger identity check: **886,949 points / 2.734375 CPU
seconds**, one metric, RHS, tangent and ledger, no accepted result. The
failure receipt recorded the phase but omitted the failed scalar identities.
Add those quantities to failure receipts without changing any mathematics
or tolerance. A single diagnostic repeat is restricted to that first case,
2,000,000 points / 15 CPU seconds, with a guard against a second metric.
Its purpose is to recover the missing failure quantities; this is not a new
four-case result or authorization for another trajectory.

That repeat, frozen at `ffcdf61`, used **886,949 points / 2.671875 CPU
seconds** and recovered the cause. Only forward-versus-reverse contraction
exceeded the bound: `5.40387e-13` versus `5.68434e-14`. The mass-ledger
difference was `4.16334e-17`, telescoping `2.77556e-17`, face decomposition
`2.33103e-18`, and primal-face difference `8.27181e-25`. The original bound
scaled by the small final terms, omitting cancellation within 128 panels.
The correction retains absolute contraction work before cancellation and
uses a path-count rounding bound, described in the mathematical owner.
No mass derivative, source, transport, wall equation or other tolerance changes.
An independent Decimal cancellation control also rejects an injected error.
The 55 current controls pass; all executions so far charge **155,213 points /
2.9375 active CPU seconds**. After independent review and a new source
freeze, resume the fixed four-case creation with the first failed attempt's
886,949 points / 2.734375 CPU seconds charged against its original
12,000,000-point / 90-second allowance. The separate one-case diagnosis
keeps its own declared receipt. Grok confirmed the accumulation bound; its
review removed the arbitrary unit floor and added a nontrivial attenuation
to the independent cancellation control. The subsequent accepted result is
recorded in the outcome above.

## Terminal cap mass-functional derivative, 16 October 2026

**Outcome:** the [completed calculation](current-result.md#terminal-cap-mass-functional-derivative-16-october-2026)
and creation validation pass with 1,711,610 points / 6.6875 process CPU seconds.
One standalone check passes with 855,805 points / 4.078125 CPU seconds.
The finite instantaneous defect remains `-0.6612504470` at the finest
backward width, compared with same-rate history `-0.6407403055`.
Its quadrature difference is `1.09497e-5`; the finest stencil has no recorded
limiter branch changes. Coarser tiny-secant changes remain disclosed.
This ends the probe set and redirects to discrete cap energy balance.
No new trajectory or grid scan is called for. The raw output stays in
ignored local storage with a [tracked summary](../results/summaries/nsc-atlas-cap-mass-direction-v1.json).
The admission below preserves the completed calculation's scope.

The completed selected-cap audit isolates a terminal temporal mass
discrepancy that is much larger than its incoming-trace flux change. The
[mass-transport identity and tangent](nsc-atlas-finite-amplitude.md#cap-mass-balance-under-field-transport)
connect the radial functional to bulk field transport and moving-interface
terms. The missing connection is its derivative under the implemented cap
transport rates. Evaluate that selected-state question before another
trajectory or a changed physical wall law.

| Item | Fixed terminal functional calculation |
| --- | --- |
| Inputs and domain | Pin the completed cap-consistency payload and its seven producer/control bindings; authenticate its predecessor payload/manifest and 104 bindings once through the existing hash-only reader. Use only cut 128, both amplitudes and quadrature orders four/eight, unchanged 64-cell raw cap states and stored profile jets. No profile fit or previous analysis replay. |
| Rates and reconstruction | Reconstruct four base cap metrics and check their saved mass/field values. Reuse the corresponding stored hybrid `r_s,T_s,Vdot_plus,w_plus`; set `R_x=r_s/Vdot_plus`, `T_x=T_s/Vdot_plus`. Evaluate the existing `background_cap_rhs` once for each base, yielding exactly four cap RHS calls. Keep these directions fixed for the local functional probes. No full coupled stage, exterior prefix/map, source solve, gauge quadrature or trajectory. |
| Directional probe | Let `dx` be the stored saved-cut spacing and `h=dx/8`. Vary `(R,T,Ip,Iq)` linearly backward along the selected rates at offsets `-2h,-h,-h/2,-h/4`. Reconstruct each cap metric with its fixed profile and order. Together with the base this is exactly 20 cap metrics per analysis, supplying backward-three-point derivatives at `h,h/2,h/4`. Each input must retain the finite profile domain, positive metric, valid radius and finite reconstructed fields. The terminal zero case has no upper-profile margin; forward probes are excluded. |
| Measurements | Retain functional derivatives, step-width and order differences, left and imposed-face fluxes, global mass balance `D_x n_W-R_x rho_n-T_x J`, and the version using each previously recorded radial derivative. Verify their difference `R_x e_r`. Compare with both saved-history time widths using the owner's `e_hist,hyb` bridge: first express the saved mass derivative with hybrid rates, then report its difference from the original measured-rate residual separately. Retain hybrid-minus-measured `R_x,T_x` differences. Record limiter secant/sign/selection ties and branch changes, radial/profile split counts and domain margins for every probe. Do not interpret three widths or a branch crossing as an established convergence order. |
| Controls | Use manufactured smooth moving functionals to check the derivative coefficients, geometric/background participation, width scaling and global-versus-radial balance identity. Exercise a one-sided branch example, immutable inputs, exact stored-rate conversion, common-budget exhaustion and first-failure retention. Synthetic dispatch controls verify 20 metrics/four RHS calls and that no forward probe or full stage is invoked. Controls share 30 process CPU seconds/300,000 points; no retained physical diagnostic is part of them. |
| Allowance | One creation including validation has 60 process CPU seconds/4,000,000 points. One standalone check has 30 CPU seconds/2,000,000 points. Static accounting from the prior metric implementation estimates fewer than one million points per analysis; this is not a measured guarantee. Freeze producers and controls before creation. |
| Decision and stop | A stable instantaneous residual comparable with the saved-history discrepancy redirects to the cap's discrete mass/energy balance. A small instantaneous residual redirects to time integration or saved-history differentiation. Material width/order/branch sensitivity leaves attribution unresolved and ends this probe set. Stop on the first invalid input/domain, nonfinite field, failed reconstruction or budget exhaustion with an ignored diagnostic receipt. No new evolution or wall law follows from this calculation alone. |

The [analyzer](../scripts/analyze_nsc_atlas_cap_mass_direction.py) and
[recorder](../scripts/diagnose_nsc_atlas_cap_mass_direction.py) pass independent
source review. All 26 controls pass before creation. The first control run
stopped after six passes at a synthetic bridge discrepancy of `5.55112e-17`
against an unscaled `5e-17` threshold. Replacing that sub-epsilon cutoff
with a rounding allowance scaled to the operands corrected the control;
the numerical implementation was unchanged. The affected orchestration
test and remaining unexecuted controls then passed. Across both attempts,
the shared allowance records **171,957 points / 1.34375 process CPU seconds**,
including exhausted child attempts. No retained physical calculation was
part of the controls. Producers and controls were frozen at `fc71709` before creation.

## Selected cap mass and wall-curve audit, 16 October 2026

**Outcome:** the [completed audit](current-result.md#selected-cap-mass-and-wall-curve-audit-16-october-2026)
and creation validation pass on the first attempt with 4,393,156 points /
15.359375 process CPU seconds. One standalone check passes with 2,196,578
points / 8.453125 CPU seconds. Pulse-time residuals are strongly sensitive
to radial derivative width; quadrature-order sensitivity is small. At the
terminal cut, radial mass/lapse residuals are small while the temporal mass
residual remains and is not explained by the incoming-trace flux mismatch.
This redirects to the cap mass-functional/discrete-transport connection
before another evolution or changed wall law. The admission below remains
the fixed scope of the completed calculation.

The saved exterior-row audit localizes the sensitivity but leaves its
origin ambiguous. Before another evolution or boundary-stencil change,
test the independent cap temporal mass equation and the actual wall
curve. The [cap consistency derivation](nsc-atlas-finite-amplitude.md#independent-cap-mass-and-wall-curve-consistency)
reuses the chain rule, polar-areal Einstein equations and full acceleration
formula already in this owner chain. The missing connection is their
evaluation on reconstructed saved cap data, with actual left field traces.

| Item | Fixed cap consistency calculation |
| --- | --- |
| Inputs and domain | Authenticate the completed operator comparison and all its 104 bindings once. Use its zero-backreaction and finite cases at target cuts 63, 68 and 128. Reconstruct both quadrature orders 4/8 at all 16 unique neighboring cuts: 61--70 and 122, 124--128, for both amplitudes. Instantiate the unchanged `CharacteristicProfile` directly from stored labels, values and two jets. No source/profile fit or prior trajectory reconstruction. |
| Metric and derivatives | Reuse `background_cap_metric` and `background_cap_profile`, with 64 cells and the saved cap states and wall `R,T`. Differentiate returned `n_faces` and `delta_faces` at the wall using backward-three-point first derivatives at radial strides one/two. Use centered-three-point temporal derivatives at strides one/two for targets 63/68; use backward-three-point first and backward-four-point second temporal derivatives at 128. All supports stay inside the saved domain. The matched order-4/order-8 neighbor sets control temporal as well as radial quadrature sensitivity. |
| Boundary inputs | At each target, amplitude and quadrature order, reconstruct the exterior prefix with the unchanged quadratic departure gradient and exact saved source/gauge rows, then evaluate the hybrid wall map: 12 selected prefix/map calls, reusing target cap metrics. Its incoming characteristic and force predictions do not require the log-lapse proper derivative, so one map call with that auxiliary input zero suffices for the selected outputs. No full coupled stage, cap RHS, source integration, gauge quadrature or time evolution. |
| Measurements | Record normalized and actual temporal/radial mass residuals, lapse residual, reconstructed left `Pi,Phi`, actual imposed hybrid incoming characteristic and a separately labelled wall-curve kinematic trace. Retain flux changes between those traces. Derive the curve proper-time factor from measured `R_x,T_x,A,delta`; compare it with the derivative of saved `s`. Use measured second time derivatives and radial/temporal mass derivatives in the full geometric acceleration, then compare with the hybrid prediction. Do not substitute a desired constraint, velocity or acceleration to define its own test. |
| Controls | Synthetic moving profiles check the fixed-radius chain rule, analytic vacuum subtraction, potential quotient, normalized alpha-zero limit, radial and temporal stencil coefficients and actual-curve acceleration. Check target/neighbour indices, input mutation, common-budget exhaustion, first-failure retention and that no full stage/RHS/evolution is invoked. Include a scalar trace mismatch that changes the flux without changing the primary left-trace residual. Controls share 30 CPU seconds/1,000,000 points; no retained physical diagnostic is part of them. |
| Allowance | One creation including validation has 120 process CPU seconds/8,000,000 points. One standalone check has 60 CPU seconds/4,000,000 points. Static helper accounting estimates about 2.6 million points per analysis for 64 cap metrics, 12 prefix/maps and selected reconstructions; that is a planning estimate, not a measured guarantee. Freeze new producers and controls before creation. |
| Decision and stop | A residual stable against these widths and integration orders identifies a cap/interface consistency question; a strongly changing residual identifies unresolved numerical differentiation or reconstruction. Trace mismatch is a discrete boundary fact, not proof of a new source law. Retain the cap contributions in the J/O identities. Stop on nonfinite data, invalid profile/domain, non-timelike reconstructed curve, failed reconstruction or budget exhaustion, with an ignored failure receipt before any revised calculation. No new evolution, physical coupling or wall equation is admitted by this diagnostic. |

The 34 synthetic mathematical, orchestration and record controls pass on
their first combined run with **62,119 aggregate points / 0.34375 process
CPU seconds**, including deliberate exhausted attempts. The acceleration
control contracts an independently calculated Christoffel acceleration
with the outward normal. Measured `D_x A` enters the curve normalization;
its difference from an algebraic chain-rule reconstruction is retained.
All controls use one cumulative allowance and no retained physical replay.
The new analyzer, recorder and controls are frozen before creation.

## Saved near-wall row audit, 16 October 2026

**Outcome:** the [fixed saved-row result](current-result.md#saved-near-wall-row-audit-16-october-2026)
and creation validation pass on the first attempt with 97,350 points and
1.59375 process CPU seconds. One standalone check passes with 48,675 points
and 1.5625 CPU seconds. Curvature/lapse mismatch is concentrated near the
wall and sensitive to derivative width. Its origin remains ambiguous:
some third differences decrease under the changed evolution even where
the residual grows. End this same-grid scan and use the independent cap
mass/lapse consistency check before another evolution or wall-law change.
The admission below remains the fixed scope of this completed calculation.

The fixed-grid comparison below localizes its incoming and outgoing maxima
at `eta=63/64`, while native and wider wall-curvature estimates disagree.
The next calculation inspects that numerical connection before changing a
boundary law or evolving another window. Polynomial interpolation supplies
the derivative identities; the missing information is their contribution on
the actual saved fields.

| Item | Fixed saved-array diagnostic |
| --- | --- |
| Inputs and domain | Authenticate the operator comparison, its manifest and all 104 bindings once. Reuse its bound 64-cell predecessor and exact 1,025-row source/gauge cache. Inspect cuts 63, 68 and 128, nodes 56--64 at the wall and 28--36 as an interior control, for both methods and both amplitudes. Cuts 63/68 are the retained incoming/raw-outgoing maximum cuts; 128 tests the terminal first-jet discrepancy. No scan selects another cut after seeing this result. |
| Measurements | Retain signed curvature, background-lapse, background-radius, departure-product and scalar-stress contributions to the incoming constraint. Compare independent centered derivatives on interior patch nodes and the existing `p3`/`p5` wall-estimate families at strides one/two (three/four nodes for the first/second derivative in `p3`, five for both in `p5`). Cross-check retained incoming constraints where available. Report the signed third-difference identity between selected and backward/centered quadratic gradients; distinguish changing a stencil on identical arrays from comparing different evolved arrays with the same diagnostic. The outgoing maximum selects a cut but no new time-derivative reconstruction is admitted here. |
| Reuse and limits | Reconstruct the physical wall endpoint from the saved wall and analytic vacuum/source offsets, with no required wall derivative substituted. Reuse the frozen vacuum and wall-estimate helpers. No stage, profile fit, source solve, gauge quadrature, complete-history diagnostic or trajectory replay. This audit can localize a sensitive boundary shape; it cannot by itself identify whether that shape was generated by the closure or by cap/interface data. |
| Controls and allowance | Controls share 30 process CPU seconds/100,000 charged points. Check constant through cubic polynomials, source offsets and scalar stress, vacuum cancellation, row/channel orientation and a penultimate-node impulse: native wall curvatures receive coefficients -5 and -26/3 times the impulse divided by physical spacing squared, while stride-two estimates omit it. Check mutation, invalid inputs and budget rejection. Creation including validation has 30 CPU seconds/150,000 points; one standalone check has the same allowance. Freeze the new producer and controls before creation. |
| Decision and stop | A localized, width-sensitive last-row contribution redirects attention to the discrete closure. A smooth, width-consistent residual supports the independent cap mass/lapse consistency audit. Ambiguity ends this same-grid diagnostic branch and also redirects to the independent cap check. Stop on invalid input, nonfinite output, failed reconstruction or exhausted budget, retaining an ignored failure receipt before designing any retry. This result supplies neither a new boundary condition nor a converged junction certificate. |

The eight synthetic controls pass on their first attempt with 60,460
aggregate charged points and 0.15625 process CPU seconds, including 195
points in the deliberately exhausted child budget. Independent source and
test review passes. No retained physical diagnostic was part of those
controls; the new producer and controls are frozen before creation.

## Fixed-grid exterior-operator comparison, 16 October 2026

**Outcome:** the [completed pair](current-result.md#fixed-grid-exterior-operator-comparison-16-october-2026)
and creation validation pass on the first attempt with 55,398,966 total
points and 256.25 process CPU seconds. One standalone check passes with
4,040,624 points and no trajectory. Final finite kinematic defects shrink
about 21--33 times, and paired wall-scalar method sensitivity is 0.524076%
of the new peak. Incoming-constraint RMS decreases, but its maximum at
`eta=63/64` rises slightly; native wall-curvature peaks worsen while
wider estimates improve. The measured next question is the near-wall
discrete closure, using saved arrays before another evolution. Cap
constraints cannot yet be inferred from these unresolved second jets.
The fixed admission and implementation controls are retained below.

This comparison may execute after the derivative, coupled-stage and callback
driver controls below pass and their independent review is resolved. It tests
whether the improved exterior derivative reduces the measured errors with
the physical candidate and both spatial grids held fixed.

| Item | Fixed same-grid comparison |
| --- | --- |
| Inputs | Authenticate the completed 64-cell spatial payload, its manifest and all 92 bound paths. Read its bound 32-cell ancestor only for the unchanged 1,025-row cache and computational profile. Reuse the original 513-row source/gauge arrays for saved-history diagnostics. No old trajectory, source solve, cache reconstruction or profile-stage extraction. |
| Changed connection | Use the new common quadratic departure gradient in both prefix jets and ALE rates through the explicit stage/SSPRK3 callbacks. All physical equation helpers, source's 129 Hermite knots, profile's 65 knots, amplitudes, cap discretization and wall laws stay fixed. |
| New cases | One zero and one original finite-amplitude window, each at 64 exterior/cap cells and 512 SSPRK3 steps, with the retained exact midpoint keys and 129 saved cuts. Two final reconstructions are included. No 128-cell case or profile repair. |
| Measurements | Compare raw zero, raw finite and finite-minus-zero fields, wall observables and total cap integrals directly at the same saved coordinates. Cap radii can differ; no adjacent-cell summation or interpolation is needed between the two 64-cell methods. Reuse old native-64 constraint grids and independent wall estimates; compute new 129/65-cut constraints and three-/five-point wall estimates at strides one/two. Retain the actual inherited outgoing-source constraint, quiet scalar values, outgoing cap trace, profile margins, wall kinematic defects and common-wall flux identities. |
| Allowance | One creation shares 450 process CPU seconds/80,000,000 points, including the two evolutions and creation validation. One standalone check has 60 CPU seconds/8,000,000 points without time evolution. Recorder controls share 30 CPU seconds/1,000,000 points; analyzer controls have a separate 30 CPU seconds/1,000,000 points. The previous 60.09-million-point creation leaves allowance for the new charged derivative while this comparison reuses old diagnostic arrays. Stage totals are measured anew. |
| Verification | Preserve exact caches, coordinates and model; reconstruct saved physical arrays from raw states; recompute two final stages and diagnostics without replaying either trajectory. Controls check changed-input rejection, the actual new stage callback, raw/paired same-grid differences, and equal-width constraint/wall comparisons. Freeze producers before creation. |
| Stop and decision | Stop at the first domain, nonfinite, metric, wall, Courant `.4`, reconstruction or budget failure, with an ignored failure receipt and no retry. If constraints improve while wall observations remain stable, investigate the remaining cap/interface or source/profile accuracy using the measured results. If they worsen or oscillate, diagnose the new closure before another run. A method comparison supplies sensitivity, not a convergence order, continuum limit or complete junction certificate. |

Implementation and independent source/control review pass before creation.
The 24 recorder controls use 26,027 shared points/2.125 CPU seconds;
the seven analyzer controls use 587,907 points/0.4375 CPU seconds.
Both pass on their first attempts without evaluating retained physical
histories. Together with the operator and coupled-stage suites below,
89 controls pass, reporting 2,047,112 aggregate points and 3.84375 CPU
seconds. Producers and controls are frozen for the admitted comparison.

## Second-order exterior departure operator, 16 October 2026

The 32/64-cell comparison improves constraints at equal diagnostic width;
the remaining finite wall kinematic defects approximately halve. The
existing exterior derivative is first order and supplies both the prefix
jets and moving-grid rate. Improve that numerical connection with the
physical source, cap laws and computational profile held fixed.

For `h=1/N` and negative grid velocity, differentiate departures by
`G_i=(-3Z_i+4Z_(i+1)-Z_(i+2))/(2h)` through `i=N-2`,
`G_(N-1)=(Z_N-Z_(N-2))/(2h)`, and
`G_N=(3Z_N-4Z_(N-1)+Z_(N-2))/(2h)`. These are derivatives of local
quadratic interpolants. Use the same `G` in the transverse prefix jets
and the ALE term. The physical wall value supplies inflow; no prescribed
wall acceleration or required derivative replaces the reconstructed jet.

| Item | Fixed operator development and controls |
| --- | --- |
| Reuse and missing connection | Reuse the exact vacuum, departure prefixes, cap metric/RHS, hybrid wall equations and SSPRK3 coefficients. A new explicit stage core accepts the departure-gradient function; a callback SSPRK3 driver supplies it at every stage and final diagnostic. Frozen stage/driver bytes remain replay inputs. |
| Alternatives | Refactoring frozen owners would invalidate retained bindings. Calling the old stage then replacing rates leaves dependent prefix, wall and cap quantities inconsistent unless the full stage is repeated. The selected explicit callbacks copy only coupling/time orchestration, keep the physical equation owners shared, and become the extension point for subsequent operators. |
| Fixed scope | At least four exterior cells, bounded by existing solver limits; uniform closed eta grid, positive `W,W_x`, source departure exactly zero. Preserve source/cache coordinates, cap representation, wall normalization, primitive conservation and all model parameters. No physical evolution or input regeneration is admitted by these controls. |
| Derivative controls | Constant/linear/quadratic exactness, analytic vacuum preservation, endpoint orientation, invalid/nonfinite shapes, shared-budget exhaustion, and smooth manufactured inflow/advection at fixed negative velocity and the actual `a(eta)=eta*b` grid flow. Compare 16/32/64 cells with prescribed analytic inflow and bounded SSPRK3 steps. Include a finite homogeneous amplification check; it is not a coupled PDE stability proof. |
| Stage and driver controls | First-order callback reproduces the complete frozen stage to arithmetic precision; the selected gradient reaches both prefix and ALE consumers. Independent active-scalar common-wall-flux identities, all three stages plus final callback, exact midpoint keys, unchanged snapshots, and failure/count/budget propagation. |
| Allowance | Derivative and manufactured-advection controls share 45 CPU seconds/1,500,000 points. Stage/driver controls share 30 CPU seconds/300,000 points. No prepared source solve, saved physical diagnostic replay or coupled research trajectory. Every attempted control is counted; an unexpected failure is diagnosed before rerunning its affected scope. |
| Decision and failure forecast | An inflow instability or nondecreasing smooth advection error redirects the boundary stencil. A failed first-order equivalence or common-wall identity redirects the adapter before numerical research. Excessive new derivative sensitivity requires a measured local comparison before any full window. Passing controls permits designing a separately recorded fixed-grid operator comparison; it does not establish a continuum solution or inherit a stability theorem from the previous stencil. |

The unchanged profile's prefront leakage and the source's curvature-trace
convention remain separate work. This development does not silently fix
either one.

The derivative's 21 focused controls pass on their first run, using
1,405,762 ordinary points plus a deliberately exhausted five-point attempt,
1,405,767 counted points in total and 0.78125 process CPU seconds. On
16/32/64 cells the smooth constant-speed maximum errors are
`0.0248193 / 0.00766079 / 0.00241073`, versus the first-order
`0.145860 / 0.076993 / 0.039708`. The corresponding variable-speed
errors are `0.0325281 / 0.00953550 / 0.00298090`, versus
`0.136196 / 0.074108 / 0.038693`. Refinement ratios exceed three for
these smooth examples; they are not an asserted coupled convergence order.
At 64 cells the wall-jet errors also decrease, although the quadratic
wall jet is worse at 16 cells. Finite homogeneous tests at Courant numbers
at most 0.2991 have maximum measured L2 gains below 1.15, approaching one
on the finer grids. These tests use three selected initial shapes and do
not establish a uniform norm bound or nonsmooth-front convergence.

The stage and callback-driver suites also pass on their first run:
37 controls, 27,411 aggregate points and 0.5 process CPU seconds.
Supplying the first-order callback reproduces every inherited stage field
and all three old SSPRK3 stages exactly. Active canonical/conformal wall
flux checks pass; the one reconstructed gradient reaches both consumers.
Nonlinear and nonautonomous manufactured time-step ratios approach eight,
consistent with the retained third-order time scheme on those examples.
No coupled research window is part of these controls.

## Fixed-input spatial refinement, 16 October 2026

**Outcome:** the [two new evolutions](current-result.md#fixed-input-spatial-refinement-16-october-2026)
complete on the first attempt with 60,087,926 total points and 256.8125
CPU seconds, including 6,684,624 points for creation validation. One
standalone check passes with 6,684,624 points and no trajectory replay.
The paired wall-scalar spatial difference is `2.69431e-7`, or `0.774698%`
of the finer response peak; geometry and clock differences are below
`0.071471%` of their respective peaks. The incoming constraint maximum
decreases from `5.54765e-4` to `3.53518e-4`; at identical diagnostic
width it is `3.31042e-4` on the new solution. The source-subtracted
outgoing maximum decreases from `7.78011e-3` to `4.93833e-3` (native)
or `5.15143e-3` (equal width).

This is an actual spatial improvement, not merely a changed diagnostic.
Wall curvature still changes substantially with width and evolution:
equal-width weighted wall differences reach `1.84444e-4` and
`2.22042e-4` for the three- and five-point alternatives. The profile's
quiet leakage and source-trace sampling sensitivity remain. No automatic
third grid or junction-failure conclusion follows. The admission fixed
before the calculation is retained below.

The completed wall audit leaves a spatially unresolved curvature estimate.
Another first-derivative balance on the same saved grid could avoid an
explicit second derivative, but would still share its evolution error.
Obtain independent spatial data before interpreting the wall constraint.

| Item | Fixed spatial comparison |
| --- | --- |
| Inputs and model | Hash-authenticate the completed 32-cell/512-step record, its manifest and all 86 bound inputs. Reuse its 1,025 exact source/gauge cache rows and original computational profile inputs. Hold the action, amplitudes, source's 129 Hermite knots, profile's 65 knots, physical window, gauge and boundary laws fixed. No source solve, cache reconstruction, profile-stage extraction or old trajectory replay. |
| New cases | Zero and original finite amplitude at 64 exterior and cap cells, each with the same 512 SSPRK3 steps and cached explicit midpoints. Retain 129 accepted cuts: 1,024 steps and 3,074 stages including two final reconstructions. Existing operators support this grid without modification. |
| Comparisons | Common saved time and exterior eta coordinates; sum adjacent total physical cap integrals at common normalized radial cells. Retain raw and finite-minus-zero wall geometry, clocks and scalar histories, with their full spatial differences and peak locations. Physical cap radii can differ and are identified as such. |
| Constraint controls | Retain native 129/65-cut diagnostics and a 32-node-interval restriction of the 64-cell exterior. Compare the restricted diagnostic with the old 32-cell result using identical derivative widths, then compare native and restricted derivatives on the new solution. Separate inherited `r_source E_source/r` from the outgoing residual. Retain bottom-origin, wall-origin and separator masks. Incoming constraints are not presumed zero at a sourced wall. |
| Wall curvature | Compare three- and five-point one-sided estimates on the new native grid and its every-other-node restriction. The latter uses the same physical eta widths as the old 32-cell stencil. Compare weighted wall values at identical saved cuts; no new birth-transport interpolation is needed to test their spatial stability. |
| Other diagnostics | Preserve the own-trace zero-cap comparison, outgoing cap signal, quiet scalar controls, profile-domain margins, final wall kinematic defects and constructed flux identities. The known profile-support defect is not corrected within this comparison. |
| Allowance | One creation has a shared 450 process CPU seconds/80,000,000-point ceiling, including input authentication, unchanged-profile construction, both trajectories, diagnostics and creation validation. The prior pair's 30.83 million trajectory points and cap-panel growth suggest 45--65 million points and 250--350 CPU seconds, not an execution guarantee. One standalone check has 60 CPU seconds/8,000,000 points and no trajectory. Focused new controls share 60 CPU seconds/2,000,000 points. |
| Verification and stop | Independent restriction/derivative examples, raw-to-physical reconstruction of all saved states, two final stages, exact cache/time/grid preservation, zero and known-limit checks, changed-input/output rejection and shared-budget accounting. Stop at the first domain, nonfinite, metric, wall, Courant `.4`, consistency or budget failure and preserve an ignored `.build` receipt. No automatic retry or resolution increase. |
| Decision | Determine whether the field and incoming-constraint differences decrease under actual spatial refinement and whether wall curvature becomes consistent at equal widths. If not, identify the remaining spatial operator, source regularity or cap/interface mismatch before changing one of them. Two grids do not establish an order, a continuum limit or a propagated error enclosure. |

Implementation and independent source review are complete. The 24 recorder
controls pass on their first run with 21,311 points and 2.171875 CPU
seconds. The 24 diagnostic controls pass with aggregate 943,303 points
and 0.75 CPU seconds, including an isolated rerun after correcting an
expected-array shape in the quartic test and the final 129-cut/64-versus-32
orchestration control. The latter stubs repeated summaries and profile
work while exercising actual restrictions, comparisons and wall curvature;
it is not a physical evolution. Combined measured control work is
964,614 points and 2.921875 CPU seconds, plus a deliberate 99-point
exhaustion attempt. New sources were frozen before the evolutions;
existing results and their bound producers remain unchanged.

## Incoming wall curvature and transported constraint, 16 October 2026

The preceding saved-data audit finds a much larger incoming constraint on
wall-origin stencils than on bottom-origin stencils. Its missing boundary
datum is the wall's pure transverse curvature. Calculate that datum from
saved departures, without inserting the hybrid solver's prescribed wall
curvature or assuming a zero constraint.

| Item | Fixed wall diagnostic |
| --- | --- |
| Inputs | Hash-authenticate the 512-step refinement record and its 86 bound paths, plus the saved constraint-transport audit and its producer/control bindings. Use the original two amplitudes, 129 cuts and 32-cell eta grid. No inherited scientific validator or source/cache/profile/stage/trajectory replay. |
| Wall reconstruction | Evaluate the existing analytic Schwarzschild vacuum interface only at the 129 reference wall coordinates `(W(x),x)`, sharing these queries between amplitudes. Recover endpoint departures from physical wall values minus vacuum values and the source offset. Interior departures are already stored. |
| Independent derivatives | Add analytic vacuum first/second jets to one-sided derivatives of departures, with `d_Ubar=d_eta/W`. Compare four specified stencils: three-point first/four-point second derivatives, and five-point first/second derivatives, each at eta strides one and two. Form `E_Ubar,W=r_Ubar_Ubar-L_Ubar*r_Ubar+alpha*r*u_Ubar^2/(2D)` and its separate terms. No prescribed acceleration or wall rate enters. |
| Transport | Tabulate `B(W(x))=r_W E_Ubar,W`. Interpolate it piecewise linearly to each wall-origin interior `Ubar`, using 129 wall samples or the 65 alternating samples. An interpolation bracket may use only cuts at or before the observation cut. Keep unsupported points explicit and compare only common supported points; no extrapolation or future wall datum. |
| Independent subsets | Retain all wall-origin stencils, and separately those whose three-point eta stencil is disjoint from the union of all four boundary stencils. The latter boundary union starts at eta index 24 of 32. Report overlapping-stencil results separately, because shared field samples can correlate derivative errors. |
| Comparisons | Keep raw interior `r E_Ubar`, transported wall prediction and signed remainder. Vary the wall derivative stencil with interpolation held fixed, and vary interpolation sampling with the wall stencil held fixed. Retain maxima, RMS and locations. These are diagnostic sensitivities, not evolution refinements or an error enclosure. |
| Allowance | Creation and one standalone check each use at most 30 process CPU seconds and 2,000,000 arithmetic points, including authentication and the 129 analytic vacuum evaluations. Focused controls share 30 CPU seconds and 200,000 points. No numerical source or time evolution. |
| Checks and stop | Independent polynomial derivative oracles, nonzero transported boundary data, interpolation causality/coverage, disjoint-stencil masks and changed-input/budget controls. Stop on first inconsistent reconstruction, domain/nonfinite failure, invalid bracket or exhaustion; preserve an ignored `.build` receipt. Output is creation-only. |
| Decision | Determine whether an independently measured wall datum accounts for the interior incoming constraint, and whether that conclusion survives derivative/interpolation changes away from overlapping stencils. If it does, evaluate its inherited outgoing-source contribution and the corresponding cap/junction residuals. Nonzero wall constraint alone is not a junction failure: with zero cap constraints the smooth identities give `J_+=(A_-/H^2)O_+`. If transport fails, diagnose spatial derivatives or mixed transport before changing the physical model. |

**Outcome:** the [saved wall audit](current-result.md#incoming-wall-curvature-and-transported-constraint-16-october-2026)
completes with 1,391,719 points and 1.046875 CPU seconds; its one
independent standalone check passes with the same point count. Each uses
129 analytic vacuum evaluations and no source, cache, profile, stage or
trajectory replay. All 19 focused controls pass on their first run with
71,773 points and 0.0625 CPU seconds. Producer and controls remain frozen.

All 916 wall-origin points have causal interpolation brackets at both
samplings. On the 162 disjoint-stencil points, raw weighted incoming
constraint RMS is `7.16444e-5`; every wall-subtracted alternative increases
it, to `8.14450e-5`--`1.38386e-4`. Spatial-stencil changes in transported
wall data reach `1.08923e-4`--`2.67968e-4`, whereas 65/129 interpolation
changes reach only `2.53441e-6`--`3.86554e-6`. The disjoint subset covers
only `x=0.45157688..0.48086746`, `eta=0.53125..0.6875`, excluding the
earlier near-wall maximum. Zero controls remain at about `1e-12`.

The wall datum is not spatially resolved by this comparison and does not
reliably account for interior transport. An aggregate improvement on
overlapping stencils cannot settle that result. The next calculation
addresses spatial differentiation or bulk transport before a geometric
junction-acceleration interpretation. The conditional outgoing-source
baseline above still requires control of the actual cap constraints and
acceleration.

## Saved constraint transport and sampling, 16 October 2026

The 512-step comparison reduces temporal scalar sensitivity while leaving
the finite sampled null constraints. Before changing the evolution, separate
the source constraint already supplied by the Hermite representation from
the residual inside the domain. The diagnostic uses the existing exact
transport identity and only retained arrays.

| Item | Fixed saved-data calculation |
| --- | --- |
| Inputs | Hash-authenticate the completed `nsc-atlas-ssprk3-refinement-v1` payload, manifest and all 86 bound paths. Use its two saved histories, source/gauge cache and existing 129/65-cut departure diagnostics. No recursive physical validator, source evaluation, profile fit, stage or trajectory. |
| Missing connection | Compare the measured `E_x` to `r_source E_source/r`, retaining both the inherited source contribution and the remainder. Calculate `r E_x-r_source E_source` separately. Reconstruct radius from saved physical increments and the saved source radius, and independently form the actual source constraint from its cached first/second jets. |
| Domain and sampling | Same zero/finite cases, 32 cells and original finite interval. Compare 65-cut diagnostic values with the identical interior cuts and eta nodes of the 129-cut diagnostic. Keep signed arrays, maxima and locations; use the union of their front-crossing time-stencil masks. This is diagnostic differentiation sensitivity, not time-evolution refinement. |
| Incoming constraint | Partition the three-point eta stencil by `Ubar=eta W(x)` relative to `W(.4)`: wholly bottom-origin, wholly wall-origin, or crossing their separator. A quiet bottom supplies zero incoming constraint. The wall-origin value remains unspecified; the saved first-jet wall defect cannot supply pure transverse curvature. Temporal coarsening leaves this spatial stencil unchanged. |
| Regularity convention | Every saved cut is a source Hermite knot. At the interior cuts the cached second derivative uses the right trace; centred differences of evolved departures mix neighboring pieces. The subtracted remainder retains this sampling issue and is not identified with continuous evolution error. |
| Allowance | One creation uses at most 30 process CPU seconds and 1,000,000 arithmetic points including authentication and saved-array validation. One standalone check has the same allowance. Focused manufactured controls share 30 CPU seconds and 300,000 points. No physical source or evolution is repeated. |
| Checks and stop | Independent nonzero-source transport examples, coordinate/index/sign checks, stencil partitions, source-jet identity, changed-input/output rejection and shared-budget controls. Stop on the first mismatch, nonfinite value, invalid radius/grid or exhaustion; preserve a diagnostic receipt under ignored `.build`. Creation requires a fresh path. |
| Decision | Determine how much of the sampled constraint is inherited and how much changes with the diagnostic stencil. Then choose a source representation, spatial derivative or independent incoming-wall diagnostic according to the measured remainder. No automatic time doubling or altered physical assumptions. |

**Outcome:** the [saved audit](current-result.md#saved-constraint-transport-and-sampling-16-october-2026)
completes with 689,270 points and 0.921875 CPU seconds. All 21 controls
pass with aggregate 72,316 points and 0.359375 CPU seconds, including
four rerun mutation controls after an isolated fixture-alias correction.
Independent review passes. The finite inherited `E_x` maximum is
`1.29723e-3`; the subtracted maximum remains `7.78011e-3`. Matched
129/65 diagnostic sensitivity reaches `1.34189e-3`. The incoming
constraint is largest in the wall-origin region (`5.54765e-4`, versus
`1.15980e-5` on bottom-origin stencils). This selects an independent
one-sided wall-curvature and transport diagnostic before a change to the
evolution. The source and diagnostic trace conventions remain explicit.
One standalone saved-arithmetic check passes with 689,270 points and
no source/cache/profile/stage or trajectory replay.

## Midpoint-complete SSPRK3 comparison, 16 October 2026

**Outcome:** the [512-step pair](current-result.md#midpoint-complete-ssprk3-comparison-16-october-2026)
completes with 174.90625 process CPU seconds and 33,604,170 points.
Standalone authentication passes with 1,418,330 points. The paired scalar
256/512 difference is `5.84100e-7`, or 1.67616% of the finer peak, versus
the earlier `7.57345e-6`. This reduction does not establish an order or
uniform bound. Finite constraints and the quiet cap-profile defect persist.
The old SSPRK2 extrapolant is farther from the new 512-step history than
from the 256-step history. No further automatic time refinement follows;
the next connection is the identified profile support and spatial/source
constraint work. The original fixed admission follows.

The saved four-case histories localize 90% of the squared paired scalar
time difference to cuts 89--98. On the zero reference, direct Simpson
quadrature of the known background chain has defects `2.06210e-6` and
`1.28889e-7` at 128/256 steps; their difference closely follows the
actual zero-wall time difference. Along the finite histories, profile
knots lie inside 45 coarse intervals, which contain 99.9886% of the
squared increments of the paired scalar time difference. This is
saved-data localization, not a causal isolation or an assumed smooth
convergence order. It motivates one additional time resolution with the
profile and every spatial/physical input held fixed.

| Item | Fixed midpoint-complete comparison |
| --- | --- |
| Borrowed operators | Existing SSPRK3 coupled operator, the same 129-knot `SourceSnapshot`, and the analytic reference-wall gauge with fixed Gauss quadrature. No field sample is differentiated to obtain a jet. |
| Missing connection | The 513-row cache lacks the midpoint values needed for 512 SSPRK3 steps. Evaluate the existing source interface there; do not change its source ODE solution or interpolation. Evaluate the reference gauge directly at those labels. |
| Cache | Interleave the original 513 coordinates, preserved exactly, with their 512 explicit midpoints. For each amplitude reconstruct the same source snapshot from the saved 257 exact rows restricted to the original 129 knots; require its values and jets at all old coordinates to match every original cache row. Sample only the new midpoints afterward. The scalar preparation remains analytic; source geometry/transverse jets retain their existing Hermite approximation and residuals. |
| Gauge control | Preserve all original gauge rows. Compute the 512 new rows at order 64 and independently at order 32, charging every query, inversion and quadrature node. Compare all physical values/jets to absolute tolerance `5e-12`; iteration counts need not agree. Check clock orientation, normalization and original inverse-residual bounds. These differences are a quadrature sensitivity check, not an enclosure. |
| Cases and domain | Exactly two new cases: zero and the original finite amplitude, 32 cells/512 steps over the original interval; 129 saved slices, 1,024 steps and 3,074 stages including two final reconstructions. Reuse the original computational profile without support edits. |
| Comparison | Raw and paired 256/512 histories, with the retained 128-step differences and the old SSPRK2 extrapolant kept distinct. Locate maxima and inspect profile-crossing intervals. Retain constraints, wall compatibility, outgoing cap signal, quiet controls and own-trace cap pullback. A small final endpoint difference cannot stand in for the history maximum. |
| Shared creation allowance | 300 process CPU seconds and 50,000,000 points for authentication, new cache and both gauge orders, profile construction, both trajectories, diagnostics and validation together. Source-snapshot construction and old/new queries charge 2,314 points before later validation. No new source ODE integration or old trajectory replay. |
| Authentication | Verify preserved old keys/rows and reproduce the midpoint cache with the same source interface and both gauge quadratures; reconstruct all saved physical views, two final stages and diagnostics. One standalone check has 60 CPU seconds/5,000,000 points and performs no trajectory. Producer/operator inputs and original results remain frozen. |
| Controls and stop | New cache/recorder controls have one shared 60-CPU-second/2,000,000-point allowance. Include exact-key preservation, independent scalar/source-jet and gauge identities, mismatched old rows, bad midpoints, output tampering and budget failures. The first domain, arithmetic, profile, metric, wall, Courant `.4`, comparison-tolerance or budget failure stops the calculation with an ignored `.build` receipt; no automatic retry. |
| Decision | Determine the remaining scalar time sensitivity on the finer grid and whether its response can be separated from spatial constraints. A poor comparison prompts a profile-aware temporal discretization question, not an automatic larger grid. The confirmed early profile-support defect stays a separate representation correction. |

The cache and recorder implementations have passed 36 focused controls and
independent review. Combined control work uses 540,731 points and 2.890625
process CPU seconds, plus an intentional 513-point exhaustion attempt.
This includes the isolated correction of a test-fixture alias. Controls
include actual analytic compact-scalar jets, manufactured
source curvature, independent gauge mass/clock identities, exact old-row
preservation, paired-history indexing, tampering and shared-budget checks.
The producer and its 86 input/source bindings remain frozen with the
completed comparison. No new source integration was involved.

## Fixed-input SSPRK3 coupled window, 16 October 2026

**Outcome:** the [four-case comparison](current-result.md#fixed-input-ssprk3-coupled-window-16-october-2026)
completes with 130.734375 process CPU seconds and 24,939,592 points.
One standalone authentication passes with 956,575 points, four final
radial stages and no trajectory replay. Geometry and clock differences
are small, but the paired scalar 128/256 difference is `7.57345e-6`,
or 21.54% of the fine peak. The fine scalar agrees with the old SSPRK2
pointwise extrapolant within `5.60381e-8`; this does not establish an
error enclosure or convergence order. Finite constraints and small
pre-support scalar leakage remain measured. Saved-array and source/profile
inspection precede a decision about further time resolution. The fixed
admission below records the calculation chosen before its outcome.

The six-stage local comparison below preserves the exact scalar balance
and moves the paired increment substantially closer to the independent
SSPRK2 asymptotic extrapolant. This supports a full-window temporal
comparison while holding the physical model and all spatial operators
fixed. It does not certify either local endpoint as exact.

| Item | Fixed SSPRK3 window |
| --- | --- |
| Inputs | Hash-authenticate the original six-case background record, its 66 bound inputs and the new local-step/operator evidence. Reuse its original 513 source/gauge rows and stored computational profile inputs. No old trajectory validation replay, profile-stage extraction, source solve or reference quadrature. |
| Model and domain | Same action, source, amplitudes `0` and `0.00012983127142551443`, quiet corner, spatial grids, wall-label gauge and original finite interval. The characteristic background is unchanged. |
| Cases | Zero/finite amplitude, each at 32 cells with 128 and 256 SSPRK3 steps. Supply exact cached midpoint nodes. Save 129 accepted slices per case: 768 steps and 2,308 coupled stages including four final reconstructions. |
| Operators | Existing first-order exterior departure derivative, order-two prefix quadrature, characteristic cap departures, order-four positive partial cap quadrature and unchanged hybrid wall equations. Only the time update changes. |
| Shared creation allowance | 300 process CPU seconds and 50,000,000 charged points for input authentication, profile construction, all four cases, diagnostics and creation validation together. Preserve the caller's accumulated budget. |
| Diagnostics | Raw and finite-minus-zero temporal histories in wall geometry, clocks, scalar, physical cap integrals and exterior fields. Compare the new cases with the existing 32-cell SSPRK2 256/512-step histories at common saved cuts. Retain raw/departure constraint diagnostics at 129/65 cuts, own-trace zero-cap pullback, outgoing cap history, profile margins and final wall compatibility defects. |
| Authentication | Reconstruct every saved physical state from raw departures, source, reference and profile; reconstruct four final coupled stages and all saved-array diagnostics. Verify hash bindings, shapes, clocks, accounting and finite domains. Do not replay a trajectory to check it. A standalone authentication has 60 CPU seconds/5,000,000 points. |
| Stop | First profile or exterior domain exit, nonfinite arithmetic, nonpositive cap metric, invalid hybrid incoming sum, non-timelike wall, stage Courant value above `.4`, or shared-budget exhaustion. Preserve the failed state and phase under ignored `.build`; diagnose before a new calculation. No automatic retries or changed physical input. |
| Decision | Determine whether temporal sensitivity decreases enough to separate it from the remaining spatial constraint and cap errors. Two SSPRK3 resolutions do not establish an observed order or a propagated error enclosure. Unchanged spatial deficiencies remain part of the result. |

Recorder/authentication controls have 30 process CPU seconds and 100,000
charged points. The physical comparison begins only after these controls
and independent review. Historical results and their producers stay fixed.

The recorder and independent authentication controls are complete: 27 tests
pass with 37,821 shared points and 1.796875 process CPU seconds. They cover
both exact time schedules, real raw-to-physical reconstruction, four final
stage checks, paired/RK2 comparison indexing, changed inputs and accounting,
and creation/standalone clocks beginning before input authentication. The
producer and its 80 input/source bindings remain frozen with the result.

## SSPRK3 with fixed physical inputs, 16 October 2026

**Outcome:** all 20 operator controls pass with 27,607 points and 0.5625
CPU seconds, plus two intentional exhaustion charges; the 17 local
producer controls pass with 3,808 points and 0.140625 CPU seconds.
The [physical local comparison](../results/development/nsc-atlas-ssprk3-step-audit-v1.json)
uses six stages, 70,118 points and 1.1875 CPU seconds; its standalone
check passes. Its paired discrepancy from the SSPRK2 asymptotic
extrapolant is `-2.71979e-9`, compared with `-2.67869e-8` for the two
SSPRK2 half steps. This is a controlled local improvement against that
comparison, not an exact-error result. The fixed admission follows.

The exact local ledger below identifies a temporal cancellation to improve.
Keep the current action, amplitude, source/gauge cache, characteristic
profile, cap reconstruction and exterior derivatives fixed. Borrow the
three-stage third-order method in [Gottlieb, Shu and Tadmor, equation 4.2](https://math.umd.edu/~tadmor/pub/linear-stability/Gottlieb-Shu-Tadmor.SIREV-01.pdf).
Apply it to the same raw departures at stage times `x`, `x+h`, `x+h/2`.
Its formal order does not establish stability or third-order convergence
of this spatially reconstructed, piecewise-profile coupled system.

Implement an explicit cached-midpoint interface: each stage reads the
existing source and reference data at its exact supplied key. Validate
the midpoint to floating-point accuracy without replacing it by a new
lookup coordinate. Restore physical geometry independently at every
stage and accepted endpoint. The original 513-node cache supports
128/256-step SSPRK3 comparisons; it does not supply the midpoints for
512 steps.

The operator controls have **30 process CPU seconds and 150,000 points**.
Use independent Runge--Kutta order identities, nonautonomous polynomial
forcing, nonlinear convergence, curved reference reconstruction, matched
profile cancellation, actual quiet-vacuum stages and endpoint/domain/CFL/
shared-budget failures. Preserve the `.4` stage Courant guard and the
first-failure context. Profile knots and limiter changes must remain
visible rather than assuming smooth convergence across them.

After controls and independent review, admit one local diagnostic at the
same saved cut 82: one SSPRK3 step for each amplitude, **six coupled
stages total**, using cache indices 328, 330 and 329 in stage order.
Its **30-process-CPU-second/150,000-point** budget includes hash
authentication, profile reconstruction and the weighted five-part ledger.
The diagnostic controls separately have 30 CPU seconds/100,000 points.
Authenticate and reuse the saved SSPRK2 local results without replaying
their twelve stages. Compare endpoints and paired increments against
their coarse, two-half-step and `(4*fine-coarse)/3` values. That last value
is an asymptotic extrapolant, not an independent exact solution or an
error certificate.

The exact ledger uses stage weights `(1/6,1/6,2/3)`. Retain all five
contributions, the common wall primitive and reference-geometry defects.
Stop on the first input, domain, arithmetic, CFL or budget failure; no
automatic retry, source interpolation, profile replacement or full
trajectory is admitted by this local calculation. One standalone check
may reproduce the same six stages. Its result determines whether a
subsequent fixed 128/256 trajectory comparison is useful; that comparison
needs its own finite shared budget and controls before it runs.

## Actual scalar-step balance, 16 October 2026

**Outcome:** the [six local steps](current-result.md#actual-scalar-step-balance-16-october-2026)
finish with twelve stages, 138,891 points and 1.546875 process CPU seconds.
Standalone checking reproduces the result. The paired local scalar change
is `8.03607e-8`, chiefly the difference of the background-chain contributions
after substantial cancellation in each amplitude. The five-part sum closes
within `5.11e-17`. All 22 controls pass with 5,312 points and 0.125 CPU
seconds. This locates the local contributions without attributing the
whole trajectory's error or changing the physical model. The admission
below records the fixed calculation chosen beforehand.

The saved-profile audit below exposes large cancelling terms. Its frozen
rates and straight paths cannot identify the error of an accepted coupled
step. The next calculation reuses the existing SSPRK2 update, conformal
field map and original 513-node source/gauge cache. It changes no physical
assumption or numerical evolution operator.

From the same saved starting state at cut 82, compare one original
256-grid interval with its two 512-grid halves, separately for the zero
and finite 32-cell/256-step cases. This is two local experiments, reset to
their saved initial states; it is not another trajectory campaign.
Use exactly twelve coupled stage evaluations in total. Accepted endpoints
need profile and frame evaluation, not another radial stage.

For each accepted step retain canonical and invariant field increments,
predictor and accepted geometry, and the exact decomposition into integrated
exterior flux, integrated returning-cap flux, nonlinear frame remainder,
reference-geometry remainder and background-chain remainder. The radius
and clock evolved by the code are departures from the reference gauge:
the accepted physical radius is not assumed equal to the trapezoid of
the total radius rate. Retain that difference explicitly. Sum the two
half-step ledgers before comparing them with the coarse ledger, then form
finite-minus-zero comparisons. Check the independent scalar balance and
every reconstructed endpoint.

Manufactured controls have a shared **30-process-CPU-second/100,000-point**
allowance. After those controls and independent review, the physical local
diagnostic has **30 process CPU seconds and 250,000 charged points**,
including hash authentication, reconstruction and all twelve stages.
A standalone check may repeat this same local calculation. No full history,
source solve, cache extension or altered profile is admitted. Preserve all
existing domain, positive-metric, timelike-wall and finite-arithmetic checks;
stop on the first failure with its stage and budget context. Diagnose any
failure before designing another calculation.

The decision is which terms carry the local one-step/two-half-step
difference and how they cancel. A local result does not establish the
propagated error over the retained window or justify removing a background
term from the evolution.

## Saved-profile time-integration diagnosis, 16 October 2026

The [diagnostic](../scripts/diagnose_nsc_atlas_background_time.py) is
implemented and independently reviewed. Its 35 controls pass with 10,956
points and 0.265625 process CPU seconds, plus the intentional three-point
exhaustion check. The [completed audit](../results/development/nsc-atlas-background-time-audit-v1.json)
uses 57,152 charged points and 1.0625 process CPU seconds. Its standalone
check passes. It authenticates all 66 input bindings without a source solve
or time evolution. The known reference-chain time difference reaches
`1.24476e-4`, while the actual zero-case scalar time difference reaches
`2.87645e-6`: substantial cancellation remains. At cut 82 the local paired
rate change is carried mainly by the exterior derivative, but it does not
explain the accumulated scalar difference. The following admission records
the fixed scope chosen before the audit.

The six-case result below removes dominant cap pulse smearing, but its
paired wall-scalar time change is `5.71668e-6`, or `15.83%` of the fine
peak. Radial quadrature is much smaller. The missing connection is to
separate integration of the known profile from changes in the physical
exterior/cap wall flux before selecting another temporal method.

Write the canonical wall primitive as `chi=B(T,R)+Z`, where `B=chi0`
and `Z=sum(departure Iq)`. The existing stage gives
`Z_x=J_chi-B_x`, with `B_x=T_x*(1-R^2)*p0+R_x*q0`. Therefore one
SSPRK2 step satisfies

$$
\chi_{n+1}-\chi_n=\tfrac h2(J_{\chi,n}+J_{\chi,*})
+B_{n+1}-B_n-\tfrac h2(B_{x,n}+B_{x,*}).
$$

The last three terms are the known-profile chain quadrature contribution.
They are not automatically an error in the method or the cause of the
observed scalar difference. The conformal wall rate also has the exact
split `J_c=Q_c+C_minus/(2*x_s)` from the existing hybrid junction map.
These identities supply the independent algebraic controls.

After manufactured controls and review, admit one saved-data audit with
**30 process CPU seconds and 150,000 charged points** for the whole audit:

- Authenticate the six-case payload, manifest and all 66 source/input
  bindings by hashes. Reuse the bound profile and original 513-node input
  cache; do not recursively replay scientific authentication or solve a source.
- On all 513 reference-gauge rows, evaluate `B,B_x` and compare exact
  `B(x)-B(x0)` with cumulative trapezoids on the existing 256/512 grids.
  Compare their difference at the 129 saved cuts with the observed zero
  scalar time difference. Separate signed contributions of steps crossing
  profile knots and steps wholly within a polynomial piece.
- Reconstruct exactly four frozen stages at saved cut 82: zero/finite
  amplitude, 32 cells, 256/512 steps. Retain `B_x`, `J_chi`, exterior `Q_c`,
  outgoing cap contribution, and the exact chain/split residuals. No stage
  advances the stored state.
- At each frozen `T,R,T_x,R_x`, evaluate the profile along the linear
  clock/radius path over one original 256-step interval and its two halves.
  Use analytic endpoint primitives and derivatives to measure local
  trapezoid defects and identify profile-knot crossings. This prescribed
  linear path is not the actual accelerated SSPRK2 trajectory.
- Preserve profile coverage, finite arithmetic, all existing stage domains
  and the shared budget. Stop on the first failure and diagnose it before
  any new calculation. No extrapolation, replacement source, time-cache
  extension, physical time evolution or automatic retry is admitted.

Manufactured polynomial, knot-partition, chain/split and accounting controls
have a separate 30-CPU-second/100,000-point allowance. A recorded check may
repeat this same bounded audit, never the trajectories. The decision is
whether known-profile integration can generate the observed scale, or
whether exterior scalar jets and wall coupling require the next correction;
the audit must preserve any unresolved cancellation or ambiguity.

## Coupled characteristic-background comparison, 16 October 2026

**Outcome:** all six fixed cases and creation authentication completed in
228.9375 process CPU seconds with 43,551,872 charged points. The
[saved result](current-result.md#characteristic-background-coupled-evolution-16-october-2026)
removes the large incoming cap transport error and substantially improves
paired wall geometry/clock agreement. Scalar time sensitivity and exterior
pointwise constraints remain unresolved. No source integration or change
of physical preparation was used. The fixed admission below records the
conditions chosen before that comparison.

The new profile, radial cap, input adapter and common stage have passed
72 manufactured/independent controls, including exact endpoint coverage.
The input adapter uses 8,129 points
plus 180 intentional budget-exhaustion charges; the coupled stage uses
26,661 points. The existing source, exterior and wall algorithms stay
fixed. The changed cap variables carry the known characteristic waveform
analytically, while departures retain the nonlinear response and every
incoming discrepancy. The next run answers whether this removes the
measured leading cap transport error in the original coupled preparation.

The following fixed comparison is admitted after its recorder and
authentication controls pass:

| Item | Characteristic-background comparison |
| --- | --- |
| Parent inputs | Authenticate the compatible record and its bound inputs. Reuse the original 513-point source and gauge caches. No new source integration or reference quadrature. |
| Computational profile | Reconstruct same-stage rates at the 129 saved zero-amplitude, 32-cell cuts. Select every second cut for 65 profile knots; retain the other 64 as a value/first-jet interpolation comparison. Define `F,F',F''` through the documented adapter. Reuse this one profile in every case. Its finite label interval is a query guard, not permission to extend the physical source. |
| Physical model and domain | Same action, compact preparation, `alpha=0` or `0.00012983127142551443`, original incoming interval, quiet corner and reference wall-label gauge. Actual radius, lapse and both clocks evolve. |
| Six cases | Each amplitude uses 16 cells/256 steps, 32 cells/256 steps and 32 cells/512 steps. Save 129 accepted slices per case. This adds the previously missing zero-amplitude time comparison. |
| Numerical methods | Unchanged first-order exterior departure derivative and order-two prefix quadrature; cap characteristic departures and positive partial quadrature of order four, split at cell and profile knots; fixed SSPRK2. No second-order exterior evolution is inferred from the frozen audit. |
| Shared allowance | 600 process CPU seconds and 100,000,000 charged points for parent authentication, profile extraction, all six evolutions, diagnostics and creation authentication together. The larger finite allowance covers profile-knot panels and positive partial quadrature; no child reset or uncharged source replacement. |
| Independent comparisons | Paired finite-minus-zero spatial and temporal responses; raw physical histories and departure constraints; cap incoming-characteristic pullback in the zero case; actual total cap integral restriction; predecessor overlap; profile coverage and interpolation sensitivity. |
| Radial quadrature control | At saved cuts 69, 79 and 128 of the finite 32-cell/256-step case, hold the state/profile fixed and compare cap orders four/eight. Retain wall mass coefficient, lapse, metric, outgoing trace and internal gamma consistency. This is a local quadrature comparison, not a time-evolution error enclosure. |
| Stop | First profile-domain exit, original exterior tube exit, nonfinite arithmetic, loss of positive cap metric or hybrid incoming sum, non-timelike wall, Courant value above `.4`, or cumulative budget failure. Preserve the first state/context and diagnose before designing another run. |
| Decision | Measure whether the leading pulse is transported accurately and whether the finite response gains paired space/time agreement. Remaining derivative, source and nonlinear residuals stay visible. This window still precedes central arrival and does not establish formation, stability, a global error bound or observational agreement. |

The variable number of radial panels makes stage charges depend on the
actual state. Record measured shared-budget deltas and validate their
arithmetic and ceilings. Authentication reconstructs saved physical
states, final radial/wall states and diagnostics without time evolution;
it does not claim to recompute every unsaved stage's charge.

## Characteristic background and frozen-stage transport diagnosis, 16 October 2026

The next connection is the common coupled stage. Its new input adapter
builds the computational profile from supplied zero-amplitude wall values
and same-stage rates: `zeta=T+atanh(R)` and
`F'=c_x/(T_x+R_x/(1-R^2))`. Finite differences of these supplied first
jets select interpolation curvature `F''`; they are not reported as a
measured physical second derivative. The resulting single quintic
interpolant supplies all runtime background jets. No source integration
or time evolution occurs inside this adapter.

The coupled adapter keeps the original exterior discretization, source,
wall equations and fixed SSPRK2, replacing only the cap representation.
It restores the full background wall rate and saves both raw departure
integrals and reconstructed physical integrals. Input-adapter controls
have 30 CPU seconds/20,000 points; common-stage controls have 30 CPU
seconds/100,000 points. They test manufactured input chains, actual moving
vacuum evolution, finite common-wall conversion, free incoming correction,
domain failure and cumulative-budget failure. A new original-source run
still requires its concrete profile coverage, cases and shared creation
budget to be specified after these controls.

The compatible record removes the vacuum and field-primitive defects but
still spreads the known incoming cap pulse. Reuse the exact de Sitter
reduced-wave formula and the same nonlinear radial constraints. Define one
finite C2 quintic Hermite profile `F` from explicit values and first/second
jets. Its analytic derivatives belong to that same interpolant. The
background `chi0=F(T+atanh r)-F(T-atanh r)` has exact physical-cell
integrals; evolve departures in those integrals and retain the full
finite incoming boundary discrepancy. The profile is a computational
subtraction, not an extra prescribed finite wall value.

Implement the background, conservative departure balances, and total-field
cap constraints with quadrature split at cell boundaries and mapped
profile knots. Include all baseline/departure cross terms in the energy
source. Operator controls use manufactured fields, known regular-axis
limits, exact incoming propagation, moving-wall chains and independent
radial integration; the profile suite has a cumulative 30-process-CPU-second
and 20,000-point limit. The first cap controls exposed negative nodal mass
from the Gauss antiderivative on a narrow positive-energy pulse. The
corrected cap suite uses direct positive-weight partial integrals, including
the independent order-4/8/16 comparison. Its additional partial quadrature
has a fixed 30-process-CPU-second/100,000-point allowance. No clipping or
physical run is used to repair the failed numerical representation.
A failed control identifies a correction before any
physical successor is designed. No new prepared-source evolution is
admitted here.

Separately, authenticate the compatible payload and its 51 source/input
hashes, without replaying its source or time solve. A frozen-data audit has
one 30-process-CPU-second/100,000-point allowance. It compares both saved
zero-case cap primitives with incoming-characteristic pullback of their own
wall histories, keeping linear/PCHIP interpolation sensitivity. At four
fixed finite 32-cell saved cuts (indices 16, 69, 79 and 128), cross first-
and second-order departure derivatives with prefix quadrature orders two
and four. Reconstruct the cap once per cut and hold it fixed in the four
crossed evaluations. Record outgoing jets, inferred incoming jets and wall
kinematic defects. The closed polynomial prefix control at 9/17/33 nodes
uses a separate 30-second/20,000-point operator allowance. Stop at the
first domain, finite-arithmetic or cumulative-budget failure.

This audit decides whether the exterior discrepancy chiefly changes with
the transverse derivative or prefix quadrature on the same state. It does
not measure the error in that state or in the frozen cap. The background
controls decide whether the cap can carry the known wave without numerical
spreading while leaving the nonlinear correction free. These decisions
precede another coupled run or a central-arrival extension.

**Outcome:** the frozen-state audit completes with 48,910 points and
0.859375 process CPU seconds. Its [record](../results/development/nsc-atlas-compatible-transport-audit-v1.json)
retains 51 authenticated input bindings, both cap pullbacks and the 16
crossed exterior evaluations. Geometric wall mismatches chiefly respond
to the transverse derivative choice; scalar mismatches do not uniformly
improve. This does not select a new exterior evolution scheme.

The characteristic profile and corrected cap operators pass 31 controls:
861 profile points and 41,332 cap points, using .046875 and .328125 process
CPU seconds respectively in the combined check. The narrow manufactured
pulse is smaller than a cap cell. Its wall-mass coefficient errors against
independent radial ODE integration are `.0291139`, `6.43643e-10` and
`5.57581e-11` at quadrature orders 4, 8 and 16. The corresponding internal
gamma consistency differences are `1.29328e-4`, `1.29928e-6` and
`3.28205e-14`. Both the exact matched background and a distinct nonzero
incoming correction are tested. These controls validate the new regional
representation; they are not a new result for the original preparation.

## Compatible reduced cap and vacuum departures, 16 October 2026

The [first coupled window](current-result.md#first-coupled-exterior-wall-and-cap-window-16-october-2026)
locates two numerical errors: a quiet vacuum layer beside the wall and a
cap field/gradient mismatch that corrupts the leading outgoing-wave
cancellation. The successor keeps the original action, source, amplitude,
incoming interval and wall-label gauge. It changes the numerical variables
and integrates known background terms analytically.

For the exterior, write the physical fields as exact Schwarzschild fields,
the unchanged source departure, and evolved spatial departures. Differentiate
the vacuum analytically and integrate the nonlinear prefix differences.
Evolve wall departures in radius, lapse and clocks. This preserves the exact
vacuum during time integration, while the original source and finite
responses remain dynamical. Retain raw independently differenced constraints
and the finite-difference residual of the exact background separately.

For the cap, use the canonical reduced field `chi=r*f`, `p=r*Pi` and
`q=chi_r=f+r*Phi`. Their exact equations are
`p_T=(B*q)_r-B_r*f-r*exp(delta)*U_alpha` and `q_T=(B*p)_r`.
The radial constraints give `B_r=exp(delta)*(2*m/r^2-r*P)`.
At zero amplitude, `m=r^3/2`, `P=3` and `U_0=2*f`, so the reaction
vanishes identically and the outgoing combination `p-q` obeys a homogeneous
transport equation. Reconstruct the field and its potential from the same
conservative gradient, with the regular centre and common wall value.
The canonical reduced field differs from the conformal `c` at finite
amplitude; the existing exact conversion and hybrid wall law remain in use.

First implement and independently verify these regional operators and
their common stage. Manufactured identities, regular-axis limits, analytic
incoming waves and exact vacuum controls have an aggregate ceiling of
30 process CPU seconds and 20,000 charged points per operator test suite.
Check source preservation, moving-grid fluxes, the actual wall-value chain
rule, and nonzero manufactured departures. Exact preservation of an imposed
identity alone is insufficient: compare propagated fields and outgoing
characteristics with independent analytic solutions. Preserve the first
run and its bound producers. No new physical trajectory is admitted by
these implementation tests; its fixed comparison and cumulative allowance
will be recorded after the operators and boundary data are checked.

Stop or revise a discrete representation when its axis, boundary data count
or conservation identities fail. The next decision is whether these
equivalent equations remove the diagnosed numerical errors with a measured
refinement improvement. Increasing the grid alone is not the proposed fix.

The four new operator suites now pass together: **94 controls** cover the
vacuum prefixes, reduced cap, coupled stage and independent departure
diagnostics. The manufactured cap propagation reduces its primitive error
from `.00311851` to `.00142199` on 16 and 32 cells. Its outgoing integral
defect is at most `2.53e-19`; its actual wall time-integration error,
`6.69e-6`, is retained separately. These are manufactured controls rather
than a new result for the prepared physical source. They admit the following
fixed successor comparison once the recorder's authentication controls pass.

The recorder's 30 synthetic checks now pass, including raw-to-physical
reconstruction, modified-data rejection, fixed-cache use and cumulative
creation/authentication limits. The bounded physical comparison below is
ready; no prepared-source trajectory was used as an implementation test.

| Item | Fixed compatible-evolution comparison |
| --- | --- |
| Inputs | Authenticate the saved first coupled record and reuse its 513 exact-index source and reference-gauge rows. No new source integration, reference quadrature, interpolation choice, amplitude or preparation. |
| Domain | The identical original `.4<=x<=x_last`, `0<=Ubar<=W_ref(x)` and `0<=r<=R(x)` window, quiet corner and clock origins. |
| State | Exact-vacuum exterior departures and four wall geometry/clock departures; cap `Ip,Iq` physical integrals. The common scalar wall value is derived from `sum(Iq)` and converted between canonical and conformal fields by the same action. |
| Cases | Zero and existing nonzero amplitude on 16 and 32 cells, each with 256 steps, plus the nonzero 32-cell case with 512 steps. Retain 129 accepted slices per case. Use the same first-order exterior upwinding, prefix order two, cap quadrature order four and fixed SSPRK2. The cap characteristic reconstruction and exact vacuum prefixes are the declared changes. |
| Cumulative budget | One shared 180-process-CPU-second and 8,000,000-point allowance, including authentication of the parent, every stage, diagnostics and creation authentication. No budget reset or uncharged replacement source. |
| Independent checks | Preserve raw saved-field pure constraints and add analytic-vacuum-plus-differenced-departure constraints. Retain source curvature defects, numerical wall-jet mismatch, final independent radial constraints, paired spatial and separate temporal comparisons, and the original rectangle overlap. Compare the predecessor's wall and exterior histories at identical coordinates. |
| Cap controls | Measure the actual outgoing reduced characteristic and independently tested propagation. Label primitive and global field-jump identities as consequences of the representation. Restrict fine conserved integrals by summing adjacent cells; do not average them as point values. |
| Stop | Same exterior tube (`.2<=r<=.8`, `-2<=log C<=2`, `-4<=c<=4`), positive cap `A,B` and hybrid incoming sum, timelike wall, maximum stage Courant `.4`, finite arithmetic and cumulative limits. Preserve the first failure and state; no retry, source replacement, enlarged tube or adaptive step. |
| Decision | Determine whether the vacuum layer and leading outgoing-wave error are removed in the coupled calculation, and whether the finite response improves under the declared spatial/time comparisons. A failed or unresolved comparison produces a diagnosis. A successful comparison supplies the next central-arrival-domain construction; it is not a uniform error enclosure, stability proof or observation of nature. |

**Outcome:** all five compatible cases finish at their original endpoint,
using 64.40625 process CPU seconds and 5,371,737 charged points including
creation authentication. The [record](../results/development/nsc-atlas-compatible-evolution-v1/compatible.json)
binds 51 source/input files. Independent standalone authentication passes
with 1,202,991 points and no source integration, reference quadrature or
time evolution. The vacuum geometric departures stay at roundoff and the
zero-amplitude outgoing reduced characteristic is exactly zero. Both
diagnosed predecessor errors are removed.

The nonzero response still needs a transport accuracy correction. The
paired spatial scalar change falls to `2.49582e-6`, about 6.9% of its fine
response maximum, but the separate finite-case time change `1.24653e-4`
exceeds that response. There is no zero-amplitude 512-step case, so
cancellation of the leading time error has not been measured. Nonzero
departure constraints improve spatially but remain largest beside the wall
during the pulse. The cap's nonlinear outgoing trace changes from
`1.60392e-4` to `3.44027e-4` with spatial refinement. A saved-history
incoming-characteristic comparison identifies substantial leading cap
pulse spreading, rather than an additional algebraic junction condition,
as the next connection. See the [numerical outcome](nsc-atlas-finite-amplitude.md#compatible-window-outcome-and-the-remaining-cap-transport).
No central-arrival extension is inferred from this finite-window result.

## Coupled characteristic exterior and polar cap, 16 October 2026

The completed rectangle supplies an independent overlap test for the next
coupled evolution. Its terminal leg is not the moving wall. Construct the
coupled numerical stage on the full exterior interval
`0<=Ubar<=W(x)`, where `W(.4)=Ubar_b>0`, and the regular cap
`0<=r<=r_wall(T)`. Moving coordinates keep both boundaries on grid endpoints
without starting a strip of zero width. The original source, central clock,
action and hybrid transmission law remain fixed.

For the exterior, the mixed equations determine the longitudinal first jets
from the source by ordered radial integrals. With `P=X_Ubar` supplied by
the same stage fields, first integrate `(r*r_x)_Ubar=-C/4`, then
`(r*f_x)_Ubar=-r_x*f_Ubar` for the normalized canonical scalar, then the
existing lapse equation. The moving-coordinate update includes its grid
advection. A source-plus-increment representation preserves the known compact
profile analytically. The wall's values and incoming derivatives must come
from the hybrid map at its actual events; discrepancy with the numerical
one-sided regional jets is a monitored coupling residual.

Implement and independently check the source snapshot, exterior prefix
operator and regular cap constraint/evolution operators before a physical
coupled run. A snapshot copies bounded initial data while the original
source budget is live; later polynomial evaluation debits the caller's
cumulative budget and performs no reintegration or clock reset. Derivatives
of an interpolated source must be derivatives of that same interpolant.
Any source-ODE target or constraint defect is reported separately.

Manufactured, analytic-vacuum and saved-data controls have a 30-process-CPU-
second and 20,000-evaluation ceiling for each operator's complete test run.
No new physical source integration, amplitude choice or spacetime run is
included at this implementation stage. The cap scheme must specify its axis,
radial constraints, one incoming scalar characteristic, moving-grid advection
and common-field compatibility. The exterior gauge must match the quiet
corner and use only the remaining wall-born null-label freedom. Resolve a
failed identity or boundary count before assembling the coupled step. Its
physical calculation will receive an explicit domain, shared budget and
refinement design after these connections are concrete.

The source snapshot, prefix and cap operators now pass their manufactured
and vacuum controls. The assembled simultaneous stage passes 14 controls,
including the actual wall acceleration, common-field chain rule, quiet
scalar, zero-amplitude geometric independence and fixed-step stopping.
Independent review confirms the interface data count. This admits the
following first coupled calculation; it does not assume its outcome.

| Item | Fixed coupled calculation |
| --- | --- |
| Borrowed mathematics and connection | The exact conformal action, original characteristic source, regular polar constraints and hybrid wall law derived in the [nonlinear owner](nsc-atlas-finite-amplitude.md). Connect these through one simultaneous exterior/cap/wall stage. |
| Domain | Original `x=.4..x_last`, exterior `0<=Ubar<=W_ref(x)`, cap `0<=r<=r_wall(x)`. The vacuum reference has wall radii `.45..55`; the sourced radius evolves. The original incoming labels and quiet corner remain fixed. |
| Gauge and state | Prescribe only `W_ref(x)`. Evolve the wall radius, common scalar, exterior log lapse, central clock and elapsed proper time. Source-plus-increment exterior state and weighted cap cells use the same stage. |
| Source and amplitude | Original compact pulse and action, `alpha=0` and `alpha=0.00012983127142551443`. Two source integrations retain `rtol=2e-10`, `atol=1e-12`, maximum step `source_length/32`; sample 257 points and the 33-node bottom while each child budget is live. A 129-knot snapshot retains the analytic pulse and has an independent withheld-point comparison against those source rows. |
| Spatial and temporal comparison | Exterior and cap both use 16 or 32 cells. Four zero/finite cases use 256 fixed SSPRK2 steps; a fifth finite 32-cell case uses 512 steps. Store 129 accepted slices per case. Exterior forward upwinding and the cap finite-volume reconstruction are first order; prefix quadrature is trapezoidal, cap quadrature order four. |
| Reference evaluation | Evaluate the reference wall gauge on the common 513-point time grid before the marches. Reuse exact-index rows; no interpolation of its derivatives. Keep inversion and quadrature accounting and endpoint defects. |
| Budget | One cumulative 180-process-CPU-second and 8,000,000-evaluation allowance for creation, including inputs, stages and checks. Each source child retains its own 30-second/20,000-evaluation limit and debits the parent. No child clock is reset. |
| Controls | Exact vacuum geometry and wall clocks; quiet field; finite-minus-zero spatial changes with raw errors retained; the separate 256/512-step comparison; source interpolation defects; independent pure null constraints on saved histories; one-sided wall compatibility, regular cap constraints and field compatibility. Compare the original rectangle only at common physical null coordinates, with interpolation explicitly identified. |
| Stop | First nonfinite state, declared exterior tube exit (`0.2<=r<=0.8`, `-2<=log C<=2`, `-4<=c<=4`), loss of positive cap `A` or incoming sum `H`, wall timelikeness, characteristic Courant value above `.4`, or budget failure. Preserve the case, stage and reason. Do not enlarge a tube, alter a step, amplitude or source, or retry silently. |
| Decision | A finite coupled numerical trajectory with measured changes and residuals, or an explicit incompatibility/numerical diagnosis. This interval does not reach the first central arrival; it cannot by itself complete the full observer response, formation/stability analysis or a uniform error enclosure. |

The five cases test one source law and one action. Zero amplitude retains
the normalized scalar equation while its stress vanishes. Comparison with
the source-only amplitude threshold does not turn that threshold into a
coupled existence bound. Completion of a recorded numerical trajectory
must be distinguished from its convergence and physical interpretation.

**Outcome:** all five cases finish on their declared grids, with 1,536
time steps, 3,077 simultaneous stages and no restart or guard failure.
Creation plus authentication uses 29.234375 process CPU seconds and
3,582,985 evaluations. The [source-bound record](../results/development/nsc-atlas-coupled-evolution-v1/coupled.json)
authenticates without source or time replay. Vacuum geometry and the
independent rectangle comparison improve under spatial refinement; the
separate time-step changes in geometric response are much smaller.
The maximum `E_Ubar` persists in an early quiet wall layer, cap
`Phi-f_r` compatibility is not resolved, and the scalar wall response
changes sign between spatial grids. The [measured results](current-result.md#first-coupled-exterior-wall-and-cap-window-16-october-2026)
retain these limitations. Diagnose the spatial/interface and cap
compatibility errors from these saved data before designing a successor;
this outcome does not admit an unbounded refinement campaign or a new
physical source.

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
