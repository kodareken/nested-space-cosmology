# Cosmological atlas: regions, boundaries and observer measurements

The programme asks how regions, their boundaries and their fields can form a
nested geometry, and what an observer inside one region can infer about the
rest. Spherical symmetry supplies the first tractable building block. The
paper's connecting argument is geometry, admissible attachment, specified
state, transmission, then a local measurement. Particle, cosmological and
quantum descriptions are questions to connect through explicit models; their
shared vocabulary alone does not identify their physical objects.

This is the programme owner. The [living paper](../paper/nested-space-cosmology.md)
presents its argument, the [inventory](research-inventory.md) records what has
been read and preserved, and the [existing admission rule](prior-art-and-open-claim.md)
controls new calculation. Those owners retain their separate roles. Reuse a
solved equation when available; identify whether the missing ingredient is its
source, state, domain, boundary data, or connection to another result before
starting a derivation.

## Current map and independent approaches

The [v1 motif](../results/development/nsc-cosmological-atlas-v1/atlas.json)
contains two potentially unequal children and a potential grandchild, with
finite descent depth two. Its separate Schwarzschild charts and closed-FLRW
child ansatz are historical starting assumptions. The [v2 parent](nsc-atlas-parent-initial-data.md)
adds shared Brill-Lindquist spatial data. The [v4 attachment](nsc-cosmological-atlas-v4.md)
places one sourced interior and a vacuum bridge in a published many-body
existence construction through selected smooth exterior ends. It does not
supply an explicit connector metric, a second established child or the common
future. The [v5 observer map](nsc-cosmological-atlas-v5.md) and
[v6 conditional waveform](nsc-cosmological-atlas-v6.md) use a separate isolated
spherical reference. The [finite exterior preparation](nsc-atlas-field-preparation.md)
supplies a conditional data-to-wall-to-observer operator in that reference.
The [v8 prepared response](nsc-atlas-prepared-response.md) evaluates it for the
independent compact input, with an analytically positive post-input tail and
numerical refinement checks over a finite arrival window. These causal and field results
do not automatically extend to the glued many-body spacetime.

The [conformal stress check](nsc-atlas-scalar-wall-consistency.md) now fixes
the field normalization and checks its quadratic surface stress and
tangential flux at that same wall. Both scalar corrections vanish, so no
extra material surface source is needed for those conditions. The sourced
bulk evolution, normal-force balance and stability remain open. The
[characteristic initialization](nsc-atlas-characteristic-initialization.md)
now fixes one initial null-leg radius correction from the same compact
probe with the corner geometry and affine normalization held fixed. Its
post-input radial slope differs from the reference. This is not a mass
change or an evolved wall response by itself. The
[source-leg mass calculation](nsc-atlas-source-mass.md) supplies the
transverse-radius equation and geometric mass on the same null surface.
Full frame conversion and improved-stress transport give the same
negative constant coefficient after the compact input. Later bulk and
wall evolution remain open. The [normal-force successor](nsc-atlas-scalar-normal-force.md)
now establishes zero quadratic regular normal-stress jump at the same
conformal wall. Fixed vacuum/tension sources then require zero variation
of mean extrinsic-curvature trace, including metric and wall displacement.
The [wall-mass closure](nsc-atlas-wall-mass-closure.md) now supplies an
analytic exterior-side geometric-mass functional from the fixed leading
wall trace. It retains both null fluxes and the full frame correction.
Numerical evaluation and coupled cap/wall geometry are still unresolved;
the wall tail does not inherit the initial leg's mass plateau or late sign.

Keep three approaches available:

- **Geometric atlas:** assemble compatible published regions and junctions;
  keep topology, source and initial-data assumptions on each attachment.
- **Exterior-response atlas:** reuse the [finite response and memory machinery](nsc-evolving-reduction.md)
  with a physically justified Hamiltonian and state. An atlas of locally
  indistinguishable exteriors is useful when reconstruction is nonunique.
- **Connected mouths:** retain the charged semiclassical construction in the
  [arXiv library](arxiv-studies.md) as a separate communication route. Its
  sources and causal connections cannot be imported into the neutral scalar
  reference without a compatibility argument.

One-hole multiple daughters and separate holes producing different children
are different branching hypotheses. The archived conformal cyclic cosmology
citation, arXiv:2503.24263, remains a separate comparison question in the
[library](arxiv-studies.md), with no matching equations absorbed by this
programme update. Review its conformal matching and state assumptions before
connecting it to any black-hole branch or energy-recycling claim. None of
these alternatives becomes compulsory when another one fails.

<a id="programme-dependencies"></a>

## Research dependencies and useful decisions

The stable calculation IDs below identify missing connections, not scheduled
campaigns. A calculation earns its place through **prerequisites -> decision
-> further calculation unlocked**. A finite obstruction is valuable when it
rules out a proposed connection under fixed assumptions. No method can promise
that every open branch succeeds.

| Connection ID | Prerequisites and decision | Further calculation unlocked |
|---|---|---|
| `geometry-baseline` | Fix a borrowed region, topology and validity domain; decide which overlaps can be compared. | `boundary-classification`, alternative region comparisons. |
| `boundary-classification` | Specify whether the boundary is material, causal, effective or a geometric junction; check the appropriate source and matching conditions. | `wall-transmission`, `clock-entropy`, `field-orientation`. |
| `wall-transmission` | Supply a field/action, wall source and state class; decide the admissible transmission law. The v6 probe law is already available. | `wall-trace` after compatible preparation data exist. |
| `parent-preparation` | Locate a data surface in the retained causal domain and specify admissible exterior-side field data; decide whether they can supply the wall input. The finite region-IV problem now has a conditional existence certificate; remote generation remains open. | `wall-trace`, or a recorded exclusion that redirects the source route. |
| `wall-trace` | Combine compatible preparation, propagation and wall transmission; determine the incident trace rather than assuming its origin. The region-IV trace is evaluated for the fixed independent input in v8. | `observer-response` for the independently specified finite exterior data. |
| `observer-response` | Use the observer, clocks and fixed transfer operator with a justified trace; v8 resolves a finite central waveform and a positive post-input field tail. | `observer-inference`, comparisons of sources and attachments. |
| `observer-inference` | The [fixed-model inverse](nsc-atlas-observer-inference.md) recovers a compatible preparation from a complete signed central trace and gives positive compact-source exclusions. A second geometry and a physical preparation law remain unresolved. | A preparation-restricted response-range comparison or a source-linked independent channel. |
| `probe-stress` | Fix physical field normalization and sufficient regularity; borrow the full Hilbert stress and coupled shell laws. The [conformal shell check](nsc-atlas-scalar-wall-consistency.md) establishes zero additional surface stress and tangential flux jump for this probe. | `gravitational-initialization`, `wall-normal-force`, a later preparation-energy criterion; no amplitude bound is supplied by surface cancellation. |
| `gravitational-initialization` | Specify compatible corner geometry and a null normalization. The [source-leg constraint](nsc-atlas-characteristic-initialization.md) fixes the quadratic radius correction for the original compact field trace; the zero bottom leg retains its reference geometry. | `source-leg-mass` after the transverse radius equation is supplied; no mass shift is inferred from a radial slope alone. |
| `source-leg-mass` | The [mixed equation and mass law](nsc-atlas-source-mass.md) transport the fixed corner datum and determine geometric mass on the same initial leg. Full frame conversion, geometric and improved-stress checks agree. | `wall-mass-transport`; a finite initial-surface mass is not an asymptotic energy or evolved wall response. |
| `wall-normal-force` | The [full normal projection and Gauss identity](nsc-atlas-scalar-normal-force.md) establish zero additional normal-stress jump for the conformal probe. Fixed vacuum/tension sources require zero variation of mean curvature. | `bulk-backreaction` with this necessary boundary condition; it includes metric and wall displacement rather than setting displacement to zero. |
| `wall-mass-transport` | The [analytic closure](nsc-atlas-wall-mass-closure.md) uses the fixed leading wall trace, its first derivative, quiet mass anchor and reference geometry. Both null fluxes plus the full frame term agree with original-frame improved-stress transport. | `wall-mass-evaluation` or `cap-wall-closure`; the quadratic exterior functional is closed, without numerical mass values or a solved coupled boundary. |
| `wall-mass-evaluation` | Decide whether a finite mass value resolves the next junction uncertainty. Specify quadrature and independent field/geometry error control; the waveform's old samples do not certify the mass integral. | A controlled numerical exterior input to `cap-wall-closure`, if that value changes the matching decision. |
| `cap-wall-closure` | Identify the remaining cap-side mass and geometric junction data needed alongside the exterior functional and local stress conditions. Retain the cap vacuum and displacement effects on its nonconstant background mass. | A compatible finite metric/wall response or a precise missing datum; no new bulk solver is presumed. |
| `bulk-backreaction` | Supply any propagation equations and remaining junction conditions actually required after the local closures. Initial null geometry, mass functionals and the normal-force condition do not themselves prove a sourced bulk future. | A finite validity estimate or a sourced metric/wall response, then an independently shared gravitational response channel. |
| `inheritance-map` | Specify a physical transition/state map compatible with the regional laws; decide what is inherited rather than merely relabelled. | `regional-scales`, `branching-recurrence`. |
| `regional-scales` | Fix regional rulers and dimensionless observables; decide whether a change is units, state, dressing or a different law. | Conditional coupling or propagation comparisons. |
| `field-orientation` | Specify field lifts, conserved current and attachment orientation; decide which charge identifications are allowed. | A charge-transport or asymmetry calculation with a stated initial state. |
| `clock-entropy` | Supply proper clocks, clock domains and a state evolution; decide which durations and entropy changes are comparable. | A finite arrow-of-time discriminator, without a cross-universe age assumption. |
| `branching-recurrence` | Supply an offspring rule, physical transition and state/energy bookkeeping; decide what a finite descendant ensemble predicts. | Selection or recurrence comparisons with their probability assumptions exposed. |

The main response dependency chain is
`parent-preparation + wall-transmission -> wall-trace -> observer-response -> observer-inference`.
The v6 waveform establishes `observer-response` conditional on a prescribed
trace. The [finite preparation certificate](nsc-atlas-field-preparation.md)
connects independent exterior-side data to a unique trace and response operator;
the [v8 evaluation](nsc-atlas-prepared-response.md) now closes the numerical
trace and response for those finite data. Observer identifiability under
preparation freedom and the remote origin of the input remain distinct
questions. Inheritance, charge, clock and alternative-geometry questions
may proceed independently when their
own prerequisites are met. A different exterior model needs its own dependency
edges, not a forced passage through this scalar example.

The parallel physical-response branch is
`wall-transmission + parent-preparation + physical normalization -> probe-stress -> gravitational-initialization -> source-leg-mass`.
The local branch `probe-stress -> wall-normal-force` needs no mass quadrature.
The existing `wall-trace`, canonical normalization and quiet mass anchor
now close `wall-mass-transport` analytically; no numerical source-leg mass
curve or normal-force result is needed for the identity itself. The next
decision is `wall-mass-evaluation` or a named missing `cap-wall-closure`
input. These branches can supply inputs to `bulk-backreaction` when needed.
They supply prerequisites for interpreting the field as a gravitational source;
they do not supply another proof of the waveform or a formation history.

## Eight linked research questions

These themes organize the existing research census; they do not declare a new
audit or eight solved connections. The stable IDs used by the programme record
are `sphere`, `insideoutside`, `quantum`, `barriers`, `scales`, `charge`, `time`
and `branching`.

<a id="question-spheres"></a>

### Spheres

**Owner and current result:** [parent initial data](nsc-atlas-parent-initial-data.md)
and [sourced slice](nsc-atlas-connection-slice.md); the
[v4 record](../results/development/nsc-cosmological-atlas-v4/atlas.json)
combines locally spherical bodies through an existence construction.
**Borrowed foundation:** spherical symmetry and the Euclidean fixed-volume,
isotropic surface-energy result recorded in the [Wulff reading card](arxiv-studies.md).
**Missing connection:** a specified dynamical selection or stability argument
for the relevant material or gravitational system. **Next question:** which
departure from symmetry changes the attachment or an observer measurement,
and what source permits it? **Knowledge boundary:** a sphere's simplicity is
motivation; neither two competing forces nor the cited Euclidean minimization
proves a universal most-probable cosmological shape. S2 orbits do not fix global
S3 topology, and a multi-hole parent is not globally spherical.

<a id="question-inside-outside"></a>

### Inside and outside

**Owner and current result:** [connection gluing](nsc-atlas-connection-gluing.md),
[causal classification](nsc-atlas-causal-classification.md) and
[observer access](../results/development/nsc-cosmological-atlas-v5/observer-access.json).
The [prepared response](nsc-atlas-prepared-response.md) adds a resolved finite
scalar waveform and tail. The sourced initial slice and isolated reference access verdicts have different
domains. **Borrowed foundation:** Israel matching, published thin-wall causal
completions and vacuum end gluing. **Missing connection:** globally compatible
preparation/evolution data beyond the isolated reference and protected domains.
**Next question:** what exterior properties can that observer identify when
preparation freedom is included, and what wider history could supply the input?
**Knowledge boundary:** geometric attachment, descent, causal
influence and correlation remain distinct; the v4 theorem supplies no generic
sibling signal channel or many-body future.

<a id="question-particle-wave"></a>

### Particle and wave descriptions

**Owner and current result:** [Dirac tetrads](nsc-dirac-tetrad.md),
[local-observer correspondence](nsc-local-observer-correspondence.md) and
[scalar transfer](nsc-atlas-field-transfer.md). The
[v6 field record](../results/development/nsc-cosmological-atlas-v6/field-transfer.json)
is a conditional classical response; the [v8 evaluation](nsc-atlas-prepared-response.md)
propagates independent exterior data to that observer. Finite fermionic observer constructions
belong to their separately specified operator/state models.
**Borrowed foundation:** field equations, mode decompositions, regular-centre
conditions and finite retained-region reduction. **Missing connection:** a
physical identification of a localized particle sector, its quantum state and
the proposed regional geometry. **Next question:** which shared observable and
state map would make a field/mode description comparable with a localized
regional excitation? **Knowledge boundary:** a classical scalar pulse does not
quantize the atlas, and wave behaviour does not establish literal nested
particle interiors. No universal particle-spectrum calculation is required
before answering the finite signal question.

<a id="question-gradient-boundaries"></a>

### Gradient zones and boundaries

**Owner and current result:** [junction sources](nsc-atlas-junction-sources.md),
[moving boundaries](nsc-atlas-moving-boundaries.md) and
[field matching](nsc-atlas-field-matching.md); the
[v6 bundle](../results/development/nsc-cosmological-atlas-v6/narrative.md)
uses a specified pure-tension wall and scalar jump.
**Borrowed foundation:** junction conservation, characteristics and
distributional field matching. **Missing connection:** the source and transport
law corresponding to any newly proposed energy-gradient zone. **Next question:**
does that zone act as a material interface, a finite field-propagation barrier,
an observer horizon, or an attachment of charts, and what observable separates
these possibilities? **Knowledge boundary:** those categories are not
interchangeable. The [RN equality crossings](nsc-rn-gradient-bubble.md) remain
diagnostics; a pressure crossing alone establishes neither mechanical balance
nor centre regularity.

<a id="question-regional-scales"></a>

### Regional scales and constants

**Owner and current result:** [finite inheritance](nsc-nested-qualities.md),
[scale binding](nsc-scale-binding.md) and [spectral endpoint](nsc-spectral-endpoint.md).
The atlas's [unit controls](../tests/test_nsc_cosmological_atlas.py) separate
rescaling from dimensionless changes; the spectral and charged-throat results
retain their own models. **Borrowed foundation:** dimensional analysis,
normalized finite operators, known spectra and the imported charged-throat
scale relation. **Missing connection:** a physical transition that selects
dimensionless parameters or an additional propagation coupling.
**Next question:** which dimensionless observable could distinguish inherited
laws, changed state and environmental dressing in a specified attachment?
**Knowledge boundary:** the finite inheritance identity does not derive
regional constants. Planck length is a derived scale; a shortest allowed
wavelength requires another hypothesis. Absolute zero is a thermodynamic
state limit, not the other end of a wavelength interval. Coordinate light
delays alone do not change locally measured vacuum light speed.

<a id="question-charge"></a>

### Charge and orientation

**Owner and current result:** [charged-source route](nsc-charged-self-sourcing-route.md)
and the exploratory [swap analysis](child-neck/swap-symmetry-note.md), compared
with the [CPT and black-mirror reading cards](arxiv-studies.md).
They retain particular field, topology and symmetry assumptions; the active
v6 real scalar is neutral. **Borrowed foundation:** conserved field currents,
charge conjugation and the cited CPT/boundary constructions.
**Missing connection:** an attachment-selected field identification and
admissible state that transports or generates a charge asymmetry.
**Next question:** which transformation preserves the chosen action and
current across the actual boundary? **Knowledge boundary:** reversing an
inside/outside orientation does not by itself turn matter into antimatter.
The swap note's algebra and memory-sourced comparison table have not been
independently reverified by this programme inventory.

<a id="question-time"></a>

### Time and its direction

**Owner and current result:** [causal observers](nsc-atlas-causal-observers.md),
[clock inheritance](nsc-clock-horizon.md) and
[the access table](../results/development/nsc-cosmological-atlas-v5/observer-access.json).
The reference assigns separate centre, wall and parent proper clocks with
explicit origins and finite access domains. **Borrowed foundation:** proper
time, null propagation, conformal coordinates and finite state evolution.
**Missing connection:** a stated entropy/coarse-graining choice and physical
state dynamics across an attachment. **Next question:** what clock-domain
and entropy measurements would distinguish an inherited time orientation
from a chosen low-entropy initial condition? **Knowledge boundary:** a
coordinate divergence is not an infinite elapsed lifetime, and finite clock
maps establish neither a universal age nor an inevitable arrow. The retained
w2 clock-factor discrepancy remains in the [inventory](research-inventory.md).

<a id="question-branching"></a>

### Branching, inheritance and recycling

**Owner and current result:** the [finite atlas](../results/development/nsc-cosmological-atlas-v1/atlas.json),
[finite inheritance theorem](nsc-nested-qualities.md) and
[branching literature cards](arxiv-studies.md). The graph permits unequal
children and finite descendants; it does not determine offspring rates.
**Borrowed foundation:** Barrabes-Frolov multiple-daughter geometry and
Altenberg multitype offspring/inheritance mathematics, in their recorded
reading domains. **Missing connection:** formation, a transition kernel,
admissible state transfer and regional energy bookkeeping.
**Next question:** which physical transition could supply an offspring or
inheritance rule, and which finite population observable would discriminate
it? **Knowledge boundary:** recycling and recurrence remain hypotheses.
Existence of observers supplies no numerical probability proof; conservation
does not by itself imply repeated personal identity or infinite cosmic
recurrence. A finite graph or preset matter fraction cannot establish them.

<a id="signal-origin-cycle"></a>

## Signal-origin cycle: certified preparation and measured response

The existing [v5 access table](../results/development/nsc-cosmological-atlas-v5/observer-access.json)
already excludes the right static parent observer P as a source reaching C
through this attachment: future-directed causal V cannot decrease from P's
V>0 to the wall's V<0. Repeating that exclusion is not the next discovery.
The [v6 transfer](nsc-atlas-field-transfer.md) already computes a centre
waveform from prescribed exterior-side wall data. Its uncomputed origin is
the missing connection between these two results.

The [finite preparation derivation](nsc-atlas-field-preparation.md) now closes a
conditional data-to-wall-to-observer operator on the early region-IV side of
the retained wall. It keeps the isolated geometry, neutral conformal scalar
and fixed matching law. Independent characteristic data are prescribed on
U=U0, with zero bottom data and zero cap Cauchy data at t_a=t_W(R/ell=0.45).
The retained wedge has 0.4 <= R/ell <= 0.55 and never reaches a horizon or
singularity. An analytic no-return bound inherits the outgoing wall condition;
a supremum-norm contraction bound q_* = 0.138820992 < 1 proves a unique finite
solution and determines its wall trace. The existing cap transfer then defines
the central response operator. Borrowed field and contraction mathematics are
being applied to these specified data and this wedge, not rediscovered.

The v7 certificate is a derived conditional existence and transfer result.
Its input is not chosen to reproduce the v6 wall trace. Comparing with v6's
earlier zero cap at R/ell=0.05 requires the additional
no-earlier-incoming-history assumption stated in the derivation. "Exterior
side" means the retained Schwarzschild side of the material wall, not the
right asymptotically flat exterior I. The known P exclusion remains unchanged.

The [v8 evaluation](nsc-atlas-prepared-response.md) now computes the operator
and composes it with the existing observer transfer. A bipolar central field
is followed by a positive tail over the admitted post-input arrival interval
-0.340772595<=T/ell<=-0.250230701. Positivity has an independent analytic
proof; the measured waveform passes spatial and geometry refinement and
three positive-integral bracket checks. The first unsegmented prototype
stopped at its call cap; event-aligned integration succeeded with the same
theory, source and tolerances. These are numerical measurements supported by
analytic bounds, not rigorous enclosures of every waveform value or a claim
of persistence after the window.

This closes `wall-trace` and `observer-response` for this preparation. Remote generation of those
data, a globally prepared quantum state, regular wall formation and
backreaction remain separate connections. Alternative attachments,
connected-mouth models and response-equivalent exteriors stay available.
Each further calculation still requires its own finite domain, decision,
budget and redirect outcome under the [existing admission rule](prior-art-and-open-claim.md).

<a id="observer-identifiability"></a>

### Observer identifiability under preparation freedom

The next useful question is whether an interior measurement distinguishes
exterior structure once the allowed input preparations are stated. The v8
calculation fixes one exterior and one source; it does not establish that its
tail uniquely identifies that geometry. A change in source could imitate a
change in propagation. The v6 prescribed-wall pulse is a calibration, and
the v8 zero-kernel control is a numerical omission, so neither automatically
supplies a competing physical exterior.

The [analytic inverse](nsc-atlas-observer-inference.md) now closes a bounded
part of this question. With the same geometry, known clock, zero cap and
no-return window, a complete compatible signed central trace uniquely
determines its characteristic preparation. C2 source and wall histories
that vanish initially correspond to C1 observer histories that vanish on
the corresponding clock interval. The inverse has a uniform amplitude
bound and an explicit corner check. No new wave solve is required.

That enlarged class does not preserve the original source's positivity,
late cutoff or shape. Nonzero positive sources ending at the prescribed
cutoff must give a positive post-input field and negative accumulated
signed field. An entirely positive observer bump ending before that cutoff
requires a source that continues afterwards, and is excluded by the
compact-source class. The [record](../results/development/nsc-atlas-observer-inference-v1/inference.json)
states the distinction without claiming equivalence of two geometries.

The next connection is therefore an independently justified preparation
restriction and its response-range comparison, or an additional channel
with a shared physical source law. Simply adding a detector in the same
unchanged cap can be redundant once the full wall history is fixed.
The source class, common clocks and measured quantity must precede that
comparison; its budget and stop rule need their own admission. The
geometric atlas, connected mouths, inheritance and clock routes remain
available under their own prerequisites.

The selected next step was to check the existing probe's own stress before
introducing a new source or channel. The [surface check](nsc-atlas-scalar-wall-consistency.md)
finds compatibility with the fixed material wall, using its canonical
normalization and original C3 preparation. It does not yet restrict the
source amplitude or produce a second observable. The
[gravitational initialization](nsc-atlas-characteristic-initialization.md)
then applies the borrowed null constraint to the same compact field trace.
The specified corner geometry and affine gauge force a radius correction.
A generic bulk-stress bound was deferred because it would not provide
these missing gravitational data. The [mass successor](nsc-atlas-source-mass.md)
then closes the mixed-radius equation on the same surface and supplies a
coordinate-invariant geometric response. Its conformal-frame difference
must be retained during the input. The [normal-force check](nsc-atlas-scalar-normal-force.md)
now closes the local stress jump without computing individual pressure
curves. The required mean-curvature variation combines metric and wall
motion. The [wall-mass closure](nsc-atlas-wall-mass-closure.md) now provides
the exterior-side analytic functional directly from those leading traces,
using both null fluxes and the full conformal endpoint term. Its first
derivatives suffice for a future evaluation; no sampled waveform is
differentiated in the proof. The next choice is controlled evaluation if
a mass value resolves the matching question, or the named missing cap/wall
input first. The cap vacuum and radius-dependent reference mass prevent
silently extending the exterior shortcut through the wall. Coupled
evolution, stability and finite-amplitude validity remain open.

## Owners and reproducible output

- [Exterior wall-mass closure](nsc-atlas-wall-mass-closure.md) and its
  [exact-algebra recorder](../scripts/record_nsc_atlas_wall_mass_closure.py)
  bind the conditional finite functional to known leading field traces.
  Numerical evaluation and remaining cap/wall junction data stay separate.

- [Normal-force compatibility](nsc-atlas-scalar-normal-force.md) and its
  [exact-algebra recorder](../scripts/record_nsc_atlas_scalar_normal_force.py)
  supply the remaining local scalar stress jump and a necessary geometric
  wall condition. They leave individual pressures and coupled evolution open.

- [Source-leg mass](nsc-atlas-source-mass.md), its
  [adapter](../src/recursive_horizons/nsc_atlas_source_mass.py) and
  [recorder](../scripts/record_nsc_atlas_source_mass.py) connect the fixed initial
  preparation to an invariant geometric response through borrowed mass
  laws, mixed-radius transport and full frame conversion. Later bulk and
  wall evolution remain unresolved.

- [Characteristic initialization](nsc-atlas-characteristic-initialization.md),
  [constraint adapter](../src/recursive_horizons/nsc_atlas_characteristic_initialization.py)
  and [recorder](../scripts/record_nsc_atlas_characteristic_initialization.py)
  connect the original preparation to its intrinsic gravitational source-leg
  constraint. Coupled evolution and quantitative amplitude validity remain open.

- [Scalar wall consistency](nsc-atlas-scalar-wall-consistency.md) and its
  [exact-algebra recorder](../scripts/record_nsc_atlas_scalar_wall_consistency.py)
  separate two checked surface conditions from open bulk backreaction.

- [Observer-inverse proof](nsc-atlas-observer-inference.md) and
  [declaration checker](../scripts/record_nsc_atlas_observer_inference.py)
  bind the conditional analytic result to unchanged v7/v8 evidence, with
  no numerical reconstruction or second-model equivalence claimed.

- [Model](../src/recursive_horizons/nsc_cosmological_atlas.py): region/link
  records, finite validation, imported junction/clock constraints, fixed-theory
  check, causal graph bounds, and explicitly supplied Hamiltonian adapter.
- [Producer](../scripts/derive_nsc_cosmological_atlas.py): creation-only JSON,
  SVG patch maps, PNG previews through the existing matplotlib, illustrated
  narrative and portable source/output hashes. Previews use the SVG's geometry
  and text, with a flat equivalent of its decorative background gradient.
- [Tests](../tests/test_nsc_cosmological_atlas.py): incompatible junctions,
  unequal children, clock exits, causal delays/disconnection, unit changes,
  finite depth/cycles, fixed-theory perturbations and invisible exteriors.
- [Illustrated narrative](../results/development/nsc-cosmological-atlas-v1/narrative.md)
  and [machine-readable atlas](../results/development/nsc-cosmological-atlas-v1/atlas.json).
- [Prepared-response owner](nsc-atlas-prepared-response.md),
  [numerical adapter](../src/recursive_horizons/nsc_atlas_prepared_response.py),
  [v8 producer](../scripts/derive_nsc_cosmological_atlas_v8.py) and
  [independent controls](../tests/test_nsc_atlas_prepared_response.py) own the
  latest finite response and its analytic tail checks. The
  [v8 response record](../results/development/nsc-cosmological-atlas-v8/prepared-response.json)
  keeps its samples, numerical diagnostics and finite compute counts.

The [reproduction guide](reproducing.md) owns the current root commands,
Windows interpreter fallback and separate historical renderer profile. Use
`scripts/replay_history.py` in that explicitly selected profile for the frozen
v1-v6 bundles; v7 has its separate explicit `atlas-v7` replay command and
recorded profile. Current producers intentionally have different source
bindings. They may create a fresh bundle at a supplied `--output` path, but cannot make an
old bundle pass by silently loading old sources. The living paper's overview
is supplied by [the overview producer](../scripts/build_geometry_overview.py).

The producer refuses an existing output directory. A failure after claiming
the directory leaves an inspectable partial bundle, without a success manifest;
use a new versioned output directory after correction. Verification is read-only.
Source hashes include transitive local Python imports, the producer, test,
this note and local notes used as evidence. Paths are repository-relative;
no host path, clock timestamp or
unverified historical revision supplies authority. Byte hashes authenticate
outputs; numerical replay uses relative 1e-10 and absolute 1e-12 tolerances.

## Interfaces and physical meaning

`build_atlas()` returns independent scenario records. Every link has one of
four kinds, a source and an evidence statement: descent, attachment, causal,
or correlation. An unresolved link stays unresolved. The validated descent
graph is acyclic and bounded. Explicit `contained_in` interfaces retain their
host's declared generation, so moving through an interface cannot reset the
lineage count. Arbitrary attachments and causal networks may have cycles.

`angular_junction()` checks equal induced areal radius and real timelike
normal, and evaluates the imported angular Israel condition. It does not
solve temporal extrinsic curvature, shell matter or conservation. The
small-neck orientation exposes its negative surface-density requirement.
`null_collision_residual()` is only the imported collision-sphere condition.

`static_clock()` is a static Schwarzschild observer's elapsed proper time
within a finite declared interval. It rejects horizons and domain exits.
`earliest_arrival()` composes only explicitly declared positive-delay causal
edges sharing a named clock and finite interval. It gives a conditional graph
bound, not a geodesic solution. Physical candidate links presently lack the
complete numeric clock maps, so they do not acquire numerical arrival times.
Synthetic delay fixtures test this distinction. Entanglement, common ancestry
and correlation never become signal routes by themselves.

`require_fixed_theory()` permits state changes while rejecting changes to a
candidate's theory data and declared regional bindings of fixed couplings.
This is a record-level guard, not proof that a given state perturbation
satisfies its field equations.

`evolve_supplied_hamiltonian()` accepts a finite Hermitian matrix or callable
`H(t)` (dimension inferred at the first finite sample and validated on every
call), uniform finite time grid, retained indices, covariance and explicit origin. It reuses
`nsc_evolving_reduction.evolve_retained_region`; it neither infers coefficients
from the diagram nor derives a relativistic signal cone. The saved four-mode
integration is an adapter fixture. Its unequal occupations exercise different
children without assigning cosmological values to the coefficients.

The response-equivalence control uses known Schur algebra: changing an exactly
uncoupled exterior mode changes the exterior spectrum but not the local
resolvent. The covariance control uses the same visible state. This supplies
a concrete limit on exterior reconstruction, not a new inverse theorem. Its
resolvent lives in the upper half-plane, with poles excluded in boundary
values. A direct hidden-mode spectral probe would distinguish the examples;
adding a coupling is a separate model and changes the local response problem.

## Borrowed sources and remaining connections

The complete versioned URLs and reading boundaries live in the generated
atlas source records. Key anchors are [Barrabes–Frolov](https://arxiv.org/html/hep-th/9511136v1)
for one-hole branching and null shells, [Emparan et al.](https://arxiv.org/html/2012.07821v2)
for a controlled connected-mouth alternative, and [Poplawski v2](https://arxiv.org/html/2007.11556v2)
for the anisotropic spin-fluid hypothesis and its explicit geometric limitation.
The exact Kantowski–Sachs metric is not a finished route to the required S3
geometry. The connected-mouth source domain cannot be silently replaced with
the exploratory one-Dirac-field neutral-child source.
The one-hole l>a disconnection threshold specifically uses the
near-singularity approximation, equal nucleation radius, and a<<R0<<2m in
equations 9-12; it is not a universal separation rule for sibling universes.

Beyond the v4 existence application and the isolated v5-v8 reference, missing
work includes an explicit common connector/evolution, compatible interface
sources, an admissible joint state, interpatch clocks and exterior responses
outside the computed finite IV domain. Dimensionless
inherited laws, environmental dressing and state choices must stay distinct.
Only after those connections exist does a local dark-sector, antimatter or
time discriminator become a meaningful target. There is no abundance fit here.

## Admission and limits

The source, domain and state compatibility map is the admission category.
Each calculation records six fields: question, borrowed mathematics, missing
connection, decision and unlock, finite budget, and stop condition.
Junction identities, clocks, finite Schur elimination and the numerical
adapter are imported mathematics or software controls. Neither their replay
nor their tests are a new physical result. This atlas does not establish a
globally solved multi-black-hole spacetime, stable child formation, inherited
constant values, sibling signaling, or identification of our observed cosmos
with any displayed child. No existing paused numerical campaign is resumed.

## Retained parent, wall and clock calculations

The v2–v3 additions have separate domains. Brill–Lindquist supplies a shared
time-symmetric vacuum slice,

$$
h_{ab}=\psi^4\delta_{ab},\qquad K_{ab}=0,\qquad
\psi=1+\sum_i\frac{m_i}{2|x-x_i|}.
$$

Away from the punctures, $R[h]=-8\psi^{-5}\Delta\psi=0$.
The outer ADM charge is $\sum_i m_i$; puncture-end charges are
$M_i=m_i[1+\sum_{j\ne i}m_j/(2d_{ij})]$. For $m_A=1,m_B=1.5,d=12$
these are $2.5,1.0625,1.5625$. They are asymptotic charges, not horizon
masses. The three ends of $\mathbb R^3$ minus two punctures do not establish
formed children or causal channels. The exact companion contribution
$m_j/[2\sqrt{d^2+r^2-2dr\mu}]$ gives angular variation in
$q_{AB}=r^2\psi^4\Omega_{AB}$; its sphere curvature is
$\mathcal K=(1-2\Delta_{S^2}\log\psi)/(r^2\psi^4)$. The recorded 17-point
profiles at coordinate radii $1,1.5$ are probes, not located horizons.

The [N-body scalar-flat theorem](https://arxiv.org/html/0909.1101v1)
instead supplies an existence family retaining prescribed spherical collars
at sufficiently large separation, with positive-mass smooth AE ends. It
supplies neither a numerical separation threshold nor a connector metric
for the BL parameters. The [non-time-symmetric extension](https://arxiv.org/html/1004.1355v1)
needs smooth qualifying $(g,K)$ bodies and timelike ADM energy-momenta.
An ideal material delta-layer is outside its automatically verified scope.
BL's axial Killing field is a nonzero KID, so the separate trivial-KID
localized-neck theorem does not apply directly; this does not exclude other
gluing routes. The 2002 generic positive-mass-bridge argument was withdrawn
by its erratum, while its separate many-Schwarzschild construction remained.
A cutoff interpolation alone does not solve the annular constraints.

The comoving radiation screen fixes $a=A\sin\eta$,
$t=A(1-\cos\eta)$, $s=\sin\chi_0$, $c_0=\cos\chi_0$ and constant
Schwarzschild mass $M$. Constant wall tension $\sigma_0$, with
$P=-\sigma_0$, implies the necessary angular condition

$$
k^2R^4-2c_0kR^3+2MR-A^2s^4=0,\qquad k=4\pi\sigma_0.
$$

It cannot hold on an interval of changing $R$ because $A^2s^4\ne0$.
This excludes either cap orientation and either sign of constant tension,
without excluding a separate de Sitter/vacuum-wall solution. The 33-point
witness uses $A=1,M/A=0.4,\chi_0=\pi/3,1.25\le\eta\le1.45$ and $R>2M$;
samples support the algebraic argument rather than proving it. With vacuum
outside and constant pure tension, conservation also forces zero relative
speed through radiation: boundary motion alone does not repair the source.
The v3 vacuum-wall example and retained-collar family still need compatible
complete Cauchy data and a protected evolution domain before being identified
as one spacetime. Neither edition supplies a multi-region Hamiltonian.

The fixed Bronnikov observer calculation uses
$r=\sqrt{1+\rho^2}$, $\beta=\sqrt{1-A}$ and
$A=1+3\rho+3(1+\rho^2)(\arctan\rho-\pi/2)$, in units $c=L_{\rm throat}=1$.
Its null expansions are
$\theta_\pm=2\rho(-\beta\pm1)/(1+\rho^2)$: exterior, trapped,
minimum-sphere and anti-trapped regions remain distinct. The horizon is
$\rho_h\simeq1.90069160547$, $r_h\simeq2.14770309380$.
The $\mathbb R\times S^2$ interior approaches de Sitter with $\Lambda=9\pi$;
it is not finite-radius FLRW or a compact spatial universe. Its horizon
tangent has $\kappa_h\simeq0.238325799634$ and affine parameter
$\lambda=\lambda_0+C e^{\kappa_h\tau}$, $C>0$. Infinite negative PG time
therefore has finite affine length in that patch; no singularity or global
beginning follows. An assigned ancestry clock requires the convergence test
$T_0=\sum_{m\ge0}(\prod_{j=1}^m q_j)\Delta\widehat\tau_m$.
Infinitely many rooms do not guarantee infinite proper or affine duration.
The geometry's quantum source, physical clock transport and global extension
remain unresolved.
