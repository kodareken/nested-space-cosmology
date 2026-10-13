# Atlas v8: from a finite exterior preparation to an interior signal

The [response owner](nsc-atlas-prepared-response.md) owns the equations,
analytic bounds and numerical admission. This evidence edition evaluates the
operator certified by [finite preparation](nsc-atlas-field-preparation.md). The
[canonical paper](../paper/nested-space-cosmology.md) remains the living
narrative; the [programme owner](nsc-cosmological-atlas.md) connects the
calculation to the wider map.

## Result and limits

The retained v7 step is a conditional existence certificate for independent
scalar data in a finite characteristic strip on the white-hole side of the
fixed wall, away from the singularity. Its integral-operator norm bound is
$0.138821<1$. It defines a unique trace and the composed interior response
operator; v7 did not evaluate the waveform. The research dependency graph
records prerequisites and stopping conditions separately from causal links.
Neither its nodes nor this certificate creates a route from the right parent
or connects the nine earlier scenarios.

An independently specified compact scalar input in the retained Schwarzschild
sector IV now has an evaluated route through the fixed moving wall to the
interior centre. The finite central waveform is bipolar, followed by a
positive field tail after the input ends. The tail satisfies an independent
analytic comparison and is resolved under spatial and geometric refinement.
The analytic sign result and the numerical refinement diagnostics are distinct;
the latter do not enclose the exact continuum solution rigorously.

The geometry, classical scalar, coupling, source profile and zero cap data
remain those of the v7 preparation. An initial unsegmented integration stopped
at its call budget. Splitting the same integration at known source and moving
wall/grid events passed the unchanged accuracy thresholds and finite budget.
The [response owner](nsc-atlas-prepared-response.md) records both attempts and
the independent checks. Established integration methods are borrowed; this
is a calculated connection in one specified geometry, not a new generic solver.

All ten v7 scenarios remain unchanged. The eleventh scenario advances only
the `wall-trace` and `observer-response` research nodes for this preparation.
The right static parent remains causally excluded. Prescribing data in IV
does not identify their remote source or establish a quantum state, formation
history, backreaction, stability, energy flux or permanent detector memory.
The computed tail is confined to the declared observer window. No inference
about its infinite-time limit or a unique exterior follows from it.

## Evidence and figures

The [manifest](../results/development/nsc-cosmological-atlas-v8/manifest.json)
binds source bytes and these outputs:

| Output | Purpose |
| --- | --- |
| [atlas.json](../results/development/nsc-cosmological-atlas-v8/atlas.json) | Eleven finite scenarios with separate attachment, causal and research annotations. |
| [programme.json](../results/development/nsc-cosmological-atlas-v8/programme.json) | Prerequisites, decisions and downstream questions; no research edge creates a causal path. |
| [preparation-access.json](../results/development/nsc-cosmological-atlas-v8/preparation-access.json) | Retained IV domain, independent data, analytic operator certificate and observer access. |
| [prepared-response.json](../results/development/nsc-cosmological-atlas-v8/prepared-response.json) | Sampled waveform, analytic tail comparison, independent controls, refinement and execution budget. |
| [narrative.md](../results/development/nsc-cosmological-atlas-v8/narrative.md) | Readable connection, access limits and measured response. |
| [prepared-response.png](../results/development/nsc-cosmological-atlas-v8/prepared-response.png), [SVG](../results/development/nsc-cosmological-atlas-v8/prepared-response.svg) | Source and wall reduced fields, central response and zero-kernel diagnostic, then the magnified finite tail. |

The zero-kernel curve is a diagnostic of the calculation; it is not another
physical spacetime. Tail comparison points at normalized source coordinates
0.8125, 0.875 and 1 were declared before execution. Refinement and these
independent analytic comparisons test different possible errors.

## Reproduce

Use the explicit original root context in the
[reproduction guide](reproducing.md#atlas-v8-original-prepared-response-context):

```sh
python scripts/replay_history.py prepare-v8
python scripts/replay_history.py atlas-v8
```

The first command authenticates the retained evidence and materializes commit
`1154baf290d901e73058392e143d0acb37d22951` without a solve. The second runs its
unchanged original checker with the recorded Windows scientific/renderer profile.
The [producer](../scripts/derive_nsc_cosmological_atlas_v8.py) binds the
[response implementation](../src/recursive_horizons/nsc_atlas_prepared_response.py),
its scientific owners and upstream evidence; it excludes the canonical
manuscript to avoid a source/evidence cycle. Because current owners have evolved,
the current producer's direct `--check` of this frozen bundle fails strictly.
The replay command selects history explicitly; it does not hide that mismatch.
Changed sources, altered outputs and attempted overwrites fail.
The [producer tests](../tests/test_nsc_cosmological_atlas_v8.py)
cover those failures and preservation of the preceding scenarios.

Replay preserved v7 with `python scripts/replay_history.py atlas-v7`, selecting
its original sources and recorded Windows runtime as described in the
[historical replay instructions](reproducing.md#atlas-v7-original-root-context).
The v1–v6 selectors and frozen evidence remain unchanged.

The [analytic successor](nsc-atlas-observer-inference.md) now derives the
fixed-model preparation inverse and tests restrictions on the allowed source
class. It supplies no second geometry. Whether different admissible
exterior/preparation pairs yield the same response still requires a named
compatible comparison and an independent preparation constraint or discriminator.
The stored v8 evidence remains unchanged.
