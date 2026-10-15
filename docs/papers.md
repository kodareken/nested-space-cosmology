# Paper guide

The one living manuscript is **Nested-Space Cosmology: The Geometry of a Universe**:
[Markdown source](../paper/nested-space-cosmology.md) and
[illustrated PDF](../paper/nested-space-cosmology.pdf). Its revision is independent
of the software release. Future writing updates these files. The
[catalog](../paper/catalog.json) identifies current, supporting and historical editions.

## The common research thread

The current working manuscript includes the
[smooth exterior response](nsc-atlas-bulk-response.md): the prepared field
now drives a calculated radius and lapse on the original characteristic
rectangle. Independent mass transport and null constraints check that
propagation. The next connection is its continuation to the moving wall
and through the cap, needed for corrected observer measurements.

Regions, boundaries and fields define possible nested geometries. Attachment and
initial state then determine what an interior observer can access and measure.
Spherical symmetry supplies the first tractable charts; it does not assert global
multi-hole topology or a literal shape for every particle.

The paper moves from motivation and mathematical language to the map, causal
access, calculated responses and missing physical connections. The scalar pulse
and finite transfer/memory results retain their different domains. A mode-space
complement is not automatically an exterior spacetime. Descent, attachment,
influence and correlation remain distinct.

The [observer-inverse proof](nsc-atlas-observer-inference.md) connects the
calculated signal back to a preparation within the fixed geometry. Its
positive compact-source exclusions explain why the next comparison must
state what preparations are physically allowed before inferring an exterior.

The [conformal stress check](nsc-atlas-scalar-wall-consistency.md) connects
the same prepared probe to the wall's quadratic stress conditions. Its
surface coefficient and tangential flux jump vanish; the smooth bulk
gravitational evolution requires compatible initial data. The
[characteristic initialization](nsc-atlas-characteristic-initialization.md)
now computes the required source-leg radius correction with the original
compact preparation, corner geometry and null normalization fixed. This
connects the field to one gravitational constraint; transverse evolution,
mass change and the coupled wall response require additional equations.
The [source-leg mass successor](nsc-atlas-source-mass.md) now supplies
the mixed-radius connection on that same initial surface. The geometric
mass has a negative constant quadratic coefficient after the compact
trace, checked through full frame conversion and improved-stress transport.
The paper distinguishes this local initial-data result from bulk evolution,
asymptotic energy and the still-open dynamical wall response.
The [normal-force check](nsc-atlas-scalar-normal-force.md) now adds zero
regular normal-stress jump for the same conformal probe. It supplies a
necessary mean-curvature condition for a later sourced wall, with bulk
metric and wall displacement both included in that variation.
The [wall-mass closure](nsc-atlas-wall-mass-closure.md) now shows that the
fixed leading trace supplies an analytic exterior-side mass functional,
using both null fluxes and the full frame correction. Its numerical
evaluation remains separate from the unsolved cap/wall geometry. The next
decision is a controlled finite evaluation or a named remaining junction
input, according to which unlocks the next physical connection.

## Supporting and historical editions

- [Foundation v0.26.0 PDF](../paper/archive/foundation-v0.26.0/paper/nested-space-cosmology.pdf),
  [source](../paper/archive/foundation-v0.26.0/paper/nested-space-cosmology.md) and
  [archive receipt](../paper/archive/foundation-v0.26.0/archive-manifest.json),
  retaining metadata, claim bindings and the renderer.
- [Regional Transfer, Geometric Feedback, and Local Memory](../paper/finite-regeneration.pdf):
  the finite realization article, [LaTeX](../paper/finite-regeneration/main.tex)
  and [evidence manifest](../paper/finite-regeneration-manifest.json).
- [Local incoming-gate companion](../paper/local-incoming-gate-draft.pdf):
  the preserved OPEN application checkpoint; its campaign remains paused.
- Atlas [v1 evidence](https://doi.org/10.5281/zenodo.23090015),
  [v2 manuscript](https://doi.org/10.5281/zenodo.23090015),
  [v3](https://doi.org/10.5281/zenodo.23090015), [v4](https://doi.org/10.5281/zenodo.23090015),
  [v5](https://doi.org/10.5281/zenodo.23090015), [v6](https://doi.org/10.5281/zenodo.23090015):
  source-bound historical stages, not competing current drafts.

Historical text retains original logical paths. Materialize its original layout
with `python scripts/replay_history.py prepare` when following old `lab/` links.
The living paper uses root paths. The [reproduction guide](reproducing.md)
distinguishes authentication, numerical replay and publication builds.
No arXiv submission occurred in reassembly.
