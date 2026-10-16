# Paper guide

The one living manuscript is **Nested-Space Cosmology: The Geometry of a Universe**:
[Markdown source](../paper/nested-space-cosmology.md) and
[illustrated PDF](../paper/nested-space-cosmology.pdf). Its revision is independent
of the software release. Future writing updates these files. The
[catalog](../paper/catalog.json) identifies current, supporting and historical editions.

## The common research thread

The current working manuscript carries the prepared field through the
[smooth exterior response](nsc-atlas-bulk-response.md),
[moving wall](nsc-atlas-wall-evolution.md) and
[cap metric and central null clock](nsc-atlas-cap-response.md).
Independent radial and timelike mass paths check this quadratic geometry.
The central arrival map now includes both travel through the changed cap
and the moved wall event at a fixed exterior label. The higher-order scalar
calculation uses that same source and clock convention.
Section 13.13 now derives and evaluates the cap's outgoing scalar backscatter
and its exterior wall derivative, with five numerical comparisons. This
boundary input is used by Section 13.14, which now computes the exterior
scalar correction with separate forcing-cache, spatial and ODE comparisons.
The [central reconstruction](nsc-atlas-central-scalar-response.md) now
evaluates the formal cubic coefficient at five interior arrival cuts. Its
[radial](../results/development/nsc-atlas-central-scalar-v1/central.json) and
[ray](../results/development/nsc-atlas-central-ray-v1/ray.json) comparisons
retain separate incoming-quadrature and exterior cache/grid/ODE sensitivities.
The [exact nonlinear model](nsc-atlas-finite-amplitude.md) specifies the
source, smooth equations and material interface. A full finite-amplitude
solution and uniform remainder are not established by these coefficients.

The [characteristic-background coupled window](current-result.md#characteristic-background-coupled-evolution-16-october-2026)
now supplies six finite-amplitude numerical histories with a common moving
wall and cap. Carrying the known incoming pulse analytically removes the
dominant cap transport error. Paired wall geometry and clock responses
agree closely under the declared spatial/time refinements, while the
scalar time sensitivity, cap outgoing response and exterior equation
residuals remain explicit. The [comparison figure](figures/nsc-atlas-background-evolution.png)
and saved fields connect this numerical advance to the manuscript without
claiming a uniform error enclosure or central-arrival solution.

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
The paper distinguishes this local initial-data result from asymptotic
energy and from the subsequent bulk and wall-response calculations.
The [normal-force check](nsc-atlas-scalar-normal-force.md) now adds zero
regular normal-stress jump for the same conformal probe. It supplies a
necessary mean-curvature condition for a later sourced wall, with bulk
metric and wall displacement both included in that variation.
The [wall-mass closure](nsc-atlas-wall-mass-closure.md) now shows that the
fixed leading trace supplies an analytic exterior-side mass functional,
using both null fluxes and the full frame correction. Its numerical
evaluation and compatible cap/wall response are now recorded in the
successor calculations above. Those formal coefficients retain their
perturbative domain; formation, stability and a physical preparation
mechanism remain later questions.

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
