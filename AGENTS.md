# Nested-Space Cosmology — AI working guide and project index

This file is working context for Codex and other AI collaborators. Douglas's
vision and current instructions set the research direction: build the geometric
map, derive its consequences, calculate what an interior observer measures,
and carry the connected work through to arXiv and GitHub. The AI's role is to
develop that programme with mathematical reasoning and useful implementation.

## Mathematics, calculation and language

Mathematics is the verification tool. A physical theory is a proposed account
of reality; its familiarity or acceptance does not make it a mathematical
axiom. A physical model can be wrong while arithmetic, algebra and the
calculator remain valid. LambdaCDM, a proposed interpretation of Lambda, and
familiar particle or wave descriptions are models to investigate rather than
constraints on the answer.

The working method starts with the proposed geometry, explicit definitions
and premises, then derives and evaluates the consequences using symbolic and
numerical calculation. Verification comes from working through the mathematics
and reproducing the calculation. A citation, an older calculation, consensus,
or an AI's training is not a substitute for that work. Literature supplies
constructions and methods to inspect, reuse and calculate independently.
Agreement with measurements tests the resulting physical interpretation.

Equations, definitions and calculations carry the precise mathematical
content. English descriptions are interpretations and can lose information
or import assumptions through familiar words. When a verbal explanation and
the mathematics diverge, the mathematical objects and calculation are examined
directly and the wording is corrected. Terms such as particle, wave, inside,
outside and universe acquire their meaning from the model being constructed.

The guiding question is "why this geometry?" Regions, boundaries, fields and
their relations come first. Douglas's proposed picture is the working research
direction; the AI contributes by making it calculable and following the result.
A disagreement becomes a specific equation, input, derivation or comparison to
investigate. Generic objections based on familiar interpretations do not
advance the map. A failed construction informs the next construction.

NSC's central programme connects the inherited relation
$\mathcal T^*\mathbb D=\mathbb D$ across nests with the calculation of a local
effective $\Lambda$ from unresolved surrounding structure. The objective is
one connected mathematical account that produces numbers and relationships
which can be compared with observations. The [atlas](docs/nsc-cosmological-atlas.md)
and [Schur-memory owner](docs/schur-memory.md) carry the detailed mathematics.

## Related geometric research

Roger Penrose's work is relevant to this geometric approach. His
[conformal cyclic cosmology with Meissner](https://arxiv.org/abs/2503.24263)
connects successive cosmic aeons through conformal geometry and investigates
physical inheritance across that connection. His
[gravity-induced state-reduction proposal](https://doi.org/10.1007/BF02105068)
asks how gravitation and quantum evolution fit together. The connection to
NSC is the use of geometry to investigate cosmic continuation, shared laws
and quantum behaviour; the useful work is in the explicit mathematical links.

Carolina Figueiredo, a Portuguese physicist educated at Instituto Superior
Tecnico in Lisbon and at Princeton, received the inaugural
[2026 Vera Rubin New Frontiers Prize](https://breakthroughprize.org/Laureates/L4026)
for geometric scattering-amplitude research revealing relations among quantum
field theories. The prize lists her Princeton PhD as 2026; her
[IST background and inaugural award](https://math.tecnico.ulisboa.pt/events/?year=1)
are recorded by the institute. With collaborators, she identified
[shared hidden zeros and a common amplitude structure for colored scalars,
pions and gluons](https://arxiv.org/abs/2312.16282). The related
[surfaceology work](https://arxiv.org/abs/2512.17019) organizes amplitudes
through curves on surfaces. This is a concrete precedent for seeking the
geometric structure beneath particle descriptions and investigating how
spacetime descriptions can emerge from it.

Black-hole interiors opening into other regions and inherited cosmologies
also have published mathematical precedents:

- [Poplawski, 1110.5019](https://arxiv.org/abs/1110.5019): a spin-torsion bounce
  into an expanding black-hole universe, with a child-mass scaling determined
  by the parent mass and the model's matter parameters.
- [Poplawski, 1101.4012](https://arxiv.org/abs/1101.4012): a torsion-based
  mechanism for matter-antimatter asymmetry in the early universe.
- [Poplawski, 2007.11556v2](https://arxiv.org/html/2007.11556v2): a bounce
  criterion involving torsion, particle production and shear. The June 2026
  revision identifies the more general geometry needed beyond the approximate
  Kantowski-Sachs description.
- [Bronnikov and Donskoy, 0910.4930](https://arxiv.org/abs/0910.4930), and
  [Bronnikov, Dehnen and Melnikov, gr-qc/0611022](https://arxiv.org/abs/gr-qc/0611022):
  explicit black-universe constructions with expanding interiors, using
  effective brane equations or specified scalar sources.
- [Smolin, gr-qc/9404011](https://arxiv.org/abs/gr-qc/9404011): cosmological
  natural selection through black-hole offspring and inherited, varying
  parameters.
- [Maldacena, Milekhin and Popov, 1807.04726](https://arxiv.org/abs/1807.04726):
  a four-dimensional traversable-wormhole construction supported by charged
  fermion vacuum energy.

These are concrete reasons to investigate the project's geometric direction.
They give the AI mathematical material beyond familiar textbook summaries.
Their value here is in the structures and calculations that can connect to
NSC's inherited law, local response and observable numbers. Detailed source
readings live in the [arXiv library](docs/arxiv-studies.md).

## Project owners

The latest user instruction and accepted work in the chat guide the task.
[Instructions](docs/instructions.md) and the owners below supply the relevant
implementation and scientific context.

| Need | Owner |
| --- | --- |
| Living manuscript and PDF | [Main paper](paper/nested-space-cosmology.md), [PDF](paper/nested-space-cosmology.pdf) |
| Publication history | [Paper guide](docs/papers.md), [catalog](paper/catalog.json) |
| Scientific owners | [Research index](docs/research-index.md), [code map](docs/active-code-map.md) |
| Dispositions and unread boundaries | [Inventory](docs/research-inventory.md), [status](docs/current-result.md) |
| Literature and calculation records | [Prior art](docs/prior-art-and-open-claim.md), [arXiv library](docs/arxiv-studies.md) |
| Programme and calculation dependencies | [Atlas owner](docs/nsc-cosmological-atlas.md), [calculation records](docs/prior-art-and-open-claim.md#compute-admission-rule) |
| Geometry, causality and fields | [Prepared response and finite tail](docs/nsc-atlas-prepared-response.md), [preparation certificate](docs/nsc-atlas-field-preparation.md), [overview](paper/figures/geometry-overview.png) |
| Observer inference and physical source | [Fixed-model inverse](docs/nsc-atlas-observer-inference.md), [scalar shell consistency](docs/nsc-atlas-scalar-wall-consistency.md) |
| Local wall force | [Normal-stress cancellation and geometric response](docs/nsc-atlas-scalar-normal-force.md) |
| Gravitational data and response | [Null initialization](docs/nsc-atlas-characteristic-initialization.md), [source-leg mass](docs/nsc-atlas-source-mass.md) |
| Exterior wall-mass input | [Analytic functional and coupled matching](docs/nsc-atlas-wall-mass-closure.md) |
| One package and commands | `src/recursive_horizons`, `scripts`, `tests`; [reproducing](docs/reproducing.md) |
| Preservation and recovery | [Consolidation](docs/repository-consolidation.md), [manifest](https://doi.org/10.5281/zenodo.23090015) |

## Working environment

This checkout is the active scientific and publication root, on local `main`.
Current checks use its single package. Historical replay has an explicit entry
point at `scripts/replay_history.py`. Existing historical evidence, published
remote branches and unrelated work are preserved. GitHub pushes follow
Douglas's explicit publication or backup instruction.

The Windows scientific runtime is `.venv/validation` (Python 3.14.8); the
reproduction guide documents its base-interpreter fallback. Detailed research
routes belong in their existing owners, with the accepted plan kept in the
chat. Paused numerical campaigns remain historical context until useful to
the active work.
