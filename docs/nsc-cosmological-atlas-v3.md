# Moving boundaries and a parent with spherical collars

This successor advances two distinct connections in the exterior-first atlas:
the source law of a moving spherical wall, and a published way to retain
spherical regions inside a common nonspherical parent. The
[v2 atlas](nsc-cosmological-atlas-v2.md), its manuscript and all bound sources
remain frozen. This note owns the successor's assembly, not a new gluing or
thin-shell theorem.

## Decision and implementation map

The previous unresolved label, "moving boundary", hid a physical dependency:
a boundary cannot acquire arbitrary motion while its material law and bulk
flux stay fixed. We compare three routes. A pure constant-tension wall with
radiation inside and vacuum outside is constrained by conservation. A
published de Sitter/vacuum wall changes the interior and supplies a separate
moving solution. A published many-body construction changes the parent family
and supplies exact Schwarzschild collars under its existence assumptions.
None of these changes silently repairs the original Brill-Lindquist data.

- [Moving-boundary owner](nsc-atlas-moving-boundaries.md): imports conservation
  and signed shell dynamics; its adapter owns fixed parameters, finite samples,
  proper-time quadrature and checks.
- [Parent-gluing owner](nsc-atlas-parent-gluing.md): owns the literature reading,
  theorem hypotheses, corrections and unverified interface with the wall.
- [Successor producer](../scripts/derive_nsc_cosmological_atlas_v3.py): reuses
  v2's parent and source screen, adds the independent candidates and emits
  the map, numerical figure, machine-readable records and narrative.
- [Working manuscript v3](../../paper/atlas-working-draft-v3.md): explains the
  combined story; the predecessor manuscript stays bound to v2 evidence.

No new dependency, fitted coupling, offspring probability or observer clock
is introduced. The wall adapter owns its finite proper-time interval. The
gluing route owns symbolic existence conditions rather than invented numeric
thresholds. Failure stays local to the affected attachment.

## Calculation admission

**Borrowed:** the existing v2 baseline and source obstruction; general
Israel/Codazzi matching and published signed de Sitter/Schwarzschild dynamics;
published many-body vacuum initial-data gluing. Exact versions, equations and
reading boundaries are in the two owners above.

**Missing connection:** whether boundary motion alone rescues the radiation
candidate, and whether published geometry already supplies a common parent
with exact spherical collars. Joining the moving matter wall's Cauchy data to
such a collar is a further question, not an output of either borrowed result.

**Decision and unlock:** apply the conservation obstruction without a new
radiation evolution; admit a separate finite vacuum-wall patch; retain a
theorem-backed parent family whose interface can be checked next. Do not
spend a solver campaign rediscovering the known shell or gluing construction.

**Finite domain and budget:** one fixed wall model, 33 radius samples and its
bounded quadrature/checks, as recorded by the adapter; one two-body gluing
scenario represented by hypotheses; inherited descent cutoff two. Generation
and replay target less than 60 CPU seconds. There is no global spacetime
simulation or search over source coefficients.

**Stop or redirect:** stop this pure-wall/radiation route on the conservation
constraint. If the signed vacuum patch fails the unsquared junction, do not
draw it as admissible. If the many-body theorem's hypotheses or the Cauchy
interface are unverified, retain those exact missing conditions; do not infer
an attachment from matching mass labels.

**Verification and attribution:** direct tensor contraction, unsquared
junction and temporal-curvature checks support the wall adapter. The atlas
checks distinguish imported local matching, existence of a parent family,
and an actual shared spacetime. Neither this assembly nor evaluating published
formulas is a claim of new-to-world mathematics. Bounded searches do not
certify priority.

## Results and replay

The [saved map](../results/development/nsc-cosmological-atlas-v3/map.svg),
[moving-boundary record](../results/development/nsc-cosmological-atlas-v3/moving-boundary.json),
[parent-family card](../results/development/nsc-cosmological-atlas-v3/parent-gluing.json)
and [illustrated narrative](../results/development/nsc-cosmological-atlas-v3/narrative.md)
are created together. The manifest binds all outputs and transitive source
files, including the manuscript. Creation refuses an existing directory;
`--check` replays read-only. The v1 and v2 producers remain independently
replayable.

```sh
python scripts/lab.py scripts/derive_nsc_cosmological_atlas_v3.py --check
python scripts/lab.py -m pytest tests/test_nsc_atlas_moving_boundaries.py tests/test_nsc_cosmological_atlas_v3.py -q
```

The shell's proper time is local to its worldtube. It is not a child's age or
a shared parent clock. The alternate parent is an existence family, not an
explicit replacement metric for the sampled Brill-Lindquist slice. No sibling
signal path, formation history, inherited constants or local dark-sector
prediction is supplied by these records.

## Finish-line proof and failure forecast

Completion is falsifiable by an unsquared-junction or flux discrepancy, a
map that promotes a theorem's unverified hypothesis to a solved connection,
a changed v1/v2 source hash, or a failed replay of the bound successor. The
new manuscript must state the same domains as the adapter and map.

Three likely obstacles remain. First, the moving wall and a time-symmetric
vacuum collar may need incompatible initial slices; the next interface check
must compare complete Cauchy data and retained sides. Second, a gluing theorem
may guarantee existence without a usable separation bound; preserve that
limit instead of assigning the Brill-Lindquist coordinate separation to it.
Third, reproducing global shell or gluing simulations could consume the work
without closing either gap; stop at the finite imported control and return to
the interface question.
