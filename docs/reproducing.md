# Reproducing the root project

There is one package, `nested-space-cosmology`, importing `recursive_horizons`
from root `src`. Use the pinned versions in [pyproject.toml](../pyproject.toml).
Scientific conclusions retain the domains in their owning notes.
The current paper and test dependencies use Matplotlib 3.11.2; the scientific
stack includes python-flint 0.9.0.

## Current work

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

`check_lab_snapshot.py` authenticates the current paths, hashes and sizes in the
v2 [retained-file receipts](lab-snapshot.json), plus the unchanged objects and
pins in [source history](../archive/source-history/manifest.json).
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
geometry and the corrected observer signal remain uncomputed.

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
installing the current project package there:

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
