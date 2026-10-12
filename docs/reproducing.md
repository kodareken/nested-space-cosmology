# Reproducing the root project

There is one package, `nested-space-cosmology`, importing `recursive_horizons`
from root `src`. Use the pinned versions in [pyproject.toml](../pyproject.toml).
Scientific conclusions retain the domains in their owning notes.

## Current work

```sh
python -m pip install -e '.[paper,dev]'
python scripts/check_lab_snapshot.py
python scripts/check_publication.py
python scripts/check_claim_bindings.py
python scripts/test_reassembly.py
python scripts/build_geometry_overview.py
python scripts/build_characteristic_initialization_figure.py
python scripts/build_source_mass_figure.py
python scripts/build_paper.py --check
```

The affected test suite covers atlas attachment incompatibility, disconnected
siblings, causal delays, unequal children, clock-domain exits, unit rescaling,
finite-depth dependence, field transfer, preservation failures and publication
contracts. `make check test paper-check` uses the current environment. `make verify`
also adds explicit historical replay and requires the separate `HISTORY_PYTHON`
profile described below. Running every historical simulation is neither required
nor authorized by these convenience commands.

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

The original v1–v6 atlas manifests record Python 3.14.8 and Matplotlib 3.11.2;
byte-identical historical image replay is verified on Windows.
Their unchanged checkers authenticate recorded files and require byte-identical
regenerated images. This rendering profile is separate from the current paper's
Matplotlib 3.11.1 profile. [requirements-history.txt](../requirements-history.txt)
pins the recorded scientific versions, verified replay dependency Pillow 12.3.0
and the historical checking tools. Pillow and pytest versions were not recorded in
the atlas manifests. The version pins retain the recorded producer environment;
Windows selects the platform on which exact historical replay is verified.
Install it in an isolated environment, without installing the current
package or its paper extras there. From Windows with Python 3.14.8 available:

```powershell
py -3.14 -m venv .venv/history
.\.venv\history\Scripts\python.exe -m pip install -r requirements-history.txt
.\.venv\history\Scripts\python.exe scripts/replay_history.py all
```

The following selectors use that historical interpreter (activate the environment
first if using the short `python` command):

```sh
python scripts/replay_history.py prepare
python scripts/replay_history.py atlas
python scripts/replay_history.py finite
python scripts/replay_history.py companion
python scripts/replay_history.py all
```

These commands select the pre-reassembly source commit explicitly, materialize its
original layout in ignored `.build/replay`, check current preserved artifacts
against those bytes, and run the original checkers. All six atlas bundles are
recomputed within their finite domains. Supporting-paper authentication checks
their frozen source, evidence and archives; it does not rerun their full campaigns.
The original context contains the former foundation paper at its historical path.

With Make installed, select the historical interpreter explicitly:

```sh
make history-check HISTORY_PYTHON=.venv/history/Scripts/python.exe
make verify PYTHON=python HISTORY_PYTHON=.venv/history/Scripts/python.exe
```

`HISTORY_PYTHON` defaults to `PYTHON` for an already selected historical environment;
it never chooses another interpreter after a failed check. GitHub Actions runs the
current paper checks on Ubuntu, the legacy historical replay in a Windows job
with Matplotlib 3.11.2, and the explicit v7/v8 replays in a separate Windows job
with Matplotlib 3.11.1. Both Windows jobs use Python 3.14.8 and retain the original
evidence comparisons. Other historical rendering profiles have not been
established as byte-identical.

Current atlas producers use root sources. Their source bindings differ from the
historical bundles, so passing an old bundle to a current producer must fail.
To create a new bundle, supply a fresh `--output` directory; creation is exclusive.
A changed historical context, changed retained evidence, missing Git object or
corrupt archive is an error, never a reason to silently select alternative bytes.

### Atlas v7: original root context

Atlas v7 predates the evaluated-response successor. Authenticate its retained
evidence and materialize its original sources with `prepare-v7`; use `atlas-v7`
for the original checker's full replay, including exact figure bytes:

```sh
python scripts/replay_history.py prepare-v7
python scripts/replay_history.py atlas-v7
```

These selectors use commit `7fee894ad5a4da2ba378f724e45e926401f82670` under
`.build/replay/7fee894ad5a4`. `prepare-v7` performs no scientific replay and
selects no renderer. Full `atlas-v7` replay requires its recorded Windows
profile: Python 3.14.8, NumPy 2.5.1, SciPy 1.17.1, SymPy 1.14.0 and
Matplotlib 3.11.1. The existing Windows validation environment supplies this
profile; its fallback invocation is
`python scripts/run_python.py scripts/replay_history.py atlas-v7`.
The v1–v6 history environment above has a different Matplotlib pin.

Neither v7 selector changes the legacy `all` selector or v1–v6 replay.
An incompatible runtime, modified preserved output or missing source object
fails explicitly. Current changes to scientific owners are not substituted
for v7's original source bytes.

### Atlas v8: original prepared-response context

Atlas v8 predates the analytic inverse and the updated literature/atlas owners.
Use its original commit explicitly to authenticate or replay that result:

```sh
python scripts/replay_history.py prepare-v8
python scripts/replay_history.py atlas-v8
```

These selectors use commit `1154baf290d901e73058392e143d0acb37d22951` under
`.build/replay/1154baf290d9`. `prepare-v8` authenticates retained evidence and
source bindings against that commit without running a field solve. `atlas-v8`
runs the unchanged original checker and its finite numerical budget, including
exact figure comparisons, in the same Windows/Python 3.14.8/Matplotlib 3.11.1
profile specified for v7 above. The Windows fallback invocation is
`python scripts/run_python.py scripts/replay_history.py atlas-v8`.

These commands do not add v8 to legacy `all`, loosen numerical/image checks,
or silently recover from a current-source mismatch. Calling the current v8
producer's `--check` on the frozen canonical bundle after a bound owner changes
must fail. Neither the historical result nor its hash is rewritten to make the
current source fit. Use the explicit context to replay the original result.

Hosted-runner boundary: on 12 October 2026, the first Windows v8 CI replay
failed its `prepared-response.json` comparison; a subsequent runner passed
the unchanged original producer and tolerances. The differing fields from
that first run were not retained, so its cause remains unresolved. Package
versions alone do not bind CPU dispatch or BLAS configuration. Each v7/v8 CI
run now logs that runtime, and a failure preserves original-producer outputs
and field-level comparison diagnostics. A successful later run does not
certify identical replay on every machine or excuse a future mismatch.

A subsequent [v5 diagnostic run](https://github.com/kodareken/nested-space-cosmology/actions/runs/38313295989/job/115024003729)
retained a PNG mismatch while its SVG matched exactly and its JSON passed
the original numerical tolerance. Only two pixels differed, each by one
red/green level. That runner selected OpenBLAS SkylakeX; the known passing
local replay selects Haswell in both NumPy and SciPy. The original atlas
manifests did not record CPU dispatch. The Windows
replay jobs now explicitly select and verify `OPENBLAS_CORETYPE=Haswell`
before running the unchanged checkers. This binds a previously automatic
runtime choice using [OpenBLAS's documented selector](https://github.com/OpenMathLib/OpenBLAS/blob/develop/USAGE.md),
without changing evidence or tolerances. It does not retrospectively identify
the uncaptured first v8 mismatch. For local PowerShell replay in a compatible
x86-64 environment, set `$env:OPENBLAS_CORETYPE='Haswell'` before the commands
above; other processor architectures require their own established profile.

The former `make reproduce`, `make reproduce-development` and `make reproduce-exact`
campaign routes now stop with guidance. Their frozen release inputs differ from
the merged current sources; running the old producer against the current package
would fail authentication. Full campaign recomputation is no longer a GitHub Actions
option. A deliberate rerun first needs a named question, a fresh admission and a
finite budget, then selection of the source context produced by
`python scripts/replay_history.py prepare`. Run its historical producer with that
context's isolated source environment; never change pinned hashes or substitute
the current package to make an old record pass. Routine authentication remains
`make history-check` or `python scripts/replay_history.py all`.

The [migration checker](../scripts/reassemble_repository.py) authenticates every
destination and explicit successor. [Preservation and recovery](repository-consolidation.md)
describe the original-path manifest and the named implementation variants.

## Windows runtime

The existing validation environment uses Python 3.14.8 and the project's pinned
scientific stack. If its moved Windows redirector cannot locate `pyvenv.cfg`, use
`python scripts/run_python.py <script-or-python-arguments>` with the base interpreter;
for example, `python scripts/run_python.py scripts/test_reassembly.py`. Do not install a
second project package or restore the old directory. Normal installations use
their activated interpreter. Root `scripts/lab.py` is only a compatibility launcher:
it now runs in this root and does not switch scientific packages.

## Scope

The main [paper](../paper/nested-space-cosmology.md) is a living atlas manuscript.
The [paper guide](papers.md) routes archived and supporting editions. Spherical
symmetry, specified states and fixed parameters are assumptions of individual
tests. No new constants, propagation or recurrence campaign was admitted during
repository reassembly. New calculations use the [existing admission rule](prior-art-and-open-claim.md).
