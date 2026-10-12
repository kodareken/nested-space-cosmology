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
python scripts/build_paper.py --check
```

The affected test suite covers atlas attachment incompatibility, disconnected
siblings, causal delays, unequal children, clock-domain exits, unit rescaling,
finite-depth dependence, field transfer, preservation failures and publication
contracts. `make check test paper-check` uses the current environment. `make verify`
also adds explicit historical replay and requires the separate `HISTORY_PYTHON`
profile described below. Running every historical simulation is neither required
nor authorized by these convenience commands.

The current atlas v8 evaluates the admitted finite exterior preparation through
the moving wall to the interior observer:

```sh
python scripts/derive_nsc_cosmological_atlas_v8.py --output .build/atlas-v8-review
python scripts/derive_nsc_cosmological_atlas_v8.py --output .build/atlas-v8-review --check
```

Creation requires a fresh output directory; `--check` authenticates the bound
sources and outputs before regeneration. Omit `--output` with `--check` to check
the canonical v8 bundle. The fixed input, finite numerical budget, refinement
diagnostics and independent analytic tail bound are recorded in the
[response owner](nsc-atlas-prepared-response.md). Refinement diagnostics are not
a rigorous numerical enclosure. The [v8 guide](nsc-cosmological-atlas-v8.md)
describes the outputs and the unchanged causal limits.

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
current paper checks on Ubuntu and the exact historical replay in a separate Windows
job with Python 3.14.8. Both jobs retain their original evidence comparisons. Other
historical rendering profiles have not been established as byte-identical.

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
