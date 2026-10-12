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

The current finite preparation bundle can be regenerated without a PDE solve:

```sh
python scripts/derive_nsc_cosmological_atlas_v7.py --output .build/atlas-v7-review
python scripts/derive_nsc_cosmological_atlas_v7.py --output .build/atlas-v7-review --check
```

Creation requires a fresh output directory. This evaluates the admitted finite
geometry and analytic bound, not a new numerical observer waveform. The
[preparation owner](nsc-atlas-field-preparation.md) states the domain and limits.

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

The original atlas manifests record Python 3.14.8 and Matplotlib 3.11.2;
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
