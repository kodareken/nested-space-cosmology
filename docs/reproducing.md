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
contracts. `make verify` adds explicit historical replay. Running every historical
simulation is neither required nor authorized by these convenience commands.

The paper builder renders local figures, authenticates their bytes and compares
two builds in check mode. Visual inspection is still required for a manuscript
change. The software release version and the manuscript revision are independent.

## Preserved editions and source contexts

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

Current atlas producers use root sources. Their source bindings differ from the
historical bundles, so passing an old bundle to a current producer must fail.
To create a new bundle, supply a fresh `--output` directory; creation is exclusive.
A changed historical context, changed retained evidence, missing Git object or
corrupt archive is an error, never a reason to silently select alternative bytes.

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
