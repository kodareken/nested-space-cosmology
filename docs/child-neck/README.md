# Child and neck: an exploratory chain, 7 to 9 October 2026

The [10 October scope correction](child-and-neck-status-note.md) distinguishes the
approximate density-cancellation scale from a real closed-FLRW turning point.
The original records remain unchanged and exploratory. Review the
[status and corrections](child-and-neck-status-note.md) before reuse.

Status: exploratory and hypothesis-first. The working picture is a space that makes black holes, each hosting a space that makes more, with the two sides of a throat treated as one space. The notes assume that picture and compute what follows, with one Dirac field, a 10 solar-mass parent and order-of-magnitude accuracy unless a note says otherwise. Nothing here is entered in the [claim ledger](../claim-ledger.md), and none of it is a measured repository result or a verdict on the picture. Each note names its own unaccounted pieces. A proposal is not an observation.

Drafted with AI assistance under Douglas Ek's direction. Human review has not been done.

Start with [child-and-neck-status-note.md](child-and-neck-status-note.md): one table of current results, a log of six corrections to earlier notes, the named missing ingredients, and the open forks in recommended order.

## Notes, scripts and records

| Note | What it computes | Script | Record and figure |
|---|---|---|---|
| [kasner-shear-production-note.md](child-and-neck-status-note.md) | Pair production by shear in a stretching (Kasner) collapse, and where the stretch stops | [kasner_shear_production.py](https://doi.org/10.5281/zenodo.23090015), [kasner_stopping_scales.py](https://doi.org/10.5281/zenodo.23090015), [kasner_shear_production_figure.py](https://doi.org/10.5281/zenodo.23090015) | [kasner_shear_production.json](../../results/development/child-neck-v1/kasner_shear_production.json) [png](https://doi.org/10.5281/zenodo.23090015), [kasner_stopping_scales.json](../../results/development/child-neck-v1/kasner_stopping_scales.json) [png](https://doi.org/10.5281/zenodo.23090015) |
| [eck-collapse-pair-production-note.md](child-and-neck-status-note.md) | Pair production in the Einstein-Cartan black-hole universe: a self-consistency check | [eck_collapse_consistency_check.py](https://doi.org/10.5281/zenodo.23090015) | [json](../../results/development/child-neck-v1/eck_collapse_consistency_check.json) |
| [child-energy-floor-and-thread-length-note.md](child-and-neck-status-note.md) | Child energy, the no-production baseline, and what sets the thread length | [child_energy_floor.py](https://doi.org/10.5281/zenodo.23090015) | [json](../../results/development/child-neck-v1/child_energy_floor.json) [png](https://doi.org/10.5281/zenodo.23090015) |
| [child-expansion-and-matter-fraction-note.md](child-and-neck-status-note.md) | After the bounce: how far the child expands, what it weighs, and what the antimatter fraction does | [child_expansion_ladder.py](https://doi.org/10.5281/zenodo.23090015), [child_matter_fraction.py](https://doi.org/10.5281/zenodo.23090015) | [child_expansion_ladder.json](../../results/development/child-neck-v1/child_expansion_ladder.json) [png](../../results/development/child-neck-v1/child_expansion_ladder.png), [child_matter_fraction.json](../../results/development/child-neck-v1/child_matter_fraction.json) |
| [matter-antimatter-and-time-note.md](child-and-neck-status-note.md) | Matter and antimatter in the stretched thread, the parent's books, and what the imbalance does for time | [child_lifetime_entropy.py](https://doi.org/10.5281/zenodo.23090015) | [json](../../results/development/child-neck-v1/child_lifetime_entropy.json) [png](../../results/development/child-neck-v1/child_lifetime_entropy.png) |
| [torsion-splitting-and-inherited-spin-note.md](child-and-neck-status-note.md) | The Hehl-Datta splitting at the bounce: how much spin the child needs, and what the parent can give | [torsion_splitting_asymmetry.py](https://doi.org/10.5281/zenodo.23090015) | [json](../../results/development/child-neck-v1/torsion_splitting_asymmetry.json) |
| [second-quantised-splitting-note.md](child-and-neck-status-note.md) | Does the torsion splitting survive second quantisation? The vector-axial term | [second_quantised_splitting.py](https://doi.org/10.5281/zenodo.23090015) | [json](https://doi.org/10.5281/zenodo.23090015) |
| [nesting-thread-length-note.md](child-and-neck-status-note.md) | What sets the child's thread length, and why it is not the lever for nesting | [nesting_thread_length.py](https://doi.org/10.5281/zenodo.23090015) | [json](https://doi.org/10.5281/zenodo.23090015) [png](https://doi.org/10.5281/zenodo.23090015) |
| [collapse-threshold-and-expansion-note.md](child-and-neck-status-note.md) | What sets the collapse threshold, and why the Sun does not expand | none saved | none |
| [junction-neck-note.md](child-and-neck-status-note.md) | The child glued to the parent's exterior: the neck, its opening, and what caps the mass ratio | [junction_neck.py](https://doi.org/10.5281/zenodo.23090015) | [json](../../results/development/child-neck-v1/junction_neck.json) [png](../../results/development/child-neck-v1/junction_neck.png) |
| [israel-matching-note.md](child-and-neck-status-note.md) | Full Israel matching of the closed child to the parent's Schwarzschild exterior | [israel_matching.py](https://doi.org/10.5281/zenodo.23090015) | [json](../../results/development/child-neck-v1/israel_matching.json) [png](../../results/development/child-neck-v1/israel_matching.png) |
| [torsion-layer-note.md](child-and-neck-status-note.md) | Can torsion matter above the bounce density supply the neck's negative-energy layer? | [torsion_layer_eos.py](https://doi.org/10.5281/zenodo.23090015) | [json](../../results/development/child-neck-v1/torsion_layer_eos.json) [png](../../results/development/child-neck-v1/torsion_layer_eos.png) |
| [casimir-layer-note.md](child-and-neck-status-note.md) | Dirac vacuum layer at the neck, and the relation between a layer's pressure/energy ratio and its energy growth | [casimir_layer_check.py](https://doi.org/10.5281/zenodo.23090015) | [json](../../results/development/child-neck-v1/casimir_layer_check.json) [png](../../results/development/child-neck-v1/casimir_layer_check.png) |
| [composite-layer-note.md](https://doi.org/10.5281/zenodo.23090015) | Composite neck layer: Dirac vacuum layer plus torsion matter, matched to the thread path (superseded; the header names the successor) | [composite_layer.py](https://doi.org/10.5281/zenodo.23090015) | [json](https://doi.org/10.5281/zenodo.23090015) [png](https://doi.org/10.5281/zenodo.23090015) |
| [balanced-layers-note.md](child-and-neck-status-note.md) | Normal-stress balance of the neck layer, and the layer whose own ratio is -1/2 | [balanced_layers.py](https://doi.org/10.5281/zenodo.23090015) | [json](https://doi.org/10.5281/zenodo.23090015) [png](https://doi.org/10.5281/zenodo.23090015) |
| [neck-flux-note.md](child-and-neck-status-note.md) | What the child gives to the neck: energy, quanta and net baryon number swept up | [neck_flux.py](https://doi.org/10.5281/zenodo.23090015) | [json](../../results/development/child-neck-v1/neck_flux.json) [png](../../results/development/child-neck-v1/neck_flux.png) |
| [swap-symmetry-note.md](child-and-neck-status-note.md) | Which swaps of the throat keep the Dirac action, and what they do to matter and helicity | [swap_symmetry.py](https://doi.org/10.5281/zenodo.23090015) | [json](https://doi.org/10.5281/zenodo.23090015) [png](https://doi.org/10.5281/zenodo.23090015) |

Historical child-neck scripts are preserved in the [archived evidence](https://doi.org/10.5281/zenodo.23090015). Retained records and figures are in [results/development/child-neck-v1/](../../results/development/child-neck-v1/).

## Historical script layout

The following command describes the archived layout. Each script reads and writes its JSON in the working directory and imports its siblings, using the records folder and its matching scripts folder:

```sh
cd lab/results/development/child-neck-v1
PYTHONPATH=../../../scripts/child_neck python -B ../../../scripts/child_neck/swap_symmetry.py
```

On PowerShell set `$env:PYTHONPATH` first. The shipped JSON files are the inputs the scripts read, so any one script runs on its own. All 19 scripts ran from this layout in about 100 s total on 16 cores, and 18 of 18 JSON records regenerate to within a relative `1e-9` of the saved ones. They were run with Python 3.13 and numpy, scipy, sympy, mpmath, numba and matplotlib. They have not been run under `.venv/validation`. The launcher `scripts/lab.py` sets the working directory to `lab/`, so use the command above rather than the launcher.

## Notes on provenance

- The figures are the versions made during the analysis session. Five scripts also draw plainer figures of their own (`child_expansion_ladder`, `child_lifetime_entropy`, `junction_neck`, `kasner_shear_production`, `nesting_thread_length`), which overwrite the shipped PNG of the same name when run.
- `swap_symmetry.py` gained a `__main__` block when it was installed here, so it writes its own record. The analysis the block runs is unchanged.
- File names beginning `eck` are as first saved; the scripts refer to them by those names.
- Line endings are LF. Notes that were first saved on Windows had CRLF.
- `collapse-threshold-and-expansion-note.md` has no saved script.

## Studies used

Recorded in [docs/arxiv-studies.md](../arxiv-studies.md): arXiv 1110.5019 and 1101.4012 (Popławski) and hep-th/0507253 (Freidel, Minic, Takeuchi).
