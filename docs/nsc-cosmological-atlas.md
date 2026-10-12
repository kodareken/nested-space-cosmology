# Finite exterior-first cosmological atlas

The atlas assembles possible exteriors, interfaces and descendants before
attempting a local dark-sector fit. It is a map of mathematical candidates
and missing physical connections, not a claimed solution of the universe.

The default motif has one unresolved parent, two off-centre black-hole
interfaces, two potentially unequal children and one potential grandchild.
Its finite descent cutoff is two. Each interface has a separate imported
Schwarzschild exterior chart. No common two-hole parent metric has been
supplied. A closed FLRW child is an ansatz with unknown source and scale
factor. S2 symmetry orbits do not imply a global S3 spatial slice, and the
illustration is a patch diagram rather than a three-dimensional embedding.

The three alternatives are one-hole multiple daughters, a conditional
three-mouth traversable construction, and a finite response-equivalence
control. Their assumptions remain separate. Multiple daughters inside one
black hole are distinct from separate holes producing different children.

## Owners and reproducible output

- [Model](../src/recursive_horizons/nsc_cosmological_atlas.py): region/link
  records, finite validation, imported junction/clock constraints, fixed-theory
  check, causal graph bounds, and explicitly supplied Hamiltonian adapter.
- [Producer](../scripts/derive_nsc_cosmological_atlas.py): creation-only JSON,
  SVG patch maps, PNG previews through the existing matplotlib, illustrated
  narrative and portable source/output hashes. Previews use the SVG's geometry
  and text, with a flat equivalent of its decorative background gradient.
- [Tests](../tests/test_nsc_cosmological_atlas.py): incompatible junctions,
  unequal children, clock exits, causal delays/disconnection, unit changes,
  finite depth/cycles, fixed-theory perturbations and invisible exteriors.
- [Illustrated narrative](../results/development/nsc-cosmological-atlas-v1/narrative.md)
  and [machine-readable atlas](../results/development/nsc-cosmological-atlas-v1/atlas.json).

From the repository root, with the existing validation interpreter:

```powershell
.venv/validation/Scripts/python.exe scripts/lab.py scripts/derive_nsc_cosmological_atlas.py
.venv/validation/Scripts/python.exe scripts/lab.py scripts/derive_nsc_cosmological_atlas.py --check
.venv/validation/Scripts/python.exe scripts/lab.py -m pytest tests/test_nsc_cosmological_atlas.py -q
```

At the initial atlas verification on Windows, the validation redirector
reported `failed to locate pyvenv.cfg: Access is denied`. The installed base
Python 3.14.8 ran with the same validation site-packages and explicit lab
source paths, without installing anything. Since `scripts/lab.py` replaces
`PYTHONPATH`, that fallback invokes this standalone producer directly:

```powershell
$atlasRoot = (Get-Location).Path
$env:PYTHONPATH = "$atlasRoot/lab/src;$atlasRoot/lab/scripts;$atlasRoot/.venv/validation/Lib/site-packages"
& "$env:LOCALAPPDATA/Programs/Python/Python314/python.exe" lab/scripts/derive_nsc_cosmological_atlas.py --check
& "$env:LOCALAPPDATA/Programs/Python/Python314/python.exe" -m pytest lab/tests/test_nsc_cosmological_atlas.py -q
```

The producer refuses an existing output directory. A failure after claiming
the directory leaves an inspectable partial bundle, without a success manifest;
use a new versioned output directory after correction. Verification is read-only.
Source hashes include transitive local Python imports, the producer, test,
this note and local notes used as evidence. Paths are repository-relative;
no host path, clock timestamp or
unverified historical revision supplies authority. Byte hashes authenticate
outputs; numerical replay uses relative 1e-10 and absolute 1e-12 tolerances.

## Interfaces and physical meaning

`build_atlas()` returns independent scenario records. Every link has one of
four kinds, a source and an evidence statement: descent, attachment, causal,
or correlation. An unresolved link stays unresolved. The validated descent
graph is acyclic and bounded. Explicit `contained_in` interfaces retain their
host's declared generation, so moving through an interface cannot reset the
lineage count. Arbitrary attachments and causal networks may have cycles.

`angular_junction()` checks equal induced areal radius and real timelike
normal, and evaluates the imported angular Israel condition. It does not
solve temporal extrinsic curvature, shell matter or conservation. The
small-neck orientation exposes its negative surface-density requirement.
`null_collision_residual()` is only the imported collision-sphere condition.

`static_clock()` is a static Schwarzschild observer's elapsed proper time
within a finite declared interval. It rejects horizons and domain exits.
`earliest_arrival()` composes only explicitly declared positive-delay causal
edges sharing a named clock and finite interval. It gives a conditional graph
bound, not a geodesic solution. Physical candidate links presently lack the
complete numeric clock maps, so they do not acquire numerical arrival times.
Synthetic delay fixtures test this distinction. Entanglement, common ancestry
and correlation never become signal routes by themselves.

`require_fixed_theory()` permits state changes while rejecting changes to a
candidate's theory data and declared regional bindings of fixed couplings.
This is a record-level guard, not proof that a given state perturbation
satisfies its field equations.

`evolve_supplied_hamiltonian()` accepts a finite Hermitian matrix or callable
`H(t)` (dimension inferred at the first finite sample and validated on every
call), uniform finite time grid, retained indices, covariance and explicit origin. It reuses
`nsc_evolving_reduction.evolve_retained_region`; it neither infers coefficients
from the diagram nor derives a relativistic signal cone. The saved four-mode
integration is an adapter fixture. Its unequal occupations exercise different
children without assigning cosmological values to the coefficients.

The response-equivalence control uses known Schur algebra: changing an exactly
uncoupled exterior mode changes the exterior spectrum but not the local
resolvent. The covariance control uses the same visible state. This supplies
a concrete limit on exterior reconstruction, not a new inverse theorem. Its
resolvent lives in the upper half-plane, with poles excluded in boundary
values. A direct hidden-mode spectral probe would distinguish the examples;
adding a coupling is a separate model and changes the local response problem.

## Borrowed sources and remaining connections

The complete versioned URLs and reading boundaries live in the generated
atlas source records. Key anchors are [Barrabes–Frolov](https://arxiv.org/html/hep-th/9511136v1)
for one-hole branching and null shells, [Emparan et al.](https://arxiv.org/html/2012.07821v2)
for a controlled connected-mouth alternative, and [Poplawski v2](https://arxiv.org/html/2007.11556v2)
for the anisotropic spin-fluid hypothesis and its explicit geometric limitation.
The exact Kantowski–Sachs metric is not a finished route to the required S3
geometry. The connected-mouth source domain cannot be silently replaced with
the exploratory one-Dirac-field neutral-child source.
The one-hole l>a disconnection threshold specifically uses the
near-singularity approximation, equal nucleation radius, and a<<R0<<2m in
equations 9-12; it is not a universal separation rule for sibling universes.

The new scientific work would supply a common parent geometry; physical
interface sources and full matching; an admissible joint state; interpatch
clock conversions; and a causal exterior-response prediction. Dimensionless
inherited laws, environmental dressing and state choices must stay distinct.
Only after those connections exist does a local dark-sector, antimatter or
time discriminator become a meaningful target. There is no abundance fit here.

## Admission and limits

The source, domain and state compatibility map is the admission category.
Each calculation records six fields: question, borrowed mathematics, missing
connection, decision and unlock, finite budget, and stop condition.
Junction identities, clocks, finite Schur elimination and the numerical
adapter are imported mathematics or software controls. Neither their replay
nor their tests are a new physical result. This atlas does not establish a
globally solved multi-black-hole spacetime, stable child formation, inherited
constant values, sibling signaling, or identification of our observed cosmos
with any displayed child. No existing paused numerical campaign is resumed.
