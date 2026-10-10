# Atlas v4: the first sourced interior and a finite observer domain

The [v4 producer](../scripts/derive_nsc_cosmological_atlas_v4.py) preserves the
six prior scenarios and adds an explicit cap/wall/bridge/AE initial slice.
[Slice mathematics](nsc-atlas-connection-slice.md),
[end-only gluing](nsc-atlas-connection-gluing.md) and
[observer domain](nsc-atlas-connection-observer.md) own the new assumptions,
admissions, equations, tests and reading boundaries. No prior source-bound
file is changed. The [working paper v4](../../paper/atlas-working-draft-v4.md)
and the [bundle](../results/development/nsc-cosmological-atlas-v4/narrative.md)
provide the illustrated narrative.

The concrete two-body existence application uses the sourced body A with
M_A/ell=0.4 and a complete vacuum Schwarzschild bridge B with M_B/ell=0.6.
Only their selected smooth vacuum ends are inputs to the published many-body
gluing theorem. Exact open overlaps reattach their unmodified remainders.
There is one designated common outer parent end and B's original other end;
there are two total AE ends. B's other end is not a second established child.
The common connector is an existence family at sufficiently large separation,
not an explicitly solved metric or a certified value of the old BL distance.

The observer at R_obs/ell=1.6 has a protected local vacuum domain if the
assembly retains the right initial exterior through R_cut/ell=4 and an open
overlap beyond it. Its proper-time bound is an open endpoint. Equivalent
initial data imply the same local vacuum geometry in that domain even for
different distant assemblies. No global shell evolution theorem is invoked.

## Reproduction and evidence

From the repository root with the scientific environment:

```sh
python scripts/lab.py scripts/derive_nsc_cosmological_atlas_v4.py --check
python scripts/lab.py -m pytest tests/test_nsc_atlas_connection_slice.py tests/test_nsc_atlas_connection_observer.py tests/test_nsc_cosmological_atlas_v4.py -q
```

Without --check, generation is creation-only. The manifest binds all inherited
v1/v2/v3 sources plus this version's adapters, tests, notes, producer and paper.
The atlas record labels historical assessments superseded by this application;
historical scenario records remain unchanged. Imported branching, communicating
mouths, response equivalence and the radiation-wall exclusion remain separate.

Checks compare independent surface stress projections, metric derivatives,
bulk constraints, bridge regularity, unit changes, null quadrature and Kruskal
coordinates. The plot R(l) is an intrinsic radial profile, not a Euclidean
embedding: |R'| exceeds one on part of the slice. Finite samples end at R=6ell
while the exact analytic end extends to infinity. No finite numerical sample
certifies arbitrary asymptotic behaviour; that continuation is the imported
Schwarzschild solution.

## What this unlocks

The previously missing body-to-vacuum-end interface is supplied in one fixed
model. Its parent assembly is justified at initial-data existence level. The
next physical question is the causal completion and stability of this sourced
body, or a discriminating field/state observable beyond the demonstrated
vacuum degeneracy. First borrow the published completion and obstruction
results, then admit only a named unresolved finite check. Parent formation,
offspring production, parameter inheritance, antimatter and local dark-sector
signatures remain separate open questions. No new-to-world claim is made.
