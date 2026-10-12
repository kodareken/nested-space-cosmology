# Shared parent data and a fixed-source attachment test

The second atlas adds two independent pieces of mathematics to the
[first candidate map](nsc-cosmological-atlas.md): a published two-hole spatial
initial-data set and a finite source-compatibility screen. Their domains remain
separate. A spherical junction test is not a gluing construction inside a
nonspherical two-hole spacetime.

## What the map now contains

The [parent adapter](nsc-atlas-parent-initial-data.md) imports time-symmetric
vacuum Brill–Lindquist data with two unequal puncture parameters. The three
asymptotically flat ends, their ADM masses, and the spatial metric are known.
No lapse, spacetime development, horizon locations, child formation or
cross-region clock map follows merely from plotting that slice. The plot
shows conformal coordinates and a conformal factor, not a physical embedding.
Finite angular variation near each puncture shows why two exact spherical
charts cannot simply be pasted into the common parent.

The [source test](nsc-atlas-junction-sources.md) takes an exact closed-radiation
FLRW history, a comoving spherical boundary, a constant-mass Schwarzschild
exterior and a prescribed constant-tension wall. An algebraic obstruction
excludes this combination on any open interval of changing radius. Samples
of the full junction stress and conservation identities illustrate and
independently check the result; samples are not the proof. Allowing a moving
boundary, changing either bulk, or supplying a different source opens a
different candidate. No wall law is fitted to the desired trajectory.

The original one-hole multiple-daughter, communicating-mouth and
response-equivalence alternatives remain available. Their failure conditions
are independent of this source screen. The finite descent cutoff remains two.
No causal delay is inferred from the spatial parent data or from correlation.

## Owners and evidence

- [Producer](../scripts/derive_nsc_cosmological_atlas_v2.py) combines the
  unchanged v1 scenario model with the two new adapters.
- [Integration tests](../tests/test_nsc_cosmological_atlas_v2.py) check that
  the rejected source is not promoted to a general exclusion and that the
  spatial initial data create neither offspring nor signal routes.
- [Atlas and illustrated narrative](../results/development/nsc-cosmological-atlas-v2/narrative.md)
  link the machine-readable inputs, parent record and junction record.
- [Working manuscript](../../paper/atlas-working-draft.md) explains the
  construction, conditional obstruction, existing response machinery and
  open physical connections. It is an unpublished draft in the companion
  research direction; the released manuscripts remain frozen.

## Admission and limits

The adapters each record borrowed results, the missing connection, the
decision and unlock, a finite domain and budget, stopping conditions, and
independent checks. The parent is imported integration work. The junction
result applies established matching theory to one previously unresolved
atlas candidate; no new-to-world theorem or new physical mechanism is claimed.
The response engine from v1 remains usable only with a supplied Hamiltonian.
These additions supply no physical multi-region Hamiltonian or local
dark-sector prediction.

The figure producer adds a fixed 181 by 101 coordinate grid to the parent
adapter's small controls, masking disks of coordinate radius 0.28 around the
punctures. It renders four figures and the 33-point junction witness with the
existing plotting stack. This is a bounded illustration budget of one CPU
minute per production, not an evolution or an adaptive parameter search.

The evidence is reproducible with the existing scientific Python stack:

```powershell
python scripts/lab.py scripts/derive_nsc_cosmological_atlas_v2.py
python scripts/lab.py scripts/derive_nsc_cosmological_atlas_v2.py --check
```

The producer creates a new directory exclusively and refuses overwriting.
The read-only check authenticates source/output hashes and replays the bounded
calculation. The Windows base-interpreter fallback is documented in the
[v1 owner](nsc-cosmological-atlas.md#owners-and-reproducible-output).
The manifest records dependency versions. JSON numerical replay uses relative
`1e-10` and absolute `1e-12` tolerances; those tolerances are integration
checks, not error bars for a cosmological prediction. Figure hashes bind the
rendered outputs; exact figure replay requires the recorded rendering stack.

The next missing connection is source-supported moving-boundary dynamics or
a different published interior, followed by a justified response and clock
map. The constant-tension/comoving-radiation combination is retired from that
queue. A known de Sitter/domain-wall solution may be imported as a separate
controlled attachment; it does not repair the rejected radiation attachment.
