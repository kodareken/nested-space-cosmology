# Physical rho=1 source-bound successor, v2

This successor starts from the immutable v1 OPEN coverage owner. It traces
the stored `rho=1` three-column fibers in the authenticated source inventory,
splits the 59 straddling low/subgap panels at each V5 joined threshold, and
builds row and panel manifests for the complete 164-panel / 49372-row
column-accuracy universe.

Stored construction is not restated as a count. Low-mode panels are the
accepted `mode_at_one` recoveries restricted by `restrict_resolved_modes`.
Middle panels are the original archived order-16 producer used by the
incoming baseline. Subgap panels are the eight inherited real rows at
`rho=1`. Axial source values on this slice are the PDE restriction
`F_z = G(1) F` of those stored columns. The Pauli 2-norm of `G(1)` supplies
a directed Lipschitz map from a coherent Frobenius fiber error, original
quadrature weights, signed multiplicities, channel degeneracies and the
owned E/-E antiunitary pairing into N,beta majorants for
`baseline_low_subgap` and the physical upstream slot.

The 62 V5 action-covered middle panels remain source-accuracy OPEN. Their
action certificate is not imported as column error. Analytic ell-zero
radius-history response on groups 0, 13 and 23 is reused; those groups'
baseline allocations stay. A finite N,beta component is emitted only when
every inventory identity carries an enclosure payload.

The missing primitive is a directed a-posteriori enclosure of the stored
`rho=1` fibers, and of `G(1)F`, against the exact physical incoming modes
for each construction method. SVD amplification, ODE tolerances and the
historical `4.998e-9` aggregate are not that enclosure. This slice does not
launch those generators. A resource forecast records the 29052 low-mode,
20016 middle and 304 subgap rows that still need construction enclosures.
A partial or resource stop remains OPEN and does not close Gate 2.

```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python \
  scripts/derive_nsc_ks_rho1_source_bound_v2.py --check
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python \
  -m pytest -q tests/test_nsc_ks_rho1_source_bound_v2.py
```
