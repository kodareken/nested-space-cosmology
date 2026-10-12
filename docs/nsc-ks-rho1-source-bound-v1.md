# Physical rho=1 source-bound core, v1

This slice owns the coverage ledger and validator for physical reconstruction
error of the stored `rho=1` source columns. It traces the authenticated V5
inventory, the retained source payload, the low/subgap successor, and the
upstream continuation owner whose physical `rho=1` input remains missing.

The low/subgap classifier selects 102 panel-level windows and 29356 rows:
64 low-mode panels and 38 subgap panels. That count is not the count of rows
below the V5 joined thresholds. Authenticated row inspection finds 27828 rows
below, 1528 rows at or above, and 59 straddling panels. The 62 middle panels /
20016 rows are V5 action-covered, but an action certificate is not a physical
source-column reconstruction bound. The physical column-accuracy universe is
therefore still all 164 panels / 49372 rows until a successor supplies a proved
bridge. A production campaign must split the straddling windows explicitly.

Group-12 and group-32 order-24 middle windows remain reusable action successors;
they are not silently promoted to source accuracy. Analytic ell-zero
radius-history response on groups 0, 13 and 23 is reused while their baseline
stays. E/-E pairing is the owned antiunitary map with original weights; it is
not an extra multiplicity or a signed-energy fold.

The v1 recorder inspects two remaining panels and at most three coherent
fibers. It does not run a source campaign and cannot close a directed N,beta
budget from a contribution count, samples, zero-filled rows, V5 action values,
or the superseded coarse aggregate near `4.998e-9`. Missing bounds stay `null`.
The advertised `--check` mode rebuilds and compares the stable record.

```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python \
  scripts/derive_nsc_ks_rho1_source_bound_v1.py --check
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python \
  -m pytest -q tests/test_nsc_ks_rho1_source_bound.py
```
