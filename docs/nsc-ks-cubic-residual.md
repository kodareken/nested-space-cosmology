# Cubic leading term in the saved value residual

The [cubic inventory](nsc-ks-cubic-inventory.md) assembles formal moments and
sets `physical_residual_modified` false. This successor is the executable
addition of those **leading** moments to the saved finite-cutoff value
residual. It is not another metadata sibling, not a new field or source
solver, and not a UV certificate.

Owner: `src/recursive_horizons/nsc_ks_cubic_residual.py`.
Driver: `scripts/derive_nsc_ks_cubic_residual.py`.

## What is added

Both characteristic signs stay explicit. For the shared unit-angular channels
\((N_2,N_4,N_M)\) and \((J_2,J_4,J_M)\),

$$
N=-\sum_s\sum_{i=1}^{3} w_i B_{s,i},\qquad
\beta=+\sum_s\sum_{i=1}^{3} w_i C_{s,i}.
$$

The weights \(w_i\) are the inventory **leading** moments, not the raw
moments. Each nonzero-angular family contributes \(\mu/(2\pi)\), and that
family's leading piece divides by \(2E_c^2\) before the sum. Groups 10, 11,
12, 31 and 32 use \(E_c=320\); every other nonzero-angular group uses
\(E_c=160\). The energy folding factor is 1. \(N_M\) already contains the
surface mass term, so that term is not added again.

The saved value residual is the replayed nodal sum of baseline, geometry,
matter and edge from `nsc-ks-gate-value-v2-current-dop-g1024-d48-metadata2`,
history `0b0e4ced`, 129 exact target nodes. Its own baseline is still the
original length-2 center pair inside `assemble_residual`. The cubic term is
applied only after that nodal residual exists. A length-2 cubic pair, a
one-entry channel, or one `z` domain repeated across the nodes is rejected.
Missing nodes are left empty; they are not filled from the center.

Approximate residual changes are binary64 midpoints. Enclosure endpoints are
stored separately. A midpoint is not a bound.

## What blocks certification

`C_M`, uniform \(C_4\) on \(I\), and the complete UV tail stay null. The
higher remainder was not constructed. `certify` therefore raises. The
physical local gate stays OPEN. This file does not change source-row coverage
and does not run a field or source evolution.

The demonstration marches the existing shared channel basis at 1024 cells and
Taylor order 8, on exact saved target nodes, both characteristic signs, until
the 180-second CPU budget. A smaller cell count is not substituted. Nodes
that do not finish inside the budget are omitted.

```sh
/Users/admin/Documents/BlackHoles-Infinity/.venv/validation/bin/python \
  /Users/admin/Documents/BlackHoles-Infinity-Worktrees/grok-uv-residual-20260929/scripts/lab.py \
  -m pytest -q tests/test_nsc_ks_cubic_residual.py tests/test_nsc_ks_cubic_inventory.py
/Users/admin/Documents/BlackHoles-Infinity/.venv/validation/bin/python \
  /Users/admin/Documents/BlackHoles-Infinity-Worktrees/grok-uv-residual-20260929/scripts/lab.py \
  scripts/derive_nsc_ks_cubic_residual.py --check
```
