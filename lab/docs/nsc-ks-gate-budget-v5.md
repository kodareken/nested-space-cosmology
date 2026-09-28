# Gate error budget v5

The 2026-09-28 Mac integration writes the corrected v5-schema successor to
`results/development/nsc-ks-gate-budget-v5-review1.json`, using the reviewed
continuous-enclosure-v2 record. The original v5 and v4 bytes remain unchanged.
The narrower per-component allocations are provisional; the known covered-region
bound already exceeds its slot. Reconcile using measured pilots and the fixed
total reserve, rather than silently replacing previously owned uncertainties.


Successor of the [v4 nine-component budget](nsc-ks-gate-budget-independent.md)
for the switched law `C_Sigma[g]=U_g C_up U_g^dagger` on `I=S(1)+[.12,.18]`.
The nine named components and their allocations are frozen:

| Component | Allocation |
|---|---|
| field_space_time | `5e-12` |
| changed_history_UV_tail | `5e-12` |
| baseline_low_subgap | `3e-12` |
| upstream | `2e-12` |
| energy_interpolation | `1e-12` |
| covered_regions | `1e-12` |
| phase_value | `1e-12` |
| between_node | `1e-12` |
| arithmetic | `1e-12` |

The allocation total is `2e-11`. At least `1e-11` remains for the physical
residual. No allocation may be moved silently. A component is filled only from
an authenticated applicable PASS enclosure for the current profile, switched
state law, and 257-node grid. Evidence path bytes must match `sha256`.
Diagnostic, sample, ENCLOSED, REUSE, OPEN, 129-node phase, missing evidence,
or a fabricated digest do not fill v5.

Predecessor v4 bytes stay immutable. Its three enclosed rows are recorded as
predecessor evidence. They are not silent v5 PASS fills. The current record
therefore remains OPEN with null component bounds. The integers 60 and 120
are the declared family requirement; they are not measured coverage. Covered
family counts stay null until the source inventory is authenticated.
Admission of a closed v5 budget additionally requires all nine finite, a
directed total at most `2e-11` componentwise, authenticated 60/120 source
coverage, and the residual reserve. Caller-written coverage flags are not
coverage. An OPEN or incomplete budget cannot authorize the search driver and
is not EXISTENCE or NON-EXISTENCE.

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_ks_gate_budget_v5.py
python3 scripts/derive_nsc_ks_gate_budget_v5.py --check
```
