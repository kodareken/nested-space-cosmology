# Incoming source with adjacent-window refinements

The updated incoming matter source adds four explicit disjoint increments
to the previous source-update record: group12 LOW`[32,40]`, group22 middle
`[40,160]`, and group32 LOW`[32,40]` plus middle`[40,320]`. The original
receipts remain intact. The group32 union is a budget check, not a third
source contribution.

Reuse the existing source-update normalization, local/reference contributions,
certified energy coefficients and new source-window records. The missing
connection is their composition into the same incoming constraint equation.
Check interval coverage and explicit before/after identities, add each
increment once, and propagate the partial budget. Stop at this composition;
do not rerun source, radial, reference or local generators.

The main risks are double counting a shared boundary/union, applying a bound
without its matching source correction, and promoting a component budget to
complete source accuracy. The owner rejects overlaps of positive measure,
requires source-order and certificate agreement, and leaves the complete
source error null. Endpoints shared by adjacent integration cells have zero
measure and introduce no physical cutoff or boundary term.

The [record](../results/development/nsc-incoming-source-update-v2.json)
contains the updated source and incoming lapse/shift approximants, explicit
increments and residuals. The physical thermal contributions are bounded,
not assumed absent. The remaining middle, low/subgap and rigorous quadrature
errors still precede a certified numerical constraint root. No physical IV,
parent history, coupling, seed, action term or metric timestep is selected.

```sh
python3 scripts/derive_nsc_incoming_source_update_v2.py --write
python3 scripts/derive_nsc_incoming_source_update_v2.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_source_update_v2.py
```
