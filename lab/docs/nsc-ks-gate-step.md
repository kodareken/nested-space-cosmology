# Gate step ledger

One owner records steps in the declared `(w, U)` class. The merit definition
is written before any row is appended. Existing rows are not rewritten.
A linear prediction is not a measured merit.

```sh
python3 scripts/derive_nsc_ks_gate_step.py --declare
python3 scripts/derive_nsc_ks_gate_step.py --check
```
