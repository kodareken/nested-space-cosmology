# Local incoming-gate status after L(E) and n=32 geometry leftover

`L(E)` flips sign on the five accepted histories. Baseline `L` is a
nearly constant `-0.05712`, but geometry, matter and edge overwhelm it.
n=32 geometry leftover cannot take a usable step (`condition ~ 2.36e16`)
and truncated SVD leftovers stay about `10^{-4}`. No family was evolved.

State law, class and interval are unchanged. OPEN is not Done.

```sh
python scripts/derive_nsc_local_incoming_gate_status_after_n32.py --check
```
