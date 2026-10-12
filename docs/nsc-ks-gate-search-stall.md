# Value-only search stall

The one-knob floor moves with the solver by about `8e-4`. That is
numerical error, so the trust-region search is not started. No family
is evolved. The stall is not NON-EXISTENCE. The next owner is the same
declared `(w, U)` class on `I`, after the ruler is fixed.

```sh
python3 scripts/derive_nsc_ks_gate_search_stall.py --record
python3 scripts/derive_nsc_ks_gate_search_stall.py --check
```
