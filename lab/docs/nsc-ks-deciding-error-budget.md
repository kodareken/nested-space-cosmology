# Deciding local-gate error leftovers at the fifth clipped history

Geometry interpolation on `I` is enclosed at `1e-13` for the accepted
fifth clipped history. The four leftovers that can still decide the gate
stay `None`:

- changed-history UV tail at the actual 160/320 cutoffs
- full between-node residual remainder
- final-history field / evaluation accuracy
- remaining low/subgap source bound

The Pauli/subgap matrix identity is not a source-error bound. The residual
Chebyshev tail is a resolution diagnostic, not a matter remainder. Missing
bounds stay `None`. They cannot create EXISTENCE while the n=16 leftover
stays about `10^{-4}`.

```sh
python scripts/derive_nsc_ks_deciding_error_budget.py --check
```
