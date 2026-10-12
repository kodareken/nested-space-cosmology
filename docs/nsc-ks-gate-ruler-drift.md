# Ruler drift of one history under two tangent sets

`b4b53f0b…` was evolved twice from the same upstream source. Iterate6
carried 64 tangent directions. The n=64 columns froze that history and
carried 64 further modes at amplitude zero. The summed matter changes at
the 47 nodes differ by about `5.78e-5` in N and `1.28e-10` in beta.

That difference is an indicator. It is not an error bound, and it is not
subtracted from the residual. N reported below this scale on the present
grid, degree and step control is not reproducible.

```sh
python3 scripts/derive_nsc_ks_gate_ruler_drift.py --record
python3 scripts/derive_nsc_ks_gate_ruler_drift.py --check
```
