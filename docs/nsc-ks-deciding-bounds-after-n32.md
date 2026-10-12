# Deciding bounds stay None after the n=32 leftover floor

`L(E)` is not a locked gap. The n=32 geometry leftover stays on the current
residual scale. UV, full between-node, field and low/subgap remain `None`.
Geometry between-node remainder stays `1e-13`. No new field or source
evolution is started.

```sh
python scripts/derive_nsc_ks_deciding_bounds_after_n32.py --check
```
