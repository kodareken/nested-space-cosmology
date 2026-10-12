# Deciding bounds are not spent against a 1e-4 leftover

UV tail, full between-node remainder, field accuracy and low/subgap remain
`None`. The freed-jet unclipped leftover is still about `1.5e-4`. Those
missing bounds cannot create EXISTENCE and are not large enough, on any
owned estimate, to turn the leftover into a necessary gap. Geometry
between-node remainder stays `1e-13`. No new field or source evolution is
started.

```sh
python scripts/derive_nsc_ks_deciding_bounds_after_leftover.py --check
```
