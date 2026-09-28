# Geometry between-node remainder on I

The local incoming residual is known on the owned 47-node solve/verification
union. A nodal maximum is not `sup_I`. This owner evaluates the
local/reference geometry term on 129 collocation nodes inside `I` and compares
it to the barycentric interpolant of the same term on the owned 47 nodes.

That difference encloses the geometry interpolation remainder. The
matter/source remainder is not evaluated on the dense nodes, so the full
between-node bound, UV tail, field error and low/subgap entries stay `None`.
The physical local incoming gate remains OPEN.

```sh
python scripts/derive_nsc_ks_between_node_geometry.py --check
```
