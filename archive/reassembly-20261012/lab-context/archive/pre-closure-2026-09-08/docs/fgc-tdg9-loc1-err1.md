# FGC-1-TDG9-LOC1-ERR1: independent-root enclosure diagnosis

`FGC-1-TDG9-LOC1-ERR1` is a compact, premise-only diagnosis bound to sealed
LOC1 authority commit `14e26c43edbc88003d68cbcf24681bb090d98534`.

LOC1 stopped before publication at
`localization/independent_disagreement` for retry `3`, channel `u:alpha`,
component `value_V`, level `D01`. The primary and independent routes found the
same two co-maximizers and exact maximum, but counted `8,176` versus `8,544`
candidates. Exactly `184` tiny value-Hermite cubics acquired `368` false
interior candidates in the independent route.

The first witness factors as

```text
V(t)  = delta t^2 (2t - 3)
V'(t) = 6 delta t (t - 1).
```

Its derivative roots are exactly the excluded endpoints `0` and `1`. The v1
independent evaluator used a fixed absolute square-root enclosure whose scale
was wider than this tiny exact-square discriminant; the transformed intervals
were then clipped to `[0,1/2]` and `[1/2,1]` and misclassified as interior.

The thrown `independent_disagreement` was observed by the operator but never
serialized. No LOC1 manifest or terminal exists. At ERR1 construction, the
LOC1 output and staging namespaces were absent and the current 115-leaf store
snapshot had digest
`5917d70de3080dcc6b7abeb7fdb7cc3370c6140cbeb37bde0c142d6a2d36a445`.
An independent regeneration of immutable PREF2 owner result `7e676cb4...`
matched its sealed store manifest `46b57bc3...`; this binds the construction
snapshot to the authenticated PREF lineage. It is not an unavailable LOC1
before/after comparison. Ordinary ERR1 verification is compact and
store-blind; only its one-time construction consumed live absence, store, and
exact retry-3 reconstruction evidence.

ERR1 diagnoses a numerical evaluator. It is not a LOC1 success, convergence
result, GR-0 calibration, action/model obstruction, candidate result,
mechanism result, or physical result.
