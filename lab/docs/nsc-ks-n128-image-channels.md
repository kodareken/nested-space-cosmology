# Sup-norm image and channels at the n=128 truncation

The best measured point remains `ad759424…`. Its zero-pad lift is
`3af275b0…`. This register does not evolve a history.

Geometry on the declared interval is compared with the barycentric
interpolant of the same term on the owned 47-node union. The weights are
the explicit products in node order. SciPy's barycentric evaluator is not
repeatable on this union, so it is not the recorded remainder. The raw
maximum and the `nextafter(raw+1e-14)` enclosure are the same for the
64-coefficient history and the zero-pad. Both components exceed `1e-13`,
so the enclosure is not the older `1e-13` quote.

UV, the full residual between-node bound, field error and low/subgap stay
`None`. The fine-matter owner is bound to a different axial profile. The
linear UV contraction has no numerical tail on `I`. The group13 subgap
panel is not a constraint error at this history. Dense matter is not stored.

The geometry principal determinant on these slots excludes zero and keeps
the certified sign. That is invertibility, not EXISTENCE. The equilibrated
coupled Jacobian has 94 positive singular values. Newton rank is 90 because
the absolute floor is `1e-14` while the leading value is below 1. The absent
step is the fat shape, 256 unknowns against 94 equations, not a certified
left null. The stable block has ratio `1e6` and whitens to 66 columns.

Inside singular ratio `1e6`, a beta-primary sup-norm screen holds N at its
current maximum. The predicted image still misses the `3e-11` ball. A
joint N gain remains positive when beta must fall by `1e-8`, and it does
not close the N gap. `L(E)` and the flux combination are scored on this
residual only. A sign change, or a gap at one history, is not scoped
NON-EXISTENCE. n=256 is not authorized.

```sh
python3 scripts/derive_nsc_ks_n128_image_channels.py --write
python3 scripts/derive_nsc_ks_n128_image_channels.py --check
```
