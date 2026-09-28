# Freed-jet n=16 linear proposal, not evolved

Dropping the declared solver jet makes the 32-direction rectangular system
square and kills the solve-node leftover. The same unclipped step still
leaves all-node maxima about `(1.52e-4, 7.77e-6)`, and the unit clip stays
on the same leftover scale as the current residual. That leftover cannot
reach `3e-11`, so no family is evolved. The physical local incoming gate
remains OPEN.

```sh
python scripts/derive_nsc_ks_n16_freed_jet_proposal.py --check
```
