# Low/subgap source decision

The existing V5 action bound is reused only on its authenticated spectral
regions. The exact zero radius-history response of groups 0, 13 and 23 is
also reused, while their nonzero baseline remains. Neither fact bounds the
reconstruction error of the retained physical source columns.

The prior conditional coarse lapse enclosure is about `4.998e-9`, dominated
by groups 12 and 32 with group 14 still incomplete. It is far above the
approximately `1.1001e-11` N headroom remaining after the three known budget
components. The approved next method is therefore a window-specific source
reconstruction bound for groups 12/32 and completion of group 14. The coarse
enclosure is not scaled unchanged, and Pauli, signed-source or nested-node
identities are not relabelled as source-error bounds.

```sh
python3 scripts/derive_nsc_ks_low_subgap_decision.py --check
```
