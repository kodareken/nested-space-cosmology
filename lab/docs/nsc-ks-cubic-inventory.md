# Original-inventory assembly of the shared cubic basis

The loader checks the original 33-channel census and 60 nonzero-angular signed
families against the Cauchy, cutoff-bridge and source-inventory records.
Groups 0, 13, 23 have zero angular label and are absent from this cubic response;
that does not erase their baseline uncertainty. Groups 10, 11, 12, 31, 32 use cutoff
320; the other nonzero-angular groups use 160. Each angular family contributes
its recorded multiplicity divided by 2*pi. The two characteristic signs remain
explicit: no additional energy folding factor is applied.

For channel powers (ell^2,ell^4,m^2 ell^2), raw moments sum those powers with
mu/(2*pi). Leading-integral moments additionally divide each family's term by
2*E_c^2 before summation. Using one common cutoff would be wrong. The N action
coefficient has the existing minus sign; beta has the plus sign. NM already
contains its surface mass term, which is not added again.

Assembly requires both characteristic bases to share the profile, upstream
slice and incoming-coordinate box. Missing values or changed labels fail.
The result is a formal cubic coefficient and its leading-integral expression,
not a complete UV tail: C_M remains null and no physical residual is modified.
Full incoming-interval coverage remains separate.

Three tests check the actual census, the group14 normalization against
(60/pi,300/pi,15*pi), cutoff/multiplicity mutations, both signs and matching
domains. Independent review is pending. No production record has been emitted.

Owners: `src/recursive_horizons/nsc_ks_cubic_inventory.py` and
`tests/test_nsc_ks_cubic_inventory.py`.
