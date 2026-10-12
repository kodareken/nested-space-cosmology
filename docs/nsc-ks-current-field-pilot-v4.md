# Current-history local-Fourier field pilot, v4

This successor tests the missing continuous field proof method on the immutable
family `14_1`, DOP853 cell 122 capture for history `0b0e4ced...`. It does not
evolve a new field or source.

The v3 proof retained Fourier modes obtained from a sampled transform and used
a global derivative-L1 inequality for every omitted coefficient. That inequality
is valid, but the flat cutoff ramps make it many orders too large. Version 4
integrates each retained coefficient directly with Arb and encloses the
exponentially flat endpoint strips. It then bounds the analytic profile minus
that retained Fourier polynomial on local panels using midpoint Taylor jets and
a directed derivative remainder. No fitted decay or sampled maximum enters the
bound.

The resulting coefficients and uniform profile-tail bounds feed the unchanged
`DifferenceResidualPolynomial`. The history-minus-reference residual is formed
before norms. Time-Taylor and reciprocal-radius remainders remain separately
added through `difference_operator_remainder_bounds`.

The executable owner is
`scripts/validate_nsc_ks_current_field_pilot_v4.py`; the reusable method is
`src/recursive_horizons/nsc_ks_local_profile_fourier_bound.py`. The record
reports the directed one-cell values. Even a small one-cell value is not a
whole-cone field error: all accepted time cells, representative source families,
the Gronwall transport and the source error remain required before the method
can enter the nine-component gate budget.

On the declared family `14_1` cell, the local-Fourier method encloses the
continuous normalized residual at approximately
`(3.00849e-16, 1.36374e-12)` for axial orders zero and one. The corresponding
v3 enclosure was approximately `(5.481e-4, 6.429e-1)`. The finite retained
band contributes approximately `(1.34233e-16, 6.09639e-13)`; the separate
time/radius remainder is approximately `(1.13638e-23, 2.71047e-20)`.

This resolves the named *one-cell method* obstruction caused by the global
cutoff derivative L1. It does not show that the sum over every accepted time
cell and every source family fits the field allocation. Whole-cone capture,
Gronwall propagation, contraction and source error remain the next decision.

Reproduce with:

```sh
PYTHONPATH=src .venv/validation/bin/python scripts/validate_nsc_ks_current_field_pilot_v4.py --check
```

`--check` verifies the immutable record and its dependencies. `--replay`
recomputes the directed integrals and local panels before comparing the stable
scientific fields.
