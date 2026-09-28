# Current-history whole-cone field core, v1

This successor owns the reusable whole-cone field API. It restores the
immutable v4 local-Fourier coefficient payload by hash, streams
difference-reconstruction cells, and sums directed residual bounds. The
history-minus-reference polynomial is formed before norms. It does not
replace the v4 one-cell record and it does not run family `14_1`.

The bounded polynomial is the normalized defect `S=X_xi-h L X`. Its cell
bound is already the contribution to `integral |R| d rho`; the accumulator
therefore sums it directly and never multiplies it by the cell width again.

The executable owner is `scripts/validate_nsc_ks_current_field_cone_v1.py`.
The reusable API is `src/recursive_horizons/nsc_ks_current_field_cone.py`.

## Bindings

| Object | Role |
|---|---|
| `WholeConeFieldConfig` | Family, profile identity, source columns, period, settings, radius tail, and the v4 payload digest |
| `WholeConeAccumulator` | Streamed `accumulate_segment`; segments and field balls are not retained |
| `WholeConeCheckpoint` | Exact prefix totals and dense-segment hashes; every prefix segment must replay in order before new cells are accepted |
| `WholeConeFieldResult` | Directed polynomial, remainder and residual-integral sums |
| `finalize_family` | Adapters to `propagate_difference_error` and `difference_matter_error` |

Physical `rho=1` source error remains unresolved. Resource exhaustion,
method failure and missing source error stay failures; they are not
converted into `NON_EXISTENCE` or `PASS`. Runtime and memory observations
are excluded from the scientific digest.

On the saved family-`14_1` cell 122, restoring the v4 payload and
evaluating the same difference residual keeps the recorded v4 enclosure
inside the restored directed total. That comparison is one cell. Coverage
of every accepted time cell and every source family remains open.

Reproduce with:

```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/validate_nsc_ks_current_field_cone_v1.py --check
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_nsc_ks_current_field_cone.py tests/test_nsc_ks_local_profile_fourier_bound.py
```
