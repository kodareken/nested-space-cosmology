# Replayable witness for one subgap row

This records the already bounded original `group14/low16_1` row 15 and its
archived opposite-energy, opposite-angular partner. It does not derive a new
action, source law, or cosmological identity. The local gate stays **OPEN**.

The phase transport is the completed mixed-transport owner. The selection is
the one-row pilot: $E=0.24867511687395624$, $m=\pi/2$, $\ell=+\sqrt{5}$,
order 8, collar $10^{-10}$, `rtol=2e-13`, `atol=2e-15`, degree 16, 48 metric
terms, 4 defect subdivisions, tube $10^{-5}$, and inner request
`horizon_rho+1.01e-4`. The Bloch comparison uses the existing horizon frame,
`reflection_from_phase`, `initial_bloch`, and `validate_bloch` at the archived
`rho_up`. Covariance errors are unweighted. No $N/\beta$ aggregate is formed.

Owners:

- `src/recursive_horizons/nsc_subgap_row_witness.py`
- `scripts/derive_nsc_subgap_row15_upstream.py`

Planned outputs, not written by the tests:

- `results/development/nsc-subgap-row15-upstream-v1.json`
- `results/development/artifacts/nsc-subgap-row15-upstream-v1.npz`

## Phase prefix

Each transported DOP853 cell is stored as eleven float64 columns:
`y0`, `y1`, `theta0`, `theta1`, and all seven phase coefficients `F0` through
`F6`. The prefix must be finite, strictly inward, contiguous in both node and
angle, and satisfy `theta0+F0 == theta1` exactly. The original log-amplitude
component is omitted. The subgap phase proof does not read it.

Replay rebuilds `Dop853DenseOutput` interpolants and a SciPy `OdeSolution`,
then a `run.sol` closure whose only free variable is `dense`.
`original_dense_solution` already requires that contract, so the mixed-transport
owner runs unchanged and does not solve again. A direct list of phase segments
would duplicate the owner's cell loop. The closure is the smaller adapter:
capture checks that the rebuilt segments match the original 843 anchored cells
exactly, including the last kept cell, whose endpoint agrees with `theta0+F0`
and the following solver node.

The mode namespace carries the archived energy, mass, angular momentum,
background, and outer radius. A substituted label is rejected before transport.

## Signed source

The NPZ stores the exact positive and negative `2x3` columns and `3x3`
covariances. The negative partner is the archived complement, not a second
Jost solve. Its error is recomputed with the sign flip
$(n_x,n_y,n_z)\mapsto(n_x,-n_y,-n_z)$. Copying the positive scalar is rejected.
The quadrature weight is recorded and not applied.

`--check` revalidates the saved phase prefix and Bloch trace. It does not call
`solve_jost` or `capture_bloch`. CPU and wall time are printed only; they are
not part of the equality record. The frame errors are included in the initial
Bloch enclosure; no separate horizon-sewing budget number is added. Physical
$\rho=1$ and full upstream-budget entries remain null because this record
covers only this pair at the archived upstream radius.

| Quantity | Observed upper |
|---|---|
| Phase cells (contracting / expanding) | 843 (820 / 23) |
| Inner phase error | $1.814146078\times 10^{-9}$ |
| Bloch cells | 671 |
| Positive covariance error | $1.047331588\times 10^{-10}$ |
| Negative covariance error | $1.046705871\times 10^{-10}$ |

Zero and straddling rates did not occur on this row. One bounded row is not a
source campaign. The other rows, the other positive angular family, quadrature,
and contraction into the changed-minus-reference source remain open.

```sh
.venv/validation/bin/python scripts/lab.py -m pytest -q tests/test_nsc_subgap_row_witness.py
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_subgap_row15_upstream.py --check
```

`--check` reads the planned files after `--record`. The tests publish only
under a temporary directory.
