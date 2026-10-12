# Original upstream covariance and its signed partner

This successor carries the [source covariance proof](nsc-subgap-source-covariance.md)
to the exact upstream slice used by the production evaluator. The saved upstream
columns are restored through `RetainedUpstreamArchive`, with their original
source and preparation digests. No preparation ODE is repeated and no saved
columns are replaced. The separate Bloch curve is a proof witness.

For the same original group14/low16_1 row0, the exact endpoint is

$$y_{\rm up}=\log\bigl(q_h-\pi/2-\arctan\rho_{\rm up}\bigr).$$

The earlier continuous defect method encloses the complete positive covariance
there and compares it with $A_{\rm up}C_{\rm src}A_{\rm up}^\dagger$ from the
actual archived columns. Thus this comparison includes their numerical
continuation error; it is not a bound on a newly substituted preparation.
The quadrature weight is retained in the record but is not applied to this
single-energy covariance norm. Physical source contraction applies it later.

## Negative energy requires its actual covariance

The exact homogeneous interior fundamental matrix is unitary. The complete
subgap sewing has Gram identity because its reflection has unit modulus.
Consequently, on this preparation slice only,

$$Q_- = I-\sigma_3\overline{Q_+}\sigma_3,$$

with the opposite angular sign. In Bloch coordinates this is
$(n_x,n_y,n_z)\mapsto(n_x,-n_y,-n_z)$, which preserves the error norm.
The same validated ball can therefore be compared directly with the actual
negative numerical $A_-C_-A_-^\dagger$. This comparison includes the numerical
Gram and covariance-complement discrepancies. It does not simply copy the
positive scalar error bound.

For example, take $A=\operatorname{diag}(1,2)$ with a third zero column and
$C_+=\operatorname{diag}(1,0,0)$. The positive covariance can be exact while
the complementary negative covariance has operator error 3. A regression
requires the signed bound to detect this. The exact Gram identity above is
not assumed pointwise after inhomogeneous evolution.

## Scope and replay

The record covers one positive energy/angular mode and its original
opposite-energy/opposite-angular partner. The other positive angular family,
remaining rows, quadrature remainder and effective-source contraction are not
covered. The physical upstream budget component remains null and the local
gate remains OPEN.

Owners are `src/recursive_horizons/nsc_subgap_upstream_covariance.py` and
`scripts/derive_nsc_subgap_upstream_covariance.py`. The v1 JSON binds the
authenticated original archive and a compact NPZ witness.
Replay consumes that saved witness and revalidates its continuous defect;
it does not ask the numerical integrator to reproduce the same curve.
The initial time must match the enclosed horizon initializer, and the source
arrays in the witness must equal the authenticated archived arrays. Payload
mutation and shifted-start controls exercise these bindings.

```sh
.venv/validation/bin/python scripts/lab.py -m pytest -q tests/test_nsc_subgap_upstream_covariance.py
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_subgap_upstream_covariance.py --check
```
