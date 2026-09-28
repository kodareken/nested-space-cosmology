# Continuous endpoint norms and signed insertion error

The field campaign now connects its completed numerical field error to the
two matter contractions. This fills a missing method, not the physical source
error, full-family coverage or the gate certificate.

For each source column let the captured endpoint difference interpolant be

\[
D_{a j}(z)=\sum_k d_{a j k}e^{ik(z-z_0)},\qquad
\Delta F_{a j}=e^{-iE_jz}D_{a j}.
\]

Arb's DFT encloses the Fourier coefficients of the supplied binary endpoint
values. With the original weights already multiplying those coefficients,
define

\[
b_a=\sum_k\left(\sum_j|d_{a j k}|^2\right)^{1/2},\qquad
c_a=\sum_k\left(\sum_j|(k-E_j)d_{a j k}|^2\right)^{1/2}.
\]

Then uniformly over the entire numerical period,
`||delta F||_F <= sqrt(sum_a b_a^2)` and
`||delta F_z||_F <= sqrt(sum_a c_a^2)`. This is the same Fourier triangle
estimate as the existing `endpoint_norms` owner, applied directly to the
anchored endpoint without rebuilding eight time coefficients. The actual
envelope derivative is retained. A nodal maximum cannot replace this bound.

`complete_endpoint_inputs` joins those norms and the homogeneous reference
norms to the four already-propagated field errors. Geometry denominators use
lower dyadic endpoints of the owned exact rational lower bounds, rather than
the older generic upper-slot conversion. The campaign uses the
authenticated terminal segment only after the complete time interval has
been enclosed. Missing covariance or multiplicity keeps the contraction open.

Positive and negative energy contributions use their own source covariances.
The existing opposite-angular antiunitary map preserves field/error norms;
it does not make the covariance norms equal. `contract_signed_pair` bounds
the sectors separately and adds their N and beta uppers with directed
arithmetic. Each sector retains the original per-angular-sign multiplicity.
No extra fold or duplicated source weight is introduced.

The recorded pilot evaluates the norms at the endpoint of the archived
family-14 cell 122. That endpoint is not rho=1. The pilot is a check of the
continuous norm method on real source data, not an incoming constraint or a
whole-time error bound. The tests additionally exercise a peak between grid
nodes, the Nyquist derivative, unequal coherent signed covariances and direct
perturbations of the existing constraint contraction.

```sh
PYTHONPATH=src .venv/validation/bin/python -m pytest -q tests/test_nsc_ks_endpoint_contraction.py tests/test_nsc_ks_current_field_campaign.py
PYTHONPATH=src .venv/validation/bin/python scripts/derive_nsc_ks_endpoint_norm_pilot.py --check
```
