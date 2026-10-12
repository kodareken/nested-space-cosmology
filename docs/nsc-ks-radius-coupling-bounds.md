# Continuous radius-coupling integrals

This closes two previously unowned inputs of the reference/difference field
error method: `M_integral` and `Mz_integral`. It uses the existing radius
history, generator and exact rational profile bounds. It introduces no source,
new physical law or final gate claim.

Write `r = r0 + delta_r`, with `r0` homogeneous in z. The actual generator
in `nsc_ks_difference_envelope._radius_multipliers` gives

\[
M=-\frac{i\ell}{a}\frac{\delta r}{r_0r}\sigma_2,
\qquad M_z=-\frac{i\ell}{a}\frac{\delta r_z}{r^2}\sigma_2.
\]

The row-sum norm of the Pauli matrix is one. In the upstream half slab,
`-sigma <= s <= 0`, the existing relation `|d rho| = a |ds|` removes the
factor `1/a`. Since `|chi| <= 1`, if `W_j,U_j` are global derivative bounds,

\[
J_j=\frac{\sigma^2W_j}{2}+\frac{\sigma^4U_j}{24},\qquad
\int\|M\|\,|d\rho|\le\frac{|\ell|J_0}{r_{0,\min}r_{\min}},
\qquad
\int\|M_z\|\,|d\rho|\le\frac{|\ell|J_1}{r_{\min}^2}.
\]

These are continuous analytical bounds over the full compact upstream support,
including cutoff ramps, not sampled maxima or convergence estimates. The
reference and actual radius lowers come from the existing owned slab proof.
All polynomial arithmetic is rational; conversion to Arb supplies an outward
upper for the existing propagation interface. The small symbolic control
independently differentiates the exact coupling and integrates the polynomial
majorant. Zero history and the massless zero-angular sector give exact zero.

The whole-cone campaign now fills these two slots, checking its geometry and
angular channel before doing so. Physical preparation error, reference residual,
full-time difference residual and continuous contraction norms remain separate
obligations. The historical cone-v1 record keeps its original bytes and bindings.

```sh
PYTHONPATH=src .venv/validation/bin/python -m pytest -q tests/test_nsc_ks_radius_coupling_bounds.py tests/test_nsc_ks_current_field_campaign.py
PYTHONPATH=src .venv/validation/bin/python scripts/derive_nsc_ks_radius_coupling_bounds.py --check
```
