# Uniform-energy polynomial preflight for the vacuum field

The scalar-phase change `G -> G+(E/B)I` gives diagonal entries `(0,2E/B)`
and leaves the vacuum projector unchanged. With `E=28+4x`, represent the
new approximate column by a Chebyshev polynomial of degree24. The physical
equation's residual has degree25; its last coefficient is retained.

## Reuse and fixed stopping decision

Reuse the [certified initializer](nsc-incoming-matched-horizon-initializer.md)
and [fixed-energy Taylor residual machinery](nsc-incoming-validated-vacuum-propagation.md).
The missing connection is a whole-real-energy error estimate without many
independent field solves. First certify the initial degree24 interpolant,
verify exact phase/energy identities, and take **one** degree36 time step
of size1/32 for both signs. Share the geometry coefficients and keep analytic
time radius1/8. Stop against the local uniform-column gate `1e-20`; record
actual runtime and extrapolated work before any full-panel decision.

Risks are silently discarding the degree25 energy residual, normalizing an
endpoint polynomial into a rational function, and using complex-energy
growth assumptions for an unevaluated propagated field. None is needed.

## Initial interpolation, with a separate projector error

The explicit matched initial ratio is rational in energy. Its normalization
continues using coefficient conjugates at the same complex E. On the
Bernstein ellipse `rho=8`, the enclosing energy rectangle is
`Re(E) in [11.75,44.25]`, `abs(Im(E))<=15.75`. A directed bound on the ratio
excludes zeros of `1+w(E)*wbar_coeff(E)` and fixes the normalization branch.
Only the initial formula is continued to complex energy.

The degree24 interpolant is built at25 Chebyshev-Lobatto nodes using a
directed DCT. The exact constant column `e1` is separated before bounding
`v-e1`. For each complex component, the imported bound is

$$
\|f-p_{24}\|_\infty\le\frac{4M\rho^{-24}}{\rho-1}.
$$

This is [Trefethen, ATAP Theorem8.2, Eq8.3](https://raw.githubusercontent.com/chebfun/ATAP/development/chap8.m).
The vector error is the Euclidean combination of the two component bounds.
Directed coefficient centering adds the sum of coefficient-radius norms.
The inherited mathematical affine-projector error stays **separate** from
this column interpolation error. Frozen binary mass/angular labels are used
in both the initializer and propagated operator; source kappa is unchanged.

## Continuous residual over real energy

For `vhat=sum c_k(y) T_k(x)`, the exact Chebyshev multiplication rules include
`x T0=T1` and `x Tk=(T(k-1)+T(k+1))/2` for k>=1. Projecting time-recursion
coefficients to degree24 defines the approximation; it does not remove the
degree25 contribution from the actual equation. Every time-residual
coefficient retains all26 energy coefficients. Since `abs(Tk)<=1`,

$$
\sup_{x\in[-1,1]}\|R(y,x)\|\le\sum_{k=0}^{25}\|R_k(y)\|.
$$

Pointwise unitarity in real E integrates this uniform residual without a
Gronwall factor or a propagated complex-energy bound. The existing directed
Cauchy coefficients enclose the geometry. Factoring them into `1/B`,
`r sqrt(delta/B)` and `sqrt(delta/B)` lets both angular signs share one
geometry calculation. The analytic time-generator tail and all coefficient
and step-centering errors remain charged separately.

## Endpoint algebra and unexecuted future work

The raw polynomial covariance has degree48; a linear metric vertex gives
degree49. Gauss48 would integrate that polynomial part exactly. Renormalizing
the endpoint polynomial is prohibited in this argument because it would
make the covariance rational. Its norm defect must remain within the
certified covariance error. Ad4 source subtraction needs its own directed
spectral quadrature/arithmetic account.

The [preflight receipt](../results/development/nsc-incoming-energy-panel-preflight.json)
reports measured initializer/shared-geometry/step times, the finite operation
count for699 future steps, a linear runtime extrapolation and a2x planning
allowance. Those timings are planning estimates, not wall-clock guarantees.
No full-panel propagation or physical source integral has been performed.
All earlier field receipts, modes and source values remain intact.

```sh
python3 scripts/derive_nsc_incoming_energy_panel_preflight.py --write
python3 scripts/derive_nsc_incoming_energy_panel_preflight.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_energy_panel.py
```
