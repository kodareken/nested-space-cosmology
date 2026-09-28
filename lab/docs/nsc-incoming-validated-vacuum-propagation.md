# One-frequency vacuum propagation with a continuous residual bound

For the fixed canonical field equation `v_y=-i G(y) v`, exact unitarity gives

$$
\|v_{\rm exact}(y_1)-\widehat v(y_1)\|
\le e_{\rm initial}+\sum_j\int_{I_j}
\|\widehat v'_j+iG\widehat v_j\|\,dy+e_{\rm centering}+e_{\rm endpoints}.
$$

The [receipt](../results/development/nsc-incoming-validated-vacuum-propagation.json)
applies this reused unitary-defect estimate to group32 at **E=24 only**, for
both angular signs. It includes initial column representation, continuous
polynomial residuals, coefficient arithmetic, step centering and the exact
incoming target coordinate. No spectral integral follows from a one-frequency
field result.

## Reuse, missing connection and fixed gate

Reuse the [certified matched vacuum initializer](nsc-incoming-matched-horizon-initializer.md),
the unchanged fixed profile and the unitary residual theorem already used
by the incoming vacuum-tail owner. The missing connection is a bounded field
solution error between the finite initializer and the incoming rho1 slice.
Perform one fixed run at Taylor degree36, real step at most1/32, analytic
radius1/8 and80-digit directed arithmetic. The pointwise weighted lapse-kernel
allocation is `1e-13`. Stop with PASS or an explicit obstruction at this
resolution; do not infer an energy-window result or escalate resolution
automatically. At most two CPU workers handle the two angular signs.

The likely obstacles are geometry cancellation near the horizon, interval
state wrapping, and hiding numerical residual behind solver tolerances.
Stable exact-profile formulas, centered approximate columns with a separate
norm-error radius, and a whole-step residual enclosure address these directly.

## Exact geometry and physical coordinates

Use `y=log(delta)` with true exact-profile horizon distance. The initial
distance is the inherited `1e-10`; the target is
`delta1=q_h_exact-3*pi/4`, the existing rho1 surface. The coordinate range is
fixed by these geometric endpoints, not an invented propagation duration.
Stored binary mass and signed angular labels are used explicitly in the
Hamiltonian. The source occupation kappa stays unchanged.

The generator is evaluated through `B=W(q_h-delta)/delta`. The exact identity

$$
B=b_0+(-3\sin2q_h-\cos2q_h)\delta\,\operatorname{sinc}(\delta)^2
+(-3\cos2q_h+\sin2q_h)\,[\operatorname{sinc}(2\delta)-1]
$$

removes the small-delta subtraction. Directed sinc series retain degree36
with an explicit geometric tail bound. Each complex y disk must have
positive real lower bounds for B and `sin(q_h-delta)`. These margins exclude
zeros and fix the square-root branch throughout the disk. The old truncated
`interior_W` is not substituted for the exact profile.

## Directed coefficients and continuous residuals

At each real step center,128 directed circle evaluations on radius1/16
produce generator Taylor coefficients through degree36. A directed radix-two
DFT implements the Cauchy coefficient formula. A norm bound M on the larger
radius R=1/8 disk encloses its analytic alias tail:

$$
|a_j-\widehat a_j|\le\frac{M}{R^j}
\frac{(r/R)^{128}}{1-(r/R)^{128}}.
$$

The approximate column coefficients are midpoints of the directed recurrence
`(n+1)c[n+1]=-i sum G[j] c[n-j]`. They are recorded mathematical point
coefficients; their uncertainty is charged to the residual instead of
propagating an interval state. The actual polynomial residual coefficients,
including all low orders, are enclosed through degree72. No low-order
coefficient is silently set to zero.

The omitted analytic generator tail contributes at most

$$
\frac{M\,\sup\|\widehat v\|\,h}{38}
\frac{(h/R)^{37}}{1-h/R}
$$

to the step residual integral. Polynomial coefficient-norm bounds integrate
as `sum R_n h^(n+1)/(n+1)`. Directed Horner evaluation encloses each step
endpoint; choosing its next center adds an explicit rounding radius. Exact
unitarity sums these independent error radii without Gronwall growth or
interval state wrapping. Checkpoints retain completed centers and residual
certificates so a resumed run evaluates only missing steps.

## Initial state, units and scope

The new initializer's binary ratio defines a normalized mathematical column.
Directed normalization supplies its representation radius. Its separately
certified affine-vacuum projector discrepancy is carried once. For final
column error e, the raw outer-product covariance error is bounded by
`(2+e)e`, plus that initial projector discrepancy. The raw column norm defect
is retained and checked; no archived column is normalized or modified.

The final error is converted with the inherited signed-family factor to a
lapse action-gradient **spectral kernel per positive-energy dE**. The two
actual signs are summed once. This is a pointwise bound at E24, not the
error in a source integral, and includes no thermal source approximation.
There is no metric timestep, physical IV choice or stationarity claim.

Replay contracts saved continuous residual envelopes, geometry margins and
endpoint accounting without running field steps. The fixed-resolution proof
can support a subsequent energy-panel calculation only when its own uniform
field and quadrature errors are established.

```sh
python3 scripts/derive_nsc_incoming_validated_vacuum_propagation.py --prepare --workers 2
python3 scripts/derive_nsc_incoming_validated_vacuum_propagation.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_validated_vacuum_propagation.py
```
