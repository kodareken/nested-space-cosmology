# Directed low-energy affine-vacuum projector

## Reuse, decision and stopping condition

The target is `Re Pv01 > 1/4` at rho1 for group14, positive angular sign,
and the fixed binary energy0.5562120090641313. Pv is the dominant affine
vacuum component of the existing coherent source, not a replacement state.
The separate thermal/coherent remainder retains the original covariance and
reflection-phase dependence. No archived column is altered or normalized.

Reuse the exact profile root enclosure, the stable sinc-minus-one identity,
lossless directed interval serialization and the unitary defect estimate.
The guarded group32 high-energy propagator is unchanged; only its generic
sinc helper is imported. One2048-step run at60 directed decimal digits is
allowed, with120 CPU seconds and no adaptive refinement. Stop PASS only if
the final lower enclosure exceeds1/4; otherwise report OPEN. This is a
coarse sign predicate, not a3e-11 source-precision calculation.

## Generator and omitted affine tail

Use delta=qh_exact−q, y=log(delta), with qh enclosed from the existing exact
profile root. For real delta,

\[
B=W/\delta=b_0+\mathrm{linear}\,\delta(1+\operatorname{sincm1}\delta)^2
 +\mathrm{sincCoeff}\operatorname{sincm1}(2\delta),
\]

where `b0=2*kappa_geometry`, `linear=−3*sin(2qh)−cos(2qh)`, and
`sincCoeff=−3*cos(2qh)+sin(2qh)`. The reused complex sinc enclosure is projected
to its real part because the exact function is real here. Geometry positivity
is checked before division. Let `r=1/sin(qh−delta)` and

\[
H=\begin{pmatrix}-D&\bar p\\p&D\end{pmatrix},\qquad
D=E/B,\quad p=\sqrt\delta(-mr+i\ell)/\sqrt B.
\]

The mass/angular intervals contain both their stored binary labels and
pi/2,sqrt5. E and source occupation parameters are unchanged binary values.
Every frozen matrix is built from three exact real dyadic numbers, so it
is exactly Hermitian before its exponential is evaluated.

The numerical initial column is exactly e1. Its y coordinate is the lower
endpoint of directed log(1e−8). The omitted projector tail is bounded by
`2*A*sqrt(dmax)`, where dmax covers both1e−8 and exp(actual y0), and
`|Hoff|<=A*sqrt(delta)` throughout(0,dmax]. This follows from the constant
projector residual `norm([H,diag(1,0)])=|Hoff|` and the unitary defect integral.
It includes starting-coordinate rounding instead of declaring the finite
initializer exact. Source phases and Jost coherence are not changed.

## Point-column algorithm and directed recurrence

For each of2048 consecutive real cells, enclose the actual H on the whole
cell. Select an exact Hermitian Hmid from its midpoint enclosure. If h is
the actual cell length, use

\[
\eta=h\sup\|H-H_{\rm mid}\|,
\qquad
\epsilon_{j+1}=\epsilon_j+\eta\|u_j\|+r_j.
\]

For a traceless Hermitian difference its operator norm is exactly the
Euclidean norm of the three real coordinates D,Re p,Im p; interval ranges
therefore give a conservative supremum. Directed sine, cosine and square
root evaluate `exp(−i*h*Hmid)u_j`. Recenter that rectangle to an exact point
column; its full Euclidean radius is r_j and enters the recurrence. Only
point columns propagate. Unitarity prevents amplification of the previous
error. No wrapped interval state or normalization is used.

The exact solution starting from e1 has norm1, so its projector difference
from the final numerical outer product is at most
`(1+norm(u))*epsilon`. Add the omitted affine tail. The physical rho1 endpoint
`log(qh−3*pi/4)` is an interval; the last numerical point is its midpoint.
A further projector bridge `2*sup(norm(H))*abs(endpoint−point)` encloses that
uncertainty. Form Re(u0*conj(u1)) with interval arithmetic first, then add
all projector error radii and outward-round its final real interval.

The proof artifact stores each cell, real H range, frozen point, column
chain, interval update image, recenter radius and all additive bounds.
Replay checks coverage and recalculates error sums/final predicate without
another propagation. Independent review precedes production. The final
record binds this specific channel/energy/configuration and all inherited
owners; the original coherent source and all old receipts stay unchanged.

The numerical result and all error components are in the
[development receipt](../results/development/nsc-incoming-low-energy-vacuum-enclosure.json).
Its compressed proof preserves lossless interval endpoints and point-column
bits. The verification command replays those error sums only:

```sh
python3 scripts/derive_nsc_incoming_low_energy_vacuum_enclosure.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_low_energy_vacuum_enclosure.py
```
