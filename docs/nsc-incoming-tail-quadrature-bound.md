# Directed numerical certificate for the archived group22 infinite tail

Reuse the existing order16/ad4 paired tail, its exact `c1/c2/ad1` endpoint
cancellation, the directed coefficient and Gauss-root owners, and the Cauchy
quadrature estimate. The missing quantity is numerical integration error
relative to the archived order13 float, including its truncation and rounding.
The existing radial vacuum-error certificate remains unchanged.

One pilot covers group22 above its actual endpoint160. Compactify with
`x=1/E`, giving the same integral of `source(x)/x^2` on `[0,1/160]`.
Use one48-node,80-digit directed Gauss rule, with circle center `1/320`,
radius `1/160` and32 directed arcs. Stop after this pilot and its replay:
PASS if the lapse-error bound is below `3e-11`, otherwise retain OPEN and
identify the analytic or numerical obstruction. No old Taylor, radial,
mode, scattering or source-baseline calculation is rerun.

## Exact cancellation and analytic continuation

The vacuum-minus-ad0 difference reuses the stable normalized-projector
difference formula. Its base series is `S0=U*x/(1+sqrt(1+y))`, where
`U=lambda/r+i*m` and `y=a^2*(m^2+lambda^2/r^2)*x^2`. The difference from
the16-term polynomial begins at degree2. An exact algebraic identity checks
this rearrangement before any directed evaluation. The owned `c1/c2/ad1`
identity then supplies `source=x^2*g`, `g(0)=0`; no numerical coefficient
is clipped or inferred from near-zero samples.

The remaining adiabatic terms cannot add a pole. A small integer-valuation
certificate propagates the existing cross-product/normalization recurrence:
`h_transverse=O(1)`, `h_z=O(x^-1)`, and inverse gap squared `O(x^2)`.
The leading `b0_z=1+O(x^2)` has a q-independent constant. The resulting
transverse/longitudinal powers for ad1 through ad4 are `(2,3)`, `(3,4)`,
`(4,5)`, `(5,6)`. Thus ad2–ad4 contribute powers at least1 to the reduced
source `source/x^2`. Only the ad1 constant requires the imported exact
cancellation. No homogeneous source integration is repeated.

The existing generated ad4 function is reused. A fail-closed AST whitelist
replaces only its energy-gap square root and four odd half powers
with the root `E/sqrt(C)*sqrt(1+C*B/E^2)`, where the original generated
assignments define `C` and `B`. This is the branch analytic at inverse energy
zero and equal to the original positive-energy branch. The original
assignments and replacement counts are checked exactly. All other generated
functions remain unchanged; exact half-integer constants retain their guard.
Energy-dependent integer powers use exact multiplication, with the reciprocal
base taken before exponentiation. Division by a power is rewritten as
multiplication by that inverse power. This avoids a diagnosed interval-wrapper
failure in the original circle pilot: `1/(gap_root^9)` enclosed zero after
raising a wide rectangle, despite valid disk pole margins. A localized arc
control verifies this algebra-preserving correction at the same resolution.

The resulting expression is meromorphic with at most a finite-order pole
at `x=0`. The imported exact leading cancellations remove that pole.
Directed triangle bounds exclude both normalization zeros and gap-root zeros
on the whole disk. Complex continuation conjugates coefficients at the same
x, never conjugates x itself. The circle bounds and positive exact Gauss
rule give the already-owned `4*h*M*(h/R)^(2*n)/(1-h/R)` error estimate.

The [record](../results/development/nsc-incoming-tail-quadrature-bound.json)
preserves the archived order13 value and encloses the same exact numerical
recipe directly. Its distance to the directed integral bounds archived
truncation, coefficient rounding, summation and integration together.
The existing elementary thermal formula is evaluated outward and stored with
lossless interval endpoints so its nonzero exponentially small bound survives
floating underflow. The physical order16 defect is imported, not recomputed.

```sh
python3 scripts/derive_nsc_incoming_tail_quadrature_bound.py --prepare
python3 scripts/derive_nsc_incoming_tail_quadrature_bound.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_tail_quadrature_bound.py
```
