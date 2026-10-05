# One propagation direction: the state has no transverse degree of freedom

Moving a wave to the side to get a sphere needs a second dimension for the
oscillation to run in. That is what was claimed about this construction and
it has now been checked rather than inferred from the line element.

Three independent confirmations, in the code rather than in a note.

**The angular index is required to be exactly one.** The response tangent
raises unless the state is `(2, points, 1)`:

```
if array.ndim != 3 or array.shape[0] != 2 or array.shape[2] != 1:
    raise ValueError("response tangent expects the rank-one column (2, points, 1)")
```

**The angle is a hardcoded constant.** The connection is evaluated once at
`theta = 0.5 * math.pi` and never varied:

```
387:    theta = 0.5 * math.pi
```

**The recorded state has no angular axis at all.** Every dynamical array is
one-dimensional in radius:

```
coupled_phi  (2, 869, 1)      two spinor components, one angular point
coupled_mass (869,)
radius       (869,)
frozen_lapse (869,)
frozen_shift (869,)
R4           (869,)
```

## What this settles

There is no transverse mode because there is **no slot for one**. The
dynamics are radial: `lambda_pm = -beta +- N` along `partial_r`, and the
angles enter only through the passive factor `r^2 dOmega^2` in the line
element.

The arrays are radial. That is a fact about the state vector. It is not a
fact about the geometry.

## Withdrawn: "this is not a sphere"

The sentence that the radial interior, surface and exterior "are not a sphere
with an interior", and that the chart "asserts a sphere through `r^2 dOmega^2`
and supplies only a radius", is withdrawn. It read a symmetry reduction as a
missing dimension.

Two dimensions give the circle. The inside of that circle is the third
dimension. The circle plus that inside is the ball. In this chart the circle
at areal radius `r` is the factor `4 pi r^2` already in the force integral,
and the inside is smaller `r`. Front and back are the same circle, because
the symmetry sets them equal. They are not two calculations that were never
run.

The gradient balance `B(r)` is that walk into the ball. It matches the
sphere. The earlier reading that it fails to match the sphere is the
withdrawn sentence.

What the missing angular slot still forbids is a front that differs from the
back. That is a different object. It is not required for the ball.

## What this does not say

It does not say transverse structure is impossible in the physics being
modelled. It says this state vector has one radial degree of freedom, and
that this degree of freedom is the interior of the circle.

## Scope

A construction check on the code as written, with one withdrawn reading.
**It does not measure the actual universe.** The arrays are radial. The
geometry they reduce is the ball: circle, plus the inside as the third
dimension.

## Provenance

- `lab/src/recursive_horizons/nsc_rn_response.py` line 106
- `lab/src/recursive_horizons/nsc_rn_observables.py` line 387
- station payloads under `lab/results/development/`
