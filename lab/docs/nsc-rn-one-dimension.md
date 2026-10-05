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

So everything measured in this line — the shell, its amplitude, its
wavelength, the factor-of-two compression, the three-sided structure, the
floor — is **one-dimensional radial structure**. The three sides are a radial
interior, surface and exterior. They are not a sphere with an interior, and
the picture cannot be made one by moving a wave to the side, because there is
no side to move it into.

## What this does not say

It does not say transverse structure is impossible in the physics being
modelled. It says this state vector, this operator and this chart have one
degree of freedom, and a construction with more would be a different one
rather than a wider version of this one.

That is the same kind of boundary as the parent/child cut: the construction
asserts a nesting and does not supply the physics for it. Here it asserts a
sphere through `r^2 dOmega^2` and supplies only a radius.

## Scope

A construction check on the code as written. **It does not measure the
actual universe** and it does not test the idea of transverse propagation.
It establishes what the instrument can and cannot do, which is what any
result read off it has to be bounded by.

## Provenance

- `lab/src/recursive_horizons/nsc_rn_response.py` line 106
- `lab/src/recursive_horizons/nsc_rn_observables.py` line 387
- station payloads under `lab/results/development/`
