# The obstruction is a pole, not stiffness

The previous commit called the stop at the inner horizon "stiffness". That
was an assertion. It is measured here, and the assertion was wrong.

## What the equation actually does

The characteristic right-hand side divides by the two speeds:

```
d_plus  = ( ... ) / speed_plus
d_minus = ( ... ) / speed_minus
```

and `speed_plus` vanishes wherever the shift equals the lapse — which is at
**both** horizons, not only at `r_plus`. So `d_plus` carries a simple pole at
`r_minus` as well.

The code clears the pole at `r_plus` properly: `_horizon_jet` evaluates the
spinor and its derivative *on* the horizon without dividing by the speed that
vanishes there. It does not do the same at `r_minus`, and so the inward
branch integrates straight into the pole.

## Both poles are on the outgoing leg

| mu | `r_minus` | `r_plus` | outgoing zeroes | ingoing zeroes |
|---|---|---|---|---|
| 1.002 | 0.938723 | 1.065277 | `0.938723`, `1.065277` | **0** |
| 1.010 | 0.868226 | 1.151774 | `0.868226`, `1.151774` | **0** |
| 1.040 | 0.754343 | 1.325657 | `0.754343`, `1.325657` | **0** |

`speed_minus` never vanishes — it is `-(N + β)`, which is `−2` at both
horizons and nowhere near zero. So the entire obstruction lives on the
outgoing leg, and there are exactly two of them.

## The pole is simple, and it is measured

For `mu = 1.04`, approaching `r_minus = 0.754343` from above:

| distance above | `speed_plus` | `1 / speed_plus` |
|---|---|---|
| `1e-1` | `-3.178e-02` | `-3.15e+01` |
| `1e-2` | `-4.792e-03` | `-2.09e+02` |
| `1e-3` | `-4.996e-04` | `-2.00e+03` |
| `1e-4` | `-5.014e-05` | `-1.99e+04` |

Each decade closer gives a decade larger reciprocal: a simple pole, with a
coefficient of order two for `mu = 1.04`. It is stronger when the horizons
are closer together — about `-1.4` for `mu = 1.002` and `-0.53` for
`mu = 1.010`.

That is why the solver reports a step-size failure. The solution is
unbounded at `r_minus`, so an adaptive integrator asks for ever smaller
steps and finally runs out of float64 room.

## What this corrects

"Stiffness" was the wrong word and I should not have used it without
checking. The solver is not struggling with a badly conditioned but bounded
right-hand side. It is walking into an unbounded one, at a radius I can name.

## What this suggests, without claiming it

The pole at `r_plus` was cleared by evaluating the jet on the horizon instead
of dividing by the vanishing speed. **The same construction applies at
`r_minus`**, and nothing in the geometry says it cannot be done. Whether it
is worth doing depends on what is on the other side, which is still unknown.

That is a concrete proposal rather than a vague hope, and it is a smaller
piece of work than anything else currently on the table.

## Scope

A construction check on the reference chart. **This does not measure the
actual universe.** Locating the pole says where the obstruction is. It says
nothing about what is on the other side of it.

Three values of the mass ratio, `40000` samples per entry, no station
payloads, no marching.

## Provenance

- reduction: `lab/scripts/derive_nsc_rn_pole_structure.py`
- record: `lab/results/development/nsc-rn-pole-structure-v1.json`
- tests: `lab/tests/test_nsc_rn_pole_structure.py`
