# How far in the machinery actually reaches

Three gates stood between the marching code and the region below the inner
horizon. All three were refusals by function, not by geometry. All three are
now overridable and **all three are off by default**, so a run that does not
ask for the change behaves exactly as before.

| gate | before | now |
|---|---|---|
| `vacuum_chart` | refuses radii at or below `r_minus` | `allow_below_inner_horizon` extends to `Q^2/(2M)` |
| `build_grid` | requires `r_minus < excision` | same flag moves the lower bound to `Q^2/(2M)` |
| `prepare_radial_packet` | packet must sit above `r_minus` | same flag |
| `scattering_mode` | grid must stay above `r_minus` | same flag |

With all four open, how far in does the packet actually integrate?

| lower edge | result |
|---|---|
| `1.1000` | ok |
| `0.9500` | ok |
| `0.8500` | ok |
| `0.8000` | ok |
| `0.7700` | ok |
| `0.7550` | ok |
| `0.7000` | **fails: required step size is less than spacing between numbers** |
| `0.6000` | **fails, same** |

## The boundary is the inner horizon, and the reason is numerical

For `mu = 1.04` the inner horizon is at `r_minus = 0.754343`. The
integrator reaches `0.7550`, which is `0.0007` **above** it, and fails at
`0.70`, which is below it.

So the construction's machinery reaches precisely to the inner horizon and
no further — and the reason is not a convention. The characteristic ODE
stiffens there, and float64 cannot take the step the solver asks for. That
is a different kind of obstruction from the three guards, and a more
interesting one: the guards said "we did not look," and this says "the
integrator stops."

This is also the fourth correction of the same shape in this line. Three
"infinite" walls turned out to be unmeasured conventions — memory scaling as
`N` squared, resolution, the chart's own guard. This one is real. It is a
stiff boundary, not a choice.

## What it does not say

Whether the third region contains anything. Whether it is static in a useful
sense. Whether matter crossing the inner horizon outward does anything. A
stiffer integrator, a different variable, or a different chart might cross;
nothing here shows that any of them would.

That is now a well-posed question rather than a vague one: the obstruction
is localised, it is at a known radius, and it has a known cause.

## Scope

Construction checks on finite marched collars. **This does not measure the
actual universe.**

Four default-off guards were added to two production modules. The curated
suite is unchanged at `368 passed, 1 skipped`, which is the evidence that
nothing moved for a caller that does not ask for the change. Tests pin both
that the guards are still present and that the defaults still refuse.

## Provenance

- `lab/src/recursive_horizons/nsc_rn_source.py`: `vacuum_chart`,
  `scattering_mode`, `prepare_radial_packet`
- `lab/src/recursive_horizons/nsc_rn_pg.py`: `build_grid`
