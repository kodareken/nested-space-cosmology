# What the inner-horizon fix actually requires

The pole structure is measured and the fix is identified. This records what
the fix costs, because the cost is the thing that decides whether it is
worth doing.

## The fix

The outgoing speed vanishes at both horizons. The code clears the pole at
`r_plus` by evaluating the spinor jet on the horizon rather than dividing
by the vanishing speed. **The same construction applies at `r_minus`**, and
`_horizon_jet` has been generalised to take a radius so it can be built
there. That part is done and is harmless: the default still builds it at
`r_plus` and the curated suite is unchanged.

Cutting the inward branch at `r_minus` and restarting it from that jet is
the part that is not done.

## Why it was not done

Implementing the split surfaced a second requirement that the first hides.
The characteristic right-hand side calls `_point_chart`, and `_point_chart`
calls `vacuum_chart` **without the flag**. So even with a correct branch
split, the integrator's right-hand side still refuses radii below
`r_minus`, and it does so from inside the solver's step loop.

The chain that has to be threaded is:

```
scattering_mode  ->  _integrate_branch  ->  rhs closure
                                           ->  _rhs_characteristic
                                           ->  _point_chart
                                           ->  vacuum_chart
```

Four functions, the last two of them on the hot path of an adaptive ODE
solver. That is a real refactor of core numerical code rather than a local
edit.

An earlier attempt at the split was **reverted** rather than left in place.
It had a wrong filled-node count and still failed for want of the flag, and
half-applied changes to the solver that produces every recorded family are
not worth carrying. The curated suite confirms the revert: `368 passed,
1 skipped`.

## What is left, precisely

1. Thread `allow_below_inner_horizon` through `_point_chart` and
   `_rhs_characteristic`, defaulting to off so the solver's behaviour is
   unchanged for every caller that does not ask.
2. Cut the inward branch at `r_minus` and restart it from the inner jet.
3. Fill the nodes near `r_minus` by Taylor, as is already done near
   `r_plus`, so no branch integrates into its own starting point.

Steps 1 and 3 are mechanical. Step 2 is where the work is, and it is where a
bug would be expensive, because every recorded family in this programme was
produced by the path being changed.

## What this does not establish

That the fix works. That the third region can be marched. That anything in
it is a nest. None of that is claimed here, and the failure mode if the
refactor goes wrong is a silently corrupted recorded family rather than a
loud error.

That is the argument for doing it as its own piece of work, with its own
record, rather than folded into the end of a session.

## Scope

A design note. **It does not measure the actual universe** and it changes
no physics. One production function was generalised in a backward-compatible
way and nothing else was altered.

## Provenance

- `lab/src/recursive_horizons/nsc_rn_source.py`: `_horizon_jet` takes an
  optional radius, defaulting to `r_plus`
- the reverted split is described here rather than committed
- tests: `lab/tests/test_nsc_rn_third_region_reach.py` continues to pin that
  the defaults still refuse
