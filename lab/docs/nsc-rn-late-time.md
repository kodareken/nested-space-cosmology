# The shell does not empty

The question is whether the compressed shell decays to nothing or leaves a
remainder. It is worth asking because the two answers mean different things,
and only one of them is observable: a shell that empties has emptied, and a
shell that bottoms out and then grows again has a remainder that persists
and returns.

| | edge `0.77` | edge `0.76` |
|---|---|---|
| peak, at `t = 1.7` | `3.6180e-04` | `3.6729e-04` |
| trough, at `t = 4.9` | `8.9552e-06` | `8.8983e-06` |
| trough as a fraction of peak | **`0.0248`** | **`0.0242`** |
| final, at `t = 6.0` | `1.5597e-05` | `1.6298e-05` |
| final as a fraction of peak | `0.0431` | `0.0444` |
| turns over | yes | yes |
| rise is monotone | yes | yes |
| **ever reaches zero** | **no** | **no** |

The shell compresses, decays by more than a factor of thirty, reaches a
floor at about two and a half percent of its peak, **turns over, and grows
again** for the remainder of the run. It does not approach zero at any point,
and the two edges agree on the floor to three figures.

So the answer to the question is the one that leaves something: the decay
terminates on a remainder, and the remainder returns.

The floor is at `t = 4.9` in both runs, which is roughly three times the
rise time after the peak. Whether it has settled by `t = 6.0` is not
established; the rise is still monotone at the last recorded station and had
not begun to turn.

The horizon drifts with the same shape over the same interval, rising,
dipping and rising again. That coincidence is recorded and not explained.

## Scope

Construction checks on finite marched collars. **This does not measure the
actual universe.** A remainder that persists here is a property of this
construction on these grids. Nothing here says it is physical, and the runs
stop at `t = 6.0` with the amplitude still climbing, so the long-time state
is not known.

Two runs, two edge positions, one release radius, one mass ratio, one source
strength, one grid resolution.

## Provenance

- reduction: `lab/scripts/derive_nsc_rn_late_time.py`
- record: `lab/results/development/nsc-rn-late-time-v1.json`
- tests: `lab/tests/test_nsc_rn_late_time.py`
