# The shape of the compressed layer

The mass-shell work measured how much the density peak grows against the
inner edge: `x1.20` at the horizon midpoint, saturating at `x1.58`. That is an
amplitude. It says nothing about what the compressed layer *looks like*.

This reads the shape, at the compression peak near `t = 1.7`, across all six
edge positions. Each profile is divided by its own value at the inner edge, so
the question asked is whether the shape moves with the boundary or only the
amplitude does.

## The layer narrows as the boundary deepens

| `r_inner` | peak sampled | edge density | half-width | cells to the horizon |
|---|---|---|---|---|
| `1.0400` | `1.4` | `3.144e-04` | `1.7121` | `8.009` |
| `0.8500` | `1.5` | `3.761e-04` | `1.3637` | `2.666` |
| `0.8000` | `1.7` | `4.202e-04` | `0.9346` | `1.270` |
| `0.7700` | `1.7` | `4.504e-04` | `0.7915` | `0.435` |
| `0.7600` | `1.7` | `4.592e-04` | `0.7558` | `0.157` |
| `0.7555` | `1.7` | `4.410e-04` | `0.8279` | `0.032` |

The density at the boundary rises by `46%` from the midpoint to `0.76`, and
the layer **narrows monotonically** over the same range, from `1.712` to
`0.756` in units of `r_m`, a factor of `2.27`. Amplitude up, width down: the
layer is genuinely compressing, not merely scaling.

### The peak is not at one time across the sweep

The compression peak **moves with the edge**: `1.4` at the midpoint, `1.5` at
`0.85`, `1.7` from `0.80` down. The first version of this reduction took a
single shared peak time and sampled every run at `1.7`, which put the two
shallowest edges *past* their peak. Their half-widths came out at `1.3554`
and `1.1125` instead of `1.7121` and `1.3637`.

The trend survives — it still narrows monotonically and by a larger factor
than reported — but the first two numbers in the table above were wrong and
are corrected here. The driver now takes a peak time per run and refuses a
run with no station near its own peak rather than falling back to a
neighbouring run's.

That is the fourth time in this line that a statistic read off the wrong
clock has produced a wrong number: the shell fit that ran through a phase
change, the curvature reading that took its mass peak from a different march,
the trajectory comparison that matched stations by file name, and now this.

For five of the six the density falls monotonically from the boundary with no
interface and no second maximum. That is a gradient, in the literal sense.

## And at the deepest boundary it stops being one

At `0.7555` the half-width turns back up to `0.8279` and the density at the
boundary falls from `4.592e-04` to `4.410e-04`. The same place where the
compression gain turned over.

More sharply: the **matter density maximum leaves the boundary**. At every
other edge the maximum sits exactly on node `0`. At `0.7555` it sits one cell
in, at `r = 0.7915`, standing **`22%` above** the boundary value.

That is a real, resolved feature and not a tie. The margin is stated because
the companion quantity is marginal: the **mass gradient** at that same node is
only `0.6%` above the boundary value, which is not a resolved detachment. The
record reports the two maxima separately and claims the detachment only for
the matter density.

So the picture has three stages, not one:

- from `1.04` to `0.76`: a monotone gradient against the boundary, tightening
  as the boundary deepens;
- at `0.76`: the tightest layer, still attached;
- at `0.7555`: the maximum detaches and becomes an interior peak at `r = 0.79`.

## What this does not establish

`0.7555` is **`0.032` cells** from the inner horizon at `0.7543`. On this grid
that is less than a tenth of a cell of separation. The detachment is a `22%`
margin at the *first* node, which is resolved, but its radius is resolved only
to one grid cell and it sits at the edge of what the grid can distinguish from
the horizon.

So the sequence above is a description of **this construction on this grid**,
not a statement that the compressed layer forms an interior shell near a
horizon. The stage boundaries would move under a different grid, and the
`8000`-point runs that could test that were only run to the cpu cap at
`t = 1.5`, before the peak.

## Scope

Construction checks on finite marched collars. **This does not measure the
actual universe.** A gradient zone that tightens under compression in one
numerical construction is not evidence that any real object has one.

Six edge positions, one release radius, one mass ratio, one source strength,
one grid resolution, one moment in time.

## Provenance

- reduction: `lab/scripts/derive_nsc_rn_layer_shape.py`
- record: `lab/results/development/nsc-rn-layer-shape-v3.json`; `v1` and `v2`
  are retained, `v2` having tested detachment against a threshold below `1`
  and so listed every run
- tests: `lab/tests/test_nsc_rn_layer_shape.py`

## What would settle it

The `8000`-point runs at the deep edges. They resolve the edge region by six
times and would put the detachment's radius, and the half-width sequence, on a
much finer grid. They reached the cpu cap at `t = 1.5` and need to be carried
past the peak at `1.7`, which is roughly another `1200 s` of CPU each.

That is now known to be affordable, which is why it is the obvious next step
rather than a resource question.