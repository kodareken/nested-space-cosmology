# The 8000-point verdict: the compression is x2.05, not x1.58

The previous commit said the wide-grid results were converged, on the
strength of the inner mass stock agreeing between `869` and `8000` points to
about `1e-6`. That was true and it did not imply what I claimed.

Six `8000`-point marches at the six edge positions, `21` stations each, past
the compression peak at `t` near `1.7`.

| `r_inner` | gain at 869 | **gain at 8000** | width at 869 | **width at 8000** |
|---|---|---|---|---|
| `1.0400` | `1.1972` | **`1.2571`** | `1.7121` | **`1.6798`** |
| `0.8500` | `1.3905` | **`1.6039`** | `1.3637` | **`1.2539`** |
| `0.8000` | `1.4328` | **`1.8391`** | `0.9346` | **`0.8074`** |
| `0.7700` | `1.5758` | **`2.0481`** | `0.7915` | **`0.6364`** |
| `0.7600` | `1.5825` | **`2.0295`** | `0.7558` | **`0.5858`** |
| `0.7555` | `1.5849` | **`1.9952`** | `0.8279` | **`0.6172`** |

**The saturation is near `x2.0`, not `x1.58`.** The coarse grid under-resolved
it by five to thirty percent, most at the deep edges. The narrowing of the
layer is likewise stronger: `2.265` on the coarse grid, `2.867` at `8000`
points.

The turn-over also moves: it appears at `0.76` rather than at `0.7555`, with
the maximum `x2.0481` at `0.77`.

## Why the trajectory converged and the gain did not

The inner mass stock agreed between resolutions to `5e-09` at `t = 0.1`
rising to about `1.2e-06` at the compression peak. That agreement is real.

But the gain is a ratio of two **peak amplitudes of `dM/dr`**, and `dM/dr` is
a spatial derivative of the very quantity that converged. A finite difference
of a converged field is not converged to the same accuracy: it amplifies
whatever remains between grid points.

So the check that passed verified the trajectory and I generalised it to a
quantity derived from it. That was my error, and the numbers now say so
plainly: `x1.58` was a coarse-grid artifact.

## What survived refinement

**The three-sided shell.** Side counts are `[2, 2, 2, 2, 2, 3]` at both
resolutions. Five edges carry a boundary-attached layer with only two sides;
only the deepest, `0.7555`, has a genuine interior, a surface and an
exterior. That is the one structural finding in this session that does not
move when the ruler changes, and it is the direct answer to the observation
that a sphere has three sides where a two-sided quantity would have two.

**The monotonicity of the narrowing.** Still monotone from `1.04` down to
`0.76`, still turning back up at `0.7555`, at both resolutions.

## What has to be withdrawn

The claim that the wide-grid gains were converged. They were not. `x1.58` is
retained as a record of what the coarse grid produced and is superseded by
`x2.05`.

The convergence claim itself was narrower than stated: the trajectory
converged, and nothing downstream of a derivative of it was tested until now.

## Scope

Construction checks on finite marched collars. **This does not measure the
actual universe.**

Six edges, one release radius, one mass ratio, one source strength, two
resolutions. The `8000`-point runs reach `t = 2.0`, which is past the peak but
not far past the relaxation, so the relaxation fraction has not been compared
across resolutions and is not claimed to be.

## Provenance

- marches: `nsc-rn-inner-packet-n8000d-{exc080,exc077,exc076,exc07555}` and
  `nsc-rn-inner-packet-n8000-{def1040,exc085}`, `21` stations each
- reductions: `derive_nsc_rn_inner_packet.py`, `derive_nsc_rn_layer_shape.py`
- records: `nsc-rn-layer-shape-n8000-v1.json`, and the six
  `nsc-rn-inner-packet-shell-n8000*-v1.json`

## What would settle it

A third resolution. Two points bracket nothing, and the gain moved by up to
thirty percent between the two already tried, so there is no reason to
believe either is converged. `16000` points would make the gap about four
cells and cost about `230 s` of CPU per station, which is now affordable but
is roughly a day's marching for a full sweep.