# The three-sided shell is withdrawn: its interior is sub-cell

A claim is withdrawn here rather than qualified, because the measurement
does not support it at any resolution tried.

The claim was that the compressed layer gains a **third side** at the
deepest edge — an interior, a surface and an exterior — against two sides
everywhere else. It was offered as a direct answer to the observation that a
sphere has three sides where a two-sided quantity has two.

## What the interior actually is

The third side appears only at `r_inner = 0.7555`. Its interior is the
material between the boundary and `r_minus`:

| `r_inner` | to `r_minus` | cells at 869 | cells at 8000 |
|---|---|---|---|
| `1.0400` | `0.28566` | `7.94` | `73.1` |
| `0.8500` | `0.09566` | `2.66` | `24.5` |
| `0.8000` | `0.04566` | `1.27` | `11.7` |
| `0.7700` | `0.01566` | `0.44` | `4.01` |
| `0.7600` | `0.00566` | `0.16` | `1.45` |
| **`0.7555`** | **`0.00116`** | **`0.032`** | **`0.296`** |

**The interior is 0.032 cells wide at `869` points and 0.296 cells at `8000`.**
It is smaller than one cell at every resolution this line has run.

The matter maximum that was read as the "surface" sits `0.00391` above the
edge — exactly **one cell** at `8000` points.

## Why it is most likely a boundary-layer effect

The third side appears only where the excision is nearly on `r_minus`. And
there the outgoing speed goes to zero: it changes sign at `r_minus`, so within
`0.0012` of it the chart is nearly null. A sourced density measured on a
nearly-null chart, one cell out, is exactly the configuration in which a
boundary layer produces a dip at the boundary node and a maximum just beyond
it.

That is what the numbers look like. A one-cell dip plus a one-cell maximum,
appearing only where the boundary condition is about to reverse, is not the
same thing as a resolved interior.

## Withdrawn

**The three-sided shell.** The side counts `[2, 2, 2, 2, 2, 3]` stand as a
statement about the recorded arrays, and nothing more. They do not establish
a third side, an interior, or a sphere with an inside.

The `22%` overshoot and the `1.19` overshoot at `8000` points are likewise
statements about two adjacent nodes at a boundary that is `0.3` cells from
where the chart turns over.

## What survives

Everything radial and thick. The layer narrows from `1.680` to `0.586`,
about `2.87` times, across a region tens of cells wide at both resolutions.
The amplitude rises by `1.57` over the same range. Those are measured over
many cells and are not boundary artefacts.

The wave reading survives with that correction: the layer is **louder and
shorter**, and it is a one-dimensional radial layer. It is not a sphere with
an interior, and this commit is why that can be said rather than assumed.

## What would settle it

A grid fine enough that the region between the boundary and `r_minus` spans
several cells, which at this geometry means roughly `40000` points, or an
excision placed far enough from `r_minus` that a three-sided structure could
be measured away from the boundary condition. Neither is a small run.

## Scope

A withdrawal. **It does not measure the actual universe.** It removes a
claim that was made on this line and does not survive the check applied to
it.

## Provenance

- the interior width is the distance from the edge to `r_minus`, from
  `nsc-rn-geometry-map-v6.json` and the edge positions in
  `nsc-rn-layer-shape-n8000-v1.json`
- superseded claims: `lab/docs/nsc-rn-layer-shape.md` and the three-sided
  tests in `lab/tests/test_nsc_rn_layer_shape.py`, both updated
