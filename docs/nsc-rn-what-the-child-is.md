# What the child is, read off the construction

This has been asked three times and left open. It turns out the construction
declares the answer in plain text, and it is not the answer the numbers
suggested.

From `nsc_rn_reference`:

```
"excision": "pure outflow on an open interval between r_minus and r_plus",
"parent":   "large-r exterior",
"child":    "near the outer horizon, different areal interval",
```

So:

| region | interval |
|---|---|
| **trapped / excision** | `(0.7543, 1.3257)` |
| **child** | `[1.3614, 1.6113]` |
| **parent** | `(1.6113, 32]` |

## The child is not the interior

The child begins one gap-sixteenth *above* the outer horizon, not below it.
It sits in the parent, adjacent to the horizon, separated from it by a band
of `0.0357`.

The hierarchy the code declares is therefore

```
parent  ⊃  child  ⊃  [ horizon ]  ⊃  trapped region
```

and **the horizon is the membrane between the child and the trapped region**,
not between the parent and the child.

## What that means for this stretch's measurements

Every compression measurement in this line — the `x2.05`, the `2.87` layer
narrowing, the three-sided shell, the `2.5%` floor, the `94%` retention —
was taken in the **trapped region**, which the construction calls the
*excision*, not the child.

That is not a complaint about the work. It is a correction to what the work
was about. The trapped region is where matter that cannot escape ends up,
and that turned out to be a real and well-behaved object. It is simply not
the region the construction names "child".

## What this does not settle

**Whether the child is or contains a universe.** It is an areal interval on
one chart with `N = 1` and a known characteristic structure. That is a
description, not an occupancy.

**What separates the child from the parent.** The construction draws the
cut at `r_plus + gap/2` and nothing else. There is no condition, no
membrane and no measurement there — the earlier check found the metric is
smooth across it to `3.1e-04` on a field that varies by `0.5`. So the
parent/child boundary is a **label with no physical signature in this
geometry**.

That is the load-bearing gap now, and it is sharp: the construction nests
parent and child, and the only physical interface it contains is at
`r_plus`, which is *below* the child, not at its boundary.

## Scope

A reading of the construction's own declarations plus the recorded
measurements. **It does not measure the actual universe.**

## Provenance

- declarations: `lab/src/recursive_horizons/nsc_rn_reference.py`
- intervals: `nsc_rn_geometry_map-v6.json`
- smoothness across the parent/child cut: measured in the nesting-chain
  document
